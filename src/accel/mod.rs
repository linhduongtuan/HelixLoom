//! Candle on the GPUs of the workstation: Apple's GPU (Metal) on a Mac, every CUDA card on a Linux or Windows machine, shared
//! out by what each can do (`pool`). The work that is worth sending there: distances between many sequences and the
//! resampled columns of a bootstrap (`distance`), and the likelihood of a tree over many alignment columns (`seq::likelihood::gpu`),
//! which maximum-likelihood search, model fitting and Bayesian sampling all call. The plain Rust code stays the reference:
//! the tests compare the two, and a person can compare them on their own cards with `helixloom gpu-check`.
//!
//! Not for a browser (candle does not build for it) and not without the `gpu` feature.

pub mod check;
pub mod distance;
pub mod pool;
pub mod torch;

pub use pool::{mode, Mode, Pool};

/// `candle` or `torch` (or `auto`, the start: PyTorch when `HELIXLOOM_TORCH_PYTHON` names a Python, else candle): which library
/// does the work that goes to the GPUs. Both give the same answers (the tests and `helixloom gpu-check` compare them with the
/// processor's); candle is built into the program, PyTorch runs in Python workers (`torch`).
pub const BACKEND_VARIABLE: &str = "HELIXLOOM_ACCEL_BACKEND";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Backend {
    Candle,
    Torch,
}

/// The backend that the environment asks for.
pub fn backend() -> Backend {
    parse_backend(
        &std::env::var(BACKEND_VARIABLE).unwrap_or_default(),
        torch::configured_python().is_some(),
    )
}

fn parse_backend(text: &str, python_named: bool) -> Backend {
    match text.trim().to_ascii_lowercase().as_str() {
        "torch" | "pytorch" => Backend::Torch,
        "candle" => Backend::Candle,
        _ if python_named => Backend::Torch,
        _ => Backend::Candle,
    }
}

/// The joins of neighbour-joining over `dist` when the GPUs are to make them (PyTorch, and enough sequences to be worth it): see
/// `seq::tree::neighbor_joining`. `None` where the processor does it, and when the GPUs fail (said on the standard error).
pub fn neighbor_joins(dist: &[Vec<f64>]) -> Option<(Vec<crate::seq::tree::Join>, Vec<Vec<f64>>)> {
    if backend() != Backend::Torch || dist.len() <= 3 {
        return None;
    }
    let n = dist.len() as f64;
    if !worth_it(n * n * n, threshold::join_cells()) {
        return None;
    }
    let pool = torch::TorchPool::global()?;
    match torch::nj::joins(pool, dist) {
        Ok(joins) => {
            // once, in the job's log: nothing else says which of the two made the tree
            static SAID: std::sync::Once = std::sync::Once::new();
            SAID.call_once(|| {
                eprintln!(
                    "HelixLoom: neighbour-joining of {} sequences made on {}",
                    dist.len(),
                    pool.workers()[pool.fastest()].name
                )
            });
            Some(joins)
        }
        Err(why) => {
            note_failure(
                "the GPUs failed on neighbour-joining",
                &format!("{why}; the processor does it"),
            );
            None
        }
    }
}

/// The alignments of `bs` with the sequence `a` when the PyTorch workers are to make them (see `torch::align`): **only where the person
/// asked** (`HELIXLOOM_ACCEL=force`, or `HELIXLOOM_ACCEL_ALIGN_CELLS` set to a size, in cells of all the tables together, from which they
/// are), because the processor is faster everywhere it has been measured. `None` where the processor does it, and when the workers fail (said
/// on the standard error, once). Each `b` and `a` have a base and no pair is too large for the processor's `align`.
pub fn align_pairs(
    a: &[u8],
    bs: &[&[u8]],
    sc: &crate::seq::align::Scoring,
    free_ends: bool,
) -> Option<Vec<crate::seq::align::Alignment>> {
    if backend() != Backend::Torch || bs.is_empty() {
        return None;
    }
    let cells: f64 = bs
        .iter()
        .map(|b| (a.len() + 1) as f64 * (b.len() + 1) as f64)
        .sum();
    if !worth_it(cells, threshold::align_cells()) {
        return None;
    }
    let pool = torch::TorchPool::global()?;
    match torch::align::alignments(pool, a, bs, sc, free_ends) {
        Ok(done) => {
            // once, in the job's log: nothing else says which of the two made the alignments
            static SAID: std::sync::Once = std::sync::Once::new();
            SAID.call_once(|| {
                eprintln!(
                    "HelixLoom: alignments of {} sequences made on {} PyTorch worker(s), as asked",
                    bs.len(),
                    pool.len()
                )
            });
            Some(done)
        }
        Err(why) => {
            note_failure(
                "the GPUs failed on alignments",
                &format!("{why}; the processor does them"),
            );
            None
        }
    }
}

/// Whether the backend in use has a device to compute on (which, for PyTorch, starts its workers the first time).
pub fn available() -> bool {
    match backend() {
        Backend::Candle => Pool::global().is_some(),
        Backend::Torch => torch::may_have_devices() && torch::TorchPool::global().is_some(),
    }
}

/// The devices of whichever backend is in use, as `helixloom gpu-check` and its tests try them: all together, or one at a time.
pub enum Target {
    Candle(Pool),
    Torch(std::sync::Arc<torch::TorchPool>),
}

impl Target {
    /// The devices of the backend in use, found afresh (not the pool that jobs of this process share): `Ok(None)` where it has none,
    /// `Err` where it cannot start.
    pub fn discover() -> Result<Option<Target>, String> {
        match backend() {
            Backend::Candle => Ok(Pool::discover().map(Target::Candle)),
            Backend::Torch => Ok(torch::TorchPool::from_environment()?
                .map(|pool| Target::Torch(std::sync::Arc::new(pool)))),
        }
    }

    pub fn len(&self) -> usize {
        match self {
            Target::Candle(p) => p.len(),
            Target::Torch(p) => p.len(),
        }
    }

    /// Device `k` alone, with all the work.
    pub fn only(&self, k: usize) -> Target {
        match self {
            Target::Candle(p) => {
                let m = &p.members()[k];
                Target::Candle(Pool::new(vec![pool::Member::new(
                    m.device.clone(),
                    m.name.clone(),
                    1.0,
                )]))
            }
            Target::Torch(p) => Target::Torch(std::sync::Arc::new(p.only(k))),
        }
    }

    /// For each device, its name and the shares of the likelihood and of the distances it is given (the same for candle).
    pub fn cards(&self) -> Vec<(String, f64, f64)> {
        match self {
            Target::Candle(p) => p
                .members()
                .iter()
                .map(|m| (m.name.clone(), m.weight, m.weight))
                .collect(),
            Target::Torch(p) => p.cards(),
        }
    }

    /// What the backend is, for the report.
    pub fn backend_name(&self) -> String {
        match self {
            Target::Candle(_) => "candle".into(),
            Target::Torch(p) => p.versions.clone(),
        }
    }

    /// The joins of neighbour-joining over `dist` on the devices; `None` for a backend that cannot (candle).
    #[allow(clippy::type_complexity)]
    pub fn joins(
        &self,
        dist: &[Vec<f64>],
    ) -> Option<Result<(Vec<crate::seq::tree::Join>, Vec<Vec<f64>>), String>> {
        match self {
            Target::Candle(_) => None,
            Target::Torch(p) => Some(torch::nj::joins(p, dist)),
        }
    }

    /// The alignments of `bs` with `a` on the devices (opt-in work, see `align_pairs`); `None` for a backend that cannot (candle).
    pub fn alignments(
        &self,
        a: &[u8],
        bs: &[&[u8]],
        sc: &crate::seq::align::Scoring,
        free_ends: bool,
    ) -> Option<Result<Vec<crate::seq::align::Alignment>, String>> {
        match self {
            Target::Candle(_) => None,
            Target::Torch(p) => Some(torch::align::alignments(p, a, bs, sc, free_ends)),
        }
    }

    /// The rows of an alignment on all the devices.
    pub fn distances(&self, rows: &[Vec<u8>]) -> Result<Distances, String> {
        match self {
            Target::Candle(p) => distance::Engine::new(p, rows).map(Distances::Candle),
            Target::Torch(p) => torch::distance::Engine::new(p, rows).map(Distances::Torch),
        }
    }
}

/// The columns of an alignment on the devices of whichever backend is in use, for the counts that make distances.
pub enum Distances {
    Candle(distance::Engine),
    Torch(torch::distance::Engine),
}

impl Distances {
    /// Put the rows on the devices: `Ok(None)` where there is no device, `Err` with the reason where there is one and it cannot take them.
    pub fn open(rows: &[Vec<u8>]) -> Result<Option<Distances>, String> {
        match backend() {
            Backend::Candle => Pool::global()
                .map(|pool| distance::Engine::new(pool, rows).map(Distances::Candle))
                .transpose(),
            Backend::Torch => torch::TorchPool::global()
                .map(|pool| torch::distance::Engine::new(pool, rows).map(Distances::Torch))
                .transpose(),
        }
    }

    /// See `distance::Engine::counts`.
    pub fn counts(&self, weights: Option<&[u32]>) -> Result<distance::Counts, String> {
        match self {
            Distances::Candle(e) => e.counts(weights),
            Distances::Torch(e) => e.counts(weights),
        }
    }

    /// Which devices, and what share of the rows each does, in words.
    pub fn describe(&self) -> String {
        match self {
            Distances::Candle(e) => e.describe(),
            Distances::Torch(e) => e.describe(),
        }
    }
}

/// How big a piece of work must be before the GPUs are worth it: the numbers of cells (rows times columns of what the cards would
/// compute over) from which they beat the processor. **They come from one machine, the Apple GPU of a Mac** (`helixloom gpu-check
/// --bench --quick`, November 2026: distances won 5.5 times at 64 sequences by 20 000 columns and lost at 64 by 2 000; a likelihood
/// lost up to about 6 million cells and won 1.2 times at 6.4 million), so on an NVIDIA card they are a place to start from:
/// `gpu-check --bench` measures where the cards of a machine overtake the processor, and the environment variables below
/// override them.
pub mod threshold {
    /// Patterns times tips of a likelihood.
    pub const LIKELIHOOD_CELLS: f64 = 1.0e7;
    pub const LIKELIHOOD_VARIABLE: &str = "HELIXLOOM_ACCEL_LIKELIHOOD_CELLS";
    /// Sequences squared times columns of a set of distances, for all the resamples of a bootstrap together.
    pub const DISTANCE_CELLS: f64 = 1.0e8;
    pub const DISTANCE_VARIABLE: &str = "HELIXLOOM_ACCEL_DISTANCE_CELLS";

    /// With PyTorch the distances need more to be worth it: its workers start in every job process that uses them (importing PyTorch and
    /// opening the cards take seconds that candle, built in, does not), so a set of distances that the processor does in a second is
    /// not worth sending. Roughly: ten seconds of one core's time (a guess from a Mac; `gpu-check --bench` on the workstation decides).
    pub const DISTANCE_CELLS_TORCH: f64 = 5.0e9;

    /// Sequences cubed, from which neighbour-joining goes to the GPUs (PyTorch only): the processor takes about n cubed steps, which is
    /// twenty seconds or so of one core at this size (about 2,700 sequences; on a Mac, 600 sequences took 0.15 s and the time grows as
    /// the cube). A guess, as the others are: the cards' side has been tried only on the processor's PyTorch, which loses to the Rust
    /// code at these sizes, so `gpu-check --bench` on the cards is what decides.
    pub const JOIN_CELLS: f64 = 2.0e10;
    pub const JOIN_VARIABLE: &str = "HELIXLOOM_ACCEL_JOIN_CELLS";

    /// Cells of all the tables of a set of alignments, from which they go to the PyTorch workers: **never, to start with** (infinite).
    /// The processor does them 100 times or more faster than a Mac's GPU or a CPU worker did; a card of another kind may win, which
    /// `gpu-check --bench` says, and then this variable is how to use it.
    pub const ALIGN_CELLS: f64 = f64::INFINITY;
    pub const ALIGN_VARIABLE: &str = "HELIXLOOM_ACCEL_ALIGN_CELLS";

    /// The size from which a likelihood goes to the GPUs: the person's variable, else the start for the backend in use.
    pub fn likelihood_cells() -> f64 {
        from_environment(LIKELIHOOD_VARIABLE, LIKELIHOOD_CELLS)
    }

    /// The same for neighbour-joining.
    pub fn join_cells() -> f64 {
        from_environment(JOIN_VARIABLE, JOIN_CELLS)
    }

    /// The same for alignments.
    pub fn align_cells() -> f64 {
        from_environment(ALIGN_VARIABLE, ALIGN_CELLS)
    }

    /// The same for distances.
    pub fn distance_cells() -> f64 {
        let default = match super::backend() {
            super::Backend::Candle => DISTANCE_CELLS,
            super::Backend::Torch => DISTANCE_CELLS_TORCH,
        };
        from_environment(DISTANCE_VARIABLE, default)
    }

    /// The threshold in `variable` if it holds a number, else `default`.
    pub fn from_environment(variable: &str, default: f64) -> f64 {
        std::env::var(variable)
            .ok()
            .and_then(|v| v.trim().parse::<f64>().ok())
            .filter(|v| v.is_finite() && *v >= 0.0)
            .unwrap_or(default)
    }
}

/// Whether a tree job of this size should be given GPUs, by its size alone (what the server asks: it must not touch a GPU itself, and
/// a card is found only inside the process that runs the job). A likelihood (maximum likelihood, Bayesian) costs patterns times tips,
/// of which the number of columns is the most there can be; distances (and the bootstrap of them) cost sequences squared times
/// columns for each resample; neighbour-joining, with PyTorch, costs sequences cubed.
pub fn tree_job_is_big(
    tips: usize,
    columns: usize,
    likelihood: bool,
    joins: bool,
    replicates: usize,
) -> bool {
    match mode() {
        Mode::Off => false,
        Mode::Force => true,
        Mode::Auto if likelihood => tips as f64 * columns as f64 >= threshold::likelihood_cells(),
        Mode::Auto => {
            // neighbour-joining itself (only PyTorch does it): sequences cubed
            let joining = joins
                && backend() == Backend::Torch
                && (tips as f64).powi(3) >= threshold::join_cells();
            distances_are_big(tips, columns, replicates) || joining
        }
    }
}

/// Whether the distances of a tree job (`replicates` sets of them, for a bootstrap) are big enough to be given GPUs.
pub fn distances_are_big(tips: usize, columns: usize, replicates: usize) -> bool {
    match mode() {
        Mode::Off => false,
        Mode::Force => true,
        Mode::Auto => {
            replicates.max(1) as f64 * tips as f64 * tips as f64 * columns as f64
                >= threshold::distance_cells()
        }
    }
}

/// Whether work of `cells` cells should go to the GPUs, given the mode and the `threshold`: never when they are off or there
/// are none (of the backend in use), always when forced, else when it is big enough.
pub fn worth_it(cells: f64, threshold: f64) -> bool {
    match mode() {
        Mode::Off => false,
        Mode::Force => available(),
        Mode::Auto => cells >= threshold && available(),
    }
}

/// Say once, on the standard error (which a job's log keeps), that something about the GPUs failed (`what`: the same words every time
/// it happens, which is how the repeats are told) and why. The processor does the work then, which is slower, and without this line
/// nobody would know. Each `what` is said once in a process.
pub fn note_failure(what: &str, why: &str) {
    use std::sync::{Mutex, OnceLock, PoisonError};
    static SAID: OnceLock<Mutex<Vec<String>>> = OnceLock::new();
    let mut said = SAID
        .get_or_init(Mutex::default)
        .lock()
        .unwrap_or_else(PoisonError::into_inner);
    if !said.iter().any(|w| w == what) {
        said.push(what.to_string());
        eprintln!("HelixLoom: {what}: {why}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_backend_is_what_the_person_names_else_pytorch_only_where_a_python_is_given() {
        assert_eq!(parse_backend("", false), Backend::Candle);
        assert_eq!(parse_backend("auto", false), Backend::Candle);
        assert_eq!(parse_backend("", true), Backend::Torch);
        assert_eq!(parse_backend("auto", true), Backend::Torch);
        assert_eq!(parse_backend(" Torch ", false), Backend::Torch);
        assert_eq!(parse_backend("pytorch", false), Backend::Torch);
        // a name for candle beats the python being given
        assert_eq!(parse_backend("CANDLE", true), Backend::Candle);
        assert_eq!(parse_backend("tensorflow", false), Backend::Candle);
    }
}
