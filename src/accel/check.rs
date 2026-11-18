//! `helixloom gpu-check`: finds the GPUs this program can use, tries each of them and then all together on the two kinds of work
//! that go there (distances and the likelihood of a tree) against the processor's own answers, and, with `--bench`, measures how big
//! the work must be before the cards are faster than the processor. It is how a person finds out, on their own cards, what could
//! not be tried when this was written: whether a card computes the same numbers, how the work is cut among cards that are not
//! alike, and where the crossover lies.

use super::distance::Counts;
use super::pool::{mode, MODE_VARIABLE, WEIGHTS_VARIABLE};
use super::threshold::{DISTANCE_VARIABLE, LIKELIHOOD_VARIABLE};
use super::torch::{DEVICES_VARIABLE, PYTHON_VARIABLE};
use super::{Target, BACKEND_VARIABLE};
use crate::seq::align::Scoring;
use crate::seq::likelihood::check as likelihood;
use crate::seq::tree::{aligned_distance, cpu_joins, join_tree, splits, DistanceModel, Rng};
use std::time::Instant;

pub const USAGE: &str = "\
helixloom gpu-check [--bench] [--quick]

Finds the GPUs that the backend in use can reach, tries each one and then all together against the processor's answers, and says how
the work is shared among them. The backend is candle (built into the program: CUDA cards on a Linux or Windows machine built with
`--features cuda`, the Apple GPU on a Mac) or, where HELIXLOOM_TORCH_PYTHON names a Python that has PyTorch (or
HELIXLOOM_ACCEL_BACKEND=torch), PyTorch in worker processes, one for each card (HELIXLOOM_TORCH_DEVICES=cpu tries it without a GPU).

  --bench   also time the processor against the GPUs on alignments of growing size, and say from what size the GPUs win
  --quick   with --bench, smaller sizes (a minute or two instead of several)

Exit status: 0 when every check passed, 1 when one failed, 2 when there is no GPU to try.
";

/// A likelihood or a distance that differs from the processor's by more than this, as a share, is a failure: 32-bit arithmetic agrees
/// to about the seventh digit, so a card that is out in the sixth is not computing the same thing.
const LIKELIHOOD_TOLERANCE: f64 = 1e-5;

struct Report {
    failures: usize,
}

impl Report {
    fn say(&mut self, ok: bool, text: String) {
        println!("  [{}] {text}", if ok { "ok" } else { "FAIL" });
        if !ok {
            self.failures += 1;
        }
    }
}

/// `n` sequences of `columns` columns from a common ancestor, with gaps and ambiguity codes here and there.
fn alignment(n: usize, columns: usize, seed: u64) -> Vec<Vec<u8>> {
    let mut rng = Rng::new(seed);
    let letters = b"ACGT";
    let ancestor: Vec<u8> = (0..columns).map(|_| letters[rng.below(4)]).collect();
    (0..n)
        .map(|_| {
            ancestor
                .iter()
                .map(|&c| {
                    let x = if rng.below(100) < 15 {
                        letters[rng.below(4)]
                    } else {
                        c
                    };
                    match rng.below(40) {
                        0 => b'-',
                        1 => b'N',
                        2 => x.to_ascii_lowercase(),
                        _ => x,
                    }
                })
                .collect()
        })
        .collect()
}

/// The distances of every pair on the processor, under one model (`None` where a pair cannot be compared).
fn cpu_distances(rows: &[Vec<u8>], model: DistanceModel) -> Vec<Vec<Option<f64>>> {
    let n = rows.len();
    let mut d = vec![vec![Some(0.0); n]; n];
    for i in 0..n {
        for j in i + 1..n {
            let v = aligned_distance(&rows[i], &rows[j], model);
            (d[i][j], d[j][i]) = (v, v);
        }
    }
    d
}

/// Whether `counts` give the processor's distances for `rows` under `model`, to the last digit (and the same pairs that cannot be
/// compared).
fn distances_equal(counts: &Counts, rows: &[Vec<u8>], model: DistanceModel) -> bool {
    let want = cpu_distances(rows, model);
    match counts.distances(model) {
        Ok(got) => (0..rows.len()).all(|i| (0..rows.len()).all(|j| want[i][j] == Some(got[i][j]))),
        Err(_) => false,
    }
}

/// Whether `pool` gives the processor's distances to the last digit, for the three models; and how long the cards took.
fn check_distances(pool: &Target) -> Result<(bool, f64), String> {
    let rows = alignment(40, 3000, 5);
    let start = Instant::now();
    let counts = pool.distances(&rows)?.counts(None)?;
    let seconds = start.elapsed().as_secs_f64();
    let same = DistanceModel::ALL
        .into_iter()
        .all(|model| distances_equal(&counts, &rows, model));
    Ok((same, seconds))
}

/// Whether the columns taken as often as a bootstrap resample takes them (a weight for each, on the cards) give the distances of the
/// resampled alignment on the processor. The weighted products are another path through the cards than the plain ones.
fn check_bootstrap(pool: &Target) -> Result<bool, String> {
    let rows = alignment(30, 2000, 7);
    let width = rows[0].len();
    let mut rng = Rng::new(11);
    let picks: Vec<usize> = (0..width).map(|_| rng.below(width)).collect();
    let mut weights = vec![0u32; width];
    picks.iter().for_each(|&p| weights[p] += 1);
    let resampled: Vec<Vec<u8>> = rows
        .iter()
        .map(|r| picks.iter().map(|&c| r[c]).collect())
        .collect();
    let counts = pool.distances(&rows)?.counts(Some(&weights))?;
    Ok(DistanceModel::ALL
        .into_iter()
        .all(|model| distances_equal(&counts, &resampled, model)))
}

/// Several threads at once on the same cards, each getting the processor's distances. candle's Metal backend gave wrong answers when
/// they did not take turns; the program takes turns on a card (`Member::exclusive`) whatever the card, so this checks that taking
/// turns works on these cards, not whether CUDA would have needed it.
fn check_threads(pool: &Target) -> Result<bool, String> {
    let rows = alignment(24, 1200, 9);
    std::thread::scope(|s| {
        let jobs: Vec<_> = (0..6)
            .map(|_| {
                s.spawn(|| -> Result<bool, String> {
                    let mut same = true;
                    for _ in 0..3 {
                        let counts = pool.distances(&rows)?.counts(None)?;
                        same &= distances_equal(&counts, &rows, DistanceModel::P);
                    }
                    Ok(same)
                })
            })
            .collect();
        jobs.into_iter().try_fold(true, |all, job| {
            let one = job.join().map_err(|_| "a thread panicked".to_string())??;
            Ok::<bool, String>(all && one)
        })
    })
}

/// Distances between `n` random points in the plane: a matrix with no ties, which two neighbour-joinings must agree on.
fn plane_distances(n: usize, seed: u64) -> Vec<Vec<f64>> {
    let mut rng = Rng::new(seed);
    let points: Vec<(f64, f64)> = (0..n)
        .map(|_| (rng.unit() * 10.0, rng.unit() * 10.0))
        .collect();
    points
        .iter()
        .map(|a| {
            points
                .iter()
                .map(|b| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt())
                .collect()
        })
        .collect()
}

/// Whether the devices make the joins of neighbour-joining as the processor does, over 120 sequences: every join but the last (which
/// ties for four clusters, whatever the distances) the same, and the same tree. `None` where the backend does not do it.
fn check_joins(pool: &Target) -> Option<Result<bool, String>> {
    let dist = plane_distances(120, 3);
    let got = pool.joins(&dist)?;
    Some(got.map(|(joins, last)| {
        let (want, want_last) = cpu_joins(&dist);
        let names: Vec<String> = (0..dist.len()).map(|k| format!("s{k}")).collect();
        let strict = want.len().saturating_sub(1);
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-9 * b.abs().max(1.0);
        joins.len() == want.len()
            && joins
                .iter()
                .zip(&want)
                .take(strict)
                .all(|(g, w)| (g.i, g.j) == (w.i, w.j) && close(g.li, w.li) && close(g.lj, w.lj))
            && splits(&join_tree(&names, &joins, &last), &names)
                == splits(&join_tree(&names, &want, &want_last), &names)
    }))
}

/// Everything that is tried on one pool, with its name for the report.
fn check_pool(report: &mut Report, name: &str, pool: &Target) {
    println!("{name}:");
    match check_distances(pool) {
        Ok((same, seconds)) => report.say(
            same,
            format!("distances between 40 sequences of 3000 columns equal the processor's to the last digit ({:.0} ms)", 1000.0 * seconds),
        ),
        Err(why) => report.say(false, format!("distances: {why}")),
    }
    match check_bootstrap(pool) {
        Ok(same) => report.say(
            same,
            "distances of a bootstrap resample (the columns weighted) equal the processor's to the last digit".into(),
        ),
        Err(why) => report.say(false, format!("bootstrap distances: {why}")),
    }
    if let Some(result) = check_joins(pool) {
        match result {
            Ok(same) => report.say(
                same,
                "neighbour-joining over 120 sequences makes the processor's joins and tree".into(),
            ),
            // a device with no 64-bit floats (a Mac's GPU) cannot do it, and the processor does it: not a fault
            Err(why) if why.contains("64-bit") => {
                println!("  [--] neighbour-joining: {why}; the processor does it")
            }
            Err(why) => report.say(false, format!("neighbour-joining: {why}")),
        }
    }
    match likelihood::agree(pool, 24, 1500, 1) {
        Ok(a) => report.say(
            a.worst_relative <= LIKELIHOOD_TOLERANCE,
            format!(
                "the likelihood of a 24-tip tree over {} patterns, under {} models, from one thread and from par_map threads, differs from the processor's by at most {:.1e} of its size ({:.1e} log-likelihood units; allowed {:.0e} of its size)",
                a.patterns, a.models, a.worst_relative, a.worst_absolute, LIKELIHOOD_TOLERANCE
            ),
        ),
        Err(why) => report.say(false, format!("likelihood: {why}")),
    }
    match check_threads(pool) {
        Ok(same) => report.say(
            same,
            "six threads at once on these cards all get the processor's distances".into(),
        ),
        Err(why) => report.say(false, format!("threads: {why}")),
    }
}

/// The least work (cells) from which the cards won every size tried: the smallest `cells` such that every row of `rows` with that
/// many cells or more is a win (`true`). The time of the cards depends on the shape of the work as well as its size, so one
/// win among losses, by the cell count alone, says little; and a size at which they lost bars every smaller one from the answer.
fn crossover(mut rows: Vec<(f64, bool)>) -> Option<f64> {
    // equal sizes with different results: the loss comes last in the walk from the largest, and ends it
    rows.sort_by(|a, b| a.0.total_cmp(&b.0).then(b.1.cmp(&a.1)));
    let mut least = None;
    for &(cells, won) in rows.iter().rev() {
        if !won {
            break;
        }
        least = Some(cells);
    }
    least
}

/// The time of the processor and of `pool` for likelihoods of growing size, and the least work (patterns times tips) from which the
/// cards were at least half as fast again every time.
fn bench_likelihood(pool: &Target, quick: bool) -> Option<f64> {
    println!("Likelihood of one tree (one evaluation):");
    println!(
        "  {:>5} {:>9} {:>11} {:>11} {:>8}",
        "tips", "patterns", "processor", "GPUs", "faster"
    );
    let (tips, columns): (&[usize], &[usize]) = if quick {
        (&[16, 64], &[1_000, 20_000, 100_000])
    } else {
        (&[16, 64, 256], &[1_000, 10_000, 100_000, 400_000])
    };
    let mut tried = Vec::new();
    for &t in tips {
        for &c in columns {
            match likelihood::time(pool, t, c) {
                Ok(m) => {
                    let ratio = m.cpu_seconds / m.cards_seconds.max(1e-9);
                    println!(
                        "  {:>5} {:>9} {:>8.1} ms {:>8.1} ms {:>7.1}x",
                        m.tips,
                        m.patterns,
                        1000.0 * m.cpu_seconds,
                        1000.0 * m.cards_seconds,
                        ratio
                    );
                    tried.push(((m.patterns * m.tips) as f64, ratio >= 1.5));
                }
                Err(why) => println!("  {t:>5} {c:>9}  {why}"),
            }
        }
    }
    crossover(tried)
}

/// The same for the distances of every pair (the columns are sent to the cards as part of the time).
fn bench_distances(pool: &Target, quick: bool) -> Option<f64> {
    println!("Distances between all pairs (the processor on one core; the GPUs including sending the columns):");
    println!(
        "  {:>5} {:>8} {:>11} {:>11} {:>8}",
        "seqs", "columns", "processor", "GPUs", "faster"
    );
    let (seqs, columns): (&[usize], &[usize]) = if quick {
        (&[64, 256], &[2_000, 20_000])
    } else {
        (&[64, 256, 1024], &[2_000, 20_000])
    };
    let mut tried = Vec::new();
    for &n in seqs {
        for &l in columns {
            let rows = alignment(n, l, 3);
            let start = Instant::now();
            let _ = cpu_distances(&rows, DistanceModel::JukesCantor);
            let cpu = start.elapsed().as_secs_f64();
            // once to warm the cards up (a card's first run of a size compiles its kernels), then the timed one
            let _ = pool.distances(&rows).and_then(|e| e.counts(None));
            let start = Instant::now();
            let cards = pool.distances(&rows).and_then(|e| e.counts(None));
            let gpu = start.elapsed().as_secs_f64();
            match cards {
                Ok(_) => {
                    let ratio = cpu / gpu.max(1e-9);
                    println!(
                        "  {n:>5} {l:>8} {:>8.1} ms {:>8.1} ms {ratio:>7.1}x",
                        1000.0 * cpu,
                        1000.0 * gpu
                    );
                    tried.push(((n * n * l) as f64, ratio >= 1.5));
                }
                Err(why) => println!("  {n:>5} {l:>8}  {why}"),
            }
        }
    }
    crossover(tried)
}

/// What to look at when no device can be used, for the backend asked for.
fn no_gpu(why: Option<&str>) {
    println!("  No GPU can be used by this program.");
    if let Some(why) = why {
        println!("  - {why}");
    }
    if super::backend() == super::Backend::Torch {
        println!("  - PyTorch is used through a Python that has it installed: `{PYTHON_VARIABLE}=/path/to/python` (try `/path/to/python -c 'import torch; print(torch.cuda.device_count())'`). `{DEVICES_VARIABLE}=cpu` tries the code without a GPU.");
        println!("  - An RTX 50 card needs a PyTorch built for CUDA 12.8 or later (`python -c 'import torch; print(torch.cuda.get_arch_list())'` lists `sm_120` for it).");
    } else {
        println!("  - On Linux or Windows the program must be built with `cargo build --release --features cuda` (it needs the CUDA toolkit; see DEPLOY.md), or use PyTorch instead: `{PYTHON_VARIABLE}=/path/to/python`.");
        println!("  - `nvidia-smi` should list the cards, and the driver should be new enough for them (an RTX 50 card needs a driver of 570 or later, the open kernel modules, and CUDA toolkit 12.8 or later).");
        println!("  - Build for the oldest card of the machine, `CUDA_COMPUTE_CAP=86 cargo build --release --features cuda` for an RTX 3090 with an RTX 50 card beside it: see DEPLOY.md.");
    }
    println!("  - `{MODE_VARIABLE}=off` switches the GPUs off, on purpose.");
}

/// The same for neighbour-joining (sequences cubed). `None` for a backend that does not do it or a device that cannot (no 64-bit floats).
fn bench_joins(pool: &Target, quick: bool) -> Option<Option<f64>> {
    let probe = plane_distances(5, 1);
    match pool.joins(&probe)? {
        Ok(_) => {}
        Err(why) => {
            println!("Neighbour-joining: {why}; the processor does it.");
            return None;
        }
    }
    println!("Neighbour-joining (the processor, sharing the work among its cores when there are many sequences):");
    println!(
        "  {:>5} {:>11} {:>11} {:>8}",
        "seqs", "processor", "GPUs", "faster"
    );
    let sizes: &[usize] = if quick {
        &[100, 300, 600]
    } else {
        &[100, 300, 600, 1200]
    };
    let mut tried = Vec::new();
    for &n in sizes {
        let dist = plane_distances(n, 9);
        let start = Instant::now();
        let _ = cpu_joins(&dist);
        let cpu = start.elapsed().as_secs_f64();
        let _ = pool.joins(&dist);
        let start = Instant::now();
        match pool.joins(&dist)? {
            Ok(_) => {
                let gpu = start.elapsed().as_secs_f64();
                let ratio = cpu / gpu.max(1e-9);
                println!(
                    "  {n:>5} {:>8.1} ms {:>8.1} ms {ratio:>7.1}x",
                    1000.0 * cpu,
                    1000.0 * gpu
                );
                tried.push(((n * n * n) as f64, ratio >= 1.5));
            }
            Err(why) => println!("  {n:>5}  {why}"),
        }
    }
    Some(crossover(tried))
}

/// The same for pairwise alignments, which are **opt-in work**: the processor, with eight pairs to a vector instruction on every core, was
/// 60 to 1,000 times faster than a worker on the CPU device or a Mac's GPU, so they go to the cards only where a person sets
/// `HELIXLOOM_ACCEL_ALIGN_CELLS`; this says whether the cards of this machine are ever worth it (cells are those of all the tables
/// together). `None` for a backend that does not do it or a device that cannot.
fn bench_alignments(pool: &Target, quick: bool) -> Option<Option<f64>> {
    let sc = Scoring::default();
    let mut rng = Rng::new(21);
    let mut sequence =
        |len: usize| -> Vec<u8> { (0..len).map(|_| b"ACGT"[rng.below(4)]).collect() };
    let probe = (sequence(20), sequence(20));
    match pool.alignments(&probe.0, &[probe.1.as_slice()], &sc, true)? {
        Ok(_) => {}
        Err(why) => {
            println!("Alignments: {why}; the processor does them.");
            return None;
        }
    }
    println!("Pairwise alignments (the processor with all its cores, and the cards, which are used for them only on request):");
    println!(
        "  {:>13} {:>11} {:>11} {:>8}",
        "pairs x bases", "processor", "GPUs", "faster"
    );
    let sizes: &[(usize, usize)] = if quick {
        &[(32, 300), (64, 500)]
    } else {
        &[(32, 300), (64, 500), (128, 1000)]
    };
    let mut tried = Vec::new();
    for &(count, len) in sizes {
        let a = sequence(len);
        let bs: Vec<Vec<u8>> = (0..count).map(|_| sequence(len)).collect();
        let refs: Vec<&[u8]> = bs.iter().map(Vec::as_slice).collect();
        let _ = pool.alignments(&a, &refs, &sc, true);
        let start = Instant::now();
        let _ = crate::seq::align::align_many_on_processor(&a, &refs, &sc, true);
        let cpu = start.elapsed().as_secs_f64();
        let start = Instant::now();
        match pool.alignments(&a, &refs, &sc, true)? {
            Ok(_) => {
                let gpu = start.elapsed().as_secs_f64();
                let ratio = cpu / gpu.max(1e-9);
                println!(
                    "  {:>13} {:>8.1} ms {:>8.1} ms {ratio:>7.2}x",
                    format!("{count} x {len}"),
                    1000.0 * cpu,
                    1000.0 * gpu
                );
                tried.push((((len + 1) * (len + 1) * count) as f64, ratio >= 1.5));
            }
            Err(why) => println!("  {count} x {len}  {why}"),
        }
    }
    Some(crossover(tried))
}

pub fn run(args: &[String]) -> i32 {
    if args.iter().any(|a| a == "--help" || a == "-h") {
        print!("{USAGE}");
        return 0;
    }
    if let Some(other) = args
        .iter()
        .find(|a| !matches!(a.as_str(), "--bench" | "--quick"))
    {
        eprintln!("Unknown option {other}.\n\n{USAGE}");
        return 2;
    }
    let (bench, quick) = (
        args.iter().any(|a| a == "--bench"),
        args.iter().any(|a| a == "--quick"),
    );
    // cards numbered as `nvidia-smi` numbers them (by their place on the bus), as the server numbers them for its jobs: CUDA's own
    // order is fastest first, and then "card 1" here could be another card than the one `--gpus 1` means. Before any card is opened.
    if std::env::var_os("CUDA_DEVICE_ORDER").is_none() {
        std::env::set_var("CUDA_DEVICE_ORDER", "PCI_BUS_ID");
    }
    println!("HelixLoom GPU check");
    println!(
        "  {MODE_VARIABLE}: {:?}   {WEIGHTS_VARIABLE}: {}",
        mode(),
        std::env::var(WEIGHTS_VARIABLE).unwrap_or_else(|_| "(measured)".into())
    );
    println!(
        "  {BACKEND_VARIABLE}: {:?} (candle, or PyTorch where {PYTHON_VARIABLE} names a Python)",
        super::backend()
    );
    let pool = match Target::discover() {
        Ok(Some(pool)) => pool,
        Ok(None) => {
            no_gpu(None);
            return 2;
        }
        Err(why) => {
            no_gpu(Some(&why));
            return 2;
        }
    };
    println!("  Backend: {}", pool.backend_name());
    let names = crate::remote::gpus::detect();
    println!("  {} GPU(s) found:", pool.len());
    let cards = pool.cards();
    for (k, (name, likelihood_share, distance_share)) in cards.iter().enumerate() {
        // the name `nvidia-smi` gives, when it is not what the device is called already
        let smi = names
            .get(k)
            .filter(|g| names.len() == pool.len() && !name.contains(&g.name))
            .map_or(String::new(), |g| format!(" ({})", g.name));
        println!(
            "    {k}: {name}{smi}, share of the likelihood {:.0} %, of the distances {:.0} %",
            100.0 * likelihood_share,
            100.0 * distance_share
        );
    }
    let mut report = Report { failures: 0 };
    for (k, (name, _, _)) in cards.iter().enumerate() {
        check_pool(
            &mut report,
            &format!("Card {k} alone ({name})"),
            &pool.only(k),
        );
    }
    if pool.len() > 1 {
        check_pool(
            &mut report,
            "All cards together (the work cut by their shares)",
            &pool,
        );
    }
    if bench {
        println!("Benchmark, all cards together:");
        let likelihood_cells = bench_likelihood(&pool, quick);
        let distance_cells = bench_distances(&pool, quick);
        let join_cells = bench_joins(&pool, quick);
        let align_cells = bench_alignments(&pool, quick);
        println!("Suggested from this run, the least work from which the GPUs won every size tried (the program's own guesses are in `accel/mod.rs`):");
        match likelihood_cells {
            Some(c) => println!("  {LIKELIHOOD_VARIABLE}={c:.0}"),
            None => println!("  the GPUs were not 1.5 times faster for the largest likelihood tried: leave {LIKELIHOOD_VARIABLE} high, or use {MODE_VARIABLE}=off"),
        }
        match distance_cells {
            Some(c) => println!("  {DISTANCE_VARIABLE}={c:.0}"),
            None => println!("  the GPUs were not 1.5 times faster for the largest distances tried: leave {DISTANCE_VARIABLE} high"),
        }
        // these two start PyTorch workers, which take seconds: the suggestion does not count that, so take more than it says
        if let Some(found) = join_cells {
            match found {
                Some(c) => println!("  {}={c:.0}   (the workers' start not counted)", super::threshold::JOIN_VARIABLE),
                None => println!("  the GPUs were not 1.5 times faster for the largest neighbour-joining tried: leave {} high", super::threshold::JOIN_VARIABLE),
            }
        }
        if let Some(found) = align_cells {
            match found {
                Some(c) => println!("  {}={c:.0}   (the workers' start not counted; alignments go to the cards only with this set)", super::threshold::ALIGN_VARIABLE),
                None => println!("  the GPUs were not 1.5 times faster for the largest alignments tried: leave {} unset (the processor does them)", super::threshold::ALIGN_VARIABLE),
            }
        }
    }
    if report.failures == 0 {
        println!("All checks passed.");
        0
    } else {
        println!("{} check(s) FAILED: do not trust the GPUs on this machine until that is understood (HELIXLOOM_ACCEL=off uses the processor only).", report.failures);
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::accel::pool::{Member, Pool};
    use candle_core::Device;

    #[test]
    fn the_crossover_is_where_the_cards_won_every_size_from_there_up() {
        // lost, lost, won, won: from the third
        assert_eq!(
            crossover(vec![
                (1.0, false),
                (10.0, false),
                (100.0, true),
                (1000.0, true)
            ]),
            Some(100.0)
        );
        // a win among losses (the shape of the work, not its size, made it) does not count, and a loss at a bigger size bars it
        assert_eq!(
            crossover(vec![
                (1.0, true),
                (10.0, false),
                (100.0, true),
                (1000.0, true)
            ]),
            Some(100.0)
        );
        assert_eq!(crossover(vec![(1.0, true), (100.0, false)]), None);
        // the order they were tried in does not matter, and equal sizes with different results count as a loss
        assert_eq!(
            crossover(vec![(1000.0, true), (10.0, false), (100.0, true)]),
            Some(100.0)
        );
        assert_eq!(crossover(vec![(5.0, true), (5.0, false)]), None);
        assert_eq!(crossover(vec![(5.0, false), (5.0, true)]), None);
        assert_eq!(crossover(Vec::new()), None);
    }

    #[test]
    fn the_checks_pass_on_a_pool_that_computes_the_processors_numbers() {
        // a CPU device runs the very same tensor code as a card: the checks are what is tried here, and must pass
        let pool = Target::Candle(Pool::new(vec![
            Member::new(Device::Cpu, "CPU".into(), 2.0),
            Member::new(Device::Cpu, "CPU".into(), 1.0),
        ]));
        assert_eq!(check_distances(&pool).map(|(same, _)| same), Ok(true));
        assert_eq!(check_bootstrap(&pool), Ok(true));
        assert_eq!(check_threads(&pool), Ok(true));
        let a = likelihood::agree(&pool, 12, 300, 1).expect("an agreement");
        assert!(a.worst_relative <= LIKELIHOOD_TOLERANCE, "{a:?}");
        assert!(a.worst_absolute >= 0.0 && a.worst_absolute.is_finite());
    }

    #[test]
    fn the_checks_pass_on_pytorch_workers_too_and_one_card_alone_is_a_pool_of_its_own() {
        let Some(pool) = crate::accel::torch::tests::cpu_pool(&[2.0, 1.0]) else {
            return;
        };
        let pool = Target::Torch(pool);
        assert_eq!(pool.len(), 2);
        let cards = pool.cards();
        assert_eq!(cards.len(), 2);
        assert!(cards[0].1 > cards[1].1, "{cards:?}");
        let alone = pool.only(1);
        assert_eq!(alone.len(), 1);
        for target in [&pool, &alone] {
            assert_eq!(check_distances(target).map(|(same, _)| same), Ok(true));
            assert_eq!(check_bootstrap(target), Ok(true));
            assert_eq!(check_threads(target), Ok(true));
            let a = likelihood::agree(target, 12, 300, 1).expect("an agreement");
            assert!(a.worst_relative <= LIKELIHOOD_TOLERANCE, "{a:?}");
        }
    }
}
