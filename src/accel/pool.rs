//! The GPUs of this machine as one pool, for work that is shared out among cards that are not alike (an RTX 3090 beside an
//! RTX 5060 Ti, say): each member has a weight, the share of the work it is given, which is measured (how fast the card runs
//! the operations of this program) or set by the person (`HELIXLOOM_GPU_WEIGHTS=0.6,0.4`).

use crate::seq::phylocnn::gpu::devices;
use candle_core::{Device, Tensor};
use std::ops::Range;
use std::sync::{Arc, Mutex, MutexGuard, OnceLock, PoisonError};
use std::time::Instant;

/// `off` never uses the GPUs, `force` uses them whatever the size of the work (to try them), `auto` (the default) only when
/// the work is big enough to be worth it.
pub const MODE_VARIABLE: &str = "HELIXLOOM_ACCEL";
/// The weights of the cards, in the order of their numbers: `0.6,0.4`.
pub const WEIGHTS_VARIABLE: &str = "HELIXLOOM_GPU_WEIGHTS";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Off,
    Auto,
    Force,
}

/// What the person asked for with `HELIXLOOM_ACCEL`.
pub fn mode() -> Mode {
    parse_mode(&std::env::var(MODE_VARIABLE).unwrap_or_default())
}

fn parse_mode(text: &str) -> Mode {
    match text.trim().to_ascii_lowercase().as_str() {
        "off" | "0" | "no" | "false" | "cpu" => Mode::Off,
        "force" | "on" | "always" => Mode::Force,
        _ => Mode::Auto,
    }
}

#[derive(Clone)]
pub struct Member {
    pub device: Device,
    pub name: String,
    /// The share of the work, relative to the others (the weights of a pool add up to 1).
    pub weight: f64,
    /// Held while a thread works on the card. candle's Metal backend shares one command buffer among all the threads that
    /// use a device, and several threads at once got wrong answers, silently (the tests caught it): so one thread at a time
    /// works on a card, whatever the card. (Cards of the pool still work side by side, and a thread waits for a card only
    /// while another uses it.) Whether CUDA needs it too has not been tried; the cost where it does not is some overlap of
    /// the host's launching of kernels.
    busy: Arc<Mutex<()>>,
}

impl Member {
    pub fn new(device: Device, name: String, weight: f64) -> Member {
        Member {
            device,
            name,
            weight,
            busy: Arc::new(Mutex::new(())),
        }
    }

    /// Wait for the card and hold it until the guard goes out of scope. Not to be taken again by the thread that holds it.
    pub fn exclusive(&self) -> MutexGuard<'_, ()> {
        self.busy.lock().unwrap_or_else(PoisonError::into_inner)
    }
}

#[derive(Clone)]
pub struct Pool {
    members: Vec<Member>,
}

impl Pool {
    /// A pool of `members`, their weights made to add up to 1 (a weight that is not a positive number counts as the mean).
    pub fn new(mut members: Vec<Member>) -> Pool {
        let mut weights: Vec<f64> = members.iter().map(|m| m.weight).collect();
        normalize(&mut weights);
        for (m, w) in members.iter_mut().zip(weights) {
            m.weight = w;
        }
        Pool { members }
    }

    pub fn members(&self) -> &[Member] {
        &self.members
    }

    pub fn len(&self) -> usize {
        self.members.len()
    }

    #[cfg_attr(not(test), allow(dead_code))]
    pub fn is_empty(&self) -> bool {
        self.members.is_empty()
    }

    /// Every GPU this process can use, weighted by how fast each runs the operations of this program (or as
    /// `HELIXLOOM_GPU_WEIGHTS` says); `None` when there is none, or `HELIXLOOM_ACCEL=off`.
    pub fn discover() -> Option<Pool> {
        if mode() == Mode::Off {
            return None;
        }
        let found = devices();
        if found.is_empty() {
            return None;
        }
        let given = given_weights(found.len());
        let members = members_from(found, given.as_deref(), speed);
        (!members.is_empty()).then(|| Pool::new(members))
    }

    /// The pool of this process, found once: `None` when there is no GPU to use.
    pub fn global() -> Option<&'static Pool> {
        static POOL: OnceLock<Option<Pool>> = OnceLock::new();
        POOL.get_or_init(Pool::discover).as_ref()
    }

    /// Cut `total` units of work into one contiguous range for each member, in proportion to the weights (see `shares`).
    pub fn shares(&self, total: usize) -> Vec<Range<usize>> {
        let weights: Vec<f64> = self.members.iter().map(|m| m.weight).collect();
        shares(&weights, total)
    }

    /// Run `f` for every member at once, each on a thread of its own, and give the answers in the order of the members.
    pub fn run<T: Send>(&self, f: impl Fn(usize, &Member) -> T + Sync) -> Vec<T> {
        let f = &f;
        std::thread::scope(|s| {
            let jobs: Vec<_> = self
                .members
                .iter()
                .enumerate()
                .map(|(k, m)| s.spawn(move || f(k, m)))
                .collect();
            jobs.into_iter()
                .map(|j| j.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                .collect()
        })
    }
}

/// Make `weights` add up to 1: a weight that is not a positive number counts as the mean of the good ones (all equal when none is).
pub(super) fn normalize(weights: &mut [f64]) {
    let good: Vec<f64> = weights
        .iter()
        .copied()
        .filter(|w| w.is_finite() && *w > 0.0)
        .collect();
    let mean = if good.is_empty() {
        1.0
    } else {
        good.iter().sum::<f64>() / good.len() as f64
    };
    for w in weights.iter_mut() {
        if !(w.is_finite() && *w > 0.0) {
            *w = mean;
        }
    }
    let total: f64 = weights.iter().sum();
    for w in weights.iter_mut() {
        *w /= total.max(f64::MIN_POSITIVE);
    }
}

/// Cut `total` units of work into one contiguous range for each of the (normalised) `weights`, in proportion to them; what the
/// rounding down leaves over goes to those with the largest fractions, so that the ranges cover `0..total` exactly.
pub(super) fn shares(weights: &[f64], total: usize) -> Vec<Range<usize>> {
    let ideal: Vec<f64> = weights.iter().map(|w| w * total as f64).collect();
    let mut sizes: Vec<usize> = ideal.iter().map(|x| x.floor() as usize).collect();
    // what the rounding down left over is fewer units than there are members: one more each to those with the largest fractions
    let left = total - sizes.iter().sum::<usize>().min(total);
    let mut by_fraction: Vec<usize> = (0..sizes.len()).collect();
    by_fraction.sort_by(|&a, &b| {
        (ideal[b] - ideal[b].floor())
            .total_cmp(&(ideal[a] - ideal[a].floor()))
            .then(a.cmp(&b))
    });
    for &k in by_fraction.iter().take(left) {
        sizes[k] += 1;
    }
    let mut start = 0;
    sizes
        .into_iter()
        .map(|n| {
            let range = start..start + n;
            start += n;
            range
        })
        .collect()
}

/// The members for `found` cards: weighted as `given` says, or by what `measure` makes of each. A card that `measure` cannot get a
/// number for cannot run the operations of this program (a toolkit too old for it, say), and is left out, with a line that says so:
/// kept, it would fail every share of the work it was given and take the other cards' work down with it.
fn members_from(
    found: Vec<Device>,
    given: Option<&[f64]>,
    measure: impl Fn(&Device) -> Option<f64>,
) -> Vec<Member> {
    found
        .into_iter()
        .enumerate()
        .filter_map(|(k, device)| {
            let weight = match given {
                Some(w) => w[k],
                None => match measure(&device) {
                    Some(w) => w,
                    None => {
                        super::note_failure(
                            &format!("{} is left out of the pool", label(&device, k)),
                            "it could not run the operations of the speed test (an older CUDA toolkit than the card needs, say)",
                        );
                        return None;
                    }
                },
            };
            Some(Member::new(device.clone(), label(&device, k), weight))
        })
        .collect()
}

/// What a card is called: its number, and the kind of device.
fn label(device: &Device, k: usize) -> String {
    match device {
        Device::Cuda(_) => format!("CUDA {k}"),
        Device::Metal(_) => "Apple GPU (Metal)".to_string(),
        Device::Cpu => "CPU".to_string(),
    }
}

/// The weights written in `HELIXLOOM_GPU_WEIGHTS`, when there is one for each of `count` cards and all are positive.
pub(super) fn given_weights(count: usize) -> Option<Vec<f64>> {
    parse_weights(&std::env::var(WEIGHTS_VARIABLE).ok()?, count)
}

fn parse_weights(text: &str, count: usize) -> Option<Vec<f64>> {
    let weights: Vec<f64> = text
        .split(',')
        .map(|w| w.trim().parse::<f64>())
        .collect::<Result<_, _>>()
        .ok()?;
    (weights.len() == count && weights.iter().all(|w| w.is_finite() && *w > 0.0)).then_some(weights)
}

/// How fast `device` runs the operations of the likelihood and the distances: one over the seconds that a few of their kinds
/// (elementwise arithmetic over a big array, a maximum over its rows, a matrix product) take, after one run to warm it up
/// (the first run of a card compiles its kernels). `None` if the device cannot run them.
pub fn speed(device: &Device) -> Option<f64> {
    const ROWS: usize = 1 << 20;
    let run = || -> candle_core::Result<f64> {
        let v = Tensor::ones((ROWS, 4), candle_core::DType::F32, device)?;
        let p = Tensor::ones((1, 4, 4), candle_core::DType::F32, device)?;
        let a = Tensor::ones((2048, 2048), candle_core::DType::F32, device)?;
        let start = Instant::now();
        for _ in 0..4 {
            let mixed = v.unsqueeze(1)?.broadcast_mul(&p)?.sum(2)?;
            let top = mixed.max_keepdim(1)?;
            let scaled = mixed.broadcast_div(&top)?;
            let _ = scaled.log()?;
            let _ = a.matmul(&a)?;
        }
        // the answer must be read back, or the card may not have started
        let _ = v.sum_all()?.to_scalar::<f32>()?;
        Ok(start.elapsed().as_secs_f64())
    };
    run().ok()?;
    let seconds = run().ok()?;
    (seconds > 0.0).then(|| 1.0 / seconds)
}

#[cfg(test)]
mod tests;
