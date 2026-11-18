//! The PyTorch backend: the same work as the candle backend (distances between sequences, the likelihood of a tree), done by Python
//! worker processes that run PyTorch, one for every GPU the job was given. The program keeps everything else (the search, the
//! queue, the sums); a worker keeps the columns of the alignment on its card and answers for the counts of a matrix of distances
//! and for the likelihood of a tree over a range of the patterns.
//!
//! Where candle needs the CUDA toolkit to build the program (and a build for the oldest of two different cards), this needs only a
//! Python with PyTorch on the workstation: `HELIXLOOM_TORCH_PYTHON=/path/to/python`. The path is the owner's, in the server's
//! environment, as `--fold-command` is: nothing from a request is ever a program.
//!
//! A worker serves one request at a time, so a card is used by one thread at a time (as `Member::exclusive` makes candle's), and
//! cards work side by side because each has its own worker.

pub mod align;
pub mod distance;
pub mod nj;
pub mod wire;

use super::note_failure;
use super::pool::{given_weights, normalize, shares};
use std::io::BufReader;
use std::ops::Range;
use std::process::{Child, ChildStdin, Command, Stdio};
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::mpsc::{channel, Receiver, RecvTimeoutError};
use std::sync::{Arc, Mutex, OnceLock, PoisonError};
use std::time::Duration;
use wire::Message;

/// The Python (with PyTorch installed) that the workers run in.
pub const PYTHON_VARIABLE: &str = "HELIXLOOM_TORCH_PYTHON";
/// The devices to use, instead of every card PyTorch finds: `cuda:0,cuda:1`, `mps`, or `cpu,cpu` (to try the code without a GPU).
pub const DEVICES_VARIABLE: &str = "HELIXLOOM_TORCH_DEVICES";
/// What a worker runs: `worker.py`, which is part of the program.
const SCRIPT: &str = include_str!("worker.py");
/// How long a worker may take to start (importing PyTorch cold, opening the card, a first kernel): then it is killed.
const START_TIMEOUT: Duration = Duration::from_secs(300);
/// How long one request may take. A GPU that has stopped answering would otherwise hold the job until its time limit; the longest
/// real request (the distances of a thousand sequences by a million columns on one card) is far shorter.
const CALL_TIMEOUT: Duration = Duration::from_secs(1800);
/// The least patterns each worker is given when a likelihood is shared out: below that, a card is not the limit, the launching of
/// its kernels is (about ten operations for each inner node of the tree whatever the number of patterns), and cutting the patterns
/// among workers adds a trip over a pipe for each and no speed. An estimate for a 3090 and a 5060 Ti, not measured on them.
const LEAST_PATTERNS_EACH: usize = 200_000;

/// The python the owner named, if any.
pub fn configured_python() -> Option<String> {
    std::env::var(PYTHON_VARIABLE)
        .ok()
        .map(|p| p.trim().to_string())
        .filter(|p| !p.is_empty())
}

/// Whether a job process could have a device for PyTorch: not when it was shown no GPU (`CUDA_VISIBLE_DEVICES` set to nothing, as the
/// server does for a job that was given none) and no device was named. Starting workers there would only import PyTorch, for seconds,
/// to find none.
pub fn may_have_devices() -> bool {
    can_have_devices(
        std::env::var("CUDA_VISIBLE_DEVICES").ok().as_deref(),
        std::env::var(DEVICES_VARIABLE).is_ok_and(|d| !d.trim().is_empty()),
    )
}

fn can_have_devices(visible: Option<&str>, devices_named: bool) -> bool {
    devices_named || visible.is_none_or(|v| !v.trim().is_empty())
}

/// What a worker said about itself as it started.
#[derive(Clone, Debug)]
pub struct Hello {
    /// How many devices of its kind PyTorch found (0: nothing to compute on).
    pub count: usize,
    pub name: String,
    pub kind: String,
    pub torch: String,
    pub cuda: String,
    pub arches: String,
}

struct Io {
    stdin: ChildStdin,
    /// What the worker sends, read by a thread of its own, so that waiting for it can time out.
    replies: Receiver<std::io::Result<Message>>,
}

pub struct Worker {
    /// The device and what it is: `cuda:1 NVIDIA GeForce RTX 5060 Ti (16311 MB, sm_120)`.
    pub name: String,
    /// `None` once the worker has failed or stopped: it is not asked again.
    io: Mutex<Option<Io>>,
    child: Mutex<Child>,
    /// How long one request may take before the worker is taken to be hung and is killed.
    timeout: Duration,
}

impl Worker {
    /// Start a worker on `device` (`auto`: the first CUDA card, else the Apple GPU), and read what it says about itself.
    fn spawn(python: &str, device: &str) -> Result<(Worker, Hello), String> {
        Worker::spawn_script(python, SCRIPT, &[], device, START_TIMEOUT, CALL_TIMEOUT)
    }

    /// The same running `script` as `python -c` runs it, with `env` added to its environment (the tests wrap `worker.py`, which they
    /// hand over in a variable, in code of their own).
    fn spawn_script(
        python: &str,
        script: &str,
        env: &[(&str, &str)],
        device: &str,
        start_timeout: Duration,
        timeout: Duration,
    ) -> Result<(Worker, Hello), String> {
        let mut child = Command::new(python)
            .args(["-I", "-u", "-c", script, device])
            .envs(env.iter().copied())
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            // what PyTorch says goes to the job's log
            .stderr(Stdio::inherit())
            .spawn()
            .map_err(|e| format!("Cannot start {python}: {e}"))?;
        let (Some(stdin), Some(stdout)) = (child.stdin.take(), child.stdout.take()) else {
            return Err(format!("Cannot talk to {python}."));
        };
        // a thread reads the worker's messages for as long as it lives, and ends when the pipe does
        let (send, replies) = channel();
        std::thread::spawn(move || {
            let mut out = BufReader::new(stdout);
            loop {
                let next = Message::read_from(&mut out);
                let failed = next.is_err();
                if send.send(next).is_err() || failed {
                    return;
                }
            }
        });
        let fail = |child: &mut Child, why: String| -> Result<(Worker, Hello), String> {
            let _ = child.kill();
            let _ = child.wait();
            Err(why)
        };
        let said = match replies.recv_timeout(start_timeout) {
            Ok(Ok(said)) => said,
            Ok(Err(e)) => {
                return fail(
                    &mut child,
                    format!("The PyTorch worker did not start ({e}): is PyTorch (and NumPy) installed for {python}, as `-I` runs it, without the user's own packages? Its own message is in the log."),
                )
            }
            Err(_) => {
                return fail(
                    &mut child,
                    format!("The PyTorch worker did not start within {} seconds.", start_timeout.as_secs()),
                )
            }
        };
        let hello = Hello {
            count: said
                .u32s("count")
                .and_then(|c| c.first().copied())
                .unwrap_or(0) as usize,
            name: said.text("name").unwrap_or_default(),
            kind: said.text("kind").unwrap_or_default(),
            torch: said.text("torch").unwrap_or_default(),
            cuda: said.text("cuda").unwrap_or_default(),
            arches: said.text("arches").unwrap_or_default(),
        };
        // the first worker is started on `auto`: what that came to is its name
        let place = match (device, hello.kind.as_str()) {
            ("auto", "cuda") => "cuda:0".to_string(),
            ("auto", kind) => kind.to_string(),
            (device, _) => device.to_string(),
        };
        let worker = Worker {
            name: format!("{place} {}", hello.name),
            io: Mutex::new(Some(Io { stdin, replies })),
            child: Mutex::new(child),
            timeout,
        };
        Ok((worker, hello))
    }

    /// Send `message` and wait for the answer (one request at a time on a worker). An error says what failed; a worker whose pipe
    /// broke, or that did not answer in time (and is killed), is not used again, and one that only reported an error (out of memory,
    /// say) stays up.
    pub fn call(&self, message: &Message) -> Result<Message, String> {
        let mut guard = self.io.lock().unwrap_or_else(PoisonError::into_inner);
        let Some(io) = guard.as_mut() else {
            return Err(format!("The PyTorch worker for {} has stopped.", self.name));
        };
        let sent = message
            .write_to(&mut io.stdin)
            .map_err(|e| e.to_string())
            .and_then(|()| match io.replies.recv_timeout(self.timeout) {
                Ok(Ok(reply)) => Ok(reply),
                Ok(Err(e)) => Err(e.to_string()),
                Err(RecvTimeoutError::Timeout) => Err(format!(
                    "it did not answer within {} seconds, and was stopped",
                    self.timeout.as_secs()
                )),
                Err(RecvTimeoutError::Disconnected) => Err("its output ended".to_string()),
            });
        match sent {
            Ok(reply) if reply.op == "ok" => Ok(reply),
            Ok(reply) => Err(format!(
                "{}: {}",
                self.name,
                reply.text("text").unwrap_or_else(|| "an error".into())
            )),
            Err(e) => {
                *guard = None;
                // a worker that hung is not left running
                let _ = self
                    .child
                    .lock()
                    .unwrap_or_else(PoisonError::into_inner)
                    .kill();
                Err(format!("The PyTorch worker for {} stopped: {e}", self.name))
            }
        }
    }
}

impl Drop for Worker {
    fn drop(&mut self) {
        // closing its input ends it; a worker stuck in a kernel is killed
        *self.io.lock().unwrap_or_else(PoisonError::into_inner) = None;
        let mut child = self.child.lock().unwrap_or_else(PoisonError::into_inner);
        let _ = child.kill();
        let _ = child.wait();
    }
}

/// The workers of a job, one for each card, with the share of the work each is given.
pub struct TorchPool {
    workers: Vec<Arc<Worker>>,
    /// For the likelihood (limited by memory speed) and for the distances (by matrix products): the shares differ between cards that
    /// are unlike, and add up to 1 each.
    likelihood_weights: Vec<f64>,
    distance_weights: Vec<f64>,
    /// What PyTorch and its CUDA are, for the report.
    pub versions: String,
}

impl TorchPool {
    /// Start workers and measure the cards. `Ok(None)` when PyTorch finds no device to compute on. `devices` is what
    /// `HELIXLOOM_TORCH_DEVICES` said (none: every card PyTorch finds).
    pub fn start(python: &str, devices: Option<Vec<String>>) -> Result<Option<TorchPool>, String> {
        TorchPool::start_weighted(python, devices, None)
    }

    /// The same with the shares of the work given (one for each device, as `HELIXLOOM_GPU_WEIGHTS` gives them) instead of measured.
    pub fn start_weighted(
        python: &str,
        devices: Option<Vec<String>>,
        weights: Option<Vec<f64>>,
    ) -> Result<Option<TorchPool>, String> {
        let (first, hello) = Worker::spawn(python, devices.as_ref().map_or("auto", |d| &d[0]))?;
        if hello.count == 0 {
            note_failure(
                "PyTorch found no device to compute on",
                &format!(
                    "PyTorch {}, CUDA {}, built for [{}], CUDA_VISIBLE_DEVICES={:?}: a PyTorch built without CUDA, a driver it cannot use, or no GPU given to this job",
                    hello.torch,
                    hello.cuda,
                    hello.arches,
                    std::env::var("CUDA_VISIBLE_DEVICES").unwrap_or_else(|_| "(not set)".into())
                ),
            );
            return Ok(None);
        }
        let specs: Vec<String> = match devices {
            Some(list) => list,
            None if hello.kind == "cuda" => (0..hello.count).map(|k| format!("cuda:{k}")).collect(),
            None => vec![hello.kind.clone()],
        };
        let versions = format!(
            "PyTorch {}{}{}",
            hello.torch,
            if hello.cuda == "None" {
                String::new()
            } else {
                format!(", CUDA {}", hello.cuda)
            },
            if hello.arches.is_empty() {
                String::new()
            } else {
                format!(" (built for {})", hello.arches)
            }
        );
        // the others start at the same time: importing PyTorch takes seconds
        let mut workers = vec![first];
        let rest: Vec<Result<Worker, String>> = std::thread::scope(|s| {
            let jobs: Vec<_> = specs[1..]
                .iter()
                .map(|spec| s.spawn(move || Worker::spawn(python, spec).map(|(w, _)| w)))
                .collect();
            jobs.into_iter()
                .map(|j| j.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                .collect()
        });
        for started in rest {
            match started {
                Ok(w) => workers.push(w),
                Err(why) => note_failure("a PyTorch worker did not start", &why),
            }
        }
        let given = weights
            .filter(|w| w.len() == workers.len())
            .or_else(|| given_weights(workers.len()));
        let (likelihood_weights, distance_weights, workers) = match given {
            Some(w) => (w.clone(), w, workers),
            None => measure(workers),
        };
        if workers.is_empty() {
            return Ok(None);
        }
        let (mut likelihood_weights, mut distance_weights) = (likelihood_weights, distance_weights);
        normalize(&mut likelihood_weights);
        normalize(&mut distance_weights);
        let pool = TorchPool {
            workers: workers.into_iter().map(Arc::new).collect(),
            likelihood_weights,
            distance_weights,
            versions,
        };
        // what was found and how the work will be cut, in the job's log: the first thing to look at when a card seems idle
        eprintln!(
            "HelixLoom: {}: share of the likelihood {}; of the distances {}",
            pool.versions,
            pool.likelihood_description(),
            pool.distance_description()
        );
        Ok(Some(pool))
    }

    /// The workers that the environment asks for: the python named by `HELIXLOOM_TORCH_PYTHON` (else `python3`) on the devices that
    /// `HELIXLOOM_TORCH_DEVICES` lists (else every card PyTorch finds).
    pub fn from_environment() -> Result<Option<TorchPool>, String> {
        let python = configured_python().unwrap_or_else(|| "python3".into());
        let devices: Vec<String> = std::env::var(DEVICES_VARIABLE)
            .unwrap_or_default()
            .split(',')
            .map(|d| d.trim().to_string())
            .filter(|d| !d.is_empty())
            .collect();
        TorchPool::start(&python, (!devices.is_empty()).then_some(devices))
    }

    /// The pool of this process, started the first time it is asked for: `None` when it has no device, or could not be started (said on
    /// the standard error). Work that is not big enough never asks, so a small job pays nothing for it.
    pub fn global() -> Option<&'static Arc<TorchPool>> {
        static POOL: OnceLock<Option<Arc<TorchPool>>> = OnceLock::new();
        POOL.get_or_init(|| match TorchPool::from_environment() {
            Ok(pool) => pool.map(Arc::new),
            Err(why) => {
                note_failure("PyTorch could not be started", &why);
                None
            }
        })
        .as_ref()
    }

    /// Only worker `k`, with all the work: to try one card alone (it shares the worker, which goes on until the last pool that has it ends).
    pub fn only(&self, k: usize) -> TorchPool {
        TorchPool {
            workers: vec![Arc::clone(&self.workers[k])],
            likelihood_weights: vec![1.0],
            distance_weights: vec![1.0],
            versions: self.versions.clone(),
        }
    }

    pub fn len(&self) -> usize {
        self.workers.len()
    }

    /// The worker that is fastest at the likelihood (limited by memory speed, as neighbour-joining is): for work that is not cut among
    /// workers.
    pub fn fastest(&self) -> usize {
        (0..self.workers.len())
            .max_by(|&a, &b| {
                self.likelihood_weights[a]
                    .total_cmp(&self.likelihood_weights[b])
                    .then(b.cmp(&a))
            })
            .unwrap_or(0)
    }

    pub fn workers(&self) -> &[Arc<Worker>] {
        &self.workers
    }

    /// The ranges of `total` patterns that each worker does, by how fast each runs the likelihood, but no worker is given fewer than
    /// `LEAST_PATTERNS_EACH`: a likelihood of fewer goes whole to the fastest worker, or to as many of the fastest as it has enough for.
    pub fn likelihood_shares(&self, total: usize) -> Vec<Range<usize>> {
        split(&self.likelihood_weights, total, LEAST_PATTERNS_EACH)
    }

    /// The same by how fast each runs matrix products, for the distances.
    pub fn distance_shares(&self, total: usize) -> Vec<Range<usize>> {
        shares(&self.distance_weights, total)
    }

    /// For each worker, its name and its shares of the likelihood and of the distances.
    pub fn cards(&self) -> Vec<(String, f64, f64)> {
        self.workers
            .iter()
            .zip(self.likelihood_weights.iter().zip(&self.distance_weights))
            .map(|(w, (l, d))| (w.name.clone(), *l, *d))
            .collect()
    }

    /// Which devices, and what share of the work each does, in words.
    pub fn describe(&self, weights: &[f64]) -> String {
        self.workers
            .iter()
            .zip(weights)
            .map(|(w, share)| format!("{} {:.0} %", w.name, 100.0 * share))
            .collect::<Vec<_>>()
            .join(", ")
    }

    pub fn likelihood_description(&self) -> String {
        self.describe(&self.likelihood_weights)
    }

    pub fn distance_description(&self) -> String {
        self.describe(&self.distance_weights)
    }

    /// Run `f` for every worker at once, each on a thread of its own, and give the answers in the order of the workers.
    pub fn run<T: Send>(&self, f: impl Fn(usize, &Worker) -> T + Sync) -> Vec<T> {
        let f = &f;
        std::thread::scope(|s| {
            let jobs: Vec<_> = self
                .workers
                .iter()
                .enumerate()
                .map(|(k, w)| s.spawn(move || f(k, w.as_ref())))
                .collect();
            jobs.into_iter()
                .map(|j| j.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                .collect()
        })
    }
}

/// `total` units cut among the `weights`, as `shares` cuts them, but among only as many of the heaviest as `least` units each allows
/// (at least one): the others get nothing.
fn split(weights: &[f64], total: usize, least: usize) -> Vec<Range<usize>> {
    let mut weights = weights.to_vec();
    normalize(&mut weights);
    let used = (total / least.max(1)).clamp(1, weights.len().max(1));
    if used >= weights.len() {
        return shares(&weights, total);
    }
    let mut order: Vec<usize> = (0..weights.len()).collect();
    order.sort_by(|&a, &b| weights[b].total_cmp(&weights[a]).then(a.cmp(&b)));
    let chosen: Vec<usize> = order.into_iter().take(used).collect();
    let mut among: Vec<f64> = chosen.iter().map(|&k| weights[k]).collect();
    normalize(&mut among);
    let mut cut = shares(&among, total).into_iter();
    let mut out = vec![0..0; weights.len()];
    // contiguous ranges, in the order of the workers: those that are not used get empty ones at the end of what came before
    let mut taken: Vec<(usize, Range<usize>)> = chosen.iter().copied().zip(cut.by_ref()).collect();
    taken.sort_by_key(|(k, _)| *k);
    let mut start = 0;
    for (k, range) in taken {
        out[k] = start..start + range.len();
        start += range.len();
    }
    out
}

/// How fast each worker's device is, for the likelihood and for the distances. A worker that cannot run the test (a CUDA library too
/// old for the card, say) is left out, with a line that says so: kept, it would fail every share of the work it was given.
fn measure(workers: Vec<Worker>) -> (Vec<f64>, Vec<f64>, Vec<Worker>) {
    let speeds: Vec<Result<(f64, f64), String>> = std::thread::scope(|s| {
        let jobs: Vec<_> = workers
            .iter()
            .map(|w| {
                s.spawn(move || {
                    let reply = w.call(&Message::new("speed"))?;
                    match (reply.f64s("memory"), reply.f64s("matmul")) {
                        (Some([m, ..]), Some([x, ..])) => Ok((*m, *x)),
                        _ => Err("it gave no speed".to_string()),
                    }
                })
            })
            .collect();
        jobs.into_iter()
            .map(|j| j.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    });
    let (mut likelihood, mut distance, mut kept) = (Vec::new(), Vec::new(), Vec::new());
    for (worker, speed) in workers.into_iter().zip(speeds) {
        match speed {
            Ok((m, x)) => {
                likelihood.push(m);
                distance.push(x);
                kept.push(worker);
            }
            Err(why) => note_failure(
                &format!("{} is left out of the pool", worker.name),
                &format!("it could not run the speed test: {why}"),
            ),
        }
    }
    (likelihood, distance, kept)
}

static NEXT_SET: AtomicU32 = AtomicU32::new(1);

/// A number for a set of columns that the workers keep, not used before in this process.
pub fn next_set() -> u32 {
    NEXT_SET.fetch_add(1, Ordering::Relaxed)
}

/// Tell every worker to forget a set of columns (their memory comes back). A worker that has stopped does not matter.
pub(crate) fn release(pool: &TorchPool, id: u32) {
    for w in pool.workers() {
        let _ = w.call(&Message::new("release").with_number("id", id as usize));
    }
}

#[cfg(test)]
pub(crate) mod tests;
