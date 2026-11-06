//! The server's jobs: the queue, the GPUs handed out one job at a time, the process that runs each job, and everything
//! kept on disk so that results wait for the person after the job is done, and survive a restart of the server.
//!
//! A job lives in a folder of its own, `<data>/jobs/<id>/`: `meta.txt` (one line: what the list shows), `request.bin`
//! (what was sent), `result.bin` (what came back), `log.txt` (what the job's process printed). Each job runs in a new
//! process, `helixloom run-job`, so that a job that crashes, runs too long or is cancelled cannot harm the server, and
//! a GPU job can be shown exactly one GPU.

mod config;
mod info;

use super::gpus::Gpu;
use super::one_line;
use super::tasks::Kind;
use info::{counter_of, made_order};
use std::collections::{BTreeMap, HashSet, VecDeque};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, ExitStatus, Stdio};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Condvar, Mutex, MutexGuard, PoisonError};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

pub use config::Config;
pub use info::{valid_id, JobInfo, Status};

/// How often a running job's process is looked at.
const WATCH: Duration = Duration::from_millis(200);
/// How often finished jobs older than the retention time are cleared away.
const SWEEP: Duration = Duration::from_secs(3600);
const REQUEST_FILE: &str = "request.bin";
const RESULT_FILE: &str = "result.bin";
const ERROR_FILE: &str = "error.txt";
const META_FILE: &str = "meta.txt";
const LOG_FILE: &str = "log.txt";

/// What a worker may use for a job.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Slot {
    Cpu,
    Gpu(usize),
}

/// Hands out the GPUs one job at a time, and a limited number of CPU jobs at once.
#[derive(Debug)]
pub struct Scheduler {
    free_gpus: Vec<usize>,
    free_cpu: usize,
}

impl Scheduler {
    pub fn new(gpus: &[usize], cpu_jobs: usize) -> Scheduler {
        let mut free_gpus = gpus.to_vec();
        free_gpus.sort_unstable();
        free_gpus.dedup();
        Scheduler {
            free_gpus,
            free_cpu: cpu_jobs.max(1),
        }
    }

    pub fn can_take(&self, needs_gpu: bool) -> bool {
        if needs_gpu {
            !self.free_gpus.is_empty()
        } else {
            self.free_cpu > 0
        }
    }

    /// The lowest free GPU for a job that needs one, or one of the CPU places.
    pub fn take(&mut self, needs_gpu: bool) -> Option<Slot> {
        if !self.can_take(needs_gpu) {
            return None;
        }
        if needs_gpu {
            Some(Slot::Gpu(self.free_gpus.remove(0)))
        } else {
            self.free_cpu -= 1;
            Some(Slot::Cpu)
        }
    }

    pub fn give_back(&mut self, slot: Slot) {
        match slot {
            Slot::Cpu => self.free_cpu += 1,
            Slot::Gpu(id) => {
                self.free_gpus.push(id);
                self.free_gpus.sort_unstable();
            }
        }
    }
}

struct State {
    jobs: BTreeMap<String, JobInfo>,
    queue: VecDeque<String>,
    scheduler: Scheduler,
    /// Running jobs somebody asked to stop.
    cancel: HashSet<String>,
}

/// How a job's process ended.
enum Outcome {
    Done,
    Failed(String),
    Cancelled,
}

/// Why a job could not be handed over.
#[derive(Debug, PartialEq, Eq)]
pub enum Lookup {
    NotFound,
    /// The job exists but its result does not yet (or never will).
    NotReady(Status),
    Failed(String),
}

/// What the status page shows.
pub struct Overview {
    pub running: usize,
    pub queued: usize,
    /// Each GPU and the job using it, if any.
    pub gpus: Vec<(Gpu, Option<String>)>,
}

pub struct Shared {
    pub config: Config,
    // ponytail: one lock for the whole table of jobs, which is plenty for the hundreds of jobs one workstation keeps;
    // the table is read whole for every list, so thousands of jobs would want an index
    state: Mutex<State>,
    wake: Condvar,
    counter: AtomicU64,
    /// Bytes of uploads being taken in now: a body is held in memory until the job it makes is stored.
    uploading: AtomicUsize,
}

/// Room for uploads is never less than this, whatever the largest upload is set to.
const MIN_UPLOAD_ROOM: usize = 3 << 20;
/// How many times over an upload is in memory while it is being taken in: the body, the container made from it, and the
/// sequences read out of that.
const UPLOAD_COPIES: usize = 3;

/// Room taken for an upload; given back when this goes out of scope.
pub struct UploadRoom<'a> {
    room: &'a AtomicUsize,
    bytes: usize,
}

impl Drop for UploadRoom<'_> {
    fn drop(&mut self) {
        self.room.fetch_sub(self.bytes, Ordering::SeqCst);
    }
}

fn now() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs())
}

/// Write a file so that it is either all there or not there: to a temporary name, then renamed.
fn write_atomic(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    let tmp = path.with_extension("part");
    std::fs::write(&tmp, bytes)?;
    std::fs::rename(&tmp, path)
}

/// Make the folder for the owner only. A folder that is already there is left as it is (it may be one the person chose, a
/// shared one even); only the folders made now are closed to everybody else.
#[cfg(unix)]
fn private_dir(path: &Path) -> std::io::Result<()> {
    use std::os::unix::fs::DirBuilderExt;
    std::fs::DirBuilder::new()
        .recursive(true)
        .mode(0o700)
        .create(path)
}

#[cfg(not(unix))]
fn private_dir(path: &Path) -> std::io::Result<()> {
    std::fs::create_dir_all(path)
}

/// What the folders of jobs held when the server started.
struct Found {
    jobs: BTreeMap<String, JobInfo>,
    /// The jobs that wait, oldest first.
    queue: VecDeque<String>,
    /// The highest number in any folder's name, readable or not: a new job never takes the number of a folder that is
    /// already there.
    highest: u64,
}

/// Read the jobs the folders under `jobs_dir` hold, and settle those that were not over when the server stopped.
fn read_jobs(jobs_dir: &Path, keep: Duration) -> Result<Found, String> {
    let entries = std::fs::read_dir(jobs_dir)
        .map_err(|e| format!("Could not read {}: {e}", jobs_dir.display()))?;
    let mut jobs = BTreeMap::new();
    let mut waiting: Vec<((u64, u64), String)> = Vec::new();
    let mut highest = 0;
    for entry in entries.flatten() {
        let dir = entry.path();
        let id = entry.file_name().to_string_lossy().into_owned();
        if !valid_id(&id) {
            continue;
        }
        highest = highest.max(counter_of(&id));
        let Some(job) = read_job(&dir, &id) else {
            clear_if_old(&dir, keep);
            continue;
        };
        let job = settled(job, &dir);
        if job.status == Status::Queued {
            waiting.push((made_order(job.created, &id), id.clone()));
        }
        jobs.insert(id, job);
    }
    waiting.sort();
    Ok(Found {
        jobs,
        queue: waiting.into_iter().map(|(_, id)| id).collect(),
        highest,
    })
}

/// The job a folder holds, if its `meta.txt` can be read and is the one for the folder's own id.
fn read_job(dir: &Path, id: &str) -> Option<JobInfo> {
    std::fs::read_to_string(dir.join(META_FILE))
        .ok()
        .and_then(|t| JobInfo::from_line(&t))
        .filter(|j| j.id == id)
}

/// A folder of this server that holds no readable job (a crash left it half made) is cleared away once it is as old as a
/// result may get.
fn clear_if_old(dir: &Path, keep: Duration) {
    let age = std::fs::metadata(dir)
        .and_then(|m| m.modified())
        .ok()
        .and_then(|t| t.elapsed().ok());
    if age.is_some_and(|age| age > keep) {
        let _ = std::fs::remove_dir_all(dir);
    }
}

/// A job that was running when the server stopped did not finish, unless its result is there (a result only appears whole,
/// by a rename, so one that is there is complete: the server stopped just before it said so); one that was waiting still
/// waits (unless its request is gone).
fn settled(job: JobInfo, dir: &Path) -> JobInfo {
    let lost = job.status == Status::Running
        || (job.status == Status::Queued && !dir.join(REQUEST_FILE).is_file());
    if !lost {
        return job;
    }
    let finished = now();
    let job = if dir.join(RESULT_FILE).is_file() {
        JobInfo {
            status: Status::Done,
            finished,
            ..job
        }
    } else {
        JobInfo {
            status: Status::Failed,
            finished,
            error: "The server stopped while this job was running.".into(),
            ..job
        }
    };
    let _ = write_atomic(&dir.join(META_FILE), job.to_line().as_bytes());
    job
}

/// Stop a job's process and everything it started (its process group on Unix, its tree on Windows).
// ponytail: on Unix this asks the shell for its own `kill` (there is no killpg in std, and a minimal system, a container
// say, may have no `kill` program); a crate such as libc would do it without starting a process.
#[cfg_attr(not(any(unix, windows)), allow(unused_variables))]
pub(super) fn kill_tree(pid: u32) {
    #[cfg(unix)]
    let _ = Command::new("sh")
        .args(["-c", "kill -KILL -- \"-$1\"", "sh", &pid.to_string()])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();
    #[cfg(windows)]
    let _ = Command::new("taskkill")
        .args(["/PID", &pid.to_string(), "/T", "/F"])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status();
}

/// Stop a job's process and what it started, and wait until it is gone.
fn stop_process(child: &mut Child) {
    kill_tree(child.id());
    let _ = child.kill();
    let _ = child.wait();
}

/// How a job's process ended by itself: with a whole result, or with the reason the job wrote down.
fn ended(dir: &Path, status: ExitStatus) -> Outcome {
    if dir.join(RESULT_FILE).is_file() && status.success() {
        return Outcome::Done;
    }
    let said = std::fs::read_to_string(dir.join(ERROR_FILE)).unwrap_or_default();
    Outcome::Failed(if said.trim().is_empty() {
        format!("The job's process stopped with {status}.")
    } else {
        said
    })
}

impl Shared {
    /// Open the job store in `config.data_dir` (reading the jobs left there) and start the workers.
    /// With `run_workers` false jobs only queue: for tests of everything but running.
    pub fn open(config: Config, run_workers: bool) -> Result<Arc<Shared>, String> {
        let jobs_dir = config.data_dir.join("jobs");
        private_dir(&config.data_dir)
            .and_then(|()| private_dir(&jobs_dir))
            .map_err(|e| format!("Could not prepare {}: {e}", jobs_dir.display()))?;
        let found = read_jobs(&jobs_dir, config.keep)?;
        let scheduler = Scheduler::new(
            &config.gpus.iter().map(|g| g.id).collect::<Vec<_>>(),
            config.cpu_jobs,
        );
        let workers = config.cpu_jobs.max(1) + config.gpus.len();
        let shared = Arc::new(Shared {
            config,
            state: Mutex::new(State {
                jobs: found.jobs,
                queue: found.queue,
                scheduler,
                cancel: HashSet::new(),
            }),
            wake: Condvar::new(),
            counter: AtomicU64::new(found.highest),
            uploading: AtomicUsize::new(0),
        });
        if run_workers {
            shared.start_workers(workers);
        }
        Ok(shared)
    }

    /// `count` workers, and one thread that clears away the jobs that have been kept long enough.
    fn start_workers(self: &Arc<Shared>, count: usize) {
        for _ in 0..count {
            let shared = Arc::clone(self);
            std::thread::spawn(move || shared.work());
        }
        let sweeper = Arc::clone(self);
        std::thread::spawn(move || loop {
            sweeper.sweep();
            std::thread::sleep(SWEEP);
        });
    }

    fn lock(&self) -> MutexGuard<'_, State> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Room in memory for an upload of `bytes`, or `None` when the uploads under way use it all: the sender is asked to
    /// try again, so that a few big uploads cannot fill the memory of the workstation. An upload is counted three times
    /// over, for the body, the container made from it and the sequences read out of that, and two of the largest fit at
    /// once.
    pub fn reserve_upload(&self, bytes: usize) -> Option<UploadRoom<'_>> {
        let bytes = bytes.saturating_mul(UPLOAD_COPIES);
        let most = self
            .config
            .max_upload
            .saturating_mul(2 * UPLOAD_COPIES)
            .max(MIN_UPLOAD_ROOM);
        let mut taken = self.uploading.load(Ordering::SeqCst);
        loop {
            let after = taken.checked_add(bytes).filter(|&a| a <= most)?;
            match self
                .uploading
                .compare_exchange(taken, after, Ordering::SeqCst, Ordering::SeqCst)
            {
                Ok(_) => {
                    return Some(UploadRoom {
                        room: &self.uploading,
                        bytes,
                    })
                }
                Err(seen) => taken = seen,
            }
        }
    }

    pub fn job_dir(&self, id: &str) -> PathBuf {
        self.config.data_dir.join("jobs").join(id)
    }

    /// The bytes of upload room in use now (each upload counted as it is held in memory): for tests that wait for it.
    #[cfg(test)]
    pub fn uploading_bytes(&self) -> usize {
        self.uploading.load(Ordering::SeqCst)
    }

    /// Put a job in the table as it is: for tests that stand in for a worker.
    #[cfg(test)]
    pub fn replace_for_tests(&self, job: JobInfo) {
        self.lock().jobs.insert(job.id.clone(), job);
    }

    fn save(&self, job: &JobInfo) {
        // a status that cannot be written (the disk is full, say) would be forgotten at the next start: say so where the
        // person who runs the server looks
        if let Err(e) = write_atomic(
            &self.job_dir(&job.id).join(META_FILE),
            job.to_line().as_bytes(),
        ) {
            eprintln!("Could not record the state of job {}: {e}", job.id);
        }
    }

    /// Why no more jobs are taken: too many wait, or too many are kept altogether.
    fn why_full(&self, st: &State) -> Option<String> {
        if st.queue.len() >= self.config.max_queued {
            return Some(format!(
                "The queue is full ({} jobs are waiting).",
                self.config.max_queued
            ));
        }
        if st.jobs.len() >= self.config.max_kept {
            return Some(format!(
                "The workstation keeps {} jobs already: remove some of the finished ones.",
                self.config.max_kept
            ));
        }
        None
    }

    /// Put a request in the queue; the id it will be known by.
    pub fn submit(
        &self,
        kind: Kind,
        label: &str,
        request: &[u8],
        needs_gpu: bool,
    ) -> Result<String, String> {
        let id = format!(
            "j{:x}-{:x}",
            now(),
            self.counter.fetch_add(1, Ordering::SeqCst) + 1
        );
        // when there is no room, say so before the request is written to the disk (and again below, for two requests that
        // came together)
        if let Some(why) = self.why_full(&self.lock()) {
            return Err(why);
        }
        let dir = self.job_dir(&id);
        let job = JobInfo {
            id: id.clone(),
            kind,
            label: one_line(label),
            status: Status::Queued,
            created: now(),
            started: 0,
            finished: 0,
            gpu: None,
            needs_gpu,
            error: String::new(),
        };
        private_dir(&dir)
            .and_then(|()| write_atomic(&dir.join(REQUEST_FILE), request))
            .and_then(|()| write_atomic(&dir.join(META_FILE), job.to_line().as_bytes()))
            .map_err(|e| {
                let _ = std::fs::remove_dir_all(&dir);
                format!("The workstation could not store the job: {e}")
            })?;
        let mut st = self.lock();
        if let Some(why) = self.why_full(&st) {
            drop(st);
            let _ = std::fs::remove_dir_all(&dir);
            return Err(why);
        }
        st.jobs.insert(id.clone(), job);
        st.queue.push_back(id.clone());
        self.wake.notify_all();
        Ok(id)
    }

    /// Every job, oldest first.
    pub fn list(&self) -> Vec<JobInfo> {
        let mut all: Vec<JobInfo> = self.lock().jobs.values().cloned().collect();
        all.sort_by_key(|j| made_order(j.created, &j.id));
        all
    }

    pub fn get(&self, id: &str) -> Option<JobInfo> {
        self.lock().jobs.get(id).cloned()
    }

    /// The bytes of a finished job's result.
    pub fn result(&self, id: &str) -> Result<Vec<u8>, Lookup> {
        match self.get(id) {
            None => Err(Lookup::NotFound),
            Some(job) if job.status != Status::Done => Err(Lookup::NotReady(job.status)),
            Some(_) => std::fs::read(self.job_dir(id).join(RESULT_FILE))
                .map_err(|e| Lookup::Failed(format!("The result could not be read: {e}"))),
        }
    }

    /// Remove a job: one that waits leaves the queue, one that runs is stopped (it stays in the list as cancelled
    /// until it is removed again), one that is over is forgotten and its files deleted.
    pub fn delete(&self, id: &str) -> Result<(), Lookup> {
        let mut st = self.lock();
        let status = st.jobs.get(id).ok_or(Lookup::NotFound)?.status;
        if status == Status::Running {
            // the worker that owns the process sees this within a moment and stops it: only a worker that still has the
            // process to wait for may kill it, so a process number that has been handed to something else is never hit
            st.cancel.insert(id.to_string());
            return Ok(());
        }
        st.queue.retain(|q| q != id);
        st.jobs.remove(id);
        drop(st);
        let _ = std::fs::remove_dir_all(self.job_dir(id));
        Ok(())
    }

    pub fn overview(&self) -> Overview {
        let st = self.lock();
        let busy = |gpu: usize| {
            st.jobs
                .values()
                .find(|j| j.status == Status::Running && j.gpu == Some(gpu))
                .map(|j| j.id.clone())
        };
        Overview {
            running: st
                .jobs
                .values()
                .filter(|j| j.status == Status::Running)
                .count(),
            queued: st.queue.len(),
            gpus: self
                .config
                .gpus
                .iter()
                .map(|g| (g.clone(), busy(g.id)))
                .collect(),
        }
    }

    /// Forget finished jobs that have been kept long enough.
    pub fn sweep(&self) {
        let now = now();
        let old: Vec<String> = self
            .lock()
            .jobs
            .values()
            .filter(|j| j.expired(now, self.config.keep))
            .map(|j| j.id.clone())
            .collect();
        for id in old {
            let _ = self.delete(&id);
        }
    }

    /// A worker: take the next job that can run, run it, give back what it used.
    fn work(&self) {
        loop {
            let (id, slot) = self.next_job();
            let outcome = self.run_job(&id, slot);
            self.finish(&id, slot, outcome);
        }
    }

    /// Wait for a job that can start (a GPU job needs a free GPU, any other a free CPU place).
    fn next_job(&self) -> (String, Slot) {
        let mut st = self.lock();
        loop {
            let pick = st.queue.iter().position(|id| {
                st.jobs
                    .get(id)
                    .is_some_and(|j| st.scheduler.can_take(j.needs_gpu))
            });
            if let Some(slot_at) = pick {
                let id = st.queue.remove(slot_at).unwrap_or_default();
                let needs = st.jobs.get(&id).is_some_and(|j| j.needs_gpu);
                if let (Some(slot), Some(job)) = (st.scheduler.take(needs), st.jobs.get_mut(&id)) {
                    job.status = Status::Running;
                    job.started = now();
                    job.gpu = match slot {
                        Slot::Gpu(g) => Some(g),
                        Slot::Cpu => None,
                    };
                    let job = job.clone();
                    // the disk is written without holding the table, so a slow disk does not hold up the others
                    drop(st);
                    self.save(&job);
                    return (id, slot);
                }
            }
            st = self.wake.wait(st).unwrap_or_else(PoisonError::into_inner);
        }
    }

    /// Run a job's process and wait for it, stopping it if it is cancelled or runs too long.
    fn run_job(&self, id: &str, slot: Slot) -> Outcome {
        let dir = self.job_dir(id);
        let tmp = dir.join("tmp");
        let outcome = match self.start_process(&dir, &tmp, slot) {
            Ok(mut child) => self.watch(id, &mut child, &dir),
            Err(why) => Outcome::Failed(why),
        };
        let _ = std::fs::remove_dir_all(&tmp);
        outcome
    }

    /// Start `helixloom run-job` for the job in `dir`, in a process group of its own so that all of it can be stopped at once.
    fn start_process(&self, dir: &Path, tmp: &Path, slot: Slot) -> Result<Child, String> {
        std::fs::create_dir_all(tmp).map_err(|e| format!("Could not make a work folder: {e}"))?;
        let log =
            File::create(dir.join(LOG_FILE)).map_err(|e| format!("Could not open the log: {e}"))?;
        let err = log.try_clone().map_err(|e| e.to_string())?;
        let exe =
            std::env::current_exe().map_err(|e| format!("Could not find this program: {e}"))?;
        let mut cmd = Command::new(exe);
        cmd.arg("run-job")
            .arg(dir)
            .arg("--fold-seconds")
            .arg(self.config.job_limit.as_secs().to_string())
            // the job ends by itself if this process goes (it is told which process that is, so that it is not fooled
            // by a server that went before the job had looked)
            .arg("--server-pid")
            .arg(std::process::id().to_string());
        if let Some(words) = &self.config.fold_command {
            cmd.arg("--").args(words);
        }
        // the job sees the GPU it was given and no other, keeps its temporary files in its own folder, and is never
        // handed the token (a server started with HELIXLOOM_TOKEN has it in its environment). The GPUs are numbered as
        // `nvidia-smi` numbers them (by place on the bus), not in CUDA's own order (fastest first), so that the number
        // on the status page and in `--gpus` is the card the job gets.
        cmd.env(
            "CUDA_VISIBLE_DEVICES",
            match slot {
                Slot::Gpu(g) => g.to_string(),
                Slot::Cpu => String::new(),
            },
        )
        .env("CUDA_DEVICE_ORDER", "PCI_BUS_ID")
        .env_remove("HELIXLOOM_TOKEN")
        .env("TMPDIR", tmp)
        .env("TEMP", tmp)
        .env("TMP", tmp)
        .stdin(Stdio::null())
        .stdout(log)
        .stderr(err);
        #[cfg(unix)]
        {
            use std::os::unix::process::CommandExt;
            cmd.process_group(0);
        }
        cmd.spawn()
            .map_err(|e| format!("Could not start the job's process: {e}"))
    }

    /// Wait for a job's process, stopping it if the job is cancelled or runs longer than it may.
    fn watch(&self, id: &str, child: &mut Child, dir: &Path) -> Outcome {
        let started = Instant::now();
        loop {
            if self.lock().cancel.contains(id) {
                stop_process(child);
                return Outcome::Cancelled;
            }
            match child.try_wait() {
                Ok(Some(status)) => {
                    // a program the job started (MAFFT, a folding program) may have outlived it, if it crashed or was killed:
                    // it would go on holding the GPU that is about to be given to the next job
                    kill_tree(child.id());
                    return ended(dir, status);
                }
                Ok(None) if started.elapsed() > self.config.job_limit => {
                    stop_process(child);
                    return Outcome::Failed(format!(
                        "The job ran longer than {} hours and was stopped.",
                        self.config.job_limit.as_secs() / 3600
                    ));
                }
                Ok(None) => std::thread::sleep(WATCH),
                Err(e) => return Outcome::Failed(format!("Could not wait for the job: {e}")),
            }
        }
    }

    fn finish(&self, id: &str, slot: Slot, outcome: Outcome) {
        // what was sent is not needed any more (the result, or the reason, is what is kept), and it may be hundreds of
        // megabytes. It goes before the job is shown as over, so that nobody sees "done" with the request still on the disk.
        let _ = std::fs::remove_file(self.job_dir(id).join(REQUEST_FILE));
        let mut st = self.lock();
        st.scheduler.give_back(slot);
        // a job somebody asked to stop is cancelled, unless it had already finished its work
        let outcome = if st.cancel.remove(id) && !matches!(outcome, Outcome::Done) {
            Outcome::Cancelled
        } else {
            outcome
        };
        let over = st.jobs.get_mut(id).map(|job| {
            job.finished = now();
            match outcome {
                Outcome::Done => job.status = Status::Done,
                Outcome::Cancelled => job.status = Status::Cancelled,
                Outcome::Failed(why) => {
                    job.status = Status::Failed;
                    job.error = one_line(&why);
                }
            }
            job.clone()
        });
        drop(st);
        if let Some(job) = over {
            self.save(&job);
        }
        self.wake.notify_all();
    }
}

#[cfg(test)]
mod tests;
