//! What the server is set up to do: the options of `helixloom serve` as one value, which the rest of the server reads.

use crate::remote::gpus::Gpu;
use crate::remote::one_line;
use std::path::PathBuf;
use std::process::Command;
use std::time::Duration;

#[derive(Clone)]
pub struct Config {
    pub name: String,
    pub data_dir: PathBuf,
    /// Anyone who sends this token may use the server.
    pub token: String,
    pub gpus: Vec<Gpu>,
    /// How many jobs that use no GPU run at once.
    pub cpu_jobs: usize,
    /// The largest request, in bytes.
    pub max_upload: usize,
    /// How long a finished job's result is kept.
    pub keep: Duration,
    /// How long one job may run.
    pub job_limit: Duration,
    /// The most jobs waiting at once.
    pub max_queued: usize,
    /// The most jobs kept in all, finished ones included, so that a client that never removes its jobs cannot fill the disk.
    pub max_kept: usize,
    /// The program that folds proteins on a GPU (its words), if the workstation has one.
    pub fold_command: Option<Vec<String>>,
    /// How many GPUs that program uses for one structure: a fold job is given this many at once (and no more, whatever the
    /// app asks for).
    pub fold_gpus: usize,
    /// The programs found on this machine and their versions, for the status page.
    pub tools: Vec<(String, String)>,
}

/// Everything but the token, which a log or a failed test must never print.
impl std::fmt::Debug for Config {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Config")
            .field("name", &self.name)
            .field("data_dir", &self.data_dir)
            .field("token", &"<hidden>")
            .field("gpus", &self.gpus)
            .field("cpu_jobs", &self.cpu_jobs)
            .field("max_upload", &self.max_upload)
            .field("keep", &self.keep)
            .field("job_limit", &self.job_limit)
            .field("max_queued", &self.max_queued)
            .field("max_kept", &self.max_kept)
            .field("fold_command", &self.fold_command)
            .field("fold_gpus", &self.fold_gpus)
            .field("tools", &self.tools)
            .finish()
    }
}

impl Config {
    pub fn new(token: String, data_dir: PathBuf) -> Config {
        Config {
            name: machine_name(),
            data_dir,
            token,
            gpus: Vec::new(),
            cpu_jobs: 2,
            max_upload: 512 << 20,
            keep: Duration::from_secs(7 * 24 * 3600),
            job_limit: Duration::from_secs(24 * 3600),
            max_queued: 1000,
            max_kept: 2000,
            fold_command: None,
            fold_gpus: 1,
            tools: Vec::new(),
        }
    }
}

/// This machine's name, for the status page: the environment's, else what `hostname` says.
fn machine_name() -> String {
    std::env::var("HOSTNAME")
        .or_else(|_| std::env::var("COMPUTERNAME"))
        .ok()
        .or_else(|| {
            Command::new("hostname")
                .output()
                .ok()
                .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        })
        .map(|n| one_line(&n))
        .filter(|n| !n.is_empty())
        .unwrap_or_else(|| "workstation".into())
}
