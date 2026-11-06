//! What the list of jobs shows about a job: its status and its line of text (which is also what `meta.txt` holds), and the
//! numbers that make an id.

use crate::remote::one_line;
use crate::remote::tasks::Kind;
use std::time::Duration;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Queued,
    Running,
    Done,
    Failed,
    Cancelled,
}

impl Status {
    pub fn name(self) -> &'static str {
        match self {
            Status::Queued => "queued",
            Status::Running => "running",
            Status::Done => "done",
            Status::Failed => "failed",
            Status::Cancelled => "cancelled",
        }
    }

    pub fn parse(name: &str) -> Option<Status> {
        [
            Status::Queued,
            Status::Running,
            Status::Done,
            Status::Failed,
            Status::Cancelled,
        ]
        .into_iter()
        .find(|s| s.name() == name)
    }

    /// Finished, one way or another: nothing more will happen to the job.
    pub fn is_over(self) -> bool {
        matches!(self, Status::Done | Status::Failed | Status::Cancelled)
    }
}

/// What the list of jobs shows about one job (and what `meta.txt` holds).
#[derive(Clone, Debug, PartialEq)]
pub struct JobInfo {
    pub id: String,
    pub kind: Kind,
    pub label: String,
    pub status: Status,
    /// Seconds since 1970; 0 where it has not happened yet.
    pub created: u64,
    pub started: u64,
    pub finished: u64,
    pub gpu: Option<usize>,
    pub needs_gpu: bool,
    pub error: String,
}

impl JobInfo {
    /// The job as one line of tab-separated fields.
    pub fn to_line(&self) -> String {
        [
            self.id.clone(),
            self.kind.name().to_string(),
            self.status.name().to_string(),
            self.created.to_string(),
            self.started.to_string(),
            self.finished.to_string(),
            self.gpu.map_or(String::new(), |g| g.to_string()),
            u8::from(self.needs_gpu).to_string(),
            one_line(&self.label),
            one_line(&self.error),
        ]
        .join("\t")
    }

    pub fn from_line(line: &str) -> Option<JobInfo> {
        let f: Vec<&str> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
        let [id, kind, status, created, started, finished, gpu, needs_gpu, label, error] =
            f.as_slice()
        else {
            return None;
        };
        // a line may come from a server this program does not know: its id goes into addresses, its text onto the screen
        if !valid_id(id) {
            return None;
        }
        Some(JobInfo {
            id: (*id).to_string(),
            kind: Kind::parse(kind)?,
            status: Status::parse(status)?,
            created: created.parse().ok()?,
            started: started.parse().ok()?,
            finished: finished.parse().ok()?,
            gpu: if gpu.is_empty() {
                None
            } else {
                Some(gpu.parse().ok()?)
            },
            needs_gpu: *needs_gpu == "1",
            label: one_line(label),
            error: one_line(error),
        })
    }

    /// Whether a finished job has been kept as long as `keep` allows.
    pub(super) fn expired(&self, now: u64, keep: Duration) -> bool {
        self.status.is_over() && now.saturating_sub(self.finished) > keep.as_secs()
    }
}

/// The order in which jobs were made: by time, then by the counter in the id as a number (as text, `…-10` would come
/// before `…-9`, and a dozen jobs sent in one second would be listed out of order).
pub(super) fn made_order(created: u64, id: &str) -> (u64, u64) {
    (created, counter_of(id))
}

/// The counter in an id (the part after the dash, in hex).
pub(super) fn counter_of(id: &str) -> u64 {
    id.split_once('-')
        .and_then(|(_, c)| u64::from_str_radix(c, 16).ok())
        .unwrap_or(0)
}

/// An id is `j`, the creation time and a counter in hex: nothing a path could be made of.
pub fn valid_id(id: &str) -> bool {
    let Some(rest) = id.strip_prefix('j') else {
        return false;
    };
    let Some((time, count)) = rest.split_once('-') else {
        return false;
    };
    let hex = |s: &str, most: usize| {
        (1..=most).contains(&s.len())
            && s.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    };
    hex(time, 16) && hex(count, 16)
}
