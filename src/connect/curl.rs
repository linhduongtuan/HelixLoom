//! The fetcher that reaches the web with `curl`.

use super::{Fetcher, Reply};
use std::path::PathBuf;
use std::process::Command;
use std::sync::Mutex;
use std::time::Duration;
use web_time::Instant;

/// The longest a request may take, and the largest answer taken (the largest PDB entries are tens of megabytes).
const SECONDS: &str = "60";
const MOST_BYTES: &str = "200000000";
/// NCBI asks for no more than three requests a second without a key of one's own.
const NCBI_GAP: Duration = Duration::from_millis(350);
const NCBI_HOST: &str = "https://eutils.ncbi.nlm.nih.gov/";

static LAST_NCBI: Mutex<Option<Instant>> = Mutex::new(None);

/// How long to wait before a request, given when the last one to the same service was made.
pub fn wait_needed(last: Option<Instant>, now: Instant, gap: Duration) -> Duration {
    last.map_or(Duration::ZERO, |l| {
        gap.saturating_sub(now.saturating_duration_since(l))
    })
}

pub struct Curl {
    program: PathBuf,
    /// Only the tests, which ask a server on this computer, allow a plain `http` address.
    allow_http: bool,
}

impl Curl {
    pub fn new(program: PathBuf) -> Curl {
        Curl {
            program,
            allow_http: false,
        }
    }

    #[cfg(test)]
    pub fn for_tests(program: PathBuf) -> Curl {
        Curl {
            program,
            allow_http: true,
        }
    }
}

impl Fetcher for Curl {
    fn get(&self, url: &str) -> Result<Reply, String> {
        if url.starts_with(NCBI_HOST) {
            // one request at a time to NCBI, spaced, whichever thread asks
            let mut last = LAST_NCBI.lock().unwrap_or_else(|e| e.into_inner());
            std::thread::sleep(wait_needed(*last, Instant::now(), NCBI_GAP));
            *last = Some(Instant::now());
        }
        let protocols = if self.allow_http {
            "=http,https"
        } else {
            "=https"
        };
        let output = Command::new(&self.program)
            // `-q` must come first: it leaves out a `.curlrc` of the person's, which could change anything below. `--globoff`:
            // no `{}` or `[]` expansion of the address. At most five redirects.
            .arg("-q")
            .args([
                "--silent",
                "--show-error",
                "--globoff",
                "--location",
                "--max-redirs",
                "5",
                "--max-time",
                SECONDS,
            ])
            .args([
                "--max-filesize",
                MOST_BYTES,
                "--proto",
                protocols,
                "--proto-redir",
                protocols,
            ])
            .args(["--write-out", "\n%{http_code}", "--user-agent", "HelixLoom"])
            .arg("--")
            .arg(url)
            .output()
            .map_err(|e| format!("could not start curl: {e}"))?;
        if !output.status.success() {
            let why = String::from_utf8_lossy(&output.stderr);
            let why = why.trim();
            return Err(if why.is_empty() {
                "curl gave no reason".to_string()
            } else {
                why.to_string()
            });
        }
        let answer = String::from_utf8_lossy(&output.stdout);
        let (body, status) = answer.rsplit_once('\n').ok_or("curl gave no status")?;
        let status = status
            .trim()
            .parse()
            .map_err(|_| format!("curl gave a status that is not a number: {status}"))?;
        Ok(Reply {
            status,
            body: body.to_string(),
        })
    }
}
