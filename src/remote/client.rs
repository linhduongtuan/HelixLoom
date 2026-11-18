//! The app's side of the workstation: send a job, wait for it, fetch what came back, look at the list of jobs.

use super::http::{Allowed, Reply, RESULT_REPLY, SMALL_REPLY};
use super::jobs::{valid_id, JobInfo, Status};
use super::one_line;
use super::settings::GpusPerJob;
use super::tasks::{Kind, Request, Response};
use super::wire::{meta, Container};
use super::DEFAULT_PORT;
use std::time::Duration;

#[cfg(not(target_arch = "wasm32"))]
use super::block::block_on;

/// How often a waiting job is asked about.
const POLL: Duration = Duration::from_secs(1);
/// So many answers in a row that fail (a network that dropped out) end the waiting; the job itself goes on.
const MAX_FAILED_POLLS: usize = 30;
/// The most jobs, GPUs and programs that are taken from what a workstation says about itself.
const MAX_LISTED: usize = 10_000;
const MAX_GPUS: usize = 64;

/// The most that is taken as the result of a tree, PhyloCNN estimates or a structure: far more than any of them is. (A tree
/// is read into a node for every few characters, and a structure is worked over residue against residue, so a reply
/// many times larger would take the program's memory or time.)
const SMALL_RESULT: usize = 8 << 20;

/// How large a result of each kind may be taken: only an alignment, and a pipeline that holds one, can be large.
fn reply_for(kind: Kind) -> Allowed {
    match kind {
        // a pipeline file may end with the alignments of many files
        Kind::Align | Kind::Pipeline | Kind::Flow => RESULT_REPLY,
        Kind::Tree | Kind::PhyloCnn | Kind::Fold | Kind::Species | Kind::Motif | Kind::Domains => {
            Allowed {
                bytes: SMALL_RESULT,
                time: Duration::from_secs(600),
            }
        }
    }
}

/// `host_port` as it is, when it has a port, else with `default_port` added (`None`: left without one).
fn with_port(host_port: &str, default_port: Option<u16>) -> Result<String, String> {
    // a port, unless the text is a bare IPv6 address in brackets or a name
    let has_port = match host_port.rsplit_once(':') {
        Some((host, port)) if !host.is_empty() && (!host.contains(':') || host.ends_with(']')) => {
            port.parse::<u16>()
                .map_err(|_| format!("\"{port}\" is not a port number."))?;
            true
        }
        _ => false,
    };
    if has_port {
        Ok(host_port.to_string())
    } else if host_port.contains(':') && !host_port.starts_with('[') {
        Err("An IPv6 address needs brackets: [::1]:7878.".into())
    } else {
        Ok(default_port.map_or_else(|| host_port.to_string(), |p| format!("{host_port}:{p}")))
    }
}

/// The text after the scheme (if any) of what was typed as an address: no path, no user, no query.
fn host_part(text: &str) -> Result<&str, String> {
    let host_port = text
        .strip_prefix("http://")
        .or_else(|| text.strip_prefix("https://"))
        .unwrap_or(text)
        .trim_end_matches('/');
    if host_port.is_empty() || host_port.contains(['/', ' ', '\t', '@', '?', '#']) {
        return Err("The workstation's address should look like host:7878.".into());
    }
    Ok(host_port)
}

/// What was typed as the workstation's address, as `host:port`. The port is 7878 when none is given. `https://` is
/// refused: the server does not speak TLS, and the way to a safe connection is an SSH tunnel or a private network.
pub fn parse_address(text: &str) -> Result<String, String> {
    let text = text.trim();
    if text.starts_with("https://") {
        return Err("The workstation server does not use TLS: use http://, through an SSH tunnel or a private network.".into());
    }
    with_port(host_part(text)?, Some(DEFAULT_PORT))
}

/// The same for a page, which reaches the workstation through the browser's `fetch`: the address keeps its scheme
/// (`http://host:7878`), and `https://` is allowed because the browser, not this program, does the TLS (a server behind a
/// proxy that has a certificate). A page served over https cannot call a plain http address: the browser refuses. A bare
/// `host` is `http://host:7878`; an address that has a scheme is a URL, and without a port it has the scheme's own (so the
/// page's own address, which a browser writes without the default port, is taken as it is).
pub fn parse_web_address(text: &str) -> Result<String, String> {
    let text = text.trim();
    let secure = text.starts_with("https://");
    let has_scheme = secure || text.starts_with("http://");
    let default_port = (!has_scheme).then_some(DEFAULT_PORT);
    let scheme = if secure { "https" } else { "http" };
    Ok(format!(
        "{scheme}://{}",
        with_port(host_part(text)?, default_port)?
    ))
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GpuState {
    pub id: usize,
    pub name: String,
    /// The job using it, if any.
    pub job: Option<String>,
}

/// What a workstation says about itself.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ServerStatus {
    pub name: String,
    pub version: String,
    pub protocol: u32,
    pub running: usize,
    pub queued: usize,
    pub cpu_jobs: usize,
    pub max_upload: usize,
    pub keep_hours: u64,
    /// It has a program that folds proteins on a GPU.
    pub fold: bool,
    pub gpus: Vec<GpuState>,
    pub tools: Vec<(String, String)>,
}

impl ServerStatus {
    pub fn parse(c: &Container) -> Result<ServerStatus, String> {
        let m = meta(c.require("meta")?);
        // what a server says about itself is shown as short lines of text, whatever it sent
        let text = |k: &str| one_line(m.get(k).copied().unwrap_or(""));
        let number = |k: &str| -> Result<usize, String> {
            m.get(k)
                .and_then(|v| v.parse().ok())
                .ok_or_else(|| format!("The workstation's status has no \"{k}\"."))
        };
        let gpus = c
            .text("gpus")
            .unwrap_or("")
            .lines()
            .take(MAX_GPUS)
            .filter_map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                let [id, name, job] = f.as_slice() else {
                    return None;
                };
                Some(GpuState {
                    id: id.parse().ok()?,
                    name: one_line(name),
                    job: valid_id(job).then(|| (*job).to_string()),
                })
            })
            .collect();
        let tools = c
            .text("tools")
            .unwrap_or("")
            .lines()
            .take(MAX_GPUS)
            .filter_map(|l| l.split_once('\t'))
            .map(|(n, v)| (one_line(n), one_line(v)))
            .collect();
        Ok(ServerStatus {
            name: text("name"),
            version: text("version"),
            protocol: u32::try_from(number("protocol")?).unwrap_or(0),
            running: number("running")?,
            queued: number("queued")?,
            cpu_jobs: number("cpu_jobs")?,
            max_upload: number("max_upload")?,
            keep_hours: number("keep_hours")? as u64,
            fold: number("fold")? == 1,
            gpus,
            tools,
        })
    }
}

#[derive(Clone, PartialEq, Eq)]
pub struct Remote {
    address: String,
    token: String,
    /// How many GPUs one job is to have at once, when it can use several.
    gpus: GpusPerJob,
}

/// The address only: the token is a password, and a log line or a failed test must never print it.
impl std::fmt::Debug for Remote {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Remote")
            .field("address", &self.address)
            .field("token", &"<hidden>")
            .field("gpus", &self.gpus)
            .finish()
    }
}

impl Remote {
    /// A connection to the workstation at `address` (as typed), with its token.
    pub fn new(address: &str, token: &str) -> Result<Remote, String> {
        let token = token.trim();
        if token.is_empty() {
            return Err(
                "There is no token: paste the one made on the workstation with `helixloom token`."
                    .into(),
            );
        }
        // what `helixloom token` makes is letters and digits: anything else is a mistake in pasting, and a browser's message about
        // a header it refuses may quote it
        if !token.bytes().all(|b| b.is_ascii_graphic()) {
            return Err("The token should be the letters and digits that `helixloom token` makes, with no spaces or other characters.".into());
        }
        Ok(Remote {
            address: if cfg!(target_arch = "wasm32") {
                parse_web_address(address)?
            } else {
                parse_address(address)?
            },
            token: token.to_string(),
            gpus: GpusPerJob::default(),
        })
    }

    /// The same connection, asking for `gpus` GPUs for each job that can use several.
    pub fn with_gpus(self, gpus: GpusPerJob) -> Remote {
        Remote { gpus, ..self }
    }

    pub fn address(&self) -> &str {
        &self.address
    }

    /// The token, for saving the connection.
    pub fn token(&self) -> &str {
        &self.token
    }

    /// One request. When it fails: the status the server answered with (0 when there was no answer) and why, in words.
    async fn attempt(
        &self,
        method: &str,
        target: &str,
        body: &[u8],
        allowed: Allowed,
    ) -> Result<Reply, (u16, String)> {
        let reply = self
            .send(method, target, body, allowed)
            .await
            .map_err(|why| (0, why))?;
        match reply.status {
            200..=299 => Ok(reply),
            401 => Err((401, "The workstation refused the token: check that it is the one in the file the server was started with.".into())),
            status => {
                // what a server says is shown as one short line, whatever it sent
                let said = super::one_line(&String::from_utf8_lossy(&reply.body));
                Err((
                    status,
                    if said.is_empty() {
                        format!("The workstation answered {status}.")
                    } else {
                        said
                    },
                ))
            }
        }
    }

    /// The request itself: on a plain connection in a program, through the browser's `fetch` in a page.
    async fn send(
        &self,
        method: &str,
        target: &str,
        body: &[u8],
        allowed: Allowed,
    ) -> Result<Reply, String> {
        #[cfg(not(target_arch = "wasm32"))]
        {
            super::http::call(&self.address, method, target, &self.token, body, allowed)
        }
        #[cfg(target_arch = "wasm32")]
        {
            super::fetch::call(&self.address, method, target, &self.token, body, allowed).await
        }
    }

    async fn call(
        &self,
        method: &str,
        target: &str,
        body: &[u8],
        allowed: Allowed,
    ) -> Result<Reply, String> {
        self.attempt(method, target, body, allowed)
            .await
            .map_err(|(_, why)| why)
    }

    async fn container(
        &self,
        method: &str,
        target: &str,
        body: &[u8],
    ) -> Result<Container, String> {
        Container::decode(&self.call(method, target, body, SMALL_REPLY).await?.body)
    }

    pub async fn status_async(&self) -> Result<ServerStatus, String> {
        ServerStatus::parse(&self.container("GET", "/v1/status", b"").await?)
    }

    /// Every job on the workstation, oldest first.
    pub async fn jobs_async(&self) -> Result<Vec<JobInfo>, String> {
        let c = self.container("GET", "/v1/jobs", b"").await?;
        Ok(c.require("jobs")?
            .lines()
            .take(MAX_LISTED)
            .filter_map(JobInfo::from_line)
            .collect())
    }

    /// One job's line. When that fails: the status the server answered with (0 when there was no answer) and why.
    pub async fn job_async(&self, id: &str) -> Result<JobInfo, (u16, String)> {
        let reply = self
            .attempt("GET", &format!("/v1/jobs/{id}"), b"", SMALL_REPLY)
            .await?;
        Container::decode(&reply.body)
            .and_then(|c| {
                JobInfo::from_line(c.require("job")?).ok_or_else(|| {
                    "The workstation sent a job this program cannot read.".to_string()
                })
            })
            .map_err(|why| (0, why))
    }

    /// Send a request to be run, with a label for the list of jobs; the job's id.
    pub async fn submit_async(&self, request: &Request, label: &str) -> Result<String, String> {
        let kind = request.kind();
        let job = format!(
            "kind={}\nlabel={}\ngpus={}\n",
            kind.name(),
            super::one_line(label),
            self.gpus.wire()
        );
        let body = request.encode().with("job", job).encode();
        let reply = self.container("POST", "/v1/jobs", &body).await?;
        // the id goes into the address of every later request: it is only taken if it is one this server makes
        let id = reply.require("id")?;
        if !valid_id(id) {
            return Err("The workstation sent a job number this program does not accept.".into());
        }
        Ok(id.to_string())
    }

    /// A finished job's result.
    pub async fn result_async(&self, kind: Kind, id: &str) -> Result<Response, String> {
        let reply = self
            .call(
                "GET",
                &format!("/v1/jobs/{id}/result"),
                b"",
                reply_for(kind),
            )
            .await?;
        Response::decode(kind, &Container::decode(&reply.body)?)
    }

    /// Remove a job, or stop it if it is running.
    pub async fn delete_async(&self, id: &str) -> Result<(), String> {
        self.call("DELETE", &format!("/v1/jobs/{id}"), b"", SMALL_REPLY)
            .await
            .map(|_| ())
    }

    /// Send a job, wait until it is over and give its result. If the connection drops while waiting, the job goes on
    /// on the workstation and can be fetched from its list.
    pub async fn run_async(&self, request: &Request, label: &str) -> Result<Response, String> {
        let id = self.submit_async(request, label).await?;
        let mut failures = 0;
        loop {
            // ponytail: one new connection per question, a second apart (the server answers one request per connection);
            // a held-open connection or a long poll would save the set-up if a great many jobs were watched at once
            sleep(POLL).await;
            match self.job_async(&id).await {
                // the job is gone (removed from the Workstation page, say) or the token no longer fits: no use asking again
                Err((404, _)) => {
                    return Err("The job is no longer on the workstation: it was removed there before it was over.".into())
                }
                Err((401, why)) => return Err(why),
                Ok(job) => {
                    failures = 0;
                    match job.status {
                        Status::Done => {
                            // the work is done and kept there: say so if only the fetching went wrong
                            return self.result_async(request.kind(), &id).await.map_err(|why| {
                                format!(
                                    "The job finished on the workstation, but its result could not be fetched ({why}). \
                                     Open it from the Workstation page."
                                )
                            });
                        }
                        Status::Failed => {
                            return Err(if job.error.is_empty() {
                                "The job failed on the workstation.".to_string()
                            } else {
                                job.error
                            })
                        }
                        Status::Cancelled => {
                            return Err("The job was stopped on the workstation.".into())
                        }
                        Status::Queued | Status::Running => {}
                    }
                }
                Err((_, why)) => {
                    failures += 1;
                    if failures >= MAX_FAILED_POLLS {
                        return Err(format!(
                            "Lost touch with the workstation ({why}). The job ({id}) goes on there: fetch it from the Workstation page."
                        ));
                    }
                }
            }
        }
    }
}

/// A wait that does not stop the thread's other work in a page (a timer of the browser), and waits in a program.
#[cfg(not(target_arch = "wasm32"))]
async fn sleep(time: Duration) {
    std::thread::sleep(time);
}
#[cfg(target_arch = "wasm32")]
use super::fetch::sleep;

/// A program has threads to wait on: the same calls, which block until they are answered. (A page cannot block, and asks
/// with the `_async` ones.)
#[cfg(not(target_arch = "wasm32"))]
impl Remote {
    pub fn status(&self) -> Result<ServerStatus, String> {
        block_on(self.status_async())
    }
}

/// The rest of the blocking calls, for the tests that stand in for a program.
#[cfg(test)]
impl Remote {
    pub fn jobs(&self) -> Result<Vec<JobInfo>, String> {
        block_on(self.jobs_async())
    }

    pub fn job(&self, id: &str) -> Result<JobInfo, (u16, String)> {
        block_on(self.job_async(id))
    }

    pub fn submit(&self, request: &Request, label: &str) -> Result<String, String> {
        block_on(self.submit_async(request, label))
    }

    pub fn result(&self, kind: Kind, id: &str) -> Result<Response, String> {
        block_on(self.result_async(kind, id))
    }

    pub fn delete(&self, id: &str) -> Result<(), String> {
        block_on(self.delete_async(id))
    }

    pub fn run(&self, request: &Request, label: &str) -> Result<Response, String> {
        block_on(self.run_async(request, label))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_address_is_made_into_host_and_port() {
        for (typed, host_port) in [
            ("workstation.local", "workstation.local:7878"),
            ("  workstation.local:9000 ", "workstation.local:9000"),
            ("http://10.0.0.5", "10.0.0.5:7878"),
            ("http://10.0.0.5:8080/", "10.0.0.5:8080"),
            ("localhost:7878", "localhost:7878"),
            ("[::1]:7878", "[::1]:7878"),
            ("[::1]", "[::1]:7878"),
        ] {
            assert_eq!(parse_address(typed).as_deref(), Ok(host_port), "{typed}");
        }
        for bad in [
            "",
            "   ",
            "host:notaport",
            "host:99999",
            "a b:1",
            "http://host/path",
            "user@host",
            "::1",
            "host:7878?x",
        ] {
            assert!(parse_address(bad).is_err(), "{bad:?}");
        }
        assert!(parse_address("https://host").unwrap_err().contains("TLS"));
    }

    #[test]
    fn a_page_keeps_the_scheme_of_the_address_it_is_given() {
        for (typed, url) in [
            ("workstation.local", "http://workstation.local:7878"),
            ("  localhost:7878 ", "http://localhost:7878"),
            ("http://10.0.0.5:8080/", "http://10.0.0.5:8080"),
            // an address with a scheme is a URL: no port typed means the scheme's own (a page served on port 80 writes its
            // origin without one)
            ("http://10.0.0.5", "http://10.0.0.5"),
            // the browser does the TLS, so a name behind a proxy with a certificate is reached as it is
            ("https://ws.example.org", "https://ws.example.org"),
            (
                "https://ws.example.org:8443/",
                "https://ws.example.org:8443",
            ),
            ("[::1]", "http://[::1]:7878"),
            ("http://[::1]:9", "http://[::1]:9"),
        ] {
            assert_eq!(parse_web_address(typed).as_deref(), Ok(url), "{typed}");
        }
        for bad in [
            "",
            "  ",
            "http://",
            "https://",
            "ftp://host",
            "https://host/path",
            "user@host",
            "host:notaport",
            "https://host:99999",
            "::1",
        ] {
            assert!(parse_web_address(bad).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn a_token_is_plain_text_or_it_is_refused_before_it_is_sent() {
        for good in ["abc123", "a-token_of.sixteen+chars/="] {
            assert!(Remote::new("host", good).is_ok(), "{good}");
        }
        // a pasting mistake, and what a header cannot hold: a browser's message about it may quote it
        for bad in [
            "to ken",
            "tok\nen",
            "tok\ren: x",
            "tökén",
            "tok\ten",
            "to\u{0}ken",
        ] {
            let why = Remote::new("host", bad).unwrap_err();
            assert!(why.contains("letters and digits"), "{bad:?}: {why}");
            assert!(!why.contains(bad), "the message must not repeat the token");
        }
    }

    #[test]
    fn a_connection_needs_an_address_and_a_token() {
        assert!(Remote::new("host", "").is_err());
        assert!(Remote::new("host", "   ").is_err());
        assert!(Remote::new("", "token").is_err());
        let r = Remote::new("host", " token ").expect("a connection");
        assert_eq!(r.address(), "host:7878");
    }

    #[test]
    fn a_status_is_read_from_what_the_server_sends() {
        let c = Container::new()
            .with(
                "meta",
                "name=lab-gpu\nversion=0.1.0\nprotocol=1\nrunning=1\nqueued=2\ncpu_jobs=2\nmax_upload=1000\nkeep_hours=168\nfold=1\n",
            )
            .with("gpus", "0\tNVIDIA A\tj1-1\n1\tNVIDIA B\t\n")
            .with("tools", "MAFFT\tv7\nIQ-TREE\t3.1\n");
        let s = ServerStatus::parse(&c).expect("a status");
        assert_eq!(
            (
                s.name.as_str(),
                s.protocol,
                s.running,
                s.queued,
                s.keep_hours
            ),
            ("lab-gpu", 1, 1, 2, 168)
        );
        assert!(s.fold);
        assert_eq!(
            s.gpus[0],
            GpuState {
                id: 0,
                name: "NVIDIA A".into(),
                job: Some("j1-1".into())
            }
        );
        assert_eq!(s.gpus[1].job, None);
        assert_eq!(
            s.tools,
            [
                ("MAFFT".to_string(), "v7".to_string()),
                ("IQ-TREE".to_string(), "3.1".to_string())
            ]
        );
        assert!(ServerStatus::parse(&Container::new().with("meta", "name=x\n")).is_err());
    }

    /// A server on this computer that reads one request and then does what `send` says: what a server that is not what it
    /// should be does. Its address.
    fn hostile(send: impl FnOnce(&mut std::net::TcpStream) + Send + 'static) -> String {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").expect("a port");
        let address = listener.local_addr().expect("an address").to_string();
        std::thread::spawn(move || {
            let (mut stream, _) = listener.accept().expect("a connection");
            read_request(&mut stream);
            send(&mut stream);
        });
        address
    }

    /// Read one request whole: its header, and as many bytes of body as the header says. A server that answers and closes
    /// with some of the request still unread turns its answer into a reset on the way (the header and the body are sent as
    /// two writes, so one read may get only the header), and the test would then see "connection reset", not the answer.
    fn read_request(stream: &mut std::net::TcpStream) {
        use std::io::Read;
        let mut got = Vec::new();
        let mut chunk = [0u8; 4096];
        loop {
            let n = stream.read(&mut chunk).unwrap_or(0);
            if n == 0 {
                return;
            }
            got.extend_from_slice(&chunk[..n]);
            let Some(end) = got.windows(4).position(|w| w == b"\r\n\r\n") else {
                continue;
            };
            let length = String::from_utf8_lossy(&got[..end])
                .lines()
                .filter_map(|l| l.split_once(':'))
                .find(|(name, _)| name.eq_ignore_ascii_case("content-length"))
                .and_then(|(_, value)| value.trim().parse::<usize>().ok())
                .unwrap_or(0);
            if got.len() >= end + 4 + length {
                return;
            }
        }
    }

    fn remote(at: &str) -> Remote {
        Remote::new(at, "a-token-of-sixteen-chars").expect("a connection")
    }

    #[test]
    fn a_reply_larger_than_a_result_can_be_is_refused_not_read_into_memory() {
        use std::io::Write;
        let at = hostile(|s| {
            let _ = s.write_all(b"HTTP/1.1 200 OK\r\nContent-Length: 20000000\r\n\r\n");
            let piece = vec![b'x'; 1 << 16];
            for _ in 0..400 {
                if s.write_all(&piece).is_err() {
                    break;
                }
            }
        });
        let why = remote(&at)
            .result(Kind::Tree, "j1-1")
            .err()
            .expect("refused");
        assert!(
            why.contains("larger than this program takes (8 MB)"),
            "{why}"
        );
    }

    #[test]
    fn a_server_that_sends_a_byte_now_and_then_cannot_keep_the_program_waiting() {
        use std::io::Write;
        let at = hostile(|s| {
            let _ = s.write_all(b"HTTP/1.1 200 OK\r\nContent-Length: 100\r\n\r\n");
            for _ in 0..100 {
                if s.write_all(b"x").is_err() {
                    break;
                }
                std::thread::sleep(Duration::from_millis(100));
            }
        });
        let allowed = Allowed {
            bytes: 1 << 20,
            time: Duration::from_millis(400),
        };
        let started = std::time::Instant::now();
        let why = crate::remote::http::call(&at, "GET", "/", "token", b"", allowed)
            .err()
            .expect("too slow");
        assert!(why.contains("took too long"), "{why}");
        assert!(
            started.elapsed() < Duration::from_secs(5),
            "{:?}",
            started.elapsed()
        );
    }

    #[test]
    fn what_a_server_says_about_itself_and_its_jobs_is_cleaned_before_it_is_shown_or_used() {
        // a name with a line break and a terminal escape, a job number that is a path, an overlong label
        let dirty = Container::new()
            .with("meta", "name=lab\u{1b}[31m\nx\nversion=1\nprotocol=1\nrunning=0\nqueued=0\ncpu_jobs=2\nmax_upload=1\nkeep_hours=1\nfold=0\n")
            .with("gpus", "0\tGPU\u{1b}[2J\tj1-1\n1\tB\t../../etc\n")
            .with("tools", "MAFFT\t7\u{7}\n");
        let s = ServerStatus::parse(&dirty).expect("a status");
        assert!(!s.name.contains(char::is_control) && !s.gpus[0].name.contains(char::is_control));
        assert!(!s.tools[0].1.contains(char::is_control));
        assert_eq!(s.gpus[0].job.as_deref(), Some("j1-1"));
        assert_eq!(
            s.gpus[1].job, None,
            "a job number that is not one of the server's is dropped"
        );
        let line = |id: &str, label: &str| format!("{id}\talign\tdone\t1\t2\t3\t\t0\t{label}\t");
        assert!(JobInfo::from_line(&line("../../etc/passwd", "x")).is_none());
        assert!(JobInfo::from_line(&line("j1-1", "ok")).is_some());
        let long = JobInfo::from_line(&line("j1-1", &"y".repeat(5000))).expect("a job");
        assert_eq!(long.label.chars().count(), 200);
    }

    #[test]
    fn a_job_number_that_is_not_one_the_server_makes_is_not_taken_from_a_reply() {
        use std::io::Write;
        let at = hostile(|s| {
            let body = Container::new().with("id", "../../etc/passwd").encode();
            let head = format!(
                "HTTP/1.1 201 Created\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                body.len()
            );
            let _ = s.write_all(head.as_bytes());
            let _ = s.write_all(&body);
        });
        let request = Request::Fold {
            sequence: b"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ".to_vec(),
            local: false,
        };
        let why = remote(&at).submit(&request, "x").expect_err("refused");
        assert!(why.contains("job number"), "{why}");
    }

    /// A server on this computer that takes one request whole, answers that the job is queued as `j1-1`, and hands what it was
    /// sent to the receiver: its address.
    fn recording() -> (String, std::sync::mpsc::Receiver<Vec<u8>>) {
        use std::io::{Read, Write};
        let listener = std::net::TcpListener::bind("127.0.0.1:0").expect("a port");
        let address = listener.local_addr().expect("an address").to_string();
        let (send, receive) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let (mut stream, _) = listener.accept().expect("a connection");
            let mut got = Vec::new();
            let mut chunk = [0u8; 4096];
            // until the header is in and as many bytes of the body as it says
            let whole = |got: &[u8]| {
                let text = String::from_utf8_lossy(got);
                let Some(end) = text.find("\r\n\r\n") else {
                    return false;
                };
                let length: usize = text[..end]
                    .lines()
                    .find_map(|l| l.strip_prefix("Content-Length: "))
                    .and_then(|n| n.trim().parse().ok())
                    .unwrap_or(0);
                got.len() >= end + 4 + length
            };
            while !whole(&got) {
                let n = stream.read(&mut chunk).unwrap_or(0);
                if n == 0 {
                    break;
                }
                got.extend_from_slice(&chunk[..n]);
            }
            let body = Container::new().with("id", "j1-1").encode();
            let head = format!(
                "HTTP/1.1 201 Created\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                body.len()
            );
            let _ = stream.write_all(head.as_bytes());
            let _ = stream.write_all(&body);
            let _ = send.send(got);
        });
        (address, receive)
    }

    #[test]
    fn a_job_is_sent_asking_for_the_gpus_the_connection_was_set_to() {
        let request = Request::Fold {
            sequence: b"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ".to_vec(),
            local: true,
        };
        // a connection that has not been told asks for one, as before jobs could have several
        for (gpus, line) in [
            (None, "gpus=1\n"),
            (Some(GpusPerJob::Count(1)), "gpus=1\n"),
            (Some(GpusPerJob::Count(3)), "gpus=3\n"),
            (Some(GpusPerJob::All), "gpus=0\n"),
        ] {
            let (at, sent) = recording();
            let connection = gpus.map_or_else(|| remote(&at), |g| remote(&at).with_gpus(g));
            assert_eq!(connection.submit(&request, "x").expect("queued"), "j1-1");
            let body = sent
                .recv_timeout(Duration::from_secs(10))
                .expect("what was sent");
            let text = String::from_utf8_lossy(&body);
            assert!(text.contains(line), "{gpus:?}: {text}");
            assert!(text.contains("kind=fold\n"), "{text}");
        }
    }

    #[test]
    fn two_connections_to_one_workstation_are_the_same_only_when_they_ask_for_the_same_gpus() {
        let a = Remote::new("h:1", "token-token-token-1").expect("a connection");
        assert_eq!(a.clone(), a);
        assert_ne!(a.clone().with_gpus(GpusPerJob::All), a);
        assert_eq!(
            a.clone().with_gpus(GpusPerJob::Count(1)),
            a,
            "one is what a new connection asks for"
        );
        // the token is still never printed, and what the connection asks for is
        let shown = format!("{:?}", a.with_gpus(GpusPerJob::Count(2)));
        assert!(
            shown.contains("<hidden>") && shown.contains("Count(2)"),
            "{shown}"
        );
        assert!(!shown.contains("token-token"), "{shown}");
    }

    #[test]
    fn what_a_server_says_when_it_refuses_is_shown_as_one_short_line() {
        use std::io::Write;
        let at = hostile(|s| {
            let body = format!("bad\u{1b}[31m\nthing {}", "z".repeat(1000));
            let head = format!(
                "HTTP/1.1 400 Bad Request\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
                body.len()
            );
            let _ = s.write_all(head.as_bytes());
            let _ = s.write_all(body.as_bytes());
        });
        let why = remote(&at).status().expect_err("refused");
        assert!(
            !why.contains(char::is_control) && why.chars().count() <= 200,
            "{why:?}"
        );
    }
}
