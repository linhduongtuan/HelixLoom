//! The app's side of the workstation: send a job, wait for it, fetch what came back, look at the list of jobs.

use super::http::{self, Allowed, Reply, RESULT_REPLY, SMALL_REPLY};
use super::jobs::{valid_id, JobInfo, Status};
use super::one_line;
use super::tasks::{Kind, Request, Response};
use super::wire::{meta, Container};
use super::DEFAULT_PORT;
use std::time::Duration;

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
        Kind::Align | Kind::Pipeline => RESULT_REPLY,
        Kind::Tree | Kind::PhyloCnn | Kind::Fold | Kind::Species => Allowed {
            bytes: SMALL_RESULT,
            time: Duration::from_secs(600),
        },
    }
}

/// What was typed as the workstation's address, as `host:port`. The port is 7878 when none is given. `https://` is
/// refused: the server does not speak TLS, and the way to a safe connection is an SSH tunnel or a private network.
pub fn parse_address(text: &str) -> Result<String, String> {
    let text = text.trim();
    if text.starts_with("https://") {
        return Err("The workstation server does not use TLS: use http://, through an SSH tunnel or a private network.".into());
    }
    let host_port = text
        .strip_prefix("http://")
        .unwrap_or(text)
        .trim_end_matches('/');
    if host_port.is_empty() || host_port.contains(['/', ' ', '\t', '@', '?', '#']) {
        return Err("The workstation's address should look like host:7878.".into());
    }
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
        Ok(format!("{host_port}:{DEFAULT_PORT}"))
    }
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
}

/// The address only: the token is a password, and a log line or a failed test must never print it.
impl std::fmt::Debug for Remote {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Remote")
            .field("address", &self.address)
            .field("token", &"<hidden>")
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
        Ok(Remote {
            address: parse_address(address)?,
            token: token.to_string(),
        })
    }

    pub fn address(&self) -> &str {
        &self.address
    }

    /// The token, for saving the connection.
    pub fn token(&self) -> &str {
        &self.token
    }

    /// One request. When it fails: the status the server answered with (0 when there was no answer) and why, in words.
    fn attempt(
        &self,
        method: &str,
        target: &str,
        body: &[u8],
        allowed: Allowed,
    ) -> Result<Reply, (u16, String)> {
        let reply = http::call(&self.address, method, target, &self.token, body, allowed)
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

    fn call(
        &self,
        method: &str,
        target: &str,
        body: &[u8],
        allowed: Allowed,
    ) -> Result<Reply, String> {
        self.attempt(method, target, body, allowed)
            .map_err(|(_, why)| why)
    }

    fn container(&self, method: &str, target: &str, body: &[u8]) -> Result<Container, String> {
        Container::decode(&self.call(method, target, body, SMALL_REPLY)?.body)
    }

    pub fn status(&self) -> Result<ServerStatus, String> {
        ServerStatus::parse(&self.container("GET", "/v1/status", b"")?)
    }

    /// Every job on the workstation, oldest first.
    pub fn jobs(&self) -> Result<Vec<JobInfo>, String> {
        let c = self.container("GET", "/v1/jobs", b"")?;
        Ok(c.require("jobs")?
            .lines()
            .take(MAX_LISTED)
            .filter_map(JobInfo::from_line)
            .collect())
    }

    /// One job's line. When that fails: the status the server answered with (0 when there was no answer) and why.
    pub fn job(&self, id: &str) -> Result<JobInfo, (u16, String)> {
        let reply = self.attempt("GET", &format!("/v1/jobs/{id}"), b"", SMALL_REPLY)?;
        Container::decode(&reply.body)
            .and_then(|c| {
                JobInfo::from_line(c.require("job")?).ok_or_else(|| {
                    "The workstation sent a job this program cannot read.".to_string()
                })
            })
            .map_err(|why| (0, why))
    }

    /// Send a request to be run, with a label for the list of jobs; the job's id.
    pub fn submit(&self, request: &Request, label: &str) -> Result<String, String> {
        let kind = request.kind();
        let job = format!("kind={}\nlabel={}\n", kind.name(), super::one_line(label));
        let body = request.encode().with("job", job).encode();
        let reply = self.container("POST", "/v1/jobs", &body)?;
        // the id goes into the address of every later request: it is only taken if it is one this server makes
        let id = reply.require("id")?;
        if !valid_id(id) {
            return Err("The workstation sent a job number this program does not accept.".into());
        }
        Ok(id.to_string())
    }

    /// A finished job's result.
    pub fn result(&self, kind: Kind, id: &str) -> Result<Response, String> {
        let reply = self.call(
            "GET",
            &format!("/v1/jobs/{id}/result"),
            b"",
            reply_for(kind),
        )?;
        Response::decode(kind, &Container::decode(&reply.body)?)
    }

    /// Remove a job, or stop it if it is running.
    pub fn delete(&self, id: &str) -> Result<(), String> {
        self.call("DELETE", &format!("/v1/jobs/{id}"), b"", SMALL_REPLY)
            .map(|_| ())
    }

    /// Send a job, wait until it is over and give its result. If the connection drops while waiting, the job goes on
    /// on the workstation and can be fetched from its list.
    pub fn run(&self, request: &Request, label: &str) -> Result<Response, String> {
        let id = self.submit(request, label)?;
        let mut failures = 0;
        loop {
            // ponytail: one new connection per question, a second apart (the server answers one request per connection);
            // a held-open connection or a long poll would save the set-up if a great many jobs were watched at once
            std::thread::sleep(POLL);
            match self.job(&id) {
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
                            return self.result(request.kind(), &id).map_err(|why| {
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
        use std::io::Read;
        let listener = std::net::TcpListener::bind("127.0.0.1:0").expect("a port");
        let address = listener.local_addr().expect("an address").to_string();
        std::thread::spawn(move || {
            let (mut stream, _) = listener.accept().expect("a connection");
            let mut request = [0u8; 4096];
            let _ = stream.read(&mut request);
            send(&mut stream);
        });
        address
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
        let why = http::call(&at, "GET", "/", "token", b"", allowed)
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
