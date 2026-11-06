//! The server: it listens for connections, checks the token before it reads a body, answers the few requests the app
//! makes and hands jobs to the store. Nothing a request says is used as a path, a command or a program: a job is a
//! kind from a short list, its parts are checked when it arrives, and it runs in a process the server starts itself.
//!
//! The messages are not encrypted. Reach the server over an SSH tunnel or a private network (see DEPLOY.md), and never
//! put its port on the open internet.

use super::http::{self, Failure, Within, IO_TIMEOUT};
use super::jobs::{valid_id, Config, JobInfo, Lookup, Shared, Status};
use super::tasks::{Kind, Request};
use super::wire::{meta, Container};
use super::PROTOCOL;
use std::collections::HashMap;
use std::net::{IpAddr, SocketAddr, TcpListener, TcpStream};
use std::sync::{Arc, Mutex, PoisonError};
use std::time::Duration;

/// Connections handled at once; more are told to try again.
const MAX_CONNECTIONS: usize = 64;
/// Connections handled at once from one address, so that one computer cannot take every place. (Through an SSH tunnel
/// every connection comes from the workstation itself, so this is also what all the computers of one person share.)
const MAX_PER_PEER: usize = 8;
/// How long a connection may take to send its whole header. A real client sends it in one piece.
const HEAD_TIMEOUT: Duration = Duration::from_secs(5);
/// Requests other than a job upload carry no body; this is room for a small one anyway.
const SMALL_BODY: usize = 64 * 1024;
/// The slowest upload that is waited for, in bytes a second (the largest upload, at this rate, takes about three hours).
const MIN_UPLOAD_RATE: usize = 50 * 1024;
const CONTAINER: &str = "application/x-helixloom";

/// Who is connected: how many connections there are from each address, in all and one by one.
#[derive(Default)]
struct Places {
    from: Mutex<HashMap<IpAddr, usize>>,
}

/// A place taken by a connection; given back when this goes out of scope, however the connection ends.
struct Place {
    places: Arc<Places>,
    address: IpAddr,
}

impl Places {
    /// A place for a connection from `address`, or `None` when there are too many connections, or too many from there.
    fn take(self: &Arc<Self>, address: IpAddr) -> Option<Place> {
        let mut from = self.from.lock().unwrap_or_else(PoisonError::into_inner);
        let all: usize = from.values().sum();
        let here = from.get(&address).copied().unwrap_or(0);
        if all >= MAX_CONNECTIONS || here >= MAX_PER_PEER {
            return None;
        }
        from.insert(address, here + 1);
        Some(Place {
            places: Arc::clone(self),
            address,
        })
    }
}

impl Drop for Place {
    fn drop(&mut self) {
        let mut from = self
            .places
            .from
            .lock()
            .unwrap_or_else(PoisonError::into_inner);
        match from.get_mut(&self.address) {
            Some(n) if *n > 1 => *n -= 1,
            _ => {
                from.remove(&self.address);
            }
        }
    }
}

pub struct Server {
    pub addr: SocketAddr,
    pub shared: Arc<Shared>,
}

/// An answer: a status, what kind of body, and the body.
pub struct Reply {
    pub status: u16,
    pub content_type: &'static str,
    pub body: Vec<u8>,
}

impl Reply {
    fn text(status: u16, message: &str) -> Reply {
        Reply {
            status,
            content_type: "text/plain; charset=utf-8",
            body: message.as_bytes().to_vec(),
        }
    }

    fn container(status: u16, c: &Container) -> Reply {
        Reply {
            status,
            content_type: CONTAINER,
            body: c.encode(),
        }
    }
}

/// Whether two tokens are the same, taking as long as the longer one whatever the answer.
fn same_token(given: &str, expected: &str) -> bool {
    let (a, b) = (given.as_bytes(), expected.as_bytes());
    let mut diff = a.len() ^ b.len();
    for i in 0..a.len().max(b.len()) {
        diff |= usize::from(a.get(i).copied().unwrap_or(0) ^ b.get(i).copied().unwrap_or(0));
    }
    diff == 0
}

/// The token in an `Authorization: Bearer <token>` header.
fn bearer(header: Option<&str>) -> &str {
    header
        .and_then(|h| h.split_once(' '))
        .filter(|(scheme, _)| scheme.eq_ignore_ascii_case("bearer"))
        .map_or("", |(_, token)| token.trim())
}

fn status_reply(shared: &Shared) -> Reply {
    let c = &shared.config;
    let o = shared.overview();
    let tools: String = c.tools.iter().map(|(n, v)| format!("{n}\t{v}\n")).collect();
    let gpus: String = o
        .gpus
        .iter()
        .map(|(g, job)| {
            format!(
                "{}\t{}\t{}\n",
                g.id,
                super::one_line(&g.name),
                job.as_deref().unwrap_or("")
            )
        })
        .collect();
    let facts = [
        ("name", super::one_line(&c.name)),
        ("version", env!("CARGO_PKG_VERSION").to_string()),
        ("protocol", PROTOCOL.to_string()),
        ("running", o.running.to_string()),
        ("queued", o.queued.to_string()),
        ("cpu_jobs", c.cpu_jobs.to_string()),
        ("max_upload", c.max_upload.to_string()),
        ("keep_hours", (c.keep.as_secs() / 3600).to_string()),
        ("fold", u8::from(c.fold_command.is_some()).to_string()),
    ];
    let meta: String = facts.iter().map(|(k, v)| format!("{k}={v}\n")).collect();
    Reply::container(
        200,
        &Container::new()
            .with("meta", meta)
            .with("gpus", gpus)
            .with("tools", tools),
    )
}

fn job_lines(jobs: &[JobInfo]) -> String {
    jobs.iter().map(|j| format!("{}\n", j.to_line())).collect()
}

fn submit(shared: &Shared, body: &[u8]) -> Reply {
    let bad = |why: String| Reply::text(400, &why);
    let container = match Container::decode(body) {
        Ok(c) => c,
        Err(why) => return bad(why),
    };
    let job = match container.require("job") {
        Ok(t) => meta(t),
        Err(why) => return bad(why),
    };
    let Some(kind) = job.get("kind").and_then(|k| Kind::parse(k)) else {
        return bad("The job's kind is not one this server runs.".into());
    };
    let request = match Request::decode(kind, &container) {
        Ok(r) => r,
        Err(why) => return bad(why),
    };
    let needs_gpu = request.needs_gpu();
    if needs_gpu && shared.config.gpus.is_empty() {
        return bad("This workstation has no GPU set up for that job (see --gpus).".into());
    }
    if matches!(request, Request::Fold { local: true, .. }) && shared.config.fold_command.is_none()
    {
        return bad(
            "This workstation has no folding program configured (see --fold-command).".into(),
        );
    }
    match shared.submit(
        kind,
        job.get("label").copied().unwrap_or(""),
        body,
        needs_gpu,
    ) {
        Ok(id) => Reply::container(201, &Container::new().with("id", id)),
        Err(why) => Reply::text(503, &why),
    }
}

fn lookup_failure(why: &Lookup) -> Reply {
    match why {
        Lookup::NotFound => Reply::text(404, "There is no such job."),
        Lookup::NotReady(Status::Failed) => {
            Reply::text(409, "The job failed, so it has no result.")
        }
        Lookup::NotReady(Status::Cancelled) => {
            Reply::text(409, "The job was stopped, so it has no result.")
        }
        Lookup::NotReady(s) => Reply::text(
            409,
            &format!("The job is {}: its result is not ready.", s.name()),
        ),
        Lookup::Failed(m) => Reply::text(500, m),
    }
}

/// Answer a request that has been accepted (its token checked, its body read).
fn route(shared: &Shared, method: &str, target: &str, body: &[u8]) -> Reply {
    let path = target.split_once('?').map_or(target, |(p, _)| p);
    let parts: Vec<&str> = path.trim_matches('/').split('/').collect();
    match (method, parts.as_slice()) {
        ("GET", ["v1", "status"]) => status_reply(shared),
        ("GET", ["v1", "jobs"]) => Reply::container(
            200,
            &Container::new().with("jobs", job_lines(&shared.list())),
        ),
        ("POST", ["v1", "jobs"]) => submit(shared, body),
        ("GET", ["v1", "jobs", id]) if valid_id(id) => shared.get(id).map_or_else(
            || lookup_failure(&Lookup::NotFound),
            |j| Reply::container(200, &Container::new().with("job", j.to_line())),
        ),
        ("GET", ["v1", "jobs", id, "result"]) if valid_id(id) => match shared.result(id) {
            Ok(bytes) => Reply {
                status: 200,
                content_type: CONTAINER,
                body: bytes,
            },
            Err(why) => lookup_failure(&why),
        },
        ("DELETE", ["v1", "jobs", id]) if valid_id(id) => match shared.delete(id) {
            Ok(()) => Reply {
                status: 204,
                content_type: "text/plain",
                body: Vec::new(),
            },
            Err(why) => lookup_failure(&why),
        },
        (_, ["v1", "status"] | ["v1", "jobs"]) => {
            Reply::text(405, "That method is not allowed here.")
        }
        (_, ["v1", "jobs", id] | ["v1", "jobs", id, "result"]) if valid_id(id) => {
            Reply::text(405, "That method is not allowed here.")
        }
        (_, ["v1", "jobs", ..]) => Reply::text(404, "There is no such job."),
        _ => Reply::text(404, "There is nothing at that address."),
    }
}

/// Read one request from a connection, answer it, and close.
fn serve_connection(stream: &mut TcpStream, shared: &Shared, head_time: Duration) {
    let _ = stream.set_write_timeout(Some(IO_TIMEOUT));
    // the whole header within the time, however slowly it is sent: a stranger cannot hold a place by dribbling
    let head = http::read_head(&mut Within::new(stream, head_time));
    let reply = match head {
        Err(Failure::Closed(_)) => return,
        Err(why) => Reply::text(why.status(), &why.message()),
        Ok(head) if !same_token(bearer(head.header("authorization")), &shared.config.token) => {
            Reply::text(401, "The token is missing or wrong.")
        }
        Ok(head) => match take_in(stream, shared, head) {
            Some(reply) => reply,
            None => return,
        },
    };
    let _ = http::respond(stream, reply.status, reply.content_type, &reply.body);
    http::finish(stream);
}

/// How long an upload of `length` bytes may take in all: a minute, and as long as it needs at the slowest rate that is
/// waited for. Without that a sender that gives a byte now and then could keep its share of the memory for ever.
fn upload_time(length: usize) -> Duration {
    Duration::from_secs(60 + (length / MIN_UPLOAD_RATE) as u64)
}

/// Read the body of a request that has the right token and answer the request; `None` when the connection is gone and
/// there is nobody to answer.
fn take_in(stream: &mut TcpStream, shared: &Shared, head: http::Head) -> Option<Reply> {
    // only a job upload may be large
    let upload = head.method == "POST" && head.target.split('?').next() == Some("/v1/jobs");
    let limit = if upload {
        shared.config.max_upload
    } else {
        SMALL_BODY
    };
    let length = match http::announced_length(&head) {
        Ok(n) => n,
        Err(why) => return Some(Reply::text(why.status(), &why.message())),
    };
    // a body is held in memory, so only so many bytes of bodies are taken in at once (a body that is too large is
    // refused without being read, and takes no room)
    let _room = if length <= limit {
        match shared.reserve_upload(length) {
            Some(room) => Some(room),
            None => {
                return Some(Reply::text(
                    503,
                    "The workstation is taking in other uploads: try again in a moment.",
                ))
            }
        }
    } else {
        None
    };
    // the whole body within the time that its length allows, however slowly it comes
    let mut body = Within::new(stream, upload_time(length)).idle(IO_TIMEOUT);
    match http::read_body(&mut body, head, limit) {
        Ok(request) => Some(route(
            shared,
            &request.head.method,
            &request.head.target,
            &request.body,
        )),
        Err(Failure::Closed(_)) => None,
        Err(why) => Some(Reply::text(why.status(), &why.message())),
    }
}

/// Open the store, bind `listen` and start answering. With `run_workers` false jobs queue but never run.
pub fn start(config: Config, listen: SocketAddr, run_workers: bool) -> Result<Server, String> {
    if config.token.len() < super::MIN_TOKEN {
        return Err(format!(
            "The token must be at least {} characters (make one with `helixloom token`).",
            super::MIN_TOKEN
        ));
    }
    let listener =
        TcpListener::bind(listen).map_err(|e| format!("Cannot listen on {listen}: {e}"))?;
    let addr = listener.local_addr().map_err(|e| e.to_string())?;
    let shared = Shared::open(config, run_workers)?;
    // ponytail: a thread for each connection, at most MAX_CONNECTIONS; an async server only if many people share one
    // workstation
    let places = Arc::new(Places::default());
    let accept = Arc::clone(&shared);
    std::thread::spawn(move || {
        for stream in listener.incoming() {
            let Ok(mut stream) = stream else {
                // an error that goes on (too many open files, say) would otherwise be tried again as fast as the processor
                // allows
                std::thread::sleep(Duration::from_millis(50));
                continue;
            };
            let place = stream
                .peer_addr()
                .ok()
                .and_then(|peer| places.take(peer.ip()));
            let Some(place) = place else {
                http::turn_away(&mut stream, 503, "Too many connections: try again soon.");
                continue;
            };
            let shared = Arc::clone(&accept);
            // if no thread can be started the closure, and the place in it, is dropped and the connection closes
            let _ = std::thread::Builder::new()
                .name("connection".into())
                .spawn(move || {
                    // given back when the thread ends, even if it panics
                    let _place = place;
                    serve_connection(&mut stream, &shared, HEAD_TIMEOUT);
                });
        }
    });
    Ok(Server { addr, shared })
}

#[cfg(test)]
mod tests;
