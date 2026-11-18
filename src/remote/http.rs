//! A small HTTP/1.1 for the job server and its client, over plain TCP: one request and one reply on each connection
//! (`Connection: close`), bodies of a stated length only (no chunked transfer), strict limits on everything a stranger
//! can make the server read. The server reads the header first and decides whether to read the body at all, so someone
//! with no token never gets a large upload accepted.

use std::io::{ErrorKind, Read, Write};
use std::net::{Shutdown, TcpStream, ToSocketAddrs};
use std::time::{Duration, Instant};

/// The longest request header the server reads.
pub const HEADER_LIMIT: usize = 16 * 1024;
const MAX_HEADERS: usize = 64;
const MAX_TARGET: usize = 512;
const MAX_METHOD: usize = 8;
pub const CONNECT_TIMEOUT: Duration = Duration::from_secs(10);
/// How long one read or write may wait.
pub const IO_TIMEOUT: Duration = Duration::from_secs(60);
/// The slowest upload that is waited for, in bytes a second (the largest upload, at this rate, takes about three hours).
pub const MIN_UPLOAD_RATE: usize = 50 * 1024;
/// The most a client reads back in one reply: as much as the largest upload, since a result is not much larger than what
/// made it.
pub const MAX_RESPONSE: usize = 512 << 20;

/// What a reply may be: how many bytes of body, and how long it may take to arrive in all (a server that sends a byte now
/// and then cannot keep the client waiting for ever).
#[derive(Clone, Copy)]
pub struct Allowed {
    pub bytes: usize,
    pub time: Duration,
}

/// A status, a list of jobs, a job: small, and quick.
pub const SMALL_REPLY: Allowed = Allowed {
    bytes: 4 << 20,
    time: Duration::from_secs(30),
};
/// A result: as large as a result can be, and it may take a while over a slow link.
pub const RESULT_REPLY: Allowed = Allowed {
    bytes: MAX_RESPONSE,
    time: Duration::from_secs(3600),
};
/// After a reply that left a body unread, this much is read and thrown away so the reply is not lost to a reset...
const DRAIN_LIMIT: usize = 1 << 16;
/// ...for this long at most, however slowly it comes.
const DRAIN_TIME: Duration = Duration::from_secs(2);
/// The body buffer grows by what arrives, so a claimed length alone costs no memory.
const FIRST_BODY_ALLOC: usize = 1 << 20;

#[derive(Debug, PartialEq, Eq)]
pub enum Failure {
    /// The request is not valid HTTP, or breaks a limit: 400.
    Bad(String),
    /// The body is longer than the limit: 413.
    TooLarge,
    /// Something HTTP allows that this server does not: 501.
    Unsupported(String),
    /// The connection failed or ended; no reply can be sent.
    Closed(String),
}

impl Failure {
    pub fn status(&self) -> u16 {
        match self {
            Failure::Bad(_) | Failure::Closed(_) => 400,
            Failure::TooLarge => 413,
            Failure::Unsupported(_) => 501,
        }
    }

    pub fn message(&self) -> String {
        match self {
            Failure::Bad(m) | Failure::Unsupported(m) | Failure::Closed(m) => m.clone(),
            Failure::TooLarge => "The upload is larger than this server accepts.".into(),
        }
    }
}

/// A request's first line and headers; the body is not read yet.
#[derive(Debug)]
pub struct Head {
    pub method: String,
    pub target: String,
    headers: Vec<(String, String)>,
    /// Bytes of the body that arrived with the header.
    leftover: Vec<u8>,
}

impl Head {
    /// The value of header `name` (compared without regard to case).
    pub fn header(&self, name: &str) -> Option<&str> {
        self.headers
            .iter()
            .find(|(n, _)| n.eq_ignore_ascii_case(name))
            .map(|(_, v)| v.as_str())
    }
}

pub struct Request {
    pub head: Head,
    pub body: Vec<u8>,
}

fn find(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack.windows(needle.len()).position(|w| w == needle)
}

fn token_char(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b"!#$%&'*+-.^_`|~".contains(&b)
}

/// Read from the stream up to the blank line that ends a header: what has come (the header, and perhaps the start of the
/// body) and where the header ends.
fn read_to_blank_line<S: Read>(stream: &mut S) -> Result<(Vec<u8>, usize), Failure> {
    let mut buf: Vec<u8> = Vec::new();
    let mut chunk = [0u8; 4096];
    let end = loop {
        if let Some(p) = find(&buf, b"\r\n\r\n") {
            break p;
        }
        if buf.len() > HEADER_LIMIT {
            return Err(Failure::Bad("The request header is too long.".into()));
        }
        let n = stream
            .read(&mut chunk)
            .map_err(|e| Failure::Closed(e.to_string()))?;
        if n == 0 {
            return Err(Failure::Closed(
                "The connection ended inside the request.".into(),
            ));
        }
        buf.extend(&chunk[..n]);
    };
    if end > HEADER_LIMIT {
        return Err(Failure::Bad("The request header is too long.".into()));
    }
    Ok((buf, end))
}

/// Whether a first line is one this server reads: a method in capitals, a target of plain characters that starts at the
/// root, and HTTP/1.1 (or 1.0).
fn reads_request_line(method: &str, target: &str, version: &str) -> bool {
    !method.is_empty()
        && method.len() <= MAX_METHOD
        && method.bytes().all(|b| b.is_ascii_uppercase())
        && (version == "HTTP/1.1" || version == "HTTP/1.0")
        && target.starts_with('/')
        && target.len() <= MAX_TARGET
        && !target.bytes().any(|b| b.is_ascii_control() || b > 0x7e)
}

/// Read a request's header: the first line and the headers, checked against the limits.
pub fn read_head<S: Read>(stream: &mut S) -> Result<Head, Failure> {
    let (buf, end) = read_to_blank_line(stream)?;
    let text = std::str::from_utf8(&buf[..end])
        .map_err(|_| Failure::Bad("The request header is not text.".into()))?;
    let mut lines = text.split("\r\n");
    let first = lines.next().unwrap_or_default();
    let parts: Vec<&str> = first.split(' ').collect();
    let [method, target, version] = parts.as_slice() else {
        return Err(Failure::Bad(
            "The first line is not \"METHOD target HTTP/1.1\".".into(),
        ));
    };
    if !reads_request_line(method, target, version) {
        return Err(Failure::Bad(
            "The first line is not a request this server reads.".into(),
        ));
    }
    let mut headers = Vec::new();
    for line in lines {
        let (name, value) = line
            .split_once(':')
            .filter(|(n, _)| !n.is_empty() && n.bytes().all(token_char))
            .ok_or_else(|| Failure::Bad("A header line is not \"name: value\".".into()))?;
        if headers.len() >= MAX_HEADERS {
            return Err(Failure::Bad("The request has too many headers.".into()));
        }
        headers.push((name.to_string(), value.trim().to_string()));
    }
    Ok(Head {
        method: (*method).to_string(),
        target: (*target).to_string(),
        headers,
        leftover: buf[end + 4..].to_vec(),
    })
}

/// The length of the body that a header announces: none, or one clear number. A body sent in chunks is not accepted.
pub fn announced_length(head: &Head) -> Result<usize, Failure> {
    if head.header("transfer-encoding").is_some() {
        return Err(Failure::Unsupported(
            "Chunked uploads are not supported: send a Content-Length.".into(),
        ));
    }
    let lengths: Vec<&str> = head
        .headers
        .iter()
        .filter(|(n, _)| n.eq_ignore_ascii_case("content-length"))
        .map(|(_, v)| v.as_str())
        .collect();
    match lengths.as_slice() {
        [] => Ok(0),
        [one] if one.bytes().all(|b| b.is_ascii_digit()) && !one.is_empty() => one
            .parse::<usize>()
            .map_err(|_| Failure::Bad("The Content-Length is not a number.".into())),
        _ => Err(Failure::Bad(
            "The Content-Length is not one clear number.".into(),
        )),
    }
}

/// Read the body the header announces, as long as it is no longer than `limit`.
pub fn read_body<S: Read + Write>(
    stream: &mut S,
    mut head: Head,
    limit: usize,
) -> Result<Request, Failure> {
    let length = announced_length(&head)?;
    if length > limit {
        return Err(Failure::TooLarge);
    }
    let mut body = std::mem::take(&mut head.leftover);
    body.truncate(length);
    if body.len() < length {
        let expects = head
            .header("expect")
            .is_some_and(|e| e.eq_ignore_ascii_case("100-continue"));
        if expects {
            stream
                .write_all(b"HTTP/1.1 100 Continue\r\n\r\n")
                .and_then(|()| stream.flush())
                .map_err(|e| Failure::Closed(e.to_string()))?;
        }
        body.reserve(length.min(FIRST_BODY_ALLOC).saturating_sub(body.len()));
        let mut chunk = vec![0u8; 1 << 16];
        while body.len() < length {
            let want = chunk.len().min(length - body.len());
            let n = stream
                .read(&mut chunk[..want])
                .map_err(|e| Failure::Closed(e.to_string()))?;
            if n == 0 {
                return Err(Failure::Closed(
                    "The connection ended inside the body.".into(),
                ));
            }
            body.extend(&chunk[..n]);
        }
    }
    Ok(Request { head, body })
}

pub fn reason(status: u16) -> &'static str {
    match status {
        200 => "OK",
        201 => "Created",
        204 => "No Content",
        304 => "Not Modified",
        400 => "Bad Request",
        401 => "Unauthorized",
        404 => "Not Found",
        405 => "Method Not Allowed",
        409 => "Conflict",
        413 => "Payload Too Large",
        429 => "Too Many Requests",
        501 => "Not Implemented",
        503 => "Service Unavailable",
        _ => "Internal Server Error",
    }
}

/// Write a reply: the status, a body of the stated length, and the end of the connection.
pub fn respond<W: Write>(
    out: &mut W,
    status: u16,
    content_type: &str,
    body: &[u8],
) -> std::io::Result<()> {
    respond_with(out, status, content_type, &[], body)
}

/// The same with more headers. One that has the name of a usual header (`Cache-Control`) takes its place. A header that
/// holds a line break is left out: a value must never be able to end the header and start another reply.
pub fn respond_with<W: Write>(
    out: &mut W,
    status: u16,
    content_type: &str,
    extra: &[(&str, String)],
    body: &[u8],
) -> std::io::Result<()> {
    let extra: Vec<&(&str, String)> = extra
        .iter()
        .filter(|(n, v)| !n.contains(['\r', '\n']) && !v.contains(['\r', '\n']))
        .collect();
    write!(
        out,
        "HTTP/1.1 {status} {}\r\nContent-Type: {content_type}\r\n",
        reason(status)
    )?;
    // a "not modified" ends with its header, and has no length
    if status != 304 {
        write!(out, "Content-Length: {}\r\n", body.len())?;
    }
    out.write_all(b"Connection: close\r\n")?;
    if !extra
        .iter()
        .any(|(n, _)| n.eq_ignore_ascii_case("cache-control"))
    {
        out.write_all(b"Cache-Control: no-store\r\n")?;
    }
    out.write_all(b"X-Content-Type-Options: nosniff\r\n")?;
    for (name, value) in extra {
        write!(out, "{name}: {value}\r\n")?;
    }
    out.write_all(b"\r\n")?;
    out.write_all(body)?;
    out.flush()
}

/// A connection that can be read or written only until a deadline. A timeout on each read is not enough: a sender that gives one
/// byte just before every timeout would hold the connection, and the thread and place it takes, for hours.
pub struct Within<'a> {
    stream: &'a TcpStream,
    until: Instant,
    /// The longest one read may wait, whatever time is left.
    idle: Duration,
}

impl<'a> Within<'a> {
    pub fn new(stream: &'a TcpStream, time: Duration) -> Within<'a> {
        Within {
            stream,
            until: Instant::now() + time,
            idle: time,
        }
    }

    /// The same, but no one read may wait longer than `idle`.
    pub fn idle(mut self, idle: Duration) -> Within<'a> {
        self.idle = idle;
        self
    }
}

impl Read for Within<'_> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let left = self.until.saturating_duration_since(Instant::now());
        if left.is_zero() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::TimedOut,
                "The sender took too long.",
            ));
        }
        // setting the wait can fail on a connection the other side has already reset (macOS says "invalid argument"), and
        // what it sent before that may still be there to read: the read tells the truth about the connection
        let _ = self.stream.set_read_timeout(Some(left.min(self.idle)));
        match (&*self.stream).read(buf) {
            // a read that waited out its timeout says "would block" on some systems and "timed out" on others
            Err(e) if matches!(e.kind(), ErrorKind::WouldBlock | ErrorKind::TimedOut) => Err(
                std::io::Error::new(ErrorKind::TimedOut, "The sender took too long."),
            ),
            other => other,
        }
    }
}

/// Turn a connection away with a reply, without waiting for anything: what the other side has sent already is read and
/// thrown away (once, if there is any), so that closing the connection is less likely to turn the reply into a reset.
pub fn turn_away(stream: &mut TcpStream, status: u16, message: &str) {
    let _ = respond(stream, status, "text/plain", message.as_bytes());
    let _ = stream.shutdown(Shutdown::Write);
    if stream.set_nonblocking(true).is_ok() {
        let mut sink = [0u8; 4096];
        let _ = stream.read(&mut sink);
    }
}

impl Write for Within<'_> {
    /// Written within the same deadline: a receiver that takes a few bytes before every timeout would otherwise hold the
    /// connection, and the place it takes, for as long as the reply is long.
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let left = self.until.saturating_duration_since(Instant::now());
        if left.is_zero() {
            return Err(std::io::Error::new(
                ErrorKind::TimedOut,
                "The receiver took too long.",
            ));
        }
        let _ = self.stream.set_write_timeout(Some(left.min(self.idle)));
        match (&*self.stream).write(bytes) {
            Err(e) if matches!(e.kind(), ErrorKind::WouldBlock | ErrorKind::TimedOut) => Err(
                std::io::Error::new(ErrorKind::TimedOut, "The receiver took too long."),
            ),
            other => other,
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        (&*self.stream).flush()
    }
}

/// After a reply, read and throw away what the other side is still sending (a little of it, for a short time), so
/// closing the connection does not turn the reply into a reset on its way.
pub fn finish(stream: &mut TcpStream) {
    let _ = stream.shutdown(Shutdown::Write);
    let mut sink = [0u8; 4096];
    let mut total = 0;
    let mut other = Within::new(stream, DRAIN_TIME);
    while total < DRAIN_LIMIT {
        match other.read(&mut sink) {
            Ok(0) | Err(_) => break,
            Ok(n) => total += n,
        }
    }
}

pub struct Reply {
    pub status: u16,
    pub body: Vec<u8>,
}

/// A reply read from the bytes of a whole connection.
pub fn parse_reply(bytes: &[u8]) -> Result<Reply, String> {
    let end =
        find(bytes, b"\r\n\r\n").ok_or("The workstation's reply ended before its header did.")?;
    let head = std::str::from_utf8(&bytes[..end]).map_err(|_| "The reply header is not text.")?;
    let mut lines = head.split("\r\n");
    let status: u16 = lines
        .next()
        .and_then(|l| l.split(' ').nth(1))
        .and_then(|s| s.parse().ok())
        .ok_or("The reply does not start with an HTTP status.")?;
    let length = lines
        .filter_map(|l| l.split_once(':'))
        .find(|(n, _)| n.eq_ignore_ascii_case("content-length"))
        .and_then(|(_, v)| v.trim().parse::<usize>().ok());
    let mut body = bytes[end + 4..].to_vec();
    if let Some(n) = length {
        if body.len() < n {
            return Err("The workstation's reply was cut short.".into());
        }
        body.truncate(n);
    }
    Ok(Reply { status, body })
}

/// A connection to `address` (`host:port`): to the first of its addresses that answers.
fn connect(address: &str) -> Result<TcpStream, String> {
    // a page cannot open a TCP connection; it reaches the workstation through the browser's fetch (`fetch.rs`)
    if cfg!(target_arch = "wasm32") {
        return Err(
            "A page reaches the workstation through fetch, not through a TCP connection.".into(),
        );
    }
    let addrs = address
        .to_socket_addrs()
        .map_err(|e| format!("Cannot find {address}: {e}"))?;
    let mut last = format!("Cannot reach {address}.");
    for a in addrs {
        match TcpStream::connect_timeout(&a, CONNECT_TIMEOUT) {
            Ok(stream) => return Ok(stream),
            Err(e) => last = format!("Cannot reach {address}: {e}"),
        }
    }
    Err(last)
}

/// Send one request to `address` (`host:port`) and read the whole reply, as long as it is no larger than `allowed` says
/// and arrives in the time it says.
pub fn call(
    address: &str,
    method: &str,
    target: &str,
    token: &str,
    body: &[u8],
    allowed: Allowed,
) -> Result<Reply, String> {
    let mut stream = connect(address)?;
    let _ = stream.set_read_timeout(Some(IO_TIMEOUT));
    let _ = stream.set_write_timeout(Some(IO_TIMEOUT));
    let head = format!(
        "{method} {target} HTTP/1.1\r\nHost: {address}\r\nAuthorization: Bearer {token}\r\nContent-Type: application/x-helixloom\r\nContent-Length: {}\r\nConnection: close\r\n\r\n",
        body.len()
    );
    // a server may answer (a refusal, say) before it has read the whole body: read the answer even if the write failed
    let wrote = stream
        .write_all(head.as_bytes())
        .and_then(|()| stream.write_all(body))
        .and_then(|()| stream.flush());
    let mut got = Vec::new();
    let most = (allowed.bytes + HEADER_LIMIT) as u64;
    let read = Within::new(&stream, allowed.time)
        .idle(IO_TIMEOUT)
        .take(most)
        .read_to_end(&mut got);
    if got.len() as u64 >= most {
        return Err(format!(
            "The reply of {address} is larger than this program takes ({} MB).",
            allowed.bytes >> 20
        ));
    }
    match parse_reply(&got) {
        Ok(reply) => Ok(reply),
        Err(why) => Err(match (wrote, read) {
            (Err(e), _) | (_, Err(e)) => format!("The connection to {address} failed: {e}"),
            _ => why,
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// The text of a reply written by `respond_with`.
    fn written(status: u16, extra: &[(&str, String)], body: &[u8]) -> String {
        let mut out = Vec::new();
        respond_with(&mut out, status, "text/plain", extra, body).expect("written");
        String::from_utf8(out).expect("text")
    }

    #[test]
    fn extra_headers_are_added_and_one_that_ends_the_header_is_left_out() {
        let plain = written(200, &[], b"hi");
        assert!(plain
            .starts_with("HTTP/1.1 200 OK\r\nContent-Type: text/plain\r\nContent-Length: 2\r\n"));
        assert!(plain.contains("Cache-Control: no-store\r\n") && plain.ends_with("\r\n\r\nhi"));
        let extra = written(200, &[("ETag", "\"abc\"".into())], b"hi");
        assert!(extra.contains("ETag: \"abc\"\r\n"));
        // a value or a name that holds a line break would end the header and start another reply
        for (name, value) in [
            ("X-A", "1\r\nSet-Cookie: stolen=1"),
            ("X-A", "1\nX-B: 2"),
            ("X-A\r\nX-B", "2"),
        ] {
            let split = written(200, &[(name, value.to_string())], b"hi");
            assert!(
                !split.contains("Set-Cookie") && !split.contains("X-B"),
                "{split:?}"
            );
            assert_eq!(split.matches("\r\n\r\n").count(), 1, "{split:?}");
        }
    }

    #[test]
    fn a_reply_that_the_other_side_does_not_take_is_given_up_when_its_time_is_up() {
        // a reader that connects and then reads nothing: each write alone would wait its whole timeout, for ever
        let listener = std::net::TcpListener::bind("127.0.0.1:0").expect("a port");
        let _silent =
            TcpStream::connect(listener.local_addr().expect("an address")).expect("a connection");
        let (stream, _) = listener.accept().expect("a connection");
        // the stream's own timeout is long, so that a test that fails does not hang
        stream
            .set_write_timeout(Some(Duration::from_secs(5)))
            .expect("a timeout");
        let started = Instant::now();
        let mut out =
            Within::new(&stream, Duration::from_millis(300)).idle(Duration::from_secs(60));
        let big = vec![0u8; 256 << 20];
        let why = out.write_all(&big).expect_err("given up");
        assert_eq!(why.kind(), ErrorKind::TimedOut);
        assert!(
            started.elapsed() < Duration::from_secs(3),
            "took {:?}",
            started.elapsed()
        );
    }

    #[test]
    fn a_header_that_names_a_usual_one_takes_its_place_and_a_not_modified_has_no_length() {
        let kept = written(200, &[("cache-control", "no-cache".into())], b"hi");
        assert!(kept.contains("cache-control: no-cache\r\n"));
        assert!(!kept.contains("no-store"), "{kept}");
        let unchanged = written(304, &[], b"");
        assert!(unchanged.starts_with("HTTP/1.1 304 Not Modified\r\n"));
        assert!(!unchanged.contains("Content-Length"), "{unchanged}");
    }

    /// A connection made of bytes: what is read, and what was written.
    struct Duplex {
        input: Cursor<Vec<u8>>,
        output: Vec<u8>,
    }

    impl Duplex {
        fn new(bytes: &[u8]) -> Duplex {
            Duplex {
                input: Cursor::new(bytes.to_vec()),
                output: Vec::new(),
            }
        }
    }

    impl Read for Duplex {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            self.input.read(buf)
        }
    }

    impl Write for Duplex {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.output.extend(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    fn request(bytes: &[u8], limit: usize) -> Result<Request, Failure> {
        let mut d = Duplex::new(bytes);
        let head = read_head(&mut d)?;
        read_body(&mut d, head, limit)
    }

    #[test]
    fn a_get_and_a_post_are_read_with_their_headers_and_body() {
        let get = request(
            b"GET /v1/status HTTP/1.1\r\nHost: x\r\nAuthorization: Bearer abc\r\n\r\n",
            100,
        )
        .expect("a request");
        assert_eq!(
            (get.head.method.as_str(), get.head.target.as_str()),
            ("GET", "/v1/status")
        );
        assert_eq!(get.head.header("authorization"), Some("Bearer abc"));
        assert_eq!(get.head.header("AUTHORIZATION"), Some("Bearer abc"));
        assert!(get.body.is_empty());
        let post = request(
            b"POST /v1/jobs HTTP/1.1\r\nContent-Length: 5\r\n\r\nhelloEXTRA",
            100,
        )
        .expect("a request");
        assert_eq!(post.body, b"hello");
        // a body that arrives in pieces
        let big = vec![b'x'; 200_000];
        let mut bytes = format!(
            "POST /v1/jobs HTTP/1.1\r\nContent-Length: {}\r\n\r\n",
            big.len()
        )
        .into_bytes();
        bytes.extend(&big);
        assert_eq!(request(&bytes, 300_000).expect("a request").body, big);
    }

    #[test]
    fn requests_that_break_a_rule_are_refused_with_the_right_status() {
        let status = |bytes: &[u8], limit| request(bytes, limit).err().map(|f| f.status());
        assert_eq!(
            status(b"get / HTTP/1.1\r\n\r\n", 10),
            Some(400),
            "lower case method"
        );
        assert_eq!(status(b"GET / HTTP/2\r\n\r\n", 10), Some(400));
        assert_eq!(status(b"GET relative HTTP/1.1\r\n\r\n", 10), Some(400));
        assert_eq!(status(b"GET /a b HTTP/1.1\r\n\r\n", 10), Some(400));
        assert_eq!(status(b"GET\r\n\r\n", 10), Some(400));
        assert_eq!(
            status(b"GET / HTTP/1.1\r\nno colon here\r\n\r\n", 10),
            Some(400)
        );
        assert_eq!(
            status(b"GET / HTTP/1.1\r\n folded: header\r\n\r\n", 10),
            Some(400)
        );
        assert_eq!(
            status(b"GET / HTTP/1.1\r\nBad Name: x\r\n\r\n", 10),
            Some(400)
        );
        assert_eq!(
            status(b"POST / HTTP/1.1\r\nContent-Length: abc\r\n\r\n", 10),
            Some(400)
        );
        assert_eq!(
            status(b"POST / HTTP/1.1\r\nContent-Length: -5\r\n\r\n", 10),
            Some(400)
        );
        assert_eq!(
            status(
                b"POST / HTTP/1.1\r\nContent-Length: 1\r\nContent-Length: 2\r\n\r\nxx",
                10
            ),
            Some(400)
        );
        assert_eq!(
            status(b"POST / HTTP/1.1\r\nTransfer-Encoding: chunked\r\n\r\n", 10),
            Some(501)
        );
        assert_eq!(
            status(b"POST / HTTP/1.1\r\nContent-Length: 11\r\n\r\n", 10),
            Some(413)
        );
        // ends inside the header, or inside the body
        assert!(matches!(
            request(b"GET / HTTP/1.1\r\nHost:", 10),
            Err(Failure::Closed(_))
        ));
        assert!(matches!(
            request(b"POST / HTTP/1.1\r\nContent-Length: 9\r\n\r\nabc", 10),
            Err(Failure::Closed(_))
        ));
        // a header that never ends, and one with too many lines
        assert_eq!(status(&vec![b'a'; HEADER_LIMIT * 2], 10), Some(400));
        let many: String = (0..=MAX_HEADERS).map(|i| format!("H{i}: v\r\n")).collect();
        assert_eq!(
            status(format!("GET / HTTP/1.1\r\n{many}\r\n").as_bytes(), 10),
            Some(400)
        );
    }

    #[test]
    fn a_body_too_large_is_refused_before_it_is_read() {
        let mut d = Duplex::new(b"POST / HTTP/1.1\r\nContent-Length: 999999999999\r\n\r\nabc");
        let head = read_head(&mut d).expect("a header");
        assert_eq!(read_body(&mut d, head, 1000).err(), Some(Failure::TooLarge));
        // nothing was asked of the sender
        assert!(d.output.is_empty());
    }

    #[test]
    fn expect_continue_is_answered_before_the_body_is_read() {
        let mut d =
            Duplex::new(b"POST / HTTP/1.1\r\nExpect: 100-continue\r\nContent-Length: 3\r\n\r\nabc");
        let head = read_head(&mut d).expect("a header");
        // the body arrived with the header here, so no answer is needed; with an empty rest it is
        let r = read_body(&mut d, head, 10).expect("a request");
        assert_eq!(r.body, b"abc");
        let mut d =
            Duplex::new(b"POST / HTTP/1.1\r\nExpect: 100-continue\r\nContent-Length: 3\r\n\r\n");
        let head = read_head(&mut d).expect("a header");
        d.input = Cursor::new(b"abc".to_vec());
        let r = read_body(&mut d, head, 10).expect("a request");
        assert_eq!(r.body, b"abc");
        assert!(String::from_utf8_lossy(&d.output).starts_with("HTTP/1.1 100 Continue"));
    }

    #[test]
    fn a_reply_is_written_with_its_length_and_read_back() {
        let mut out = Vec::new();
        respond(&mut out, 404, "text/plain", b"no such job").expect("written");
        let text = String::from_utf8_lossy(&out).to_string();
        assert!(text.starts_with("HTTP/1.1 404 Not Found\r\n"));
        assert!(text.contains("Content-Length: 11\r\n") && text.contains("Connection: close"));
        let reply = parse_reply(&out).expect("a reply");
        assert_eq!(
            (reply.status, reply.body.as_slice()),
            (404, &b"no such job"[..])
        );
        // a reply cut short, or not a reply
        assert!(parse_reply(&out[..out.len() - 3]).is_err());
        assert!(parse_reply(b"HTTP/1.1 200 OK\r\nContent-Length: 4").is_err());
        assert!(parse_reply(b"garbage\r\n\r\n").is_err());
        assert!(parse_reply(b"").is_err());
        // a body longer than its length is trimmed to it
        let long =
            parse_reply(b"HTTP/1.1 200 OK\r\nContent-Length: 2\r\n\r\nabcd").expect("a reply");
        assert_eq!(long.body, b"ab");
    }

    /// A sender on this computer that gives `byte` every `pause`, `times` times, and a connection from it.
    fn dribbler(
        byte: u8,
        pause: Duration,
        times: usize,
    ) -> (TcpStream, std::thread::JoinHandle<()>) {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").expect("a port");
        let address = listener.local_addr().expect("an address");
        let sender = std::thread::spawn(move || {
            let mut out = TcpStream::connect(address).expect("a connection");
            for _ in 0..times {
                if out.write_all(&[byte]).is_err() {
                    break;
                }
                std::thread::sleep(pause);
            }
        });
        let (stream, _) = listener.accept().expect("an accepted connection");
        (stream, sender)
    }

    #[test]
    fn a_header_sent_a_byte_at_a_time_is_dropped_when_its_time_is_up() {
        // a byte every 40 ms would pass any timeout set on each read; the time for the whole header is what counts
        let (stream, sender) = dribbler(b'a', Duration::from_millis(40), 200);
        let started = Instant::now();
        let outcome = read_head(&mut Within::new(&stream, Duration::from_millis(300)));
        assert!(matches!(outcome, Err(Failure::Closed(_))), "{outcome:?}");
        assert!(
            started.elapsed() < Duration::from_secs(3),
            "{:?}",
            started.elapsed()
        );
        drop(stream);
        sender.join().expect("the sender ends");
    }

    #[test]
    fn a_header_that_comes_in_time_is_read_whole() {
        let listener = std::net::TcpListener::bind("127.0.0.1:0").expect("a port");
        let address = listener.local_addr().expect("an address");
        let sender = std::thread::spawn(move || {
            let mut out = TcpStream::connect(address).expect("a connection");
            out.write_all(
                b"GET /v1/status HTTP/1.1\r\nHost: x\r\nAuthorization: Bearer abc\r\n\r\n",
            )
            .expect("sent");
        });
        let (stream, _) = listener.accept().expect("an accepted connection");
        let head = read_head(&mut Within::new(&stream, Duration::from_secs(5))).expect("a header");
        assert_eq!(
            (head.method.as_str(), head.target.as_str()),
            ("GET", "/v1/status")
        );
        assert_eq!(head.header("authorization"), Some("Bearer abc"));
        sender.join().expect("the sender ends");
    }

    #[test]
    fn what_is_thrown_away_after_a_reply_is_thrown_away_for_a_short_time_only() {
        // a sender that goes on for ten seconds is let go after DRAIN_TIME
        let (mut stream, sender) = dribbler(b'z', Duration::from_millis(100), 100);
        let started = Instant::now();
        finish(&mut stream);
        assert!(
            started.elapsed() < DRAIN_TIME + Duration::from_secs(2),
            "{:?}",
            started.elapsed()
        );
        drop(stream);
        sender.join().expect("the sender ends");
    }

    #[test]
    fn an_unreachable_server_is_reported_in_words() {
        // nothing listens on port 1 of this computer
        let why = call(
            "127.0.0.1:1",
            "GET",
            "/v1/status",
            "token",
            b"",
            SMALL_REPLY,
        )
        .err()
        .expect("an error");
        assert!(why.contains("127.0.0.1:1"), "{why}");
        assert!(call(
            "not a host name at all:99999",
            "GET",
            "/",
            "t",
            b"",
            SMALL_REPLY
        )
        .is_err());
    }
}
