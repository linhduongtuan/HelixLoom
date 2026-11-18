//! What a browser needs from the server besides the jobs: the page itself (`--web-dir`) and permission to call the server
//! from a page that another server hosts (`--allow-origin`).
//!
//! The page is read into memory once, when the server starts, and a request is only ever looked up in that table: nothing
//! a request says is used as a path, so no request can name a file that was not in the folder at the start. The page is
//! public (it is the program, not data) and is sent without a token; every `/v1/` request still needs one.
//!
//! A file that has a packed copy beside it (`name.gz`, as `scripts/build-web.sh` makes them) is sent packed to a browser
//! that says it takes gzip: the program is 16 MB and its packed copy a third of that. The server packs nothing itself.

use super::http::Head;
use super::one_line;
use super::server::{Body, Reply};
use std::collections::HashMap;
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Duration;

/// A page may be at most this many files, of at most `MAX_FILE` bytes each and `MAX_TOTAL` in all: it is held in memory.
const MAX_FILES: usize = 1000;
const MAX_FILE: u64 = 64 << 20;
const MAX_TOTAL: u64 = 128 << 20;
/// How many folders down the page's files are looked for.
const MAX_DEPTH: usize = 4;
/// How long a browser may keep the answer to a question about what it may send, in seconds.
const PREFLIGHT_KEPT: &str = "600";
/// The page goes to anyone who asks, with no token, so it is sent to only so many at a time (the others are told to try
/// again), which leaves places for the person who has the token (a server serves 64 connections at once)...
const MAX_PAGE_SENDS: usize = 16;
/// ...and in a time that even a slow link meets: this long, and as long as the file needs at this many bytes a second (the
/// 13 MB of the page: about a minute, which a link of 2 Mbit/s still meets). The slower floor of a token holder's results
/// (`http::MIN_UPLOAD_RATE`) would let a stranger who reads a little at a time keep a place for five minutes.
const PAGE_FIRST: Duration = Duration::from_secs(15);
const PAGE_MIN_RATE: usize = 256 * 1024;

/// What the name of a packed copy ends with, and the two bytes every gzip file begins with.
const PACKED: &str = ".gz";
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

/// The checksum a gzip file ends with (CRC-32, as zlib and PNG compute it) of `bytes`.
fn crc32(bytes: &[u8]) -> u32 {
    // the table of the 256 bytes, each through the eight steps of the sum
    let table: [u32; 256] = std::array::from_fn(|byte| {
        (0..8).fold(byte as u32, |sum, _| {
            if sum & 1 == 1 {
                0xedb8_8320 ^ (sum >> 1)
            } else {
                sum >> 1
            }
        })
    });
    !bytes.iter().fold(!0u32, |sum, &byte| {
        table[((sum ^ u32::from(byte)) & 0xff) as usize] ^ (sum >> 8)
    })
}

/// Whether `packed` is the gzip of `file`, as far as can be told without unpacking it (the server has nothing to unpack
/// with): it begins as a gzip file does, and its last eight bytes, which hold the checksum and the length of what was packed,
/// are the checksum and the length of `file`. A copy left over from another build is of other bytes and fails this, though it
/// be of the same length (the page is, in every build: the names in it are as long each time).
fn is_packed_copy(packed: &[u8], file: &[u8]) -> bool {
    let Some(tail) = packed.len().checked_sub(8) else {
        return false;
    };
    let word = |at: usize| {
        u32::from_le_bytes([packed[at], packed[at + 1], packed[at + 2], packed[at + 3]])
    };
    // the length is kept in 32 bits: of a longer file, what is left of it after that many
    packed.starts_with(&GZIP_MAGIC)
        && word(tail + 4) == file.len() as u32
        && word(tail) == crc32(file)
}

/// Whether a request says that gzip may be sent: `gzip` among what `Accept-Encoding` lists, and not with a weight of nothing.
fn takes_gzip(head: &Head) -> bool {
    head.header("accept-encoding").is_some_and(|listed| {
        listed.split(',').any(|item| {
            let mut parts = item.split(';').map(str::trim);
            let name = parts.next().unwrap_or_default();
            let refused = parts.any(|p| {
                p.strip_prefix("q=")
                    .and_then(|q| q.parse::<f32>().ok())
                    .is_some_and(|q| q <= 0.0)
            });
            name.eq_ignore_ascii_case("gzip") && !refused
        })
    })
}

/// Extra reply headers, name and value.
pub type Headers = Vec<(&'static str, String)>;

/// One file of the page, as it is sent.
pub struct File {
    pub content_type: &'static str,
    /// A tag of its bytes, so that a browser that has the file is told so instead of being sent it again.
    pub etag: String,
    pub bytes: Vec<u8>,
}

/// The page and the origins that may call the server.
#[derive(Default)]
pub struct Site {
    files: HashMap<String, Arc<File>>,
    /// The packed copies, under the path of the file each is a copy of.
    packed: HashMap<String, Arc<File>>,
    origins: Vec<String>,
    /// How many files of the page are being sent at this moment.
    sending: AtomicUsize,
}

/// A place among those that the page is sent to; given back when this goes out of scope, however the sending ends.
pub struct SendPlace<'a>(&'a AtomicUsize);

impl Drop for SendPlace<'_> {
    fn drop(&mut self) {
        self.0.fetch_sub(1, Ordering::SeqCst);
    }
}

/// How long a file of `length` bytes may take to be sent to someone who asked for the page.
pub fn page_time(length: usize) -> Duration {
    PAGE_FIRST + Duration::from_secs((length / PAGE_MIN_RATE) as u64)
}

/// The type a browser is told for a file name, by its extension.
pub fn content_type(name: &str) -> &'static str {
    let extension = Path::new(name)
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_ascii_lowercase);
    match extension.as_deref() {
        Some("html") => "text/html; charset=utf-8",
        Some("js" | "mjs") => "text/javascript; charset=utf-8",
        Some("css") => "text/css; charset=utf-8",
        // a browser refuses to compile a program that is not sent with this exact type
        Some("wasm") => "application/wasm",
        Some("json") => "application/json",
        // the manifest a browser installs the page by: the type the standard gives it
        Some("webmanifest") => "application/manifest+json",
        Some("png") => "image/png",
        Some("svg") => "image/svg+xml",
        Some("ico") => "image/x-icon",
        Some("txt") => "text/plain; charset=utf-8",
        Some("woff2") => "font/woff2",
        Some("gz") => "application/gzip",
        _ => "application/octet-stream",
    }
}

/// An origin as `--allow-origin` takes it: `http://host[:port]` or `https://host[:port]`, in lower case, or why not.
pub fn parse_origin(text: &str) -> Result<String, String> {
    let bad = || {
        format!(
            "\"{}\" is not an origin: write it as https://host or http://host:port, with no path.",
            one_line(text)
        )
    };
    let lower = text.trim().to_ascii_lowercase();
    let rest = lower
        .strip_prefix("http://")
        .or_else(|| lower.strip_prefix("https://"))
        .ok_or_else(bad)?;
    let plain = rest
        .bytes()
        .all(|b| b.is_ascii_alphanumeric() || b"-._:[]".contains(&b));
    // a host, and a port after it (an IPv6 address is in brackets, and has colons of its own)
    let (host, port) = match rest.rsplit_once(':') {
        Some((host, port)) if !host.contains(':') || host.ends_with(']') => (host, Some(port)),
        _ => (rest, None),
    };
    let host_ok = !host.is_empty() && (!host.contains(':') || host.starts_with('['));
    let port_ok = port.is_none_or(|p| p.parse::<u16>().is_ok_and(|n| n > 0));
    if !(plain && host_ok && port_ok) {
        return Err(bad());
    }
    // a browser leaves the default port out of the Origin it sends: spelt with it, the origin could never match
    let default = if lower.starts_with("https://") {
        ":443"
    } else {
        ":80"
    };
    Ok(lower.strip_suffix(default).unwrap_or(&lower).to_string())
}

/// A name of a file that can be served: plain characters only, and not hidden, so that what a request can ask for is
/// exactly what it can spell.
fn servable(name: &str) -> bool {
    !name.is_empty()
        && !name.starts_with('.')
        && name
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b))
}

/// A tag for some bytes (FNV-1a): to tell whether a browser's copy is the file, not to guard against anyone.
fn etag(bytes: &[u8]) -> String {
    let hash = bytes.iter().fold(0xcbf2_9ce4_8422_2325_u64, |h, &b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    });
    format!("\"{hash:016x}\"")
}

/// Read the files of `dir` into `files`, under the path a request would use (`prefix` and the name). Links and anything that
/// is not a plain file or folder are left out.
fn collect(
    dir: &Path,
    prefix: &str,
    depth: usize,
    files: &mut HashMap<String, Arc<File>>,
    total: &mut u64,
) -> Result<(), String> {
    let read = |e: std::io::Error| format!("Cannot read {}: {e}", dir.display());
    for entry in std::fs::read_dir(dir).map_err(read)? {
        let entry = entry.map_err(read)?;
        let file_name = entry.file_name();
        let Some(name) = file_name.to_str().filter(|n| servable(n)) else {
            continue;
        };
        let path = format!("{prefix}/{name}");
        // the type of the entry itself: a link is a link, whatever it points at
        let kind = entry.file_type().map_err(read)?;
        if kind.is_dir() {
            if depth < MAX_DEPTH {
                collect(&entry.path(), &path, depth + 1, files, total)?;
            }
        } else if kind.is_file() {
            if files.len() >= MAX_FILES {
                return Err(format!("The web folder has more than {MAX_FILES} files."));
            }
            let length = entry.metadata().map_err(read)?.len();
            *total += length;
            if length > MAX_FILE || *total > MAX_TOTAL {
                return Err(format!(
                    "The web folder is too large to serve from memory ({path} makes it more than {} MB).",
                    MAX_TOTAL >> 20
                ));
            }
            let bytes = std::fs::read(entry.path()).map_err(read)?;
            files.insert(
                path,
                Arc::new(File {
                    content_type: content_type(name),
                    etag: etag(&bytes),
                    bytes,
                }),
            );
        }
    }
    Ok(())
}

impl Site {
    /// The files of `dir` (and its folders, a few levels down), served as the page; `origins` may call the server.
    pub fn load(dir: &Path, origins: &[String]) -> Result<Site, String> {
        let meta = std::fs::metadata(dir)
            .map_err(|e| format!("Cannot read the web folder {}: {e}", dir.display()))?;
        if !meta.is_dir() {
            return Err(format!("{} is not a folder.", dir.display()));
        }
        let mut files = HashMap::new();
        collect(dir, "", 0, &mut files, &mut 0)?;
        if files.keys().any(|path| path.starts_with("/v1/")) {
            return Err(format!(
                "{} has a folder called v1: the addresses under /v1/ are the server's, and need a token.",
                dir.display()
            ));
        }
        let index = files.get("/index.html").cloned().ok_or_else(|| {
            format!(
                "{} has no index.html: build the web version first (trunk build --release, see docs/WEB.md).",
                dir.display()
            )
        })?;
        files.insert("/".to_string(), index);
        // a file with its packed copy beside it: the copy is sent in its place to a browser that takes it, and is no longer
        // a file of the page by its own name (it is held once). A `.gz` that is not the copy of the file beside it stays
        // what it is: a file of the folder.
        let mut packed: HashMap<String, Arc<File>> = HashMap::new();
        let paired: Vec<String> = files
            .iter()
            .filter(|(path, file)| {
                files
                    .get(&format!("{path}{PACKED}"))
                    .is_some_and(|copy| is_packed_copy(&copy.bytes, &file.bytes))
            })
            .map(|(path, _)| path.clone())
            .collect();
        for path in paired {
            let Some(content_type) = files.get(&path).map(|file| file.content_type) else {
                continue;
            };
            let Some(copy) = files.remove(&format!("{path}{PACKED}")) else {
                continue;
            };
            let copy = Arc::try_unwrap(copy).unwrap_or_else(|shared| File {
                content_type: shared.content_type,
                etag: shared.etag.clone(),
                bytes: shared.bytes.clone(),
            });
            packed.insert(
                path,
                Arc::new(File {
                    content_type,
                    ..copy
                }),
            );
        }
        // the root is the page
        if let Some(copy) = packed.get("/index.html").cloned() {
            packed.insert("/".to_string(), copy);
        }
        Ok(Site {
            files,
            packed,
            origins: origins.to_vec(),
            ..Site::default()
        })
    }

    /// A server that serves no page and lets `origins` call it.
    pub fn with_origins(origins: &[String]) -> Site {
        Site {
            origins: origins.to_vec(),
            ..Site::default()
        }
    }

    /// A place to send a file of the page in, or `None` when it is being sent to as many as it may be.
    pub fn take_send_place(&self) -> Option<SendPlace<'_>> {
        let mut sending = self.sending.load(Ordering::SeqCst);
        loop {
            if sending >= MAX_PAGE_SENDS {
                return None;
            }
            match self.sending.compare_exchange(
                sending,
                sending + 1,
                Ordering::SeqCst,
                Ordering::SeqCst,
            ) {
                Ok(_) => return Some(SendPlace(&self.sending)),
                Err(now) => sending = now,
            }
        }
    }

    /// The path of a file of the page that holds `needle`, if one does: the page is public, so nothing secret may be in it.
    pub fn file_containing(&self, needle: &[u8]) -> Option<String> {
        if needle.is_empty() {
            return None;
        }
        let mut paths: Vec<&String> = self.files.keys().filter(|p| p.as_str() != "/").collect();
        paths.sort();
        paths
            .into_iter()
            .find(|p| {
                self.files[*p]
                    .bytes
                    .windows(needle.len())
                    .any(|w| w == needle)
            })
            .cloned()
    }

    /// The headers an answer to a request from `origin` carries about who may read it: the permission when the origin was
    /// named, and, whenever any origin was, that the answer depends on who asked (so that a cache in front, as a proxy is,
    /// does not hand one origin the copy made for another).
    pub fn origin_headers(&self, origin: Option<&str>) -> Headers {
        match self.allowed_origin(origin) {
            Some(allowed) => cors_headers(allowed),
            None if !self.origins.is_empty() => vec![("Vary", "Origin".to_string())],
            None => Vec::new(),
        }
    }

    /// Whether there is a page to serve.
    pub fn has_page(&self) -> bool {
        !self.files.is_empty()
    }

    /// The origins that may call the server.
    pub fn origins(&self) -> &[String] {
        &self.origins
    }

    /// The paths that are served, sorted.
    #[cfg(test)]
    pub fn paths(&self) -> Vec<String> {
        let mut paths: Vec<String> = self.files.keys().cloned().collect();
        paths.sort();
        paths
    }

    /// The origin, in the spelling it was allowed with, when a request from it may be answered across origins.
    pub fn allowed_origin(&self, origin: Option<&str>) -> Option<&str> {
        let origin = origin?;
        self.origins
            .iter()
            .find(|allowed| allowed.as_str() == origin)
            .map(String::as_str)
    }

    /// An answer that needs no token: a file of the page, a "not modified", or the permission a browser asks for before it
    /// calls from another origin. `None` for everything else.
    pub fn answer(&self, head: &Head) -> Option<(Reply, Headers)> {
        match head.method.as_str() {
            "GET" => self.page_file(head),
            "OPTIONS" => self.preflight(head),
            _ => None,
        }
    }

    fn page_file(&self, head: &Head) -> Option<(Reply, Headers)> {
        let path = head.target.split('?').next().unwrap_or_default();
        let plain = self.files.get(path)?;
        let copy = self.packed.get(path);
        let packed = copy.filter(|_| takes_gzip(head));
        let file = packed.unwrap_or(plain);
        let mut headers = vec![
            ("ETag", file.etag.clone()),
            // the browser may keep the file but must ask whether it is still the one: a new build replaces it
            ("Cache-Control", "no-cache".to_string()),
            // the page takes a token as typed: no other page may frame it, and a link out of it says where it came from to no one
            ("X-Frame-Options", "DENY".to_string()),
            ("Referrer-Policy", "no-referrer".to_string()),
        ];
        if packed.is_some() {
            headers.push(("Content-Encoding", "gzip".to_string()));
        }
        if copy.is_some() {
            // which of the two is sent depends on what the browser takes: a cache in front must not hand one the other's
            headers.push(("Vary", "Accept-Encoding".to_string()));
        }
        let has_it = head.header("if-none-match").is_some_and(|tags| {
            tags.split(',')
                .map(|t| t.trim().trim_start_matches("W/"))
                .any(|t| t == "*" || t == file.etag)
        });
        let (status, body) = if has_it {
            (304, Vec::new().into())
        } else {
            (200, Body::Page(Arc::clone(file)))
        };
        let reply = Reply {
            status,
            content_type: file.content_type,
            body,
        };
        Some((reply, headers))
    }

    /// The answer to a browser that asks what a page of an allowed origin may send.
    fn preflight(&self, head: &Head) -> Option<(Reply, Headers)> {
        self.allowed_origin(head.header("origin"))?;
        head.header("access-control-request-method")?;
        let mut headers = vec![
            (
                "Access-Control-Allow-Methods",
                "GET, POST, DELETE".to_string(),
            ),
            (
                "Access-Control-Allow-Headers",
                "authorization, content-type".to_string(),
            ),
            ("Access-Control-Max-Age", PREFLIGHT_KEPT.to_string()),
        ];
        // a page on the internet that reaches for a computer on a private network is asked about as well
        if head
            .header("access-control-request-private-network")
            .is_some_and(|v| v.eq_ignore_ascii_case("true"))
        {
            headers.push(("Access-Control-Allow-Private-Network", "true".to_string()));
        }
        let reply = Reply {
            status: 204,
            content_type: "text/plain",
            body: Vec::new().into(),
        };
        Some((reply, headers))
    }
}

/// The headers that tell a browser it may read the answer, for a request from `origin` (one that was allowed).
pub fn cors_headers(origin: &str) -> Headers {
    vec![
        ("Access-Control-Allow-Origin", origin.to_string()),
        // the answer depends on who asked, so a cache must not hand it to somebody else
        ("Vary", "Origin".to_string()),
    ]
}

#[cfg(test)]
mod tests;
