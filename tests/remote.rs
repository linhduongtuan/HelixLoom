//! The workstation server as a person runs it: the real program, started with `serve` and spoken to over plain HTTP with
//! none of the program's own code in between. Each job runs in a process of its own, as it does on a workstation, so this
//! checks the wire format, the token, the limits, the job processes and the results kept across a restart from outside.

use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpStream;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::mpsc;
use std::time::{Duration, Instant};

const TOKEN: &str = "an-end-to-end-token-of-some-length";
const WAIT: Duration = Duration::from_secs(60);
const FASTA: &str = ">a\nACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGATGC\n\
                     >b\nACGTACGTAGCTAGCTTGCTAGCATCGATCGATCGATGC\n\
                     >c\nACGTACGTAGCTAGCTAGCTAGCATCGATCGTTCGATGC\n";

/// A folder for one test, removed when the test is over.
struct Folder(PathBuf);

impl Folder {
    fn new(label: &str) -> Folder {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let path = std::env::temp_dir().join(format!(
            "helixloom-e2e-{label}-{}-{nanos}",
            std::process::id()
        ));
        std::fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }

    fn join(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    /// A file holding `text`, readable by this user only.
    fn file(&self, name: &str, text: &str) -> PathBuf {
        let path = self.join(name);
        std::fs::write(&path, text).expect("a file");
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o600))
                .expect("permissions");
        }
        path
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A server process that is stopped when this goes out of scope.
struct Running {
    child: Child,
    address: String,
}

impl Drop for Running {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

/// Start `helixloom serve` on a free port of this computer, with no GPU, and wait until it says where it listens.
fn start(folder: &Folder, token_file: &Path, extra: &[&str]) -> Running {
    let mut command = serve_command(folder, extra);
    command.arg("--token-file").arg(token_file);
    launch(command)
}

/// The command of a server on a free port of this computer, with no GPU unless `extra` says otherwise.
fn serve_command(folder: &Folder, extra: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_helixloom"));
    command
        .args([
            "serve",
            "--listen",
            "127.0.0.1:0",
            "--gpus",
            "none",
            "--name",
            "end-to-end",
        ])
        .arg("--data-dir")
        .arg(folder.join("data"))
        .args(extra);
    command
}

/// Run a server command and wait until it says where it listens.
fn launch(mut command: Command) -> Running {
    let child = command
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .expect("the server starts");
    let mut running = Running {
        child,
        address: String::new(),
    };
    let stdout = running.child.stdout.take().expect("the server's output");
    let (send, receive) = mpsc::channel();
    std::thread::spawn(move || {
        for line in BufReader::new(stdout).lines().map_while(Result::ok) {
            if let Some(at) = line.strip_prefix("HelixLoom server listening on ") {
                let _ = send.send(at.trim().to_string());
            }
        }
    });
    running.address = receive
        .recv_timeout(WAIT)
        .expect("the server says where it listens");
    running
}

/// One request, one reply: the status and the body.
fn call(
    address: &str,
    method: &str,
    path: &str,
    token: Option<&str>,
    body: &[u8],
) -> (u16, Vec<u8>) {
    let mut stream = TcpStream::connect(address).expect("a connection");
    stream.set_read_timeout(Some(WAIT)).expect("a timeout");
    let auth = token.map_or(String::new(), |t| format!("Authorization: Bearer {t}\r\n"));
    let head = format!(
        "{method} {path} HTTP/1.1\r\nHost: {address}\r\n{auth}Content-Length: {}\r\nConnection: close\r\n\r\n",
        body.len()
    );
    stream.write_all(head.as_bytes()).expect("sent");
    stream.write_all(body).expect("sent");
    let mut reply = Vec::new();
    stream.read_to_end(&mut reply).expect("a reply");
    split_reply(&reply)
}

fn split_reply(reply: &[u8]) -> (u16, Vec<u8>) {
    let end = reply
        .windows(4)
        .position(|w| w == b"\r\n\r\n")
        .expect("a header in the reply");
    let head = String::from_utf8_lossy(&reply[..end]).into_owned();
    let status = head
        .split(' ')
        .nth(1)
        .and_then(|s| s.parse().ok())
        .expect("a status line");
    (status, reply[end + 4..].to_vec())
}

/// Send only a header that announces `length` bytes of body, and read the answer: a server that refuses must do so
/// from the header, without waiting for bytes that are never coming.
fn call_header_only(address: &str, path: &str, token: Option<&str>, length: usize) -> u16 {
    let mut stream = TcpStream::connect(address).expect("a connection");
    stream.set_read_timeout(Some(WAIT)).expect("a timeout");
    let auth = token.map_or(String::new(), |t| format!("Authorization: Bearer {t}\r\n"));
    let head = format!(
        "POST {path} HTTP/1.1\r\n{auth}Content-Length: {length}\r\nConnection: close\r\n\r\n"
    );
    stream.write_all(head.as_bytes()).expect("sent");
    let mut reply = Vec::new();
    let _ = stream.read_to_end(&mut reply);
    split_reply(&reply).0
}

/// A container: for each section `@name length`, the bytes, and a line break.
fn container(sections: &[(&str, &[u8])]) -> Vec<u8> {
    let mut out = Vec::new();
    for (name, data) in sections {
        out.extend(format!("@{name} {}\n", data.len()).as_bytes());
        out.extend(*data);
        out.push(b'\n');
    }
    out
}

fn parse(mut rest: &[u8]) -> Vec<(String, Vec<u8>)> {
    let mut out = Vec::new();
    while !rest.is_empty() {
        let end = rest
            .iter()
            .position(|&b| b == b'\n')
            .expect("a header line");
        let head = std::str::from_utf8(&rest[..end]).expect("a header in text");
        let (name, length) = head
            .strip_prefix('@')
            .and_then(|h| h.split_once(' '))
            .expect("\"@name length\"");
        let length: usize = length.parse().expect("a length");
        out.push((name.to_string(), rest[end + 1..end + 1 + length].to_vec()));
        assert_eq!(rest[end + 1 + length], b'\n', "a line break after {name}");
        rest = &rest[end + 2 + length..];
    }
    out
}

fn section(parsed: &[(String, Vec<u8>)], name: &str) -> String {
    let found = parsed.iter().find(|(n, _)| n == name);
    String::from_utf8_lossy(&found.unwrap_or_else(|| panic!("no \"{name}\" section")).1)
        .into_owned()
}

fn align_job(label: &str, fasta: &str) -> Vec<u8> {
    let job = format!("kind=align\nlabel={label}\n");
    container(&[
        ("meta", b"engine=builtin\n"),
        ("sequences", fasta.as_bytes()),
        ("job", job.as_bytes()),
    ])
}

/// Send a job; its id.
fn submit(address: &str, body: &[u8]) -> String {
    let (status, reply) = call(address, "POST", "/v1/jobs", Some(TOKEN), body);
    assert_eq!(status, 201, "{}", String::from_utf8_lossy(&reply));
    section(&parse(&reply), "id")
}

/// The fields of a job's line in the list: id, kind, status, created, started, finished, gpu, needs_gpu, label, error.
fn job_line(address: &str, id: &str) -> Vec<String> {
    let (status, reply) = call(address, "GET", &format!("/v1/jobs/{id}"), Some(TOKEN), b"");
    assert_eq!(status, 200, "{}", String::from_utf8_lossy(&reply));
    let line = section(&parse(&reply), "job");
    let fields: Vec<String> = line
        .trim_end_matches(['\r', '\n'])
        .split('\t')
        .map(str::to_string)
        .collect();
    assert_eq!(fields.len(), 10, "{line:?}");
    fields
}

/// Wait until the job is over (done, failed or cancelled) and give its line.
fn wait_over(address: &str, id: &str) -> Vec<String> {
    let end = Instant::now() + WAIT;
    loop {
        let fields = job_line(address, id);
        if matches!(fields[2].as_str(), "done" | "failed" | "cancelled") {
            return fields;
        }
        assert!(Instant::now() < end, "the job is still {}", fields[2]);
        std::thread::sleep(Duration::from_millis(100));
    }
}

fn result(address: &str, id: &str) -> Vec<u8> {
    let (status, reply) = call(
        address,
        "GET",
        &format!("/v1/jobs/{id}/result"),
        Some(TOKEN),
        b"",
    );
    assert_eq!(status, 200, "{}", String::from_utf8_lossy(&reply));
    reply
}

#[test]
fn a_job_runs_in_its_own_process_and_its_result_waits_for_the_person() {
    let folder = Folder::new("jobs");
    let token = folder.file("token.txt", TOKEN);
    let server = start(&folder, &token, &[]);
    let at = server.address.as_str();

    // nothing is told to a visitor without the token, or with another one, or with part of it
    for header in [
        None,
        Some("wrong-token-of-the-same-length-x"),
        Some(&TOKEN[..TOKEN.len() - 1]),
    ] {
        assert_eq!(
            call(at, "GET", "/v1/status", header, b"").0,
            401,
            "{header:?}"
        );
        assert_eq!(
            call(at, "GET", "/v1/jobs", header, b"").0,
            401,
            "{header:?}"
        );
    }
    let (code, reply) = call(at, "GET", "/v1/status", Some(TOKEN), b"");
    assert_eq!(code, 200);
    let status = parse(&reply);
    let meta = section(&status, "meta");
    assert!(
        meta.contains("name=end-to-end\n")
            && meta.contains("protocol=1\n")
            && meta.contains("running=0\n"),
        "{meta}"
    );
    assert_eq!(
        section(&status, "gpus"),
        "",
        "the server was told to use no GPU"
    );

    // a job goes in, runs, and its answer is the alignment, kept until it is fetched
    let id = submit(at, &align_job("three isolates", FASTA));
    let line = wait_over(at, &id);
    assert_eq!(
        (line[1].as_str(), line[2].as_str(), line[8].as_str()),
        ("align", "done", "three isolates"),
        "{line:?}"
    );
    // what was sent is not kept once the job is over: the result is
    let kept = folder.join("data").join("jobs").join(&id);
    assert!(kept.join("result.bin").is_file() && !kept.join("request.bin").exists());
    let answer = parse(&result(at, &id));
    assert!(section(&answer, "meta").contains("engine="), "{answer:?}");
    let rows = section(&answer, "rows");
    let aligned: Vec<&str> = rows.lines().filter(|l| !l.starts_with('>')).collect();
    assert_eq!(aligned.len(), 3, "{rows}");
    assert!(
        aligned.iter().all(|r| r.len() == aligned[0].len()),
        "an alignment's rows are as long as each other"
    );
    assert!(rows.contains(">a\n") && rows.contains(">c\n"));
    let (_, listed) = call(at, "GET", "/v1/jobs", Some(TOKEN), b"");
    assert!(section(&parse(&listed), "jobs").contains(&id));

    // a job that cannot be done says why instead of vanishing
    let one = ">only\nACGT\n";
    let bad = submit(at, &align_job("a lonely sequence", one));
    let line = wait_over(at, &bad);
    assert_eq!(line[2], "failed", "{line:?}");
    assert!(line[9].contains("at least two sequences"), "{line:?}");
    let (code, _) = call(
        at,
        "GET",
        &format!("/v1/jobs/{bad}/result"),
        Some(TOKEN),
        b"",
    );
    assert_eq!(code, 409, "a failed job has no result");

    // requests that make no sense are refused in words, and the names of jobs never reach the disk
    let nonsense = container(&[("job", b"kind=mine-bitcoin\nlabel=x\n")]);
    assert_eq!(call(at, "POST", "/v1/jobs", Some(TOKEN), &nonsense).0, 400);
    assert_eq!(
        call(at, "POST", "/v1/jobs", Some(TOKEN), b"not a container").0,
        400
    );
    for path in [
        "/v1/jobs/j1-1",
        "/v1/jobs/../../etc/passwd",
        "/v1/jobs/..%2f..%2fetc",
        "/v1/jobs/j1-1/result",
        "/v1/nothing",
    ] {
        assert_eq!(call(at, "GET", path, Some(TOKEN), b"").0, 404, "{path}");
    }
    assert_eq!(call(at, "PUT", "/v1/jobs", Some(TOKEN), b"").0, 405);
    assert_eq!(
        call(at, "GET", "/v1/status", Some(TOKEN), b"").0,
        200,
        "the server is still well"
    );

    // removing a job removes its result
    assert_eq!(
        call(at, "DELETE", &format!("/v1/jobs/{id}"), Some(TOKEN), b"").0,
        204
    );
    assert_eq!(
        call(at, "GET", &format!("/v1/jobs/{id}"), Some(TOKEN), b"").0,
        404
    );
    assert_eq!(
        call(at, "DELETE", &format!("/v1/jobs/{id}"), Some(TOKEN), b"").0,
        404
    );
}

/// Six aligned sequences of 60 columns, each a little different from the first.
fn alignment() -> String {
    let base = b"ACGTACGTAGCTAGCTAGCTAGCATCGATCGATCGATGCTAGCTAGCATCGATCGATCGA";
    (0..6)
        .map(|r| {
            let mut row = base.to_vec();
            for k in 0..r * 3 {
                row[(k * 7 + r) % base.len()] = b"ACGT"[(k + r) % 4];
            }
            format!(">r{r}\n{}\n", String::from_utf8_lossy(&row))
        })
        .collect()
}

/// A job that takes far longer than a test does: the longest Bayesian run the server accepts.
fn long_tree_job() -> Vec<u8> {
    let alignment = alignment();
    container(&[
        (
            "meta",
            b"analysis=bayes\nmodel=p\nbootstrap=0\ngenerations=10000000\n",
        ),
        ("alignment", alignment.as_bytes()),
        ("job", b"kind=tree\nlabel=a long tree\n"),
    ])
}

/// The job processes of this folder that are alive, from the list of processes.
#[cfg(unix)]
fn job_processes(folder: &Folder) -> usize {
    let listing = Command::new("ps")
        .args(["-axww", "-o", "command="])
        .output()
        .expect("ps runs");
    let here = folder.0.to_string_lossy().into_owned();
    String::from_utf8_lossy(&listing.stdout)
        .lines()
        .filter(|l| l.contains("run-job") && l.contains(&here))
        .count()
}

#[cfg(unix)]
#[test]
fn a_running_job_is_stopped_when_asked_and_its_process_ends() {
    let folder = Folder::new("stop");
    let token = folder.file("token.txt", TOKEN);
    let server = start(&folder, &token, &[]);
    let at = server.address.as_str();
    let id = submit(at, &long_tree_job());
    let end = Instant::now() + WAIT;
    while job_line(at, &id)[2] != "running" {
        assert!(Instant::now() < end, "the job never started");
        std::thread::sleep(Duration::from_millis(50));
    }
    assert_eq!(job_processes(&folder), 1, "one process runs the job");

    assert_eq!(
        call(at, "DELETE", &format!("/v1/jobs/{id}"), Some(TOKEN), b"").0,
        204
    );
    let line = wait_over(at, &id);
    assert_eq!(line[2], "cancelled", "{line:?}");
    while job_processes(&folder) > 0 {
        assert!(Instant::now() < end, "the job's process is still there");
        std::thread::sleep(Duration::from_millis(50));
    }
    assert_eq!(
        call(
            at,
            "GET",
            &format!("/v1/jobs/{id}/result"),
            Some(TOKEN),
            b""
        )
        .0,
        409
    );
    // the stopped job stays in the list until it is removed again
    assert_eq!(
        call(at, "DELETE", &format!("/v1/jobs/{id}"), Some(TOKEN), b"").0,
        204
    );
    assert_eq!(
        call(at, "GET", &format!("/v1/jobs/{id}"), Some(TOKEN), b"").0,
        404
    );
    // and the server goes on working
    let next = submit(at, &align_job("after the stop", FASTA));
    assert_eq!(wait_over(at, &next)[2], "done");
}

#[cfg(unix)]
#[test]
fn a_job_ends_by_itself_when_its_server_is_gone() {
    let folder = Folder::new("orphan");
    let token = folder.file("token.txt", TOKEN);
    let mut server = start(&folder, &token, &[]);
    let at = server.address.clone();
    let id = submit(&at, &long_tree_job());
    let end = Instant::now() + WAIT;
    while job_line(&at, &id)[2] != "running" {
        assert!(Instant::now() < end, "the job never started");
        std::thread::sleep(Duration::from_millis(50));
    }
    assert_eq!(job_processes(&folder), 1);
    // the server stops without a word, as when it crashes or is stopped with Ctrl-C
    server.child.kill().expect("the server is stopped");
    server.child.wait().expect("the server is gone");
    while job_processes(&folder) > 0 {
        assert!(Instant::now() < end, "the job went on without its server");
        std::thread::sleep(Duration::from_millis(100));
    }
}

fn fold_job(label: &str) -> Vec<u8> {
    fold_job_saying(&format!("kind=fold\nlabel={label}\n"))
}

/// A folding job whose sender asks for `gpus` GPUs (0 asks for every GPU there is).
#[cfg(unix)]
fn fold_job_asking(label: &str, gpus: usize) -> Vec<u8> {
    fold_job_saying(&format!("kind=fold\nlabel={label}\ngpus={gpus}\n"))
}

fn fold_job_saying(job: &str) -> Vec<u8> {
    container(&[
        ("meta", b"local=1\n"),
        ("sequence", b"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ"),
        ("job", job.as_bytes()),
    ])
}

/// What GPUs are given to a job's process: "0", "0,1", ..., or nothing for a job that uses none.
fn gpu_of(line: &[String]) -> &str {
    &line[6]
}

/// The command of a stand-in for a folding program: it takes a moment, then writes the GPUs it was shown.
#[cfg(unix)]
fn fold_command(folder: &Folder) -> String {
    let script = folder.file(
        "fold.sh",
        "#!/bin/sh\nsleep 1.5\nprintf 'REMARK GPU=%s\\nATOM      1  N   MET A   1       0.000   0.000   0.000  1.00 90.00           N\\nEND\\n' \"$CUDA_VISIBLE_DEVICES\" > \"$2\"\n",
    );
    format!("sh {}", script.display())
}

/// The GPUs that the stand-in program of a finished job says it was shown.
#[cfg(unix)]
fn gpus_seen(address: &str, id: &str) -> String {
    let pdb = section(&parse(&result(address, id)), "pdb");
    pdb.lines()
        .find_map(|l| l.strip_prefix("REMARK GPU="))
        .expect("the program wrote the GPUs it saw")
        .to_string()
}

/// Wait until a job is running; its line.
#[cfg(unix)]
fn wait_running(address: &str, id: &str) -> Vec<String> {
    let end = Instant::now() + WAIT;
    loop {
        let fields = job_line(address, id);
        if fields[2] == "running" {
            return fields;
        }
        assert!(Instant::now() < end, "the job is still {}", fields[2]);
        std::thread::sleep(Duration::from_millis(20));
    }
}

/// Who holds each GPU of the server, by the GPU's place in the list: a job's id, or nothing for a free one.
#[cfg(unix)]
fn gpu_holders(address: &str) -> Vec<String> {
    let (_, reply) = call(address, "GET", "/v1/status", Some(TOKEN), b"");
    let rows = section(&parse(&reply), "gpus");
    rows.lines()
        .map(|row| row.split('\t').nth(2).unwrap_or("").to_string())
        .collect()
}

#[cfg(unix)]
#[test]
fn each_gpu_job_gets_a_gpu_of_its_own_sees_only_that_one_and_the_next_waits_for_a_free_one() {
    let folder = Folder::new("gpus");
    let token = folder.file("token.txt", TOKEN);
    let command = fold_command(&folder);
    let server = start(
        &folder,
        &token,
        &["--gpus", "0,1", "--fold-command", &command],
    );
    let at = server.address.as_str();

    let (_, reply) = call(at, "GET", "/v1/status", Some(TOKEN), b"");
    let status = parse(&reply);
    assert!(section(&status, "meta").contains("fold=1\n"));
    assert_eq!(section(&status, "gpus").lines().count(), 2, "{status:?}");

    let ids: Vec<String> = (1..=3)
        .map(|n| submit(at, &fold_job(&format!("fold {n}"))))
        .collect();
    // two jobs run at once, one on each GPU; the third waits for a GPU to be free
    let end = Instant::now() + WAIT;
    let running = loop {
        let lines: Vec<Vec<String>> = ids.iter().map(|id| job_line(at, id)).collect();
        if lines.iter().filter(|l| l[2] == "running").count() == 2 {
            break lines;
        }
        assert!(
            Instant::now() < end,
            "two jobs never ran together: {lines:?}"
        );
        std::thread::sleep(Duration::from_millis(30));
    };
    let mut on: Vec<&str> = running
        .iter()
        .filter(|l| l[2] == "running")
        .map(|l| gpu_of(l))
        .collect();
    on.sort_unstable();
    assert_eq!(on, ["0", "1"], "{running:?}");
    let waiting = running
        .iter()
        .find(|l| l[2] == "queued")
        .expect("a job that waits");
    assert_eq!(waiting[7], "1", "it is waiting for a GPU");
    let (_, reply) = call(at, "GET", "/v1/status", Some(TOKEN), b"");
    let busy = section(&parse(&reply), "gpus");
    assert!(
        ids.iter().filter(|id| busy.contains(id.as_str())).count() == 2,
        "both GPUs are busy: {busy}"
    );

    // each job's program saw exactly one GPU, and the two that ran together saw different ones
    let mut seen = Vec::new();
    for id in &ids {
        let line = wait_over(at, id);
        assert_eq!(line[2], "done", "{line:?}");
        let gpu = gpus_seen(at, id);
        assert!(
            gpu == gpu_of(&line),
            "the job's line says GPU {} but the program saw {gpu}",
            gpu_of(&line)
        );
        seen.push(gpu);
    }
    let first_two: std::collections::BTreeSet<&String> = seen.iter().take(2).collect();
    assert_eq!(first_two.len(), 2, "{seen:?}");
    assert!(seen[2] == "0" || seen[2] == "1", "{seen:?}");
}

#[cfg(unix)]
#[test]
fn a_job_is_given_several_gpus_together_and_the_jobs_behind_it_keep_their_places() {
    let folder = Folder::new("several-gpus");
    let token = folder.file("token.txt", TOKEN);
    let command = fold_command(&folder);
    // three GPUs, and the folding program may use two of them at once
    let server = start(
        &folder,
        &token,
        &[
            "--gpus",
            "0,1,2",
            "--fold-gpus",
            "2",
            "--fold-command",
            &command,
        ],
    );
    let at = server.address.as_str();

    // two jobs that ask for two GPUs each, and behind them one that asks for a single GPU
    let first = submit(at, &fold_job_asking("two", 2));
    let second = submit(at, &fold_job_asking("two more", 2));
    let third = submit(at, &fold_job_asking("one", 1));
    let line = wait_running(at, &first);
    assert_eq!((gpu_of(&line), line[7].as_str()), ("0,1", "2"), "{line:?}");
    // GPU 2 is free, but the second job cannot use it alone, and the third, which could, stays behind the second
    for (id, wanted) in [(&second, "2"), (&third, "1")] {
        let waiting = job_line(at, id);
        let (status, gpus) = (waiting[2].as_str(), gpu_of(&waiting));
        assert_eq!((status, gpus, waiting[7].as_str()), ("queued", "", wanted));
    }
    // both of the first job's GPUs show it
    assert_eq!(gpu_holders(at), [first.as_str(), first.as_str(), ""]);

    // when the first job is over the other two run together, and each program was shown exactly its job's GPUs
    let mut seen = Vec::new();
    for id in [&first, &second, &third] {
        let line = wait_over(at, id);
        assert_eq!(line[2], "done", "{line:?}");
        let gpus = gpus_seen(at, id);
        assert_eq!(gpus, gpu_of(&line), "the job's line says what it was shown");
        seen.push(gpus);
    }
    assert_eq!(seen, ["0,1", "0,1", "2"]);
}

#[cfg(unix)]
#[test]
fn a_job_is_given_one_gpu_when_the_server_has_not_said_its_program_uses_more() {
    let folder = Folder::new("one-gpu-each");
    let token = folder.file("token.txt", TOKEN);
    let command = fold_command(&folder);
    // three GPUs; no `--fold-gpus`, so the folding program is given one GPU at a time
    let server = start(
        &folder,
        &token,
        &["--gpus", "0,1,2", "--fold-command", &command],
    );
    let at = server.address.as_str();

    // two GPUs asked for, every GPU asked for, and more GPUs than there are
    let ids: Vec<String> = [2, 0, 9]
        .iter()
        .map(|&asked| submit(at, &fold_job_asking(&format!("asks {asked}"), asked)))
        .collect();
    for id in &ids {
        assert_eq!(job_line(at, id)[7], "1", "one GPU is wanted");
    }
    let mut seen = Vec::new();
    for id in &ids {
        let line = wait_over(at, id);
        assert_eq!(line[2], "done", "{line:?}");
        let gpus = gpus_seen(at, id);
        assert_eq!(gpus, gpu_of(&line), "the job's line says what it was shown");
        assert!(!gpus.contains(','), "one GPU, not {gpus}");
        seen.push(gpus);
    }
    // the three ran together, a GPU each
    seen.sort_unstable();
    assert_eq!(seen, ["0", "1", "2"]);
}

#[test]
fn results_and_waiting_jobs_survive_a_restart_of_the_server() {
    let folder = Folder::new("restart");
    let token = folder.file("token.txt", TOKEN);
    let (id, kept) = {
        let server = start(&folder, &token, &[]);
        let id = submit(&server.address, &align_job("kept across a restart", FASTA));
        assert_eq!(wait_over(&server.address, &id)[2], "done");
        let kept = result(&server.address, &id);
        (id, kept)
    };
    // the first server is gone (stopped when it went out of scope); a second one uses the same folder
    let again = start(&folder, &token, &[]);
    let line = job_line(&again.address, &id);
    assert_eq!(
        (line[2].as_str(), line[8].as_str()),
        ("done", "kept across a restart")
    );
    assert_eq!(
        result(&again.address, &id),
        kept,
        "the result is exactly what was stored"
    );
    // and it still takes new work
    let next = submit(&again.address, &align_job("after the restart", FASTA));
    assert_ne!(next, id);
    assert_eq!(wait_over(&again.address, &next)[2], "done");
}

#[test]
fn an_upload_over_the_limit_is_refused_from_its_header_and_a_stranger_gets_nothing_read() {
    let folder = Folder::new("limits");
    let token = folder.file("token.txt", TOKEN);
    let server = start(&folder, &token, &["--max-upload-mb", "1"]);
    let at = server.address.as_str();
    let two_mebibytes = 2 << 20;
    assert_eq!(
        call_header_only(at, "/v1/jobs", Some(TOKEN), two_mebibytes),
        413
    );
    // without the right token the server does not even ask for the body
    assert_eq!(call_header_only(at, "/v1/jobs", None, two_mebibytes), 401);
    assert_eq!(
        call_header_only(
            at,
            "/v1/jobs",
            Some("not-the-token-at-all-not-at-all-no"),
            two_mebibytes
        ),
        401
    );
    // a body that is not announced by a length is not accepted
    let mut stream = TcpStream::connect(at).expect("a connection");
    stream.set_read_timeout(Some(WAIT)).expect("a timeout");
    let head = format!("POST /v1/jobs HTTP/1.1\r\nAuthorization: Bearer {TOKEN}\r\nTransfer-Encoding: chunked\r\n\r\n");
    stream.write_all(head.as_bytes()).expect("sent");
    let mut reply = Vec::new();
    let _ = stream.read_to_end(&mut reply);
    assert_eq!(split_reply(&reply).0, 501);
    assert_eq!(
        call(at, "GET", "/v1/status", Some(TOKEN), b"").0,
        200,
        "the server is still well"
    );
}

/// Run the program with `args` and give its exit code, what it printed and what it complained of.
fn run(args: &[&str]) -> (i32, String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_helixloom"))
        .args(args)
        .stdin(Stdio::null())
        .output()
        .expect("the program runs");
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).into_owned(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
    )
}

#[test]
fn the_command_line_makes_tokens_asks_a_server_and_refuses_a_weak_one() {
    let (code, out, _) = run(&["token"]);
    let token = out.trim();
    assert!(
        code == 0 && token.len() >= 32 && token.bytes().all(|b| b.is_ascii_hexdigit()),
        "{out:?}"
    );
    assert_ne!(run(&["token"]).1, out, "a token is never made twice");

    let folder = Folder::new("cli");
    let good = folder.file("token.txt", TOKEN);
    let server = start(&folder, &good, &[]);
    let (code, out, err) = run(&[
        "status",
        "--server",
        &server.address,
        "--token-file",
        good.to_str().expect("a path"),
    ]);
    assert_eq!(code, 0, "{err}");
    assert!(
        out.contains("end-to-end") && out.contains("GPUs: none"),
        "{out}"
    );
    let wrong = folder.file("wrong.txt", "another-token-of-some-length-here");
    let (code, _, err) = run(&[
        "status",
        "--server",
        &server.address,
        "--token-file",
        wrong.to_str().expect("a path"),
    ]);
    assert_eq!(code, 1);
    assert!(err.contains("refused"), "{err}");

    // a server with a weak token does not start at all, and one with no token does not either
    let weak = folder.file("weak.txt", "short");
    let (code, _, err) = run(&[
        "serve",
        "--listen",
        "127.0.0.1:0",
        "--gpus",
        "none",
        "--token-file",
        weak.to_str().expect("a path"),
        "--data-dir",
        folder.join("weak-data").to_str().expect("a path"),
    ]);
    assert_ne!(code, 0);
    assert!(err.contains("32"), "{err}");
    let out = Command::new(env!("CARGO_BIN_EXE_helixloom"))
        .args(["serve", "--listen", "127.0.0.1:0", "--gpus", "none"])
        .env_remove("HELIXLOOM_TOKEN")
        .stdin(Stdio::null())
        .output()
        .expect("the program runs");
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("no token"));
}

#[cfg(unix)]
#[test]
fn the_token_in_the_environment_is_taken_out_before_any_program_is_run() {
    use std::os::unix::fs::PermissionsExt;
    let folder = Folder::new("env");
    let bin = folder.join("bin");
    std::fs::create_dir_all(&bin).expect("a folder");
    let dump = folder.join("nvidia-smi-environment.txt");
    // a stand-in nvidia-smi that records the environment it is started in, and a folding program that says what it saw
    let smi = bin.join("nvidia-smi");
    std::fs::write(
        &smi,
        format!(
            "#!/bin/sh\nenv > '{}'\necho 'GPU 0: Fake GPU (UUID: GPU-1234)'\n",
            dump.display()
        ),
    )
    .expect("a script");
    std::fs::set_permissions(&smi, std::fs::Permissions::from_mode(0o755)).expect("permissions");
    let fold = folder.file(
        "fold.sh",
        "#!/bin/sh\nprintf 'REMARK TOKEN=%s ORDER=%s GPU=%s\\nATOM      1  N   MET A   1       0.000   0.000   0.000  1.00 90.00           N\\nEND\\n' \"${HELIXLOOM_TOKEN:-none}\" \"$CUDA_DEVICE_ORDER\" \"$CUDA_VISIBLE_DEVICES\" > \"$2\"\n",
    );
    let path = format!(
        "{}:{}",
        bin.display(),
        std::env::var("PATH").unwrap_or_default()
    );
    // the token comes from the environment, as it may on a workstation, and there is no token file at all
    let mut command = serve_command(
        &folder,
        &[
            "--gpus",
            "auto",
            "--fold-command",
            &format!("sh {}", fold.display()),
        ],
    );
    command.env("HELIXLOOM_TOKEN", TOKEN).env("PATH", path);
    let server = launch(command);
    let at = server.address.as_str();

    let (_, reply) = call(at, "GET", "/v1/status", Some(TOKEN), b"");
    let gpus = section(&parse(&reply), "gpus");
    assert!(gpus.starts_with("0\tFake GPU\t"), "{gpus:?}");
    let seen = std::fs::read_to_string(&dump).expect("the stand-in ran");
    assert!(
        seen.contains("PATH="),
        "it was given an environment: {seen}"
    );
    assert!(
        !seen.contains("HELIXLOOM_TOKEN"),
        "nvidia-smi was handed the token"
    );

    // a job's own program is not handed it either, and is shown the GPU as nvidia-smi numbers it
    let id = submit(at, &fold_job("what the program sees"));
    assert_eq!(wait_over(at, &id)[2], "done");
    let pdb = section(&parse(&result(at, &id)), "pdb");
    assert!(
        pdb.contains("REMARK TOKEN=none ORDER=PCI_BUS_ID GPU=0"),
        "{pdb}"
    );
}
