use super::*;
use crate::pipeline::tools::work_folder;
use crate::remote::client::{Remote, ServerStatus};
use crate::remote::gpus::Gpu;
use crate::remote::tasks::{Request, Response};
use std::io::{Read, Write};

const TOKEN: &str = "a-token-of-at-least-sixteen-chars";

fn config(label: &str) -> Config {
    let dir = work_folder(label).expect("a folder");
    let mut c = Config::new(TOKEN.into(), dir);
    c.gpus = vec![Gpu {
        id: 0,
        name: "GPU A".into(),
    }];
    c.max_upload = 4096;
    c.tools = vec![("MAFFT".into(), "7.5".into())];
    c
}

/// A server that queues jobs and never runs them, on a free port of this computer.
fn serving(label: &str) -> Server {
    start(
        config(label),
        "127.0.0.1:0".parse().expect("an address"),
        false,
    )
    .expect("a server")
}

fn align_request() -> Request {
    Request::Align {
        seqs: vec![
            ("a".into(), b"ACGTACGT".to_vec()),
            ("b".into(), b"ACGTTCGT".to_vec()),
        ],
        mafft: false,
    }
}

fn upload(request: &Request) -> Vec<u8> {
    request
        .encode()
        .with(
            "job",
            format!("kind={}\nlabel=a test\n", request.kind().name()),
        )
        .encode()
}

fn text(reply: &Reply) -> String {
    String::from_utf8_lossy(&reply.body).to_string()
}

#[test]
fn tokens_are_compared_whole_and_taken_from_the_header_only_as_a_bearer() {
    assert!(same_token("abc", "abc"));
    assert!(!same_token("abc", "abd"));
    assert!(!same_token("abc", "abcd"));
    assert!(!same_token("", "abc"));
    assert!(!same_token("abcd", "abc"));
    assert!(same_token("", ""));
    assert_eq!(bearer(Some("Bearer secret")), "secret");
    assert_eq!(bearer(Some("bearer  secret ")), "secret");
    for none in [None, Some("Basic secret"), Some("Bearer"), Some("secret")] {
        assert_eq!(bearer(none), "", "{none:?}");
    }
}

#[test]
fn the_routes_answer_status_jobs_results_and_removal() {
    let shared = Shared::open(config("server-routes"), false).expect("a store");
    let get = |target: &str| route(&shared, "GET", target, b"");

    let status = get("/v1/status");
    assert_eq!(status.status, 200);
    let s = ServerStatus::parse(&Container::decode(&status.body).expect("a container"))
        .expect("a status");
    assert_eq!(
        (s.protocol, s.queued, s.running),
        (crate::remote::PROTOCOL, 0, 0)
    );
    assert_eq!(s.gpus.len(), 1);
    assert_eq!(s.tools, [("MAFFT".to_string(), "7.5".to_string())]);

    // a job goes in, shows in the list, has no result yet, and can be removed
    let sent = route(&shared, "POST", "/v1/jobs", &upload(&align_request()));
    assert_eq!(sent.status, 201, "{}", text(&sent));
    let id = Container::decode(&sent.body)
        .expect("a container")
        .require("id")
        .expect("an id")
        .to_string();
    let list = Container::decode(&get("/v1/jobs").body).expect("a container");
    assert_eq!(list.require("jobs").expect("jobs").lines().count(), 1);
    let one = Container::decode(&get(&format!("/v1/jobs/{id}")).body).expect("a container");
    assert!(one.require("job").expect("a job").contains("a test"));
    assert_eq!(get(&format!("/v1/jobs/{id}/result")).status, 409);
    assert_eq!(
        route(&shared, "DELETE", &format!("/v1/jobs/{id}"), b"").status,
        204
    );
    assert_eq!(get(&format!("/v1/jobs/{id}")).status, 404);
    assert_eq!(
        route(&shared, "DELETE", &format!("/v1/jobs/{id}"), b"").status,
        404
    );
    // a query string changes nothing
    assert_eq!(get("/v1/status?x=1").status, 200);
}

#[test]
fn a_finished_jobs_result_is_handed_over_as_stored() {
    let shared = Shared::open(config("server-result"), false).expect("a store");
    let id = shared
        .submit(Kind::Align, "x", b"request", 0)
        .expect("queued");
    let mut job = shared.get(&id).expect("a job");
    job.status = Status::Done;
    std::fs::write(shared.job_dir(&id).join("result.bin"), b"@meta 0\n\n").expect("a result");
    shared.replace_for_tests(job);
    let r = route(&shared, "GET", &format!("/v1/jobs/{id}/result"), b"");
    assert_eq!((r.status, r.body.as_slice()), (200, &b"@meta 0\n\n"[..]));
    assert_eq!(r.content_type, "application/x-helixloom");
}

#[test]
fn requests_that_are_wrong_are_refused_with_the_right_status() {
    let shared = Shared::open(config("server-wrong"), false).expect("a store");
    let status =
        |method: &str, target: &str, body: &[u8]| route(&shared, method, target, body).status;
    assert_eq!(status("GET", "/", b""), 404);
    assert_eq!(status("GET", "/v2/status", b""), 404);
    assert_eq!(status("GET", "/v1/nothing", b""), 404);
    assert_eq!(status("POST", "/v1/status", b""), 405);
    assert_eq!(status("PUT", "/v1/jobs", b""), 405);
    assert_eq!(status("GET", "/v1/jobs/j1-1", b""), 404, "no such job");
    assert_eq!(status("PUT", "/v1/jobs/j1-1", b""), 405);
    // ids that are not ids never reach the store
    for id in [
        "..",
        "../..",
        "j1-1/../..",
        "%2e%2e",
        "j",
        "J1-1",
        "<script>",
    ] {
        assert_eq!(status("GET", &format!("/v1/jobs/{id}"), b""), 404, "{id}");
        assert_eq!(
            status("DELETE", &format!("/v1/jobs/{id}"), b""),
            404,
            "{id}"
        );
        assert_eq!(
            status("GET", &format!("/v1/jobs/{id}/result"), b""),
            404,
            "{id}"
        );
    }
    // uploads that are not jobs
    assert_eq!(status("POST", "/v1/jobs", b"not a container"), 400);
    assert_eq!(status("POST", "/v1/jobs", b""), 400);
    let no_job = align_request().encode().encode();
    assert_eq!(status("POST", "/v1/jobs", &no_job), 400, "no job section");
    let unknown = align_request()
        .encode()
        .with("job", "kind=bitcoin\n")
        .encode();
    assert_eq!(status("POST", "/v1/jobs", &unknown), 400);
    let broken = Container::new()
        .with("job", "kind=align\n")
        .with("meta", "engine=builtin\n")
        .with("sequences", "")
        .encode();
    assert_eq!(status("POST", "/v1/jobs", &broken), 400);
    assert!(shared.list().is_empty(), "nothing wrong was queued");
}

#[test]
fn a_job_that_needs_what_the_workstation_lacks_is_refused_at_once() {
    let mut c = config("server-lacks");
    c.gpus.clear();
    let shared = Shared::open(c, false).expect("a store");
    let cnn = |newick: &str| Request::PhyloCnn {
        newick: newick.into(),
        sampling: 0.5,
        task: crate::seq::phylocnn::Task::Bd,
        priors: crate::seq::phylocnn::Priors::Simulated,
        gpu: true,
    };
    // many trees on a GPU need one (in a build that can use it)
    let many = route(
        &shared,
        "POST",
        "/v1/jobs",
        &upload(&cnn("((A:1,B:1):1,C:2);\n((A:1,B:2):1,C:2);")),
    );
    if cfg!(feature = "gpu") {
        assert_eq!(many.status, 400);
        assert!(text(&many).contains("no GPU"), "{}", text(&many));
    } else {
        assert_eq!(many.status, 201, "{}", text(&many));
    }
    // one tree runs on the processor whatever was asked, so it needs no GPU
    let one = route(
        &shared,
        "POST",
        "/v1/jobs",
        &upload(&cnn("((A:1,B:1):1,C:2);")),
    );
    assert_eq!(one.status, 201, "{}", text(&one));
    // a local fold needs a folding program, a fold through the web service does not
    let fold = |local| Request::Fold {
        sequence: b"MQIFVKTLTGKTITLEVEP".to_vec(),
        local,
    };
    let r = route(&shared, "POST", "/v1/jobs", &upload(&fold(true)));
    assert_eq!(r.status, 400);
    let web = route(&shared, "POST", "/v1/jobs", &upload(&fold(false)));
    assert_eq!(web.status, 201, "{}", text(&web));
}

/// A fold job sent with the given `gpus=` line (none: as a program from before jobs could have several GPUs sends it).
fn fold_upload(gpus: Option<&str>) -> Vec<u8> {
    let line = gpus.map_or(String::new(), |g| format!("gpus={g}\n"));
    let request = Request::Fold {
        sequence: b"MQIFVKTLTGKTITLEVEP".to_vec(),
        local: true,
    };
    request
        .encode()
        .with("job", format!("kind=fold\nlabel=a fold\n{line}"))
        .encode()
}

#[test]
fn a_job_is_queued_for_the_gpus_it_asked_for_within_what_the_workstation_has_and_the_program_uses()
{
    let mut c = config("server-gpus");
    c.gpus = (0..4)
        .map(|id| Gpu {
            id,
            name: format!("GPU {id}"),
        })
        .collect();
    c.fold_command = Some(vec!["fold".into()]);
    // the folding program uses two GPUs for a structure
    c.fold_gpus = 2;
    let shared = Shared::open(c, false).expect("a store");
    let wanted = |gpus: Option<&str>| {
        let sent = route(&shared, "POST", "/v1/jobs", &fold_upload(gpus));
        assert_eq!(sent.status, 201, "{}", text(&sent));
        let id = Container::decode(&sent.body)
            .expect("a container")
            .require("id")
            .expect("an id")
            .to_string();
        shared.get(&id).expect("a job").gpus_wanted
    };
    assert_eq!(wanted(None), 1, "a program that says nothing asks for one");
    assert_eq!(wanted(Some("1")), 1);
    assert_eq!(wanted(Some("2")), 2);
    assert_eq!(wanted(Some("4")), 2, "the program uses two, so two");
    assert_eq!(
        wanted(Some("0")),
        2,
        "every GPU, within what the program uses"
    );
    assert_eq!(wanted(Some(" 2 ")), 2);
    assert_eq!(
        wanted(Some("many")),
        1,
        "words nobody could mean are one GPU"
    );
    assert_eq!(wanted(Some("-1")), 1);
    // a job that needs no GPU is queued for none, whatever it asked for
    let align = request_with_gpus(&align_request(), "3");
    let sent = route(&shared, "POST", "/v1/jobs", &align);
    assert_eq!(sent.status, 201, "{}", text(&sent));
    let id = Container::decode(&sent.body)
        .expect("a container")
        .require("id")
        .expect("an id")
        .to_string();
    assert_eq!(shared.get(&id).expect("a job").gpus_wanted, 0);
}

/// `request` as an upload that asks for `gpus` GPUs.
fn request_with_gpus(request: &Request, gpus: &str) -> Vec<u8> {
    request
        .encode()
        .with(
            "job",
            format!(
                "kind={}\nlabel=a test\ngpus={gpus}\n",
                request.kind().name()
            ),
        )
        .encode()
}

#[test]
fn a_phylocnn_job_is_given_no_more_gpus_than_it_has_trees() {
    if !cfg!(feature = "gpu") {
        return;
    }
    let mut c = config("server-trees");
    c.gpus = (0..4)
        .map(|id| Gpu {
            id,
            name: format!("GPU {id}"),
        })
        .collect();
    let shared = Shared::open(c, false).expect("a store");
    let cnn = Request::PhyloCnn {
        newick: "((A:1,B:1):1,C:2);\n((A:1,B:2):1,C:2);\n((A:2,B:1):1,C:2);".into(),
        sampling: 0.5,
        task: crate::seq::phylocnn::Task::Bd,
        priors: crate::seq::phylocnn::Priors::Simulated,
        gpu: true,
    };
    // three trees: three GPUs of the four at most, whatever is asked for
    for (asked, given) in [("1", 1), ("2", 2), ("0", 3), ("4", 3)] {
        let sent = route(&shared, "POST", "/v1/jobs", &request_with_gpus(&cnn, asked));
        assert_eq!(sent.status, 201, "{}", text(&sent));
        let id = Container::decode(&sent.body)
            .expect("a container")
            .require("id")
            .expect("an id")
            .to_string();
        assert_eq!(
            shared.get(&id).expect("a job").gpus_wanted,
            given,
            "asked for {asked}"
        );
    }
}

#[test]
fn over_the_network_a_wrong_token_is_refused_and_the_right_one_is_served() {
    let server = serving("server-net");
    let address = server.addr.to_string();
    let good = Remote::new(&address, TOKEN).expect("a connection");
    let status = good.status().expect("a status");
    assert_eq!(status.gpus.len(), 1);
    let bad = Remote::new(&address, "another-token-of-sixteen-chars").expect("a connection");
    assert!(bad.status().unwrap_err().contains("refused the token"));
    assert!(bad.jobs().is_err() && bad.submit(&align_request(), "x").is_err());
    assert!(server.shared.list().is_empty());
    // no token at all, and a token that is a prefix of the right one
    for header in [
        "",
        "Authorization: Bearer \r\n",
        &format!("Authorization: Bearer {}\r\n", &TOKEN[..TOKEN.len() - 1]),
    ] {
        let reply = raw(&server, &format!("GET /v1/status HTTP/1.1\r\n{header}\r\n"));
        assert!(reply.starts_with("HTTP/1.1 401"), "{header:?} gave {reply}");
    }
}

/// Send raw text to the server and read what it answers.
fn raw(server: &Server, request: &str) -> String {
    let mut stream = TcpStream::connect(server.addr).expect("a connection");
    stream
        .set_read_timeout(Some(Duration::from_secs(10)))
        .expect("a timeout");
    stream.write_all(request.as_bytes()).expect("sent");
    // nothing more is coming: a server that waits for more sees the end
    let _ = stream.shutdown(std::net::Shutdown::Write);
    let mut got = String::new();
    let _ = stream.read_to_string(&mut got);
    got
}

#[test]
fn a_job_goes_over_the_network_into_the_queue_and_comes_out_of_it() {
    let server = serving("server-queue");
    let remote = Remote::new(&server.addr.to_string(), TOKEN).expect("a connection");
    let id = remote
        .submit(&align_request(), "two short sequences\nand a second line")
        .expect("sent");
    let jobs = remote.jobs().expect("jobs");
    assert_eq!(jobs.len(), 1);
    assert_eq!(
        (jobs[0].id.as_str(), jobs[0].status),
        (id.as_str(), Status::Queued)
    );
    assert_eq!(jobs[0].label, "two short sequences and a second line");
    assert_eq!(remote.job(&id).expect("a job").kind, Kind::Align);
    assert!(remote
        .result(Kind::Align, &id)
        .err()
        .expect("not ready")
        .contains("not ready"));
    remote.delete(&id).expect("removed");
    assert!(remote.jobs().expect("jobs").is_empty());
    assert!(remote.job(&id).is_err());
    // the result of a job comes back through the same path as the server stored it
    let id = remote.submit(&align_request(), "again").expect("sent");
    let response = Response::Align(crate::pipeline::steps::Aligned {
        rows: vec![
            ("a".into(), b"AC-GT".to_vec()),
            ("b".into(), b"ACCGT".to_vec()),
        ],
        engine: "stand-in".into(),
        seconds: 0.5,
    });
    let mut job = server.shared.get(&id).expect("a job");
    job.status = Status::Done;
    std::fs::write(
        server.shared.job_dir(&id).join("result.bin"),
        response.encode().encode(),
    )
    .expect("a result");
    server.shared.replace_for_tests(job);
    let Response::Align(got) = remote.result(Kind::Align, &id).expect("a result") else {
        panic!("an alignment");
    };
    assert_eq!((got.rows.len(), got.engine.as_str()), (2, "stand-in"));
}

#[test]
fn nothing_unauthorised_is_read_and_nothing_too_large_is_accepted() {
    let server = serving("server-limits");
    // a huge upload with the wrong token is refused from its header: the body is never asked for
    let started = std::time::Instant::now();
    let reply = raw(&server, "POST /v1/jobs HTTP/1.1\r\nAuthorization: Bearer wrong\r\nContent-Length: 999999999999\r\n\r\n");
    assert!(reply.starts_with("HTTP/1.1 401"), "{reply}");
    assert!(started.elapsed() < Duration::from_secs(5));
    // with the right token, a body over the limit (4096 in these tests) is refused as too large
    let head = format!(
        "POST /v1/jobs HTTP/1.1\r\nAuthorization: Bearer {TOKEN}\r\nContent-Length: 8192\r\n\r\n"
    );
    let reply = raw(&server, &format!("{head}{}", "x".repeat(8192)));
    assert!(
        reply.starts_with("HTTP/1.1 413"),
        "{}",
        &reply[..reply.len().min(60)]
    );
    // a body that is not chunked-and-sized, and one with an impossible length
    let chunked = format!("POST /v1/jobs HTTP/1.1\r\nAuthorization: Bearer {TOKEN}\r\nTransfer-Encoding: chunked\r\n\r\n");
    assert!(raw(&server, &chunked).starts_with("HTTP/1.1 501"));
    let nonsense = format!(
        "POST /v1/jobs HTTP/1.1\r\nAuthorization: Bearer {TOKEN}\r\nContent-Length: nope\r\n\r\n"
    );
    assert!(raw(&server, &nonsense).starts_with("HTTP/1.1 400"));
    // garbage, and a request that stops short
    assert!(raw(&server, "GARBAGE\r\n\r\n").starts_with("HTTP/1.1 400"));
    assert!(raw(&server, "").is_empty());
    // small requests that carry a body they should not are only allowed a small one
    let sneaky = format!("GET /v1/status HTTP/1.1\r\nAuthorization: Bearer {TOKEN}\r\nContent-Length: 100000\r\n\r\n");
    assert!(raw(&server, &sneaky).starts_with("HTTP/1.1 413"));
    // and after all that the server still answers
    assert!(Remote::new(&server.addr.to_string(), TOKEN)
        .expect("a connection")
        .status()
        .is_ok());
}

#[test]
fn a_token_that_is_too_short_stops_the_server_from_starting() {
    let mut c = config("server-short");
    c.token = "short".into();
    let why = start(c, "127.0.0.1:0".parse().expect("an address"), false)
        .err()
        .expect("an error");
    assert!(why.contains("at least"), "{why}");
}

#[test]
fn a_connection_that_dribbles_its_header_is_let_go_when_the_time_is_up() {
    let shared = Shared::open(config("server-dribble"), false).expect("a store");
    let listener = std::net::TcpListener::bind("127.0.0.1:0").expect("a port");
    let address = listener.local_addr().expect("an address");
    // a byte every 40 ms: no timeout set on each read would ever fire, and the header never ends
    let sender = std::thread::spawn(move || {
        let mut out = TcpStream::connect(address).expect("a connection");
        for byte in b"GET /v1/status HTTP/1.1\r\nX-Slow: "
            .iter()
            .cycle()
            .take(300)
        {
            if out.write_all(&[*byte]).is_err() {
                break;
            }
            std::thread::sleep(Duration::from_millis(40));
        }
    });
    let (mut stream, _) = listener.accept().expect("an accepted connection");
    let started = std::time::Instant::now();
    serve_connection(&mut stream, &shared, Duration::from_millis(300));
    assert!(
        started.elapsed() < Duration::from_secs(3),
        "{:?}",
        started.elapsed()
    );
    drop(stream);
    sender.join().expect("the sender ends");
}

#[test]
fn only_so_many_bytes_of_uploads_are_taken_in_at_once() {
    let mut c = config("server-room");
    // room for two uploads of this size in memory at once, and no more
    c.max_upload = 600_000;
    let server = start(c, "127.0.0.1:0".parse().expect("an address"), false).expect("a server");
    let announce = |length: usize| {
        format!("POST /v1/jobs HTTP/1.1\r\nAuthorization: Bearer {TOKEN}\r\nContent-Length: {length}\r\n\r\n")
    };
    // two uploads that have said how big they are and are waiting for the rest: they hold their room
    let hold = || {
        let mut s = TcpStream::connect(server.addr).expect("a connection");
        s.write_all(announce(600_000).as_bytes()).expect("sent");
        s
    };
    let held = (hold(), hold());
    // the server takes a moment to see them: wait until both hold their room (each counted three times over, as it is in
    // memory), so that a third finds it taken
    let end = std::time::Instant::now() + Duration::from_secs(20);
    while server.shared.uploading_bytes() < 2 * 3 * 600_000 {
        assert!(
            std::time::Instant::now() < end,
            "the two uploads never arrived"
        );
        std::thread::sleep(Duration::from_millis(10));
    }
    // a third is told to try again
    let end = std::time::Instant::now() + Duration::from_secs(20);
    loop {
        let reply = raw(&server, &announce(600_000));
        if reply.starts_with("HTTP/1.1 503") {
            assert!(reply.contains("try again"), "{reply}");
            break;
        }
        assert!(std::time::Instant::now() < end, "never refused: {reply:?}");
        std::thread::sleep(Duration::from_millis(50));
    }
    // a request with no body needs no room, so the server goes on answering
    assert!(Remote::new(&server.addr.to_string(), TOKEN)
        .expect("a connection")
        .status()
        .is_ok());
    // when the first two are gone their room comes back
    drop(held);
    let end = std::time::Instant::now() + Duration::from_secs(20);
    while raw(&server, &announce(600_000)).starts_with("HTTP/1.1 503") {
        assert!(std::time::Instant::now() < end, "the room never came back");
        std::thread::sleep(Duration::from_millis(50));
    }
}

#[test]
fn the_places_for_connections_are_shared_out_and_given_back() {
    let places = Arc::new(Places::default());
    let address = |text: &str| -> IpAddr { text.parse().expect("an address") };
    // one address can have MAX_PER_PEER places and no more, while another is not held up by it
    let mine: Vec<Place> = (0..MAX_PER_PEER)
        .map(|_| places.take(address("10.0.0.1")).expect("a place"))
        .collect();
    assert!(places.take(address("10.0.0.1")).is_none());
    assert!(places.take(address("10.0.0.2")).is_some());
    // a place given back is free again
    drop(mine);
    assert!(places.take(address("10.0.0.1")).is_some());
    // and in all, whoever asks, there are MAX_CONNECTIONS places
    let mut held = Vec::new();
    for k in 0..MAX_CONNECTIONS / MAX_PER_PEER {
        for _ in 0..MAX_PER_PEER {
            held.push(
                places
                    .take(address(&format!("10.1.0.{k}")))
                    .expect("a place"),
            );
        }
    }
    assert!(
        places.take(address("10.2.0.1")).is_none(),
        "every place is taken"
    );
    drop(held);
    assert!(places.take(address("10.2.0.1")).is_some());
}

#[test]
fn one_computer_cannot_keep_the_server_from_everybody_else() {
    let server = serving("server-peer");
    // as many connections as one address may have, saying nothing and staying there
    let silent: Vec<TcpStream> = (0..MAX_PER_PEER)
        .map(|_| TcpStream::connect(server.addr).expect("a connection"))
        .collect();
    // the server sees them within a moment; then the next connection from this address is turned away, with the reply
    // sent before a word is said (a connection that is let in says nothing for a moment, and is not turned away)
    let probe = || {
        let mut s = TcpStream::connect(server.addr).expect("a connection");
        s.set_read_timeout(Some(Duration::from_millis(500)))
            .expect("a timeout");
        let mut got = String::new();
        let _ = s.read_to_string(&mut got);
        got
    };
    let end = std::time::Instant::now() + Duration::from_secs(20);
    loop {
        let reply = probe();
        if reply.starts_with("HTTP/1.1 503") {
            assert!(reply.contains("Too many connections"), "{reply}");
            break;
        }
        assert!(std::time::Instant::now() < end, "never refused: {reply:?}");
        std::thread::sleep(Duration::from_millis(30));
    }
    // when they go, the server answers again
    drop(silent);
    let remote = Remote::new(&server.addr.to_string(), TOKEN).expect("a connection");
    let end = std::time::Instant::now() + Duration::from_secs(20);
    while remote.status().is_err() {
        assert!(
            std::time::Instant::now() < end,
            "the places never came back"
        );
        std::thread::sleep(Duration::from_millis(30));
    }
}

#[test]
fn an_upload_may_take_as_long_as_its_size_needs_at_the_slowest_rate_that_is_waited_for() {
    assert_eq!(upload_time(0), Duration::from_secs(60));
    // 50 KiB a second: a megabyte takes 20 seconds more
    assert_eq!(upload_time(1 << 20), Duration::from_secs(60 + 20));
    // the largest default upload is given hours, a slow link is not turned away, but not for ever
    let largest = upload_time(512 << 20);
    assert!(
        largest > Duration::from_secs(2 * 3600) && largest < Duration::from_secs(4 * 3600),
        "{largest:?}"
    );
}

#[test]
fn a_job_that_is_removed_while_it_is_waited_for_ends_the_waiting_at_once() {
    let server = serving("server-gone");
    let remote = Remote::new(&server.addr.to_string(), TOKEN).expect("a connection");
    let waiting = remote.clone();
    let waiter = std::thread::spawn(move || waiting.run(&align_request(), "to be removed").err());
    // the job is queued (the workers are off); removing it is what the Workstation page's button does
    let end = std::time::Instant::now() + Duration::from_secs(20);
    let id = loop {
        if let Some(job) = server.shared.list().first() {
            break job.id.clone();
        }
        assert!(std::time::Instant::now() < end, "the job never arrived");
        std::thread::sleep(Duration::from_millis(20));
    };
    remote.delete(&id).expect("removed");
    // not thirty failed polls later: the answer is "no such job", which is final
    let started = std::time::Instant::now();
    let why = waiter.join().expect("the waiter ends").expect("an error");
    assert!(why.contains("no longer on the workstation"), "{why}");
    assert!(
        started.elapsed() < Duration::from_secs(10),
        "{:?}",
        started.elapsed()
    );
}
