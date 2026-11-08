use super::*;
use crate::pipeline::tools::work_folder;
use crate::remote::tasks::Response;

fn args(list: &[&str]) -> Vec<String> {
    list.iter().map(|s| (*s).to_string()).collect()
}

fn no_env(_: &str) -> Option<String> {
    None
}

fn token_file(label: &str, text: &str) -> PathBuf {
    let path = work_folder(label).expect("a folder").join("token.txt");
    std::fs::write(&path, text).expect("a token file");
    path
}

#[test]
fn the_server_starts_with_safe_defaults_given_only_a_token() {
    let file = token_file(
        "cli-defaults",
        "  a-token-of-sixteen-chars-or-more \nsecond line ignored\n",
    );
    let o = serve_options(
        &args(&[
            "--token-file",
            &file.display().to_string(),
            "--gpus",
            "none",
        ]),
        no_env,
    )
    .expect("options");
    // listens on this computer only until told otherwise
    assert_eq!(o.listen, "127.0.0.1:7878".parse().expect("an address"));
    assert!(o.listen.ip().is_loopback());
    let c = o.config;
    assert_eq!(c.token, "a-token-of-sixteen-chars-or-more");
    assert_eq!(
        (c.cpu_jobs, c.max_upload, c.max_queued),
        (2, 512 << 20, 1000)
    );
    assert_eq!(
        (c.keep.as_secs(), c.job_limit.as_secs()),
        (7 * 86400, 24 * 3600)
    );
    assert!(c.gpus.is_empty() && c.fold_command.is_none());
    // a structure is folded on one GPU unless the person who runs the server says the program uses more
    assert_eq!(c.fold_gpus, 1);
    // the folder is named as an absolute path, so that what the jobs are told about it holds wherever they run
    assert!(c.data_dir.is_absolute() && c.data_dir.ends_with("helixloom-data"));
}

#[test]
fn every_option_is_read() {
    let file = token_file("cli-all", "a-token-of-sixteen-chars-or-more");
    let o = serve_options(
        &args(&[
            "--listen",
            "0.0.0.0:9000",
            "--token-file",
            &file.display().to_string(),
            "--data-dir",
            "/var/lib/helixloom",
            "--gpus",
            "0,2-3",
            "--cpu-jobs",
            "4",
            "--max-upload-mb",
            "100",
            "--keep-days",
            "3",
            "--job-hours",
            "12",
            "--fold-command",
            "/opt/esm/run.sh --fast",
            "--fold-gpus",
            "2",
            "--name",
            "lab\nGPU",
        ]),
        no_env,
    )
    .expect("options");
    assert_eq!(o.listen.port(), 9000);
    assert_eq!(o.config.fold_gpus, 2);
    assert!(!o.listen.ip().is_loopback());
    let c = o.config;
    assert_eq!(
        c.gpus.iter().map(|g| g.id).collect::<Vec<_>>(),
        vec![0, 2, 3]
    );
    assert_eq!(
        (
            c.cpu_jobs,
            c.max_upload,
            c.keep.as_secs(),
            c.job_limit.as_secs()
        ),
        (4, 100 << 20, 3 * 86400, 12 * 3600)
    );
    assert_eq!(
        c.fold_command,
        Some(vec!["/opt/esm/run.sh".to_string(), "--fast".to_string()])
    );
    assert_eq!(c.name, "lab GPU");
    assert_eq!(c.data_dir, PathBuf::from("/var/lib/helixloom"));
}

#[test]
fn the_token_may_come_from_the_environment_but_never_from_the_command_line() {
    let env =
        |k: &str| (k == "HELIXLOOM_TOKEN").then(|| "  a-token-from-the-environment ".to_string());
    let o = serve_options(&args(&["--gpus", "none"]), env).expect("options");
    assert_eq!(o.config.token, "a-token-from-the-environment");
    // there is no --token option: a token typed on the command line would show in `ps`
    assert!(serve_options(&args(&["--token", "secret"]), env).is_err());
    let why = serve_options(&args(&["--gpus", "none"]), no_env)
        .err()
        .expect("an error");
    assert!(why.contains("no token"), "{why}");
}

#[test]
fn wrong_options_are_refused_in_words() {
    let file = token_file("cli-wrong", "a-token-of-sixteen-chars-or-more");
    let t = file.display().to_string();
    for (what, list) in [
        ("unknown", vec!["--token-file", &t, "--frobnicate", "1"]),
        ("no value", vec!["--token-file", &t, "--cpu-jobs"]),
        (
            "not a number",
            vec!["--token-file", &t, "--cpu-jobs", "many"],
        ),
        ("zero", vec!["--token-file", &t, "--cpu-jobs", "0"]),
        (
            "zero upload",
            vec!["--token-file", &t, "--max-upload-mb", "0"],
        ),
        (
            "absurd upload",
            vec!["--token-file", &t, "--max-upload-mb", "99999999999"],
        ),
        (
            "absurd days",
            vec!["--token-file", &t, "--keep-days", "99999999999999999"],
        ),
        (
            "absurd jobs",
            vec!["--token-file", &t, "--cpu-jobs", "100000"],
        ),
        ("bad gpu", vec!["--token-file", &t, "--gpus", "x"]),
        (
            "zero fold gpus",
            vec!["--token-file", &t, "--fold-gpus", "0"],
        ),
        (
            "absurd fold gpus",
            vec!["--token-file", &t, "--fold-gpus", "1000"],
        ),
        (
            "fold gpus in words",
            vec!["--token-file", &t, "--fold-gpus", "two"],
        ),
        (
            "bad address",
            vec!["--token-file", &t, "--listen", "not an address"],
        ),
        ("missing file", vec!["--token-file", "/no/such/token"]),
    ] {
        assert!(serve_options(&args(&list), no_env).is_err(), "{what}");
    }
    let empty = token_file("cli-empty", "  \n");
    assert!(serve_options(
        &args(&["--token-file", &empty.display().to_string()]),
        no_env
    )
    .err()
    .expect("an error")
    .contains("empty"));
}

#[test]
fn a_new_token_is_long_hex_and_never_the_same_twice() {
    let a = random_token();
    let b = random_token();
    assert_eq!(a.len(), 64);
    assert!(a.bytes().all(|c| c.is_ascii_hexdigit()));
    assert_ne!(a, b);
    assert!(a.len() >= crate::remote::MIN_TOKEN);
}

#[test]
fn the_process_of_a_job_reads_its_request_and_writes_its_result() {
    let dir = work_folder("cli-run-job").expect("a folder");
    let request = Request::Align {
        seqs: crate::sample::tree_demo()
            .into_iter()
            .map(|d| (d.name, d.seq))
            .collect(),
        mafft: false,
    };
    let bytes = request
        .encode()
        .with("job", "kind=align\nlabel=test\n")
        .encode();
    std::fs::write(dir.join("request.bin"), bytes).expect("a request");
    assert_eq!(run_job(&[dir.display().to_string()]), 0);
    let result = std::fs::read(dir.join("result.bin")).expect("a result");
    let Response::Align(a) = crate::remote::tasks::Response::decode(
        Kind::Align,
        &Container::decode(&result).expect("a container"),
    )
    .expect("a response") else {
        panic!("an alignment");
    };
    assert_eq!(a.rows.len(), 6);
    assert!(!dir.join("error.txt").exists());
}

#[test]
fn a_job_that_cannot_run_says_why_in_a_file_for_the_server() {
    let dir = work_folder("cli-run-bad").expect("a folder");
    // no request at all
    assert_eq!(run_job(&[dir.display().to_string()]), 1);
    assert!(std::fs::read_to_string(dir.join("error.txt"))
        .expect("an error")
        .contains("request"));
    // a request of a kind nobody runs, and one that is damaged
    for body in [&b"@job 15\nkind=nonsense\n\n"[..], b"garbage"] {
        std::fs::write(dir.join("request.bin"), body).expect("a request");
        assert_eq!(run_job(&[dir.display().to_string()]), 1);
        assert!(!dir.join("result.bin").exists());
    }
    assert_eq!(run_job(&[]), 2);
}

#[test]
fn the_commands_of_the_server_are_told_from_the_apps_own_arguments() {
    assert_eq!(dispatch(&args(&["helixloom"])), None);
    assert_eq!(
        dispatch(&args(&["helixloom", "a.gb", "--view", "home"])),
        None
    );
    assert_eq!(
        dispatch(&args(&["helixloom", "--screenshot", "x.png"])),
        None
    );
    assert_eq!(dispatch(&args(&["helixloom", "token"])), Some(0));
    assert_eq!(dispatch(&args(&["helixloom", "help"])), Some(0));
    assert_eq!(dispatch(&args(&["helixloom", "serve", "--help"])), Some(0));
    assert_eq!(dispatch(&args(&["helixloom", "status"])), Some(2));
}
