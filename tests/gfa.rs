//! Inspect graphs through the command's file and standard-input interfaces.

use std::io::Write;
use std::process::{Command, Output, Stdio};

const BIN: &str = env!("CARGO_BIN_EXE_helixloom");
const TABLE_HEADER: &str = "segment\tlength_bp\tdepth\tkmer_count\tlink_records\n";

fn run_stdin(command: &str, gfa: &str) -> Output {
    let mut child = Command::new(BIN)
        .args(["gfa", command, "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(gfa.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

fn assert_table(out: Output, rows: &str) {
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        format!("{TABLE_HEADER}{rows}")
    );
}

fn assert_summary(out: Output, rows: &[&str]) {
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = String::from_utf8(out.stdout).unwrap();
    assert_eq!(
        text.lines().count(),
        rows.len(),
        "unexpected rows in {text}"
    );
    for row in rows {
        assert!(
            text.lines().any(|line| line == *row),
            "missing {row:?} in {text}"
        );
    }
}

#[test]
fn help_lists_the_command_and_explains_its_inputs() {
    let run = |args: &[&str]| {
        let out = Command::new(BIN).args(args).output().unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        String::from_utf8(out.stdout).unwrap()
    };
    let general = run(&["--help"]);
    let gfa_help = run(&["gfa", "--help"]);
    for command in ["inspect", "segments"] {
        let file_form = format!("gfa {command} FILE");
        let stdin_form = format!("gfa {command} -");
        assert!(general.contains(&file_form), "{general}");
        assert!(
            gfa_help.contains(&file_form) && gfa_help.contains(&stdin_form),
            "{gfa_help}"
        );
        let text = run(&["gfa", command, "--help"]);
        assert!(
            text.contains(&file_form) && text.contains(&stdin_form),
            "{text}"
        );
    }
}

#[test]
fn missing_files_and_incorrect_arguments_fail_without_a_report() {
    let missing =
        std::env::temp_dir().join(format!("helixloom-missing-{}.gfa", std::process::id()));
    for (args, diagnostic) in [
        (
            vec!["gfa", "inspect", missing.to_str().unwrap()],
            missing.to_str().unwrap(),
        ),
        (vec!["gfa", "inspect"], "gfa inspect FILE"),
        (vec!["gfa", "inspect", "-", "extra"], "gfa inspect FILE"),
    ] {
        let out = Command::new(BIN).args(args).output().unwrap();
        assert!(!out.status.success());
        assert!(out.stdout.is_empty());
        assert!(String::from_utf8_lossy(&out.stderr).contains(diagnostic));
    }
}

#[test]
fn inspects_the_graph_and_reports_segment_tags_written_by_dbg() {
    let dir = std::env::temp_dir().join(format!("helixloom-gfa-dbg-{}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    std::fs::write(dir.join("reads.fa"), ">tiny\nAGCCATGGT\n").unwrap();
    let dbg = Command::new(BIN)
        .current_dir(&dir)
        .args([
            "dbg",
            "-k",
            "5",
            "--cov-cutoff",
            "off",
            "--isolated-length",
            "0",
            "-o",
            ".",
            "reads.fa",
        ])
        .output()
        .unwrap();
    let inspect = Command::new(BIN)
        .current_dir(&dir)
        .args(["gfa", "inspect", "assembly_graph.gfa"])
        .output()
        .unwrap();
    let segments = Command::new(BIN)
        .current_dir(&dir)
        .args(["gfa", "segments", "assembly_graph.gfa"])
        .output()
        .unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    assert!(
        dbg.status.success(),
        "{}",
        String::from_utf8_lossy(&dbg.stderr)
    );
    assert_summary(
        inspect,
        &[
            "segments\t1",
            "links\t1",
            "components\t1",
            "isolated_segments\t0",
            "known_segment_bases\t7",
            "segments_without_length\t0",
        ],
    );
    assert_table(segments, "1\t7\t1.333333\t4\t1\n");
}

#[test]
fn resolves_forward_links_and_distinguishes_declared_and_unknown_lengths_on_stdin() {
    let out = run_stdin(
        "inspect",
        "L\ta\t+\tb\t-\t1M\nS\ta\tACGT\nS\tb\t*\tLN:i:6\nS\tx\t*\n",
    );
    assert_summary(
        out,
        &[
            "segments\t3",
            "links\t1",
            "components\t2",
            "isolated_segments\t1",
            "known_segment_bases\t10",
            "segments_without_length\t1",
        ],
    );
}

#[test]
fn an_unresolved_link_fails_without_a_partial_report() {
    for command in ["inspect", "segments"] {
        let out = run_stdin(command, "S\ta\tACGT\nL\ta\t+\tmissing\t+\t0M\n");
        assert!(!out.status.success());
        assert!(out.stdout.is_empty());
        assert!(String::from_utf8_lossy(&out.stderr).contains("missing"));
    }
}

#[test]
fn lists_segments_in_declaration_order_with_incident_records_and_unknown_fields() {
    let gfa = "L\tz\t-\ta\t+\t*\n\
               S\tb\t*\tLN:i:0\tDP:f:0\n\
               S\ta\t*\tKC:i:0\n\
               L\tz\t-\ta\t+\t*\n\
               L\tz\t+\tz\t-\t1M\n\
               S\tz\tACTG\tDP:f:1.23456789\tKC:i:19";
    assert_table(
        run_stdin("segments", gfa),
        "b\t0\t0\t.\t0\na\t.\t.\t0\t2\nz\t4\t1.23456789\t19\t3\n",
    );
    assert_table(run_stdin("segments", "H\tVN:Z:1.0\n"), "");
}

#[test]
fn a_late_malformed_record_emits_no_table_header_or_rows() {
    let out = run_stdin("segments", "S\ta\tAC\nS\tb\t*\tDP:f:2\nS\tc\tA\tKC:i:-1");
    assert!(!out.status.success());
    assert!(out.stdout.is_empty());
    assert!(String::from_utf8_lossy(&out.stderr).contains("line 3"));
}

#[test]
fn a_pipeline_folder_writes_the_same_tsv_as_inspect_and_reuses_unchanged_files() {
    let dir = std::env::temp_dir().join(format!("helixloom-gfa-pipeline-{}", std::process::id()));
    std::fs::create_dir_all(dir.join("inputs")).unwrap();
    let cmd = |args: &[&str]| {
        Command::new(BIN)
            .args(args)
            .current_dir(&dir)
            .output()
            .unwrap()
    };
    let example = cmd(&["pipeline", "example", "gfa"]);
    assert!(example.status.success());
    std::fs::write(dir.join("graphs.txt"), example.stdout).unwrap();
    std::fs::write(dir.join("inputs/sample one.gfa"), "S\ta\tACGT\n").unwrap();
    std::fs::write(dir.join("inputs/sample.2.GFA"), "H\tVN:Z:1.0\n").unwrap();
    let run = || {
        cmd(&[
            "pipeline",
            "run",
            "graphs.txt",
            "inputs",
            "--out-dir",
            "results",
            "--cache",
            "cache",
        ])
    };
    let first = run();
    assert!(
        first.status.success(),
        "{}",
        String::from_utf8_lossy(&first.stderr)
    );
    for (input, output) in [
        ("inputs/sample one.gfa", "results/sample one.tsv"),
        ("inputs/sample.2.GFA", "results/sample.2.tsv"),
    ] {
        let direct = cmd(&["gfa", "inspect", input]);
        assert!(direct.status.success());
        assert_eq!(std::fs::read(dir.join(output)).unwrap(), direct.stdout);
    }
    std::fs::write(dir.join("inputs/sample.2.GFA"), "S\tb\t*\tLN:i:5\n").unwrap();
    std::fs::write(dir.join("inputs/bad.gfa"), "L\tmissing\t+\tmissing\t+\t*\n").unwrap();
    let second = run();
    assert!(!second.status.success());
    let log = String::from_utf8(second.stdout).unwrap();
    for expected in [
        "sample one: gfa-inspect: found in the cache",
        "sample.2: gfa-inspect: ran",
        "bad: gfa-inspect: FAILED",
        "line 1",
    ] {
        assert!(log.contains(expected), "{expected:?} missing from {log}");
    }
    assert!(!dir.join("results/bad.tsv").exists());
    assert_eq!(
        std::fs::read(dir.join("results/sample.2.tsv")).unwrap(),
        cmd(&["gfa", "inspect", "inputs/sample.2.GFA"]).stdout
    );
    std::fs::remove_dir_all(dir).unwrap();
}
