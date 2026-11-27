use super::*;
use crate::seq::minimap::sim::{genome, read, Rng};

fn args(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| s.to_string()).collect()
}

#[test]
fn every_option_is_read_and_the_defaults_are_the_documented_ones() {
    let d = parse(&args(&["a.fq"])).unwrap();
    assert_eq!(d.params, Params::default());
    assert!(d.output.is_none());
    let a = parse(&args(&[
        "-o",
        "out.fa",
        "--rounds",
        "3",
        "--coverage",
        "20",
        "--candidates",
        "90",
        "--min-support",
        "6",
        "--min-length",
        "100",
        "-k",
        "13",
        "-w",
        "4",
        "a.fq",
        "b.fa",
    ]))
    .unwrap();
    assert_eq!(a.output.as_deref(), Some("out.fa"));
    assert_eq!(a.inputs, ["a.fq", "b.fa"]);
    let p = a.params;
    assert_eq!(
        (
            p.rounds,
            p.coverage,
            p.max_candidates,
            p.min_support,
            p.min_piece,
            p.k,
            p.w
        ),
        (3, 20, 90, 6, 100, 13, 4)
    );
}

#[test]
fn options_that_are_wrong_say_what_is_wrong() {
    let bad = |a: &[&str]| parse(&args(a)).err().unwrap();
    assert!(bad(&["-o", "x"]).contains("Name the reads"));
    assert!(bad(&["--nope", "a"]).contains("--nope is not an option"));
    assert!(bad(&["a", "--rounds"]).contains("--rounds needs a value"));
    assert!(bad(&["--coverage", "x", "a"]).contains("--coverage needs a whole number"));
    assert!(bad(&["--rounds", "0", "a"]).contains("--rounds must be from 1 to 5"));
    assert!(bad(&["--rounds", "6", "a"]).contains("--rounds must be from 1 to 5"));
    assert!(bad(&["--candidates", "0", "a"]).contains("at least 1"));
    assert!(bad(&["--min-support", "0", "a"]).contains("at least 1"));
    assert!(bad(&["-k", "33", "a"]).contains("k must be from 1 to 32"));
    assert!(bad(&["-k", "0", "a"]).contains("k must be from 1 to 32"));
    assert!(bad(&["-w", "0", "a"]).contains("w must be 1 or more"));
    for ok in [
        &["--rounds", "5", "a"][..],
        &["--rounds", "1", "a"],
        &["-k", "32", "a"],
        &["--coverage", "0", "--min-length", "0", "a"],
    ] {
        assert!(parse(&args(ok)).is_ok(), "{ok:?}");
    }
}

#[test]
fn fasta_is_written_in_lines_of_eighty() {
    let r = Read {
        name: "a b".into(),
        seq: vec![b'A'; 170],
        qual: Vec::new(),
    };
    let text = fasta(&[r]);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 4);
    assert_eq!(lines[0], ">a b");
    assert_eq!(
        (lines[1].len(), lines[2].len(), lines[3].len()),
        (80, 80, 10)
    );
}

#[test]
fn the_summary_counts_reads_and_bases() {
    let r = |n: usize| Read {
        name: "r".into(),
        seq: vec![b'C'; n],
        qual: Vec::new(),
    };
    assert_eq!(
        summary(&[r(10), r(20)], &[r(25)]),
        "2 reads of 30 bases in, 1 of 25 bases out.\n"
    );
}

#[test]
fn a_run_corrects_the_reads_of_a_file_and_writes_them_where_it_is_told() {
    let dir = std::env::temp_dir().join(format!("helixloom-correct-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let g = genome(5000, 105);
    let mut r = Rng::new(9);
    let mut text = String::new();
    for i in 0..60 {
        let start = r.below(5000 - 1500);
        let s = read(&g, start, 1500, i % 2 == 0, 0.03, &mut r);
        text += &format!(
            "@r{i}\n{}\n+\n{}\n",
            String::from_utf8_lossy(&s),
            "0".repeat(s.len())
        );
    }
    let p = |n: &str| dir.join(n).display().to_string();
    std::fs::write(p("in.fq"), text).unwrap();
    assert_eq!(run(&args(&["-o", &p("out.fa"), &p("in.fq")])), 0);
    let out = std::fs::read_to_string(p("out.fa")).unwrap();
    assert!(out.starts_with(">r"));
    assert!(out.lines().filter(|l| l.starts_with('>')).count() > 40);
    // without -o the answer goes to the standard output
    assert_eq!(run(&args(&["--rounds", "1", &p("in.fq")])), 0);
    assert_eq!(run(&args(&["-o", "/no/such/folder/x", &p("in.fq")])), 1);
    assert_eq!(run(&args(&[&p("none.fq")])), 1);
    assert_eq!(run(&[]), 0);
    assert_eq!(run(&args(&["--help"])), 0);
    std::fs::remove_dir_all(&dir).ok();
}
