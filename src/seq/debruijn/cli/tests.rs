use super::*;
use crate::seq::debruijn::sim::{random, reads};
use crate::seq::dna::revcom;

fn args(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| s.to_string()).collect()
}

fn folder(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("helixloom-dbg-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn fastq(reads: &[Vec<u8>]) -> String {
    reads
        .iter()
        .enumerate()
        .map(|(i, r)| {
            format!(
                "@r{i}\n{}\n+\n{}\n",
                String::from_utf8_lossy(r),
                "I".repeat(r.len())
            )
        })
        .collect()
}

#[test]
fn the_defaults_are_those_of_the_help() {
    let a = parse(&args(&["reads.fq"])).unwrap();
    assert_eq!(a.files, ["reads.fq"]);
    assert_eq!(a.options, Options::default());
    assert_eq!((a.min_qual, a.min_length), (10, 200));
    assert_eq!((a.fasta, a.gfa, a.quiet), (None, None, false));
}

#[test]
fn every_option_is_read() {
    let a = parse(&args(&[
        "-1",
        "a_1.fq",
        "-2",
        "a_2.fq",
        "-s",
        "u.fq",
        "x.fa",
        "-k",
        "21, 33,55",
        "--cov-cutoff",
        "3.5",
        "--min-count",
        "2",
        "--tip-length",
        "80",
        "--tip-ratio",
        "4",
        "--bulge-length",
        "90",
        "--bulge-diff",
        "0.2",
        "--min-qual",
        "20",
        "--min-length",
        "500",
        "--isolated-length",
        "50",
        "--fasta",
        "c.fa",
        "--gfa",
        "g.gfa",
        "-q",
    ]))
    .unwrap();
    assert_eq!(a.files, ["a_1.fq", "a_2.fq", "u.fq", "x.fa"]);
    assert_eq!(a.options.ks, [21, 33, 55]);
    assert_eq!(a.options.cutoff, Cutoff::Value(3.5));
    assert_eq!(
        (
            a.options.min_count,
            a.options.tip_length,
            a.options.bulge_length
        ),
        (2, Some(80), Some(90))
    );
    assert_eq!((a.options.tip_ratio, a.options.bulge_diff), (4.0, 0.2));
    assert_eq!((a.min_qual, a.min_length, a.quiet), (20, 500, true));
    assert_eq!(a.options.isolated_length, 50);
    assert_eq!(
        (a.fasta.as_deref(), a.gfa.as_deref()),
        (Some("c.fa"), Some("g.gfa"))
    );
    let b = parse(&args(&["r.fq", "--cov-cutoff", "off"])).unwrap();
    assert_eq!(b.options.cutoff, Cutoff::Off);
    let c = parse(&args(&["r.fq", "--cov-cutoff", "auto"])).unwrap();
    assert_eq!(c.options.cutoff, Cutoff::Auto);
    let d = parse(&args(&["r.fq", "-o", "out"])).unwrap();
    assert_eq!(
        (d.fasta.as_deref(), d.gfa.as_deref()),
        (Some("out/contigs.fasta"), Some("out/assembly_graph.gfa"))
    );
}

#[test]
fn a_wrong_option_says_so() {
    let why = |a: &[&str]| parse(&args(a)).err().unwrap();
    assert!(why(&[]).contains("Name the FASTA or FASTQ"));
    assert!(why(&["-k", "21"]).contains("Name the FASTA or FASTQ"));
    assert!(why(&["r.fq", "--blast"]).contains("--blast is not an option"));
    assert!(why(&["r.fq", "-k", "x"]).contains("-k needs a number, not x"));
    assert!(why(&["r.fq", "-k", "21,,33"]).contains("-k needs a number"));
    assert!(why(&["r.fq", "-k"]).contains("-k needs a value"));
    assert!(why(&["r.fq", "--cov-cutoff", "-2"]).contains("0 or more"));
    assert!(why(&["r.fq", "--cov-cutoff", "lots"]).contains("--cov-cutoff needs a number"));
    assert!(why(&["r.fq", "--tip-ratio", "nan"]).contains("0 or more"));
    assert!(why(&["r.fq", "--min-count", "0"]).contains("1 or more"));
    assert!(why(&["r.fq", "--min-qual", "300"]).contains("--min-qual needs a number"));
}

#[test]
fn reads_are_assembled_to_a_fasta_and_a_gfa_in_a_folder() {
    let d = folder("run");
    let genome = random(3000, 21);
    let r = reads(&genome, false, 100, 40.0, 0.002, 22);
    let input = d.join("reads.fq");
    std::fs::write(&input, fastq(&r)).unwrap();
    let out = d.join("out");
    let code = run(&args(&[
        "-o",
        out.to_str().unwrap(),
        "-q",
        input.to_str().unwrap(),
    ]));
    assert_eq!(code, 0);
    let fa = std::fs::read_to_string(out.join("contigs.fasta")).unwrap();
    let mut lines = fa.lines();
    let header = lines.next().unwrap();
    assert!(header.starts_with(">NODE_1_length_3"), "{header}");
    let seq: String = lines.collect();
    assert!(
        seq.as_bytes() == genome || seq.as_bytes() == revcom(&genome),
        "contig is not the genome"
    );
    let gfa = std::fs::read_to_string(out.join("assembly_graph.gfa")).unwrap();
    assert!(gfa.starts_with("H\tVN:Z:1.0\nS\t1\t"));
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn the_two_files_of_a_pair_and_a_k_that_is_asked_for_are_used() {
    let d = folder("pair");
    let genome = random(2000, 23);
    let r = reads(&genome, false, 100, 40.0, 0.0, 24);
    let (a, b) = r.split_at(r.len() / 2);
    let (f1, f2) = (d.join("a_1.fq"), d.join("a_2.fq"));
    std::fs::write(&f1, fastq(a)).unwrap();
    std::fs::write(&f2, fastq(b)).unwrap();
    let out = d.join("contigs.fa");
    let code = run(&args(&[
        "-1",
        f1.to_str().unwrap(),
        "-2",
        f2.to_str().unwrap(),
        "-k",
        "31",
        "--min-length",
        "1000",
        "--fasta",
        out.to_str().unwrap(),
        "-q",
    ]));
    assert_eq!(code, 0);
    let fa = std::fs::read_to_string(&out).unwrap();
    assert_eq!(fa.matches('>').count(), 1);
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn what_cannot_be_run_says_why_and_the_help_is_the_help() {
    let d = folder("errors");
    assert_eq!(run(&args(&["--help"])), 0);
    assert_eq!(run(&[]), 0);
    let why = |a: &[&str]| run_inner(&args(a)).err().unwrap();
    assert!(why(&["nothing.fq"]).contains("Could not read nothing.fq"));
    // reads that are all unknown bases
    let n = d.join("n.fa");
    std::fs::write(&n, format!(">r\n{}\n", "N".repeat(100))).unwrap();
    assert!(why(&[n.to_str().unwrap()]).contains("No run of 21 bases or more is left"));
    // reads shorter than every default k
    let short = d.join("s.fa");
    std::fs::write(
        &short,
        format!(">r\n{}\n", String::from_utf8_lossy(&random(40, 1))),
    )
    .unwrap();
    assert_eq!(run(&args(&["-q", short.to_str().unwrap(), "-k", "45"])), 1);
    // a k that is even
    assert!(why(&[short.to_str().unwrap(), "-k", "22"]).contains("odd"));
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn a_cutoff_and_a_ratio_of_zero_are_allowed() {
    let a = parse(&args(&[
        "r.fq",
        "--cov-cutoff",
        "0",
        "--tip-ratio",
        "0",
        "--bulge-diff",
        "0",
    ]))
    .unwrap();
    assert_eq!(a.options.cutoff, Cutoff::Value(0.0));
    assert_eq!((a.options.tip_ratio, a.options.bulge_diff), (0.0, 0.0));
}
