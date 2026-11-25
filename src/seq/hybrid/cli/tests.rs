use super::*;
use crate::seq::debruijn::output::{gfa, number};
use crate::seq::hybrid::sim::{is_circle_of, long_reads, ring, short_graph};

fn args(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| s.to_string()).collect()
}

fn folder(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("helixloom-hybrid-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn fasta_of(reads: &[(String, Vec<u8>)]) -> String {
    reads
        .iter()
        .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
        .collect()
}

#[test]
fn the_defaults_are_those_of_the_help() {
    let a = parse(&args(&["g.gfa", "long.fq"])).unwrap();
    assert_eq!(
        (a.graph.as_str(), a.reads.as_slice()),
        ("g.gfa", ["long.fq".to_string()].as_slice())
    );
    assert_eq!(a.options, Options::default());
    assert_eq!(
        (a.genes, a.overlap, a.fasta, a.gfa, a.quiet),
        (None, 0, None, None, false)
    );
}

#[test]
fn every_option_is_read() {
    let a = parse(&args(&[
        "--mode",
        "bold",
        "--min-quality",
        "40",
        "--min-reads",
        "3",
        "--single-depth",
        "12.5",
        "--depth-tolerance",
        "0.1",
        "--min-align",
        "500",
        "--end-slack",
        "100",
        "--start-genes",
        "genes.fa",
        "--overlap",
        "54",
        "--fasta",
        "a.fa",
        "--gfa",
        "a.gfa",
        "-q",
        "g.gfa",
        "x.fa",
        "y.fq",
    ]))
    .unwrap();
    assert_eq!(a.options.mode, Mode::Bold);
    assert_eq!(
        (a.options.min_quality, a.options.min_reads),
        (Some(40.0), 3)
    );
    assert_eq!(
        (a.options.copies.single_depth, a.options.copies.tolerance),
        (Some(12.5), 0.1)
    );
    assert_eq!(
        (a.options.bridges.min_align, a.options.bridges.end_slack),
        (500, 100)
    );
    assert_eq!((a.genes.as_deref(), a.overlap), (Some("genes.fa"), 54));
    assert_eq!(
        (a.fasta.as_deref(), a.gfa.as_deref(), a.quiet),
        (Some("a.fa"), Some("a.gfa"), true)
    );
    assert_eq!(a.reads, ["x.fa", "y.fq"]);
    let b = parse(&args(&[
        "-o",
        "out",
        "g.gfa",
        "r.fq",
        "--mode",
        "conservative",
    ]))
    .unwrap();
    assert_eq!(b.options.mode, Mode::Conservative);
    assert_eq!(
        (b.fasta.as_deref(), b.gfa.as_deref()),
        (Some("out/assembly.fasta"), Some("out/assembly.gfa"))
    );
    assert_eq!(
        parse(&args(&["g", "r", "--mode", "normal"]))
            .unwrap()
            .options
            .mode,
        Mode::Normal
    );
}

#[test]
fn a_wrong_option_says_so() {
    let why = |a: &[&str]| parse(&args(a)).err().unwrap();
    assert!(why(&[]).contains("Name the GFA"));
    assert!(why(&["g.gfa"]).contains("Name the FASTA or FASTQ file of the long reads"));
    assert!(why(&["g", "r", "--blast"]).contains("--blast is not an option"));
    assert!(why(&["g", "r", "--mode", "fast"]).contains("conservative, normal or bold, not fast"));
    assert!(why(&["g", "r", "--min-quality", "101"]).contains("from 0 to 100"));
    assert!(why(&["g", "r", "--min-quality", "x"]).contains("--min-quality needs a number, not x"));
    assert!(why(&["g", "r", "--min-quality", "-1"]).contains("0 or more"));
    assert!(why(&["g", "r", "--min-reads", "0"]).contains("1 or more"));
    assert!(why(&["g", "r", "--single-depth", "0"]).contains("more than 0"));
    assert!(why(&["g", "r", "--min-align"]).contains("--min-align needs a value"));
}

#[test]
fn zero_and_a_hundred_are_allowed_where_they_mean_something() {
    let a = parse(&args(&["g", "r", "--min-quality", "100"])).unwrap();
    assert_eq!(a.options.min_quality, Some(100.0));
    let a = parse(&args(&[
        "g",
        "r",
        "--min-quality",
        "0",
        "--depth-tolerance",
        "0",
        "--min-align",
        "0",
        "--end-slack",
        "0",
    ]))
    .unwrap();
    assert_eq!(a.options.min_quality, Some(0.0));
    assert_eq!(a.options.copies.tolerance, 0.0);
    assert_eq!(
        (a.options.bridges.min_align, a.options.bridges.end_slack),
        (0, 0)
    );
}

#[test]
fn a_graph_and_long_reads_are_joined_to_a_circle_in_two_files() {
    let d = folder("run");
    let r = ring(3000, 1200, 1);
    let g = short_graph(&r.genome, true, 2);
    let graph_path = d.join("graph.gfa");
    std::fs::write(&graph_path, gfa(&g, &number(&g))).unwrap();
    let reads = long_reads(&r.genome, true, 5000, 40.0, 0.02, 3);
    let reads_path = d.join("long.fa");
    std::fs::write(&reads_path, fasta_of(&reads)).unwrap();
    let out = d.join("out");
    let code = run(&args(&[
        "-o",
        out.to_str().unwrap(),
        "-q",
        graph_path.to_str().unwrap(),
        reads_path.to_str().unwrap(),
    ]));
    assert_eq!(code, 0);
    let fa = std::fs::read_to_string(out.join("assembly.fasta")).unwrap();
    let mut lines = fa.lines();
    let header = lines.next().unwrap();
    assert!(
        header.starts_with(">1 length=11400 depth=") && header.ends_with("x circular=true"),
        "{header}"
    );
    let seq: String = lines.collect();
    assert!(is_circle_of(seq.as_bytes(), &r.genome));
    let out_gfa = std::fs::read_to_string(out.join("assembly.gfa")).unwrap();
    assert!(out_gfa.starts_with("H\tVN:Z:1.0\nS\t1\t"));
    assert_eq!(out_gfa.lines().filter(|l| l.starts_with('S')).count(), 1);
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn a_start_gene_turns_the_circle() {
    let d = folder("gene");
    let r = ring(3000, 1200, 1);
    let g = short_graph(&r.genome, true, 2);
    let graph_path = d.join("graph.gfa");
    std::fs::write(&graph_path, gfa(&g, &number(&g))).unwrap();
    let reads_path = d.join("long.fa");
    std::fs::write(
        &reads_path,
        fasta_of(&long_reads(&r.genome, true, 5000, 40.0, 0.02, 3)),
    )
    .unwrap();
    let gene = r.genome[5000..5400].to_vec();
    let genes_path = d.join("genes.fa");
    std::fs::write(
        &genes_path,
        format!(">dnaA\n{}\n", String::from_utf8_lossy(&gene)),
    )
    .unwrap();
    let fa_path = d.join("a.fa");
    let code = run(&args(&[
        "--start-genes",
        genes_path.to_str().unwrap(),
        "--fasta",
        fa_path.to_str().unwrap(),
        "-q",
        graph_path.to_str().unwrap(),
        reads_path.to_str().unwrap(),
    ]));
    assert_eq!(code, 0);
    let fa = std::fs::read_to_string(&fa_path).unwrap();
    let seq: String = fa.lines().skip(1).collect();
    let start = String::from_utf8_lossy(&gene).to_string();
    // the circle is either strand of the genome: the gene is first, on its own strand
    assert!(
        seq.starts_with(&start)
            || seq
                .starts_with(&String::from_utf8_lossy(&crate::seq::dna::revcom(&gene)).to_string())
    );
    assert!(
        seq.starts_with(&start),
        "the gene must be first on the strand it lies on"
    );
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn what_cannot_be_run_says_why_and_the_help_is_the_help() {
    let d = folder("errors");
    assert_eq!(run(&args(&["--help"])), 0);
    assert_eq!(run(&[]), 0);
    let why = |a: &[&str]| run_inner(&args(a)).err().unwrap();
    assert!(why(&["nothing.gfa", "r.fa"]).contains("Could not read nothing.gfa"));
    let g = d.join("g.gfa");
    std::fs::write(&g, "S\tA\t*\n").unwrap();
    assert!(why(&[g.to_str().unwrap(), "r.fa"]).contains("no bases"));
    std::fs::write(&g, "S\tA\tAAACCCGGG\tDP:f:5\n").unwrap();
    assert!(why(&[g.to_str().unwrap(), "missing.fa"]).contains("Could not read missing.fa"));
    std::fs::remove_dir_all(&d).ok();
}
