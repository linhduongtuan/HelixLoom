use super::*;
use crate::seq::fastq::{to_fastq, FastqRead};
use crate::seq::minimap::sim::genome;
use crate::seq::variants::sim::{plain, plant, reads as simulated};

fn args(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| s.to_string()).collect()
}

fn folder(name: &str) -> std::path::PathBuf {
    let dir =
        std::env::temp_dir().join(format!("helixloom-variants-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn path(dir: &std::path::Path, name: &str) -> String {
    dir.join(name).display().to_string()
}

fn fasta(name: &str, seq: &[u8]) -> String {
    format!(">{name}\n{}\n", String::from_utf8_lossy(seq))
}

#[test]
fn every_option_is_read_and_the_defaults_are_those_of_the_usage() {
    let d = parse(&args(&["ref.fa", "reads.fq"])).unwrap();
    assert_eq!(
        (d.reference.as_str(), d.reads.as_str()),
        ("ref.fa", "reads.fq")
    );
    assert_eq!(d.filter, Filter::default());
    assert_eq!(d.params, Params::default());
    assert_eq!((d.mapping.k, d.mapping.w, d.mapping.cigar), (15, 10, true));
    assert_eq!(d.sample, "sample");
    assert!(d.output.is_none() && d.consensus.is_none() && d.depth.is_none() && d.sam.is_none());
    let a = parse(&args(&[
        "-o",
        "v.vcf",
        "--consensus",
        "c.fa",
        "--depth",
        "d.tsv",
        "--sam",
        "m.sam",
        "--sample",
        "S1",
        "-q",
        "7",
        "-Q",
        "20",
        "-d",
        "100",
        "-D",
        "8",
        "--min-qual",
        "30",
        "--prior",
        "0.01",
        "--indel-prior",
        "0.001",
        "--indel-error",
        "0.2",
        "--default-quality",
        "12",
        "-k",
        "17",
        "-w",
        "5",
        "ref.fa",
        "reads.fq",
    ]))
    .unwrap();
    assert_eq!(
        (
            a.output.as_deref(),
            a.consensus.as_deref(),
            a.depth.as_deref(),
            a.sam.as_deref(),
            a.sample.as_str()
        ),
        (
            Some("v.vcf"),
            Some("c.fa"),
            Some("d.tsv"),
            Some("m.sam"),
            "S1"
        )
    );
    assert_eq!(
        (
            a.filter.min_mapq,
            a.filter.min_base_quality,
            a.filter.max_depth,
            a.filter.default_quality
        ),
        (7, 20, 100, 12)
    );
    assert_eq!(
        (
            a.params.min_depth,
            a.params.min_quality,
            a.params.snp_prior,
            a.params.indel_prior,
            a.params.indel_error
        ),
        (8, 30.0, 0.01, 0.001, 0.2)
    );
    assert_eq!((a.mapping.k, a.mapping.w), (17, 5));
}

#[test]
fn options_that_are_wrong_say_what_is_wrong() {
    let bad = |a: &[&str]| parse(&args(a)).err().unwrap();
    assert!(bad(&["ref.fa"]).contains("Name a reference and the reads"));
    assert!(bad(&["a", "b", "c"]).contains("Name a reference and the reads"));
    assert!(bad(&["--nope", "a", "b"]).contains("--nope is not an option"));
    assert!(bad(&["-Q", "x", "a", "b"]).contains("-Q needs a number, not x"));
    assert!(bad(&["a", "b", "-o"]).contains("-o needs a value"));
    assert!(bad(&["-k", "0", "a", "b"]).contains("k must be from 1 to 32"));
    assert!(bad(&["-k", "33", "a", "b"]).contains("k must be from 1 to 32"));
    assert!(bad(&["-w", "0", "a", "b"]).contains("w must be 1 or more"));
    assert!(bad(&["-d", "0", "a", "b"]).contains("-d must be from 1 to"));
    assert!(bad(&["-d", "65536", "a", "b"]).contains("-d must be from 1 to"));
    assert!(bad(&["--default-quality", "94", "a", "b"]).contains("default quality"));
    assert!(bad(&["--prior", "0", "a", "b"]).contains("priors"));
    assert!(bad(&["-D", "0", "a", "b"]).contains("depth"));
    // the ends that are fine
    assert!(parse(&args(&["-d", "65535", "--default-quality", "93", "a", "b"])).is_ok());
    assert_eq!(parse(&args(&["-k", "32", "a", "b"])).unwrap().mapping.k, 32);
    assert_eq!(parse(&args(&["-k", "1", "a", "b"])).unwrap().mapping.k, 1);
}

#[test]
fn a_sam_file_is_told_by_its_name_or_by_a_header_line_and_not_by_a_fastq_name() {
    assert!(is_sam("x.SAM", "anything"));
    assert!(is_sam("x.txt", "@HD\tVN:1.6\n"));
    assert!(is_sam("x.txt", "@SQ\tSN:a\tLN:4\n"));
    assert!(is_sam("x.txt", "@PG\tID:x\n"));
    assert!(is_sam("x.txt", "@RG\tID:x\n"));
    assert!(is_sam("x.txt", "@CO\tnote\n"));
    assert!(!is_sam("x.fq", "@HD\nACGT\n+\nIIII\n"));
    assert!(!is_sam("x.fq", "@read1\nACGT\n+\nIIII\n"));
    assert!(!is_sam("x.fa", ">a\nACGT\n"));
}

#[test]
fn a_reference_with_the_same_name_twice_or_an_empty_sequence_is_refused() {
    let d = folder("reference");
    std::fs::write(d.join("dup.fa"), ">a\nACGT\n>a\nTTTT\n").unwrap();
    std::fs::write(d.join("empty.fa"), ">a\nACGT\n>b\n").unwrap();
    std::fs::write(d.join("ok.fa"), ">a x\nACGT\n>b\nTTTT\n").unwrap();
    let e = reference(&path(&d, "dup.fa")).unwrap_err();
    assert!(e.contains("named a"), "{e}");
    let e = reference(&path(&d, "empty.fa")).unwrap_err();
    assert!(e.contains("b has no bases"), "{e}");
    let ok = reference(&path(&d, "ok.fa")).unwrap();
    assert_eq!(
        ok,
        [
            ("a".to_string(), b"ACGT".to_vec()),
            ("b".to_string(), b"TTTT".to_vec())
        ]
    );
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn with_no_arguments_or_help_the_usage_is_written_and_the_answer_is_nothing_wrong() {
    assert_eq!(run(&[]), 0);
    assert_eq!(run(&args(&["--help"])), 0);
    assert_eq!(run(&args(&["-h"])), 0);
    assert!(USAGE.contains("helixloom variants"));
    assert_eq!(run(&args(&["--nope"])), 1);
    assert_eq!(run(&args(&["/no/such/ref.fa", "/no/such/reads.fq"])), 1);
}

/// A sample with a substitution, a deletion of two bases and an insertion of one, and reads from it with the errors of nanopore reads.
struct Case {
    dir: std::path::PathBuf,
    sample: Vec<u8>,
    snp: usize,
    deleted: usize,
    inserted: usize,
}

fn case(name: &str, seed: u32) -> Case {
    let dir = folder(name);
    let reference = genome(3000, seed);
    let deleted = plain(&reference, 1500);
    let inserted = plain(&reference, 2200);
    let base = *b"ACGT"
        .iter()
        .find(|&&c| c != reference[inserted] && c != reference[inserted + 1])
        .unwrap();
    let snp = 800;
    let sample = plant(&reference, &[snp], &[(deleted, 2)], &[(inserted, &[base])]);
    let reads = simulated(&sample, 90, 1500, 0.02, 15, seed + 1);
    let fastq: Vec<FastqRead> = reads
        .into_iter()
        .map(|r| FastqRead {
            name: r.name,
            seq: r.seq,
            qual: r.qual,
        })
        .collect();
    std::fs::write(dir.join("ref.fa"), fasta("ref", &reference)).unwrap();
    std::fs::write(dir.join("reads.fq"), to_fastq(&fastq)).unwrap();
    Case {
        dir,
        sample,
        snp,
        deleted,
        inserted,
    }
}

#[test]
fn reads_of_a_sample_with_three_planted_differences_give_those_and_the_consensus_of_the_sample() {
    let c = case("planted", 51);
    let p = |n: &str| path(&c.dir, n);
    let code = run(&args(&[
        "-o",
        &p("v.vcf"),
        "--consensus",
        &p("c.fa"),
        "--depth",
        &p("d.tsv"),
        "--sam",
        &p("m.sam"),
        &p("ref.fa"),
        &p("reads.fq"),
    ]));
    assert_eq!(code, 0);
    let vcf = std::fs::read_to_string(p("v.vcf")).unwrap();
    let rows: Vec<Vec<&str>> = vcf
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect())
        .collect();
    let found: Vec<(usize, usize, usize)> = rows
        .iter()
        .map(|r| (r[1].parse().unwrap(), r[3].len(), r[4].len()))
        .collect();
    assert_eq!(
        found,
        [(c.snp + 1, 1, 1), (c.deleted, 3, 1), (c.inserted + 1, 1, 2)],
        "{vcf}"
    );
    // the consensus is the sample, but for the ends where few reads reach
    let cons = std::fs::read_to_string(p("c.fa")).unwrap();
    let bases: Vec<u8> = cons.lines().skip(1).flat_map(|l| l.bytes()).collect();
    assert_eq!(bases.len(), c.sample.len());
    assert_eq!(
        &bases[100..bases.len() - 100],
        &c.sample[100..c.sample.len() - 100]
    );
    // the depth table has a row to a position, and the middle of the genome is well covered
    let depth = std::fs::read_to_string(p("d.tsv")).unwrap();
    assert_eq!(depth.lines().count(), 3001);
    let mid: u32 = depth
        .lines()
        .nth(1500)
        .unwrap()
        .split('\t')
        .nth(2)
        .unwrap()
        .parse()
        .unwrap();
    assert!(mid >= 20, "{mid}");
    // and the SAM that was written is one that this program can read
    let sam = std::fs::read_to_string(p("m.sam")).unwrap();
    assert!(sam.starts_with("@HD\tVN:1.6"));
    assert!(sam.contains("@SQ\tSN:ref\tLN:3000"));
    assert_eq!(
        run(&args(&["-o", &p("again.vcf"), &p("ref.fa"), &p("m.sam")])),
        0
    );
    assert_eq!(std::fs::read_to_string(p("again.vcf")).unwrap(), vcf);
    std::fs::remove_dir_all(&c.dir).ok();
}

#[test]
fn a_sam_file_that_does_not_fit_the_reference_is_refused_and_so_is_a_second_mapping_asked_of_it() {
    let d = folder("badsam");
    let bases = genome(100, 5);
    std::fs::write(d.join("ref.fa"), fasta("ref", &bases)).unwrap();
    let write = |name: &str, text: &str| std::fs::write(d.join(name), text).unwrap();
    write("other.sam", "@SQ\tSN:zzz\tLN:100\n");
    write("short.sam", "@SQ\tSN:ref\tLN:99\n");
    write("lost.sam", "a\t0\tzzz\t1\t60\t4M\t*\t0\t0\tACGT\t*\n");
    write(
        "ok.sam",
        "@SQ\tSN:ref\tLN:100\na\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t*\n",
    );
    let run_with = |reads: &str, extra: &[&str]| {
        let mut a = vec!["-o", "/dev/null"];
        a.extend_from_slice(extra);
        let (r, s) = (path(&d, "ref.fa"), path(&d, reads));
        a.push(&r);
        a.push(&s);
        let args = args(&a);
        parse(&args).and_then(|a| {
            let targets = reference(&a.reference)?;
            alignments(&a, &targets).map(|_| ())
        })
    };
    let e = run_with("other.sam", &[]).unwrap_err();
    assert!(e.contains("target zzz") && e.contains("does not"), "{e}");
    let e = run_with("short.sam", &[]).unwrap_err();
    assert!(e.contains("length of 99") && e.contains("100"), "{e}");
    let e = run_with("lost.sam", &[]).unwrap_err();
    assert!(e.contains("lies on zzz"), "{e}");
    assert!(run_with("ok.sam", &[]).is_ok());
    let e = run_with("ok.sam", &["--sam", "x.sam"]).unwrap_err();
    assert!(e.contains("SAM file already"), "{e}");
    let e = run_with("missing.sam", &[]).unwrap_err();
    assert!(e.contains("Could not read"), "{e}");
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn the_vcf_goes_to_the_standard_output_when_no_file_is_named_and_a_file_that_cannot_be_written_is_an_error(
) {
    let c = case("stdout", 52);
    let p = |n: &str| path(&c.dir, n);
    assert_eq!(run(&args(&[&p("ref.fa"), &p("reads.fq")])), 0);
    assert_eq!(
        run(&args(&[
            "-o",
            "/no/such/folder/v.vcf",
            &p("ref.fa"),
            &p("reads.fq")
        ])),
        1
    );
    std::fs::remove_dir_all(&c.dir).ok();
}
