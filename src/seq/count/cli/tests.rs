use super::*;
use crate::seq::minimap::sim::genome;
use crate::seq::minimap::Options;
use crate::seq::variants::mapped::map_to_sam;
use crate::seq::variants::reads::Read;
use crate::seq::variants::sam as samfile;

fn args(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| s.to_string()).collect()
}

fn folder(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("helixloom-count-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn every_option_is_read() {
    let d = parse(&args(&["-a", "g.gtf", "x.sam"])).unwrap();
    assert_eq!(
        (
            d.annotation.as_str(),
            d.types.as_slice(),
            d.attribute.as_str()
        ),
        ("g.gtf", &["exon".to_string()][..], "gene_id")
    );
    assert_eq!(d.rules, Rules::default());
    assert!(d.output.is_none());
    let a = parse(&args(&[
        "-a",
        "g",
        "-o",
        "out",
        "-t",
        "CDS,gene,",
        "-g",
        "ID",
        "-f",
        "-O",
        "-M",
        "--fraction",
        "--primary",
        "-Q",
        "20",
        "-s",
        "2",
        "--minOverlap",
        "5",
        "--fracOverlap",
        "0.5",
        "--fracOverlapFeature",
        "0.25",
        "--nonOverlap",
        "7",
        "--nonOverlapFeature",
        "9",
        "--largestOverlap",
        "--ignoreDup",
        "--splitOnly",
        "--each-read",
        "a.sam",
        "b.sam",
    ]))
    .unwrap();
    assert_eq!(a.output.as_deref(), Some("out"));
    assert_eq!(a.types, ["CDS", "gene"]);
    assert_eq!(a.attribute, "ID");
    assert_eq!(a.files, ["a.sam", "b.sam"]);
    let r = a.rules;
    assert!(r.feature_level && r.multi_overlap && r.multi_mapping && r.fraction && r.primary_only);
    assert!(r.largest_overlap && r.ignore_duplicates && r.each_read);
    assert_eq!(
        (
            r.min_mapq,
            r.strandedness,
            r.min_overlap,
            r.frac_overlap,
            r.frac_overlap_feature
        ),
        (20, 2, 5, 0.5, 0.25)
    );
    assert_eq!(
        (r.non_overlap, r.non_overlap_feature, r.splice),
        (Some(7), Some(9), Splice::Split)
    );
    let w = parse(&args(&["-a", "g", "--nonSplitOnly", "x"])).unwrap();
    assert_eq!(w.rules.splice, Splice::Whole);
}

#[test]
fn options_that_are_wrong_say_what_is_wrong() {
    let bad = |a: &[&str]| parse(&args(a)).err().unwrap();
    assert!(bad(&["x.sam"]).contains("annotation with -a"));
    assert!(bad(&["-a", "g"]).contains("at least one SAM file"));
    assert!(bad(&["-a", "g", "--nope", "x"]).contains("--nope is not an option"));
    assert!(bad(&["-a"]).contains("-a needs a value"));
    assert!(bad(&["-a", "g", "-Q", "x", "x"]).contains("-Q needs a whole number"));
    assert!(bad(&["-a", "g", "-s", "3", "x"]).contains("-s must be 0, 1 or 2"));
    assert!(bad(&["-a", "g", "--fracOverlap", "1.1", "x"]).contains("share between 0 and 1"));
    assert!(
        bad(&["-a", "g", "--fracOverlapFeature", "-0.1", "x"]).contains("share between 0 and 1")
    );
    assert!(bad(&["-a", "g", "-t", ",", "x"]).contains("at least one type"));
    assert!(bad(&["-a", "g", "--fraction", "x"]).contains("--fraction needs -M, -O or both"));
    assert!(bad(&["-a", "g", "--splitOnly", "--nonSplitOnly", "x"]).contains("nothing to count"));
    for ok in [
        &["-a", "g", "--fraction", "-M", "x"][..],
        &["-a", "g", "--fraction", "-O", "x"],
        &[
            "-a",
            "g",
            "-s",
            "2",
            "--fracOverlap",
            "1",
            "--fracOverlapFeature",
            "0",
            "x",
        ],
    ] {
        assert!(parse(&args(ok)).is_ok(), "{ok:?}");
    }
}

#[test]
fn with_no_arguments_or_help_the_usage_is_written() {
    assert_eq!(run(&[]), 0);
    assert_eq!(run(&args(&["--help"])), 0);
    assert_eq!(run(&args(&["--nope"])), 1);
    assert!(USAGE.contains("helixloom count"));
}

#[test]
fn reads_simulated_from_genes_and_mapped_are_counted_on_the_genes_they_came_from() {
    let d = folder("genes");
    let g = genome(20_000, 201);
    // three genes, apart from one another; reads come from inside one of them
    let genes = [
        (1000usize, 3000usize, '+'),
        (6000, 8000, '-'),
        (12_000, 15_000, '+'),
    ];
    let wanted = [7usize, 3, 5];
    let mut reads = Vec::new();
    for (k, (&(from, to, _), &n)) in genes.iter().zip(&wanted).enumerate() {
        for i in 0..n {
            let start = from + (i * 97) % (to - from - 600);
            reads.push(Read {
                name: format!("g{k}-{i}"),
                seq: g[start..start + 600].to_vec(),
                qual: Vec::new(),
            });
        }
    }
    // one read from between genes
    reads.push(Read {
        name: "between".into(),
        seq: g[4000..4600].to_vec(),
        qual: Vec::new(),
    });
    let o = Options {
        cigar: true,
        ..Options::default()
    };
    let sam = map_to_sam(&[("ref".to_string(), g)], &reads, &o).unwrap();
    std::fs::write(d.join("m.sam"), samfile::write(&sam)).unwrap();
    let gtf: String = genes
        .iter()
        .enumerate()
        .map(|(k, (from, to, strand))| {
            format!(
                "ref\ttest\texon\t{}\t{to}\t.\t{strand}\t.\tgene_id \"gene{k}\"; transcript_id \"t{k}\";\n",
                from + 1
            )
        })
        .collect();
    std::fs::write(d.join("g.gtf"), gtf).unwrap();
    let p = |n: &str| d.join(n).display().to_string();
    assert_eq!(
        run(&args(&[
            "-a",
            &p("g.gtf"),
            "-o",
            &p("counts.tsv"),
            &p("m.sam"),
            &p("m.sam")
        ])),
        0
    );
    let table = std::fs::read_to_string(p("counts.tsv")).unwrap();
    let rows: Vec<Vec<&str>> = table
        .lines()
        .skip(2)
        .map(|l| l.split('\t').collect())
        .collect();
    for (k, n) in wanted.iter().enumerate() {
        assert_eq!(rows[k][0], format!("gene{k}"));
        assert_eq!(
            (rows[k][6], rows[k][7]),
            (n.to_string().as_str(), n.to_string().as_str())
        );
    }
    let summary = std::fs::read_to_string(p("counts.tsv.summary")).unwrap();
    assert!(summary.contains("Assigned\t15\t15\n"), "{summary}");
    assert!(
        summary.contains("Unassigned_NoFeatures\t1\t1\n"),
        "{summary}"
    );
    // without -o the table goes to the standard output
    assert_eq!(run(&args(&["-a", &p("g.gtf"), &p("m.sam")])), 0);
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn a_file_that_cannot_be_read_or_is_wrong_says_which_one() {
    let d = folder("errors");
    let p = |n: &str| d.join(n).display().to_string();
    std::fs::write(
        d.join("g.gtf"),
        "c\t.\texon\t1\t9\t.\t+\t.\tgene_id \"x\";\n",
    )
    .unwrap();
    std::fs::write(d.join("pair.sam"), "r\t1\tc\t1\t60\t4M\t*\t0\t0\tACGT\t*\n").unwrap();
    std::fs::write(d.join("bad.sam"), "r\t0\tc\n").unwrap();
    std::fs::write(
        d.join("other.sam"),
        "r\t0\tzz\t1\t60\t4M\t*\t0\t0\tACGT\t*\n",
    )
    .unwrap();
    let inner = |a: &[&str]| run_inner(&args(a)).unwrap_err();
    assert!(inner(&["-a", &p("none.gtf"), &p("pair.sam")]).contains("Could not read"));
    assert!(inner(&["-a", &p("g.gtf"), &p("none.sam")]).contains("Could not read"));
    let e = inner(&["-a", &p("g.gtf"), &p("pair.sam")]);
    assert!(e.contains("pair.sam") && e.contains("--each-read"), "{e}");
    let e = inner(&["-a", &p("g.gtf"), &p("bad.sam")]);
    assert!(e.contains("bad.sam") && e.contains("columns"), "{e}");
    let e = inner(&["-a", &p("g.gtf"), "-t", "mRNA", &p("other.sam")]);
    assert!(e.contains("g.gtf") && e.contains("mRNA"), "{e}");
    // names that do not match are not an error: the read overlaps nothing
    assert_eq!(
        run(&args(&[
            "-a",
            &p("g.gtf"),
            "-o",
            &p("o.tsv"),
            &p("other.sam")
        ])),
        0
    );
    assert!(std::fs::read_to_string(p("o.tsv.summary"))
        .unwrap()
        .contains("Unassigned_NoFeatures\t1"));
    assert_eq!(
        run(&args(&[
            "-a",
            &p("g.gtf"),
            "-o",
            "/no/such/folder/o",
            &p("other.sam")
        ])),
        1
    );
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn whether_a_read_lies_on_a_sequence_of_the_annotation_is_noticed() {
    let a = annotation::parse(
        "c\t.\texon\t1\t9\t.\t+\t.\tgene_id \"x\";\n",
        &["exon".to_string()],
        "gene_id",
    )
    .unwrap();
    let d = folder("seen");
    let p = |n: &str| d.join(n).display().to_string();
    std::fs::write(
        d.join("yes.sam"),
        "r\t0\tc\t100\t60\t4M\t*\t0\t0\tACGT\t*\n",
    )
    .unwrap();
    std::fs::write(d.join("no.sam"), "r\t0\tzz\t1\t60\t4M\t*\t0\t0\tACGT\t*\n").unwrap();
    let seen = |file: &str| {
        let a_args = parse(&args(&["-a", "g", &p(file)])).unwrap();
        let index = Index::new(&a, false);
        counts(&a_args, &index).unwrap().1
    };
    assert!(seen("yes.sam"));
    assert!(!seen("no.sam"));
    std::fs::remove_dir_all(&d).ok();
}
