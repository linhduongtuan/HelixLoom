use super::*;
use crate::seq::minimap::sim::{genome, read, Rng};

fn args(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| s.to_string()).collect()
}

fn folder(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("helixloom-map-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn options_are_read_and_a_preset_is_what_the_flags_after_it_change() {
    let a = parse(&args(&[
        "-k", "17", "-c", "-p", "0.5", "ref.fa", "reads.fq",
    ]))
    .unwrap();
    assert_eq!(
        (a.options.k, a.options.cigar, a.options.secondary_ratio),
        (17, true, 0.5)
    );
    assert_eq!(a.inputs, ["ref.fa", "reads.fq"]);
    // the preset comes first whatever the place of -x: the flags change what it sets
    let b = parse(&args(&["-w", "7", "-x", "ava-ont", "reads.fq"])).unwrap();
    assert_eq!(
        (b.options.w, b.options.min_score, b.options.overlaps),
        (7, 100.0, true)
    );
    let o = parse(&args(&[
        "-o", "out.paf", "-g", "5000", "-r", "100", "-n", "2", "-m", "20", "-N", "1", "-M", "0.3",
        "-f", "0.01", "t", "q",
    ]))
    .unwrap();
    assert_eq!(o.output.as_deref(), Some("out.paf"));
    assert_eq!(
        (
            o.options.max_gap,
            o.options.bandwidth,
            o.options.min_anchors,
            o.options.min_score
        ),
        (5000, 100, 2, 20.0)
    );
    assert_eq!(
        (
            o.options.max_secondary,
            o.options.mask_level,
            o.options.occ_fraction
        ),
        (1, 0.3, 0.01)
    );
}

#[test]
fn a_flag_that_is_wrong_says_so() {
    let why = |a: &[&str]| parse(&args(a)).err().unwrap();
    assert!(why(&["-k", "x", "t", "q"]).contains("-k needs a number, not x"));
    assert!(why(&["-k", "0", "t", "q"]).contains("k must be"));
    assert!(why(&["-k", "33", "t", "q"]).contains("k must be"));
    assert!(why(&["-w", "0", "t", "q"]).contains("w must be"));
    assert!(why(&["-z", "t", "q"]).contains("-z is not an option"));
    assert!(why(&["-k"]).contains("-k needs a value"));
    assert!(why(&["-x", "splice", "t"]).contains("splice is not a preset"));
    assert!(why(&["-x"]).contains("-x needs a preset"));
}

#[test]
fn a_run_maps_the_reads_in_a_file_and_writes_paf_where_it_is_told() {
    let dir = folder("run");
    let g = genome(15_000, 71);
    let mut r = Rng::new(3);
    let fasta = format!(">ref a reference\n{}\n", String::from_utf8_lossy(&g));
    let mut fastq = String::new();
    for n in 0..4 {
        let q = read(&g, 1000 + 3000 * n, 2000, n % 2 == 1, 0.03, &mut r);
        let bases = String::from_utf8_lossy(&q);
        fastq += &format!("@read{n}\n{bases}\n+\n{}\n", "I".repeat(q.len()));
    }
    let (t, q, out) = (
        dir.join("ref.fa"),
        dir.join("reads.fq"),
        dir.join("out.paf"),
    );
    std::fs::write(&t, fasta).unwrap();
    std::fs::write(&q, fastq).unwrap();
    let a = args(&[
        "-c",
        "-o",
        out.to_str().unwrap(),
        t.to_str().unwrap(),
        q.to_str().unwrap(),
    ]);
    assert_eq!(run(&a), 0);
    let paf = std::fs::read_to_string(&out).unwrap();
    let rows: Vec<Vec<&str>> = paf.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.iter().filter(|f| f.contains(&"tp:A:P")).count(), 4);
    // the name is the first word of the header, and with -c the lines have a CIGAR
    assert!(rows
        .iter()
        .all(|f| f[5] == "ref" && f.iter().any(|x| x.starts_with("cg:Z:"))));
    // overlaps from one file
    let ava = dir.join("ava.paf");
    let a = args(&[
        "-x",
        "ava-ont",
        "-o",
        ava.to_str().unwrap(),
        q.to_str().unwrap(),
    ]);
    assert_eq!(run(&a), 0);
    assert!(std::fs::read_to_string(&ava)
        .unwrap()
        .lines()
        .all(|l| l.split('\t').nth(11) == Some("255")));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn what_cannot_be_run_says_why_and_the_help_is_the_help() {
    let dir = folder("errors");
    let f = dir.join("x.fa");
    std::fs::write(&f, ">a\nACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    let why = |a: &[&str]| run_inner(&args(a)).err().unwrap();
    assert!(why(&["/nonexistent/ref.fa", "/nonexistent/reads.fq"])
        .contains("Could not read /nonexistent/ref.fa"));
    assert!(why(&[f.to_str().unwrap()]).contains("Name the reads after the target"));
    assert!(why(&[]).contains("Name a target"));
    assert!(
        why(&["-x", "ava-ont", f.to_str().unwrap(), f.to_str().unwrap()])
            .contains("one file of reads")
    );
    let junk = dir.join("junk.txt");
    std::fs::write(&junk, "hello").unwrap();
    assert!(why(&[junk.to_str().unwrap(), f.to_str().unwrap()]).contains("neither FASTA nor FASTQ"));
    assert_eq!(run(&args(&["--help"])), 0);
    assert_eq!(run(&[]), 0);
    assert_eq!(run(&args(&["/nonexistent", "/nonexistent"])), 1);
    assert!(USAGE.contains("-x PRESET") && USAGE.contains("docs/MAP.md"));
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn with_a_the_answer_is_sam_with_the_header_the_qualities_and_a_line_for_a_read_that_did_not_map() {
    let dir = folder("sam");
    let g = genome(9_000, 72);
    let mut r = Rng::new(4);
    let good = read(&g, 2000, 2500, true, 0.03, &mut r);
    let junk = genome(1200, 999);
    let text = |name: &str, bases: &[u8], q: char| {
        format!(
            "@{name}\n{}\n+\n{}\n",
            String::from_utf8_lossy(bases),
            q.to_string().repeat(bases.len())
        )
    };
    let (t, q, out) = (
        dir.join("ref.fa"),
        dir.join("reads.fq"),
        dir.join("out.sam"),
    );
    std::fs::write(&t, format!(">ref\n{}\n", String::from_utf8_lossy(&g))).unwrap();
    std::fs::write(&q, text("good", &good, '5') + &text("junk", &junk, 'I')).unwrap();
    let a = args(&[
        "-a",
        "-o",
        out.to_str().unwrap(),
        t.to_str().unwrap(),
        q.to_str().unwrap(),
    ]);
    assert!(parse(&a).unwrap().sam);
    assert!(parse(&a).unwrap().options.cigar);
    assert_eq!(run(&a), 0);
    let sam = std::fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = sam.lines().collect();
    assert_eq!(lines[0], "@HD\tVN:1.6\tSO:unsorted");
    assert_eq!(lines[1], "@SQ\tSN:ref\tLN:9000");
    let rows: Vec<Vec<&str>> = lines
        .iter()
        .filter(|l| !l.starts_with('@'))
        .map(|l| l.split('\t').collect())
        .collect();
    let on = rows
        .iter()
        .find(|f| f[0] == "good" && f[1] == "16")
        .expect("the read maps on the other strand");
    assert_eq!(on[2], "ref");
    assert!(on[5].contains('M'));
    // the qualities are the file's, '5' is 20
    assert_eq!(on[9].len(), on[10].len());
    assert!(on[10].bytes().all(|b| b == b'5'));
    let lost = rows.iter().find(|f| f[0] == "junk").unwrap();
    assert_eq!(&lost[1..6], ["4", "*", "0", "0", "*"]);
    // the others: overlaps are PAF and a single file is not enough
    let ava = args(&["-a", "-x", "ava-ont", q.to_str().unwrap()]);
    assert_eq!(run(&ava), 1);
    assert_eq!(run(&args(&["-a", q.to_str().unwrap()])), 1);
    std::fs::remove_dir_all(&dir).ok();
}
