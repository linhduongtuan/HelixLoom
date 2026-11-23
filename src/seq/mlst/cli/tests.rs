use super::*;
use crate::seq::mlst::tests::{
    allele, allele_of, contigs, fasta, put, scheme, scheme_from, write_scheme,
};

fn args(a: &[&str]) -> Vec<String> {
    a.iter().map(|s| s.to_string()).collect()
}

fn folder(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("helixloom-mlst-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn the_options_are_read_and_the_defaults_are_those_of_the_module() {
    let a = parse(&args(&["--scheme", "s", "g.fa"])).unwrap();
    assert_eq!(
        (a.params, a.full, a.csv, a.nopath),
        (Params::default(), false, false, false)
    );
    assert_eq!(a.schemes, [PathBuf::from("s")]);
    let b = parse(&args(&[
        "--schemes",
        "d",
        "--minid",
        "90",
        "--mincov",
        "20",
        "--minscore",
        "80",
        "--full",
        "--csv",
        "--nopath",
        "--json",
        "j.json",
        "--novel",
        "n.fa",
        "--label",
        "x",
        "-o",
        "o.tsv",
        "g.fa",
    ]))
    .unwrap();
    assert_eq!(
        (
            b.params.min_identity,
            b.params.min_coverage,
            b.params.min_score
        ),
        (90.0, 20.0, 80.0)
    );
    assert!(b.full && b.csv && b.nopath);
    assert_eq!(
        (
            b.json.as_deref(),
            b.novel.as_deref(),
            b.label.as_deref(),
            b.output.as_deref()
        ),
        (Some("j.json"), Some("n.fa"), Some("x"), Some("o.tsv"))
    );
    assert_eq!(b.scan, Some(PathBuf::from("d")));
}

#[test]
fn a_wrong_option_says_so() {
    let why = |a: &[&str]| parse(&args(a)).err().unwrap();
    assert!(why(&["g.fa"]).contains("Name a scheme"));
    assert!(why(&["--scheme", "s"]).contains("FASTA file of a genome"));
    assert!(
        why(&["--scheme", "s", "--minid", "x", "g.fa"]).contains("--minid needs a number, not x")
    );
    assert!(why(&["--scheme", "s", "--minid", "120", "g.fa"]).contains("from 0 to 100"));
    assert!(why(&["--scheme", "s", "--minscore", "-5", "g.fa"]).contains("from 0 to 100"));
    assert!(why(&["--scheme", "s", "--label", "x", "a.fa", "b.fa"]).contains("one genome"));
    assert!(why(&["--scheme", "s", "--blast", "g.fa"]).contains("--blast is not an option"));
    assert!(why(&["--scheme"]).contains("--scheme needs a value"));
}

#[test]
fn the_limits_of_zero_and_one_hundred_are_allowed() {
    let a = parse(&args(&[
        "--scheme",
        "s",
        "--minid",
        "100",
        "--mincov",
        "0",
        "--minscore",
        "100",
        "g.fa",
    ]))
    .unwrap();
    assert_eq!(
        (
            a.params.min_identity,
            a.params.min_coverage,
            a.params.min_score
        ),
        (100.0, 0.0, 100.0)
    );
}

#[test]
fn a_scheme_whose_score_is_the_limit_is_chosen_and_one_just_under_is_not() {
    let s = scheme("synth", 1, 6);
    // one locus with a novel allele: 63 points
    let mut novel = allele(0, 3);
    novel[5] = if novel[5] == b'A' { b'C' } else { b'A' };
    let c = Contigs::new(contigs(31, &[put(novel, false, 0)])).unwrap();
    let at = |min_score: f64| {
        best(
            std::slice::from_ref(&s),
            &c,
            &Params {
                min_score,
                ..Params::default()
            },
        )
    };
    assert_eq!(at(63.0).unwrap().score, 63.0);
    assert!(at(63.5).is_none());
}

#[test]
fn a_folder_of_schemes_is_read_in_alphabetical_order_and_only_its_folders_count() {
    let dir = folder("scan");
    for n in ["m", "c", "x", "a", "q", "f", "k", "b", "z", "e"] {
        std::fs::create_dir_all(dir.join(n)).unwrap();
    }
    std::fs::write(dir.join("stray.txt"), "").unwrap();
    let names: Vec<String> = scan(&dir)
        .unwrap()
        .iter()
        .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
        .collect();
    assert_eq!(names, ["a", "b", "c", "e", "f", "k", "m", "q", "x", "z"]);
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn a_genome_is_typed_against_a_scheme_folder_and_written_in_each_form() {
    let dir = folder("run");
    let s = scheme("synth", 7, 6);
    write_scheme(&dir, &s);
    // type 3, and one allele of it replaced by a novel one: ~ in the calls
    let ids: Vec<usize> = (0..7).map(|l| (3 + l) % 6 + 1).collect();
    let mut puts: Vec<_> = ids
        .iter()
        .enumerate()
        .map(|(l, &id)| put(allele(l, id), l % 2 == 1, l % 3))
        .collect();
    let genome_path = dir.join("g.fa");
    std::fs::write(&genome_path, fasta(&contigs(1, &puts))).unwrap();
    let out = dir.join("out.tsv");
    let a = args(&[
        "--scheme",
        dir.join("synth").to_str().unwrap(),
        "--full",
        "--nopath",
        "-o",
        out.to_str().unwrap(),
        genome_path.to_str().unwrap(),
    ]);
    assert_eq!(run(&a), 0);
    let text = std::fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "FILE\tSCHEME\tST\tSTATUS\tSCORE\tALLELES");
    let want: Vec<String> = ids
        .iter()
        .enumerate()
        .map(|(l, id)| format!("L{l}({id})"))
        .collect();
    assert_eq!(
        lines[1],
        format!("g.fa\tsynth\t3\tPERFECT\t100\t{}", want.join(";"))
    );
    // a novel allele, with the other forms
    let mut changed = allele(2, ids[2]);
    for at in [100, 200, 300] {
        changed[at] = if changed[at] == b'G' { b'T' } else { b'G' };
    }
    puts[2] = put(changed.clone(), false, 2);
    let novel_genome = dir.join("n.fa");
    std::fs::write(&novel_genome, fasta(&contigs(2, &puts))).unwrap();
    let (json, novel, csv) = (dir.join("o.json"), dir.join("novel.fa"), dir.join("o.csv"));
    let a = args(&[
        "--scheme",
        dir.join("synth").to_str().unwrap(),
        "--csv",
        "--label",
        "sample1",
        "--json",
        json.to_str().unwrap(),
        "--novel",
        novel.to_str().unwrap(),
        "-o",
        csv.to_str().unwrap(),
        novel_genome.to_str().unwrap(),
    ]);
    assert_eq!(run(&a), 0);
    let line = std::fs::read_to_string(&csv).unwrap();
    assert!(line.starts_with("sample1,synth,-,L0("), "{line}");
    assert!(line.contains(&format!("L2(~{})", ids[2])), "{line}");
    assert!(std::fs::read_to_string(&json)
        .unwrap()
        .contains("\"id\": \"sample1\""));
    let fa = std::fs::read_to_string(&novel).unwrap();
    assert!(
        fa.starts_with(">synth.L2-")
            && fa.contains(" sample1\n")
            && fa.contains(&String::from_utf8(changed).unwrap()),
        "{fa}"
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn of_several_schemes_the_one_that_scores_best_is_used_and_none_under_the_limit_gives_a_line_of_dashes(
) {
    let dir = folder("auto");
    let (first, second) = (
        scheme_from("aaa", 2000, 7, 6),
        scheme_from("bbb", 3000, 7, 6),
    );
    write_scheme(&dir, &first);
    write_scheme(&dir, &second);
    // a genome of the second scheme's type 4
    let ids: Vec<usize> = (0..7).map(|l| (4 + l) % 6 + 1).collect();
    let puts: Vec<_> = ids
        .iter()
        .enumerate()
        .map(|(l, &id)| put(allele_of(3000, l, id), false, l % 3))
        .collect();
    let g = dir.join("g.fa");
    std::fs::write(&g, fasta(&contigs(3, &puts))).unwrap();
    let schemes = vec![
        Scheme::load(&dir.join("aaa")).unwrap(),
        Scheme::load(&dir.join("bbb")).unwrap(),
    ];
    let c = Contigs::new(crate::pipeline::steps::read_fasta(
        &std::fs::read_to_string(&g).unwrap(),
    ))
    .unwrap();
    let t = best(&schemes, &c, &Params::default()).expect("a scheme");
    assert_eq!(
        (t.scheme.as_str(), t.sequence_type.as_deref()),
        ("bbb", Some("4"))
    );
    // the same through the command line, with a folder of schemes
    let out = dir.join("out.tsv");
    assert_eq!(
        run(&args(&[
            "--schemes",
            dir.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            g.to_str().unwrap()
        ])),
        0
    );
    assert!(std::fs::read_to_string(&out)
        .unwrap()
        .contains("\tbbb\t4\t"));
    // a genome of neither: no scheme
    let other = dir.join("other.fa");
    std::fs::write(&other, fasta(&contigs(4, &[]))).unwrap();
    assert_eq!(
        run(&args(&[
            "--scheme",
            dir.join("aaa").to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            other.to_str().unwrap()
        ])),
        0
    );
    assert_eq!(
        std::fs::read_to_string(&out).unwrap().trim_end(),
        format!("{}\t-\t-", other.display())
    );
    // two that score alike: the first wins (a genome of both: the same alleles are not shared, so make the scheme twice)
    let twin = scheme_from("ccc", 2000, 7, 6);
    write_scheme(&dir, &twin);
    let tied = vec![
        Scheme::load(&dir.join("ccc")).unwrap(),
        Scheme::load(&dir.join("aaa")).unwrap(),
    ];
    let ids: Vec<usize> = (0..7).map(|l| (2 + l) % 6 + 1).collect();
    let puts: Vec<_> = ids
        .iter()
        .enumerate()
        .map(|(l, &id)| put(allele_of(2000, l, id), false, l % 3))
        .collect();
    let c = Contigs::new(contigs(5, &puts)).unwrap();
    assert_eq!(best(&tied, &c, &Params::default()).unwrap().scheme, "ccc");
    assert_eq!(
        best(
            &tied,
            &c,
            &Params {
                min_score: 100.0,
                ..Params::default()
            }
        )
        .unwrap()
        .scheme,
        "ccc"
    );
    assert!(best(&tied[..0], &c, &Params::default()).is_none());
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn what_cannot_be_run_says_why_and_the_help_is_the_help() {
    let dir = folder("errors");
    let why = |a: &[&str]| run_inner(&args(a)).err().unwrap();
    assert!(
        why(&["--scheme", "/nonexistent/scheme", "g.fa"]).contains("no profile table")
            || why(&["--scheme", "/nonexistent/scheme", "g.fa"]).contains("Could not read")
    );
    write_scheme(&dir, &scheme("synth", 2, 3));
    let s = dir.join("synth");
    assert!(why(&["--scheme", s.to_str().unwrap(), "/nonexistent/g.fa"])
        .contains("Could not read /nonexistent/g.fa"));
    let empty = dir.join("empty.fa");
    std::fs::write(&empty, "no fasta here").unwrap();
    assert!(
        why(&["--scheme", s.to_str().unwrap(), empty.to_str().unwrap()])
            .contains("holds no sequences")
    );
    let nothing = dir.join("nothing");
    std::fs::create_dir_all(&nothing).unwrap();
    assert!(
        why(&["--schemes", nothing.to_str().unwrap(), "g.fa"]).contains("holds no scheme folder")
    );
    assert_eq!(run(&args(&["--help"])), 0);
    assert_eq!(run(&[]), 0);
    assert_eq!(run(&args(&["--scheme", "/nonexistent", "g.fa"])), 1);
    assert!(USAGE.contains("--schemes FOLDER") && USAGE.contains("docs/MLST.md"));
    std::fs::remove_dir_all(&dir).ok();
}
