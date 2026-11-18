//! `helixloom motif` as a pipeline runs it: the real program, started with files and options. Sequences are made here with a
//! motif put in each at a known place and strand, so that what it writes can be checked against what was planted.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const WORD: &[u8] = b"ACCGTTAGCA";

struct Folder(PathBuf);

impl Folder {
    fn new(label: &str) -> Folder {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let path = std::env::temp_dir().join(format!(
            "helixloom-motif-{label}-{}-{nanos}",
            std::process::id()
        ));
        std::fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }

    fn file(&self, name: &str, text: &str) -> PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, text).expect("a file");
        path
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn motif(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_helixloom"))
        .arg("motif")
        .args(args)
        .output()
        .expect("the program starts")
}

fn path(p: &Path) -> &str {
    p.to_str().expect("a path that is text")
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn stdout(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).into_owned()
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'A',
    }
}

/// `count` sequences of 150 random bases with WORD in each, on a random strand, as FASTA (`s0`, `s1`, ...), with where and on which
/// strand: a simple generator that is the same on every run.
fn planted(count: usize) -> (String, Vec<(usize, bool)>) {
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = move |n: usize| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 33) % n as u64) as usize
    };
    let (mut fasta, mut places) = (String::new(), Vec::new());
    for i in 0..count {
        let mut seq: Vec<u8> = (0..150).map(|_| b"ACGT"[next(4)]).collect();
        let reverse = next(2) == 1;
        let copy: Vec<u8> = if reverse {
            WORD.iter().rev().map(|&b| complement(b)).collect()
        } else {
            WORD.to_vec()
        };
        let at = next(150 - WORD.len() + 1);
        seq[at..at + WORD.len()].copy_from_slice(&copy);
        fasta.push_str(&format!(
            ">s{i} planted\n{}\n",
            String::from_utf8_lossy(&seq)
        ));
        places.push((at, reverse));
    }
    (fasta, places)
}

/// The most likely base of each column of the first motif of a MEME text.
fn consensus(meme: &str) -> String {
    meme.lines()
        .skip_while(|l| !l.starts_with("letter-probability"))
        .skip(1)
        .take_while(|l| l.starts_with(' '))
        .map(|l| {
            let p: Vec<f64> = l
                .split_whitespace()
                .map(|x| x.parse().expect("a number"))
                .collect();
            let best = (0..4).max_by(|&a, &b| p[a].total_cmp(&p[b])).unwrap_or(0);
            "ACGT".chars().nth(best).unwrap_or('A')
        })
        .collect()
}

fn differences(a: &str, b: &str) -> usize {
    a.chars().zip(b.chars()).filter(|(x, y)| x != y).count()
}

#[test]
fn a_planted_motif_comes_out_as_meme_text_and_a_table_of_the_places() {
    let folder = Folder::new("planted");
    let (fasta, places) = planted(30);
    let input = folder.file("seqs.fa", &fasta);
    let (meme_path, sites_path) = (folder.0.join("m.meme"), folder.0.join("sites.tsv"));
    let out = motif(&[
        "-i",
        path(&input),
        "-o",
        path(&meme_path),
        "--sites",
        path(&sites_path),
        "-w",
        "10",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stdout(&out).is_empty(), "the motifs went to the file");
    let meme = std::fs::read_to_string(&meme_path).expect("the motifs");
    assert!(meme.starts_with("MEME version 4\n"), "{meme}");
    assert_eq!(meme.matches("MOTIF ").count(), 1);
    assert!(meme.contains("w= 10 nsites="), "{meme}");
    // the consensus is the word, or the word of the other strand
    let found = consensus(&meme);
    let other: String = WORD.iter().rev().map(|&b| complement(b) as char).collect();
    let word = String::from_utf8_lossy(WORD).into_owned();
    assert!(
        differences(&found, &word).min(differences(&found, &other)) <= 1,
        "{found} against {word} / {other}"
    );
    // the places: almost every sequence has one, at the planted place on the strand it was planted (or all on the other)
    let table = std::fs::read_to_string(&sites_path).expect("the places");
    let mut lines = table.lines();
    assert_eq!(
        lines.next(),
        Some("motif\tsequence\tposition\tstrand\tscore\tprobability")
    );
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
    assert!(rows.len() >= 27, "{} of 30 sequences", rows.len());
    let (mut same, mut flipped, mut at_place) = (0, 0, 0);
    for row in &rows {
        let k: usize = row[1][1..].parse().expect("s<number>");
        let (at, reverse) = places[k];
        at_place += usize::from(row[2].parse::<usize>().expect("a position") == at + 1);
        match (row[3], reverse) {
            ("+", false) | ("-", true) => same += 1,
            _ => flipped += 1,
        }
    }
    assert!(
        at_place * 10 >= rows.len() * 9,
        "{at_place} of {} at the place",
        rows.len()
    );
    assert!(
        same.max(flipped) * 10 >= rows.len() * 9,
        "strands: {same} as planted, {flipped} the other way"
    );
    // the log says what was found
    let log = stderr(&out);
    assert!(log.contains("30 sequences, 4500 bases"), "{log}");
    assert!(log.contains("motif 1:"), "{log}");
}

#[test]
fn several_motifs_go_to_the_standard_output_as_one_file_and_more_effort_is_taken() {
    let folder = Folder::new("several");
    let (fasta, _) = planted(20);
    let input = folder.file("seqs.fa", &fasta);
    let out = motif(&[
        "-i",
        path(&input),
        "-w",
        "10",
        "-n",
        "2",
        "--one-strand",
        "--effort",
        "2",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    let meme = stdout(&out);
    assert_eq!(meme.matches("MEME version 4").count(), 1, "{meme}");
    assert_eq!(meme.matches("MOTIF ").count(), 2, "{meme}");
    // the last argument that is no option is the input, as in ASTER
    let bare = motif(&["-w", "10", path(&input)]);
    assert!(bare.status.success(), "{}", stderr(&bare));
    assert_eq!(stdout(&bare).matches("MOTIF ").count(), 1);
}

#[test]
fn what_cannot_be_done_is_said_on_the_standard_error_with_a_failing_exit_code() {
    let folder = Folder::new("refused");
    let (fasta, _) = planted(4);
    let input = folder.file("seqs.fa", &fasta);
    let one = folder.file("one.fa", ">only\nACGTACGTACGT\n");
    for (args, why) in [
        (vec![], "No sequences given"),
        (vec!["-i", "no-such-file.fa"], "Could not read"),
        (vec!["-i", path(&one)], "at least two sequences"),
        (vec!["-i", path(&input), "-w", "2"], "between 3 and 40"),
        (vec!["-i", path(&input), "-w", "wide"], "needs a number"),
        (vec!["-i", path(&input), "--effort", "0"], "--effort"),
        (vec!["-i", path(&input), "--nonsense", "1"], "not known"),
    ] {
        let out = motif(&args);
        assert!(!out.status.success(), "{args:?}");
        assert!(stderr(&out).contains(why), "{args:?}: {}", stderr(&out));
        assert!(stdout(&out).is_empty(), "{args:?}: nothing on the output");
    }
    let help = motif(&["--help"]);
    assert!(help.status.success());
    let text = stdout(&help);
    for part in ["--sites", "--effort", "--one-strand", "MEME"] {
        assert!(text.contains(part), "{part}");
    }
}
