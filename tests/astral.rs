//! `helixloom astral` as a pipeline runs it: the real program, started with files and options, with none of the program's own
//! code in between. The gene trees are the example of the species-tree window (40 genes on 8 species simulated from the tree
//! in `assets/species/truth.tre`); the score 1733 is what this engine and ASTRAL 5.7.8 both give for them.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const GENES: &str = include_str!("../assets/species/genes.tre");
const TRUTH: &str = include_str!("../assets/species/truth.tre");
const SCORE: &str = "Final quartet score is: 1733";

/// A folder for one test, removed when the test is over.
struct Folder(PathBuf);

impl Folder {
    fn new(label: &str) -> Folder {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let path = std::env::temp_dir().join(format!(
            "helixloom-astral-{label}-{}-{nanos}",
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

fn astral(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_helixloom"))
        .arg("astral")
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

#[test]
fn it_writes_the_species_tree_to_a_file_and_the_log_to_the_standard_error() {
    let dir = Folder::new("search");
    let genes = dir.file("genes.tre", GENES);
    let found = dir.0.join("species.tre");
    let out = astral(&["-i", path(&genes), "-o", path(&found)]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stdout(&out).is_empty(), "the tree goes to the file");
    let log = stderr(&out);
    assert!(
        log.contains("40 gene trees on 8 taxa") && log.contains(SCORE),
        "{log}"
    );
    let tree = std::fs::read_to_string(&found).expect("the tree");
    assert!(tree.ends_with(";\n"), "{tree}");
    for species in "ABCDEFGH".chars() {
        assert!(tree.contains(species), "{species} in {tree}");
    }
    // the local posterior probability and a length on every inner branch, as ASTRAL writes them
    assert!(tree.contains(")0.6:") || tree.contains(")1:"), "{tree}");
}

#[test]
fn the_same_input_gives_the_same_bytes_every_time_whatever_the_option_spelling() {
    let dir = Folder::new("same");
    let genes = dir.file("genes.tre", GENES);
    let first = astral(&["-i", path(&genes)]);
    assert!(first.status.success(), "{}", stderr(&first));
    assert!(stdout(&first).trim_end().ends_with(';'));
    for _ in 0..4 {
        // a last argument that is no option is the input, as in ASTER
        let again = astral(&[path(&genes)]);
        assert_eq!(again.stdout, first.stdout);
    }
    let long = astral(&["--input", path(&genes), "--branch-annotate", "3"]);
    assert_eq!(long.stdout, first.stdout);
}

#[test]
fn a_given_tree_is_scored_and_annotated_with_what_each_option_asks_for() {
    let dir = Folder::new("score");
    let genes = dir.file("genes.tre", GENES);
    let truth = dir.file("truth.tre", TRUTH);
    let out = astral(&["-i", path(&genes), "-q", path(&truth), "-t", "2"]);
    assert!(out.status.success(), "{}", stderr(&out));
    let log = stderr(&out);
    assert!(
        log.contains(SCORE),
        "the tree the genes came from has the best score: {log}"
    );
    assert!(!log.contains("Searched"), "nothing was searched: {log}");
    let tree = stdout(&out);
    for part in ["pp1=", "q1=", "f1=", "QC=", "EN="] {
        assert!(tree.contains(part), "{part} in {tree}");
    }
    let quartets = astral(&["-i", path(&genes), "-q", path(&truth), "-t", "8"]);
    assert!(stdout(&quartets).contains("'[q1="), "{}", stdout(&quartets));
    let none = astral(&["-i", path(&genes), "-q", path(&truth), "-t", "0"]);
    assert!(!stdout(&none).contains("'["), "{}", stdout(&none));
}

#[test]
fn the_exact_search_and_the_search_of_clusters_find_the_same_best_score_on_eight_species() {
    let dir = Folder::new("exact");
    let genes = dir.file("genes.tre", GENES);
    let exact = astral(&["-i", path(&genes), "-x"]);
    assert!(exact.status.success(), "{}", stderr(&exact));
    assert!(stderr(&exact).contains(SCORE), "{}", stderr(&exact));
    let plain = astral(&["-i", path(&genes), "-p", "0", "--no-polish"]);
    assert!(stderr(&plain).contains(SCORE), "{}", stderr(&plain));
}

#[test]
fn what_cannot_be_done_is_said_on_the_standard_error_with_a_failing_exit_code() {
    let dir = Folder::new("refuse");
    let genes = dir.file("genes.tre", GENES);
    let family = dir.file("family.tre", "((A,B),(A,C),D);\n");
    let few = dir.file("few.tre", "((A,B),C);\n((A,C),B);\n");
    let cases: [(Vec<&str>, &str); 6] = [
        (vec![], "No gene trees given"),
        (
            vec!["-i", "/nonexistent/genes.tre"],
            "/nonexistent/genes.tre",
        ),
        (vec!["-i", path(&family)], "ASTRAL-Pro"),
        (vec!["-i", path(&few)], "at least four taxa"),
        (vec!["-i", path(&genes), "--nonsense", "1"], "--nonsense"),
        (vec!["-i", path(&genes), "-t", "5"], "-t"),
    ];
    for (args, said) in cases {
        let out = astral(&args);
        assert_eq!(out.status.code(), Some(1), "{args:?}");
        assert!(stderr(&out).contains(said), "{args:?}: {}", stderr(&out));
        assert!(stdout(&out).is_empty(), "{args:?}");
    }
}

#[test]
fn help_lists_the_options_and_succeeds_and_the_other_name_of_the_command_works() {
    let help = astral(&["--help"]);
    assert!(help.status.success());
    let text = stdout(&help);
    for part in [
        "Usage: helixloom astral",
        "--score-tree",
        "--branch-annotate",
        "--exact",
        "--lambda",
    ] {
        assert!(text.contains(part), "{part}");
    }
    let dir = Folder::new("alias");
    let genes = dir.file("genes.tre", GENES);
    let out = Command::new(env!("CARGO_BIN_EXE_helixloom"))
        .args(["species-tree", "-i", path(&genes)])
        .output()
        .expect("the program starts");
    assert!(out.status.success() && String::from_utf8_lossy(&out.stderr).contains(SCORE));
}

#[test]
fn names_with_spaces_and_brackets_are_written_in_quotes_and_the_tree_reads_back() {
    let dir = Folder::new("quotes");
    let genes = dir.file(
        "genes.tre",
        "(('Homo sapiens','Pan troglodytes'),('Nomascus (gibbon)','Pongo pygmaeus'));\n\
         (('Homo sapiens','Nomascus (gibbon)'),('Pan troglodytes','Pongo pygmaeus'));\n\
         (('Homo sapiens','Pan troglodytes'),('Nomascus (gibbon)','Pongo pygmaeus'));\n",
    );
    let out = astral(&["-i", path(&genes)]);
    assert!(out.status.success(), "{}", stderr(&out));
    let tree = stdout(&out);
    for name in [
        "'Homo sapiens'",
        "'Nomascus (gibbon)'",
        "'Pan troglodytes'",
        "'Pongo pygmaeus'",
    ] {
        assert!(tree.contains(name), "{name} in {tree}");
    }
    // the text that was written is a tree that can be given back: scored against the same genes it is the same tree
    let species = dir.file("species.tre", &tree);
    let again = astral(&["-i", path(&genes), "-q", path(&species)]);
    assert!(again.status.success(), "{}", stderr(&again));
    let score = |out: &Output| {
        stderr(out)
            .lines()
            .find(|l| l.starts_with("Final quartet score"))
            .map(str::to_string)
    };
    assert!(score(&out).is_some() && score(&out) == score(&again));
}

#[test]
fn a_given_tree_is_hung_from_the_species_asked_for_and_one_that_is_not_there_is_refused() {
    let dir = Folder::new("root");
    let genes = dir.file("genes.tre", GENES);
    let truth = dir.file("truth.tre", TRUTH);
    let out = astral(&["-i", path(&genes), "-q", path(&truth), "--root", "H"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stdout(&out).starts_with("(H,"), "{}", stdout(&out));
    let nobody = astral(&["-i", path(&genes), "-q", path(&truth), "--root", "NOPE"]);
    assert_eq!(nobody.status.code(), Some(1));
    assert!(
        stderr(&nobody).contains("root taxon NOPE"),
        "{}",
        stderr(&nobody)
    );
}
