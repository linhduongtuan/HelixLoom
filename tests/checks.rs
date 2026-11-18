//! `helixloom checks` as a person runs it: the real program, started with words or none.

use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_helixloom"))
        .arg("checks")
        .args(args)
        .output()
        .expect("the program starts")
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn stdout(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).into_owned()
}

#[test]
fn the_checks_command_prints_every_analysis_with_its_grade_and_what_was_not_checked() {
    let out = run(&[]);
    assert!(out.status.success(), "{}", stderr(&out));
    let text = stdout(&out);
    assert!(
        text.starts_with("How each analysis was checked: "),
        "{text}"
    );
    assert!(
        text.contains("\nSpecies tree from gene trees [reference program]\n"),
        "{text}"
    );
    assert!(
        text.contains("\nStructure prediction (a folding program) [not checked]\n"),
        "{text}"
    );
    // every analysis has its four lines
    let analyses = text
        .lines()
        .filter(|l| l.ends_with(']') && !l.starts_with(' '))
        .count();
    assert!(analyses >= 30, "{analyses} analyses in {text}");
    for line in [
        "  compared with: ",
        "  result: ",
        "  not checked: ",
        "  tests: ",
    ] {
        assert_eq!(
            text.lines().filter(|l| l.starts_with(line)).count(),
            analyses,
            "{line}"
        );
    }
}

#[test]
fn words_leave_the_checks_that_hold_them_all_and_words_that_find_none_fail() {
    let out = run(&["ASTRAL", "5.7.8"]);
    assert!(out.status.success(), "{}", stderr(&out));
    let text = stdout(&out);
    assert!(
        text.starts_with("How each analysis was checked: 1 analysis (1 reference program,"),
        "{text}"
    );
    assert!(!text.contains("Primer design"), "{text}");

    let out = run(&["nosuchword"]);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(stdout(&out), "");
    assert_eq!(stderr(&out).trim(), "No analysis matches nosuchword.");
}
