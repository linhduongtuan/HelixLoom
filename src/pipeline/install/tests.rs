use super::*;
use crate::pipeline::tools::{CURL, HMMSCAN, MAFFT};

fn found(managers: &[Manager]) -> Vec<Found> {
    managers
        .iter()
        .map(|m| Found {
            manager: *m,
            path: PathBuf::from(format!("/somewhere/{}", m.command())),
        })
        .collect()
}

#[test]
fn each_manager_is_told_the_name_of_the_programs_own_package_and_nothing_a_person_typed() {
    assert_eq!(
        arguments(&HMMSCAN, Manager::Homebrew),
        Some(vec!["install".to_string(), "hmmer".to_string()])
    );
    for manager in [Manager::Micromamba, Manager::Mamba, Manager::Conda] {
        assert_eq!(
            arguments(&MAFFT, manager).expect("a package").join(" "),
            "install --yes --channel conda-forge --channel bioconda mafft"
        );
    }
    // curl has a formula and no package in Bioconda
    assert!(arguments(&CURL, Manager::Homebrew).is_some());
    assert_eq!(arguments(&CURL, Manager::Conda), None);
}

#[test]
fn a_manager_is_named_by_its_label_or_its_program_in_any_case() {
    assert_eq!(Manager::parse("brew"), Some(Manager::Homebrew));
    assert_eq!(Manager::parse("Homebrew"), Some(Manager::Homebrew));
    assert_eq!(Manager::parse("MICROMAMBA"), Some(Manager::Micromamba));
    assert_eq!(Manager::parse("mamba"), Some(Manager::Mamba));
    assert_eq!(Manager::parse("conda"), Some(Manager::Conda));
    assert_eq!(Manager::parse("apt"), None);
    assert_eq!(Manager::parse(""), None);
}

#[test]
fn the_first_manager_that_can_install_the_program_is_the_plan_and_the_command_is_what_a_person_would_type(
) {
    let all = found(&[Manager::Homebrew, Manager::Micromamba, Manager::Conda]);
    let plan_of = |wanted, environment| plan(&HMMSCAN, &all, wanted, environment);
    let first = plan_of(None, None).expect("a plan");
    assert_eq!(
        (first.program, first.manager, first.line().as_str()),
        ("HMMER", Manager::Homebrew, "brew install hmmer")
    );
    assert_eq!(first.path, PathBuf::from("/somewhere/brew"));
    // the one asked for, when it is there and can
    let asked = plan_of(Some(Manager::Conda), Some("/envs/work")).expect("a plan");
    assert_eq!(
        asked.line(),
        "conda install --yes --channel conda-forge --channel bioconda hmmer"
    );
    assert_eq!(asked.path, PathBuf::from("/somewhere/conda"));
    // without Homebrew the first of conda's kin that is there, when an environment is active
    let kin = found(&[Manager::Micromamba, Manager::Conda]);
    assert_eq!(
        plan(&HMMSCAN, &kin, None, Some("/envs/work"))
            .expect("a plan")
            .manager,
        Manager::Micromamba
    );
}

#[test]
fn what_cannot_be_installed_says_why_and_how_to_do_it_by_hand() {
    // no package manager at all
    let none = plan(&MAFFT, &[], None, None).expect_err("no manager");
    assert!(
        none.starts_with("No package manager was found on this computer"),
        "{none}"
    );
    assert!(
        none.contains(
            "By hand: brew install mafft   or   conda install -c bioconda -c conda-forge mafft."
        ),
        "{none}"
    );
    // the one asked for is not there, and what is there is said
    let missing = plan(
        &MAFFT,
        &found(&[Manager::Homebrew]),
        Some(Manager::Mamba),
        None,
    )
    .expect_err("not found");
    assert_eq!(
        missing,
        "mamba was not found on this computer; it has Homebrew."
    );
    // conda and its kin need an environment that is active, and an empty name is none
    for environment in [None, Some(""), Some("  ")] {
        let no_env = plan(&MAFFT, &found(&[Manager::Micromamba]), None, environment)
            .expect_err("no environment");
        assert!(
            no_env.contains("micromamba installs into the environment that is active, and none is: activate one (`micromamba activate NAME`)"),
            "{no_env}"
        );
    }
    // a program the manager has no package of: said for each manager that was tried
    let no_package = plan(
        &CURL,
        &found(&[Manager::Conda, Manager::Mamba]),
        None,
        Some("/envs/work"),
    )
    .expect_err("no package");
    assert_eq!(
        no_package,
        "HelixLoom cannot install curl here: conda has no package of curl; mamba has no package of curl. By hand: brew install curl."
    );
    // Homebrew needs no environment
    assert!(plan(&MAFFT, &found(&[Manager::Homebrew]), None, None).is_ok());
}

#[test]
fn the_notice_shows_the_command_where_the_program_comes_from_and_its_licence_before_anything_runs()
{
    let brew = plan(&HMMSCAN, &found(&[Manager::Homebrew]), None, None).expect("a plan");
    let notice = brew.notice(&HMMSCAN);
    assert!(notice.starts_with("To install HMMER ("), "{notice}");
    assert!(
        notice.contains("\n\n    brew install hmmer\n\n"),
        "{notice}"
    );
    assert!(
        notice.contains("for this computer, from Homebrew's servers"),
        "{notice}"
    );
    // the very program that will be started is named, not only its short name
    assert!(
        notice.contains("with the brew at /somewhere/brew."),
        "{notice}"
    );
    assert_eq!(brew.environment, None);
    assert!(
        notice.contains("Its licence is BSD-3-Clause as Homebrew lists it"),
        "{notice}"
    );
    let conda = plan(
        &HMMSCAN,
        &found(&[Manager::Conda]),
        None,
        Some("/envs/work"),
    )
    .expect("a plan");
    assert!(
        conda
            .notice(&HMMSCAN)
            .contains("with the conda at /somewhere/conda. This downloads the program into the environment that is active (/envs/work), from the conda-forge and Bioconda channels."),
        "{}",
        conda.notice(&HMMSCAN)
    );
    assert_eq!(conda.environment.as_deref(), Some("/envs/work"));
}

#[test]
fn only_a_plain_yes_is_a_yes() {
    for yes in ["y", "Y", "yes", " YES \n"] {
        assert!(agreed(yes), "{yes:?}");
    }
    for no in ["", "\n", "n", "no", "yep", "ok", "y es"] {
        assert!(!agreed(no), "{no:?}");
    }
}

/// A program that stands in for a package manager: it prints what it was given on both streams, and ends as `code` says.
#[cfg(unix)]
fn stand_in(dir: &std::path::Path, code: i32) -> PathBuf {
    use std::os::unix::fs::PermissionsExt;
    let path = dir.join("brew");
    std::fs::write(
        &path,
        format!("#!/bin/sh\necho \"asked: $*\"\necho \"a line of progress\" 1>&2\necho\necho \"the last line\"\nexit {code}\n"),
    )
    .expect("the stand-in");
    std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).expect("made runnable");
    path
}

#[cfg(unix)]
#[test]
fn the_manager_is_run_with_the_plan_and_every_line_it_prints_is_passed_on_as_it_comes() {
    let dir = tools::work_folder("install-test").expect("a folder");
    let plan = Plan {
        program: "HMMER",
        manager: Manager::Homebrew,
        path: stand_in(&dir, 0),
        args: arguments(&HMMSCAN, Manager::Homebrew).expect("a package"),
        environment: None,
    };
    let mut said = Vec::new();
    assert_eq!(run(&plan, |line| said.push(line.to_string())), Ok(()));
    said.sort();
    assert_eq!(
        said,
        [
            "",
            "a line of progress",
            "asked: install hmmer",
            "the last line"
        ]
    );
    let _ = std::fs::remove_dir_all(dir);
}

#[cfg(unix)]
#[test]
fn a_manager_that_fails_or_cannot_be_started_is_an_error_with_its_last_lines() {
    let dir = tools::work_folder("install-test").expect("a folder");
    let mut plan = Plan {
        program: "HMMER",
        manager: Manager::Homebrew,
        path: stand_in(&dir, 3),
        args: vec!["install".into(), "hmmer".into()],
        environment: None,
    };
    let failed = run(&plan, |_| {}).expect_err("it ends with 3");
    assert!(failed.starts_with("brew stopped ("), "{failed}");
    for line in [
        "asked: install hmmer",
        "a line of progress",
        "the last line",
    ] {
        assert!(failed.contains(line), "{line} in {failed}");
    }
    plan.path = dir.join("no such program");
    let missing = run(&plan, |_| {}).expect_err("nothing to start");
    assert!(missing.starts_with("Could not start "), "{missing}");
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn what_a_terminal_would_obey_is_taken_out_of_a_line_and_the_rest_is_kept() {
    // colours and cursor moves
    assert_eq!(
        clean("\u{1b}[1;32m==>\u{1b}[0m Pouring \u{1b}[2Khmmer"),
        "==> Pouring hmmer"
    );
    // a title, and a write to the clipboard, to their bell or to their end code
    assert_eq!(clean("\u{1b}]0;a title\u{7}after"), "after");
    assert_eq!(clean("before\u{1b}]52;c;QUJD\u{1b}\\after"), "beforeafter");
    // any other code is two characters; a bell, a backspace and a carriage return go, a tab stays
    assert_eq!(clean("a\u{1b}Mb\u{7}c\u{8}d\re\tf"), "abcde\tf");
    // a code that the line ends in the middle of takes nothing more
    assert_eq!(clean("end\u{1b}"), "end");
    assert_eq!(clean("end\u{1b}[31"), "end");
    assert_eq!(clean("plain é 中 text"), "plain é 中 text");
}

/// A stand-in that runs the lines of `script` as a shell would.
#[cfg(unix)]
fn script(dir: &std::path::Path, script: &str) -> Plan {
    use std::os::unix::fs::PermissionsExt;
    let path = dir.join("brew");
    std::fs::write(&path, format!("#!/bin/sh\n{script}\n")).expect("the stand-in");
    std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).expect("made runnable");
    Plan {
        program: "HMMER",
        manager: Manager::Homebrew,
        path,
        args: Vec::new(),
        environment: None,
    }
}

#[cfg(unix)]
#[test]
fn a_line_that_never_ends_comes_in_pieces_and_colour_codes_do_not_reach_the_person() {
    let dir = tools::work_folder("install-test").expect("a folder");
    // 10,000 letters with no end of line, then a coloured line
    let plan = script(
        &dir,
        "head -c 10000 /dev/zero | tr '\\0' 'x'\nprintf '\\n\\033[32mdone\\033[0m\\n'",
    );
    let mut said = Vec::new();
    assert_eq!(run(&plan, |line| said.push(line.to_string())), Ok(()));
    assert!(
        said.iter().all(|line| line.len() <= 4_096),
        "{:?}",
        said.iter().map(String::len).collect::<Vec<_>>()
    );
    let letters: usize = said.iter().map(|line| line.matches('x').count()).sum();
    assert_eq!(letters, 10_000);
    assert_eq!(said.last().map(String::as_str), Some("done"));
    let _ = std::fs::remove_dir_all(dir);
}

#[cfg(unix)]
#[test]
fn a_process_the_manager_leaves_behind_holding_its_output_is_not_waited_for() {
    let dir = tools::work_folder("install-test").expect("a folder");
    // the sleeper keeps the manager's output open for twenty seconds after the manager has ended
    let plan = script(&dir, "sleep 20 &\necho installed\nexit 0");
    let started = std::time::Instant::now();
    let mut said = Vec::new();
    assert_eq!(run(&plan, |line| said.push(line.to_string())), Ok(()));
    assert!(
        started.elapsed() < Duration::from_secs(15),
        "{:?}",
        started.elapsed()
    );
    assert_eq!(said, ["installed"]);
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn strings_meant_for_the_terminal_itself_and_marks_that_turn_the_text_are_taken_out_too() {
    assert_eq!(clean("a\u{1b}Pa string for the terminal\u{1b}\\b"), "ab");
    assert_eq!(clean("a\u{1b}_another\u{7}b"), "ab");
    // a mark that makes what follows read from right to left would show a name as another
    assert_eq!(clean("safe\u{202e}txt.exe"), "safetxt.exe");
    assert_eq!(clean("a\u{2066}b\u{2069}c"), "abc");
    // a control character of the second set
    assert_eq!(clean("a\u{9b}b"), "ab");
}

#[cfg(unix)]
#[test]
fn a_process_left_behind_that_goes_on_printing_is_listened_to_for_a_moment_only() {
    let dir = tools::work_folder("install-test").expect("a folder");
    // twenty seconds of a line every tenth of a second, from a process that outlives the manager
    let plan = script(
        &dir,
        "(i=0; while [ $i -lt 200 ]; do echo tick; sleep 0.1; i=$((i+1)); done) &\necho installed\nexit 0",
    );
    let started = std::time::Instant::now();
    let mut said = Vec::new();
    assert_eq!(run(&plan, |line| said.push(line.to_string())), Ok(()));
    assert!(
        started.elapsed() < Duration::from_secs(10),
        "{:?}",
        started.elapsed()
    );
    assert!(said.iter().any(|line| line == "installed"), "{said:?}");
    let _ = std::fs::remove_dir_all(dir);
}

#[cfg(unix)]
#[test]
fn a_bar_of_progress_is_its_last_state_and_a_character_that_a_long_line_is_cut_in_stays_whole() {
    let dir = tools::work_folder("install-test").expect("a folder");
    // a bar drawn three times over on one line; then a letter and 3,000 letters of two bytes each with no end of line
    // before the last, so that the cut after 4,096 bytes falls in the middle of one
    let plan = script(
        &dir,
        "printf '10%%\\r50%%\\r100%%\\n'\nprintf 'x'\ni=0; while [ $i -lt 3000 ]; do printf 'é'; i=$((i+1)); done\nprintf '\\n'",
    );
    let mut said = Vec::new();
    assert_eq!(run(&plan, |line| said.push(line.to_string())), Ok(()));
    assert_eq!(
        said.first().map(String::as_str),
        Some("100%"),
        "{:?}",
        &said[..1]
    );
    let long: String = said[1..].concat();
    assert_eq!(long.chars().count(), 3_001);
    assert_eq!(long, format!("x{}", "é".repeat(3_000)));
    assert!(
        said.len() >= 3,
        "the long line came in pieces: {}",
        said.len()
    );
    let _ = std::fs::remove_dir_all(dir);
}

#[cfg(unix)]
#[test]
fn a_failure_is_told_by_the_last_three_lines_the_manager_printed() {
    let dir = tools::work_folder("install-test").expect("a folder");
    let plan = script(
        &dir,
        "echo one\necho two\necho\necho three\necho four\nexit 2",
    );
    let failed = run(&plan, |_| {}).expect_err("it ends with 2");
    assert!(failed.ends_with(": two | three | four"), "{failed}");
    assert!(!failed.contains("one"), "{failed}");
    let _ = std::fs::remove_dir_all(dir);
}
