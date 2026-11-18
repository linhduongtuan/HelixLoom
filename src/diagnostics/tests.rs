use super::*;

fn facts() -> Facts {
    Facts {
        version: "0.2.0".to_string(),
        os: "macos".to_string(),
        arch: "aarch64".to_string(),
        features: vec!["gpu".to_string()],
        programs: Some(vec![
            Looked {
                label: "MAFFT".to_string(),
                found: Some((
                    "/Users/ada/tools/bin/mafft".to_string(),
                    "v7.526".to_string(),
                )),
            },
            Looked {
                label: "curl".to_string(),
                found: Some(("/usr/bin/curl".to_string(), String::new())),
            },
            Looked {
                label: "Snakemake".to_string(),
                found: None,
            },
        ]),
        managers: vec!["Homebrew".to_string(), "conda".to_string()],
        models: vec![("phylocnn".to_string(), false), ("pfam".to_string(), true)],
        settings: "the usual folder".to_string(),
        workstation: "connected (HelixLoom 0.2.0, 2 GPUs)".to_string(),
        appearance: "dark, 110 %".to_string(),
        last_crash: Some(
            "2026-11-13, HelixLoom 0.2.0, thread main\nat /Users/ada/src/helixloom/src/x.rs:9\n"
                .to_string(),
        ),
        home: Some("/Users/ada".to_string()),
    }
}

#[test]
fn the_text_says_the_version_the_computer_and_what_was_found_a_line_for_each() {
    assert_eq!(
        text(&facts()),
        "HelixLoom 0.2.0\n\
         System: macos on aarch64\n\
         Built with: gpu\n\
         Appearance: dark, 110 %\n\
         Settings folder: the usual folder\n\
         Workstation: connected (HelixLoom 0.2.0, 2 GPUs)\n\
         Programs:\n\
         \x20 MAFFT: ~/tools/bin/mafft (v7.526)\n\
         \x20 curl: /usr/bin/curl\n\
         \x20 Snakemake: not found\n\
         Package managers: Homebrew, conda\n\
         Model files:\n\
         \x20 phylocnn: not installed\n\
         \x20 pfam: installed\n\
         Last fault:\n\
         \x20 2026-11-13, HelixLoom 0.2.0, thread main\n\
         \x20 at ~/src/helixloom/src/x.rs:9\n"
    );
}

#[test]
fn what_is_not_there_is_said_plainly() {
    let bare = Facts {
        features: Vec::new(),
        programs: None,
        managers: Vec::new(),
        settings: "none".to_string(),
        last_crash: None,
        ..facts()
    };
    let said = text(&bare);
    for line in [
        "Built with: no optional features\n",
        "Settings folder: none\n",
        "Programs:\n  still being looked for\n",
        "Package managers: none found\n",
        "Last fault: none noted\n",
    ] {
        assert!(said.contains(line), "{line} in {said}");
    }
}

#[test]
fn the_home_folder_is_written_as_a_tilde_and_nothing_of_the_persons_name_is_left_in_the_text() {
    let home = Some("/Users/ada");
    assert_eq!(tilde("/Users/ada/x/y", home), "~/x/y");
    assert_eq!(tilde("/Users/ada", home), "~");
    assert_eq!(tilde("/Users/ada/", Some("/Users/ada/")), "~/");
    // another folder whose name begins the same is not the home folder
    assert_eq!(tilde("/Users/adam/x", home), "/Users/adam/x");
    assert_eq!(tilde("/usr/bin/curl", home), "/usr/bin/curl");
    assert_eq!(tilde("/usr/bin/curl", None), "/usr/bin/curl");
    assert_eq!(tilde("/usr/bin/curl", Some("")), "/usr/bin/curl");
    // in the middle of a line too, as the place of a fault is written
    assert_eq!(tilde("at /Users/ada/src/x.rs:9", home), "at ~/src/x.rs:9");
    assert_eq!(tilde("in /Users/ada", home), "in ~");
    assert_eq!(tilde("/Users/ada/a and /Users/ada/b", home), "~/a and ~/b");
    // whatever follows it that is no part of a name: a colon, a bracket, a quote, a comma
    assert_eq!(tilde("/Users/ada: no such file", home), "~: no such file");
    assert_eq!(tilde("(/Users/ada), '/Users/ada'", home), "(~), '~'");
    // a name that goes on is another folder
    assert_eq!(
        tilde("/Users/ada.old/x and /Users/ada-2/y", home),
        "/Users/ada.old/x and /Users/ada-2/y"
    );
    // the other kind of slash, as another system writes its paths
    assert_eq!(
        tilde("C:\\Users\\ada\\tools\\mafft.exe", Some("C:\\Users\\ada")),
        "~\\tools\\mafft.exe"
    );
    // a home folder that is the top of everything hides nothing
    assert_eq!(tilde("/usr/bin/curl", Some("/")), "/usr/bin/curl");
    assert!(!text(&facts()).contains("ada"), "{}", text(&facts()));
}

#[test]
fn the_note_of_a_fault_says_when_and_where_and_never_what_the_program_said() {
    assert_eq!(
        crash_note("2026-11-13", "0.2.0", "main", "src/x.rs:9"),
        "2026-11-13, HelixLoom 0.2.0, thread main\nat src/x.rs:9\n"
    );
    // each thing is one line of a bounded length, whatever it holds
    let note = crash_note("d", "v", "a\nthread", &"p".repeat(5_000));
    assert_eq!(note.lines().count(), 2, "{note}");
    assert!(
        note.starts_with("d, HelixLoom v, thread a thread\nat ppp"),
        "{note}"
    );
    assert!(note.len() < 300, "{}", note.len());
}

#[test]
fn the_place_of_a_fault_is_told_without_the_folders_of_the_computer_that_built_the_program() {
    // the program's own files are named from `src/` as they are
    assert_eq!(place_of("src/seq/fold.rs", 512), "src/seq/fold.rs:512");
    // a file of another package by the package, not by the home folder of whoever built the program
    assert_eq!(
        place_of("/Users/builder/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/egui-0.36.2/src/ui.rs", 7),
        "egui-0.36.2/src/ui.rs:7"
    );
    assert_eq!(
        place_of(
            "C:\\Users\\builder\\.cargo\\registry\\src\\index.crates.io-1\\png-0.18.0\\src\\lib.rs",
            3
        ),
        "png-0.18.0/src/lib.rs:3"
    );
    // any other whole path from its last `src/`, and one with none says no more than that
    assert_eq!(
        place_of("/Users/builder/work/helixloom/src/main.rs", 9),
        "src/main.rs:9"
    );
    assert_eq!(
        place_of("/rustc/abc123/library/core/option.rs", 2),
        "a file of the program:2"
    );
}

#[test]
fn a_version_is_kept_only_when_it_looks_like_one() {
    for version in [
        "v7.526",
        "3.1.4",
        "0.23.4",
        "8.7.1",
        "24.10.0-edge",
        "1.2_3+b",
    ] {
        assert_eq!(plain_version(version), version);
    }
    // a path, a sentence, something very long: no version rather than that in a report
    for other in [
        "",
        "/opt/x/envs/my.project/bin/python",
        "version 3.4",
        "a.b:c",
        &"1.".repeat(20),
    ] {
        assert_eq!(plain_version(other), "", "{other}");
    }
}

#[test]
fn the_last_fault_is_read_from_the_settings_folder_and_a_folder_without_one_has_none() {
    let dir = crate::pipeline::tools::work_folder("diagnostics").expect("a folder");
    assert_eq!(last_crash(&dir), None);
    std::fs::write(dir.join(CRASH_FILE), "   \n").expect("written");
    assert_eq!(last_crash(&dir), None, "an empty note is none");
    let note = crash_note("2026-11-13", "0.2.0", "main", "src/x.rs:9");
    std::fs::write(dir.join(CRASH_FILE), &note).expect("written");
    assert_eq!(last_crash(&dir), Some(note));
    // a file that grew large is read by its start only
    std::fs::write(dir.join(CRASH_FILE), "y".repeat(100_000)).expect("written");
    assert_eq!(last_crash(&dir).map(|n| n.len()), Some(4_000));
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn the_note_a_fault_leaves_names_its_place_without_the_builders_folders_and_its_day() {
    let note = fault_note(
        Some((
            "/Users/builder/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/egui-0.36.2/src/ui.rs",
            7,
        )),
        "main",
        // the fourteenth of November 2026, at noon
        1_794_657_600,
    );
    assert_eq!(
        note,
        format!(
            "2026-11-14, HelixLoom {}, thread main\nat egui-0.36.2/src/ui.rs:7\n",
            env!("CARGO_PKG_VERSION")
        )
    );
    assert!(fault_note(None, "t", 0).ends_with("\nat an unknown place\n"));
}

#[test]
fn the_settings_folder_is_said_in_words_and_a_program_keeps_its_version_only_when_it_is_one() {
    assert_eq!(settings_in_words(false, false), "none");
    assert_eq!(settings_in_words(false, true), "none");
    assert_eq!(settings_in_words(true, false), "the usual folder");
    assert_eq!(
        settings_in_words(true, true),
        "a folder named by HELIXLOOM_CONFIG_DIR"
    );

    let path = Path::new("/opt/homebrew/bin/mafft");
    let found = looked("MAFFT", Some((path, "v7.526")));
    assert_eq!(found.label, "MAFFT");
    assert_eq!(
        found.found,
        Some(("/opt/homebrew/bin/mafft".to_string(), "v7.526".to_string()))
    );
    // a first line that is a path with a number in it is no version
    let odd = looked("tool", Some((path, "/opt/py3.9/bin/python")));
    assert_eq!(odd.found.map(|(_, version)| version), Some(String::new()));
    assert_eq!(looked("tool", None).found, None);
}
