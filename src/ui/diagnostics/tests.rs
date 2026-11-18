use super::*;
use crate::ui::app::{Action, Startup};
use crate::ui::commands;
use crate::ui::testing::Screen;

#[test]
fn the_command_copies_the_version_the_system_and_what_was_found_and_nothing_of_the_persons_name() {
    let dir = tools::work_folder("diagnostics-ui").expect("a folder");
    // the place of a fault is a path of the person who built the program: under the home folder here
    let home = std::env::var("HOME").ok().filter(|h| h.len() > 1);
    let place = format!("{}/src/x.rs:9", home.as_deref().unwrap_or("/somewhere"));
    let note = diagnostics::crash_note("2026-11-13", "0.2.0", "main", &place);
    std::fs::write(dir.join(diagnostics::CRASH_FILE), &note).expect("written");
    let mut app = App::new(Startup {
        settings_dir: Some(dir.clone()),
        ..Startup::default()
    });
    let mut screen = Screen::new();
    app.platform.no_programs_here();
    screen.settle(&mut app);
    assert_eq!(screen.copied(), None);
    app.run(Action::Diagnostics, screen.ctx());
    screen.settle(&mut app);
    assert_eq!(app.status, COPIED);
    let copied = screen.copied().expect("the diagnostics").to_string();
    assert_eq!(copied, app.diagnostics_text());
    let version = format!("HelixLoom {}\n", env!("CARGO_PKG_VERSION"));
    assert!(copied.starts_with(&version), "{copied}");
    for line in [
        format!(
            "System: {} on {}\n",
            std::env::consts::OS,
            std::env::consts::ARCH
        ),
        "Workstation: not connected\n".to_string(),
        // the folder the settings are in is said in words: its path is not in the text
        "Settings folder: the usual folder\n".to_string(),
        "Appearance: light, 100 %\n".to_string(),
        // on a computer where no program was found, each is said to be missing
        "  MAFFT: not found\n".to_string(),
        "  Nextflow: not found\n".to_string(),
        "Model files:\n  phylocnn: ".to_string(),
        "\n  pfam: ".to_string(),
        format!(
            "Last fault:\n  2026-11-13, HelixLoom 0.2.0, thread main\n  at {}/src/x.rs:9\n",
            if home.is_some() { "~" } else { "/somewhere" }
        ),
    ] {
        assert!(copied.contains(&line), "{line} in {copied}");
    }
    // the home folder is written as a tilde wherever a path has it
    if let Some(home) = &home {
        assert!(!copied.contains(home.as_str()), "{home} in {copied}");
    }
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn the_command_is_in_the_help_group_of_the_palette_and_needs_no_sequence() {
    let command = commands::all()
        .iter()
        .find(|c| matches!(c.action, Action::Diagnostics))
        .expect("a command");
    assert_eq!((command.group, command.title), ("Help", "Copy diagnostics"));
    assert!(!commands::acts_on_the_sequence(&Action::Diagnostics));
}

#[test]
fn a_fault_on_any_thread_leaves_its_note_in_the_folder_it_was_told() {
    let dir = tools::work_folder("fault-note").expect("a folder");
    diagnostics::note_faults_in(dir.clone());
    // other tests may fault on purpose at the same moment and leave their note instead: tried a few times
    let mut noted = None;
    for attempt in 0..5 {
        let marker = format!("a fault made by the test, number {attempt}");
        let said = marker.clone();
        let joined = std::thread::Builder::new()
            .name("faulty".to_string())
            .spawn(move || panic!("{said}"))
            .expect("a thread")
            .join();
        assert!(joined.is_err());
        // the note of this fault is the one of this thread, at this file
        if let Some(note) = diagnostics::last_crash(&dir)
            .filter(|n| n.contains("thread faulty") && n.contains("src/ui/diagnostics/tests.rs:"))
        {
            noted = Some(note);
            break;
        }
    }
    let note = noted.expect("the note of the fault");
    let lines: Vec<&str> = note.lines().collect();
    // the day and the thread, and the place in the code: not what the fault said
    assert_eq!(lines.len(), 2, "{note}");
    assert!(!note.contains("a fault made by the test"), "{note}");
    assert!(
        lines[0].ends_with(&format!(
            ", HelixLoom {}, thread faulty",
            env!("CARGO_PKG_VERSION")
        )),
        "{note}"
    );
    assert!(
        lines[1].starts_with("at src/ui/diagnostics/tests.rs:"),
        "{note}"
    );
    let _ = std::fs::remove_dir_all(dir);
}
