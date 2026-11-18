use super::*;
use crate::pipeline::install::Manager;
use crate::ui::app::{Action, Startup};
use crate::ui::commands;
use crate::ui::testing::Screen;
use std::path::PathBuf;
use std::time::{Duration, Instant};

/// A computer with Homebrew alone, which is not really there.
fn only_homebrew() -> Vec<Found> {
    vec![Found {
        manager: Manager::Homebrew,
        path: PathBuf::from("/nowhere/brew"),
    }]
}

fn no_manager() -> Vec<Found> {
    Vec::new()
}

/// An installer that must not be started.
fn never(_: &Plan, _: &mut dyn FnMut(&str)) -> Result<(), String> {
    panic!("nothing may run before the person says so")
}

/// What the stand-in for a package manager was asked to run.
static ASKED: Mutex<Vec<String>> = Mutex::new(Vec::new());

fn pretend(plan: &Plan, say: &mut dyn FnMut(&str)) -> Result<(), String> {
    ASKED.lock().expect("the list").push(plan.line());
    say("==> Fetching the program");
    say("==> Pouring it");
    Ok(())
}

fn failing(_: &Plan, say: &mut dyn FnMut(&str)) -> Result<(), String> {
    say("Error: no formula of that name");
    Err("brew stopped (exit status: 1): Error: no formula of that name".to_string())
}

/// The window open on a computer where no program was found, with stand-ins for the managers and the installer.
fn opened(managers: fn() -> Vec<Found>, installer: Installer) -> (App, Screen) {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    app.platform.no_programs_here();
    app.programs_win.managers = managers;
    app.programs_win.installer = installer;
    app.run(Action::Programs, screen.ctx());
    screen.settle(&mut app);
    (app, screen)
}

#[test]
fn the_window_lists_every_program_with_what_it_is_for_its_licence_and_that_it_is_not_found() {
    let (app, screen) = opened(only_homebrew, never);
    let text = screen.text();
    assert!(text.contains(TITLE), "{text}");
    // the first programs of the registry are in view, each with what it is for, how to install it by hand, and its licence
    for program in &tools::ALL[..2] {
        assert!(text.contains(program.label), "{} in {text}", program.label);
        assert!(text.contains(program.role), "{} in {text}", program.role);
        assert!(
            text.contains(&format!(
                "By hand: {}",
                tools::install_ways(program).join("   or   ")
            )),
            "{text}"
        );
        assert!(
            text.contains(&format!("Licence: {}", program.licence)),
            "{text}"
        );
    }
    assert!(text.matches("not found").count() >= 2, "{text}");
    assert!(!app.programs_win.is_installing());
    // the command is in the palette, and needs no sequence in view
    let command = commands::all()
        .iter()
        .find(|c| matches!(c.action, Action::Programs))
        .expect("a command");
    assert_eq!(
        (command.group, command.title),
        ("Tools", "Programs HelixLoom calls…")
    );
    assert!(!commands::acts_on_the_sequence(&Action::Programs));
}

#[test]
fn install_shows_the_command_first_and_nothing_runs_until_the_person_says_so() {
    let (mut app, mut screen) = opened(only_homebrew, never);
    assert!(!screen.text().contains(RUN));
    // the first program of the registry is MAFFT
    screen.click(&mut app, INSTALL);
    screen.settle(&mut app);
    let text = screen.text();
    let plan = app.programs_win.asked.clone().expect("a plan that waits");
    assert_eq!(plan.line(), "brew install mafft");
    assert!(text.contains(&plan.notice(&tools::MAFFT)), "{text}");
    for button in [RUN, COPY, CANCEL] {
        assert!(text.contains(button), "{button} in {text}");
    }
    // while a plan waits, no other can be asked for: the Install button of the next program does nothing
    assert!(!app.programs_win.is_installing());
    let buttons = screen.find_all(INSTALL);
    assert!(
        buttons.len() >= 2,
        "an Install button for each program in view"
    );
    screen.click_at(&mut app, buttons[1]);
    screen.settle(&mut app);
    assert_eq!(
        app.programs_win.asked.as_ref().map(Plan::line).as_deref(),
        Some("brew install mafft")
    );

    // the command can be copied, to be run by hand
    screen.click(&mut app, COPY);
    screen.settle(&mut app);
    assert_eq!(screen.copied(), Some("brew install mafft"));
    assert!(app.programs_win.asked.is_some(), "copying is no answer");

    // and it can be left: nothing ran
    screen.click(&mut app, CANCEL);
    screen.settle(&mut app);
    assert!(app.programs_win.asked.is_none());
    assert!(!screen.text().contains(RUN));
}

#[test]
fn run_it_runs_the_plan_shows_what_the_manager_prints_and_looks_for_the_programs_again() {
    let (mut app, mut screen) = opened(only_homebrew, pretend);
    screen.click(&mut app, INSTALL);
    screen.settle(&mut app);
    assert!(app.platform.programs_known_at(&Place::Here));
    screen.click(&mut app, RUN);
    let end = Instant::now() + Duration::from_secs(30);
    loop {
        screen.frame(&mut app, vec![]);
        if !app.programs_win.is_installing() && app.programs_win.asked.is_none() {
            break;
        }
        assert!(Instant::now() < end, "gave up waiting");
        std::thread::sleep(Duration::from_millis(10));
    }
    assert_eq!(*ASKED.lock().expect("the list"), ["brew install mafft"]);
    assert_eq!(
        app.programs_win.message,
        "brew install mafft ended well: looking for MAFFT again."
    );
    // what was found before is forgotten, so that the programs are looked for again
    assert!(!app.platform.programs_known_at(&Place::Here));
    screen.frame(&mut app, vec![]);
    let text = screen.text();
    for line in ["==> Fetching the program", "==> Pouring it"] {
        assert!(text.contains(line), "{line} in {text}");
    }
}

#[test]
fn an_install_that_fails_says_why_and_a_computer_with_no_package_manager_says_what_to_do_by_hand() {
    let (mut app, mut screen) = opened(only_homebrew, failing);
    screen.click(&mut app, INSTALL);
    screen.settle(&mut app);
    screen.click(&mut app, RUN);
    let end = Instant::now() + Duration::from_secs(30);
    while app.programs_win.is_installing() {
        screen.frame(&mut app, vec![]);
        assert!(Instant::now() < end, "gave up waiting");
        std::thread::sleep(Duration::from_millis(10));
    }
    assert_eq!(
        app.programs_win.message,
        "brew stopped (exit status: 1): Error: no formula of that name"
    );
    // what was found before is still what is shown: nothing was looked for again
    assert!(app.platform.programs_known_at(&Place::Here));

    let (mut app, mut screen) = opened(no_manager, never);
    screen.click(&mut app, INSTALL);
    screen.settle(&mut app);
    assert!(app.programs_win.asked.is_none());
    assert!(
        app.programs_win
            .message
            .starts_with("No package manager was found on this computer"),
        "{}",
        app.programs_win.message
    );
    assert!(!screen.text().contains(RUN));
}

#[test]
fn a_program_that_is_found_shows_where_and_its_version_and_has_no_install_button() {
    let found = Installed {
        path: PathBuf::from("/opt/homebrew/bin/hmmscan"),
        version: "# HMMER 3.4 (Aug 2023); http://hmmer.org/".to_string(),
    };
    assert_eq!(
        found_line(&found),
        "/opt/homebrew/bin/hmmscan  (version 3.4)"
    );
    let plain = Installed {
        path: PathBuf::from("/usr/bin/tool"),
        version: "no number here".to_string(),
    };
    assert_eq!(found_line(&plain), "/usr/bin/tool");
}

#[test]
fn the_view_named_programs_starts_with_the_window_open() {
    let mut app = App::new(Startup {
        view: "programs".into(),
        ..Startup::default()
    });
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert!(app.programs_win.open);
    assert!(screen.text().contains(TITLE));
}
