use super::update::{Msg, State};
use crate::ui::app::{App, Startup};
use crate::ui::testing::Screen;
use crate::update::{Asset, Downloaded, Outcome, Release};
use eframe::egui;
use std::path::PathBuf;

fn release() -> Release {
    Release {
        tag: "v9.9.9".into(),
        notes: "What is new in this one".into(),
        page: "https://github.com/linhduongtuan/HelixLoom/releases/tag/v9.9.9".into(),
        assets: vec![
            Asset {
                name: "HelixLoom-v9.9.9-macos-arm64.dmg".into(),
                url: "https://github.com/linhduongtuan/HelixLoom/releases/download/v9.9.9/HelixLoom-v9.9.9-macos-arm64.dmg".into(),
                size: 9_000_000,
            },
            Asset {
                name: "HelixLoom-v9.9.9-windows-x64-NOT-TRIED.zip".into(),
                url: "https://github.com/linhduongtuan/HelixLoom/releases/download/v9.9.9/HelixLoom-v9.9.9-windows-x64-NOT-TRIED.zip".into(),
                size: 9_000_000,
            },
        ],
    }
}

fn saved() -> Downloaded {
    Downloaded {
        path: PathBuf::from("/Users/x/Downloads/HelixLoom-v9.9.9-macos-arm64.dmg"),
        verified: true,
    }
}

#[test]
fn the_state_follows_what_the_work_comes_back_with() {
    let check = |state: State, msg| state.after(msg);
    assert_eq!(
        check(
            State::Checking,
            Msg::Checked(Ok(Outcome::UpToDate {
                latest: "v0.2.0".into()
            }))
        ),
        State::UpToDate("v0.2.0".into())
    );
    assert_eq!(
        check(
            State::Checking,
            Msg::Checked(Ok(Outcome::Available(release())))
        ),
        State::Available(release())
    );
    assert_eq!(
        check(State::Checking, Msg::Checked(Err("no network".into()))),
        State::Failed("no network".into())
    );
    assert_eq!(
        check(State::Downloading(release()), Msg::Downloaded(Ok(saved()))),
        State::Downloaded(release(), saved())
    );
    assert_eq!(
        check(
            State::Downloading(release()),
            Msg::Downloaded(Err("cut off".into()))
        ),
        State::Failed("cut off".into())
    );
    // an answer that nothing was waiting for changes nothing
    assert_eq!(
        check(State::Idle, Msg::Downloaded(Ok(saved()))),
        State::Idle
    );
    assert_eq!(
        check(
            State::Available(release()),
            Msg::Checked(Err("late".into()))
        ),
        State::Available(release())
    );
    assert_eq!(
        check(
            State::UpToDate("v1".into()),
            Msg::Downloaded(Err("late".into()))
        ),
        State::UpToDate("v1".into())
    );
}

#[test]
fn the_window_says_what_each_state_is_and_asks_nothing_until_a_button_is_pressed() {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    app.update_win.open = true;
    screen.settle(&mut app);
    let text = screen.text();
    assert!(
        text.contains(&format!("This is HelixLoom {}", env!("CARGO_PKG_VERSION"))),
        "{text}"
    );
    assert!(
        text.contains("Check for updates")
            && text.contains("Nothing is asked of any server until you press"),
        "{text}"
    );
    assert_eq!(
        app.update_win.state,
        State::Idle,
        "opening the window checks nothing"
    );
    for (state, words) in [
        (
            State::UpToDate("v0.2.0".into()),
            "This is the latest version",
        ),
        (
            State::Failed("no network".into()),
            "Could not do that: no network",
        ),
        (State::Checking, "Asking GitHub"),
        (
            State::Available(release()),
            "A later version is out: v9.9.9",
        ),
        (State::Downloading(release()), "Downloading"),
        (State::Downloaded(release(), saved()), "checksum matches"),
    ] {
        app.update_win.state = state;
        screen.settle(&mut app);
        let text = screen.text();
        assert!(text.contains(words), "{words}: {text}");
    }
    // with a later version out: the notes, the release page, and a file for this computer or the reason there is none
    app.update_win.state = State::Available(release());
    screen.settle(&mut app);
    let text = screen.text();
    assert!(
        text.contains("What is new in this one") && text.contains("Open the release page"),
        "{text}"
    );
    assert!(
        text.contains("Download HelixLoom-v9.9.9")
            || text.contains("no file for this computer's system"),
        "{text}"
    );
}

#[test]
fn the_command_and_the_menu_item_open_the_window() {
    let mut app = App::new(Startup::default());
    assert!(!app.update_win.open);
    let found = crate::ui::commands::all()
        .iter()
        .find(|c| c.title.to_lowercase().contains("check for updates"))
        .expect("the palette has the command");
    // run as the palette runs it: the window is open, and nothing has been asked of GitHub
    app.run(found.action.clone(), &egui::Context::default());
    assert!(app.update_win.open);
    assert_eq!(app.update_win.state, State::Idle);
}
