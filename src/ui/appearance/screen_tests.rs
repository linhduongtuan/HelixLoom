//! The Appearance window and what it does to the whole window, drawn frame by frame as a person would meet it.

use super::*;
use crate::ui::app::{Action, App, Startup};
use crate::ui::testing::Screen;
use eframe::egui::{Event, Key, Modifiers};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, Ordering};
use web_time::Instant;

static NEXT: AtomicU32 = AtomicU32::new(0);

/// A settings folder of its own for one test, removed when the test ends.
struct Folder(PathBuf);

impl Folder {
    fn new() -> Folder {
        let path = std::env::temp_dir().join(format!(
            "helixloom-appearance-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }

    fn file(&self) -> PathBuf {
        self.0.join("appearance.txt")
    }

    /// A start with the Appearance window open (a view that is named also keeps the tour away) that keeps its settings here.
    fn start(&self) -> Startup {
        Startup {
            view: "appearance".to_string(),
            settings_dir: Some(self.0.clone()),
            ..Startup::default()
        }
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn open_window() -> (App, Screen) {
    let mut app = App::new(Startup {
        view: "appearance".to_string(),
        ..Startup::default()
    });
    let mut screen = Screen::new();
    screen.settle(&mut app);
    (app, screen)
}

fn key(key: Key, modifiers: Modifiers) -> Event {
    Event::Key {
        key,
        physical_key: Some(key),
        pressed: true,
        repeat: false,
        modifiers,
    }
}

#[test]
fn the_view_opens_the_window_with_every_setting_in_it() {
    let (_, screen) = open_window();
    let text = screen.text();
    for word in [
        "Appearance",
        "THEME",
        "Same as the system",
        "Light",
        "Dark",
        "SIZE",
        "80%",
        "100%",
        "200%",
        "SPACING",
        "Comfortable",
        "Compact",
        "ACCENT COLOUR",
        "Indigo",
        "CONTRAST",
        "High contrast",
        "COLOURS OF BASES AND RESIDUES",
        "Colour-blind safe",
        "Back to the defaults",
    ] {
        assert!(text.contains(word), "no \"{word}\" in\n{text}");
    }
    // the sample of colours shows the letters of bases and of one amino acid of each group
    for letter in ["A", "C", "G", "T", "L", "K", "D", "S", "P", "Y"] {
        assert!(screen.find(letter).is_some(), "no {letter} in the sample");
    }
}

#[test]
fn nothing_is_said_to_be_kept_when_there_is_no_settings_folder_and_it_is_when_there_is() {
    let (_, screen) = open_window();
    assert!(screen.text().contains("not kept"), "{}", screen.text());
    let folder = Folder::new();
    let mut app = App::new(folder.start());
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert!(
        screen.text().contains("kept for the next time"),
        "{}",
        screen.text()
    );
}

#[test]
fn each_choice_is_used_by_the_window_at_once() {
    let (mut app, mut screen) = open_window();
    assert_eq!(screen.ctx().zoom_factor(), 1.0);
    assert_eq!(screen.ctx().global_style().spacing.interact_size.y, 26.0);

    screen.click(&mut app, "150%");
    screen.settle(&mut app);
    assert_eq!(app.prefs.appearance.scale, 150);
    assert_eq!(screen.ctx().zoom_factor(), 1.5);

    screen.click(&mut app, "Compact");
    screen.settle(&mut app);
    assert_eq!(app.prefs.appearance.density, Density::Compact);
    assert_eq!(
        screen.ctx().global_style().spacing.interact_size.y,
        MIN_TARGET
    );

    screen.click(&mut app, "High contrast");
    screen.settle(&mut app);
    assert_eq!(app.prefs.appearance.contrast, Contrast::High);
    assert_eq!(screen.ctx().global_style().visuals.window_stroke.width, 1.5);
    assert_eq!(
        app.prefs.look().text,
        Color32::BLACK,
        "the chrome is drawn with it too"
    );

    screen.click(&mut app, "Colour-blind safe");
    screen.settle(&mut app);
    assert_eq!(app.prefs.appearance.colours, Colours::Safe);
    assert!(app.prefs.look().safe_colours);

    screen.click(&mut app, "Dark");
    screen.settle(&mut app);
    assert_eq!(app.prefs.theme, ThemeMode::Dark);
    assert!(app.prefs.dark);
}

#[test]
fn the_accent_is_chosen_by_its_colour_and_the_chrome_takes_it() {
    let (mut app, mut screen) = open_window();
    // the swatches have no text: they are found by where the program drew them (the name of the one chosen is beside them), once the
    // window, which fades in, is drawn in full colour
    for _ in 0..12 {
        screen.frame(&mut app, vec![]);
    }
    assert!(screen.text().contains("Indigo"));
    let swatches: Vec<_> = Accent::ALL.iter().map(|a| a.swatch(false)).collect();
    let spots = screen.circles_of(&swatches);
    assert_eq!(spots.len(), Accent::ALL.len(), "{spots:?}");
    let at = spots[2];
    screen.click_at(&mut app, at);
    screen.settle(&mut app);
    assert_eq!(app.prefs.appearance.accent, Accent::Teal);
    assert!(screen.text().contains("Teal"));
    assert_eq!(app.prefs.look().accent, Accent::Teal.swatch(false));
}

#[test]
fn back_to_the_defaults_puts_every_setting_back_but_not_the_theme() {
    let (mut app, mut screen) = open_window();
    app.prefs.appearance =
        Appearance::from_spec("scale=110,density=compact,accent=purple,contrast=high,colours=safe");
    app.prefs.theme = ThemeMode::Dark;
    screen.settle(&mut app);
    screen.click(&mut app, "Back to the defaults");
    screen.settle(&mut app);
    assert_eq!(app.prefs.appearance, Appearance::default());
    assert_eq!(screen.ctx().zoom_factor(), 1.0);
    assert_eq!(app.prefs.theme, ThemeMode::Dark);
}

#[test]
fn the_choices_are_kept_in_the_settings_folder_and_read_at_the_next_start() {
    let folder = Folder::new();
    {
        let mut app = App::new(folder.start());
        let mut screen = Screen::new();
        screen.settle(&mut app);
        assert!(!folder.file().exists(), "a start writes nothing");
        // the window takes a moment to lay itself out again at a new size
        for choice in ["125%", "Compact", "Dark"] {
            screen.click(&mut app, choice);
            screen.settle(&mut app);
        }
    }
    let text = std::fs::read_to_string(folder.file()).expect("a file");
    assert_eq!(
        text,
        "theme dark\nscale 125\ndensity compact\naccent indigo\ncontrast standard\ncolours standard\n"
    );
    let mut again = App::new(Startup {
        settings_dir: Some(folder.0.clone()),
        ..Startup::default()
    });
    // reading it back writes nothing again: the time of the file is taken before the first frame, and the time is let pass
    let before = std::fs::metadata(folder.file())
        .and_then(|m| m.modified())
        .expect("a time");
    std::thread::sleep(std::time::Duration::from_millis(30));
    let mut screen = Screen::new();
    screen.settle(&mut again);
    assert_eq!(again.prefs.theme, ThemeMode::Dark);
    assert_eq!(again.prefs.appearance.scale, 125);
    assert_eq!(again.prefs.appearance.density, Density::Compact);
    assert_eq!(
        screen.ctx().zoom_factor(),
        1.25,
        "and used by the window from its first frames"
    );
    screen.settle(&mut again);
    let after = std::fs::metadata(folder.file())
        .and_then(|m| m.modified())
        .expect("a time");
    assert_eq!(before, after);
}

#[test]
fn the_command_line_is_put_on_what_the_folder_holds_and_a_run_given_it_keeps_nothing() {
    let folder = Folder::new();
    let held = "theme light\nscale 90\naccent teal\n";
    std::fs::write(folder.file(), held).expect("a file");
    let mut app = App::new(Startup {
        dark: true,
        appearance: Some("scale=150".to_string()),
        settings_dir: Some(folder.0.clone()),
        view: "appearance".to_string(),
        ..Startup::default()
    });
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert_eq!(app.prefs.theme, ThemeMode::Dark);
    assert_eq!(
        (app.prefs.appearance.scale, app.prefs.appearance.accent),
        (150, Accent::Teal),
        "what the command line does not name stays as it was kept"
    );
    assert!(
        screen.text().contains("none of them is kept"),
        "{}",
        screen.text()
    );
    // a choice made in this run is used, but neither it nor what was forced is written over what was chosen before
    screen.click(&mut app, "Compact");
    screen.settle(&mut app);
    assert_eq!(app.prefs.appearance.density, Density::Compact);
    assert_eq!(
        std::fs::read_to_string(folder.file()).expect("a file"),
        held
    );
}

#[test]
fn a_write_that_fails_is_tried_again_later_and_the_window_says_that_it_is_not_kept() {
    use std::time::Duration;
    let folder = Folder::new();
    // a settings folder that cannot be written: where it should be there is a file
    let blocked = folder.0.join("settings");
    std::fs::write(&blocked, "in the way").expect("a file");
    let mut app = App::new(Startup {
        view: "appearance".to_string(),
        settings_dir: Some(blocked.clone()),
        ..Startup::default()
    });
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert!(
        screen.text().contains("kept for the next time"),
        "{}",
        screen.text()
    );
    screen.click(&mut app, "125%");
    screen.settle(&mut app);
    assert!(app.keeper.appearance_failed);
    assert!(
        screen.text().contains("could not be kept"),
        "{}",
        screen.text()
    );
    // the next try is not at once: asked a moment later it says how long there is to wait, and writes nothing
    let (theme, now) = (app.prefs.theme, app.prefs.appearance);
    let soon = Instant::now() + Duration::from_secs(1);
    let wait = app
        .keeper
        .keep_appearance(theme, &now, soon)
        .expect("a try later");
    assert!(
        wait < Duration::from_millis(4500),
        "{wait:?}: counted from the failure, not from now"
    );
    assert!(blocked.is_file());
    // the folder can be written now: the try after the wait goes through, and the window says it is kept again
    std::fs::remove_file(&blocked).expect("removed");
    std::fs::create_dir(&blocked).expect("a folder");
    let later = Instant::now() + Duration::from_secs(60);
    assert_eq!(app.keeper.keep_appearance(theme, &now, later), None);
    assert!(!app.keeper.appearance_failed);
    assert!(blocked.join("appearance.txt").is_file());
    screen.settle(&mut app);
    assert!(
        screen.text().contains("kept for the next time"),
        "{}",
        screen.text()
    );
}

#[test]
fn a_file_that_is_damaged_is_started_from_the_defaults_and_replaced_by_the_first_change() {
    let folder = Folder::new();
    std::fs::write(folder.file(), "scale 7000\nnonsense\n\u{0}\u{1}").expect("a file");
    let mut app = App::new(folder.start());
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert_eq!(app.prefs.appearance, Appearance::default());
    screen.click(&mut app, "110%");
    screen.settle(&mut app);
    let text = std::fs::read_to_string(folder.file()).expect("a file");
    assert!(text.contains("scale 110"), "{text}");
}

#[test]
fn the_keys_step_through_the_sizes_one_at_a_time_and_the_windows_own_zoom_keys_are_off() {
    let (mut app, mut screen) = open_window();
    assert!(
        !screen.ctx().options(|o| o.zoom_with_keyboard),
        "egui's own steps of a tenth are not used"
    );
    for (press, expected) in [
        (key(Key::Equals, Modifiers::COMMAND), 1.1),
        (key(Key::Equals, Modifiers::COMMAND), 1.25),
        (key(Key::Plus, Modifiers::COMMAND | Modifiers::SHIFT), 1.5),
        (key(Key::Minus, Modifiers::COMMAND), 1.25),
        (key(Key::Num0, Modifiers::COMMAND), 1.0),
    ] {
        screen.frame(&mut app, vec![press]);
        screen.settle(&mut app);
        assert!(
            (screen.ctx().zoom_factor() - expected).abs() < 1e-6,
            "{expected}: {}",
            screen.ctx().zoom_factor()
        );
    }
    assert!(app.status.contains("100%"), "{}", app.status);
}

#[test]
fn the_commands_of_the_palette_and_the_menu_do_the_same_as_the_keys() {
    let (mut app, mut screen) = open_window();
    let ctx = screen.ctx().clone();
    app.run(Action::Zoom(1), &ctx);
    app.run(Action::Zoom(1), &ctx);
    assert_eq!(app.prefs.appearance.scale, 125);
    app.run(Action::Zoom(-1), &ctx);
    assert_eq!(app.prefs.appearance.scale, 110);
    app.run(Action::Zoom(0), &ctx);
    assert_eq!(app.prefs.appearance.scale, 100);
    app.appearance_open = false;
    app.run(Action::Appearance, &ctx);
    assert!(app.appearance_open);
    // from another page too, without being taken to the sequence
    app.workspace = crate::ui::platform::Workspace::Trees;
    app.run(Action::Zoom(1), &ctx);
    assert_eq!(app.workspace, crate::ui::platform::Workspace::Trees);
    screen.settle(&mut app);
}

#[test]
fn every_page_and_the_window_draw_at_every_size_and_setting_without_a_panic() {
    use crate::ui::platform::Workspace;
    for scale in SCALES {
        let spec =
            format!("scale={scale},density=compact,accent=graphite,contrast=high,colours=safe");
        for dark in [false, true] {
            let mut app = App::new(Startup {
                view: "appearance".to_string(),
                dark,
                appearance: Some(spec.clone()),
                ..Startup::default()
            });
            let mut screen = Screen::new();
            for ws in Workspace::ALL {
                app.workspace = ws;
                screen.settle(&mut app);
            }
            assert!((screen.ctx().zoom_factor() - f32::from(scale) / 100.0).abs() < 1e-6);
        }
    }
}

/// The longest piece of text in the last frame that is nothing but base letters: a row of the alignment, or a stretch of one.
fn longest_run_of_bases(screen: &Screen) -> usize {
    screen
        .text()
        .lines()
        .filter(|line| !line.is_empty() && line.bytes().all(|b| b"ACGT".contains(&b)))
        .map(str::len)
        .max()
        .unwrap_or(0)
}

/// An app showing the example alignment with `spec`, drawn until the alignment (which is made in the background) is on the screen.
fn alignment_with(spec: &str) -> Screen {
    let mut app = App::new(Startup {
        view: "alignment".to_string(),
        dark: true,
        appearance: Some(spec.to_string()),
        ..Startup::default()
    });
    let mut screen = Screen::new();
    for _ in 0..600 {
        screen.frame(&mut app, vec![]);
        if longest_run_of_bases(&screen) > 0 {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(5));
    }
    screen.settle(&mut app);
    screen
}

#[test]
fn the_alignment_letters_keep_the_text_colour_unless_the_safe_colours_need_black_and_white() {
    let standard = longest_run_of_bases(&alignment_with(""));
    let safe = longest_run_of_bases(&alignment_with("colours=safe"));
    // a row is one piece of text in one colour, as it always was
    assert!(standard >= 60, "{standard}");
    // and with the safe colours each stretch of a colour that reads is a piece of its own
    assert!(safe > 0 && safe < 20, "{safe}");
    // high contrast chooses black or white by the cell as well, whichever colours the cells have
    let high = longest_run_of_bases(&alignment_with("contrast=high"));
    assert!(high > 0 && high < 40, "{high}");
}

#[test]
fn a_swatch_is_reached_with_the_tab_key_has_a_ring_and_is_chosen_with_enter() {
    let (mut app, mut screen) = open_window();
    for _ in 0..12 {
        screen.frame(&mut app, vec![]);
    }
    assert!(
        screen.circles_of_radius(window::FOCUS_RING).is_empty(),
        "nothing has the keyboard yet"
    );
    let swatches: Vec<_> = Accent::ALL.iter().map(|a| a.swatch(false)).collect();
    assert_eq!(screen.circles_of(&swatches).len(), Accent::ALL.len());
    // Tab goes from one control to the next; the ring shows when it comes to a swatch (the window scrolls to show it, so the places
    // of the swatches are looked up again in every frame)
    let mut reached = None;
    for _ in 0..400 {
        screen.frame(&mut app, vec![key(Key::Tab, Modifiers::NONE)]);
        let (spots, rings) = (
            screen.circles_of(&swatches),
            screen.circles_of_radius(window::FOCUS_RING),
        );
        reached = spots
            .iter()
            .position(|spot| rings.iter().any(|ring| (*ring - *spot).length() < 1.0));
        if reached.is_some() {
            break;
        }
    }
    let reached = reached.expect("no swatch was reached by Tab");
    // and Enter chooses it, as a click would
    let wanted = Accent::ALL[reached];
    screen.frame(&mut app, vec![key(Key::Enter, Modifiers::NONE)]);
    screen.settle(&mut app);
    assert_eq!(app.prefs.appearance.accent, wanted);
    // only one swatch has the ring at a time
    assert_eq!(screen.circles_of_radius(window::FOCUS_RING).len(), 1);
}

#[test]
fn the_window_says_where_the_choices_are_kept_and_what_could_not_be() {
    use crate::ui::autosave::Keeper;
    use crate::ui::session::{Memory, Store};
    // no settings at all
    let mut bare = App::new(Startup::default());
    assert!(bare
        .appearance_note()
        .contains("There is no settings folder"));
    // a folder of the computer
    let folder = Folder::new();
    let mut desktop = App::new(folder.start());
    assert_eq!(
        desktop.appearance_note(),
        "Each choice is used at once and kept for the next time."
    );
    desktop.keeper.appearance_failed = true;
    assert!(desktop
        .appearance_note()
        .contains("the settings folder cannot be written"));
    // a page of a browser
    let table = std::rc::Rc::new(Memory::default());
    bare.keeper = Keeper::with(Some(Store::on(table)));
    assert_eq!(
        bare.appearance_note(),
        "Each choice is used at once and kept in this browser for the next visit."
    );
    bare.keeper.appearance_failed = true;
    let failed = bare.appearance_note();
    assert!(
        failed.contains("the browser's storage for this site cannot be written"),
        "{failed}"
    );
    assert!(!failed.contains("settings folder"), "{failed}");
    // a run that was given settings says that first, wherever they would be kept
    bare.appearance_locked = true;
    assert!(bare.appearance_note().contains("none of them is kept"));
}
