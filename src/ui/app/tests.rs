//! The keyboard shortcuts of the window.

use super::*;
use crate::ui::headless::frame_input;
use eframe::egui::{Event, Key, Modifiers};

/// What one press of `key` with `modifiers` asks the window to do.
fn pressed(key: Key, modifiers: Modifiers) -> Vec<Action> {
    let ctx = egui::Context::default();
    let app = App::new(Startup::default());
    let mut input = frame_input(0);
    input.events.push(Event::Key {
        key,
        physical_key: Some(key),
        pressed: true,
        repeat: false,
        modifiers,
    });
    let mut actions = Vec::new();
    let mut output = ctx.run_ui(input, |ui| app.shortcuts(ui.ctx(), &mut actions));
    output.textures_delta.clear();
    actions
}

#[test]
fn command_one_and_two_open_the_two_calculators() {
    let one = pressed(Key::Num1, Modifiers::COMMAND);
    assert!(matches!(one.as_slice(), [Action::Calculator]), "{one:?}");
    let two = pressed(Key::Num2, Modifiers::COMMAND);
    assert!(
        matches!(two.as_slice(), [Action::MolecularCalculator]),
        "{two:?}"
    );
}

#[test]
fn command_alt_and_a_number_goes_to_a_workspace_and_opens_no_calculator() {
    for (key, workspace) in [
        (Key::Num1, Workspace::ALL[0]),
        (Key::Num2, Workspace::ALL[1]),
        (Key::Num9, Workspace::ALL[8]),
    ] {
        let actions = pressed(key, Modifiers::COMMAND | Modifiers::ALT);
        assert!(
            matches!(actions.as_slice(), [Action::Go(w)] if *w == workspace),
            "{key:?}: {actions:?}"
        );
    }
}

#[test]
fn command_k_opens_the_palette_and_annotating_moved_to_shift_command_k() {
    let k = pressed(Key::K, Modifiers::COMMAND);
    assert!(matches!(k.as_slice(), [Action::Palette]), "{k:?}");
    let shifted = pressed(Key::K, Modifiers::COMMAND | Modifiers::SHIFT);
    assert!(
        matches!(shifted.as_slice(), [Action::Annotate]),
        "{shifted:?}"
    );
}

#[test]
fn what_works_on_the_sequence_takes_the_person_to_it_from_any_page() {
    let ctx = egui::Context::default();
    let mut app = App::new(Startup::default());
    assert_eq!(app.workspace, Workspace::Home);
    app.run(Action::Upper, &ctx);
    assert_eq!(app.workspace, Workspace::Sequences);
    // a tool that does not need the sequence leaves the page where it is
    app.workspace = Workspace::Trees;
    app.run(Action::Calculator, &ctx);
    assert_eq!(app.workspace, Workspace::Trees);
}

#[test]
fn the_colours_follow_the_system_until_light_or_dark_is_chosen_and_follow_it_again_after() {
    use crate::ui::state::ThemeMode;
    let ctx = egui::Context::default();
    let mut app = App::new(Startup::default());
    let frame = |app: &mut App, n: u32, system: Option<egui::Theme>| {
        let mut input = frame_input(n);
        input.system_theme = system;
        let mut output = ctx.run_ui(input, |ui| app.draw(ui));
        output.textures_delta.clear();
        (app.prefs.dark, ctx.global_style().visuals.dark_mode)
    };
    assert_eq!(
        frame(&mut app, 0, None),
        (false, false),
        "no word from the system: light"
    );
    assert_eq!(frame(&mut app, 1, Some(egui::Theme::Light)), (false, false));
    // the system turns dark, and the window follows without being asked (a Mac does this in the evening)
    assert_eq!(frame(&mut app, 2, Some(egui::Theme::Dark)), (true, true));
    assert_eq!(frame(&mut app, 3, Some(egui::Theme::Light)), (false, false));
    // a choice of the person's is kept whatever the system does, and giving it up follows the system again
    app.prefs.theme = ThemeMode::Dark;
    assert_eq!(frame(&mut app, 4, Some(egui::Theme::Light)), (true, true));
    app.prefs.theme = ThemeMode::Light;
    assert_eq!(frame(&mut app, 5, Some(egui::Theme::Dark)), (false, false));
    app.prefs.theme = ThemeMode::System;
    assert_eq!(frame(&mut app, 6, Some(egui::Theme::Dark)), (true, true));
}
