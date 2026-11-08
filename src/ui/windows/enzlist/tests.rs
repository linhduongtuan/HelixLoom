use super::*;
use crate::ui::app::{Action, Startup};
use crate::ui::commands::acts_on_the_sequence;
use crate::ui::headless::frame_input;
use eframe::egui::{
    epaint::Shape, pos2, Context, Event, FullOutput, Key, Modifiers, OutputCommand, PointerButton,
    Rect,
};
use std::cell::Cell;

thread_local! {
    /// The frames this test has drawn: the time of a frame must never go back.
    static FRAMES: Cell<u32> = const { Cell::new(0) };
}

fn app() -> App {
    App::new(Startup::default())
}

fn frame(app: &mut App, ctx: &Context, events: Vec<Event>) -> FullOutput {
    let n = FRAMES.with(|f| f.replace(f.get() + 1));
    let mut input = frame_input(n);
    input.events = events;
    let mut output = ctx.run_ui(input, |ui| app.draw(ui));
    output.textures_delta.clear();
    output
}

/// Where `text` was drawn in a frame: the first piece of text that holds it.
fn text_rect(shapes: &[Shape], text: &str) -> Option<Rect> {
    shapes.iter().find_map(|shape| match shape {
        Shape::Text(t) if t.galley.text().contains(text) => Some(t.visual_bounding_rect()),
        Shape::Vec(inner) => text_rect(inner, text),
        _ => None,
    })
}

/// A real click on the middle of the piece of text `text` that a frame drew (the button or choice that holds it); the frame
/// that took the click.
fn click_on(app: &mut App, ctx: &Context, shown: &FullOutput, text: &str) -> FullOutput {
    let shapes: Vec<Shape> = shown.shapes.iter().map(|s| s.shape.clone()).collect();
    let rect = text_rect(&shapes, text).unwrap_or_else(|| panic!("{text} was not drawn"));
    let at = pos2(rect.center().x, rect.center().y);
    let button = |pressed| Event::PointerButton {
        pos: at,
        button: PointerButton::Primary,
        pressed,
        modifiers: Modifiers::NONE,
    };
    frame(app, ctx, vec![Event::PointerMoved(at)]);
    frame(app, ctx, vec![button(true), button(false)])
}

/// Two frames, so that a window that has just opened is laid out and can be clicked; the last one.
fn settle(app: &mut App, ctx: &Context) -> FullOutput {
    frame(app, ctx, vec![]);
    frame(app, ctx, vec![])
}

fn run_of(app: &App) -> &ListRun {
    app.enz_list.run.as_ref().expect("a list was made")
}

/// What follows the heading in a list's text: the lines the original's `list_enzymes_text` made.
fn body_of(app: &App) -> String {
    let (name, len) = (&app.docs[0].doc.name, app.docs[0].doc.len());
    let run = run_of(app);
    let head = enzlist::report(PROGRAM, name, run.choice.which, None, len, "");
    run.report
        .strip_prefix(&head)
        .unwrap_or(&run.report)
        .to_string()
}

/// The enzymes of the demo plasmid with how many sites each has.
fn site_counts(app: &mut App) -> Vec<(String, usize)> {
    let sites = app.docs[0].sites(&app.enzymes);
    let names = app.enzymes.iter().map(|e| e.name.clone());
    names.zip(sites.iter().map(Vec::len)).collect()
}

#[test]
fn site_usage_lists_every_enzyme_with_its_number_of_sites_in_columns() {
    let mut app = app();
    app.run_quick_list(Quick::SiteUsage);

    let counts = site_counts(&mut app);
    let run = run_of(&app);
    let (name, len) = (app.docs[0].doc.name.clone(), app.docs[0].doc.len());
    let heading = format!("{PROGRAM}\n{name} \nFrom 0 to {len}. \n\nAll Enzymes\n\n");

    assert!(run.report.starts_with(&heading));
    assert_eq!(run.shown, app.enzymes.len());
    assert_eq!(run.title, format!("{name} All Enzymes"));
    let body = body_of(&app);
    let first = body.lines().next().expect("a first row");
    assert_eq!(first.len(), 6 * 20, "six columns of twelve and eight");
    let (eco, n) = counts.iter().find(|(e, _)| e == "EcoRI").expect("EcoRI");
    assert!(body.contains(&format!("{eco:<12}{n:<8}")));
}

#[test]
fn unique_and_absent_lists_hold_the_enzymes_that_cut_once_and_never() {
    let mut app = app();
    let counts = site_counts(&mut app);
    let named = |wanted: usize| {
        let mut names: Vec<&str> = counts
            .iter()
            .filter(|(_, n)| *n == wanted)
            .map(|(e, _)| e.as_str())
            .collect();
        names.sort_unstable();
        names
    };

    app.run_quick_list(Quick::UniqueSites);
    let body = body_of(&app);
    let mut unique: Vec<&str> = body.split_whitespace().collect();
    unique.sort_unstable();
    assert_eq!(unique, named(1));
    assert!(run_of(&app).title.ends_with(" Unique Sites"));

    app.run_quick_list(Quick::AbsentSites);
    let body = body_of(&app);
    let mut absent: Vec<&str> = body.split_whitespace().collect();
    absent.sort_unstable();
    assert_eq!(absent, named(0));
    assert!(run_of(&app).title.ends_with(" Absent Sites"));
}

#[test]
fn the_dialog_lists_with_its_choices_closes_and_keeps_them_for_next_time() {
    let mut app = app();
    app.open_list_dialog();
    app.enz_list.format = Format::BySite;
    app.enz_list.which = Which::Unique;
    app.enz_list.columns = 3;
    app.list_dialog_ok();

    assert!(!app.enz_list.dialog_open && app.enz_list.result_open);
    let choice = run_of(&app).choice;
    assert_eq!(
        (choice.format, choice.which, choice.columns),
        (Format::BySite, Which::Unique, 3)
    );

    app.enz_list.selection_only = true;
    app.open_list_dialog();
    assert!(app.enz_list.dialog_open);
    assert_eq!(app.enz_list.format, Format::BySite, "the format is kept");
    assert_eq!(app.enz_list.columns, 3, "so are the columns");
    assert!(!app.enz_list.selection_only, "DNA starts at All");
}

#[test]
fn a_quick_list_uses_the_columns_the_dialog_has() {
    let mut app = app();
    app.enz_list.columns = 2;

    app.run_quick_list(Quick::AbsentSites);

    assert_eq!(run_of(&app).choice.columns, 2);
    let first = body_of(&app).lines().next().expect("a row").len();
    assert_eq!(first, 2 * 15);
}

#[test]
fn selected_dna_only_lists_the_sites_on_the_selection_and_the_heading_says_where() {
    let mut app = app();
    app.docs[0].select(326, 446);
    app.enz_list.format = Format::Positions;
    app.enz_list.selection_only = true;
    app.list_dialog_ok();

    let report = &run_of(&app).report;
    assert!(
        report.contains("\nFrom 326 to 446. \n\nAll Enzymes\n\n"),
        "{report}"
    );
    for line in body_of(&app).lines() {
        let places = line.get(37..).unwrap_or("");
        for p in places
            .split_whitespace()
            .map(|p| p.parse::<usize>().expect("a place"))
        {
            assert!((326..=446).contains(&p), "{line}");
        }
    }
}

#[test]
fn selected_dna_only_with_nothing_selected_is_all_of_the_dna() {
    let mut app = app();
    app.docs[0].place_caret(10);
    let len = app.docs[0].doc.len();
    app.enz_list.selection_only = true;

    app.list_dialog_ok();

    assert!(run_of(&app)
        .report
        .contains(&format!("\nFrom 0 to {len}. \n")));
}

#[test]
fn selected_enzymes_follow_the_enzyme_selector_and_none_selected_lists_nothing() {
    let mut app = app();
    app.enz_list.which = Which::Selected;
    app.enz_list.format = Format::Names;
    app.selected.clear();
    app.list_dialog_ok();
    assert_eq!(run_of(&app).shown, 0);
    assert_eq!(body_of(&app), "");

    app.selected = ["KpnI", "EcoRI"].iter().map(|s| s.to_string()).collect();
    app.list_dialog_ok();
    assert_eq!(run_of(&app).shown, 2);
    let names: Vec<String> = body_of(&app).split_whitespace().map(String::from).collect();
    assert_eq!(names, ["EcoRI", "KpnI"]);
}

const STALE: &str = "The sequence has changed since this list was made";

/// Whether a frame drew the note that the list is out of date.
fn says_stale(shown: &FullOutput) -> bool {
    let shapes: Vec<Shape> = shown.shapes.iter().map(|s| s.shape.clone()).collect();
    text_rect(&shapes, STALE).is_some()
}

#[test]
fn a_list_goes_stale_when_the_sequence_is_edited_and_listing_again_makes_a_new_one() {
    let (ctx, mut app) = (Context::default(), app());
    app.run_quick_list(Quick::UniqueSites);
    let choice = run_of(&app).choice;
    let before = settle(&mut app, &ctx);
    assert!(!says_stale(&before), "a list that is current has no note");

    app.docs[0].rev += 1;
    let after = frame(&mut app, &ctx, vec![]);
    assert!(says_stale(&after), "the note is drawn");

    app.make_list(choice);
    let again = frame(&mut app, &ctx, vec![]);
    assert!(!says_stale(&again), "listing again takes it away");
    assert_eq!(run_of(&app).basis.rev, app.docs[0].rev);
}

#[test]
fn a_list_goes_stale_when_the_topology_or_the_dam_dcm_switch_change_the_sites() {
    let (ctx, mut app) = (Context::default(), app());
    app.run_quick_list(Quick::UniqueSites);
    settle(&mut app, &ctx);

    app.docs[0].doc.methylated = !app.docs[0].doc.methylated;
    assert!(says_stale(&frame(&mut app, &ctx, vec![])), "Dam/Dcm");

    app.run_quick_list(Quick::UniqueSites);
    settle(&mut app, &ctx);
    app.docs[0].doc.circular = !app.docs[0].doc.circular;
    assert!(
        says_stale(&frame(&mut app, &ctx, vec![])),
        "circular or linear"
    );
}

#[test]
fn a_list_of_the_selected_enzymes_goes_stale_when_the_selector_changes_and_the_others_do_not() {
    let (ctx, mut app) = (Context::default(), app());
    app.selected = ["EcoRI"].iter().map(|s| s.to_string()).collect();
    app.enz_list.which = Which::Selected;
    app.list_dialog_ok();
    settle(&mut app, &ctx);

    app.selected.insert("KpnI".to_string());
    assert!(
        says_stale(&frame(&mut app, &ctx, vec![])),
        "the selected enzymes"
    );

    app.run_quick_list(Quick::UniqueSites);
    settle(&mut app, &ctx);
    app.selected.insert("SacI".to_string());
    assert!(
        !says_stale(&frame(&mut app, &ctx, vec![])),
        "a list of all the enzymes does not depend on it"
    );
}

#[test]
fn both_windows_draw_in_every_format() {
    let (ctx, mut app) = (Context::default(), app());
    for format in Format::ALL {
        app.enz_list.format = format;
        app.open_list_dialog();
        app.list_dialog_ok();
        app.enz_list.dialog_open = true;
        settle(&mut app, &ctx);
    }
    assert!(app.enz_list.result_open && app.enz_list.dialog_open);
}

#[test]
fn the_menu_and_the_palette_reach_the_lists() {
    assert!(acts_on_the_sequence(&Action::ListEnzymes));
    assert!(acts_on_the_sequence(&Action::QuickList(Quick::SiteUsage)));
    let labels: Vec<&str> = Quick::ALL.iter().map(|q| q.label()).collect();
    assert_eq!(labels, ["Site Usage", "Unique Sites", "Absent Sites"]);
    let commands = crate::ui::commands::all();
    let found = |typed: &str| {
        let hits = crate::ui::commands::search(commands, typed, 8);
        hits.into_iter()
            .map(|i| commands[i].action.clone())
            .collect::<Vec<Action>>()
    };
    assert!(found("list enzymes")
        .iter()
        .any(|a| matches!(a, Action::ListEnzymes)));
    assert!(found("quick list unique")
        .iter()
        .any(|a| matches!(a, Action::QuickList(Quick::UniqueSites))));
    assert!(found("absent")
        .iter()
        .any(|a| matches!(a, Action::QuickList(Quick::AbsentSites))));
    assert!(found("site usage")
        .iter()
        .any(|a| matches!(a, Action::QuickList(Quick::SiteUsage))));
}

#[test]
fn without_a_sequence_a_list_says_so() {
    let mut app = app();
    app.docs.clear();

    app.run_quick_list(Quick::SiteUsage);

    assert!(app.enz_list.run.is_none());
    assert_eq!(app.status, "Open a sequence first.");
}

#[test]
fn clicking_the_choices_and_ok_makes_the_list_that_was_asked_for_and_closes_the_dialog() {
    let (ctx, mut app) = (Context::default(), app());
    app.open_list_dialog();
    let shown = settle(&mut app, &ctx);

    let shown = click_on(&mut app, &ctx, &shown, "Unique Enzymes");
    assert_eq!(
        app.enz_list.which,
        Which::Unique,
        "the radio button took the click"
    );
    let shown = click_on(&mut app, &ctx, &shown, "Selected DNA only");
    assert!(app.enz_list.selection_only);
    click_on(&mut app, &ctx, &shown, "OK");

    assert!(!app.enz_list.dialog_open && app.enz_list.result_open);
    let choice = run_of(&app).choice;
    assert_eq!((choice.which, choice.selection_only), (Which::Unique, true));
}

#[test]
fn return_in_the_dialog_is_ok() {
    let (ctx, mut app) = (Context::default(), app());
    app.open_list_dialog();
    settle(&mut app, &ctx);
    let key = |pressed| Event::Key {
        key: Key::Enter,
        physical_key: Some(Key::Enter),
        pressed,
        repeat: false,
        modifiers: Modifiers::NONE,
    };

    frame(&mut app, &ctx, vec![key(true), key(false)]);

    assert!(!app.enz_list.dialog_open && app.enz_list.result_open);
}

#[test]
fn closing_the_dialog_lists_nothing() {
    let (ctx, mut app) = (Context::default(), app());
    app.open_list_dialog();
    app.enz_list.dialog_open = false;

    settle(&mut app, &ctx);

    assert!(app.enz_list.run.is_none() && !app.enz_list.result_open);
}

#[test]
fn list_again_and_copy_as_text_work_from_the_result_window() {
    let (ctx, mut app) = (Context::default(), app());
    app.run_quick_list(Quick::SiteUsage);
    settle(&mut app, &ctx);
    app.docs[0].rev += 1;
    let shown = frame(&mut app, &ctx, vec![]);
    assert!(says_stale(&shown), "the list is stale");

    let shown = click_on(&mut app, &ctx, &shown, "List again");
    assert_eq!(run_of(&app).basis.rev, app.docs[0].rev, "it was made again");
    let after = click_on(&mut app, &ctx, &shown, "Copy as text");

    let copied: Vec<&String> = after
        .platform_output
        .commands
        .iter()
        .filter_map(|c| match c {
            OutputCommand::CopyText(text) => Some(text),
            _ => None,
        })
        .collect();
    assert_eq!(copied, [&run_of(&app).report]);
}

#[test]
fn a_count_says_enzyme_for_one_and_enzymes_for_any_other_number() {
    assert_eq!(count_of(0), "0 enzymes");
    assert_eq!(count_of(1), "1 enzyme");
    assert_eq!(count_of(116), "116 enzymes");
    assert_eq!(count_of(1_200), "1,200 enzymes");
}

#[test]
fn the_message_after_a_list_counts_the_enzymes_in_it() {
    let mut app = app();
    app.enz_list.which = Which::Selected;
    app.selected = ["EcoRI"].iter().map(|s| s.to_string()).collect();

    app.list_dialog_ok();

    assert_eq!(app.status, "Selected Sites: 1 enzyme.");
}

#[test]
fn a_result_window_opens_as_wide_as_its_longest_line_but_not_wider_than_the_screen() {
    let (ctx, mut app) = (Context::default(), app());
    settle(&mut app, &ctx);
    let screen = ctx.content_rect().size();

    let short = natural_size(&ctx, shape_of("ab\ncd\n"));
    let wide = natural_size(&ctx, shape_of(&"x".repeat(100)));
    let huge = natural_size(&ctx, shape_of(&"x".repeat(5_000)));
    let tall = natural_size(&ctx, shape_of(&"x\n".repeat(5_000)));

    assert_eq!(short, MIN_WINDOW, "a small text opens at the least size");
    assert!(
        wide[0] > short[0] && wide[0] < screen.x - SCREEN_ROOM,
        "{wide:?}"
    );
    assert!(
        (huge[0] - (screen.x - SCREEN_ROOM)).abs() < 0.01,
        "{huge:?} on {screen:?}"
    );
    assert!(
        (tall[1] - (screen.y - SCREEN_ROOM)).abs() < 0.01,
        "{tall:?} on {screen:?}"
    );
}

#[test]
fn the_shape_of_a_text_is_its_lines_and_the_characters_of_the_longest() {
    assert_eq!(shape_of(""), (0, 0));
    assert_eq!(shape_of("ab\ncdef\n"), (2, 4));
    assert_eq!(shape_of("a\n\nbcd"), (3, 3));
}

fn key_press(key: Key) -> Vec<Event> {
    let event = |pressed| Event::Key {
        key,
        physical_key: Some(key),
        pressed,
        repeat: false,
        modifiers: Modifiers::NONE,
    };
    vec![event(true), event(false)]
}

#[test]
fn escape_and_the_cancel_button_close_the_dialog_without_listing() {
    let (ctx, mut app) = (Context::default(), app());
    app.open_list_dialog();
    settle(&mut app, &ctx);
    frame(&mut app, &ctx, key_press(Key::Escape));
    assert!(
        !app.enz_list.dialog_open && app.enz_list.run.is_none(),
        "Escape"
    );

    app.open_list_dialog();
    let shown = settle(&mut app, &ctx);
    click_on(&mut app, &ctx, &shown, "Cancel");
    assert!(
        !app.enz_list.dialog_open && app.enz_list.run.is_none(),
        "Cancel"
    );
}

#[test]
fn return_that_runs_a_palette_command_is_not_also_ok_in_the_dialog() {
    let (ctx, mut app) = (Context::default(), app());
    app.open_list_dialog();
    settle(&mut app, &ctx);
    app.palette.toggle();
    settle(&mut app, &ctx);

    frame(&mut app, &ctx, key_press(Key::Enter));

    assert!(app.enz_list.dialog_open, "the dialog is still there");
    assert!(app.enz_list.run.is_none(), "and listed nothing");
}

#[test]
fn return_in_a_text_field_is_not_for_the_dialog() {
    let (ctx, mut app) = (Context::default(), app());
    app.find_win.open = true;
    settle(&mut app, &ctx);
    ctx.memory_mut(|m| m.request_focus(Id::new("find_text")));
    settle(&mut app, &ctx);
    app.open_list_dialog();
    settle(&mut app, &ctx);
    assert!(
        ctx.text_edit_focused(),
        "Find's text field has the keyboard"
    );

    // Find takes the key and lets go of the field in the same frame; the dialog still must not press OK
    frame(&mut app, &ctx, key_press(Key::Enter));

    assert!(app.enz_list.run.is_none(), "Return did not press OK");
    assert!(app.enz_list.dialog_open, "and the dialog is still there");
}

#[test]
fn the_dialog_does_not_open_without_a_sequence() {
    let mut app = app();
    app.docs.clear();

    app.open_list_dialog();

    assert!(!app.enz_list.dialog_open);
    assert_eq!(app.status, "Open a sequence first.");
}
