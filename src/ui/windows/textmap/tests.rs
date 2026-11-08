use super::view::base_at;
use super::*;
use crate::seq::textmap::{Frame, MARGIN};
use crate::ui::app::{Action, Startup};
use crate::ui::headless::frame_input;
use eframe::egui::{
    epaint::Shape, pos2, text::CCursor, Context, Event, FullOutput, Key, Modifiers, PointerButton,
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
    frame_holding(app, ctx, events, Modifiers::NONE)
}

/// A frame with these modifier keys held down.
fn frame_holding(app: &mut App, ctx: &Context, events: Vec<Event>, held: Modifiers) -> FullOutput {
    let n = FRAMES.with(|f| f.replace(f.get() + 1));
    let mut input = frame_input(n);
    // which modifier keys are held is an event of its own, and holds until the next one says otherwise
    input.events = std::iter::once(Event::ModifiersChanged(held))
        .chain(events)
        .collect();
    let mut output = ctx.run_ui(input, |ui| app.draw(ui));
    output.textures_delta.clear();
    output
}

/// Frames enough for a window that has just opened to be laid out (a window with a scroll area is a frame or two slower to find
/// its size than a small one) and clicked; the last one.
fn settle(app: &mut App, ctx: &Context) -> FullOutput {
    (0..5).fold(frame(app, ctx, vec![]), |_, _| frame(app, ctx, vec![]))
}

fn shapes_of(shown: &FullOutput) -> Vec<Shape> {
    shown.shapes.iter().map(|s| s.shape.clone()).collect()
}

/// Where `text` was drawn in a frame: the first piece of text that holds it.
fn text_rect(shapes: &[Shape], text: &str) -> Option<Rect> {
    shapes.iter().find_map(|shape| match shape {
        Shape::Text(t) if t.galley.text().contains(text) => Some(t.visual_bounding_rect()),
        Shape::Vec(inner) => text_rect(inner, text),
        _ => None,
    })
}

/// A real click at a place of the window.
fn click_at(app: &mut App, ctx: &Context, at: eframe::egui::Pos2, modifiers: Modifiers) {
    let button = |pressed| Event::PointerButton {
        pos: at,
        button: PointerButton::Primary,
        pressed,
        modifiers,
    };
    frame_holding(app, ctx, vec![Event::PointerMoved(at)], modifiers);
    frame_holding(app, ctx, vec![button(true), button(false)], modifiers);
}

/// A real click on the middle of the piece of text `text` that a frame drew.
fn click_on(app: &mut App, ctx: &Context, shown: &FullOutput, text: &str) {
    let rect = text_rect(&shapes_of(shown), text).unwrap_or_else(|| panic!("{text} was not drawn"));
    click_at(
        app,
        ctx,
        pos2(rect.center().x, rect.center().y),
        Modifiers::NONE,
    );
}

fn run_of(app: &App) -> &MapRun {
    app.text_map.run.as_ref().expect("a page was made")
}

fn key_press(key: Key, modifiers: Modifiers) -> Vec<Event> {
    let event = |pressed| Event::Key {
        key,
        physical_key: Some(key),
        pressed,
        repeat: false,
        modifiers,
    };
    vec![event(true), event(false)]
}

#[test]
fn a_box_moves_one_place_and_stays_where_it_is_at_the_ends() {
    let order = Frame::ALL.to_vec();
    let moved = move_frame(&order, Frame::Dna, true);
    assert_eq!(moved[1], Frame::Dna);
    assert_eq!(moved[2], Frame::Enzymes);
    assert_eq!(move_frame(&order, Frame::Index, true), order);
    assert_eq!(move_frame(&order, Frame::Space, false), order);
    // the order that was given is not changed, and nothing is lost
    assert_eq!(order, Frame::ALL.to_vec());
    let mut sorted = moved.clone();
    sorted.sort_by_key(|f| f.label());
    let mut all = Frame::ALL.to_vec();
    all.sort_by_key(|f| f.label());
    assert_eq!(sorted, all);
}

#[test]
fn the_page_of_the_demo_plasmid_has_its_date_name_and_rows() {
    let mut app = app();
    app.make_text_map();
    let run = run_of(&app);
    let (name, len) = (app.docs[0].doc.name.clone(), app.docs[0].doc.len());
    let lines: Vec<&str> = run.page.text.lines().collect();
    // the date and the time (UTC), then the name of the sequence
    assert!(
        lines[0].ends_with("UTC") && lines[0].contains(", 20"),
        "{}",
        lines[0]
    );
    assert_eq!(lines[1], name);
    // a line of DNA for each hundred bases (the start of the dialog), numbered on both sides
    assert_eq!(run.page.dna_rows.len(), len.div_ceil(100));
    assert_eq!(run.page.dna_rows[0].number, 1);
    assert_eq!(run.title, format!("{name} Text Map"));
    assert!(app.text_map.result_open);
}

#[test]
fn the_dialog_opens_with_the_start_of_the_original_and_its_ok_makes_the_page() {
    let (ctx, mut app) = (Context::default(), app());
    assert_eq!(app.text_map.settings, Settings::default());
    app.open_text_map_dialog();
    let shown = settle(&mut app, &ctx);
    for label in [
        "Index Line",
        "Enzymes",
        "DNA",
        "Translation",
        "2nd Strand",
        "Features",
    ] {
        assert!(text_rect(&shapes_of(&shown), label).is_some(), "{label}");
    }
    click_on(&mut app, &ctx, &shown, "OK");
    assert!(app.text_map.run.is_some(), "OK made the page");
    assert!(!app.text_map.dialog_open);
}

#[test]
fn cancel_and_escape_close_the_dialog_without_making_a_page() {
    let (ctx, mut app) = (Context::default(), app());
    app.open_text_map_dialog();
    let shown = settle(&mut app, &ctx);
    click_on(&mut app, &ctx, &shown, "Cancel");
    assert!(!app.text_map.dialog_open && app.text_map.run.is_none());

    app.open_text_map_dialog();
    settle(&mut app, &ctx);
    frame(&mut app, &ctx, key_press(Key::Escape, Modifiers::NONE));
    assert!(!app.text_map.dialog_open && app.text_map.run.is_none());
}

#[test]
fn the_shortcut_and_the_palette_open_the_dialog() {
    let (ctx, mut app) = (Context::default(), app());
    settle(&mut app, &ctx);
    let shift_command = Modifiers::COMMAND | Modifiers::SHIFT;
    frame(&mut app, &ctx, key_press(Key::T, shift_command));
    assert!(app.text_map.dialog_open);
    let known = crate::ui::commands::all()
        .iter()
        .any(|c| matches!(c.action, Action::TextMap));
    assert!(known, "the command palette lists Text Map");
}

#[test]
fn a_dialog_without_a_sequence_says_so_and_does_not_open() {
    let mut app = app();
    app.docs.clear();
    app.open_text_map_dialog();
    assert!(!app.text_map.dialog_open);
    assert_eq!(app.status, "Open a sequence first.");
}

#[test]
fn the_dna_of_the_selection_only_when_the_dialog_says_so() {
    let mut app = app();
    app.docs[0].select(200, 260);
    // "All" is the choice of a start with the Selection Only switch off
    app.open_text_map_dialog();
    assert!(!app.text_map.selection_only);
    app.make_text_map();
    assert_eq!(run_of(&app).page.dna_rows[0].number, 1);

    app.text_map.selection_only = true;
    app.make_text_map();
    let rows = &run_of(&app).page.dna_rows;
    assert_eq!((rows[0].number, rows.len()), (201, 1));
    assert_eq!(rows[0].bases, 60);

    // the switch of the Enzymes menu decides what the dialog starts from, as in the original
    app.prefs.selection_only = true;
    app.open_text_map_dialog();
    assert!(app.text_map.selection_only);
}

#[test]
fn the_features_are_coloured_on_the_lines_of_dna_and_under_them() {
    let mut app = app();
    app.make_text_map();
    let run = run_of(&app);
    let coloured: Vec<usize> = run
        .painted
        .pieces_in_rows()
        .iter()
        .enumerate()
        .filter(|(_, pieces)| !pieces.is_empty())
        .map(|(n, _)| n + 1)
        .collect();
    let dna_lines: Vec<usize> = run.page.dna_rows.iter().map(|r| r.line).collect();
    assert!(
        dna_lines.iter().any(|l| coloured.contains(l)),
        "a line of DNA has a colour behind some of its bases"
    );
    assert!(
        coloured.iter().any(|l| !dna_lines.contains(l)),
        "so do the lines of the features"
    );
}

#[test]
fn a_click_on_a_line_of_dna_puts_the_caret_on_that_base_and_a_shift_click_selects_up_to_it() {
    let (ctx, mut app) = (Context::default(), app());
    app.make_text_map();
    let shown = settle(&mut app, &ctx);
    // the first line of DNA: the 10th base, then the 25th with Shift
    let first = &run_of(&app).page.dna_rows[0];
    let line = app.docs[0].doc.seq[..first.bases]
        .iter()
        .map(|&b| char::from(b))
        .collect::<String>();
    let galley = shapes_of(&shown)
        .iter()
        .find_map(|s| match s {
            Shape::Text(t) if t.galley.text().contains(&line) => Some((t.clone(), t.pos)),
            _ => None,
        })
        .expect("the line of DNA was drawn");
    let at = |col: usize| {
        let rect = galley.0.galley.pos_from_cursor(CCursor::new(col));
        pos2(
            galley.1.x + rect.center().x + 1.0,
            galley.1.y + rect.center().y,
        )
    };
    click_at(&mut app, &ctx, at(MARGIN + 9), Modifiers::NONE);
    assert_eq!((app.docs[0].caret, app.docs[0].anchor), (9, None));
    click_at(&mut app, &ctx, at(MARGIN + 24), Modifiers::SHIFT);
    assert_eq!(app.docs[0].selection(), Some((9, 24)));
}

#[test]
fn only_the_lines_of_dna_answer_a_click() {
    let mut app = app();
    app.make_text_map();
    let run = run_of(&app);
    let row = &run.page.dna_rows[1];
    assert_eq!(base_at(run, row.line, MARGIN), Some(row.number - 1));
    assert_eq!(
        base_at(run, row.line, MARGIN + row.bases - 1),
        Some(row.number - 2 + row.bases)
    );
    // the margin, past the end of the line, the title and the blank line under the page
    assert_eq!(base_at(run, row.line, MARGIN - 1), None);
    assert_eq!(base_at(run, row.line, MARGIN + row.bases), None);
    assert_eq!(base_at(run, 2, MARGIN), None);
}

#[test]
fn a_page_says_when_the_sequence_has_changed_since_it_was_made() {
    let (ctx, mut app) = (Context::default(), app());
    app.make_text_map();
    let shown = settle(&mut app, &ctx);
    let note = "The sequence has changed since this page was made";
    assert!(text_rect(&shapes_of(&shown), note).is_none());
    app.docs[0].edit(|d| d.seq.push(b'A'));
    let shown = settle(&mut app, &ctx);
    assert!(text_rect(&shapes_of(&shown), note).is_some());
    // making it again makes it new
    click_on(&mut app, &ctx, &shown, "Make again");
    let shown = settle(&mut app, &ctx);
    assert!(text_rect(&shapes_of(&shown), note).is_none());
}

#[test]
fn what_the_dialog_holds_is_what_the_page_is_made_with() {
    let mut app = app();
    app.text_map.settings.width = 40;
    app.text_map.settings.numbers = crate::seq::textmap::LineNumbers::None;
    app.text_map.settings.show_index = false;
    app.text_map.settings.show_features = false;
    app.make_text_map();
    let run = run_of(&app);
    let len = app.docs[0].doc.len();
    assert_eq!(run.page.dna_rows.len(), len.div_ceil(40));
    let dna_line = run
        .page
        .text
        .lines()
        .nth(run.page.dna_rows[0].line - 1)
        .expect("a line of DNA");
    assert!(
        dna_line.starts_with("       ") && dna_line.ends_with(" "),
        "{dna_line:?}"
    );
}

#[test]
fn the_clock_line_is_a_date_and_a_time() {
    let line = clock_line();
    let words: Vec<&str> = line.split(' ').collect();
    assert!(
        ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"].contains(&words[0]),
        "{line}"
    );
    assert!(line.ends_with(" UTC") && line.contains(':'), "{line}");
}

#[test]
fn the_second_strand_above_the_dna_is_the_error_the_original_stops_with_and_says_what_to_do() {
    let mut app = app();
    app.text_map.settings.second_strand = true;
    app.text_map.settings.order = move_frame(
        &move_frame(&Frame::ALL, Frame::SecondStrand, true),
        Frame::SecondStrand,
        true,
    );
    // index, enzymes, dna, translation, second strand: two places up puts it above the DNA
    let order = &app.text_map.settings.order;
    assert!(
        order.iter().position(|&f| f == Frame::SecondStrand)
            < order.iter().position(|&f| f == Frame::Dna)
    );
    app.make_text_map();
    assert!(app.text_map.run.is_none(), "no page was made");
    assert!(
        app.status
            .starts_with("Text Map: can't read \"textline\": no such variable"),
        "{}",
        app.status
    );
    assert!(app
        .status
        .contains("put the 2nd Strand box below the DNA box"));
    // below the DNA it is made
    app.text_map.settings.order = Frame::ALL.to_vec();
    app.make_text_map();
    assert!(app.text_map.run.is_some());
}
