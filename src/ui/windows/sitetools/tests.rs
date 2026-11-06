use super::*;
use crate::ui::app::Startup;
use crate::ui::headless::frame_input;
use eframe::egui::{pos2, Context, Event, Modifiers, PointerButton};

/// The coding stretch of the demo plasmid and enzymes that make a few places in it (the `--view silent` screen).
fn app_on_the_coding_stretch() -> App {
    let mut app = App::new(Startup::default());
    app.docs[0].select(326, 446);
    app.selected = [
        "BamHI", "EcoRV", "HindIII", "KpnI", "NheI", "PstI", "SacI", "XhoI",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    app
}

fn frame(app: &mut App, ctx: &Context, n: u32, events: Vec<Event>) {
    let mut input = frame_input(n);
    input.events = events;
    let mut output = ctx.run_ui(input, |ui| app.draw(ui));
    output.textures_delta.clear();
}

/// A click in the middle of the widget `id`, which must have been drawn in the last frame.
fn click(app: &mut App, ctx: &Context, n: u32, id: Id) {
    let rect = ctx.read_response(id).expect("the row was drawn").rect;
    let at = pos2(rect.left() + 20.0, rect.center().y);
    let button = |pressed| Event::PointerButton {
        pos: at,
        button: PointerButton::Primary,
        pressed,
        modifiers: Modifiers::NONE,
    };
    frame(app, ctx, n, vec![Event::PointerMoved(at)]);
    frame(app, ctx, n + 1, vec![button(true), button(false)]);
}

#[test]
fn silent_sites_lists_places_in_the_selected_stretch_and_a_click_selects_the_codons() {
    let (ctx, mut app) = (Context::default(), app_on_the_coding_stretch());
    app.run_silent_sites();
    assert!(app.site_win.silent_open);
    let run = app.site_win.silent.as_ref().expect("a run");
    assert_eq!(
        (run.at.from, run.at.to),
        (326, 446),
        "only the selected stretch"
    );
    assert!(!run.all && run.enzymes == 8);
    let hits = run.result.hits.clone();
    assert!(!hits.is_empty(), "the demo's coding stretch has places");
    assert!(hits
        .iter()
        .all(|h| h.codons_at.0 >= 326 && h.codons_at.1 <= 446));
    // the places come in order of the site's start
    assert!(hits.windows(2).all(|w| w[0].site_start <= w[1].site_start));

    frame(&mut app, &ctx, 0, vec![]);
    frame(&mut app, &ctx, 1, vec![]);
    app.workspace = Workspace::Home;
    click(&mut app, &ctx, 2, silent_row_id(0));
    assert_eq!(app.docs[0].selection(), Some(hits[0].codons_at));
    assert_eq!(
        app.workspace,
        Workspace::Sequences,
        "the click goes to the sequence"
    );
}

#[test]
fn a_search_that_the_sequence_has_changed_since_cannot_be_clicked() {
    let (ctx, mut app) = (Context::default(), app_on_the_coding_stretch());
    app.run_silent_sites();
    frame(&mut app, &ctx, 0, vec![]);
    frame(&mut app, &ctx, 1, vec![]);
    let rect = ctx.read_response(silent_row_id(0)).expect("drawn").rect;
    app.docs[0].select(0, 5);
    app.docs[0].rev += 1;
    frame(&mut app, &ctx, 2, vec![]);
    let at = pos2(rect.left() + 20.0, rect.center().y);
    let button = |pressed| Event::PointerButton {
        pos: at,
        button: PointerButton::Primary,
        pressed,
        modifiers: Modifiers::NONE,
    };
    frame(&mut app, &ctx, 3, vec![Event::PointerMoved(at)]);
    frame(&mut app, &ctx, 4, vec![button(true), button(false)]);
    assert_eq!(
        app.docs[0].selection(),
        Some((0, 5)),
        "the stale row did nothing"
    );
    assert!(!app.is_current(app.site_win.silent.as_ref().map(|r| &r.at)));
}

#[test]
fn with_no_enzyme_chosen_every_enzyme_is_searched_and_the_window_says_so() {
    let (ctx, mut app) = (Context::default(), App::new(Startup::default()));
    app.selected.clear();
    app.run_silent_sites();
    let run = app.site_win.silent.as_ref().expect("a run");
    assert!(run.all);
    assert_eq!(run.enzymes, app.enzymes.len());
    assert_eq!((run.at.from, run.at.to), (0, app.docs[0].doc.len()));
    frame(&mut app, &ctx, 0, vec![]);
    frame(&mut app, &ctx, 1, vec![]);
}

#[test]
fn the_copied_text_is_the_originals_report() {
    let mut app = app_on_the_coding_stretch();
    app.run_silent_sites();
    let text = app.silent_text();
    let name = app.docs[0].doc.name.clone();
    assert!(text.starts_with(&format!(
        "{PROGRAM}\n{name} \nFrom 327 to 446. \nSilent Sites\n\n"
    )));
    let hits = &app.site_win.silent.as_ref().expect("a run").result.hits;
    assert_eq!(
        text.matches("\n\t\t ").count(),
        hits.len(),
        "one residue line for each place"
    );
}

#[test]
fn a_search_without_a_sequence_says_so_and_opens_nothing() {
    let mut app = App::new(Startup::default());
    app.docs.clear();
    app.run_silent_sites();
    app.run_diagnostic_sites();
    assert!(!app.site_win.silent_open && app.site_win.silent.is_none());
    assert!(app.site_win.diagnostic.is_none());
    assert_eq!(app.status, "Open a sequence first.");
}

#[test]
fn diagnostic_sites_search_the_chosen_enzymes_and_a_click_selects_the_bases() {
    let (ctx, mut app) = (Context::default(), app_on_the_coding_stretch());
    app.site_win.diagnostic_open = true;
    app.site_win.mismatches = 2;
    app.run_diagnostic_sites();
    let run = app.site_win.diagnostic.as_ref().expect("a run");
    let found = run.result.found.clone();
    assert!(!found.is_empty());
    assert!(found.iter().all(|f| f.mismatches >= 1 && f.mismatches <= 2));
    assert_eq!(run.dna.len(), 120);

    frame(&mut app, &ctx, 0, vec![]);
    frame(&mut app, &ctx, 1, vec![]);
    app.docs[0].select(0, 3);
    click(&mut app, &ctx, 2, diagnostic_row_id(0));
    let (from, len) = (326 + found[0].index, found[0].pattern.len());
    assert_eq!(app.docs[0].selection(), Some((from, from + len)));
}

#[test]
fn diagnostic_sites_can_search_the_enzymes_that_are_absent() {
    let mut app = App::new(Startup::default());
    app.site_win.absent = true;
    app.site_win.mismatches = 1;
    app.run_diagnostic_sites();
    // the sites the sequence has, with Dam and Dcm methylation blocking some, as the Enzyme Selector counts them
    let counts: Vec<usize> = app.docs[0]
        .sites(&app.enzymes)
        .iter()
        .map(Vec::len)
        .collect();
    let absent = counts.iter().filter(|&&n| n == 0).count();
    let run = app.site_win.diagnostic.as_ref().expect("a run");
    assert_eq!(run.enzymes, absent);
    assert!(absent > 0 && absent < app.enzymes.len());
    for f in &run.result.found {
        let k = app
            .enzymes
            .iter()
            .position(|e| e.name == f.enzyme)
            .expect("an enzyme");
        assert_eq!(counts[k], 0, "{} has a site", f.enzyme);
    }
}

#[test]
fn the_sites_that_are_already_there_follow_the_methylation_switch_of_the_sequence() {
    let mut app = App::new(Startup::default());
    app.site_win.absent = true;
    for methylated in [true, false] {
        app.docs[0].doc.methylated = methylated;
        let absent = app.docs[0]
            .sites(&app.enzymes)
            .iter()
            .filter(|sites| sites.is_empty())
            .count();
        app.run_diagnostic_sites();
        let run = app.site_win.diagnostic.as_ref().expect("a run");
        assert_eq!(run.enzymes, absent, "methylated: {methylated}");
    }
}

#[test]
fn the_chosen_enzymes_keep_the_order_of_the_enzyme_list() {
    let mut app = App::new(Startup::default());
    app.selected = ["XhoI", "BamHI", "EcoRI"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let names: Vec<&str> = app
        .chosen_enzymes()
        .iter()
        .map(|&k| app.enzymes[k].name.as_str())
        .collect();
    assert_eq!(names, vec!["BamHI", "EcoRI", "XhoI"]);
}

#[test]
fn both_tools_are_in_the_palette_and_take_the_person_to_the_sequence() {
    use crate::ui::app::Action;
    use crate::ui::commands::{acts_on_the_sequence, all};
    let answers: Vec<(&str, bool)> = all()
        .iter()
        .map(|c| {
            (
                c.title,
                matches!(c.action, Action::SilentSites | Action::DiagnosticSites),
            )
        })
        .filter(|&(_, ours)| ours)
        .collect();
    let titles: Vec<&str> = answers.iter().map(|&(t, _)| t).collect();
    assert_eq!(
        titles,
        vec!["Silent restriction sites…", "Diagnostic restriction sites…"]
    );
    assert!(acts_on_the_sequence(&Action::SilentSites));
    assert!(acts_on_the_sequence(&Action::DiagnosticSites));
}

/// A long sequence with every enzyme gives far more places than the window draws; a frame must still be quick.
#[test]
fn a_list_longer_than_the_window_draws_stays_quick() {
    use crate::seq::doc::Doc;
    use crate::ui::state::DocView;
    use std::time::Instant;
    // a repeatable mixed sequence of 30,000 bases
    let mut state = 12_345u32;
    let seq: Vec<u8> = (0..30_000)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            b"ACGT"[(state >> 16) as usize % 4]
        })
        .collect();
    let (ctx, mut app) = (Context::default(), App::new(Startup::default()));
    app.docs[0] = DocView::new(Doc::new("Long", seq));
    app.selected.clear();
    app.run_silent_sites();
    let total = app
        .site_win
        .silent
        .as_ref()
        .expect("a run")
        .result
        .hits
        .len();
    assert!(total > ROWS_DRAWN, "{total} places");
    frame(&mut app, &ctx, 0, vec![]);
    let started = Instant::now();
    for n in 1..6 {
        frame(&mut app, &ctx, n, vec![]);
    }
    let each = started.elapsed() / 5;
    eprintln!("{total} places; a frame with {ROWS_DRAWN} rows takes {each:?}");
    assert!(each.as_millis() < 500, "a frame took {each:?}");
}

#[test]
fn searching_again_after_a_click_on_a_row_keeps_the_stretch_but_a_new_selection_changes_it() {
    let (ctx, mut app) = (Context::default(), app_on_the_coding_stretch());
    app.run_silent_sites();
    frame(&mut app, &ctx, 0, vec![]);
    frame(&mut app, &ctx, 1, vec![]);
    click(&mut app, &ctx, 2, silent_row_id(0));
    assert_ne!(
        app.docs[0].selection(),
        Some((326, 446)),
        "the click selected the codons"
    );
    app.run_silent_sites();
    let at = &app.site_win.silent.as_ref().expect("a run").at;
    assert_eq!(
        (at.from, at.to),
        (326, 446),
        "the stretch the list was made from"
    );
    // a stretch the person selects is the one searched
    app.docs[0].select(330, 400);
    app.run_silent_sites();
    let at = &app.site_win.silent.as_ref().expect("a run").at;
    assert_eq!((at.from, at.to), (330, 400));
}

#[test]
fn finding_again_after_a_click_on_a_row_keeps_the_stretch_too() {
    let (ctx, mut app) = (Context::default(), app_on_the_coding_stretch());
    app.site_win.diagnostic_open = true;
    app.run_diagnostic_sites();
    frame(&mut app, &ctx, 0, vec![]);
    frame(&mut app, &ctx, 1, vec![]);
    click(&mut app, &ctx, 2, diagnostic_row_id(0));
    app.site_win.mismatches = 3;
    app.run_diagnostic_sites();
    let run = app.site_win.diagnostic.as_ref().expect("a run");
    assert_eq!((run.at.from, run.at.to), (326, 446));
    assert_eq!(run.mismatches, 3, "the run says what it was made with");
    assert!(!run.absent);
}

#[test]
fn a_click_anywhere_in_a_row_selects_it() {
    let (ctx, mut app) = (Context::default(), app_on_the_coding_stretch());
    app.run_silent_sites();
    frame(&mut app, &ctx, 0, vec![]);
    frame(&mut app, &ctx, 1, vec![]);
    let hit = app.site_win.silent.as_ref().expect("a run").result.hits[0].clone();
    let rect = ctx.read_response(silent_row_id(0)).expect("drawn").rect;
    let at = pos2(rect.left() + 30.0, rect.bottom() - 3.0);
    let press = |pressed| Event::PointerButton {
        pos: at,
        button: PointerButton::Primary,
        pressed,
        modifiers: Modifiers::NONE,
    };
    frame(&mut app, &ctx, 2, vec![Event::PointerMoved(at)]);
    frame(&mut app, &ctx, 3, vec![press(true), press(false)]);
    assert_eq!(
        app.docs[0].selection(),
        Some(hit.codons_at),
        "the last line of the row"
    );
}

#[test]
fn two_sequences_with_the_same_name_are_told_apart() {
    use crate::seq::doc::Doc;
    let mut app = app_on_the_coding_stretch();
    app.run_silent_sites();
    assert!(app.is_current(app.site_win.silent.as_ref().map(|r| &r.at)));
    let twin = Doc::new(app.docs[0].doc.name.clone(), app.docs[0].doc.seq.clone());
    app.add_doc(twin);
    assert_ne!(app.docs[0].id, app.docs[1].id);
    assert!(
        !app.is_current(app.site_win.silent.as_ref().map(|r| &r.at)),
        "the twin in front is not the sequence the list was made from"
    );
}

#[test]
fn a_search_with_no_enzymes_says_so_and_draws() {
    let (ctx, mut app) = (Context::default(), app_on_the_coding_stretch());
    app.selected.clear();
    app.site_win.diagnostic_open = true;
    app.run_diagnostic_sites();
    assert_eq!(app.site_win.diagnostic.as_ref().expect("a run").enzymes, 0);
    frame(&mut app, &ctx, 0, vec![]);
    frame(&mut app, &ctx, 1, vec![]);
}
