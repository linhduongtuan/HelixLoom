use super::*;
use crate::seq::astral::BranchOut;
use crate::seq::treetools::robinson_foulds;
use crate::ui::app::Startup;
use crate::ui::headless::frame_input;
use eframe::egui::{
    epaint::Shape, pos2, Context, Event, FullOutput, Modifiers, PointerButton, Rect,
};
use std::cell::Cell;
use std::time::Instant;

thread_local! {
    /// The frames this test has drawn: the time of a frame must never go back, whatever helper draws it.
    static FRAMES: Cell<u32> = const { Cell::new(0) };
}

fn frame(app: &mut App, ctx: &Context, events: Vec<Event>) -> FullOutput {
    let n = FRAMES.with(|f| f.replace(f.get() + 1));
    let mut input = frame_input(n);
    input.events = events;
    let mut output = ctx.run_ui(input, |ui| app.draw(ui));
    output.textures_delta.clear();
    output
}

/// Draw frames until `done` says that what was waited for has happened.
fn draw_until(app: &mut App, ctx: &Context, mut done: impl FnMut(&mut App) -> bool) {
    let end = Instant::now() + Duration::from_secs(60);
    loop {
        frame(app, ctx, vec![]);
        if done(app) {
            return;
        }
        assert!(Instant::now() < end, "gave up waiting");
        std::thread::sleep(Duration::from_millis(5));
    }
}

/// Where `text` was drawn in a frame: the first piece of text that holds it.
fn text_rect(shapes: &[Shape], text: &str) -> Option<Rect> {
    shapes.iter().find_map(|shape| match shape {
        Shape::Text(t) if t.galley.text().contains(text) => Some(t.visual_bounding_rect()),
        Shape::Vec(inner) => text_rect(inner, text),
        _ => None,
    })
}

/// Click the middle of the piece of text `text` that a frame drew: a real click on the button that holds it.
fn click_on(app: &mut App, ctx: &Context, shown: &FullOutput, text: &str) {
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
    frame(app, ctx, vec![button(true), button(false)]);
}

fn open_window() -> (App, Context) {
    let mut app = App::new(Startup::default());
    app.species_win.open = true;
    (app, Context::default())
}

/// The first tree of `text`, as the tools that compare trees read it.
fn tree(text: &str) -> NwTree {
    NwTree::parse_all(text).expect("a tree").remove(0)
}

fn distance_from_truth(newick: &str) -> usize {
    robinson_foulds(&tree(newick), &tree(SPECIES_TRUTH))
        .expect("the same tips")
        .distance
}

#[test]
fn the_example_finds_the_species_tree_it_was_simulated_from() {
    let (mut app, ctx) = open_window();
    app.load_species_example();
    assert!(app.species_win.result.is_none());
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let r = app.species_win.result.clone().expect("a result");
    assert_eq!((r.genes, r.taxa.len(), r.branches.len()), (40, 8, 5));
    assert!(r.clusters.is_some());
    assert_eq!(distance_from_truth(&r.newick), 0, "{}", r.newick);
    let said = &app.species_win.message;
    assert!(
        said.contains("40 gene trees on 8 species") && said.contains("on this computer"),
        "{said}"
    );
    assert!(!app.species_win.is_running() && app.running_analyses() == 0);
}

#[test]
fn pressing_the_button_runs_the_search_and_pressing_open_in_trees_draws_the_tree() {
    let (mut app, ctx) = open_window();
    app.load_species_example();
    frame(&mut app, &ctx, vec![]);
    let shown = frame(&mut app, &ctx, vec![]);
    click_on(&mut app, &ctx, &shown, "Find the species tree");
    let w = &app.species_win;
    assert!(
        w.run_requested || w.is_running() || w.result.is_some(),
        "the click did not start a run"
    );
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());

    let shown = frame(&mut app, &ctx, vec![]);
    click_on(&mut app, &ctx, &shown, "Open in Trees");
    assert_eq!(app.workspace, Workspace::Trees);
    let drawn = app
        .tree_win
        .tree
        .as_ref()
        .expect("a tree on the Trees page");
    assert_eq!(drawn.leaf_count(), 8);
    let said = &app.platform.trees.message;
    assert!(said.contains("40 gene trees"), "{said}");
}

#[test]
fn a_given_tree_is_scored_not_searched_and_has_the_score_the_search_found_for_that_tree() {
    let (mut app, ctx) = open_window();
    app.load_species_example();
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let searched = app.species_win.result.clone().expect("a result");

    app.species_win.score_given = true;
    app.run_species();
    draw_until(&mut app, &ctx, |a| {
        a.species_win
            .result
            .as_ref()
            .is_some_and(|r| r.clusters.is_none())
    });
    let scored = app.species_win.result.clone().expect("a result");
    // the example's own species tree is the tree the search found
    assert_eq!(
        (scored.score, scored.branches.len()),
        (searched.score, searched.branches.len())
    );
    let said = &app.species_win.message;
    assert!(said.contains("Scored the tree you gave"), "{said}");

    // a tree that does not have the species of the genes is refused, with what is wrong
    app.species_win.given = "((A,B),(C,Z));".into();
    app.run_species();
    draw_until(&mut app, &ctx, |a| {
        !a.species_win.is_running() && !a.species_win.run_requested
    });
    assert!(app.species_win.result.is_none());
    let said = &app.species_win.message;
    assert!(!said.is_empty() && !said.starts_with("Scored"), "{said}");
}

/// Press Find the species tree with these gene trees and settle: what the window says, and whether it shows a result.
fn run_with(app: &mut App, ctx: &Context, genes: &str) -> (String, bool) {
    app.species_win.genes = genes.into();
    app.run_species();
    draw_until(app, ctx, |a| {
        !a.species_win.run_requested && !a.species_win.is_running()
    });
    (
        app.species_win.message.clone(),
        app.species_win.result.is_some(),
    )
}

#[test]
fn nothing_to_run_and_trees_that_cannot_be_used_are_said_and_nothing_is_shown() {
    let (mut app, ctx) = open_window();
    let (said, shown) = run_with(&mut app, &ctx, "");
    assert!(said.contains("Paste the gene trees"), "{said}");
    assert!(!shown && !app.species_win.is_running());

    // a tree to score that was not given is said, not quietly replaced by a search
    app.species_win.score_given = true;
    app.species_win.given.clear();
    let (said, shown) = run_with(&mut app, &ctx, SPECIES_GENE_TREES);
    assert!(!shown && said.contains("species tree to score"), "{said}");
    app.species_win.score_given = false;

    // three species are too few for a quartet
    let (said, shown) = run_with(&mut app, &ctx, "((A,B),C);\n((A,C),B);");
    assert!(
        !shown && !said.is_empty() && !said.contains("Working"),
        "{said}"
    );

    // a species twice in a gene is a gene family: the other program of the family is for that
    let (said, shown) = run_with(&mut app, &ctx, "((A,B),(A,C),D);");
    assert!(!shown && said.contains("ASTRAL-Pro"), "{said}");

    // text that is no tree at all
    let (said, shown) = run_with(&mut app, &ctx, "this is not a tree");
    assert!(
        !shown && !said.is_empty() && !said.contains("Working"),
        "{said}"
    );

    // the exact search takes few species: thirteen are too many
    app.species_win.exact = true;
    let thirteen = "((((((((((((a,b),c),d),e),f),g),h),i),j),k),l),m);\n".repeat(3);
    let (said, shown) = run_with(&mut app, &ctx, &thirteen);
    assert!(!shown && said.contains("at most 12"), "{said}");
}

#[test]
fn a_file_of_gene_trees_stays_out_of_the_box_and_is_the_text_that_runs() {
    let (mut app, ctx) = open_window();
    app.species_win.genes = "typed text that the file replaces".into();
    let text = SPECIES_GENE_TREES.repeat(3);
    app.load_species_file(Source::Bytes {
        name: "genes.tre",
        bytes: text.as_bytes(),
    });
    let w = &app.species_win;
    assert!(w.genes.is_empty(), "the box is not given the file");
    assert_eq!(
        w.loaded.as_ref().map(|l| (l.name.as_str(), l.trees)),
        Some(("genes.tre", 120))
    );
    assert_eq!(w.weight(), Weight::of_species(&text));
    let said = &w.message;
    assert!(
        said.contains("genes.tre") && said.contains("120 gene trees"),
        "{said}"
    );
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    assert_eq!(app.species_win.result.as_ref().map(|r| r.genes), Some(120));

    // a file that cannot be read says so and keeps what was there
    app.load_species_file(Source::Bytes {
        name: "bad.tre",
        bytes: &[0xff, 0xfe, 0x00],
    });
    assert!(app.species_win.loaded.is_some());
    assert!(!app.species_win.message.is_empty());
    // Clear empties the box, the file, the result and the message
    app.species_win.clear();
    let w = &app.species_win;
    assert!(w.loaded.is_none() && w.genes.is_empty() && w.result.is_none() && w.message.is_empty());
}

/// A search that takes a little while, as one over many genes does, for two gene trees (so that its answer is told from the
/// example's forty).
fn slow_search(ctx: &Context) -> Job<Result<SpeciesResult, String>> {
    Job::spawn(ctx, || {
        std::thread::sleep(Duration::from_millis(200));
        species::run(&SpeciesJob::search(
            "((A,B),(C,(D,E)));\n((A,C),(B,(D,E)));",
        ))
    })
}

#[test]
fn clearing_or_loading_other_trees_lets_a_run_that_was_going_end_but_never_shows_its_answer() {
    let (mut app, ctx) = open_window();
    app.species_win.job = Some(slow_search(&ctx));
    app.species_win.genes = "to be cleared".into();
    app.species_win.clear();
    // a search cannot be stopped: it still counts as running (it is using the cores) and a new run waits for it
    assert!(app.species_win.is_running() && app.running_analyses() == 1);
    app.load_species_example();
    app.run_species();
    frame(&mut app, &ctx, vec![]);
    assert!(app.species_win.result.is_none(), "the new run waits");
    // when the old search ends its answer is let go, and the new run gives its own
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let r = app.species_win.result.clone().expect("a result");
    assert_eq!(r.genes, 40, "the example's answer, not the old search's");
    let said = &app.species_win.message;
    assert!(said.contains("40 gene trees on 8 species"), "{said}");
    assert!(!app.species_win.discard && app.running_analyses() == 0);

    // loading a file while a search is going does the same, and leaves the file's message
    app.species_win.job = Some(slow_search(&ctx));
    app.load_species_file(Source::Bytes {
        name: "other.tre",
        bytes: SPECIES_GENE_TREES.as_bytes(),
    });
    assert!(app.species_win.discard);
    draw_until(&mut app, &ctx, |a| !a.species_win.is_running());
    assert!(app.species_win.result.is_none());
    let said = &app.species_win.message;
    assert!(said.contains("other.tre"), "{said}");
}

#[test]
fn the_view_that_pictures_the_window_runs_the_example() {
    let mut app = App::new(Startup {
        view: "species".into(),
        ..Startup::default()
    });
    assert!(app.species_win.open && app.species_win.run_requested);
    let ctx = Context::default();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
}

#[test]
fn the_window_draws_in_every_state_in_both_themes() {
    let ctx = Context::default();
    for dark in [false, true] {
        let mut app = App::new(Startup {
            dark,
            ..Startup::default()
        });
        app.species_win.open = true;
        frame(&mut app, &ctx, vec![]); // nothing yet
        app.load_species_example();
        frame(&mut app, &ctx, vec![]); // the example in the box
        app.run_species();
        draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
        frame(&mut app, &ctx, vec![]); // a result
        app.species_win.score_given = true;
        frame(&mut app, &ctx, vec![]); // a tree to score
        app.species_win.score_given = false;
        app.load_species_file(Source::Bytes {
            name: "g.tre",
            bytes: SPECIES_GENE_TREES.as_bytes(),
        });
        frame(&mut app, &ctx, vec![]); // a file
        app.species_win.message = "x".repeat(500);
        frame(&mut app, &ctx, vec![]); // a long message
        for mode in crate::remote::settings::Compute::ALL {
            app.platform.workstation.compute = mode;
            frame(&mut app, &ctx, vec![]);
        }
    }
}

fn branch(taxa: usize) -> BranchOut {
    BranchOut {
        taxa: (0..taxa).map(|i| format!("sp{i}")).collect(),
        f: [0.5, 0.3, 0.2],
        q: [0.5, 0.3, 0.2],
        pp: [0.97, 0.02, 0.01],
        en: 10.0,
        quartets: 100,
        length: 1.25,
    }
}

#[test]
fn long_clades_are_abbreviated_and_a_table_that_is_cut_says_so() {
    assert_eq!(side(&branch(2).taxa), "sp0, sp1");
    assert_eq!(side(&branch(NAMED).taxa), "sp0, sp1, sp2, sp3");
    assert_eq!(
        side(&branch(NAMED + 3).taxa),
        "sp0, sp1, sp2, sp3, … (3 more)"
    );
    assert_eq!(length_text(1.25), "1.250");
    assert_eq!(length_text(f64::INFINITY), "∞");
    assert_eq!(length_text(f64::NAN), "–");
    assert_eq!(cut_note(MAX_LISTED), None);
    let said = cut_note(MAX_LISTED + 50).expect("a note");
    assert!(
        said.contains(&format!("first {MAX_LISTED} of {}", MAX_LISTED + 50)),
        "{said}"
    );
}

#[test]
fn a_result_with_many_branches_is_summed_up_and_drawn_whatever_its_size() {
    let result = SpeciesResult {
        newick: "(sp0,sp1,(sp2,sp3));".into(),
        score: 1,
        normalized: 0.5,
        taxa: (0..8).map(|i| format!("sp{i}")).collect(),
        genes: 3,
        skipped: 2,
        clusters: Some(9),
        branches: (0..MAX_LISTED + 50).map(|i| branch(2 + i % 6)).collect(),
        branch_count: MAX_LISTED + 50,
        seconds: 0.5,
    };
    let said = summary(&result, "the workstation");
    for part in [
        "3 gene trees on 8 species",
        "2 with fewer than 4",
        "50.0 %",
        "9 groups",
        "on the workstation",
    ] {
        assert!(said.contains(part), "{part} in {said}");
    }
    let (mut app, ctx) = open_window();
    app.show_species_result(result);
    assert!(app.species_win.open);
    assert_eq!(app.species_win.message, said);
    frame(&mut app, &ctx, vec![]);
    frame(&mut app, &ctx, vec![]);
}

#[test]
fn a_tree_with_infinite_lengths_is_drawn_without_them_and_one_that_cannot_be_read_says_so() {
    // a prior of 0 gives an infinite length where every gene agrees and none where no gene covers a branch
    let (mut app, ctx) = open_window();
    let mut result = species::run(&SpeciesJob::search(SPECIES_GENE_TREES)).expect("a species tree");
    result.newick = "(A,((C,E)0.95:Infinity,(G,((B,F)1:NaN,(D,H)1:2.5)0.86:0.2)1:0.5):0.1);".into();
    app.show_species_result(result);
    app.show_species_in_trees();
    assert_eq!(app.workspace, Workspace::Trees);
    assert_eq!(app.tree_win.tree.as_ref().map(|t| t.leaf_count()), Some(8));

    // text that is no tree is said, and nothing is drawn over the tree that was there
    let mut broken = species::run(&SpeciesJob::search(SPECIES_GENE_TREES)).expect("a species tree");
    broken.newick = "((a,b".into();
    app.workspace = Workspace::Sequences;
    app.show_species_result(broken);
    app.show_species_in_trees();
    assert_eq!(app.workspace, Workspace::Sequences);
    let said = &app.species_win.message;
    assert!(said.contains("cannot be drawn"), "{said}");
    frame(&mut app, &ctx, vec![]);
}
