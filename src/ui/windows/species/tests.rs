use super::*;
use crate::seq::astral::BranchOut;
use crate::seq::treetools::robinson_foulds;
use crate::ui::app::Startup;
use crate::ui::headless::frame_input;
use crate::ui::provenance as made_ui;
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

/// Whether a frame drew a piece of text that is exactly `text` (a column header, not a sentence that has it).
fn drew_exactly(shown: &FullOutput, text: &str) -> bool {
    fn in_shapes(shapes: &[Shape], text: &str) -> bool {
        shapes.iter().any(|shape| match shape {
            Shape::Text(t) => t.galley.text() == text,
            Shape::Vec(inner) => in_shapes(inner, text),
            _ => false,
        })
    }
    let shapes: Vec<Shape> = shown.shapes.iter().map(|s| s.shape.clone()).collect();
    in_shapes(&shapes, text)
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

    // a species twice in a gene is a gene family: the message says which option reads those
    let (said, shown) = run_with(&mut app, &ctx, "((A,B),(A,C),D);");
    assert!(
        !shown && said.contains("ASTRAL-Pro") && said.contains("gene families option"),
        "{said}"
    );

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
        bootstrap: None,
        polytomy: None,
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
        weights: crate::seq::astral::Weights::None,
        weighted_score: None,
        families: false,
        bootstrap: None,
        polytomy: false,
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

// --- weights ----------------------------------------------------------------------------------------------------------------

/// Three gene trees sure that A goes with B, and four that join A with C and are hardly sure of it.
fn disagreeing() -> String {
    let ab = "(A:0.1,B:0.1,(C:0.1,D:0.1)0.95:0.1);\n".repeat(3);
    let ac = "(A:0.1,C:0.1,(B:0.1,D:0.1)0.05:0.1);\n".repeat(4);
    format!("{ab}{ac}")
}

#[test]
fn the_weights_are_chosen_in_the_settings_and_a_weighted_run_says_what_it_weighed() {
    let (mut app, ctx) = open_window();
    app.species_win.genes = disagreeing();
    // counted one each, the four win
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let plain = app.species_win.result.clone().expect("a result");
    assert_eq!((plain.weights, plain.weighted_score), (Weights::None, None));
    assert!(
        !app.species_win.message.contains("weighted"),
        "{}",
        app.species_win.message
    );
    // weighted by support, the three that are sure
    app.species_win.weights = Weights::Support;
    app.run_species();
    draw_until(&mut app, &ctx, |a| {
        a.species_win
            .result
            .as_ref()
            .is_some_and(|r| r.weights == Weights::Support)
    });
    let weighted = app.species_win.result.clone().expect("a result");
    assert_ne!(weighted.newick, plain.newick);
    assert!((weighted.weighted_score.expect("a weighted score") - 2.85).abs() < 1e-9);
    let said = &app.species_win.message;
    assert!(
        said.contains("weighted by the support of the gene trees' branches: 2.85"),
        "{said}"
    );
}

#[test]
fn gene_trees_without_what_the_weights_read_are_refused_in_the_window_with_the_reason() {
    let (mut app, ctx) = open_window();
    app.species_win.genes = "((A,B),(C,D));\n((A,C),(B,D));".into();
    app.species_win.weights = Weights::Hybrid;
    app.run_species();
    draw_until(&mut app, &ctx, |a| {
        !a.species_win.message.is_empty() && !a.species_win.is_running()
    });
    assert!(app.species_win.result.is_none());
    assert!(
        app.species_win.message.contains("no branch support"),
        "{}",
        app.species_win.message
    );
}

#[test]
fn weighted_trees_count_as_a_bigger_job_for_the_choice_of_where_to_run() {
    let (mut app, _) = open_window();
    app.species_win.genes = crate::sample::SPECIES_GENE_TREES.to_string();
    let plain = app.species_win.weight().leaves;
    app.species_win.weights = Weights::Hybrid;
    assert_eq!(
        app.species_win.weight().leaves,
        plain * WEIGHTED_WORDS as usize
    );
    // a file of trees too: it was counted when it was opened
    app.species_win.weights = Weights::None;
    app.load_species_file(Source::Bytes {
        name: "genes.tre",
        bytes: crate::sample::SPECIES_GENE_TREES.as_bytes(),
    });
    let file = app.species_win.weight().leaves;
    app.species_win.weights = Weights::Length;
    assert_eq!(
        app.species_win.weight().leaves,
        file * WEIGHTED_WORDS as usize
    );
}

#[test]
fn gene_families_are_chosen_in_the_settings_and_their_result_says_so_with_no_table_of_branches() {
    let (mut app, ctx) = open_window();
    app.species_win.genes =
        "(((A,B),(C,(D,E))),((A,B),(C,(D,E))));\n(((A,B),(C,(D,E))),((A,B),(D,E)));\n".into();
    // as single-copy trees they are refused for the repeated species
    app.run_species();
    draw_until(&mut app, &ctx, |a| {
        !a.species_win.message.is_empty() && !a.species_win.is_running()
    });
    assert!(
        app.species_win.message.contains("more than once"),
        "{}",
        app.species_win.message
    );
    // as families they are read; the weights are not used with them, whatever was chosen
    app.species_win.families = true;
    app.species_win.weights = Weights::Hybrid;
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let result = app.species_win.result.clone().expect("a result");
    assert!(result.families && result.branches.is_empty());
    assert_eq!(result.weights, Weights::None);
    let said = &app.species_win.message;
    assert!(
        said.contains("quartet classes the gene family trees have"),
        "{said}"
    );
}

#[test]
fn a_bootstrap_is_asked_for_in_the_settings_and_its_support_is_in_the_table_and_the_message() {
    let (mut app, ctx) = open_window();
    app.species_win.genes = "((A,B),(C,D));\n".repeat(6);
    // none by default: no support
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    assert!(app
        .species_win
        .result
        .as_ref()
        .is_some_and(|r| r.bootstrap.is_none()));
    // with replicates: the support of every branch, and what was done
    app.species_win.replicates = 10;
    app.run_species();
    draw_until(&mut app, &ctx, |a| {
        a.species_win
            .result
            .as_ref()
            .is_some_and(|r| r.bootstrap.is_some())
    });
    let result = app.species_win.result.clone().expect("a result");
    assert_eq!(result.branches[0].bootstrap, Some(100.0));
    let said = &app.species_win.message;
    assert!(
        said.contains("Bootstrap support from 10 of 10 replicates resampling the genes"),
        "{said}"
    );
    assert!(result.newick.contains(")100:"), "{}", result.newick);
    // the sites need the bootstrap trees of the genes, and say so when there are none
    app.species_win.resample = Resample::Sites;
    app.run_species();
    draw_until(&mut app, &ctx, |a| {
        !a.species_win.is_running() && a.species_win.result.is_none()
    });
    assert!(
        app.species_win.message.contains("bootstrap trees"),
        "{}",
        app.species_win.message
    );
}

#[test]
fn the_polytomy_test_is_asked_for_in_the_settings_and_its_p_values_are_in_the_table_and_the_message(
) {
    let (mut app, ctx) = open_window();
    // thirty gene trees, a third for each way to resolve the branch: the test cannot tell the branch from a polytomy
    app.species_win.genes = ["((A,B),(C,D));\n", "((A,C),(B,D));\n", "((A,D),(B,C));\n"]
        .concat()
        .repeat(10);
    // none by default
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let plain = app.species_win.result.clone().expect("a result");
    assert!(!plain.polytomy && plain.branches.iter().all(|b| b.polytomy.is_none()));
    assert!(!app.species_win.message.contains("Polytomy test"));
    app.species_win.polytomy = true;
    assert!(app.species_win.request().polytomy);
    app.run_species();
    draw_until(&mut app, &ctx, |a| {
        a.species_win.result.as_ref().is_some_and(|r| r.polytomy)
    });
    let result = app.species_win.result.clone().expect("a result");
    let test = result.branches[0]
        .polytomy
        .expect("thirty gene trees are enough");
    assert_eq!((test.genes, test.chi2, test.p), (30, 0.0, 1.0));
    let said = &app.species_win.message;
    assert!(
        said.contains("Polytomy test: 1 of 1 branches tested")
            && said.contains("1 not rejected as polytomies"),
        "{said}"
    );
    // the p-value is the label of the branch in the tree
    assert!(result.newick.contains(")1:"), "{}", result.newick);
    // gene families have nothing to test: the setting is not sent with them
    app.species_win.families = true;
    assert!(!app.species_win.request().polytomy);
    app.species_win.families = false;
    // the table has the column when the test was run, and not when it was not, a branch that was not tested has a dash
    let shown = frame(&mut app, &ctx, vec![]);
    assert!(drew_exactly(&shown, "Polytomy p"));
    let mut few = (*result).clone();
    few.branches[0].polytomy = None;
    app.species_win.result = Some(std::sync::Arc::new(few));
    let shown = frame(&mut app, &ctx, vec![]);
    assert!(
        drew_exactly(&shown, "Polytomy p"),
        "still tested, so the column stays"
    );
    app.species_win.result = Some(plain);
    let shown = frame(&mut app, &ctx, vec![]);
    assert!(!drew_exactly(&shown, "Polytomy p"));
}

#[test]
fn the_bootstrap_trees_of_the_genes_are_sent_only_when_the_sites_are_resampled() {
    let (mut app, _) = open_window();
    app.species_win.genes = "((A,B),(C,D));\n".repeat(2);
    app.species_win.site_trees = Some(("list.txt".into(), vec!["((A,B),(C,D));\n".into(); 2]));
    app.species_win.replicates = 5;
    app.species_win.seed = 11;
    // genes alone: none of them
    let job = app.species_win.request();
    assert!(job.site_trees.is_empty());
    assert_eq!(
        (job.replicates, job.seed, job.resample),
        (5, 11, Resample::Genes)
    );
    app.species_win.resample = Resample::Both;
    assert_eq!(app.species_win.request().site_trees.len(), 2);
    // no replicates: no bootstrap, and nothing carried for it
    app.species_win.replicates = 0;
    assert!(app.species_win.request().site_trees.is_empty());
}

#[test]
fn replicates_make_the_same_trees_a_bigger_job_for_the_choice_of_where_to_run() {
    let (mut app, _) = open_window();
    app.species_win.genes = crate::sample::SPECIES_GENE_TREES.to_string();
    let plain = app.species_win.weight().leaves;
    app.species_win.replicates = 100;
    assert_eq!(app.species_win.weight().leaves, plain * 101);
    app.species_win.weights = Weights::Hybrid;
    assert_eq!(
        app.species_win.weight().leaves,
        plain * 101 * WEIGHTED_WORDS as usize
    );
}

#[test]
fn a_result_that_was_not_weighted_as_asked_is_refused() {
    let found =
        species::run(&SpeciesJob::search(crate::sample::SPECIES_GENE_TREES)).expect("a result");
    let asked = |weights| crate::remote::species::Asked {
        weights,
        families: false,
        bootstrap: None,
        polytomy: false,
    };
    assert!(as_asked(asked(Weights::None), found.clone()).is_ok());
    let why = as_asked(asked(Weights::Hybrid), found).expect_err("not weighted");
    assert!(why.contains("older version"), "{why}");
}

#[test]
fn a_species_tree_says_what_made_it_and_running_again_brings_its_settings_back() {
    let (mut app, ctx) = open_window();
    app.load_species_example();
    app.species_win.polytomy = true;
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let (made, _) = app.species_win.made.clone().expect("a record");
    assert_eq!(
        (made.analysis.as_str(), made.ran_on.as_str()),
        ("Species tree", "this computer")
    );
    assert_eq!(made.inputs.len(), 1);
    assert_eq!(made.inputs[0].name, "gene trees");
    assert_eq!(made.inputs[0].size, app.species_win.genes.len());
    assert!(made
        .parameters
        .contains(&("polytomy test".to_string(), "yes".to_string())));
    assert!(
        !made
            .parameters
            .iter()
            .any(|(k, _)| k == "seed" || k == "resampling"),
        "no bootstrap was asked for, so its settings are not recorded"
    );
    assert_eq!(
        made.programs,
        [format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))]
    );
    let (inputs, settings) = app.species_win.as_it_is();
    assert!(made.changed_since(&inputs, &settings).is_empty());
    // two settings are changed, and a tree is added to the gene trees: all three are said
    app.species_win.polytomy = false;
    app.species_win.lambda = 0.25;
    app.species_win.genes.push_str("\n((A,B),(C,D));");
    let (inputs, settings) = app.species_win.as_it_is();
    assert_eq!(
        made.changed_since(&inputs, &settings),
        [
            "gene trees has changed",
            "polytomy test was yes, is now no",
            "branch length prior was 0.5, is now 0.25"
        ]
    );
    // every other setting is changed as well: Run again must bring each one back
    let kept = app.species_win.made.clone().expect("a record").1;
    let w = &mut app.species_win;
    (w.score_given, w.families, w.exact, w.extra) = (true, true, true, false);
    (w.weights, w.resample) = (Weights::Hybrid, Resample::Both);
    (w.replicates, w.seed, w.root) = (9, 4, "E".to_string());
    assert_ne!(w.kept(), kept);
    // the section is opened and Run again pressed: the settings of the result are back, and the trees are as they are now
    let shown = frame(&mut app, &ctx, vec![]);
    click_on(&mut app, &ctx, &shown, made_ui::TITLE);
    for _ in 0..8 {
        frame(&mut app, &ctx, vec![]);
    }
    let shown = frame(&mut app, &ctx, vec![]);
    click_on(&mut app, &ctx, &shown, made_ui::AGAIN);
    assert!(app.species_win.polytomy && app.species_win.lambda == 0.5);
    assert_eq!(
        app.species_win.kept(),
        kept,
        "every setting of the result is back"
    );
    draw_until(&mut app, &ctx, |a| {
        a.species_win.result.is_some() && a.species_win.made.is_some()
    });
    let (again, _) = app.species_win.made.clone().expect("a new record");
    assert_eq!(again.parameters, made.parameters);
    assert_ne!(again.inputs, made.inputs, "the gene trees had changed");
}

#[test]
fn a_species_job_is_described_with_the_inputs_and_settings_it_uses() {
    let mut job = SpeciesJob::search("((A,B),(C,D));\n");
    let (inputs, settings) = described(&job);
    assert_eq!(inputs.len(), 1);
    let keys: Vec<&str> = settings.iter().map(|(k, _)| k.as_str()).collect();
    assert_eq!(
        keys,
        [
            "weights",
            "gene families",
            "polytomy test",
            "exact search",
            "polytomy resolutions",
            "branch length prior",
            "hung from",
            "bootstrap replicates"
        ]
    );
    assert!(settings.contains(&("hung from".to_string(), "the first species".to_string())));
    assert!(settings.contains(&("bootstrap replicates".to_string(), "0".to_string())));
    // a tree to score, bootstrap trees and a bootstrap add an input each and the settings of the bootstrap
    job.score_tree = "((A,B),(C,D));".into();
    job.site_trees = vec!["((A,B),(C,D));\n".into(), "((A,C),(B,D));\n".into()];
    job.replicates = 20;
    job.root = "C".into();
    let (inputs, settings) = described(&job);
    let names: Vec<&str> = inputs.iter().map(|i| i.name.as_str()).collect();
    assert_eq!(
        names,
        [
            "gene trees",
            "species tree to score",
            "bootstrap trees of the genes"
        ]
    );
    assert_eq!(
        inputs[2].size, 30,
        "the bootstrap trees of all the genes together"
    );
    assert!(settings.contains(&("resampling".to_string(), "genes".to_string())));
    assert!(settings.contains(&("seed".to_string(), "692".to_string())));
    assert!(settings.contains(&("hung from".to_string(), "C".to_string())));
}

#[test]
fn a_result_from_the_workstations_list_or_a_cleared_window_has_no_record() {
    let (mut app, ctx) = open_window();
    app.load_species_example();
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    assert!(app.species_win.made.is_some());
    let found = (*app.species_win.result.clone().expect("a result")).clone();
    app.show_species_result(found);
    assert!(
        app.species_win.made.is_none(),
        "what it was made from is not known here"
    );
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.made.is_some());
    // a run that cannot start (no gene trees) takes the result away, and its record with it
    let kept_text = std::mem::take(&mut app.species_win.genes);
    app.run_species();
    draw_until(&mut app, &ctx, |a| !a.species_win.run_requested);
    assert!(app.species_win.result.is_none());
    assert!(
        app.species_win.made.is_none(),
        "no result, so no record of one"
    );
    app.species_win.genes = kept_text;
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.made.is_some());
    app.species_win.clear();
    assert!(app.species_win.made.is_none());
}

#[test]
fn the_record_of_a_species_tree_is_what_its_copy_button_copies() {
    let (mut app, ctx) = open_window();
    app.load_species_example();
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.made.is_some());
    let record = app.species_win.made.as_ref().expect("a record").0.to_text();
    assert!(record.starts_with("Species tree\nMade on "), "{record}");
    let mut screen = crate::ui::testing::Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::TITLE);
    screen.settle(&mut app);
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::COPY);
    assert_eq!(screen.copied(), Some(record.as_str()));
    assert!(
        app.species_win.message.contains("Copied the record"),
        "{}",
        app.species_win.message
    );
}

#[test]
fn the_link_under_the_record_of_a_species_tree_opens_the_page_of_checks() {
    let (mut app, ctx) = open_window();
    app.load_species_example();
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.made.is_some());
    let mut screen = crate::ui::testing::Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::TITLE);
    screen.settle(&mut app);
    screen.settle(&mut app);
    assert!(!app.checks_win.open);
    screen.click(&mut app, made_ui::CHECKED);
    screen.settle(&mut app);
    assert!(app.checks_win.open);
}

#[test]
fn a_species_tree_is_a_part_of_the_report_with_its_drawing_its_scores_and_the_check_against_astral()
{
    let (mut app, ctx) = open_window();
    app.load_species_example();
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.made.is_some());
    let parts = app.report_parts();
    let part = parts
        .iter()
        .find(|p| p.title == "Species tree")
        .expect("a part for the species tree");
    assert_eq!(
        part.checked.map(|c| c.analysis),
        Some("Species tree from gene trees")
    );
    let result = app.species_win.result.clone().expect("a result");
    let (svg, caption) = part.figure.clone().expect("the tree drawn");
    assert_eq!(
        caption,
        format!("The species tree from {} gene trees.", result.genes)
    );
    assert!(
        svg.matches("</text>").count() >= 4,
        "the names of the species in {svg}"
    );
    assert_eq!(part.tables[0].caption, "Scores");
    assert_eq!(
        part.tables[0].rows,
        [[
            result.genes.to_string(),
            result.score.to_string(),
            format!("{:.4}", result.normalized)
        ]]
    );
    assert_eq!(part.tables[1].rows, [[result.newick.trim().to_string()]]);
    assert!(part.notes.is_empty());
    // a tree whose text cannot be read is not drawn, and the part says so instead of saying nothing
    let mut broken = (*result).clone();
    broken.newick = "((A,B".to_string();
    app.species_win.result = Some(std::sync::Arc::new(broken));
    let parts = app.report_parts();
    let part = parts
        .iter()
        .find(|p| p.title == "Species tree")
        .expect("a part for the species tree");
    assert!(part.figure.is_none());
    assert_eq!(part.notes, [crate::ui::platform::report_texts::NOT_DRAWN]);
    // without its result the window has nothing to report
    app.species_win.clear();
    assert!(app.report_parts().iter().all(|p| p.title != "Species tree"));
}
