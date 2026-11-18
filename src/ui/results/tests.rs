use super::*;
use crate::provenance::{Asked, Input};
use crate::seq::tree::Node;
use crate::ui::app::Startup;
use crate::ui::platform::Workspace;
use crate::ui::testing::Screen;

fn app() -> App {
    let mut app = App::new(Startup::default());
    app.platform.no_programs_here();
    app
}

/// A tree whose names, lengths and supports are the ones a plain Newick text would lose.
fn tree() -> Node {
    let leaf = |name: &str, length: f64| Node {
        name: name.to_string(),
        length,
        support: None,
        children: Vec::new(),
    };
    let pair = Node {
        name: String::new(),
        length: 0.1 + 0.2,
        support: Some(87.5),
        children: vec![
            leaf("E. coli (K-12); 'lab'", 1e-9),
            leaf("two words", 0.123456789012),
        ],
    };
    Node {
        name: String::new(),
        length: 0.0,
        support: None,
        children: vec![pair, leaf("third:one", 2.0)],
    }
}

fn record() -> Provenance {
    Asked {
        analysis: "tree estimation",
        inputs: vec![Input::of("a row", b"ACGT-")],
        parameters: vec![("method".to_string(), "Bayesian".to_string())],
        ran_on: "this computer",
    }
    .made(Vec::new(), 61.5)
}

const SETTINGS: EstimateOptions = EstimateOptions {
    analysis: Analysis::Bayesian,
    model: DistanceModel::JukesCantor,
    bootstrap: 250,
    generations: 20_000,
};

fn stamp_of(app: &App, kind: &str) -> Option<u64> {
    let found = app.result_stamps().into_iter().find(|(k, _)| *k == kind);
    found.map(|(_, stamp)| stamp)
}

#[test]
fn the_tree_comes_back_exactly_with_the_lines_about_it_and_the_record_and_settings_that_made_it() {
    let mut shown = app();
    shown.show_tree(tree(), "shown");
    shown.tree_win.details = vec!["2 runs of 20000 generations".into(), "a\ttab".into()];
    shown.platform.phylogeny.made = Some((record(), SETTINGS));
    let text = shown.result_text(TREE).expect("a tree");
    let mut back = app();
    back.put_result_back(TREE, &text).expect("put back");
    assert_eq!(back.tree_win.tree, Some(tree()));
    assert_eq!(back.tree_win.details, shown.tree_win.details);
    let (made, o) = back.platform.phylogeny.made.as_ref().expect("the record");
    assert_eq!(*made, record());
    assert!(o.analysis == Analysis::Bayesian && o.model == DistanceModel::JukesCantor);
    assert_eq!((o.bootstrap, o.generations), (250, 20_000));
    assert_eq!(back.platform.trees.message, TREE_IS_BACK);
    assert_eq!(back.result_text(TREE), Some(text));
    // a tree that was opened or pasted, which nothing here estimated, comes back with no record
    shown.show_tree(tree(), "opened");
    let text = shown.result_text(TREE).expect("a tree");
    assert!(!text.contains("made"), "{text}");
    back.put_result_back(TREE, &text).expect("put back");
    assert!(back.tree_win.tree.is_some() && back.platform.phylogeny.made.is_none());
}

#[test]
fn a_kept_tree_that_cannot_be_read_changes_nothing_and_settings_the_page_does_not_allow_are_no_record(
) {
    let mut shown = app();
    shown.show_tree(tree(), "shown");
    let before = (shown.result_stamps(), shown.result_text(TREE));
    for text in [
        "",
        "detail\tonly a line\n",
        "newick\t((a,b);\n",
        "newick\t\n",
    ] {
        assert!(shown.put_result_back(TREE, text).is_err(), "{text:?}");
        assert_eq!((shown.result_stamps(), shown.result_text(TREE)), before);
    }
    assert!(shown.put_result_back("pipeline", "anything").is_err());
    assert_eq!(shown.result_text("pipeline"), None);
    // "Run again" runs with the settings: numbers that the page would not let a person choose are not taken from a file
    let kept = record().kept();
    for (analysis, model, bootstrap, generations) in [
        ("bayes", "jc", "1001", "20000"),
        ("bayes", "jc", "250", "999"),
        ("bayes", "jc", "250", "1000001"),
        ("bayes", "jc", "-1", "20000"),
        ("bayes", "f84", "250", "20000"),
        ("parsimony", "jc", "250", "20000"),
    ] {
        let made = line(&["made", analysis, model, bootstrap, generations, &kept]);
        let text = format!("newick\t(a:1,b:1);\nlater\tline\n{made}");
        shown
            .put_result_back(TREE, &text)
            .expect("the tree is read");
        assert!(shown.platform.phylogeny.made.is_none(), "{made}");
    }
    let made = line(&["made", "bayes", "jc", "1000", "1000000", &kept]);
    shown
        .put_result_back(TREE, &format!("newick\t(a:1,b:1);\n{made}"))
        .expect("the tree is read");
    assert!(shown.platform.phylogeny.made.is_some());
}

#[test]
fn every_way_a_tree_arrives_or_changes_gives_it_another_stamp_and_looking_at_it_does_not() {
    let mut shown = app();
    assert!(
        shown.result_stamps().is_empty(),
        "nothing is shown at the start"
    );
    let mut seen = Vec::new();
    let mut changed = |app: &App, what: &str| {
        let stamp = stamp_of(app, TREE).unwrap_or_else(|| panic!("no tree after {what}"));
        assert!(!seen.contains(&stamp), "{what}");
        seen.push(stamp);
    };
    shown.show_tree(tree(), "shown");
    changed(&shown, "a tree shown");
    shown.show_tree(tree(), "shown");
    changed(&shown, "the same tree shown again");
    shown.platform.phylogeny.made = Some((record(), SETTINGS));
    changed(&shown, "a record");
    // the buttons of the Trees page
    shown.workspace = Workspace::Trees;
    let mut screen = Screen::new();
    screen.settle(&mut shown);
    let before = stamp_of(&shown, TREE);
    screen.settle(&mut shown);
    assert_eq!(
        stamp_of(&shown, TREE),
        before,
        "drawing the page changes nothing"
    );
    screen.click(&mut shown, "Ladderize");
    changed(&shown, "Ladderize");
    screen.click(&mut shown, "Root at midpoint");
    changed(&shown, "Root at midpoint");
    // the Tree window builds one from the open sequences, and takes the one before away first
    for doc in crate::sample::tree_demo() {
        shown.add_doc(doc);
    }
    shown.build_tree();
    changed(&shown, "the Tree window");
    // the Phylogeny page estimates one from the alignment
    let mut estimating = app();
    estimating.load_example(Workspace::Phylogeny);
    estimating.workspace = Workspace::Phylogeny;
    let mut screen = Screen::new();
    for _ in 0..2000 {
        screen.frame(&mut estimating, Vec::new());
        let waiting = estimating.platform.phylogeny.run_when_ready
            || estimating.platform.phylogeny.is_running()
            || estimating.platform.alignment.is_running();
        if !waiting && estimating.tree_win.tree.is_some() {
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(5));
    }
    let estimated = stamp_of(&estimating, TREE).expect("an estimated tree");
    assert!(estimating.platform.phylogeny.made.is_some());
    assert_eq!(estimated % 2, 1, "the stamp of a tree with its record");
    assert!(
        estimated > 1,
        "the estimation counted as a change of the tree"
    );
}

#[test]
fn each_page_gives_its_result_under_its_own_name_and_only_when_it_shows_one() {
    let mut shown = app();
    assert!(shown.result_stamps().is_empty());
    for kind in [ALIGNMENT, TREE, STRUCTURE, PIPELINES] {
        assert_eq!(shown.result_text(kind), None, "{kind}");
    }
    shown.platform.alignment.add("only", b"ACGT");
    assert_eq!(shown.result_stamps().len(), 1);
    assert_eq!(shown.result_stamps()[0].0, ALIGNMENT);
    shown.show_tree(tree(), "shown");
    shown.platform.structure.load_example();
    let run = crate::pipeline::Outcome {
        log: vec!["Align: 2 sequences.".to_string()],
        ..Default::default()
    };
    shown
        .platform
        .pipelines
        .show_fetched(run, "fetched".to_string());
    let kinds: Vec<&str> = shown.result_stamps().iter().map(|(k, _)| *k).collect();
    assert_eq!(kinds, [ALIGNMENT, TREE, STRUCTURE, PIPELINES]);
    // each text goes back to its own page, and to no other
    let mut back = app();
    for kind in kinds {
        let text = shown.result_text(kind).expect(kind);
        back.put_result_back(kind, &text).expect(kind);
        assert_eq!(back.result_text(kind), Some(text.clone()), "{kind}");
        for other in [ALIGNMENT, TREE, STRUCTURE, PIPELINES] {
            let taken = app().put_result_back(other, &text).is_ok();
            assert_eq!(taken, other == kind, "{kind} as {other}");
        }
    }
    assert_eq!(back.platform.alignment.inputs.len(), 1);
    assert!(back.platform.structure.loaded.is_some());
    let run = back.platform.pipelines.outcome.as_ref().expect("a run");
    assert_eq!(run.log, ["Align: 2 sequences."]);
    assert_eq!(
        named([ALIGNMENT, "later", TREE, STRUCTURE, PIPELINES]),
        "the Alignment page, the tree, the structure and the Pipelines page"
    );
    assert_eq!(named([TREE]), "the tree");
}

#[test]
fn a_tree_from_elsewhere_does_not_keep_the_lines_about_the_tree_that_was_there() {
    let mut shown = app();
    shown.show_tree(tree(), "estimated");
    shown.tree_win.details = vec!["Best model by BIC: HKY85".into()];
    shown.platform.phylogeny.made = Some((record(), SETTINGS));
    // another tree is pasted or opened: nothing here made it, so nothing is said of how it was made
    shown.show_tree(tree(), "opened");
    assert!(shown.tree_win.details.is_empty());
    let text = shown.result_text(TREE).expect("a tree");
    assert!(!text.contains("detail") && !text.contains("made"), "{text}");
}

#[test]
fn a_tree_drawn_from_a_pipeline_run_comes_with_the_log_of_that_run_and_not_with_the_lines_of_a_tree_before(
) {
    let mut shown = app();
    shown.workspace = Workspace::Pipelines;
    shown.load_example(Workspace::Pipelines);
    shown.tree_win.details = vec!["of the tree before".into()];
    let mut screen = Screen::new();
    screen.settle(&mut shown);
    let end = std::time::Instant::now() + std::time::Duration::from_secs(120);
    while shown.platform.busy_count() > 0 || shown.platform.pipelines.run_requested {
        assert!(
            std::time::Instant::now() < end,
            "gave up waiting for the run"
        );
        std::thread::sleep(std::time::Duration::from_millis(10));
        screen.frame(&mut shown, Vec::new());
    }
    screen.settle(&mut shown);
    let run = shown.platform.pipelines.outcome.as_ref().expect("a run");
    let log = run.log.clone();
    assert!(!log.is_empty());
    screen.click(&mut shown, "Draw in Trees");
    assert!(shown.tree_win.tree.is_some());
    assert_eq!(shown.tree_win.details, log);
    shown.workspace = Workspace::Pipelines;
    screen.settle(&mut shown);
    // the run is something to keep, and a second run of the same is another one
    let first = stamp_of(&shown, PIPELINES).expect("a run to keep");
    screen.click(&mut shown, "Run pipeline");
    while shown.platform.busy_count() > 0 || shown.platform.pipelines.run_requested {
        assert!(
            std::time::Instant::now() < end,
            "gave up waiting for the run"
        );
        std::thread::sleep(std::time::Duration::from_millis(10));
        screen.frame(&mut shown, Vec::new());
    }
    screen.settle(&mut shown);
    let second = stamp_of(&shown, PIPELINES).expect("a run to keep");
    assert_ne!(first, second);
}
