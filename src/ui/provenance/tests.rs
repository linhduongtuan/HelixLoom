use super::*;
use crate::ui::app::{App, Startup};
use crate::ui::platform::alignment::Engine;
use crate::ui::platform::Workspace;
use crate::ui::testing::Screen;
use std::time::{Duration, Instant};

/// Draw frames until `done` says so (an analysis runs on a thread of its own).
fn until(screen: &mut Screen, app: &mut App, mut done: impl FnMut(&App) -> bool) {
    let end = Instant::now() + Duration::from_secs(120);
    loop {
        screen.frame(app, vec![]);
        if done(app) {
            return;
        }
        assert!(Instant::now() < end, "gave up waiting");
        std::thread::sleep(Duration::from_millis(10));
    }
}

/// The Alignment page with its example aligned.
fn aligned() -> (App, Screen) {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    app.load_example(Workspace::Alignment);
    until(&mut screen, &mut app, |a| {
        a.platform.alignment.result.is_some() && !a.platform.alignment.is_running()
    });
    screen.settle(&mut app);
    (app, screen)
}

#[test]
fn an_alignment_says_what_made_it_and_that_nothing_has_changed() {
    let (mut app, mut screen) = aligned();
    let page = &app.platform.alignment;
    let (made, _) = page.made.as_ref().expect("a record of what made it");
    assert_eq!(made.analysis, "Alignment");
    assert_eq!(made.ran_on, "this computer");
    assert_eq!(
        made.parameters,
        [("engine".to_string(), "built in".to_string())]
    );
    assert_eq!(
        made.inputs.len(),
        page.inputs.iter().filter(|i| i.keep).count()
    );
    assert!(made.inputs.len() >= 2);
    assert!(
        made.programs[0].contains("centre-star"),
        "{:?}",
        made.programs
    );
    assert!(made
        .programs
        .last()
        .is_some_and(|p| p.starts_with("HelixLoom ")));
    let (inputs, parameters) = page.as_it_is();
    assert!(made.changed_since(&inputs, &parameters).is_empty());
    // the section is closed until it is clicked; open, it says the day, the programs and that nothing has changed
    assert!(!screen.text().contains(UNCHANGED));
    screen.click(&mut app, TITLE);
    screen.settle(&mut app);
    screen.settle(&mut app);
    let text = screen.text();
    assert!(text.contains(UNCHANGED), "{text}");
    assert!(text.contains("Settings: engine=built in"), "{text}");
    assert!(text.contains(&format!(
        "Inputs: {}",
        app.platform
            .alignment
            .made
            .as_ref()
            .map_or(0, |m| m.0.inputs.len())
    )));
}

#[test]
fn what_changed_since_is_listed_and_running_again_uses_the_settings_of_the_result() {
    let (mut app, mut screen) = aligned();
    screen.click(&mut app, TITLE);
    screen.settle(&mut app);
    // a sequence is left out and another engine chosen: both are said
    let name = app.platform.alignment.inputs[0].name.clone();
    app.platform.alignment.inputs[0].keep = false;
    app.platform.alignment.choose_engine(Engine::Mafft);
    screen.settle(&mut app);
    let text = screen.text();
    assert!(text.contains("Changed since:"), "{text}");
    assert!(
        text.contains(&format!("• {name} is no longer there")),
        "{text}"
    );
    assert!(
        text.contains("• engine was built in, is now MAFFT"),
        "{text}"
    );
    assert!(!text.contains(UNCHANGED));
    // running again takes the engine of the result (not the one chosen since) and the sequences as they are now
    let before = app
        .platform
        .alignment
        .made
        .as_ref()
        .expect("made")
        .0
        .inputs
        .len();
    screen.click(&mut app, AGAIN);
    until(&mut screen, &mut app, |a| {
        !a.platform.alignment.is_running()
            && a.platform
                .alignment
                .made
                .as_ref()
                .is_some_and(|m| m.0.inputs.len() == before - 1)
    });
    let (made, _) = app.platform.alignment.made.as_ref().expect("a new record");
    assert_eq!(
        made.parameters,
        [("engine".to_string(), "built in".to_string())]
    );
    assert!(!made.inputs.iter().any(|i| i.name == name));
    screen.settle(&mut app);
    assert!(
        screen.text().contains(UNCHANGED),
        "the new result is of the inputs as they are"
    );
}

#[test]
fn the_record_is_copied_as_text_and_an_alignment_from_elsewhere_has_none() {
    let (mut app, mut screen) = aligned();
    screen.click(&mut app, TITLE);
    // the section opens over a few frames, and its buttons are at its foot
    screen.settle(&mut app);
    screen.settle(&mut app);
    screen.click(&mut app, COPY);
    assert!(
        app.platform.alignment.message.contains("Copied the record"),
        "{}",
        app.platform.alignment.message
    );
    let record = app
        .platform
        .alignment
        .made
        .as_ref()
        .expect("a record")
        .0
        .to_text();
    assert_eq!(
        screen.copied(),
        Some(record.as_str()),
        "the record itself is what was copied"
    );
    // an alignment that was made elsewhere (a pipeline, a job fetched from the workstation) carries no record, and shows no section
    let rows = app.platform.alignment.result.clone().expect("a result");
    app.platform.alignment.show_aligned(rows, "From elsewhere.");
    assert!(app.platform.alignment.made.is_none());
    screen.settle(&mut app);
    assert!(!screen.text().contains(TITLE));
}

/// The Phylogeny page with the tree of its example estimated.
fn estimated() -> (App, Screen) {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    app.load_example(Workspace::Phylogeny);
    app.workspace = Workspace::Phylogeny;
    until(&mut screen, &mut app, |a| {
        a.platform.phylogeny.made.is_some() && !a.platform.phylogeny.is_running()
    });
    screen.settle(&mut app);
    (app, screen)
}

#[test]
fn a_tree_says_what_made_it_and_running_again_brings_its_settings_back() {
    use crate::seq::estimate::Analysis;
    let (mut app, mut screen) = estimated();
    let rows = app
        .platform
        .alignment
        .result
        .as_ref()
        .expect("an alignment")
        .rows
        .len();
    let (made, options) = app.platform.phylogeny.made.clone().expect("a record");
    assert_eq!(made.analysis, "Tree estimation");
    assert_eq!(made.inputs.len(), rows, "the rows of the alignment");
    assert_eq!(made.parameters[0].0, "method");
    assert_eq!(
        made.programs,
        [format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))]
    );
    let method = made.parameters[0].1.clone();
    screen.click(&mut app, TITLE);
    screen.settle(&mut app);
    screen.settle(&mut app);
    assert!(screen.text().contains(UNCHANGED), "{}", screen.text());
    // another method is chosen: it is said, and running again puts the first one back and estimates with it
    let other = if options.analysis == Analysis::Upgma {
        Analysis::NeighborJoining
    } else {
        Analysis::Upgma
    };
    app.tree_win.analysis = other;
    // the other settings are changed too, and a row is taken out of the alignment: only a new run can show that
    app.tree_win.bootstrap = options.bootstrap + 7;
    app.tree_win.generations = options.generations + 1000;
    app.tree_win.model = *crate::seq::tree::DistanceModel::ALL
        .iter()
        .find(|m| **m != options.model)
        .expect("another model");
    let gone = app
        .platform
        .alignment
        .result
        .as_mut()
        .and_then(|a| a.rows.pop())
        .expect("a row")
        .0;
    screen.settle(&mut app);
    let text = screen.text();
    assert!(
        text.contains(&format!("• method was {method}, is now")),
        "{text}"
    );
    assert!(
        text.contains(&format!("• {gone} is no longer there")),
        "{text}"
    );
    screen.click(&mut app, AGAIN);
    let w = &app.tree_win;
    assert!(
        w.analysis == options.analysis
            && w.model == options.model
            && w.bootstrap == options.bootstrap
            && w.generations == options.generations,
        "every setting of the result is back"
    );
    until(&mut screen, &mut app, |a| {
        !a.platform.phylogeny.is_running()
            && a.platform
                .phylogeny
                .made
                .as_ref()
                .is_some_and(|m| m.0.inputs.len() == rows - 1)
    });
    let (again, _) = app.platform.phylogeny.made.clone().expect("a record");
    assert_eq!(
        again.parameters, made.parameters,
        "made with the settings of the result"
    );
    assert!(
        !again.inputs.iter().any(|i| i.name == gone),
        "and with the alignment as it is now"
    );
    // the record of the tree is what its Copy button copies
    screen.settle(&mut app);
    screen.click(&mut app, COPY);
    assert_eq!(screen.copied(), Some(again.to_text().as_str()));
    assert!(app.platform.phylogeny.message.contains("Copied the record"));
}

#[test]
fn the_settings_recorded_for_a_tree_are_the_ones_its_method_uses() {
    use crate::pipeline::steps::Aligned;
    use crate::seq::estimate::{Analysis, EstimateOptions};
    use crate::seq::tree::DistanceModel;
    use crate::ui::platform::phylogeny::as_it_is;
    let aligned = Aligned {
        rows: vec![
            ("a".into(), b"AC-T".to_vec()),
            ("b".into(), b"ACGT".to_vec()),
        ],
        engine: String::new(),
        seconds: 0.0,
    };
    let o = |analysis| EstimateOptions {
        analysis,
        model: DistanceModel::ALL[0],
        bootstrap: 100,
        generations: 5000,
    };
    let keys = |analysis| -> Vec<String> {
        as_it_is(&aligned, &o(analysis))
            .1
            .into_iter()
            .map(|(k, _)| k)
            .collect()
    };
    assert_eq!(
        keys(Analysis::NeighborJoining),
        ["method", "distances", "bootstrap resamples"]
    );
    assert_eq!(
        keys(Analysis::Upgma),
        ["method", "distances", "bootstrap resamples"]
    );
    assert_eq!(keys(Analysis::Likelihood), ["method"]);
    assert_eq!(keys(Analysis::Bayesian), ["method", "generations per run"]);
    let (inputs, settings) = as_it_is(&aligned, &o(Analysis::Bayesian));
    assert_eq!(inputs, [Input::of("a", b"AC-T"), Input::of("b", b"ACGT")]);
    assert_eq!(settings[0].1, "Bayesian");
    assert_eq!(settings[1].1, "5000");
}

#[test]
fn a_tree_that_came_from_elsewhere_carries_no_record_of_the_estimation_before_it() {
    let (mut app, mut screen) = estimated();
    assert!(app.platform.phylogeny.made.is_some());
    let tree = app.tree_win.tree.clone().expect("a tree");
    app.show_tree(tree, "Opened from a file.");
    assert!(app.platform.phylogeny.made.is_none());
    app.workspace = Workspace::Phylogeny;
    screen.settle(&mut app);
    assert!(!screen.text().contains(TITLE));
}

#[test]
fn a_tree_built_in_the_tree_window_or_changed_on_the_trees_page_is_not_shown_under_the_estimations_record(
) {
    // the Phylogenetic Tree window builds its own tree of the open sequences: the record of the estimation does not describe it
    let (mut app, mut screen) = estimated();
    app.build_tree();
    assert!(
        app.platform.phylogeny.made.is_none(),
        "Build tree replaced the tree"
    );
    app.workspace = Workspace::Phylogeny;
    screen.settle(&mut app);
    assert!(!screen.text().contains(TITLE));
    // rooting at the midpoint and ladderizing change the tree that is shown: it is no longer the one the estimation made
    for button in ["Root at midpoint", "Ladderize"] {
        let (mut app, mut screen) = estimated();
        app.workspace = Workspace::Trees;
        screen.settle(&mut app);
        screen.click(&mut app, button);
        assert!(app.platform.phylogeny.made.is_none(), "{button}");
    }
}

#[test]
fn mafft_that_is_not_installed_where_the_run_goes_is_not_what_the_record_says_ran() {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    app.platform.no_programs_here();
    for doc in crate::sample::tree_demo() {
        app.platform.alignment.add(&doc.name, &doc.seq);
    }
    app.workspace = Workspace::Alignment;
    app.platform.alignment.choose_engine(Engine::Mafft);
    screen.settle(&mut app);
    screen.click(&mut app, "Align");
    until(&mut screen, &mut app, |a| {
        a.platform.alignment.made.is_some() && !a.platform.alignment.is_running()
    });
    let page = &app.platform.alignment;
    let (made, engine) = page.made.as_ref().expect("a record");
    // the built-in aligner ran, and the record, the kept engine and the page all say so
    assert!(
        made.programs[0].contains("centre-star"),
        "{:?}",
        made.programs
    );
    assert_eq!(
        made.parameters,
        [("engine".to_string(), "built in".to_string())]
    );
    assert!(*engine == Engine::BuiltIn);
    let (inputs, settings) = page.as_it_is();
    assert!(
        made.changed_since(&inputs, &settings).is_empty(),
        "the page's choice went back to what ran"
    );
}

#[test]
fn the_link_under_a_record_opens_the_page_that_says_how_each_analysis_was_checked() {
    for (mut app, mut screen) in [aligned(), estimated()] {
        assert!(!app.checks_win.open);
        screen.click(&mut app, TITLE);
        screen.settle(&mut app);
        screen.settle(&mut app);
        screen.click(&mut app, CHECKED);
        screen.settle(&mut app);
        assert!(app.checks_win.open);
        assert!(screen.text().contains(crate::ui::checks::TITLE));
    }
}
