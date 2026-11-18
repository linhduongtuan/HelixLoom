use super::*;
use crate::ui::app::Startup;
use std::time::Instant;

const TREES: &str = "((A:1,B:1):2,C:3);";

/// An answer of the networks for one tree, made up: the networks are the authors' files and are not installed everywhere.
fn answer() -> CnnResult {
    CnnResult {
        task: Task::Bd,
        per_tree: vec![vec![Estimate {
            name: "R0",
            value: 2.5,
        }]],
        skipped: Vec::new(),
        seconds: 0.75,
        engine: "CPU".into(),
        trees: 1,
        tips: 3,
    }
}

fn asked_for(w: &PhyloCnnWin) -> (Asked, Kept) {
    let kept = (w.sampling, w.priors, w.task);
    let (inputs, parameters) = described(&w.newick, kept);
    let asked = Asked {
        analysis: "PhyloCNN estimates",
        inputs,
        parameters,
        ran_on: "this computer",
    };
    (asked, kept)
}

#[test]
fn a_run_is_asked_with_its_trees_and_settings_and_its_record_is_only_of_an_answer() {
    let ctx = egui::Context::default();
    let mut app = App::new(Startup::default());
    let w = &mut app.phylocnn_win;
    (w.newick, w.sampling, w.priors, w.task) = (TREES.into(), 0.4, Priors::Simulated, Task::Bd);
    let expected = asked_for(w);
    // estimates from before, with their record: both go when another run starts
    w.take(report(answer()), Some(expected.clone()));
    assert!(w.made.is_some() && !w.table.is_empty());
    app.run_phylocnn();
    app.tick_phylocnn(&ctx);
    assert!(app.phylocnn_win.made.is_none() && app.phylocnn_win.table.is_empty());
    assert_eq!(
        app.phylocnn_win.asked.as_ref(),
        Some(&expected),
        "what the run is asked is kept when it starts"
    );
    assert_eq!(expected.0.inputs[0].name, "trees");
    assert_eq!(expected.0.inputs[0].size, TREES.len());
    assert_eq!(
        expected.0.parameters,
        [
            ("estimate".to_string(), Task::Bd.label().to_string()),
            ("trained with".to_string(), "general priors".to_string()),
            ("sampling probability".to_string(), "0.4".to_string())
        ]
    );
    let until = Instant::now() + Duration::from_secs(120);
    while app.phylocnn_win.is_running() && Instant::now() < until {
        std::thread::sleep(Duration::from_millis(10));
        app.tick_phylocnn(&ctx);
    }
    assert!(!app.phylocnn_win.is_running(), "the run did not end");
    // where the networks are not installed the run says so and leaves no record; where they are, the record is of that run
    let w = &app.phylocnn_win;
    assert_eq!(w.made.is_some(), !w.table.is_empty(), "{}", w.message);
    if let Some((made, kept)) = &w.made {
        assert_eq!((&made.inputs, kept), (&expected.0.inputs, &expected.1));
    }
}

#[test]
fn the_record_of_the_estimates_names_the_networks_and_says_what_has_changed_since() {
    let mut w = PhyloCnnWin {
        newick: TREES.into(),
        sampling: 0.4,
        task: Task::Bd,
        ..PhyloCnnWin::default()
    };
    let asked = asked_for(&w);
    w.take(report(answer()), Some(asked));
    let (made, kept) = w.made.clone().expect("a record");
    assert_eq!(
        (made.analysis.as_str(), made.ran_on.as_str(), made.seconds),
        ("PhyloCNN estimates", "this computer", 0.75)
    );
    assert_eq!(
        made.programs,
        [
            "PhyloCNN networks of Perez & Gascuel (2025), on the CPU".to_string(),
            format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))
        ]
    );
    assert_eq!(kept, (0.4, Priors::Simulated, Task::Bd));
    let (inputs, settings) = w.as_it_is();
    assert!(made.changed_since(&inputs, &settings).is_empty());
    (w.sampling, w.priors, w.task) = (0.5, Priors::Hiv, Task::Bdss);
    w.newick.push_str("\n((A:1,B:1):2,C:3);");
    let (inputs, settings) = w.as_it_is();
    assert_eq!(
        made.changed_since(&inputs, &settings),
        [
            "trees has changed".to_string(),
            format!(
                "estimate was {}, is now {}",
                Task::Bd.label(),
                Task::Bdss.label()
            ),
            "trained with was general priors, is now HIV priors".to_string(),
            "sampling probability was 0.4, is now 0.5".to_string()
        ]
    );
    // an answer that does not say what the networks ran on (an older workstation's) is not given an "on the" with nothing after
    let asked = asked_for(&w);
    let unnamed = CnnResult {
        engine: String::new(),
        ..answer()
    };
    w.take(report(unnamed), Some(asked));
    assert_eq!(
        w.made.as_ref().expect("a record").0.programs[0],
        "PhyloCNN networks of Perez & Gascuel (2025)"
    );
    // another tree in the box takes the estimates away, and their record with them
    w.set_tree(TREES.into(), String::new());
    assert!(w.table.is_empty() && w.made.is_none());
}

#[test]
fn the_record_of_the_estimates_is_copied_run_again_brings_the_settings_back_and_a_fetched_answer_has_none(
) {
    let mut app = App::new(Startup::default());
    let w = &mut app.phylocnn_win;
    (w.open, w.newick, w.sampling, w.task) = (true, TREES.into(), 0.4, Task::Bd);
    let asked = asked_for(w);
    w.take(report(answer()), Some(asked));
    let note = w.message.clone();
    let record = w.made.as_ref().expect("a record").0.to_text();
    assert!(
        record.starts_with("PhyloCNN estimates\nMade on "),
        "{record}"
    );
    let mut screen = crate::ui::testing::Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::TITLE);
    screen.settle(&mut app);
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::COPY);
    assert_eq!(screen.copied(), Some(record.as_str()));
    assert_eq!(
        app.phylocnn_win.message, note,
        "what is said of the estimates stays"
    );
    assert!(app
        .toasts
        .texts()
        .iter()
        .any(|t| t.contains("Copied the record")));
    assert!(!app.checks_win.open);
    screen.click(&mut app, made_ui::CHECKED);
    screen.settle(&mut app);
    assert!(app.checks_win.open);
    app.checks_win.open = false;
    // the settings are changed; Run again asks for a run with the ones the estimates were made with
    let w = &mut app.phylocnn_win;
    (w.sampling, w.priors, w.task) = (0.9, Priors::Hiv, Task::Bdss);
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::AGAIN);
    let w = &app.phylocnn_win;
    assert_eq!(
        (w.sampling, w.priors, w.task),
        (0.4, Priors::Simulated, Task::Bd)
    );
    assert!(
        w.run_requested || w.is_running() || w.made.is_none(),
        "a run was asked for"
    );
    // an answer fetched from the workstation's list has no record: its trees and settings are not known here
    app.show_phylocnn_result(answer());
    assert!(!app.phylocnn_win.table.is_empty() && app.phylocnn_win.made.is_none());
}
