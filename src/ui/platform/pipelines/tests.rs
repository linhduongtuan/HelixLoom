use super::*;
use crate::ui::app::Startup;
use crate::ui::testing::Screen;

#[test]
fn choosing_a_tree_aligns_first_and_the_bootstrap_has_its_minimum_size() {
    let mut page = PipelinesPage {
        align: false,
        ..PipelinesPage::default()
    };
    let c = page.choice();
    assert!(c.align && c.tree && c.primers && c.use_installed);
    assert_eq!(c.bootstrap, 0);
    page.tree = false;
    page.bootstrap = true;
    let c = page.choice();
    assert!(!c.align && !c.tree);
    assert_eq!(c.bootstrap, BOOTSTRAP_REPLICATES);
}

#[test]
fn the_target_reads_in_one_based_bases() {
    assert_eq!(target_text(None), "anywhere along it");
    assert_eq!(target_text(Some((200, 60))), "around bases 201 to 260");
}

#[test]
fn the_example_has_six_isolates_longer_than_the_primers_target() {
    let seqs = example_sequences();
    assert_eq!(seqs.len(), 6);
    let (start, length) = EXAMPLE_TARGET;
    assert!(seqs.iter().all(|(_, bases)| bases.len() >= start + length));
}

#[test]
fn the_example_runs_through_all_three_steps_with_the_built_in_methods() {
    let page = PipelinesPage {
        use_installed: false,
        ..PipelinesPage::default()
    };
    let seqs = example_sequences();
    let settings = crate::seq::primer3::Settings::default();
    let done = pipeline::run(
        &seqs,
        Some(EXAMPLE_TARGET),
        &settings,
        &page.choice(),
        &Available::default(),
    )
    .expect("a run");
    assert_eq!(done.log.len(), 3);
    assert_eq!(done.aligned.as_ref().expect("an alignment").rows.len(), 6);
    let newick = &done.tree.as_ref().expect("a tree").newick;
    assert_eq!(NwTree::parse_all(newick).expect("a tree").len(), 1);
    assert!(!done.primers.as_ref().expect("primers").pairs.is_empty());
}

/// Draw frames until a run that was asked for has started and ended.
fn wait(app: &mut App, screen: &mut Screen) {
    let until = std::time::Instant::now() + Duration::from_secs(120);
    loop {
        screen.frame(app, vec![]);
        let page = &app.platform.pipelines;
        if !page.run_requested && !page.is_running() {
            break;
        }
        assert!(std::time::Instant::now() < until, "the run did not end");
        std::thread::sleep(Duration::from_millis(15));
    }
    screen.settle(app);
}

/// The page with its example run through the built-in methods.
fn ran() -> (App, Screen) {
    let mut app = App::new(Startup::default());
    app.platform.pipelines.use_installed = false;
    app.load_example(Workspace::Pipelines);
    app.workspace = Workspace::Pipelines;
    let mut screen = Screen::new();
    wait(&mut app, &mut screen);
    assert!(
        app.platform.pipelines.outcome.is_some(),
        "{}",
        app.platform.pipelines.message
    );
    (app, screen)
}

#[test]
fn a_run_leaves_a_record_of_its_sequences_its_steps_and_the_primer_settings() {
    let (app, _) = ran();
    let page = &app.platform.pipelines;
    let (made, kept) = page.made.clone().expect("a record");
    assert_eq!(
        (made.analysis.as_str(), made.ran_on.as_str()),
        ("Pipeline", "this computer")
    );
    let seqs = example_sequences();
    let names: Vec<&str> = made.inputs.iter().map(|i| i.name.as_str()).collect();
    assert_eq!(names.len(), 7, "{names:?}");
    assert_eq!(
        (names[0], made.inputs[0].size),
        (seqs[0].0.as_str(), seqs[0].1.len())
    );
    assert_eq!(names[6], "primer settings");
    let pair = |k: &str, v: &str| (k.to_string(), v.to_string());
    assert_eq!(
        made.parameters,
        [
            pair("steps", "align, tree, primers"),
            pair("installed programs", "not used"),
            pair("bootstrap replicates", "0"),
            pair("primers go", "around bases 201 to 260")
        ]
    );
    assert_eq!(kept, (Source::Example, [true, true, true, false, false]));
    // the programs are the ones the three parts say made them, each once, and this program last
    let outcome = page.outcome.as_ref().expect("an outcome");
    let (aligned, tree, primers) = (
        outcome.aligned.as_ref().expect("an alignment"),
        outcome.tree.as_ref().expect("a tree"),
        outcome.primers.as_ref().expect("primers"),
    );
    for engine in [&aligned.engine, &tree.engine, &primers.engine] {
        assert_eq!(
            made.programs.iter().filter(|p| *p == engine).count(),
            1,
            "{:?}",
            made.programs
        );
    }
    assert_eq!(
        made.programs.last(),
        Some(&format!("HelixLoom {}", env!("CARGO_PKG_VERSION")))
    );
    assert_eq!(
        made.seconds,
        aligned.seconds + tree.seconds + primers.seconds
    );
    // nothing has changed; then the primer settings do; then the primers are left out and a sequence is closed
    let settings = app.p3_win.settings.clone();
    let choice = page.choice();
    let (inputs, now) = described(&seqs, Some(EXAMPLE_TARGET), &settings, &choice);
    assert!(made.changed_since(&inputs, &now).is_empty());
    let other = crate::seq::primer3::Settings {
        opt_size: settings.opt_size + 1,
        ..settings.clone()
    };
    let (inputs, now) = described(&seqs, Some(EXAMPLE_TARGET), &other, &choice);
    assert_eq!(
        made.changed_since(&inputs, &now),
        ["primer settings has changed"]
    );
    let fewer = Choice {
        primers: false,
        ..choice
    };
    let (inputs, now) = described(&seqs[..5], None, &settings, &fewer);
    assert_eq!(
        made.changed_since(&inputs, &now),
        [
            format!("{} is no longer there", seqs[5].0),
            "primer settings is no longer there".to_string(),
            "steps was align, tree, primers, is now align, tree".to_string(),
            "primers go is no longer a setting".to_string()
        ]
    );
    // without a tree there is no bootstrap to record
    let no_tree = Choice {
        tree: false,
        ..fewer
    };
    let (_, now) = described(&seqs, None, &settings, &no_tree);
    assert!(now.iter().all(|(k, _)| k != "bootstrap replicates"));
}

#[test]
fn the_record_of_a_run_says_on_the_page_what_changed_is_copied_and_runs_it_again_as_it_was_ticked()
{
    let (mut app, mut screen) = ran();
    let (made, kept) = app.platform.pipelines.made.clone().expect("a record");
    screen.click(&mut app, made_ui::TITLE);
    screen.settle(&mut app);
    assert!(
        screen.text().contains(made_ui::UNCHANGED),
        "{}",
        screen.text()
    );
    // a box is unticked: the open record says so
    app.platform.pipelines.primers = false;
    screen.settle(&mut app);
    let text = screen.text();
    assert!(text.contains("primers go is no longer a setting"), "{text}");
    screen.click(&mut app, made_ui::COPY);
    assert_eq!(screen.copied(), Some(made.to_text().as_str()));
    assert!(
        app.platform.pipelines.message.contains("Copied the record"),
        "{}",
        app.platform.pipelines.message
    );
    assert!(!app.checks_win.open);
    screen.click(&mut app, made_ui::CHECKED);
    screen.settle(&mut app);
    assert!(app.checks_win.open);
    app.checks_win.open = false;
    // everything is ticked another way; Run again brings back what the run was ticked with and runs
    let page = &mut app.platform.pipelines;
    (page.tree, page.bootstrap, page.use_installed) = (false, true, true);
    assert_ne!(page.kept(), kept);
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::AGAIN);
    assert_eq!(app.platform.pipelines.kept(), kept);
    wait(&mut app, &mut screen);
    assert_eq!(app.platform.pipelines.message, "Done.", "a run ran");
    let (again, _) = app.platform.pipelines.made.clone().expect("a new record");
    assert_eq!(
        (&again.parameters, &again.inputs),
        (&made.parameters, &made.inputs)
    );
}

#[test]
fn a_run_that_cannot_start_takes_the_result_away_and_its_record_with_it() {
    let (mut app, mut screen) = ran();
    assert!(app.platform.pipelines.made.is_some());
    // the sequences open in Sequences: the one a new window has, which is too few to align
    app.platform.pipelines.source = Source::Open;
    app.platform.pipelines.run_requested = true;
    wait(&mut app, &mut screen);
    let page = &app.platform.pipelines;
    assert!(
        page.message.contains("at least two sequences"),
        "{}",
        page.message
    );
    assert!(page.outcome.is_none() && page.made.is_none());
}

#[test]
fn starting_a_run_takes_the_run_that_was_shown_out_of_what_the_session_keeps() {
    let mut app = App::new(Startup::default());
    // the file mode shows a run too, so the page has something to keep all the way through
    app.show_example_flow_run();
    let run = Outcome {
        log: vec!["Align: 2 sequences.".to_string()],
        ..Outcome::default()
    };
    app.platform
        .pipelines
        .show_fetched(run, "fetched".to_string());
    let before = app.platform.pipelines.stamp().expect("a run to keep");
    // with no sequence open the run does not start, but the run that was shown is gone all the same
    app.docs.clear();
    app.platform.pipelines.source = Source::Open;
    app.start_pipeline(&egui::Context::default());
    assert!(app.platform.pipelines.outcome.is_none());
    let after = app
        .platform
        .pipelines
        .stamp()
        .expect("the file's run is still there");
    assert_ne!(before, after);
    let kept = app.platform.pipelines.kept_text().expect("a run to keep");
    let begins = |word: &str| kept.lines().any(|line| line.starts_with(word));
    assert!(begins("file\t") && !begins("route\t"), "{kept}");
}
