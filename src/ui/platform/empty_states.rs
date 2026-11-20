//! Every page that waits for something to show says what to do, with buttons: the rule these tests hold.

use super::Workspace;
use crate::ui::app::{App, Startup};
use crate::ui::testing::Screen;

/// Draw the page `ws` of a fresh app.
fn page(ws: Workspace) -> (App, Screen) {
    let mut app = App::new(Startup::default());
    app.workspace = ws;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    (app, screen)
}

/// Click the button that says `label` on the page `ws` of `app`.
fn click(app: &mut App, ws: Workspace, label: &str) {
    app.workspace = ws;
    let mut screen = Screen::new();
    screen.settle(app);
    screen.click(app, label);
}

#[test]
fn a_page_with_nothing_to_show_offers_the_buttons_that_get_it_started() {
    // (button, how many times the page says it): a page's left panel has its own "Open…" button, so the card makes two
    for (ws, buttons) in [
        (
            Workspace::Alignment,
            &[
                ("Add open sequences", 2),
                ("Open FASTA…", 2),
                ("Try the example", 1),
            ][..],
        ),
        (
            Workspace::Signals,
            &[("Open trace…", 2), ("Try the example", 1)][..],
        ),
        (
            Workspace::Trees,
            &[("Open Newick or NEXUS…", 2), ("Try the example", 1)][..],
        ),
        (
            Workspace::Structure,
            &[("Open PDB…", 2), ("Try the example", 1)][..],
        ),
        (
            Workspace::Phylogeny,
            &[("Try the example", 1), ("Add sequences in Alignment", 1)][..],
        ),
        (
            Workspace::Pipelines,
            &[("Try the example", 1), ("Open a sequence", 1)][..],
        ),
    ] {
        let (_, screen) = page(ws);
        let text = screen.text();
        for (label, times) in buttons {
            assert!(
                text.matches(label).count() >= *times,
                "{ws:?} says \"{label}\" fewer than {times} times:\n{text}"
            );
        }
    }
}

fn fresh() -> App {
    App::new(Startup::default())
}

#[test]
fn try_the_example_on_the_alignment_page_fills_it_with_the_example() {
    let mut app = fresh();
    app.platform.alignment.inputs.clear();
    click(&mut app, Workspace::Alignment, "Try the example");
    assert_eq!(
        app.platform.alignment.inputs.len(),
        crate::sample::tree_demo().len()
    );
}

#[test]
fn try_the_example_on_the_signals_trees_and_structure_pages_opens_the_example() {
    let mut app = fresh();
    click(&mut app, Workspace::Signals, "Try the example");
    assert!(app.platform.signals.loaded.is_some());

    let mut app = fresh();
    click(&mut app, Workspace::Trees, "Try the example");
    assert!(app.tree_win.tree.is_some());

    let mut app = fresh();
    click(&mut app, Workspace::Structure, "Try the example");
    assert!(app.platform.structure.loaded.is_some());
}

#[test]
fn try_the_example_on_the_phylogeny_page_loads_the_alignment() {
    let mut app = fresh();
    app.platform.alignment.inputs.clear();
    click(&mut app, Workspace::Phylogeny, "Try the example");
    assert!(!app.platform.alignment.inputs.is_empty());
}

#[test]
fn the_phylogeny_page_leads_on_to_alignment_for_sequences() {
    let mut app = fresh();
    click(&mut app, Workspace::Phylogeny, "Add sequences in Alignment");
    assert_eq!(app.workspace, Workspace::Alignment);
}

#[test]
fn the_pipelines_page_has_an_example_and_a_way_to_a_sequence() {
    let mut app = fresh();
    click(&mut app, Workspace::Pipelines, "Open a sequence");
    assert_eq!(app.workspace, Workspace::Sequences);

    let mut app = fresh();
    app.platform.pipelines.source = super::pipelines::Source::Open;
    click(&mut app, Workspace::Pipelines, "Try the example");
    assert!(matches!(
        app.platform.pipelines.source,
        super::pipelines::Source::Example
    ));
}
