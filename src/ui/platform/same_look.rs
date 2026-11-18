//! Every workspace page is built in one frame (`docs/DESIGN.md`, *The frame of a page*): its controls in one card of one
//! width at the left, its results beside them from one place on, its modes and the figure button in the header, a heading
//! over the first group of controls, and the button that runs the page's analysis filled. The rules these tests hold.

use super::pipelines::flow;
use super::{Workspace, RAIL_WIDTH};
use crate::ui::app::{App, Startup};
use crate::ui::design::{self, Look};
use crate::ui::testing::Screen;
use eframe::egui::{pos2, Color32, Rect};
use std::time::{Duration, Instant};

/// The pages that have controls and results (Home and Sequences are laid out for what they are).
const PAGES: [Workspace; 7] = [
    Workspace::Signals,
    Workspace::Alignment,
    Workspace::Phylogeny,
    Workspace::Trees,
    Workspace::Structure,
    Workspace::Pipelines,
    Workspace::Workstation,
];

/// The test screen, in points.
const SCREEN: (f32, f32) = (1280.0, 800.0);
/// Under the bar at the top of the window.
const UNDER_THE_BAR: f32 = 46.0;

/// An app on the page `ws`, with its example if asked. No program is installed for it, so what it runs is the same on
/// every computer, and it does not wait for a search.
fn app_on(ws: Workspace, example: bool) -> App {
    let mut app = App::new(Startup::default());
    app.platform.no_programs_here();
    app.workspace = ws;
    if example {
        app.load_example(ws);
    }
    app
}

/// The same page in a smaller window: the window's scale makes everything larger, so the window holds fewer points
/// (853 by 533 at 150, 640 by 400 at 200), as a small window does, or a tab of a browser on a small screen.
fn scaled(mut app: App, scale: u16) -> App {
    app.prefs.appearance.scale = scale;
    app
}

/// Draw `app` until what it was asked to run has ended.
fn drawn(mut app: App) -> (App, Screen, Look) {
    let mut screen = Screen::new();
    screen.settle(&mut app);
    let end = Instant::now() + Duration::from_secs(120);
    let waits = |app: &App| {
        let p = &app.platform;
        p.busy_count() > 0 || p.phylogeny.run_when_ready || p.pipelines.run_requested
    };
    while waits(&app) {
        assert!(Instant::now() < end, "{:?}: gave up waiting", app.workspace);
        std::thread::sleep(Duration::from_millis(10));
        screen.frame(&mut app, vec![]);
    }
    screen.settle(&mut app);
    let look = app.prefs.look();
    (app, screen, look)
}

fn page(ws: Workspace, example: bool) -> (App, Screen, Look) {
    drawn(app_on(ws, example))
}

/// The cards of the column of controls: what is filled as a card is and begins at the left edge of the page, under the
/// header. A page has one.
fn column_cards(screen: &Screen, look: &Look) -> Vec<Rect> {
    let left = RAIL_WIDTH + f32::from(design::PAGE_MARGIN);
    screen
        .rects_filled(look.surface)
        .into_iter()
        .filter(|r| (r.min.x - left).abs() < 1.0 && r.min.y > UNDER_THE_BAR && r.height() > 40.0)
        .collect()
}

/// The card of the controls.
fn controls(screen: &Screen, look: &Look) -> Rect {
    *column_cards(screen, look)
        .first()
        .expect("a card of controls at the left edge of the page")
}

/// The room beside the controls, down to the bottom of the screen: where the results are.
fn beside(card: Rect) -> Rect {
    // a few points clear of the card, so that the line a column ends with is not taken for a result
    Rect::from_min_max(
        pos2(card.max.x + 4.0, card.min.y - 1.0),
        pos2(SCREEN.0, SCREEN.1),
    )
}

/// What is filled behind the button of the controls that says `label`: the smallest filled thing its words lie in.
fn fill_behind(screen: &Screen, card: Rect, label: &str) -> Color32 {
    let at = screen
        .find_all(label)
        .into_iter()
        .find(|at| card.contains(*at))
        .unwrap_or_else(|| panic!("no {label:?} in the controls"));
    screen
        .fills_at(at)
        .into_iter()
        .min_by(|a, b| a.0.area().partial_cmp(&b.0.area()).unwrap())
        .unwrap_or_else(|| panic!("nothing filled behind {label:?}"))
        .1
}

#[test]
fn the_controls_of_every_page_are_one_card_of_one_width_and_the_results_begin_at_one_place() {
    let left = RAIL_WIDTH + f32::from(design::PAGE_MARGIN);
    let results = left + design::CONTROLS_WIDTH + f32::from(design::PAGE_MARGIN);
    let mut tops = Vec::new();
    for example in [false, true] {
        for ws in PAGES {
            let (_app, screen, look) = page(ws, example);
            let told = format!("{ws:?}, example {example}");
            let cards = column_cards(&screen, &look);
            assert_eq!(
                cards.len(),
                1,
                "{told}: one card of controls, not {cards:?}"
            );
            let card = cards[0];
            assert!(
                (card.width() - design::CONTROLS_WIDTH).abs() <= 1.0,
                "{told}: the controls are {card:?}, not {} wide",
                design::CONTROLS_WIDTH
            );
            // whatever the page shows there: a card, a view, a row of chips, a line of words
            let begins = screen.leftmost_in(beside(card));
            assert!(
                begins.is_some_and(|x| (x - results).abs() <= 1.0),
                "{told}: the results begin at {begins:?}, not at {results}"
            );
            tops.push((card.min.y, told));
        }
    }
    // and at one height: nothing of a page stands between its header and its controls
    let (first, _) = tops[0];
    for (top, told) in &tops {
        assert!(
            (top - first).abs() <= 0.5,
            "{told}: the controls begin at {top}, not at {first}"
        );
    }
}

#[test]
fn the_first_thing_in_the_controls_of_every_page_is_a_heading() {
    let mut apps: Vec<(String, App)> = [false, true]
        .into_iter()
        .flat_map(|example| {
            PAGES.map(|ws| (format!("{ws:?}, example {example}"), app_on(ws, example)))
        })
        .collect();
    let mut of_a_file = app_on(Workspace::Pipelines, false);
    of_a_file.platform.pipelines.mode = flow::Mode::File;
    apps.push(("Pipelines, of a file".into(), of_a_file));
    for (told, app) in apps {
        let (_app, screen, look) = drawn(app);
        let words = screen.words_in(controls(&screen, &look));
        let first = words.first().expect("words in the controls");
        assert!(
            first.size == design::SECTION_SIZE
                && first.text.chars().any(char::is_alphabetic)
                && first.text == first.text.to_uppercase(),
            "{told}: the controls begin with {:?} in letters of {}, which is no heading",
            first.text,
            first.size
        );
    }
}

#[test]
fn the_modes_of_a_page_and_its_figure_button_are_in_the_header_and_nowhere_else() {
    // once on the page, above the controls and right of them: where a page's actions are
    let in_header = |screen: &Screen, look: &Look, label: &str, told: &str| {
        let card = controls(screen, look);
        let found = screen.find_all(label);
        assert_eq!(found.len(), 1, "{told}: {label:?} is on the page once");
        let at = found[0];
        assert!(
            at.y > UNDER_THE_BAR && at.y < card.min.y && at.x > card.max.x,
            "{told}: {label:?} is at {at:?}, not in the header (the controls are {card:?})"
        );
    };
    for (ws, modes) in [
        (Workspace::Signals, ["Sanger traces", "Reads (FASTQ)"]),
        (Workspace::Pipelines, ["Fixed route", "Pipeline file"]),
    ] {
        for example in [false, true] {
            let (_app, screen, look) = page(ws, example);
            let told = format!("{ws:?}, example {example}");
            modes
                .iter()
                .for_each(|m| in_header(&screen, &look, m, &told));
        }
    }
    // the figure button is there before there is anything to save, pale, so that the header does not move
    for ws in [
        Workspace::Signals,
        Workspace::Alignment,
        Workspace::Phylogeny,
        Workspace::Trees,
        Workspace::Structure,
    ] {
        for example in [false, true] {
            let (_app, screen, look) = page(ws, example);
            let told = format!("{ws:?}, example {example}");
            in_header(&screen, &look, crate::ui::figure::BUTTON, &told);
            assert_eq!(
                screen.is_pale(crate::ui::figure::BUTTON),
                Some(!example),
                "{told}: the figure button is pale until there is something to save"
            );
        }
    }
}

#[test]
fn a_view_is_framed_as_a_card_is_and_begins_where_the_results_begin() {
    let results = RAIL_WIDTH
        + f32::from(design::PAGE_MARGIN)
        + design::CONTROLS_WIDTH
        + f32::from(design::PAGE_MARGIN);
    for ws in [
        Workspace::Signals,
        Workspace::Alignment,
        Workspace::Phylogeny,
        Workspace::Trees,
        Workspace::Structure,
    ] {
        let (_app, screen, look) = page(ws, true);
        let card = controls(&screen, &look);
        let outlined = screen.rects_outlined(look.border);
        // the tall things on the surface of a card beside the controls: the page's view (a card of words is not tall)
        let views: Vec<Rect> = screen
            .rects_filled(look.surface)
            .into_iter()
            .filter(|r| r.min.x > card.max.x && r.min.y >= card.min.y - 1.0 && r.height() > 150.0)
            .filter(|r| (r.min.x - results).abs() <= 1.0)
            .collect();
        assert!(!views.is_empty(), "{ws:?}: a view where the results begin");
        for view in views {
            assert!(
                outlined.contains(&view),
                "{ws:?}: the view {view:?} has the outline of a card"
            );
        }
    }
}

#[test]
fn where_an_analysis_runs_is_said_right_above_the_button_that_runs_it() {
    let mut of_a_file = app_on(Workspace::Pipelines, false);
    of_a_file.platform.pipelines.mode = flow::Mode::File;
    for (told, app, label) in [
        ("Alignment", app_on(Workspace::Alignment, false), "Align"),
        (
            "Alignment, example",
            app_on(Workspace::Alignment, true),
            "Align",
        ),
        (
            "Phylogeny",
            app_on(Workspace::Phylogeny, true),
            "Estimate tree",
        ),
        (
            "Structure",
            app_on(Workspace::Structure, false),
            "Predict structure",
        ),
        (
            "Structure, example",
            app_on(Workspace::Structure, true),
            "Predict structure",
        ),
        (
            "Pipelines",
            app_on(Workspace::Pipelines, false),
            "Run pipeline",
        ),
        ("Pipelines, of a file", of_a_file, "Run pipeline"),
    ] {
        let (_app, screen, look) = drawn(app);
        let words = screen.words_in(controls(&screen, &look));
        let button = words
            .iter()
            .position(|w| w.text == label)
            .unwrap_or_else(|| panic!("{told}: no {label:?} in the controls"));
        let above = &words[button.checked_sub(1).expect("something above the button")];
        assert!(
            above.text.starts_with("Runs on this computer."),
            "{told}: above {label:?} is {:?}, not where it runs",
            above.text
        );
        assert!(
            words[button].at.y - above.at.y < 60.0,
            "{told}: {label:?} is {} under the line that says where it runs",
            words[button].at.y - above.at.y
        );
    }
}

#[test]
fn the_button_that_runs_the_analysis_of_a_page_is_filled_and_pale_until_it_can_be_pressed() {
    let mut of_a_file = app_on(Workspace::Pipelines, false);
    of_a_file.platform.pipelines.mode = flow::Mode::File;
    let mut typed = app_on(Workspace::Workstation, false);
    typed.platform.workstation.address = "workstation.test:7878".into();
    typed.platform.workstation.token = "a token of a test".into();
    // (what, the page, the button, whether it can be pressed)
    let states = [
        (
            "Workstation, nothing typed",
            app_on(Workspace::Workstation, false),
            "Connect",
            false,
        ),
        ("Workstation, typed", typed, "Connect", true),
        (
            "Alignment, empty",
            app_on(Workspace::Alignment, false),
            "Align",
            false,
        ),
        (
            "Alignment, example",
            app_on(Workspace::Alignment, true),
            "Align",
            true,
        ),
        (
            "Phylogeny, empty",
            app_on(Workspace::Phylogeny, false),
            "Estimate tree",
            false,
        ),
        (
            "Phylogeny, example",
            app_on(Workspace::Phylogeny, true),
            "Estimate tree",
            true,
        ),
        // nothing typed to fold, and no program to ask the service with
        (
            "Structure, empty",
            app_on(Workspace::Structure, false),
            "Predict structure",
            false,
        ),
        (
            "Structure, example",
            app_on(Workspace::Structure, true),
            "Predict structure",
            false,
        ),
        (
            "Pipelines, its route",
            app_on(Workspace::Pipelines, false),
            "Run pipeline",
            true,
        ),
        (
            "Pipelines, example",
            app_on(Workspace::Pipelines, true),
            "Run pipeline",
            true,
        ),
        (
            "Pipelines, of a file, no folder",
            of_a_file,
            "Run pipeline",
            false,
        ),
    ];
    let mut pale = None;
    for (told, app, label, ready) in states {
        let (_app, screen, look) = drawn(app);
        let fill = fill_behind(&screen, controls(&screen, &look), label);
        if ready {
            assert_eq!(
                fill, look.accent,
                "{told}: {label:?} is filled with the accent"
            );
            continue;
        }
        // the accent seen through: the same pale on every page
        assert!(
            fill != look.accent && fill.a() < 255 && fill.b() > fill.r() && fill.b() > fill.g(),
            "{told}: {label:?} is on {fill:?}, not on a pale accent ({:?})",
            look.accent
        );
        assert_eq!(*pale.get_or_insert(fill), fill, "{told}: {label:?}");
    }
}

#[test]
fn in_a_small_window_the_frame_of_a_page_is_the_same_and_the_header_is_never_drawn_over_itself() {
    let left = RAIL_WIDTH + f32::from(design::PAGE_MARGIN);
    let results = left + design::CONTROLS_WIDTH + f32::from(design::PAGE_MARGIN);
    for scale in [150, 200] {
        for ws in PAGES {
            let (_app, screen, look) = drawn(scaled(app_on(ws, true), scale));
            let told = format!("{ws:?} at {scale} %");
            let card = controls(&screen, &look);
            assert!(
                (card.width() - design::CONTROLS_WIDTH).abs() <= 1.0,
                "{told}: the controls are {card:?}"
            );
            // a header is two rows at most (its actions under its title, in a window with no room beside it): the page
            // begins right under it, and has most of the window
            let window_height = SCREEN.1 / (f32::from(scale) / 100.0);
            assert!(
                card.min.y < UNDER_THE_BAR + 150.0 && card.height() > window_height / 3.0,
                "{told}: the controls are {card:?} in a window {window_height} high"
            );
            let begins = screen.leftmost_in(beside(card));
            assert!(
                begins.is_some_and(|x| (x - results).abs() <= 1.0),
                "{told}: the results begin at {begins:?}, not at {results}"
            );
            // the header: the title and the line of purpose at the left edge, the page's actions wherever there is room
            let header = Rect::from_min_max(
                pos2(RAIL_WIDTH, UNDER_THE_BAR),
                pos2(SCREEN.0, card.min.y - 1.0),
            );
            let words = screen.words_in(header);
            let (own, actions): (Vec<_>, Vec<_>) = words
                .iter()
                .partition(|w| (w.at.x - left).abs() < 1.0 && (w.size == 22.0 || w.size == 13.0));
            assert_eq!(own.len(), 2, "{told}: a title and a line of purpose");
            // the actions are in the page, not pushed out at its left
            for action in &actions {
                assert!(
                    action.at.x >= left,
                    "{told}: {:?} begins at {}, left of the page",
                    action.text,
                    action.at.x
                );
            }
            for mine in &own {
                // one row each: a line that does not fit is cut short, it does not go on to a second row
                assert!(
                    mine.room.height() < 2.0 * mine.size,
                    "{told}: {:?} is {} high",
                    mine.text,
                    mine.room.height()
                );
                assert!(
                    mine.room.max.x
                        <= SCREEN.0 / (f32::from(scale) / 100.0) - f32::from(design::PAGE_MARGIN)
                            + 1.0,
                    "{told}: {:?} runs out of the window",
                    mine.text
                );
                for action in &actions {
                    assert!(
                        !mine.room.intersects(action.room),
                        "{told}: {:?} is drawn over {:?}",
                        action.text,
                        mine.text
                    );
                }
            }
        }
    }
}

/// The outlined things that begin where the results begin, from the top down: on the Trees page the tree, and under it
/// the card that says what the tree is when the window is too narrow for that card to stand at the tree's right.
fn at_the_results_edge(screen: &Screen, look: &Look) -> Vec<Rect> {
    let results = RAIL_WIDTH
        + f32::from(design::PAGE_MARGIN)
        + design::CONTROLS_WIDTH
        + f32::from(design::PAGE_MARGIN);
    let mut found: Vec<Rect> = screen
        .rects_outlined(look.border)
        .into_iter()
        .filter(|r| (r.min.x - results).abs() <= 1.0 && r.height() > 40.0)
        .collect();
    found.sort_by(|a, b| a.min.y.partial_cmp(&b.min.y).unwrap());
    found
}

/// The card of the controls of the Trees page, its tree, and where the card that says what the tree is begins.
fn tree_and_card(screen: &Screen, look: &Look) -> (Rect, Rect, eframe::egui::Pos2) {
    let card = controls(screen, look);
    let tree = *at_the_results_edge(screen, look).first().expect("the tree");
    let at = screen.find("THIS TREE").expect("what the tree is");
    (card, tree, at)
}

#[test]
fn in_a_window_too_narrow_for_it_the_card_of_a_tree_is_under_the_tree_and_the_tree_keeps_its_room()
{
    // (how much larger everything is, which makes the window that much smaller; whether the card is under the tree).
    // At 110 the tree has 420 points beside the card and at 125 it would have 280: the least it is given there is 360
    for (scale, under) in [(100, false), (110, false), (125, true), (150, true)] {
        let (_app, screen, look) = drawn(scaled(app_on(Workspace::Trees, true), scale));
        let (card, tree, at) = tree_and_card(&screen, &look);
        let told = format!("at {scale} %: the tree is {tree:?}, its card begins at {at:?}");
        assert!(
            at.x > card.max.x,
            "{told}: beside the controls, not in them"
        );
        if under {
            let window = SCREEN.1 / (f32::from(scale) / 100.0);
            assert!(at.y > tree.max.y && at.y < window, "{told}: under the tree");
            // the whole width of the results
            assert!(
                tree.max.x > SCREEN.0 / (f32::from(scale) / 100.0) - 30.0,
                "{told}"
            );
        } else {
            assert!(
                at.x > tree.max.x && at.y < tree.max.y,
                "{told}: at its right"
            );
        }
        assert!(
            tree.width() >= 360.0,
            "{told}: the tree has too little width"
        );
        assert_eq!(
            column_cards(&screen, &look).len(),
            1,
            "{told}: still one card of controls"
        );
    }
}

#[test]
fn picking_a_clade_in_a_small_window_leaves_the_tree_the_room_it_had_and_the_card_in_the_window() {
    for (scale, least) in [(150, 120.0), (200, 50.0)] {
        let window = SCREEN.1 / (f32::from(scale) / 100.0);
        let (mut app, mut screen, look) = drawn(scaled(app_on(Workspace::Trees, true), scale));
        let (_card, before, _) = tree_and_card(&screen, &look);
        // the whole tree is picked: the card lists its 24 tips, which is far more than the card has room for
        app.platform.trees.selected = Some(0);
        screen.settle(&mut app);
        let (_card, after, at) = tree_and_card(&screen, &look);
        let told = format!("at {scale} %");
        assert_eq!(
            before, after,
            "{told}: the tree is where it was and as large"
        );
        assert!(
            after.height() >= least,
            "{told}: the tree is {} high",
            after.height()
        );
        assert!(at.y > after.max.y && at.y < window, "{told}: {at:?}");
        // what the card holds is cut where the card ends, and the card ends in the window
        let under = *at_the_results_edge(&screen, &look)
            .get(1)
            .expect("the card under the tree");
        assert!(
            under.max.y <= window,
            "{told}: the card {under:?} ends out of a window {window} high"
        );
    }
}
