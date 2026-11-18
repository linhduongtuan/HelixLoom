use super::*;
use std::sync::atomic::{AtomicU32, Ordering};
use std::{fs, path::PathBuf};

static NEXT: AtomicU32 = AtomicU32::new(0);

/// A folder of its own for one test, removed when the test ends.
struct Folder(PathBuf);

impl Folder {
    fn new() -> Folder {
        let path = std::env::temp_dir().join(format!(
            "helixloom-tour-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[test]
fn the_first_run_with_a_settings_folder_starts_at_the_first_card() {
    let f = Folder::new();
    let t = Tour::first_run(Some(&f.0), false);
    assert_eq!(t.number(), Some((1, 3)));
    assert_eq!(t.current().map(|s| s.target), Some(Target::Search));
}

#[test]
fn no_tour_without_a_folder_to_remember_it_in_or_when_something_was_asked_for() {
    let f = Folder::new();
    assert!(Tour::first_run(None, false).current().is_none());
    assert!(Tour::first_run(Some(&f.0), true).current().is_none());
}

#[test]
fn it_goes_through_the_cards_in_order_and_ends_after_the_last() {
    let f = Folder::new();
    let mut t = Tour::first_run(Some(&f.0), false);
    let targets: Vec<Target> = std::iter::from_fn(|| {
        let s = t.current()?;
        t.advance().expect("kept");
        Some(s.target)
    })
    .collect();
    assert_eq!(targets, [Target::Search, Target::Rail, Target::Jobs]);
    assert!(t.current().is_none());
}

#[test]
fn only_the_last_card_is_the_last() {
    let f = Folder::new();
    let mut t = Tour::first_run(Some(&f.0), false);
    assert!(!t.is_last());
    t.advance().unwrap();
    assert!(!t.is_last());
    t.advance().unwrap();
    assert!(t.is_last());
}

#[test]
fn a_tour_that_ended_is_not_shown_again_by_the_next_run() {
    let f = Folder::new();
    let mut t = Tour::first_run(Some(&f.0), false);
    t.finish().unwrap();
    assert!(Tour::first_run(Some(&f.0), false).current().is_none());
}

#[test]
fn skipping_counts_as_seeing_it() {
    let f = Folder::new();
    let mut t = Tour::first_run(Some(&f.0), false);
    t.advance().unwrap();
    t.finish().unwrap();
    assert!(t.current().is_none());
    assert!(f.0.join(SEEN).exists());
}

#[test]
fn asking_for_the_tour_again_shows_it_and_does_not_forget_that_it_was_seen() {
    let f = Folder::new();
    let mut t = Tour::first_run(Some(&f.0), false);
    t.finish().unwrap();
    t.replay();
    assert_eq!(t.number(), Some((1, 3)));
    assert!(f.0.join(SEEN).exists());
}

#[test]
fn a_mark_that_cannot_be_written_is_an_error_but_the_tour_still_ends() {
    let f = Folder::new();
    // a file where the folder should be: nothing can be made inside it
    let blocked = f.0.join("blocked");
    fs::write(&blocked, "x").unwrap();
    let mut t = Tour::first_run(Some(&blocked.join("settings")), false);
    assert!(t.current().is_some());
    assert!(t.finish().is_err());
    assert!(t.current().is_none());
}

#[test]
fn the_card_stays_inside_the_window() {
    let screen = Rect::from_min_size(Pos2::ZERO, vec2(1280.0, 800.0));
    let spots = [
        // the search box at the right edge, the pages at the left, a pill at the very corner
        (
            Some(Rect::from_min_size(pos2(1100.0, 8.0), vec2(170.0, 30.0))),
            Target::Search,
        ),
        (
            Some(Rect::from_min_size(pos2(0.0, 50.0), vec2(72.0, 700.0))),
            Target::Rail,
        ),
        (
            Some(Rect::from_min_size(pos2(1250.0, 770.0), vec2(30.0, 30.0))),
            Target::Jobs,
        ),
        (None, Target::Jobs),
    ];
    for (target, which) in spots {
        let p = card_position(target, which, screen);
        let card = Rect::from_min_size(p, CARD);
        assert!(screen.contains_rect(card), "{which:?} {target:?}: {card:?}");
    }
}

#[test]
fn the_card_sits_below_the_search_box_and_beside_the_pages() {
    let screen = Rect::from_min_size(Pos2::ZERO, vec2(1280.0, 800.0));
    let search = Rect::from_min_size(pos2(900.0, 8.0), vec2(290.0, 30.0));
    let p = card_position(Some(search), Target::Search, screen);
    assert!(p.y >= search.bottom());
    assert!((p.x + CARD.x - search.right()).abs() < 1.0);
    let rail = Rect::from_min_size(pos2(0.0, 50.0), vec2(72.0, 700.0));
    let p = card_position(Some(rail), Target::Rail, screen);
    assert!(p.x >= rail.right());
}

// --- the tour in the window ---

use crate::ui::app::{App, Startup};
use crate::ui::testing::Screen;
use eframe::egui::Event;

fn first_run_app(f: &Folder) -> App {
    App::new(Startup {
        settings_dir: Some(f.0.clone()),
        ..Startup::default()
    })
}

#[test]
fn the_first_run_shows_the_first_card_and_finds_what_it_points_at() {
    let f = Folder::new();
    let mut app = first_run_app(&f);
    let mut screen = Screen::new();
    screen.settle(&mut app);
    let text = screen.text();
    assert!(text.contains("Search for anything"), "{text}");
    assert!(text.contains("STEP 1 OF 3"), "{text}");
    for target in [Target::Search, Target::Rail, Target::Jobs] {
        let rect = app.tour.rect(target).expect("the window says where it is");
        assert!(rect.width() > 0.0 && rect.height() > 0.0, "{target:?}");
    }
}

#[test]
fn clicking_through_the_cards_ends_the_tour_and_the_next_run_is_without_it() {
    let f = Folder::new();
    let mut app = first_run_app(&f);
    let mut screen = Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, "Next");
    assert!(screen.text().contains("The pages"), "{}", screen.text());
    screen.click(&mut app, "Next");
    assert!(screen.text().contains("Where the jobs run"));
    screen.click(&mut app, "Done");
    assert!(app.tour.current().is_none());
    assert!(!screen.text().contains("Where the jobs run"));

    let mut again = first_run_app(&f);
    let mut screen = Screen::new();
    screen.settle(&mut again);
    assert!(again.tour.current().is_none());
}

#[test]
fn skipping_ends_it_at_once() {
    let f = Folder::new();
    let mut app = first_run_app(&f);
    let mut screen = Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, "Skip the tour");
    assert!(app.tour.current().is_none());
    assert!(f.0.join(SEEN).exists());
}

#[test]
fn no_card_shows_when_the_program_has_nowhere_to_keep_the_mark() {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert!(!screen.text().contains("STEP 1 OF 3"));
}

#[test]
fn the_last_card_opens_the_sample_project() {
    let f = Folder::new();
    let mut app = first_run_app(&f);
    app.docs.clear();
    app.tour.advance().unwrap();
    app.tour.advance().unwrap();
    let mut screen = Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, "Open the sample project");
    assert!(app.tour.current().is_none());
    assert!(app
        .docs
        .iter()
        .any(|d| d.doc.name == crate::sample::DEMO_NAME));
    assert!(app.platform.signals.loaded.is_some());
    assert!(app.tree_win.tree.is_some());
    assert_eq!(app.workspace, Workspace::Sequences);
}

#[test]
fn the_sample_project_does_not_open_the_plasmid_twice() {
    let mut app = App::new(Startup::default());
    app.open_sample_project();
    app.open_sample_project();
    let plasmids = app
        .docs
        .iter()
        .filter(|d| d.doc.name == crate::sample::DEMO_NAME)
        .count();
    assert_eq!(plasmids, 1);
}

#[test]
fn help_can_bring_the_tour_back() {
    let f = Folder::new();
    let mut app = first_run_app(&f);
    app.tour.finish().unwrap();
    let ctx = egui::Context::default();
    app.run(Action::Tour, &ctx);
    assert_eq!(app.tour.number(), Some((1, 3)));
}

#[test]
fn the_card_does_not_panic_in_a_window_smaller_than_itself() {
    let target = Some(Rect::from_min_size(pos2(10.0, 10.0), vec2(30.0, 30.0)));
    for size in [vec2(300.0, 150.0), vec2(320.0, 568.0), vec2(0.0, 0.0)] {
        let screen = Rect::from_min_size(Pos2::ZERO, size);
        for which in [Target::Search, Target::Rail, Target::Jobs] {
            for target in [target, None] {
                let p = card_position(target, which, screen);
                assert!(p.x >= screen.left() + MARGIN && p.y >= screen.top() + MARGIN);
            }
        }
    }
}

fn escape() -> Vec<Event> {
    vec![Event::Key {
        key: Key::Escape,
        physical_key: Some(Key::Escape),
        pressed: true,
        repeat: false,
        modifiers: Modifiers::NONE,
    }]
}

#[test]
fn escape_ends_the_tour() {
    let f = Folder::new();
    let mut app = first_run_app(&f);
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert!(app.tour.current().is_some());
    screen.frame(&mut app, escape());
    assert!(app.tour.current().is_none());
    assert!(f.0.join(SEEN).exists());
}

#[test]
fn escape_that_closes_the_palette_leaves_the_tour_where_it_is() {
    let f = Folder::new();
    let mut app = first_run_app(&f);
    let mut screen = Screen::new();
    screen.settle(&mut app);
    app.palette.toggle();
    screen.settle(&mut app);
    assert!(app.palette.is_open());
    screen.frame(&mut app, escape());
    assert!(!app.palette.is_open(), "the palette took the key");
    assert_eq!(app.tour.number(), Some((1, 3)));
    // the next Esc is the tour's
    screen.frame(&mut app, escape());
    assert!(app.tour.current().is_none());
}

#[test]
fn in_a_browser_the_tour_is_shown_at_a_first_visit_and_that_it_was_seen_is_kept_by_the_browser() {
    use crate::ui::session::{Memory, Shelf, Store};
    let table = std::rc::Rc::new(Memory::default());
    let visit = || Tour::with(Some(Store::on(table.clone())), false);
    let mut first = visit();
    assert_eq!(first.number(), Some((1, STEPS.len())));
    assert!(table.get(SEEN).is_none(), "not until it has ended");
    first.finish().expect("kept");
    assert!(table.get(SEEN).is_some());
    assert!(visit().current().is_none(), "not at the next visit");
    // a visit that asks for something particular has no tour, seen or not
    table.remove(SEEN);
    assert!(Tour::with(Some(Store::on(table.clone())), true)
        .current()
        .is_none());
    // a browser that keeps nothing: the tour ends for this visit, and why it was not kept is told
    table.room.set(Some(0));
    let mut full = visit();
    assert!(full.finish().is_err());
    assert!(full.current().is_none());
}
