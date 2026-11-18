//! The pages' figures: the button is clicked on each page as a person clicks it, and what is saved is looked at.

use super::*;
use crate::ui::app::{App, Startup};
use crate::ui::platform::Workspace;
use crate::ui::testing::Screen;
use eframe::egui::{pos2, vec2, Event};
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

/// The app on `page` with its example loaded.
fn on(page: Workspace) -> (App, Screen) {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    app.load_example(page);
    app.workspace = page;
    screen.settle(&mut app);
    (app, screen)
}

/// Choose `kind` in the menu of the page's figure button; the one figure it saves, as its name and its bytes.
fn save_as(app: &mut App, screen: &mut Screen, kind: &str) -> (String, Vec<u8>) {
    let _ = take_saved();
    screen.click(app, BUTTON);
    screen.settle(app);
    screen.click(app, kind);
    screen.settle(app);
    let mut saved = take_saved();
    assert_eq!(saved.len(), 1, "one figure for one click");
    saved.remove(0)
}

fn svg_of(app: &mut App, screen: &mut Screen) -> (String, String) {
    let (name, bytes) = save_as(app, screen, SVG);
    (name, String::from_utf8(bytes).expect("an SVG is text"))
}

/// The wheel of the mouse turned by `delta` points.
fn wheel(delta: Vec2) -> Event {
    Event::MouseWheel {
        unit: eframe::egui::MouseWheelUnit::Point,
        delta,
        phase: eframe::egui::TouchPhase::Move,
        modifiers: eframe::egui::Modifiers::NONE,
    }
}

/// The names of the tips of a tree.
fn tips_of(node: &crate::seq::tree::Node) -> Vec<String> {
    if node.is_leaf() {
        return vec![node.name.clone()];
    }
    node.children.iter().flat_map(tips_of).collect()
}

/// The width and the height a PNG says it has.
fn png_size(bytes: &[u8]) -> (u32, u32) {
    assert_eq!(&bytes[..8], b"\x89PNG\r\n\x1a\n");
    let number = |at: usize| u32::from_be_bytes(bytes[at..at + 4].try_into().expect("four bytes"));
    (number(16), number(20))
}

/// The width and the height an SVG says it has, in points.
fn svg_size(svg: &str) -> (u32, u32) {
    let first = svg.lines().next().expect("a first line");
    let number = |name: &str| {
        let from = first.find(&format!(" {name}=\"")).expect("the size") + name.len() + 3;
        let to = from + first[from..].find('"').expect("the end");
        first[from..to].parse::<f32>().expect("a number") as u32
    };
    (number("width"), number("height"))
}

#[test]
fn the_map_is_saved_with_its_features_and_its_enzymes_and_none_of_the_buttons_beside_it() {
    let (mut app, mut screen) = on(Workspace::Sequences);
    let (name, svg) = svg_of(&mut app, &mut screen);
    assert_eq!(name, "map.svg");
    // the name of the plasmid and its length in the middle, a feature and its filled arrow, an enzyme with its place
    for text in [">pDemo</text>", ">T7 promoter</text>", " bp</text>"] {
        assert!(svg.contains(text), "{text} in {svg}");
    }
    assert!(svg.contains("<path d=\"M"), "the features are filled areas");
    assert!(
        svg.matches("<line ").count() > 10,
        "ticks and the lines to the enzymes"
    );
    // only the map: none of the controls of the panel it is in, and none of the rest of the window
    for control in [
        "Circular",
        "Linear",
        "Enzymes…",
        BUTTON,
        "Search or run a command",
    ] {
        assert!(!svg.contains(control), "{control} in the figure");
    }
    assert!(app.status.contains("Saved map.svg."), "{}", app.status);
    // drawn in room enough for the names of all its enzymes, however small its panel in the window is
    let (width, height) = svg_size(&svg);
    assert!(width >= 900 && height >= 700, "{width} by {height}");
    for line in svg.lines().filter(|l| l.starts_with("<text ")) {
        let y: f32 = line
            .split(" y=\"")
            .nth(1)
            .and_then(|rest| rest.split('"').next())
            .and_then(|y| y.parse().ok())
            .expect("where the text stands");
        assert!(
            y > 8.0 && y < height as f32,
            "a name runs off the figure: {line}"
        );
    }

    // the linear map is another picture of the same plasmid
    app.map_linear = true;
    screen.settle(&mut app);
    let (_, linear) = svg_of(&mut app, &mut screen);
    assert_ne!(linear, svg);
    assert!(linear.contains(">T7 promoter</text>"), "{linear}");

    // and a PNG has four pixels to each point of the SVG
    let (width, height) = svg_size(&linear);
    let (name, png) = save_as(&mut app, &mut screen, PNG);
    assert_eq!(name, "map.png");
    assert_eq!(png_size(&png), (width * 4, height * 4));
}

#[test]
fn the_alignment_is_saved_as_the_window_shows_it_scrolled_where_it_is_scrolled() {
    let (mut app, mut screen) = on(Workspace::Alignment);
    until(&mut screen, &mut app, |a| {
        a.platform.alignment.result.is_some() && !a.platform.alignment.is_running()
    });
    screen.settle(&mut app);
    let (name, first) = svg_of(&mut app, &mut screen);
    assert_eq!(name, "alignment.svg");
    assert_eq!(app.platform.alignment.message, "Saved alignment.svg.");
    let rows = &app
        .platform
        .alignment
        .result
        .as_ref()
        .expect("an alignment")
        .rows;
    assert!(
        first.contains(&format!(">{}</text>", rows[0].0)),
        "the name of a row in {first}"
    );
    assert!(first.contains(">Conservation</text>"), "{first}");
    // the letters are text in the font of sequences, cut where the view ends
    assert!(first.contains("font-family=\"'JetBrains Mono'"), "{first}");
    assert!(first.contains("<g clip-path=\"url(#cut"), "{first}");
    assert!(!first.contains("Copy FASTA") && !first.contains(BUTTON));

    // scrolled to the right in the window, the figure begins further along too
    let over = screen.find("Conservation").expect("the viewer") + vec2(300.0, 60.0);
    screen.frame(
        &mut app,
        vec![Event::PointerMoved(over), wheel(vec2(-400.0, 0.0))],
    );
    screen.settle(&mut app);
    let (_, second) = svg_of(&mut app, &mut screen);
    assert_ne!(second, first);
}

#[test]
fn a_trace_is_saved_with_its_four_curves_and_the_letters_of_its_bases() {
    let (mut app, mut screen) = on(Workspace::Signals);
    let (name, svg) = svg_of(&mut app, &mut screen);
    assert_eq!(name, "trace.svg");
    assert_eq!(app.platform.signals.message, "Saved trace.svg.");
    assert!(
        svg.matches("<polyline ").count() >= 4,
        "a curve for each base in {svg}"
    );
    for base in ["A", "C", "G", "T"] {
        assert!(svg.contains(&format!(">{base}</text>")), "{base} in {svg}");
    }
    assert!(!svg.contains("Export as SCF") && !svg.contains(BUTTON));

    // scrolled along in the window, the figure shows the later part too
    let over = pos2(760.0, 420.0);
    screen.frame(
        &mut app,
        vec![Event::PointerMoved(over), wheel(vec2(-600.0, 0.0))],
    );
    screen.settle(&mut app);
    let (_, later) = svg_of(&mut app, &mut screen);
    assert_ne!(later, svg);
}

#[test]
fn a_tree_is_saved_with_the_names_of_its_tips_from_the_trees_page_and_from_the_phylogeny_page() {
    let (mut app, mut screen) = on(Workspace::Trees);
    let tree = app.tree_win.tree.clone().expect("the example tree");
    let tips = tips_of(&tree);
    let (name, svg) = svg_of(&mut app, &mut screen);
    assert_eq!(name, "tree.svg");
    assert_eq!(app.platform.trees.message, "Saved tree.svg.");
    for tip in &tips {
        assert!(svg.contains(&format!(">{tip}</text>")), "{tip} in {svg}");
    }
    assert!(
        svg.matches("<line ").count() >= tips.len(),
        "a branch for each tip"
    );
    assert!(!svg.contains("Fit to window") && !svg.contains("Save Newick"));
    let (width, height) = svg_size(&svg);
    let (_, png) = save_as(&mut app, &mut screen, PNG);
    assert_eq!(png_size(&png), (width * 4, height * 4));

    // zoomed in the window, the figure is zoomed as far
    let over = pos2(760.0, 420.0);
    screen.frame(
        &mut app,
        vec![Event::PointerMoved(over), wheel(vec2(0.0, 240.0))],
    );
    screen.settle(&mut app);
    let (_, zoomed) = svg_of(&mut app, &mut screen);
    assert_ne!(zoomed, svg);

    // the tree that the Phylogeny page estimated
    let (mut app, mut screen) = on(Workspace::Phylogeny);
    until(&mut screen, &mut app, |a| {
        a.platform.phylogeny.made.is_some() && !a.platform.phylogeny.is_running()
    });
    screen.settle(&mut app);
    let estimated = app.tree_win.tree.clone().expect("the estimated tree");
    let (name, svg) = svg_of(&mut app, &mut screen);
    assert_eq!(name, "tree.svg");
    assert_eq!(app.platform.phylogeny.message, "Saved tree.svg.");
    for tip in tips_of(&estimated) {
        assert!(svg.contains(&format!(">{tip}</text>")), "{tip} in {svg}");
    }
}

#[test]
fn a_structure_is_saved_as_it_is_turned_in_the_window() {
    let (mut app, mut screen) = on(Workspace::Structure);
    let (name, svg) = svg_of(&mut app, &mut screen);
    assert_eq!(name, "structure.svg");
    assert_eq!(app.platform.structure.message, "Saved structure.svg.");
    let drawn = svg.matches("<line ").count()
        + svg.matches("<circle ").count()
        + svg.matches("<polyline ").count()
        + svg.matches("<polygon ").count()
        + svg.matches("<path ").count();
    assert!(
        drawn > 50,
        "a chain of 76 residues is drawn with {drawn} shapes"
    );
    assert!(!svg.contains("Save PDB") && !svg.contains("Reset view"));
    // what the pointer does is said on the screen only
    assert!(screen.text().contains("drag to turn"));
    assert!(!svg.contains("drag to turn"), "{svg}");
    // turned by dragging across it, it is another picture
    let (from, to) = (pos2(800.0, 420.0), pos2(880.0, 420.0));
    let button = |pos, pressed| Event::PointerButton {
        pos,
        button: eframe::egui::PointerButton::Primary,
        pressed,
        modifiers: eframe::egui::Modifiers::NONE,
    };
    screen.frame(&mut app, vec![Event::PointerMoved(from)]);
    screen.frame(&mut app, vec![button(from, true)]);
    screen.frame(&mut app, vec![Event::PointerMoved(to)]);
    screen.frame(&mut app, vec![button(to, false)]);
    screen.settle(&mut app);
    let (_, turned) = svg_of(&mut app, &mut screen);
    assert_ne!(turned, svg);
}

#[test]
fn what_was_asked_for_is_given_once_and_forgotten_when_the_next_frame_begins() {
    let looks = Looks {
        dark: false,
        appearance: Appearance::default(),
    };
    frame_begins(looks);
    assert_eq!(wanted(), None);
    WANTED.set(Some(Format::Png));
    assert_eq!(wanted(), Some(Format::Png));
    assert_eq!(wanted(), None, "once");
    // a page whose view was not drawn in the frame of the click saves nothing in a later frame
    WANTED.set(Some(Format::Svg));
    frame_begins(looks);
    assert_eq!(wanted(), None);
}

#[test]
fn a_figure_has_the_colours_of_the_window_dark_or_light() {
    let (mut app, mut screen) = on(Workspace::Trees);
    let (_, light) = svg_of(&mut app, &mut screen);
    app.run(
        crate::ui::app::Action::Theme(crate::ui::state::ThemeMode::Dark),
        screen.ctx(),
    );
    screen.settle(&mut app);
    let (_, dark) = svg_of(&mut app, &mut screen);
    let ground = |svg: &str| svg.lines().nth(1).expect("the ground").to_string();
    assert_ne!(ground(&light), ground(&dark));
    let look = app.prefs.look();
    let [r, g, b, _] = look.bg.to_srgba_unmultiplied();
    assert!(
        ground(&dark).contains(&format!("fill=\"#{r:02x}{g:02x}{b:02x}\"")),
        "{}",
        ground(&dark)
    );
}
