use super::*;
use crate::ui::app::{App, Startup};
use crate::ui::platform::Workspace;
use crate::ui::testing::Screen;

fn cut() -> String {
    crate::sample::ubiquitin_in_two_chains()
}

fn app_with(text: &str) -> (App, Screen) {
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Structure;
    app.platform
        .structure
        .show_pdb("cut ubiquitin", text)
        .expect("a structure");
    let mut screen = Screen::new();
    screen.settle(&mut app);
    (app, screen)
}

/// Draw frames until the measuring is over, and a few more.
fn finish(app: &mut App, screen: &mut Screen) {
    let until = std::time::Instant::now() + Duration::from_secs(60);
    while app.platform.structure.interface.is_running() {
        assert!(std::time::Instant::now() < until, "no answer");
        std::thread::sleep(Duration::from_millis(10));
        screen.frame(app, vec![]);
    }
    screen.settle(app);
}

fn analyse(app: &mut App, screen: &mut Screen) {
    screen.click(app, "Analyse the interface");
    finish(app, screen);
}

#[test]
fn a_structure_of_one_chain_says_an_interface_needs_two() {
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Structure;
    app.platform.structure.load_example();
    let mut screen = Screen::new();
    screen.settle(&mut app);
    let text = screen.text();
    assert!(text.contains("INTERFACE"), "{text}");
    assert!(
        text.contains("An interface needs at least two chains; this structure has one."),
        "{text}"
    );
    assert!(!text.contains("Analyse the interface"));
}

#[test]
fn a_structure_of_two_chains_starts_with_the_first_chain_against_the_rest() {
    let (app, screen) = app_with(&cut());
    let sides = &app.platform.structure.interface.sides;
    assert_eq!(sides, &[('A', Side::First), ('B', Side::Second)]);
    let text = screen.text();
    for want in [
        "Chain A",
        "Chain B",
        "Analyse the interface",
        "Which chains meet",
    ] {
        assert!(text.contains(want), "no \"{want}\":\n{text}");
    }
}

#[test]
fn analysing_shows_the_area_the_contacts_and_a_table_of_the_residues() {
    let (mut app, mut screen) = app_with(&cut());
    analyse(&mut app, &mut screen);
    let page = &app.platform.structure.interface;
    let done = page.measured.as_ref().expect("a measurement");
    // the numbers the independent implementation gives for this cut (see seq/interface/tests.rs)
    assert!(
        (done.interface.buried_total() - 2021.15).abs() < 3.0,
        "{}",
        done.interface.buried_total()
    );
    assert_eq!(done.interface.contacts.len(), 65);
    let text = screen.text();
    for want in [
        "interface 1011 Å²",
        "65 contacts",
        "Buried: 2021 Å² in all",
        "Buried Å²",
        "A ILE 23",
        "Colour by buried area",
        "Copy the table",
    ] {
        assert!(text.contains(want), "no \"{want}\":\n{text}");
    }
    assert!(
        page.message.ends_with("residues at the interface."),
        "{}",
        page.message
    );
}

#[test]
fn the_structure_is_coloured_by_what_each_residue_buries() {
    let (mut app, mut screen) = app_with(&cut());
    analyse(&mut app, &mut screen);
    assert!(app
        .platform
        .structure
        .loaded
        .as_ref()
        .unwrap()
        .fold
        .residues
        .iter()
        .all(|r| r.buried.is_none()));
    screen.click(&mut app, "Colour by buried area");
    let page = &app.platform.structure;
    assert_eq!(page.style.scheme, Scheme::Interface);
    let fold = &page.loaded.as_ref().unwrap().fold;
    let at = |n: i32| {
        fold.residues
            .iter()
            .find(|r| r.number == n)
            .and_then(|r| r.buried)
    };
    assert!((at(23).unwrap() - 118.87).abs() < 0.6, "{:?}", at(23));
    assert_eq!(at(33), Some(0.0));
    assert!(
        screen.text().contains("buries 100 Å² or more"),
        "the legend says what the colours mean"
    );
}

#[test]
fn a_row_of_the_table_picks_its_residue_in_the_structure() {
    let (mut app, mut screen) = app_with(&cut());
    analyse(&mut app, &mut screen);
    screen.click(&mut app, "A ILE 23");
    let page = &app.platform.structure;
    let picked = page.style.picked.expect("a residue is picked");
    let r = &page.loaded.as_ref().unwrap().fold.residues[picked];
    assert_eq!((r.chain, r.number), ('A', 23));
}

#[test]
fn the_table_can_be_copied_for_a_spreadsheet() {
    let (mut app, mut screen) = app_with(&cut());
    analyse(&mut app, &mut screen);
    let rows = app
        .platform
        .structure
        .interface
        .measured
        .as_ref()
        .unwrap()
        .rows
        .len();
    screen.click(&mut app, "Copy the table");
    assert_eq!(
        app.platform.structure.interface.message,
        format!("Copied {rows} rows.")
    );
}

#[test]
fn changing_a_side_takes_back_what_was_measured_and_the_colours() {
    let (mut app, mut screen) = app_with(&cut());
    analyse(&mut app, &mut screen);
    screen.click(&mut app, "Colour by buried area");
    // the first "left out" button is chain A's
    screen.click(&mut app, "left out");
    let page = &app.platform.structure;
    assert_eq!(page.interface.sides[0], ('A', Side::Neither));
    assert!(page.interface.measured.is_none());
    assert_ne!(page.style.scheme, Scheme::Interface);
    assert!(page
        .loaded
        .as_ref()
        .unwrap()
        .fold
        .residues
        .iter()
        .all(|r| r.buried.is_none()));
}

#[test]
fn a_new_structure_forgets_the_measurement_of_the_old_one() {
    let (mut app, mut screen) = app_with(&cut());
    analyse(&mut app, &mut screen);
    app.platform.structure.show_pdb("again", &cut()).unwrap();
    assert!(app.platform.structure.interface.measured.is_none());
    assert_eq!(app.platform.structure.interface.sides.len(), 2);
    screen.settle(&mut app);
    assert!(!screen.text().contains("Colour by buried area"));
}

/// A measurement that takes a moment, in the place of the real one: to change things while it goes on.
fn slow_job(app: &mut App) {
    let ctx = egui::Context::default();
    app.platform.structure.interface.job = Some(Job::spawn(&ctx, || {
        std::thread::sleep(Duration::from_millis(300));
        Ok(Interface::default())
    }));
    app.platform.structure.interface.message = "Measuring…".into();
}

#[test]
fn changing_a_side_while_a_measurement_goes_on_drops_it_and_its_message() {
    let (mut app, mut screen) = app_with(&cut());
    slow_job(&mut app);
    screen.click(&mut app, "left out");
    let page = &app.platform.structure.interface;
    assert!(
        !page.is_running(),
        "the measuring for the old sides is dropped"
    );
    assert!(page.message.is_empty(), "{}", page.message);
    std::thread::sleep(Duration::from_millis(400));
    screen.settle(&mut app);
    assert!(
        app.platform.structure.interface.measured.is_none(),
        "its answer never arrives"
    );
}

#[test]
fn showing_another_structure_while_a_measurement_goes_on_drops_it() {
    let (mut app, mut screen) = app_with(&cut());
    slow_job(&mut app);
    app.platform.structure.show_pdb("another", &cut()).unwrap();
    assert!(!app.platform.structure.interface.is_running());
    std::thread::sleep(Duration::from_millis(400));
    screen.settle(&mut app);
    assert!(app.platform.structure.interface.measured.is_none());
}

#[test]
fn analysing_again_after_colouring_does_not_leave_the_structure_all_one_grey() {
    let (mut app, mut screen) = app_with(&cut());
    analyse(&mut app, &mut screen);
    screen.click(&mut app, "Colour by buried area");
    assert_eq!(app.platform.structure.style.scheme, Scheme::Interface);
    analyse(&mut app, &mut screen);
    let page = &app.platform.structure;
    assert_ne!(
        page.style.scheme,
        Scheme::Interface,
        "no colours of the interface until they are asked for again"
    );
    assert!(page
        .loaded
        .as_ref()
        .unwrap()
        .fold
        .residues
        .iter()
        .all(|r| r.buried.is_none()));
    assert!(page.interface.measured.is_some());
}

#[test]
fn a_row_that_is_not_a_drawn_residue_says_so_and_keeps_the_pick() {
    let (mut app, mut screen) = app_with(&cut());
    analyse(&mut app, &mut screen);
    app.platform.structure.style.picked = Some(3);
    let ctx = egui::Context::default();
    app.apply_interface_clicks(
        InterfaceClicks {
            pick: Some(usize::MAX),
            ..InterfaceClicks::default()
        },
        &ctx,
    );
    let page = &app.platform.structure;
    assert_eq!(page.style.picked, Some(3));
    assert!(
        page.interface
            .message
            .contains("not part of the drawn protein chain"),
        "{}",
        page.interface.message
    );
}

#[test]
fn an_insertion_code_shows_in_the_table() {
    // residue 23 of chain A as 23A: a file that numbers an antibody this way has such residues
    let text: String = cut()
        .lines()
        .map(
            |l| match (l.starts_with("ATOM"), l.get(22..26).map(str::trim)) {
                (true, Some("23")) => format!("{}A{}", &l[..26], &l[27..]),
                _ => l.to_string(),
            },
        )
        .collect::<Vec<_>>()
        .join("\n");
    let (mut app, mut screen) = app_with(&text);
    analyse(&mut app, &mut screen);
    let done = app.platform.structure.interface.measured.as_ref().unwrap();
    assert!(done.rows.iter().any(|r| r.number_text() == "23A"));
    assert!(interface::tsv(&done.rows).contains("\tILE\t23A\t"));
    assert!(screen.text().contains("A ILE 23A"));
}

#[test]
fn the_note_describes_a_contact_as_the_pair_of_residues_that_the_count_is_of() {
    let (mut app, mut screen) = app_with(&cut());
    analyse(&mut app, &mut screen);
    // the chip counts residue pairs (65 for this cut), so the note must not call a contact a pair of atoms
    assert!(screen.text().contains("65 contacts"));
    assert!(EXPLANATION.contains("A contact is a pair of residues, one from each side"));
    assert!(EXPLANATION.contains("no angle is measured, so some are not real bonds"));
}

#[test]
fn the_antibody_example_loads_a_real_complex_with_the_antibody_on_one_side_and_lysozyme_on_the_other(
) {
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Structure;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, "Antibody");
    let page = &app.platform.structure;
    assert!(page.loaded.as_ref().unwrap().name.contains("1VFB"));
    assert_eq!(
        page.interface.sides,
        [('A', Side::First), ('B', Side::First), ('C', Side::Second)]
    );
    analyse(&mut app, &mut screen);
    let done = app.platform.structure.interface.measured.as_ref().unwrap();
    assert!((done.interface.interface_area() - 686.0).abs() < 5.0);
    // the table lists the antibody's residues first (the first side), and scrolls: the antigen's come below what is drawn
    assert!(done
        .rows
        .iter()
        .any(|r| (r.chain, r.name.as_str(), r.number) == ('C', "GLN", 121)));
    let text = screen.text();
    assert!(
        text.contains("B TYR 101"),
        "the antibody residue that buries most"
    );
    assert!(text.contains("39 contacts"));
}
