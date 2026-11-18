use crate::pipeline::steps::Aligned;
use crate::ui::app::{App, Startup};
use crate::ui::platform::Workspace;
use crate::ui::testing::Screen;

const UBIQUITIN: &str =
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
/// So many columns stand before ubiquitin in the rows of the wide family, and after it.
const FLANK: usize = 150;

/// Ubiquitin and two relatives with a long stretch of something else before them, as a domain has in its protein: the
/// residue `j` of the structure is in column `FLANK + j`. With all 76 residues the same stretch follows them too; `upto`
/// below 76 ends the rows there, so that the last residues of the structure are in no row.
fn wide_family(upto: usize) -> Aligned {
    let flank = "GS".repeat(FLANK / 2);
    let after = if upto == UBIQUITIN.len() {
        flank.as_str()
    } else {
        ""
    };
    let row = |middle: &str| format!("{flank}{}{after}", &middle[..upto]).into_bytes();
    Aligned {
        rows: vec![
            ("ubiquitin".to_string(), row(UBIQUITIN)),
            (
                "relative 1".to_string(),
                row(&UBIQUITIN.replace("EVEPS", "EVWPS")),
            ),
            (
                "relative 2".to_string(),
                row(&UBIQUITIN.replace("KTITL", "KTWTL")),
            ),
        ],
        engine: "test".into(),
        seconds: 0.0,
    }
}

/// The example structure (ubiquitin, residues numbered from 1) and the wide family, not yet lined up.
fn app_with(family: Aligned) -> App {
    let mut app = App::new(Startup::default());
    app.platform.structure.load_example();
    app.platform.alignment.set_result(family);
    app
}

fn picked(app: &App) -> Option<usize> {
    app.platform.structure.picked()
}

fn selected(app: &App) -> Option<usize> {
    app.platform.alignment.selected_column()
}

#[test]
fn a_column_picks_its_residue_once_the_structure_is_lined_up_with_the_alignment() {
    let mut app = app_with(wide_family(76));
    // not lined up yet: the column is selected, and the structure is left as it is
    app.platform.select_column(Some(FLANK + 20));
    assert_eq!(selected(&app), Some(FLANK + 20));
    assert_eq!(picked(&app), None);
    assert!(app
        .platform
        .structure
        .residues_in_column(FLANK + 20)
        .is_empty());

    // lining them up shows the column that is selected already
    app.colour_structure_by_alignment();
    assert_eq!(picked(&app), Some(20));
    assert_eq!(selected(&app), Some(FLANK + 20));
    assert_eq!(
        app.platform.structure.residues_in_column(FLANK + 20),
        ["D 21 of chain A"]
    );

    // from now on each column picks its residue, a column of the flank picks none, and no column picks none
    app.platform.select_column(Some(FLANK + 30));
    assert_eq!(picked(&app), Some(30));
    app.platform.select_column(Some(3));
    assert_eq!(picked(&app), None);
    app.platform.select_column(Some(FLANK));
    assert_eq!(picked(&app), Some(0));
    app.platform.select_column(None);
    assert_eq!((picked(&app), selected(&app)), (None, None));
}

#[test]
fn a_residue_picked_in_the_structure_selects_its_column_and_one_without_a_column_selects_none() {
    // the rows end six residues before the structure does
    let mut app = app_with(wide_family(70));
    app.colour_structure_by_alignment();
    app.platform.structure.pick(Some(40));
    app.platform.residue_picked();
    assert_eq!(selected(&app), Some(FLANK + 40));
    assert!(
        app.platform.alignment.reveal,
        "the page scrolls to it when it is next shown"
    );
    // a residue that no row has
    app.platform.structure.pick(Some(75));
    app.platform.residue_picked();
    assert_eq!(selected(&app), None);
    // and letting go of the residue lets go of the column
    app.platform.structure.pick(Some(12));
    app.platform.residue_picked();
    assert_eq!(selected(&app), Some(FLANK + 12));
    app.platform.structure.pick(None);
    app.platform.residue_picked();
    assert_eq!(selected(&app), None);
}

#[test]
fn a_structure_that_was_not_lined_up_leaves_the_column_alone() {
    let mut app = app_with(wide_family(76));
    app.platform.select_column(Some(5));
    app.platform.structure.pick(Some(3));
    app.platform.residue_picked();
    assert_eq!(selected(&app), Some(5));
    assert!(!app.platform.alignment.reveal);
    app.platform.structure.pick(None);
    app.platform.residue_picked();
    assert_eq!(selected(&app), Some(5));
}

#[test]
fn lining_up_selects_the_column_of_a_residue_that_is_picked_when_no_column_is() {
    let mut app = app_with(wide_family(76));
    app.platform.structure.pick(Some(33));
    app.colour_structure_by_alignment();
    assert_eq!(selected(&app), Some(FLANK + 33));
    assert_eq!(picked(&app), Some(33));
    // an alignment of another protein lines nothing up, and changes neither
    let mut other = app_with(Aligned {
        rows: vec![
            (
                "a".to_string(),
                b"MSDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEY".to_vec(),
            ),
            (
                "b".to_string(),
                b"MSDWIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEY".to_vec(),
            ),
        ],
        engine: "test".into(),
        seconds: 0.0,
    });
    other.platform.select_column(Some(4));
    other.platform.structure.pick(Some(9));
    other.colour_structure_by_alignment();
    assert_eq!((selected(&other), picked(&other)), (Some(4), Some(9)));
}

#[test]
fn the_link_goes_when_the_alignment_or_its_score_changes() {
    let ctx = eframe::egui::Context::default();
    // another score of the same alignment takes the colours off the structure, and the columns go with them
    let mut scored = app_with(wide_family(76));
    scored.colour_structure_by_alignment();
    scored.platform.select_column(Some(FLANK + 7));
    assert_eq!(picked(&scored), Some(7));
    scored
        .platform
        .alignment
        .score_for_test(crate::ui::platform::alignment::Score::Entropy);
    scored.platform.alignment.refresh_scores();
    scored.structure_tick(&ctx);
    assert!(scored
        .platform
        .structure
        .residues_in_column(FLANK + 7)
        .is_empty());
    // the column is a column of the alignment still, and no longer moves the residue
    assert_eq!(selected(&scored), Some(FLANK + 7));
    scored.platform.select_column(Some(FLANK + 9));
    assert_eq!(picked(&scored), Some(7));

    let mut app = app_with(wide_family(76));
    app.colour_structure_by_alignment();
    app.platform.select_column(Some(FLANK + 7));
    assert_eq!(picked(&app), Some(7));
    // another alignment: its columns are other columns
    app.platform.alignment.set_result(wide_family(70));
    assert_eq!(selected(&app), None);
    app.structure_tick(&ctx);
    assert!(app
        .platform
        .structure
        .residues_in_column(FLANK + 7)
        .is_empty());
    // the residue stays picked (it is a residue of the structure still), and a column no longer moves it
    assert_eq!(picked(&app), Some(7));
    app.platform.select_column(Some(FLANK + 9));
    assert_eq!(picked(&app), Some(7));
    // nor does the residue move the column
    app.platform.structure.pick(Some(2));
    app.platform.residue_picked();
    assert_eq!(selected(&app), Some(FLANK + 9));
}

#[test]
fn the_pictured_example_shows_one_column_on_both_pages() {
    for on_structure in [false, true] {
        // as the command line names the two pictures
        let view = if on_structure {
            "linked-residue"
        } else {
            "linked-column"
        };
        let mut app = App::new(Startup {
            view: view.to_string(),
            ..Startup::default()
        });
        let column = selected(&app).expect("a column is selected");
        // no row has a gap, so the residue is the one with the column's number
        assert_eq!(picked(&app), Some(column));
        let want = if on_structure {
            Workspace::Structure
        } else {
            Workspace::Alignment
        };
        assert_eq!(app.workspace, want);
        // the made-up rows say that they are made up, and the page says where the alignment came from
        let rows = &app.platform.alignment.result.as_ref().expect("rows").rows;
        assert_eq!(rows[0].0, "ubiquitin (1UBQ)");
        assert!(rows[1..]
            .iter()
            .all(|(name, _)| name.starts_with("made-up variant ")));
        assert!(
            app.platform.alignment.made.is_none(),
            "nothing here made it"
        );
        // every variant is ubiquitin with some residues changed, never all of them
        for (name, row) in &rows[1..] {
            let changed = row.iter().zip(&rows[0].1).filter(|(a, b)| a != b).count();
            assert!((1..=8).contains(&changed), "{name}: {changed} changed");
            assert_eq!(row.len(), rows[0].1.len());
        }
        // and the name the structure page gives the residue is on the alignment's link
        let mut screen = Screen::new();
        app.workspace = Workspace::Alignment;
        screen.settle(&mut app);
        let residues = app.platform.structure.residues_in_column(column);
        assert_eq!(residues.len(), 1, "{residues:?}");
        let link = format!("Show {} in 3D", residues[0]);
        assert!(screen.text().contains(&link), "{}", screen.text());
    }
}

/// The example structure twice: chain A, and a copy of it called chain B (in the same place: only its sequence matters).
fn two_copies() -> String {
    let one = include_str!("../../../../assets/structure/1ubq.pdb");
    let atoms: Vec<&str> = one.lines().filter(|l| l.starts_with("ATOM")).collect();
    let copy = atoms.iter().map(|line| {
        // the chain is the twenty-second character of an ATOM line
        format!("{}B{}", &line[..21], &line[22..])
    });
    atoms
        .iter()
        .map(|line| line.to_string())
        .chain(copy)
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn the_copies_of_a_chain_share_a_column_and_the_one_that_is_picked_is_named_first() {
    let mut app = App::new(Startup::default());
    app.platform
        .structure
        .show_pdb("two copies", &two_copies())
        .expect("a structure");
    app.platform.alignment.set_result(wide_family(76));
    app.colour_structure_by_alignment();
    let named = |app: &App| app.platform.structure.residues_in_column(FLANK + 20);
    // a column picks the residue of the first chain, and both are named
    app.platform.select_column(Some(FLANK + 20));
    assert_eq!(picked(&app), Some(20));
    assert_eq!(named(&app), ["D 21 of chain A", "D 21 of chain B"]);
    // the same residue of the second chain is in the same column: the column stays, the chain that was clicked stays
    // picked, and it is the one the Alignment page names
    app.platform.structure.pick(Some(76 + 20));
    app.platform.residue_picked();
    assert_eq!(selected(&app), Some(FLANK + 20));
    assert_eq!(picked(&app), Some(76 + 20));
    assert_eq!(named(&app), ["D 21 of chain B", "D 21 of chain A"]);
}

/// The lowest text `label` on the screen (the sequence strip is under the structure).
fn lowest(screen: &Screen, label: &str) -> eframe::egui::Pos2 {
    screen
        .find_all(label)
        .into_iter()
        .max_by(|a, b| a.y.total_cmp(&b.y))
        .unwrap_or_else(|| panic!("no \"{label}\" on the screen:\n{}", screen.text()))
}

#[test]
fn a_click_on_a_residue_of_the_strip_selects_its_column_and_the_alignment_scrolls_to_it() {
    let mut app = app_with(wide_family(76));
    app.colour_structure_by_alignment();
    app.workspace = Workspace::Structure;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    // the strip writes the number of every tenth residue over its letter: the letter under "20" is residue 20, in the
    // row of letters that begins with the M of the first residue
    let number = lowest(&screen, "20");
    let letters = lowest(&screen, "M").y;
    assert!(letters > number.y, "the letters are under the numbers");
    screen.click_at(&mut app, eframe::egui::pos2(number.x, letters));
    screen.settle(&mut app);
    assert_eq!(picked(&app), Some(19));
    assert_eq!(selected(&app), Some(FLANK + 19));
    // the structure says which column it is, in the picture and as a link
    let text = screen.text();
    assert!(text.contains(" · column 170"), "{text}");
    let link = "S 20 is column 170 of the alignment";
    assert!(text.contains(link), "{text}");

    // the link goes to the Alignment page, which has scrolled to the column: without that, column 170 is far to the right
    // of the window and its number is not drawn
    screen.click(&mut app, link);
    screen.settle(&mut app);
    assert_eq!(app.workspace, Workspace::Alignment);
    let ruler = screen.find_all("170");
    assert_eq!(
        ruler.len(),
        1,
        "the ruler's 170 is in view: {}",
        screen.text()
    );
    // away from both edges of the letters (the names end about 130 points right of the middle of a name)
    let name = screen.find("ubiquitin").expect("the names");
    assert!(
        ruler[0].x > name.x + 250.0 && ruler[0].x < 1280.0 - 250.0,
        "{ruler:?}"
    );
    assert!(!app.platform.alignment.reveal, "it scrolls once");
    let text = screen.text();
    assert!(text.contains("Column 170 of 376"), "{text}");

    // another residue picked in the structure brings its column to the very same place: the page scrolls to the column
    app.platform.structure.pick(Some(49));
    app.platform.residue_picked();
    screen.settle(&mut app);
    let again = screen.find_all("200");
    assert_eq!(again.len(), 1, "the ruler's 200: {}", screen.text());
    assert!(
        (again[0].x - ruler[0].x).abs() < 1.0,
        "{again:?} where {ruler:?} was"
    );

    // scrolled as it is, a click on a letter selects the column under the pointer (ten columns left of the ruler's 200
    // is column 190), and a click on a name selects none
    let before = screen.find_all("190");
    assert_eq!(before.len(), 1, "{before:?}");
    let cell = (again[0].x - before[0].x) / 10.0;
    screen.click_at(
        &mut app,
        eframe::egui::pos2(before[0].x - cell / 2.0, name.y),
    );
    screen.settle(&mut app);
    assert_eq!(selected(&app), Some(189));
    assert_eq!(picked(&app), Some(39));
    screen.click_at(&mut app, name);
    screen.settle(&mut app);
    assert_eq!(selected(&app), Some(189), "a name is not a column");

    // the page names the residue, and that goes back to the structure
    let text = screen.text();
    let back = "Show Q 40 of chain A in 3D";
    assert!(text.contains(back), "{text}");
    screen.click(&mut app, back);
    screen.settle(&mut app);
    assert_eq!(app.workspace, Workspace::Structure);
    assert_eq!(picked(&app), Some(39));
}

#[test]
fn a_click_on_a_residue_in_the_picture_selects_its_column_and_a_click_beside_the_structure_lets_both_go(
) {
    let mut app = app_with(wide_family(76));
    app.colour_structure_by_alignment();
    app.platform.select_column(Some(FLANK + 30));
    app.workspace = Workspace::Structure;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    // the picked residue has a ring round its alpha carbon, in the picture (which is right of the options)
    let rings = |screen: &Screen| -> Vec<eframe::egui::Pos2> {
        screen
            .circles_of_radius(9.0)
            .into_iter()
            .filter(|at| at.x > 420.0)
            .collect()
    };
    let ring = rings(&screen);
    assert_eq!(ring.len(), 1, "{ring:?}");
    // let go of it, then click where it was
    app.platform.select_column(None);
    screen.settle(&mut app);
    assert!(rings(&screen).is_empty());
    screen.click_at(&mut app, ring[0]);
    screen.settle(&mut app);
    assert_eq!(picked(&app), Some(30));
    assert_eq!(selected(&app), Some(FLANK + 30));
    assert!(app.platform.alignment.reveal);
    assert_eq!(rings(&screen), ring);
    // the words in the corner of the picture are beside the structure: a click there lets go of the residue and its column
    let beside = screen
        .find(super::super::viewer3d::HINT)
        .expect("the hint of the picture");
    screen.click_at(&mut app, beside);
    screen.settle(&mut app);
    assert_eq!((picked(&app), selected(&app)), (None, None));
}

#[test]
fn beside_the_link_to_the_residue_deselect_still_lets_the_column_go_and_goes_nowhere() {
    let mut app = app_with(wide_family(76));
    app.colour_structure_by_alignment();
    app.platform.select_column(Some(FLANK + 5));
    // the page is drawn with the column scrolled to, as if it had been selected from the structure
    app.platform.structure.pick(Some(5));
    app.platform.residue_picked();
    app.workspace = Workspace::Alignment;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    let text = screen.text();
    assert!(text.contains("Show K 6 of chain A in 3D"), "{text}");
    screen.click(&mut app, "Deselect");
    screen.settle(&mut app);
    assert_eq!((selected(&app), picked(&app)), (None, None));
    assert_eq!(app.workspace, Workspace::Alignment);
}

/// The structure page with ubiquitin not yet lined up, a column selected in the alignment and another residue picked.
fn before_colouring(column: usize) -> (App, Screen) {
    let mut app = app_with(wide_family(76));
    app.platform.select_column(Some(column));
    app.platform.structure.pick(Some(9));
    app.workspace = Workspace::Structure;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    (app, screen)
}

#[test]
fn colouring_keeps_the_column_that_is_selected_whatever_residue_was_picked() {
    // a column of the flank, which no residue of the structure is in: the residue that was picked is of another column and
    // is let go, and the column stays
    let (mut app, mut screen) = before_colouring(3);
    screen.click(&mut app, "Colour by conservation");
    screen.settle(&mut app);
    assert!(
        app.platform.structure.column_of_pick().is_some(),
        "lined up"
    );
    assert_eq!((selected(&app), picked(&app)), (Some(3), None));
    assert!(
        !app.platform.alignment.reveal,
        "nothing was selected from this page"
    );
    // a column of ubiquitin: its residue is picked instead of the one that was
    let (mut app, mut screen) = before_colouring(FLANK + 12);
    screen.click(&mut app, "Colour by conservation");
    screen.settle(&mut app);
    assert_eq!((selected(&app), picked(&app)), (Some(FLANK + 12), Some(12)));
    assert!(!app.platform.alignment.reveal);
}

#[test]
fn a_click_on_a_column_shows_in_the_structure_and_resetting_the_view_lets_both_go() {
    let mut app = app_with(wide_family(76));
    app.colour_structure_by_alignment();
    // a residue is picked through its column
    app.platform.select_column(Some(FLANK + 5));
    assert_eq!(picked(&app), Some(5));
    app.workspace = Workspace::Alignment;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    // the window opens at the first columns, which are the flank: no residue of the structure is in them, so a click on
    // one lets go of the residue that was picked
    let tenth = screen.find_all("10");
    assert_eq!(tenth.len(), 1, "{tenth:?}");
    let row = screen.find("ubiquitin").expect("a row").y;
    screen.click_at(&mut app, eframe::egui::pos2(tenth[0].x - 3.0, row));
    screen.settle(&mut app);
    assert_eq!(selected(&app), Some(9));
    assert_eq!(picked(&app), None);
    let text = screen.text();
    assert!(text.contains("Column 10 of 376"), "{text}");
    assert!(!text.contains("in 3D"), "{text}");

    // a column of ubiquitin, selected as a click would
    app.platform.select_column(Some(FLANK + 5));
    app.workspace = Workspace::Structure;
    screen.settle(&mut app);
    assert_eq!(picked(&app), Some(5));
    assert!(screen.text().contains(" · column 156"), "{}", screen.text());
    // Reset view lets go of the residue, and the column goes with it
    screen.click(&mut app, "Reset view");
    screen.settle(&mut app);
    assert_eq!((picked(&app), selected(&app)), (None, None));
}
