use super::*;
use crate::pipeline::steps::Aligned;
use crate::ui::app::{App, Startup};
use crate::ui::platform::alignment::Score;
use crate::ui::platform::Workspace;
use crate::ui::testing::Screen;

fn aligned(rows: &[(&str, &str)]) -> Aligned {
    Aligned {
        rows: rows
            .iter()
            .map(|(name, row)| (name.to_string(), row.as_bytes().to_vec()))
            .collect(),
        engine: "test".into(),
        seconds: 0.0,
    }
}

/// Five proteins of ten columns: the first four are alike in every row, the fifth is K in four rows and R in one, the
/// seventh is a gap in four rows, and the last is one of three letters.
fn small() -> Aligned {
    aligned(&[
        ("one", "MKLVKA-AAG"),
        ("two", "MKLVKAWAAG"),
        ("three", "MKLVRA-AAS"),
        ("four", "MKLVKA-AAT"),
        ("five", "MKLVKA-AAT"),
    ])
}

fn page_with(alignment: Aligned) -> AlignmentPage {
    let mut page = AlignmentPage::default();
    page.set_result(alignment);
    page
}

#[test]
fn the_facts_of_a_column_are_its_letters_its_gaps_its_score_and_its_block() {
    let page = page_with(small());
    let facts = page.column_facts(4).expect("the fifth column");
    assert_eq!(
        (facts.column, facts.of, facts.consensus, facts.gaps),
        (4, 10, b'K', 0)
    );
    assert_eq!(facts.symbols, [(b'K', 4), (b'R', 1)]);
    // the score that the bars show: identity, the share of the rows with the commonest letter
    assert_eq!(facts.score_name, "Identity");
    assert!((facts.score.expect("a score") - 0.8).abs() < 1e-6);
    assert!(!facts.conserved);
    assert_eq!(facts.block, None);
    assert_eq!(
        facts.words(),
        "Column 5 of 10 · consensus K · K 4, R 1 · no gaps · Identity 0.80"
    );
    // a column of the block of four that are alike in every row
    let first = page.column_facts(0).expect("the first column");
    assert!(first.conserved);
    assert_eq!(first.block, Some((0, 3)));
    // the last column of the block is in it, and the one after is not
    assert_eq!(
        page.column_facts(3).expect("the fourth").block,
        Some((0, 3))
    );
    assert_eq!(page.column_facts(4).expect("the fifth").block, None);
    assert_eq!(
        first.words(),
        "Column 1 of 10 · consensus M · M 5 · no gaps · Identity 1.00, among the conserved · in the conserved block 1-4"
    );
    // mostly gaps: the consensus is a gap, and the gaps are counted
    let gapped = page.column_facts(6).expect("the seventh column");
    assert_eq!(
        gapped.words(),
        "Column 7 of 10 · consensus a gap · W 1 · 4 gaps · Identity 0.20"
    );
    // a column that is a gap in every row has no letters to name
    let empty = page_with(aligned(&[("a", "A-C"), ("b", "A-C"), ("c", "A-C")]));
    assert_eq!(
        empty.column_facts(1).expect("a column").words(),
        "Column 2 of 3 · consensus a gap · 3 gaps · Identity 0.00"
    );
    // one gap is one gap, and there is no column past the end
    let one = page_with(aligned(&[("a", "A-"), ("b", "AC"), ("c", "AC")]));
    assert!(one
        .column_facts(1)
        .expect("a column")
        .words()
        .contains(" · 1 gap · "));
    assert!(page.column_facts(10).is_none());
    assert!(AlignmentPage::default().column_facts(0).is_none());
}

#[test]
fn a_column_that_has_no_score_says_so_and_many_letters_are_cut_short() {
    let mut page = page_with(small());
    // the scores made for many sequences give none to a column that is mostly gaps
    page.score_for_test(Score::Jsd);
    page.refresh_scores();
    let gapped = page.column_facts(6).expect("the seventh column");
    assert_eq!(gapped.score, None);
    assert_eq!(gapped.score_name, "Jensen-Shannon divergence");
    assert!(
        gapped.words().ends_with(" · not scored (too many gaps)"),
        "{}",
        gapped.words()
    );
    // eight letters in one column: the six commonest are named and the others counted
    let many = page_with(aligned(&[
        ("a", "A"),
        ("b", "A"),
        ("c", "C"),
        ("d", "D"),
        ("e", "E"),
        ("f", "F"),
        ("g", "G"),
        ("h", "H"),
        ("i", "I"),
    ]));
    let words = many.column_facts(0).expect("a column").words();
    assert!(
        words.contains(" · A 2, C 1, D 1, E 1, F 1, G 1 and 2 more · "),
        "{words}"
    );
}

#[test]
fn a_column_is_selected_only_if_the_alignment_has_it_and_a_new_alignment_forgets_it() {
    let mut page = page_with(small());
    assert_eq!(page.selected_column(), None);
    page.select_column(Some(9), false);
    assert_eq!(page.selected_column(), Some(9));
    page.select_column(Some(10), false);
    assert_eq!(page.selected_column(), None, "there is no eleventh column");
    // the page scrolls to a column that was selected elsewhere: not to one that is not there, nor to none
    page.select_column(Some(10), true);
    assert!(!page.reveal);
    page.select_column(None, true);
    assert!(!page.reveal);
    page.select_column(Some(3), true);
    assert!(page.reveal);
    page.set_result(small());
    assert_eq!(page.selected_column(), None);
    assert!(!page.reveal, "nothing is left to scroll to");
    // and without an alignment there is no column to select
    let mut empty = AlignmentPage::default();
    empty.select_column(Some(0), false);
    assert_eq!(empty.selected_column(), None);
}

/// The window on the Alignment page, with ubiquitin and two relatives aligned (76 columns, no gaps).
fn window() -> (App, Screen) {
    let ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
    let relative = ubiquitin.replace("EVEPS", "EVWPS");
    let mut app = App::new(Startup::default());
    app.platform.alignment.set_result(aligned(&[
        ("ubiquitin", ubiquitin),
        ("relative 1", &relative),
        ("relative 2", ubiquitin),
    ]));
    app.workspace = Workspace::Alignment;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    (app, screen)
}

/// Where the letter of column `k` (counted from 0; `k + 1` must be a number of the ruler) is in the row called `row`,
/// and how wide a column is.
fn letter_at(screen: &Screen, k: usize, row: &str) -> (eframe::egui::Pos2, f32) {
    let tick = |n: usize| {
        let found = screen.find_all(&n.to_string());
        assert_eq!(found.len(), 1, "the ruler's {n}: {found:?}");
        found[0].x
    };
    // the ruler writes a number over the right edge of every tenth column
    let cell = (tick(20) - tick(10)) / 10.0;
    let y = screen.find(row).expect("the name of the row").y;
    (eframe::egui::pos2(tick(k + 1) - cell / 2.0, y), cell)
}

/// The marks round a column: what is outlined in the accent and taller than a row (an icon of the rail has a small
/// outline in that colour too).
fn marks(screen: &Screen, app: &App) -> Vec<eframe::egui::Rect> {
    screen
        .rects_outlined(app.prefs.look().accent)
        .into_iter()
        .filter(|r| r.height() > 30.0)
        .collect()
}

/// The tints one column wide (`cell`): what stands behind the top of the selected column.
fn tints(screen: &Screen, app: &App, cell: f32) -> Vec<eframe::egui::Rect> {
    screen
        .rects_filled(app.prefs.look().accent_soft)
        .into_iter()
        .filter(|r| (r.width() - cell).abs() < 0.01)
        .collect()
}

#[test]
fn the_link_to_the_structure_names_the_first_residue_and_counts_the_other_chains() {
    let residue = |chain: char| format!("L 56 of chain {chain}");
    assert_eq!(link_words(&[]), None);
    assert_eq!(
        link_words(&[residue('A')]).as_deref(),
        Some("Show L 56 of chain A in 3D")
    );
    assert_eq!(
        link_words(&[residue('A'), residue('B')]).as_deref(),
        Some("Show L 56 of chain A in 3D (and 1 more chain)")
    );
    assert_eq!(
        link_words(&[residue('A'), residue('B'), residue('C')]).as_deref(),
        Some("Show L 56 of chain A in 3D (and 2 more chains)")
    );
}

#[test]
fn a_click_on_a_letter_selects_its_column_and_the_page_says_what_it_holds() {
    let (mut app, mut screen) = window();
    assert!(screen.text().contains(HINT), "{}", screen.text());
    assert!(marks(&screen, &app).is_empty());

    // the letter of the twentieth column in the second row
    let (at, cell) = letter_at(&screen, 19, "relative 1");
    screen.click_at(&mut app, at);
    screen.settle(&mut app);
    assert_eq!(app.platform.alignment.selected_column(), Some(19));
    let text = screen.text();
    assert!(
        text.contains("Column 20 of 76 · consensus S · S 3 · no gaps · Identity 1.00"),
        "{text}"
    );
    assert!(!text.contains(HINT));
    // the column is outlined where its letters are: one column wide, round the place that was clicked, over all three rows
    let outlined = marks(&screen, &app);
    assert_eq!(outlined.len(), 1, "{outlined:?}");
    let mark = outlined[0];
    assert!(
        (mark.width() - cell).abs() < 0.01,
        "{mark:?} for a column of {cell}"
    );
    assert!(
        mark.left() < at.x && at.x < mark.right(),
        "{mark:?} round {at:?}"
    );
    let top = screen.find("ubiquitin").expect("the first row").y;
    let bottom = screen.find("relative 2").expect("the last row").y;
    assert!(mark.top() < top && bottom < mark.bottom(), "{mark:?}");
    // and above the rows its number, bar and consensus letter stand on a tint, so that it is found from the top too
    let tinted = tints(&screen, &app, cell);
    assert_eq!(tinted.len(), 1, "{tinted:?}");
    assert!(
        tinted[0].left() < at.x && at.x < tinted[0].right(),
        "{tinted:?}"
    );
    assert!(
        tinted[0].bottom() < top && tinted[0].bottom() >= mark.top() - 0.01,
        "{tinted:?}"
    );

    // another column takes its place
    let (other, _) = letter_at(&screen, 29, "ubiquitin");
    screen.click_at(&mut app, other);
    screen.settle(&mut app);
    assert_eq!(app.platform.alignment.selected_column(), Some(29));
    assert_eq!(marks(&screen, &app).len(), 1);

    // a click on the column that is selected lets it go
    screen.click_at(&mut app, other);
    screen.settle(&mut app);
    assert_eq!(app.platform.alignment.selected_column(), None);
    assert!(screen.text().contains(HINT));
    assert!(marks(&screen, &app).is_empty());
    assert!(tints(&screen, &app, cell).is_empty());
}

#[test]
fn the_ruler_selects_a_column_too_and_a_name_or_the_room_past_the_last_column_selects_none() {
    let (mut app, mut screen) = window();
    // the consensus letter of the tenth column, above the rows
    let (letter, cell) = letter_at(&screen, 9, "ubiquitin");
    let ruler = screen.find_all("10")[0];
    screen.click_at(&mut app, eframe::egui::pos2(letter.x, ruler.y));
    screen.settle(&mut app);
    assert_eq!(app.platform.alignment.selected_column(), Some(9));
    // the name of a row
    let name = screen.find("relative 2").expect("a name");
    screen.click_at(&mut app, name);
    screen.settle(&mut app);
    assert_eq!(
        app.platform.alignment.selected_column(),
        Some(9),
        "a click on a name leaves the column"
    );
    // past the seventy-sixth column
    let (seventieth, _) = letter_at(&screen, 69, "ubiquitin");
    let past = eframe::egui::pos2(seventieth.x + 7.0 * cell, seventieth.y);
    screen.click_at(&mut app, past);
    screen.settle(&mut app);
    assert_eq!(app.platform.alignment.selected_column(), Some(9));
    // the last column is still a column
    let last = eframe::egui::pos2(seventieth.x + 6.0 * cell, seventieth.y);
    screen.click_at(&mut app, last);
    screen.settle(&mut app);
    assert_eq!(app.platform.alignment.selected_column(), Some(75));
    // the spare room under the last row is under a column, and is not the column: a row is 20 points high, and its name
    // is in the middle of it
    app.platform.select_column(None);
    screen.settle(&mut app);
    let under = screen.find("relative 2").expect("the last row").y + 16.0;
    screen.click_at(&mut app, eframe::egui::pos2(letter.x, under));
    screen.settle(&mut app);
    assert_eq!(app.platform.alignment.selected_column(), None);
    // (the last row itself is the column)
    screen.click_at(&mut app, eframe::egui::pos2(letter.x, under - 8.0));
    screen.settle(&mut app);
    assert_eq!(app.platform.alignment.selected_column(), Some(9));
}

#[test]
fn clear_lets_the_column_go() {
    let (mut app, mut screen) = window();
    let (at, _) = letter_at(&screen, 39, "ubiquitin");
    screen.click_at(&mut app, at);
    screen.settle(&mut app);
    assert_eq!(app.platform.alignment.selected_column(), Some(39));
    screen.click(&mut app, CLEAR);
    screen.settle(&mut app);
    assert_eq!(app.platform.alignment.selected_column(), None);
    assert!(screen.find(CLEAR).is_none(), "nothing is left to clear");
    assert_eq!(app.workspace, Workspace::Alignment, "it goes nowhere");
}

#[test]
fn a_figure_of_the_alignment_and_the_report_do_not_scroll_and_the_report_marks_no_column() {
    let mut page = page_with(small());
    page.select_column(Some(2), true);
    let look = crate::ui::design::look(false);
    let view = page.view(&look).expect("a view for the report");
    assert_eq!(view.selected, None);
    assert!(!view.reveal);
}
