use super::*;

fn orf(start: usize) -> Orf {
    Orf {
        start,
        length: 30,
        reverse: false,
        frame: 0,
        open_end: false,
    }
}

#[test]
fn orf_features_wrap_around_the_origin() {
    let o = Orf {
        start: 90,
        length: 30,
        reverse: true,
        frame: 1,
        open_end: false,
    };
    let f = orf_feature(&o, 100, 2);
    assert_eq!(f.spans, vec![(0, 20), (90, 100)]);
    assert!(f.reverse);
    assert_eq!(f.name, "ORF 2");
    assert_eq!(orf_feature(&orf(10), 100, 1).spans, vec![(10, 40)]);
}

#[test]
fn min_choices_match_the_original_menu() {
    assert_eq!(
        ORF_MIN_CHOICES,
        [0, 30, 60, 90, 180, 300, 600, 900, 1200, 1500]
    );
}

#[test]
fn the_bases_of_an_orf_end_where_the_sequence_ends() {
    assert_eq!(orf_bases(&orf(10), 100), (10, 40));
    assert_eq!(orf_bases(&orf(70), 100), (70, 100));
    // one that runs past the end of a circular sequence: the part before the end
    assert_eq!(orf_bases(&orf(90), 100), (90, 100));
}

#[test]
fn the_heading_names_the_orf_whose_bases_are_exactly_the_selected_ones() {
    let reverse = Orf {
        reverse: true,
        ..orf(50)
    };
    let orfs = [orf(10), reverse.clone(), orf(1200)];
    let words = |selection| selected_words(&orfs, 5000, selection);
    assert_eq!(
        words(Some((10, 40))).as_deref(),
        Some("the selection is the forward ORF at 11–40")
    );
    assert_eq!(
        words(Some((50, 80))).as_deref(),
        Some("the selection is the reverse ORF at 51–80")
    );
    assert_eq!(
        words(Some((1200, 1230))).as_deref(),
        Some("the selection is the forward ORF at 1,201–1,230")
    );
    // a base more or less is not the ORF, and nothing selected is none
    assert_eq!(words(Some((10, 41))), None);
    assert_eq!(words(Some((11, 40))), None);
    assert_eq!(words(None), None);
    // the same bases on both strands
    let both = [orf(50), reverse];
    assert_eq!(
        selected_words(&both, 5000, Some((50, 80))).as_deref(),
        Some("the selection is 2 of these ORFs, at 51–80")
    );
    // one that ends on the last base of the sequence is all there: no part of it is beyond the end
    assert_eq!(
        selected_words(&[orf(70)], 100, Some((70, 100))).as_deref(),
        Some("the selection is the forward ORF at 71–100")
    );
    // one that runs past the end is the selection of its part before the end, as a click on it selects
    assert_eq!(
        selected_words(&[orf(90)], 100, Some((90, 100))).as_deref(),
        Some(
            "the selection is the forward ORF at 91–100 (its part before the end of the sequence)"
        )
    );
}

#[test]
fn the_last_base_of_an_orf_that_ends_the_sequence_is_the_last_of_the_sequence() {
    assert_eq!(end_shown(&orf(10), 100), 40);
    // from base 71 for thirty bases: it ends on base 100, which round a circle is counted 0
    assert_eq!(end_shown(&orf(70), 100), 100);
    // and one that runs past the end ends twenty bases after the start
    assert_eq!(end_shown(&orf(90), 100), 20);
}

/// Where the row of the list's one ORF says "forward": the last such word that is drawn, since the window is drawn after the
/// page behind it (whose inspector says "forward" of selected bases too).
fn row(screen: &crate::ui::testing::Screen) -> eframe::egui::Pos2 {
    *screen
        .find_all("forward")
        .last()
        .expect("the row of the ORF")
}

/// Whether the row is marked: the fill of a chosen button stands behind its word, the one that stands behind "Top" among
/// the strands (which is the strand the list of these tests is of).
fn row_is_marked(screen: &crate::ui::testing::Screen) -> bool {
    let chosen = screen
        .fills_at(screen.find("Top").expect("the strand that is chosen"))
        .last()
        .expect("a chosen button has a fill")
        .1;
    screen
        .fills_at(row(screen))
        .iter()
        .any(|(_, fill)| *fill == chosen)
}

#[test]
fn the_list_marks_the_orf_that_is_selected_wherever_it_was_selected() {
    use crate::ui::testing::Screen;
    let mut app = app_with(false, false);
    app.open_orf_list();
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert!(
        screen.text().contains("1 ORF in 600 bp"),
        "{}",
        screen.text()
    );
    assert!(!screen.text().contains("the selection is"));
    assert!(!row_is_marked(&screen));
    assert_eq!(app.orf_win.scrolled_to, None);

    // selected somewhere else, as Find Next ORF or a drag across exactly its bases does
    app.docs[0].select(ORF_START, ORF_START + 126);
    screen.settle(&mut app);
    let text = screen.text();
    assert!(
        text.contains("1 ORF in 600 bp · the selection is the forward ORF at 101–226"),
        "{text}"
    );
    assert!(row_is_marked(&screen));
    let id = app.docs[0].id;
    assert_eq!(
        app.orf_win.scrolled_to,
        Some((id, (ORF_START, ORF_START + 126)))
    );

    // other bases are no ORF: the heading and the row say so again
    app.docs[0].select(ORF_START, ORF_START + 125);
    screen.settle(&mut app);
    assert!(!screen.text().contains("the selection is"));
    assert!(!row_is_marked(&screen));
    assert_eq!(app.orf_win.scrolled_to, None);

    // and a click on the row selects the ORF, and the heading names it
    // the window takes some frames to come to its size: a click needs the row to be where it was looked for
    let mut at = row(&screen);
    for _ in 0..60 {
        screen.frame(&mut app, vec![]);
        let now = row(&screen);
        if now == at {
            break;
        }
        at = now;
    }
    screen.click_at(&mut app, at);
    screen.settle(&mut app);
    assert_eq!(app.docs[0].selection(), Some((ORF_START, ORF_START + 126)));
    assert!(screen
        .text()
        .contains("the selection is the forward ORF at 101–226"));
}

/// So many frames a window takes to fade in, and the list to scroll to a row (egui does both over a third of a second
/// at most, and a frame of these tests is a sixtieth).
const A_WHILE: usize = 30;

fn a_while(screen: &mut crate::ui::testing::Screen, app: &mut App) {
    for _ in 0..A_WHILE {
        screen.frame(app, vec![]);
    }
}

/// Draw frames until the window is there in its own colours and has come to its size (the button under the list stays
/// where it is).
fn steady(screen: &mut crate::ui::testing::Screen, app: &mut App) {
    a_while(screen, app);
    let button = |screen: &crate::ui::testing::Screen| screen.find("Add these ORFs as features");
    let mut at = button(screen);
    for _ in 0..60 {
        screen.frame(app, vec![]);
        let now = button(screen);
        if now.is_some() && now == at {
            return;
        }
        at = now;
    }
    panic!("the window never came to rest");
}

/// The colours the map of the list draws a forward ORF and a reverse one in.
const FORWARD_FILL: Color32 = Color32::from_rgb(233, 208, 36);

#[test]
fn the_bar_of_the_selected_orf_has_a_heavy_line_round_it_and_a_click_on_a_bar_selects_its_orf() {
    use crate::ui::testing::Screen;
    let mut app = app_with(false, false);
    app.open_orf_list();
    let mut screen = Screen::new();
    steady(&mut screen, &mut app);
    let ink = app.palette().fg;
    assert!(screen.paths_stroked(ink, MARKED_LINE).is_empty());
    // the one bar of the map
    let bars = screen.paths_filled(FORWARD_FILL);
    assert_eq!(bars.len(), 1, "{bars:?}");
    screen.click_at(&mut app, bars[0].center());
    screen.settle(&mut app);
    assert_eq!(app.docs[0].selection(), Some((ORF_START, ORF_START + 126)));
    // the line is round that bar
    let lines = screen.paths_stroked(ink, MARKED_LINE);
    assert_eq!(lines.len(), 1, "{lines:?}");
    assert!(
        (lines[0].center() - bars[0].center()).length() < 1.0,
        "{lines:?} round {bars:?}"
    );
    // other bases: no line
    app.docs[0].select(ORF_START, ORF_START + 30);
    screen.settle(&mut app);
    assert!(screen.paths_stroked(ink, MARKED_LINE).is_empty());
}

#[test]
fn of_an_orf_that_runs_past_the_end_the_piece_before_the_end_is_what_is_selected_and_what_is_marked(
) {
    use crate::ui::testing::Screen;
    // the ORF of 126 bases from base 541 of 600, round the circle
    let mut seq = vec![b'C'; 600];
    let orf: Vec<u8> = [&b"ATG"[..], &b"GCT".repeat(40), &b"TAA"[..]].concat();
    for (j, &base) in orf.iter().enumerate() {
        seq[(540 + j) % 600] = base;
    }
    let mut doc = crate::seq::doc::Doc::new("t", seq);
    doc.circular = true;
    let mut app = App::new(crate::ui::app::Startup::default());
    app.docs = vec![crate::ui::state::DocView::new(doc)];
    app.prefs.orf_min = 90;
    app.prefs.orf_strand = Strand::Top;
    assert_eq!(starts(&app), [540]);
    app.open_orf_list();
    let mut screen = Screen::new();
    steady(&mut screen, &mut app);
    // two pieces of one bar: before the end, and from the start on
    let bars = screen.paths_filled(FORWARD_FILL);
    assert_eq!(bars.len(), 2, "{bars:?}");
    let first = app.orfs_of(0, false)[0].clone();
    app.select_orf(0, &first);
    assert_eq!(app.docs[0].selection(), Some((540, 600)));
    screen.settle(&mut app);
    let text = screen.text();
    assert!(
        text.contains(
            "the selection is the forward ORF at 541–600 (its part before the end of the sequence)"
        ),
        "{text}"
    );
    // the line is round the piece before the end only: the right one of the two
    let lines = screen.paths_stroked(app.palette().fg, MARKED_LINE);
    assert_eq!(lines.len(), 1, "{lines:?}");
    let right = if bars[0].center().x > bars[1].center().x {
        bars[0]
    } else {
        bars[1]
    };
    assert!(
        (lines[0].center() - right.center()).length() < 1.0,
        "{lines:?}"
    );
}

/// Forty forward ORFs of 126 bases, 150 bases apart from base 21: more rows than the list shows at once.
fn app_with_forty() -> App {
    let orf: Vec<u8> = [&b"ATG"[..], &b"GCT".repeat(40), &b"TAA"[..]].concat();
    let mut seq = vec![b'C'; 40 * 150 + 40];
    for k in 0..40 {
        seq[20 + 150 * k..20 + 150 * k + orf.len()].copy_from_slice(&orf);
    }
    let mut app = App::new(crate::ui::app::Startup::default());
    app.docs = vec![crate::ui::state::DocView::new(crate::seq::doc::Doc::new(
        "many", seq,
    ))];
    app.prefs.orf_min = 90;
    app.prefs.orf_strand = Strand::Top;
    app
}

/// The wheel of the mouse turned by `delta` points.
fn wheel(delta: eframe::egui::Vec2) -> eframe::egui::Event {
    eframe::egui::Event::MouseWheel {
        unit: eframe::egui::MouseWheelUnit::Point,
        delta,
        phase: eframe::egui::TouchPhase::Move,
        modifiers: eframe::egui::Modifiers::NONE,
    }
}

#[test]
fn the_list_scrolls_to_the_selected_orf_once_and_again_when_it_is_opened_anew() {
    use crate::ui::testing::Screen;
    let mut app = app_with_forty();
    assert_eq!(starts(&app).len(), 40);
    app.open_orf_list();
    let mut screen = Screen::new();
    steady(&mut screen, &mut app);
    // the row of the last ORF (it starts at base 5,871) is far below what the list shows
    let last = "5,871";
    assert!(screen.find("171").is_some(), "the second row is shown");
    assert!(screen.find(last).is_none());
    // selected somewhere else, it is scrolled to
    app.docs[0].select(5870, 5996);
    a_while(&mut screen, &mut app);
    let at = screen
        .find(last)
        .expect("the row of the selected ORF is in view");
    assert!(at.y < 640.0, "{at:?}");

    // the list can be scrolled away from it: it is not pulled back in every frame
    screen.frame(&mut app, vec![eframe::egui::Event::PointerMoved(at)]);
    screen.frame(&mut app, vec![wheel(eframe::egui::vec2(0.0, 4000.0))]);
    a_while(&mut screen, &mut app);
    assert!(screen.find(last).is_none(), "scrolled away, it stays away");
    assert!(screen.find("171").is_some());

    // closed and opened anew on the same selection, it shows it again
    app.orf_win.list_open = false;
    screen.frame(&mut app, vec![]);
    app.open_orf_list();
    a_while(&mut screen, &mut app);
    assert!(
        screen.find(last).is_some(),
        "opened anew, it shows the selected ORF"
    );
}

#[test]
fn a_reverse_orf_and_a_list_opened_on_a_stretch_name_the_selected_orf_too() {
    use crate::ui::testing::Screen;
    let mut reverse = app_with(false, true);
    reverse.open_orf_list();
    reverse.docs[0].select(ORF_START, ORF_START + 126);
    let mut screen = Screen::new();
    screen.settle(&mut reverse);
    let text = screen.text();
    assert!(
        text.contains("1 ORF in 600 bp · the selection is the reverse ORF at 101–226"),
        "{text}"
    );
    // opened on bases 51 to 500, whose ORF is numbered as in all of the sequence
    let mut stretch = app_with(false, false);
    stretch.prefs.selection_only = true;
    stretch.docs[0].select(50, 500);
    stretch.open_orf_list();
    stretch.docs[0].select(ORF_START, ORF_START + 126);
    let mut screen = Screen::new();
    screen.settle(&mut stretch);
    let text = screen.text();
    assert!(
        text.contains("1 ORF in bases 51–500 · the selection is the forward ORF at 101–226"),
        "{text}"
    );
}

#[test]
fn the_picture_of_a_selected_orf_selects_an_orf_of_the_demo_plasmid() {
    use crate::ui::testing::Screen;
    let mut app = App::new(crate::ui::app::Startup {
        view: "orf-selected".to_string(),
        ..crate::ui::app::Startup::default()
    });
    let mut screen = Screen::new();
    screen.settle(&mut app);
    let text = screen.text();
    assert!(
        text.contains("2 ORFs in 2,724 bp · the selection is the forward ORF at 288–902"),
        "{text}"
    );
}

#[test]
fn the_window_keeps_its_size_and_its_button_stays_in_it() {
    use crate::ui::testing::Screen;
    let mut app = app_with(false, false);
    app.open_orf_list();
    let mut screen = Screen::new();
    let button = |screen: &Screen| {
        screen
            .find("Add these ORFs as features")
            .expect("the button of the window")
    };
    for _ in 0..30 {
        screen.frame(&mut app, vec![]);
    }
    let early = button(&screen);
    for _ in 0..60 {
        screen.frame(&mut app, vec![]);
    }
    // it does not grow from frame to frame, and it is the size it was given (430 points of content under its title, a
    // little below the top bar: its button is at about 590), not as tall as the screen (800) nor shorter
    assert_eq!(button(&screen), early);
    assert!((early.y - 590.0).abs() < 15.0, "{early:?}");
}

const ORF_START: usize = 100;

/// 600 bases with one forward ORF (ATG, forty GCT, TAA) from base 100, counted from 0: 126 bases.
fn app_with(circular: bool, reverse: bool) -> App {
    let mut seq = vec![b'C'; 600];
    let orf: Vec<u8> = [&b"ATG"[..], &b"GCT".repeat(40), &b"TAA"[..]].concat();
    let orf = if reverse {
        crate::seq::dna::revcom(&orf)
    } else {
        orf
    };
    seq[ORF_START..ORF_START + orf.len()].copy_from_slice(&orf);
    let mut doc = crate::seq::doc::Doc::new("t", seq);
    doc.circular = circular;
    let mut app = App::new(crate::ui::app::Startup::default());
    app.docs = vec![crate::ui::state::DocView::new(doc)];
    app.prefs.orf_min = 90;
    app.prefs.orf_strand = if reverse { Strand::Bottom } else { Strand::Top };
    app
}

fn starts(app: &App) -> Vec<usize> {
    app.orfs_of(0, false).iter().map(|o| o.start).collect()
}

#[test]
fn with_selection_only_off_the_list_looks_at_all_of_the_sequence() {
    let mut app = app_with(false, false);
    app.docs[0].select(50, 300);

    app.open_orf_list();

    assert!(app.orf_region(0).is_none());
    assert_eq!(starts(&app), [ORF_START]);
}

#[test]
fn with_selection_only_on_the_list_is_opened_on_the_selected_stretch_numbered_as_in_all_of_it() {
    let mut app = app_with(false, false);
    app.prefs.selection_only = true;
    app.docs[0].select(50, 300);

    app.open_orf_list();

    assert_eq!(app.orf_region(0), Some((50, 300)));
    assert_eq!(starts(&app), [ORF_START]);
    // a stretch that begins inside the ORF has no ATG of it: no ORF
    app.docs[0].select(200, 400);
    app.open_orf_list();
    assert!(starts(&app).is_empty());
}

#[test]
fn with_selection_only_on_and_nothing_selected_the_list_is_opened_on_all_of_the_sequence() {
    let mut app = app_with(false, false);
    app.prefs.selection_only = true;
    app.docs[0].place_caret(10);

    app.open_orf_list();

    assert!(app.orf_region(0).is_none());
    assert_eq!(starts(&app).len(), 1);
}

#[test]
fn a_click_on_a_row_does_not_make_the_list_follow_the_selection_it_makes() {
    let mut app = app_with(false, false);
    app.prefs.selection_only = true;
    app.docs[0].select(50, 500);
    app.open_orf_list();
    let before = app.orfs_of(0, false);

    app.select_orf(0, &before[0]);

    assert_eq!(
        app.docs[0].selection(),
        Some((100, 226)),
        "the ORF is selected"
    );
    assert_eq!(
        app.orf_region(0),
        Some((50, 500)),
        "the list keeps its stretch"
    );
    assert_eq!(app.orfs_of(0, false), before);
}

#[test]
fn a_stretch_is_a_linear_piece_even_in_a_circular_sequence_so_an_orf_cut_by_its_end_is_open() {
    let mut app = app_with(true, false);
    app.prefs.selection_only = true;
    app.docs[0].select(50, 200);

    app.open_orf_list();
    let orfs = app.orfs_of(0, false);

    // ATG and thirty-two GCT end at base 200: no stop, so no wrapping round the origin either
    assert_eq!(orfs.len(), 1);
    assert_eq!(
        (orfs[0].start, orfs[0].length, orfs[0].open_end),
        (ORF_START, 99, true)
    );
    // in all of the circle the ORF is whole
    app.prefs.selection_only = false;
    app.open_orf_list();
    let whole = app.orfs_of(0, false);
    assert_eq!((whole[0].length, whole[0].open_end), (126, false));
}

#[test]
fn a_reverse_strand_orf_is_numbered_from_its_left_end_in_all_of_the_sequence() {
    let mut app = app_with(false, true);
    app.prefs.selection_only = true;
    app.docs[0].select(50, 300);

    app.open_orf_list();
    let orfs = app.orfs_of(0, false);

    assert_eq!(orfs.len(), 1);
    assert_eq!(
        (orfs[0].start, orfs[0].length, orfs[0].reverse),
        (ORF_START, 126, true)
    );
}

#[test]
fn a_region_of_another_sequence_or_one_the_sequence_has_shrunk_below_is_not_used() {
    let mut app = app_with(false, false);
    app.prefs.selection_only = true;
    app.docs[0].select(50, 300);
    app.open_orf_list();

    app.docs[0].id += 1;
    assert!(app.orf_region(0).is_none(), "another sequence");
    app.docs[0].id -= 1;
    app.orf_win.region = app.orf_win.region.map(|r| OrfRegion {
        from: 700,
        to: 900,
        ..r
    });
    assert!(app.orf_region(0).is_none(), "beyond the end");
    app.orf_win.region = app.orf_win.region.map(|r| OrfRegion {
        from: 500,
        to: 900,
        ..r
    });
    assert_eq!(app.orf_region(0), Some((500, 600)), "cut to the end");
}

#[test]
fn the_heading_the_amino_acids_and_the_length_say_what_was_looked_at_and_what_is_cut_off() {
    assert_eq!(orf_heading(1, None, 2724), "1 ORF in 2,724 bp");
    assert_eq!(
        orf_heading(3, Some((50, 300)), 2724),
        "3 ORFs in bases 51–300"
    );
    let whole = Orf {
        start: 0,
        length: 126,
        reverse: false,
        frame: 0,
        open_end: false,
    };
    let cut = Orf {
        open_end: true,
        ..whole.clone()
    };
    assert_eq!(
        (amino_acids(&whole), length_text(&whole)),
        (41, "126".to_string())
    );
    assert_eq!(
        (amino_acids(&cut), length_text(&cut)),
        (42, "126+".to_string())
    );
}
