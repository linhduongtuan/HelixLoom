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
