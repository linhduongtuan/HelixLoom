use crate::pipeline::steps::Aligned;
use crate::ui::app::{Action, App, Startup};
use crate::ui::platform::alignment::Role;

/// The same 1,000 bases that primer3's and the probe's own tests use, as `n` rows of one alignment.
fn aligned(n: usize) -> Aligned {
    let mut x = 777u32;
    let bases: Vec<u8> = (0..1000)
        .map(|_| {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((x >> 16) & 3) as usize]
        })
        .collect();
    Aligned {
        rows: (0..n)
            .map(|i| (format!("strain{}", i + 1), bases.clone()))
            .collect(),
        engine: "test".into(),
        seconds: 0.0,
    }
}

fn app_with(n: usize) -> App {
    let mut app = App::new(Startup::default());
    app.platform.alignment.set_result(aligned(n));
    app
}

#[test]
fn without_an_alignment_the_window_says_to_align_first() {
    let mut app = App::new(Startup::default());
    app.run_aln_primers();
    assert!(app.aln_win.report.is_none());
    assert!(
        app.aln_win.message.contains("no alignment"),
        "{}",
        app.aln_win.message
    );
}

#[test]
fn the_button_opens_the_window_and_the_pick_fills_it() {
    let mut app = app_with(3);
    assert!(!app.aln_win.open);
    app.run_aln_primers();
    let report = app.aln_win.report.as_ref().expect("a report");
    assert_eq!((report.rows, report.columns), (3, 1000));
    assert!(!report.designs.is_empty());
    assert!(
        app.aln_win.message.contains("covers 3 of 3"),
        "{}",
        app.aln_win.message
    );
}

#[test]
fn a_share_lower_than_the_one_the_sequences_have_finds_nothing_and_says_why() {
    let mut app = app_with(2);
    // the rows differ everywhere: one is the complement of the other
    let swap = |b: u8| match b {
        b'A' => b'C',
        b'C' => b'A',
        b'G' => b'T',
        _ => b'G',
    };
    let mut a = aligned(2);
    a.rows[1].1 = a.rows[1].1.iter().map(|&b| swap(b)).collect();
    app.platform.alignment.set_result(a);
    app.run_aln_primers();
    let said = app.aln_win.message.clone();
    assert!(
        said.contains("No column is shared") && said.contains("lower the share"),
        "{said}"
    );
    assert!(app
        .aln_win
        .report
        .as_ref()
        .is_some_and(|r| r.designs.is_empty()));
}

#[test]
fn a_bad_setting_is_an_error_in_the_window_and_not_a_crash() {
    let mut app = app_with(2);
    app.aln_win.use_region = true;
    app.aln_win.region = (900, 1500);
    app.run_aln_primers();
    // the region is clamped by the form, but a region outside the alignment that got here is refused
    assert!(app.aln_win.report.is_none());
    assert!(
        app.aln_win.message.contains("not inside"),
        "{}",
        app.aln_win.message
    );
}

#[test]
fn showing_a_design_marks_its_columns_and_selects_the_first() {
    let mut app = app_with(3);
    app.run_aln_primers();
    let mut out: Vec<Action> = Vec::new();
    let asks = super::aln_primers::test_asks_show(0);
    let ctx = eframe::egui::Context::default();
    app.apply_aln_primer_asks(&ctx, asks, &mut out);
    let marks = app.platform.alignment.marks().to_vec();
    assert_eq!(marks.len(), 2);
    assert_eq!((marks[0].role, marks[1].role), (Role::Left, Role::Right));
    assert_eq!(
        app.platform.alignment.selected_column(),
        Some(marks[0].first)
    );
    assert!(matches!(out.first(), Some(Action::Go(_))));
    // a new alignment clears them
    app.platform.alignment.set_result(aligned(3));
    assert!(app.platform.alignment.marks().is_empty());
}

#[test]
fn a_design_is_not_shown_on_an_alignment_that_has_changed_since() {
    let mut app = app_with(3);
    app.run_aln_primers();
    let mut other = aligned(3);
    // a base that is not the one there (the position may well hold an A)
    other.rows[0].1[10] = if other.rows[0].1[10] == b'A' {
        b'C'
    } else {
        b'A'
    };
    app.platform.alignment.set_result(other);
    let mut out = Vec::new();
    let ctx = eframe::egui::Context::default();
    app.apply_aln_primer_asks(&ctx, super::aln_primers::test_asks_show(0), &mut out);
    assert!(app.platform.alignment.marks().is_empty() && out.is_empty());
}

#[test]
fn the_consensus_opens_with_the_primers_as_features_in_the_right_places() {
    let mut app = app_with(3);
    app.aln_win.settings.with_probe = true;
    (
        app.aln_win.settings.primers.product_min,
        app.aln_win.settings.primers.product_max,
    ) = (70, 150);
    app.run_aln_primers();
    let docs_before = app.docs.len();
    let mut out = Vec::new();
    let ctx = eframe::egui::Context::default();
    app.apply_aln_primer_asks(&ctx, super::aln_primers::test_asks_consensus(0), &mut out);
    assert_eq!(app.docs.len(), docs_before + 1);
    let doc = &app.docs.last().expect("a document").doc;
    assert_eq!(doc.seq.len(), 1000);
    let report = app.aln_win.report.as_ref().expect("a report");
    let d = &report.designs[0];
    let names: Vec<&str> = doc.features.iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names, ["Left primer 1", "Right primer 1", "Probe 1"]);
    let left = &doc.features[0];
    assert_eq!(left.spans, [(d.left.oligo.start, d.left.oligo.end)]);
    assert!(!left.reverse && doc.features[1].reverse);
    // a feature's bases are the oligo (as it is on the consensus strand)
    let (a, b) = left.spans[0];
    assert_eq!(&doc.seq[a..b], d.left.oligo.seq.as_slice());
}

#[test]
fn a_region_keeps_every_oligo_inside_it_and_the_marks_are_the_columns_of_the_oligos() {
    let mut app = app_with(3);
    app.aln_win.use_region = true;
    app.aln_win.region = (101, 700);
    app.run_aln_primers();
    let report = app.aln_win.report.clone().expect("a report");
    assert_eq!(report.usable, 600);
    assert!(!report.designs.is_empty());
    for d in &report.designs {
        for s in [&d.left, &d.right] {
            assert!(s.oligo.start >= 100 && s.oligo.end <= 700, "{s:?}");
        }
    }
    let mut out: Vec<Action> = Vec::new();
    let ctx = eframe::egui::Context::default();
    app.apply_aln_primer_asks(&ctx, super::aln_primers::test_asks_show(0), &mut out);
    let d = &report.designs[0];
    let marks = app.platform.alignment.marks();
    assert_eq!(
        (marks[0].first, marks[0].last),
        (d.left.oligo.start, d.left.oligo.end - 1)
    );
    assert_eq!(
        (marks[1].first, marks[1].last),
        (d.right.oligo.start, d.right.oligo.end - 1)
    );
}
