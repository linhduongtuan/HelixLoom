use crate::ui::app::{App, Startup};

/// The demo plasmid with the window open and the primers picked around its MCS.
fn app() -> App {
    App::new(Startup {
        view: "primer3".to_string(),
        ..Startup::default()
    })
}

#[test]
fn without_the_probe_choice_the_window_picks_primers_alone_as_it_did() {
    let app = app();
    assert!(!app.p3_win.with_probe);
    assert!(!app.p3_win.found.is_empty());
    assert!(app.p3_win.found.iter().all(|(_, probe)| probe.is_none()));
    assert!(!app.p3_win.message.contains("probe"));
}

#[test]
fn with_the_probe_choice_every_pair_comes_with_a_probe_between_its_primers() {
    let mut app = app();
    app.p3_win.with_probe = true;
    app.run_primer3();
    let found = &app.p3_win.found;
    assert!(!found.is_empty() && found.len() <= app.p3_win.settings.num_return);
    for (pair, probe) in found {
        let probe = probe.as_ref().expect("a probe");
        assert!(probe.start > pair.left.end && probe.end < pair.right.start);
        assert!(probe.seq[0] != b'G');
    }
    assert!(
        app.p3_win.message.contains("probe"),
        "{}",
        app.p3_win.message
    );
}

#[test]
fn when_no_probe_fits_the_window_says_so_and_does_not_say_that_no_pair_fits() {
    let mut app = app();
    app.p3_win.with_probe = true;
    // a probe that must melt forty degrees above the primers is none
    (app.p3_win.probe.min_above, app.p3_win.probe.max_above) = (40.0, 50.0);
    app.run_primer3();
    assert!(app.p3_win.found.is_empty());
    let said = &app.p3_win.message;
    assert!(
        said.contains("probe") && !said.starts_with("No pair fits"),
        "{said}"
    );
    // and when not even a pair fits, it is the pairs that it names
    app.p3_win.settings.product_max = app.p3_win.settings.product_min;
    app.run_primer3();
    assert!(
        app.p3_win.message.starts_with("No pair fits"),
        "{}",
        app.p3_win.message
    );
}

#[test]
fn every_assay_puts_its_probe_on_the_strand_and_the_span_it_has() {
    let mut app = app();
    app.p3_win.with_probe = true;
    app.run_primer3();
    let found = app.p3_win.found.clone();
    assert!(found.len() > 1);
    for (k, (_, probe)) in found.iter().enumerate() {
        let probe = probe.as_ref().expect("a probe");
        app.add_pair_features(k);
        let made = app.docs[0].doc.features.last().expect("a feature");
        assert_eq!(made.name, format!("Probe {}", k + 1));
        assert_eq!(made.spans, vec![(probe.start, probe.end)]);
        assert_eq!(made.reverse, probe.reverse);
    }
    // the probes are not all on one strand, or the check above would not tell a fixed strand from the right one
    let minus = found
        .iter()
        .filter(|f| f.1.as_ref().is_some_and(|p| p.reverse))
        .count();
    assert!(
        minus > 0 && minus < found.len(),
        "{minus} of {}",
        found.len()
    );
}

#[test]
fn features_are_not_put_on_a_sequence_that_the_pairs_were_not_picked_for() {
    let mut app = app();
    app.p3_win.with_probe = true;
    app.run_primer3();
    let before = app.docs[0].doc.features.len();
    // the sequence is edited after the pairs were picked
    app.docs[0].edit(|d| d.seq.insert(0, b'A'));
    app.add_pair_features(0);
    assert_eq!(app.docs[0].doc.features.len(), before);
    // picked again, they go on
    app.run_primer3();
    app.add_pair_features(0);
    assert_eq!(app.docs[0].doc.features.len(), before + 3);
}

#[test]
fn the_features_of_an_assay_are_the_two_primers_and_the_probe_on_its_own_strand() {
    let mut app = app();
    app.p3_win.with_probe = true;
    app.run_primer3();
    let (pair, probe) = app.p3_win.found[0].clone();
    let probe = probe.expect("a probe");
    let before = app.docs[0].doc.features.len();
    app.add_pair_features(0);
    let added = &app.docs[0].doc.features[before..];
    assert_eq!(added.len(), 3);
    let kinds: Vec<&str> = added.iter().map(|f| f.kind.as_str()).collect();
    assert_eq!(kinds, ["primer_bind", "primer_bind", "misc_feature"]);
    assert_eq!(added[0].spans, vec![(pair.left.start, pair.left.end)]);
    let made = &added[2];
    assert_eq!(made.name, "Probe 1");
    assert_eq!(made.spans, vec![(probe.start, probe.end)]);
    assert_eq!(made.reverse, probe.reverse);
    // what the note holds is what would be ordered, 5' to 3'
    let note = &made
        .qualifiers
        .iter()
        .find(|q| q.0 == "note")
        .expect("a note")
        .1;
    assert!(
        note.contains(&String::from_utf8_lossy(&probe.seq).into_owned()),
        "{note}"
    );
    // without a probe it is the two primers, as before
    let mut alone = self::app();
    let before = alone.docs[0].doc.features.len();
    alone.add_pair_features(0);
    assert_eq!(alone.docs[0].doc.features.len(), before + 2);
}
