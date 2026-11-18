use super::*;
use std::time::Instant;

fn finish(w: &mut MotifWin, ctx: &egui::Context) {
    let until = Instant::now() + Duration::from_secs(120);
    while w.running() && Instant::now() < until {
        std::thread::sleep(Duration::from_millis(20));
        w.poll(ctx);
    }
    assert!(!w.running(), "the search did not end");
}

#[test]
fn sequences_are_read_from_fasta_or_from_lines_and_too_few_or_too_many_are_refused() {
    let (names, seqs) = read_sequences(">a one\nACGT\nACGT\n>b\nTTTT\n").expect("sequences");
    assert_eq!(names, ["a", "b"]);
    assert_eq!(seqs, [b"ACGTACGT".to_vec(), b"TTTT".to_vec()]);
    let (names, seqs) = read_sequences("ACGT ACGT\n\nTTGG\n").expect("sequences");
    assert_eq!(names, ["sequence 1", "sequence 2"]);
    assert_eq!(seqs, [b"ACGTACGT".to_vec(), b"TTGG".to_vec()]);
    assert!(read_sequences("").is_err());
    assert!(read_sequences("ACGT").unwrap_err().contains("at least two"));
    let many = "ACGT\n".repeat(MOST_SEQUENCES + 1);
    assert!(read_sequences(&many).unwrap_err().contains("more than"));
}

#[test]
fn the_logo_draws_a_certain_column_as_one_tall_letter_and_a_flat_one_as_nothing() {
    let certain = Motif::from_sites("m", &[&b"ACGT"[..]; 30]).expect("a motif");
    let columns = logo(&certain);
    assert_eq!(columns.len(), 4);
    for (column, want) in columns.iter().zip("ACGT".chars()) {
        let top = column[3];
        assert_eq!(top.0, want);
        assert!(top.1 > 0.9, "{top:?}");
        // sorted: the tallest is last, none is higher than the column
        assert!(column.windows(2).all(|p| p[0].1 <= p[1].1));
        assert!(column.iter().map(|l| l.1).sum::<f32>() <= 1.0001);
    }
    let flat = Motif {
        name: "flat".into(),
        probs: vec![[0.25; 4]; 3],
    };
    assert!(logo(&flat).iter().all(|c| c.iter().all(|l| l.1 < 1e-6)));
}

#[test]
fn a_search_runs_in_the_background_and_its_motif_is_the_one_in_the_example_which_a_scan_then_finds_again(
) {
    let ctx = egui::Context::default();
    let mut w = MotifWin {
        text: crate::sample::motif_demo(),
        settings: Settings {
            width: 10,
            ..Settings::default()
        },
        ..MotifWin::default()
    };
    w.start(&ctx);
    assert!(w.running());
    // a second press while it goes does nothing
    w.start(&ctx);
    finish(&mut w, &ctx);
    let run = w.run.as_ref().expect("a run");
    assert_eq!(run.found.len(), 1);
    let found = &run.found[0];
    let consensus = found.motif.consensus();
    let planted = "ACCGTTAGCA";
    let other_strand: String = planted
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect();
    let off = |word: &str| {
        consensus
            .bytes()
            .zip(word.bytes())
            .filter(|(a, b)| a != b)
            .count()
    };
    assert!(off(planted).min(off(&other_strand)) <= 2, "{consensus}");
    assert!(found.sites.len() >= 9, "{} sites", found.sites.len());
    assert!(
        w.message.contains("1 motif in 12 sequences"),
        "{}",
        w.message
    );
    // a scan at a strict p-value finds the sites back, each as significant as was asked
    w.p_value = 1e-4;
    w.scan_chosen();
    assert!(w.hit_total >= 8, "{} hits", w.hit_total);
    assert!(w.hits.iter().all(|h| h.p_value <= 1e-4 + 1e-12));
    assert!(w.chosen_meme().expect("text").contains("MOTIF motif_1"));
}

#[test]
fn a_search_that_cannot_start_says_why_and_starts_nothing() {
    let ctx = egui::Context::default();
    let mut w = MotifWin::default();
    w.start(&ctx);
    assert!(!w.running());
    assert!(w.message.contains("at least two"), "{}", w.message);
    // a motif wider than every sequence: the search says so
    w.text = "ACGT\nACGT\n".into();
    w.settings.width = 12;
    w.start(&ctx);
    finish(&mut w, &ctx);
    assert!(w.run.is_none());
    assert!(w.message.contains("as long as the motif"), "{}", w.message);
}

#[test]
fn the_logo_takes_the_colours_of_the_trace_when_the_safe_colours_are_chosen_and_keeps_its_own_otherwise(
) {
    use crate::ui::appearance::{Appearance, Colours};
    for dark in [false, true] {
        let plain = Appearance::default().look(dark);
        let safe = Appearance {
            colours: Colours::Safe,
            ..Appearance::default()
        }
        .look(dark);
        for letter in ['A', 'C', 'G', 'T'] {
            assert_ne!(
                letter_colour(letter, &plain),
                letter_colour(letter, &safe),
                "{letter}"
            );
            assert_eq!(
                letter_colour(letter, &safe),
                crate::ui::design::trace_color(letter as u8, &safe)
            );
        }
        // the letters of the designed logo are the colours it always had
        assert_eq!(letter_colour('A', &plain), Color32::from_rgb(0, 160, 60));
        assert_eq!(letter_colour('T', &plain), Color32::from_rgb(210, 40, 40));
    }
}

#[test]
fn what_the_workstation_found_is_taken_only_if_it_is_about_the_sequences_that_were_sent() {
    let (names, seqs) = read_sequences(&crate::sample::motif_demo()).expect("the example");
    let settings = Settings {
        width: 10,
        ..Settings::default()
    };
    let result = remote_motif::run(&MotifJob {
        seqs: seqs.clone(),
        settings,
    })
    .expect("motifs");
    let run = from_workstation(result.clone(), names.clone(), seqs.clone(), settings, 1.5)
        .expect("an answer that fits");
    assert_eq!((run.found.len(), run.seconds), (result.found.len(), 1.5));
    assert_eq!(run.names, names);
    // fewer sequences than the answer has sites in, or shorter ones: refused with the reason, and no run is made
    let said = from_workstation(
        result.clone(),
        names[..1].to_vec(),
        seqs[..1].to_vec(),
        settings,
        1.0,
    )
    .err()
    .expect("refused");
    assert!(said.contains("do not fit"), "{said}");
    let short: Vec<Vec<u8>> = seqs.iter().map(|s| s[..12].to_vec()).collect();
    assert!(from_workstation(result, names, short, settings, 1.0).is_err());
}

fn example_window() -> MotifWin {
    MotifWin {
        text: crate::sample::motif_demo(),
        settings: Settings {
            width: 10,
            ..Settings::default()
        },
        ..MotifWin::default()
    }
}

#[test]
fn a_search_leaves_a_record_of_its_sequences_and_settings_which_says_what_has_changed_since() {
    let ctx = egui::Context::default();
    let mut w = example_window();
    w.start(&ctx);
    finish(&mut w, &ctx);
    let made = w.made.clone().expect("a record");
    let (names, seqs) = read_sequences(&w.text).expect("the example");
    assert_eq!(
        (made.analysis.as_str(), made.ran_on.as_str()),
        ("Motif search", "this computer")
    );
    assert_eq!(made.inputs.len(), 12);
    assert_eq!(
        (made.inputs[3].name.as_str(), made.inputs[3].size),
        (names[3].as_str(), seqs[3].len())
    );
    let keys: Vec<&str> = made.parameters.iter().map(|(k, _)| k.as_str()).collect();
    assert_eq!(
        keys,
        [
            "width",
            "motifs",
            "both strands",
            "starting words",
            "rounds",
            "effort"
        ]
    );
    assert!(made
        .parameters
        .contains(&("width".to_string(), "10".to_string())));
    assert_eq!(
        made.programs,
        [format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))]
    );
    assert_eq!(made.seconds, w.run.as_ref().expect("a run").seconds);
    let (inputs, settings) = w.as_it_is();
    assert!(made.changed_since(&inputs, &settings).is_empty());
    // a setting is changed and a base of one sequence: both are said, the sequence by its name
    w.settings.width = 8;
    w.settings.both_strands = !w.settings.both_strands;
    let at = w.text.rfind('A').expect("a base");
    w.text.replace_range(at..=at, "C");
    let (inputs, settings) = w.as_it_is();
    let changes = made.changed_since(&inputs, &settings);
    assert_eq!(changes.len(), 3, "{changes:?}");
    assert!(changes[0].ends_with("has changed"), "{changes:?}");
    assert_eq!(changes[1], "width was 10, is now 8");
    assert!(changes[2].starts_with("both strands was "), "{changes:?}");
    // text that can no longer be read as sequences: every sequence of the record is gone, and that is said, not hidden
    w.text.clear();
    let (inputs, _) = w.as_it_is();
    assert!(inputs.is_empty());
    // a search that finds nothing or cannot run takes the motifs away, and their record with them
    w.text = "ACGT\nACGT\n".into();
    w.settings.width = 12;
    w.start(&ctx);
    finish(&mut w, &ctx);
    assert!(w.run.is_none() && w.made.is_none());
}

#[test]
fn the_record_of_the_motifs_is_copied_and_running_again_brings_their_settings_back() {
    use crate::ui::provenance as made_ui;
    let mut app = App::new(crate::ui::app::Startup::default());
    app.motif_win.show_example();
    let made = app.motif_win.made.clone().expect("a record");
    let record = made.to_text();
    assert!(record.starts_with("Motif search\nMade on "), "{record}");
    let mut screen = crate::ui::testing::Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::TITLE);
    screen.settle(&mut app);
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::COPY);
    assert_eq!(screen.copied(), Some(record.as_str()));
    assert!(
        app.motif_win.message.contains("Copied the record"),
        "{}",
        app.motif_win.message
    );
    // the link under the record opens the page that says how each analysis was checked
    assert!(!app.checks_win.open);
    screen.click(&mut app, made_ui::CHECKED);
    screen.settle(&mut app);
    assert!(app.checks_win.open);
    app.checks_win.open = false;
    // the settings are changed; Run again searches with the ones the motifs were found with
    app.motif_win.settings.width = 7;
    app.motif_win.settings.motifs = 3;
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::AGAIN);
    assert_eq!(
        (app.motif_win.settings.width, app.motif_win.settings.motifs),
        (10, 1)
    );
    let until = Instant::now() + Duration::from_secs(120);
    while (app.motif_win.running() || app.motif_win.made.is_none()) && Instant::now() < until {
        std::thread::sleep(Duration::from_millis(20));
        screen.frame(&mut app, vec![]);
    }
    let again = app.motif_win.made.clone().expect("a new record");
    assert_eq!(
        (&again.parameters, &again.inputs),
        (&made.parameters, &made.inputs)
    );
    assert!(
        app.motif_win.message.contains("motif in 12 sequences"),
        "a search ran: {}",
        app.motif_win.message
    );
}

#[test]
fn motifs_fetched_from_the_workstations_list_have_no_record() {
    let mut app = App::new(crate::ui::app::Startup::default());
    app.motif_win.show_example();
    assert!(app.motif_win.made.is_some());
    let (_, seqs) = read_sequences(&app.motif_win.text).expect("the example");
    let result = remote_motif::run(&MotifJob {
        seqs,
        settings: app.motif_win.settings,
    })
    .expect("motifs");
    app.show_motif_result(result);
    assert!(app.motif_win.run.is_some(), "{}", app.motif_win.message);
    assert!(
        app.motif_win.made.is_none(),
        "the settings it was found with are not known here"
    );
    assert!(
        app.motif_win.message.ends_with("on the workstation."),
        "{}",
        app.motif_win.message
    );
    // the example is then searched on this computer, and its record and what is said of it say so
    app.motif_win.show_example();
    let made = app.motif_win.made.as_ref().expect("a record");
    assert_eq!(made.ran_on, "this computer");
    assert!(
        !app.motif_win.message.contains("workstation"),
        "{}",
        app.motif_win.message
    );
}
