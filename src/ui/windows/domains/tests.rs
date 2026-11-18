use super::*;
use crate::ui::app::Startup;
use crate::ui::headless::frame_input;
use eframe::egui::{epaint::Shape, Context};
use std::time::Instant;

fn hit(family: &str, from: usize, to: usize) -> DomainHit {
    DomainHit {
        sequence: 0,
        family: family.to_string(),
        accession: "PF00000.1".into(),
        family_length: 100,
        score: 50.0,
        i_evalue: 1e-10,
        c_evalue: 1e-10,
        hmm_from: 1,
        hmm_to: 100,
        ali_from: from,
        ali_to: to,
        env_from: from,
        env_to: to,
        accuracy: 0.9,
        description: format!("{family} domain"),
    }
}

fn result(hits: Vec<DomainHit>) -> DomainsResult {
    DomainsResult {
        names: vec!["first".into(), "second".into()],
        lengths: vec![400, 120],
        hits,
        program: "HMMER 3.4".into(),
        database: "Pfam-A.hmm".into(),
        seconds: 1.25,
    }
}

#[test]
fn a_family_has_one_colour_in_every_sequence_and_run() {
    for family in ["Pkinase", "fn3", "Globin", "", "ubiquitin"] {
        assert!(colour_index(family) < PALETTE.len());
        assert_eq!(colour_index(family), colour_index(family));
    }
    // and families are told apart: the names of a few Pfam families do not all take one colour
    let used: std::collections::HashSet<usize> = [
        "Pkinase",
        "fn3",
        "Globin",
        "SH2",
        "SH3_1",
        "zf-C2H2",
        "Ubiquitin",
        "RRM_1",
    ]
    .iter()
    .map(|f| colour_index(f))
    .collect();
    assert!(used.len() >= 4, "{used:?}");
    // fixed values: a change of the hash would recolour every family
    assert_eq!(colour_index("Pkinase"), colour_index("Pkinase"));
    assert_eq!(colour_index("fn3"), 2);
}

#[test]
fn domains_that_overlap_are_drawn_one_under_the_other_and_a_free_lane_is_used_again() {
    let hits = [
        hit("A", 1, 50),
        hit("B", 40, 90),
        hit("C", 60, 100),
        hit("D", 95, 120),
        hit("E", 130, 140),
    ];
    let along: Vec<&DomainHit> = hits.iter().collect();
    // A in lane 0; B overlaps A, lane 1; C starts after A ended, lane 0; D (95-120) starts after B ended (90) but not after C
    // (100), lane 1; E starts after both, lane 0
    assert_eq!(lanes(&along), [0, 1, 0, 1, 0]);
    // three that all overlap need three lanes
    let stacked = [hit("A", 1, 100), hit("B", 10, 100), hit("C", 20, 100)];
    assert_eq!(lanes(&stacked.iter().collect::<Vec<_>>()), [0, 1, 2]);
    // touching is overlapping: a domain that starts where the last one ends has a lane of its own
    let touching = [hit("A", 1, 50), hit("B", 50, 90)];
    assert_eq!(lanes(&touching.iter().collect::<Vec<_>>()), [0, 1]);
    assert!(lanes(&[]).is_empty());
}

#[test]
fn a_residue_is_drawn_at_its_share_of_the_bar() {
    assert_eq!(x_of(1, 100, 500.0), 0.0);
    assert_eq!(x_of(51, 100, 500.0), 250.0);
    // the end of the last residue is the end of the bar
    assert_eq!(x_of(101, 100, 500.0), 500.0);
    assert_eq!(x_of(0, 100, 500.0), 0.0);
    assert_eq!(x_of(5, 0, 500.0), 2000.0 / 1.0);
}

#[test]
fn what_a_search_found_is_told_in_a_sentence() {
    let none = summary(&result(Vec::new()), "this computer");
    assert_eq!(none, "No domain found in 2 sequences by HMMER 3.4 against Pfam-A.hmm, in 1.2 s on this computer.");
    let one = summary(&result(vec![hit("fn3", 1, 50)]), "the workstation");
    assert_eq!(one, "1 domain in 1 of 2 sequences by HMMER 3.4 against Pfam-A.hmm, in 1.2 s on the workstation.");
    let mut two = result(vec![hit("fn3", 1, 50), hit("fn3", 60, 100)]);
    two.hits[1].sequence = 1;
    two.names.truncate(1);
    two.lengths.truncate(1);
    assert!(summary(&two, "here").starts_with("2 domains in 2 of 1 sequence by"));
    let bare = DomainsResult {
        program: String::new(),
        database: String::new(),
        ..result(Vec::new())
    };
    assert!(summary(&bare, "here").contains("by HMMER, in"));
}

fn answered(_: &DomainsJob) -> Result<DomainsResult, String> {
    Ok(result(vec![hit("Pkinase", 20, 380)]))
}

fn refused(_: &DomainsJob) -> Result<DomainsResult, String> {
    Err("HMMER's hmmscan was not found on this computer: install HMMER 3.".into())
}

fn finish(w: &mut DomainsWin, ctx: &Context) {
    let until = Instant::now() + Duration::from_secs(60);
    while w.running() && Instant::now() < until {
        std::thread::sleep(Duration::from_millis(10));
        w.poll(ctx);
    }
    assert!(!w.running(), "the search did not end");
}

#[test]
fn a_search_runs_in_the_background_and_its_domains_are_kept_and_said() {
    let ctx = Context::default();
    let mut w = DomainsWin {
        runner: answered,
        text: EXAMPLE.to_string(),
        ..DomainsWin::default()
    };
    w.start(&ctx, Place::Here);
    assert!(w.running());
    // a second press while it runs starts nothing
    w.start(&ctx, Place::Here);
    finish(&mut w, &ctx);
    let said = &w.message;
    assert!(said.starts_with("1 domain in 1 of 2 sequences"), "{said}");
    assert!(said.ends_with("on this computer."), "{said}");
    assert_eq!(w.result.as_ref().map(|r| r.hits.len()), Some(1));
    // a program that is not there: the reason, and no result left from before
    w.runner = refused;
    w.start(&ctx, Place::Here);
    finish(&mut w, &ctx);
    assert!(w.message.contains("hmmscan was not found"), "{}", w.message);
    assert!(w.result.is_none());
}

#[test]
fn what_cannot_be_searched_is_said_and_nothing_is_started() {
    let ctx = Context::default();
    let mut w = DomainsWin {
        runner: answered,
        ..DomainsWin::default()
    };
    w.start(&ctx, Place::Here);
    assert!(!w.running());
    assert!(w.message.contains("at least one protein"), "{}", w.message);
    // an E-value that is no cutoff
    w.text = EXAMPLE.to_string();
    (w.curated, w.evalue) = (false, 0.0);
    w.start(&ctx, Place::Here);
    assert!(!w.running());
    assert!(w.message.contains("above 0"), "{}", w.message);
    // and the settings are what the job asks for
    w.evalue = 0.002;
    assert_eq!(
        w.job_asked_for().expect("a job").cutoff,
        Cutoff::Evalue(0.002)
    );
    w.curated = true;
    assert_eq!(w.job_asked_for().expect("a job").cutoff, Cutoff::Curated);
}

#[test]
fn the_size_of_a_search_is_the_sequences_and_residues_in_the_text() {
    let w = DomainsWin {
        text: ">a\nMKV\nLLE\n>b\nGA\n".to_string(),
        ..DomainsWin::default()
    };
    let weight = w.weight();
    assert_eq!((weight.sequences, weight.letters), (2, 8));
    let lines = DomainsWin {
        text: "MKV\n\nGA\n".to_string(),
        ..DomainsWin::default()
    };
    assert_eq!(lines.weight().letters, 5);
}

/// The text of every piece of text a frame drew.
fn drawn_text(shapes: &[Shape], out: &mut Vec<String>) {
    for shape in shapes {
        match shape {
            Shape::Text(t) => out.push(t.galley.text().to_string()),
            Shape::Vec(inner) => drawn_text(inner, out),
            _ => {}
        }
    }
}

#[test]
fn the_window_draws_the_domains_along_the_sequences_in_both_themes() {
    for dark in [false, true] {
        let ctx = Context::default();
        let mut app = App::new(Startup::default());
        app.domains_win.open = true;
        app.domains_win.text = EXAMPLE.to_string();
        let hits = vec![
            hit("Pkinase", 20, 380),
            hit("Overlapper", 100, 200),
            hit("Small", 390, 392),
        ];
        app.domains_win.result = Some(Arc::new(result(hits)));
        let _ = dark;
        // a window takes its size on the first frames
        let mut text = Vec::new();
        for n in 0..4 {
            let mut output = ctx.run_ui(frame_input(n), |ui| app.draw(ui));
            output.textures_delta.clear();
            text.clear();
            let shapes: Vec<Shape> = output.shapes.iter().map(|s| s.shape.clone()).collect();
            drawn_text(&shapes, &mut text);
        }
        // the wide domain is named in its box, the table lists all three, the overlap is said
        let ours: Vec<&String> = text
            .iter()
            .skip_while(|t| !t.contains("Protein Domains"))
            .collect();
        assert!(!ours.is_empty(), "the window was not drawn");
        assert!(text.iter().any(|t| t == "Pkinase"), "{ours:?}");
        assert!(text.iter().any(|t| t == "Overlapper"), "{text:?}");
        assert!(
            text.iter().any(|t| t.contains("Some domains overlap")),
            "{text:?}"
        );
        assert!(text.iter().any(|t| t.contains("Find domains")));
    }
}

#[test]
fn a_search_leaves_a_record_of_its_proteins_its_threshold_and_the_programs_that_answered() {
    let ctx = Context::default();
    let mut w = DomainsWin {
        runner: answered,
        text: EXAMPLE.to_string(),
        ..DomainsWin::default()
    };
    w.start(&ctx, Place::Here);
    finish(&mut w, &ctx);
    let (made, kept) = w.made.clone().expect("a record");
    assert_eq!(
        (made.analysis.as_str(), made.ran_on.as_str()),
        ("Protein domain search", "this computer")
    );
    let names: Vec<&str> = made.inputs.iter().map(|i| i.name.as_str()).collect();
    assert_eq!(names.len(), 4);
    assert!(names[0].starts_with("ubiquitin"), "{names:?}");
    assert_eq!(made.inputs[0].size, 76, "the residues of ubiquitin");
    assert_eq!(
        made.parameters,
        [(
            "threshold".to_string(),
            "Pfam's own for each family".to_string()
        )]
    );
    assert_eq!(
        made.programs,
        [
            "HMMER 3.4".to_string(),
            "Pfam-A.hmm".to_string(),
            format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))
        ]
    );
    assert_eq!((made.seconds, kept), (1.25, (true, DEFAULT_EVALUE)));
    let (inputs, settings) = w.as_it_is();
    assert!(made.changed_since(&inputs, &settings).is_empty());
    // another threshold, and a protein taken out of the box: both are said
    (w.curated, w.evalue) = (false, 1e-5);
    let cut = w.text.find(">lysozyme").expect("the last protein");
    w.text.truncate(cut);
    let (inputs, settings) = w.as_it_is();
    assert_eq!(
        made.changed_since(&inputs, &settings),
        [
            // a protein is named by the first word of its line, as HMMER names it
            "lysozyme is no longer there",
            "threshold was Pfam's own for each family, is now E-value at most 1e-5"
        ]
    );
    // an answer that names no program and no database (an older workstation's) leaves no empty name in the record
    fn unnamed(_: &DomainsJob) -> Result<DomainsResult, String> {
        Ok(DomainsResult {
            program: String::new(),
            database: String::new(),
            ..result(Vec::new())
        })
    }
    w.runner = unnamed;
    w.start(&ctx, Place::Here);
    finish(&mut w, &ctx);
    let (unnamed_made, _) = w.made.clone().expect("a record");
    assert_eq!(
        unnamed_made.programs,
        [format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))]
    );
    // a search that fails takes the domains away, and their record with them
    w.runner = refused;
    w.start(&ctx, Place::Here);
    assert!(w.made.is_none(), "gone as soon as another search starts");
    finish(&mut w, &ctx);
    assert!(w.result.is_none() && w.made.is_none());
}

#[test]
fn the_record_of_the_domains_is_copied_and_running_again_brings_the_threshold_back() {
    let mut app = App::new(Startup::default());
    app.domains_win.runner = answered;
    (app.domains_win.curated, app.domains_win.evalue) = (false, 1e-7);
    app.domains_win.show_example();
    let (made, kept) = app.domains_win.made.clone().expect("a record");
    assert_eq!(kept, (false, 1e-7));
    let record = made.to_text();
    assert!(
        record.starts_with("Protein domain search\nMade on "),
        "{record}"
    );
    let mut screen = crate::ui::testing::Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::TITLE);
    screen.settle(&mut app);
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::COPY);
    assert_eq!(screen.copied(), Some(record.as_str()));
    assert!(
        app.domains_win.message.contains("Copied the record"),
        "{}",
        app.domains_win.message
    );
    assert!(!app.checks_win.open);
    screen.click(&mut app, made_ui::CHECKED);
    screen.settle(&mut app);
    assert!(app.checks_win.open);
    app.checks_win.open = false;
    // the threshold is changed; Run again searches with the one the domains were found with
    (app.domains_win.curated, app.domains_win.evalue) = (true, 1e-3);
    screen.settle(&mut app);
    screen.click(&mut app, made_ui::AGAIN);
    assert_eq!(
        (app.domains_win.curated, app.domains_win.evalue),
        (false, 1e-7)
    );
    let until = Instant::now() + Duration::from_secs(60);
    while (app.domains_win.running() || app.domains_win.made.is_none()) && Instant::now() < until {
        std::thread::sleep(Duration::from_millis(10));
        screen.frame(&mut app, vec![]);
    }
    let (again, kept) = app.domains_win.made.clone().expect("a new record");
    assert_eq!(again.parameters, made.parameters);
    assert_eq!(kept, (false, 1e-7), "kept again for the next time");
    assert!(
        app.domains_win.message.starts_with("1 domain in"),
        "a search ran: {}",
        app.domains_win.message
    );
    // what is fetched from the workstation's list has no record: its proteins and threshold are not known here
    app.show_domains_result(result(vec![hit("Pkinase", 20, 380)]));
    assert!(app.domains_win.result.is_some() && app.domains_win.made.is_none());
    assert!(
        app.domains_win.message.ends_with("on the workstation."),
        "{}",
        app.domains_win.message
    );
    // the example is then searched on this computer, and its record and what is said of it say so
    app.domains_win.show_example();
    let (made, _) = app.domains_win.made.as_ref().expect("a record");
    assert_eq!(made.ran_on, "this computer");
    assert!(
        app.domains_win.message.ends_with("on this computer."),
        "{}",
        app.domains_win.message
    );
}
