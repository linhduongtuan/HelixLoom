use super::*;

#[test]
fn the_longest_orf_of_the_demo_plasmid_is_translated_when_nothing_is_picked() {
    let doc = crate::sample::demo_plasmid();
    let protein = protein_of(&doc, None, GENETIC_CODES[0].table).expect("a protein");
    assert!(protein.starts_with(b"M"));
    assert!(protein.len() > 100, "{}", protein.len());
    assert!(!protein.contains(&b'*'));
}

#[test]
fn a_picked_stretch_is_read_from_its_first_base_up_to_a_stop() {
    let doc = Doc::new("x", b"ATGAAATTTTAGGGGCCC".to_vec());
    let table = GENETIC_CODES[0].table;
    assert_eq!(
        protein_of(&doc, Some((0, 18)), table),
        Some(b"MKF".to_vec())
    );
    // nothing before the stop codon: nothing to fold
    assert_eq!(protein_of(&doc, Some((9, 12)), table), None);
    // nothing picked and no ORF long enough
    assert_eq!(protein_of(&doc, None, table), None);
}

#[test]
fn the_example_structure_loads_and_a_predicted_one_starts_on_confidence() {
    let mut page = StructurePage::default();
    page.show_pdb("ubiquitin", EXAMPLE_PDB)
        .expect("the example");
    assert_eq!(page.style.scheme, Scheme::Secondary);
    assert_eq!(
        page.loaded.as_ref().expect("loaded").fold.residues.len(),
        76
    );
    // a prediction: confidence in the temperature column, no EXPDTA
    let predicted: String = EXAMPLE_PDB
        .lines()
        .filter(|l| l.starts_with("ATOM"))
        .map(|l| format!("{}  0.91{}\n", &l[..60], &l[66..]))
        .collect();
    page.show_pdb("prediction", &predicted)
        .expect("a prediction");
    assert_eq!(page.style.scheme, Scheme::Confidence);
    let s = page.loaded.as_ref().expect("loaded").fold.summary();
    assert!((s.confidence.expect("a confidence") - 91.0).abs() < 0.1);
    assert!(page.show_pdb("nothing", "not a structure").is_err());
}

#[test]
fn a_new_structure_clears_the_pick_and_keeps_the_way_it_is_drawn() {
    let mut page = StructurePage::default();
    page.show_pdb("ubiquitin", EXAMPLE_PDB)
        .expect("the example");
    page.style.picked = Some(5);
    page.style.representation = Representation::Backbone;
    page.camera.zoom = 4.0;
    page.show_pdb("again", EXAMPLE_PDB).expect("the example");
    assert_eq!(page.style.picked, None);
    assert_eq!(page.style.representation, Representation::Backbone);
    assert_eq!(page.camera, Camera::default());
}

#[test]
fn an_atom_file_with_no_protein_in_it_is_refused_and_leaves_the_page_as_it_was() {
    let ligand = "HETATM    1 FE   HEM A   1       0.000   0.000   0.000  1.00 10.00          FE\n";
    let mut page = StructurePage::default();
    assert!(page
        .show_pdb("heme", ligand)
        .unwrap_err()
        .contains("no protein chain"));
    assert!(page.loaded.is_none());
}

#[test]
fn a_file_that_cannot_be_read_is_reported_not_raised() {
    let mut page = StructurePage::default();
    page.open_path(Path::new("/definitely/not/here.pdb"));
    assert!(
        page.message.starts_with("Could not read"),
        "{}",
        page.message
    );
    assert!(page.loaded.is_none());
}

/// Press "Predict structure" in a window whose heavy jobs go to a server on this computer, which has (or has not)
/// a program that folds on its GPU; give back the request that arrived and what the page made of the answer.
fn predict_on_the_workstation(gpu_fold: bool) -> (Request, StructurePage) {
    use crate::remote::gpus::Gpu;
    use crate::remote::jobs::{Config, Status};
    use crate::remote::tasks::{Folded, Kind, Response};
    use crate::remote::wire::Container;
    use crate::ui::app::Startup;
    const TOKEN: &str = "a-token-of-at-least-sixteen-chars";
    let folder = crate::pipeline::tools::work_folder("structure-remote").expect("a folder");
    let mut config = Config::new(TOKEN.into(), folder);
    if gpu_fold {
        config.gpus = vec![Gpu {
            id: 0,
            name: "GPU A".into(),
        }];
        // never started: the server's workers are off
        config.fold_command = Some(vec!["/bin/true".into()]);
    }
    let addr = "127.0.0.1:0".parse().expect("an address");
    let server = crate::remote::server::start(config, addr, false).expect("a server");

    let mut app = App::new(Startup::default());
    let ctx = egui::Context::default();
    let w = &mut app.platform.workstation;
    (w.address, w.token) = (server.addr.to_string(), TOKEN.into());
    w.connect();
    w.set_compute(crate::remote::settings::Compute::Workstation);
    // the page asks nothing of the workstation until it has said what it has
    let end = Instant::now() + Duration::from_secs(20);
    while app.platform.workstation.tools().is_none() {
        app.workstation_tick(&ctx);
        assert!(Instant::now() < end, "the server did not answer");
        std::thread::sleep(Duration::from_millis(20));
    }
    assert_eq!(app.platform.workstation.folds_on_gpu(), gpu_fold);
    let residues = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
    app.platform.structure.sequence = residues.into();
    let clicks = Clicks {
        predict: true,
        ..Clicks::default()
    };
    app.apply_structure_clicks(clicks, None, &ctx);
    assert!(app.platform.structure.is_running());

    let end = Instant::now() + Duration::from_secs(20);
    let queued = loop {
        if let Some(job) = server.shared.list().into_iter().next() {
            break job;
        }
        assert!(Instant::now() < end, "the job never arrived");
        std::thread::sleep(Duration::from_millis(20));
    };
    assert_eq!((queued.kind, queued.needs_gpu()), (Kind::Fold, gpu_fold));
    assert_eq!(queued.label, "Fold of 76 residues");
    let folder = server.shared.job_dir(&queued.id);
    let sent = std::fs::read(folder.join("request.bin")).expect("a request");
    let request = Request::decode(Kind::Fold, &Container::decode(&sent).expect("a container"))
        .expect("a request");
    let folded = Folded {
        pdb: EXAMPLE_PDB.into(),
        engine: "a stand-in".into(),
        seconds: 2.0,
    };
    std::fs::write(
        folder.join("result.bin"),
        Response::Fold(folded).encode().encode(),
    )
    .expect("a result");
    let mut done = queued;
    done.status = Status::Done;
    server.shared.replace_for_tests(done);

    while app.platform.structure.is_running() {
        app.structure_tick(&ctx);
        assert!(
            Instant::now() < end + Duration::from_secs(40),
            "no structure came back"
        );
        std::thread::sleep(Duration::from_millis(20));
    }
    (request, std::mem::take(&mut app.platform.structure))
}

#[test]
fn with_the_workstation_on_the_sequence_goes_there_and_the_structure_comes_back() {
    // a workstation with a folding program of its own is asked to use it; one without asks the web service
    for (gpu_fold, local) in [(true, true), (false, false)] {
        let (request, page) = predict_on_the_workstation(gpu_fold);
        let Request::Fold {
            sequence,
            local: asked,
        } = request
        else {
            panic!("a structure was asked for");
        };
        assert_eq!((sequence.len(), asked), (76, local));
        assert!(page.loaded.is_some(), "{}", page.message);
        assert!(page.message.starts_with("Predicted in"), "{}", page.message);
        // the prediction says what made it: the sequence that was sent, who was asked to fold it, and what answered
        let made = page.made.clone().expect("a record");
        assert_eq!(
            (made.analysis.as_str(), made.ran_on.as_str()),
            ("Structure prediction", "the workstation")
        );
        assert_eq!(made.inputs, [Input::of("protein sequence", &sequence)]);
        assert_eq!(
            made.parameters,
            [("predictor".to_string(), predictor(Some(local)).to_string())]
        );
        assert_eq!(
            made.programs,
            [
                "a stand-in".to_string(),
                format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))
            ]
        );
        let (inputs, settings) = page.as_it_is(Some(local));
        assert!(made.changed_since(&inputs, &settings).is_empty());
    }
}

/// A page that shows a prediction with its record, as one asked from this computer leaves it.
fn page_with_a_prediction() -> StructurePage {
    let mut page = StructurePage {
        sequence: "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
            .into(),
        ..StructurePage::default()
    };
    page.show_pdb("ESMFold prediction", EXAMPLE_PDB)
        .expect("a structure");
    let (inputs, parameters) = described(&page.checked().expect("a protein"), None);
    let asked = Asked {
        analysis: "Structure prediction",
        inputs,
        parameters,
        ran_on: "this computer",
    };
    page.made = Some(asked.made(vec!["the ESM Atlas web service".into()], 3.5));
    page
}

#[test]
fn the_record_of_a_prediction_says_what_has_changed_and_goes_when_another_structure_is_shown() {
    assert_ne!(predictor(Some(true)), predictor(Some(false)));
    assert_eq!(predictor(Some(false)), predictor(None));
    let mut page = page_with_a_prediction();
    let made = page.made.clone().expect("a record");
    let (inputs, settings) = page.as_it_is(None);
    assert!(made.changed_since(&inputs, &settings).is_empty());
    // another sequence in the box, and a workstation that folds with its own program: both are said
    page.sequence.push('W');
    let (inputs, settings) = page.as_it_is(Some(true));
    assert_eq!(
        made.changed_since(&inputs, &settings),
        [
            "protein sequence has changed".to_string(),
            format!(
                "predictor was {}, is now {}",
                predictor(None),
                predictor(Some(true))
            )
        ]
    );
    // Run again is offered while who would fold now is who made this one: the service, asked from here or by a workstation
    assert!(made_by(&made, None) && made_by(&made, Some(false)));
    assert!(!made_by(&made, Some(true)));
    // a structure that is opened was made by nothing here: the record of the prediction goes
    page.load_example();
    assert!(page.loaded.is_some() && page.made.is_none());
}

#[test]
fn the_record_of_a_prediction_is_copied_and_running_again_asks_for_the_prediction() {
    use crate::ui::app::Startup;
    use crate::ui::platform::Workspace;
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Structure;
    app.platform.structure = page_with_a_prediction();
    let record = app
        .platform
        .structure
        .made
        .as_ref()
        .expect("a record")
        .to_text();
    assert!(
        record.starts_with("Structure prediction\nMade on "),
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
        app.platform.structure.message.contains("Copied the record"),
        "{}",
        app.platform.structure.message
    );
    assert!(!app.checks_win.open);
    screen.click(&mut app, made_ui::CHECKED);
    screen.settle(&mut app);
    assert!(app.checks_win.open);
    // Run again is the press of Predict: here with no curl to ask the service with, so nothing is sent and that is said
    let ctx = egui::Context::default();
    let clicks = Clicks {
        made: made_ui::Asks {
            again: true,
            ..made_ui::Asks::default()
        },
        ..Clicks::default()
    };
    app.apply_structure_clicks(clicks, None, &ctx);
    assert!(!app.platform.structure.is_running());
    assert!(
        app.platform
            .structure
            .message
            .contains("curl was not found"),
        "{}",
        app.platform.structure.message
    );
}

/// An alignment (rows aligned, none gapped) of ubiquitin and two relatives that differ from it in a few places.
fn ubiquitin_family() -> crate::pipeline::steps::Aligned {
    let ubiquitin =
        b"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG".to_vec();
    let change = |at: &[usize], to: u8| {
        let mut row = ubiquitin.clone();
        at.iter().for_each(|&k| row[k] = to);
        row
    };
    crate::pipeline::steps::Aligned {
        rows: vec![
            ("ubiquitin".to_string(), ubiquitin.clone()),
            ("relative 1".to_string(), change(&[5, 14, 20], b'W')),
            ("relative 2".to_string(), change(&[6, 8, 10], b'W')),
        ],
        engine: "test".into(),
        seconds: 0.0,
    }
}

fn app_with_the_family() -> App {
    let mut app = App::new(crate::ui::app::Startup::default());
    app.platform.structure.load_example();
    app.platform.alignment.set_result(ubiquitin_family());
    app
}

#[test]
fn a_structure_is_coloured_by_the_alignment_and_says_where_the_colours_came_from() {
    let mut app = app_with_the_family();
    app.colour_structure_by_alignment();
    let page = &app.platform.structure;
    assert_eq!(page.style.scheme, Scheme::Conservation);
    let fold = &page.loaded.as_ref().expect("a structure").fold;
    assert!(fold.residues.iter().all(|r| r.conservation.is_some()));
    assert!(
        page.message.contains("chain A") && page.message.contains("ubiquitin"),
        "{}",
        page.message
    );
    assert!(page.message.contains("76 of 76"), "{}", page.message);
    assert!(page.conservation_from.is_some());
}

#[test]
fn an_alignment_of_another_protein_leaves_the_structure_as_it_was_and_says_so() {
    let mut app = App::new(crate::ui::app::Startup::default());
    app.platform.structure.load_example();
    let other = b"MSDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEY".to_vec();
    let mut similar = other.clone();
    similar[3] = b'W';
    app.platform
        .alignment
        .set_result(crate::pipeline::steps::Aligned {
            rows: vec![
                ("other".to_string(), other),
                ("similar".to_string(), similar),
            ],
            engine: "test".into(),
            seconds: 0.0,
        });
    let before = app.platform.structure.style.scheme;
    app.colour_structure_by_alignment();
    let page = &app.platform.structure;
    assert_eq!(page.style.scheme, before);
    assert!(
        page.message.contains("nothing was coloured"),
        "{}",
        page.message
    );
    assert!(page.conservation_from.is_none());
}

#[test]
fn colours_made_from_scores_that_have_since_changed_are_taken_off_and_the_view_goes_back() {
    let ctx = egui::Context::default();
    let mut app = app_with_the_family();
    app.colour_structure_by_alignment();
    // nothing has changed: a frame leaves the colours
    app.structure_tick(&ctx);
    assert_eq!(app.platform.structure.style.scheme, Scheme::Conservation);
    // a new score for the same alignment
    app.platform
        .alignment
        .score_for_test(crate::ui::platform::alignment::Score::Entropy);
    app.platform.alignment.refresh_scores();
    app.structure_tick(&ctx);
    let page = &app.platform.structure;
    assert_ne!(page.style.scheme, Scheme::Conservation);
    let fold = &page.loaded.as_ref().expect("a structure").fold;
    assert!(fold.residues.iter().all(|r| r.conservation.is_none()));
    assert!(page.message.contains("taken off"), "{}", page.message);
    assert!(page.conservation_from.is_none());
    // and again for a new alignment
    app.colour_structure_by_alignment();
    assert_eq!(app.platform.structure.style.scheme, Scheme::Conservation);
    app.platform.alignment.set_result(ubiquitin_family());
    app.structure_tick(&ctx);
    assert_ne!(app.platform.structure.style.scheme, Scheme::Conservation);
    // a new structure forgets that it was coloured
    app.colour_structure_by_alignment();
    app.platform.structure.load_example();
    assert!(app.platform.structure.conservation_from.is_none());
}

#[test]
fn in_a_small_window_the_structure_is_still_drawn_under_the_whole_row_of_choices() {
    use crate::ui::app::{App, Startup};
    use crate::ui::platform::Workspace;
    use crate::ui::testing::Screen;
    let mut app = App::new(Startup::default());
    app.load_example(Workspace::Structure);
    // a window that holds 853 by 533 points: the choices above the structure no longer fit in two rows
    app.prefs.appearance.scale = 150;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    let look = app.prefs.look();
    let window = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1280.0 / 1.5, 800.0 / 1.5));
    // every choice is in the window, and so is the structure, with at least the room it is promised
    let choices = [
        "Ribbon",
        "All bonds",
        "Confidence",
        "N to C",
        "Interface",
        "Conservation",
        "Spin",
    ];
    let lowest = choices
        .iter()
        .map(|choice| {
            let at = screen
                .find(choice)
                .unwrap_or_else(|| panic!("no {choice:?} on the screen"));
            assert!(
                window.contains(at),
                "{choice:?} is at {at:?}, out of the window"
            );
            at.y
        })
        .fold(0.0, f32::max);
    let view = screen
        .rects_filled(look.surface)
        .into_iter()
        .find(|r| r.min.x > 400.0 && r.height() >= MIN_VIEWER_HEIGHT - 1.0)
        .expect("the structure, in a frame");
    assert!(
        view.min.y > lowest,
        "under the choices: {view:?}, the lowest choice at {lowest}"
    );
    assert!(
        view.min.y < window.max.y - 100.0,
        "{view:?} begins in the window"
    );
    // the choice that is made is drawn as made, and one that is not is not: something is filled behind the first only
    let behind = |screen: &Screen, choice: &str| {
        let at = screen.find(choice).expect("the choice");
        screen
            .fills_at(at)
            .into_iter()
            .filter(|(room, _)| room.height() < 60.0)
            .count()
    };
    assert_eq!(
        (behind(&screen, "Ribbon"), behind(&screen, "Backbone")),
        (1, 0)
    );
    screen.click(&mut app, "Backbone");
    screen.settle(&mut app);
    assert_eq!(
        app.platform.structure.style.representation,
        Representation::Backbone
    );
    assert_eq!(
        (behind(&screen, "Ribbon"), behind(&screen, "Backbone")),
        (0, 1)
    );
    // one that may not be chosen yet is pale, and stays as it is when it is clicked
    assert_eq!(screen.is_pale("Conservation"), Some(true));
    assert_eq!(screen.is_pale("Chain"), Some(false));
    let before = app.platform.structure.style.scheme;
    screen.click(&mut app, "Conservation");
    screen.settle(&mut app);
    assert_eq!(app.platform.structure.style.scheme, before);
    // and one that may is chosen by a click
    screen.click(&mut app, "Chain");
    screen.settle(&mut app);
    assert_eq!(app.platform.structure.style.scheme, Scheme::Chain);
}
