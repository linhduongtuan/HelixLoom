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
    assert_eq!((queued.kind, queued.needs_gpu), (Kind::Fold, gpu_fold));
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
    }
}
