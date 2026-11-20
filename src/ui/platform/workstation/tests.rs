use super::*;
use crate::pipeline::steps::Named;
use crate::pipeline::tools::{work_folder, Installed};
use crate::pipeline::Choice;
use crate::remote::client::GpuState;
use crate::remote::gpus::Gpu;
use crate::remote::jobs::Config;
use crate::remote::run::{execute, Env};
use crate::remote::server::{start, Server};
use crate::remote::species::{self, SpeciesJob};
use crate::remote::tasks::{Folded, Kind, Request};
use crate::remote::weight::{
    HEAVY_GENERATIONS, HEAVY_LETTERS, HEAVY_REPLICATES, HEAVY_SEQUENCES, HEAVY_TREES,
};
use crate::remote::wire::Container;
use crate::seq::estimate::Analysis;
use crate::seq::phylocnn::{Priors, Task};
use crate::seq::primer3::Settings;
use crate::seq::tree::DistanceModel;
use crate::ui::app::Startup;
use crate::ui::headless::frame_input;
use crate::ui::platform::route::Source;
use crate::ui::platform::Platform;

const TOKEN: &str = "a-token-of-at-least-sixteen-chars";
const UBIQUITIN: &str = include_str!("../../../../assets/structure/1ubq.pdb");

fn status(tools: &[(&str, &str)]) -> ServerStatus {
    ServerStatus {
        name: "lab-gpu".into(),
        version: "0.1.0".into(),
        protocol: PROTOCOL,
        running: 1,
        queued: 2,
        cpu_jobs: 2,
        max_upload: 1 << 20,
        keep_hours: 168,
        fold: false,
        gpus: vec![GpuState {
            id: 0,
            name: "NVIDIA A".into(),
            job: None,
        }],
        tools: tools
            .iter()
            .map(|(n, v)| ((*n).to_string(), (*v).to_string()))
            .collect(),
    }
}

fn job(id: &str, kind: Kind, status: Status) -> JobInfo {
    JobInfo {
        id: id.into(),
        kind,
        label: format!("job {id}"),
        status,
        created: 100,
        started: 110,
        finished: 0,
        gpus: Vec::new(),
        gpus_wanted: 0,
        error: String::new(),
    }
}

/// A page that has been given an address and a token, connected and set to send every job there.
fn connected(address: &str) -> WorkstationPage {
    let mut page = WorkstationPage {
        address: address.into(),
        token: TOKEN.into(),
        ..WorkstationPage::open(None)
    };
    page.connect();
    page.set_compute(Compute::Workstation);
    page
}

/// Whether a job of no size at all goes to the workstation: true only when the choice sends every job there.
fn goes_there(page: &WorkstationPage) -> bool {
    page.place(&Weight::default()).is_there()
}

/// A job just over the size at which Automatic sends it to the workstation.
fn big() -> Weight {
    Weight::of_sequences(HEAVY_SEQUENCES + 1, 0)
}

fn isolates() -> Vec<Named> {
    crate::sample::tree_demo()
        .into_iter()
        .map(|d| (d.name, d.seq))
        .collect()
}

#[test]
fn the_page_opens_as_it_was_left_and_keeps_what_it_is_told() {
    let dir = work_folder("workstation-page").expect("a folder");
    let fresh = WorkstationPage::open(Some(dir.clone()));
    assert!(fresh.active.is_none() && fresh.compute == Compute::Local && !goes_there(&fresh));

    let mut page = WorkstationPage {
        address: "lab-gpu".into(),
        token: format!("  {TOKEN} "),
        ..fresh
    };
    page.connect();
    assert!(page.problem.is_empty(), "{}", page.problem);
    assert_eq!(
        page.active.as_ref().map(Remote::address),
        Some("lab-gpu:7878")
    );
    // being connected is not the same as sending the jobs there
    assert!(!goes_there(&page));
    page.set_compute(Compute::Workstation);
    assert!(goes_there(&page));

    let again = WorkstationPage::open(Some(dir.clone()));
    assert_eq!(
        (again.address.as_str(), again.token.as_str()),
        ("lab-gpu:7878", TOKEN)
    );
    assert!(again.compute == Compute::Workstation && goes_there(&again));

    // forgetting takes the token off the disk too, and the choice with it
    let mut page = again;
    page.forget();
    assert!(!goes_there(&page) && page.token.is_empty() && page.address.is_empty());
    let cleared = WorkstationPage::open(Some(dir));
    assert!(
        cleared.active.is_none() && cleared.token.is_empty() && cleared.compute == Compute::Local
    );
}

#[test]
fn heavy_jobs_go_only_where_the_person_said_and_only_over_a_good_connection() {
    let mut page = WorkstationPage::open(None);
    for mode in [Compute::Workstation, Compute::Automatic] {
        page.set_compute(mode);
        assert!(
            page.compute == Compute::Local && !goes_there(&page) && !page.place(&big()).is_there(),
            "no connection to send to: {mode:?}"
        );
    }

    page.token = TOKEN.into();
    for (address, said) in [
        ("lab-gpu:99999", "port"),
        ("https://lab-gpu", "TLS"),
        ("", "address"),
    ] {
        page.address = address.into();
        page.connect();
        assert!(
            page.active.is_none() && page.problem.contains(said),
            "{address}: {}",
            page.problem
        );
    }
    page.address = "lab-gpu".into();
    page.token = "   ".into();
    page.connect();
    assert!(
        page.active.is_none() && page.problem.contains("token"),
        "{}",
        page.problem
    );

    // a saved file that says "use it" without a token switches nothing on, in the old words or the new
    for (label, text) in [
        ("old", "address=lab-gpu\ntoken=\nuse=1\n"),
        (
            "workstation",
            "address=lab-gpu\ntoken=\ncompute=workstation\n",
        ),
        ("automatic", "address=lab-gpu\ntoken=\ncompute=automatic\n"),
    ] {
        let dir = work_folder(&format!("workstation-bad-saved-{label}")).expect("a folder");
        std::fs::write(dir.join(settings::FILE), text).expect("a file");
        assert_eq!(WorkstationPage::open(Some(dir)).compute, Compute::Local);
    }
}

#[test]
fn the_workstations_programs_stand_in_for_this_computers_while_jobs_go_there() {
    let mine = Installed {
        path: "/usr/bin/mafft".into(),
        version: "v7.526 (2023)".into(),
    };
    let mut platform = Platform::new(None);
    platform.programs = Some(vec![("MAFFT", Some(mine.clone())), ("IQ-TREE", None)]);
    assert_eq!(platform.installed("MAFFT"), Some(mine));
    assert!(platform.installed("IQ-TREE").is_none());
    assert!(platform.programs_known_at(&Place::Here));

    platform.workstation = connected("lab-gpu");
    let there = platform.place(&Weight::default());
    assert!(there.is_there());
    // nothing is known until it has said; if it cannot be reached, nothing is there
    assert!(
        !platform.programs_known_at(&there) && platform.installed_at("MAFFT", &there).is_none()
    );
    platform.workstation.problem = "no route to host".into();
    assert!(platform.programs_known_at(&there) && platform.installed_at("MAFFT", &there).is_none());
    platform.workstation.problem.clear();
    platform.workstation.status = Some(status(&[("IQ-TREE", "3.1.0")]));
    assert!(platform.installed_at("MAFFT", &there).is_none());
    let theirs = platform
        .installed_at("IQ-TREE", &there)
        .expect("the workstation's program");
    assert_eq!(theirs.version, "3.1.0");
    assert!(
        theirs.path.as_os_str().is_empty(),
        "only the version is known of a program elsewhere"
    );
    assert_eq!(platform.workstation.gpu_name(), Some("NVIDIA A"));
    // this computer's programs are still there for a job that stays
    assert!(platform.installed("MAFFT").is_some());

    // with the choice back on this computer a job's programs are this computer's again
    platform.workstation.set_compute(Compute::Local);
    let here = platform.place(&Weight::default());
    assert!(!here.is_there());
    assert!(
        platform.installed_at("MAFFT", &here).is_some()
            && platform.installed_at("IQ-TREE", &here).is_none()
    );
}

#[test]
fn the_choice_of_where_jobs_run_sends_each_job_by_its_size_and_is_kept() {
    let dir = work_folder("workstation-choice").expect("a folder");
    let mut page = WorkstationPage {
        address: "lab-gpu".into(),
        token: TOKEN.into(),
        ..WorkstationPage::open(Some(dir.clone()))
    };
    page.connect();
    let small = Weight::of_sequences(HEAVY_SEQUENCES, 1000);
    let on_gpu = Weight::on_gpu();
    let heavy = [
        big(),
        on_gpu,
        Weight::of_sequences(3, HEAVY_LETTERS + 1),
        Weight::of_tree(10, 100, Analysis::NeighborJoining, HEAVY_REPLICATES + 1, 0),
        Weight::of_tree(10, 100, Analysis::Bayesian, 0, HEAVY_GENERATIONS + 1),
        Weight::of_trees(&";".repeat(HEAVY_TREES + 1), false),
        Weight::of_trees(";", true),
    ];
    let light = [
        small,
        Weight::default(),
        Weight::of_tree(10, 100, Analysis::NeighborJoining, HEAVY_REPLICATES, 0),
        // settings the method does not use do not make a job big
        Weight::of_tree(10, 100, Analysis::Likelihood, 100_000, 100_000_000),
        Weight::of_tree(10, 100, Analysis::Bayesian, 100_000, HEAVY_GENERATIONS),
        Weight::of_trees(&";".repeat(HEAVY_TREES), false),
    ];

    // this computer: nothing leaves it, whatever the size
    assert_eq!(page.compute, Compute::Local);
    for weight in heavy.iter().chain(&light) {
        assert!(!page.place(weight).is_there(), "{weight:?}");
    }
    // the workstation: everything goes there, even a job of no size
    page.set_compute(Compute::Workstation);
    for weight in heavy.iter().chain(&light) {
        assert!(page.place(weight).is_there(), "{weight:?}");
    }
    // automatic: the big ones and the ones that use its GPU go there, the small ones stay
    page.set_compute(Compute::Automatic);
    for weight in &heavy {
        assert!(
            weight.is_heavy() && page.place(weight).is_there(),
            "{weight:?}"
        );
    }
    for weight in &light {
        assert!(
            !weight.is_heavy() && !page.place(weight).is_there(),
            "{weight:?}"
        );
    }
    // what a page shows says the same
    let shown = page.view(&small);
    assert!(shown.connected && !shown.there && shown.mode == Compute::Automatic);
    assert!(page.view(&big()).there);

    // the choice is kept with the connection and read again, and forgetting the workstation forgets it
    let again = WorkstationPage::open(Some(dir.clone()));
    assert_eq!(again.compute, Compute::Automatic);
    assert!(again.place(&big()).is_there() && !again.place(&small).is_there());
    page.forget();
    assert_eq!(WorkstationPage::open(Some(dir)).compute, Compute::Local);
}

#[test]
fn a_click_on_the_choice_of_a_page_is_applied_and_saved_at_the_next_frame() {
    let dir = work_folder("workstation-picker").expect("a folder");
    let mut app = App::new(Startup::default());
    let ctx = egui::Context::default();
    app.platform.workstation = WorkstationPage {
        address: "127.0.0.1:1".into(),
        token: TOKEN.into(),
        ..WorkstationPage::open(Some(dir.clone()))
    };
    app.platform.workstation.connect();
    // a page that is drawn while the workstation page is borrowed leaves its click in the context
    ctx.data_mut(|d| d.insert_temp(choice_id(), Compute::Automatic));
    app.workstation_tick(&ctx);
    assert_eq!(app.platform.workstation.compute, Compute::Automatic);
    assert_eq!(WorkstationPage::open(Some(dir)).compute, Compute::Automatic);
    // and the click is used up
    app.platform.workstation.set_compute(Compute::Local);
    app.workstation_tick(&ctx);
    assert_eq!(app.platform.workstation.compute, Compute::Local);
}

#[test]
fn looks_and_fetches_are_taken_in_and_failures_say_why() {
    let mut page = connected("lab-gpu");
    let jobs = vec![
        job("j1-1", Kind::Align, Status::Done),
        job("j1-2", Kind::Tree, Status::Queued),
    ];
    assert!(page.finish(Ok(Done::Looked(status(&[]), jobs))).is_none());
    assert_eq!(page.message, "Connected to lab-gpu.");
    assert_eq!(page.jobs.len(), 2);

    assert!(page
        .finish(Err("The workstation refused the token.".into()))
        .is_none());
    assert_eq!(page.problem, "The workstation refused the token.");
    // the last good look stays on the page beside the problem
    assert!(page.status.is_some() && page.jobs.len() == 2);

    page.finish(Ok(Done::Removed("j1-1".into(), false)));
    assert_eq!(
        page.jobs.iter().map(|j| j.id.as_str()).collect::<Vec<_>>(),
        ["j1-2"]
    );
    assert_eq!(page.message, "Removed the job.");
    // a job that was running is only asked to stop: it stays until the next look says what became of it
    page.finish(Ok(Done::Removed("j1-2".into(), true)));
    assert_eq!(page.jobs.len(), 1);
    assert!(page.message.contains("stop"), "{}", page.message);
    assert!(page.looked.is_none(), "the next look is at once");
    let folded = Folded {
        pdb: "ATOM\n".into(),
        engine: "x".into(),
        seconds: 1.0,
    };
    let fetched = page.finish(Ok(Done::Fetched(
        job("j1-2", Kind::Fold, Status::Done),
        Response::Fold(folded),
    )));
    assert!(fetched.is_some() && page.problem.is_empty());
}

#[test]
fn a_workstation_that_has_not_been_heard_from_is_tried_again_and_one_that_has_is_left_alone() {
    let ctx = egui::Context::default();
    let long_ago = |page: &mut WorkstationPage| {
        page.looked = Instant::now().checked_sub(REFRESH + Duration::from_millis(50));
    };
    let mut app = App::new(Startup::default());
    // nothing listens on port 1 of this computer
    app.platform.workstation = connected("127.0.0.1:1");
    app.workstation_tick(&ctx);
    assert!(
        app.platform.workstation.work.is_some(),
        "the first look starts at once"
    );
    let end = Instant::now() + Duration::from_secs(20);
    while app.platform.workstation.work.is_some() || app.platform.workstation.problem.is_empty() {
        app.workstation_tick(&ctx);
        assert!(Instant::now() < end, "the look never ended");
        std::thread::sleep(Duration::from_millis(10));
    }
    // too soon to ask again, then time to
    app.workstation_tick(&ctx);
    assert!(app.platform.workstation.work.is_none());
    long_ago(&mut app.platform.workstation);
    app.workstation_tick(&ctx);
    assert!(
        app.platform.workstation.work.is_some(),
        "tried again while on another page"
    );
    assert_eq!(app.workspace, Workspace::Home);

    // once it has been heard from, a page that is not showing does not ask again
    let mut app = App::new(Startup::default());
    app.platform.workstation = connected("127.0.0.1:1");
    app.platform.workstation.status = Some(status(&[]));
    long_ago(&mut app.platform.workstation);
    app.workstation_tick(&ctx);
    assert!(app.platform.workstation.work.is_none());
    // but while it shows, it keeps looking
    app.workspace = Workspace::Workstation;
    app.workstation_tick(&ctx);
    assert!(app.platform.workstation.work.is_some());
}

#[test]
fn durations_and_job_lines_are_written_for_reading() {
    assert_eq!(duration_text(5), "5 s");
    assert_eq!(duration_text(90), "1.5 min");
    assert_eq!(duration_text(7200), "2.0 h");
    assert_eq!(
        (
            keep_text(168).as_str(),
            keep_text(24).as_str(),
            keep_text(36).as_str()
        ),
        ("7 days", "24 h", "36 h")
    );
    let long = "x".repeat(100);
    assert_eq!(shorten(&long).chars().count(), LABEL_CHARS);
    assert!(shorten(&long).ends_with('…'));
    assert_eq!(shorten("short"), "short");

    let mut j = job("j1-1", Kind::Align, Status::Queued);
    assert_eq!(detail(&j), "waiting");
    j.gpus_wanted = 1;
    assert_eq!(detail(&j), "waiting for a GPU");
    j.gpus_wanted = 3;
    assert_eq!(detail(&j), "waiting for 3 GPUs");
    (j.status, j.gpus) = (Status::Running, vec![2]);
    assert_eq!(detail(&j), "on GPU 2");
    j.gpus = vec![0, 1, 3];
    assert_eq!(detail(&j), "on GPUs 0, 1, 3");
    (j.status, j.finished) = (Status::Done, 130);
    assert_eq!(detail(&j), "took 20 s");
    j.status = Status::Failed;
    assert_eq!(detail(&j), "");
}

#[test]
fn a_fetched_result_opens_in_the_page_it_belongs_to() {
    let env = Env::default();
    let mut app = App::new(Startup::default());

    let aligned = execute(
        Request::Align {
            seqs: isolates(),
            mafft: false,
            progressive: None,
        },
        &env,
    )
    .and_then(Response::into_align)
    .expect("aligned");
    app.open_remote_result(
        &job("j1-1", Kind::Align, Status::Done),
        Response::Align(aligned.clone()),
    );
    assert_eq!(app.workspace, Workspace::Alignment);
    let page = &app.platform.alignment;
    assert_eq!(
        page.result.as_ref().map(|a| a.rows.len()),
        Some(aligned.rows.len())
    );
    assert!(
        page.inputs.iter().all(|i| !i.bases.contains(&b'-')),
        "the list holds the sequences without gaps"
    );
    assert!(page.message.contains("job j1-1"), "{}", page.message);

    let tree = execute(
        Request::Tree {
            alignment: aligned.rows,
            analysis: Analysis::NeighborJoining,
            model: DistanceModel::P,
            bootstrap: 0,
            generations: 1000,
        },
        &env,
    )
    .expect("a tree");
    // the lines about a tree that was there before are not this tree's; its own come with it
    app.tree_win.details = vec!["of the tree before".to_string()];
    let Response::Tree(mut estimate) = tree else {
        panic!("a tree was asked for");
    };
    estimate.details.push("a line of its own".to_string());
    let lines = estimate.details.clone();
    let tree = Response::Tree(estimate);
    app.open_remote_result(&job("j1-2", Kind::Tree, Status::Done), tree);
    assert_eq!(app.workspace, Workspace::Trees);
    assert!(app.tree_win.tree.is_some());
    assert_eq!(app.tree_win.details, lines);

    let choice = Choice {
        align: true,
        tree: false,
        primers: false,
        use_installed: false,
        bootstrap: 0,
    };
    let run = execute(
        Request::Pipeline {
            seqs: isolates(),
            target: None,
            settings: Settings::default(),
            choice,
        },
        &env,
    )
    .expect("a pipeline");
    // a run of the page from before left its record: the fetched result is not what that record is of
    let stale = crate::provenance::Asked {
        analysis: "Pipeline",
        inputs: Vec::new(),
        parameters: Vec::new(),
        ran_on: "this computer",
    }
    .made(Vec::new(), 1.0);
    let ticked = (crate::ui::platform::pipelines::Source::Example, [true; 5]);
    app.platform.pipelines.made = Some((stale, ticked));
    app.open_remote_result(&job("j1-3", Kind::Pipeline, Status::Done), run);
    assert_eq!(app.workspace, Workspace::Pipelines);
    assert!(app.platform.pipelines.outcome.is_some());
    assert!(
        app.platform.pipelines.made.is_none(),
        "what the workstation made is not what the page's last run was made from"
    );

    // the networks are the authors' files, which this repository does not carry: without them there is no estimate to open
    if !crate::models::phylocnn_missing() {
        let newick = crate::seq::phylocnn::example_tree().0.trim().to_string();
        let estimates = execute(
            Request::PhyloCnn {
                newick,
                sampling: 0.25,
                task: Task::Bdss,
                priors: Priors::Hiv,
                gpu: false,
            },
            &env,
        )
        .expect("estimates");
        app.open_remote_result(&job("j1-4", Kind::PhyloCnn, Status::Done), estimates);
        assert!(app.phylocnn_win.open && !app.phylocnn_win.table.is_empty());
    }

    let fold = |pdb: &str| {
        Response::Fold(Folded {
            pdb: pdb.into(),
            engine: "x".into(),
            seconds: 3.0,
        })
    };
    app.open_remote_result(&job("j1-5", Kind::Fold, Status::Done), fold(UBIQUITIN));
    assert_eq!(app.workspace, Workspace::Structure);
    assert!(app.platform.structure.loaded.is_some());
    // text that is not a structure is said so, not drawn
    app.open_remote_result(&job("j1-6", Kind::Fold, Status::Done), fold("nothing"));
    assert!(
        app.platform.structure.message.contains("cannot be shown"),
        "{}",
        app.platform.structure.message
    );

    let found = species::run(&SpeciesJob::search(crate::sample::SPECIES_GENE_TREES))
        .expect("a species tree");
    app.open_remote_result(
        &job("j1-7", Kind::Species, Status::Done),
        Response::Species(found),
    );
    assert!(app.species_win.open);
    assert_eq!(app.species_win.result.as_ref().map(|r| r.genes), Some(40));
    assert!(
        app.species_win.message.contains("on the workstation"),
        "{}",
        app.species_win.message
    );
}

#[test]
fn a_domain_result_opens_in_its_window_and_says_it_came_from_the_workstation() {
    let mut app = App::new(Startup::default());
    let found = crate::remote::domains::DomainsResult {
        names: vec!["ubiquitin".into()],
        lengths: vec![76],
        hits: Vec::new(),
        program: "HMMER 3.4".into(),
        database: "Pfam-A.hmm".into(),
        seconds: 2.0,
    };
    app.open_remote_result(
        &job("j1-10", Kind::Domains, Status::Done),
        Response::Domains(found),
    );
    assert!(app.domains_win.open);
    assert_eq!(
        app.domains_win.result.as_ref().map(|r| r.names.len()),
        Some(1)
    );
    assert!(
        app.domains_win
            .message
            .starts_with("No domain found in 1 sequence by HMMER 3.4")
            && app.domains_win.message.ends_with("on the workstation."),
        "{}",
        app.domains_win.message
    );
}

#[test]
fn a_domain_search_is_sent_there_when_the_choice_says_so_and_its_answer_or_its_reason_comes_back() {
    let server = serving("workstation-domains");
    let (mut app, ctx) = window_for(&server);
    app.domains_win.open = true;
    app.domains_win.text =
        ">ubiquitin\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
            .into();
    let place = app.platform.place(&app.domains_win.weight());
    assert!(place.is_there(), "the choice is the workstation");
    app.domains_win.start(&ctx, place);
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    assert!(app.domains_win.running() && app.running_analyses() >= 1);
    // (a worker here: this computer may or may not have HMMER and Pfam, so the job is run, and what it says is kept either way)
    let queued = server
        .shared
        .list()
        .into_iter()
        .find(|j| j.status == Status::Queued)
        .expect("a queued job");
    let folder = server.shared.job_dir(&queued.id);
    let sent = std::fs::read(folder.join("request.bin")).expect("a request");
    let request = Request::decode(queued.kind, &Container::decode(&sent).expect("a container"))
        .expect("a request");
    let mut finished = queued.clone();
    match execute(request.clone(), &Env::default()) {
        Ok(response) => {
            std::fs::write(folder.join("result.bin"), response.encode().encode())
                .expect("a result");
            finished.status = Status::Done;
        }
        Err(why) => (finished.status, finished.error) = (Status::Failed, why),
    }
    server.shared.replace_for_tests(finished);
    assert_eq!(queued.kind, Kind::Domains);
    assert_eq!(queued.label, "Domain search of 1 sequences");
    assert!(matches!(
        request,
        Request::Domains(ref job) if job.queries.len() == 1 && job.queries[0].0 == "ubiquitin"
            && job.queries[0].1.len() == 76 && job.cutoff == crate::pipeline::domains::Cutoff::Curated
    ));
    draw_until(&mut app, &ctx, |a| !a.domains_win.running());
    // this computer may or may not have HMMER and Pfam: either the domains come back, or the reason does
    let said = &app.domains_win.message;
    assert!(
        said.contains("on the workstation") || said.contains("hmmscan") || said.contains("Pfam"),
        "{said}"
    );
}

#[test]
fn a_motif_result_opens_in_the_finder_only_if_it_is_about_the_sequences_the_window_holds() {
    let mut app = App::new(Startup::default());
    let text = crate::sample::motif_demo();
    let (_, seqs) = crate::seq::motif::read_sequences(&text).expect("the example");
    let found = crate::remote::motif::run(&crate::remote::motif::MotifJob {
        seqs,
        settings: crate::seq::motif::Settings {
            width: 10,
            ..crate::seq::motif::Settings::default()
        },
    })
    .expect("motifs");
    // the window holds other sequences: the result is not shown as if it were theirs
    app.motif_win.text = "ACGTACGTAC\nTTGACCAGGT\n".into();
    app.open_remote_result(
        &job("j1-9", Kind::Motif, Status::Done),
        Response::Motif(found.clone()),
    );
    assert!(app.motif_win.open && app.motif_win.run.is_none());
    assert!(
        app.motif_win.message.contains("no longer holds"),
        "{}",
        app.motif_win.message
    );
    // the sequences it is about: shown, and the message says where it was made
    app.motif_win.text = text;
    app.open_remote_result(
        &job("j1-9", Kind::Motif, Status::Done),
        Response::Motif(found),
    );
    assert_eq!(app.motif_win.run.as_ref().map(|r| r.found.len()), Some(1));
    assert!(
        app.motif_win.message.contains("on the workstation"),
        "{}",
        app.motif_win.message
    );
}

#[test]
fn a_motif_search_is_computed_there_when_the_choice_says_so_and_shown_when_it_comes_back() {
    let server = serving("workstation-motif");
    let (mut app, ctx) = window_for(&server);
    app.motif_win.open = true;
    app.motif_win.text = crate::sample::motif_demo();
    app.motif_win.settings.width = 10;
    app.motif_win.settings.effort = 3.0;
    app.start_motifs(&ctx);
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    assert!(app.motif_win.running() && app.running_analyses() >= 1);

    let (queued, request) = answer(&server);
    assert_eq!(queued.kind, Kind::Motif);
    assert_eq!(queued.label, "Motif search in 12 sequences");
    assert!(matches!(
        request,
        Request::Motif(ref job) if job.seqs.len() == 12 && job.settings.width == 10
            && job.settings.effort == 3.0
    ));
    draw_until(&mut app, &ctx, |a| a.motif_win.run.is_some());
    let said = &app.motif_win.message;
    assert!(said.contains("on the workstation"), "{said}");
    assert_eq!(app.motif_win.run.as_ref().map(|r| r.found.len()), Some(1));
    assert!(!app.motif_win.running());
}

#[test]
fn the_report_of_a_pipeline_file_that_the_workstation_ran_opens_on_the_page_that_runs_them() {
    let (mut app, _) = (App::new(Startup::default()), ());
    let done = crate::remote::flow::run(
        crate::remote::flow::FlowJob {
            pipeline: crate::pipeline::flow::Pipeline::loci_to_species(),
            loci: crate::pipeline::flow::tests::loci(4),
        },
        None,
    )
    .expect("a result");
    app.open_remote_result(&job("j1-8", Kind::Flow, Status::Done), Response::Flow(done));
    assert_eq!(app.workspace, Workspace::Pipelines);
    let page = &app.platform.pipelines;
    assert_eq!(page.mode, crate::ui::platform::pipelines::flow::Mode::File);
    assert_eq!(page.flow.ran_on, "the workstation");
    assert!(
        page.flow.message.contains("job j1-8"),
        "{}",
        page.flow.message
    );
    let report = page.flow.report.as_ref().expect("a report");
    assert!(matches!(
        report.output,
        Some(crate::pipeline::flow::Value::Tree(_))
    ));
}

// ---------------------------------------------------------------------------------------------------------------
// The whole window against a server on this computer. The server's workers are off, so a job waits in its queue until
// the test does what a worker would: runs the request here and stores the answer.
// ---------------------------------------------------------------------------------------------------------------

fn serving(label: &str) -> Server {
    let mut config = Config::new(TOKEN.into(), work_folder(label).expect("a folder"));
    config.gpus = vec![Gpu {
        id: 0,
        name: "GPU A".into(),
    }];
    config.tools = vec![("MAFFT".into(), "7.5".into())];
    start(config, "127.0.0.1:0".parse().expect("an address"), false).expect("a server")
}

/// Draw one frame of the whole window; the textures it asks for are of no use here, so they are let go.
fn frame(app: &mut App, ctx: &egui::Context, n: u32) {
    let mut output = ctx.run_ui(frame_input(n), |ui| app.draw(ui));
    output.textures_delta.clear();
}

/// Draw frames of the whole window until `done` says that what was waited for has happened.
fn draw_until(app: &mut App, ctx: &egui::Context, mut done: impl FnMut(&mut App) -> bool) {
    let end = Instant::now() + Duration::from_secs(60);
    for n in 0.. {
        frame(app, ctx, n);
        if done(app) {
            return;
        }
        assert!(Instant::now() < end, "gave up waiting at frame {n}");
        std::thread::sleep(Duration::from_millis(15));
    }
}

/// What a worker does with the oldest queued job: run its request here, keep the answer, and mark the job done.
fn answer(server: &Server) -> (JobInfo, Request) {
    let queued = server
        .shared
        .list()
        .into_iter()
        .find(|j| j.status == Status::Queued)
        .expect("a queued job");
    let folder = server.shared.job_dir(&queued.id);
    let sent = std::fs::read(folder.join("request.bin")).expect("a request");
    let request = Request::decode(queued.kind, &Container::decode(&sent).expect("a container"))
        .expect("a request");
    let response = execute(request.clone(), &Env::default()).expect("a result");
    std::fs::write(folder.join("result.bin"), response.encode().encode()).expect("a result");
    let mut done = queued.clone();
    done.status = Status::Done;
    server.shared.replace_for_tests(done);
    (queued, request)
}

/// A window set to send its heavy jobs to `server`.
fn window_for(server: &Server) -> (App, egui::Context) {
    let mut app = App::new(Startup::default());
    app.platform.workstation = connected(&server.addr.to_string());
    (app, egui::Context::default())
}

#[test]
fn an_alignment_asked_for_with_the_workstation_on_is_computed_there_and_shown_when_it_comes_back() {
    let server = serving("workstation-align");
    let (mut app, ctx) = window_for(&server);
    app.load_example(Workspace::Alignment);
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    assert!(app.platform.alignment.is_running());
    assert!(
        app.platform.alignment.message.contains("workstation"),
        "{}",
        app.platform.alignment.message
    );

    let (queued, request) = answer(&server);
    assert_eq!(queued.label, "Alignment of 6 sequences (built-in)");
    let Request::Align { seqs, mafft, .. } = request else {
        panic!("an alignment was asked for");
    };
    assert_eq!((seqs.len(), mafft), (6, false));
    draw_until(&mut app, &ctx, |a| a.platform.alignment.result.is_some());
    assert!(
        app.platform
            .alignment
            .message
            .starts_with("Aligned 6 sequences"),
        "{}",
        app.platform.alignment.message
    );
    // the workstation said what it has, and the page showed that while it waited
    assert_eq!(app.platform.workstation.tools().map(<[_]>::len), Some(1));
}

#[test]
fn the_progressive_method_is_asked_for_there_with_the_settings_the_page_shows() {
    let server = serving("workstation-align-progressive");
    let (mut app, ctx) = window_for(&server);
    let options = crate::seq::clustal::Options {
        guide: crate::seq::clustal::Guide::Upgma,
        free_end_gaps: false,
        refine: 2,
        fast: true,
    };
    app.platform.alignment.choose_progressive(options);
    app.load_example(Workspace::Alignment);
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    let (queued, request) = answer(&server);
    assert_eq!(queued.label, "Alignment of 6 sequences (progressive)");
    let Request::Align {
        mafft, progressive, ..
    } = request
    else {
        panic!("an alignment was asked for");
    };
    assert!(!mafft);
    assert_eq!(progressive, Some(options));
    draw_until(&mut app, &ctx, |a| a.platform.alignment.result.is_some());
    let result = app
        .platform
        .alignment
        .result
        .as_ref()
        .expect("an alignment");
    assert!(
        result.engine.contains("UPGMA guide tree") && result.engine.contains("refined"),
        "{}",
        result.engine
    );
}

#[test]
fn a_workstation_that_does_not_know_the_progressive_method_is_told_from_one_that_does() {
    let server = serving("workstation-align-older");
    let (mut app, ctx) = window_for(&server);
    app.platform
        .alignment
        .choose_progressive(crate::seq::clustal::Options::default());
    app.load_example(Workspace::Alignment);
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    // an earlier version reads the request without the setting it does not know, and runs the built-in method
    let queued = server
        .shared
        .list()
        .into_iter()
        .find(|j| j.status == Status::Queued)
        .expect("a queued job");
    let folder = server.shared.job_dir(&queued.id);
    let sent = std::fs::read(folder.join("request.bin")).expect("a request");
    let Request::Align { seqs, .. } =
        Request::decode(queued.kind, &Container::decode(&sent).expect("a container"))
            .expect("a request")
    else {
        panic!("an alignment was asked for");
    };
    let older = Request::Align {
        seqs,
        mafft: false,
        progressive: None,
    };
    let response = execute(older, &Env::default()).expect("a result");
    std::fs::write(folder.join("result.bin"), response.encode().encode()).expect("a result");
    let mut done = queued.clone();
    done.status = Status::Done;
    server.shared.replace_for_tests(done);
    draw_until(&mut app, &ctx, |a| a.platform.alignment.result.is_some());
    let said = &app.platform.alignment.message;
    assert!(
        said.contains("progressive method was asked for") && said.contains("earlier version"),
        "{said}"
    );
}

#[test]
fn a_tree_is_estimated_there_when_the_alignment_is_ready() {
    let server = serving("workstation-tree");
    let (mut app, ctx) = window_for(&server);
    let aligned = execute(
        Request::Align {
            seqs: isolates(),
            mafft: false,
            progressive: None,
        },
        &Env::default(),
    )
    .and_then(Response::into_align)
    .expect("aligned");
    app.platform.alignment.set_result(aligned);
    app.platform.phylogeny.run_when_ready = true;
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());

    let (queued, request) = answer(&server);
    assert_eq!(queued.kind, Kind::Tree);
    assert!(
        queued.label.starts_with("Tree (") && queued.label.ends_with("of 6 sequences"),
        "{}",
        queued.label
    );
    assert!(matches!(request, Request::Tree { ref alignment, .. } if alignment.len() == 6));
    draw_until(&mut app, &ctx, |a| a.tree_win.tree.is_some());
    assert!(
        app.platform
            .phylogeny
            .message
            .starts_with("Estimated a tree of 6 tips"),
        "{}",
        app.platform.phylogeny.message
    );
}

#[test]
fn a_pipeline_waits_for_the_workstation_to_say_what_it_has_and_then_runs_there() {
    let server = serving("workstation-pipeline");
    let (mut app, ctx) = window_for(&server);
    app.platform.pipelines.use_installed = false;
    app.load_example(Workspace::Pipelines);
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    assert!(
        app.platform.workstation.tools().is_some(),
        "the run started only once the programs were known"
    );

    let (queued, request) = answer(&server);
    assert_eq!(queued.label, "Pipeline on 6 sequences");
    assert!(
        matches!(request, Request::Pipeline { ref choice, target: Some((200, 60)), .. } if choice.align && choice.tree && choice.primers)
    );
    draw_until(&mut app, &ctx, |a| a.platform.pipelines.outcome.is_some());
    let outcome = app.platform.pipelines.outcome.as_ref().expect("an outcome");
    assert!(outcome.aligned.is_some() && outcome.tree.is_some() && outcome.primers.is_some());
}

#[test]
fn phylocnn_estimates_are_computed_there_and_its_gpu_is_the_one_offered() {
    // the networks are the authors' files, which this repository does not carry
    if crate::models::phylocnn_missing() {
        return;
    }
    let server = serving("workstation-cnn");
    let (mut app, ctx) = window_for(&server);
    app.phylocnn_win.open = true;
    app.load_phylocnn_example();
    app.run_phylocnn();
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());

    let (queued, request) = answer(&server);
    assert_eq!(queued.kind, Kind::PhyloCnn);
    assert!(matches!(
        request,
        Request::PhyloCnn {
            task: Task::Bdss,
            priors: Priors::Hiv,
            gpu: false,
            ..
        }
    ));
    draw_until(&mut app, &ctx, |a| !a.phylocnn_win.table.is_empty());
    assert_eq!(app.platform.workstation.gpu_name(), Some("GPU A"));
}

#[test]
fn a_job_that_failed_there_says_why_and_a_missing_workstation_is_said_too() {
    let server = serving("workstation-fail");
    let (mut app, ctx) = window_for(&server);
    app.load_example(Workspace::Alignment);
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    let mut failed = server.shared.list().remove(0);
    failed.status = Status::Failed;
    failed.error = "MAFFT stopped with exit status: 1.".into();
    server.shared.replace_for_tests(failed);
    draw_until(&mut app, &ctx, |a| !a.platform.alignment.is_running());
    assert_eq!(
        app.platform.alignment.message,
        "MAFFT stopped with exit status: 1."
    );
    assert!(app.platform.alignment.result.is_none());

    // a workstation that is not there: the run says so rather than falling back to this computer
    let mut app = App::new(Startup::default());
    app.platform.workstation = connected("127.0.0.1:1");
    let ctx = egui::Context::default();
    app.load_example(Workspace::Alignment);
    draw_until(&mut app, &ctx, |a| !a.platform.alignment.is_running());
    let said = &app.platform.alignment.message;
    assert!(
        said.contains("127.0.0.1:1") || said.contains("connect"),
        "{said}"
    );
    assert!(app.platform.alignment.result.is_none());
}

#[test]
fn with_the_choice_on_automatic_a_small_alignment_stays_here_and_a_big_one_goes_to_the_workstation()
{
    let server = serving("workstation-automatic");
    let (mut app, ctx) = window_for(&server);
    app.platform.workstation.set_compute(Compute::Automatic);

    // six isolates are a small job: they are aligned on this computer and nothing is sent
    app.load_example(Workspace::Alignment);
    draw_until(&mut app, &ctx, |a| a.platform.alignment.result.is_some());
    assert!(server.shared.list().is_empty(), "a small job was sent");
    assert!(
        !app.platform.alignment.message.contains("workstation"),
        "{}",
        app.platform.alignment.message
    );

    // one more than the size at which a job is big
    let page = &mut app.platform.alignment;
    page.inputs.clear();
    for k in 0..=HEAVY_SEQUENCES {
        page.add(&format!("s{k}"), b"ACGTACGTACGTACGT");
    }
    assert!(page.weight().is_heavy());
    app.start_alignment(&ctx);
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    let (queued, request) = answer(&server);
    assert_eq!(queued.kind, Kind::Align);
    assert_eq!(queued.label, "Alignment of 101 sequences (built-in)");
    assert!(matches!(request, Request::Align { ref seqs, .. } if seqs.len() == 101));
    draw_until(&mut app, &ctx, |a| {
        a.platform
            .alignment
            .result
            .as_ref()
            .is_some_and(|r| r.rows.len() == 101)
    });
}

#[test]
fn with_the_choice_on_automatic_ticking_the_gpu_sends_even_one_tree_to_the_workstation() {
    // the networks are the authors' files, which this repository does not carry
    if crate::models::phylocnn_missing() {
        return;
    }
    let server = serving("workstation-automatic-gpu");
    let (mut app, ctx) = window_for(&server);
    app.platform.workstation.set_compute(Compute::Automatic);
    app.platform.workstation.status = Some(status(&[]));
    app.phylocnn_win.open = true;
    app.load_phylocnn_example();

    // one tree on the processor is a small job: it runs here
    app.run_phylocnn();
    draw_until(&mut app, &ctx, |a| !a.phylocnn_win.table.is_empty());
    assert!(server.shared.list().is_empty(), "a small job was sent");

    // asking for the workstation's GPU is what sends it there
    app.phylocnn_win.use_gpu = true;
    app.run_phylocnn();
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    let (queued, request) = answer(&server);
    assert_eq!(queued.kind, Kind::PhyloCnn);
    assert!(matches!(request, Request::PhyloCnn { gpu: true, .. }));
}

#[test]
fn a_species_tree_is_computed_there_when_the_choice_says_so_and_shown_when_it_comes_back() {
    let server = serving("workstation-species");
    let (mut app, ctx) = window_for(&server);
    app.species_win.open = true;
    app.load_species_example();
    app.run_species();
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    assert!(app.species_win.is_running() && app.running_analyses() == 1);

    let (queued, request) = answer(&server);
    assert_eq!(queued.kind, Kind::Species);
    assert_eq!(queued.label, "Species tree from 40 gene trees");
    assert!(matches!(
        request,
        Request::Species(ref job) if job.genes.matches(';').count() == 40
            && job.score_tree.is_empty() && job.extra == 1 && !job.exact
    ));
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let said = &app.species_win.message;
    assert!(said.contains("on the workstation"), "{said}");
    assert_eq!(
        app.species_win.result.as_ref().map(|r| r.branches.len()),
        Some(5)
    );
    assert!(!app.species_win.is_running());
}

/// A folder of `n` loci of six taxa each, as the files of an alignment, removed when the test ends.
struct Loci(std::path::PathBuf);

impl Loci {
    fn new(label: &str, n: u64) -> Loci {
        let dir = work_folder(label).expect("a folder");
        for (name, value) in crate::pipeline::flow::tests::loci(n) {
            std::fs::write(dir.join(format!("{name}.fa")), value.to_text()).expect("a file");
        }
        Loci(dir)
    }
}

impl Drop for Loci {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

#[test]
fn a_pipeline_file_is_run_there_when_the_choice_says_so_and_its_report_comes_back_to_its_page() {
    let server = serving("workstation-flow");
    let (mut app, ctx) = window_for(&server);
    let loci = Loci::new("workstation-flow-loci", 4);
    app.platform.pipelines.flow.folder = Some(loci.0.clone());
    app.start_flow(&ctx);
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    assert!(app.platform.pipelines.flow.is_running() && app.running_analyses() == 1);

    let (queued, request) = answer(&server);
    assert_eq!(queued.kind, Kind::Flow);
    assert!(queued.label.ends_with("on 4 files"), "{}", queued.label);
    assert!(matches!(
        request,
        Request::Flow(ref job) if job.loci.len() == 4
            && job.pipeline == crate::pipeline::flow::Pipeline::loci_to_species()
    ));
    draw_until(&mut app, &ctx, |a| {
        a.platform.pipelines.flow.report.is_some()
    });
    let page = &app.platform.pipelines.flow;
    assert_eq!(page.ran_on, "the workstation");
    assert!(page.message.starts_with("5 steps run"), "{}", page.message);
    let report = page.report.as_ref().expect("a report");
    assert!(matches!(
        report.output,
        Some(crate::pipeline::flow::Value::Tree(_))
    ));
    assert!(!page.is_running());
}

#[test]
fn with_the_choice_on_automatic_a_small_pipeline_file_stays_here_and_one_over_many_files_goes_there(
) {
    let server = serving("workstation-flow-automatic");
    let (mut app, ctx) = window_for(&server);
    app.platform.workstation.set_compute(Compute::Automatic);
    let few = Loci::new("workstation-flow-few", 3);
    app.platform.pipelines.flow.folder = Some(few.0.clone());
    app.start_flow(&ctx);
    draw_until(&mut app, &ctx, |a| {
        a.platform.pipelines.flow.report.is_some()
    });
    assert!(server.shared.list().is_empty(), "a small run was sent");
    assert_eq!(app.platform.pipelines.flow.ran_on, "this computer");

    // more files than the trees at which a run is big
    let many = Loci::new("workstation-flow-many", HEAVY_TREES as u64 + 1);
    app.platform.pipelines.flow.folder = Some(many.0.clone());
    app.start_flow(&ctx);
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    let (queued, _) = answer(&server);
    assert_eq!(queued.kind, Kind::Flow);
    draw_until(&mut app, &ctx, |a| {
        a.platform.pipelines.flow.report.is_some()
    });
    assert_eq!(app.platform.pipelines.flow.ran_on, "the workstation");
}

#[test]
fn a_weighted_species_search_is_computed_there_with_its_weights_and_shown_as_weighted() {
    let server = serving("workstation-species-weighted");
    let (mut app, ctx) = window_for(&server);
    app.species_win.open = true;
    app.species_win.genes = "(A:0.1,B:0.1,(C:0.1,D:0.1)0.95:0.1);\n".repeat(3)
        + &"(A:0.1,C:0.1,(B:0.1,D:0.1)0.05:0.1);\n".repeat(4);
    app.species_win.weights = crate::seq::astral::Weights::Hybrid;
    app.run_species();
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    let (queued, request) = answer(&server);
    assert_eq!(queued.kind, Kind::Species);
    assert!(matches!(
        request,
        Request::Species(ref job) if job.weights == crate::seq::astral::Weights::Hybrid
    ));
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let result = app.species_win.result.clone().expect("a result");
    assert_eq!(result.weights, crate::seq::astral::Weights::Hybrid);
    assert!(result.weighted_score.is_some());
    let said = &app.species_win.message;
    assert!(
        said.contains("on the workstation") && said.contains("weighted by both"),
        "{said}"
    );
}

#[test]
fn a_bootstrap_is_computed_there_with_its_replicates_and_shown_with_its_support() {
    let server = serving("workstation-species-bootstrap");
    let (mut app, ctx) = window_for(&server);
    app.species_win.open = true;
    app.species_win.genes = "((A,B),(C,D));\n".repeat(6);
    app.species_win.replicates = 6;
    app.species_win.seed = 3;
    app.run_species();
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    let (queued, request) = answer(&server);
    assert_eq!(queued.kind, Kind::Species);
    assert!(matches!(request, Request::Species(ref job) if job.replicates == 6 && job.seed == 3));
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let result = app.species_win.result.clone().expect("a result");
    assert_eq!(result.branches[0].bootstrap, Some(100.0));
    let said = &app.species_win.message;
    assert!(
        said.contains("on the workstation") && said.contains("Bootstrap support from 6 of 6"),
        "{said}"
    );
}

#[test]
fn gene_families_are_computed_there_when_asked_for_and_shown_as_families() {
    let server = serving("workstation-species-families");
    let (mut app, ctx) = window_for(&server);
    app.species_win.open = true;
    app.species_win.genes =
        "(((A,B),(C,(D,E))),((A,B),(C,(D,E))));\n(((A,B),(C,(D,E))),((A,B),(D,E)));\n".into();
    app.species_win.families = true;
    app.run_species();
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    let (queued, request) = answer(&server);
    assert_eq!(queued.kind, Kind::Species);
    assert!(matches!(request, Request::Species(ref job) if job.families));
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    let result = app.species_win.result.clone().expect("a result");
    assert!(result.families && result.score > 0);
    assert!(
        app.species_win.message.contains("on the workstation"),
        "{}",
        app.species_win.message
    );
}

#[test]
fn with_the_choice_on_automatic_a_small_species_search_stays_here_and_a_big_one_goes_there() {
    let server = serving("workstation-species-automatic");
    let (mut app, ctx) = window_for(&server);
    app.platform.workstation.set_compute(Compute::Automatic);
    app.species_win.open = true;
    app.load_species_example();
    app.run_species();
    draw_until(&mut app, &ctx, |a| a.species_win.result.is_some());
    assert!(server.shared.list().is_empty(), "a small job was sent");
    assert!(app.species_win.message.contains("on this computer"));

    // 170 copies of the example have just over the leaves at which a search is big
    let text = crate::sample::SPECIES_GENE_TREES.repeat(170);
    assert!(Weight::of_species(&text).is_heavy());
    app.load_species_file(Source::Bytes {
        name: "many.tre",
        bytes: text.as_bytes(),
    });
    app.run_species();
    draw_until(&mut app, &ctx, |_| !server.shared.list().is_empty());
    let (queued, _) = answer(&server);
    assert_eq!(queued.kind, Kind::Species);
    draw_until(&mut app, &ctx, |a| {
        a.species_win
            .result
            .as_ref()
            .is_some_and(|r| r.genes == 6_800)
    });
    assert!(app.species_win.message.contains("on the workstation"));
}

#[test]
fn every_page_that_starts_a_job_draws_the_choice_in_every_mode_with_and_without_a_workstation() {
    let ctx = egui::Context::default();
    for dark in [false, true] {
        for connected_to in [None, Some("127.0.0.1:1")] {
            for mode in Compute::ALL {
                let mut app = App::new(Startup {
                    dark,
                    ..Startup::default()
                });
                if let Some(address) = connected_to {
                    app.platform.workstation = connected(address);
                    app.platform.workstation.set_compute(mode);
                }
                app.phylocnn_win.open = true;
                for (n, workspace) in [
                    Workspace::Alignment,
                    Workspace::Phylogeny,
                    Workspace::Pipelines,
                    Workspace::Structure,
                    Workspace::Workstation,
                    Workspace::Home,
                ]
                .into_iter()
                .enumerate()
                {
                    app.workspace = workspace;
                    frame(&mut app, &ctx, n as u32);
                }
            }
        }
    }
}

#[test]
fn the_page_draws_in_every_state_in_both_themes() {
    let ctx = egui::Context::default();
    for dark in [false, true] {
        let mut app = App::new(Startup {
            dark,
            ..Startup::default()
        });
        app.workspace = Workspace::Workstation;
        // not connected: the set-up steps
        for n in 0..2 {
            frame(&mut app, &ctx, n);
        }
        // connected, with a status, a busy GPU and one job of every kind of status
        let mut page = connected("127.0.0.1:1");
        let mut s = status(&[("MAFFT", "7.5"), ("curl", "")]);
        s.fold = true;
        s.gpus.push(GpuState {
            id: 1,
            name: "NVIDIA B".into(),
            job: Some("j1-2".into()),
        });
        let mut running = job("j1-2", Kind::PhyloCnn, Status::Running);
        (running.gpus, running.gpus_wanted) = (vec![1], 1);
        let mut failed = job("j1-4", Kind::Fold, Status::Failed);
        failed.error = "the program stopped".into();
        page.status = Some(s);
        page.jobs = vec![
            job("j1-1", Kind::Align, Status::Done),
            running,
            job("j1-3", Kind::Tree, Status::Queued),
            failed,
            job("j1-5", Kind::Pipeline, Status::Cancelled),
        ];
        page.problem = "Could not refresh: timed out".into();
        page.looked = Some(Instant::now());
        app.platform.workstation = page;
        for n in 2..5 {
            frame(&mut app, &ctx, n);
        }
        // a workstation that speaks another protocol is warned about
        if let Some(s) = app.platform.workstation.status.as_mut() {
            s.protocol = PROTOCOL + 1;
        }
        frame(&mut app, &ctx, 5);
    }
}

/// Where `text` was drawn in a frame: the first piece of text that holds it.
fn text_rect(shapes: &[egui::epaint::Shape], text: &str) -> Option<egui::Rect> {
    shapes.iter().find_map(|shape| match shape {
        egui::epaint::Shape::Text(t) if t.galley.text().contains(text) => {
            Some(t.visual_bounding_rect())
        }
        egui::epaint::Shape::Vec(inner) => text_rect(inner, text),
        _ => None,
    })
}

/// One frame of the whole window with `events`, numbered by a counter that goes up with each (the time of a frame must never
/// go back).
fn frame_with(
    app: &mut App,
    ctx: &egui::Context,
    n: &mut u32,
    events: Vec<egui::Event>,
) -> egui::FullOutput {
    let mut input = frame_input(*n);
    input.events = events;
    *n += 1;
    let mut output = ctx.run_ui(input, |ui| app.draw(ui));
    output.textures_delta.clear();
    output
}

/// A page connected to nothing in particular, set to send every job to the workstation, which has said that it has `gpus`
/// GPUs.
fn page_with_gpus(dir: PathBuf, gpus: u8) -> WorkstationPage {
    let mut page = WorkstationPage {
        address: "127.0.0.1:1".into(),
        token: TOKEN.into(),
        ..WorkstationPage::open(Some(dir))
    };
    page.connect();
    page.set_compute(Compute::Workstation);
    let mut s = status(&[]);
    for id in 1..gpus {
        s.gpus.push(GpuState {
            id: usize::from(id),
            name: format!("NVIDIA {id}"),
            job: None,
        });
    }
    page.status = Some(s);
    // it has just looked: it does not try the closed port while the test draws
    page.looked = Some(Instant::now());
    page
}

#[test]
fn the_choice_of_gpus_for_one_job_is_on_the_page_only_when_the_workstation_has_several() {
    let ctx = egui::Context::default();
    let dir = work_folder("workstation-gpus-shown").expect("a folder");
    for (gpus, shown) in [(1, false), (2, true), (8, true)] {
        let mut app = App::new(Startup::default());
        app.workspace = Workspace::Workstation;
        app.platform.workstation = page_with_gpus(dir.clone(), gpus);
        frame(&mut app, &ctx, u32::from(gpus) * 10);
        let mut input = frame_input(u32::from(gpus) * 10 + 1);
        input.events.clear();
        let drawn = ctx.run_ui(input, |ui| app.draw(ui));
        let shapes: Vec<egui::epaint::Shape> =
            drawn.shapes.iter().map(|s| s.shape.clone()).collect();
        // a heading is drawn in capitals
        let found = text_rect(&shapes, "GPUS FOR ONE JOB").is_some();
        assert_eq!(found, shown, "{gpus} GPUs");
    }
}

#[test]
fn a_click_on_the_gpus_choice_is_kept_and_the_jobs_sent_afterwards_ask_for_that_many() {
    use egui::{pos2, Event, Modifiers, PointerButton};
    let dir = work_folder("workstation-gpus-click").expect("a folder");
    let ctx = egui::Context::default();
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Workstation;
    app.platform.workstation = page_with_gpus(dir.clone(), 4);
    assert_eq!(app.platform.workstation.gpus, GpusPerJob::Count(1));

    let mut n = 0;
    let centre_of = |shown: &egui::FullOutput, text: &str| {
        let shapes: Vec<egui::epaint::Shape> =
            shown.shapes.iter().map(|s| s.shape.clone()).collect();
        let rect = text_rect(&shapes, text).unwrap_or_else(|| panic!("{text} was not drawn"));
        pos2(rect.center().x, rect.center().y)
    };
    let click = |app: &mut App, n: &mut u32, at: egui::Pos2| {
        let button = |pressed| Event::PointerButton {
            pos: at,
            button: PointerButton::Primary,
            pressed,
            modifiers: Modifiers::NONE,
        };
        frame_with(app, &ctx, n, vec![Event::PointerMoved(at)]);
        frame_with(app, &ctx, n, vec![button(true), button(false)])
    };

    frame_with(&mut app, &ctx, &mut n, vec![]);
    let shown = frame_with(&mut app, &ctx, &mut n, vec![]);
    // the choice shows what it is now; a click opens its list, and a click on a number in the list chooses it
    let at = centre_of(&shown, "1 GPU");
    click(&mut app, &mut n, at);
    let opened = frame_with(&mut app, &ctx, &mut n, vec![]);
    let at = centre_of(&opened, "3 GPUs");
    click(&mut app, &mut n, at);
    frame_with(&mut app, &ctx, &mut n, vec![]);
    assert_eq!(app.platform.workstation.gpus, GpusPerJob::Count(3));

    // it is kept with the connection, and the jobs that go there ask for it
    assert_eq!(WorkstationPage::open(Some(dir)).gpus, GpusPerJob::Count(3));
    let Place::There(remote) = app.platform.workstation.place(&Weight::default()) else {
        panic!("every job goes to the workstation");
    };
    let plain = Remote::new("127.0.0.1:1", TOKEN).expect("a connection");
    assert_eq!(remote, plain.with_gpus(GpusPerJob::Count(3)));
}

#[test]
fn the_choices_are_one_gpu_up_to_all_the_workstations_and_then_every_gpu_it_has() {
    assert_eq!(
        compute::gpu_choices(3),
        [
            GpusPerJob::Count(1),
            GpusPerJob::Count(2),
            GpusPerJob::Count(3),
            GpusPerJob::All
        ]
    );
    assert_eq!(compute::gpu_choices(0), [GpusPerJob::All]);
    // never more than a job can be asked for
    assert_eq!(
        compute::gpu_choices(1000).len(),
        crate::remote::settings::MOST_GPUS_PER_JOB + 1
    );
    // each choice says what it does
    for choice in compute::gpu_choices(4) {
        assert!(!compute::explain_gpus(choice).is_empty(), "{choice:?}");
    }
    assert!(compute::explain_gpus(GpusPerJob::All).contains("every GPU"));
    assert!(compute::explain_gpus(GpusPerJob::Count(1)).contains("side by side"));
    assert!(compute::explain_gpus(GpusPerJob::Count(4)).contains("4 at once"));
}

#[test]
fn forgetting_the_workstation_forgets_the_gpus_choice_with_it() {
    let dir = work_folder("workstation-gpus-forget").expect("a folder");
    let mut page = page_with_gpus(dir.clone(), 2);
    page.set_gpus(GpusPerJob::All);
    assert_eq!(
        WorkstationPage::open(Some(dir.clone())).gpus,
        GpusPerJob::All
    );
    page.forget();
    assert_eq!(page.gpus, GpusPerJob::Count(1));
    assert_eq!(WorkstationPage::open(Some(dir)).gpus, GpusPerJob::Count(1));
}

#[test]
fn the_button_of_an_analysis_is_pressed_once_ready_never_while_it_runs_and_a_wheel_turns_beside_it()
{
    use eframe::egui::epaint::Shape;
    use eframe::egui::{Event, Modifiers, PointerButton, Pos2};
    fn walk(shape: &Shape, see: &mut dyn FnMut(&Shape)) {
        match shape {
            Shape::Vec(v) => v.iter().for_each(|s| walk(s, see)),
            other => see(other),
        }
    }
    let look = design::look(false);
    let view = ComputeView {
        mode: Compute::Local,
        connected: false,
        there: false,
    };
    // (where the words of the button are, whether a wheel is drawn, whether a click on the button pressed it)
    let tried = |state: (bool, bool)| -> (bool, bool) {
        let ctx = egui::Context::default();
        let mut pressed = false;
        let mut n = 0;
        let mut frame = |events: Vec<Event>| {
            let mut input = crate::ui::headless::frame_input(n);
            n += 1;
            input.events = events;
            let mut out = ctx.run_ui(input, |ui| {
                pressed |= super::run_button(ui, &look, view, "Go on", state);
            });
            out.textures_delta.clear();
            out.shapes
        };
        let shapes = frame(vec![]);
        let (mut at, mut wheel) = (None::<Pos2>, false);
        for clipped in &shapes {
            walk(&clipped.shape, &mut |s| match s {
                Shape::Text(t) if t.galley.text() == "Go on" => {
                    at = Some(t.pos + t.galley.rect.center().to_vec2());
                }
                Shape::Path(_) | Shape::LineSegment { .. } => wheel = true,
                _ => {}
            });
        }
        let at = at.expect("the button says what it does");
        let button = |down| Event::PointerButton {
            pos: at,
            button: PointerButton::Primary,
            pressed: down,
            modifiers: Modifiers::NONE,
        };
        frame(vec![Event::PointerMoved(at)]);
        frame(vec![button(true)]);
        frame(vec![button(false)]);
        frame(vec![]);
        (wheel, pressed)
    };
    assert_eq!(
        tried((true, false)),
        (false, true),
        "ready: no wheel, and a click runs it"
    );
    assert_eq!(
        tried((false, false)),
        (false, false),
        "not ready: a click does nothing"
    );
    assert_eq!(
        tried((true, true)),
        (true, false),
        "running: a wheel, and a click does nothing"
    );
}

/// A table that stands where the storage a browser gives one tab does: it is that tab's alone, and is gone with it.
fn a_tab() -> std::rc::Rc<crate::ui::session::Memory> {
    std::rc::Rc::new(crate::ui::session::Memory::default())
}

/// The names a tab's table holds.
fn kept_in(tab: &crate::ui::session::Memory) -> Vec<String> {
    crate::ui::session::Shelf::names(tab)
}

/// What a tab's table holds of the connection.
fn entry_of(tab: &crate::ui::session::Memory) -> Option<String> {
    crate::ui::session::Shelf::get(tab, settings::FILE)
}

/// A page of the tab with `address` and the token typed, and Connect pressed.
fn connected_in(tab: &std::rc::Rc<crate::ui::session::Memory>, address: &str) -> WorkstationPage {
    let mut page = WorkstationPage {
        address: address.into(),
        token: TOKEN.into(),
        ..WorkstationPage::on_tab(tab.clone())
    };
    page.connect();
    assert!(page.problem.is_empty(), "{}", page.problem);
    page
}

#[test]
fn in_a_tab_the_connection_is_kept_for_the_tab_and_comes_back_when_the_page_is_loaded_again() {
    let tab = a_tab();
    let fresh = WorkstationPage::on_tab(tab.clone());
    assert!(fresh.active.is_none() && fresh.token.is_empty() && kept_in(&tab).is_empty());
    let mut page = connected_in(&tab, "lab-gpu");
    assert_eq!(
        page.message, "Connecting…",
        "it was kept: nothing to say of that"
    );
    page.set_compute(Compute::Workstation);
    // the tab's table holds it, under one name of its own
    assert_eq!(kept_in(&tab), [settings::FILE]);
    let kept = entry_of(&tab).expect("kept");
    assert!(
        kept.contains("address=lab-gpu:7878") && kept.contains("compute=workstation"),
        "the address and the choice are kept"
    );
    assert!(kept.contains(&format!("token={TOKEN}")), "and the token");

    // the page is loaded again in the same tab: it is connected as it was, and sends its jobs where it did
    let mut again = WorkstationPage::on_tab(tab.clone());
    assert_eq!(
        (again.address.as_str(), again.token.as_str()),
        ("lab-gpu:7878", TOKEN)
    );
    assert!(again.active.is_some() && again.compute == Compute::Workstation && goes_there(&again));
    // the other choices are kept by the tab as they are made, each of them
    again.set_compute(Compute::Automatic);
    again.set_gpus(GpusPerJob::Count(2));
    let third = WorkstationPage::on_tab(tab.clone());
    assert_eq!(
        (third.compute, third.gpus),
        (Compute::Automatic, GpusPerJob::Count(2))
    );

    // another tab, or this page once its tab was closed, has a table of its own, and that one is empty
    let other = WorkstationPage::on_tab(a_tab());
    assert!(other.active.is_none() && other.token.is_empty() && other.compute == Compute::Local);

    // Forget takes it out of the tab's table at once: nothing is left to load
    let mut page = third;
    page.forget();
    assert!(kept_in(&tab).is_empty(), "{:?}", kept_in(&tab));
    let cleared = WorkstationPage::on_tab(tab.clone());
    assert!(
        cleared.active.is_none() && cleared.token.is_empty() && cleared.compute == Compute::Local
    );
    // and the page that forgot still keeps what it is told next
    page.address = "another-gpu".into();
    page.token = TOKEN.into();
    page.connect();
    assert!(entry_of(&tab).is_some_and(|kept| kept.contains("address=another-gpu:7878")));
}

#[test]
fn a_tab_whose_storage_takes_nothing_still_connects_and_says_that_it_will_ask_again() {
    let tab = a_tab();
    tab.room.set(Some(0));
    let page = connected_in(&tab, "lab-gpu");
    assert!(page.active.is_some());
    assert_eq!(page.message, NOT_KEPT_FOR_THE_TAB);
    assert!(
        !NOT_KEPT_FOR_THE_TAB.contains("Connected"),
        "the workstation has not answered yet: the page does not say that it has"
    );
    assert!(kept_in(&tab).is_empty());
}

#[test]
fn a_connection_that_could_not_be_kept_does_not_leave_the_one_before_it_to_come_back() {
    let tab = a_tab();
    let _first = connected_in(&tab, "the-first-gpu");
    assert!(entry_of(&tab).is_some_and(|kept| kept.contains("the-first-gpu")));
    // the table takes no more of this entry, and the person connects to another workstation
    *tab.refused.borrow_mut() = Some(settings::FILE.to_string());
    let second = connected_in(&tab, "the-second-gpu");
    assert_eq!(second.message, NOT_KEPT_FOR_THE_TAB);
    assert_eq!(
        second.active.as_ref().map(Remote::address),
        Some("the-second-gpu:7878")
    );
    // the first is not what the page comes back to
    assert_eq!(entry_of(&tab), None);
    let again = WorkstationPage::on_tab(tab);
    assert!(again.active.is_none() && again.token.is_empty());
}

#[test]
fn an_entry_of_a_tab_that_makes_no_connection_is_removed_and_nothing_of_it_is_shown() {
    // not what the page kept: it keeps only what made a connection. A token with a space in it makes none
    for unusable in [
        "address=lab-gpu\ntoken=two words\ncompute=workstation\n",
        "address=\ntoken=abcdefghijklmnopqrstuvwx\ncompute=automatic\n",
        "nothing that can be read",
    ] {
        let tab = a_tab();
        crate::ui::session::Shelf::set(tab.as_ref(), settings::FILE, unusable).expect("set");
        let page = WorkstationPage::on_tab(tab.clone());
        assert!(page.active.is_none(), "{unusable:?}");
        assert!(
            page.token.is_empty() && page.address.is_empty(),
            "{unusable:?}"
        );
        assert_eq!(page.compute, Compute::Local, "{unusable:?}");
        assert_eq!(entry_of(&tab), None, "{unusable:?}");
    }
    // a folder's file is the person's own, and is left as they wrote it
    let dir = work_folder("workstation-unusable").expect("a folder");
    let text = "address=lab-gpu\ntoken=two words\ncompute=workstation\n";
    std::fs::write(dir.join(settings::FILE), text).expect("a file");
    let page = WorkstationPage::open(Some(dir.clone()));
    assert!(page.active.is_none() && page.compute == Compute::Local);
    assert_eq!(
        std::fs::read_to_string(dir.join(settings::FILE)).expect("still there"),
        text
    );
}

#[test]
fn the_page_says_what_is_true_of_where_the_token_is_kept() {
    // a tab: kept for the tab, and how to be rid of it; nothing that a browser does not hold to
    let tab = privacy_text(true, true);
    assert!(
        tab.contains("for this tab") && tab.contains("Forget removes them"),
        "{tab}"
    );
    assert!(
        tab.contains("duplicated") && tab.contains("brings the tab back"),
        "{tab}"
    );
    assert!(
        !tab.contains("until this tab is closed"),
        "a browser may bring a closed tab back"
    );
    // whether the page is in a browser does not decide it: what the page keeps does
    assert_eq!(privacy_text(true, false), tab);
    // a browser that gives the page no storage, and the program on a computer
    assert!(privacy_text(false, true).contains("it is not kept"));
    assert!(privacy_text(false, false).contains("not encrypted"));
    assert!(!privacy_text(false, false).contains("tab"));
    for words in [tab, privacy_text(false, true), privacy_text(false, false)] {
        assert!(words.contains("SSH tunnel"), "{words}");
    }
    assert!(forget_hint(true).contains("this tab") && forget_hint(false).contains("this computer"));
}

#[test]
fn the_page_of_a_tab_shows_those_words_and_never_the_token() {
    use crate::ui::app::App;
    use crate::ui::platform::Workspace;
    use crate::ui::testing::Screen;
    let said_of = |page: WorkstationPage| {
        let mut app = App::new(Startup::default());
        app.platform.workstation = WorkstationPage {
            address: "lab-gpu".into(),
            token: TOKEN.into(),
            ..page
        };
        app.workspace = Workspace::Workstation;
        let mut screen = Screen::new();
        screen.settle(&mut app);
        screen.text()
    };
    // the words of a paragraph are drawn in rows: one that is short enough to stand in one row is looked for
    let in_a_tab = said_of(WorkstationPage::on_tab(a_tab()));
    assert!(in_a_tab.contains("for this tab"), "{in_a_tab}");
    assert!(!in_a_tab.contains(TOKEN), "the token is typed unseen");
    // where nothing is kept, the page does not say that anything is
    let nowhere = said_of(WorkstationPage::open(None));
    assert!(!nowhere.contains("for this tab"), "{nowhere}");
    assert!(!nowhere.contains(TOKEN));
}
