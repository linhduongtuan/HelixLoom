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
        std::fs::write(dir.join("workstation.txt"), text).expect("a file");
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
    app.open_remote_result(&job("j1-2", Kind::Tree, Status::Done), tree);
    assert_eq!(app.workspace, Workspace::Trees);
    assert!(app.tree_win.tree.is_some());

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
    app.open_remote_result(&job("j1-3", Kind::Pipeline, Status::Done), run);
    assert_eq!(app.workspace, Workspace::Pipelines);
    assert!(app.platform.pipelines.outcome.is_some());

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
    let Request::Align { seqs, mafft } = request else {
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
fn a_tree_is_estimated_there_when_the_alignment_is_ready() {
    let server = serving("workstation-tree");
    let (mut app, ctx) = window_for(&server);
    let aligned = execute(
        Request::Align {
            seqs: isolates(),
            mafft: false,
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
