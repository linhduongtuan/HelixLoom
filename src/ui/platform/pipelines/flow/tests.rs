use super::*;
use crate::pipeline::flow::tests::locus;
use crate::pipeline::flow::Value;
use crate::ui::app::{App, Startup};
use crate::ui::testing::Screen;
use std::sync::atomic::{AtomicU32, Ordering};

static NEXT: AtomicU32 = AtomicU32::new(0);

struct Folder(PathBuf);

impl Folder {
    fn new() -> Folder {
        let path = std::env::temp_dir().join(format!(
            "helixloom-flowpage-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(path.join("loci")).expect("a folder");
        Folder(path)
    }

    fn loci(&self, n: u64) {
        for k in 1..=n {
            let text = Value::Alignment(locus(k * 17, 150)).to_text();
            std::fs::write(self.0.join(format!("loci/gene{k}.fasta")), text).unwrap();
        }
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// An app on the Pipelines page in its file mode, with a folder of three loci and a cache.
fn page(f: &Folder) -> (App, Screen) {
    f.loci(3);
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Pipelines;
    let p = &mut app.platform.pipelines;
    p.mode = Mode::File;
    p.flow.folder = Some(f.0.join("loci"));
    p.flow.cache = Some(f.0.join("cache"));
    let mut screen = Screen::new();
    screen.settle(&mut app);
    (app, screen)
}

/// Draw frames until the run is over, and one more to show what it gave.
fn finish(app: &mut App, screen: &mut Screen) {
    let until = std::time::Instant::now() + Duration::from_secs(120);
    while app.platform.pipelines.flow.is_running() {
        assert!(std::time::Instant::now() < until, "the run did not end");
        std::thread::sleep(Duration::from_millis(20));
        screen.frame(app, vec![]);
    }
    screen.settle(app);
}

#[test]
fn the_file_mode_is_chosen_on_the_page_and_shows_the_pipeline_and_what_it_needs() {
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Pipelines;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert_eq!(app.platform.pipelines.mode, Mode::Route);
    screen.click(&mut app, "Pipeline file");
    assert_eq!(app.platform.pipelines.mode, Mode::File);
    let text = screen.text();
    for want in [
        "Gene trees, then a species tree",
        "1. tree engine=builtin bootstrap=0",
        "2. species",
        "Open a pipeline file…",
        "Choose the folder…",
        "No folder chosen.",
        "None: every step runs every time.",
        "No run yet",
    ] {
        assert!(text.contains(want), "no \"{want}\":\n{text}");
    }
    screen.click(&mut app, "Fixed route");
    assert_eq!(app.platform.pipelines.mode, Mode::Route);
}

#[test]
fn a_pipeline_to_start_from_replaces_the_one_in_use() {
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Pipelines;
    app.platform.pipelines.mode = Mode::File;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, "Align, gene trees, then a species tree");
    let p = &app.platform.pipelines.flow.pipeline;
    assert_eq!(*p, Pipeline::unaligned_loci_to_species());
    assert!(screen.text().contains("1. align engine=builtin"));
}

#[test]
fn a_file_that_cannot_be_read_leaves_the_pipeline_as_it_was_and_says_why() {
    let mut page = FlowPage::default();
    let before = page.pipeline.clone();
    page.load_text(
        "helixloom-pipeline 1\ninput: alignment\nstep align\n",
        Path::new("x.txt"),
    );
    assert_eq!(page.pipeline, before);
    assert!(
        page.message.starts_with("x.txt: ") && page.message.contains("Step 1 (align)"),
        "{}",
        page.message
    );
    page.load_text(
        &Pipeline::unaligned_loci_to_species().to_text(),
        Path::new("y.txt"),
    );
    assert_eq!(page.pipeline, Pipeline::unaligned_loci_to_species());
    assert_eq!(page.message, "Read y.txt.");
}

#[test]
fn nothing_runs_until_a_folder_is_chosen() {
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Pipelines;
    app.platform.pipelines.mode = Mode::File;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, "Run pipeline");
    assert!(!app.platform.pipelines.flow.is_running());
    assert!(app.platform.pipelines.flow.report.is_none());
}

#[test]
fn a_run_shows_what_each_step_did_and_the_second_finds_it_all_in_the_cache() {
    let f = Folder::new();
    let (mut app, mut screen) = page(&f);
    screen.click(&mut app, "Run pipeline");
    // the run is on a thread of its own, and the job button counts it; a quick run may be over by the time the click is drawn
    if app.platform.pipelines.flow.is_running() {
        assert_eq!(app.platform.busy_count(), 1);
    }
    finish(&mut app, &mut screen);
    assert!(
        screen.text().contains("4 answers kept"),
        "counted when the run ended"
    );
    let report = app
        .platform
        .pipelines
        .flow
        .report
        .clone()
        .expect("a report");
    assert_eq!((report.ran(), report.cached(), report.failed()), (4, 0, 0));
    let text = screen.text();
    for want in [
        "What each step did",
        "gene1",
        "gene3",
        "all files",
        "The species tree",
        "4 run",
        "0 from the cache",
    ] {
        assert!(text.contains(want), "no \"{want}\":\n{text}");
    }

    screen.click(&mut app, "Run pipeline");
    finish(&mut app, &mut screen);
    let again = app.platform.pipelines.flow.report.clone().unwrap();
    assert_eq!((again.ran(), again.cached()), (0, 4));
    assert_eq!(again.output, report.output);
    assert!(screen.text().contains("4 from the cache"));
    assert_eq!(app.platform.busy_count(), 0);
}

#[test]
fn a_run_of_a_file_is_something_for_the_session_to_keep_and_each_run_and_each_other_pipeline_is_another(
) {
    let f = Folder::new();
    let (mut app, mut screen) = page(&f);
    let stamp = |app: &App| app.platform.pipelines.stamp();
    assert_eq!(stamp(&app), None, "no run yet");
    screen.click(&mut app, super::super::RUN);
    finish(&mut app, &mut screen);
    let first = stamp(&app).expect("a run to keep");
    // the same run again is another run (its steps come from the cache now)
    screen.click(&mut app, super::super::RUN);
    finish(&mut app, &mut screen);
    let second = stamp(&app).expect("a run to keep");
    assert_ne!(first, second);
    // another pipeline is chosen: the report is no longer the answer to it, and nothing is left to keep
    let chosen = app.platform.pipelines.flow.pipeline.clone();
    app.platform.pipelines.flow.choose(chosen);
    assert_eq!(stamp(&app), None);
    // a result fetched from the workstation's list, and the example of the pictures, are runs to keep too
    let result = crate::remote::flow::FlowResult {
        pipeline: Pipeline::loci_to_species(),
        report: Report::default(),
        seconds: 0.0,
    };
    app.show_flow_result(result.clone(), "fetched");
    let fetched = stamp(&app).expect("a run to keep");
    // the same one fetched again is another showing, and so is the example shown twice
    app.show_flow_result(result, "fetched");
    let again = stamp(&app).expect("a run to keep");
    assert_ne!(again, fetched);
    app.show_example_flow_run();
    let example = stamp(&app).expect("a run to keep");
    app.show_example_flow_run();
    let twice = stamp(&app).expect("a run to keep");
    assert!(example != again && twice != example);
}

#[test]
fn a_locus_that_fails_is_named_and_the_others_go_on() {
    let f = Folder::new();
    let (mut app, mut screen) = page(&f);
    std::fs::write(
        f.0.join("loci/broken.fasta"),
        ">a\nACGTACGT\n>b\nACGTACGA\n",
    )
    .unwrap();
    screen.click(&mut app, "Run pipeline");
    finish(&mut app, &mut screen);
    let text = screen.text();
    assert!(text.contains("FAILED"), "{text}");
    assert!(
        text.contains("Left out of what was gathered: broken"),
        "{text}"
    );
    assert!(text.contains("1 failed"), "{text}");
    let page = &app.platform.pipelines.flow;
    let report = page.report.as_ref().unwrap();
    assert!(report.output.is_some());
    // the run has its record all the same: the file that failed is one of its inputs, and what its step said is not
    // among the programs
    let (made, _) = page.made.as_ref().expect("a record");
    assert!(made.inputs.iter().any(|i| i.name == "broken"));
    let failed: Vec<&String> = report
        .items
        .iter()
        .flat_map(|i| &i.steps)
        .filter(|s| matches!(s.status, Status::Failed(_)))
        .map(|s| &s.note)
        .collect();
    assert_eq!(failed.len(), 1);
    assert!(
        failed[0].is_empty() || !made.programs.contains(failed[0]),
        "{:?}",
        made.programs
    );
}

#[test]
fn the_program_of_a_step_is_what_its_note_names_before_the_model_and_before_why_it_was_not_kept() {
    use crate::pipeline::flow::StepReport;
    let step = |step: Step, status: Status, note: &str| StepReport {
        step,
        status,
        seconds: 0.0,
        note: note.to_string(),
    };
    let tree = Pipeline::loci_to_species().steps[0].clone();
    assert!(matches!(tree, Step::Tree { .. }));
    let align = Pipeline::unaligned_loci_to_species().steps[0].clone();
    assert!(matches!(align, Step::Align { .. }));
    for (report, program) in [
        (
            step(tree.clone(), Status::Ran, "IQ-TREE 3.0.1, model HKY+F+G4"),
            Some("IQ-TREE 3.0.1"),
        ),
        (
            step(
                tree.clone(),
                Status::Cached,
                "IQ-TREE 3.0.1, model HKY+F+G4",
            ),
            Some("IQ-TREE 3.0.1"),
        ),
        (
            step(
                align.clone(),
                Status::Ran,
                "MAFFT v7.526; not kept in the cache: the folder is full",
            ),
            Some("MAFFT v7.526"),
        ),
        (
            step(
                tree.clone(),
                Status::Failed("no tree".into()),
                "IQ-TREE 3.0.1",
            ),
            None,
        ),
        (
            step(Step::Species, Status::Ran, "4 gene trees, 6 species"),
            None,
        ),
        (step(align.clone(), Status::Ran, ""), None),
    ] {
        assert_eq!(report.program(), program, "{report:?}");
    }
}

#[test]
fn a_folder_says_how_many_files_it_holds_their_bytes_and_when_the_newest_was_changed() {
    let f = Folder::new();
    f.loci(2);
    let dir = f.0.join("loci");
    let kind = Pipeline::loci_to_species().input;
    let first = engine::folder_stamp(&dir, kind).expect("a listing");
    let on_disk: u64 = ["gene1.fasta", "gene2.fasta"]
        .iter()
        .map(|name| std::fs::metadata(dir.join(name)).unwrap().len())
        .sum();
    assert_eq!((first.files, first.bytes), (2, on_disk));
    assert_eq!(engine::folder_size(&dir, kind), Some((2, on_disk)));
    assert_eq!(engine::folder_stamp(&dir, kind), Some(first.clone()));
    let before = folder_words(Some(&dir), kind);
    assert!(before.starts_with(&format!(
        "loci: 2 files, {on_disk} bytes, newest changed 20"
    )));
    assert!(
        before.ends_with(&format!(" UTC, listing {}", &first.listing[..8])),
        "{before}"
    );
    let set_time = |name: &str, seconds: u64| {
        let file = std::fs::File::options()
            .write(true)
            .open(dir.join(name))
            .expect("the file");
        file.set_modified(std::time::UNIX_EPOCH + Duration::from_secs(seconds))
            .expect("a time");
    };
    // a file is written again a while later with as many bytes: the folder no longer says the same
    set_time("gene1.fasta", first.newest + 90);
    let later = engine::folder_stamp(&dir, kind).expect("a listing");
    assert_eq!(
        (later.files, later.bytes, later.newest),
        (2, on_disk, first.newest + 90)
    );
    assert_ne!(later.listing, first.listing);
    assert_ne!(folder_words(Some(&dir), kind), before);
    // the file that is not the newest is put back from an older copy: the newest time is the same, the listing is not
    set_time("gene2.fasta", first.newest - 500);
    let older = engine::folder_stamp(&dir, kind).expect("a listing");
    assert_eq!(
        (older.files, older.bytes, older.newest),
        (later.files, later.bytes, later.newest)
    );
    assert_ne!(older.listing, later.listing);
    // a file gets another name: as many files, as many bytes, the same times, and another listing
    std::fs::rename(dir.join("gene2.fasta"), dir.join("geneX.fasta")).expect("another name");
    set_time("geneX.fasta", first.newest - 500);
    let renamed = engine::folder_stamp(&dir, kind).expect("a listing");
    assert_eq!(
        (renamed.files, renamed.bytes, renamed.newest),
        (older.files, older.bytes, older.newest)
    );
    assert_ne!(renamed.listing, older.listing);
    // no folder, and one that is not there
    assert_eq!(folder_words(None, kind), "none chosen");
    assert_eq!(
        folder_words(Some(&f.0.join("nowhere")), kind),
        "nowhere: cannot be listed"
    );
}

#[test]
fn the_species_tree_can_be_drawn_in_trees() {
    let f = Folder::new();
    let (mut app, mut screen) = page(&f);
    screen.click(&mut app, "Run pipeline");
    finish(&mut app, &mut screen);
    screen.click(&mut app, "Draw in Trees");
    assert_eq!(app.workspace, Workspace::Trees);
    assert!(app.tree_win.tree.is_some());
}

#[test]
fn a_folder_with_nothing_to_run_on_says_so_in_the_message() {
    let f = Folder::new();
    let (mut app, mut screen) = page(&f);
    let empty = f.0.join("empty");
    std::fs::create_dir_all(&empty).unwrap();
    app.platform.pipelines.flow.folder = Some(empty);
    screen.click(&mut app, "Run pipeline");
    finish(&mut app, &mut screen);
    let message = app.platform.pipelines.flow.message.clone();
    assert!(message.contains("holds no file of alignment"), "{message}");
    assert!(app.platform.pipelines.flow.report.is_none());
}

#[test]
fn what_a_run_is_made_of_is_not_changed_while_it_goes_on() {
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Pipelines;
    app.platform.pipelines.mode = Mode::File;
    let ctx = egui::Context::default();
    app.platform.pipelines.flow.job = Some(Job::spawn(&ctx, || {
        std::thread::sleep(Duration::from_millis(400));
        Ok((Report::default(), Vec::new()))
    }));
    let mut screen = Screen::new();
    screen.settle(&mut app);
    screen.click(&mut app, "Align, gene trees, then a species tree");
    assert_eq!(
        app.platform.pipelines.flow.pipeline,
        Pipeline::loci_to_species(),
        "ignored while running"
    );
    finish(&mut app, &mut screen);
    screen.click(&mut app, "Align, gene trees, then a species tree");
    assert_eq!(
        app.platform.pipelines.flow.pipeline,
        Pipeline::unaligned_loci_to_species()
    );
}

#[test]
fn the_cache_is_counted_when_it_changes_and_not_in_every_frame() {
    let f = Folder::new();
    let (mut app, mut screen) = page(&f);
    app.platform.pipelines.flow.count_kept();
    assert_eq!(app.platform.pipelines.flow.kept, Some(0));
    screen.settle(&mut app);
    assert!(screen.text().contains("0 answers kept"));
    // two answers appear in the folder behind the page's back
    let cache = Cache::open(&f.0.join("cache")).unwrap();
    cache.put("a", &Value::Tree("(a,b);".into()), "").unwrap();
    cache.put("b", &Value::Tree("(a,b);".into()), "").unwrap();
    screen.settle(&mut app);
    assert!(
        screen.text().contains("0 answers kept"),
        "the page does not look into the folder in every frame"
    );
    app.platform.pipelines.flow.count_kept();
    screen.settle(&mut app);
    assert!(screen.text().contains("2 answers kept"));
}

#[test]
fn a_run_leaves_a_record_of_the_files_it_read_and_the_steps_and_says_when_the_folder_has_changed() {
    let f = Folder::new();
    let (mut app, mut screen) = page(&f);
    screen.click(&mut app, "Run pipeline");
    finish(&mut app, &mut screen);
    let page = &app.platform.pipelines.flow;
    let (made, pipeline) = page.made.clone().expect("a record");
    assert_eq!(
        (made.analysis.as_str(), made.ran_on.as_str()),
        ("Pipeline file", "this computer")
    );
    assert_eq!(pipeline, page.pipeline);
    // every file by its name, with the size and the checksum of what was read from it
    let on_disk = std::fs::read(f.0.join("loci/gene2.fasta")).expect("the file");
    let names: Vec<&str> = made.inputs.iter().map(|i| i.name.as_str()).collect();
    assert_eq!(names, ["gene1", "gene2", "gene3"]);
    assert_eq!(made.inputs[1], Input::of("gene2", &on_disk));
    // the pipeline step by step, and the folder as it was listed
    let keys: Vec<&str> = made.parameters.iter().map(|(k, _)| k.as_str()).collect();
    assert_eq!(keys, ["pipeline", "input", "step 1", "step 2", "folder"]);
    assert_eq!(made.parameters[0].1, page.pipeline.name);
    assert!(
        made.parameters[2].1.starts_with("tree "),
        "{:?}",
        made.parameters
    );
    assert_eq!(made.parameters[3].1, "species");
    let folder = &made.parameters[4].1;
    assert!(folder.starts_with("loci: 3 files, "), "{folder}");
    assert!(folder.contains(" UTC, listing "), "{folder}");
    // what made the answers: the tree builder of the three files, once, without its model, and this program last. What
    // the species step found (how many trees, how the quartets agree) is not a program
    let report = page.report.as_ref().expect("a report");
    let tree_note = &report.items[0].steps[0].note;
    let builder = tree_note.split(", model ").next().expect("a name");
    assert!(
        tree_note.contains(", model ") && !builder.is_empty(),
        "{tree_note}"
    );
    assert_eq!(
        made.programs,
        [
            builder.to_string(),
            format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))
        ]
    );
    assert!(report.gathered[0].note.contains("gene trees"));
    assert_eq!(report.gathered[0].program(), None);
    assert_eq!(report.items[0].steps[0].program(), Some(builder));
    let time = screen.ctx().input(|i| i.time);
    let (inputs, now) = page.as_it_is(time);
    assert!(made.changed_since(&inputs, &now).is_empty());

    // a file is added to the folder. The folder is not listed again at once (a frame must not list it every time)...
    std::fs::write(f.0.join("loci/gene4.fasta"), &on_disk).expect("a fourth file");
    let page = &app.platform.pipelines.flow;
    let (inputs, now) = page.as_it_is(time + LOOK_AGAIN * 0.5);
    assert!(made.changed_since(&inputs, &now).is_empty());
    // ...but a moment later it is, and the record says that the folder is no longer what was read
    let (inputs, now) = page.as_it_is(time + LOOK_AGAIN);
    let changes = made.changed_since(&inputs, &now);
    assert_eq!(changes.len(), 1, "{changes:?}");
    assert!(
        changes[0].starts_with("folder was loci: 3 files, "),
        "{changes:?}"
    );
    assert!(
        changes[0].contains(", is now loci: 4 files, "),
        "{changes:?}"
    );
}

#[test]
fn the_record_of_a_run_is_copied_and_running_again_brings_its_pipeline_back() {
    let f = Folder::new();
    let (mut app, mut screen) = page(&f);
    screen.click(&mut app, "Run pipeline");
    finish(&mut app, &mut screen);
    let (made, pipeline) = app.platform.pipelines.flow.made.clone().expect("a record");
    screen.click(&mut app, made_ui::TITLE);
    screen.settle(&mut app);
    assert!(
        screen.text().contains(made_ui::UNCHANGED),
        "{}",
        screen.text()
    );
    screen.click(&mut app, made_ui::COPY);
    assert_eq!(screen.copied(), Some(made.to_text().as_str()));
    assert!(
        app.platform
            .pipelines
            .flow
            .message
            .contains("Copied the record"),
        "{}",
        app.platform.pipelines.flow.message
    );
    assert!(!app.checks_win.open);
    screen.click(&mut app, made_ui::CHECKED);
    screen.settle(&mut app);
    assert!(app.checks_win.open);
    app.checks_win.open = false;
    // a file is added while the record is open: a moment later the page itself says that the folder changed
    let one = std::fs::read(f.0.join("loci/gene1.fasta")).expect("a file");
    std::fs::write(f.0.join("loci/gene9.fasta"), one).expect("another file");
    for _ in 0..(LOOK_AGAIN * 60.0) as usize + 5 {
        screen.frame(&mut app, vec![]);
    }
    let text = screen.text();
    assert!(text.contains("folder was loci: 3 files, "), "{text}");
    std::fs::remove_file(f.0.join("loci/gene9.fasta")).expect("gone again");
    // another pipeline is put in use under the result: the open record says what differs, and Run again brings back its own
    app.platform.pipelines.flow.pipeline = Pipeline::unaligned_loci_to_species();
    screen.settle(&mut app);
    let text = screen.text();
    assert!(text.contains("pipeline was "), "{text}");
    screen.click(&mut app, made_ui::AGAIN);
    assert_eq!(app.platform.pipelines.flow.pipeline, pipeline);
    finish(&mut app, &mut screen);
    let said = &app.platform.pipelines.flow.message;
    assert!(!said.contains("Copied"), "a run ran: {said}");
    let (again, _) = app
        .platform
        .pipelines
        .flow
        .made
        .clone()
        .expect("a new record");
    assert_eq!(
        (&again.parameters, &again.inputs),
        (&made.parameters, &made.inputs)
    );
    // a pipeline to start from, or a report fetched from the workstation's list, is not what that record is of
    let report = app
        .platform
        .pipelines
        .flow
        .report
        .clone()
        .expect("a report");
    app.platform
        .pipelines
        .flow
        .choose(Pipeline::loci_to_species());
    assert!(app.platform.pipelines.flow.made.is_none());
    screen.click(&mut app, "Run pipeline");
    finish(&mut app, &mut screen);
    assert!(app.platform.pipelines.flow.made.is_some());
    let result = FlowResult {
        pipeline: Pipeline::loci_to_species(),
        report,
        seconds: 1.0,
    };
    app.show_flow_result(result, "From the workstation.");
    let page = &app.platform.pipelines.flow;
    assert!(page.report.is_some() && page.made.is_none());
    assert_eq!(page.ran_on, "the workstation");
    // the example for a picture runs here, whatever ran last, and has no record either
    app.show_example_flow_run();
    let page = &app.platform.pipelines.flow;
    assert!(page.report.is_some() && page.made.is_none());
    assert_eq!(page.ran_on, "this computer");
}
