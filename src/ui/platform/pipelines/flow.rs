//! The Pipelines page's second mode: a pipeline that is a file of steps, run on a folder of files (a locus each), with the
//! answers of every step kept so that a run that stopped, or that is made again after a change, does only what is new.

use crate::pipeline::flow::{self as engine, Cache, Pipeline, Report, Status, Step};
use crate::provenance::{Asked, Input, Provenance};
use crate::remote::flow::{FlowJob, FlowResult};
use crate::remote::tasks::{Request, Response};
use crate::remote::weight::Weight;
use crate::seq::newick::NwTree;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::dialog::{self, Wanted};
use crate::ui::fonts::FontWeight;
use crate::ui::job::{Job, Poll};
use crate::ui::platform::workstation::{run_button, ComputeView, Place};
use crate::ui::platform::Workspace;
use crate::ui::provenance as made_ui;
use eframe::egui::{self, RichText, Ui};
use std::cell::RefCell;
use std::path::{Path, PathBuf};
use std::time::Duration;

const WORKING_REDRAW: Duration = Duration::from_millis(200);
/// While the record of a run is open, the folder is listed again this often (in seconds), to say whether its files changed.
const LOOK_AGAIN: f64 = 2.0;
/// Rows of the answer's Newick text that are shown.
const NEWICK_ROWS: usize = 3;

/// The route of the page that the person is on.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum Mode {
    /// The fixed route: align, tree, primers.
    #[default]
    Route,
    /// A pipeline file run on a folder.
    File,
}

/// What a run gives: its report, and the files it read, each by its name, its size and a checksum; or why there is none.
type Ran = Result<(Report, Vec<Input>), String>;

pub struct FlowPage {
    pub pipeline: Pipeline,
    pub folder: Option<PathBuf>,
    pub cache: Option<PathBuf>,
    pub message: String,
    pub report: Option<Report>,
    /// Counted up whenever `report` is set (`set_report`): the session keeps the run, and tells by this that it is another one.
    pub(super) shown: u64,
    /// How many answers the cache folder keeps, counted when the folder is chosen and when a run or a forgetting changes it:
    /// counting is a pass over the folder, which a frame must not make.
    pub kept: Option<usize>,
    /// How many files the folder holds and how many bytes, counted when it is chosen and when a run starts (a listing of the
    /// folder, which a frame must not make): it is how big the run is, for the choice of where it goes.
    pub size: Option<(usize, u64)>,
    /// Where the run that made `report` was done: "this computer" or "the workstation"; empty before one.
    pub ran_on: &'static str,
    job: Option<Job<Ran>>,
    /// What the run that is going was asked, with its pipeline to bring back for *Run again*.
    asked: Option<(Asked, Pipeline)>,
    /// The record of the run that `report` is of; none for a report fetched from the workstation's list, whose files are not
    /// known here.
    pub made: Option<(Provenance, Pipeline)>,
    /// When the folder was last listed for the record (in the seconds of the window), and what it was then, in words
    /// (`folder_words`). In a cell, because it is listed again where the record is drawn, which only reads the page.
    looked: RefCell<(f64, String)>,
}

/// How much of the checksum of a folder's listing is shown.
const LISTING_SHOWN: usize = 8;

/// A folder of loci in words, from a listing and not a reading of it: its name, how many files it holds of the kind the
/// pipeline takes, their bytes, when the newest was changed, and the start of a checksum of the listing (every file's name,
/// size and time of change). Two listings that say the same are taken as the same files.
// ponytail: a file written again within one second to the same size is not seen as changed; compare checksums of the files
// (a reading of every one) if that ever matters
pub fn folder_words(dir: Option<&Path>, kind: engine::Kind) -> String {
    let Some(dir) = dir else {
        return "none chosen".to_string();
    };
    let Some(stamp) = engine::folder_stamp(dir, kind) else {
        return format!("{}: cannot be listed", short(dir));
    };
    format!(
        "{}: {} files, {} bytes, newest changed {} {:02}:{:02}:{:02} UTC, listing {}",
        short(dir),
        stamp.files,
        stamp.bytes,
        crate::connect::utc_date(stamp.newest),
        (stamp.newest / 3600) % 24,
        (stamp.newest / 60) % 60,
        stamp.newest % 60,
        &stamp.listing[..LISTING_SHOWN.min(stamp.listing.len())]
    )
}

/// The settings of a run, as its record says them: the pipeline step by step, and the folder as it was listed.
pub fn settings_of(pipeline: &Pipeline, folder: &str) -> Vec<(String, String)> {
    let mut all = vec![
        ("pipeline".to_string(), pipeline.name.clone()),
        ("input".to_string(), pipeline.input.word().to_string()),
    ];
    for (k, step) in pipeline.steps.iter().enumerate() {
        let said: Vec<String> = std::iter::once(step.name().to_string())
            .chain(step.params().iter().map(|(k, v)| format!("{k}={v}")))
            .collect();
        all.push((format!("step {}", k + 1), said.join(" ")));
    }
    all.push(("folder".to_string(), folder.to_string()));
    all
}

/// The files of a run as its record names them: each by its name, with the size and the checksum of what was read from it
/// (its sequences or trees as text, not the bytes on the disk).
fn inputs_of(items: &[(String, engine::Value)]) -> Vec<Input> {
    items
        .iter()
        .map(|(name, value)| Input::of(name, value.as_text().as_bytes()))
        .collect()
}

/// The programs that made the answers of a report (`StepReport::program`: the aligner and the tree builder of the steps that
/// did not fail), each once, and the time the steps took together.
fn made_by(report: &Report) -> (Vec<String>, f64) {
    let steps = report
        .items
        .iter()
        .flat_map(|i| &i.steps)
        .chain(&report.gathered);
    let mut programs: Vec<String> = Vec::new();
    let mut seconds = 0.0;
    for step in steps {
        seconds += step.seconds;
        if let Some(program) = step.program().filter(|p| !programs.iter().any(|q| q == p)) {
            programs.push(program.to_string());
        }
    }
    (programs, seconds)
}

impl Default for FlowPage {
    fn default() -> Self {
        FlowPage {
            pipeline: Pipeline::loci_to_species(),
            folder: None,
            cache: None,
            message: String::new(),
            report: None,
            shown: 0,
            kept: None,
            size: None,
            ran_on: "",
            job: None,
            asked: None,
            made: None,
            looked: RefCell::new((0.0, String::new())),
        }
    }
}

impl FlowPage {
    pub fn is_running(&self) -> bool {
        self.job.is_some()
    }

    /// Take the pipeline of a file's text; say why when it cannot be read, and keep the one that is there.
    pub fn load_text(&mut self, text: &str, from: &Path) {
        match Pipeline::parse(text) {
            Ok(p) => {
                self.message = format!("Read {}.", from.display());
                self.pipeline = p;
                self.forget_run();
                self.count_files();
            }
            Err(why) => self.message = format!("{}: {why}", from.display()),
        }
    }

    /// Count the files of the folder again.
    pub fn count_files(&mut self) {
        self.size = self
            .folder
            .as_deref()
            .and_then(|dir| engine::folder_size(dir, self.pipeline.input));
    }

    /// How big the run is: its files, and the most bootstrap replicates a step of it asks for.
    pub fn weight(&self) -> Weight {
        let (files, bytes) = self.size.unwrap_or_default();
        let bootstrap = self
            .pipeline
            .steps
            .iter()
            .map(|s| match s {
                Step::Tree { bootstrap, .. } => *bootstrap,
                _ => 0,
            })
            .max()
            .unwrap_or(0);
        Weight::of_flow(
            files,
            usize::try_from(bytes).unwrap_or(usize::MAX),
            bootstrap,
        )
    }

    /// Count the answers in the cache folder again.
    pub fn count_kept(&mut self) {
        self.kept = self.cache.as_deref().map(Cache::kept_in);
    }

    /// A pipeline to start from; what was run before is no longer the answer to it.
    pub fn choose(&mut self, p: Pipeline) {
        self.pipeline = p;
        self.forget_run();
        self.message.clear();
        self.count_files();
    }

    /// The report of the last run is no longer shown, and its record goes with it.
    pub(super) fn forget_run(&mut self) {
        self.set_report(None);
        self.made = None;
    }

    /// Show `report`, or none. Every change of the report that is shown goes through here, so that it is counted.
    pub(super) fn set_report(&mut self, report: Option<Report>) {
        self.report = report;
        self.shown += 1;
    }

    /// The folder in words (`folder_words`), listed now: `now` is the time of the window.
    fn look_at_folder(&self, now: f64) -> String {
        let words = folder_words(self.folder.as_deref(), self.pipeline.input);
        *self.looked.borrow_mut() = (now, words.clone());
        words
    }

    /// The files and the settings as they are at `now` (the time of the window), to say what has changed since the run that
    /// is shown. The files are not read again for this (a frame must not): they are the ones of the record, and the folder's
    /// line says, from a listing, whether they are still what was read. The listing is made again when the last one is
    /// `LOOK_AGAIN` old, and only here: where the record is open.
    // ponytail: the listing is made on the thread that draws; give it to a job if a folder on a slow drive makes the page
    // stutter while the record is open
    pub fn as_it_is(&self, now: f64) -> (Vec<Input>, Vec<(String, String)>) {
        let inputs = self
            .made
            .as_ref()
            .map_or_else(Vec::new, |(made, _)| made.inputs.clone());
        let (looked_at, words) = self.looked.borrow().clone();
        let folder = if now - looked_at >= LOOK_AGAIN {
            self.look_at_folder(now)
        } else {
            words
        };
        (inputs, settings_of(&self.pipeline, &folder))
    }
}

#[derive(Default)]
pub struct FlowClicks {
    load: bool,
    save: bool,
    preset: Option<usize>,
    folder: bool,
    cache: bool,
    no_cache: bool,
    forget: bool,
    run: bool,
    copy: bool,
    to_trees: bool,
    /// What was pressed under the record of the run.
    made: made_ui::Asks,
}

/// The pipelines the page offers to start from.
fn presets() -> [Pipeline; 3] {
    [
        Pipeline::loci_to_species(),
        Pipeline::unaligned_loci_to_species(),
        Pipeline::inspect_graphs(),
    ]
}

fn step_line(k: usize, step: &Step) -> String {
    let settings: Vec<String> = step
        .params()
        .iter()
        .map(|(k, v)| format!("{k}={v}"))
        .collect();
    format!(
        "{}. {} {}  ({} → {})",
        k + 1,
        step.name(),
        settings.join(" "),
        step.takes().word(),
        step.gives().word()
    )
}

fn short(path: &Path) -> String {
    path.file_name().map_or_else(
        || path.display().to_string(),
        |n| n.to_string_lossy().into_owned(),
    )
}

pub(super) fn options(
    ui: &mut Ui,
    look: &Look,
    page: &mut FlowPage,
    view: ComputeView,
    c: &mut FlowClicks,
) {
    design::section(ui, look, "Pipeline");
    ui.label(
        RichText::new(&page.pipeline.name)
            .semibold()
            .color(look.text),
    );
    for (k, step) in page.pipeline.steps.iter().enumerate() {
        ui.label(
            RichText::new(step_line(k, step))
                .monospace()
                .size(12.0)
                .color(look.dim),
        );
    }
    if let Some(text) = page.pipeline.split_text() {
        ui.label(RichText::new(text).size(12.0).color(look.dim));
    }
    ui.horizontal_wrapped(|ui| {
        c.load = ui.button("Open a pipeline file…").clicked();
        c.save = ui.button("Save…").clicked();
    });
    ui.label(RichText::new("Or start from:").size(12.0).color(look.dim));
    for (k, p) in presets().iter().enumerate() {
        if ui.small_button(&p.name).clicked() {
            c.preset = Some(k);
        }
    }

    design::section(ui, look, "Input");
    match &page.folder {
        Some(f) => {
            ui.label(RichText::new(short(f)).semibold().color(look.text))
                .on_hover_text(f.display().to_string());
        }
        None => {
            ui.label(RichText::new("No folder chosen.").color(look.dim));
        }
    }
    ui.label(
        RichText::new(format!(
            "Input files: {} ({}).",
            page.pipeline.input.word(),
            engine::extensions(page.pipeline.input).join(", ")
        ))
        .size(12.0)
        .color(look.dim),
    );
    c.folder = ui.button("Choose the folder…").clicked();

    design::section(ui, look, "Cache");
    match &page.cache {
        Some(f) => {
            ui.label(RichText::new(short(f)).semibold().color(look.text))
                .on_hover_text(f.display().to_string());
            let kept = page
                .kept
                .map_or(String::new(), |n| format!("{n} answers kept: "));
            ui.label(
                RichText::new(format!(
                    "{kept}a step that has its answer here is not run again."
                ))
                .size(12.0)
                .color(look.dim),
            );
        }
        None => {
            ui.label(RichText::new("None: every step runs every time.").color(look.dim));
        }
    }
    ui.horizontal_wrapped(|ui| {
        c.cache = ui.button("Choose the cache folder…").clicked();
        if page.cache.is_some() {
            c.no_cache = ui.button("Keep nothing").clicked();
            c.forget = ui.button("Forget the answers").clicked();
        }
    });

    let state = (page.folder.is_some(), page.is_running());
    c.run = run_button(ui, look, view, super::RUN, state);
    if view.there && page.cache.is_some() {
        ui.label(
            RichText::new(
                "On the workstation the answers are kept there, in a folder of its own: the cache folder above is not used.",
            )
            .size(12.0)
            .color(look.dim),
        );
    }
    if !page.message.is_empty() {
        ui.add_space(4.0);
        ui.label(RichText::new(&page.message).size(12.0).color(look.dim));
    }
}

fn status_text(look: &Look, s: &Status, seconds: f64) -> RichText {
    match s {
        Status::Ran => RichText::new(format!("ran {seconds:.2} s")).color(look.text),
        Status::Cached => RichText::new("from the cache").color(look.good),
        Status::Failed(_) => RichText::new("FAILED").color(look.bad).semibold(),
    }
}

fn graph_summaries(ui: &mut Ui, look: &Look, report: &Report) {
    let mut graphs = report
        .items
        .iter()
        .filter_map(|item| match &item.output {
            Some(engine::Value::GraphSummary(summary)) => Some((&item.name, summary)),
            _ => None,
        })
        .peekable();
    if graphs.peek().is_none() {
        return;
    }
    design::titled_card(ui, look, "Graph summaries", |ui| {
        egui::ScrollArea::horizontal().show(ui, |ui| {
            egui::Grid::new("flow_graph_summaries")
                .striped(true)
                .show(ui, |ui| {
                    for label in [
                        "File",
                        "Segments",
                        "Links",
                        "Components",
                        "Isolated",
                        "Known bases",
                        "Unknown lengths",
                    ] {
                        ui.strong(label);
                    }
                    ui.end_row();
                    for (name, graph) in graphs {
                        ui.label(name);
                        ui.label(graph.segments.to_string());
                        ui.label(graph.links.to_string());
                        ui.label(graph.components.to_string());
                        ui.label(graph.isolated.to_string());
                        ui.label(graph.known_bases.to_string());
                        ui.label(graph.unknown_lengths.to_string());
                        if ui
                            .small_button("Copy TSV")
                            .on_hover_text(format!("Copy the summary for {name}"))
                            .clicked()
                        {
                            ui.ctx().copy_text(graph.to_tsv());
                        }
                        ui.end_row();
                    }
                });
        });
    });
    ui.add_space(design::CARD_GAP);
}

pub(super) fn results(ui: &mut Ui, look: &Look, page: &FlowPage, c: &mut FlowClicks) {
    let Some(report) = &page.report else {
        let message = if page.is_running() {
            "Working… the page stays free while the pipeline runs."
        } else if page.pipeline.input == engine::Kind::Gfa {
            "Choose a folder of .gfa files and press Run pipeline to see a summary for each graph."
        } else {
            "Choose the folder of loci and press Run pipeline. The gene trees and the species tree appear here, \
             with what each step did: run now, or found in the cache from an earlier run."
        };
        if design::empty_state(ui, look, "No run yet", message, &["Choose the folder…"]) == Some(0)
        {
            c.folder = true;
        }
        return;
    };
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            design::card(ui, look, |ui| {
                if !page.ran_on.is_empty() {
                    ui.label(
                        RichText::new(format!("Ran on {}.", page.ran_on))
                            .size(12.0)
                            .color(look.dim),
                    );
                }
                if let Some((made, _)) = &page.made {
                    let can_run = page.folder.is_some() && !page.is_running();
                    let time = ui.input(|i| i.time);
                    c.made = made_ui::section(ui, look, made, || page.as_it_is(time), can_run);
                }
                ui.horizontal_wrapped(|ui| {
                    design::chip(ui, look, &format!("{} run", report.ran()), Tone::Accent);
                    design::chip(
                        ui,
                        look,
                        &format!("{} from the cache", report.cached()),
                        Tone::Good,
                    );
                    if report.failed() > 0 {
                        design::chip(ui, look, &format!("{} failed", report.failed()), Tone::Bad);
                    }
                });
                if !report.left_out.is_empty() {
                    let label = if report.gathered.is_empty() {
                        "Failed inputs"
                    } else {
                        "Left out of what was gathered"
                    };
                    ui.label(
                        RichText::new(format!("{label}: {}", report.left_out.join(", ")))
                            .color(look.bad),
                    );
                }
                for warning in &report.warnings {
                    ui.label(RichText::new(warning).size(12.0).color(look.bad));
                }
            });
            ui.add_space(design::CARD_GAP);
            graph_summaries(ui, look, report);
            design::titled_card(ui, look, "What each step did", |ui| {
                egui::Grid::new("flow_steps")
                    .striped(true)
                    .spacing([14.0, 3.0])
                    .show(ui, |ui| {
                        let rows = report
                            .items
                            .iter()
                            .map(|i| (i.name.as_str(), &i.steps))
                            .chain(std::iter::once(("all files", &report.gathered)));
                        for (name, steps) in rows {
                            for r in steps {
                                ui.label(name);
                                ui.label(r.step.name());
                                let text = status_text(look, &r.status, r.seconds);
                                let shown = ui.label(text);
                                match &r.status {
                                    Status::Failed(why) => {
                                        shown.on_hover_text(why);
                                        ui.label(RichText::new(why).size(12.0).color(look.bad));
                                    }
                                    _ => {
                                        ui.label(RichText::new(&r.note).size(12.0).color(look.dim));
                                    }
                                }
                                ui.end_row();
                            }
                        }
                    });
            });
            if let Some(engine::Value::Tree(newick)) = &report.output {
                ui.add_space(design::CARD_GAP);
                design::titled_card(ui, look, "The species tree", |ui| {
                    let mut text = newick.clone();
                    ui.add(
                        egui::TextEdit::multiline(&mut text)
                            .font(egui::TextStyle::Monospace)
                            .desired_rows(NEWICK_ROWS)
                            .desired_width(f32::INFINITY),
                    );
                    ui.horizontal_wrapped(|ui| {
                        c.copy = ui.button("Copy Newick").clicked();
                        c.to_trees = design::primary_button(ui, look, "Draw in Trees").clicked();
                    });
                });
            }
        });
}

impl App {
    pub(super) fn apply_flow_clicks(
        &mut self,
        c: FlowClicks,
        ctx: &egui::Context,
        out: &mut Vec<Action>,
    ) {
        self.checks_win.open |= c.made.checks;
        let page = &mut self.platform.pipelines.flow;
        // what a run is made of is not changed under it: its answer would then be shown beside what it was not made from
        let busy = page.is_running();
        if let (true, Some((made, _))) = (c.made.copy, &page.made) {
            ctx.copy_text(made.to_text());
            page.message = "Copied the record of the run.".into();
        }
        let again = c.made.again && !busy;
        if let (true, Some(pipeline)) = (again, page.made.as_ref().map(|(_, p)| p.clone())) {
            // with the pipeline the run that is shown was made with, on the folder as it is now
            page.pipeline = pipeline;
        }
        if c.load && !busy {
            let ask = dialog::FileDialog::new().set_title("Open a pipeline file");
            if let Some(path) = dialog::pick_one(Wanted::Anything, ask) {
                match std::fs::read_to_string(&path) {
                    Ok(text) => page.load_text(&text, &path),
                    Err(e) => page.message = format!("Could not read {}: {e}", path.display()),
                }
            }
        }
        if c.save {
            let ask = dialog::FileDialog::new().set_title("Save the pipeline file");
            if let Some(saved) =
                dialog::save_bytes(ask, "pipeline.txt", page.pipeline.to_text().as_bytes())
            {
                page.message = saved.message;
            }
        }
        if let Some(k) = c.preset.filter(|_| !busy) {
            if let Some(p) = presets().into_iter().nth(k) {
                page.choose(p);
            }
        }
        if c.folder && !busy {
            if let Some(f) = choose_folder("Choose the input folder") {
                page.folder = Some(f);
                page.forget_run();
                page.count_files();
            }
        }
        if c.cache && !busy {
            if let Some(f) = choose_folder("Choose the cache folder") {
                page.cache = Some(f);
                page.count_kept();
            }
        }
        if c.no_cache && !busy {
            page.cache = None;
            page.count_kept();
        }
        if c.forget && !busy {
            if let Some(dir) = &page.cache {
                page.message = match Cache::open(dir).and_then(|cache| cache.clear()) {
                    Ok(n) => format!("Forgot {n} answers."),
                    Err(why) => why,
                };
                page.count_kept();
            }
        }
        if c.run || again {
            self.start_flow(ctx);
        }
        let page = &mut self.platform.pipelines.flow;
        if c.copy {
            if let Some(engine::Value::Tree(newick)) =
                page.report.as_ref().and_then(|r| r.output.clone())
            {
                ctx.copy_text(newick);
                page.message = "Copied the species tree as Newick.".into();
            }
        }
        if c.to_trees {
            let node = match page.report.as_ref().and_then(|r| r.output.as_ref()) {
                Some(engine::Value::Tree(newick)) => NwTree::parse_all(newick)
                    .ok()
                    .and_then(|trees| trees.first().map(NwTree::to_node)),
                _ => None,
            };
            if let Some(node) = node {
                self.show_tree(node, "The species tree of the pipeline run.");
                out.push(Action::Go(Workspace::Trees));
            }
        }
    }

    pub(in crate::ui) fn start_flow(&mut self, ctx: &egui::Context) {
        let page = &mut self.platform.pipelines.flow;
        let Some(folder) = page.folder.clone() else {
            page.message = "Choose an input folder first.".into();
            return;
        };
        // the folder may have changed since it was chosen, and how big it is says where the run goes
        page.count_files();
        let weight = page.weight();
        let place = self.platform.place(&weight);
        let page = &mut self.platform.pipelines.flow;
        let (pipeline, cache) = (page.pipeline.clone(), page.cache.clone());
        page.forget_run();
        page.ran_on = match place {
            Place::There(_) => "the workstation",
            Place::Here => "this computer",
        };
        // what the run is asked: the files are named when they have been read, on the thread of the run
        let folder_now = page.look_at_folder(ctx.input(|i| i.time));
        let asked = Asked {
            analysis: "Pipeline file",
            inputs: Vec::new(),
            parameters: settings_of(&pipeline, &folder_now),
            ran_on: page.ran_on,
        };
        page.asked = Some((asked, pipeline.clone()));
        let Place::There(remote) = place else {
            page.message = "Working…".into();
            page.job = Some(Job::spawn(ctx, move || {
                let items = engine::read_folder(&folder, pipeline.input)?;
                let read = inputs_of(&items);
                let cache = cache.as_deref().map(Cache::open).transpose()?;
                let report = engine::run(
                    &pipeline,
                    items,
                    cache.as_ref(),
                    &engine::available_for(&pipeline),
                )?;
                Ok((report, read))
            }));
            return;
        };
        // the files are read here and sent; the workstation keeps the answers of the steps in a folder of its own
        page.message = "Working on the workstation…".into();
        page.job = Some(Job::spawn_async(ctx, async move {
            let loci = engine::read_folder(&folder, pipeline.input)?;
            let read = inputs_of(&loci);
            let label =
                crate::remote::one_line(&format!("{} on {} files", pipeline.name, loci.len()));
            let request = Request::Flow(FlowJob { pipeline, loci });
            remote
                .run_async(&request, &label)
                .await
                .and_then(Response::into_flow)
                .map(|result| (result.report, read))
        }));
    }

    /// Show what the workstation made of a pipeline file, in the page that runs them. A run that is going on this page is not
    /// replaced by it: its answer would then come over the one that was asked for.
    pub(in crate::ui) fn show_flow_result(&mut self, result: FlowResult, note: &str) {
        let pipelines = &mut self.platform.pipelines;
        if pipelines.flow.is_running() {
            pipelines.flow.message =
                "A run is going on this page: fetch the workstation's result again when it is over."
                    .into();
        } else {
            pipelines.mode = Mode::File;
            let page = &mut pipelines.flow;
            page.pipeline = result.pipeline;
            page.ran_on = "the workstation";
            page.message = format!("{note} {}", result.report.summary());
            // fetched from the workstation's list: the files it was made from are not known here
            page.set_report(Some(result.report));
            page.made = None;
            page.count_files();
        }
        self.workspace = Workspace::Pipelines;
    }

    /// For a picture: the example loci (windows of the demo isolates) run at once, as if a folder of them had been chosen.
    pub(in crate::ui) fn show_example_flow_run(&mut self) {
        let isolates = crate::sample::tree_demo();
        let items = (0..4)
            .map(|k| {
                let rows = isolates
                    .iter()
                    .map(|d| (d.name.clone(), d.seq[k * 120..(k + 1) * 120 + 60].to_vec()))
                    .collect();
                (format!("locus{}", k + 1), engine::Value::Alignment(rows))
            })
            .collect();
        let page = &mut self.platform.pipelines.flow;
        // the example runs here, wherever the last run was
        page.ran_on = "this computer";
        match engine::run(&page.pipeline, items, None, &Default::default()) {
            Ok(report) => {
                page.message = report.summary();
                page.set_report(Some(report));
                page.made = None;
            }
            Err(why) => page.message = why,
        }
        page.folder = Some(PathBuf::from("example-loci"));
    }

    /// Collect a finished run; once a frame.
    pub(super) fn flow_tick(&mut self, ctx: &egui::Context) {
        let page = &mut self.platform.pipelines.flow;
        let Some(polled) = page.job.as_ref().map(Job::poll) else {
            return;
        };
        match polled {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(done) => {
                page.job = None;
                let asked = page.asked.take();
                match done {
                    Ok((report, read)) => {
                        page.message = report.summary();
                        page.made = asked.map(|(asked, pipeline)| {
                            let (programs, seconds) = made_by(&report);
                            let asked = Asked {
                                inputs: read,
                                ..asked
                            };
                            (asked.made(programs, seconds), pipeline)
                        });
                        page.set_report(Some(report));
                    }
                    Err(why) => page.message = why,
                }
                page.count_kept();
            }
            Poll::Died => {
                page.job = None;
                page.message = "The run stopped without an answer.".into();
            }
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
fn choose_folder(title: &str) -> Option<PathBuf> {
    dialog::FileDialog::new().set_title(title).pick_folder()
}

/// A browser has no folders to choose.
#[cfg(target_arch = "wasm32")]
fn choose_folder(_title: &str) -> Option<PathBuf> {
    None
}

#[cfg(test)]
mod tests;
