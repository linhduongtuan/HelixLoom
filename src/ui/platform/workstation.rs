//! Workstation: connect this computer to the one that does the heavy work. Jobs sent there (alignments, trees, pipelines,
//! PhyloCNN estimates, structure predictions) run on its processors and GPUs and wait there until they are fetched; this
//! page says how the workstation is, lists its jobs and opens a result in the page it belongs to.

mod browser;
mod compute;

use super::Workspace;
use crate::pipeline::tools;
use crate::remote::client::{Remote, ServerStatus};
use crate::remote::jobs::{JobInfo, Status};
use crate::remote::settings::{self, Compute, Saved};
use crate::remote::tasks::Response;
use crate::remote::weight::Weight;
use crate::remote::PROTOCOL;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::job::{Job, Poll};
use eframe::egui::{self, Align, Layout, RichText, ScrollArea, Ui};
use std::path::PathBuf;
use std::time::Duration;
use web_time::Instant;

pub use compute::Place;
use compute::{choice_id, explain};
pub(in crate::ui) use compute::{compute_picker, ComputeView};

const OPTIONS_WIDTH: f32 = 340.0;
const CARD_GAP: f32 = 10.0;
/// How often the page asks the workstation how it is while the page is showing.
const REFRESH: Duration = Duration::from_secs(3);
const WORKING_REDRAW: Duration = Duration::from_millis(200);
/// The newest jobs the list shows; older ones stay on the workstation until they expire.
const LISTED: usize = 50;
const LABEL_CHARS: usize = 54;

/// What a piece of work with the workstation brought back.
enum Done {
    Looked(ServerStatus, Vec<JobInfo>),
    /// A job was removed, or (the flag) asked to stop because it was running.
    Removed(String, bool),
    Fetched(JobInfo, Response),
}

/// What a click on a job asked for, waiting for the page to be free: one request to the workstation is on its way at a time,
/// and a second would otherwise push the first out.
enum Pending {
    Fetch(JobInfo),
    Remove(String, bool),
}

#[derive(Default)]
pub struct WorkstationPage {
    pub address: String,
    pub token: String,
    show_token: bool,
    /// The connection the person made last; jobs go there as `compute` says.
    active: Option<Remote>,
    /// Where the heavy jobs are computed. It is `Local` while no workstation is connected.
    pub compute: Compute,
    status: Option<ServerStatus>,
    jobs: Vec<JobInfo>,
    /// Why the last request failed; empty when it worked.
    problem: String,
    message: String,
    work: Option<Job<Result<Done, String>>>,
    pending: Option<Pending>,
    looked: Option<Instant>,
    /// Where the connection is kept between runs; none keeps nothing.
    dir: Option<PathBuf>,
}

impl WorkstationPage {
    /// Whether the workstation that was connected has answered, and not failed since: its dot on the top bar is green.
    pub fn answering(&self) -> bool {
        self.status.is_some() && self.problem.is_empty()
    }

    /// The page as it was left last time: connected again, and with the choice of where jobs run, if there was one. With
    /// no folder nothing is read or kept.
    pub fn open(dir: Option<PathBuf>) -> WorkstationPage {
        let saved = dir.as_deref().map(settings::load_from).unwrap_or_default();
        let active = Remote::new(&saved.address, &saved.token).ok();
        WorkstationPage {
            address: saved.address,
            token: saved.token,
            compute: if active.is_some() {
                saved.compute
            } else {
                Compute::Local
            },
            active,
            dir,
            ..WorkstationPage::default()
        }
    }

    /// Whether a job of this weight goes to the workstation: when one is connected and the choice says so (Automatic: when
    /// the job is big or can use a GPU there).
    fn goes_there(&self, weight: &Weight) -> bool {
        self.active.is_some()
            && match self.compute {
                Compute::Local => false,
                Compute::Workstation => true,
                Compute::Automatic => weight.is_heavy(),
            }
    }

    /// Where a job of this weight runs.
    pub fn place(&self, weight: &Weight) -> Place {
        match &self.active {
            Some(remote) if self.goes_there(weight) => Place::There(remote.clone()),
            _ => Place::Here,
        }
    }

    /// What a page shows about where a job of this weight would run.
    pub(in crate::ui) fn view(&self, weight: &Weight) -> ComputeView {
        ComputeView {
            mode: self.compute,
            connected: self.active.is_some(),
            there: self.goes_there(weight),
        }
    }

    /// The programs the workstation has, as (name, version): known once it has said, empty when it could not.
    pub fn tools(&self) -> Option<&[(String, String)]> {
        match &self.status {
            Some(s) => Some(&s.tools),
            None if !self.problem.is_empty() => Some(&[]),
            None => None,
        }
    }

    /// The name of its first GPU, when it has one.
    pub fn gpu_name(&self) -> Option<&str> {
        self.status
            .as_ref()
            .and_then(|s| s.gpus.first())
            .map(|g| g.name.as_str())
    }

    /// Whether it has a program that folds proteins on its GPUs.
    pub fn folds_on_gpu(&self) -> bool {
        self.status.as_ref().is_some_and(|s| s.fold)
    }

    fn save(&mut self) {
        let Some(dir) = &self.dir else {
            return;
        };
        let saved = self.active.as_ref().map_or_else(Saved::default, |r| Saved {
            address: r.address().to_string(),
            token: r.token().to_string(),
            compute: self.compute,
        });
        if let Err(why) = settings::save_to(dir, &saved) {
            self.message = why;
        }
    }

    /// Use the address and token as typed, and look at what is there.
    pub(super) fn connect(&mut self) {
        match Remote::new(&self.address, &self.token) {
            Ok(remote) => {
                self.active = Some(remote);
                (self.status, self.jobs, self.looked, self.work) = (None, Vec::new(), None, None);
                self.problem.clear();
                self.message = "Connecting…".into();
                self.save();
            }
            Err(why) => self.problem = why,
        }
    }

    /// Forget the workstation: the connection, and the token with it.
    fn forget(&mut self) {
        *self = WorkstationPage {
            dir: self.dir.take(),
            ..WorkstationPage::default()
        };
        // saving may say that it could not (and then that is what the person should read)
        self.message = "Forgot the workstation.".into();
        self.save();
    }

    /// Choose where the heavy jobs are computed. Without a connected workstation only this computer can be chosen.
    pub(super) fn set_compute(&mut self, mode: Compute) {
        self.compute = if self.active.is_some() {
            mode
        } else {
            Compute::Local
        };
        if self.compute != Compute::Local && self.status.is_none() {
            self.looked = None;
        }
        self.save();
    }

    /// Ask the workstation how it is and what its jobs are doing.
    fn look(&mut self, ctx: &egui::Context) {
        let Some(remote) = self.active.clone() else {
            return;
        };
        self.looked = Some(Instant::now());
        self.work = Some(Job::spawn(ctx, move || {
            Ok(Done::Looked(remote.status()?, remote.jobs()?))
        }));
    }

    fn fetch(&mut self, ctx: &egui::Context, job: JobInfo) {
        let Some(remote) = self.active.clone() else {
            return;
        };
        self.message = format!("Fetching \"{}\"…", job.label);
        self.work = Some(Job::spawn(ctx, move || {
            let response = remote.result(job.kind, &job.id)?;
            Ok(Done::Fetched(job, response))
        }));
    }

    fn remove(&mut self, ctx: &egui::Context, id: String, stopping: bool) {
        let Some(remote) = self.active.clone() else {
            return;
        };
        self.work = Some(Job::spawn(ctx, move || {
            remote.delete(&id)?;
            Ok(Done::Removed(id, stopping))
        }));
    }

    /// Start what a click asked for, if nothing else is on its way.
    fn start_pending(&mut self, ctx: &egui::Context) {
        if self.work.is_some() {
            return;
        }
        match self.pending.take() {
            Some(Pending::Fetch(job)) => self.fetch(ctx, job),
            Some(Pending::Remove(id, stopping)) => self.remove(ctx, id, stopping),
            None => {}
        }
    }

    /// Take in what a piece of work brought back; a fetched result is for the app to open.
    fn finish(&mut self, done: Result<Done, String>) -> Option<(JobInfo, Response)> {
        match done {
            Ok(Done::Looked(status, jobs)) => {
                if self.status.is_none() {
                    self.message = format!("Connected to {}.", status.name);
                }
                (self.status, self.jobs) = (Some(status), jobs);
                self.problem.clear();
                None
            }
            Ok(Done::Removed(id, stopping)) => {
                // a job that was running is only asked to stop: it stays in the list, as stopped, at the next look
                if stopping {
                    self.message = "Asked the workstation to stop the job.".into();
                } else {
                    self.jobs.retain(|j| j.id != id);
                    self.message = "Removed the job.".into();
                }
                self.looked = None;
                None
            }
            Ok(Done::Fetched(job, response)) => {
                self.problem.clear();
                Some((job, response))
            }
            Err(why) => {
                self.problem = why;
                None
            }
        }
    }
}

#[derive(Default)]
struct Clicks {
    connect: bool,
    forget: bool,
    refresh: bool,
    /// Where the heavy jobs are to run from now on.
    choose: Option<Compute>,
    open: Option<JobInfo>,
    /// The job to remove, and whether it is running (so that it is stopped, not removed).
    remove: Option<(String, bool)>,
}

impl App {
    /// Keep the page's work going, whichever page is showing: the first look at a workstation that was used last time,
    /// a look every few seconds while the page shows, and the result of a fetch.
    pub(super) fn workstation_tick(&mut self, ctx: &egui::Context) {
        let showing = self.workspace == Workspace::Workstation;
        let page = &mut self.platform.workstation;
        let mut fetched = None;
        if let Some(polled) = page.work.as_ref().map(Job::poll) {
            match polled {
                Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
                Poll::Done(done) => {
                    page.work = None;
                    fetched = page.finish(done);
                }
                Poll::Died => {
                    page.work = None;
                    page.problem =
                        "The request to the workstation stopped without an answer.".into();
                }
            }
        }
        // a choice made on another page waits here: that page is drawn while this one is borrowed
        if let Some(mode) = ctx.data_mut(|d| d.remove_temp::<Compute>(choice_id())) {
            page.set_compute(mode);
        }
        page.start_pending(ctx);
        let stale = page.looked.is_none_or(|t| t.elapsed() >= REFRESH);
        // jobs may go there but it is not heard from yet (the network was not up when the program started, say): try
        // again, whichever page shows, so that the pages offer what the workstation has as soon as it can be reached
        let unheard = page.compute != Compute::Local && page.status.is_none();
        let wanted = stale && (showing || unheard);
        if page.active.is_some() && page.work.is_none() && wanted {
            page.look(ctx);
        }
        if showing && page.active.is_some() {
            ctx.request_repaint_after(REFRESH);
        }
        if let Some((job, response)) = fetched {
            self.open_remote_result(&job, response);
        }
    }

    /// Show what a job on the workstation made in the page it belongs to, and go there.
    fn open_remote_result(&mut self, job: &JobInfo, response: Response) {
        let note = format!("From the workstation's job \"{}\".", job.label);
        match response {
            Response::Align(aligned) => {
                self.platform.alignment.show_aligned(aligned, &note);
                self.workspace = Workspace::Alignment;
            }
            Response::Tree(estimate) => {
                self.tree_win.details = estimate.details;
                self.show_tree(estimate.tree, &note);
                self.workspace = Workspace::Trees;
            }
            Response::Pipeline(outcome) => {
                let page = &mut self.platform.pipelines;
                (page.outcome, page.message) = (Some(outcome), note);
                self.workspace = Workspace::Pipelines;
            }
            Response::PhyloCnn(result) => self.show_phylocnn_result(result),
            Response::Species(result) => self.show_species_result(result),
            Response::Fold(folded) => {
                let page = &mut self.platform.structure;
                match page.show_pdb("ESMFold prediction", &folded.pdb) {
                    Ok(()) => page.message = format!("{note} {}", page.message),
                    Err(why) => {
                        page.message = format!("The workstation's structure cannot be shown: {why}")
                    }
                }
                self.workspace = Workspace::Structure;
            }
        }
    }

    pub(super) fn workstation_page(&mut self, ui: &mut Ui, _out: &mut Vec<Action>) {
        if cfg!(target_arch = "wasm32") {
            return browser::not_available(ui);
        }
        let look = design::look(self.prefs.dark);
        let mut c = Clicks::default();
        let page = &mut self.platform.workstation;
        design::page(ui, |ui| {
            design::page_header(
                ui,
                &look,
                "Workstation",
                "Send the heavy work to a computer with GPUs and read the results here",
                |ui| {
                    if page.active.is_some() {
                        c.refresh = ui.button("Refresh").clicked();
                    }
                },
            );
            egui::Panel::left("workstation_options")
                .exact_size(OPTIONS_WIDTH)
                .resizable(false)
                .frame(egui::Frame::new())
                .show(ui, |ui| {
                    ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            design::card(ui, &look, |ui| connection_panel(ui, &look, page, &mut c));
                        });
                });
            egui::CentralPanel::default()
                .frame(egui::Frame::new().inner_margin(egui::Margin {
                    left: design::PAGE_MARGIN,
                    ..egui::Margin::ZERO
                }))
                .show(ui, |ui| overview(ui, &look, page, &mut c));
        });
        let ctx = ui.ctx().clone();
        self.apply_workstation_clicks(c, &ctx);
    }

    fn apply_workstation_clicks(&mut self, c: Clicks, ctx: &egui::Context) {
        let page = &mut self.platform.workstation;
        if c.connect {
            page.connect();
        }
        if c.forget {
            page.forget();
        }
        if let Some(mode) = c.choose {
            page.set_compute(mode);
        }
        if c.refresh {
            page.looked = None;
        }
        if let Some(job) = c.open {
            page.pending = Some(Pending::Fetch(job));
        }
        if let Some((id, stopping)) = c.remove {
            page.pending = Some(Pending::Remove(id, stopping));
        }
        page.start_pending(ctx);
    }
}

fn dim(look: &Look, text: impl Into<String>) -> RichText {
    RichText::new(text).size(12.0).color(look.dim)
}

fn connection_panel(ui: &mut Ui, look: &Look, page: &mut WorkstationPage, c: &mut Clicks) {
    design::section(ui, look, "The workstation");
    ui.label("Address");
    ui.add(
        egui::TextEdit::singleline(&mut page.address)
            .hint_text("workstation.local:7878")
            .desired_width(f32::INFINITY),
    );
    ui.label("Token");
    ui.add(
        egui::TextEdit::singleline(&mut page.token)
            .password(!page.show_token)
            .hint_text("made with helixloom token")
            .desired_width(f32::INFINITY),
    );
    ui.checkbox(&mut page.show_token, "Show the token");
    ui.horizontal(|ui| {
        let typed = !page.address.trim().is_empty() && !page.token.trim().is_empty();
        c.connect = ui
            .add_enabled_ui(typed, |ui| design::primary_button(ui, look, "Connect"))
            .inner
            .clicked();
        if page.active.is_some() {
            c.forget = ui
                .button("Forget")
                .on_hover_text("Remove the address and the token from this computer")
                .clicked();
        }
    });
    if !page.problem.is_empty() {
        ui.label(RichText::new(&page.problem).size(12.0).color(look.bad));
    } else if !page.message.is_empty() {
        ui.label(dim(look, &page.message));
    }

    heavy_jobs_section(ui, look, page, c);

    design::section(ui, look, "Keep it private");
    ui.label(
        RichText::new(
            "The connection is not encrypted: the token and your sequences travel as they are. Use it on a private \
             network or through an SSH tunnel (DEPLOY.md shows how).",
        )
        .size(11.5)
        .color(look.warn),
    );
}

/// Where the heavy jobs are computed: this computer, the workstation, or each job by its size. The same choice is on every
/// page that starts a job.
fn heavy_jobs_section(ui: &mut Ui, look: &Look, page: &WorkstationPage, c: &mut Clicks) {
    design::section(ui, look, "Where jobs run");
    let connected = page.active.is_some();
    let mut mode = page.compute;
    for choice in Compute::ALL {
        // this computer can always be chosen; the others need a workstation to send to
        let usable = connected || choice == Compute::Local;
        ui.add_enabled_ui(usable, |ui| {
            ui.radio_value(&mut mode, choice, choice.label());
        });
        ui.label(dim(look, explain(choice)));
    }
    if mode != page.compute {
        c.choose = Some(mode);
    }
    if !connected {
        ui.label(dim(look, "Connect a workstation to send jobs there."));
    }
}

fn overview(ui: &mut Ui, look: &Look, page: &WorkstationPage, c: &mut Clicks) {
    ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            if page.active.is_none() {
                return setup_card(ui, look);
            }
            status_card(ui, look, page);
            ui.add_space(CARD_GAP);
            jobs_card(ui, look, page, c);
        });
}

/// A card with a title.
fn titled(ui: &mut Ui, look: &Look, title: &str, add: impl FnOnce(&mut Ui)) {
    design::card(ui, look, |ui| {
        ui.set_min_width(ui.available_width());
        ui.label(RichText::new(title).size(15.0).strong().color(look.text));
        ui.add_space(4.0);
        add(ui);
    });
}

fn setup_card(ui: &mut Ui, look: &Look) {
    titled(ui, look, "Set up a workstation", |ui| {
        ui.label(RichText::new("Run the heavy work on a computer with GPUs and read the results here, on this one.").color(look.dim));
        ui.add_space(6.0);
        let step = |ui: &mut Ui, text: &str, command: &str| {
            ui.label(RichText::new(text).color(look.text));
            ui.label(
                RichText::new(command)
                    .monospace()
                    .size(12.0)
                    .color(look.accent),
            );
        };
        step(
            ui,
            "1. On the workstation, make a token and keep it in a file:",
            "helixloom token > token.txt",
        );
        step(
            ui,
            "2. Start the server there:",
            "helixloom serve --listen 0.0.0.0:7878 --token-file token.txt",
        );
        ui.label(
            RichText::new(
                "3. Type the workstation's address and the token on the left, and press Connect.",
            )
            .color(look.text),
        );
        ui.add_space(6.0);
        ui.label(dim(look, "DEPLOY.md has the whole setup: building for CUDA, choosing GPUs, starting the server at boot and reaching it safely."));
    });
}

fn status_card(ui: &mut Ui, look: &Look, page: &WorkstationPage) {
    let title = page
        .status
        .as_ref()
        .map(|s| s.name.as_str())
        .filter(|n| !n.is_empty())
        .unwrap_or("The workstation");
    titled(ui, look, title, |ui| {
        let Some(s) = &page.status else {
            let (text, colour) = if page.problem.is_empty() {
                ("Asking the workstation…", look.dim)
            } else {
                (page.problem.as_str(), look.bad)
            };
            ui.label(RichText::new(text).color(colour));
            return;
        };
        summary_chips(ui, look, s);
        if s.protocol != PROTOCOL {
            ui.label(
                RichText::new(format!(
                    "This workstation speaks version {} of the protocol and this program version {PROTOCOL}: update the older one.",
                    s.protocol
                ))
                .color(look.bad),
            );
        }
        if !s.gpus.is_empty() {
            ui.add_space(4.0);
            gpu_grid(ui, look, s, &page.jobs);
        }
        ui.add_space(4.0);
        tool_chips(ui, look, s);
        if !page.problem.is_empty() {
            ui.label(
                RichText::new(format!("Could not refresh: {}", page.problem))
                    .size(12.0)
                    .color(look.bad),
            );
        }
    });
}

/// How the workstation is, in a few chips: its version, what it is doing, its GPUs and CPU places, how long results stay.
fn summary_chips(ui: &mut Ui, look: &Look, s: &ServerStatus) {
    ui.horizontal_wrapped(|ui| {
        design::chip(ui, look, &format!("version {}", s.version), Tone::Neutral);
        let busy = if s.running > 0 {
            Tone::Accent
        } else {
            Tone::Neutral
        };
        let running = format!("{} running · {} waiting", s.running, s.queued);
        design::chip(ui, look, &running, busy);
        let gpus = s.gpus.len();
        let tone = if gpus > 0 { Tone::Good } else { Tone::Neutral };
        let plural = if gpus == 1 { "" } else { "s" };
        design::chip(ui, look, &format!("{gpus} GPU{plural}"), tone);
        let cpu = format!("{} jobs at once without a GPU", s.cpu_jobs);
        design::chip(ui, look, &cpu, Tone::Neutral);
        let kept = format!("results kept {}", keep_text(s.keep_hours));
        design::chip(ui, look, &kept, Tone::Neutral);
    });
}

/// Each GPU, with the job it is busy with if there is one (named as the list of jobs names it).
fn gpu_grid(ui: &mut Ui, look: &Look, s: &ServerStatus, jobs: &[JobInfo]) {
    egui::Grid::new("workstation_gpus")
        .num_columns(3)
        .spacing([14.0, 4.0])
        .show(ui, |ui| {
            for g in &s.gpus {
                let number = RichText::new(format!("GPU {}", g.id)).strong();
                ui.label(number.color(look.text));
                ui.label(RichText::new(&g.name).color(look.text));
                match &g.job {
                    Some(id) => {
                        let what = jobs
                            .iter()
                            .find(|j| &j.id == id)
                            .map_or(id.as_str(), |j| j.label.as_str());
                        let busy = format!("busy: {}", shorten(what));
                        design::chip(ui, look, &busy, Tone::Accent);
                    }
                    None => design::chip(ui, look, "free", Tone::Good),
                }
                ui.end_row();
            }
        });
}

/// The programs the workstation has found, and how it folds proteins.
fn tool_chips(ui: &mut Ui, look: &Look, s: &ServerStatus) {
    ui.horizontal_wrapped(|ui| {
        for program in tools::ALL {
            match s.tools.iter().find(|(name, _)| name == program.label) {
                Some((_, version)) => {
                    let found = format!("{} {version}", program.label);
                    design::chip(ui, look, found.trim(), Tone::Good);
                }
                None => {
                    let missing = format!("{} not found", program.label);
                    design::chip(ui, look, &missing, Tone::Neutral);
                }
            }
        }
        if s.fold {
            design::chip(ui, look, "folds proteins on its GPUs", Tone::Good);
        } else {
            design::chip(
                ui,
                look,
                "folds through the ESM Atlas service",
                Tone::Neutral,
            );
        }
    });
}

fn jobs_card(ui: &mut Ui, look: &Look, page: &WorkstationPage, c: &mut Clicks) {
    titled(ui, look, "Jobs", |ui| {
        if page.jobs.is_empty() {
            ui.label(dim(
                look,
                "No jobs yet. Under \"Where jobs run\" choose Workstation or Automatic, then run an alignment, a tree or a pipeline.",
            ));
            return;
        }
        for job in page.jobs.iter().rev().take(LISTED) {
            job_row(ui, look, job, c);
        }
        if page.jobs.len() > LISTED {
            ui.label(dim(
                look,
                format!("and {} older ones", page.jobs.len() - LISTED),
            ));
        }
    });
}

fn job_row(ui: &mut Ui, look: &Look, job: &JobInfo, c: &mut Clicks) {
    ui.horizontal(|ui| {
        design::chip(ui, look, job.status.name(), status_tone(job.status));
        ui.label(RichText::new(shorten(&job.label)).color(look.text))
            .on_hover_text(format!("{} · {}", job.kind.label(), job.id));
        ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
            let remove = if job.status.is_over() {
                "Remove"
            } else {
                "Stop"
            };
            let hover = if job.status.is_over() {
                "Delete the job and its result from the workstation"
            } else {
                "Stop the job"
            };
            if ui.button(remove).on_hover_text(hover).clicked() {
                c.remove = Some((job.id.clone(), !job.status.is_over()));
            }
            if job.status == Status::Done
                && ui
                    .button("Open")
                    .on_hover_text("Fetch the result and show it")
                    .clicked()
            {
                c.open = Some(job.clone());
            }
            ui.label(dim(look, detail(job)));
        });
    });
    if job.status == Status::Failed && !job.error.is_empty() {
        ui.label(RichText::new(&job.error).size(12.0).color(look.bad));
    }
}

fn status_tone(status: Status) -> Tone {
    match status {
        Status::Queued | Status::Cancelled => Tone::Neutral,
        Status::Running => Tone::Accent,
        Status::Done => Tone::Good,
        Status::Failed => Tone::Bad,
    }
}

/// What the list says about a job beside its name.
fn detail(job: &JobInfo) -> String {
    match (job.status, job.gpu) {
        (Status::Queued, _) if job.needs_gpu => "waiting for a GPU".into(),
        (Status::Queued, _) => "waiting".into(),
        (Status::Running, Some(gpu)) => format!("on GPU {gpu}"),
        (Status::Running, None) => "running".into(),
        (Status::Done, _) if job.finished >= job.started && job.started > 0 => {
            format!("took {}", duration_text(job.finished - job.started))
        }
        _ => String::new(),
    }
}

fn duration_text(seconds: u64) -> String {
    match seconds {
        0..=59 => format!("{seconds} s"),
        60..=3599 => format!("{:.1} min", seconds as f64 / 60.0),
        _ => format!("{:.1} h", seconds as f64 / 3600.0),
    }
}

fn keep_text(hours: u64) -> String {
    if hours >= 48 && hours.is_multiple_of(24) {
        format!("{} days", hours / 24)
    } else {
        format!("{hours} h")
    }
}

fn shorten(text: &str) -> String {
    if text.chars().count() <= LABEL_CHARS {
        return text.to_string();
    }
    let cut: String = text.chars().take(LABEL_CHARS - 1).collect();
    format!("{cut}…")
}

#[cfg(test)]
mod tests;
