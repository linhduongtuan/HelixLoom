//! Workstation: connect this computer to the one that does the heavy work. Jobs sent there (alignments, trees, pipelines,
//! PhyloCNN estimates, structure predictions) run on its processors and GPUs and wait there until they are fetched; this
//! page says how the workstation is, lists its jobs and opens a result in the page it belongs to.

mod compute;

use super::Workspace;
use crate::pipeline::tools;
use crate::remote::client::{Remote, ServerStatus};
use crate::remote::jobs::{JobInfo, Status};
use crate::remote::settings::{self, Compute, GpusPerJob, Saved};
use crate::remote::tasks::Response;
use crate::remote::weight::Weight;
use crate::remote::PROTOCOL;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::dialog::IN_BROWSER;
use crate::ui::fonts::FontWeight;
use crate::ui::job::{Job, Poll};
use crate::ui::session::Shelf;
use eframe::egui::{self, Align, Layout, RichText, ScrollArea, Ui};
use std::path::PathBuf;
use std::time::Duration;
use web_time::Instant;

pub use compute::Place;
use compute::{choice_id, heavy_jobs_section};
pub(in crate::ui) use compute::{compute_picker, run_button, ComputeView};

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
    /// How many GPUs of the workstation one job is given at once, when it can use several.
    pub gpus: GpusPerJob,
    status: Option<ServerStatus>,
    jobs: Vec<JobInfo>,
    /// Why the last request failed; empty when it worked.
    problem: String,
    message: String,
    work: Option<Job<Result<Done, String>>>,
    pending: Option<Pending>,
    looked: Option<Instant>,
    /// Where the connection is kept between runs.
    keep: Keep,
}

/// What a page says when the browser would not keep the connection for its tab.
const NOT_KEPT_FOR_THE_TAB: &str =
    "This browser would not keep the connection for the tab: it is asked for again when the page is loaded again.";

/// What the page says about who can see the token and where it is kept. A tab's table is not as private as it sounds, and
/// the words say what is true of it: a tab that was duplicated has a copy of its own, and a browser can bring a closed tab
/// back with what it held.
fn privacy_text(kept_by_the_tab: bool, in_a_browser: bool) -> &'static str {
    if kept_by_the_tab {
        "The address and the token are kept by this browser for this tab, so that loading the page again does not ask \
         for them. Forget removes them; so does closing the tab, unless the browser brings the tab back (a tab that was \
         duplicated has a copy of its own). On a computer that others use, press Forget before you leave. Unless the \
         address starts with https://, the token and your sequences travel as they are: use a private network or an \
         SSH tunnel (DEPLOY.md shows how)."
    } else if in_a_browser {
        "The token stays in this page only: it is not kept, and loading the page again forgets it. Unless the address \
         starts with https://, the token and your sequences travel as they are: use a private network or an SSH \
         tunnel (DEPLOY.md shows how)."
    } else {
        "The connection is not encrypted: the token and your sequences travel as they are. Use it on a private \
         network or through an SSH tunnel (DEPLOY.md shows how)."
    }
}

/// What the pointer is told over the Forget button.
fn forget_hint(kept_by_the_tab: bool) -> &'static str {
    if kept_by_the_tab {
        "Remove the address and the token from this tab, now"
    } else {
        "Remove the address and the token from this computer"
    }
}

/// Where the connection (the address, the token, where jobs run) is kept from one run of the program to the next.
#[derive(Default)]
enum Keep {
    /// Nowhere: a picture, a test, a page whose browser gives it no storage.
    #[default]
    Nowhere,
    /// In the settings folder of the computer, in a file that only the person may read.
    Folder(PathBuf),
    /// In the table a browser gives one tab: the tab's alone, there when the page is loaded again in it, and gone when
    /// the tab is closed. The token is a password, and the owner of the program said (17 November 2026) that a page may
    /// keep it that long and no longer: it is never put in the table of the site, which stays.
    // only a page, or a test, has a tab
    #[cfg_attr(not(any(target_arch = "wasm32", test)), allow(dead_code))]
    Tab(std::rc::Rc<dyn Shelf>),
}

impl Keep {
    /// The place a page has when it is given no settings folder: its tab's table in a browser, none anywhere else.
    #[cfg(target_arch = "wasm32")]
    fn of_a_page() -> Keep {
        crate::ui::session::tab_shelf().map_or(Keep::Nowhere, Keep::Tab)
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn of_a_page() -> Keep {
        Keep::Nowhere
    }

    /// What was kept; nothing where nothing is.
    fn load(&self) -> Saved {
        match self {
            Keep::Nowhere => Saved::default(),
            Keep::Folder(dir) => settings::load_from(dir),
            Keep::Tab(tab) => tab
                .get(settings::FILE)
                .map(|text| Saved::parse(&text))
                .unwrap_or_default(),
        }
    }

    /// Keep `saved`, or nothing of a connection when it is none (the workstation was forgotten): a tab's table then has
    /// no entry at all, and a file is written empty of an address and a token.
    fn save(&self, saved: &Saved) -> Result<(), String> {
        match self {
            Keep::Nowhere => Ok(()),
            Keep::Folder(dir) => settings::save_to(dir, saved),
            Keep::Tab(tab) if *saved == Saved::default() => {
                tab.remove(settings::FILE);
                Ok(())
            }
            // what the table says of why (it speaks of the work that is kept in the site's table) is not what is true
            // here: only the keeping of the connection failed. An entry from before is then not what the person has
            // just typed: it is removed, or the page would come back connected to where it was before
            Keep::Tab(tab) => tab.set(settings::FILE, &saved.to_text()).map_err(|_| {
                tab.remove(settings::FILE);
                NOT_KEPT_FOR_THE_TAB.to_string()
            }),
        }
    }

    fn is_a_tab(&self) -> bool {
        matches!(self, Keep::Tab(_))
    }
}

/// The address the page itself was loaded from, in a browser; `None` in a program, and for a page opened from a file.
#[cfg(target_arch = "wasm32")]
fn page_origin() -> Option<String> {
    web_sys::window()?
        .location()
        .origin()
        .ok()
        .filter(|origin| origin.starts_with("http"))
}

#[cfg(not(target_arch = "wasm32"))]
fn page_origin() -> Option<String> {
    None
}

impl WorkstationPage {
    /// Whether the workstation that was connected has answered, and not failed since: its dot on the top bar is green.
    pub fn answering(&self) -> bool {
        self.status.is_some() && self.problem.is_empty()
    }

    /// The page as it was left last time: connected again, and with the choice of where jobs run, if there was one. With
    /// no folder it is a page of a browser, which keeps the connection for its tab, or a picture, which keeps nothing.
    pub fn open(dir: Option<PathBuf>) -> WorkstationPage {
        WorkstationPage::opened(dir.map_or_else(Keep::of_a_page, Keep::Folder))
    }

    /// A page of a browser whose tab has the table `tab`.
    #[cfg(test)]
    pub(in crate::ui) fn on_tab(tab: std::rc::Rc<dyn Shelf>) -> WorkstationPage {
        WorkstationPage::opened(Keep::Tab(tab))
    }

    fn opened(keep: Keep) -> WorkstationPage {
        let mut saved = keep.load();
        let active = Remote::new(&saved.address, &saved.token)
            .ok()
            .map(|r| r.with_gpus(saved.gpus));
        // what a tab keeps was a connection when the page kept it (`connect` checks what is typed first). An entry that
        // makes none is not this page's doing (any page of the site that is loaded in the tab can write there): it is
        // removed, and nothing of it is shown, so that no token is left where the page has no button to forget it with
        if let (None, Keep::Tab(tab)) = (&active, &keep) {
            tab.remove(settings::FILE);
            saved = Saved::default();
        }
        WorkstationPage {
            // in a page that the workstation's own server handed out, the workstation is where the page came from
            address: if saved.address.is_empty() {
                page_origin().unwrap_or_default()
            } else {
                saved.address
            },
            token: saved.token,
            compute: if active.is_some() {
                saved.compute
            } else {
                Compute::Local
            },
            gpus: saved.gpus,
            active,
            keep,
            ..WorkstationPage::default()
        }
    }

    /// How many GPUs the workstation has, as far as it has said: none until it has answered.
    pub fn gpu_count(&self) -> usize {
        self.status.as_ref().map_or(0, |s| s.gpus.len())
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
    /// Whether a workstation answers, and how many GPUs it has: what a report of a problem may hold of it. Never its address.
    pub fn for_diagnostics(&self) -> String {
        if !self.answering() {
            return "not connected".to_string();
        }
        match self.gpu_count() {
            0 => "connected (no GPU)".to_string(),
            1 => "connected (1 GPU)".to_string(),
            n => format!("connected ({n} GPUs)"),
        }
    }

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
        let saved = self.active.as_ref().map_or_else(Saved::default, |r| Saved {
            address: r.address().to_string(),
            token: r.token().to_string(),
            compute: self.compute,
            gpus: self.gpus,
        });
        if let Err(why) = self.keep.save(&saved) {
            self.message = why;
        }
    }

    /// Use the address and token as typed, and look at what is there.
    pub(super) fn connect(&mut self) {
        match Remote::new(&self.address, &self.token) {
            Ok(remote) => {
                self.active = Some(remote.with_gpus(self.gpus));
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
            keep: std::mem::take(&mut self.keep),
            // in a page that the workstation handed out, its own address is still the one to use
            address: page_origin().unwrap_or_default(),
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

    /// Choose how many GPUs one job is given at once (when it can use several): the jobs sent from now on ask for that many.
    pub(super) fn set_gpus(&mut self, gpus: GpusPerJob) {
        self.gpus = gpus;
        self.active = self.active.take().map(|r| r.with_gpus(gpus));
        self.save();
    }

    /// Ask the workstation how it is and what its jobs are doing.
    fn look(&mut self, ctx: &egui::Context) {
        let Some(remote) = self.active.clone() else {
            return;
        };
        self.looked = Some(Instant::now());
        self.work = Some(Job::spawn_async(ctx, async move {
            Ok(Done::Looked(
                remote.status_async().await?,
                remote.jobs_async().await?,
            ))
        }));
    }

    fn fetch(&mut self, ctx: &egui::Context, job: JobInfo) {
        let Some(remote) = self.active.clone() else {
            return;
        };
        self.message = format!("Fetching \"{}\"…", job.label);
        self.work = Some(Job::spawn_async(ctx, async move {
            let response = remote.result_async(job.kind, &job.id).await?;
            Ok(Done::Fetched(job, response))
        }));
    }

    fn remove(&mut self, ctx: &egui::Context, id: String, stopping: bool) {
        let Some(remote) = self.active.clone() else {
            return;
        };
        self.work = Some(Job::spawn_async(ctx, async move {
            remote.delete_async(&id).await?;
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
    /// How many GPUs one job is to be given from now on.
    gpus: Option<GpusPerJob>,
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
                self.show_tree(estimate.tree, &note);
                self.tree_win.details = estimate.details;
                self.workspace = Workspace::Trees;
            }
            Response::Pipeline(outcome) => {
                self.platform.pipelines.show_fetched(outcome, note);
                self.workspace = Workspace::Pipelines;
            }
            Response::PhyloCnn(result) => self.show_phylocnn_result(result),
            Response::Species(result) => self.show_species_result(result),
            Response::Motif(result) => self.show_motif_result(result),
            Response::Domains(result) => self.show_domains_result(result),
            Response::Flow(result) => self.show_flow_result(result, &note),
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
        let look = self.prefs.look();
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
            design::controls_panel(ui, &look, "workstation_options", |ui| {
                connection_panel(ui, &look, page, &mut c)
            });
            design::results_area(ui, |ui| overview(ui, &look, page, &mut c));
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
        if let Some(gpus) = c.gpus {
            page.set_gpus(gpus);
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
            .hint_text(if IN_BROWSER {
                "http://workstation.local:7878"
            } else {
                "workstation.local:7878"
            })
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
            let from = forget_hint(page.keep.is_a_tab());
            c.forget = ui.button("Forget").on_hover_text(from).clicked();
        }
    });
    if !page.problem.is_empty() {
        ui.label(RichText::new(&page.problem).size(12.0).color(look.bad));
    } else if !page.message.is_empty() {
        ui.label(dim(look, &page.message));
    }

    heavy_jobs_section(ui, look, page, c);

    design::section(ui, look, "Keep it private");
    let privacy = privacy_text(page.keep.is_a_tab(), IN_BROWSER);
    ui.label(RichText::new(privacy).size(11.5).color(look.warn));
}

fn overview(ui: &mut Ui, look: &Look, page: &WorkstationPage, c: &mut Clicks) {
    ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            if page.active.is_none() {
                return setup_card(ui, look);
            }
            status_card(ui, look, page);
            ui.add_space(design::CARD_GAP);
            jobs_card(ui, look, page, c);
        });
}

/// A card with a title.
fn setup_card(ui: &mut Ui, look: &Look) {
    design::titled_card(ui, look, "Set up a workstation", |ui| {
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
        if IN_BROWSER {
            step(
                ui,
                "2. Start the server there, and let it hand out this page:",
                "helixloom serve --listen 0.0.0.0:7878 --token-file token.txt --web-dir dist",
            );
            ui.label(
                RichText::new(
                    "3. Open http://HOST:7878/ in this browser: the address on the left is then filled in. Type the token \
                     and press Connect.",
                )
                .color(look.text),
            );
            ui.label(dim(look, "A page that is hosted somewhere else can call the server too, if the server is started with --allow-origin and that page's address: then replace the address on the left, which is that of the page, with the workstation's."));
        } else {
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
        }
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
    design::titled_card(ui, look, title, |ui| {
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
                let number = RichText::new(format!("GPU {}", g.id)).semibold();
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
    design::titled_card(ui, look, "Jobs", |ui| {
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
    match (job.status, job.gpus.as_slice()) {
        (Status::Queued, _) if job.gpus_wanted > 1 => {
            format!("waiting for {} GPUs", job.gpus_wanted)
        }
        (Status::Queued, _) if job.needs_gpu() => "waiting for a GPU".into(),
        (Status::Queued, _) => "waiting".into(),
        (Status::Running, [one]) => format!("on GPU {one}"),
        (Status::Running, [_, ..]) => format!("on GPUs {}", numbers(&job.gpus)),
        (Status::Running, []) => "running".into(),
        (Status::Done, _) if job.finished >= job.started && job.started > 0 => {
            format!("took {}", duration_text(job.finished - job.started))
        }
        _ => String::new(),
    }
}

/// GPU numbers for a sentence: `0, 1, 3`.
fn numbers(ids: &[usize]) -> String {
    let words: Vec<String> = ids.iter().map(usize::to_string).collect();
    words.join(", ")
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
