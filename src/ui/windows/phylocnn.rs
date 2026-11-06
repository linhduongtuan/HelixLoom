//! Tools → PhyloCNN Epidemic Estimates: epidemic parameters, or which transmission model fits,
//! from a tree whose branch lengths are times, with the trained networks of Perez & Gascuel (2025).
//! Several trees (a posterior sample, say) are estimated together and summarised, on the CPU or the GPU.

use super::treetools::{newick_text, read_newick_file};
use crate::models;
use crate::remote::run::run_phylocnn;
use crate::remote::settings::Compute;
use crate::remote::tasks::{CnnResult, Request, Response};
use crate::remote::weight::Weight;
use crate::seq::phylocnn::{
    example_tree, summarize, Estimate, Priors, Spread, Task, MAX_TIPS, MIN_TRAINED_TIPS,
};
use crate::ui::app::App;
use crate::ui::design;
use crate::ui::dialog::{self, Wanted};
use crate::ui::job::{Job, Poll};
use crate::ui::platform::route::Source;
use crate::ui::platform::workstation::{compute_picker, Place};
use eframe::egui::{self, RichText, ScrollArea};
use std::time::Duration;

/// The authors analyse the Zurich tree with a quarter of the infected people sampled.
const EXAMPLE_SAMPLING: f64 = 0.25;
/// Smallest and largest share of infected people the networks were trained on.
const SAMPLING_RANGE: std::ops::RangeInclusive<f64> = 0.01..=1.0;
const TREE_BOX_MAX_HEIGHT: f32 = 110.0;
/// While a run is going the window looks again this often to move its spinner.
const WORKING_REDRAW: Duration = Duration::from_millis(200);
/// While the networks are missing the window looks again for them this often, in seconds.
const MODELS_RECHECK: f64 = 2.0;

/// What a run of the networks found.
struct Report {
    task: Task,
    trees: usize,
    used: usize,
    tips: usize,
    skipped: Vec<(usize, String)>,
    /// The estimates when one tree was used.
    one: Vec<Estimate>,
    /// Their spread over the trees when several were.
    spread: Vec<Spread>,
    seconds: f64,
    engine: String,
}

/// What the window shows of a run, wherever it ran.
fn report(r: CnnResult) -> Report {
    let used = r.per_tree.len();
    Report {
        task: r.task,
        trees: r.trees,
        used,
        tips: r.tips,
        skipped: r.skipped,
        one: if used == 1 {
            r.per_tree[0].clone()
        } else {
            Vec::new()
        },
        spread: if used > 1 {
            summarize(&r.per_tree)
        } else {
            Vec::new()
        },
        seconds: r.seconds,
        engine: r.engine,
    }
}

pub struct PhyloCnnWin {
    pub open: bool,
    pub newick: String,
    pub sampling: f64,
    pub priors: Priors,
    pub task: Task,
    /// Run many trees on the GPU (only offered when there is one).
    pub use_gpu: bool,
    /// The GPU's name once looked for: `Some(None)` when there is no GPU.
    gpu: Option<Option<String>>,
    run_requested: bool,
    job: Option<Job<Result<Report, String>>>,
    /// Rows of the last answer, as text.
    pub table: Vec<Vec<String>>,
    pub header: Vec<&'static str>,
    pub message: String,
    /// Whether the trained networks are installed, and when that was looked at (they are files of the authors', not of this
    /// program; they may be copied in by hand while the window is open).
    models: Option<bool>,
    models_checked: f64,
}

impl Default for PhyloCnnWin {
    fn default() -> Self {
        PhyloCnnWin {
            open: false,
            newick: String::new(),
            sampling: EXAMPLE_SAMPLING,
            priors: Priors::Simulated,
            task: Task::ChooseModel,
            use_gpu: false,
            gpu: None,
            run_requested: false,
            job: None,
            table: Vec::new(),
            header: Vec::new(),
            message: String::new(),
            models: None,
            models_checked: 0.0,
        }
    }
}

impl PhyloCnnWin {
    /// Whether the networks are installed, looked at again every few seconds.
    fn models_ready(&mut self, ctx: &egui::Context) -> bool {
        let now = ctx.input(|i| i.time);
        if self.models.is_none() || now - self.models_checked > MODELS_RECHECK {
            self.models = Some(models::installed(&models::PHYLOCNN));
            self.models_checked = now;
        }
        let ready = self.models == Some(true);
        if !ready {
            ctx.request_repaint_after(Duration::from_secs_f64(MODELS_RECHECK));
        }
        ready
    }

    /// Whether a run is going on.
    pub fn is_running(&self) -> bool {
        self.job.is_some()
    }

    fn set_tree(&mut self, text: String, message: String) {
        self.newick = text;
        self.table.clear();
        self.message = message;
    }

    /// The GPU's name, looked for the first time it is asked.
    fn gpu_name(&mut self) -> Option<&str> {
        #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
        let found = || crate::seq::phylocnn::gpu::device().map(|(_, name)| name);
        #[cfg(not(all(feature = "gpu", not(target_arch = "wasm32"))))]
        let found = || None;
        self.gpu.get_or_insert_with(found).as_deref()
    }

    fn show(task: Task, value: f64) -> String {
        match task {
            Task::ChooseModel => format!("{:.1} %", value * 100.0),
            _ => format!("{value:.3}"),
        }
    }

    fn take(&mut self, report: Report) {
        let task = report.task;
        if report.one.is_empty() {
            self.header = vec!["", "mean", "median", "95 % interval"];
            self.table = report
                .spread
                .iter()
                .map(|s| {
                    vec![
                        s.name.to_string(),
                        Self::show(task, s.mean),
                        Self::show(task, s.median),
                        format!(
                            "{} to {}",
                            Self::show(task, s.low),
                            Self::show(task, s.high)
                        ),
                    ]
                })
                .collect();
        } else {
            self.header = vec!["", "estimate"];
            self.table = report
                .one
                .iter()
                .map(|e| vec![e.name.to_string(), Self::show(task, e.value)])
                .collect();
        }
        let mut notes = vec![if report.trees == 1 {
            format!("{} tips.", report.tips)
        } else {
            format!(
                "{} of {} trees used, {} tips in the first.",
                report.used, report.trees, report.tips
            )
        }];
        if !report.skipped.is_empty() {
            let (first, why) = &report.skipped[0];
            notes.push(format!(
                "{} left out, the first being tree {}: {why}",
                report.skipped.len(),
                first + 1
            ));
        }
        if !(MIN_TRAINED_TIPS..=MAX_TIPS).contains(&report.tips) {
            notes.push(format!(
                "The networks learned from trees of {MIN_TRAINED_TIPS} to {MAX_TIPS} tips, so read this as a rough guide."
            ));
        }
        if report.used > 1 {
            notes.push(format!(
                "Networks ran on the {} in {:.2} s.",
                report.engine, report.seconds
            ));
        }
        self.message = notes.join(" ");
    }
}

impl App {
    /// Show what the networks found on the workstation, and open the window to show it.
    pub fn show_phylocnn_result(&mut self, result: CnnResult) {
        self.phylocnn_win.open = true;
        self.phylocnn_win.take(report(result));
    }

    /// The Zurich HIV tree with the settings the authors use for it; where that file is not installed, a made-up tree of
    /// as many tips.
    pub fn load_phylocnn_example(&mut self) {
        let (tree, zurich) = example_tree();
        let w = &mut self.phylocnn_win;
        (w.sampling, w.priors, w.task) = (EXAMPLE_SAMPLING, Priors::Hiv, Task::Bdss);
        let note = if zurich {
            "Loaded the Zurich HIV tree (200 sequences, Rasmussen et al. 2017). Press Estimate."
        } else {
            "The Zurich HIV tree is not installed, so this is a made-up tree of 200 tips. Press Estimate."
        };
        w.set_tree(tree.trim().to_string(), note.into());
    }

    /// Ask for the folder that holds the authors' files, and copy them into the models folder.
    fn install_phylocnn_models(&mut self) {
        let w = &mut self.phylocnn_win;
        if dialog::IN_BROWSER {
            w.message = "The networks are not available in the browser yet.".into();
            return;
        }
        let Some(source) = dialog::FileDialog::new()
            .set_title("The folder with the PhyloCNN network files")
            .pick_folder()
        else {
            return;
        };
        let Some(dest) = models::install_dir() else {
            w.message = "There is no folder to keep the files in on this computer.".into();
            return;
        };
        w.message = match models::install_from(&source, &dest, &models::PHYLOCNN) {
            Ok(count) => format!("Installed {count} files into {}.", dest.display()),
            Err(why) => why,
        };
        w.models = None;
    }

    fn open_phylocnn_tree(&mut self) {
        if let Some(read) = read_newick_file(Wanted::PhyloCnnTree) {
            self.use_phylocnn_tree(read);
        }
    }

    /// Put a Newick file that was read into the window's box, or say why it could not be read.
    fn use_phylocnn_tree(&mut self, read: Result<(String, String), String>) {
        let w = &mut self.phylocnn_win;
        match read {
            Ok((name, text)) => w.set_tree(text, format!("Loaded {name}. Press Estimate.")),
            Err(e) => w.message = e,
        }
    }

    /// A Newick file for the window, from a source the browser brought.
    pub(in crate::ui) fn load_phylocnn_tree(&mut self, file: Source) {
        self.use_phylocnn_tree(newick_text(file));
    }

    /// Ask for the networks to be run on the trees in the box; the window starts the run on its next draw.
    pub fn run_phylocnn(&mut self) {
        self.phylocnn_win.run_requested = true;
    }

    fn start_phylocnn(&mut self, ctx: &egui::Context) {
        let remote_gpu = self.platform.workstation.gpu_name().is_some();
        let w = &self.phylocnn_win;
        let place = self
            .platform
            .place(&Weight::of_trees(&w.newick, w.use_gpu && remote_gpu));
        let w = &mut self.phylocnn_win;
        w.run_requested = false;
        w.table.clear();
        let (newick, sampling, task, priors) = (w.newick.clone(), w.sampling, w.task, w.priors);
        let Place::There(remote) = place else {
            let gpu = w.use_gpu && w.gpu_name().is_some();
            w.message = "Working…".into();
            w.job = Some(Job::spawn(ctx, move || {
                run_phylocnn(&newick, sampling, task, priors, gpu).map(report)
            }));
            return;
        };
        w.message = "Working on the workstation…".into();
        let request = Request::PhyloCnn {
            newick,
            sampling,
            task,
            priors,
            gpu: w.use_gpu && remote_gpu,
        };
        w.job = Some(Job::spawn(ctx, move || {
            remote
                .run(&request, "PhyloCNN estimates")
                .and_then(Response::into_phylocnn)
                .map(report)
        }));
    }

    /// Start the run that was asked for and take its answer when it comes, whether the window is showing or not: a run
    /// whose window was closed still ends (the answer waits in the window), and does not go on for ever as "running".
    pub fn tick_phylocnn(&mut self, ctx: &egui::Context) {
        if self.phylocnn_win.run_requested && self.phylocnn_win.job.is_none() {
            self.start_phylocnn(ctx);
        }
        if let Some(job) = &self.phylocnn_win.job {
            match job.poll() {
                Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
                Poll::Done(outcome) => {
                    self.phylocnn_win.job = None;
                    match outcome {
                        Ok(report) => self.phylocnn_win.take(report),
                        Err(why) => self.phylocnn_win.message = why,
                    }
                }
                Poll::Died => {
                    self.phylocnn_win.job = None;
                    self.phylocnn_win.message = "The run stopped without an answer.".into();
                }
            }
        }
    }

    pub fn phylocnn_window(&mut self, ctx: &egui::Context) {
        self.tick_phylocnn(ctx);
        if !self.phylocnn_win.open {
            return;
        }
        let (mut open, mut run, mut load, mut example, mut install) =
            (true, false, false, false, false);
        let models_ready = self.phylocnn_win.models_ready(ctx);
        // where the run would go with the GPU as ticked: it is the workstation's GPU that is on offer when jobs may go there
        let there_gpu = self.platform.workstation.gpu_name().map(str::to_string);
        let weight = Weight::of_trees(
            &self.phylocnn_win.newick,
            self.phylocnn_win.use_gpu && there_gpu.is_some(),
        );
        let view = self.platform.workstation.view(&weight);
        let look = design::look(self.prefs.dark);
        crate::ui::design::window("PhyloCNN Epidemic Estimates")
            .open(&mut open)
            .default_width(560.0)
            .show(ctx, |ui| {
                let w = &mut self.phylocnn_win;
                ui.label(
                    "A rooted tree whose branch lengths are times (years, say), with one tip per \
                     sampled infection. Networks trained on simulated outbreaks estimate how fast \
                     the infection spreads, or which transmission model fits best. Several trees, \
                     one after the other, are estimated together and summarised.",
                );
                if !models_ready {
                    ui.horizontal_wrapped(|ui| {
                        ui.colored_label(
                            ui.visuals().warn_fg_color,
                            "The trained networks are not installed: they belong to the authors of PhyloCNN and are not \
                             part of HelixLoom.",
                        );
                        install = ui.button("Install from folder…").clicked();
                    });
                }
                ui.horizontal(|ui| {
                    load = ui.button("Open tree file…").clicked();
                    example = ui
                        .button("Zurich HIV example")
                        .on_hover_text("200 HIV-1 sequences from the Swiss HIV Cohort Study")
                        .clicked();
                    if ui.button("Clear").clicked() {
                        w.set_tree(String::new(), String::new());
                    }
                });
                ScrollArea::vertical()
                    .max_height(TREE_BOX_MAX_HEIGHT)
                    .show(ui, |ui| {
                        ui.add(
                            egui::TextEdit::multiline(&mut w.newick)
                                .font(egui::TextStyle::Monospace)
                                .hint_text("Paste a Newick tree such as ((A:1,B:1):2,C:3);")
                                .desired_rows(4)
                                .desired_width(f32::INFINITY),
                        );
                    });
                egui::Grid::new("phylocnn_settings")
                    .num_columns(2)
                    .show(ui, |ui| {
                        ui.label("Sampling probability");
                        ui.add(
                            egui::DragValue::new(&mut w.sampling)
                                .range(SAMPLING_RANGE)
                                .speed(0.01)
                                .fixed_decimals(2),
                        )
                        .on_hover_text("The share of all infected people that the tree contains");
                        ui.end_row();
                        ui.label("Trained with");
                        egui::ComboBox::from_id_salt("phylocnn_priors")
                            .selected_text(w.priors.label())
                            .show_ui(ui, |ui| {
                                for p in Priors::ALL {
                                    ui.selectable_value(&mut w.priors, p, p.label());
                                }
                            });
                        ui.end_row();
                    });
                if !w.task.trained_with(w.priors) {
                    w.task = Task::Bdss; // trained with every set of priors
                }
                for t in Task::ALL {
                    let usable = t.trained_with(w.priors);
                    let radio = egui::RadioButton::new(w.task == t, t.label());
                    if ui.add_enabled(usable, radio).clicked() {
                        w.task = t;
                    }
                }
                compute_picker(ui, &look, view);
                // ticking the GPU is what sends a small job to the workstation when the choice is Automatic
                let gpu = match view.mode {
                    Compute::Workstation => there_gpu,
                    Compute::Automatic => there_gpu.or_else(|| w.gpu_name().map(str::to_string)),
                    Compute::Local => w.gpu_name().map(str::to_string),
                };
                if let Some(gpu) = gpu {
                    ui.horizontal(|ui| {
                        ui.label("Many trees run on:");
                        ui.radio_value(&mut w.use_gpu, false, "CPU, every core");
                        ui.radio_value(&mut w.use_gpu, true, gpu).on_hover_text(
                            "Same answers. The GPU builds its programs the first time (about 15 s \
                             here), and for these small networks the CPU is still faster.",
                        );
                    });
                }
                ui.horizontal(|ui| {
                    run = ui
                        .add_enabled(
                            w.job.is_none(),
                            egui::Button::new(RichText::new("Estimate").strong()),
                        )
                        .clicked();
                    if w.job.is_some() {
                        ui.spinner();
                    }
                    ui.label(RichText::new(w.priors.note()).small());
                });
                if !w.message.is_empty() {
                    ui.label(&w.message);
                }
                if !w.table.is_empty() {
                    egui::Grid::new("phylocnn_results")
                        .striped(true)
                        .show(ui, |ui| {
                            for h in &w.header {
                                ui.label(RichText::new(*h).strong());
                            }
                            ui.end_row();
                            for row in &w.table {
                                for (k, cell) in row.iter().enumerate() {
                                    if k == 0 {
                                        ui.label(cell);
                                    } else {
                                        ui.label(RichText::new(cell).monospace().strong());
                                    }
                                }
                                ui.end_row();
                            }
                        });
                    ui.label(
                        RichText::new("Times come back in the units of your branch lengths.")
                            .small(),
                    );
                }
            });
        if load {
            self.open_phylocnn_tree();
        }
        if install {
            self.install_phylocnn_models();
        }
        if example {
            self.load_phylocnn_example();
        }
        if run {
            self.run_phylocnn();
        }
        self.phylocnn_win.open = open;
    }
}
