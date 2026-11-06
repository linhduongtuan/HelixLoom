//! Phylogeny: estimate a tree from the alignment made on the Alignment page, with neighbour-joining, UPGMA,
//! maximum likelihood or Bayesian sampling, on a worker thread so the window stays free; the result is drawn
//! beside the log and becomes the tree the Trees page shows.

use super::treeview::{self, Camera, Options};
use super::workstation::{compute_picker, ComputeView, Place};
use super::Workspace;
use crate::remote::tasks::{Request, Response};
use crate::remote::weight::Weight;
use crate::seq::estimate::{estimate, Analysis, Estimate, EstimateOptions};
use crate::seq::tree::DistanceModel;
use crate::seq::treelayout::lay_out;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::job::{Job, Poll};
use eframe::egui::{self, RichText, ScrollArea, Ui};
use std::time::Duration;
use web_time::Instant;

const OPTIONS_WIDTH: f32 = 320.0;
const WORKING_REDRAW: Duration = Duration::from_millis(150);
const MAX_REPLICATES: usize = 1000;
const GENERATION_RANGE: std::ops::RangeInclusive<usize> = 1_000..=1_000_000;

#[derive(Default)]
pub struct PhylogenyPage {
    job: Option<(Job<Result<Estimate, String>>, Instant)>,
    pub message: String,
    camera: Camera,
    /// Estimate as soon as an alignment is there (an example was loaded).
    pub run_when_ready: bool,
}

impl PhylogenyPage {
    pub fn is_running(&self) -> bool {
        self.job.is_some()
    }
}

#[derive(Default)]
struct Clicks {
    run: bool,
    to_alignment: bool,
    to_trees: bool,
    phylocnn: bool,
    tree_tools: bool,
}

/// The methods the page offers, with what each is good for.
const METHODS: [(Analysis, &str, &str); 4] = [
    (
        Analysis::NeighborJoining,
        "Neighbour-joining",
        "Fast; from distances between the sequences.",
    ),
    (
        Analysis::Upgma,
        "UPGMA",
        "Distances, with the assumption of a steady molecular clock.",
    ),
    (
        Analysis::Likelihood,
        "Maximum likelihood",
        "Models JC69, K80, F81 and HKY85 compared by BIC; the tree improved by swaps.",
    ),
    (
        Analysis::Bayesian,
        "Bayesian",
        "HKY85, two independent runs; branch labels are posterior probabilities.",
    ),
];

impl App {
    /// Start a waiting estimation once the alignment is there, and collect a finished one.
    pub(super) fn phylogeny_tick(&mut self, ctx: &egui::Context) {
        let ready =
            self.platform.alignment.result.is_some() && !self.platform.alignment.is_running();
        if self.platform.phylogeny.run_when_ready && ready && !self.platform.phylogeny.is_running()
        {
            self.platform.phylogeny.run_when_ready = false;
            self.start_estimation(ctx);
        }
        let page = &mut self.platform.phylogeny;
        let Some(polled) = page.job.as_ref().map(|(job, _)| job.poll()) else {
            return;
        };
        match polled {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(Ok(done)) => {
                let seconds = page
                    .job
                    .take()
                    .map_or(0.0, |(_, started)| started.elapsed().as_secs_f64());
                let tips = done.tree.leaf_count();
                page.message = format!("Estimated a tree of {tips} tips in {seconds:.2} s.");
                page.camera = Camera::default();
                self.tree_win.details = done.details;
                self.tree_win.tree = Some(done.tree);
                self.platform.trees.reset_view();
            }
            Poll::Done(Err(why)) => (page.message, page.job) = (why, None),
            Poll::Died => {
                (page.message, page.job) =
                    ("The estimation stopped without an answer.".into(), None)
            }
        }
    }

    fn start_estimation(&mut self, ctx: &egui::Context) {
        let Some(aligned) = &self.platform.alignment.result else {
            self.platform.phylogeny.message = "Align some sequences first.".into();
            return;
        };
        let w = &self.tree_win;
        let options = EstimateOptions {
            analysis: w.analysis,
            model: w.model,
            bootstrap: w.bootstrap,
            generations: w.generations,
        };
        let weight = Weight::of_tree(
            aligned.rows.len(),
            aligned.columns(),
            options.analysis,
            options.bootstrap,
            options.generations,
        );
        let place = self.platform.place(&weight);
        let page = &mut self.platform.phylogeny;
        let Place::There(remote) = place else {
            let msa = aligned.msa();
            page.message = "Estimating…".into();
            page.job = Some((
                Job::spawn(ctx, move || estimate(&msa, &options)),
                Instant::now(),
            ));
            return;
        };
        page.message = "Estimating on the workstation…".into();
        let label = format!(
            "Tree ({}) of {} sequences",
            options.analysis.name(),
            aligned.rows.len()
        );
        let request = Request::Tree {
            alignment: aligned.rows.clone(),
            analysis: options.analysis,
            model: options.model,
            bootstrap: options.bootstrap,
            generations: options.generations,
        };
        page.job = Some((
            Job::spawn(ctx, move || {
                remote.run(&request, &label).and_then(Response::into_tree)
            }),
            Instant::now(),
        ));
    }

    pub(super) fn phylogeny_page(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let look = design::look(self.prefs.dark);
        let mut c = Clicks::default();
        let have_tree = self.tree_win.tree.is_some();
        let summary = self
            .platform
            .alignment
            .result
            .as_ref()
            .map(|a| (a.rows.len(), a.columns(), a.engine.clone()));
        let tree = self.tree_win.tree.clone();
        let details = self.tree_win.details.clone();
        // where an estimate with the settings as they are would run
        let weight = summary
            .as_ref()
            .map_or_else(Weight::default, |(rows, columns, _)| {
                let w = &self.tree_win;
                Weight::of_tree(*rows, *columns, w.analysis, w.bootstrap, w.generations)
            });
        let view = self.platform.workstation.view(&weight);
        let (page, options) = (&mut self.platform.phylogeny, &mut self.tree_win);
        design::page(ui, |ui| {
            design::page_header(
                ui,
                &look,
                "Phylogenetic estimation",
                "Build a tree from an alignment; it opens in Trees",
                |ui| {
                    ui.add_enabled_ui(have_tree, |ui| {
                        c.to_trees = design::primary_button(ui, &look, "Open in Trees").clicked();
                    });
                },
            );
            egui::Panel::left("phylo_options")
                .exact_size(OPTIONS_WIDTH)
                .resizable(false)
                .frame(egui::Frame::new())
                .show(ui, |ui| {
                    ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            design::card(ui, &look, |ui| {
                                let input = (summary.as_ref(), view);
                                options_panel(ui, &look, page, options, input, &mut c)
                            });
                        });
                });
            egui::CentralPanel::default()
                .frame(egui::Frame::new().inner_margin(egui::Margin {
                    left: design::PAGE_MARGIN,
                    ..egui::Margin::ZERO
                }))
                .show(ui, |ui| {
                    result_panel(ui, &look, page, tree.as_ref(), &details)
                });
        });
        if c.run {
            self.start_estimation(ui.ctx());
        }
        if c.to_alignment {
            out.push(Action::Go(Workspace::Alignment));
        }
        if c.to_trees {
            out.push(Action::Go(Workspace::Trees));
        }
        if c.phylocnn {
            self.phylocnn_win.open = true;
        }
        if c.tree_tools {
            self.treetools_win.open = true;
        }
    }
}

fn options_panel(
    ui: &mut Ui,
    look: &Look,
    page: &mut PhylogenyPage,
    w: &mut crate::ui::windows::tree::TreeWin,
    (alignment, view): (Option<&(usize, usize, String)>, ComputeView),
    c: &mut Clicks,
) {
    design::section(ui, look, "Alignment");
    match alignment {
        Some((rows, columns, engine)) => {
            ui.horizontal_wrapped(|ui| {
                design::chip(ui, look, &format!("{rows} sequences"), Tone::Accent);
                design::chip(ui, look, &format!("{columns} columns"), Tone::Neutral);
            });
            ui.label(
                RichText::new(format!("Aligned with {engine}."))
                    .size(12.0)
                    .color(look.dim),
            );
        }
        None => {
            ui.label(RichText::new("There is no alignment yet.").color(look.dim));
        }
    }
    c.to_alignment = ui
        .button(if alignment.is_some() {
            "Change the alignment…"
        } else {
            "Go to Alignment…"
        })
        .clicked();

    design::section(ui, look, "Method");
    for (analysis, label, blurb) in METHODS {
        ui.radio_value(&mut w.analysis, analysis, label);
        ui.label(RichText::new(blurb).size(12.0).color(look.dim));
    }

    design::section(ui, look, "Settings");
    match w.analysis {
        Analysis::Likelihood => {
            ui.label(
                RichText::new("No settings: the model is chosen for you.")
                    .size(12.5)
                    .color(look.dim),
            );
        }
        Analysis::Bayesian => {
            ui.horizontal(|ui| {
                ui.label("Generations per run");
                ui.add(
                    egui::DragValue::new(&mut w.generations)
                        .range(GENERATION_RANGE)
                        .speed(100.0),
                );
            });
        }
        Analysis::NeighborJoining | Analysis::Upgma => {
            ui.horizontal_wrapped(|ui| {
                for m in DistanceModel::ALL {
                    ui.radio_value(&mut w.model, m, m.label());
                }
            });
            ui.horizontal(|ui| {
                ui.label("Bootstrap resamples");
                ui.add(egui::DragValue::new(&mut w.bootstrap).range(0..=MAX_REPLICATES));
            });
        }
    }
    ui.add_space(8.0);
    compute_picker(ui, look, view);
    ui.horizontal(|ui| {
        let ready = alignment.is_some() && !page.is_running();
        c.run = ui
            .add_enabled(
                ready,
                egui::Button::new(RichText::new("Estimate tree").strong()),
            )
            .clicked();
        if page.is_running() {
            ui.spinner();
        }
    });
    if !page.message.is_empty() {
        ui.label(RichText::new(&page.message).size(12.0).color(look.dim));
    }

    design::section(ui, look, "More");
    c.phylocnn = ui
        .button("PhyloCNN epidemic estimates…")
        .on_hover_text("Neural networks that read epidemic parameters from a tree")
        .clicked();
    c.tree_tools = ui
        .button("Tree tools…")
        .on_hover_text("Distances, Robinson-Foulds, midpoint root")
        .clicked();
}

fn result_panel(
    ui: &mut Ui,
    look: &Look,
    page: &mut PhylogenyPage,
    tree: Option<&crate::seq::tree::Node>,
    details: &[String],
) {
    let Some(tree) = tree else {
        design::card(ui, look, |ui| {
            ui.label(
                RichText::new("No tree yet")
                    .size(17.0)
                    .strong()
                    .color(look.text),
            );
            ui.label(
                RichText::new(
                    "Choose a method on the left and press Estimate tree. A preview appears here, with how the tree \
                     was made; the full viewer is on the Trees page.",
                )
                .color(look.dim),
            );
        });
        return;
    };
    if !details.is_empty() {
        design::card(ui, look, |ui| {
            for line in details {
                ui.label(RichText::new(line).size(12.5).color(look.dim));
            }
        });
        ui.add_space(8.0);
    }
    egui::Frame::new()
        .fill(look.surface)
        .stroke(egui::Stroke::new(1.0, look.border))
        .corner_radius(egui::CornerRadius::same(design::CARD_RADIUS))
        .show(ui, |ui| {
            let laid = lay_out(tree, crate::seq::treelayout::Layout::Phylogram);
            let options = Options::default();
            treeview::show(ui, look, &laid, &options, &mut page.camera, &mut None, true);
        });
}
