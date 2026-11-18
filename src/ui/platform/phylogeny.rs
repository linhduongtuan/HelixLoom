//! Phylogeny: estimate a tree from the alignment made on the Alignment page, with neighbour-joining, UPGMA,
//! maximum likelihood or Bayesian sampling, on a worker thread so the window stays free; the result is drawn
//! beside the log and becomes the tree the Trees page shows.

use super::treeview::{self, Camera, Options};
use super::workstation::{run_button, ComputeView, Place};
use super::Workspace;
use crate::pipeline::steps::Aligned;
use crate::provenance::{Asked, Input, Provenance};
use crate::remote::tasks::{Request, Response};
use crate::remote::weight::Weight;
use crate::seq::estimate::{estimate, Analysis, Estimate, EstimateOptions};
use crate::seq::tree::DistanceModel;
use crate::seq::treelayout::lay_out;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::job::{Job, Poll};
use crate::ui::provenance as made_ui;
use eframe::egui::{self, RichText, Ui};
use std::time::Duration;
use web_time::Instant;

const WORKING_REDRAW: Duration = Duration::from_millis(150);
pub(in crate::ui) const MAX_REPLICATES: usize = 1000;
pub(in crate::ui) const GENERATION_RANGE: std::ops::RangeInclusive<usize> = 1_000..=1_000_000;

#[derive(Default)]
pub struct PhylogenyPage {
    job: Option<(Job<Result<Estimate, String>>, Instant)>,
    /// What the estimation that is going was asked, with its settings; and, for the tree that it made (while that is the tree
    /// shown), the record of what made it.
    asked: Option<(Asked, EstimateOptions)>,
    pub made: Option<(Provenance, EstimateOptions)>,
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
    made: made_ui::Asks,
    run: bool,
    example: bool,
    to_alignment: bool,
    to_trees: bool,
    phylocnn: bool,
    species: bool,
    tree_tools: bool,
}

/// The inputs and the settings of an estimation: the rows of the alignment, the method, and the settings that the method uses.
pub fn as_it_is(aligned: &Aligned, o: &EstimateOptions) -> (Vec<Input>, Vec<(String, String)>) {
    let inputs = aligned
        .rows
        .iter()
        .map(|(name, row)| Input::of(name, row))
        .collect();
    let method = METHODS
        .iter()
        .find(|m| m.0 == o.analysis)
        .map_or("", |m| m.1);
    let mut settings = vec![("method".to_string(), method.to_string())];
    match o.analysis {
        Analysis::NeighborJoining | Analysis::Upgma => {
            settings.push(("distances".to_string(), o.model.label().to_string()));
            settings.push(("bootstrap resamples".to_string(), o.bootstrap.to_string()));
        }
        Analysis::Bayesian => {
            settings.push(("generations per run".to_string(), o.generations.to_string()));
        }
        Analysis::Likelihood => {}
    }
    (inputs, settings)
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
                // the methods are this program's own, so it is the only program named
                page.made = page
                    .asked
                    .take()
                    .map(|(asked, options)| (asked.made(Vec::new(), seconds), options));
                page.camera = Camera::default();
                self.tree_win.details = done.details;
                self.tree_win.set_tree(Some(done.tree));
                self.platform.trees.reset_view();
            }
            Poll::Done(Err(why)) => (page.message, page.job, page.asked) = (why, None, None),
            Poll::Died => {
                page.asked = None;
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
        let (inputs, parameters) = as_it_is(aligned, &options);
        let page = &mut self.platform.phylogeny;
        page.asked = Some((
            Asked {
                analysis: "Tree estimation",
                inputs,
                parameters,
                ran_on: if matches!(place, Place::There(_)) {
                    "the workstation"
                } else {
                    "this computer"
                },
            },
            options,
        ));
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
            Job::spawn_async(ctx, async move {
                remote
                    .run_async(&request, &label)
                    .await
                    .and_then(Response::into_tree)
            }),
            Instant::now(),
        ));
    }

    pub(super) fn phylogeny_page(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let look = self.prefs.look();
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
        let aligned = self.platform.alignment.result.as_ref();
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
                        crate::ui::figure::button(ui);
                    });
                },
            );
            design::controls_panel(ui, &look, "phylo_options", |ui| {
                let input = (summary.as_ref(), view);
                options_panel(ui, &look, page, options, input, &mut c)
            });
            design::results_area(ui, |ui| {
                let current = EstimateOptions {
                    analysis: options.analysis,
                    model: options.model,
                    bootstrap: options.bootstrap,
                    generations: options.generations,
                };
                let now = || aligned.map_or_else(Default::default, |a| as_it_is(a, &current));
                result_panel(ui, &look, page, (tree.as_ref(), &details), now, &mut c)
            });
        });
        if c.made.again {
            // with the settings the tree that is shown was made with
            if let Some((_, o)) = &self.platform.phylogeny.made {
                let w = &mut self.tree_win;
                (w.analysis, w.model, w.bootstrap, w.generations) =
                    (o.analysis, o.model, o.bootstrap, o.generations);
            }
        }
        if let (true, Some((made, _))) = (c.made.copy, &self.platform.phylogeny.made) {
            ui.ctx().copy_text(made.to_text());
            self.platform.phylogeny.message = "Copied the record of how the tree was made.".into();
        }
        self.checks_win.open |= c.made.checks;
        if c.run || c.made.again {
            self.start_estimation(ui.ctx());
        }
        if c.to_alignment {
            out.push(Action::Go(Workspace::Alignment));
        }
        if c.to_trees {
            out.push(Action::Go(Workspace::Trees));
        }
        if c.example {
            out.push(Action::Example(Workspace::Phylogeny));
        }
        if c.phylocnn {
            self.phylocnn_win.open = true;
        }
        if c.species {
            self.species_win.open = true;
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
    let state = (alignment.is_some(), page.is_running());
    c.run = run_button(ui, look, view, "Estimate tree", state);
    if !page.message.is_empty() {
        ui.label(RichText::new(&page.message).size(12.0).color(look.dim));
    }

    design::section(ui, look, "More");
    c.phylocnn = ui
        .button("PhyloCNN epidemic estimates…")
        .on_hover_text("Neural networks that read epidemic parameters from a tree")
        .clicked();
    c.species = ui
        .button("Species tree from gene trees…")
        .on_hover_text("The species tree that agrees with most gene trees (the ASTRAL method)")
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
    (tree, details): (Option<&crate::seq::tree::Node>, &[String]),
    now: impl FnOnce() -> (Vec<Input>, Vec<(String, String)>),
    c: &mut Clicks,
) {
    let Some(tree) = tree else {
        let pick = design::empty_state(
            ui,
            look,
            "No tree yet",
            "Choose a method on the left and press Estimate tree. A preview appears here, with how the tree \
             was made; the full viewer is on the Trees page.",
            &["Try the example", "Add sequences in Alignment"],
        );
        match pick {
            Some(0) => c.example = true,
            Some(1) => c.to_alignment = true,
            _ => {}
        }
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
    if let Some((made, _)) = &page.made {
        c.made = made_ui::section(ui, look, made, now, !page.is_running());
        ui.add_space(8.0);
    }
    let laid = lay_out(tree, crate::seq::treelayout::Layout::Phylogram);
    let options = Options::default();
    if let Some(format) = crate::ui::figure::wanted() {
        let mut camera = page.camera;
        page.message = crate::ui::figure::save(format, "tree", ui.available_size(), |ui| {
            treeview::show(ui, look, &laid, &options, &mut camera, &mut None, false);
        });
    }
    design::view_frame(ui, look, |ui| {
        treeview::show(ui, look, &laid, &options, &mut page.camera, &mut None, true);
    });
}
