//! Pipelines: the sequences through alignment, a tree and primer pairs in one run, with MAFFT, IQ-TREE and primer3
//! when they are installed and the built-in methods when they are not; and the same steps saved as a Snakemake
//! and a Nextflow workflow, for running anywhere those programs are.

use super::home::program_chip;
use super::workstation::{compute_picker, ComputeView, Place};
use super::Workspace;
use crate::pipeline::export::{write_folder, Plan, PrimerSpec};
use crate::pipeline::steps::{fasta_text, Named, PrimerPair};
use crate::pipeline::tools::Installed;
use crate::pipeline::{self, Available, Choice, Outcome};
use crate::remote::tasks::{Request, Response};
use crate::remote::weight::Weight;
use crate::seq::newick::NwTree;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::dialog;
use crate::ui::job::{Job, Poll};
use crate::ui::windows::tree::open_sequences;
use eframe::egui::{self, RichText, ScrollArea, Ui};
use std::time::Duration;

const OPTIONS_WIDTH: f32 = 330.0;
const RESULT_GAP: f32 = 10.0;
/// IQ-TREE's smallest number of ultrafast bootstrap replicates.
const BOOTSTRAP_REPLICATES: usize = 1000;
const WORKING_REDRAW: Duration = Duration::from_millis(200);
/// The example's primers go around these bases of its first isolate: where they start, and how many.
const EXAMPLE_TARGET: (usize, usize) = (200, 60);
/// How many of the sequences in use are listed by name.
const LISTED: usize = 3;
/// Rows of the tree's Newick text that are shown.
const NEWICK_ROWS: usize = 3;
/// The programs that run the steps, and the ones that run an exported workflow.
const ENGINES: [(&str, &str); 3] = [
    ("MAFFT", "MAFFT: built-in aligner"),
    ("IQ-TREE", "IQ-TREE: built-in tree search"),
    ("primer3", "primer3: built-in port"),
];
const RUNNERS: [&str; 2] = ["Snakemake", "Nextflow"];

/// Where the sequences of a run come from.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum Source {
    /// Every sequence open in Sequences.
    #[default]
    Open,
    /// Six made-up isolates, so the page can be tried without data.
    Example,
}

pub struct PipelinesPage {
    pub source: Source,
    pub align: bool,
    pub tree: bool,
    pub primers: bool,
    pub use_installed: bool,
    pub bootstrap: bool,
    /// Run at the next frame (the example was asked for).
    pub run_requested: bool,
    job: Option<Job<Result<Outcome, String>>>,
    pub outcome: Option<Outcome>,
    pub message: String,
}

impl Default for PipelinesPage {
    fn default() -> Self {
        PipelinesPage {
            source: Source::default(),
            align: true,
            tree: true,
            primers: true,
            use_installed: true,
            bootstrap: false,
            run_requested: false,
            job: None,
            outcome: None,
            message: String::new(),
        }
    }
}

impl PipelinesPage {
    pub fn is_running(&self) -> bool {
        self.job.is_some()
    }

    /// The steps to run. A tree needs an alignment, so choosing a tree aligns too.
    fn choice(&self) -> Choice {
        Choice {
            align: self.align || self.tree,
            tree: self.tree,
            primers: self.primers,
            use_installed: self.use_installed,
            bootstrap: if self.bootstrap {
                BOOTSTRAP_REPLICATES
            } else {
                0
            },
        }
    }
}

/// The six isolates of the example.
fn example_sequences() -> Vec<Named> {
    crate::sample::tree_demo()
        .into_iter()
        .map(|d| (d.name, d.seq))
        .collect()
}

/// Where the primers go, in words.
fn target_text(target: Option<(usize, usize)>) -> String {
    target.map_or("anywhere along it".to_string(), |(start, length)| {
        format!("around bases {} to {}", start + 1, start + length)
    })
}

#[derive(Default)]
struct Clicks {
    run: bool,
    export: bool,
    to_alignment: bool,
    to_trees: bool,
    copy_alignment: bool,
    copy_newick: bool,
}

impl App {
    /// The sequences of a run, and the bases of the first one the primers go around (start and length).
    fn pipeline_input(&self) -> (Vec<Named>, Option<(usize, usize)>) {
        match self.platform.pipelines.source {
            Source::Open => (open_sequences(&self.docs), self.pipeline_target()),
            Source::Example => (example_sequences(), self.pipeline_target()),
        }
    }

    fn pipeline_target(&self) -> Option<(usize, usize)> {
        match self.platform.pipelines.source {
            Source::Open => self
                .docs
                .first()
                .and_then(|d| d.selection())
                .map(|(a, b)| (a, b - a)),
            Source::Example => Some(EXAMPLE_TARGET),
        }
    }

    /// The name and length of each sequence a run would use.
    fn pipeline_listing(&self) -> Vec<(String, usize)> {
        match self.platform.pipelines.source {
            Source::Open => self
                .docs
                .iter()
                .map(|d| (d.doc.name.clone(), d.doc.len()))
                .collect(),
            Source::Example => example_sequences()
                .into_iter()
                .map(|(name, bases)| (name, bases.len()))
                .collect(),
        }
    }

    /// How big a run of the steps ticked over the sequences in `listing` is.
    fn pipeline_weight(&self, listing: &[(String, usize)]) -> Weight {
        let letters = listing.iter().map(|(_, length)| length).sum();
        Weight::of_pipeline(listing.len(), letters, &self.platform.pipelines.choice())
    }

    /// The programs on this computer, for a run that stays here (a run on the workstation uses its own).
    fn available(&self) -> Available {
        let p = &self.platform;
        Available {
            mafft: p.installed("MAFFT"),
            iqtree: p.installed("IQ-TREE"),
            primer3: p.installed("primer3"),
        }
    }

    /// The same steps as a workflow file would name them: the program's own name when it was found.
    fn plan(&self) -> Plan {
        let command = |label: &str, fallback: &str| {
            self.platform
                .installed(label)
                .and_then(|p| p.path.file_name().map(|n| n.to_string_lossy().into_owned()))
                .unwrap_or_else(|| fallback.to_string())
        };
        let choice = self.platform.pipelines.choice();
        Plan {
            align: choice.align,
            tree: choice.tree,
            primers: choice.primers,
            bootstrap: choice.bootstrap,
            mafft: command("MAFFT", "mafft"),
            iqtree: command("IQ-TREE", "iqtree3"),
            primer3: command("primer3", "primer3_core"),
        }
    }

    /// Start a requested run once the programs have been looked for, and collect a finished one.
    pub(super) fn pipelines_tick(&mut self, ctx: &egui::Context) {
        let p = &self.platform;
        // wait until the programs are known where the run would go
        let ready = p.pipelines.run_requested
            && !p.pipelines.is_running()
            && p.programs_known_at(&p.place(&self.pipeline_weight(&self.pipeline_listing())));
        if ready {
            self.start_pipeline(ctx);
        }
        let page = &mut self.platform.pipelines;
        let Some(polled) = page.job.as_ref().map(Job::poll) else {
            return;
        };
        match polled {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(done) => {
                page.job = None;
                match done {
                    Ok(outcome) => (page.message, page.outcome) = ("Done.".into(), Some(outcome)),
                    Err(why) => page.message = why,
                }
            }
            Poll::Died => {
                page.job = None;
                page.message = "The run stopped without an answer.".into();
            }
        }
    }

    fn start_pipeline(&mut self, ctx: &egui::Context) {
        let (seqs, target) = self.pipeline_input();
        let available = self.available();
        let settings = self.p3_win.settings.clone();
        let choice = self.platform.pipelines.choice();
        let letters = seqs.iter().map(|s| s.1.len()).sum();
        let place = self
            .platform
            .place(&Weight::of_pipeline(seqs.len(), letters, &choice));
        let page = &mut self.platform.pipelines;
        page.run_requested = false;
        page.outcome = None;
        if seqs.is_empty() {
            page.message =
                "There is no sequence: open one in Sequences, or use the example.".into();
            return;
        }
        if choice.align && seqs.len() < 2 {
            page.message = "Aligning needs at least two sequences.".into();
            return;
        }
        let Place::There(remote) = place else {
            page.message = "Working…".into();
            page.job = Some(Job::spawn(ctx, move || {
                pipeline::run(&seqs, target, &settings, &choice, &available)
            }));
            return;
        };
        page.message = "Working on the workstation…".into();
        let label = format!("Pipeline on {} sequences", seqs.len());
        let request = Request::Pipeline {
            seqs,
            target,
            settings,
            choice,
        };
        page.job = Some(Job::spawn(ctx, move || {
            remote
                .run(&request, &label)
                .and_then(Response::into_pipeline)
        }));
    }

    /// Save the Snakemake and Nextflow files for the chosen steps in a folder the person picks.
    fn export_pipeline(&mut self) {
        let (seqs, target) = self.pipeline_input();
        if seqs.is_empty() {
            self.platform.pipelines.message =
                "There is no sequence: open one in Sequences, or use the example.".into();
            return;
        }
        if dialog::IN_BROWSER {
            self.platform.pipelines.message =
                "Saving the workflow files needs a folder: not available in the browser yet."
                    .into();
            return;
        }
        let Some(dir) = dialog::FileDialog::new()
            .set_title("Folder for the workflow files")
            .pick_folder()
        else {
            return;
        };
        let settings = self.p3_win.settings.clone();
        let plan = self.plan();
        let first = seqs.first();
        let made = write_folder(
            &dir,
            &seqs,
            &plan,
            first.map(|sequence| PrimerSpec {
                sequence,
                target,
                settings: &settings,
            }),
        );
        self.platform.pipelines.message = match made {
            Ok(files) => format!(
                "Saved {} files in {}. Run `snakemake --cores 1` or `nextflow run main.nf` there.",
                files.len(),
                dir.display()
            ),
            Err(e) => e,
        };
    }

    pub(super) fn pipelines_page(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let look = design::look(self.prefs.dark);
        let listing = self.pipeline_listing();
        let (open, target) = (self.docs.len(), self.pipeline_target());
        // where a run with the steps as ticked would go, and so whose programs it would use
        let weight = self.pipeline_weight(&listing);
        let programs = self.platform.programs_at(&self.platform.place(&weight));
        let view = self.platform.workstation.view(&weight);
        let mut c = Clicks::default();
        let page = &mut self.platform.pipelines;
        let input = Input {
            listing: &listing,
            open,
            target,
            known: programs.is_some(),
            programs: programs.as_deref().unwrap_or_default(),
            view,
        };
        design::page(ui, |ui| {
            design::page_header(
                ui,
                &look,
                "Pipelines",
                "Align, build a tree and pick primers in one run, or export the same steps as a workflow",
                |_| {},
            );
            egui::Panel::left("pipelines_options")
                .exact_size(OPTIONS_WIDTH)
                .resizable(false)
                .frame(egui::Frame::new())
                .show(ui, |ui| {
                    ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            design::card(ui, &look, |ui| {
                                options_panel(ui, &look, page, &input, &mut c);
                            });
                        });
                });
            egui::CentralPanel::default()
                .frame(egui::Frame::new().inner_margin(egui::Margin {
                    left: design::PAGE_MARGIN,
                    ..egui::Margin::ZERO
                }))
                .show(ui, |ui| results(ui, &look, page, &mut c));
        });
        self.apply_pipeline_clicks(c, ui.ctx(), out);
    }

    fn apply_pipeline_clicks(&mut self, c: Clicks, ctx: &egui::Context, out: &mut Vec<Action>) {
        if c.run {
            self.platform.pipelines.run_requested = true;
        }
        if c.export {
            self.export_pipeline();
        }
        let Some(done) = &self.platform.pipelines.outcome else {
            return;
        };
        if c.copy_alignment {
            if let Some(a) = &done.aligned {
                ctx.copy_text(fasta_text(&a.rows));
                self.platform.pipelines.message = "Copied the alignment as FASTA.".into();
            }
        }
        if c.copy_newick {
            if let Some(t) = &done.tree {
                ctx.copy_text(t.newick.clone());
                self.platform.pipelines.message = "Copied the tree as Newick.".into();
            }
        }
        if c.to_alignment {
            if let Some(aligned) = done.aligned.clone() {
                let page = &mut self.platform.alignment;
                page.show_aligned(aligned, "The alignment of the Pipelines run.");
                out.push(Action::Go(Workspace::Alignment));
            }
        }
        if c.to_trees {
            let parsed = done
                .tree
                .as_ref()
                .and_then(|t| NwTree::parse_all(&t.newick).ok())
                .and_then(|trees| trees.first().map(NwTree::to_node));
            if let Some(node) = parsed {
                self.tree_win.details = done.log.clone();
                self.show_tree(node, "The tree of the Pipelines run.");
                out.push(Action::Go(Workspace::Trees));
            }
        }
    }
}

/// What the options panel shows about the run to come and the programs around.
struct Input<'a> {
    /// The name and length of each sequence the run would use.
    listing: &'a [(String, usize)],
    /// How many sequences are open in Sequences, whichever are in use.
    open: usize,
    target: Option<(usize, usize)>,
    programs: &'a [(&'static str, Option<Installed>)],
    /// The search for the programs is over, so one that is not in `programs` is missing.
    known: bool,
    /// Where the run would go: `programs` are the ones that computer has.
    view: ComputeView,
}

impl Input<'_> {
    fn found(&self, label: &str) -> Option<&Installed> {
        self.programs
            .iter()
            .find(|(l, _)| *l == label)
            .and_then(|(_, f)| f.as_ref())
    }
}

fn options_panel(
    ui: &mut Ui,
    look: &Look,
    page: &mut PipelinesPage,
    input: &Input,
    c: &mut Clicks,
) {
    design::section(ui, look, "Sequences");
    ui.radio_value(
        &mut page.source,
        Source::Open,
        format!("Those open in Sequences ({})", input.open),
    );
    ui.radio_value(
        &mut page.source,
        Source::Example,
        "The example: six made-up isolates",
    );
    for (name, length) in input.listing.iter().take(LISTED) {
        ui.label(
            RichText::new(format!("{name} · {length} bases"))
                .size(12.0)
                .color(look.dim),
        );
    }
    if input.listing.len() > LISTED {
        ui.label(
            RichText::new(format!("and {} more", input.listing.len() - LISTED))
                .size(12.0)
                .color(look.dim),
        );
    }
    if input.listing.is_empty() {
        ui.label(
            RichText::new("Nothing is open in Sequences.")
                .size(12.0)
                .color(look.dim),
        );
    }

    design::section(ui, look, "Steps");
    ui.checkbox(&mut page.align, "Align the sequences");
    ui.checkbox(&mut page.tree, "Build a tree")
        .on_hover_text("Needs the alignment, so it is made too");
    ui.checkbox(&mut page.primers, "Pick primer pairs");
    if page.primers {
        ui.label(
            RichText::new(format!(
                "On the first sequence, {}.",
                target_text(input.target)
            ))
            .size(12.0)
            .color(look.dim),
        );
    }

    design::section(ui, look, "Programs");
    compute_picker(ui, look, input.view);
    ui.horizontal_wrapped(|ui| {
        if !input.known {
            design::chip(ui, look, "Looking for installed programs…", Tone::Neutral);
            return;
        }
        for (label, instead) in ENGINES {
            program_chip(ui, look, label, input.found(label), instead);
        }
    });
    ui.checkbox(&mut page.use_installed, "Use the installed programs")
        .on_hover_text("When off, the built-in methods are used even where a program is installed");
    ui.checkbox(&mut page.bootstrap, "Ultrafast bootstrap (IQ-TREE)");

    ui.add_space(8.0);
    ui.horizontal(|ui| {
        let ready = input.known && !page.is_running();
        c.run = ui
            .add_enabled_ui(ready, |ui| design::primary_button(ui, look, "Run pipeline"))
            .inner
            .clicked();
        if page.is_running() {
            ui.spinner();
        }
    });
    // right under the button, where a failure (the workstation cannot be reached, say) is seen
    if !page.message.is_empty() {
        ui.label(RichText::new(&page.message).size(12.0).color(look.dim));
    }

    design::section(ui, look, "As a workflow");
    ui.label(
        RichText::new("The same steps as files for Snakemake and Nextflow, to run anywhere those are installed.")
            .size(12.0)
            .color(look.dim),
    );
    c.export = ui.button("Save workflow files…").clicked();
    if input.known {
        ui.horizontal_wrapped(|ui| {
            for label in RUNNERS {
                let missing = format!("{label} not found");
                program_chip(ui, look, label, input.found(label), &missing);
            }
        });
    }
}

/// A result's card: a title, and what the step made.
fn result_card(ui: &mut Ui, look: &Look, title: &str, add: impl FnOnce(&mut Ui)) {
    design::card(ui, look, |ui| {
        ui.set_min_width(ui.available_width());
        ui.label(RichText::new(title).size(15.0).strong().color(look.text));
        ui.add_space(4.0);
        add(ui);
    });
    ui.add_space(RESULT_GAP);
}

fn results(ui: &mut Ui, look: &Look, page: &PipelinesPage, c: &mut Clicks) {
    let Some(done) = &page.outcome else {
        design::card(ui, look, |ui| {
            ui.label(
                RichText::new("No run yet")
                    .size(17.0)
                    .strong()
                    .color(look.text),
            );
            ui.label(
                RichText::new(
                    "Choose the sequences and the steps on the left and press Run pipeline. The alignment, the tree and the \
                     primer pairs appear here, each with a button that opens it in its own workspace.",
                )
                .color(look.dim),
            );
        });
        return;
    };
    ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            result_card(ui, look, "What ran", |ui| {
                for line in &done.log {
                    ui.label(RichText::new(line).color(look.dim));
                }
            });
            if let Some(a) = &done.aligned {
                result_card(ui, look, "Alignment", |ui| {
                    ui.horizontal_wrapped(|ui| {
                        design::chip(
                            ui,
                            look,
                            &format!("{} sequences", a.rows.len()),
                            Tone::Neutral,
                        );
                        design::chip(ui, look, &format!("{} columns", a.columns()), Tone::Neutral);
                        design::chip(
                            ui,
                            look,
                            &format!("{} · {:.2} s", a.engine, a.seconds),
                            Tone::Accent,
                        );
                    });
                    ui.horizontal_wrapped(|ui| {
                        c.to_alignment =
                            design::primary_button(ui, look, "Open in Alignment").clicked();
                        c.copy_alignment = ui.button("Copy FASTA").clicked();
                    });
                });
            }
            if let Some(t) = &done.tree {
                result_card(ui, look, "Tree", |ui| {
                    ui.horizontal_wrapped(|ui| {
                        design::chip(ui, look, &t.engine, Tone::Accent);
                        if !t.model.is_empty() {
                            design::chip(ui, look, &format!("model {}", t.model), Tone::Neutral);
                        }
                        if let Some(l) = t.ln_likelihood {
                            design::chip(ui, look, &format!("ln L {l:.2}"), Tone::Neutral);
                        }
                        design::chip(ui, look, &format!("{:.2} s", t.seconds), Tone::Neutral);
                    });
                    // read-only: a text box on a borrowed string cannot be edited
                    let mut shown: &str = &t.newick;
                    ui.add(
                        egui::TextEdit::multiline(&mut shown)
                            .font(egui::TextStyle::Monospace)
                            .desired_rows(NEWICK_ROWS)
                            .desired_width(f32::INFINITY),
                    );
                    ui.horizontal_wrapped(|ui| {
                        c.to_trees = design::primary_button(ui, look, "Draw in Trees").clicked();
                        c.copy_newick = ui.button("Copy Newick").clicked();
                    });
                });
            }
            if let Some(p) = &done.primers {
                result_card(ui, look, "Primer pairs", |ui| {
                    ui.horizontal_wrapped(|ui| {
                        design::chip(ui, look, &format!("{} pairs", p.pairs.len()), Tone::Neutral);
                        design::chip(
                            ui,
                            look,
                            &format!("{} · {:.2} s", p.engine, p.seconds),
                            Tone::Accent,
                        );
                    });
                    if p.pairs.is_empty() {
                        ui.label(
                            RichText::new(
                                "No pair fits the primer settings; widen them in Tools → Primer3.",
                            )
                            .color(look.dim),
                        );
                    } else {
                        primer_rows(ui, look, &p.pairs);
                    }
                });
            }
        });
}

fn primer_rows(ui: &mut Ui, look: &Look, pairs: &[PrimerPair]) {
    egui::Grid::new("pipeline_primers")
        .striped(true)
        .spacing([14.0, 6.0])
        .show(ui, |ui| {
            for h in [
                "#",
                "Left primer",
                "Tm",
                "Right primer",
                "Tm",
                "Product",
                "Penalty",
            ] {
                ui.label(RichText::new(h).size(11.5).strong().color(look.dim));
            }
            ui.end_row();
            for (k, p) in pairs.iter().enumerate() {
                ui.label((k + 1).to_string());
                ui.label(RichText::new(&p.left).monospace());
                ui.label(format!("{:.1}", p.left_tm));
                ui.label(RichText::new(&p.right).monospace());
                ui.label(format!("{:.1}", p.right_tm));
                ui.label(p.product_size.to_string());
                ui.label(format!("{:.3}", p.penalty));
                ui.end_row();
            }
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn choosing_a_tree_aligns_first_and_the_bootstrap_has_its_minimum_size() {
        let mut page = PipelinesPage {
            align: false,
            ..PipelinesPage::default()
        };
        let c = page.choice();
        assert!(c.align && c.tree && c.primers && c.use_installed);
        assert_eq!(c.bootstrap, 0);
        page.tree = false;
        page.bootstrap = true;
        let c = page.choice();
        assert!(!c.align && !c.tree);
        assert_eq!(c.bootstrap, BOOTSTRAP_REPLICATES);
    }

    #[test]
    fn the_target_reads_in_one_based_bases() {
        assert_eq!(target_text(None), "anywhere along it");
        assert_eq!(target_text(Some((200, 60))), "around bases 201 to 260");
    }

    #[test]
    fn the_example_has_six_isolates_longer_than_the_primers_target() {
        let seqs = example_sequences();
        assert_eq!(seqs.len(), 6);
        let (start, length) = EXAMPLE_TARGET;
        assert!(seqs.iter().all(|(_, bases)| bases.len() >= start + length));
    }

    #[test]
    fn the_example_runs_through_all_three_steps_with_the_built_in_methods() {
        let page = PipelinesPage {
            use_installed: false,
            ..PipelinesPage::default()
        };
        let seqs = example_sequences();
        let settings = crate::seq::primer3::Settings::default();
        let done = pipeline::run(
            &seqs,
            Some(EXAMPLE_TARGET),
            &settings,
            &page.choice(),
            &Available::default(),
        )
        .expect("a run");
        assert_eq!(done.log.len(), 3);
        assert_eq!(done.aligned.as_ref().expect("an alignment").rows.len(), 6);
        let newick = &done.tree.as_ref().expect("a tree").newick;
        assert_eq!(NwTree::parse_all(newick).expect("a tree").len(), 1);
        assert!(!done.primers.as_ref().expect("primers").pairs.is_empty());
    }
}
