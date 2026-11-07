//! Tools → Species Tree from Gene Trees: the species tree that shares the most quartets with a set of gene trees, in the
//! way of ASTRAL (`seq/astral`), with a posterior probability and a length on every inner branch. It runs where the compute
//! choice says: on this computer, or as a job of the workstation.

use super::treetools::{newick_text, read_newick_file};
use crate::remote::species::{self, SpeciesJob, SpeciesResult};
use crate::remote::tasks::{Request, Response};
use crate::remote::weight::Weight;
use crate::sample::{SPECIES_GENE_TREES, SPECIES_TRUTH};
use crate::seq::astral::search::EXACT_LIMIT;
use crate::seq::astral::Options;
use crate::seq::newick::NwTree;
use crate::ui::app::App;
use crate::ui::design::{self, Look};
use crate::ui::dialog::Wanted;
use crate::ui::job::{Job, Poll};
use crate::ui::platform::route::Source;
use crate::ui::platform::trees::save_text;
use crate::ui::platform::workstation::{compute_picker, ComputeView, Place};
use crate::ui::platform::Workspace;
use eframe::egui::{self, RichText, ScrollArea, Ui};
use std::sync::Arc;
use std::time::Duration;

const GENES_BOX_MAX_HEIGHT: f32 = 90.0;
const BRANCHES_MAX_HEIGHT: f32 = 240.0;
const NEWICK_BOX_MAX_HEIGHT: f32 = 70.0;
/// The most branches the table lists: the tree has them all.
const MAX_LISTED: usize = 300;
/// The most species named for a branch in the table.
const NAMED: usize = 4;
/// The prior of the branch lengths: ASTRAL's default, and the most the window lets it be.
const LAMBDA_MAX: f64 = 10.0;
/// While a run is going the window looks again this often to move its spinner.
const WORKING_REDRAW: Duration = Duration::from_millis(200);
/// What the window leaves of the screen's height for the bars and for where it opens, and how tall it opens at least and at
/// most: the contents (a table of branches, the tree as text) need about 600 points.
const WINDOW_ROOM: f32 = 170.0;
const MIN_WINDOW_HEIGHT: f32 = 420.0;
const MAX_WINDOW_HEIGHT: f32 = 680.0;

/// A file of gene trees, kept out of the text box: a file of thousands of trees would make the box very slow to draw.
struct Loaded {
    name: String,
    text: String,
    trees: usize,
    weight: Weight,
}

pub struct SpeciesWin {
    pub open: bool,
    /// Gene trees typed or pasted, one after the other.
    pub genes: String,
    loaded: Option<Loaded>,
    /// Score the tree in `given` instead of searching.
    pub score_given: bool,
    pub given: String,
    /// Look at every group of species instead of the groups the gene trees suggest.
    pub exact: bool,
    /// Also try the ways to resolve the polytomies of the gene trees.
    pub extra: bool,
    pub lambda: f64,
    /// The species to hang the tree from when it is written; empty for the first.
    pub root: String,
    run_requested: bool,
    job: Option<Job<Result<SpeciesResult, String>>>,
    /// The run that is going was for trees that are gone: it cannot be stopped, so its answer is let go when it comes.
    discard: bool,
    /// Where the running job goes, for the message that says where it ran.
    ran_on: &'static str,
    pub result: Option<Arc<SpeciesResult>>,
    pub message: String,
}

impl Default for SpeciesWin {
    fn default() -> Self {
        let o = Options::default();
        SpeciesWin {
            open: false,
            genes: String::new(),
            loaded: None,
            score_given: false,
            given: String::new(),
            exact: o.exact,
            extra: o.extra > 0,
            lambda: o.lambda,
            root: String::new(),
            run_requested: false,
            job: None,
            discard: false,
            ran_on: "this computer",
            result: None,
            message: String::new(),
        }
    }
}

impl SpeciesWin {
    pub fn is_running(&self) -> bool {
        self.job.is_some()
    }

    /// The gene trees a run would use: the file that was opened, else what is in the box.
    fn text(&self) -> &str {
        self.loaded.as_ref().map_or(&self.genes, |l| &l.text)
    }

    /// How big the run would be, counted once for a file.
    pub fn weight(&self) -> Weight {
        match &self.loaded {
            Some(l) => l.weight,
            None => Weight::of_species(&self.genes),
        }
    }

    /// The job the settings ask for.
    fn request(&self) -> SpeciesJob {
        SpeciesJob {
            exact: self.exact,
            extra: u8::from(self.extra),
            lambda: self.lambda,
            root: self.root.trim().to_string(),
            score_tree: if self.score_given {
                self.given.clone()
            } else {
                String::new()
            },
            ..SpeciesJob::search(self.text())
        }
    }

    /// Start again with nothing. A run that was going was for trees that are gone, so its answer is not wanted; but a search
    /// cannot be stopped, so it counts as running (and keeps a new run waiting) until it ends, and then its answer is let go.
    /// (A job on the workstation goes on there and can be fetched from the Workstation page.)
    fn clear(&mut self) {
        self.genes.clear();
        self.loaded = None;
        self.result = None;
        self.message.clear();
        self.run_requested = false;
        self.discard = self.job.is_some();
    }

    fn take(&mut self, result: SpeciesResult, ran_on: &str) {
        self.message = summary(&result, ran_on);
        self.result = Some(Arc::new(result));
    }
}

/// What a run found, in a few sentences.
fn summary(r: &SpeciesResult, ran_on: &str) -> String {
    let mut notes = vec![format!(
        "{} gene trees on {} species.",
        r.genes,
        r.taxa.len()
    )];
    if r.skipped > 0 {
        notes.push(format!("{} with fewer than 4 species left out.", r.skipped));
    }
    notes.push(format!(
        "Quartet score {} ({:.1} % of the quartets the gene trees resolve).",
        r.score,
        r.normalized * 100.0
    ));
    notes.push(match r.clusters {
        Some(n) => format!(
            "Chose among {n} groups of species in {:.2} s on {ran_on}.",
            r.seconds
        ),
        None => format!(
            "Scored the tree you gave in {:.2} s on {ran_on}.",
            r.seconds
        ),
    });
    notes.join(" ")
}

/// The species on one side of a branch, as many as fit in a table cell.
fn side(taxa: &[String]) -> String {
    match taxa.len() {
        0..=NAMED => taxa.join(", "),
        n => format!("{}, … ({} more)", taxa[..NAMED].join(", "), n - NAMED),
    }
}

/// What to say when the table lists fewer branches than the tree has.
fn cut_note(branches: usize) -> Option<String> {
    (branches > MAX_LISTED).then(|| {
        format!("The first {MAX_LISTED} of {branches} inner branches; the tree has them all.")
    })
}

fn length_text(length: f64) -> String {
    if length.is_nan() {
        "–".into()
    } else if length.is_infinite() {
        "∞".into()
    } else {
        format!("{length:.3}")
    }
}

#[derive(Default)]
struct Clicks {
    open: bool,
    example: bool,
    clear: bool,
    run: bool,
    to_trees: bool,
    copy: bool,
    save: bool,
}

impl App {
    /// Show what the workstation made, and open the window to show it.
    pub fn show_species_result(&mut self, result: SpeciesResult) {
        self.species_win.open = true;
        self.species_win.take(result, "the workstation");
    }

    /// The example: genes simulated from a known species tree, which the window can also be asked to score.
    pub fn load_species_example(&mut self) {
        let w = &mut self.species_win;
        w.clear();
        w.genes = SPECIES_GENE_TREES.trim().to_string();
        w.given = SPECIES_TRUTH.trim().to_string();
        w.message =
            "Loaded 40 gene trees on 8 species, simulated under the multispecies coalescent. \
                     Press Find the species tree."
                .into();
    }

    fn open_species_file(&mut self) {
        if let Some(read) = read_newick_file(Wanted::SpeciesTrees) {
            self.use_species_file(read);
        }
    }

    /// Keep a file of gene trees that was read, or say why it could not be read.
    fn use_species_file(&mut self, read: Result<(String, String), String>) {
        let w = &mut self.species_win;
        match read {
            Ok((name, text)) => {
                let trees = text.matches(';').count();
                w.clear();
                w.message =
                    format!("Loaded {name}: {trees} gene trees. Press Find the species tree.");
                w.loaded = Some(Loaded {
                    weight: Weight::of_species(&text),
                    name,
                    text,
                    trees,
                });
            }
            Err(e) => w.message = e,
        }
    }

    /// A file of gene trees for the window, from a source the browser brought.
    pub(in crate::ui) fn load_species_file(&mut self, file: Source) {
        self.use_species_file(newick_text(file));
    }

    /// Ask for the run; the window starts it on its next draw.
    pub fn run_species(&mut self) {
        self.species_win.run_requested = true;
    }

    fn start_species(&mut self, ctx: &egui::Context) {
        let job = self.species_win.request();
        let place = self.platform.place(&self.species_win.weight());
        let w = &mut self.species_win;
        w.run_requested = false;
        w.result = None;
        if job.genes.trim().is_empty() {
            w.message = "Paste the gene trees or open a file of them first.".into();
            return;
        }
        if w.score_given && job.score_tree.trim().is_empty() {
            w.message = "Paste the species tree to score, or untick the box in Settings.".into();
            return;
        }
        let Place::There(remote) = place else {
            (w.ran_on, w.message) = ("this computer", "Working…".into());
            w.job = Some(Job::spawn(ctx, move || species::run(&job)));
            return;
        };
        (w.ran_on, w.message) = ("the workstation", "Working on the workstation…".into());
        let label = format!(
            "Species tree from {} gene trees",
            job.genes.matches(';').count()
        );
        let request = Request::Species(job);
        w.job = Some(Job::spawn(ctx, move || {
            remote
                .run(&request, &label)
                .and_then(Response::into_species)
        }));
    }

    /// Start the run that was asked for and take its answer when it comes, whether the window is showing or not.
    pub fn tick_species(&mut self, ctx: &egui::Context) {
        if self.species_win.run_requested && self.species_win.job.is_none() {
            self.start_species(ctx);
        }
        let w = &mut self.species_win;
        let Some(job) = &w.job else {
            return;
        };
        match job.poll() {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(outcome) => {
                w.job = None;
                let discarded = std::mem::take(&mut w.discard);
                match outcome {
                    // the trees it was for are gone: nothing of it is kept
                    _ if discarded => {}
                    Ok(result) => w.take(result, w.ran_on),
                    Err(why) => w.message = why,
                }
            }
            Poll::Died => {
                w.job = None;
                if !std::mem::take(&mut w.discard) {
                    w.message = "The run stopped without an answer.".into();
                }
            }
        }
    }

    /// Draw the species tree on the Trees page.
    fn show_species_in_trees(&mut self) {
        let Some(result) = self.species_win.result.clone() else {
            return;
        };
        // a prior of 0 makes some lengths infinite or not a number, which ASTRAL writes out and the reader here does not take:
        // such a branch is drawn without a length
        let drawable = result.newick.replace(":Infinity", "").replace(":NaN", "");
        match NwTree::parse_all(&drawable) {
            Ok(mut trees) if !trees.is_empty() => {
                let note = format!(
                    "The species tree from {} gene trees; branches carry the local posterior probability.",
                    result.genes
                );
                self.show_tree(trees.remove(0).to_node(), &note);
                self.workspace = Workspace::Trees;
            }
            Ok(_) => self.species_win.message = "The result holds no tree.".into(),
            Err(e) => self.species_win.message = format!("The tree cannot be drawn: {e}"),
        }
    }

    fn apply_species_clicks(&mut self, c: Clicks, ctx: &egui::Context) {
        if c.open {
            self.open_species_file();
        }
        if c.example {
            self.load_species_example();
        }
        if c.clear {
            self.species_win.clear();
        }
        if c.run {
            self.run_species();
        }
        let Some(result) = self.species_win.result.clone() else {
            return;
        };
        if c.copy {
            ctx.copy_text(result.newick.clone());
            self.species_win.message = "Copied the species tree as Newick.".into();
        }
        if c.save {
            let said = save_text("Save the species tree", "species-tree.nwk", &result.newick);
            if !said.is_empty() {
                self.species_win.message = said;
            }
        }
        if c.to_trees {
            self.show_species_in_trees();
        }
    }

    pub fn species_window(&mut self, ctx: &egui::Context) {
        self.tick_species(ctx);
        if !self.species_win.open {
            return;
        }
        let (mut open, mut c) = (true, Clicks::default());
        let view = self.platform.workstation.view(&self.species_win.weight());
        let look = design::look(self.prefs.dark);
        // a scrolling window keeps the height it first has, so it opens tall enough for a result and no taller than the screen
        let height =
            (ctx.content_rect().height() - WINDOW_ROOM).clamp(MIN_WINDOW_HEIGHT, MAX_WINDOW_HEIGHT);
        design::window("Species Tree from Gene Trees")
            .open(&mut open)
            .default_width(600.0)
            .default_height(height)
            .vscroll(true)
            .show(ctx, |ui| {
                draw(ui, &look, view, &mut self.species_win, &mut c)
            });
        self.apply_species_clicks(c, ctx);
        self.species_win.open = open;
    }
}

fn draw(ui: &mut Ui, look: &Look, view: ComputeView, w: &mut SpeciesWin, c: &mut Clicks) {
    ui.label(
        "The species tree that agrees with the most quartets of your gene trees (the method of ASTRAL), which stays \
         right when genes disagree through incomplete lineage sorting. One unrooted tree per gene, one after the \
         other; a gene may lack species.",
    );
    inputs(ui, w, c);
    ui.collapsing("Settings", |ui| settings(ui, w));
    compute_picker(ui, look, view);
    ui.horizontal(|ui| {
        let label = if w.score_given {
            "Score the tree"
        } else {
            "Find the species tree"
        };
        c.run = ui
            .add_enabled(
                !w.is_running(),
                egui::Button::new(RichText::new(label).strong()),
            )
            .clicked();
        if w.is_running() {
            ui.spinner();
        }
    });
    if !w.message.is_empty() {
        ui.label(&w.message);
    }
    if let Some(result) = w.result.clone() {
        results(ui, &result, c);
    }
}

fn inputs(ui: &mut Ui, w: &mut SpeciesWin, c: &mut Clicks) {
    ui.horizontal(|ui| {
        c.open = ui.button("Open gene trees…").clicked();
        c.example = ui
            .button("Example")
            .on_hover_text("40 gene trees on 8 species, simulated from a known species tree")
            .clicked();
        c.clear = ui.button("Clear").clicked();
    });
    match &w.loaded {
        Some(l) => {
            ui.label(RichText::new(format!("{}: {} gene trees", l.name, l.trees)).strong());
        }
        None => {
            ScrollArea::vertical()
                .max_height(GENES_BOX_MAX_HEIGHT)
                .show(ui, |ui| {
                    ui.add(
                        egui::TextEdit::multiline(&mut w.genes)
                            .font(egui::TextStyle::Monospace)
                            .hint_text("Paste the gene trees: ((A,B),(C,D)); on a line each")
                            .desired_rows(3)
                            .desired_width(f32::INFINITY),
                    );
                });
        }
    }
}

fn settings(ui: &mut Ui, w: &mut SpeciesWin) {
    ui.checkbox(
        &mut w.score_given,
        "Score a species tree that I give, instead of searching",
    );
    if w.score_given {
        ui.add(
            egui::TextEdit::multiline(&mut w.given)
                .font(egui::TextStyle::Monospace)
                .hint_text("A binary species tree in Newick")
                .desired_rows(2)
                .desired_width(f32::INFINITY),
        );
    } else {
        ui.checkbox(
            &mut w.extra,
            "Also try the ways to resolve the gene trees' polytomies",
        );
        ui.checkbox(
            &mut w.exact,
            format!("Look at every group of species (exact; at most {EXACT_LIMIT} species)"),
        );
    }
    egui::Grid::new("species_settings")
        .num_columns(2)
        .show(ui, |ui| {
            ui.label("Branch length prior");
            ui.add(
                egui::DragValue::new(&mut w.lambda)
                    .range(0.0..=LAMBDA_MAX)
                    .speed(0.01)
                    .fixed_decimals(2),
            )
            .on_hover_text(
                "The Yule prior of the branch lengths, 0.5 as in ASTRAL; 0 gives maximum likelihood lengths",
            );
            ui.end_row();
            ui.label("Hang the tree from");
            ui.add(
                egui::TextEdit::singleline(&mut w.root)
                    .hint_text("the first species")
                    .desired_width(160.0),
            );
            ui.end_row();
        });
}

fn results(ui: &mut Ui, r: &SpeciesResult, c: &mut Clicks) {
    ui.separator();
    ui.horizontal(|ui| {
        c.to_trees = ui.button("Open in Trees").clicked();
        c.copy = ui.button("Copy Newick").clicked();
        c.save = ui.button("Save…").clicked();
    });
    ScrollArea::vertical()
        .id_salt("species_branches_scroll")
        .max_height(BRANCHES_MAX_HEIGHT)
        .show(ui, |ui| branch_table(ui, r));
    if let Some(note) = cut_note(r.branch_count) {
        ui.label(note);
    }
    ui.label(
        RichText::new(
            "Posterior: the local posterior probability that the branch is right. Quartet support: the share \
             of the gene trees' quartets around the branch that agree with it. Lengths are in coalescent units.",
        )
        .small(),
    );
    let mut newick = r.newick.as_str();
    ScrollArea::vertical()
        .id_salt("species_newick_scroll")
        .max_height(NEWICK_BOX_MAX_HEIGHT)
        .show(ui, |ui| {
            ui.add(
                egui::TextEdit::multiline(&mut newick)
                    .font(egui::TextStyle::Monospace)
                    .desired_rows(2)
                    .desired_width(f32::INFINITY),
            );
        });
}

fn branch_table(ui: &mut Ui, r: &SpeciesResult) {
    egui::Grid::new("species_branches")
        .striped(true)
        .show(ui, |ui| {
            for h in [
                "Species on one side",
                "Posterior",
                "Quartet support",
                "Length",
            ] {
                ui.label(RichText::new(h).strong());
            }
            ui.end_row();
            for b in r.branches.iter().take(MAX_LISTED) {
                let names = ui.label(side(&b.taxa));
                if b.taxa.len() > NAMED {
                    names.on_hover_text(b.taxa.join(", "));
                }
                ui.label(
                    RichText::new(format!("{:.2}", b.pp[0]))
                        .monospace()
                        .strong(),
                );
                ui.label(RichText::new(format!("{:.0} %", b.q[0] * 100.0)).monospace());
                ui.label(RichText::new(length_text(b.length)).monospace());
                ui.end_row();
            }
        });
}

#[cfg(test)]
mod tests;
