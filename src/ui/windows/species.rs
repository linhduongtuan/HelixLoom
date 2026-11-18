//! Tools → Species Tree from Gene Trees: the species tree that shares the most quartets with a set of gene trees, in the
//! way of ASTRAL (`seq/astral`), with a posterior probability and a length on every inner branch. It runs where the compute
//! choice says: on this computer, or as a job of the workstation.

use super::treetools::{newick_text, read_newick_file};
use crate::provenance::{Asked, Input, Provenance};
use crate::remote::species::{self, as_asked, SpeciesJob, SpeciesResult};
use crate::remote::tasks::{Request, Response};
use crate::remote::weight::Weight;
use crate::sample::{SPECIES_GENE_TREES, SPECIES_TRUTH};
use crate::seq::astral::bootstrap::{read_list, MAX_REPLICATES};
use crate::seq::astral::search::{EXACT_LIMIT, WEIGHTED_WORDS};
use crate::seq::astral::{polytomy, BranchOut};
use crate::seq::astral::{Options, Resample, Weights};
use crate::seq::newick::NwTree;
use crate::ui::app::App;
use crate::ui::design::{self, Look};
use crate::ui::dialog::{self, Wanted};
use crate::ui::fonts::FontWeight;
use crate::ui::job::{Job, Poll};
use crate::ui::platform::route::Source;
use crate::ui::platform::trees::save_text;
use crate::ui::platform::workstation::{compute_picker, ComputeView, Place};
use crate::ui::platform::Workspace;
use crate::ui::provenance as made_ui;
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
    /// How the quartets of the gene trees are weighted.
    pub weights: Weights,
    /// The gene trees are gene family trees (ASTRAL-Pro): a species may occur in a tree several times.
    pub families: bool,
    /// Bootstrap replicates (none when 0), what they resample, and the seed that draws them.
    pub replicates: usize,
    pub resample: Resample,
    pub seed: u64,
    /// Test every branch for being a polytomy (Sayyari and Mirarab 2018).
    pub polytomy: bool,
    /// The bootstrap trees of the genes, from a list of files: its name, and the text of each gene's file.
    site_trees: Option<(String, Vec<String>)>,
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
    /// What the run that is going was asked, with its settings; and, for the result that is shown, the record of what made it
    /// (none for a result that was fetched from the workstation's list).
    asked: Option<(Asked, Kept)>,
    pub made: Option<(Provenance, Kept)>,
    pub message: String,
}

/// The settings of a run, to put back when it is run again.
#[derive(Clone, Debug, PartialEq)]
pub struct Kept {
    score_given: bool,
    weights: Weights,
    families: bool,
    replicates: usize,
    resample: Resample,
    seed: u64,
    polytomy: bool,
    exact: bool,
    extra: bool,
    lambda: f64,
    root: String,
}

/// The inputs and the settings of a species-tree job, as a record names them. Settings that the run does not use (those of a
/// bootstrap that is not asked for) are left out.
pub fn described(job: &SpeciesJob) -> (Vec<Input>, Vec<(String, String)>) {
    let mut inputs = vec![Input::of("gene trees", job.genes.as_bytes())];
    if !job.score_tree.trim().is_empty() {
        inputs.push(Input::of(
            "species tree to score",
            job.score_tree.as_bytes(),
        ));
    }
    if !job.site_trees.is_empty() {
        inputs.push(Input::of(
            "bootstrap trees of the genes",
            job.site_trees.concat().as_bytes(),
        ));
    }
    let yes_no = |b: bool| if b { "yes" } else { "no" }.to_string();
    let mut settings = vec![
        ("weights".to_string(), job.weights.name().to_string()),
        ("gene families".to_string(), yes_no(job.families)),
        ("polytomy test".to_string(), yes_no(job.polytomy)),
        ("exact search".to_string(), yes_no(job.exact)),
        ("polytomy resolutions".to_string(), yes_no(job.extra > 0)),
        ("branch length prior".to_string(), job.lambda.to_string()),
        (
            "hung from".to_string(),
            if job.root.is_empty() {
                "the first species".to_string()
            } else {
                job.root.clone()
            },
        ),
        (
            "bootstrap replicates".to_string(),
            job.replicates.to_string(),
        ),
    ];
    if job.replicates > 0 {
        settings.push(("resampling".to_string(), job.resample.name().to_string()));
        settings.push(("seed".to_string(), job.seed.to_string()));
    }
    (inputs, settings)
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
            weights: o.weights,
            families: o.families,
            replicates: o.bootstrap.replicates,
            resample: o.bootstrap.resample,
            seed: o.bootstrap.seed,
            polytomy: o.polytomy,
            site_trees: None,
            exact: o.exact,
            extra: o.extra > 0,
            lambda: o.lambda,
            root: String::new(),
            run_requested: false,
            job: None,
            discard: false,
            ran_on: "this computer",
            result: None,
            asked: None,
            made: None,
            message: String::new(),
        }
    }
}

impl SpeciesWin {
    /// The settings as they are, to keep with a run.
    fn kept(&self) -> Kept {
        Kept {
            score_given: self.score_given,
            weights: self.weights,
            families: self.families,
            replicates: self.replicates,
            resample: self.resample,
            seed: self.seed,
            polytomy: self.polytomy,
            exact: self.exact,
            extra: self.extra,
            lambda: self.lambda,
            root: self.root.clone(),
        }
    }

    /// Put the settings of a run back.
    fn restore(&mut self, k: &Kept) {
        self.score_given = k.score_given;
        self.weights = k.weights;
        self.families = k.families;
        self.replicates = k.replicates;
        self.resample = k.resample;
        self.seed = k.seed;
        self.polytomy = k.polytomy;
        self.exact = k.exact;
        self.extra = k.extra;
        self.lambda = k.lambda;
        self.root = k.root.clone();
    }

    /// The inputs and the settings of a run as they are now.
    pub fn as_it_is(&self) -> (Vec<Input>, Vec<(String, String)>) {
        described(&self.request())
    }

    pub fn is_running(&self) -> bool {
        self.job.is_some()
    }

    /// The gene trees a run would use: the file that was opened, else what is in the box.
    fn text(&self) -> &str {
        self.loaded.as_ref().map_or(&self.genes, |l| &l.text)
    }

    /// How big the run would be, counted once for a file.
    pub fn weight(&self) -> Weight {
        let plain = match &self.loaded {
            Some(l) => l.weight,
            None => Weight::of_species(&self.genes),
        };
        // a weighted count is dearer, and every replicate of a bootstrap is a search of its own: the same trees are a bigger job
        let times = match self.weights {
            Weights::None => 1,
            _ => WEIGHTED_WORDS as usize,
        } * (1 + self.replicates);
        Weight {
            leaves: plain.leaves.saturating_mul(times),
            ..plain
        }
    }

    /// The job the settings ask for.
    fn request(&self) -> SpeciesJob {
        SpeciesJob {
            weights: if self.families {
                Weights::None
            } else {
                self.weights
            },
            families: self.families,
            replicates: self.replicates,
            seed: self.seed,
            resample: self.resample,
            // gene families have no table of quartet counts to test
            polytomy: self.polytomy && !self.families,
            site_trees: match (&self.site_trees, self.resample.uses_site_trees()) {
                (Some((_, texts)), true) if self.replicates > 0 => texts.clone(),
                _ => Vec::new(),
            },
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
        self.made = None;
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
    notes.push(match r.weighted_score {
        _ if r.families => format!(
            "Quartet score {} ({:.1} % of the quartet classes the gene family trees have).",
            r.score,
            r.normalized * 100.0
        ),
        Some(w) => format!(
            "Quartet score {} (weighted by {}: {w:.2}, {:.1} % of the quartets the gene trees resolve).",
            r.score,
            weights_phrase(r.weights),
            r.normalized * 100.0
        ),
        None => format!(
            "Quartet score {} ({:.1} % of the quartets the gene trees resolve).",
            r.score,
            r.normalized * 100.0
        ),
    });
    if let Some(b) = r.bootstrap {
        notes.push(format!(
            "Bootstrap support from {} of {} replicates resampling {}{}.",
            b.done,
            b.asked,
            match b.resample {
                Resample::Genes => "the genes",
                Resample::Sites => "the sites",
                Resample::Both => "the genes and their sites",
            },
            if b.failed > 0 {
                format!(" ({} gave no tree)", b.failed)
            } else {
                String::new()
            }
        ));
    }
    if r.polytomy {
        notes.push(polytomy_phrase(&r.branches));
    }
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

/// What the test for a polytomy found, in a sentence: the branches tested, and those that it cannot tell from a polytomy.
fn polytomy_phrase(branches: &[BranchOut]) -> String {
    let tested: Vec<_> = branches.iter().filter_map(|b| b.polytomy).collect();
    let kept = tested.iter().filter(|t| !t.rejected()).count();
    format!(
        "Polytomy test: {} of {} branches tested ({} or more gene trees needed), {kept} not rejected as polytomies at p of {} or more.",
        tested.len(),
        branches.len(),
        polytomy::MIN_GENES,
        polytomy::ALPHA
    )
}

/// What a way of weighting the quartets is called in a sentence.
fn weights_phrase(w: Weights) -> &'static str {
    match w {
        Weights::None => "nothing",
        Weights::Support => "the support of the gene trees' branches",
        Weights::Length => "the lengths of the gene trees' branches",
        Weights::Hybrid => "both the support and the lengths of the gene trees' branches",
    }
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
    made: made_ui::Asks,
    boot_list: bool,
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
        // fetched from the workstation's list: what it was made from is not known here
        self.species_win.made = None;
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
        (w.made, w.asked) = (None, None);
        if job.genes.trim().is_empty() {
            w.message = "Paste the gene trees or open a file of them first.".into();
            return;
        }
        if w.score_given && job.score_tree.trim().is_empty() {
            w.message = "Paste the species tree to score, or untick the box in Settings.".into();
            return;
        }
        let (inputs, parameters) = described(&job);
        let ran_on = if matches!(place, Place::There(_)) {
            "the workstation"
        } else {
            "this computer"
        };
        w.asked = Some((
            Asked {
                analysis: "Species tree",
                inputs,
                parameters,
                ran_on,
            },
            w.kept(),
        ));
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
        let asked = job.asked();
        let request = Request::Species(job);
        w.job = Some(Job::spawn_async(ctx, async move {
            remote
                .run_async(&request, &label)
                .await
                .and_then(Response::into_species)
                .and_then(|result| as_asked(asked, result))
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
                let asked = w.asked.take();
                match outcome {
                    // the trees it was for are gone: nothing of it is kept
                    _ if discarded => {}
                    Ok(result) => {
                        // the search is this program's own, so it is the only program named
                        w.made = asked.map(|(a, k)| (a.made(Vec::new(), result.seconds), k));
                        w.take(result, w.ran_on);
                    }
                    Err(why) => w.message = why,
                }
            }
            Poll::Died => {
                w.job = None;
                w.asked = None;
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

    /// Ask for the list of the files of bootstrap trees, and read the files it names.
    fn load_bootstrap_list(&mut self) {
        let ask = dialog::FileDialog::new().set_title(
            "A list of the files of bootstrap trees (one path to a line, one file for each gene)",
        );
        let Some(path) = dialog::pick_one(Wanted::Anything, ask) else {
            return;
        };
        let name = path.file_name().map_or_else(
            || path.display().to_string(),
            |n| n.to_string_lossy().into_owned(),
        );
        let w = &mut self.species_win;
        match read_list(&path) {
            Ok(texts) => {
                w.message = format!(
                    "Read the bootstrap trees of {} genes from {name}.",
                    texts.len()
                );
                w.site_trees = Some((name, texts));
            }
            Err(why) => w.message = why,
        }
    }

    fn apply_species_clicks(&mut self, c: Clicks, ctx: &egui::Context) {
        if c.boot_list {
            self.load_bootstrap_list();
        }
        if c.open {
            self.open_species_file();
        }
        if c.example {
            self.load_species_example();
        }
        if c.clear {
            self.species_win.clear();
        }
        if c.made.again {
            // with the settings the result that is shown was made with
            if let Some((_, kept)) = self.species_win.made.clone() {
                self.species_win.restore(&kept);
            }
        }
        self.checks_win.open |= c.made.checks;
        if c.run || c.made.again {
            self.run_species();
        }
        if let (true, Some((made, _))) = (c.made.copy, &self.species_win.made) {
            ctx.copy_text(made.to_text());
            self.species_win.message = "Copied the record of how the species tree was made.".into();
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
        let look = self.prefs.look();
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
    ui.collapsing("Settings", |ui| settings(ui, w, c));
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
                egui::Button::new(RichText::new(label).semibold()),
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
        if let Some((made, _)) = &w.made {
            c.made = made_ui::section(ui, look, made, || w.as_it_is(), !w.is_running());
        }
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
            ui.label(RichText::new(format!("{}: {} gene trees", l.name, l.trees)).semibold());
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

fn weights_label(w: Weights) -> &'static str {
    match w {
        Weights::None => "None: every quartet counts one",
        Weights::Support => "By support",
        Weights::Length => "By length",
        Weights::Hybrid => "By support and length",
    }
}

fn resample_label(r: Resample) -> &'static str {
    match r {
        Resample::Genes => "The genes",
        Resample::Sites => "The sites (a bootstrap tree of each gene)",
        Resample::Both => "The genes and their sites",
    }
}

fn settings(ui: &mut Ui, w: &mut SpeciesWin, c: &mut Clicks) {
    ui.checkbox(
        &mut w.score_given,
        "Score a species tree that I give, instead of searching",
    );
    ui.checkbox(
        &mut w.families,
        "The gene trees are gene families: a species may be in a tree several times",
    )
    .on_hover_text(
        "ASTRAL-Pro (Zhang et al. 2020). The leaf names are the species names, repeated for each copy. A tree with a root of \
         two branches is rooted as written; one with three is rooted where duplications and losses cost least. Binary trees \
         only. The species tree has no branch support or length, and the weights are not used.",
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
            ui.label("Weights");
            ui.add_enabled_ui(!w.families, |ui| {
                egui::ComboBox::from_id_salt("species_weights")
                    .selected_text(weights_label(w.weights))
                    .show_ui(ui, |ui| {
                        for choice in Weights::ALL {
                            ui.selectable_value(&mut w.weights, choice, weights_label(choice));
                        }
                    })
                    .response
            })
            .inner
                .on_hover_text(
                    "Weighted ASTRAL (Zhang and Mirarab 2022): a gene tree sure of its quartets counts more than one that is not. \
                     By support: the labels of the gene trees' inner branches are read as supports (0 to 1, or percentages). \
                     By length: the branch lengths are used. Branches that the gene trees give no support or length for are \
                     taken as fully supported, or of length 0.",
                );
            ui.end_row();
            ui.label("Polytomy test");
            ui.add_enabled(
                !w.families,
                egui::Checkbox::new(&mut w.polytomy, "Test each branch"),
            )
            .on_hover_text(
                "Sayyari and Mirarab 2018: a chi-squared test of whether the gene trees show the three ways to resolve a branch \
                 equally often, which they would if the branch had no length. A small p-value says the branch is real; a large \
                 one does not say it is not. Only for a branch that at least 10 gene trees resolve a quartet around. Not with the \
                 gene families.",
            );
            ui.end_row();
            ui.label("Bootstrap replicates");
            ui.add(
                egui::DragValue::new(&mut w.replicates)
                    .range(0..=MAX_REPLICATES)
                    .speed(1),
            )
            .on_hover_text(
                "0: none. Otherwise the search is run again on this many resamplings, and the number on each branch is the \
                 share of those trees (in percent) that have it. Not with the gene families; with the weights only when the genes are resampled.",
            );
            ui.end_row();
            if w.replicates > 0 {
                ui.label("Resample");
                egui::ComboBox::from_id_salt("species_resample")
                    .selected_text(resample_label(w.resample))
                    .show_ui(ui, |ui| {
                        for choice in Resample::ALL {
                            ui.selectable_value(&mut w.resample, choice, resample_label(choice));
                        }
                    });
                ui.end_row();
                ui.label("Seed");
                ui.add(egui::DragValue::new(&mut w.seed))
                    .on_hover_text("The same seed gives the same support");
                ui.end_row();
            }
        });
    // the sites are resampled from the bootstrap trees of the genes, from files that a list names (a browser has no files)
    if w.replicates > 0 && w.resample.uses_site_trees() && !dialog::IN_BROWSER {
        ui.horizontal_wrapped(|ui| {
            c.boot_list = ui
                .button("Bootstrap trees of the genes…")
                .on_hover_text(
                    "A text file that lists, one path to a line, the files of bootstrap trees of the genes, in the order of the \
                     gene trees (a relative path is taken from the folder of the list)",
                )
                .clicked();
            match &w.site_trees {
                Some((name, texts)) => {
                    ui.label(format!("{name}: {} genes", texts.len()));
                }
                None => {
                    ui.label("none chosen yet");
                }
            }
        });
    }
}

fn results(ui: &mut Ui, r: &SpeciesResult, c: &mut Clicks) {
    ui.separator();
    ui.horizontal(|ui| {
        c.to_trees = ui.button("Open in Trees").clicked();
        c.copy = ui.button("Copy Newick").clicked();
        c.save = ui.button("Save…").clicked();
    });
    if r.families {
        ui.label(
            RichText::new(
                "Gene families: the species tree has no branch support or length (ASTRAL-Pro's are not built), \
                 so there is no table of branches.",
            )
            .small(),
        );
    } else {
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
        if r.polytomy {
            ui.label(
                RichText::new(
                    "Polytomy p: the p-value of the test that the branch has no length (under 0.05 in bold: the branch is \
                     real). A dash: fewer than 10 gene trees resolve a quartet around the branch. The branches are \
                     tested one by one, with no correction for how many there are.",
                )
                .small(),
            );
        }
        if r.weights != Weights::None {
            ui.label(
                RichText::new(
                    "These supports and lengths are ASTRAL's, counted on the quartets without their weights: the \
                     weighted support of wastral is not built.",
                )
                .small(),
            );
        }
    }
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
    let boot = r.branches.iter().any(|b| b.bootstrap.is_some());
    let test = r.polytomy;
    egui::Grid::new("species_branches")
        .striped(true)
        .show(ui, |ui| {
            for h in [
                "Species on one side",
                "Posterior",
                "Quartet support",
                "Length",
            ]
            .into_iter()
            .chain(boot.then_some("Bootstrap"))
            .chain(test.then_some("Polytomy p"))
            {
                ui.label(RichText::new(h).semibold());
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
                        .semibold(),
                );
                ui.label(RichText::new(format!("{:.0} %", b.q[0] * 100.0)).monospace());
                ui.label(RichText::new(length_text(b.length)).monospace());
                if boot {
                    let text = b.bootstrap.map_or("–".to_string(), |s| format!("{s:.0} %"));
                    ui.label(RichText::new(text).monospace().semibold());
                }
                if test {
                    match b.polytomy {
                        Some(t) => {
                            let p = RichText::new(polytomy::p_text(t.p)).monospace();
                            ui.label(if t.rejected() { p.semibold() } else { p })
                                .on_hover_text(format!(
                                    "chi-squared {:.2} from {} gene trees",
                                    t.chi2, t.genes
                                ));
                        }
                        None => {
                            ui.label(RichText::new("–").monospace());
                        }
                    }
                }
                ui.end_row();
            }
        });
}

#[cfg(test)]
mod tests;
