//! Alignment: gather sequences, align them with the built-in aligner or MAFFT, and read the result in a
//! colour-coded viewer with a consensus, conservation and identity.

use super::alnview::{self, Colouring};
use super::route::Source;
use super::workstation::{run_button, ComputeView, Place};
use super::Workspace;
use crate::pipeline::steps::{align, fasta_text, read_fasta, Aligned, Named};
use crate::pipeline::tools::{self, Installed};
use crate::provenance::{self, Asked, Provenance};
use crate::remote::tasks::{Request, Response};
use crate::remote::weight::Weight;
use crate::seq::alnstats::{self, Column, Summary};
use crate::seq::conservation::{self, Block, Method, Settings};
use crate::seq::io;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::dialog::{self, Wanted};
use crate::ui::job::{Job, Poll};
use crate::ui::provenance as made_ui;
use crate::ui::windows::tree::open_sequences;
use eframe::egui::{self, RichText, ScrollArea, Ui};
use std::time::Duration;

const LIST_HEIGHT: f32 = 230.0;
const NAME_CHARS: usize = 26;
const WORKING_REDRAW: Duration = Duration::from_millis(150);

#[derive(Clone)]
pub struct Input {
    pub name: String,
    pub bases: Vec<u8>,
    pub keep: bool,
}

#[derive(Clone, Copy, PartialEq, Eq, Default)]
pub enum Engine {
    #[default]
    BuiltIn,
    Mafft,
}

/// The rows of an alignment, the score of each column as the bars show it, and which score that is.
pub(in crate::ui) type ConservationView<'a> = (&'a [(String, Vec<u8>)], &'a [Option<f32>], Score);

/// The score that the bars above the columns show.
#[derive(Clone, Copy, PartialEq, Eq, Default, Debug)]
pub enum Score {
    /// The share of the rows that carry the commonest symbol.
    #[default]
    Identity,
    Entropy,
    /// The Jensen-Shannon divergence (Capra and Singh 2007), with sequence weights, a window and a gap penalty.
    Jsd,
}

impl Score {
    pub const ALL: [Score; 3] = [Score::Identity, Score::Entropy, Score::Jsd];

    pub(in crate::ui) fn label(self) -> &'static str {
        match self {
            Score::Identity => "Identity",
            Score::Entropy => "Shannon entropy",
            Score::Jsd => "Jensen-Shannon divergence",
        }
    }

    fn explain(self) -> &'static str {
        match self {
            Score::Identity => "The share of the sequences that have the commonest residue in the column. Simple, and steady with few sequences.",
            Score::Entropy => "One minus the entropy of the column's residues over its largest value: 1 when every sequence has the same residue. Sequences are weighted so that near copies count less, and gaps lower the score.",
            Score::Jsd => "How far the column's residues are from what chance gives (the BLOSUM62 background for proteins): the score of Capra and Singh (2007). In their tests on known catalytic and ligand-binding residues it ranked those residues better, on average over many proteins, than the other conservation scores they compared, Shannon entropy included. Sequences are weighted, gaps lower the score, and 3 columns on each side are mixed in. It is made for alignments of many sequences, and is a guide to where to look, not a verdict. For DNA the same method is used with an even background: that is this program's own extension.",
        }
    }
}

/// Fewer sequences than this and a conservation score says little (a note says so; the paper's tests used alignments of
/// many more, so this is a rule of thumb of this program, not a figure of the paper).
const FEW_SEQUENCES: usize = 10;
/// A conserved block is at least this many columns in a row.
const MIN_BLOCK: usize = 4;
/// So many blocks are listed under the chips; the rest are counted.
const MAX_BLOCKS_SHOWN: usize = 8;

/// The bars of the viewer for one score, and the conserved blocks that follow from them.
pub struct Bars {
    score: Score,
    values: Vec<Option<f32>>,
    /// The score from which a column counts as conserved: all of the rows alike for identity, else the top quarter.
    highlight: f32,
    blocks: Vec<Block>,
}

impl Bars {
    fn new(rows: &[Vec<u8>], columns: &[Column], protein: bool, score: Score) -> Bars {
        let (values, highlight) = match score {
            // every row alike is exactly 1 (the share is a count over itself)
            Score::Identity => (columns.iter().map(|c| Some(c.conservation)).collect(), 1.0),
            Score::Entropy | Score::Jsd => {
                let method = if score == Score::Jsd {
                    Method::Jsd
                } else {
                    Method::Entropy
                };
                let values = conservation::scores(
                    rows,
                    protein,
                    &Settings {
                        method,
                        ..Settings::default()
                    },
                );
                let highlight = conservation::top_quarter(&values).unwrap_or(1.0);
                (values, highlight)
            }
        };
        let blocks = conservation::blocks(&values, highlight, MIN_BLOCK);
        Bars {
            score,
            values,
            highlight,
            blocks,
        }
    }
}

/// What a reader wants to know about the alignment, worked out once when it arrives.
pub struct Stats {
    columns: Vec<Column>,
    summary: Summary,
    protein: bool,
    bars: Bars,
}

pub struct AlignmentPage {
    pub inputs: Vec<Input>,
    /// Counted up whenever a sequence is put on the list or taken off it: the session keeps the list, and tells by this that it is
    /// another list (two lists can have the same names and lengths).
    listed: u64,
    paste: String,
    engine: Engine,
    colouring: Colouring,
    score: Score,
    pub result: Option<Aligned>,
    stats: Option<Stats>,
    /// Counted up whenever the scores that the bars show change (a new alignment, a new choice): what was made from the
    /// scores earlier (the colours of a structure) is out of date when this has moved on.
    revision: u64,
    pub message: String,
    /// What the run that is going was asked, with the engine chosen for it; and, for the alignment that is shown, the record of what
    /// made it (none for an alignment that came from elsewhere).
    asked: Option<(Asked, Engine)>,
    pub made: Option<(Provenance, Engine)>,
    job: Option<Job<Result<Aligned, String>>>,
    /// Start aligning at the next frame (an example was loaded).
    pub run_when_shown: bool,
    /// The column that is selected (counted from 0): here and, once a structure is lined up with the alignment, there.
    selected: Option<usize>,
    /// The viewer scrolls to the selected column when it is next drawn: the column was selected on another page.
    pub(in crate::ui) reveal: bool,
}

mod column;
mod kept;

impl Default for AlignmentPage {
    fn default() -> Self {
        AlignmentPage {
            inputs: Vec::new(),
            listed: 0,
            paste: String::new(),
            engine: Engine::default(),
            colouring: Colouring::ByType,
            score: Score::default(),
            result: None,
            stats: None,
            revision: 0,
            message: String::new(),
            asked: None,
            made: None,
            job: None,
            run_when_shown: false,
            selected: None,
            reveal: false,
        }
    }
}

impl AlignmentPage {
    /// How big an alignment of the ticked sequences is.
    pub fn weight(&self) -> Weight {
        let (count, letters) = self
            .inputs
            .iter()
            .filter(|i| i.keep)
            .fold((0, 0), |(n, l), i| (n + 1, l + i.bases.len()));
        Weight::of_sequences(count, letters)
    }

    /// Add a sequence to the list (ticked), under a name no other one has.
    pub fn add(&mut self, name: &str, bases: &[u8]) {
        if self
            .inputs
            .iter()
            .any(|i| i.name == name && i.bases == bases)
        {
            return;
        }
        let taken = |candidate: &str| self.inputs.iter().any(|i| i.name == candidate);
        let mut unique = name.to_string();
        let mut n = 1;
        while taken(&unique) {
            n += 1;
            unique = format!("{name} #{n}");
        }
        self.inputs.push(Input {
            name: unique,
            bases: bases.to_vec(),
            keep: true,
        });
        self.listed += 1;
    }

    /// Take every sequence off the list.
    pub fn clear_inputs(&mut self) {
        self.inputs.clear();
        self.listed += 1;
    }

    /// Add the sequences of a file (FASTA, GenBank or bare letters).
    pub fn add_source(&mut self, file: Source) {
        let text = match file.text() {
            Ok(text) => text,
            Err(e) => {
                self.message = format!("Could not read {}: {e}", file.shown());
                return;
            }
        };
        let stem = file
            .path()
            .file_stem()
            .map_or("sequence".to_string(), |s| s.to_string_lossy().into_owned());
        let found = sequences_in(&text, &stem);
        self.message = format!("Added {} sequences from {}.", found.len(), file.shown());
        found.iter().for_each(|(name, bases)| self.add(name, bases));
    }

    /// Choose the engine, as the radio buttons of the page do.
    #[cfg(test)]
    pub fn choose_engine(&mut self, engine: Engine) {
        self.engine = engine;
    }

    /// The inputs and the settings of an alignment as they are now: the ticked sequences and the engine asked for.
    pub fn as_it_is(&self) -> (Vec<provenance::Input>, Vec<(String, String)>) {
        let inputs = self
            .inputs
            .iter()
            .filter(|i| i.keep)
            .map(|i| provenance::Input::of(&i.name, &i.bases))
            .collect();
        let engine = match self.engine {
            Engine::BuiltIn => "built in",
            Engine::Mafft => "MAFFT",
        };
        (inputs, vec![("engine".to_string(), engine.to_string())])
    }

    /// What the viewer draws of the alignment, for drawing it elsewhere than on the page (a report); none without an alignment.
    pub(in crate::ui) fn view<'a>(&'a self, look: &'a Look) -> Option<alnview::View<'a>> {
        let (a, s) = (self.result.as_ref()?, self.stats.as_ref()?);
        Some(alnview::View {
            rows: &a.rows,
            columns: &s.columns,
            bars: &s.bars.values,
            highlight: s.bars.highlight,
            protein: s.protein,
            colouring: self.colouring,
            look,
            scroll: None,
            // a report is of the alignment, not of what was clicked in it
            selected: None,
            reveal: false,
        })
    }

    /// Show `aligned` and work out what the viewer needs.
    pub fn set_result(&mut self, aligned: Aligned) {
        let rows: Vec<Vec<u8>> = aligned.rows.iter().map(|r| r.1.clone()).collect();
        let columns = alnstats::columns(&rows);
        let protein = alnstats::is_protein(&rows);
        // a column of the alignment that was here is not a column of this one
        (self.selected, self.reveal) = (None, false);
        self.revision += 1;
        self.stats = Some(Stats {
            bars: Bars::new(&rows, &columns, protein, self.score),
            columns,
            summary: alnstats::summary(&rows),
            protein,
        });
        self.result = Some(aligned);
    }

    /// Show an alignment that was made elsewhere (a pipeline, the workstation): its rows without the gaps become the
    /// list of sequences, as if they had been aligned here.
    pub fn show_aligned(&mut self, aligned: Aligned, message: &str) {
        // it was not made from the sequences of this page, so there is no record of what made it
        self.made = None;
        self.clear_inputs();
        for (name, row) in &aligned.rows {
            let plain: Vec<u8> = row.iter().copied().filter(|&b| b != b'-').collect();
            self.add(name, &plain);
        }
        self.set_result(aligned);
        self.message = message.to_string();
    }

    pub fn is_running(&self) -> bool {
        self.job.is_some()
    }

    /// Choose the score, as a click on its button does (for the tests of what depends on it).
    #[cfg(test)]
    pub(in crate::ui) fn score_for_test(&mut self, score: Score) {
        self.score = score;
    }

    /// Which version of the scores the bars show: it changes when they do.
    pub(in crate::ui) fn revision(&self) -> u64 {
        self.revision
    }

    /// What the structure page needs of the alignment: its rows, the score of each column as the bars show it, and which
    /// score that is.
    pub(in crate::ui) fn conservation(&self) -> Option<ConservationView<'_>> {
        let (result, stats) = (self.result.as_ref()?, self.stats.as_ref()?);
        Some((&result.rows, &stats.bars.values, stats.bars.score))
    }

    /// The bars of the viewer, for the score chosen now: worked out again when the choice has changed since they were.
    pub fn refresh_scores(&mut self) {
        let (Some(result), Some(stats)) = (&self.result, &mut self.stats) else {
            return;
        };
        if stats.bars.score != self.score {
            let rows: Vec<Vec<u8>> = result.rows.iter().map(|r| r.1.clone()).collect();
            stats.bars = Bars::new(&rows, &stats.columns, stats.protein, self.score);
            self.revision += 1;
        }
    }
}

#[derive(Default)]
struct Clicks {
    made: made_ui::Asks,
    /// What saving a figure of the view said.
    said: String,
    add_open: bool,
    open_fasta: bool,
    add_pasted: bool,
    clear: bool,
    example: bool,
    align: bool,
    save: bool,
    copy: bool,
    to_phylogeny: bool,
    remove: Option<usize>,
    /// The column that was clicked in the viewer.
    column: Option<usize>,
    /// What was clicked on the line of the selected column.
    of_column: column::Asks,
}

/// The sequences of a GenBank or FASTA text; any other text is taken as one sequence.
fn sequences_in(text: &str, fallback_name: &str) -> Vec<Named> {
    if text.trim_start().starts_with("LOCUS") {
        return io::parse_any(text, fallback_name)
            .map(|docs| docs.into_iter().map(|d| (d.name, d.seq)).collect())
            .unwrap_or_default();
    }
    if text.contains('>') {
        return read_fasta(text);
    }
    let bases: Vec<u8> = text
        .bytes()
        .filter(|c| c.is_ascii_alphabetic() || *c == b'-')
        .collect();
    if bases.is_empty() {
        Vec::new()
    } else {
        vec![(fallback_name.to_string(), bases)]
    }
}

impl App {
    pub(super) fn alignment_page(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let look = self.prefs.look();
        let mut c = Clicks::default();
        let have = self.platform.alignment.result.is_some();
        // where this alignment would run, and so whose MAFFT it would use
        let weight = self.platform.alignment.weight();
        let place = self.platform.place(&weight);
        let mafft = self.platform.installed_at("MAFFT", &place);
        let view = self.platform.workstation.view(&weight);
        // the residues of the structure that the selected column was lined up with, as that page names them
        let residues = self
            .platform
            .alignment
            .selected_column()
            .map_or_else(Vec::new, |k| self.platform.structure.residues_in_column(k));
        let page = &mut self.platform.alignment;
        // a column that was selected on another page is scrolled to once, now that this page is drawn
        let reveal = std::mem::take(&mut page.reveal);
        design::page(ui, |ui| {
            design::page_header(
                ui,
                &look,
                "Multiple alignment",
                "Align sequences and read the result with its conservation and consensus",
                |ui| {
                    ui.add_enabled_ui(have, |ui| {
                        c.to_phylogeny =
                            design::primary_button(ui, &look, "Estimate a tree").clicked();
                        c.copy = ui.button("Copy FASTA").clicked();
                        c.save = ui.button("Save FASTA…").clicked();
                        crate::ui::figure::button(ui);
                    });
                },
            );
            design::controls_panel(ui, &look, "aln_inputs", |ui| {
                inputs_panel(ui, &look, page, (mafft.as_ref(), view), &mut c)
            });
            design::results_area(ui, |ui| {
                viewer_panel(ui, &look, page, (&residues, reveal), &mut c)
            });
        });

        self.apply_alignment_clicks(c, ui.ctx(), out);
    }

    fn apply_alignment_clicks(&mut self, c: Clicks, ctx: &egui::Context, out: &mut Vec<Action>) {
        if c.example {
            out.push(Action::Example(Workspace::Alignment));
        }
        if c.add_open {
            let docs = open_sequences(&self.docs);
            let page = &mut self.platform.alignment;
            docs.iter().for_each(|(name, bases)| page.add(name, bases));
            if docs.is_empty() {
                page.message = "No sequence is open: open a file in Sequences first.".into();
            }
        }
        let page = &mut self.platform.alignment;
        if c.open_fasta {
            let ask = dialog::FileDialog::new().set_title("Open sequences");
            if let Some(path) = dialog::pick_one(Wanted::Alignment, ask) {
                page.add_source(Source::Path(&path));
            }
        }
        if c.add_pasted {
            let found = sequences_in(&page.paste, "Pasted sequence");
            found.iter().for_each(|(name, bases)| page.add(name, bases));
            page.message = format!("Added {} pasted sequences.", found.len());
            page.paste.clear();
        }
        if c.clear {
            page.clear_inputs();
        }
        if let Some(k) = c.remove {
            if k < page.inputs.len() {
                page.inputs.remove(k);
                page.listed += 1;
            }
        }
        if c.made.again {
            // with the engine the alignment that is shown was made with
            if let Some((_, engine)) = &page.made {
                page.engine = *engine;
            }
        }
        self.checks_win.open |= c.made.checks;
        if !c.said.is_empty() {
            self.platform.alignment.message = c.said.clone();
        }
        if c.align || c.made.again {
            self.start_alignment(ctx);
        }
        let page = &mut self.platform.alignment;
        if let (true, Some((made, _))) = (c.made.copy, &page.made) {
            ctx.copy_text(made.to_text());
            page.message = "Copied the record of how the alignment was made.".into();
        }
        if let Some(a) = &page.result {
            if c.copy {
                ctx.copy_text(fasta_text(&a.rows));
                page.message = "Copied the alignment as FASTA.".into();
            }
            if c.save {
                let ask = dialog::FileDialog::new().set_title("Save the alignment");
                let fasta = fasta_text(&a.rows);
                if let Some(saved) = dialog::save_bytes(ask, "alignment.fasta", fasta.as_bytes()) {
                    page.message = saved.message;
                }
            }
        }
        if c.to_phylogeny {
            out.push(Action::Go(Workspace::Phylogeny));
        }
        if let Some(k) = c.column {
            // a click on the column that is selected lets it go
            let next = (self.platform.alignment.selected_column() != Some(k)).then_some(k);
            self.platform.select_column(next);
        }
        if c.of_column.clear {
            self.platform.select_column(None);
        }
        if c.of_column.to_structure {
            out.push(Action::Go(Workspace::Structure));
        }
    }

    /// Align the ticked sequences on a worker, so the window stays free.
    pub(super) fn start_alignment(&mut self, ctx: &egui::Context) {
        let place = self.platform.place(&self.platform.alignment.weight());
        let mafft = self.platform.installed_at("MAFFT", &place);
        let page = &mut self.platform.alignment;
        let seqs: Vec<Named> = page
            .inputs
            .iter()
            .filter(|i| i.keep)
            .map(|i| (i.name.clone(), i.bases.clone()))
            .collect();
        if seqs.len() < 2 {
            page.message = "Tick at least two sequences to align.".into();
            return;
        }
        let mafft = if page.engine == Engine::Mafft {
            mafft
        } else {
            None
        };
        // MAFFT that is not installed where this runs cannot be what aligns (its button is greyed, but a choice made earlier
        // stays): the choice goes back to the built-in aligner, so that the page, the record and a later "Run again" all say
        // what really ran
        if mafft.is_none() {
            page.engine = Engine::BuiltIn;
        }
        let (inputs, parameters) = page.as_it_is();
        let ran_on = if matches!(place, Place::There(_)) {
            "the workstation"
        } else {
            "this computer"
        };
        page.asked = Some((
            Asked {
                analysis: "Alignment",
                inputs,
                parameters,
                ran_on,
            },
            page.engine,
        ));
        let Place::There(remote) = place else {
            page.message = "Aligning…".into();
            page.job = Some(Job::spawn(ctx, move || {
                let dir = tools::WorkFolder::new("align")?;
                align(&seqs, mafft.as_ref(), &dir)
            }));
            return;
        };
        page.message = "Aligning on the workstation…".into();
        let label = format!(
            "Alignment of {} sequences ({})",
            seqs.len(),
            if mafft.is_some() { "MAFFT" } else { "built-in" }
        );
        let request = Request::Align {
            seqs,
            mafft: mafft.is_some(),
        };
        page.job = Some(Job::spawn_async(ctx, async move {
            remote
                .run_async(&request, &label)
                .await
                .and_then(Response::into_align)
        }));
    }

    /// Start an alignment that was asked for (an example), and collect a finished one.
    pub(super) fn alignment_tick(&mut self, ctx: &egui::Context) {
        if self.platform.alignment.run_when_shown && !self.platform.alignment.is_running() {
            self.platform.alignment.run_when_shown = false;
            self.start_alignment(ctx);
        }
        self.poll_alignment(ctx);
        self.platform.alignment.refresh_scores();
    }

    fn poll_alignment(&mut self, ctx: &egui::Context) {
        let page = &mut self.platform.alignment;
        let Some(polled) = page.job.as_ref().map(Job::poll) else {
            return;
        };
        match polled {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(Ok(aligned)) => {
                page.message = format!(
                    "Aligned {} sequences into {} columns with {} in {:.2} s.",
                    aligned.rows.len(),
                    aligned.columns(),
                    aligned.engine,
                    aligned.seconds
                );
                page.made = page.asked.take().map(|(asked, engine)| {
                    (
                        asked.made(vec![aligned.engine.clone()], aligned.seconds),
                        engine,
                    )
                });
                page.set_result(aligned);
                page.job = None;
            }
            Poll::Done(Err(why)) => (page.message, page.job, page.asked) = (why, None, None),
            Poll::Died => {
                page.asked = None;
                (page.message, page.job) =
                    ("The alignment stopped without an answer.".into(), None);
            }
        }
    }
}

/// `mafft` is the program that would run the alignment, and `view` says where that is.
fn inputs_panel(
    ui: &mut Ui,
    look: &Look,
    page: &mut AlignmentPage,
    (mafft, view): (Option<&Installed>, ComputeView),
    c: &mut Clicks,
) {
    design::section(ui, look, "Sequences");
    ui.horizontal_wrapped(|ui| {
        c.add_open = ui.button("Add open sequences").clicked();
        c.open_fasta = ui.button("Open FASTA…").clicked();
        c.clear = ui.button("Clear").clicked();
    });
    ui.add(
        egui::TextEdit::multiline(&mut page.paste)
            .hint_text("Paste FASTA, or one sequence")
            .font(egui::TextStyle::Monospace)
            .desired_rows(2)
            .desired_width(f32::INFINITY),
    );
    c.add_pasted = ui
        .add_enabled(
            !page.paste.trim().is_empty(),
            egui::Button::new("Add pasted"),
        )
        .clicked();
    ui.add_space(4.0);
    let kept = page.inputs.iter().filter(|i| i.keep).count();
    ui.label(
        RichText::new(format!("{} sequences, {kept} ticked", page.inputs.len()))
            .size(12.0)
            .color(look.dim),
    );
    ScrollArea::vertical()
        .max_height(LIST_HEIGHT)
        .auto_shrink([false, true])
        .show(ui, |ui| {
            for (k, input) in page.inputs.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.checkbox(&mut input.keep, "");
                    let name: String = input.name.chars().take(NAME_CHARS).collect();
                    ui.label(name).on_hover_text(&input.name);
                    ui.label(
                        RichText::new(input.bases.len().to_string())
                            .size(11.5)
                            .color(look.dim),
                    );
                    if ui.small_button("×").on_hover_text("Remove").clicked() {
                        c.remove = Some(k);
                    }
                });
            }
        });

    design::section(ui, look, "Method");
    ui.radio_value(&mut page.engine, Engine::BuiltIn, "Built-in (centre-star)")
        .on_hover_text("For closely related DNA sequences; nothing to install");
    let label = match (mafft, view.there) {
        (Some(m), _) => format!(
            "MAFFT {}",
            m.version.split_whitespace().next().unwrap_or("")
        ),
        (None, false) => "MAFFT (not installed)".to_string(),
        (None, true) => "MAFFT (not on the workstation)".to_string(),
    };
    ui.add_enabled_ui(mafft.is_some(), |ui| {
        ui.radio_value(&mut page.engine, Engine::Mafft, label)
            .on_hover_text("Handles DNA and proteins, and distant sequences");
    });
    let state = (kept >= 2, page.is_running());
    c.align = run_button(ui, look, view, "Align", state);
    // right under the button, where a failure (the workstation cannot be reached, say) is seen
    if !page.message.is_empty() {
        ui.label(RichText::new(&page.message).size(12.0).color(look.dim));
    }

    design::section(ui, look, "Conservation");
    for score in Score::ALL {
        ui.radio_value(&mut page.score, score, score.label())
            .on_hover_text(score.explain());
    }
    ui.label(
        RichText::new(page.score.explain())
            .size(11.5)
            .color(look.dim),
    );

    design::section(ui, look, "Colours");
    ui.radio_value(&mut page.colouring, Colouring::ByType, "By residue type");
    ui.radio_value(
        &mut page.colouring,
        Colouring::Differences,
        "Differences from the consensus",
    );
    ui.radio_value(&mut page.colouring, Colouring::Plain, "Plain");
}

/// `linked` is the residues of the structure that the selected column was lined up with, and whether the viewer scrolls
/// to that column in this frame.
fn viewer_panel(
    ui: &mut Ui,
    look: &Look,
    page: &AlignmentPage,
    linked: (&[String], bool),
    c: &mut Clicks,
) {
    let (residues, reveal) = linked;
    let (Some(a), Some(s)) = (&page.result, &page.stats) else {
        let pick = design::empty_state(
            ui,
            look,
            "No alignment yet",
            "Add the open sequences, a FASTA file or pasted text, then press Align. \
             The result appears here with a consensus row and a bar for how conserved each column is.",
            &["Add open sequences", "Open FASTA…", "Try the example"],
        );
        match pick {
            Some(0) => c.add_open = true,
            Some(1) => c.open_fasta = true,
            Some(2) => c.example = true,
            _ => {}
        }
        return;
    };
    let m = &s.summary;
    ui.horizontal_wrapped(|ui| {
        design::chip(ui, look, &format!("{} sequences", m.rows), Tone::Neutral);
        design::chip(ui, look, &format!("{} columns", m.columns), Tone::Neutral);
        design::chip(
            ui,
            look,
            &format!("identity {:.1} %", m.identity),
            Tone::Accent,
        );
        design::chip(ui, look, &format!("gaps {:.1} %", m.gaps), Tone::Neutral);
        design::chip(ui, look, &format!("{} variable", m.variable), Tone::Neutral);
        let blocks = s.bars.blocks.len();
        let plural = if blocks == 1 { "" } else { "s" };
        design::chip(
            ui,
            look,
            &format!("{blocks} conserved block{plural}"),
            Tone::Accent,
        );
        design::chip(
            ui,
            look,
            &format!("{} informative", m.informative),
            Tone::Neutral,
        );
        design::chip(
            ui,
            look,
            if s.protein { "protein" } else { "DNA" },
            Tone::Good,
        );
        design::chip(
            ui,
            look,
            &format!("{} · {:.2} s", a.engine, a.seconds),
            Tone::Neutral,
        );
    });
    if let Some((made, _)) = &page.made {
        c.made = made_ui::section(ui, look, made, || page.as_it_is(), !page.is_running());
    }
    if !s.bars.blocks.is_empty() {
        let shown: Vec<String> = s
            .bars
            .blocks
            .iter()
            .take(MAX_BLOCKS_SHOWN)
            .map(|b| format!("{}-{}", b.from + 1, b.to + 1))
            .collect();
        let more = s.bars.blocks.len().saturating_sub(MAX_BLOCKS_SHOWN);
        let what = if s.bars.score == Score::Identity {
            format!("Runs of {MIN_BLOCK} or more identical columns")
        } else {
            format!(
                "Runs of {MIN_BLOCK} or more columns in the top quarter of the scores (from {:.2})",
                s.bars.highlight
            )
        };
        let tail = if more > 0 {
            format!(" and {more} more")
        } else {
            String::new()
        };
        ui.label(
            RichText::new(format!("{what}: columns {}{tail}.", shown.join(", ")))
                .size(12.0)
                .color(look.dim),
        );
    }
    if s.bars.score != Score::Identity && m.rows < FEW_SEQUENCES {
        ui.label(
            RichText::new(format!(
                "With {} sequences a conservation score is rough: these scores are made for alignments of many sequences (fewer than {FEW_SEQUENCES} say little).",
                m.rows
            ))
            .size(12.0)
            .color(look.warn),
        );
    }
    let facts = page.selected.and_then(|k| page.column_facts(k));
    c.of_column = column::line(ui, look, facts.as_ref(), residues);
    ui.add_space(6.0);
    let view = |scroll, reveal| alnview::View {
        rows: &a.rows,
        columns: &s.columns,
        bars: &s.bars.values,
        highlight: s.bars.highlight,
        protein: s.protein,
        colouring: page.colouring,
        look,
        scroll,
        selected: page.selected,
        reveal,
    };
    if let Some(format) = crate::ui::figure::wanted() {
        // what the window shows: the same room, scrolled to the same place
        let at = alnview::scrolled(ui.ctx());
        c.said = crate::ui::figure::save(format, "alignment", ui.available_size(), |ui| {
            alnview::show(ui, &view(Some(at), false));
        });
    }
    design::view_frame(ui, look, |ui| {
        c.column = alnview::show(ui, &view(None, reveal))
    });
}

#[cfg(test)]
mod tests;
