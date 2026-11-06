//! Alignment: gather sequences, align them with the built-in aligner or MAFFT, and read the result in a
//! colour-coded viewer with a consensus, conservation and identity.

use super::alnview::{self, Colouring};
use super::route::Source;
use super::workstation::{compute_picker, ComputeView, Place};
use super::Workspace;
use crate::pipeline::steps::{align, fasta_text, read_fasta, Aligned, Named};
use crate::pipeline::tools::{self, Installed};
use crate::remote::tasks::{Request, Response};
use crate::remote::weight::Weight;
use crate::seq::alnstats::{self, Column, Summary};
use crate::seq::io;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::dialog::{self, Wanted};
use crate::ui::job::{Job, Poll};
use crate::ui::windows::tree::open_sequences;
use eframe::egui::{self, RichText, ScrollArea, Ui};
use std::time::Duration;

const INPUTS_WIDTH: f32 = 300.0;
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

/// What a reader wants to know about the alignment, worked out once when it arrives.
pub struct Stats {
    columns: Vec<Column>,
    summary: Summary,
    protein: bool,
}

pub struct AlignmentPage {
    pub inputs: Vec<Input>,
    paste: String,
    engine: Engine,
    colouring: Colouring,
    pub result: Option<Aligned>,
    stats: Option<Stats>,
    pub message: String,
    job: Option<Job<Result<Aligned, String>>>,
    /// Start aligning at the next frame (an example was loaded).
    pub run_when_shown: bool,
}

impl Default for AlignmentPage {
    fn default() -> Self {
        AlignmentPage {
            inputs: Vec::new(),
            paste: String::new(),
            engine: Engine::default(),
            colouring: Colouring::ByType,
            result: None,
            stats: None,
            message: String::new(),
            job: None,
            run_when_shown: false,
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

    /// Show `aligned` and work out what the viewer needs.
    pub fn set_result(&mut self, aligned: Aligned) {
        let rows: Vec<Vec<u8>> = aligned.rows.iter().map(|r| r.1.clone()).collect();
        self.stats = Some(Stats {
            columns: alnstats::columns(&rows),
            summary: alnstats::summary(&rows),
            protein: alnstats::is_protein(&rows),
        });
        self.result = Some(aligned);
    }

    /// Show an alignment that was made elsewhere (a pipeline, the workstation): its rows without the gaps become the
    /// list of sequences, as if they had been aligned here.
    pub fn show_aligned(&mut self, aligned: Aligned, message: &str) {
        self.inputs.clear();
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
}

#[derive(Default)]
struct Clicks {
    add_open: bool,
    open_fasta: bool,
    add_pasted: bool,
    clear: bool,
    align: bool,
    save: bool,
    copy: bool,
    to_phylogeny: bool,
    remove: Option<usize>,
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
        let look = design::look(self.prefs.dark);
        let mut c = Clicks::default();
        let have = self.platform.alignment.result.is_some();
        // where this alignment would run, and so whose MAFFT it would use
        let weight = self.platform.alignment.weight();
        let place = self.platform.place(&weight);
        let mafft = self.platform.installed_at("MAFFT", &place);
        let view = self.platform.workstation.view(&weight);
        let page = &mut self.platform.alignment;
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
                    });
                },
            );
            egui::Panel::left("aln_inputs")
                .exact_size(INPUTS_WIDTH)
                .resizable(false)
                .frame(egui::Frame::new())
                .show(ui, |ui| {
                    ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            design::card(ui, &look, |ui| {
                                inputs_panel(ui, &look, page, (mafft.as_ref(), view), &mut c)
                            });
                        });
                });
            egui::CentralPanel::default()
                .frame(egui::Frame::new().inner_margin(egui::Margin {
                    left: design::PAGE_MARGIN,
                    ..egui::Margin::ZERO
                }))
                .show(ui, |ui| viewer_panel(ui, &look, page));
        });

        self.apply_alignment_clicks(c, ui.ctx(), out);
    }

    fn apply_alignment_clicks(&mut self, c: Clicks, ctx: &egui::Context, out: &mut Vec<Action>) {
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
            page.inputs.clear();
        }
        if let Some(k) = c.remove {
            if k < page.inputs.len() {
                page.inputs.remove(k);
            }
        }
        if c.align {
            self.start_alignment(ctx);
        }
        let page = &mut self.platform.alignment;
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
        let Place::There(remote) = place else {
            page.message = "Aligning…".into();
            page.job = Some(Job::spawn(ctx, move || {
                let dir = tools::work_folder("align")?;
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
        page.job = Some(Job::spawn(ctx, move || {
            remote.run(&request, &label).and_then(Response::into_align)
        }));
    }

    /// Start an alignment that was asked for (an example), and collect a finished one.
    pub(super) fn alignment_tick(&mut self, ctx: &egui::Context) {
        if self.platform.alignment.run_when_shown && !self.platform.alignment.is_running() {
            self.platform.alignment.run_when_shown = false;
            self.start_alignment(ctx);
        }
        self.poll_alignment(ctx);
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
                page.set_result(aligned);
                page.job = None;
            }
            Poll::Done(Err(why)) => (page.message, page.job) = (why, None),
            Poll::Died => {
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
    compute_picker(ui, look, view);
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
    ui.add_space(6.0);
    ui.horizontal(|ui| {
        let ready = kept >= 2 && !page.is_running();
        c.align = ui
            .add_enabled(ready, egui::Button::new(RichText::new("Align").strong()))
            .clicked();
        if page.is_running() {
            ui.spinner();
        }
    });
    // right under the button, where a failure (the workstation cannot be reached, say) is seen
    if !page.message.is_empty() {
        ui.label(RichText::new(&page.message).size(12.0).color(look.dim));
    }

    design::section(ui, look, "Colours");
    ui.radio_value(&mut page.colouring, Colouring::ByType, "By residue type");
    ui.radio_value(
        &mut page.colouring,
        Colouring::Differences,
        "Differences from the consensus",
    );
    ui.radio_value(&mut page.colouring, Colouring::Plain, "Plain");
}

fn viewer_panel(ui: &mut Ui, look: &Look, page: &AlignmentPage) {
    let (Some(a), Some(s)) = (&page.result, &page.stats) else {
        design::card(ui, look, |ui| {
            ui.label(
                RichText::new("No alignment yet")
                    .size(17.0)
                    .strong()
                    .color(look.text),
            );
            ui.label(
                RichText::new(
                    "Add the open sequences, a FASTA file or pasted text on the left, then press Align. \
                     The result appears here with a consensus row and a bar for how conserved each column is.",
                )
                .color(look.dim),
            );
        });
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
    ui.add_space(6.0);
    egui::Frame::new()
        .fill(look.surface)
        .stroke(egui::Stroke::new(1.0, look.border))
        .corner_radius(egui::CornerRadius::same(design::CARD_RADIUS))
        .show(ui, |ui| {
            alnview::show(
                ui,
                &alnview::View {
                    rows: &a.rows,
                    columns: &s.columns,
                    protein: s.protein,
                    colouring: page.colouring,
                    look,
                },
            );
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_sequence_with_the_same_name_gets_a_number() {
        let mut page = AlignmentPage::default();
        page.add("a", b"ACGT");
        page.add("a", b"ACGT"); // the very same sequence is not added twice
        page.add("a", b"TTTT");
        page.add("a", b"GGGG");
        let names: Vec<&str> = page.inputs.iter().map(|i| i.name.as_str()).collect();
        assert_eq!(names, ["a", "a #2", "a #3"]);
        assert!(page.inputs.iter().all(|i| i.keep));
    }

    #[test]
    fn pasted_text_is_fasta_or_one_sequence() {
        let fasta = sequences_in(">one\nACGT\nAC\n>two\nGGGG\n", "x");
        assert_eq!(
            fasta,
            vec![
                ("one".into(), b"ACGTAC".to_vec()),
                ("two".into(), b"GGGG".to_vec())
            ]
        );
        let bare = sequences_in("acgt acgt\n123\n", "Pasted");
        assert_eq!(bare, vec![("Pasted".to_string(), b"acgtacgt".to_vec())]);
        assert!(sequences_in("  \n 12 \n", "x").is_empty());
    }
}
