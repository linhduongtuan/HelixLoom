//! The Contig Assembler tab: choose sequences, set the overlap rules, assemble, inspect the
//! contigs and send a consensus to the editor.

use super::contig_view::{self, compare_to_consensus, ContigView};
use crate::seq::assemble::{assemble, normalize, Assembly, Contig, Params, Read};
use crate::seq::doc::Doc;
use crate::seq::features::{Feature, Rgb};
use crate::seq::io;
use crate::ui::app::App;
use crate::ui::dialog::{self, Wanted};
use crate::ui::map::commify;
use crate::ui::platform::route::Source;
use crate::ui::theme::Palette;
use eframe::egui::{self, RichText, ScrollArea, Ui};

const READ_FWD_COLOR: Rgb = [153, 193, 255];
const READ_REV_COLOR: Rgb = [255, 200, 140];

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Preset {
    Reads,
    Gibson,
    Custom,
}

pub struct ReadItem {
    pub name: String,
    pub seq: Vec<u8>,
    pub on: bool,
}

/// An assembly result together with the (cleaned) sequences it was built from.
pub struct Outcome {
    pub asm: Assembly,
    pub inputs: Vec<(String, Vec<u8>)>,
    /// For each input, its index in `AssemblyWin::reads`.
    pub origin: Vec<usize>,
}

pub struct AssemblyWin {
    pub reads: Vec<ReadItem>,
    pub params: Params,
    pub preset: Preset,
    pub outcome: Option<Outcome>,
    pub contig: usize,
    pub view: ContigView,
    pub paste_open: bool,
    pub paste_text: String,
    pub message: String,
}

impl Default for AssemblyWin {
    fn default() -> Self {
        AssemblyWin {
            reads: Vec::new(),
            params: Params::reads(),
            preset: Preset::Reads,
            outcome: None,
            contig: 0,
            view: ContigView::default(),
            paste_open: false,
            paste_text: String::new(),
            message: String::new(),
        }
    }
}

impl AssemblyWin {
    pub fn add_reads(&mut self, items: impl IntoIterator<Item = (String, Vec<u8>)>) -> usize {
        let before = self.reads.len();
        self.reads.extend(
            items
                .into_iter()
                .filter(|(_, s)| !s.is_empty())
                .map(|(name, seq)| ReadItem {
                    name,
                    seq,
                    on: true,
                }),
        );
        self.reads.len() - before
    }

    pub fn load_demo(&mut self) {
        self.reads.clear();
        self.outcome = None;
        let n = self.add_reads(
            crate::sample::demo_reads()
                .into_iter()
                .map(|r| (r.name, r.seq)),
        );
        self.set_preset(Preset::Reads);
        self.message =
            format!("Loaded {n} demo reads sampled from a synthetic plasmid. Press Assemble.");
    }

    pub fn set_preset(&mut self, p: Preset) {
        self.preset = p;
        match p {
            Preset::Reads => self.params = Params::reads(),
            Preset::Gibson => self.params = Params::gibson(),
            Preset::Custom => {}
        }
    }

    fn sync_preset(&mut self) {
        self.preset = if self.params == Params::reads() {
            Preset::Reads
        } else if self.params == Params::gibson() {
            Preset::Gibson
        } else {
            Preset::Custom
        };
    }

    pub fn run(&mut self) {
        let used: Vec<usize> = (0..self.reads.len())
            .filter(|&i| self.reads[i].on)
            .collect();
        if used.len() < 2 && !(used.len() == 1 && self.params.circularize) {
            self.message =
                "Add at least two sequences to assemble (or one that closes into a circle).".into();
            return;
        }
        let inputs: Vec<(String, Vec<u8>)> = used
            .iter()
            .map(|&i| (self.reads[i].name.clone(), normalize(&self.reads[i].seq)))
            .collect();
        let reads: Vec<Read> = inputs
            .iter()
            .map(|(n, s)| Read {
                name: n.clone(),
                seq: s.clone(),
            })
            .collect();
        let asm = assemble(&reads, &self.params);
        self.message = summary(&asm, inputs.len());
        self.contig = 0;
        self.view.selected = None;
        self.view.fit_pending = true;
        self.outcome = Some(Outcome {
            asm,
            inputs,
            origin: used,
        });
    }
}

fn summary(asm: &Assembly, total: usize) -> String {
    let in_contigs: usize = asm.contigs.iter().map(|c| c.placed.len()).sum();
    if asm.contigs.is_empty() {
        return format!("No contigs: none of the {total} sequences overlap enough. Try a shorter minimum overlap or a lower identity.");
    }
    format!(
        "{} contig{} from {in_contigs} of {total} sequences; {} left unassembled. Total {} bp, N50 {} bp.",
        asm.contigs.len(),
        if asm.contigs.len() == 1 { "" } else { "s" },
        asm.unassembled.len(),
        commify(asm.total_length),
        commify(asm.n50)
    )
}

/// The consensus as an annotated sequence: every read becomes a feature on its strand.
pub fn contig_to_doc(contig: &Contig, inputs: &[(String, Vec<u8>)]) -> Doc {
    let n = contig.consensus.len();
    let mut doc = Doc::new(contig.name.clone(), contig.consensus.clone());
    doc.circular = contig.circular;
    doc.definition = format!(
        "Consensus of {} reads assembled by HelixLoom ({} bp, mean coverage {:.1}x{}).",
        contig.placed.len(),
        n,
        contig.mean_depth(),
        if contig.circular { ", circular" } else { "" }
    );
    for pl in &contig.placed {
        let len = inputs[pl.read].1.len();
        let start = if contig.circular {
            pl.offset % n.max(1)
        } else {
            pl.offset
        };
        let spans = if contig.circular && start + len > n {
            vec![(0, start + len - n), (start, n)]
        } else {
            vec![(start, (start + len).min(n))]
        };
        let mut f = Feature::new(inputs[pl.read].0.clone(), "misc_feature", 0, 1, pl.reverse);
        f.spans = spans;
        f.fwd_color = READ_FWD_COLOR;
        f.rev_color = READ_REV_COLOR;
        let how = if pl.reverse {
            "reverse complement"
        } else {
            "forward"
        };
        f.qualifiers.push((
            "note".into(),
            format!(
                "assembled read, {how}{}",
                if pl.contained {
                    ", inside a longer read"
                } else {
                    ""
                }
            ),
        ));
        doc.features.push(f);
    }
    doc.features.sort_by_key(|f| (f.start(), f.name.clone()));
    doc.dirty = true;
    doc
}

fn fasta_of(contigs: &[&Contig]) -> String {
    contigs
        .iter()
        .map(|c| {
            let d = Doc::new(c.name.clone(), c.consensus.clone());
            let note = format!(
                "{} reads, {}",
                c.placed.len(),
                if c.circular { "circular" } else { "linear" }
            );
            io::to_fasta(&d, &format!("{} {note}", c.name))
        })
        .collect()
}

impl App {
    pub fn assembly_tab(&mut self, ui: &mut Ui) {
        let pal = self.palette();
        egui::Panel::left("asm_reads")
            .resizable(true)
            .default_size(320.0)
            .min_size(240.0)
            .show(ui, |ui| self.reads_panel(ui));
        egui::CentralPanel::default().show(ui, |ui| self.results_panel(ui, &pal));
    }

    /// Everything open in the editor, as assembly input.
    pub fn add_open_sequences(&mut self) {
        let items: Vec<(String, Vec<u8>)> = self
            .docs
            .iter()
            .filter(|d| !d.doc.is_empty())
            .map(|d| (d.doc.name.clone(), d.doc.seq.clone()))
            .collect();
        let n = self.assembly.add_reads(items);
        self.assembly.message =
            format!("Added {n} open sequence{}.", if n == 1 { "" } else { "s" });
    }

    pub fn add_read_files(&mut self, paths: &[std::path::PathBuf]) {
        let files: Vec<Source> = paths.iter().map(|p| Source::Path(p)).collect();
        self.add_read_sources(&files);
    }

    /// Add the sequences in some files: the ones on this computer that were named, or the ones a browser brought.
    pub fn add_read_sources(&mut self, files: &[Source]) {
        let mut items = Vec::new();
        for file in files {
            let stem = file
                .path()
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("read")
                .to_string();
            match file
                .bytes()
                .map_err(|e| e.to_string())
                .and_then(|b| io::parse_any(&String::from_utf8_lossy(&b), &stem))
            {
                Ok(docs) => items.extend(docs.into_iter().map(|d| (d.name, d.seq))),
                Err(e) => {
                    self.assembly.message = format!("{}: {e}", file.shown());
                    return;
                }
            }
        }
        let n = self.assembly.add_reads(items);
        self.assembly.message = format!(
            "Added {n} sequence{} from {} file{}.",
            if n == 1 { "" } else { "s" },
            files.len(),
            if files.len() == 1 { "" } else { "s" }
        );
    }

    fn reads_panel(&mut self, ui: &mut Ui) {
        ui.heading("Sequences to assemble");
        ui.horizontal_wrapped(|ui| {
            if ui
                .button("Add open sequences")
                .on_hover_text("Use every sequence that is open in the editor")
                .clicked()
            {
                self.add_open_sequences();
            }
            if ui.button("Add from file…").clicked() {
                let ask = dialog::FileDialog::new().set_title("Add sequences");
                let files = dialog::pick_many(Wanted::AssemblyReads, ask);
                self.add_read_files(&files);
            }
            if ui.button("Paste FASTA…").clicked() {
                self.assembly.paste_open = true;
            }
            if ui
                .button("Demo reads")
                .on_hover_text("Reads sampled from a made-up plasmid")
                .clicked()
            {
                self.assembly.load_demo();
            }
            if ui.button("Clear").clicked() {
                self.assembly.reads.clear();
                self.assembly.outcome = None;
            }
        });
        ui.label(
            RichText::new(format!(
                "{} sequences, {} selected. You can also drop files onto the window.",
                self.assembly.reads.len(),
                self.assembly.reads.iter().filter(|r| r.on).count()
            ))
            .small(),
        );
        ui.separator();
        let placement: Vec<Option<(usize, bool)>> = {
            let mut v = vec![None; self.assembly.reads.len()];
            if let Some(o) = &self.assembly.outcome {
                for (ci, c) in o.asm.contigs.iter().enumerate() {
                    for pl in &c.placed {
                        v[o.origin[pl.read]] = Some((ci, pl.reverse));
                    }
                }
            }
            v
        };
        let mut remove: Option<usize> = None;
        ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                for (i, r) in self.assembly.reads.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        ui.checkbox(&mut r.on, "");
                        ui.label(RichText::new(&r.name).monospace());
                        ui.label(
                            RichText::new(format!("{} bp", commify(r.seq.len())))
                                .small()
                                .weak(),
                        );
                        if let Some(Some((ci, rev))) = placement.get(i) {
                            ui.label(
                                RichText::new(format!(
                                    "{} {}",
                                    if *rev { "◀" } else { "▶" },
                                    self.assembly
                                        .outcome
                                        .as_ref()
                                        .map_or("", |o| o.asm.contigs[*ci].name.as_str())
                                ))
                                .small(),
                            );
                        } else if self.assembly.outcome.is_some() && r.on {
                            ui.label(RichText::new("no overlap").small().weak());
                        }
                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            if ui
                                .small_button("×")
                                .on_hover_text("Remove from the list")
                                .clicked()
                            {
                                remove = Some(i);
                            }
                        });
                    });
                }
            });
        if let Some(i) = remove {
            self.assembly.reads.remove(i);
            self.assembly.outcome = None;
        }
    }

    pub fn assembly_windows(&mut self, ctx: &egui::Context) {
        if !self.assembly.paste_open {
            return;
        }
        let mut open = true;
        crate::ui::design::window("Paste FASTA or sequences")
            .open(&mut open)
            .default_size([460.0, 300.0])
            .show(ctx, |ui| {
                ui.label(
                    "One or more records starting with >name. Plain sequence text is accepted too.",
                );
                ui.add(
                    egui::TextEdit::multiline(&mut self.assembly.paste_text)
                        .font(egui::TextStyle::Monospace)
                        .desired_rows(10)
                        .desired_width(f32::INFINITY),
                );
                if ui.button("Add").clicked() {
                    match io::parse_any(&self.assembly.paste_text, "pasted") {
                        Ok(docs) => {
                            let n = self
                                .assembly
                                .add_reads(docs.into_iter().map(|d| (d.name, d.seq)));
                            self.assembly.message = format!(
                                "Added {n} pasted sequence{}.",
                                if n == 1 { "" } else { "s" }
                            );
                            self.assembly.paste_text.clear();
                            self.assembly.paste_open = false;
                        }
                        Err(e) => self.assembly.message = e,
                    }
                }
            });
        if !open {
            self.assembly.paste_open = false;
        }
    }

    fn parameter_panel(&mut self, ui: &mut Ui) {
        let w = &mut self.assembly;
        ui.horizontal(|ui| {
            ui.label("Rules");
            let mut preset = w.preset;
            ui.selectable_value(&mut preset, Preset::Reads, "Sequencing reads").on_hover_text("Overlaps of 20+ bases that match at least 90%; contigs stay linear unless their ends overlap");
            ui.selectable_value(&mut preset, Preset::Gibson, "Gibson-style fragments")
                .on_hover_text("The original's Gibson assembler defaults: exact 12–100 base end overlaps, circular product only");
            if preset != w.preset && preset != Preset::Custom {
                w.set_preset(preset);
            }
            if w.preset == Preset::Custom {
                ui.label(RichText::new("custom").small().italics());
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.label("Minimum overlap");
            ui.add(
                egui::DragValue::new(&mut w.params.min_overlap)
                    .range(4..=2000)
                    .suffix(" bp"),
            );
            ui.label("Minimum identity");
            let mut pct = w.params.min_identity * 100.0;
            ui.add(
                egui::DragValue::new(&mut pct)
                    .range(50.0..=100.0)
                    .speed(0.5)
                    .suffix(" %")
                    .max_decimals(1),
            );
            w.params.min_identity = pct / 100.0;
            let mut limited = w.params.max_overlap.is_some();
            ui.checkbox(&mut limited, "Limit overlap to");
            let mut max = w.params.max_overlap.unwrap_or(100);
            ui.add_enabled(
                limited,
                egui::DragValue::new(&mut max)
                    .range(w.params.min_overlap..=5000)
                    .suffix(" bp"),
            );
            w.params.max_overlap = limited.then_some(max);
        });
        ui.horizontal_wrapped(|ui| {
            ui.checkbox(&mut w.params.both_strands, "Sequences can be reversed");
            ui.checkbox(
                &mut w.params.allow_containment,
                "Allow sequences inside longer ones",
            );
            ui.checkbox(
                &mut w.params.circularize,
                "Close into a circle when the ends overlap",
            );
            ui.checkbox(&mut w.params.only_circular, "Keep circular products only");
        });
        w.sync_preset();
    }

    fn results_panel(&mut self, ui: &mut Ui, pal: &Palette) {
        self.parameter_panel(ui);
        ui.horizontal(|ui| {
            if ui
                .add(
                    egui::Button::new(RichText::new("  Assemble  ").strong())
                        .fill(pal.accent.gamma_multiply(0.35)),
                )
                .clicked()
            {
                self.assembly.run();
            }
            ui.label(RichText::new(&self.assembly.message).small());
        });
        ui.separator();
        let Some(outcome) = &self.assembly.outcome else {
            ui.add_space(24.0);
            ui.vertical_centered(|ui| {
                ui.label(RichText::new("Contig Assembler").heading());
                ui.label("Add sequences on the left (open sequences, files, pasted FASTA, or the demo reads),");
                ui.label("choose the overlap rules above, then press Assemble.");
            });
            return;
        };
        let (mut wanted_doc, mut copy, mut save_one, mut save_all) = (None, false, false, false);
        let count = outcome.asm.contigs.len();
        if count == 0 {
            ui.label("Nothing joined. The sequences are listed on the left.");
            return;
        }
        self.assembly.contig = self.assembly.contig.min(count - 1);
        egui::Grid::new("contig_table")
            .striped(true)
            .spacing([14.0, 3.0])
            .show(ui, |ui| {
                for h in [
                    "Contig",
                    "Length",
                    "Reads",
                    "Coverage",
                    "Shape",
                    "Disagreements",
                ] {
                    ui.label(RichText::new(h).strong());
                }
                ui.end_row();
                for (i, c) in outcome.asm.contigs.iter().enumerate() {
                    if ui
                        .selectable_label(self.assembly.contig == i, &c.name)
                        .clicked()
                    {
                        self.assembly.contig = i;
                        self.assembly.view.selected = None;
                        self.assembly.view.fit_pending = true;
                    }
                    ui.label(
                        RichText::new(format!("{} bp", commify(c.consensus.len()))).monospace(),
                    );
                    ui.label(c.placed.len().to_string());
                    ui.label(format!("{:.1}×", c.mean_depth()));
                    ui.label(if c.circular { "circular" } else { "linear" });
                    ui.label(if c.conflicts == 0 {
                        "none".to_string()
                    } else {
                        c.conflicts.to_string()
                    });
                    ui.end_row();
                }
            });
        if !outcome.asm.unassembled.is_empty() {
            let names: Vec<&str> = outcome
                .asm
                .unassembled
                .iter()
                .map(|&i| outcome.inputs[i].0.as_str())
                .collect();
            ui.label(
                RichText::new(format!("Unassembled: {}", names.join(", ")))
                    .small()
                    .weak(),
            );
        }
        ui.separator();
        let contig = &outcome.asm.contigs[self.assembly.contig];
        ui.horizontal(|ui| {
            let shape = if contig.circular {
                format!(
                    "circular, closed by a {} bp overlap",
                    contig.closing_overlap
                )
            } else {
                "linear".to_string()
            };
            ui.label(
                RichText::new(format!(
                    "{} — {} bp, {shape}",
                    contig.name,
                    commify(contig.consensus.len())
                ))
                .strong(),
            );
            ui.add_space(12.0);
            ui.label("Zoom");
            ui.add(
                egui::Slider::new(
                    &mut self.assembly.view.zoom,
                    contig_view::MIN_ZOOM..=contig_view::MAX_ZOOM,
                )
                .logarithmic(true)
                .show_value(false),
            );
            if ui
                .small_button("Fit")
                .on_hover_text("Show the whole contig")
                .clicked()
            {
                self.assembly.view.fit_pending = true;
            }
            ui.label(
                RichText::new(
                    "blue = forward, orange = reverse, red = disagrees with the consensus",
                )
                .small()
                .weak(),
            );
        });
        ui.horizontal(|ui| {
            if ui
                .button("Open consensus in editor")
                .on_hover_text("Creates a new sequence; every read becomes a feature")
                .clicked()
            {
                wanted_doc = Some(contig_to_doc(contig, &outcome.inputs));
            }
            copy = ui.button("Copy consensus").clicked();
            save_one = ui.button("Save consensus FASTA…").clicked();
            save_all = count > 1 && ui.button("Save all contigs FASTA…").clicked();
        });
        if let Some(k) = self
            .assembly
            .view
            .selected
            .filter(|&k| k < contig.placed.len())
        {
            let pl = &contig.placed[k];
            let (name, seq) = &outcome.inputs[pl.read];
            let (mism, compared) = compare_to_consensus(contig, pl, seq);
            ui.label(RichText::new(format!(
                "{name}: {} • bases {}–{} of the contig • {} bp • {mism} disagreement{} in {} bases{}",
                if pl.reverse { "reverse complement" } else { "forward" },
                commify(pl.offset + 1),
                commify(pl.offset + seq.len()),
                commify(seq.len()),
                if mism == 1 { "" } else { "s" },
                commify(compared),
                if pl.contained { " • inside a longer read" } else { "" }
            )).small());
        }
        let contig_ref = &outcome.asm.contigs[self.assembly.contig];
        let inputs_ref = &outcome.inputs;
        if copy {
            ui.ctx()
                .copy_text(String::from_utf8_lossy(&contig_ref.consensus).into_owned());
        }
        if save_one || save_all {
            let all: Vec<&Contig> = if save_all {
                outcome.asm.contigs.iter().collect()
            } else {
                vec![contig_ref]
            };
            let text = fasta_of(&all);
            let name = if save_all {
                "contigs.fasta".to_string()
            } else {
                format!("{}.fasta", contig_ref.name)
            };
            let ask = dialog::FileDialog::new().set_title("Save FASTA");
            if let Some(saved) = dialog::save_bytes(ask, &name, text.as_bytes()) {
                // the answer of the contig view borrows the assembly, so the line and the toast are set field by field
                self.toasts
                    .push(&saved.message, crate::ui::toasts::tone_of(&saved.message));
                self.status = saved.message;
            }
        }
        contig_view::show(ui, &mut self.assembly.view, contig_ref, inputs_ref, pal);
        if let Some(doc) = wanted_doc {
            self.add_doc(doc);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn demo_outcome() -> AssemblyWin {
        let mut w = AssemblyWin::default();
        w.load_demo();
        w.run();
        w
    }

    #[test]
    fn the_demo_reads_assemble_into_one_circular_contig() {
        let w = demo_outcome();
        let o = w.outcome.as_ref().unwrap();
        assert_eq!(o.asm.contigs.len(), 1);
        assert!(o.asm.contigs[0].circular);
        assert_eq!(o.asm.unassembled.len(), 1);
        assert!(w.message.starts_with("1 contig from"), "{}", w.message);
    }

    #[test]
    fn the_consensus_document_has_one_feature_per_read() {
        let w = demo_outcome();
        let o = w.outcome.as_ref().unwrap();
        let c = &o.asm.contigs[0];
        let doc = contig_to_doc(c, &o.inputs);
        assert_eq!(doc.seq, c.consensus);
        assert!(doc.circular);
        assert_eq!(doc.features.len(), c.placed.len());
        assert!(doc
            .features
            .iter()
            .all(|f| f.spans.iter().all(|&(s, e)| s < e && e <= doc.len())));
        // each read's feature, read off the consensus, equals the read (allowing the odd error)
        for pl in &c.placed {
            let (name, seq) = &o.inputs[pl.read];
            let f = doc.features.iter().find(|f| &f.name == name).unwrap();
            // reading the feature on its own strand gives back the read as it was supplied
            let bases = doc.feature_bases(f);
            assert_eq!(bases.len(), seq.len(), "{name}");
            let diff = bases.iter().zip(seq).filter(|(a, b)| a != b).count();
            assert!(diff <= 2, "{name} differs in {diff} places");
        }
    }

    #[test]
    fn presets_switch_with_the_parameters() {
        let mut w = AssemblyWin::default();
        w.set_preset(Preset::Gibson);
        assert_eq!(
            (
                w.params.min_overlap,
                w.params.max_overlap,
                w.params.only_circular
            ),
            (12, Some(100), true)
        );
        w.params.min_overlap = 15;
        w.sync_preset();
        assert_eq!(w.preset, Preset::Custom);
        w.set_preset(Preset::Reads);
        assert_eq!(w.params.min_overlap, 20);
    }

    #[test]
    fn running_with_too_few_sequences_explains_itself() {
        let mut w = AssemblyWin::default();
        w.add_reads([("only".to_string(), b"ACGT".repeat(20))]);
        w.params.circularize = false;
        w.run();
        assert!(w.outcome.is_none());
        assert!(w.message.contains("at least two"));
    }

    #[test]
    fn gibson_style_fragments_make_a_circle_with_all_fragments_as_features() {
        let g = crate::sample::demo_plasmid().seq;
        let n = g.len();
        let frag = |a: usize, b: usize| -> Vec<u8> { (a..b).map(|k| g[k % n]).collect() };
        let mut w = AssemblyWin::default();
        w.add_reads([
            ("f1".to_string(), frag(0, 1000)),
            ("f2".to_string(), crate::seq::dna::revcom(&frag(960, 1900))),
            ("f3".to_string(), frag(1860, n + 40)),
        ]);
        w.set_preset(Preset::Gibson);
        w.run();
        let o = w.outcome.as_ref().unwrap();
        assert_eq!(o.asm.contigs.len(), 1, "{}", w.message);
        let c = &o.asm.contigs[0];
        assert!(c.circular);
        assert_eq!(c.consensus.len(), n);
        assert_eq!(c.placed.len(), 3);
    }
}
