//! ORFs: next / previous ORF, the ORF list with a six-frame map, and the translation window.

use crate::seq::dna::{find_orfs, Orf, StartRule, Strand, GENETIC_CODES};
use crate::seq::features::Feature;
use crate::seq::nextorf::{find_next, NoStopCodon, Settings};
use crate::seq::translation::{build, DnaLine, Letters, Numbers, Options, Source};
use crate::ui::app::App;
use crate::ui::map::commify;
use crate::ui::theme::{self, Palette};
use eframe::egui::{
    self, pos2, vec2, Align2, Color32, FontId, RichText, ScrollArea, Sense, Shape, Stroke, Ui,
};
use std::borrow::Cow;

/// The minimum-length choices of the original's ORF menu.
pub const ORF_MIN_CHOICES: [usize; 10] = [0, 30, 60, 90, 180, 300, 600, 900, 1200, 1500];
/// Most residues per line the width box takes.
const MAX_TRANSLATION_WIDTH: usize = 1000;

/// What the Translate window translates: the original's "Translate:" menu.
#[derive(Clone, Copy, PartialEq, Eq, Default, Debug)]
pub enum TranslateSource {
    #[default]
    All,
    Selection,
    /// A CDS feature, by its place in the document's feature list.
    Feature(usize),
}

/// The stretch the ORF list was opened on, taken once when it was opened (as the original's `orf_map` takes its `region`): a click
/// on a row selects that ORF, and the list must not follow the selection its own clicks make.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OrfRegion {
    /// The sequence it is of (`DocView::id`), the first base from 0 and one past the last.
    pub doc: u64,
    pub from: usize,
    pub to: usize,
}

#[derive(Default)]
pub struct OrfWin {
    pub list_open: bool,
    pub region: Option<OrfRegion>,
    pub translate_open: bool,
    pub translate: Options,
    pub source: TranslateSource,
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum OrfStep {
    Next,
    Prev,
}

pub fn rule_label(r: StartRule) -> &'static str {
    match r {
        StartRule::Met => "Met (ATG)",
        StartRule::AfterStop => "After a stop codon",
    }
}

impl App {
    /// Open the Translate window the way the original's dialog starts: on the selection when there is
    /// one, otherwise on the whole sequence.
    pub fn open_translate(&mut self) {
        let selected = self
            .current_doc()
            .is_some_and(|i| self.docs[i].selection().is_some());
        self.orf_win.translate_open = true;
        self.orf_win.source = if selected {
            TranslateSource::Selection
        } else {
            TranslateSource::All
        };
    }

    /// Open the ORF list on the selected bases when the Selection Only switch is on and some are selected (`orf_map` takes its
    /// region as "all" when nothing is, L36205), on all of the sequence otherwise.
    pub fn open_orf_list(&mut self) {
        let region = self.current_doc().and_then(|i| {
            let dv = &self.docs[i];
            let selected = dv.selection().filter(|_| self.prefs.selection_only);
            selected.map(|(from, to)| OrfRegion {
                doc: dv.id,
                from,
                to,
            })
        });
        self.orf_win.region = region;
        self.orf_win.list_open = true;
    }

    /// The stretch of sequence `i` that the list was opened on, if it was opened on one and the sequence is the same and still
    /// holds it.
    fn orf_region(&self, i: usize) -> Option<(usize, usize)> {
        let dv = &self.docs[i];
        let region = self.orf_win.region.filter(|r| r.doc == dv.id)?;
        let n = dv.doc.len();
        let (a, b) = (region.from.min(n), region.to.min(n));
        (a < b).then_some((a, b))
    }

    /// The ORFs of sequence `i`, numbered as in all of it. A selected stretch is a linear piece of its own: `orf_map` reads its
    /// text alone (L36206-36219), so an ORF that runs past its ends is not one.
    fn orfs_of(&self, i: usize, nested: bool) -> Vec<Orf> {
        let d = &self.docs[i].doc;
        let (seq, circular, from) = match self.orf_region(i) {
            Some((a, b)) => (d.slice(a, b), false, a),
            None => (&d.seq[..], d.circular, 0),
        };
        let mut orfs = find_orfs(
            seq,
            circular,
            self.prefs.orf_strand,
            self.prefs.orf_min,
            self.prefs.orf_rule,
            nested,
        );
        for o in &mut orfs {
            o.start += from;
        }
        orfs
    }

    pub fn select_orf(&mut self, i: usize, o: &Orf) {
        let n = self.docs[i].doc.len();
        let end = if o.start + o.length > n {
            n
        } else {
            o.start + o.length
        };
        let dv = &mut self.docs[i];
        dv.select(o.start, end);
        dv.scroll_to = Some(o.start);
    }

    /// Select the next / previous ORF exactly as the original's Find Next / Find Previous does
    /// (`next_ORF`, AppMain.tcl line 15025, offset 0xa9b56): every start is an ORF, the list is taken in order
    /// of start, and the search begins at the selection (one base on) or at the caret.
    pub fn orf_step(&mut self, dir: OrfStep) {
        let Some(i) = self.current_doc() else { return };
        let dv = &self.docs[i];
        let settings = Settings {
            circular: dv.doc.circular,
            strand: self.prefs.orf_strand,
            min: self.prefs.orf_min,
            rule: self.prefs.orf_rule,
        };
        let found = find_next(
            &dv.doc.seq,
            &settings,
            dv.caret,
            dv.selection().map(|(a, _)| a),
            dir == OrfStep::Next,
        );
        match found {
            Ok(Some(f)) => {
                let (a, b) = (f.select.0, f.select.1.min(dv.doc.len()));
                let dv = &mut self.docs[i];
                dv.select(a, b);
                dv.scroll_to = Some(a);
                self.say(format!(
                    "ORF {}–{} ({} bp, {} aa){}",
                    a + 1,
                    b,
                    b - a,
                    (b - a) / 3 - 1,
                    if f.entry.reverse {
                        ", reverse strand"
                    } else {
                        ""
                    }
                ));
            }
            Ok(None) => self.say("No ORFs with these settings"),
            Err(NoStopCodon) => {
                self.say("This sequence has no stop codon, so there are no ORFs to step through")
            }
        }
    }

    pub fn orf_window(&mut self, ctx: &egui::Context) {
        if self.orf_win.list_open {
            let Some(i) = self.current_doc() else {
                self.orf_win.list_open = false;
                return;
            };
            let pal = self.palette();
            let mut open = true;
            crate::ui::design::window("ORFs")
                .open(&mut open)
                .default_size([620.0, 430.0])
                .show(ctx, |ui| {
                    self.orf_controls(ui);
                    let orfs = self.orfs_of(i, false);
                    let n = self.docs[i].doc.len();
                    let region = self.orf_region(i);
                    ui.label(RichText::new(orf_heading(orfs.len(), region, n)).small());
                    let mut pick = draw_orf_map(ui, &orfs, n, region.unwrap_or((0, n)), &pal);
                    ui.separator();
                    ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            egui::Grid::new("orf_grid")
                                .striped(true)
                                .spacing([14.0, 2.0])
                                .show(ui, |ui| {
                                    for h in
                                        ["Strand", "Start", "End", "Length", "Amino acids", "Frame"]
                                    {
                                        ui.label(RichText::new(h).strong());
                                    }
                                    ui.end_row();
                                    for (k, o) in orfs.iter().enumerate() {
                                        let label = if o.reverse { "reverse" } else { "forward" };
                                        if ui.selectable_label(false, label).clicked() {
                                            pick = Some(k);
                                        }
                                        ui.label(RichText::new(commify(o.start + 1)).monospace());
                                        ui.label(
                                            RichText::new(commify(if o.end(n) == 0 {
                                                n
                                            } else {
                                                o.end(n)
                                            }))
                                            .monospace(),
                                        );
                                        ui.label(RichText::new(length_text(o)).monospace());
                                        ui.label(
                                            RichText::new(amino_acids(o).to_string()).monospace(),
                                        );
                                        ui.label(format!(
                                            "{}{}",
                                            if o.reverse { "−" } else { "+" },
                                            o.frame + 1
                                        ));
                                        ui.end_row();
                                    }
                                });
                        });
                    if let Some(k) = pick {
                        self.select_orf(i, &orfs[k]);
                    }
                    if ui.button("Add these ORFs as features").clicked() {
                        let added = orfs.len();
                        let feats: Vec<Feature> = orfs
                            .iter()
                            .enumerate()
                            .map(|(k, o)| orf_feature(o, n, k + 1))
                            .collect();
                        self.docs[i].edit(|d| d.features.extend(feats));
                        self.say(format!("Added {added} ORF features"));
                    }
                });
            self.orf_win.list_open = open;
        }
        self.translate_window(ctx);
    }

    fn orf_controls(&mut self, ui: &mut Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.label("Starts with");
            egui::ComboBox::from_id_salt("orf_rule")
                .selected_text(rule_label(self.prefs.orf_rule))
                .show_ui(ui, |ui| {
                    for r in [StartRule::Met, StartRule::AfterStop] {
                        ui.selectable_value(&mut self.prefs.orf_rule, r, rule_label(r));
                    }
                });
            ui.label("Strand");
            for (s, name) in [
                (Strand::Top, "Top"),
                (Strand::Bottom, "Bottom"),
                (Strand::Both, "Both"),
            ] {
                ui.selectable_value(&mut self.prefs.orf_strand, s, name);
            }
            ui.label("Minimum");
            egui::ComboBox::from_id_salt("orf_min")
                .selected_text(format!("{} bp", self.prefs.orf_min))
                .show_ui(ui, |ui| {
                    for m in ORF_MIN_CHOICES {
                        ui.selectable_value(&mut self.prefs.orf_min, m, m.to_string());
                    }
                });
        });
    }

    // ponytail: the text is rebuilt on every frame the window is open, which is instant for plasmids;
    // cache it by (sequence, options) if a sequence of hundreds of kilobases ever feels slow.
    fn translate_window(&mut self, ctx: &egui::Context) {
        if !self.orf_win.translate_open {
            return;
        }
        let Some(i) = self.current_doc() else {
            self.orf_win.translate_open = false;
            return;
        };
        let mut open = true;
        crate::ui::design::window("Translate")
            .open(&mut open)
            .default_size([600.0, 460.0])
            .show(ctx, |ui| {
                let dv = &self.docs[i];
                let selection = dv.selection();
                let cds: Vec<(usize, &Feature)> = dv
                    .doc
                    .features
                    .iter()
                    .enumerate()
                    .filter(|(_, f)| f.kind.eq_ignore_ascii_case("CDS"))
                    .collect();
                let w = &mut self.orf_win;
                // a source that is gone (selection cleared, feature deleted) falls back to everything
                let gone = match w.source {
                    TranslateSource::All => false,
                    TranslateSource::Selection => selection.is_none(),
                    TranslateSource::Feature(k) => !cds.iter().any(|&(j, _)| j == k),
                };
                if gone {
                    w.source = TranslateSource::All;
                }
                ui.horizontal(|ui| {
                    ui.label("Translate:");
                    let shown = match w.source {
                        TranslateSource::All => "All".to_string(),
                        TranslateSource::Selection => "Selection".to_string(),
                        TranslateSource::Feature(k) => dv.doc.features[k].name.clone(),
                    };
                    egui::ComboBox::from_id_salt("translate_source")
                        .selected_text(shown)
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut w.source, TranslateSource::All, "All");
                            if selection.is_some() {
                                ui.selectable_value(
                                    &mut w.source,
                                    TranslateSource::Selection,
                                    "Selection",
                                );
                            }
                            for &(k, f) in &cds {
                                ui.selectable_value(
                                    &mut w.source,
                                    TranslateSource::Feature(k),
                                    f.name.as_str(),
                                );
                            }
                        });
                });
                let from_feature = matches!(w.source, TranslateSource::Feature(_));
                let o = &mut w.translate;
                ui.horizontal_wrapped(|ui| {
                    ui.label("Code:");
                    ui.radio_value(&mut o.letters, Letters::One, "1 Letter");
                    ui.radio_value(&mut o.letters, Letters::Three, "3 Letter");
                    ui.checkbox(&mut o.codon_spacing, "Codon Spacing");
                    ui.add_enabled(
                        !from_feature,
                        egui::Checkbox::new(&mut o.reverse_complement, "Reverse Complement"),
                    );
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Line Width:");
                    ui.add(egui::DragValue::new(&mut o.width).range(0..=MAX_TRANSLATION_WIDTH))
                        .on_hover_text("Residues per line; 0 for one long line");
                    ui.label("Line Numbers:");
                    for (n, label) in [
                        (Numbers::None, "None"),
                        (Numbers::Left, "Left"),
                        (Numbers::Right, "Right"),
                        (Numbers::Both, "Both"),
                    ] {
                        ui.radio_value(&mut o.numbers, n, label);
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("DNA:");
                    for (d, label) in [
                        (DnaLine::None, "None"),
                        (DnaLine::Above, "Above"),
                        (DnaLine::Below, "Below"),
                    ] {
                        ui.radio_value(&mut o.dna_line, d, label);
                    }
                    ui.add_space(16.0);
                    ui.label("Genetic code");
                    egui::ComboBox::from_id_salt("code")
                        .selected_text(format!(
                            "{} ({})",
                            GENETIC_CODES[self.prefs.code].name, GENETIC_CODES[self.prefs.code].id
                        ))
                        .show_ui(ui, |ui| {
                            for (k, c) in GENETIC_CODES.iter().enumerate() {
                                ui.selectable_value(
                                    &mut self.prefs.code,
                                    k,
                                    format!("{} ({})", c.name, c.id),
                                );
                            }
                        });
                });
                let (bases, source): (Cow<[u8]>, Source) = match (w.source, selection) {
                    (TranslateSource::Selection, Some((a, b))) => (
                        Cow::Borrowed(dv.doc.slice(a, b)),
                        Source::Range { first: a + 1 },
                    ),
                    (TranslateSource::Feature(k), _) => (
                        Cow::Owned(dv.doc.feature_bases(&dv.doc.features[k])),
                        Source::Feature {
                            name: &dv.doc.features[k].name,
                        },
                    ),
                    _ => (Cow::Borrowed(&dv.doc.seq[..]), Source::Range { first: 1 }),
                };
                let out = build(&bases, GENETIC_CODES[self.prefs.code].table, source, o);
                // the original lists the file's name and place above the range
                let mut text: Vec<String> = vec![dv.doc.name.clone()];
                text.extend(dv.doc.path.iter().map(|p| p.display().to_string()));
                text.extend(out.header);
                text.extend(out.lines);
                ui.separator();
                let row = ui.text_style_height(&egui::TextStyle::Monospace);
                ScrollArea::both()
                    .auto_shrink([false, false])
                    .max_height(ui.available_height() - 30.0)
                    .show_rows(ui, row, text.len(), |ui, rows| {
                        for line in &text[rows] {
                            ui.add(
                                egui::Label::new(RichText::new(line).monospace())
                                    .wrap_mode(egui::TextWrapMode::Extend),
                            );
                        }
                    });
                if ui.button("Copy").clicked() {
                    ui.ctx()
                        .copy_text(text.iter().map(|l| format!("{l}\n")).collect::<String>());
                }
            });
        self.orf_win.translate_open = open;
    }
}

fn orf_feature(o: &Orf, n: usize, number: usize) -> Feature {
    let mut f = Feature::new(
        format!("ORF {number}"),
        "CDS",
        o.start,
        (o.start + o.length).min(n),
        o.reverse,
    );
    if o.start + o.length > n {
        f.spans = vec![(0, o.start + o.length - n), (o.start, n)];
    }
    f
}

/// "3 ORFs in 2,724 bp", or "in bases 51–300" when the list was opened on a selected stretch.
fn orf_heading(count: usize, region: Option<(usize, usize)>, n: usize) -> String {
    let place = match region {
        Some((a, b)) => format!("bases {}–{}", commify(a + 1), commify(b)),
        None => format!("{} bp", commify(n)),
    };
    let plural = if count == 1 { "" } else { "s" };
    format!("{count} ORF{plural} in {place}")
}

/// The residues of an ORF: its codons less the stop codon, which an ORF that ran out of sequence before one has not got.
fn amino_acids(o: &Orf) -> usize {
    let codons = o.length / 3;
    if o.open_end {
        codons
    } else {
        codons.saturating_sub(1)
    }
}

/// The length of an ORF; a plus after it when it is cut off by the end of the sequence (or of the stretch) with no stop codon.
fn length_text(o: &Orf) -> String {
    format!("{}{}", commify(o.length), if o.open_end { "+" } else { "" })
}

/// Six rows (three frames on each strand) with every ORF as a bar; returns the clicked ORF. `span` is the stretch the map
/// covers, all of the sequence or the selected stretch the list was opened on.
fn draw_orf_map(
    ui: &mut Ui,
    orfs: &[Orf],
    n: usize,
    span: (usize, usize),
    pal: &Palette,
) -> Option<usize> {
    let row_h = 14.0;
    let (rect, resp) = ui.allocate_exact_size(
        vec2(ui.available_width(), row_h * 6.0 + 24.0),
        Sense::click(),
    );
    let p = ui.painter_at(rect);
    p.rect_filled(rect, 2.0, pal.text_bg);
    let (left, right) = (rect.left() + 34.0, rect.right() - 10.0);
    let (lo, hi) = span;
    let x_of = |pos: usize| {
        left + pos.saturating_sub(lo) as f32 / (hi - lo).max(1) as f32 * (right - left)
    };
    for r in 0..6 {
        let y = rect.top() + 18.0 + r as f32 * row_h;
        let label = format!("{}{}", if r < 3 { "+" } else { "−" }, r % 3 + 1);
        p.text(
            pos2(rect.left() + 6.0, y + row_h / 2.0),
            Align2::LEFT_CENTER,
            label,
            FontId::proportional(10.0),
            pal.dim,
        );
        p.line_segment(
            [pos2(left, y + row_h - 1.0), pos2(right, y + row_h - 1.0)],
            Stroke::new(0.5, pal.line),
        );
    }
    let hover = resp.hover_pos();
    let mut hit = None;
    for (k, o) in orfs.iter().enumerate() {
        let row = if o.reverse {
            3 + o.frame as usize
        } else {
            o.frame as usize
        };
        let y = rect.top() + 18.0 + row as f32 * row_h;
        let fill = if o.reverse {
            Color32::from_rgb(235, 96, 108)
        } else {
            Color32::from_rgb(233, 208, 36)
        };
        let segments: Vec<(usize, usize)> = if o.start + o.length > n {
            vec![(o.start, n), (0, o.start + o.length - n)]
        } else {
            vec![(o.start, o.start + o.length)]
        };
        for (s, e) in segments {
            let bar = egui::Rect::from_min_max(
                pos2(x_of(s), y + 2.0),
                pos2(x_of(e).max(x_of(s) + 2.0), y + row_h - 3.0),
            );
            let head = 4.0f32.min(bar.width());
            let pts = if o.reverse {
                vec![
                    pos2(bar.left() + head, bar.top()),
                    pos2(bar.right(), bar.top()),
                    pos2(bar.right(), bar.bottom()),
                    pos2(bar.left() + head, bar.bottom()),
                    pos2(bar.left(), bar.center().y),
                ]
            } else {
                vec![
                    pos2(bar.left(), bar.top()),
                    pos2(bar.right() - head, bar.top()),
                    pos2(bar.right(), bar.center().y),
                    pos2(bar.right() - head, bar.bottom()),
                    pos2(bar.left(), bar.bottom()),
                ]
            };
            p.add(Shape::convex_polygon(
                pts,
                fill,
                Stroke::new(1.0, theme::mix(Color32::BLACK, fill, 30.0)),
            ));
            if hover.is_some_and(|h| bar.expand(1.0).contains(h)) {
                hit = Some(k);
            }
        }
    }
    p.text(
        pos2(left, rect.top() + 3.0),
        Align2::LEFT_TOP,
        commify(lo + 1),
        FontId::proportional(10.0),
        pal.dim,
    );
    p.text(
        pos2(right, rect.top() + 3.0),
        Align2::RIGHT_TOP,
        commify(hi),
        FontId::proportional(10.0),
        pal.dim,
    );
    if let Some(k) = hit {
        let o = &orfs[k];
        let tip = format!(
            "{}–{} • {} bp • {} aa • frame {}{}",
            commify(o.start + 1),
            commify(o.end(n).max(1)),
            length_text(o),
            amino_acids(o),
            if o.reverse { "−" } else { "+" },
            o.frame + 1
        );
        let clicked = resp.clicked();
        resp.on_hover_text(tip);
        if clicked {
            return Some(k);
        }
    }
    None
}

#[cfg(test)]
mod tests;
