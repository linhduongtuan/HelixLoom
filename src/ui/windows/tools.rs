//! Tools: Align Sequences, PCR products and primer design, a Tm / GC calculator, About.

use crate::seq::align::{align, Alignment, Scoring};
use crate::seq::dna::{extinction, gc_percent, melting_temp, revcom, translate, GENETIC_CODES};
use crate::seq::doc::Doc;
use crate::seq::features::Feature;
use crate::seq::pcr::{pcr_products, PcrPrimer, PcrScan, Product};
use crate::seq::primers::{find_primers, self_scores, Primer};
use crate::ui::app::App;
use crate::ui::map::commify;
use eframe::egui::{self, text::LayoutJob, Color32, FontId, RichText, ScrollArea, TextFormat, Ui};

/// Alignment columns per line (align_linewidth).
const ALIGN_LINE_WIDTH: usize = 100;
/// Longest list of primer candidates shown.
const MAX_CANDIDATES: usize = 40;

// ---------------------------------------------------------------------------
// Align
// ---------------------------------------------------------------------------

pub struct AlignOut {
    pub aln: Alignment,
    pub name_a: String,
    pub name_b: String,
    pub offset_a: usize,
    pub offset_b: usize,
}

pub struct AlignWin {
    pub open: bool,
    pub a: usize,
    pub b: usize,
    pub sel_a: bool,
    pub sel_b: bool,
    pub reverse_b: bool,
    pub free_ends: bool,
    pub scoring: Scoring,
    pub out: Option<AlignOut>,
    pub message: String,
}

impl Default for AlignWin {
    fn default() -> Self {
        AlignWin {
            open: false,
            a: 0,
            b: 1,
            sel_a: false,
            sel_b: false,
            reverse_b: false,
            free_ends: true,
            scoring: Scoring::default(),
            out: None,
            message: String::new(),
        }
    }
}

pub struct Block {
    pub a_start: usize,
    pub b_start: usize,
    pub a: String,
    pub mid: String,
    pub b: String,
}

/// Split an alignment into display blocks with the position (1-based) where each line starts.
pub fn alignment_blocks(
    aln: &Alignment,
    width: usize,
    offset_a: usize,
    offset_b: usize,
) -> Vec<Block> {
    let (mut pa, mut pb) = (offset_a, offset_b);
    let mut blocks = Vec::new();
    for (ca, cb) in aln.a.chunks(width).zip(aln.b.chunks(width)) {
        let mid: String = ca
            .iter()
            .zip(cb)
            .map(|(x, y)| {
                if *x != b'-' && *y != b'-' && x.eq_ignore_ascii_case(y) {
                    '|'
                } else {
                    ' '
                }
            })
            .collect();
        blocks.push(Block {
            a_start: pa + 1,
            b_start: pb + 1,
            a: String::from_utf8_lossy(ca).into_owned(),
            mid,
            b: String::from_utf8_lossy(cb).into_owned(),
        });
        pa += ca.iter().filter(|&&c| c != b'-').count();
        pb += cb.iter().filter(|&&c| c != b'-').count();
    }
    blocks
}

/// One alignment line: the sequence's first letter, the position it starts at, then the bases.
fn alignment_line(ui: &mut Ui, name: &str, start: usize, job: LayoutJob) {
    ui.horizontal(|ui| {
        ui.label(
            RichText::new(format!("{} {start:>7}", name.chars().next().unwrap_or('?')))
                .monospace()
                .weak(),
        );
        ui.label(job);
    });
}

fn colored(text: &str, mid: &str, plain: Color32) -> LayoutJob {
    let mut job = LayoutJob::default();
    let font = FontId::monospace(13.0);
    let mut run = String::new();
    let mut run_color = plain;
    for (c, m) in text.chars().zip(mid.chars().chain(std::iter::repeat(' '))) {
        let color = if c == '-' {
            Color32::from_gray(150)
        } else if m == '|' {
            plain
        } else {
            Color32::from_rgb(214, 48, 48)
        };
        if color != run_color && !run.is_empty() {
            job.append(
                &run,
                0.0,
                TextFormat {
                    font_id: font.clone(),
                    color: run_color,
                    ..Default::default()
                },
            );
            run.clear();
        }
        run_color = color;
        run.push(c);
    }
    job.append(
        &run,
        0.0,
        TextFormat {
            font_id: font,
            color: run_color,
            ..Default::default()
        },
    );
    job
}

impl App {
    fn run_align(&mut self) {
        let w = &self.align_win;
        let (Some(da), Some(db)) = (self.docs.get(w.a), self.docs.get(w.b)) else {
            self.align_win.message = "Open two sequences to align.".into();
            return;
        };
        let part = |d: &crate::ui::state::DocView, sel: bool| match d.selection().filter(|_| sel) {
            Some((a, b)) => (d.doc.slice(a, b).to_vec(), a),
            None => (d.doc.seq.clone(), 0),
        };
        let (sa, offset_a) = part(da, w.sel_a);
        let (mut sb, offset_b) = part(db, w.sel_b);
        if w.reverse_b {
            sb = revcom(&sb);
        }
        let names = (
            da.doc.name.clone(),
            format!(
                "{}{}",
                db.doc.name,
                if w.reverse_b {
                    " (reverse complement)"
                } else {
                    ""
                }
            ),
        );
        match align(&sa, &sb, &w.scoring, w.free_ends) {
            Ok(aln) => {
                self.align_win.message = format!(
                    "Score {:.1} • {:.1}% identical over {} columns",
                    aln.score,
                    aln.identity() * 100.0,
                    commify(aln.columns)
                );
                self.align_win.out = Some(AlignOut {
                    aln,
                    name_a: names.0,
                    name_b: names.1,
                    offset_a,
                    offset_b,
                });
            }
            Err(e) => {
                self.align_win.message = e;
                self.align_win.out = None;
            }
        }
    }

    pub fn align_window(&mut self, ctx: &egui::Context) {
        if !self.align_win.open {
            return;
        }
        let mut open = true;
        let mut run = false;
        crate::ui::design::window("Align Sequences")
            .open(&mut open)
            .default_size([860.0, 480.0])
            .show(ctx, |ui| {
                let names: Vec<String> = self.docs.iter().map(|d| d.doc.name.clone()).collect();
                let last = names.len().saturating_sub(1);
                self.align_win.a = self.align_win.a.min(last);
                self.align_win.b = self.align_win.b.min(last);
                let w = &mut self.align_win;
                egui::Grid::new("align_pick")
                    .num_columns(3)
                    .spacing([10.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("First");
                        egui::ComboBox::from_id_salt("align_first")
                            .selected_text(names.get(w.a).cloned().unwrap_or_default())
                            .show_ui(ui, |ui| {
                                for (k, n) in names.iter().enumerate() {
                                    ui.selectable_value(&mut w.a, k, n);
                                }
                            });
                        ui.checkbox(&mut w.sel_a, "selection only");
                        ui.end_row();
                        ui.label("Second");
                        egui::ComboBox::from_id_salt("align_second")
                            .selected_text(names.get(w.b).cloned().unwrap_or_default())
                            .show_ui(ui, |ui| {
                                for (k, n) in names.iter().enumerate() {
                                    ui.selectable_value(&mut w.b, k, n);
                                }
                            });
                        ui.checkbox(&mut w.sel_b, "selection only");
                        ui.end_row();
                    });
                ui.horizontal_wrapped(|ui| {
                    ui.checkbox(&mut w.reverse_b, "Use the reverse complement of the second");
                    ui.checkbox(&mut w.free_ends, "Ignore overhanging ends");
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Scores: mismatch");
                    ui.add(
                        egui::DragValue::new(&mut w.scoring.mismatch)
                            .speed(0.05)
                            .range(-5.0..=0.0),
                    );
                    ui.label("gap");
                    ui.add(
                        egui::DragValue::new(&mut w.scoring.gap_open)
                            .speed(0.05)
                            .range(-20.0..=0.0),
                    );
                    ui.label("further gap");
                    ui.add(
                        egui::DragValue::new(&mut w.scoring.gap_extend)
                            .speed(0.05)
                            .range(-10.0..=0.0),
                    );
                    if ui.button("Original defaults").clicked() {
                        w.scoring = Scoring::default();
                    }
                });
                ui.horizontal(|ui| {
                    run = ui.button(RichText::new("  Align  ").strong()).clicked();
                    ui.label(RichText::new(&w.message).small());
                });
                ui.separator();
                if let Some(out) = &w.out {
                    ScrollArea::both()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            let plain = ui.visuals().text_color();
                            for b in alignment_blocks(
                                &out.aln,
                                ALIGN_LINE_WIDTH,
                                out.offset_a,
                                out.offset_b,
                            ) {
                                alignment_line(
                                    ui,
                                    &out.name_a,
                                    b.a_start,
                                    colored(&b.a, &b.mid, plain),
                                );
                                ui.horizontal(|ui| {
                                    ui.label(RichText::new("          ").monospace());
                                    ui.label(
                                        RichText::new(&b.mid)
                                            .monospace()
                                            .color(Color32::from_rgb(70, 140, 70)),
                                    );
                                });
                                alignment_line(
                                    ui,
                                    &out.name_b,
                                    b.b_start,
                                    colored(&b.b, &b.mid, plain),
                                );
                                ui.add_space(6.0);
                            }
                        });
                }
            });
        if run {
            self.run_align();
        }
        self.align_win.open = open;
    }
}

// ---------------------------------------------------------------------------
// PCR and primers
// ---------------------------------------------------------------------------

/// A product and the template's features that it copies (`pcr_dialog_do_pcr` takes them with `get_tags2`).
pub struct PcrResult {
    pub product: Product,
    pub carried: Vec<Feature>,
}

#[derive(Default)]
pub struct PcrWin {
    pub open: bool,
    pub fwd: String,
    pub rev: String,
    /// How a primer is matched to the template (the original's PCR dialog settings).
    pub scan: PcrScan,
    pub products: Vec<PcrResult>,
    pub message: String,
    pub fwd_cands: Vec<Primer>,
    pub rev_cands: Vec<Primer>,
}

fn clean(s: &str) -> Vec<u8> {
    s.bytes()
        .filter(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase())
        .collect()
}

/// A short report on a typed primer: length, Tm, GC and self-complementarity.
pub fn primer_summary(seq: &[u8], dna_nm: f64, salt_mm: f64) -> String {
    if seq.is_empty() {
        return String::new();
    }
    let tm = melting_temp(seq, dna_nm, salt_mm).map_or("---".to_string(), |t| format!("{t:.1}°C"));
    let (end3, total) = self_scores(seq);
    format!(
        "{} bases • Tm {tm} • GC {:.0}% • self-complementarity {total} (3' end {end3})",
        seq.len(),
        gc_percent(seq)
    )
}

impl App {
    pub(crate) fn find_products(&mut self) {
        let Some(i) = self.current_doc() else { return };
        let primer = |name: &str, text: &str| PcrPrimer {
            name: name.to_string(),
            seq: clean(text),
        };
        let (fwd, rev) = (
            primer("Forward", &self.pcr_win.fwd),
            primer("Reverse", &self.pcr_win.rev),
        );
        let d = &self.docs[i].doc;
        let found = pcr_products(&d.seq, d.circular, &d.name, &fwd, &rev, &self.pcr_win.scan);
        self.pcr_win.products = found
            .into_iter()
            .map(|product| {
                let carried = product
                    .template_span
                    .map(|(a, b)| d.features_in(a, b, false))
                    .unwrap_or_default()
                    .into_iter()
                    .map(|f| f.shifted(product.template_offset))
                    .collect();
                PcrResult { product, carried }
            })
            .collect();
        let (n, scan) = (self.pcr_win.products.len(), self.pcr_win.scan);
        self.pcr_win.message = match n {
            0 => format!(
                "No product: each primer must match the template from its 3' end for at least {} bases, with up to {} mismatches.",
                scan.min_match, scan.max_mismatch
            ),
            1 => "1 product".to_string(),
            n => format!("{n} products"),
        };
    }

    /// Mark where the typed primers bind as primer_bind features.
    fn mark_primers(&mut self, p: &Product) {
        let Some(i) = self.current_doc() else { return };
        let n = self.docs[i].doc.len();
        let marked = |name: &str, (a, b): (usize, usize), reverse: bool| {
            let mut f = Feature::new(name, "primer_bind", a, b.min(n), reverse);
            if a >= b {
                // an annealed part that runs through the origin of a circle
                f.spans = vec![(a, n), (0, b)];
            }
            f
        };
        let feats = [
            marked("PCR forward", p.fwd, false),
            marked("PCR reverse", p.rev, true),
        ];
        self.docs[i].edit(|d| d.features.extend(feats));
    }

    /// The new document the original's PCR makes: the product with its name, its comment, the template's
    /// features that it copies and the two primer features.
    pub(crate) fn open_product(&mut self, r: &PcrResult) {
        let p = &r.product;
        let mut d = Doc::new(p.name.clone(), p.seq.clone());
        d.methylated = false;
        d.comment = p.comment.join("\n");
        d.features = r.carried.iter().chain(&p.features).cloned().collect();
        d.dirty = true;
        self.add_doc(d);
    }

    fn search_primers(&mut self) {
        let Some(i) = self.current_doc() else { return };
        let dv = &self.docs[i];
        let Some((a, b)) = dv.selection() else {
            self.pcr_win.message =
                "Select a region of the sequence first; primers are searched inside it.".into();
            return;
        };
        let (fwd, rev) = (
            find_primers(&dv.doc.seq, a, b, false, &self.prefs.primer),
            find_primers(&dv.doc.seq, a, b, true, &self.prefs.primer),
        );
        self.pcr_win.message = format!(
            "{} forward and {} reverse candidates inside the selection.",
            fwd.len(),
            rev.len()
        );
        self.pcr_win.fwd_cands = fwd;
        self.pcr_win.rev_cands = rev;
    }

    pub fn pcr_window(&mut self, ctx: &egui::Context) {
        if !self.pcr_win.open {
            return;
        }
        let mut open = true;
        let (mut find, mut search) = (false, false);
        let mut open_this: Option<usize> = None;
        let mut mark: Option<Product> = None;
        let mut select: Option<(usize, usize)> = None;
        crate::ui::design::window("PCR and Primers")
            .open(&mut open)
            .default_size([700.0, 560.0])
            .show(ctx, |ui| {
                ui.label(RichText::new("PCR product").strong());
                ui.label("Forward primer (5' to 3')");
                ui.add(
                    egui::TextEdit::singleline(&mut self.pcr_win.fwd)
                        .font(egui::TextStyle::Monospace)
                        .desired_width(f32::INFINITY),
                );
                ui.label(
                    RichText::new(primer_summary(
                        &clean(&self.pcr_win.fwd),
                        self.prefs.dna_nm,
                        self.prefs.salt_mm,
                    ))
                    .small(),
                );
                ui.label("Reverse primer (5' to 3')");
                ui.add(
                    egui::TextEdit::singleline(&mut self.pcr_win.rev)
                        .font(egui::TextStyle::Monospace)
                        .desired_width(f32::INFINITY),
                );
                ui.label(
                    RichText::new(primer_summary(
                        &clean(&self.pcr_win.rev),
                        self.prefs.dna_nm,
                        self.prefs.salt_mm,
                    ))
                    .small(),
                );
                ui.horizontal_wrapped(|ui| {
                    let scan = &mut self.pcr_win.scan;
                    ui.label("3' match at least");
                    ui.add(
                        egui::DragValue::new(&mut scan.min_match)
                            .range(5..=30)
                            .suffix(" bases"),
                    );
                    ui.label("mismatches allowed");
                    ui.add(egui::DragValue::new(&mut scan.max_mismatch).range(0..=10));
                    ui.add_enabled_ui(scan.max_mismatch > 0, |ui| {
                        ui.label("exact 3' bases");
                        ui.add(egui::DragValue::new(&mut scan.min_3_match).range(0..=4));
                    });
                    find = ui.button(RichText::new("Find products").strong()).clicked();
                });
                for (k, r) in self.pcr_win.products.iter().enumerate() {
                    let p = &r.product;
                    ui.horizontal(|ui| {
                        ui.label(
                            RichText::new(format!("{} bp", commify(p.length())))
                                .monospace()
                                .strong(),
                        );
                        ui.label(format!(
                            "forward at {}–{}, reverse at {}–{}",
                            commify(p.fwd.0 + 1),
                            commify(p.fwd.1),
                            commify(p.rev.0 + 1),
                            commify(p.rev.1)
                        ));
                        if ui.small_button("Open product").clicked() {
                            open_this = Some(k);
                        }
                        if ui.small_button("Mark primers").clicked() {
                            mark = Some(p.clone());
                        }
                    });
                    for warning in &p.warnings {
                        ui.label(
                            RichText::new(warning)
                                .small()
                                .color(Color32::from_rgb(190, 110, 0)),
                        );
                    }
                }
                ui.separator();
                ui.label(RichText::new("Find PCR primers inside the selection").strong());
                ui.collapsing("Limits (the original's defaults)", |ui| {
                    let p = &mut self.prefs.primer;
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Length");
                        ui.add(egui::DragValue::new(&mut p.len_min).range(10..=60));
                        ui.add(egui::DragValue::new(&mut p.len_max).range(10..=60));
                        ui.label("Tm");
                        ui.add(egui::DragValue::new(&mut p.tm_min).range(30.0..=90.0));
                        ui.add(egui::DragValue::new(&mut p.tm_max).range(30.0..=90.0));
                        ui.label("GC %");
                        ui.add(egui::DragValue::new(&mut p.gc_min).range(0.0..=100.0));
                        ui.add(egui::DragValue::new(&mut p.gc_max).range(0.0..=100.0));
                        ui.label("3' G/C clamp");
                        ui.add(egui::DragValue::new(&mut p.clamp_min).range(0..=5));
                        ui.add(egui::DragValue::new(&mut p.clamp_max).range(0..=5));
                        ui.label("Self-complementarity");
                        ui.add(egui::DragValue::new(&mut p.self_total_max).range(0..=40));
                        ui.add(egui::DragValue::new(&mut p.self_3prime_max).range(0..=40));
                    });
                    if ui.button("Restore defaults").clicked() {
                        *p = Default::default();
                    }
                });
                ui.horizontal(|ui| {
                    search = ui.button("Search selection").clicked();
                    ui.label(RichText::new(&self.pcr_win.message).small());
                });
                ScrollArea::vertical().max_height(260.0).show(ui, |ui| {
                    for (title, cands, is_fwd) in [
                        ("Forward", &self.pcr_win.fwd_cands, true),
                        ("Reverse", &self.pcr_win.rev_cands, false),
                    ] {
                        if cands.is_empty() {
                            continue;
                        }
                        ui.label(RichText::new(title).strong());
                        egui::Grid::new(("cands", is_fwd))
                            .striped(true)
                            .spacing([10.0, 2.0])
                            .show(ui, |ui| {
                                for h in ["Tm", "GC %", "3' GC", "Length", "Sequence", ""] {
                                    ui.label(RichText::new(h).small().strong());
                                }
                                ui.end_row();
                                for c in cands.iter().take(MAX_CANDIDATES) {
                                    ui.label(format!("{:.1}", c.tm));
                                    ui.label(format!("{:.0}", c.gc));
                                    ui.label(c.clamp.to_string());
                                    ui.label(c.seq.len().to_string());
                                    if ui
                                        .selectable_label(
                                            false,
                                            RichText::new(String::from_utf8_lossy(&c.seq))
                                                .monospace(),
                                        )
                                        .on_hover_text("Select where it binds")
                                        .clicked()
                                    {
                                        select = Some((c.start, c.end));
                                    }
                                    if ui.small_button("Use").clicked() {
                                        let text = String::from_utf8_lossy(&c.seq).into_owned();
                                        if is_fwd {
                                            self.pcr_win.fwd = text;
                                        } else {
                                            self.pcr_win.rev = text;
                                        }
                                    }
                                    ui.end_row();
                                }
                            });
                    }
                });
            });
        if find {
            self.find_products();
        }
        if search {
            self.search_primers();
        }
        if let Some(p) = mark {
            self.mark_primers(&p);
        }
        if let Some(r) = open_this.and_then(|k| self.pcr_win.products.get(k)) {
            let r = PcrResult {
                product: r.product.clone(),
                carried: r.carried.clone(),
            };
            self.open_product(&r);
        }
        if let (Some((a, b)), Some(i)) = (select, self.current_doc()) {
            self.docs[i].select(a, b);
            self.docs[i].scroll_to = Some(a);
        }
        self.pcr_win.open = open;
    }

    // -----------------------------------------------------------------------
    // Calculator and About
    // -----------------------------------------------------------------------

    pub fn calc_window(&mut self, ctx: &egui::Context) {
        if !self.calc_win.open {
            return;
        }
        let mut open = true;
        crate::ui::design::window("Tm and GC Calculator")
            .open(&mut open)
            .default_size([460.0, 360.0])
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label("Sequence");
                    if ui.small_button("Use selection").clicked() {
                        if let Some(i) = self.current_doc() {
                            self.calc_win.text =
                                String::from_utf8_lossy(self.docs[i].selected_or_all())
                                    .into_owned();
                        }
                    }
                });
                ui.add(
                    egui::TextEdit::multiline(&mut self.calc_win.text)
                        .font(egui::TextStyle::Monospace)
                        .desired_rows(4)
                        .desired_width(f32::INFINITY),
                );
                ui.horizontal(|ui| {
                    ui.label("DNA");
                    ui.add(
                        egui::DragValue::new(&mut self.prefs.dna_nm)
                            .range(1.0..=100_000.0)
                            .suffix(" nM"),
                    );
                    ui.label("Salt");
                    ui.add(
                        egui::DragValue::new(&mut self.prefs.salt_mm)
                            .range(1.0..=2000.0)
                            .suffix(" mM"),
                    );
                });
                let seq = clean(&self.calc_win.text);
                ui.separator();
                egui::Grid::new("calc")
                    .num_columns(2)
                    .spacing([14.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("Length");
                        ui.label(format!("{} bases", commify(seq.len())));
                        ui.end_row();
                        ui.label("GC content");
                        ui.label(format!("{:.1} %", gc_percent(&seq)));
                        ui.end_row();
                        ui.label("Melting temperature");
                        ui.label(
                            melting_temp(&seq, self.prefs.dna_nm, self.prefs.salt_mm).map_or(
                                "--- (only for sequences up to 100 bases)".to_string(),
                                |t| format!("{t:.1} °C"),
                            ),
                        );
                        ui.end_row();
                        ui.label("Extinction coefficient");
                        ui.label(
                            extinction(&seq).map_or("--- (A, C, G and T only)".to_string(), |e| {
                                format!("{e:.1} L/mmol/cm")
                            }),
                        );
                        ui.end_row();
                        ui.label("Reverse complement");
                        ui.label(
                            RichText::new(String::from_utf8_lossy(&revcom(&seq)).into_owned())
                                .monospace(),
                        );
                        ui.end_row();
                        ui.label("Translation");
                        ui.label(
                            RichText::new(
                                String::from_utf8_lossy(&translate(
                                    &seq,
                                    GENETIC_CODES[self.prefs.code].table,
                                ))
                                .into_owned(),
                            )
                            .monospace(),
                        );
                        ui.end_row();
                    });
            });
        self.calc_win.open = open;
    }

    pub fn about_window(&mut self, ctx: &egui::Context) {
        if !self.about_open {
            return;
        }
        let mut open = true;
        crate::ui::design::window("About HelixLoom").open(&mut open).collapsible(false).resizable(false).show(ctx, |ui| {
            ui.heading(format!("HelixLoom {}", env!("CARGO_PKG_VERSION")));
            ui.label("A bioinformatics workbench: plasmids, Sanger traces, alignment, phylogeny, trees, protein structure and pipelines, on this computer or on a workstation with GPUs.");
            ui.add_space(6.0);
            ui.label("It began as a re-creation of the behaviour of ApE (A plasmid Editor) by M. Wayne Davis; the interface is its own now, and none of ApE's code or artwork is used.");
            ui.label("Restriction enzyme names and recognition sequences are public REBASE data.");
            ui.label("Built with Rust, egui/eframe (MIT / Apache-2.0), png and rfd.");
            ui.add_space(6.0);
            ui.label(RichText::new("The contig assembler is new: the original program has no read assembler.").small());
        });
        self.about_open = open;
    }
}

#[derive(Default)]
pub struct CalcWin {
    pub open: bool,
    pub text: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alignment_blocks_track_positions_across_gaps() {
        let aln = align(b"ACGTACGTAC", b"ACGTCGTAC", &Scoring::default(), false).unwrap();
        let blocks = alignment_blocks(&aln, 4, 100, 0);
        assert_eq!(blocks[0].a_start, 101);
        assert_eq!(blocks[0].b_start, 1);
        assert_eq!(blocks.iter().map(|b| b.a.len()).sum::<usize>(), aln.a.len());
        assert_eq!(blocks[0].mid.chars().filter(|&c| c == '|').count(), 4);
        // the second block starts after the bases (not columns) already used
        let used_a = aln.a[..4].iter().filter(|&&c| c != b'-').count();
        assert_eq!(blocks[1].a_start, 100 + used_a + 1);
    }

    #[test]
    fn primer_summaries_use_the_original_tm() {
        let s = primer_summary(b"ATGGCGTTTGAACGGCTCAA", 250.0, 50.0);
        assert!(s.contains("Tm 59.3°C"), "{s}");
        assert!(s.contains("GC 50%"));
        assert!(s.contains("self-complementarity 7 (3' end 7)"));
        assert_eq!(primer_summary(b"", 250.0, 50.0), "");
    }

    #[test]
    fn pasted_primers_are_cleaned() {
        assert_eq!(clean(" 5'-atg gcc-3'\n"), b"ATGGCC".to_vec());
    }
}
