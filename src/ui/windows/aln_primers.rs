//! Alignment → Pick primers and probes: assays picked from the columns that the aligned sequences share (`seq::alnprimers`), shown on
//! the alignment, and put on the consensus as features.

use crate::seq::alnprimers::{self, Design, Report, Settings, Site};
use crate::seq::doc::Doc;
use crate::seq::features::Feature;
use crate::ui::app::{Action, App};
use crate::ui::dialog;
use crate::ui::platform::alignment::{Mark, Role};
use crate::ui::platform::Workspace;
use eframe::egui::{self, RichText};
use std::hash::{Hash, Hasher};

/// The product sizes that real-time PCR uses, as in the Primer3 window.
const TAQMAN_PRODUCT: (usize, usize) = (70, 150);

pub struct AlnPrimersWin {
    pub open: bool,
    pub settings: Settings,
    /// The share of the sequences that must agree in a column, in percent (the settings keep it as a share).
    pub percent: f64,
    pub use_region: bool,
    /// The region searched, in columns counted from 1, both ends included.
    pub region: (usize, usize),
    pub report: Option<Report>,
    /// The names of the rows the report is for, and a fingerprint of the alignment it was made from: what is put on the alignment or
    /// on a consensus is the report's only while the alignment is that one.
    names: Vec<String>,
    made_for: Option<u64>,
    /// The design whose details are shown.
    pub chosen: usize,
    pub message: String,
}

impl Default for AlnPrimersWin {
    fn default() -> Self {
        AlnPrimersWin {
            open: false,
            settings: Settings::default(),
            percent: 100.0,
            use_region: false,
            region: (1, 1),
            report: None,
            names: Vec::new(),
            made_for: None,
            chosen: 0,
            message: String::new(),
        }
    }
}

fn fingerprint(rows: &[(String, Vec<u8>)]) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    rows.hash(&mut h);
    h.finish()
}

fn marks_of(d: &Design) -> Vec<Mark> {
    let mark = |s: &Site, role| Mark {
        first: s.oligo.start,
        last: s.oligo.end - 1,
        role,
    };
    let mut marks = vec![mark(&d.left, Role::Left), mark(&d.right, Role::Right)];
    marks.extend(d.probe.as_ref().map(|p| mark(p, Role::Probe)));
    marks
}

/// The consensus without its gaps, with the primers as `primer_bind` features and the probe as a `misc_feature` with its sequence.
fn consensus_doc(report: &Report, d: &Design, n: usize) -> Doc {
    let (bases, at) = alnprimers::ungapped(&report.consensus);
    let mut doc = Doc::new("Alignment consensus", bases);
    doc.definition = format!(
        "Consensus of {} aligned sequences ({} columns), with the primers and probe of design {n}.",
        report.rows, report.columns
    );
    let span = |s: &Site| Some((at[s.oligo.start]?, at[s.oligo.end - 1]? + 1));
    let mut feature = |name: String, kind: &str, s: &Site, note: Option<String>| {
        if let Some((a, b)) = span(s) {
            let mut f = Feature::new(name, kind, a, b, s.oligo.reverse);
            f.qualifiers.extend(note.map(|n| ("note".to_string(), n)));
            doc.features.push(f);
        }
    };
    feature(format!("Left primer {n}"), "primer_bind", &d.left, None);
    feature(format!("Right primer {n}"), "primer_bind", &d.right, None);
    if let Some(p) = &d.probe {
        let note = format!(
            "TaqMan probe, 5'-{}-3'",
            String::from_utf8_lossy(&p.oligo.seq)
        );
        feature(format!("Probe {n}"), "misc_feature", p, Some(note));
    }
    doc.dirty = true;
    doc
}

fn product_text(d: &Design) -> String {
    match d.product {
        Some((a, b)) if a == b => a.to_string(),
        Some((a, b)) => format!("{a}–{b}"),
        None => "–".into(),
    }
}

#[derive(Default)]
pub(super) struct Asks {
    pick: bool,
    show: Option<usize>,
    consensus: Option<usize>,
    copy: bool,
    save: bool,
}

impl App {
    /// Pick designs for the alignment on the Alignment page.
    pub fn run_aln_primers(&mut self) {
        let w = &mut self.aln_win;
        let Some(a) = self.platform.alignment.result.as_ref() else {
            w.report = None;
            w.message =
                "There is no alignment: align sequences on the Alignment page first.".into();
            return;
        };
        let mut s = w.settings.clone();
        s.min_identity = w.percent / 100.0;
        s.window = w
            .use_region
            .then(|| (w.region.0.saturating_sub(1), w.region.1));
        match alnprimers::design(&a.rows, &s) {
            Ok(report) => {
                w.message = report_message(&report, &s, w.percent);
                w.names = a.rows.iter().map(|r| r.0.clone()).collect();
                w.made_for = Some(fingerprint(&a.rows));
                w.chosen = 0;
                w.report = Some(report);
            }
            Err(why) => {
                w.report = None;
                w.message = why;
            }
        }
    }

    /// The design `k` of the report, when the alignment is still the one it was made for.
    fn design_for_now(&mut self, k: usize) -> Option<Design> {
        let same = self
            .platform
            .alignment
            .result
            .as_ref()
            .map(|a| fingerprint(&a.rows))
            == self.aln_win.made_for;
        if !same {
            self.say("The alignment has changed since these were picked: pick them again.");
            return None;
        }
        self.aln_win.report.as_ref()?.designs.get(k).cloned()
    }

    pub(super) fn apply_aln_primer_asks(
        &mut self,
        ctx: &egui::Context,
        asks: Asks,
        out: &mut Vec<Action>,
    ) {
        if asks.pick {
            self.run_aln_primers();
        }
        if let Some(k) = asks.show {
            if let Some(d) = self.design_for_now(k) {
                self.aln_win.chosen = k;
                self.platform.alignment.set_marks(marks_of(&d));
                out.push(Action::Go(Workspace::Alignment));
            }
        }
        if let Some(k) = asks.consensus {
            if let Some(d) = self.design_for_now(k) {
                let doc = self
                    .aln_win
                    .report
                    .as_ref()
                    .map(|r| consensus_doc(r, &d, k + 1));
                if let Some(doc) = doc {
                    self.add_doc(doc);
                    out.push(Action::Go(Workspace::Sequences));
                    self.say("Opened the consensus with the primers and probe as features.");
                }
            }
        }
        let w = &self.aln_win;
        let table = w
            .report
            .as_ref()
            .map(|r| alnprimers::table_text(r, &w.names));
        if let (Some(text), true) = (&table, asks.copy) {
            ctx.copy_text(text.clone());
            self.aln_win.message = "Copied the table.".into();
        }
        if let (Some(text), true) = (table, asks.save) {
            let ask = dialog::FileDialog::new().set_title("Save the primers and probes");
            if let Some(saved) = dialog::save_bytes(ask, "primers.tsv", text.as_bytes()) {
                self.aln_win.message = saved.message;
            }
        }
    }

    pub fn aln_primers_window(&mut self, ctx: &egui::Context, out: &mut Vec<Action>) {
        if !self.aln_win.open {
            return;
        }
        let columns = self
            .platform
            .alignment
            .result
            .as_ref()
            .map_or(0, |a| a.columns());
        let (mut open, mut asks) = (true, Asks::default());
        let look = self.prefs.look();
        crate::ui::design::window("Primers and Probes from the Alignment")
            .open(&mut open)
            .default_width(760.0)
            .show(ctx, |ui| {
                settings_form(ui, &mut self.aln_win, columns, &mut asks);
                if !self.aln_win.message.is_empty() {
                    ui.label(&self.aln_win.message);
                }
                if let Some(report) = &self.aln_win.report {
                    results(ui, &look, &self.aln_win, report, &mut asks);
                }
            });
        self.aln_win.open = open;
        self.apply_aln_primer_asks(ctx, asks, out);
    }
}

fn report_message(r: &Report, s: &Settings, percent: f64) -> String {
    let all = format!("{} of {} columns are usable", r.usable, r.columns);
    if r.usable == 0 {
        return format!(
            "No column is shared by {percent:.0} % of the {} sequences: lower the share.",
            r.rows
        );
    }
    if r.designs.is_empty() {
        let need = s.primers.min_size;
        return if r.longest_run < need {
            format!(
                "{all}, but the longest stretch of them is {} columns, shorter than a primer ({need}): lower the share, or the primer length.",
                r.longest_run
            )
        } else {
            format!(
                "{all} (the longest stretch {}), but no pair fits the limits: widen the Tm, size or product range{}.",
                r.longest_run,
                if s.with_probe { ", or the probe's" } else { "" }
            )
        };
    }
    let best = &r.designs[0];
    let what = if s.with_probe {
        "assays (pair and probe)"
    } else {
        "pairs"
    };
    format!(
        "{} {what}, best first ({all}). The best covers {} of {} sequences.",
        r.designs.len(),
        best.tolerated,
        r.rows
    )
}

fn count_box(ui: &mut egui::Ui, value: &mut usize, range: std::ops::RangeInclusive<usize>) {
    ui.add(egui::DragValue::new(value).range(range));
}

fn settings_form(ui: &mut egui::Ui, w: &mut AlnPrimersWin, columns: usize, asks: &mut Asks) {
    ui.label(format!(
        "Primers are picked only from columns that enough of the {columns}-column alignment's sequences share, so that one order fits them."
    ));
    let s = &mut w.settings.primers;
    egui::Grid::new("aln_p3_settings")
        .num_columns(7)
        .show(ui, |ui| {
            ui.label("Size");
            count_box(ui, &mut s.min_size, 10..=36);
            count_box(ui, &mut s.opt_size, 10..=36);
            count_box(ui, &mut s.max_size, 10..=36);
            ui.label("Product");
            count_box(ui, &mut s.product_min, 30..=10_000);
            count_box(ui, &mut s.product_max, 30..=10_000);
            ui.end_row();
            ui.label("Tm");
            ui.add(egui::DragValue::new(&mut s.min_tm).range(30.0..=90.0));
            ui.add(egui::DragValue::new(&mut s.opt_tm).range(30.0..=90.0));
            ui.add(egui::DragValue::new(&mut s.max_tm).range(30.0..=90.0));
            ui.label("Pairs");
            count_box(ui, &mut s.num_return, 1..=20);
            ui.end_row();
        });
    ui.label(
        RichText::new(
            "min / optimum / max for size and Tm; the product is measured in alignment columns",
        )
        .small(),
    );
    ui.horizontal(|ui| {
        ui.label("Columns must agree in at least");
        ui.add(
            egui::DragValue::new(&mut w.percent)
                .range(1.0..=100.0)
                .suffix(" %"),
        );
        ui.label("of the sequences; a sequence is covered with up to");
        count_box(ui, &mut w.settings.max_mismatches, 0..=3);
        ui.label("mismatches, none at a primer's 3' end");
    });
    ui.horizontal(|ui| {
        ui.checkbox(&mut w.use_region, "Only columns");
        let top = columns.max(1);
        w.region.1 = w.region.1.clamp(1, top);
        w.region.0 = w.region.0.clamp(1, w.region.1);
        count_box(ui, &mut w.region.0, 1..=top);
        ui.label("to");
        count_box(ui, &mut w.region.1, 1..=top);
    });
    ui.horizontal(|ui| {
        ui.checkbox(
            &mut w.settings.with_probe,
            "Also pick a TaqMan probe between the primers",
        );
        if w.settings.with_probe
            && ui
                .small_button("Product 70–150")
                .on_hover_text(
                    "Real-time PCR products are short; a probe needs room between the primers",
                )
                .clicked()
        {
            (
                w.settings.primers.product_min,
                w.settings.primers.product_max,
            ) = TAQMAN_PRODUCT;
        }
        asks.pick = ui
            .button(RichText::new("Pick primers and probes").strong())
            .clicked();
    });
}

fn results(
    ui: &mut egui::Ui,
    look: &crate::ui::design::Look,
    w: &AlnPrimersWin,
    report: &Report,
    asks: &mut Asks,
) {
    if report.designs.is_empty() {
        return;
    }
    ui.horizontal(|ui| {
        for role in [Role::Left, Role::Right, Role::Probe] {
            ui.label(
                RichText::new(format!("■ {}", role.name()))
                    .color(role.colour(look))
                    .small(),
            );
        }
        ui.label(RichText::new("(the colours of the marks on the alignment)").small());
    });
    let probes = report.designs.iter().any(|d| d.probe.is_some());
    egui::ScrollArea::horizontal().show(ui, |ui| {
        egui::Grid::new("aln_pairs").striped(true).show(ui, |ui| {
            let mut heads = vec!["#", "Left primer", "Tm", "Right primer", "Tm"];
            if probes {
                heads.extend(["Probe", "Tm"]);
            }
            heads.extend(["Product", "Covers", "Penalty", ""]);
            for h in heads {
                ui.label(RichText::new(h).strong());
            }
            ui.end_row();
            for (k, d) in report.designs.iter().enumerate() {
                design_row(ui, k, d, report.rows, probes, asks);
            }
        });
    });
    if let Some(d) = report.designs.get(w.chosen) {
        details(ui, w, report, d);
    }
    ui.horizontal(|ui| {
        asks.copy = ui.button("Copy table").clicked();
        asks.save = ui.button("Save table…").clicked();
    });
}

fn design_row(ui: &mut egui::Ui, k: usize, d: &Design, rows: usize, probes: bool, asks: &mut Asks) {
    let seq =
        |s: &Site| RichText::new(String::from_utf8_lossy(&s.oligo.seq).into_owned()).monospace();
    ui.label((k + 1).to_string());
    ui.label(seq(&d.left));
    ui.label(format!("{:.1}", d.left.oligo.tm));
    ui.label(seq(&d.right));
    ui.label(format!("{:.1}", d.right.oligo.tm));
    if probes {
        match &d.probe {
            Some(p) => {
                ui.label(seq(p));
                ui.label(format!("{:.1}", p.oligo.tm));
            }
            None => {
                ui.label("");
                ui.label("");
            }
        }
    }
    ui.label(product_text(d));
    ui.label(format!("{} / {rows}", d.tolerated))
        .on_hover_text(format!("{} match every oligo perfectly", d.exact));
    ui.label(format!("{:.3}", d.penalty));
    ui.horizontal(|ui| {
        if ui
            .small_button("Show")
            .on_hover_text("Mark the oligos on the alignment")
            .clicked()
        {
            asks.show = Some(k);
        }
        if ui
            .small_button("Consensus")
            .on_hover_text("Open the consensus with the oligos as features")
            .clicked()
        {
            asks.consensus = Some(k);
        }
    });
    ui.end_row();
}

fn details(ui: &mut egui::Ui, w: &AlnPrimersWin, report: &Report, d: &Design) {
    let place = |s: &Site| {
        let reference = s.in_reference.map_or(String::new(), |(a, b)| {
            format!(
                ", {a}–{b} of {}",
                w.names.first().map_or("the first sequence", String::as_str)
            )
        });
        format!("columns {}–{}{reference}", s.oligo.start + 1, s.oligo.end)
    };
    ui.separator();
    ui.label(RichText::new(format!("Design {}", w.chosen + 1)).strong());
    ui.label(format!("Left primer: {}.", place(&d.left)));
    ui.label(format!("Right primer: {}.", place(&d.right)));
    if let Some(p) = &d.probe {
        ui.label(format!("Probe: {}.", place(p)));
    }
    if d.uncovered.is_empty() {
        ui.label(format!("Covers all {} sequences.", report.rows));
    } else {
        let names: Vec<&str> = d
            .uncovered
            .iter()
            .filter_map(|&r| w.names.get(r))
            .map(String::as_str)
            .collect();
        ui.label(format!(
            "Does not cover {} of {}: {}.",
            names.len(),
            report.rows,
            names.join(", ")
        ));
    }
}

#[cfg(test)]
pub(super) fn test_asks_show(k: usize) -> Asks {
    Asks {
        show: Some(k),
        ..Asks::default()
    }
}

#[cfg(test)]
pub(super) fn test_asks_consensus(k: usize) -> Asks {
    Asks {
        consensus: Some(k),
        ..Asks::default()
    }
}
