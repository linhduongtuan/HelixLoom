//! Tools → Primer3 Pair Design: primer pairs picked with primer3's rules around the selection.

use crate::seq::features::Feature;
use crate::seq::primer3::{design, Pair, Settings};
use crate::ui::app::App;
use eframe::egui::{self, RichText};

#[derive(Default)]
pub struct Primer3Win {
    pub open: bool,
    pub settings: Settings,
    pub pairs: Vec<Pair>,
    pub message: String,
}

fn base_text(seq: &[u8]) -> String {
    String::from_utf8_lossy(seq).into_owned()
}

impl App {
    /// Pick pairs around the selection, or anywhere in the sequence when nothing is selected.
    // ponytail: linear template, on the UI thread; a product across the origin of a circle is not
    // found, and a template of many kilobases makes the window wait a few seconds.
    pub fn run_primer3(&mut self) {
        let Some(i) = self.current_doc() else { return };
        let dv = &self.docs[i];
        let target = dv.selection().map(|(a, b)| (a, b - a));
        let pairs = design(&dv.doc.seq, target, &self.p3_win.settings);
        self.p3_win.message = if pairs.is_empty() {
            "No pair fits these limits; widen the Tm, size or product range.".into()
        } else {
            let around = if target.is_some() {
                "around the selection"
            } else {
                "in the whole sequence"
            };
            format!("{} best pairs {around}, lowest penalty first.", pairs.len())
        };
        self.p3_win.pairs = pairs;
    }

    /// Mark both primers of pair `k` on the sequence as `primer_bind` features.
    fn add_pair_features(&mut self, k: usize) {
        let (Some(i), Some(pair)) = (self.current_doc(), self.p3_win.pairs.get(k).cloned()) else {
            return;
        };
        let n = k + 1;
        let made = [
            (format!("Left primer {n}"), &pair.left),
            (format!("Right primer {n}"), &pair.right),
        ]
        .map(|(name, o)| Feature::new(name, "primer_bind", o.start, o.end, o.reverse));
        self.docs[i].edit(|d| d.features.extend(made));
        self.say(format!("Added the two primers of pair {n} as features."));
    }

    pub fn primer3_window(&mut self, ctx: &egui::Context) {
        if !self.p3_win.open {
            return;
        }
        let target = self
            .current_doc()
            .and_then(|i| self.docs[i].selection())
            .map(|(a, b)| format!("the selected bases {}–{}", a + 1, b))
            .unwrap_or_else(|| "none selected: the whole sequence".into());
        let (mut open, mut pick, mut add) = (true, false, None);
        crate::ui::design::window("Primer3 Pair Design")
            .open(&mut open)
            .default_width(640.0)
            .show(ctx, |ui| {
                ui.label(format!("Product must contain: {target}."));
                let s = &mut self.p3_win.settings;
                egui::Grid::new("p3_settings")
                    .num_columns(7)
                    .show(ui, |ui| {
                        ui.label("Size");
                        ui.add(egui::DragValue::new(&mut s.min_size).range(10..=36));
                        ui.add(egui::DragValue::new(&mut s.opt_size).range(10..=36));
                        ui.add(egui::DragValue::new(&mut s.max_size).range(10..=36));
                        ui.label("Product");
                        ui.add(egui::DragValue::new(&mut s.product_min).range(30..=10_000));
                        ui.add(egui::DragValue::new(&mut s.product_max).range(30..=10_000));
                        ui.end_row();
                        ui.label("Tm");
                        ui.add(egui::DragValue::new(&mut s.min_tm).range(30.0..=90.0));
                        ui.add(egui::DragValue::new(&mut s.opt_tm).range(30.0..=90.0));
                        ui.add(egui::DragValue::new(&mut s.max_tm).range(30.0..=90.0));
                        ui.label("Pairs");
                        ui.add(egui::DragValue::new(&mut s.num_return).range(1..=20));
                        pick = ui.button(RichText::new("Pick primers").strong()).clicked();
                        ui.end_row();
                    });
                ui.label(RichText::new("min / optimum / max for size and Tm").small());
                if !self.p3_win.message.is_empty() {
                    ui.label(&self.p3_win.message);
                }
                egui::Grid::new("p3_pairs").striped(true).show(ui, |ui| {
                    for h in [
                        "#",
                        "Left primer",
                        "Tm",
                        "Right primer",
                        "Tm",
                        "Size",
                        "Penalty",
                        "",
                    ] {
                        ui.label(RichText::new(h).strong());
                    }
                    ui.end_row();
                    for (k, p) in self.p3_win.pairs.iter().enumerate() {
                        ui.label((k + 1).to_string());
                        ui.label(RichText::new(base_text(&p.left.seq)).monospace());
                        ui.label(format!("{:.1}", p.left.tm));
                        ui.label(RichText::new(base_text(&p.right.seq)).monospace());
                        ui.label(format!("{:.1}", p.right.tm));
                        ui.label(p.product_size.to_string());
                        ui.label(format!("{:.3}", p.penalty));
                        if ui.small_button("Add as features").clicked() {
                            add = Some(k);
                        }
                        ui.end_row();
                    }
                });
            });
        if pick {
            self.run_primer3();
        }
        if let Some(k) = add {
            self.add_pair_features(k);
        }
        self.p3_win.open = open;
    }
}
