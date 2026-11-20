//! Tools → Primer3 Pair Design: primer pairs picked with primer3's rules around the selection, each with a TaqMan probe between
//! its primers when that is asked for.

use crate::seq::features::Feature;
use crate::seq::primer3::{design, Oligo, Pair, Settings};
use crate::seq::probe::{design_assays, Settings as ProbeSettings};
use crate::ui::app::App;
use crate::ui::fonts::FontWeight;
use eframe::egui::{self, RichText};

#[derive(Default)]
pub struct Primer3Win {
    pub open: bool,
    pub settings: Settings,
    /// Also pick a probe between the primers of each pair.
    pub with_probe: bool,
    pub probe: ProbeSettings,
    /// The pairs found, each with its probe when one was asked for.
    pub found: Vec<(Pair, Option<Oligo>)>,
    /// A checksum of the sequence that `found` was picked for: features are put on the sequence that is open only if it is that one.
    found_for: Option<u64>,
    pub message: String,
}

/// The product sizes that real-time PCR assays use (Applied Biosystems recommends 50 to 150 bases; a probe needs room between the
/// primers).
const TAQMAN_PRODUCT: (usize, usize) = (70, 150);

fn base_text(seq: &[u8]) -> String {
    String::from_utf8_lossy(seq).into_owned()
}

fn sequence_hash(seq: &[u8]) -> u64 {
    crate::kept::hash(&String::from_utf8_lossy(seq))
}

impl App {
    /// Pick pairs around the selection, or anywhere in the sequence when nothing is selected.
    // ponytail: linear template, on the UI thread; a product across the origin of a circle is not
    // found, and a template of many kilobases makes the window wait a few seconds.
    pub fn run_primer3(&mut self) {
        let Some(i) = self.current_doc() else { return };
        let dv = &self.docs[i];
        let target = dv.selection().map(|(a, b)| (a, b - a));
        let (settings, probe) = (&self.p3_win.settings, &self.p3_win.probe);
        let with_probe = self.p3_win.with_probe;
        let found: Vec<(Pair, Option<Oligo>)> = if with_probe {
            design_assays(&dv.doc.seq, target, settings, probe)
                .into_iter()
                .map(|a| (a.pair, Some(a.probe)))
                .collect()
        } else {
            design(&dv.doc.seq, target, settings)
                .into_iter()
                .map(|pair| (pair, None))
                .collect()
        };
        self.p3_win.message = if !found.is_empty() {
            let around = if target.is_some() {
                "around the selection"
            } else {
                "in the whole sequence"
            };
            let (what, order) = if with_probe {
                ("assays (pair and probe)", "lowest total penalty first")
            } else {
                ("best pairs", "lowest penalty first")
            };
            format!("{} {what} {around}, {order}.", found.len())
        } else if with_probe && !design(&dv.doc.seq, target, settings).is_empty() {
            // the primers are not what is missing
            "Pairs fit, but none of the best has a probe between its primers that fits the probe limits; widen the probe's Tm range, size or G+C, or the product range.".into()
        } else {
            "No pair fits these limits; widen the Tm, size or product range.".into()
        };
        self.p3_win.found_for = Some(sequence_hash(&dv.doc.seq));
        self.p3_win.found = found;
    }

    /// Mark both primers of pair `k` on the sequence as `primer_bind` features, and its probe, when it has one, as a `misc_feature`
    /// on the strand it is ordered on, with its sequence in a note.
    pub(super) fn add_pair_features(&mut self, k: usize) {
        let (Some(i), Some((pair, probe))) =
            (self.current_doc(), self.p3_win.found.get(k).cloned())
        else {
            return;
        };
        // the pairs were picked for one sequence: another tab, or the same one edited, is not where they belong
        if self.p3_win.found_for != Some(sequence_hash(&self.docs[i].doc.seq)) {
            self.say("This is not the sequence the pairs were picked for; pick them again.");
            return;
        }
        let n = k + 1;
        let mut made: Vec<Feature> = [
            (format!("Left primer {n}"), &pair.left),
            (format!("Right primer {n}"), &pair.right),
        ]
        .map(|(name, o)| Feature::new(name, "primer_bind", o.start, o.end, o.reverse))
        .into();
        let said = if let Some(p) = &probe {
            let mut feature = Feature::new(
                format!("Probe {n}"),
                "misc_feature",
                p.start,
                p.end,
                p.reverse,
            );
            feature.qualifiers.push((
                "note".to_string(),
                format!("TaqMan probe, 5'-{}-3'", base_text(&p.seq)),
            ));
            made.push(feature);
            format!("Added the two primers and the probe of assay {n} as features.")
        } else {
            format!("Added the two primers of pair {n} as features.")
        };
        self.docs[i].edit(|d| d.features.extend(made));
        self.say(said);
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
                        pick = ui
                            .button(RichText::new("Pick primers").semibold())
                            .clicked();
                        ui.end_row();
                    });
                ui.label(RichText::new("min / optimum / max for size and Tm").small());
                ui.horizontal(|ui| {
                    ui.checkbox(
                        &mut self.p3_win.with_probe,
                        "Also pick a TaqMan probe between the primers",
                    );
                    if self.p3_win.with_probe
                        && ui
                            .small_button("Product 70–150")
                            .on_hover_text("Real-time PCR products are short; a probe needs room between the primers")
                            .clicked()
                    {
                        (self.p3_win.settings.product_min, self.p3_win.settings.product_max) =
                            TAQMAN_PRODUCT;
                    }
                });
                if self.p3_win.with_probe {
                    probe_settings(ui, &mut self.p3_win.probe);
                }
                if !self.p3_win.message.is_empty() {
                    ui.label(&self.p3_win.message);
                }
                let probes = self.p3_win.found.iter().any(|(_, probe)| probe.is_some());
                egui::Grid::new("p3_pairs").striped(true).show(ui, |ui| {
                    let probe_heads = ["Probe", "Tm", "G+C", "Strand"];
                    // with a probe, the penalty is the pair's and the probe's together
                    let penalty = if probes { "Total penalty" } else { "Penalty" };
                    let heads = ["#", "Left primer", "Tm", "Right primer", "Tm", "Size", penalty];
                    let shown = heads.iter().chain(probes.then_some(&probe_heads).into_iter().flatten());
                    for h in shown.chain(&[""]) {
                        ui.label(RichText::new(*h).semibold());
                    }
                    ui.end_row();
                    for (k, (p, probe)) in self.p3_win.found.iter().enumerate() {
                        ui.label((k + 1).to_string());
                        ui.label(RichText::new(base_text(&p.left.seq)).monospace());
                        ui.label(format!("{:.1}", p.left.tm));
                        ui.label(RichText::new(base_text(&p.right.seq)).monospace());
                        ui.label(format!("{:.1}", p.right.tm));
                        ui.label(p.product_size.to_string());
                        let own = probe.as_ref().map_or(0.0, |o| o.penalty);
                        ui.label(format!("{:.3}", p.penalty + own));
                        if let Some(o) = probe {
                            ui.label(RichText::new(base_text(&o.seq)).monospace());
                            ui.label(format!("{:.1}", o.tm));
                            ui.label(format!("{:.0} %", o.gc));
                            ui.label(if o.reverse { "minus" } else { "plus" });
                        } else if probes {
                            // the columns of a row without a probe stay under their heads
                            (0..4).for_each(|_| {
                                ui.label("");
                            });
                        }
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

/// The probe's limits: size, how far above the primers it melts, G+C, and the rule about its 5' end.
fn probe_settings(ui: &mut egui::Ui, s: &mut ProbeSettings) {
    egui::Grid::new("p3_probe").num_columns(5).show(ui, |ui| {
        ui.label("Probe size");
        ui.add(egui::DragValue::new(&mut s.min_size).range(10..=40));
        ui.add(egui::DragValue::new(&mut s.opt_size).range(10..=40));
        ui.add(egui::DragValue::new(&mut s.max_size).range(10..=40));
        ui.end_row();
        ui.label("Tm above primers");
        ui.add(egui::DragValue::new(&mut s.min_above).range(-10.0..=40.0));
        ui.add(egui::DragValue::new(&mut s.opt_above).range(-10.0..=40.0));
        ui.add(egui::DragValue::new(&mut s.max_above).range(-10.0..=40.0));
        ui.end_row();
        ui.label("G+C %");
        ui.add(egui::DragValue::new(&mut s.min_gc).range(0.0..=100.0));
        ui.add(egui::DragValue::new(&mut s.max_gc).range(0.0..=100.0));
        ui.checkbox(&mut s.no_five_prime_g, "no G at the 5' end");
        ui.end_row();
    });
    ui.label(
        RichText::new(
            "min / optimum / max. The probe melts above the higher of the two primers; a G next to the reporter dye quenches it. \
             Probes of both strands are tried; minor-groove-binder probes are not modelled.",
        )
        .small(),
    );
}
