//! The Molar Ratio Reaction tab: how much of each fragment to pipette (`calc_molar_reaction`).

use super::{combo, number_field, result};
use crate::seq::molcalc::presets::{enzyme_label, DEFAULTS, ENZYME_STOCKS, MOLAR_PRESETS};
use crate::seq::molcalc::reaction::{
    copy_text, molar_reaction, DnaUnits, Fragment, MolarInput, MolarOutput, FRAGMENTS,
};
use eframe::egui::{self, RichText, Ui};

pub struct MolarTab {
    preset: usize,
    fragments: Vec<Fragment>,
    volume: String,
    final_amount: String,
    units: DnaUnits,
    min_pipette: String,
    enzyme_x: String,
    /// The enzyme stock is typed (the menu's Custom…).
    custom_enzyme: bool,
    add_water: bool,
    notes: String,
}

impl Default for MolarTab {
    fn default() -> MolarTab {
        let mut tab = MolarTab {
            preset: 0,
            fragments: vec![Fragment::default(); FRAGMENTS],
            volume: String::new(),
            final_amount: String::new(),
            units: DnaUnits::Fmole,
            min_pipette: String::new(),
            enzyme_x: String::new(),
            custom_enzyme: false,
            add_water: true,
            notes: String::new(),
        };
        // the window opens on the reaction the original opens on
        let first = MOLAR_PRESETS
            .iter()
            .position(|p| p.name == DEFAULTS.molar_reaction)
            .unwrap_or(0);
        tab.choose(first);
        tab
    }
}

impl MolarTab {
    /// Take the conditions of a reaction. The fragments stay as they are: the original keeps them between reactions.
    fn choose(&mut self, index: usize) {
        let p = &MOLAR_PRESETS[index];
        self.preset = index;
        self.volume = p.volume.into();
        self.final_amount = p.final_amount.into();
        self.units = p.units;
        self.min_pipette = p.min_pipette.into();
        self.enzyme_x = p.enzyme_x.into();
        self.custom_enzyme = enzyme_label(p.enzyme_x).is_none();
        self.notes = p.notes.into();
    }

    /// Two fragments to see a result with: a 3.5 kb backbone and a 1.2 kb insert, three insert molecules to one backbone.
    pub fn example(&mut self) {
        let row = |name: &str, size: &str, conc: &str, ratio: &str| Fragment {
            name: name.into(),
            size_kb: size.into(),
            conc: conc.into(),
            ratio: ratio.into(),
        };
        self.fragments[0] = row("backbone", "3.5", "80", "1");
        self.fragments[1] = row("insert", "1.2", "45", "3");
    }

    fn input(&self) -> MolarInput<'_> {
        MolarInput {
            fragments: &self.fragments,
            reaction_volume: &self.volume,
            final_amount: &self.final_amount,
            units: self.units,
            min_pipette: &self.min_pipette,
            enzyme_x: &self.enzyme_x,
            add_water: self.add_water,
        }
    }

    pub fn show(&mut self, ui: &mut Ui) {
        self.reaction_row(ui);
        let out = molar_reaction(&self.input());
        // a ratio left empty is filled in with 1 once the row has a size and a concentration
        for (f, row) in self.fragments.iter_mut().zip(&out.rows) {
            if f.ratio.trim().is_empty() && !row.ratio.is_empty() {
                f.ratio = row.ratio.clone();
            }
        }
        self.fragment_rows(ui, &out);
        ui.separator();
        self.volumes(ui, &out);
        ui.separator();
        self.conditions(ui, out.too_dilute);
        let name = MOLAR_PRESETS[self.preset].name;
        if ui.button("Copy").clicked() {
            ui.ctx()
                .copy_text(copy_text(name, &self.fragments, &out, &self.notes));
        }
    }

    fn reaction_row(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.label("Reaction");
            egui::ComboBox::from_id_salt("molar_reaction")
                .selected_text(MOLAR_PRESETS[self.preset].name)
                .show_ui(ui, |ui| {
                    for (i, p) in MOLAR_PRESETS.iter().enumerate() {
                        if ui.selectable_label(self.preset == i, p.name).clicked() {
                            self.choose(i);
                        }
                    }
                });
            if ui
                .button("Reset Values")
                .on_hover_text("Clear the fragments")
                .clicked()
            {
                self.fragments = vec![Fragment::default(); FRAGMENTS];
            }
        });
    }

    fn fragment_rows(&mut self, ui: &mut Ui, out: &MolarOutput) {
        egui::Grid::new("molar_fragments")
            .num_columns(6)
            .spacing([6.0, 3.0])
            .show(ui, |ui| {
                let amount = if out.amount_label.is_empty() {
                    "fmoles"
                } else {
                    out.amount_label
                };
                for heading in [
                    "Fragment Name",
                    "Size (kb)",
                    "ng/ul",
                    "ratio",
                    amount,
                    "volume (ul)",
                ] {
                    ui.label(RichText::new(heading).weak());
                }
                ui.end_row();
                for (f, row) in self.fragments.iter_mut().zip(&out.rows) {
                    ui.add(egui::TextEdit::singleline(&mut f.name).desired_width(130.0));
                    number_field(ui, &mut f.size_kb, 56.0);
                    number_field(ui, &mut f.conc, 56.0);
                    number_field(ui, &mut f.ratio, 36.0);
                    ui.label(&row.amount);
                    result(ui, &row.volume);
                    ui.end_row();
                }
            });
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.add_water, "Water");
            result(ui, &out.water);
            ui.label("ul");
            ui.add_space(12.0);
            ui.label("DNA Mix Intermediate Volume");
            result(ui, &out.intermediate);
        });
    }

    fn volumes(&mut self, ui: &mut Ui, out: &MolarOutput) {
        egui::Grid::new("molar_volumes")
            .num_columns(2)
            .show(ui, |ui| {
                ui.label("Reaction Volume (ul)");
                number_field(ui, &mut self.volume, 56.0);
                ui.end_row();
                ui.label("DNA Mix Volume (ul)");
                result(ui, &out.dna_volume);
                ui.end_row();
                ui.label("Enzyme+Buffer Mix Volume");
                result(ui, &out.enzyme_volume);
                ui.end_row();
            });
    }

    /// The conditions of the reaction. The label of the final amount is red while the DNA is too dilute for the volume
    /// (the original colours that label, `conditions_frame.l1`, and not the line that says so).
    fn conditions(&mut self, ui: &mut Ui, too_dilute: bool) {
        ui.label(RichText::new("Reaction Conditions").weak());
        egui::Grid::new("molar_conditions")
            .num_columns(2)
            .show(ui, |ui| {
                let label = RichText::new("Final DNA amount/concentration");
                ui.label(if too_dilute {
                    label.color(ui.visuals().error_fg_color)
                } else {
                    label
                });
                ui.horizontal(|ui| {
                    number_field(ui, &mut self.final_amount, 56.0);
                    combo(
                        ui,
                        "molar_units",
                        &mut self.units,
                        &DnaUnits::ALL,
                        DnaUnits::label,
                    );
                    ui.label(RichText::new("1 nM= 1 fmole/ul = 0.001 pmole/ul").weak());
                });
                ui.end_row();
                ui.label("Minimum pipette volume (ul)");
                number_field(ui, &mut self.min_pipette, 56.0);
                ui.end_row();
                ui.label("Enzyme Stock concentration (X)");
                self.enzyme_stock(ui);
                ui.end_row();
                ui.label("Reaction Notes");
                ui.add(egui::TextEdit::singleline(&mut self.notes).desired_width(380.0));
                ui.end_row();
            });
    }

    /// The enzyme mix is bought as 10x, 5x, 2x…: the menu of the original, or any number typed.
    fn enzyme_stock(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            let shown = enzyme_label(&self.enzyme_x).filter(|_| !self.custom_enzyme);
            egui::ComboBox::from_id_salt("molar_enzyme")
                .selected_text(shown.unwrap_or("Custom..."))
                .show_ui(ui, |ui| {
                    for (label, x) in ENZYME_STOCKS {
                        if ui.selectable_label(shown == Some(label), label).clicked() {
                            self.enzyme_x = x.into();
                            self.custom_enzyme = false;
                        }
                    }
                    if ui.selectable_label(shown.is_none(), "Custom...").clicked() {
                        // the original starts a typed stock at 10
                        self.enzyme_x = "10".into();
                        self.custom_enzyme = true;
                    }
                });
            if shown.is_none() {
                number_field(ui, &mut self.enzyme_x, 46.0);
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_tab_opens_on_the_first_reaction_of_the_original_and_keeps_the_fragments_between_reactions(
    ) {
        let mut tab = MolarTab::default();
        assert_eq!(
            (
                tab.volume.as_str(),
                tab.final_amount.as_str(),
                tab.units,
                tab.min_pipette.as_str(),
                tab.enzyme_x.as_str()
            ),
            ("4", "40", DnaUnits::Fmole, "0.3", "2")
        );
        assert!(!tab.custom_enzyme && tab.add_water && tab.notes.starts_with("6 to 40 fmole"));
        tab.example();
        let golden_gate = MOLAR_PRESETS.iter().position(|p| p.name == "Golden Gate");
        tab.choose(golden_gate.expect("a reaction"));
        assert_eq!(
            (
                tab.volume.as_str(),
                tab.final_amount.as_str(),
                tab.units,
                tab.enzyme_x.as_str()
            ),
            ("10", "6", DnaUnits::Nm, "5")
        );
        assert_eq!(tab.fragments[0].name, "backbone", "the fragments are kept");
    }

    #[test]
    fn the_example_works_out_volumes_and_a_ratio_left_empty_is_filled_in() {
        let mut tab = MolarTab::default();
        tab.example();
        tab.fragments[0].ratio.clear();
        let out = molar_reaction(&tab.input());
        assert_eq!(out.rows[0].ratio, "1");
        assert!(!out.rows[0].volume.is_empty() && !out.rows[1].volume.is_empty());
        // three insert molecules to each backbone molecule
        let (backbone, insert) = (&out.rows[0], &out.rows[1]);
        let molecules = |amount: &str| amount.parse::<f64>().expect("a number");
        assert!((molecules(&insert.amount) / molecules(&backbone.amount) - 3.0).abs() < 0.01);
    }
}
