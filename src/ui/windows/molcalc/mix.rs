//! The Master Mix Multiplier tab: the volume of each component of a master mix for a number of reactions
//! (`calc_mm_reaction`).

use super::{combo, number_field, result};
use crate::seq::molcalc::format::{number, tcl_double};
use crate::seq::molcalc::mix::{copy_text, master_mix, Extra, MixEdit, MixInput, COMPONENTS};
use crate::seq::molcalc::presets::{DEFAULTS, MIX_PRESETS};
use eframe::egui::{self, RichText, Ui};

pub struct MixTab {
    preset: usize,
    names: Vec<String>,
    volumes: Vec<String>,
    single_volume: String,
    reactions: String,
    extra: String,
    extra_unit: Extra,
}

impl Default for MixTab {
    fn default() -> MixTab {
        let mut tab = MixTab {
            preset: 0,
            names: vec![String::new(); COMPONENTS],
            volumes: vec![String::new(); COMPONENTS],
            single_volume: String::new(),
            reactions: DEFAULTS.reactions.into(),
            extra: DEFAULTS.extra.into(),
            extra_unit: Extra::Percent,
        };
        let first = MIX_PRESETS
            .iter()
            .position(|p| p.name == DEFAULTS.mix_reaction)
            .unwrap_or(0);
        tab.choose(first);
        tab
    }
}

impl MixTab {
    /// Take the components of a master mix; the single reaction volume is their sum.
    fn choose(&mut self, index: usize) {
        self.preset = index;
        self.names = vec![String::new(); COMPONENTS];
        self.volumes = vec![String::new(); COMPONENTS];
        for (i, (name, volume)) in MIX_PRESETS[index].components.iter().enumerate() {
            self.names[i] = (*name).into();
            self.volumes[i] = (*volume).into();
        }
        let sum: f64 = self.volumes.iter().filter_map(|v| number(v)).sum();
        self.single_volume = tcl_double(sum);
    }

    fn input(&self) -> MixInput<'_> {
        MixInput {
            components: &self.volumes,
            single_volume: &self.single_volume,
            reactions: &self.reactions,
            extra: &self.extra,
            extra_unit: self.extra_unit,
        }
    }

    pub fn show(&mut self, ui: &mut Ui) {
        self.reaction_row(ui);
        // what is shown beside the fields; what was typed decides what is worked out again afterwards
        let shown = master_mix(&self.input(), MixEdit::Other);
        let mut edit = self.settings(ui);
        ui.separator();
        edit = self.component_rows(ui, edit, &shown.volumes);
        if edit != MixEdit::Other {
            let out = master_mix(&self.input(), edit);
            self.volumes = out.components;
            self.single_volume = out.single_volume;
        }
        ui.horizontal(|ui| {
            ui.label("Total Mix Volume");
            result(ui, &shown.total);
        });
        if ui.button("Copy").clicked() {
            let name = MIX_PRESETS[self.preset].name;
            let text = copy_text(
                name,
                &self.names,
                &self.single_volume,
                &self.reactions,
                &shown,
            );
            ui.ctx().copy_text(text);
        }
    }

    fn reaction_row(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.label("Reaction");
            egui::ComboBox::from_id_salt("mix_reaction")
                .selected_text(MIX_PRESETS[self.preset].name)
                .show_ui(ui, |ui| {
                    for (i, p) in MIX_PRESETS.iter().enumerate() {
                        if ui.selectable_label(self.preset == i, p.name).clicked() {
                            self.choose(i);
                        }
                    }
                });
        });
    }

    /// The single reaction volume, the number of reactions and the extra; typing the volume rescales the components.
    fn settings(&mut self, ui: &mut Ui) -> MixEdit {
        let mut edit = MixEdit::Other;
        egui::Grid::new("mix_settings")
            .num_columns(2)
            .show(ui, |ui| {
                ui.label("Single Reaction Volume");
                if number_field(ui, &mut self.single_volume, 64.0) {
                    edit = MixEdit::SingleVolume;
                }
                ui.end_row();
                ui.label("# Reactions");
                number_field(ui, &mut self.reactions, 56.0);
                ui.end_row();
                ui.label("plus");
                ui.horizontal(|ui| {
                    number_field(ui, &mut self.extra, 46.0);
                    let units = [Extra::Percent, Extra::Ul];
                    combo(ui, "mix_extra", &mut self.extra_unit, &units, Extra::label);
                });
                ui.end_row();
            });
        edit
    }

    /// The names and volumes in one reaction, and each volume in the mix; typing a volume means the single reaction volume
    /// follows.
    fn component_rows(&mut self, ui: &mut Ui, edit: MixEdit, mixed: &[String]) -> MixEdit {
        let mut edit = edit;
        egui::Grid::new("mix_components")
            .num_columns(3)
            .show(ui, |ui| {
                for heading in ["Component Name", "Volume/ Reaction", "Volume"] {
                    ui.label(RichText::new(heading).weak());
                }
                ui.end_row();
                let rows = self.names.iter_mut().zip(&mut self.volumes).zip(mixed);
                for ((name, volume), in_mix) in rows {
                    ui.add(egui::TextEdit::singleline(name).desired_width(150.0));
                    if number_field(ui, volume, 64.0) {
                        edit = MixEdit::Component;
                    }
                    result(ui, in_mix);
                    ui.end_row();
                }
            });
        edit
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_tab_opens_on_phusion_and_its_components_add_up_to_the_single_volume() {
        let tab = MixTab::default();
        assert_eq!(tab.single_volume, "50.0");
        assert_eq!(tab.names[..3], ["primer", "primer", "enzyme"]);
        assert_eq!((tab.reactions.as_str(), tab.extra.as_str()), ("1", "0"));
        let out = master_mix(&tab.input(), MixEdit::Other);
        assert_eq!(out.volumes[0], "0.25");
        assert_eq!(out.total, "50.0");
        // for 8 reactions with 10 % extra every volume is 8.8 times as much
        let tab = MixTab {
            reactions: "8".into(),
            extra: "10".into(),
            ..MixTab::default()
        };
        let out = master_mix(&tab.input(), MixEdit::Other);
        assert_eq!(
            (out.volumes[4].as_str(), out.total.as_str()),
            ("334.4", "440.0")
        );
    }
}
