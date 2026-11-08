//! Enzyme selector: pick the enzymes used by the map, highlighting and digestion.

use crate::seq::enzymes::Enzyme;
use crate::ui::app::App;
use eframe::egui::{self, RichText, ScrollArea};

#[derive(Clone, Copy, PartialEq, Eq, Default, Debug)]
pub enum Group {
    #[default]
    All,
    Blunt,
    FivePrime,
    ThreePrime,
    Four,
    Five,
    Six,
    LongerThanSix,
    Degenerate,
    Unique,
    Absent,
    Selected,
}

impl Group {
    pub const ALL: [Group; 12] = [
        Group::All,
        Group::Blunt,
        Group::FivePrime,
        Group::ThreePrime,
        Group::Four,
        Group::Five,
        Group::Six,
        Group::LongerThanSix,
        Group::Degenerate,
        Group::Unique,
        Group::Absent,
        Group::Selected,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Group::All => "All enzymes",
            Group::Blunt => "Blunt",
            Group::FivePrime => "5' overhang",
            Group::ThreePrime => "3' overhang",
            Group::Four => "Four cutters",
            Group::Five => "Five cutters",
            Group::Six => "Six cutters",
            Group::LongerThanSix => "Longer than six",
            Group::Degenerate => "Degenerate sites",
            Group::Unique => "Unique cutters",
            Group::Absent => "Not in this sequence",
            Group::Selected => "Selected",
        }
    }

    pub fn includes(self, e: &Enzyme, sites: usize, selected: bool) -> bool {
        match self {
            Group::All => true,
            Group::Blunt => e.overhang() == 0,
            Group::FivePrime => e.overhang() > 0,
            Group::ThreePrime => e.overhang() < 0,
            Group::Four => e.cutter_size() == 4,
            Group::Five => e.cutter_size() == 5,
            Group::Six => e.cutter_size() == 6,
            Group::LongerThanSix => e.cutter_size() > 6,
            Group::Degenerate => e.is_degenerate(),
            Group::Unique => sites == 1,
            Group::Absent => sites == 0,
            Group::Selected => selected,
        }
    }
}

#[derive(Default)]
pub struct EnzymeWin {
    pub open: bool,
    pub filter: String,
    pub group: Group,
}

impl App {
    /// The Selection Only switch of the Enzymes menu (`info(analyze_selection)`, L1170; Cmd+Shift+E, L11491; Shift-click on the
    /// toolbar's enzyme button, L4370): the digest, the lists of enzymes and the ORF list look at the selected bases, or at all
    /// of the sequence when nothing is selected.
    pub fn toggle_selection_only(&mut self) {
        self.prefs.selection_only = !self.prefs.selection_only;
        let text = if self.prefs.selection_only {
            "Selection Only: on. The digest, the lists of enzymes and the ORF list look at the selected bases."
        } else {
            "Selection Only: off. The digest, the lists of enzymes and the ORF list look at all of the sequence."
        };
        self.say(text);
    }

    pub fn enzyme_window(&mut self, ctx: &egui::Context) {
        if !self.enzyme_win.open {
            return;
        }
        let mut open = true;
        crate::ui::design::window("Enzyme Selector")
            .open(&mut open)
            .default_size([430.0, 480.0])
            .show(ctx, |ui| {
                let counts: Vec<usize> = match self.current_doc() {
                    Some(i) => self.docs[i]
                        .sites(&self.enzymes)
                        .iter()
                        .map(Vec::len)
                        .collect(),
                    None => vec![0; self.enzymes.len()],
                };
                ui.horizontal(|ui| {
                    ui.label("Show");
                    egui::ComboBox::from_id_salt("enzyme_group")
                        .selected_text(self.enzyme_win.group.label())
                        .show_ui(ui, |ui| {
                            for g in Group::ALL {
                                ui.selectable_value(&mut self.enzyme_win.group, g, g.label());
                            }
                        });
                    ui.label("Name");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.enzyme_win.filter).desired_width(90.0),
                    );
                });
                let shown: Vec<usize> = (0..self.enzymes.len())
                    .filter(|&i| {
                        let e = &self.enzymes[i];
                        let needle = self.enzyme_win.filter.to_ascii_lowercase();
                        (needle.is_empty() || e.name.to_ascii_lowercase().contains(&needle))
                            && self.enzyme_win.group.includes(
                                e,
                                counts[i],
                                self.selected.contains(&e.name),
                            )
                    })
                    .collect();
                ui.horizontal(|ui| {
                    if ui.button("Select shown").clicked() {
                        self.selected
                            .extend(shown.iter().map(|&i| self.enzymes[i].name.clone()));
                    }
                    if ui.button("Unselect shown").clicked() {
                        for &i in &shown {
                            self.selected.remove(&self.enzymes[i].name);
                        }
                    }
                    if ui.button("Clear all").clicked() {
                        self.selected.clear();
                    }
                    ui.checkbox(&mut self.plus_unique, "Map also shows unique cutters");
                });
                ui.label(
                    RichText::new(format!(
                        "{} selected, {} shown",
                        self.selected.len(),
                        shown.len()
                    ))
                    .small(),
                );
                ui.separator();
                ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        egui::Grid::new("enzyme_grid")
                            .num_columns(4)
                            .striped(true)
                            .spacing([14.0, 2.0])
                            .show(ui, |ui| {
                                for &i in &shown {
                                    let e = &self.enzymes[i];
                                    let mut on = self.selected.contains(&e.name);
                                    if ui.checkbox(&mut on, &e.name).changed() {
                                        if on {
                                            self.selected.insert(e.name.clone());
                                        } else {
                                            self.selected.remove(&e.name);
                                        }
                                    }
                                    ui.label(RichText::new(e.pattern_text()).monospace());
                                    ui.label(match e.overhang() {
                                        0 => "blunt".to_string(),
                                        n if n > 0 => format!("5' {n}"),
                                        n => format!("3' {}", -n),
                                    });
                                    ui.label(
                                        RichText::new(format!(
                                            "{} site{}",
                                            counts[i],
                                            if counts[i] == 1 { "" } else { "s" }
                                        ))
                                        .small(),
                                    );
                                    ui.end_row();
                                }
                            });
                    });
            });
        self.enzyme_win.open = open;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::enzymes::default_enzymes;

    #[test]
    fn groups_filter_by_overhang_and_cutter_size() {
        let all = default_enzymes();
        let by = |n: &str| all.iter().find(|e| e.name == n).unwrap();
        assert!(Group::Blunt.includes(by("EcoRV"), 0, false));
        assert!(!Group::Blunt.includes(by("EcoRI"), 0, false));
        assert!(Group::FivePrime.includes(by("EcoRI"), 0, false));
        assert!(Group::ThreePrime.includes(by("PstI"), 0, false));
        assert!(Group::Four.includes(by("AluI"), 0, false));
        assert!(Group::Six.includes(by("EcoRI"), 0, false));
        assert!(Group::LongerThanSix.includes(by("NotI"), 0, false));
        assert!(Group::Degenerate.includes(by("ApoI"), 0, false));
        assert!(Group::Unique.includes(by("EcoRI"), 1, false));
        assert!(!Group::Unique.includes(by("EcoRI"), 2, false));
        assert!(Group::Absent.includes(by("EcoRI"), 0, false));
        assert!(Group::Selected.includes(by("EcoRI"), 3, true));
    }

    #[test]
    fn the_selection_only_switch_starts_off_and_turns_over_with_a_message() {
        let mut app = App::new(crate::ui::app::Startup::default());
        assert!(!app.prefs.selection_only);

        app.toggle_selection_only();
        assert!(app.prefs.selection_only);
        assert!(
            app.status.starts_with("Selection Only: on"),
            "{}",
            app.status
        );

        app.toggle_selection_only();
        assert!(!app.prefs.selection_only);
        assert!(
            app.status.starts_with("Selection Only: off"),
            "{}",
            app.status
        );
    }

    #[test]
    fn the_palette_finds_the_switch_by_its_name_and_by_orf_and_it_works_from_any_page() {
        use crate::ui::app::Action;
        use crate::ui::commands::{acts_on_the_sequence, all, search};
        let commands = all();
        let finds = |typed: &str| {
            let hits = search(commands, typed, 8);
            hits.iter()
                .any(|&i| matches!(commands[i].action, Action::ToggleSelectionOnly))
        };
        assert!(finds("selection only"));
        assert!(finds("orf selection"));
        // a preference: it turns over wherever the person is, without a trip to the sequence
        assert!(!acts_on_the_sequence(&Action::ToggleSelectionOnly));
    }
}
