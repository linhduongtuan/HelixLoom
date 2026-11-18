//! The document view: sequence editor with its map and feature table.

use super::app::{App, Tab};
use super::editor;
use super::infobar;
use super::map::{self, MapData, MapSite};
use super::theme::{self, Palette};
use crate::seq::dna::GENETIC_CODES;
use crate::ui::fonts::FontWeight;
use eframe::egui::{self, RichText, ScrollArea, Ui};

impl App {
    pub(super) fn content(&mut self, ui: &mut Ui) {
        if self.tab == Tab::Assembly {
            self.assembly_tab(ui);
            return;
        }
        let Some(i) = self.current_doc() else { return };
        self.last_doc = i;
        let pal = self.palette();
        egui::Panel::bottom("infobar").show(ui, |ui| {
            infobar::show(ui, &mut self.docs[i], &self.prefs, &pal)
        });
        if self.prefs.show_features && !self.docs[i].doc.features.is_empty() {
            egui::Panel::bottom("feature_table")
                .resizable(true)
                .default_size(130.0)
                .show(ui, |ui| self.feature_table(ui, i));
        }
        if self.prefs.show_map || self.prefs.show_inspector {
            egui::Panel::right("map")
                .resizable(true)
                .default_size(480.0)
                .min_size(260.0)
                .show(ui, |ui| {
                    match (self.prefs.show_map, self.prefs.show_inspector) {
                        (true, true) => {
                            // tall enough, at the start, for the table, the protein and a few overlaps; it can be dragged
                            let tall = (ui.available_height() * 0.45).clamp(170.0, 300.0);
                            egui::Panel::bottom("inspector")
                                .resizable(true)
                                .default_size(tall)
                                .min_size(80.0)
                                .show(ui, |ui| self.inspector(ui, i));
                            self.map_panel(ui, i, &pal);
                        }
                        (true, false) => self.map_panel(ui, i, &pal),
                        // without the map the inspector has the whole column
                        (false, _) => self.inspector(ui, i),
                    }
                });
        }
        egui::CentralPanel::default()
            .frame(egui::Frame::new().fill(pal.window_bg))
            .show(ui, |ui| {
                let highlight = self.highlight_ranges(i);
                let code = GENETIC_CODES[self.prefs.code].table;
                let env = editor::Env {
                    pal: &pal,
                    prefs: &self.prefs,
                    highlight: &highlight,
                    code,
                };
                editor::show(ui, &mut self.docs[i], &env);
            });
    }

    /// Enzymes shown on the map and highlighted in the text: the selected ones, plus unique cutters if asked.
    pub fn shown_enzymes(&mut self, i: usize) -> Vec<usize> {
        let counts: Vec<usize> = self.docs[i]
            .sites(&self.enzymes)
            .iter()
            .map(Vec::len)
            .collect();
        (0..self.enzymes.len())
            .filter(|&e| {
                let pick = self.selected.contains(&self.enzymes[e].name)
                    || (self.plus_unique && counts[e] == 1);
                pick && counts[e] > 0
            })
            .collect()
    }

    fn highlight_ranges(&mut self, i: usize) -> Vec<(usize, usize)> {
        if !self.docs[i].highlight_sites {
            return Vec::new();
        }
        let shown = self.shown_enzymes(i);
        let n = self.docs[i].doc.len();
        let enzymes = &self.enzymes;
        let sites = self.docs[i].sites(enzymes);
        let mut out = Vec::new();
        for e in shown {
            let len = enzymes[e].site.len();
            for s in &sites[e] {
                if s.start + len > n {
                    out.push((s.start, n));
                    out.push((0, s.start + len - n));
                } else {
                    out.push((s.start, s.start + len));
                }
            }
        }
        out
    }

    fn map_panel(&mut self, ui: &mut Ui, i: usize, pal: &Palette) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.map_linear, false, "Circular");
            ui.selectable_value(&mut self.map_linear, true, "Linear");
            ui.separator();
            ui.checkbox(&mut self.plus_unique, "+ unique cutters");
            if ui.small_button("Enzymes…").clicked() {
                self.enzyme_win.open = true;
            }
            super::figure::button(ui);
        });
        let shown = self.shown_enzymes(i);
        let enzymes = &self.enzymes;
        let dv = &mut self.docs[i];
        let all = dv.sites(enzymes).to_vec();
        let sites: Vec<MapSite> = shown
            .iter()
            .flat_map(|&e| {
                // a site with no cut (`Site::cuts_inside`) has no place on the map
                all[e]
                    .iter()
                    .filter(|s| s.cuts_inside)
                    .map(move |s| MapSite {
                        label: enzymes[e].name.clone(),
                        pos: s.top_cut,
                    })
            })
            .collect();
        let selection = dv.selection();
        let selected_feature = dv.current_feature();
        let size = ui.available_size();
        let linear = self.map_linear || !dv.doc.circular;
        let mut saved = String::new();
        let clicked = {
            let data = MapData {
                doc: &dv.doc,
                sites,
                selection,
                selected_feature,
                pal,
            };
            let draw = |ui: &mut Ui| {
                if linear {
                    map::linear(ui, &data, size)
                } else {
                    map::circular(ui, &data, size)
                }
            };
            if let Some(format) = super::figure::wanted() {
                let room = size.max(super::figure::MAP_ROOM);
                let draw_in = |ui: &mut Ui| {
                    if linear {
                        map::linear(ui, &data, room)
                    } else {
                        map::circular(ui, &data, room)
                    }
                };
                saved = super::figure::save(format, "map", room, |ui| {
                    draw_in(ui);
                });
            }
            draw(ui)
        };
        if let Some(hit) = clicked {
            dv.map_click(hit);
        }
        if !saved.is_empty() {
            self.say(saved);
        }
    }

    fn feature_table(&mut self, ui: &mut Ui, i: usize) {
        let dv = &mut self.docs[i];
        let mut pick: Option<usize> = None;
        ScrollArea::both()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                egui::Grid::new("feature_grid")
                    .striped(true)
                    .spacing([14.0, 2.0])
                    .show(ui, |ui| {
                        for h in ["", "Name", "Type", "Start", "End", "Length", "Strand"] {
                            ui.label(RichText::new(h).semibold());
                        }
                        ui.end_row();
                        for (fi, f) in dv.doc.features.iter().enumerate() {
                            let (rect, _) = ui
                                .allocate_exact_size(egui::vec2(14.0, 14.0), egui::Sense::hover());
                            ui.painter().rect_filled(rect, 2.0, theme::rgb(f.color()));
                            if ui
                                .selectable_label(dv.selected_feature == Some(fi), &f.name)
                                .clicked()
                            {
                                pick = Some(fi);
                            }
                            ui.label(&f.kind);
                            ui.label(RichText::new((f.start() + 1).to_string()).monospace());
                            ui.label(RichText::new(f.end().to_string()).monospace());
                            ui.label(RichText::new(f.len().to_string()).monospace());
                            ui.label(if f.reverse {
                                "◀ reverse"
                            } else {
                                "forward ▶"
                            });
                            ui.end_row();
                        }
                    });
            });
        if let Some(fi) = pick {
            dv.select_feature(fi);
        }
    }
}
