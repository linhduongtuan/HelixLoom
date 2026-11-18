//! The inspector: a panel under the map that says what the selected feature, or the selected bases, are.
//! The map, the feature table and the sequence view all select the same thing; this is where it is read.

use super::app::App;
use super::design;
use super::theme;
use crate::seq::dna::GENETIC_CODES;
use crate::seq::inspect::{inspect, Inspection, Subject, MOST_PROTEIN};
use crate::ui::fonts::FontWeight;
use eframe::egui::{self, RichText, ScrollArea, Ui};

/// The most characters of a value that may share its line with another: the panel is narrow, and short values are paired to
/// leave room for the protein and the overlaps below.
const SHORT: usize = 12;

/// The table's rows put in lines: two short ones side by side, anything longer alone.
fn in_lines(rows: Vec<(&'static str, String)>) -> Vec<Vec<(&'static str, String)>> {
    let mut lines: Vec<Vec<(&'static str, String)>> = Vec::new();
    for row in rows {
        let short = row.1.chars().count() <= SHORT;
        match lines.last_mut() {
            Some(last) if short && last.len() == 1 && last[0].1.chars().count() <= SHORT => {
                last.push(row);
            }
            _ => lines.push(vec![row]),
        }
    }
    lines
}

/// The least height of a line of the panel's tables: the panel is small, and the default is tall enough for a button.
const ROW_HEIGHT: f32 = 15.0;

impl App {
    /// What the inspector shows for document `i`, as of now.
    pub(super) fn inspection(&self, i: usize) -> Option<Inspection> {
        let dv = &self.docs[i];
        inspect(
            &dv.doc,
            dv.current_feature(),
            dv.selection(),
            self.prefs.dna_nm,
            self.prefs.salt_mm,
            GENETIC_CODES[self.prefs.code].table,
        )
    }

    pub(super) fn inspector(&mut self, ui: &mut Ui, i: usize) {
        let look = self.prefs.look();
        let found = self.inspection(i);
        let mut pick: Option<usize> = None;
        ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                design::section(ui, &look, "Inspector");
                let Some(ins) = &found else {
                    ui.label(RichText::new("Nothing selected.").color(look.text));
                    ui.label(
                        RichText::new(
                            "Click a feature on the map, in the table or in the sequence, \
                             or drag across some bases, and what they are appears here.",
                        )
                        .color(look.dim),
                    );
                    return;
                };
                self.inspector_head(ui, i, ins);
                egui::Grid::new("inspector_rows")
                    .spacing([12.0, 2.0])
                    .min_row_height(ROW_HEIGHT)
                    .show(ui, |ui| {
                        for line in in_lines(ins.rows()) {
                            for (label, value) in line {
                                ui.label(RichText::new(label).small().color(look.dim));
                                ui.label(RichText::new(value).monospace());
                            }
                            ui.end_row();
                        }
                    });
                if let Some(p) = &ins.protein {
                    design::section(ui, &look, "Protein");
                    let more = if p.residues > MOST_PROTEIN { "…" } else { "" };
                    ui.add(
                        egui::Label::new(RichText::new(format!("{}{more}", p.shown)).monospace())
                            .wrap(),
                    );
                }
                if !ins.qualifiers.is_empty() {
                    design::section(ui, &look, "Notes");
                    egui::Grid::new("inspector_notes")
                        .spacing([12.0, 2.0])
                        .min_row_height(ROW_HEIGHT)
                        .show(ui, |ui| {
                            for (key, value) in &ins.qualifiers {
                                ui.label(RichText::new(key).small().color(look.dim));
                                ui.add(egui::Label::new(value.as_str()).wrap());
                                ui.end_row();
                            }
                        });
                }
                if !ins.overlaps.is_empty() {
                    design::section(ui, &look, "Overlaps");
                    for &oi in &ins.overlaps {
                        let f = &self.docs[i].doc.features[oi];
                        ui.horizontal(|ui| {
                            let (rect, _) = ui
                                .allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
                            ui.painter().rect_filled(rect, 2.0, theme::rgb(f.color()));
                            if ui
                                .link(format!("{} ({})", f.name, f.kind))
                                .on_hover_text("Select it")
                                .clicked()
                            {
                                pick = Some(oi);
                            }
                        });
                    }
                }
            });
        if let Some(fi) = pick {
            self.docs[i].select_feature(fi);
        }
    }

    /// The name, with the feature's own colour beside it.
    fn inspector_head(&self, ui: &mut Ui, i: usize, ins: &Inspection) {
        ui.horizontal(|ui| {
            if let Subject::Feature(fi) = ins.subject {
                let (rect, _) =
                    ui.allocate_exact_size(egui::vec2(14.0, 14.0), egui::Sense::hover());
                let color = self.docs[i].doc.features[fi].color();
                ui.painter().rect_filled(rect, 3.0, theme::rgb(color));
            }
            ui.label(RichText::new(&ins.title).semibold().size(15.0));
        });
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod layout_tests {
    use super::*;

    fn row(label: &'static str, value: &str) -> (&'static str, String) {
        (label, value.to_string())
    }

    #[test]
    fn short_values_share_a_line_and_long_ones_have_their_own() {
        let lines = in_lines(vec![
            row("Type", "CDS"),
            row("Location", "327..902"),
            row("Length", "576 bp · 192 codons"),
            row("Strand", "forward"),
            row("GC", "52%"),
            row("Tm", "4.7° C"),
            row("ORF", "M>191*/-"),
        ]);
        let labels: Vec<Vec<&str>> = lines
            .iter()
            .map(|l| l.iter().map(|(k, _)| *k).collect())
            .collect();
        assert_eq!(
            labels,
            [
                vec!["Type", "Location"],
                vec!["Length"],
                vec!["Strand", "GC"],
                vec!["Tm", "ORF"]
            ]
        );
    }

    #[test]
    fn a_long_value_is_never_squeezed_beside_another() {
        let lines = in_lines(vec![
            row("Type", "CDS"),
            row("Location", "complement(join(1..6,11..13))"),
            row("Length", "9 bp"),
        ]);
        let labels: Vec<Vec<&str>> = lines
            .iter()
            .map(|l| l.iter().map(|(k, _)| *k).collect())
            .collect();
        assert_eq!(labels, [vec!["Type"], vec!["Location"], vec!["Length"]]);
    }

    #[test]
    fn no_rows_make_no_lines() {
        assert!(in_lines(Vec::new()).is_empty());
    }
}
