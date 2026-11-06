//! The strip under the sequence: mouse position, length, selection start / length / end with
//! codon frame, ORF summary, Tm and GC, and the circular / linear and Dam/Dcm toggles.

use super::state::{DocView, Prefs};
use super::theme::Palette;
use crate::seq::dna::{gc_percent, melting_temp, orf_summary};
use eframe::egui::{self, RichText, Ui};

/// `12<0>`: a position or length followed by its codon frame, as the original prints it.
pub fn with_frame(value: usize, frame: usize) -> String {
    format!("{value}<{frame}>")
}

pub struct Readout {
    pub start_label: &'static str,
    pub start: String,
    pub length: String,
    pub end: String,
    pub orf: String,
    pub tm: String,
    pub gc: String,
}

pub fn readout(dv: &DocView, prefs: &Prefs) -> Readout {
    match dv.selection() {
        None => Readout {
            start_label: "insert@",
            start: with_frame(dv.caret + 1, dv.caret % 3),
            length: String::new(),
            end: String::new(),
            orf: String::new(),
            tm: String::new(),
            gc: String::new(),
        },
        Some((a, b)) => {
            let sel = dv.doc.slice(a, b);
            let tm = match melting_temp(sel, prefs.dna_nm, prefs.salt_mm) {
                Some(t) => format!("{t:.1}° C"),
                None => "---".to_string(),
            };
            Readout {
                start_label: "Start",
                start: with_frame(a + 1, a % 3),
                length: with_frame(b - a, (b - a) % 3),
                end: with_frame(b, (b - 1) % 3),
                orf: orf_summary(sel),
                tm,
                gc: format!("{}%", gc_percent(sel).round()),
            }
        }
    }
}

fn field(ui: &mut Ui, pal: &Palette, label: &str, value: &str, width: f32) {
    ui.allocate_ui(egui::vec2(width, 34.0), |ui| {
        ui.vertical_centered(|ui| {
            ui.spacing_mut().item_spacing.y = 0.0;
            ui.label(RichText::new(label).small().color(pal.dim));
            ui.label(RichText::new(value).monospace());
        });
    });
}

pub fn show(ui: &mut Ui, dv: &mut DocView, prefs: &Prefs, pal: &Palette) {
    let r = readout(dv, prefs);
    let hover = dv
        .hover
        .map(|h| with_frame(h + 1, h % 3))
        .unwrap_or_default();
    let total = dv.doc.len().to_string();
    ui.horizontal(|ui| {
        field(ui, pal, "Mouse", &hover, 70.0);
        field(ui, pal, "Sequence", &total, 80.0);
        field(ui, pal, r.start_label, &r.start, 80.0);
        field(
            ui,
            pal,
            if r.length.is_empty() { "" } else { "Length" },
            &r.length,
            80.0,
        );
        field(
            ui,
            pal,
            if r.end.is_empty() { "" } else { "End" },
            &r.end,
            80.0,
        );
        field(
            ui,
            pal,
            if r.orf.is_empty() { "" } else { "ORF" },
            &r.orf,
            110.0,
        );
        field(
            ui,
            pal,
            if r.tm.is_empty() { "" } else { "Tm" },
            &r.tm,
            80.0,
        );
        field(
            ui,
            pal,
            if r.gc.is_empty() { "" } else { "%GC" },
            &r.gc,
            56.0,
        );
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            let text = if dv.doc.circular {
                "circular"
            } else {
                "linear"
            };
            if ui
                .button(text)
                .on_hover_text("Click to switch between circular and linear")
                .clicked()
            {
                dv.doc.circular = !dv.doc.circular;
                dv.doc.dirty = true;
            }
            ui.checkbox(&mut dv.doc.methylated, "Dam/Dcm")
                .on_hover_text("Methylation state of the DNA.\nSet for plasmids, unset for PCR.");
        });
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::doc::Doc;

    #[test]
    fn nothing_selected_shows_only_the_insert_position() {
        let mut dv = DocView::new(Doc::new("t", b"ATGAAATAAC".to_vec()));
        dv.place_caret(4);
        let r = readout(&dv, &Prefs::default());
        assert_eq!((r.start_label, r.start.as_str()), ("insert@", "5<1>"));
        assert!(r.length.is_empty() && r.tm.is_empty() && r.orf.is_empty());
    }

    #[test]
    fn a_selection_reports_frames_orf_tm_and_gc() {
        let mut dv = DocView::new(Doc::new("t", b"CCATGGCGTTTGAACGGCTCAAGG".to_vec()));
        dv.select(2, 22);
        let r = readout(&dv, &Prefs::default());
        assert_eq!(r.start, "3<2>");
        assert_eq!(r.length, "20<2>");
        assert_eq!(r.end, "22<0>");
        assert_eq!(r.tm, "59.3° C");
        assert_eq!(r.gc, "50%");
        assert!(r.orf.starts_with('M'));
    }

    #[test]
    fn long_selections_have_no_melting_temperature() {
        let mut dv = DocView::new(Doc::new("t", b"ACGT".repeat(40)));
        dv.select(0, 120);
        assert_eq!(readout(&dv, &Prefs::default()).tm, "---");
    }
}
