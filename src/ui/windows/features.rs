//! Features: add / edit a feature, delete or clear, annotate from a feature library.

use crate::seq::dna::Hit;
use crate::seq::features::{
    annotate, parse_library, type_colors, Feature, Rgb, DEFAULT_FWD, DEFAULT_REV,
};
use crate::ui::app::App;
use crate::ui::dialog::{self, Wanted};
use crate::ui::platform::route::Source;
use eframe::egui::{self, RichText};

/// Feature types offered in the editor (GenBank feature keys).
pub const FEATURE_TYPES: &[&str] = &[
    "gene",
    "CDS",
    "mRNA",
    "exon",
    "intron",
    "promoter",
    "terminator",
    "enhancer",
    "RBS",
    "rep_origin",
    "primer_bind",
    "protein_bind",
    "misc_binding",
    "misc_feature",
    "misc_recomb",
    "misc_signal",
    "polyA_signal",
    "repeat_region",
    "regulatory",
    "5'UTR",
    "3'UTR",
    "LTR",
    "source",
];

pub struct FeatureWin {
    pub open: bool,
    /// Making one feature for every highlighted find match instead of one for a chosen range.
    pub from_found: bool,
    /// Index of the feature being edited; `None` when adding a new one.
    pub editing: Option<usize>,
    pub name: String,
    pub kind: String,
    pub reverse: bool,
    pub start: usize,
    pub end: usize,
    pub fwd: Rgb,
    pub rev: Rgb,
    pub note: String,
}

impl Default for FeatureWin {
    fn default() -> Self {
        FeatureWin {
            open: false,
            from_found: false,
            editing: None,
            name: String::new(),
            kind: "misc_feature".into(),
            reverse: false,
            start: 1,
            end: 1,
            fwd: DEFAULT_FWD,
            rev: DEFAULT_REV,
            note: String::new(),
        }
    }
}

/// A feature on `[a, b)` with the name, type, colours and note typed into the window.
fn styled_feature(w: &FeatureWin, a: usize, b: usize, reverse: bool) -> Feature {
    let mut f = Feature::new(w.name.trim(), &w.kind, a, b, reverse);
    (f.fwd_color, f.rev_color) = (w.fwd, w.rev);
    if !w.note.trim().is_empty() {
        f.qualifiers
            .push(("note".into(), w.note.trim().to_string()));
    }
    f
}

/// The feature for one find match on a sequence of `n` bases (original `feature_from_found`).
fn feature_from_hit(w: &FeatureWin, hit: &Hit, n: usize) -> Feature {
    let end = hit.start + hit.len;
    let mut f = styled_feature(w, hit.start, end.min(n), hit.reverse);
    if end > n {
        // a match across the origin of a circle is one feature in two pieces
        f.spans = vec![(0, end - n), (hit.start, n)];
    }
    f
}

/// Features from `found` that the sequence doesn't already have (same name and place).
pub fn new_features(existing: &[Feature], found: Vec<Feature>) -> Vec<Feature> {
    found
        .into_iter()
        .filter(|f| {
            !existing
                .iter()
                .any(|e| e.name == f.name && e.spans == f.spans)
        })
        .collect()
}

impl App {
    /// Open the editor for a new feature (from the selection) or an existing one.
    pub fn open_feature_window(&mut self, edit: Option<usize>) {
        let Some(i) = self.current_doc() else { return };
        let dv = &self.docs[i];
        let mut w = FeatureWin {
            open: true,
            editing: edit,
            ..Default::default()
        };
        match edit.and_then(|k| dv.doc.features.get(k)) {
            Some(f) => {
                w.name = f.name.clone();
                w.kind = f.kind.clone();
                w.reverse = f.reverse;
                (w.start, w.end) = (f.start() + 1, f.end());
                (w.fwd, w.rev) = (f.fwd_color, f.rev_color);
                w.note = f
                    .qualifiers
                    .iter()
                    .find(|(k, _)| k == "note")
                    .map(|(_, v)| v.clone())
                    .unwrap_or_default();
            }
            None => {
                let (a, b) = dv
                    .selection()
                    .unwrap_or((dv.caret, (dv.caret + 1).min(dv.doc.len())));
                (w.start, w.end) = (a + 1, b.max(a + 1));
                (w.fwd, w.rev) = type_colors(&w.kind);
                w.name = format!("Feature {}", dv.doc.features.len() + 1);
            }
        }
        self.feature_win = w;
    }

    /// "New Features From Found Highlighting": ask for a name, type and colours, then make one
    /// feature for every highlighted match.
    pub fn open_found_feature_window(&mut self) {
        let Some(i) = self.current_doc() else { return };
        if self.docs[i].find_hits.is_empty() {
            self.say("There is no text highlighted from the find dialog.");
            return;
        }
        let kind = "misc_feature";
        let (fwd, rev) = type_colors(kind);
        self.feature_win = FeatureWin {
            open: true,
            from_found: true,
            name: "New Feature".into(),
            kind: kind.into(),
            fwd,
            rev,
            ..Default::default()
        };
    }

    fn apply_found_features(&mut self, i: usize) {
        let n = self.docs[i].doc.len();
        let hits = &self.docs[i].find_hits;
        // forward matches first, then reverse ones, as the original walks foundf then foundr
        let ordered = hits
            .iter()
            .filter(|h| !h.reverse)
            .chain(hits.iter().filter(|h| h.reverse));
        let found: Vec<Feature> = ordered
            .map(|h| feature_from_hit(&self.feature_win, h, n))
            .collect();
        let count = found.len();
        self.docs[i].edit(|d| d.features.extend(found));
        self.feature_win.open = false;
        self.say(format!(
            "Added {count} feature{} from the highlighted matches.",
            if count == 1 { "" } else { "s" }
        ));
    }

    fn apply_feature_window(&mut self) {
        let Some(i) = self.current_doc() else { return };
        if self.feature_win.from_found {
            return self.apply_found_features(i);
        }
        let w = &self.feature_win;
        let n = self.docs[i].doc.len();
        let (a, b) = ((w.start.max(1) - 1).min(n), w.end.min(n));
        if b <= a {
            self.say("The feature needs at least one base.");
            return;
        }
        let mut f = styled_feature(w, a, b, w.reverse);
        let editing = w.editing;
        self.docs[i].edit(|d| match editing.and_then(|k| d.features.get_mut(k)) {
            Some(slot) => {
                f.qualifiers = slot
                    .qualifiers
                    .iter()
                    .filter(|(k, _)| k != "note")
                    .cloned()
                    .chain(f.qualifiers.clone())
                    .collect();
                *slot = f;
            }
            None => d.features.push(f),
        });
        self.docs[i].selected_feature = editing;
        self.feature_win.open = false;
    }

    pub fn delete_selected_feature(&mut self) {
        let Some(i) = self.current_doc() else { return };
        match self.docs[i].selected_feature {
            Some(k) if k < self.docs[i].doc.features.len() => {
                self.docs[i].edit(|d| {
                    d.features.remove(k);
                });
                self.docs[i].selected_feature = None;
            }
            _ => self.say("Click a feature first to choose it."),
        }
    }

    pub fn clear_features(&mut self) {
        if let Some(i) = self.current_doc() {
            self.docs[i].edit(|d| d.features.clear());
            self.docs[i].selected_feature = None;
        }
    }

    pub fn annotate_with_library(&mut self) {
        let Some(i) = self.current_doc() else { return };
        let doc = &self.docs[i].doc;
        let fresh = new_features(
            &doc.features,
            annotate(&doc.seq, doc.circular, &self.library),
        );
        let n = fresh.len();
        self.docs[i].edit(|d| {
            d.features.extend(fresh);
            d.features.sort_by_key(|f| (f.start(), f.name.clone()));
        });
        self.say(format!(
            "Annotated: {n} new feature{} found using {} library entries.",
            if n == 1 { "" } else { "s" },
            self.library.len()
        ));
    }

    pub fn open_feature_library(&mut self) {
        let ask = dialog::FileDialog::new()
            .set_title("Open feature library")
            .add_filter("Feature libraries", &["txt", "tsv"]);
        if let Some(path) = dialog::pick_one(Wanted::FeatureLibrary, ask) {
            self.load_feature_library(Source::Path(&path));
        }
    }

    /// Read a feature library and annotate the open sequence with it.
    pub(in crate::ui) fn load_feature_library(&mut self, file: Source) {
        match file.text() {
            Ok(text) => {
                let lib = parse_library(&text);
                if lib.is_empty() {
                    self.say(format!(
                        "{} has no usable feature lines (name, tab, sequence, tab, type…).",
                        file.shown()
                    ));
                } else {
                    self.library = lib;
                    self.annotate_with_library();
                }
            }
            Err(e) => self.say(format!("Couldn't read {}: {e}", file.shown())),
        }
    }

    pub fn feature_window(&mut self, ctx: &egui::Context) {
        if !self.feature_win.open {
            return;
        }
        let max = self
            .current_doc()
            .map_or(1, |i| self.docs[i].doc.len().max(1));
        let mut open = true;
        let mut apply = false;
        let hits = self
            .current_doc()
            .map_or(0, |i| self.docs[i].find_hits.len());
        let title = if self.feature_win.from_found {
            "New Features From Found Highlighting"
        } else if self.feature_win.editing.is_some() {
            "Edit Feature"
        } else {
            "New Feature"
        };
        crate::ui::design::window(title)
            .open(&mut open)
            .collapsible(false)
            .default_width(380.0)
            .show(ctx, |ui| {
                let w = &mut self.feature_win;
                egui::Grid::new("feature_form")
                    .num_columns(2)
                    .spacing([10.0, 6.0])
                    .show(ui, |ui| {
                        ui.label("Name");
                        ui.add(egui::TextEdit::singleline(&mut w.name).desired_width(240.0));
                        ui.end_row();
                        ui.label("Type");
                        let before = w.kind.clone();
                        egui::ComboBox::from_id_salt("feature_type")
                            .selected_text(&w.kind)
                            .show_ui(ui, |ui| {
                                for t in FEATURE_TYPES {
                                    ui.selectable_value(&mut w.kind, t.to_string(), *t);
                                }
                            });
                        if w.kind != before {
                            (w.fwd, w.rev) = type_colors(&w.kind);
                        }
                        ui.end_row();
                        if w.from_found {
                            ui.label("Matches");
                            ui.label(format!(
                                "{hits} highlighted by Find, each on its own strand"
                            ));
                            ui.end_row();
                        } else {
                            ui.label("Strand");
                            ui.horizontal(|ui| {
                                ui.selectable_value(&mut w.reverse, false, "Forward ▶");
                                ui.selectable_value(&mut w.reverse, true, "◀ Reverse");
                            });
                            ui.end_row();
                            ui.label("From – to");
                            ui.horizontal(|ui| {
                                ui.add(egui::DragValue::new(&mut w.start).range(1..=max));
                                ui.label("–");
                                ui.add(egui::DragValue::new(&mut w.end).range(1..=max));
                            });
                            ui.end_row();
                        }
                        ui.label("Colours");
                        ui.horizontal(|ui| {
                            ui.label("forward");
                            egui::color_picker::color_edit_button_srgb(ui, &mut w.fwd);
                            ui.label("reverse");
                            egui::color_picker::color_edit_button_srgb(ui, &mut w.rev);
                        });
                        ui.end_row();
                        ui.label("Note");
                        ui.add(
                            egui::TextEdit::multiline(&mut w.note)
                                .desired_rows(2)
                                .desired_width(240.0),
                        );
                        ui.end_row();
                    });
                let label = if w.from_found {
                    format!("Add {hits} features")
                } else if w.editing.is_some() {
                    "Apply".to_string()
                } else {
                    "Add feature".to_string()
                };
                ui.horizontal(|ui| {
                    if ui.button(RichText::new(label).strong()).clicked() {
                        apply = true;
                    }
                });
            });
        if apply {
            self.apply_feature_window();
        } else if !open {
            self.feature_win.open = false;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn annotation_skips_features_already_present() {
        let existing = vec![Feature::new("T7 promoter", "promoter", 10, 30, false)];
        let found = vec![
            Feature::new("T7 promoter", "promoter", 10, 30, false),
            Feature::new("T7 promoter", "promoter", 100, 120, false),
            Feature::new("FLAG", "misc_feature", 10, 30, false),
        ];
        let fresh = new_features(&existing, found);
        assert_eq!(
            fresh.iter().map(|f| f.spans[0].0).collect::<Vec<_>>(),
            vec![100, 10]
        );
    }

    #[test]
    fn each_find_match_becomes_a_feature_on_its_own_strand() {
        let w = FeatureWin {
            name: " New Feature ".into(),
            kind: "misc_feature".into(),
            fwd: [1, 2, 3],
            rev: [4, 5, 6],
            note: "from find".into(),
            ..Default::default()
        };
        let fwd = feature_from_hit(
            &w,
            &Hit {
                start: 10,
                len: 6,
                reverse: false,
            },
            100,
        );
        assert_eq!(
            (fwd.name.as_str(), fwd.spans.clone(), fwd.reverse),
            ("New Feature", vec![(10, 16)], false)
        );
        assert_eq!((fwd.fwd_color, fwd.rev_color), ([1, 2, 3], [4, 5, 6]));
        assert_eq!(
            fwd.qualifiers,
            vec![("note".to_string(), "from find".to_string())]
        );
        let rev = feature_from_hit(
            &w,
            &Hit {
                start: 20,
                len: 6,
                reverse: true,
            },
            100,
        );
        assert!(rev.reverse);
        // a match over the origin of a circle is cut in two pieces
        let wrap = feature_from_hit(
            &w,
            &Hit {
                start: 97,
                len: 6,
                reverse: false,
            },
            100,
        );
        assert_eq!(wrap.spans, vec![(0, 3), (97, 100)]);
    }

    #[test]
    fn every_offered_type_has_colours_or_falls_back() {
        for t in FEATURE_TYPES {
            let (f, r) = type_colors(t);
            assert!(f != [0, 0, 0] || r != [0, 0, 0], "{t}");
        }
        assert!(FEATURE_TYPES.contains(&"CDS") && FEATURE_TYPES.contains(&"promoter"));
    }
}
