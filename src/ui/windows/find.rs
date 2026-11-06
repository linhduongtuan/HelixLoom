//! Find (the original's mismatch-aware search), Jump To and Select From-To.

use crate::seq::dna::Hit;
use crate::seq::findmm::{find, Action, Outcome, Query};
use crate::seq::goto::{jump_to, select_from_to};
use crate::ui::app::App;
use eframe::egui::{self, Key, RichText};

/// The most mismatches the original's menu offers.
const MAX_MISMATCHES: usize = 25;
/// The choices of "Minimum 3' bases" besides All: 10 to 60 in steps of 5.
const MIN_3_BASES: [usize; 11] = [10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60];

pub struct FindWin {
    pub open: bool,
    pub text: String,
    pub case: bool,
    /// "also find rev-com of string"
    pub both: bool,
    /// "as literal"
    pub literal: bool,
    pub mismatches: usize,
    pub min_3_bases: Option<usize>,
    pub message: String,
    pub jump_open: bool,
    pub jump_text: String,
    pub range_open: bool,
    pub from_text: String,
    pub to_text: String,
}

impl Default for FindWin {
    fn default() -> Self {
        let q = Query::default();
        FindWin {
            open: false,
            text: String::new(),
            case: q.case_sensitive,
            both: q.rev_com,
            literal: q.literal,
            mismatches: q.mismatches,
            min_3_bases: q.min_3_bases,
            message: String::new(),
            jump_open: false,
            jump_text: String::new(),
            range_open: false,
            from_text: String::new(),
            to_text: String::new(),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Dir {
    Next,
    Prev,
}

impl App {
    fn find_query(&self) -> Query {
        let w = &self.find_win;
        Query {
            pattern: w.text.clone(),
            case_sensitive: w.case,
            rev_com: w.both,
            literal: w.literal,
            mismatches: w.mismatches,
            min_3_bases: w.min_3_bases,
            wrap: true,
        }
    }

    /// Run Find Next, Find Prev or Highlight All exactly as the original's `do_find_mm` does and show
    /// what it says. An empty box does nothing.
    fn find_run(&mut self, action: Action) {
        let Some(i) = self.current_doc() else { return };
        if self.find_win.text.is_empty() {
            return;
        }
        let q = self.find_query();
        let dv = &self.docs[i];
        let outcome = find(
            &dv.doc.seq,
            dv.doc.circular,
            &q,
            action,
            dv.caret,
            dv.selection().map(|(a, _)| a),
        );
        let dv = &mut self.docs[i];
        let n = dv.doc.len();
        self.find_win.message = match outcome {
            Outcome::Selected {
                start,
                end,
                message,
            } => {
                dv.select(start, end.min(n));
                dv.scroll_to = Some(start.min(n));
                message
            }
            Outcome::Highlighted { tags, message } => {
                dv.find_hits = tags
                    .iter()
                    .map(|t| Hit {
                        start: t.start,
                        len: t.end - t.start,
                        reverse: t.reverse,
                    })
                    .collect();
                message
            }
            Outcome::NotFound { message } => {
                if action == Action::Highlight {
                    dv.find_hits.clear();
                }
                message
            }
        };
    }

    /// Find next / previous occurrence and select it.
    pub fn find_step(&mut self, dir: Dir) {
        self.find_run(match dir {
            Dir::Next => Action::Next,
            Dir::Prev => Action::Prev,
        });
    }

    pub fn find_highlight_all(&mut self) {
        self.find_run(Action::Highlight);
    }

    pub fn find_window(&mut self, ctx: &egui::Context) {
        if self.find_win.open {
            let mut open = true;
            crate::ui::design::window("Find")
                .open(&mut open)
                .default_width(420.0)
                .show(ctx, |ui| {
                    let edit = ui.add(
                        egui::TextEdit::singleline(&mut self.find_win.text)
                            .hint_text("Sequence (IUPAC letters; | ; , between several)")
                            .desired_width(f32::INFINITY),
                    );
                    if edit.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                        self.find_step(Dir::Next);
                    }
                    ui.horizontal_wrapped(|ui| {
                        ui.checkbox(&mut self.find_win.case, "case sensitive");
                        ui.checkbox(&mut self.find_win.both, "also find rev-com of string");
                        ui.checkbox(&mut self.find_win.literal, "as literal");
                    });
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Allow mismatches:");
                        ui.add(
                            egui::DragValue::new(&mut self.find_win.mismatches)
                                .range(0..=MAX_MISMATCHES),
                        );
                        ui.label("Minimum 3' bases:");
                        egui::ComboBox::from_id_salt("find_min_3")
                            .selected_text(
                                self.find_win
                                    .min_3_bases
                                    .map_or("All".to_string(), |n| n.to_string()),
                            )
                            .show_ui(ui, |ui| {
                                ui.selectable_value(&mut self.find_win.min_3_bases, None, "All");
                                for n in MIN_3_BASES {
                                    ui.selectable_value(
                                        &mut self.find_win.min_3_bases,
                                        Some(n),
                                        n.to_string(),
                                    );
                                }
                            });
                    });
                    ui.horizontal(|ui| {
                        if ui.button("Find Next").clicked() {
                            self.find_step(Dir::Next);
                        }
                        if ui.button("Find Prev").clicked() {
                            self.find_step(Dir::Prev);
                        }
                        if ui.button("Highlight All").clicked() {
                            self.find_highlight_all();
                        }
                        if ui.button("Clear").clicked() {
                            if let Some(i) = self.current_doc() {
                                self.docs[i].find_hits.clear();
                            }
                            self.find_win.message.clear();
                        }
                    });
                    ui.label(RichText::new(&self.find_win.message).small());
                });
            self.find_win.open = open;
        }
        self.jump_windows(ctx);
    }

    /// Open Jump To the way the original does: its box starts with the caret's position.
    pub fn open_jump(&mut self) {
        let caret = self.current_doc().map_or(0, |i| self.docs[i].caret);
        self.find_win.jump_text = caret.to_string();
        self.find_win.jump_open = true;
    }

    /// Open Select From-To the way the original does: the selection's first and last base, or the caret's
    /// base and an empty "to".
    pub fn open_range(&mut self) {
        let (from, to) = match self.current_doc().map(|i| &self.docs[i]) {
            Some(dv) => match dv.selection() {
                Some((a, b)) => ((a + 1).to_string(), b.to_string()),
                None => ((dv.caret + 1).to_string(), String::new()),
            },
            None => (String::new(), String::new()),
        };
        self.find_win.from_text = from;
        self.find_win.to_text = to;
        self.find_win.range_open = true;
    }

    fn jump_windows(&mut self, ctx: &egui::Context) {
        if self.find_win.jump_open {
            let mut open = true;
            crate::ui::design::window("Jump To")
                .open(&mut open)
                .collapsible(false)
                .show(ctx, |ui| {
                    ui.label("Jump to:");
                    let edit = ui.text_edit_singleline(&mut self.find_win.jump_text);
                    // like the original's number boxes (check_char), only digits stay in
                    self.find_win.jump_text.retain(|c| c.is_ascii_digit());
                    let go = ui.button("Jump").clicked()
                        || (edit.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)));
                    if go {
                        if let Some(i) = self.current_doc() {
                            let dv = &mut self.docs[i];
                            dv.place_caret(jump_to(&self.find_win.jump_text));
                            dv.scroll_to = Some(dv.caret);
                        }
                    }
                });
            self.find_win.jump_open = open;
        }
        if self.find_win.range_open {
            let mut open = true;
            crate::ui::design::window("Select From-To")
                .open(&mut open)
                .collapsible(false)
                .show(ctx, |ui| {
                    ui.horizontal(|ui| {
                        ui.label("Select from");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.find_win.from_text)
                                .desired_width(70.0),
                        );
                        ui.label("to");
                        ui.add(
                            egui::TextEdit::singleline(&mut self.find_win.to_text)
                                .desired_width(70.0),
                        );
                    });
                    self.find_win.from_text.retain(|c| c.is_ascii_digit());
                    self.find_win.to_text.retain(|c| c.is_ascii_digit());
                    if ui.button("Select").clicked() {
                        if let Some(i) = self.current_doc() {
                            let dv = &mut self.docs[i];
                            let n = dv.doc.len();
                            let (a, b) =
                                select_from_to(&self.find_win.from_text, &self.find_win.to_text, n);
                            dv.select(a.min(b), b);
                            dv.scroll_to = Some(a.min(n));
                        }
                    }
                });
            self.find_win.range_open = open;
        }
    }
}
