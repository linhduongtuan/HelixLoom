//! The Text Map dialog: a box for each kind of line the page can hold (`enz_text_map`, AppMain.tcl L19414-19579), in the order the
//! page puts them, each with a Show switch that turns its choices on or off. The original has the boxes dragged up and down; here
//! each has an Up and a Down button.

use super::TextMapWin;
use crate::seq::textmap::{Frame, LineNumbers, Reading, Settings, INDEX_CHARS};
use crate::ui::app::App;
use crate::ui::design::{self, Look};
use eframe::egui::{self, Key, Modifiers, ScrollArea, Ui};

/// The longest a line may be asked to be, and the shortest (the dialog's spin box, L19448).
const WIDTH_RANGE: std::ops::RangeInclusive<usize> = 10..=5000;
/// The blank lines between groups (L19435) and the most bases between two marks of the index line (L19501).
const SPACE_RANGE: std::ops::RangeInclusive<usize> = 1..=7;
const SPACING_RANGE: std::ops::RangeInclusive<usize> = 1..=5000;
/// What the dialog leaves free of the window's height (the bar above it, its title, the note and the buttons), and the least
/// height the boxes get.
const DIALOG_ROOM: f32 = 260.0;
const MIN_BOXES_HEIGHT: f32 = 160.0;

/// The order with `frame` one place earlier or later; the order as it was when it cannot move that way.
pub fn move_frame(order: &[Frame], frame: Frame, up: bool) -> Vec<Frame> {
    let mut next = order.to_vec();
    let Some(at) = next.iter().position(|&f| f == frame) else {
        return next;
    };
    let to = if up { at.checked_sub(1) } else { Some(at + 1) };
    match to {
        Some(to) if to < next.len() => next.swap(at, to),
        _ => {}
    }
    next
}

/// The switch of a box, for the boxes that have one (the DNA's box has none).
fn shown(s: &mut Settings, frame: Frame) -> Option<&mut bool> {
    match frame {
        Frame::Index => Some(&mut s.show_index),
        Frame::Enzymes => Some(&mut s.show_enzymes),
        Frame::Dna => None,
        Frame::Translation => Some(&mut s.show_translation),
        Frame::SecondStrand => Some(&mut s.second_strand),
        Frame::Genes => Some(&mut s.genes),
        Frame::Features => Some(&mut s.show_features),
        Frame::Space => Some(&mut s.show_space),
    }
}

impl App {
    pub(super) fn text_map_dialog_window(&mut self, ctx: &egui::Context, typing: bool) {
        if !self.text_map.dialog_open {
            return;
        }
        let look = design::look(self.prefs.dark);
        let (mut open, mut ok, mut cancel) = (true, false, false);
        design::window("Text Map")
            .open(&mut open)
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                (ok, cancel) = self.text_map_dialog_body(ui, &look, typing);
            });
        self.text_map.dialog_open = open && !cancel;
        if ok {
            self.text_map.dialog_open = false;
            self.make_text_map();
        }
    }

    /// The boxes and the buttons; (OK, Cancel) as pressed, or as Return and Escape ask for them.
    fn text_map_dialog_body(&mut self, ui: &mut Ui, look: &Look, typing: bool) -> (bool, bool) {
        let mut moved: Option<(Frame, bool)> = None;
        let w = &mut self.text_map;
        let order = w.settings.order.clone();
        ScrollArea::vertical()
            .max_height((ui.ctx().content_rect().height() - DIALOG_ROOM).max(MIN_BOXES_HEIGHT))
            .show(ui, |ui| {
                for (n, &frame) in order.iter().enumerate() {
                    let ends = (n == 0, n + 1 == order.len());
                    if let Some(request) = frame_box(ui, look, w, frame, ends) {
                        moved = Some(request);
                    }
                }
            });
        if let Some((frame, up)) = moved {
            w.settings.order = move_frame(&w.settings.order, frame, up);
        }
        let (mut ok, mut cancel) = (false, false);
        ui.label("Move a box with its Up and Down buttons to set the order of the lines.");
        ui.horizontal(|ui| {
            ok = design::primary_button(ui, look, "OK").clicked();
            cancel = design::ghost_button(ui, look, "Cancel").clicked();
        });
        if self.dialog_has_the_keys(ui, typing) {
            // the original binds Return to OK and its <<Cancel>> (Escape) to Cancel (L19417-19418, L19564)
            ok |= ui.input_mut(|i| i.consume_key(Modifiers::NONE, Key::Enter));
            cancel |= ui.input_mut(|i| i.consume_key(Modifiers::NONE, Key::Escape));
        }
        (ok, cancel)
    }
}

/// One box: its name, its Show switch and arrows, and its choices. A request to move the box is given back.
fn frame_box(
    ui: &mut Ui,
    look: &Look,
    w: &mut TextMapWin,
    frame: Frame,
    (first, last): (bool, bool),
) -> Option<(Frame, bool)> {
    let mut request = None;
    ui.group(|ui| {
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(frame.label()).strong().color(look.text));
            if let Some(on) = shown(&mut w.settings, frame) {
                ui.checkbox(on, "Show");
            }
            ui.add_enabled_ui(!first, |ui| {
                if ui
                    .small_button("Up")
                    .on_hover_text("Move this box up")
                    .clicked()
                {
                    request = Some((frame, true));
                }
            });
            ui.add_enabled_ui(!last, |ui| {
                if ui
                    .small_button("Down")
                    .on_hover_text("Move this box down")
                    .clicked()
                {
                    request = Some((frame, false));
                }
            });
        });
        let on = shown(&mut w.settings, frame).is_none_or(|s| *s);
        ui.add_enabled_ui(on, |ui| choices(ui, w, frame));
    });
    request
}

/// The choices inside a box.
fn choices(ui: &mut Ui, w: &mut TextMapWin, frame: Frame) {
    match frame {
        Frame::Enzymes => {
            ui.horizontal(|ui| {
                ui.radio_value(&mut w.selected_only, false, "All");
                ui.radio_value(&mut w.selected_only, true, "Selected Enzymes Only");
            });
        }
        Frame::Space => {
            ui.horizontal(|ui| {
                ui.label("Lines:");
                ui.add(egui::DragValue::new(&mut w.settings.space).range(SPACE_RANGE));
            });
        }
        Frame::Dna => dna_choices(ui, w),
        Frame::Translation => translation_choices(ui, &mut w.settings),
        Frame::Index => index_choices(ui, &mut w.settings),
        Frame::Features => {
            ui.checkbox(
                &mut w.settings.crop_features,
                "Crop Labels to Sequence Width",
            );
            ui.checkbox(&mut w.settings.show_hidden_features, "Show Hidden Features");
        }
        Frame::SecondStrand | Frame::Genes => {}
    }
}

fn dna_choices(ui: &mut Ui, w: &mut TextMapWin) {
    ui.horizontal(|ui| {
        ui.radio_value(&mut w.selection_only, false, "All");
        ui.radio_value(&mut w.selection_only, true, "Selected DNA only");
    });
    ui.checkbox(&mut w.settings.copy_highlight, "Copy Highlighting");
    ui.horizontal(|ui| {
        ui.label("Characters/line:");
        ui.add(egui::DragValue::new(&mut w.settings.width).range(WIDTH_RANGE));
        ui.label("Line Numbers:");
        egui::ComboBox::from_id_salt("text_map_numbers")
            .selected_text(w.settings.numbers.label())
            .show_ui(ui, |ui| {
                for n in LineNumbers::ALL {
                    ui.selectable_value(&mut w.settings.numbers, n, n.label());
                }
            });
    });
}

fn translation_choices(ui: &mut Ui, s: &mut Settings) {
    egui::Grid::new("text_map_readings").show(ui, |ui| {
        for (n, reading) in Reading::ALL.into_iter().enumerate() {
            ui.radio_value(&mut s.reading, reading, reading.label());
            if n % 3 == 2 {
                ui.end_row();
            }
        }
    });
}

fn index_choices(ui: &mut Ui, s: &mut Settings) {
    let mark = |ui: &mut Ui, id: &str, c: &mut char| {
        egui::ComboBox::from_id_salt(id)
            .selected_text(format!("'{c}'"))
            .width(48.0)
            .show_ui(ui, |ui| {
                for choice in INDEX_CHARS {
                    ui.selectable_value(c, choice, format!("'{choice}'"));
                }
            });
    };
    ui.horizontal(|ui| {
        mark(ui, "text_map_char1", &mut s.index_char1);
        ui.label("character every");
        ui.add(egui::DragValue::new(&mut s.index_spacing1).range(SPACING_RANGE));
        ui.label("bases.");
        ui.checkbox(&mut s.index_numbers, "Numbers");
    });
    ui.horizontal(|ui| {
        mark(ui, "text_map_char2", &mut s.index_char2);
        ui.label("character every");
        ui.add(egui::DragValue::new(&mut s.index_spacing2).range(SPACING_RANGE));
        ui.label("bases.");
    });
}
