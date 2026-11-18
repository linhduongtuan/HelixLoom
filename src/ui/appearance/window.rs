//! The Appearance window: the theme and every setting of `appearance`, each changing the window as soon as it is chosen.

use super::{Accent, Appearance, Colours, Contrast, Density, MIN_TARGET, SCALES};
use crate::ui::app::App;
use crate::ui::design::{self, Look};
use crate::ui::state::ThemeMode;
use eframe::egui::{
    self, Align2, Context, FontId, Response, RichText, Sense, Stroke, Ui, WidgetInfo, WidgetType,
};

/// What is kept free of the window's height, for the title and the space above the window.
const SCREEN_MARGIN: f32 = 140.0;
/// The least height that the contents of the window scroll in.
const MIN_ROOM: f32 = 120.0;
/// How wide the window is.
const WIDTH: f32 = 500.0;
/// The rings drawn round a colour sample: the round sample itself, the ring of the one chosen (or under the pointer), and the one of the
/// sample that has the keyboard.
const SAMPLE_RADIUS: f32 = 8.5;
const CHOSEN_RING: f32 = 11.0;
pub(super) const FOCUS_RING: f32 = 12.5;
/// The size of a colour sample: a click target of the least size.
const SWATCH: f32 = MIN_TARGET;
/// The cells of the preview: wide enough for a letter, as tall as a row of the alignment.
const CELL: egui::Vec2 = egui::vec2(26.0, 22.0);

/// One choice among several, drawn as a button that is lit when it is the one chosen. Returns whether it was clicked.
fn choice<T: PartialEq + Copy>(ui: &mut Ui, current: &mut T, value: T, label: &str) -> bool {
    let clicked = ui.selectable_label(*current == value, label).clicked();
    if clicked {
        *current = value;
    }
    clicked
}

/// A round sample of an accent, ringed when chosen; its name is the tooltip and what a screen reader says.
fn swatch(ui: &mut Ui, look: &Look, accent: Accent, chosen: bool) -> Response {
    let (rect, response) = ui.allocate_exact_size(egui::vec2(SWATCH, SWATCH), Sense::click());
    let painter = ui.painter();
    painter.circle_filled(rect.center(), SAMPLE_RADIUS, accent.swatch(look.dark));
    if chosen {
        painter.circle_stroke(rect.center(), CHOSEN_RING, Stroke::new(2.0, look.text));
    } else if response.hovered() {
        painter.circle_stroke(rect.center(), CHOSEN_RING, Stroke::new(1.0, look.dim));
    }
    // the swatch that has the keyboard (Tab goes through them, Enter and Space choose) is told by a ring outside the others
    if response.has_focus() {
        painter.circle_stroke(rect.center(), FOCUS_RING, Stroke::new(2.0, look.text));
    }
    response.widget_info(|| {
        WidgetInfo::selected(WidgetType::RadioButton, true, chosen, accent.label())
    });
    response.on_hover_text(accent.label())
}

/// Bases and one amino acid of each group as the alignment would colour them, so that a choice of colours can be judged at once.
fn preview(ui: &mut Ui, look: &Look) {
    let rows: [(&[u8], bool); 2] = [(b"ACGT", false), (b"LKDSCGPY", true)];
    for (letters, protein) in rows {
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 2.0;
            for &letter in letters {
                let fill = if protein {
                    design::residue_fill(letter, look)
                } else {
                    design::base_fill(letter, look)
                }
                .unwrap_or(look.surface_alt);
                let (rect, _) = ui.allocate_exact_size(CELL, Sense::hover());
                ui.painter().rect_filled(rect, 3.0, fill);
                ui.painter().text(
                    rect.center(),
                    Align2::CENTER_CENTER,
                    char::from(letter),
                    FontId::monospace(13.0),
                    design::letter_ink(look, Some(fill)),
                );
            }
        });
    }
}

impl App {
    /// What the window says about whether the choices are kept.
    pub(super) fn appearance_note(&self) -> &'static str {
        if self.appearance_locked {
            "Each choice is used at once. This run was started with settings on the command line, so none of them is kept."
        } else if !self.keeper.has_settings() {
            "Each choice is used at once. There is no settings folder, so it is not kept."
        } else if self.keeper.appearance_failed && self.keeper.kept_in_a_browser() {
            "Each choice is used at once. The last one could not be kept (the browser's storage for this site cannot be written); it is tried again every few seconds."
        } else if self.keeper.appearance_failed {
            "Each choice is used at once. The last one could not be kept (the settings folder cannot be written); it is tried again every few seconds."
        } else if self.keeper.kept_in_a_browser() {
            "Each choice is used at once and kept in this browser for the next visit."
        } else {
            "Each choice is used at once and kept for the next time."
        }
    }

    pub(crate) fn appearance_window(&mut self, ctx: &Context) {
        if !self.appearance_open {
            return;
        }
        let look = self.prefs.look();
        let note = self.appearance_note();
        let (mut theme, mut set) = (self.prefs.theme, self.prefs.appearance);
        let mut open = true;
        // all of the screen but a margin, whatever panels there are on it
        let room = (ctx.content_rect().height() - SCREEN_MARGIN).max(MIN_ROOM);
        design::window("Appearance")
            .open(&mut open)
            .default_width(WIDTH)
            // the most that the contents may take before they scroll (the window shrinks to them when they are shorter)
            .default_height(room)
            // over the panels too, as far as the screen goes (the window is shorter than the area between them)
            .constrain_to(ctx.content_rect())
            .collapsible(false)
            .show(ctx, |ui| {
              // a window zoomed far in is shorter than this one: it scrolls, so that every choice can still be reached
              egui::ScrollArea::vertical().max_height(room).auto_shrink([false, true]).show(ui, |ui| {
                ui.label(RichText::new(note).color(look.dim));

                design::section(ui, &look, "Theme");
                ui.horizontal_wrapped(|ui| {
                    for mode in ThemeMode::ALL {
                        choice(ui, &mut theme, mode, mode.label());
                    }
                });

                design::section(ui, &look, "Size");
                ui.horizontal_wrapped(|ui| {
                    for scale in SCALES {
                        choice(ui, &mut set.scale, scale, &format!("{scale}%"));
                    }
                });
                ui.label(
                    RichText::new(format!(
                        "{0} + and {0} − make everything larger or smaller from anywhere, {0} 0 sets it back.",
                        design::command_key()
                    ))
                    .size(12.0)
                    .color(look.dim),
                );

                design::section(ui, &look, "Spacing");
                ui.horizontal_wrapped(|ui| {
                    for density in Density::ALL {
                        choice(ui, &mut set.density, *density, density.label());
                    }
                });

                design::section(ui, &look, "Accent colour");
                ui.horizontal(|ui| {
                    for accent in Accent::ALL {
                        if swatch(ui, &look, *accent, set.accent == *accent).clicked() {
                            set.accent = *accent;
                        }
                    }
                    ui.label(RichText::new(set.accent.label()).color(look.dim));
                });

                design::section(ui, &look, "Contrast");
                ui.horizontal_wrapped(|ui| {
                    for contrast in Contrast::ALL {
                        choice(ui, &mut set.contrast, *contrast, contrast.label());
                    }
                });

                design::section(ui, &look, "Colours of bases and residues");
                ui.horizontal_wrapped(|ui| {
                    for colours in Colours::ALL {
                        choice(ui, &mut set.colours, *colours, colours.label());
                    }
                });
                // the sample uses what is chosen now, not what the window was drawn with
                preview(ui, &set.look(look.dark));

                ui.add_space(8.0);
                if design::ghost_button(ui, &look, "Back to the defaults").clicked() {
                    set = Appearance::default();
                }
              });
            });
        self.prefs.theme = theme;
        self.prefs.appearance = set;
        self.appearance_open = open;
    }
}
