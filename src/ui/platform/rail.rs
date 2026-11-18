//! The navigation rail: one button per workspace, down the left edge.

use super::{icons, Workspace, RAIL_WIDTH};
use crate::ui::app::{Action, App};
use crate::ui::design::{Look, CARD_RADIUS};
use crate::ui::tour::Target;
use eframe::egui::{self, vec2, Align2, CornerRadius, FontId, Rect, Sense, Ui};

const BUTTON_HEIGHT: f32 = 64.0;
/// A short window squeezes the buttons together, down to this height, rather than push the last ones out of sight.
const MIN_BUTTON_HEIGHT: f32 = 52.0;
const BUTTON_GAP: f32 = 2.0;
const TOP_SPACE: f32 = 10.0;
/// Room kept at the bottom of the rail.
const FOOT_ROOM: f32 = 12.0;
/// The bar that marks the workspace that is open: how wide and how tall.
const MARK_WIDTH: f32 = 3.0;
const MARK_HEIGHT: f32 = 26.0;
const ICON_SIZE: f32 = 28.0;

/// How tall each of `count` buttons is when `room` points of height are there to share, with `spacing` between two of them (the
/// space that the window puts between things, which a denser or a more spacious setting changes).
fn button_height(room: f32, count: usize, spacing: f32) -> f32 {
    let each = (room - TOP_SPACE - FOOT_ROOM) / count as f32 - BUTTON_GAP - spacing;
    each.clamp(MIN_BUTTON_HEIGHT, BUTTON_HEIGHT)
}

impl App {
    pub(in crate::ui) fn rail(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let look = self.prefs.look();
        let height = button_height(
            ui.available_height(),
            Workspace::ALL.len(),
            ui.spacing().item_spacing.y,
        );
        // a shorter button moves its icon up so the label does not touch it
        let lift = (BUTTON_HEIGHT - height) * 0.25;
        // when even the shortest buttons do not fit (a window zoomed far in), the rail scrolls instead of cutting the last ones off
        egui::ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| self.rail_buttons(ui, out, &look, height, lift));
        self.tour.aim(Target::Rail, ui.max_rect());
    }

    fn rail_buttons(
        &mut self,
        ui: &mut Ui,
        out: &mut Vec<Action>,
        look: &Look,
        height: f32,
        lift: f32,
    ) {
        ui.add_space(TOP_SPACE);
        for ws in Workspace::ALL {
            let selected = self.workspace == ws;
            let (rect, resp) =
                ui.allocate_exact_size(vec2(RAIL_WIDTH - 14.0, height), Sense::click());
            let p = ui.painter();
            if selected || resp.hovered() {
                let fill = if selected {
                    look.accent_soft
                } else {
                    look.surface_alt
                };
                p.rect_filled(rect, CornerRadius::same(CARD_RADIUS), fill);
            }
            if selected {
                // a bar at the left of the button, in the colour of the program
                let bar = Rect::from_center_size(
                    egui::pos2(rect.left() + 3.0, rect.center().y),
                    vec2(MARK_WIDTH, MARK_HEIGHT),
                );
                p.rect_filled(bar, CornerRadius::same(2), look.accent);
            }
            let ink = if selected { look.accent } else { look.dim };
            let icon = Rect::from_center_size(
                rect.center_top() + vec2(0.0, 10.0 - lift + ICON_SIZE / 2.0),
                vec2(ICON_SIZE, ICON_SIZE),
            );
            icons::paint(p, ws, icon, ink, look.accent, look);
            p.text(
                rect.center_bottom() - vec2(0.0, 9.0),
                Align2::CENTER_CENTER,
                ws.label(),
                FontId::proportional(11.5),
                if selected { look.accent } else { look.text },
            );
            if resp.on_hover_text(ws.blurb()).clicked() {
                out.push(Action::Go(ws));
            }
            ui.add_space(BUTTON_GAP);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The space between two buttons that the comfortable setting gives.
    const SPACING: f32 = 6.0;

    #[test]
    fn buttons_share_a_short_rail_but_never_shrink_past_their_label() {
        let n = Workspace::ALL.len();
        assert_eq!(button_height(900.0, n, SPACING), BUTTON_HEIGHT);
        let squeezed = button_height(600.0, n, SPACING);
        assert!(
            squeezed > MIN_BUTTON_HEIGHT && squeezed < BUTTON_HEIGHT,
            "{squeezed}"
        );
        assert_eq!(button_height(200.0, n, SPACING), MIN_BUTTON_HEIGHT);
        // what is chosen fits: every button, the gaps and the theme button, whenever the rail is tall enough to allow it
        let room = 600.0;
        for spacing in [SPACING, 3.0] {
            let used = TOP_SPACE
                + FOOT_ROOM
                + n as f32 * (button_height(room, n, spacing) + BUTTON_GAP + spacing);
            assert!(used <= room + 0.01, "{used} with {spacing} between");
        }
    }
}
