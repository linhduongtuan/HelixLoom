//! The bar along the top of the window: the mark and name of the program, its menus, the search box that opens the command
//! palette, where the heavy jobs go, whether something is working, the sheet of shortcuts and the theme.

use super::app::{Action, App};
use super::brand;
use super::design::{self, Look, BUTTON_RADIUS};
use super::platform::Workspace;
use super::state::ThemeMode;
use crate::remote::settings::Compute;
use eframe::egui::{
    self, text::LayoutJob, vec2, Align, Align2, Color32, CornerRadius, FontId, Layout, Painter,
    Pos2, Rect, Response, Sense, Stroke, StrokeKind, TextFormat, Ui,
};

/// The height of the bar and the side of the mark.
pub const HEIGHT: f32 = 46.0;
const MARK: f32 = 26.0;
/// The side of a square button, and the height of a pill.
const BUTTON: f32 = 30.0;
/// The search box when there is room for it, and the window width from which there is: the menus end about 600 points from the
/// left and the right-hand group with the box is about 500 wide.
const SEARCH_WIDTH: f32 = 290.0;
const WIDE: f32 = 1140.0;
/// Below this width the pills show no words.
const NARROW: f32 = 880.0;

impl App {
    pub(super) fn top_bar(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let look = design::look(self.prefs.dark);
        let width = ui.available_width();
        ui.horizontal_centered(|ui| {
            self.wordmark(ui, &look, out);
            ui.add_space(8.0);
            self.menu_bar(ui, out);
            ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                self.theme_button(ui, &look, out);
                if ui
                    .add_sized(
                        [BUTTON, BUTTON],
                        egui::Button::new("?").corner_radius(BUTTON_RADIUS),
                    )
                    .on_hover_text("Keyboard shortcuts")
                    .clicked()
                {
                    out.push(Action::Shortcuts);
                }
                self.jobs_pill(ui, &look, width < NARROW);
                self.compute_pill(ui, &look, out, width < NARROW);
                self.search_box(ui, &look, width >= WIDE);
            });
        });
    }

    /// The mark and the name: a click goes Home.
    fn wordmark(&self, ui: &mut Ui, look: &Look, out: &mut Vec<Action>) {
        let mut job = LayoutJob::default();
        let size = FontId::proportional(17.0);
        job.append("Helix", 0.0, TextFormat::simple(size.clone(), look.text));
        job.append("Loom", 0.0, TextFormat::simple(size, look.accent));
        let galley = ui.painter().layout_job(job);
        let (rect, response) =
            ui.allocate_exact_size(vec2(MARK + 8.0 + galley.size().x, MARK), Sense::click());
        brand::paint(
            ui.painter(),
            Rect::from_min_size(rect.min, vec2(MARK, MARK)),
            look,
        );
        let at = Pos2::new(
            rect.min.x + MARK + 8.0,
            rect.center().y - galley.size().y / 2.0,
        );
        ui.painter().galley(at, galley, look.text);
        if response.on_hover_text("Home").clicked() {
            out.push(Action::Go(Workspace::Home));
        }
    }

    /// The box that opens the command palette; just a magnifying glass when the window is not wide.
    fn search_box(&mut self, ui: &mut Ui, look: &Look, wide: bool) {
        let size = if wide {
            vec2(SEARCH_WIDTH, BUTTON)
        } else {
            vec2(BUTTON, BUTTON)
        };
        let field = search_field(ui, look, size, wide);
        if field
            .on_hover_text("Search every command and tool")
            .clicked()
        {
            self.palette.toggle();
        }
    }

    /// Where the heavy jobs go, with a dot that shows whether a workstation is connected; a click opens its page.
    fn compute_pill(&self, ui: &mut Ui, look: &Look, out: &mut Vec<Action>, short: bool) {
        let ws = &self.platform.workstation;
        let (text, dot) = match ws.compute {
            Compute::Local => ("This computer", look.dim),
            Compute::Workstation | Compute::Automatic => {
                let text = if ws.compute == Compute::Workstation {
                    "Workstation"
                } else {
                    "Automatic"
                };
                (text, if ws.answering() { look.good } else { look.warn })
            }
        };
        let hover = match ws.compute {
            Compute::Local => "Jobs run on this computer. Click to connect a workstation.",
            _ if ws.answering() => "The workstation is answering. Click to see it.",
            _ => "The saved workstation has not answered. Click to see why.",
        };
        if pill(ui, look, if short { "" } else { text }, Some(dot), None)
            .on_hover_text(hover)
            .clicked()
        {
            out.push(Action::Go(Workspace::Workstation));
        }
    }

    /// A spinner and a count while analyses are going on; nothing when none is.
    fn jobs_pill(&self, ui: &mut Ui, look: &Look, short: bool) {
        let running = self.running_analyses();
        if running == 0 {
            return;
        }
        let text = if short {
            String::new()
        } else {
            format!("{running} running")
        };
        let response = pill(ui, look, &text, None, Some(look.accent));
        response.on_hover_text("Analyses are running; the window stays free while they do");
    }

    /// Cycle the theme: the system's, light, dark.
    fn theme_button(&mut self, ui: &mut Ui, look: &Look, out: &mut Vec<Action>) {
        let (rect, response) = ui.allocate_exact_size(vec2(BUTTON, BUTTON), Sense::click());
        if response.hovered() {
            ui.painter()
                .rect_filled(rect, CornerRadius::same(BUTTON_RADIUS), look.hover);
        }
        theme_icon(ui.painter(), rect.center(), self.prefs.theme, look);
        let hover = format!(
            "Theme: {}. Click to change.",
            self.prefs.theme.label().to_lowercase()
        );
        if response.on_hover_text(hover).clicked() {
            out.push(Action::Theme(self.prefs.theme.next()));
        }
    }
}

/// The search field: a rounded box with a magnifying glass, a hint and, when `keys`, the key that opens the palette. It only
/// looks like a field; a click on it is what opens the palette. Just the glass when it is no wider than a button.
pub fn search_field(ui: &mut Ui, look: &Look, size: egui::Vec2, keys: bool) -> Response {
    let (rect, response) = ui.allocate_exact_size(size, Sense::click());
    let fill = if response.hovered() {
        look.hover
    } else {
        look.surface_alt
    };
    ui.painter().rect(
        rect,
        CornerRadius::same(BUTTON_RADIUS),
        fill,
        Stroke::new(1.0, look.border),
        StrokeKind::Inside,
    );
    let wide = size.x > BUTTON + 1.0;
    let glass_x = rect.left() + if wide { 16.0 } else { size.x / 2.0 };
    let glass = Rect::from_center_size(Pos2::new(glass_x, rect.center().y), vec2(14.0, 14.0));
    magnifier(ui.painter(), glass, look.dim);
    if wide {
        ui.painter().text(
            Pos2::new(rect.left() + 32.0, rect.center().y),
            Align2::LEFT_CENTER,
            "Search or run a command…",
            FontId::proportional(13.5),
            look.dim,
        );
    }
    if wide && keys {
        let mut at = ui.new_child(
            egui::UiBuilder::new()
                .max_rect(rect.shrink2(vec2(8.0, 0.0)))
                .layout(Layout::right_to_left(Align::Center)),
        );
        design::keycap(&mut at, look, &format!("{} K", design::command_key()));
    }
    response
}

/// A rounded label that can be clicked: a dot or a spinner, and words.
fn pill(
    ui: &mut Ui,
    look: &Look,
    text: &str,
    dot: Option<Color32>,
    spinner: Option<Color32>,
) -> Response {
    let galley = (!text.is_empty()).then(|| {
        ui.painter()
            .layout_no_wrap(text.to_string(), FontId::proportional(13.0), look.text)
    });
    let words = galley.as_ref().map_or(0.0, |g| g.size().x + 6.0);
    let (rect, response) =
        ui.allocate_exact_size(vec2(14.0 + 16.0 + words + 8.0, BUTTON), Sense::click());
    if response.hovered() {
        ui.painter()
            .rect_filled(rect, CornerRadius::same(BUTTON_RADIUS), look.hover);
    }
    let centre = Pos2::new(rect.left() + 18.0, rect.center().y);
    if let Some(colour) = dot {
        ui.painter().circle_filled(centre, 4.0, colour);
    }
    if let Some(colour) = spinner {
        let area = Rect::from_center_size(centre, vec2(14.0, 14.0));
        ui.put(area, egui::Spinner::new().size(14.0).color(colour));
    }
    if let Some(galley) = galley {
        let at = Pos2::new(rect.left() + 32.0, rect.center().y - galley.size().y / 2.0);
        ui.painter().galley(at, galley, look.text);
    }
    response
}

/// A magnifying glass: a ring and a handle.
fn magnifier(painter: &Painter, rect: Rect, colour: Color32) {
    let ring = Pos2::new(rect.center().x - 1.5, rect.center().y - 1.5);
    painter.circle_stroke(ring, rect.width() * 0.36, Stroke::new(1.6, colour));
    let from = ring + vec2(rect.width() * 0.26, rect.height() * 0.26);
    painter.line_segment([from, from + vec2(3.5, 3.5)], Stroke::new(1.8, colour));
}

/// The theme as a little picture: a half-filled disc for the system's, a sun for light, a moon for dark.
fn theme_icon(painter: &Painter, centre: Pos2, mode: ThemeMode, look: &Look) {
    let ink = look.text;
    match mode {
        ThemeMode::System => {
            painter.circle_stroke(centre, 7.0, Stroke::new(1.6, ink));
            // the right half filled
            let half: Vec<Pos2> = (0..=12)
                .map(|i| {
                    let a = -std::f32::consts::FRAC_PI_2 + std::f32::consts::PI * i as f32 / 12.0;
                    centre + vec2(7.0 * a.cos(), 7.0 * a.sin())
                })
                .collect();
            painter.add(egui::Shape::convex_polygon(half, ink, Stroke::NONE));
        }
        ThemeMode::Light => {
            painter.circle_filled(centre, 4.2, ink);
            for k in 0..8 {
                let a = std::f32::consts::FRAC_PI_4 * k as f32;
                let (c, s) = (a.cos(), a.sin());
                painter.line_segment(
                    [
                        centre + vec2(6.2 * c, 6.2 * s),
                        centre + vec2(8.2 * c, 8.2 * s),
                    ],
                    Stroke::new(1.6, ink),
                );
            }
        }
        ThemeMode::Dark => {
            painter.circle_filled(centre, 7.0, ink);
            // the second disc takes a bite out of the first, in the colour of the bar
            painter.circle_filled(centre + vec2(3.6, -3.0), 6.0, look.surface);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ui::app::Startup;
    use crate::ui::headless::frame_input;

    fn frame(width: f32, app: &mut App, ctx: &egui::Context, n: u32) {
        let mut input = frame_input(n);
        input.screen_rect = Some(Rect::from_min_size(Pos2::ZERO, vec2(width, 700.0)));
        let mut output = ctx.run_ui(input, |ui| app.draw(ui));
        output.textures_delta.clear();
    }

    #[test]
    fn the_bar_draws_at_every_width_in_both_themes_with_the_palette_open_or_shut() {
        let ctx = egui::Context::default();
        for dark in [false, true] {
            let mut app = App::new(Startup {
                dark,
                ..Startup::default()
            });
            for (n, width) in [1400.0, 1100.0, 900.0, 780.0].into_iter().enumerate() {
                app.palette = Default::default();
                frame(width, &mut app, &ctx, n as u32);
                app.palette.toggle();
                frame(width, &mut app, &ctx, n as u32 + 10);
            }
        }
    }
}
