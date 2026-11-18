//! A figure of a view (`docs/ROADMAP.md` 1.7): the view is drawn once more, alone and off the screen, by the very function that
//! draws it in the window, with the fonts and the colours of the window, and what it draws is saved as SVG (`svg.rs`: lines stay
//! lines and text stays text) or as a PNG with several pixels to a point. A page asks for one where it draws its view, with the
//! data it is drawing. Two lines for a page: `figure::button(ui)` among its buttons, and where the view is drawn
//! `if let Some(format) = figure::wanted() { message = figure::save(format, "map", size, |ui| { map::circular(ui, &data, size); }) }`.
//! What was chosen is kept for the frame it was chosen in (the button is drawn before the view), and so are the looks of the window.

mod svg;

use super::appearance::Appearance;
use super::raster::{Canvas, Textures};
use super::{dialog, fonts, headless, theme};
use eframe::egui::epaint::ClippedShape;
use eframe::egui::{
    self, Color32, ColorImage, Pos2, RawInput, Rect, Ui, Vec2, ViewportId, ViewportInfo,
};

pub const BUTTON: &str = "Save figure";
pub const SVG: &str = "SVG (lines and text, any size)…";
pub const PNG: &str = "PNG (4 pixels to a point)…";
/// How many frames a view is drawn for. Every view that the tests draw is right after one (egui measures and lays out a grid
/// within a frame); the second is to spare, for a view that takes a frame to settle, and no test needs it.
const FRAMES: u32 = 2;
/// The pixels a point gets in a PNG, and the most pixels a PNG has (a larger view gets fewer to a point).
pub const PNG_SCALE: f32 = 4.0;
const MOST_PIXELS: f32 = 24e6;
/// The sides a figure can have, in points.
const SIDES: (f32, f32) = (16.0, 6000.0);
const MAX_TEXTURE_SIDE: usize = 8192;
/// The least room the map is drawn in for a figure: in the small panel it has in the window, the names of many enzymes run off
/// its top, and a figure is for showing all of them.
pub const MAP_ROOM: Vec2 = Vec2::new(900.0, 700.0);

/// The kinds of picture a figure is saved as.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Svg,
    Png,
}

/// The looks of the window, which a figure is drawn with.
#[derive(Clone, Copy)]
pub struct Looks {
    pub dark: bool,
    pub appearance: Appearance,
}

thread_local! {
    /// The looks of the window in this frame, and the kind of figure chosen with a page's button in it.
    static LOOKS: std::cell::Cell<Option<Looks>> = const { std::cell::Cell::new(None) };
    static WANTED: std::cell::Cell<Option<Format>> = const { std::cell::Cell::new(None) };
}

/// A frame of the window begins: it looks like this, and what was asked for in the frame before is forgotten (a page whose
/// view was not drawn must not have its figure saved from another page later).
pub fn frame_begins(looks: Looks) {
    LOOKS.set(Some(looks));
    WANTED.set(None);
}

/// The kind of figure that was chosen with the button in this frame, once.
pub fn wanted() -> Option<Format> {
    WANTED.take()
}

fn figure_id() -> egui::Id {
    egui::Id::new("drawn for a figure")
}

/// Whether what is drawn here is for a figure: a view leaves out what only helps on the screen (a hint of what the pointer does).
pub fn is_figure(ctx: &egui::Context) -> bool {
    ctx.data(|d| d.get_temp(figure_id())).unwrap_or(false)
}

/// A view as it was drawn off the screen.
pub struct Figure {
    size: Vec2,
    ground: Color32,
    scale: f32,
    shapes: Vec<ClippedShape>,
    ctx: egui::Context,
    textures: Textures,
}

/// The size a figure is drawn at for a view of `size`: no side smaller or larger than a figure can have.
fn room(size: Vec2) -> Vec2 {
    let side = |x: f32| {
        if x.is_finite() {
            x.clamp(SIDES.0, SIDES.1).round()
        } else {
            SIDES.0
        }
    };
    Vec2::new(side(size.x), side(size.y))
}

/// The pixels to a point of a PNG of a figure of `size`: `PNG_SCALE`, or fewer for a figure so large that it would have too many
/// pixels.
pub fn png_scale(size: Vec2) -> f32 {
    let size = room(size);
    PNG_SCALE
        .min((MOST_PIXELS / (size.x * size.y)).sqrt())
        .max(1.0)
}

/// What `add` draws in a room of `size` points, with `scale` pixels to a point (which matters to a PNG: its text is made for that
/// many).
pub fn draw(looks: Looks, size: Vec2, scale: f32, mut add: impl FnMut(&mut Ui)) -> Figure {
    let size = room(size);
    let ctx = egui::Context::default();
    fonts::install(&ctx);
    // a figure has the sizes of the window at its plain size, whatever the window is zoomed to (the zoom is asked for with
    // the looks, and a zoom that was asked for cannot be taken back before the first frame)
    let plain = Appearance {
        scale: super::appearance::NORMAL_SCALE,
        ..looks.appearance
    };
    theme::apply(&ctx, looks.dark, false, plain);
    ctx.data_mut(|d| d.insert_temp(figure_id(), true));
    ctx.set_theme(if looks.dark {
        egui::Theme::Dark
    } else {
        egui::Theme::Light
    });
    let mut textures = Textures::default();
    let (mut freed, mut shapes) = (Vec::new(), Vec::new());
    for frame in 0..FRAMES {
        let input = RawInput {
            viewports: std::iter::once((
                ViewportId::ROOT,
                ViewportInfo {
                    native_pixels_per_point: Some(scale),
                    ..Default::default()
                },
            ))
            .collect(),
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, size)),
            max_texture_side: Some(MAX_TEXTURE_SIDE),
            time: Some(f64::from(frame) / 60.0),
            ..Default::default()
        };
        let mut out = ctx.run_ui(input, |ui| add(ui));
        textures.forget(&freed);
        freed = textures.apply(&mut out.textures_delta);
        shapes = out.shapes;
    }
    Figure {
        size,
        ground: looks.appearance.look(looks.dark).bg,
        scale,
        shapes,
        ctx,
        textures,
    }
}

impl Figure {
    pub fn svg(&self) -> String {
        svg::to_svg(&self.shapes, self.size, self.ground)
    }

    /// The figure in pixels: as many to a point as it was drawn for.
    pub fn image(&self) -> ColorImage {
        self.image_on(self.ground)
    }

    /// The figure in pixels on another ground than its own (none at all, for a picture with nothing behind it).
    pub fn image_on(&self, ground: Color32) -> ColorImage {
        let jobs = self.ctx.tessellate(self.shapes.clone(), self.scale);
        let pixels =
            [self.size.x, self.size.y].map(|points| (points * self.scale).round() as usize);
        let mut canvas = Canvas::new(pixels, ground);
        canvas.paint(&jobs, &self.textures, self.scale);
        canvas.into_image()
    }
}

/// The button of a page that saves its view as a figure; the kind chosen in its menu is what `wanted` gives in this frame.
pub fn button(ui: &mut Ui) {
    ui.menu_button(BUTTON, |ui| {
        for (label, format) in [(SVG, Format::Svg), (PNG, Format::Png)] {
            if ui.button(label).clicked() {
                WANTED.set(Some(format));
                ui.close();
            }
        }
    })
    .response
    .on_hover_text("The view as it is drawn here, without the buttons around it");
}

/// Draw what `add` draws as a figure of `size`, with the looks of the window, and save it under `name` with the ending of its
/// kind; what to tell the person (nothing when they cancelled).
pub fn save(format: Format, name: &str, size: Vec2, add: impl FnMut(&mut Ui)) -> String {
    let looks = LOOKS.get().unwrap_or(Looks {
        dark: false,
        appearance: Appearance::default(),
    });
    let (file, bytes) = match format {
        Format::Svg => (
            format!("{name}.svg"),
            Ok(draw(looks, size, 1.0, add).svg().into_bytes()),
        ),
        Format::Png => (
            format!("{name}.png"),
            headless::png_bytes(&draw(looks, size, png_scale(size), add).image()),
        ),
    };
    match bytes {
        Ok(bytes) => keep(&file, &bytes),
        Err(why) => format!("Could not make the picture: {why}"),
    }
}

/// Ask where to save `bytes` under `name` and write them; what to tell the person (nothing when they cancelled). In a test
/// there is no dialog to answer: what would be saved is kept for `take_saved`.
#[cfg(not(test))]
pub(in crate::ui) fn keep(name: &str, bytes: &[u8]) -> String {
    let ask = dialog::FileDialog::new().set_title(format!("Save {name}"));
    dialog::save_bytes(ask, name, bytes).map_or_else(String::new, |saved| saved.message)
}

#[cfg(test)]
thread_local! {
    /// What the tests' figures were saved as: a test has no dialog to answer.
    static SAVED: std::cell::RefCell<Vec<(String, Vec<u8>)>> = const { std::cell::RefCell::new(Vec::new()) };
}

#[cfg(test)]
pub(in crate::ui) fn keep(name: &str, bytes: &[u8]) -> String {
    let _ = dialog::IN_BROWSER;
    SAVED.with(|s| s.borrow_mut().push((name.to_string(), bytes.to_vec())));
    format!("Saved {name}.")
}

/// The figures saved on this thread since they were last taken: the name and the bytes of each.
#[cfg(test)]
pub fn take_saved() -> Vec<(String, Vec<u8>)> {
    SAVED.with(|s| std::mem::take(&mut *s.borrow_mut()))
}

#[cfg(test)]
mod pages;
#[cfg(test)]
mod tests;
