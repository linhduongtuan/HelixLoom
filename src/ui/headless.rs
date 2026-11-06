//! The `--screenshot` picture, drawn without a window.
//!
//! The system stops drawing a window while the screen is locked (and a machine may have no display
//! at all), and then the window never gets as far as taking its picture. This draws the same frames
//! with egui alone, fills the triangles in with `raster`, and saves the result. `watch` starts it
//! when the window is too slow, and `main` when the window cannot be opened.

use super::app::{App, Startup};
use super::raster::{Canvas, Textures};
use eframe::egui::{
    self, Color32, ColorImage, FullOutput, Pos2, RawInput, Rect, ViewportId, ViewportInfo,
};
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;

/// The size of the window `main` opens, in points.
pub const WINDOW_POINTS: [f32; 2] = [1280.0, 800.0];
/// Two pixels to the point, like a Retina screen; the picture from the window has the same scale there.
const PIXELS_PER_POINT: f32 = 2.0;
/// Frames drawn before the picture is taken.
pub const SCREENSHOT_FRAMES: u32 = 120;
/// One frame's time at 60 a second, so background jobs (an alignment, a tree) get as long as in the window.
const FRAME_PACE: Duration = Duration::from_millis(16);
/// How long the window has to take its picture before it is taken without it.
const WINDOW_PATIENCE: Duration = Duration::from_secs(20);
/// About what a graphics card reports as its biggest texture.
const MAX_TEXTURE_SIDE: usize = 8192;

static PICTURE_CLAIMED: AtomicBool = AtomicBool::new(false);

/// Whoever is about to write the picture asks first, so the window and the fallback never both do.
pub fn claim_picture() -> bool {
    !PICTURE_CLAIMED.swap(true, Ordering::SeqCst)
}

pub(super) fn frame_input(frame: u32) -> RawInput {
    let screen = ViewportInfo {
        native_pixels_per_point: Some(PIXELS_PER_POINT),
        ..Default::default()
    };
    RawInput {
        viewports: std::iter::once((ViewportId::ROOT, screen)).collect(),
        screen_rect: Some(Rect::from_min_size(
            Pos2::ZERO,
            egui::vec2(WINDOW_POINTS[0], WINDOW_POINTS[1]),
        )),
        max_texture_side: Some(MAX_TEXTURE_SIDE),
        time: Some(f64::from(frame) / 60.0),
        ..Default::default()
    }
}

/// The window as it looks after `frames` frames of at least `pace` each, drawn on the processor.
/// There is no pointer, so nothing is hovered or scrolled.
fn draw(start: Startup, frames: u32, pace: Duration) -> ColorImage {
    let ctx = egui::Context::default();
    let mut app = App::new(start);
    let mut textures = Textures::default();
    let mut freed = Vec::new();
    let (mut shapes, mut scale) = (Vec::new(), PIXELS_PER_POINT);
    for frame in 0..frames.max(1) {
        let FullOutput {
            mut textures_delta,
            shapes: drawn,
            pixels_per_point,
            ..
        } = ctx.run_ui(frame_input(frame), |ui| app.draw(ui));
        textures.forget(&freed);
        freed = textures.apply(&mut textures_delta);
        (shapes, scale) = (drawn, pixels_per_point);
        std::thread::sleep(pace);
    }
    let jobs = ctx.tessellate(shapes, scale);
    let mut canvas = Canvas::new(
        WINDOW_POINTS.map(|points| (points * scale).round() as usize),
        Color32::BLACK,
    );
    canvas.paint(&jobs, &textures, scale);
    canvas.into_image()
}

pub fn save_png(path: &str, img: &ColorImage) -> Result<(), String> {
    let file = std::fs::File::create(path).map_err(|e| e.to_string())?;
    let mut enc = png::Encoder::new(
        std::io::BufWriter::new(file),
        img.size[0] as u32,
        img.size[1] as u32,
    );
    enc.set_color(png::ColorType::Rgba);
    enc.set_depth(png::BitDepth::Eight);
    let mut w = enc.write_header().map_err(|e| e.to_string())?;
    w.write_image_data(img.as_raw()).map_err(|e| e.to_string())
}

/// Draw the picture with no window, save it, and end the program: 0 when it was saved, 1 when not.
// ponytail: ends the program with `exit` instead of closing the window politely, which is fine for
// a one-shot picture; a window that is still open is simply dropped with the process.
pub fn finish(start: Startup, path: &str) -> ! {
    let code = match save_png(path, &draw(start, SCREENSHOT_FRAMES, FRAME_PACE)) {
        Ok(()) => 0,
        Err(e) => {
            eprintln!("could not save the screenshot: {e}");
            1
        }
    };
    std::process::exit(code)
}

/// Keep an eye on the window: if it has not taken the picture within `WINDOW_PATIENCE` (the screen
/// is locked, or the machine is far too busy), take it without the window.
pub fn watch(start: Startup, path: String) {
    std::thread::spawn(move || {
        std::thread::sleep(WINDOW_PATIENCE);
        if claim_picture() {
            eprintln!(
                "The window drew nothing in {} s; drawing the picture without it.",
                WINDOW_PATIENCE.as_secs()
            );
            finish(start, &path);
        }
    });
}

#[cfg(test)]
mod tests;
