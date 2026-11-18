//! The mark of HelixLoom, drawn by code: two strands that twist around each other and are held together by rungs, the way a
//! loom holds its threads. One strand takes the colour of the program, the other its second colour.

use super::design::Look;
use eframe::egui::{pos2, Painter, Pos2, Rect, Shape, Stroke};
use std::f32::consts::PI;

/// How many straight pieces a strand is drawn in, and how many rungs join the strands.
const PIECES: usize = 24;
const RUNGS: usize = 5;
/// One full twist along the height of the mark, and how far from the middle a strand swings (a share of the width).
const TURNS: f32 = 1.0;
const SWING: f32 = 0.36;

/// Where strand `phase` (0 or π) is at height `t` (0 at the top, 1 at the bottom) inside `rect`.
fn strand_at(rect: Rect, t: f32, phase: f32) -> Pos2 {
    let angle = t * 2.0 * PI * TURNS + phase;
    pos2(
        rect.center().x + rect.width() * SWING * angle.sin(),
        rect.top() + rect.height() * t,
    )
}

/// Draw the mark in `rect`: first the stretches of each strand that run behind, fainter, then the rungs, then the stretches
/// in front, so that the strands read as passing over and under each other.
pub fn paint(painter: &Painter, rect: Rect, look: &Look) {
    let width = rect.width() * 0.1;
    let strands = [(0.0, look.accent), (PI, look.accent2)];
    let piece = |phase: f32, i: usize| {
        let (a, b) = (i as f32 / PIECES as f32, (i + 1) as f32 / PIECES as f32);
        let depth = ((a + b) * PI * TURNS + phase).cos();
        (
            [strand_at(rect, a, phase), strand_at(rect, b, phase)],
            depth,
        )
    };
    // the pieces on one side (behind or in front) that follow each other are one line: pieces drawn one by one would show
    // where they meet, the faint ones most
    let runs = |phase: f32, in_front: bool| -> Vec<Vec<Pos2>> {
        let mut runs: Vec<Vec<Pos2>> = Vec::new();
        let mut run: Vec<Pos2> = Vec::new();
        for i in 0..PIECES {
            let (points, depth) = piece(phase, i);
            if (depth >= 0.0) == in_front {
                if run.is_empty() {
                    run.push(points[0]);
                }
                run.push(points[1]);
            } else if !run.is_empty() {
                runs.push(std::mem::take(&mut run));
            }
        }
        if !run.is_empty() {
            runs.push(run);
        }
        runs
    };
    for &(phase, colour) in &strands {
        for run in runs(phase, false) {
            painter.add(Shape::line(
                run,
                Stroke::new(width, colour.gamma_multiply(0.45)),
            ));
        }
    }
    for k in 0..RUNGS {
        let t = (k as f32 + 0.5) / RUNGS as f32;
        let (left, right) = (strand_at(rect, t, 0.0), strand_at(rect, t, PI));
        painter.line_segment(
            [left, right],
            Stroke::new(width * 0.6, look.dim.gamma_multiply(0.55)),
        );
    }
    for &(phase, colour) in &strands {
        for run in runs(phase, true) {
            painter.add(Shape::line(run, Stroke::new(width, colour)));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_strands_start_in_the_middle_and_stay_inside_the_square() {
        let rect = Rect::from_min_max(pos2(0.0, 0.0), pos2(24.0, 24.0));
        for phase in [0.0, PI] {
            for i in 0..=PIECES {
                let p = strand_at(rect, i as f32 / PIECES as f32, phase);
                assert!(rect.expand(0.01).contains(p), "{p:?}");
            }
            let start = strand_at(rect, 0.0, phase);
            assert!((start.x - rect.center().x).abs() < 0.01);
        }
        // at the top one strand is left of the other only after the first quarter turn
        let (a, b) = (strand_at(rect, 0.125, 0.0), strand_at(rect, 0.125, PI));
        assert!(a.x > b.x);
    }
}

/// The sizes an icon can be asked in, in pixels each way.
pub const ICON_SIZES: (u32, u32) = (16, 2048);
/// The icon is drawn in a room of this many points, whatever its size in pixels.
const ICON_ROOM: f32 = 128.0;
/// How far the tile stands in from the edge of the icon, how round its corners are, and how far the mark stands in from the
/// tile (the proportions of the icons of a Mac).
const TILE_INSET: f32 = 12.5;
const TILE_RADIUS: f32 = 23.0;
const MARK_INSET: f32 = 20.0;

/// The mark on a light tile with round corners and nothing behind it: the picture of the app, `size` pixels each way.
pub fn icon(size: u32) -> eframe::egui::ColorImage {
    use super::appearance::Appearance;
    use super::figure::{self, Looks};
    use eframe::egui::{vec2, Color32, CornerRadius};
    let size = size.clamp(ICON_SIZES.0, ICON_SIZES.1);
    let looks = Looks {
        dark: false,
        appearance: Appearance::default(),
    };
    let look = looks.appearance.look(false);
    let drawn = figure::draw(
        looks,
        vec2(ICON_ROOM, ICON_ROOM),
        size as f32 / ICON_ROOM,
        |ui| {
            let tile = ui.max_rect().shrink(TILE_INSET);
            ui.painter()
                .rect_filled(tile, CornerRadius::same(TILE_RADIUS as u8), look.surface);
            paint(ui.painter(), tile.shrink(MARK_INSET), &look);
        },
    );
    drawn.image_on(Color32::TRANSPARENT)
}

#[cfg(test)]
mod icon_tests {
    use super::*;
    use eframe::egui::Color32;

    #[test]
    fn the_icon_is_the_mark_on_a_tile_with_nothing_in_its_corners_at_the_size_asked() {
        for size in [16, 64, 256] {
            let image = icon(size);
            assert_eq!(image.size, [size as usize, size as usize]);
            let at = |x: usize, y: usize| image.pixels[y * size as usize + x];
            // nothing in the corners, the tile in the middle of an edge of it, and the two colours of the mark somewhere
            assert_eq!(at(0, 0), Color32::TRANSPARENT, "{size}");
            assert_eq!(
                at(size as usize - 1, size as usize - 1),
                Color32::TRANSPARENT,
                "{size}"
            );
            let look = crate::ui::appearance::Appearance::default().look(false);
            assert_eq!(
                at(size as usize / 2, size as usize * 14 / 100),
                look.surface,
                "{size}"
            );
        }
        let image = icon(256);
        let look = crate::ui::appearance::Appearance::default().look(false);
        let near = |a: Color32, b: Color32| {
            (0..3).all(|k| (i16::from(a.to_array()[k]) - i16::from(b.to_array()[k])).abs() < 40)
        };
        assert!(
            image.pixels.iter().any(|p| near(*p, look.accent)),
            "the strand in the colour of the program"
        );
        assert!(
            image.pixels.iter().any(|p| near(*p, look.accent2)),
            "the strand in its second colour"
        );
        // a size that no icon has is the nearest that one has
        assert_eq!(icon(1).size, [16, 16]);
    }
}
