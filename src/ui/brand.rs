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
    for &(phase, colour) in &strands {
        for i in 0..PIECES {
            let (points, depth) = piece(phase, i);
            if depth < 0.0 {
                painter.line_segment(points, Stroke::new(width, colour.gamma_multiply(0.45)));
            }
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
        // the pieces in front that follow each other are one line
        let mut run: Vec<Pos2> = Vec::new();
        for i in 0..PIECES {
            let (points, depth) = piece(phase, i);
            if depth >= 0.0 {
                if run.is_empty() {
                    run.push(points[0]);
                }
                run.push(points[1]);
            } else if !run.is_empty() {
                painter.add(Shape::line(
                    std::mem::take(&mut run),
                    Stroke::new(width, colour),
                ));
            }
        }
        if !run.is_empty() {
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
