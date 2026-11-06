//! The workspace icons, drawn with the painter in the platform's two inks: the page's text colour and the accent.

use super::Workspace;
use crate::ui::design;
use eframe::egui::{pos2, vec2, Color32, Painter, Pos2, Rect, Shape, Stroke, StrokeKind};
use std::f32::consts::TAU;

/// Paint the icon of `ws` inside `r`.
pub fn paint(p: &Painter, ws: Workspace, r: Rect, ink: Color32, accent: Color32, dark: bool) {
    let at = |x: f32, y: f32| pos2(r.left() + x * r.width(), r.top() + y * r.height());
    let thin = Stroke::new(1.8, ink);
    let bold = Stroke::new(2.5, accent);
    let line = |a: Pos2, b: Pos2, s: Stroke| {
        p.line_segment([a, b], s);
    };
    let path = |points: Vec<Pos2>, s: Stroke| {
        p.add(Shape::line(points, s));
    };
    let curve = |f: &dyn Fn(f32) -> (f32, f32)| -> Vec<Pos2> {
        (0..=32)
            .map(|i| {
                let (x, y) = f(i as f32 / 32.0);
                at(x, y)
            })
            .collect()
    };
    match ws {
        Workspace::Home => {
            path(vec![at(0.08, 0.5), at(0.5, 0.12), at(0.92, 0.5)], bold);
            path(
                vec![at(0.2, 0.46), at(0.2, 0.88), at(0.8, 0.88), at(0.8, 0.46)],
                thin,
            );
            path(
                vec![
                    at(0.42, 0.88),
                    at(0.42, 0.62),
                    at(0.58, 0.62),
                    at(0.58, 0.88),
                ],
                bold,
            );
        }
        Workspace::Sequences => {
            path(
                curve(&|t| (0.08 + 0.84 * t, 0.5 + 0.3 * (t * TAU * 1.5).sin())),
                thin,
            );
            path(
                curve(&|t| (0.08 + 0.84 * t, 0.5 - 0.3 * (t * TAU * 1.5).sin())),
                bold,
            );
            for k in 0..6 {
                let t = 0.08 + 0.168 * k as f32;
                let (x, s) = (0.08 + 0.84 * t, 0.3 * (t * TAU * 1.5).sin());
                line(
                    at(x, 0.5 + s),
                    at(x, 0.5 - s),
                    Stroke::new(1.2, ink.gamma_multiply(0.6)),
                );
            }
        }
        Workspace::Signals => {
            line(
                at(0.05, 0.85),
                at(0.95, 0.85),
                Stroke::new(1.2, ink.gamma_multiply(0.6)),
            );
            for (k, (base, height)) in [(b'A', 0.7), (b'C', 0.5), (b'G', 0.8), (b'T', 0.6)]
                .into_iter()
                .enumerate()
            {
                let centre = 0.2 + 0.2 * k as f32;
                let colour = design::trace_color(base, dark);
                path(
                    curve(&|t| {
                        let u = (t - 0.5) * 4.0;
                        (centre - 0.13 + 0.26 * t, 0.85 - height * (-u * u).exp())
                    }),
                    Stroke::new(2.0, colour),
                );
            }
        }
        Workspace::Alignment => {
            for (k, (from, to)) in [(0.08, 0.92), (0.08, 0.7), (0.3, 0.92), (0.2, 0.8)]
                .into_iter()
                .enumerate()
            {
                let y = 0.16 + 0.2 * k as f32;
                let colour = if k % 2 == 0 {
                    accent
                } else {
                    ink.gamma_multiply(0.7)
                };
                p.rect_filled(
                    Rect::from_min_max(at(from, y), at(to, y + 0.1)),
                    2.0,
                    colour,
                );
            }
        }
        Workspace::Phylogeny => {
            path(vec![at(0.06, 0.5), at(0.3, 0.5)], thin);
            path(vec![at(0.3, 0.25), at(0.3, 0.75)], thin);
            path(vec![at(0.3, 0.25), at(0.58, 0.25)], bold);
            path(vec![at(0.3, 0.75), at(0.58, 0.75)], thin);
            for (y, tips) in [(0.25, [0.14, 0.36]), (0.75, [0.64, 0.86])] {
                path(vec![at(0.58, tips[0]), at(0.58, tips[1])], thin);
                for tip in tips {
                    path(
                        vec![at(0.58, tip), at(0.92, tip)],
                        if y < 0.5 { bold } else { thin },
                    );
                }
            }
        }
        Workspace::Trees => {
            let centre = at(0.5, 0.5);
            let radius = r.width() * 0.4;
            for k in 0..8 {
                let a = TAU * (k as f32 + 0.5) / 8.0;
                let dir = vec2(a.cos(), a.sin());
                line(
                    centre + dir * radius * 0.25,
                    centre + dir * radius,
                    if k % 3 == 0 { bold } else { thin },
                );
            }
            for (a0, a1) in [(0.0, 1.0), (2.0, 3.0), (5.0, 6.0)] {
                let arc: Vec<Pos2> = (0..=8)
                    .map(|i| {
                        let a = TAU * (a0 + 0.5 + (a1 - a0) * i as f32 / 8.0) / 8.0;
                        centre + vec2(a.cos(), a.sin()) * radius * 0.55
                    })
                    .collect();
                path(arc, thin);
            }
        }
        Workspace::Structure => {
            path(
                curve(&|t| (0.08 + 0.78 * t, 0.46 + 0.26 * (t * TAU * 2.0).sin())),
                bold,
            );
            path(
                curve(&|t| (0.08 + 0.78 * t, 0.58 + 0.26 * (t * TAU * 2.0).sin())),
                Stroke::new(1.4, ink.gamma_multiply(0.7)),
            );
            p.add(Shape::convex_polygon(
                vec![at(0.84, 0.3), at(0.96, 0.46), at(0.84, 0.62)],
                accent,
                Stroke::NONE,
            ));
        }
        Workspace::Pipelines => {
            for k in 0..3 {
                let x = 0.04 + 0.34 * k as f32;
                p.rect_stroke(
                    Rect::from_min_max(at(x, 0.36), at(x + 0.24, 0.64)),
                    3.0,
                    thin,
                    StrokeKind::Middle,
                );
            }
            for k in 0..2 {
                let x = 0.28 + 0.34 * k as f32;
                line(at(x, 0.5), at(x + 0.1, 0.5), bold);
                p.add(Shape::convex_polygon(
                    vec![at(x + 0.07, 0.42), at(x + 0.12, 0.5), at(x + 0.07, 0.58)],
                    accent,
                    Stroke::NONE,
                ));
            }
        }
        Workspace::Workstation => {
            // a rack of three machines, the middle one at work
            for k in 0..3 {
                let y = 0.1 + 0.3 * k as f32;
                let busy = k == 1;
                p.rect_stroke(
                    Rect::from_min_max(at(0.08, y), at(0.92, y + 0.22)),
                    3.0,
                    if busy { bold } else { thin },
                    StrokeKind::Middle,
                );
                line(
                    at(0.2, y + 0.11),
                    at(0.5, y + 0.11),
                    Stroke::new(1.2, ink.gamma_multiply(0.6)),
                );
                p.circle_filled(at(0.78, y + 0.11), 2.2, if busy { accent } else { ink });
            }
        }
    }
}
