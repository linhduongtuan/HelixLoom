//! Toolbar icons drawn with the painter in the original's style: black outlines with a blue accent.

use super::theme::Palette;
use eframe::egui::{self, pos2, vec2, Color32, Pos2, Rect, Response, Sense, Shape, Stroke, Ui};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Icon {
    New,
    Open,
    Save,
    Cut,
    Copy,
    Paste,
    Undo,
    Redo,
    Find,
    Enzymes,
    /// The same with a frame round the middle: the original's picture while Selection Only is on (`enzymes_sonly`).
    EnzymesSelectionOnly,
    Map,
    Digest,
    Translate,
    Orf,
    Features,
    Primers,
    Align,
    Assemble,
}

/// A square icon button with a tooltip; returns the click response.
pub fn button(ui: &mut Ui, icon: Icon, tip: &str, pal: &Palette) -> Response {
    let (rect, resp) = ui.allocate_exact_size(vec2(30.0, 28.0), Sense::click());
    if resp.hovered() {
        let fill = if pal.dark {
            Color32::from_gray(70)
        } else {
            Color32::from_gray(214)
        };
        ui.painter().rect_filled(rect.shrink(1.0), 5.0, fill);
    }
    paint(ui, icon, rect.shrink2(vec2(5.0, 4.0)), pal);
    resp.on_hover_text(tip)
}

fn paint(ui: &Ui, icon: Icon, r: Rect, pal: &Palette) {
    let p = ui.painter();
    let ink = if pal.dark {
        Color32::WHITE
    } else {
        Color32::BLACK
    };
    let blue = pal.accent;
    let (k, b) = (Stroke::new(1.8, ink), Stroke::new(2.2, blue));
    let at = |x: f32, y: f32| pos2(r.left() + x * r.width(), r.top() + y * r.height());
    let seg = |a: Pos2, c: Pos2, s: Stroke| {
        p.line_segment([a, c], s);
    };
    let frame = |a: Pos2, c: Pos2, s: Stroke| {
        p.rect_stroke(Rect::from_two_pos(a, c), 3.0, s, egui::StrokeKind::Middle);
    };
    match icon {
        Icon::New => {
            frame(at(0.2, 0.05), at(0.8, 0.95), k);
            seg(at(0.5, 0.35), at(0.5, 0.75), b);
            seg(at(0.3, 0.55), at(0.7, 0.55), b);
        }
        Icon::Open => {
            frame(at(0.05, 0.3), at(0.95, 0.9), k);
            seg(at(0.05, 0.3), at(0.4, 0.3), b);
            seg(at(0.4, 0.3), at(0.5, 0.15), b);
            seg(at(0.5, 0.15), at(0.8, 0.15), b);
        }
        Icon::Save => {
            frame(at(0.1, 0.1), at(0.9, 0.9), k);
            frame(at(0.3, 0.1), at(0.7, 0.35), b);
            frame(at(0.3, 0.6), at(0.7, 0.9), b);
        }
        Icon::Cut => {
            seg(at(0.25, 0.1), at(0.7, 0.7), k);
            seg(at(0.75, 0.1), at(0.3, 0.7), k);
            p.circle_stroke(at(0.25, 0.82), r.width() * 0.14, b);
            p.circle_stroke(at(0.75, 0.82), r.width() * 0.14, b);
        }
        Icon::Copy => {
            frame(at(0.1, 0.1), at(0.6, 0.65), k);
            frame(at(0.4, 0.35), at(0.9, 0.92), b);
        }
        Icon::Paste => {
            frame(at(0.12, 0.15), at(0.88, 0.95), k);
            frame(at(0.35, 0.02), at(0.65, 0.28), b);
            seg(at(0.3, 0.55), at(0.7, 0.55), b);
            seg(at(0.3, 0.75), at(0.6, 0.75), b);
        }
        Icon::Undo | Icon::Redo => {
            let flip = if icon == Icon::Redo { -1.0 } else { 1.0 };
            let x = |v: f32| if flip > 0.0 { v } else { 1.0 - v };
            seg(at(x(0.3), 0.3), at(x(0.65), 0.3), k);
            seg(at(x(0.65), 0.3), at(x(0.9), 0.55), k);
            seg(at(x(0.9), 0.55), at(x(0.65), 0.85), k);
            seg(at(x(0.65), 0.85), at(x(0.3), 0.85), k);
            seg(at(x(0.3), 0.3), at(x(0.1), 0.3), b);
            seg(at(x(0.3), 0.3), at(x(0.3), 0.1), b);
        }
        Icon::Find => {
            p.circle_stroke(at(0.42, 0.42), r.width() * 0.32, k);
            seg(at(0.66, 0.66), at(0.95, 0.95), b);
        }
        Icon::Enzymes | Icon::EnzymesSelectionOnly => {
            seg(at(0.1, 0.25), at(0.9, 0.25), k);
            seg(at(0.1, 0.75), at(0.9, 0.75), k);
            seg(at(0.5, 0.0), at(0.5, 1.0), b);
            seg(at(0.35, 0.12), at(0.5, 0.0), b);
            seg(at(0.65, 0.12), at(0.5, 0.0), b);
            if icon == Icon::EnzymesSelectionOnly {
                frame(at(0.22, 0.08), at(0.78, 0.92), b);
            }
        }
        Icon::Map => {
            p.circle_stroke(r.center(), r.width() * 0.46, k);
            seg(at(0.5, 0.35), at(0.5, 0.65), k);
            seg(at(0.35, 0.5), at(0.65, 0.5), k);
            arc(p, r.center(), r.width() * 0.46, -2.4, -1.2, b);
            arc(p, r.center(), r.width() * 0.46, -0.6, 0.6, b);
        }
        Icon::Digest => {
            frame(at(0.05, 0.05), at(0.95, 0.95), k);
            seg(at(0.2, 0.3), at(0.45, 0.3), b);
            seg(at(0.2, 0.5), at(0.45, 0.5), b);
            seg(at(0.2, 0.7), at(0.45, 0.7), b);
            seg(at(0.6, 0.4), at(0.8, 0.4), b);
            seg(at(0.6, 0.6), at(0.8, 0.6), b);
        }
        Icon::Translate => {
            p.text(
                r.center(),
                egui::Align2::CENTER_CENTER,
                "Aa",
                egui::FontId::proportional(15.0),
                blue,
            );
            seg(at(0.05, 0.95), at(0.95, 0.95), k);
        }
        Icon::Orf => {
            seg(at(0.05, 0.5), at(0.8, 0.5), k);
            seg(at(0.8, 0.5), at(0.62, 0.3), k);
            seg(at(0.8, 0.5), at(0.62, 0.7), k);
            seg(at(0.15, 0.2), at(0.45, 0.2), b);
            seg(at(0.15, 0.8), at(0.5, 0.8), b);
        }
        Icon::Features => {
            p.add(Shape::closed_line(
                vec![
                    at(0.05, 0.3),
                    at(0.7, 0.3),
                    at(0.95, 0.5),
                    at(0.7, 0.7),
                    at(0.05, 0.7),
                ],
                k,
            ));
            p.circle_filled(at(0.2, 0.5), 2.0, blue);
        }
        Icon::Primers => {
            seg(at(0.05, 0.3), at(0.7, 0.3), k);
            seg(at(0.7, 0.3), at(0.55, 0.15), k);
            seg(at(0.95, 0.7), at(0.3, 0.7), b);
            seg(at(0.3, 0.7), at(0.45, 0.85), b);
        }
        Icon::Align => {
            seg(at(0.1, 0.2), at(0.9, 0.2), k);
            seg(at(0.1, 0.4), at(0.9, 0.4), b);
            seg(at(0.1, 0.65), at(0.9, 0.65), k);
            seg(at(0.1, 0.85), at(0.9, 0.85), b);
        }
        Icon::Assemble => {
            seg(at(0.05, 0.15), at(0.55, 0.15), k);
            seg(at(0.3, 0.35), at(0.8, 0.35), k);
            seg(at(0.5, 0.55), at(0.95, 0.55), k);
            seg(at(0.05, 0.85), at(0.95, 0.85), b);
            seg(at(0.5, 0.65), at(0.5, 0.78), b);
        }
    }
}

fn arc(p: &egui::Painter, c: Pos2, radius: f32, a0: f32, a1: f32, stroke: Stroke) {
    let steps = 12;
    let pts: Vec<Pos2> = (0..=steps)
        .map(|i| {
            let a = a0 + (a1 - a0) * i as f32 / steps as f32;
            c + vec2(a.cos(), a.sin()) * radius
        })
        .collect();
    p.add(Shape::line(pts, stroke));
}
