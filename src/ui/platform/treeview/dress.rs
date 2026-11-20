//! What the viewer adds to a tree beyond its lines and names: coloured and widened branches, labels, shapes and bars from the
//! annotations, blocks of colour behind clades, triangles for collapsed clades, tips that a search found, an axis, and a legend.

use super::{Frame, MARGIN};
use crate::seq::treelayout::{Laid, Layout};
use crate::seq::treestyle::{self, Axis, Decor, Legend, Rgb};
use crate::ui::design::Look;
use eframe::egui::{pos2, vec2, Align2, Color32, FontId, Painter, Pos2, Rect, Shape, Stroke};
use std::f32::consts::TAU;

/// What goes on the tree besides the tree itself; all of it optional.
#[derive(Clone, Copy, Default)]
pub struct Dress<'a> {
    /// Colours, widths, labels, shapes and bars for each node of the layout.
    pub decor: Option<&'a Decor>,
    /// Clades (a node of the layout) with a block of colour behind them, or on a round layout drawn in that colour.
    pub highlights: &'a [(usize, Rgb)],
    /// Tips (nodes of the layout) that a search found.
    pub found: &'a [usize],
    pub axis: Option<Axis>,
    /// What each node's annotations say, for the tip of the pointer's hover text.
    pub notes: &'a [String],
    /// A tip's name is drawn in the colour of its branch.
    pub tip_colours: bool,
}

pub(super) fn colour32(c: Rgb) -> Color32 {
    Color32::from_rgb(c[0], c[1], c[2])
}

fn tinted(c: Rgb, alpha: u8) -> Color32 {
    Color32::from_rgba_unmultiplied(c[0], c[1], c[2], alpha)
}

impl Dress<'_> {
    /// The decoration, when it is for this layout.
    pub(super) fn decor_for(&self, laid: &Laid) -> Option<&Decor> {
        self.decor.filter(|d| d.colour.len() == laid.nodes.len())
    }

    /// The colour of the nearest highlighted clade that `i` is in.
    pub(super) fn block_colour(&self, laid: &Laid, i: usize) -> Option<Rgb> {
        self.highlights
            .iter()
            .filter(|(h, _)| *h < laid.nodes.len() && laid.is_within(i, *h))
            .map(|&(h, c)| (laid.nodes[h].origin, c))
            .max_by_key(|&(origin, _)| origin)
            .map(|(_, c)| c)
    }
}

/// Blocks of colour behind the highlighted clades of a layout in rows. (A round layout draws them in the strokes instead.)
pub(super) fn blocks(
    painter: &Painter,
    laid: &Laid,
    f: &Frame,
    dress: &Dress,
    name_room: &dyn Fn(usize) -> f32,
) {
    if f.circular {
        return;
    }
    for &(h, colour) in dress
        .highlights
        .iter()
        .filter(|(h, _)| *h < laid.nodes.len())
    {
        let tips = laid.tips_under(h);
        let (Some(first), Some(last)) = (tips.first(), tips.last()) else {
            continue;
        };
        let rows = |t: usize| laid.nodes[t].cartoon.map_or(0.5, |c| c.rows / 2.0);
        let top = f.at(0.0, laid.nodes[*first].y - rows(*first)).y;
        let bottom = f.at(0.0, laid.nodes[*last].y + rows(*last)).y;
        let left = f.at(laid.nodes[h].x, 0.0).x - 4.0;
        let right = tips
            .iter()
            .map(|&t| {
                f.at(
                    laid.nodes[t].cartoon.map_or(laid.nodes[t].x, |c| c.far),
                    0.0,
                )
                .x + name_room(t)
            })
            .fold(left, f32::max);
        painter.rect_filled(
            Rect::from_min_max(pos2(left, top), pos2(right, bottom)),
            3.0,
            tinted(colour, 46),
        );
    }
}

/// A triangle for each collapsed clade, with its name or its number of tips.
pub(super) fn cartoons(
    painter: &Painter,
    look: &Look,
    laid: &Laid,
    f: &Frame,
    dress: &Dress,
    font: &FontId,
) {
    let deco = dress.decor_for(laid);
    for (i, node) in laid.nodes.iter().enumerate() {
        let Some(c) = node.cartoon else { continue };
        let apex = f.at(node.x, node.y);
        let (a, b, middle) = if f.circular {
            // the wedge of the circle that its rows take, at the angle of the node
            let half = TAU * (c.rows / laid.rows_total.max(1.0)) / 2.0 * 0.9;
            let (direction, reach) = (TAU * node.angle, c.far - node.radius);
            let corner = |d: f32| -> Pos2 {
                let (cos, sin) = ((direction + d).cos(), (direction + d).sin());
                f.at(node.x + cos * reach, node.y + sin * reach)
            };
            let tip = f.at(
                node.x + direction.cos() * reach,
                node.y + direction.sin() * reach,
            );
            if laid.kind == Layout::Circular {
                // round layout: the far corners lie on the circle of radius `far`, where the tips are
                let on =
                    |d: f32| f.at(c.far * (direction + d).cos(), c.far * (direction + d).sin());
                (on(-half), on(half), on(0.0))
            } else {
                (corner(-half), corner(half), tip)
            }
        } else {
            let half = (c.rows / 2.0 - 0.1).max(0.1);
            (
                f.at(c.far, node.y - half),
                f.at(c.far, node.y + half),
                f.at(c.far, node.y),
            )
        };
        let tone = deco
            .and_then(|d| d.colour[i])
            .or_else(|| dress.block_colour(laid, i))
            .map_or(look.dim, colour32);
        painter.add(Shape::convex_polygon(
            vec![apex, a, b],
            tone.gamma_multiply(0.35),
            Stroke::new(1.2, tone),
        ));
        let name = if node.name.is_empty() {
            format!("{} tips", node.tips)
        } else {
            node.name.clone()
        };
        let (anchor, shift) = if f.circular && node.x < 0.0 {
            (Align2::RIGHT_CENTER, -6.0)
        } else {
            (Align2::LEFT_CENTER, 6.0)
        };
        painter.text(
            middle + vec2(shift, 0.0),
            anchor,
            name,
            font.clone(),
            look.text,
        );
    }
}

/// The axis under a tree in rows, with lines up from its ticks when asked for. `bottom` is the y of the line.
pub(super) fn axis(painter: &Painter, look: &Look, laid: &Laid, f: &Frame, a: &Axis, rect: Rect) {
    if laid.kind != Layout::Phylogram {
        return;
    }
    let deepest = laid.extent().1 .0;
    let y = rect.bottom() - 24.0;
    let (x0, x1) = (f.at(0.0, 0.0).x, f.at(deepest, 0.0).x);
    painter.line_segment([pos2(x0, y), pos2(x1, y)], Stroke::new(1.2, look.text));
    let wide = ((x1 - x0) / 70.0).floor().max(2.0) as usize;
    for (at, text) in a.ticks(deepest, wide) {
        let x = f.at(at, 0.0).x;
        if a.grid {
            painter.line_segment(
                [pos2(x, rect.top() + MARGIN), pos2(x, y)],
                Stroke::new(1.0, look.border.gamma_multiply(0.6)),
            );
        }
        painter.line_segment([pos2(x, y), pos2(x, y + 5.0)], Stroke::new(1.2, look.text));
        painter.text(
            pos2(x, y + 7.0),
            Align2::CENTER_TOP,
            text,
            FontId::proportional(11.0),
            look.dim,
        );
    }
}

/// The key to the colours at the top right.
pub(super) fn legend(painter: &Painter, look: &Look, legend: &Legend, rect: Rect) {
    let font = FontId::proportional(11.0);
    match legend {
        Legend::Gradient { key, low, high } => {
            let (w, h) = (140.0, 9.0);
            let left = rect.right() - w - MARGIN;
            let top = rect.top() + MARGIN + 14.0;
            painter.text(
                pos2(left, top - 3.0),
                Align2::LEFT_BOTTOM,
                key,
                font.clone(),
                look.text,
            );
            let steps = 28;
            for k in 0..steps {
                let t = f64::from(k) / f64::from(steps - 1);
                let x = left + w * k as f32 / steps as f32;
                painter.rect_filled(
                    Rect::from_min_max(pos2(x, top), pos2(x + w / steps as f32 + 0.5, top + h)),
                    0.0,
                    colour32(treestyle::gradient(t)),
                );
            }
            painter.text(
                pos2(left, top + h + 2.0),
                Align2::LEFT_TOP,
                treestyle::format_number(*low),
                font.clone(),
                look.dim,
            );
            painter.text(
                pos2(left + w, top + h + 2.0),
                Align2::RIGHT_TOP,
                treestyle::format_number(*high),
                font,
                look.dim,
            );
        }
        Legend::Categories { key, items } => {
            let left = rect.right() - 150.0 - MARGIN;
            let mut y = rect.top() + MARGIN + 14.0;
            painter.text(
                pos2(left, y - 3.0),
                Align2::LEFT_BOTTOM,
                key,
                font.clone(),
                look.text,
            );
            for (name, colour) in items.iter().take(12) {
                painter.rect_filled(
                    Rect::from_min_size(pos2(left, y + 1.0), vec2(10.0, 10.0)),
                    2.0,
                    colour32(*colour),
                );
                painter.text(
                    pos2(left + 15.0, y + 6.0),
                    Align2::LEFT_CENTER,
                    name,
                    font.clone(),
                    look.text,
                );
                y += 14.0;
            }
        }
    }
}
