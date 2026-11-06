//! Drawing a laid-out tree: strokes, tip names, support values and a scale bar, with zoom and pan and a
//! clade that can be picked by clicking on it. Shared by the Trees page and the preview on the Phylogeny page.

use crate::seq::treelayout::{Laid, Layout, Piece};
use crate::ui::design::Look;
use eframe::egui::{
    pos2, vec2, Align2, Color32, FontId, Pos2, Rect, Response, Sense, Shape, Stroke, Ui, Vec2,
};
use std::f32::consts::TAU;

const MARGIN: f32 = 14.0;
/// Fewest pixels between two rows before names are thinned out.
const MIN_LABEL_ROW: f32 = 6.0;
const HIT_RADIUS: f32 = 9.0;
const MAX_ZOOM: f32 = 60.0;
const MIN_ZOOM: f32 = 0.25;
/// How far apart, in radians, the points of a drawn arc are.
const ARC_STEP: f32 = 0.04;
/// The colours given to groups of tips, by the start of their names.
const GROUP_COLOURS: [[u8; 3]; 8] = [
    [0x2f, 0x6f, 0xed],
    [0xd6, 0x45, 0x45],
    [0x1e, 0x9e, 0x5a],
    [0xc7, 0x7c, 0x0e],
    [0x8a, 0x4f, 0xd0],
    [0x0e, 0x9a, 0xa7],
    [0xc2, 0x3c, 0x8f],
    [0x6b, 0x7a, 0x2f],
];
const MAX_GROUPS: usize = GROUP_COLOURS.len();

#[derive(Clone, Copy)]
pub struct Options {
    pub layout: Layout,
    pub labels: bool,
    pub support: bool,
    pub scale_bar: bool,
    pub groups: bool,
    pub font: f32,
    pub line: f32,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            layout: Layout::Phylogram,
            labels: true,
            support: true,
            scale_bar: true,
            groups: false,
            font: 12.5,
            line: 1.6,
        }
    }
}

/// Where the picture sits: zoomed about the middle of the area and moved by `pan` pixels.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Camera {
    pub zoom: f32,
    pub pan: Vec2,
}

impl Default for Camera {
    fn default() -> Self {
        Camera {
            zoom: 1.0,
            pan: Vec2::ZERO,
        }
    }
}

impl Camera {
    /// Zoom by `factor` keeping the point under `cursor` (measured from the picture's origin) where it is.
    pub fn zoom_about(&mut self, cursor: Vec2, factor: f32) {
        let zoom = (self.zoom * factor).clamp(MIN_ZOOM, MAX_ZOOM);
        let factor = zoom / self.zoom;
        let relative = cursor - self.pan;
        self.pan += relative * (1.0 - factor);
        self.zoom = zoom;
    }
}

/// The group of a tip: the start of its name up to the first `_`, space or `-`.
pub fn group_of(name: &str) -> &str {
    name.split(['_', ' ', '-']).next().unwrap_or(name)
}

/// A colour for each tip name when the names fall into a few groups, else none.
pub fn group_colours(laid: &Laid) -> Vec<Option<Color32>> {
    let mut groups: Vec<&str> = Vec::new();
    for n in laid.nodes.iter().filter(|n| n.is_tip()) {
        let g = group_of(&n.name);
        if !groups.contains(&g) {
            groups.push(g);
        }
    }
    let usable = (2..=MAX_GROUPS).contains(&groups.len());
    laid.nodes
        .iter()
        .map(|n| {
            let k = groups.iter().position(|&g| g == group_of(&n.name))?;
            (usable && n.is_tip()).then(|| {
                let c = GROUP_COLOURS[k];
                Color32::from_rgb(c[0], c[1], c[2])
            })
        })
        .collect()
}

/// How the layout's units become screen points.
struct Frame {
    origin: Pos2,
    scale: Vec2,
    /// What the layout's minimum maps to.
    min: (f32, f32),
    circular: bool,
}

impl Frame {
    fn at(&self, x: f32, y: f32) -> Pos2 {
        if self.circular {
            self.origin + vec2(x * self.scale.x, y * self.scale.y)
        } else {
            self.origin
                + vec2(
                    (x - self.min.0) * self.scale.x,
                    (y - self.min.1 + 0.5) * self.scale.y,
                )
        }
    }
}

fn frame(rect: Rect, laid: &Laid, o: &Options, camera: &Camera, label_room: f32) -> Frame {
    let (min, max) = laid.extent();
    if laid.kind == Layout::Circular {
        let room = (rect.width().min(rect.height()) / 2.0 - MARGIN - label_room).max(20.0);
        let side = room / (max.0 - min.0).max(f32::EPSILON) * 2.0 * camera.zoom;
        return Frame {
            origin: rect.center() + camera.pan,
            scale: vec2(side, side),
            min,
            circular: true,
        };
    }
    let width = (rect.width() - 2.0 * MARGIN - label_room).max(40.0);
    let rows = max.1 + 1.0;
    let height = (rect.height() - 2.0 * MARGIN).max(40.0);
    let _ = o;
    Frame {
        origin: rect.min + vec2(MARGIN, MARGIN) + camera.pan,
        scale: vec2(
            width / (max.0 - min.0).max(f32::EPSILON) * camera.zoom,
            height / rows * camera.zoom,
        ),
        min,
        circular: false,
    }
}

fn arc_points(centre: Pos2, radius: f32, from: f32, to: f32) -> Vec<Pos2> {
    let steps = (((to - from).abs() / ARC_STEP).ceil() as usize).max(2);
    (0..=steps)
        .map(|i| {
            let a = from + (to - from) * i as f32 / steps as f32;
            centre + vec2(radius * a.cos(), radius * a.sin())
        })
        .collect()
}

/// Draw the tree in the free space of `ui`. `selected` is a node whose clade is shown in the accent colour;
/// a click on a node or a tip name picks it (a click elsewhere puts it back). Zoom with the wheel, move by
/// dragging. `interactive` false draws only.
pub fn show(
    ui: &mut Ui,
    look: &Look,
    laid: &Laid,
    o: &Options,
    camera: &mut Camera,
    selected: &mut Option<usize>,
    interactive: bool,
) -> Response {
    let size = ui.available_size_before_wrap().max(vec2(80.0, 80.0));
    let sense = if interactive {
        Sense::click_and_drag()
    } else {
        Sense::hover()
    };
    let (response, painter) = ui.allocate_painter(size, sense);
    let rect = response.rect;
    let font = FontId::proportional(o.font);
    let longest = laid
        .nodes
        .iter()
        .filter(|n| n.is_tip())
        .map(|n| n.name.chars().count())
        .max()
        .unwrap_or(0) as f32;
    let label_room = if o.labels {
        longest * o.font * 0.56 + 14.0
    } else {
        0.0
    };

    if interactive {
        if response.dragged() {
            camera.pan += response.drag_delta();
        }
        if response.hovered() {
            let (scroll, pinch) = ui.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
            let factor = pinch * (scroll * 0.004).exp();
            if let (Some(pointer), true) = (response.hover_pos(), (factor - 1.0).abs() > 1e-4) {
                let origin = frame(rect, laid, o, &Camera::default(), label_room).origin;
                camera.zoom_about(pointer - origin, factor);
            }
        }
    }
    let f = frame(rect, laid, o, camera, label_room);
    let tips = laid.nodes.first().map_or(1, |n| n.tips).max(1) as f32;
    let row_px = if f.circular {
        f.scale.x * laid.extent().1 .0 * TAU / tips
    } else {
        f.scale.y
    };
    let colours = if o.groups {
        group_colours(laid)
    } else {
        vec![None; laid.nodes.len()]
    };
    let base = look.text.gamma_multiply(0.85);
    let in_clade = |i: usize| selected.is_some_and(|s| laid.is_within(i, s));

    for (i, piece) in laid.pieces() {
        let picked = in_clade(i);
        let stroke = Stroke::new(
            if picked { o.line + 1.0 } else { o.line },
            if picked { look.accent } else { base },
        );
        match piece {
            Piece::Line { from, to } => {
                painter.line_segment([f.at(from.0, from.1), f.at(to.0, to.1)], stroke);
            }
            Piece::Arc { radius, from, to } => {
                let points = arc_points(f.origin, radius * f.scale.x, from, to);
                painter.add(Shape::line(points, stroke));
            }
        }
    }

    // tip names, thinned out when the rows are closer than a line of text
    let every = if o.labels && row_px < o.font * 0.95 {
        (o.font * 0.95 / row_px.max(0.1)).ceil() as usize
    } else {
        1
    };
    let show_names = o.labels && row_px >= MIN_LABEL_ROW;
    let mut tip_number = 0usize;
    for (i, node) in laid.nodes.iter().enumerate() {
        let at = f.at(node.x, node.y);
        if node.is_tip() {
            tip_number += 1;
            if !show_names
                || !(tip_number - 1).is_multiple_of(every)
                || !rect.expand(40.0).contains(at)
            {
                continue;
            }
            let colour = if in_clade(i) {
                look.accent
            } else {
                colours[i].unwrap_or(look.text)
            };
            let (anchor, shift) = if f.circular && node.x < 0.0 {
                (Align2::RIGHT_CENTER, -6.0)
            } else {
                (Align2::LEFT_CENTER, 6.0)
            };
            painter.text(
                at + vec2(shift, 0.0),
                anchor,
                &node.name,
                font.clone(),
                colour,
            );
        } else if o.support && row_px >= 11.0 {
            if let Some(s) = node.support {
                // on a short branch the number crosses the line of the parent: give it a background
                let galley = painter.layout_no_wrap(
                    format!("{s:.0}"),
                    FontId::proportional((o.font - 3.0).max(8.0)),
                    look.dim,
                );
                let box_at = Align2::RIGHT_BOTTOM.anchor_size(at + vec2(-3.0, -2.0), galley.size());
                painter.rect_filled(box_at.expand(1.0), 2.0, look.surface);
                painter.galley(box_at.min, galley, look.dim);
            }
        }
    }

    if o.scale_bar && laid.kind != Layout::Circular {
        let (_, max) = laid.extent();
        let length = Laid::nice_length(max.0 * 0.25);
        if length > 0.0 {
            let y = rect.bottom() - 16.0;
            let (x1, x2) = (
                rect.left() + MARGIN,
                rect.left() + MARGIN + length * f.scale.x,
            );
            let stroke = Stroke::new(1.6, look.text);
            painter.line_segment([pos2(x1, y), pos2(x2, y)], stroke);
            painter.text(
                pos2((x1 + x2) / 2.0, y - 4.0),
                Align2::CENTER_BOTTOM,
                format!("{length}"),
                FontId::proportional(11.0),
                look.text,
            );
        }
    }

    // what is under the pointer: a node, or a tip and its name
    let hit = response.hover_pos().and_then(|p| {
        laid.nodes
            .iter()
            .enumerate()
            .filter_map(|(i, n)| {
                let at = f.at(n.x, n.y);
                let name_room = if n.is_tip() && o.labels && show_names {
                    n.name.chars().count() as f32 * o.font * 0.56
                } else {
                    0.0
                };
                let inside = (p.x >= at.x - HIT_RADIUS && p.x <= at.x + HIT_RADIUS + name_room)
                    && (p.y - at.y).abs() <= HIT_RADIUS.max(row_px / 2.0).min(14.0);
                inside.then_some((i, (p - at).length()))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(i, _)| i)
    });
    if interactive {
        if let Some(i) = hit {
            let n = &laid.nodes[i];
            let at = f.at(n.x, n.y);
            painter.circle_filled(at, 4.0, look.accent);
            let mut text = if n.is_tip() {
                n.name.clone()
            } else {
                format!("{} tips", n.tips)
            };
            if let Some(s) = n.support {
                text.push_str(&format!("   support {s:.0}"));
            }
            text.push_str(&format!("   {:.4} from the root", n.radius));
            response.clone().on_hover_text(text);
            if response.clicked() {
                *selected = if *selected == Some(i) { None } else { Some(i) };
            }
        } else if response.clicked() {
            *selected = None;
        }
        if let Some(s) = *selected {
            painter.circle_stroke(
                f.at(laid.nodes[s].x, laid.nodes[s].y),
                5.5,
                Stroke::new(2.0, look.accent),
            );
        }
    }
    response
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::newick::NwTree;
    use crate::seq::treelayout::lay_out;

    #[test]
    fn zooming_keeps_the_point_under_the_pointer_where_it_is() {
        let mut camera = Camera::default();
        let cursor = vec2(100.0, 60.0);
        // the layout point that sits under the cursor before ...
        let before = (cursor - camera.pan) / camera.zoom;
        camera.zoom_about(cursor, 2.0);
        assert!((camera.zoom - 2.0).abs() < 1e-6);
        // ... is still under it afterwards
        let after = (cursor - camera.pan) / camera.zoom;
        assert!((before - after).length() < 1e-4, "{before:?} {after:?}");
        // and zoom stays between its limits
        camera.zoom_about(cursor, 1e6);
        assert_eq!(camera.zoom, MAX_ZOOM);
        camera.zoom_about(cursor, 1e-9);
        assert_eq!(camera.zoom, MIN_ZOOM);
    }

    #[test]
    fn groups_are_the_start_of_the_names_and_only_a_few_are_coloured() {
        assert_eq!(group_of("East_12"), "East");
        assert_eq!(group_of("West 3"), "West");
        assert_eq!(group_of("plain"), "plain");
        let two = lay_out(
            &NwTree::parse_all("((E_1:1,E_2:1):1,(W_1:1,W_2:1):1);")
                .unwrap()
                .remove(0)
                .to_node(),
            Layout::Phylogram,
        );
        let colours = group_colours(&two);
        let tips: Vec<Option<Color32>> = two
            .nodes
            .iter()
            .zip(&colours)
            .filter(|(n, _)| n.is_tip())
            .map(|(_, c)| *c)
            .collect();
        assert!(tips.iter().all(Option::is_some));
        assert_eq!(tips[0], tips[1]);
        assert_ne!(tips[1], tips[2]);
        // one group, or too many, get no colours
        let one = lay_out(
            &NwTree::parse_all("(E_1:1,E_2:1);")
                .unwrap()
                .remove(0)
                .to_node(),
            Layout::Phylogram,
        );
        assert!(group_colours(&one).iter().all(Option::is_none));
    }
}
