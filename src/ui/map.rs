//! Graphic maps: circular (angle 0 at the top, clockwise) and linear, with features, enzyme
//! sites and a position scale.

use super::editor::head_of;
use super::theme::{self, Palette};
use crate::seq::doc::Doc;
use eframe::egui::{
    self, epaint::Mesh, pos2, vec2, Align2, Color32, FontId, Pos2, Rect, Sense, Shape, Stroke, Ui,
    Vec2,
};
use std::f32::consts::{FRAC_PI_2, TAU};

pub struct MapSite {
    pub label: String,
    /// Base after which the top strand is cut (1-based position shown in the label).
    pub pos: usize,
}

pub struct MapData<'a> {
    pub doc: &'a Doc,
    pub sites: Vec<MapSite>,
    pub selection: Option<(usize, usize)>,
    pub pal: &'a Palette,
}

/// "1,000" style numbers.
pub fn commify(n: usize) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// A round spacing for scale marks so that about `target` marks fit.
pub fn nice_step(n: usize, target: usize) -> usize {
    let raw = (n.max(1) as f64 / target.max(1) as f64).max(1.0);
    let mag = 10f64.powf(raw.log10().floor());
    let f = raw / mag;
    let step = if f <= 1.0 {
        1.0
    } else if f <= 2.0 {
        2.0
    } else if f <= 5.0 {
        5.0
    } else {
        10.0
    };
    (step * mag) as usize
}

fn label_font() -> FontId {
    FontId::proportional(11.0)
}

/// Height of one row of enzyme labels on the linear map.
const ENZYME_ROW_H: f32 = 14.0;

// ---------------------------------------------------------------------------
// Circular
// ---------------------------------------------------------------------------

struct Circle {
    c: Pos2,
    n: f32,
}

impl Circle {
    fn angle(&self, pos: f32) -> f32 {
        pos / self.n * TAU - FRAC_PI_2
    }

    fn point(&self, r: f32, pos: f32) -> Pos2 {
        let a = self.angle(pos);
        self.c + vec2(a.cos(), a.sin()) * r
    }
}

/// Filled arc band with an optional arrow head at one end (angles in radians, a0 < a1).
#[allow(clippy::too_many_arguments)]
fn arc_arrow(
    p: &egui::Painter,
    c: Pos2,
    r0: f32,
    r1: f32,
    a0: f32,
    a1: f32,
    head: Option<bool>,
    fill: Color32,
    edge: Stroke,
) {
    let rm = (r0 + r1) / 2.0;
    let head_len = ((r1 - r0) / rm).min((a1 - a0) * 0.45);
    let (b0, b1) = match head {
        Some(true) => (a0, a1 - head_len),
        Some(false) => (a0 + head_len, a1),
        None => (a0, a1),
    };
    let steps = (((b1 - b0) / 0.06).ceil() as usize).max(2);
    let at = |r: f32, a: f32| c + vec2(a.cos(), a.sin()) * r;
    let mut mesh = Mesh::default();
    for i in 0..=steps {
        let a = b0 + (b1 - b0) * i as f32 / steps as f32;
        mesh.colored_vertex(at(r1, a), fill);
        mesh.colored_vertex(at(r0, a), fill);
    }
    for i in 0..steps as u32 {
        let k = i * 2;
        mesh.add_triangle(k, k + 1, k + 2);
        mesh.add_triangle(k + 1, k + 3, k + 2);
    }
    let mut outline: Vec<Pos2> = (0..=steps)
        .map(|i| at(r1, b0 + (b1 - b0) * i as f32 / steps as f32))
        .collect();
    match head {
        Some(true) => {
            let base = mesh.vertices.len() as u32;
            mesh.colored_vertex(at(r1, b1), fill);
            mesh.colored_vertex(at(r0, b1), fill);
            mesh.colored_vertex(at(rm, a1), fill);
            mesh.add_triangle(base, base + 1, base + 2);
            outline.push(at(rm, a1));
        }
        Some(false) => {
            let base = mesh.vertices.len() as u32;
            mesh.colored_vertex(at(r1, b0), fill);
            mesh.colored_vertex(at(r0, b0), fill);
            mesh.colored_vertex(at(rm, a0), fill);
            mesh.add_triangle(base, base + 1, base + 2);
            outline.insert(0, at(rm, a0));
        }
        None => {}
    }
    outline.extend(
        (0..=steps)
            .rev()
            .map(|i| at(r0, b0 + (b1 - b0) * i as f32 / steps as f32)),
    );
    p.add(Shape::mesh(mesh));
    p.add(Shape::closed_line(outline, edge));
}

pub fn circular(ui: &mut Ui, data: &MapData, size: Vec2) -> Option<usize> {
    let (rect, resp) = ui.allocate_exact_size(size, Sense::click());
    let p = ui.painter_at(rect);
    let (doc, pal) = (data.doc, data.pal);
    let n = doc.len().max(1);
    let lane_count = feature_lanes(doc).1;
    let r = ((rect.height() * 0.5 - 28.0)
        .min(rect.width() * 0.5 - 105.0 - lane_count as f32 * 10.0))
    .max(50.0);
    let circ = Circle {
        c: rect.center(),
        n: n as f32,
    };

    p.circle_stroke(circ.c, r, Stroke::new(1.5, pal.fg));
    scale_ticks(&p, &circ, r, n, pal);
    if let Some((a, b)) = data.selection {
        arc_polyline(
            &p,
            &circ,
            r - 7.0,
            a as f32,
            b as f32,
            Stroke::new(3.0, pal.select_bg),
        );
    }
    let title = p.layout_no_wrap(doc.name.clone(), FontId::proportional(14.0), pal.fg);
    let size_text = p.layout_no_wrap(
        format!("{} bp", commify(doc.len())),
        FontId::proportional(12.0),
        pal.dim,
    );
    // regions the feature names must keep clear of: the centre text and the scale numbers
    let keep_clear = vec![
        Rect::from_center_size(circ.c - vec2(0.0, 10.0), title.size()).expand(3.0),
        Rect::from_center_size(circ.c + vec2(0.0, 10.0), size_text.size()).expand(3.0),
    ];
    feature_arcs(&p, &circ, r, doc, pal, keep_clear);
    enzyme_labels(&p, &circ, r + 8.0 + lane_count as f32 * 13.0, data, rect);

    p.galley(circ.c - vec2(0.0, 10.0) - title.size() / 2.0, title, pal.fg);
    p.galley(
        circ.c + vec2(0.0, 10.0) - size_text.size() / 2.0,
        size_text,
        pal.dim,
    );

    resp.interact_pointer_pos()
        .filter(|_| resp.clicked())
        .map(|pos| {
            let v = pos - circ.c;
            let turn = (v.y.atan2(v.x) + FRAC_PI_2).rem_euclid(TAU) / TAU;
            ((turn * n as f32) as usize).min(n)
        })
}

fn arc_polyline(p: &egui::Painter, c: &Circle, r: f32, a: f32, b: f32, stroke: Stroke) {
    let (a, b) = if b < a { (a, b + c.n) } else { (a, b) };
    let steps = (((b - a) / c.n * TAU / 0.05).ceil() as usize).max(2);
    let pts: Vec<Pos2> = (0..=steps)
        .map(|i| c.point(r, a + (b - a) * i as f32 / steps as f32))
        .collect();
    p.add(Shape::line(pts, stroke));
}

fn scale_ticks(p: &egui::Painter, c: &Circle, r: f32, n: usize, pal: &Palette) {
    let major = nice_step(n, 10);
    let minor = (major / 4).max(1);
    let mut pos = 0;
    while pos < n {
        let is_major = pos % major == 0;
        let len = if is_major { 8.0 } else { 4.0 };
        p.line_segment(
            [c.point(r, pos as f32), c.point(r - len, pos as f32)],
            Stroke::new(if is_major { 1.5 } else { 1.0 }, pal.fg),
        );
        if is_major && pos > 0 {
            p.text(
                c.point(r - 20.0, pos as f32),
                Align2::CENTER_CENTER,
                commify(pos),
                FontId::proportional(10.0),
                pal.dim,
            );
        }
        pos += minor;
    }
}

/// Lane of every feature (stacked outward where they overlap) and the lane count.
fn feature_lanes(doc: &Doc) -> (Vec<usize>, usize) {
    let layout_lanes = super::layout::Layout::new(doc, 1_000_000, 1.0, 1.0, false, 0).lanes;
    let count = layout_lanes.iter().max().map_or(0, |m| m + 1);
    (layout_lanes, count)
}

fn feature_arcs(
    p: &egui::Painter,
    c: &Circle,
    r: f32,
    doc: &Doc,
    pal: &Palette,
    keep_clear: Vec<Rect>,
) {
    let (lanes, _) = feature_lanes(doc);
    let n = doc.len();
    let mut labelled: Vec<Rect> = keep_clear;
    for (fi, f) in doc.features.iter().enumerate() {
        let r0 = r + 3.0 + lanes[fi] as f32 * 13.0;
        let r1 = r0 + 10.0;
        let fill = theme::rgb(f.color());
        let edge = Stroke::new(1.0, theme::mix(Color32::BLACK, fill, 35.0));
        let (head_span, head_right) = head_of(f, n);
        for (si, &(s, e)) in f.spans.iter().enumerate() {
            let (a0, a1) = (c.angle(s as f32), c.angle(e as f32));
            let head = (si == head_span).then_some(head_right);
            arc_arrow(p, c.c, r0, r1, a0, a1.max(a0 + 0.01), head, fill, edge);
        }
        // name inside the circle at the middle of the largest span: one line first, then two,
        // at the first radius where it stays inside the ring and clear of other text
        let (s, e) = f
            .spans
            .iter()
            .copied()
            .max_by_key(|&(s, e)| e - s)
            .unwrap_or((0, 0));
        let one_line = p.layout_no_wrap(f.name.clone(), label_font(), pal.fg);
        let two_lines = p.layout(f.name.clone(), label_font(), pal.fg, 64.0);
        'place: for galley in [one_line, two_lines] {
            for ring in [30.0, 46.0, 62.0, 78.0, 94.0] {
                let rect =
                    Rect::from_center_size(c.point(r - ring, (s + e) as f32 / 2.0), galley.size())
                        .expand(1.0);
                let corners = [
                    rect.left_top(),
                    rect.right_top(),
                    rect.left_bottom(),
                    rect.right_bottom(),
                ];
                let inside = corners
                    .iter()
                    .all(|corner| (*corner - c.c).length() <= r - 10.0);
                if inside && !labelled.iter().any(|q| q.intersects(rect)) {
                    p.galley(rect.min + vec2(1.0, 1.0), galley.clone(), pal.fg);
                    labelled.push(rect);
                    break 'place;
                }
            }
        }
    }
}

fn enzyme_labels(p: &egui::Painter, c: &Circle, r_label: f32, data: &MapData, area: Rect) {
    struct Label<'a> {
        site: &'a MapSite,
        tick: Pos2,
        y: f32,
    }
    let mut sides: [Vec<Label>; 2] = [Vec::new(), Vec::new()];
    for s in &data.sites {
        let tick = c.point(r_label, s.pos as f32);
        let a = c.angle(s.pos as f32);
        sides[usize::from(a.cos() < 0.0)].push(Label {
            site: s,
            tick,
            y: tick.y,
        });
    }
    let lh = 13.0;
    for (side, labels) in sides.iter_mut().enumerate() {
        labels.sort_by(|a, b| a.y.total_cmp(&b.y));
        let mut floor = area.top() + 6.0;
        for l in labels.iter_mut() {
            l.y = l.y.max(floor);
            floor = l.y + lh;
        }
        // if the column ran past the bottom, push it back up from the end
        let mut ceiling = area.bottom() - 8.0;
        for l in labels.iter_mut().rev() {
            l.y = l.y.min(ceiling);
            ceiling = l.y - lh;
        }
        let sign = if side == 0 { 1.0 } else { -1.0 };
        let column = c.c.x + sign * (r_label + 18.0);
        for l in labels.iter() {
            let from = c.point(r_label - 8.0, l.site.pos as f32);
            let col = data.pal.dim;
            p.line_segment([from, l.tick], Stroke::new(1.0, data.pal.fg));
            p.line_segment([l.tick, pos2(column, l.y)], Stroke::new(0.8, col));
            let anchor = if side == 0 {
                Align2::LEFT_CENTER
            } else {
                Align2::RIGHT_CENTER
            };
            p.text(
                pos2(column + sign * 3.0, l.y),
                anchor,
                format!("{} ({})", l.site.label, l.site.pos),
                label_font(),
                data.pal.fg,
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Linear
// ---------------------------------------------------------------------------

pub fn linear(ui: &mut Ui, data: &MapData, size: Vec2) -> Option<usize> {
    let (rect, resp) = ui.allocate_exact_size(size, Sense::click());
    let p = ui.painter_at(rect);
    let (doc, pal) = (data.doc, data.pal);
    let n = doc.len().max(1);
    let (left, right) = (rect.left() + 24.0, rect.right() - 24.0);
    let x_of = |pos: usize| left + pos as f32 / n as f32 * (right - left);

    // Enzyme labels start at their tick and stack upward in rows, as many rows as the panel holds.
    let mut sites: Vec<&MapSite> = data.sites.iter().collect();
    sites.sort_by_key(|s| s.pos);
    let max_rows = (((rect.height() - 160.0) / ENZYME_ROW_H).floor() as usize).clamp(3, 40);
    let mut row_end: Vec<f32> = Vec::new();
    let mut placed = Vec::new();
    let mut hidden = 0;
    for s in &sites {
        let galley = p.layout_no_wrap(format!("{} ({})", s.label, s.pos), label_font(), pal.fg);
        let tick = x_of(s.pos);
        let start = tick + 2.0;
        match row_end.iter().position(|&e| e + 5.0 <= start) {
            Some(row) => {
                row_end[row] = start + galley.size().x;
                placed.push((row, tick, start, galley));
            }
            None if row_end.len() < max_rows => {
                row_end.push(start + galley.size().x);
                placed.push((row_end.len() - 1, tick, start, galley));
            }
            None => hidden += 1,
        }
    }
    let axis_y = rect.top() + 18.0 + row_end.len() as f32 * ENZYME_ROW_H + 10.0;
    let row_y = |row: usize| axis_y - 12.0 - row as f32 * ENZYME_ROW_H;
    for (row, tick, _, _) in &placed {
        p.line_segment(
            [pos2(*tick, row_y(*row)), pos2(*tick, axis_y)],
            Stroke::new(0.6, pal.dim.gamma_multiply(0.7)),
        );
    }
    for (row, _, start, galley) in &placed {
        let text_rect = Rect::from_min_size(
            pos2(*start, row_y(*row) - galley.size().y / 2.0),
            galley.size(),
        );
        p.rect_filled(text_rect.expand2(vec2(1.0, 0.0)), 0.0, pal.window_bg);
        p.galley(text_rect.min, galley.clone(), pal.fg);
    }
    if hidden > 0 {
        p.text(
            pos2(right, rect.top() + 4.0),
            Align2::RIGHT_TOP,
            format!("+{hidden} more sites (select fewer enzymes)"),
            FontId::proportional(10.0),
            pal.dim,
        );
    }

    p.line_segment(
        [pos2(left, axis_y), pos2(right, axis_y)],
        Stroke::new(1.5, pal.fg),
    );
    let major = nice_step(n, 8);
    let mut pos = 0;
    while pos <= n {
        let x = x_of(pos);
        p.line_segment(
            [pos2(x, axis_y), pos2(x, axis_y + 6.0)],
            Stroke::new(1.5, pal.fg),
        );
        p.text(
            pos2(x, axis_y + 8.0),
            Align2::CENTER_TOP,
            commify(pos),
            FontId::proportional(10.0),
            pal.dim,
        );
        pos += major;
    }
    if let Some((a, b)) = data.selection {
        p.rect_filled(
            Rect::from_min_max(pos2(x_of(a), axis_y - 3.0), pos2(x_of(b), axis_y + 3.0)),
            0.0,
            pal.select_bg,
        );
    }
    feature_bars(&p, doc, axis_y + 26.0, &x_of, pal);
    p.text(
        pos2(rect.center().x, rect.bottom() - 4.0),
        Align2::CENTER_BOTTOM,
        format!("{}  {} bp", doc.name, commify(doc.len())),
        FontId::proportional(12.0),
        pal.dim,
    );

    resp.interact_pointer_pos()
        .filter(|_| resp.clicked())
        .map(|pos| (((pos.x - left) / (right - left)).clamp(0.0, 1.0) * n as f32) as usize)
}

fn feature_bars(
    p: &egui::Painter,
    doc: &Doc,
    top: f32,
    x_of: &dyn Fn(usize) -> f32,
    pal: &Palette,
) {
    let (lanes, _) = feature_lanes(doc);
    let n = doc.len();
    for (fi, f) in doc.features.iter().enumerate() {
        let y = top + lanes[fi] as f32 * 17.0;
        let fill = theme::rgb(f.color());
        let (head_span, head_right) = head_of(f, n);
        for (si, &(s, e)) in f.spans.iter().enumerate() {
            let (xa, xb) = (x_of(s), x_of(e).max(x_of(s) + 2.0));
            let (t, b) = (y, y + 14.0);
            let hw = 7.0f32.min(xb - xa);
            let has_head = si == head_span;
            let pts = match (has_head, head_right) {
                (false, _) => vec![pos2(xa, t), pos2(xb, t), pos2(xb, b), pos2(xa, b)],
                (true, true) => vec![
                    pos2(xa, t),
                    pos2(xb - hw, t),
                    pos2(xb, (t + b) / 2.0),
                    pos2(xb - hw, b),
                    pos2(xa, b),
                ],
                (true, false) => vec![
                    pos2(xa + hw, t),
                    pos2(xb, t),
                    pos2(xb, b),
                    pos2(xa + hw, b),
                    pos2(xa, (t + b) / 2.0),
                ],
            };
            p.add(Shape::convex_polygon(
                pts,
                fill,
                Stroke::new(1.0, theme::mix(Color32::BLACK, fill, 35.0)),
            ));
            if si == 0 {
                let clip = Rect::from_min_max(pos2(xa + 1.0, t), pos2(xb - 1.0, b));
                p.with_clip_rect(clip).text(
                    pos2(xa + 3.0, (t + b) / 2.0),
                    Align2::LEFT_CENTER,
                    &f.name,
                    FontId::proportional(10.0),
                    theme::text_on(fill),
                );
            }
        }
    }
    let _ = pal;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_get_thousands_separators() {
        assert_eq!(commify(0), "0");
        assert_eq!(commify(999), "999");
        assert_eq!(commify(1000), "1,000");
        assert_eq!(commify(2686), "2,686");
        assert_eq!(commify(1234567), "1,234,567");
    }

    #[test]
    fn scale_steps_are_round_numbers() {
        assert_eq!(nice_step(2686, 10), 500);
        assert_eq!(nice_step(10_000, 10), 1000);
        assert_eq!(nice_step(120, 10), 20);
        assert_eq!(nice_step(5, 10), 1);
        assert_eq!(nice_step(0, 10), 1);
    }
}
