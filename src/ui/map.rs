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
    /// The feature the inspector is showing: drawn with a heavy outline so it can be found on the map.
    pub selected_feature: Option<usize>,
    pub pal: &'a Palette,
}

/// What a click on the map landed on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MapHit {
    /// Empty map: this base (0-based, 0..=len).
    Base(usize),
    /// A feature's arrow.
    Feature(usize),
}

/// Outline of a feature: a thin dark edge, or a heavy one in the text colour when it is the selected one.
fn feature_edge(fill: Color32, selected: bool, pal: &Palette) -> Stroke {
    if selected {
        Stroke::new(2.5, pal.fg)
    } else {
        Stroke::new(1.0, theme::mix(Color32::BLACK, fill, 35.0))
    }
}

/// Distance from the circle's line to the near and far edge of lane `lane`'s ring of features.
fn ring_band(lane: usize) -> (f32, f32) {
    let near = 3.0 + lane as f32 * 13.0;
    (near, near + 10.0)
}

/// The least angle, in radians, that an arrow of a feature is drawn over on the circular map: a feature of one base is still seen.
const MIN_ARC: f32 = 0.01;

/// The feature under a click on the circular map: `out` is how far the click is beyond the circle's line, `base` the
/// base at its angle. A feature too short to be seen is hit over the width that it is drawn with.
fn ring_hit(doc: &Doc, lanes: &[usize], out: f32, base: usize) -> Option<usize> {
    let drawn = (MIN_ARC / TAU * doc.len() as f32).ceil() as usize;
    doc.features.iter().enumerate().find_map(|(fi, f)| {
        let (near, far) = ring_band(lanes[fi]);
        let over = f
            .spans
            .iter()
            .any(|&(s, e)| base >= s && base < e.max(s + drawn));
        (out >= near - 1.0 && out <= far + 1.0 && over).then_some(fi)
    })
}

/// Height of a feature's bar on the linear map, and the distance from one lane's bars to the next.
const BAR_HEIGHT: f32 = 14.0;
const LANE_PITCH: f32 = 17.0;

/// The feature under a click on the linear map: `x_of` gives the x of a base, `top` is where lane 0's bars start.
fn bar_hit(
    doc: &Doc,
    lanes: &[usize],
    x_of: &dyn Fn(usize) -> f32,
    top: f32,
    at: Pos2,
) -> Option<usize> {
    doc.features.iter().enumerate().find_map(|(fi, f)| {
        let y = top + lanes[fi] as f32 * LANE_PITCH;
        if at.y < y || at.y > y + BAR_HEIGHT {
            return None;
        }
        let on = f.spans.iter().any(|&(s, e)| {
            let (xa, xb) = (x_of(s), x_of(e).max(x_of(s) + 2.0));
            at.x >= xa - 1.0 && at.x <= xb + 1.0
        });
        on.then_some(fi)
    })
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

pub fn circular(ui: &mut Ui, data: &MapData, size: Vec2) -> Option<MapHit> {
    let (rect, resp) = ui.allocate_exact_size(size, Sense::click());
    let p = ui.painter_at(rect);
    let (doc, pal) = (data.doc, data.pal);
    let n = doc.len().max(1);
    let (lanes, lane_count) = feature_lanes(doc);
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
    feature_arcs(&p, &circ, r, data, &lanes, keep_clear);
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
            let base = ((turn * n as f32) as usize).min(n);
            match ring_hit(doc, &lanes, v.length() - r, base) {
                Some(fi) => MapHit::Feature(fi),
                None => MapHit::Base(base),
            }
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
    data: &MapData,
    lanes: &[usize],
    keep_clear: Vec<Rect>,
) {
    let (doc, pal) = (data.doc, data.pal);
    let n = doc.len();
    let mut labelled: Vec<Rect> = keep_clear;
    for (fi, f) in doc.features.iter().enumerate() {
        let (near, far) = ring_band(lanes[fi]);
        let (r0, r1) = (r + near, r + far);
        let fill = theme::rgb(f.color());
        let edge = feature_edge(fill, data.selected_feature == Some(fi), pal);
        let (head_span, head_right) = head_of(f, n);
        for (si, &(s, e)) in f.spans.iter().enumerate() {
            let (a0, a1) = (c.angle(s as f32), c.angle(e as f32));
            let head = (si == head_span).then_some(head_right);
            arc_arrow(p, c.c, r0, r1, a0, a1.max(a0 + MIN_ARC), head, fill, edge);
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

pub fn linear(ui: &mut Ui, data: &MapData, size: Vec2) -> Option<MapHit> {
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
    let (lanes, _) = feature_lanes(doc);
    let bars_top = axis_y + 26.0;
    feature_bars(&p, data, &lanes, bars_top, &x_of);
    p.text(
        pos2(rect.center().x, rect.bottom() - 4.0),
        Align2::CENTER_BOTTOM,
        format!("{}  {} bp", doc.name, commify(doc.len())),
        FontId::proportional(12.0),
        pal.dim,
    );

    resp.interact_pointer_pos()
        .filter(|_| resp.clicked())
        .map(|pos| match bar_hit(doc, &lanes, &x_of, bars_top, pos) {
            Some(fi) => MapHit::Feature(fi),
            None => MapHit::Base(
                (((pos.x - left) / (right - left)).clamp(0.0, 1.0) * n as f32) as usize,
            ),
        })
}

fn feature_bars(
    p: &egui::Painter,
    data: &MapData,
    lanes: &[usize],
    top: f32,
    x_of: &dyn Fn(usize) -> f32,
) {
    let (doc, pal) = (data.doc, data.pal);
    let n = doc.len();
    for (fi, f) in doc.features.iter().enumerate() {
        let y = top + lanes[fi] as f32 * LANE_PITCH;
        let fill = theme::rgb(f.color());
        let (head_span, head_right) = head_of(f, n);
        for (si, &(s, e)) in f.spans.iter().enumerate() {
            let (xa, xb) = (x_of(s), x_of(e).max(x_of(s) + 2.0));
            let (t, b) = (y, y + BAR_HEIGHT);
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
                feature_edge(fill, data.selected_feature == Some(fi), pal),
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

    fn two_features() -> (Doc, Vec<usize>) {
        use crate::seq::features::Feature;
        let mut d = Doc::new("t", vec![b'A'; 100]);
        d.features.push(Feature::new("a", "CDS", 10, 40, false));
        d.features
            .push(Feature::new("b", "promoter", 30, 60, false));
        let (lanes, _) = feature_lanes(&d);
        (d, lanes)
    }

    #[test]
    fn a_click_on_the_ring_picks_the_feature_in_that_lane_at_that_base() {
        let (d, lanes) = two_features();
        assert_ne!(lanes[0], lanes[1], "overlapping features are stacked");
        let mid = |lane: usize| (ring_band(lane).0 + ring_band(lane).1) / 2.0;
        assert_eq!(ring_hit(&d, &lanes, mid(lanes[0]), 20), Some(0));
        assert_eq!(ring_hit(&d, &lanes, mid(lanes[1]), 50), Some(1));
        // in the overlap, the lane decides
        assert_eq!(ring_hit(&d, &lanes, mid(lanes[0]), 35), Some(0));
        assert_eq!(ring_hit(&d, &lanes, mid(lanes[1]), 35), Some(1));
    }

    #[test]
    fn a_feature_of_one_base_is_hit_over_the_width_it_is_drawn_with() {
        use crate::seq::features::Feature;
        let mut d = Doc::new("t", vec![b'A'; 10_000]);
        d.features
            .push(Feature::new("tiny", "misc", 5_000, 5_001, false));
        let (lanes, _) = feature_lanes(&d);
        let mid = (ring_band(lanes[0]).0 + ring_band(lanes[0]).1) / 2.0;
        // drawn about 16 bases wide: its own base and the next few are the arrow, a hundred bases on is not
        assert_eq!(ring_hit(&d, &lanes, mid, 5_000), Some(0));
        assert_eq!(ring_hit(&d, &lanes, mid, 5_010), Some(0));
        assert_eq!(ring_hit(&d, &lanes, mid, 5_100), None);
        assert_eq!(ring_hit(&d, &lanes, mid, 4_999), None);
    }

    #[test]
    fn a_click_off_the_features_picks_none() {
        let (d, lanes) = two_features();
        let mid = (ring_band(lanes[0]).0 + ring_band(lanes[0]).1) / 2.0;
        assert_eq!(
            ring_hit(&d, &lanes, mid, 70),
            None,
            "no feature at that base"
        );
        assert_eq!(ring_hit(&d, &lanes, -5.0, 20), None, "inside the circle");
        assert_eq!(ring_hit(&d, &lanes, 80.0, 20), None, "far outside it");
    }

    #[test]
    fn a_click_on_a_bar_picks_the_feature_in_that_lane_at_that_x() {
        let (d, lanes) = two_features();
        let x_of = |pos: usize| 10.0 + pos as f32 * 5.0;
        let top = 100.0;
        let y = |lane: usize| top + lane as f32 * LANE_PITCH + BAR_HEIGHT / 2.0;
        assert_eq!(
            bar_hit(&d, &lanes, &x_of, top, pos2(x_of(20), y(lanes[0]))),
            Some(0)
        );
        assert_eq!(
            bar_hit(&d, &lanes, &x_of, top, pos2(x_of(50), y(lanes[1]))),
            Some(1)
        );
        // beside the bar, and between the lanes
        assert_eq!(
            bar_hit(&d, &lanes, &x_of, top, pos2(x_of(70), y(lanes[0]))),
            None
        );
        let gap = top + BAR_HEIGHT + 1.5;
        assert_eq!(bar_hit(&d, &lanes, &x_of, top, pos2(x_of(20), gap)), None);
    }

    #[test]
    fn a_feature_thinner_than_two_pixels_can_still_be_hit() {
        use crate::seq::features::Feature;
        let mut d = Doc::new("t", vec![b'A'; 10_000]);
        d.features
            .push(Feature::new("tiny", "misc", 5_000, 5_001, false));
        let (lanes, _) = feature_lanes(&d);
        let x_of = |pos: usize| 10.0 + pos as f32 * 0.05;
        let at = pos2(x_of(5_000) + 1.0, 100.0 + BAR_HEIGHT / 2.0);
        assert_eq!(bar_hit(&d, &lanes, &x_of, 100.0, at), Some(0));
    }
}
