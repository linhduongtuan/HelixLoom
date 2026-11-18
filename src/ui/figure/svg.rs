//! What a view draws, as an SVG picture (`docs/ROADMAP.md` 1.7): every shape the view gave egui becomes the SVG element that draws
//! the same thing in the same place, a point being one SVG unit, so a figure has the lines, the colours and the text that are on
//! the screen. Text stays text (a run for each stretch of one font and colour, with the width it has on screen), so it can be
//! edited in a drawing program; the fonts are named, not carried, so a computer without them draws the text in its nearest font,
//! at the same width. Not drawn: pictures (textures), soft shadows, what a view paints through a callback, and a colour that is
//! added to what lies under it (one with no opacity of its own).

use eframe::egui::epaint::text::{LayoutJob, TextFormat};
use eframe::egui::epaint::{
    ClippedShape, ColorMode, Mesh, PathStroke, RectShape, Shape, StrokeKind, TextShape,
};
use eframe::egui::{Color32, FontFamily, Pos2, Rect, Stroke, TextureId, Vec2};
use std::fmt::Write;

/// The fonts named for text in the two families, the program's own first (`ui/fonts.rs`).
const PROPORTIONAL: &str = "Inter, 'Helvetica Neue', Helvetica, Arial, sans-serif";
const MONOSPACE: &str = "'JetBrains Mono', Menlo, Consolas, monospace";
/// Triangles of one colour are given an edge of that colour this wide, so that no hairline shows between them.
const SEAM: f32 = 0.3;

/// A number as SVG takes it: two decimals at most, and no zeros after them. What is no number (a view that divided by
/// nothing) is written as 0: `NaN` or `inf` in a file would make it one that no program opens.
fn n(x: f32) -> String {
    if !x.is_finite() {
        return "0".to_string();
    }
    let text = format!("{x:.2}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    if text == "-0" {
        "0".to_string()
    } else {
        text.to_string()
    }
}

/// ` fill="#rrggbb"` for a colour (with its opacity when it is see-through), ` fill="none"` for none; `attribute` is `fill` or
/// `stroke`.
fn paint(attribute: &str, colour: Color32) -> String {
    let [r, g, b, a] = colour.to_srgba_unmultiplied();
    if a == 0 {
        return format!(" {attribute}=\"none\"");
    }
    let opacity = if a < 255 {
        format!(" {attribute}-opacity=\"{}\"", n(f32::from(a) / 255.0))
    } else {
        String::new()
    };
    format!(" {attribute}=\"#{r:02x}{g:02x}{b:02x}\"{opacity}")
}

/// The attributes of an outline, or nothing when there is none (SVG draws no outline unless it is asked to).
fn outline(stroke: Stroke) -> String {
    if stroke.width <= 0.0 || stroke.color.a() == 0 {
        return String::new();
    }
    format!(
        "{} stroke-width=\"{}\"",
        paint("stroke", stroke.color),
        n(stroke.width)
    )
}

/// The outline of a path as a plain one: its colour is the one colour it has, or the colour at its first point.
fn plain(stroke: &PathStroke, points: &[Pos2]) -> Stroke {
    let colour = match &stroke.color {
        ColorMode::Solid(colour) => *colour,
        ColorMode::UV(at) => points
            .first()
            .map_or(Color32::TRANSPARENT, |p| at(Rect::from_points(points), *p)),
    };
    Stroke::new(stroke.width, colour)
}

/// `x,y x,y …`
fn list(points: &[Pos2]) -> String {
    points
        .iter()
        .map(|p| format!("{},{}", n(p.x), n(p.y)))
        .collect::<Vec<_>>()
        .join(" ")
}

/// What must stand for a character in the text of an element. A character that XML does not allow at all (a control
/// character from a name in somebody's file) is left out: one of them makes the whole file unreadable.
fn escaped(text: &str) -> String {
    text.chars()
        .filter(|c| *c == '\t' || !(c.is_control() || matches!(*c, '\u{FFFE}' | '\u{FFFF}')))
        .collect::<String>()
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn rect(out: &mut String, r: &RectShape) {
    // a soft shadow and a picture are not drawn
    if r.blur_width > 0.0 || r.brush.is_some() {
        return;
    }
    // SVG draws an outline astride the edge; an outline that is not drawn moves nothing
    let drawn = outline(r.stroke);
    let grow = match r.stroke_kind {
        _ if drawn.is_empty() => 0.0,
        StrokeKind::Inside => -r.stroke.width / 2.0,
        StrokeKind::Middle => 0.0,
        StrokeKind::Outside => r.stroke.width / 2.0,
    };
    let b = r.rect.expand(grow);
    if !(b.width() > 0.0 && b.height() > 0.0) {
        return;
    }
    let looks = format!("{}{drawn}", paint("fill", r.fill));
    let turned = if r.angle == 0.0 {
        String::new()
    } else {
        format!(
            " transform=\"rotate({} {} {})\"",
            n(r.angle.to_degrees()),
            n(b.center().x),
            n(b.center().y)
        )
    };
    let c = r.corner_radius;
    if c.nw == c.ne && c.ne == c.se && c.se == c.sw {
        let round = if c.nw == 0 {
            String::new()
        } else {
            format!(" rx=\"{}\"", c.nw)
        };
        let _ = writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\"{round}{looks}{turned}/>",
            n(b.min.x),
            n(b.min.y),
            n(b.width()),
            n(b.height())
        );
        return;
    }
    // corners of their own: along the top, down the right, back along the bottom and up the left, an arc at each corner
    let (l, t, rt, bt) = (b.min.x, b.min.y, b.max.x, b.max.y);
    let [nw, ne, se, sw] = [c.nw, c.ne, c.se, c.sw].map(f32::from);
    let arc =
        |radius: f32, x: f32, y: f32| format!("A{0} {0} 0 0 1 {1} {2}", n(radius), n(x), n(y));
    let _ = writeln!(
        out,
        "<path d=\"M{} {}H{}{}V{}{}H{}{}V{}{}Z\"{looks}{turned}/>",
        n(l + nw),
        n(t),
        n(rt - ne),
        arc(ne, rt, t + ne),
        n(bt - se),
        arc(se, rt - se, bt),
        n(l + sw),
        arc(sw, l, bt - sw),
        n(t + nw),
        arc(nw, l + nw, t),
    );
}

/// The triangles of a mesh, those that follow each other in one colour as one path. A triangle whose corners differ in colour
/// gets their mean: what a view draws with a mesh here are areas of one colour (the features of the circular map).
fn mesh(out: &mut String, m: &Mesh) {
    // a mesh with a picture on it is not drawn
    if m.texture_id != TextureId::default() {
        return;
    }
    let mut areas: Vec<(Color32, String)> = Vec::new();
    for triangle in m.indices.as_chunks::<3>().0 {
        let corners: Vec<_> = triangle
            .iter()
            .filter_map(|&i| m.vertices.get(i as usize))
            .collect();
        let [a, b, c] = corners[..] else {
            continue;
        };
        let mean = |channel: fn(&Color32) -> u8| {
            ((u16::from(channel(&a.color))
                + u16::from(channel(&b.color))
                + u16::from(channel(&c.color)))
                / 3) as u8
        };
        let colour = Color32::from_rgba_premultiplied(
            mean(Color32::r),
            mean(Color32::g),
            mean(Color32::b),
            mean(Color32::a),
        );
        let d =
            format!("M{}L{}L{}Z", list(&[a.pos]), list(&[b.pos]), list(&[c.pos])).replace(',', " ");
        match areas.last_mut() {
            Some((last, path)) if *last == colour => path.push_str(&d),
            _ => areas.push((colour, d)),
        }
    }
    for (colour, d) in areas {
        if colour.a() == 0 {
            continue;
        }
        // an edge of the same colour closes the hairlines between triangles; a see-through colour would show it twice
        let seam = if colour.is_opaque() {
            format!(
                "{} stroke-width=\"{}\" stroke-linejoin=\"round\"",
                paint("stroke", colour),
                n(SEAM)
            )
        } else {
            String::new()
        };
        let _ = writeln!(out, "<path d=\"{d}\"{}{seam}/>", paint("fill", colour));
    }
}

/// The section of `job` that the character at byte `at` is in.
fn section_at(job: &LayoutJob, at: usize) -> usize {
    job.sections
        .iter()
        .position(|s| (s.byte_range.start.0..s.byte_range.end.0).contains(&at))
        .unwrap_or(0)
}

/// The section of the character that a glyph draws: the next of the text's characters that is the glyph's, from where the
/// glyph before it was found. A glyph that draws none of the characters that are left (the dots that end a text cut short)
/// keeps the section of the glyph before it, and the place in the text stays where it was for the glyphs that follow.
fn section_of(
    chr: char,
    characters: &mut std::str::CharIndices,
    job: &LayoutJob,
    before: usize,
) -> usize {
    let mut ahead = characters.clone();
    match ahead.find(|(_, c)| *c == chr) {
        Some((at, _)) => {
            *characters = ahead;
            section_at(job, at)
        }
        None => before,
    }
}

/// A stretch of a row in one font and colour.
struct Run {
    section: usize,
    /// Where it starts and ends along the row, and its baseline.
    from: f32,
    to: f32,
    baseline: f32,
    /// How far the font reaches above the baseline, and the height of its line.
    ascent: f32,
    height: f32,
    text: String,
}

fn family(format: &TextFormat) -> String {
    match &format.font_id.family {
        FontFamily::Proportional => PROPORTIONAL.to_string(),
        FontFamily::Monospace => MONOSPACE.to_string(),
        // a name goes between quotes in an attribute: nothing in it may end either
        FontFamily::Name(name) => format!(
            "'{}', sans-serif",
            name.replace(['"', '\'', '<', '>', '&'], "")
        ),
    }
}

fn run(out: &mut String, t: &TextShape, r: &Run) {
    let Some(format) = t.galley.job.sections.get(r.section).map(|s| &s.format) else {
        return;
    };
    let dim = |colour: Color32| {
        if t.opacity_factor < 1.0 {
            colour.gamma_multiply(t.opacity_factor)
        } else {
            colour
        }
    };
    if format.background.a() > 0 {
        let _ = writeln!(
            out,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\"{}/>",
            n(r.from),
            n(r.baseline - r.ascent),
            n(r.to - r.from),
            n(r.height),
            paint("fill", dim(format.background))
        );
    }
    let colour = dim(match t.override_text_color {
        Some(colour) => colour,
        None if format.color == Color32::PLACEHOLDER => t.fallback_color,
        None => format.color,
    });
    // spaces alone draw nothing
    if !r.text.trim().is_empty() && colour.a() > 0 {
        let weight = format
            .coords
            .as_ref()
            .iter()
            .find(|(tag, _)| tag.to_be_bytes() == *b"wght")
            .map_or_else(String::new, |(_, w)| format!(" font-weight=\"{}\"", n(*w)));
        let slant = if format.italics {
            " font-style=\"italic\""
        } else {
            ""
        };
        // the width it has on screen, so that another font does not make it run into its neighbours
        let width = if r.text.chars().count() > 1 {
            format!(" textLength=\"{}\"", n(r.to - r.from))
        } else {
            String::new()
        };
        let _ = writeln!(
            out,
            "<text x=\"{}\" y=\"{}\" font-family=\"{}\" font-size=\"{}\"{weight}{slant}{}{width}>{}</text>",
            n(r.from),
            n(r.baseline),
            family(format),
            n(format.font_id.size),
            paint("fill", colour),
            escaped(&r.text)
        );
    }
    let mut rule = |stroke: Stroke, y: f32| {
        if stroke.width > 0.0 && stroke.color.a() > 0 {
            let _ = writeln!(
                out,
                "<line x1=\"{}\" y1=\"{1}\" x2=\"{2}\" y2=\"{1}\"{3}/>",
                n(r.from),
                n(y),
                n(r.to),
                outline(stroke)
            );
        }
    };
    let under = if format.underline.width > 0.0 {
        format.underline
    } else {
        t.underline
    };
    rule(under, r.baseline + (r.height - r.ascent) / 2.0);
    rule(format.strikethrough, r.baseline - r.ascent / 3.0);
}

fn text(out: &mut String, t: &TextShape) {
    let job = &t.galley.job;
    let mut made = String::new();
    // the glyphs come in the order of the text's characters, and a character's place in the text says which section it is in
    let mut characters = job.text.char_indices();
    let mut section = 0;
    for placed in &t.galley.rows {
        let origin = t.pos + placed.pos.to_vec2();
        let mut open: Option<Run> = None;
        for glyph in &placed.row.glyphs {
            section = section_of(glyph.chr, &mut characters, job, section);
            let (from, baseline) = (origin.x + glyph.pos.x, origin.y + glyph.pos.y);
            match &mut open {
                Some(r) if r.section == section => {
                    r.text.push(glyph.chr);
                    r.to = from + glyph.advance_width;
                }
                _ => {
                    if let Some(done) = open.take() {
                        run(&mut made, t, &done);
                    }
                    open = Some(Run {
                        section,
                        from,
                        to: from + glyph.advance_width,
                        baseline,
                        ascent: glyph.font_ascent,
                        height: glyph.font_height,
                        text: glyph.chr.to_string(),
                    });
                }
            }
        }
        if let Some(done) = open {
            run(&mut made, t, &done);
        }
    }
    if made.is_empty() {
        return;
    }
    if t.angle == 0.0 {
        out.push_str(&made);
    } else {
        // egui turns a text about the place it is put at
        let _ = writeln!(
            out,
            "<g transform=\"rotate({} {} {})\">\n{made}</g>",
            n(t.angle.to_degrees()),
            n(t.pos.x),
            n(t.pos.y)
        );
    }
}

fn shape(out: &mut String, s: &Shape) {
    match s {
        Shape::Noop | Shape::Callback(_) => {}
        Shape::Vec(shapes) => shapes.iter().for_each(|s| shape(out, s)),
        Shape::Circle(c) => {
            let _ = writeln!(
                out,
                "<circle cx=\"{}\" cy=\"{}\" r=\"{}\"{}{}/>",
                n(c.center.x),
                n(c.center.y),
                n(c.radius),
                paint("fill", c.fill),
                outline(c.stroke)
            );
        }
        Shape::Ellipse(e) => {
            let turned = if e.angle == 0.0 {
                String::new()
            } else {
                format!(
                    " transform=\"rotate({} {} {})\"",
                    n(e.angle.to_degrees()),
                    n(e.center.x),
                    n(e.center.y)
                )
            };
            let _ = writeln!(
                out,
                "<ellipse cx=\"{}\" cy=\"{}\" rx=\"{}\" ry=\"{}\"{}{}{turned}/>",
                n(e.center.x),
                n(e.center.y),
                n(e.radius.x),
                n(e.radius.y),
                paint("fill", e.fill),
                outline(e.stroke)
            );
        }
        Shape::LineSegment { points, stroke } => {
            let drawn = outline(*stroke);
            if !drawn.is_empty() {
                let _ = writeln!(
                    out,
                    "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\"{drawn}/>",
                    n(points[0].x),
                    n(points[0].y),
                    n(points[1].x),
                    n(points[1].y)
                );
            }
        }
        Shape::Path(p) => {
            if p.points.len() < 2 {
                return;
            }
            // egui fills a path only when it is closed
            let (element, fill) = if p.closed {
                ("polygon", p.fill)
            } else {
                ("polyline", Color32::TRANSPARENT)
            };
            let _ = writeln!(
                out,
                "<{element} points=\"{}\"{}{}/>",
                list(&p.points),
                paint("fill", fill),
                outline(plain(&p.stroke, &p.points))
            );
        }
        Shape::Rect(r) => rect(out, r),
        Shape::Text(t) => text(out, t),
        Shape::Mesh(m) => mesh(out, m),
        Shape::QuadraticBezier(b) => curve(out, 'Q', &b.points, b.closed, b.fill, &b.stroke),
        Shape::CubicBezier(b) => curve(out, 'C', &b.points, b.closed, b.fill, &b.stroke),
    }
}

/// A curve from its first point by its control points to its last (`Q` with one control point, `C` with two).
fn curve(
    out: &mut String,
    kind: char,
    points: &[Pos2],
    closed: bool,
    fill: Color32,
    stroke: &PathStroke,
) {
    let fill = if closed { fill } else { Color32::TRANSPARENT };
    let _ = writeln!(
        out,
        "<path d=\"M{} {kind}{}{}\"{}{}/>",
        list(&points[..1]).replace(',', " "),
        list(&points[1..]).replace(',', " "),
        if closed { "Z" } else { "" },
        paint("fill", fill),
        outline(plain(stroke, points))
    );
}

/// What stands in the names of a figure's cuts for the figure's own mark, while the figure is written: no text of a figure
/// can hold it (`escaped` leaves control characters out).
const MARK: char = '\u{1}';
/// How much of the checksum of a figure is its mark.
const MARK_LENGTH: usize = 10;

/// The picture of `shapes`, drawn in a room of `size` points on a ground of `ground`. A shape that its view cut to a part of the
/// room (a list that scrolls) is cut to the same part. The cuts are named after the figure (a checksum of what it draws): in a
/// page that holds several figures a name is looked up in the whole page, and with one name for the first cut of every figure
/// the cut of the first figure would be taken for the others'.
pub fn to_svg(shapes: &[ClippedShape], size: Vec2, ground: Color32) -> String {
    let whole = Rect::from_min_size(Pos2::ZERO, size);
    let mut cuts: Vec<Rect> = Vec::new();
    let mut body = String::new();
    // the shapes that follow each other with one cut are one group
    let mut group: Option<(Option<usize>, String)> = None;
    let close = |body: &mut String, group: Option<(Option<usize>, String)>| match group {
        Some((_, drawn)) if drawn.is_empty() => {}
        Some((Some(cut), drawn)) => {
            let _ = writeln!(body, "<g clip-path=\"url(#cut{MARK}{cut})\">\n{drawn}</g>");
        }
        Some((None, drawn)) => body.push_str(&drawn),
        None => {}
    };
    for clipped in shapes {
        let cut = if clipped.clip_rect.contains_rect(whole) {
            None
        } else {
            let part = clipped.clip_rect.intersect(whole);
            let known = cuts.iter().position(|c| *c == part);
            Some(known.unwrap_or_else(|| {
                cuts.push(part);
                cuts.len() - 1
            }))
        };
        if group.as_ref().is_none_or(|(open, _)| *open != cut) {
            close(&mut body, group.take());
            group = Some((cut, String::new()));
        }
        if let Some((_, drawn)) = &mut group {
            shape(drawn, &clipped.shape);
        }
    }
    close(&mut body, group);
    let mut out = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{0}\" height=\"{1}\" viewBox=\"0 0 {0} {1}\" xml:space=\"preserve\">\n",
        n(size.x),
        n(size.y)
    );
    if !cuts.is_empty() {
        out.push_str("<defs>\n");
        for (k, cut) in cuts.iter().enumerate() {
            let _ = writeln!(
                out,
                "<clipPath id=\"cut{MARK}{k}\"><rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\"/></clipPath>",
                n(cut.min.x),
                n(cut.min.y),
                n(cut.width().max(0.0)),
                n(cut.height().max(0.0))
            );
        }
        out.push_str("</defs>\n");
    }
    if ground.a() > 0 {
        let _ = writeln!(
            out,
            "<rect width=\"{}\" height=\"{}\"{}/>",
            n(size.x),
            n(size.y),
            paint("fill", ground)
        );
    }
    out.push_str(&body);
    out.push_str("</svg>\n");
    if cuts.is_empty() {
        return out;
    }
    let own = crate::pipeline::flow::digest(out.as_bytes());
    out.replace(MARK, &format!("{}-", &own[..MARK_LENGTH.min(own.len())]))
}

#[cfg(test)]
mod tests;
