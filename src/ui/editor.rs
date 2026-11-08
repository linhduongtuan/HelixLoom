//! The sequence view: bases with their complement, feature arrows in lanes, translation rows,
//! selection, find and enzyme highlights, plus mouse and keyboard editing.

use super::clipboard;
use super::layout::Layout;
use super::state::{DocView, Prefs, TransMode};
use super::theme::{self, Palette};
use crate::seq::dna::{complement, revcom, translate};
use crate::seq::features::Feature;
use eframe::egui::{
    self, pos2, vec2, Align, Align2, Color32, Event, FontId, Key, Painter, PointerButton, Rect,
    Sense, Shape, Stroke, Ui,
};

/// Characters the original lets you type into a sequence (allowedkeys).
const ALLOWED_KEYS: &str = "ACGTNacgtn*";

/// The original's text area starts 75 characters wide (`default_textarea_width`, AppMain.tcl:66907) and keeps that width
/// until its window has a size; after that a row holds as many characters as fit (`textarea_width_changed`).
pub const TEXT_COLUMNS: usize = 75;

/// Room left for the numbers on the left and the scroll bar on the right, in points.
const SIDE_ROOM: f32 = 28.0;

/// Fewest characters in a row. The original has no lower limit; a row narrower than this cannot be read.
const MIN_COLUMNS: usize = 10;

/// How many characters a row holds in a view `avail_w` wide whose numbers take `gutter`: all that fit, with no upper limit
/// (the original sets the width again from its window's width at every resize, `width = int((winfo width - 4) /
/// character width)`, AppMain.tcl:62344 `proc textarea_width_changed`, bound to `<Configure>` at L2067). A view with no
/// size yet keeps the starting width.
fn row_columns(avail_w: f32, gutter: f32, char_width: f32) -> usize {
    let room = avail_w - gutter - SIDE_ROOM;
    if !room.is_finite() || char_width <= 0.0 {
        return TEXT_COLUMNS;
    }
    ((room.max(0.0) / char_width) as usize).max(MIN_COLUMNS)
}

pub struct Env<'a> {
    pub pal: &'a Palette,
    pub prefs: &'a Prefs,
    /// Recognition sites to highlight, already split at the origin.
    pub highlight: &'a [(usize, usize)],
    pub code: &'a str,
}

struct Metrics {
    font: FontId,
    small: FontId,
    cw: f32,
    gutter: f32,
}

pub fn show(ui: &mut Ui, dv: &mut DocView, env: &Env) {
    let font = FontId::monospace(env.prefs.font_size);
    let cw = ui.ctx().fonts_mut(|f| f.glyph_width(&font, 'A'));
    let m = Metrics {
        small: FontId::monospace(env.prefs.font_size * 0.72),
        font,
        cw,
        gutter: cw * 7.0,
    };
    let line_h = (env.prefs.font_size * 1.3).round();
    let lane_h = (env.prefs.font_size + 5.0).round();
    let avail_w = ui.available_width();
    let per_row = row_columns(avail_w, m.gutter, cw);
    dv.per_row = per_row;
    let layout = Layout::new(
        &dv.doc,
        per_row,
        line_h,
        lane_h,
        env.prefs.second_strand,
        dv.trans.rows(),
    );

    ruler(ui, &m, &layout, env, avail_w);
    let target = dv.scroll_to.take();
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show_viewport(ui, |ui, viewport| {
            let (rect, resp) = ui.allocate_exact_size(
                vec2(avail_w, layout.total_height() + 8.0),
                Sense::click_and_drag(),
            );
            let painter = ui.painter().with_clip_rect(rect.intersect(ui.clip_rect()));
            painter.rect_filled(rect, 0.0, env.pal.text_bg);
            let (first, last) = (
                layout.row_at_y(viewport.min.y),
                layout.row_at_y(viewport.max.y),
            );
            for r in first..=last {
                paint_row(&painter, rect, &m, &layout, dv, env, r);
            }
            if let Some(pos) = target {
                let r = layout.row_of_base(pos);
                let top = rect.top() + layout.row_y[r];
                ui.scroll_to_rect(
                    Rect::from_min_size(
                        pos2(rect.left(), top),
                        vec2(8.0, layout.row_y[r + 1] - layout.row_y[r]),
                    ),
                    None,
                );
            }
            interact(ui, &resp, rect, &m, &layout, dv, env);
        });
}

/// How many tens the ruler above a row of `columns` characters writes: it adds a ten while the string is shorter than the
/// row, so the last one may reach past the row's edge and be cut short.
fn ruler_tens(columns: usize) -> usize {
    columns.div_ceil(10)
}

/// The numbers above the sequence. The original writes ten characters for every ten bases (`    *` and the number, right
/// aligned) until the string is as long as a row (`textarea_width_changed`, AppMain.tcl:62359-62364), so the last ten may
/// be cut short by the edge of the row: at 75 columns the star of the eighth ten shows and its number does not.
fn ruler(ui: &mut Ui, m: &Metrics, layout: &Layout, env: &Env, width: f32) {
    let (rect, _) = ui.allocate_exact_size(vec2(width, layout.line_h), Sense::hover());
    let x0 = rect.left() + m.gutter;
    ui.painter().rect_filled(rect, 0.0, env.pal.window_bg);
    let row = Rect::from_min_max(
        pos2(x0, rect.top()),
        pos2(x0 + layout.per_row as f32 * m.cw, rect.bottom()),
    );
    let p = ui.painter().with_clip_rect(row.intersect(rect));
    for k in 1..=ruler_tens(layout.per_row) {
        let end = x0 + (k * 10) as f32 * m.cw;
        p.text(
            pos2(end, rect.top()),
            Align2::RIGHT_TOP,
            (k * 10).to_string(),
            m.small.clone(),
            env.pal.dim,
        );
        p.text(
            pos2(end - 5.5 * m.cw, rect.top()),
            Align2::CENTER_TOP,
            "*",
            m.small.clone(),
            env.pal.dim,
        );
    }
}

/// Part of `[a, b)` inside the row, as columns.
fn clip(a: usize, b: usize, rs: usize, re: usize) -> Option<(usize, usize)> {
    let (a, b) = (a.max(rs), b.min(re));
    (a < b).then(|| (a - rs, b - rs))
}

fn paint_row(
    p: &Painter,
    rect: Rect,
    m: &Metrics,
    layout: &Layout,
    dv: &DocView,
    env: &Env,
    r: usize,
) {
    let (rs, re) = layout.row_range(r);
    let x0 = rect.left() + m.gutter;
    let y = rect.top() + layout.row_y[r];
    let seq_h = layout.line_h * if layout.second_strand { 2.0 } else { 1.0 };
    p.text(
        pos2(x0 - 10.0, y + 2.0),
        Align2::RIGHT_TOP,
        (rs + 1).to_string(),
        m.small.clone(),
        env.pal.dim,
    );

    let band = |a: usize, b: usize, color: Color32| {
        if let Some((c0, c1)) = clip(a, b, rs, re) {
            p.rect_filled(
                Rect::from_min_max(
                    pos2(x0 + c0 as f32 * m.cw, y),
                    pos2(x0 + c1 as f32 * m.cw, y + seq_h),
                ),
                0.0,
                color,
            );
        }
    };
    let n = dv.doc.len();
    for h in &dv.find_hits {
        let color = if h.reverse {
            env.pal.find_rev
        } else {
            env.pal.find_fwd
        };
        if h.start + h.len > n {
            band(h.start, n, color);
            band(0, h.start + h.len - n, color);
        } else {
            band(h.start, h.start + h.len, color);
        }
    }
    for &(a, b) in env.highlight {
        band(a, b, env.pal.enzyme.gamma_multiply(0.30));
    }
    let selection = dv.selection();
    if let Some((a, b)) = selection {
        band(a, b, env.pal.select_bg);
    }

    let seq = &dv.doc.seq[rs..re];
    segments(
        p,
        m,
        pos2(x0, y),
        seq,
        selection.and_then(|(a, b)| clip(a, b, rs, re)),
        env.pal.fg,
        env.pal.select_fg,
    );
    if layout.second_strand {
        let comp: Vec<u8> = seq.iter().map(|&c| complement(c)).collect();
        let dim = env.pal.dim;
        let sel = selection.and_then(|(a, b)| clip(a, b, rs, re));
        segments(
            p,
            m,
            pos2(x0, y + layout.line_h),
            &comp,
            sel,
            dim,
            env.pal.select_fg,
        );
    }
    for (fi, f) in dv.doc.features.iter().enumerate() {
        draw_feature(
            p,
            rect,
            m,
            layout,
            &dv.doc.features,
            fi,
            f,
            n,
            dv.selected_feature == Some(fi),
            env,
            r,
        );
    }
    if dv.trans != TransMode::Off {
        draw_translation(p, rect, m, layout, dv, env, r);
    }
    let caret_row = layout.row_of_base(dv.caret);
    let at_row_end = dv.caret == re && re < n && r + 1 < layout.rows;
    if caret_row == r && !at_row_end && selection.is_none() {
        let cx = x0 + (dv.caret - rs) as f32 * m.cw;
        p.line_segment(
            [pos2(cx, y), pos2(cx, y + layout.line_h)],
            Stroke::new(1.6, env.pal.fg),
        );
    }
}

/// Draw `text` with the selected columns in the selection colour.
fn segments(
    p: &Painter,
    m: &Metrics,
    origin: egui::Pos2,
    text: &[u8],
    sel: Option<(usize, usize)>,
    fg: Color32,
    sel_fg: Color32,
) {
    let pieces: [(usize, usize, Color32); 3] = match sel {
        Some((a, b)) => [(0, a, fg), (a, b, sel_fg), (b, text.len(), fg)],
        None => [(0, text.len(), fg), (0, 0, fg), (0, 0, fg)],
    };
    for (a, b, color) in pieces {
        if a < b {
            let s = String::from_utf8_lossy(&text[a..b]).into_owned();
            p.text(
                origin + vec2(a as f32 * m.cw, 0.0),
                Align2::LEFT_TOP,
                s,
                m.font.clone(),
                color,
            );
        }
    }
}

/// Which span of a feature carries the arrow head, and whether it is on the right.
pub fn head_of(f: &Feature, n: usize) -> (usize, bool) {
    let last = f.spans.len().saturating_sub(1);
    let wraps = f.spans.len() > 1 && f.spans[0].0 == 0 && f.spans[last].1 == n;
    match (f.reverse, wraps) {
        (false, false) => (last, true),
        (false, true) => (0, true),
        (true, false) => (0, false),
        (true, true) => (last, false),
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_feature(
    p: &Painter,
    rect: Rect,
    m: &Metrics,
    layout: &Layout,
    _all: &[Feature],
    fi: usize,
    f: &Feature,
    n: usize,
    selected: bool,
    env: &Env,
    r: usize,
) {
    let (rs, re) = layout.row_range(r);
    let (head_span, head_right) = head_of(f, n);
    let y = rect.top() + layout.lanes_y(r) + layout.lanes[fi] as f32 * layout.lane_h;
    let (top, bottom) = (y + 2.0, y + layout.lane_h - 2.0);
    let fill = theme::rgb(f.color());
    let outline = if selected {
        env.pal.accent
    } else {
        theme::mix(Color32::BLACK, fill, 35.0)
    };
    for (si, &(s, e)) in f.spans.iter().enumerate() {
        let Some((c0, c1)) = clip(s, e, rs, re) else {
            continue;
        };
        let x0 = rect.left() + m.gutter;
        let (xa, xb) = (x0 + c0 as f32 * m.cw, x0 + c1 as f32 * m.cw);
        let has_head = si == head_span && ((head_right && e <= re) || (!head_right && s >= rs));
        let hw = ((bottom - top) * 0.5).min(xb - xa);
        let mid = (top + bottom) / 2.0;
        let pts = match (has_head, head_right) {
            (false, _) => vec![
                pos2(xa, top),
                pos2(xb, top),
                pos2(xb, bottom),
                pos2(xa, bottom),
            ],
            (true, true) => vec![
                pos2(xa, top),
                pos2(xb - hw, top),
                pos2(xb, mid),
                pos2(xb - hw, bottom),
                pos2(xa, bottom),
            ],
            (true, false) => vec![
                pos2(xa + hw, top),
                pos2(xb, top),
                pos2(xb, bottom),
                pos2(xa + hw, bottom),
                pos2(xa, mid),
            ],
        };
        p.add(Shape::convex_polygon(
            pts,
            fill,
            Stroke::new(if selected { 2.0 } else { 1.0 }, outline),
        ));
        let label_rect = Rect::from_min_max(pos2(xa + 2.0, top), pos2(xb - 2.0, bottom));
        p.with_clip_rect(label_rect).text(
            pos2(xa + 4.0, mid),
            Align2::LEFT_CENTER,
            &f.name,
            FontId::proportional((env.prefs.font_size * 0.8).max(10.0)),
            theme::text_on(fill),
        );
    }
}

fn draw_translation(
    p: &Painter,
    rect: Rect,
    m: &Metrics,
    layout: &Layout,
    dv: &DocView,
    env: &Env,
    r: usize,
) {
    let (rs, re) = layout.row_range(r);
    let x0 = rect.left() + m.gutter;
    let y0 = rect.top() + layout.trans_y(r);
    let line = layout.line_h * 0.9;
    let n = dv.doc.len();
    let row_font = FontId::monospace(env.prefs.font_size * 0.85);
    for k in 0..dv.trans.rows() {
        let reverse = k >= 3;
        let frame = k % 3;
        let y = y0 + k as f32 * line;
        let mut pos = rs.saturating_sub(2);
        while pos < re {
            // codon start on the top strand: forward frames count from the origin, reverse from the end
            let ok = if reverse {
                (n.wrapping_sub(pos + 3)) % 3 == frame
            } else {
                pos % 3 == frame
            };
            if ok && pos + 3 <= n && pos + 1 >= rs && pos + 1 < re {
                let codon = &dv.doc.seq[pos..pos + 3];
                let aa = if reverse {
                    translate(&revcom(codon), env.code)
                } else {
                    translate(codon, env.code)
                };
                let aa = aa.first().copied().unwrap_or(b'X');
                let color = match aa {
                    b'*' => env.pal.enzyme,
                    b'M' => Color32::from_rgb(0, 150, 60),
                    _ => env.pal.dim,
                };
                let cx = x0 + ((pos + 1 - rs) as f32 + 0.5) * m.cw;
                let label = if env.prefs.three_letter {
                    crate::seq::dna::three_letter(aa).to_string()
                } else {
                    (aa as char).to_string()
                };
                p.text(
                    pos2(cx, y),
                    Align2::CENTER_TOP,
                    label,
                    row_font.clone(),
                    color,
                );
            }
            pos += 1;
        }
    }
}

// ---------------------------------------------------------------------------
// Interaction
// ---------------------------------------------------------------------------

fn base_at(
    layout: &Layout,
    m: &Metrics,
    local: egui::Vec2,
    nearest_boundary: bool,
) -> Option<usize> {
    let r = layout.row_at_y(local.y);
    let (rs, re) = layout.row_range(r);
    let col = (local.x - m.gutter) / m.cw;
    if nearest_boundary {
        Some(rs + col.round().clamp(0.0, (re - rs) as f32) as usize)
    } else if col >= 0.0
        && (col as usize) < re - rs
        && local.y - layout.row_y[r] < layout.line_h * 2.0
    {
        Some(rs + col as usize)
    } else {
        None
    }
}

fn feature_at(layout: &Layout, m: &Metrics, dv: &DocView, local: egui::Vec2) -> Option<usize> {
    let r = layout.row_at_y(local.y);
    let (rs, re) = layout.row_range(r);
    let rel = local.y - layout.lanes_y(r);
    if rel < 0.0 {
        return None;
    }
    let lane = (rel / layout.lane_h) as usize;
    let base = rs + ((local.x - m.gutter) / m.cw).floor().max(0.0) as usize;
    dv.doc
        .features
        .iter()
        .enumerate()
        .find_map(|(i, f)| (layout.lanes[i] == lane && base < re && f.contains(base)).then_some(i))
}

fn interact(
    ui: &mut Ui,
    resp: &egui::Response,
    rect: Rect,
    m: &Metrics,
    layout: &Layout,
    dv: &mut DocView,
    env: &Env,
) {
    let shift = ui.input(|i| i.modifiers.shift);
    if let Some(pos) = resp.hover_pos() {
        dv.hover = base_at(layout, m, pos - rect.min, false);
    } else {
        dv.hover = None;
    }
    let pressed = ui.input(|i| i.pointer.button_pressed(PointerButton::Primary));
    if pressed && resp.hovered() {
        resp.request_focus();
        let local = ui.input(|i| i.pointer.interact_pos()).map(|p| p - rect.min);
        if let Some(local) = local {
            if let Some(fi) = feature_at(layout, m, dv, local) {
                let f = &dv.doc.features[fi];
                let (a, b) = (f.start(), f.end());
                dv.selected_feature = Some(fi);
                dv.select(a, b);
            } else if let Some(b) = base_at(layout, m, local, true) {
                dv.selected_feature = None;
                if shift {
                    dv.anchor.get_or_insert(dv.caret);
                } else {
                    dv.anchor = Some(b);
                }
                dv.caret = b;
            }
        }
    } else if resp.dragged_by(PointerButton::Primary) {
        if let Some(local) = resp.interact_pointer_pos().map(|p| p - rect.min) {
            if let Some(b) = base_at(layout, m, local, true) {
                dv.caret = b;
            }
        }
    }
    if !resp.has_focus() && ui.memory(|mem| mem.focused().is_none()) {
        resp.request_focus();
    }
    if resp.has_focus() {
        ui.memory_mut(|mem| {
            mem.set_focus_lock_filter(
                resp.id,
                egui::EventFilter {
                    horizontal_arrows: true,
                    vertical_arrows: true,
                    tab: false,
                    escape: false,
                },
            )
        });
        let events = ui.input(|i| i.events.clone());
        for ev in events {
            keyboard(ui, ev, layout, dv, env);
        }
    }
}

fn move_caret(dv: &mut DocView, to: usize, extend: bool) {
    let to = to.min(dv.doc.len());
    if extend {
        dv.anchor.get_or_insert(dv.caret);
    } else {
        dv.anchor = None;
    }
    dv.caret = to;
    dv.scroll_to = Some(to);
}

fn keyboard(ui: &Ui, ev: Event, layout: &Layout, dv: &mut DocView, env: &Env) {
    let n = dv.doc.len();
    match ev {
        Event::Text(t) => {
            let typed: Vec<u8> = t
                .chars()
                .filter(|c| ALLOWED_KEYS.contains(*c))
                .map(|c| c as u8)
                .collect();
            if !typed.is_empty() {
                dv.replace_selection(&typed);
            }
        }
        Event::Paste(s) => clipboard::paste(ui.ctx(), dv, &s, false),
        Event::Copy => clipboard::copy(ui.ctx(), dv, false),
        Event::Cut => {
            clipboard::copy(ui.ctx(), dv, false);
            if dv.selection().is_some() {
                dv.replace_selection(b"");
            }
        }
        Event::Key {
            key,
            pressed: true,
            modifiers: mods,
            ..
        } => {
            let extend = mods.shift;
            let step = if mods.alt { 10 } else { 1 };
            match key {
                Key::ArrowLeft => match dv.selection() {
                    Some((a, _)) if !extend => move_caret(dv, a, false),
                    _ => move_caret(dv, dv.caret.saturating_sub(step), extend),
                },
                Key::ArrowRight => match dv.selection() {
                    Some((_, b)) if !extend => move_caret(dv, b, false),
                    _ => move_caret(dv, dv.caret + step, extend),
                },
                Key::ArrowUp => move_caret(dv, dv.caret.saturating_sub(layout.per_row), extend),
                Key::ArrowDown => move_caret(dv, dv.caret + layout.per_row, extend),
                Key::Home if mods.command => move_caret(dv, 0, extend),
                Key::End if mods.command => move_caret(dv, n, extend),
                Key::Home => {
                    move_caret(dv, layout.row_range(layout.row_of_base(dv.caret)).0, extend)
                }
                Key::End => {
                    move_caret(dv, layout.row_range(layout.row_of_base(dv.caret)).1, extend)
                }
                Key::Backspace | Key::Delete => delete_key(dv, key == Key::Backspace),
                Key::A if mods.command => dv.select(0, n),
                Key::Z if mods.command && mods.shift => dv.redo(),
                Key::Z if mods.command => dv.undo(),
                _ => {}
            }
        }
        _ => {}
    }
    let _ = env;
}

fn delete_key(dv: &mut DocView, backward: bool) {
    if dv.selection().is_some() {
        dv.replace_selection(b"");
        return;
    }
    let c = dv.caret;
    let (a, b) = if backward {
        (c.saturating_sub(1), c)
    } else {
        (c, (c + 1).min(dv.doc.len()))
    };
    if a < b {
        dv.edit(|d| d.delete(a, b));
        dv.place_caret(a);
        dv.scroll_to = Some(a);
    }
}

#[allow(dead_code)]
pub fn scroll_hint(_a: Align) {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::doc::Doc;

    #[test]
    fn delete_key_removes_selection_or_one_base() {
        let mut dv = DocView::new(Doc::new("t", b"ACGTACGT".to_vec()));
        dv.place_caret(3);
        delete_key(&mut dv, true);
        assert_eq!(dv.doc.seq, b"ACTACGT");
        delete_key(&mut dv, false);
        assert_eq!(dv.doc.seq, b"ACACGT");
        dv.select(1, 4);
        delete_key(&mut dv, true);
        assert_eq!(dv.doc.seq, b"AGT");
        dv.place_caret(0);
        delete_key(&mut dv, true);
        assert_eq!(dv.doc.seq, b"AGT");
    }

    #[test]
    fn arrow_heads_sit_at_the_end_of_the_feature() {
        let mut f = Feature::new("x", "gene", 10, 20, false);
        assert_eq!(head_of(&f, 100), (0, true));
        f.reverse = true;
        assert_eq!(head_of(&f, 100), (0, false));
        let mut w = Feature::new("w", "gene", 90, 100, false);
        w.spans = vec![(0, 10), (90, 100)];
        assert_eq!(head_of(&w, 100), (0, true));
        w.reverse = true;
        assert_eq!(head_of(&w, 100), (1, false));
    }

    #[test]
    fn clipping_maps_ranges_to_columns() {
        assert_eq!(clip(5, 30, 10, 20), Some((0, 10)));
        assert_eq!(clip(12, 15, 10, 20), Some((2, 5)));
        assert_eq!(clip(0, 10, 10, 20), None);
    }

    #[test]
    fn a_row_holds_as_many_characters_as_fit_with_no_upper_limit() {
        // numbers take 56 points and the sides 28: (1000 - 56 - 28) / 8 = 114 characters, more than the 75 it starts with
        assert_eq!(row_columns(1000.0, 56.0, 8.0), 114);
        assert_eq!(row_columns(2560.0, 56.0, 8.0), 309);
        // a window too narrow to read keeps the fewest
        assert_eq!(row_columns(100.0, 56.0, 8.0), MIN_COLUMNS);
        assert_eq!(row_columns(10.0, 56.0, 8.0), MIN_COLUMNS);
        // a view with no size yet starts at the original's 75
        assert_eq!(row_columns(f32::INFINITY, 56.0, 8.0), TEXT_COLUMNS);
        assert_eq!(row_columns(f32::NAN, 56.0, 8.0), TEXT_COLUMNS);
        assert_eq!(row_columns(1000.0, 56.0, 0.0), TEXT_COLUMNS);
    }

    #[test]
    fn the_ruler_runs_to_the_end_of_the_row_and_may_cut_its_last_ten_short() {
        // the original adds tens until the string is as long as the row
        assert_eq!(ruler_tens(75), 8);
        assert_eq!(ruler_tens(70), 7);
        assert_eq!(ruler_tens(71), 8);
        assert_eq!(ruler_tens(10), 1);
        assert_eq!(ruler_tens(0), 0);
    }
}
