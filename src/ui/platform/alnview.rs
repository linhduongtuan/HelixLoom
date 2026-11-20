//! The alignment viewer: names kept at the left and a ruler, conservation bars and the consensus kept at the
//! top while the letters scroll, with only the visible part painted so a long alignment stays quick.

use super::alignment::Mark;
use crate::seq::alnstats::Column;
use crate::ui::design::{self, Look};
use eframe::egui::{
    pos2, vec2, Align2, Color32, FontId, Rect, ScrollArea, Sense, Stroke, StrokeKind, Ui,
};

const FONT_SIZE: f32 = 13.0;
const ROW_HEIGHT: f32 = 20.0;
const NAME_WIDTH: f32 = 170.0;
const RULER_HEIGHT: f32 = 18.0;
const BARS_HEIGHT: f32 = 30.0;
const CONSENSUS_HEIGHT: f32 = 22.0;
const HEADER_HEIGHT: f32 = RULER_HEIGHT + BARS_HEIGHT + CONSENSUS_HEIGHT + 6.0;
/// A tick and a number on the ruler for every this many columns.
const RULER_STEP: usize = 10;
/// How heavy the line round the selected column is.
const SELECTED_LINE: f32 = 1.5;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Colouring {
    /// The colours of nucleotides or of amino-acid groups.
    ByType,
    /// Only the letters that differ from the consensus are coloured.
    Differences,
    Plain,
}

pub struct View<'a> {
    pub rows: &'a [(String, Vec<u8>)],
    pub columns: &'a [Column],
    /// The score that each column's bar shows, 0 to 1 (none: the column has no score and gets a mark instead).
    pub bars: &'a [Option<f32>],
    /// The bars of this score and above are drawn in full colour.
    pub highlight: f32,
    pub protein: bool,
    pub colouring: Colouring,
    pub look: &'a Look,
    /// Where the view is scrolled to as it begins: for a figure of what the window shows. None in the window itself.
    pub scroll: Option<eframe::egui::Vec2>,
    /// The column that is selected (counted from 0), which is outlined.
    pub selected: Option<usize>,
    /// Scroll so that the selected column is in the middle: it was selected somewhere else, and may be out of view.
    pub reveal: bool,
    /// Stretches of columns to mark (the oligos of a design).
    pub marks: &'a [Mark],
}

fn scroll_id() -> eframe::egui::Id {
    eframe::egui::Id::new("alignment viewer scroll")
}

/// The room that the first `rows` sequences and `columns` columns take, with a little to spare: a letter is counted a shade
/// wider than it is drawn, so that the last column is never cut.
pub fn room(rows: usize, columns: usize) -> eframe::egui::Vec2 {
    vec2(
        NAME_WIDTH + columns as f32 * FONT_SIZE * 0.62 + 24.0,
        HEADER_HEIGHT + rows as f32 * ROW_HEIGHT + 12.0,
    )
}

/// Where the viewer in the window was last scrolled to.
pub fn scrolled(ctx: &eframe::egui::Context) -> eframe::egui::Vec2 {
    ctx.data(|d| d.get_temp(scroll_id())).unwrap_or_default()
}

impl View<'_> {
    /// The background of the cell holding `c` in column `k`, if it has one.
    fn fill(&self, c: u8, k: usize) -> Option<Color32> {
        match self.colouring {
            Colouring::Plain => None,
            Colouring::ByType if self.protein => design::residue_fill(c, self.look),
            Colouring::ByType => design::base_fill(c, self.look),
            Colouring::Differences => {
                let same = c.eq_ignore_ascii_case(&self.columns[k].consensus);
                (!same && !matches!(c, b'-' | b'.')).then_some(self.look.accent_soft)
            }
        }
    }
}

/// How many letters are laid out to measure how wide one is.
const MEASURED_RUN: usize = 64;

/// How wide a letter is in the viewer's monospaced font, as the letters of a row are laid out: each a whole number of pixels from the
/// last, which is not always the number of points that one letter measures (at 80 per cent, or 110, or where the screen has 1.1 pixels
/// to the point). Cells as wide as one letter rounded up would drift away from the letters of a long row, so the width of a long run
/// is shared out instead.
fn letter_width(ui: &Ui, font: &FontId) -> f32 {
    let run = ui
        .painter()
        .layout_no_wrap("M".repeat(MEASURED_RUN), font.clone(), Color32::WHITE)
        .size()
        .x;
    run / MEASURED_RUN as f32
}

/// The range of items of `size` points whose positions (counted from `offset`) fall in `lo..hi`.
fn visible(lo: f32, hi: f32, offset: f32, size: f32, count: usize) -> std::ops::Range<usize> {
    let first = ((lo - offset) / size).floor().max(0.0) as usize;
    let last = (((hi - offset) / size).ceil().max(0.0) as usize + 1).min(count);
    first.min(last)..last
}

/// How far to the right a view that is `room` wide is scrolled to have column `k` in the middle of its letters (the names
/// take the left of the view), for columns `cell` wide in a content `content` wide: never past either end.
fn offset_to_show(k: usize, cell: f32, room: f32, content: f32) -> f32 {
    let middle = (k as f32 + 0.5) * cell - (room - NAME_WIDTH) / 2.0;
    middle.clamp(0.0, (content - room).max(0.0))
}

/// The column under a pointer at `x` (counted from the left of the content) in a view scrolled `left` to the right, if it is
/// over the letters of one: not over the names, which stay at the left, and not past the last column.
fn column_at(x: f32, left: f32, cell: f32, width: usize) -> Option<usize> {
    if x < left + NAME_WIDTH {
        return None;
    }
    let k = ((x - NAME_WIDTH) / cell).floor();
    (k >= 0.0 && (k as usize) < width).then_some(k as usize)
}

/// Whether a pointer at `y` (counted from the top of the content) in a view scrolled `top` down is over the columns: on the
/// ruler, the bars and the consensus, which stay at the top, or on one of the `rows`; not in the spare room under the last.
fn over_columns(y: f32, top: f32, rows: usize) -> bool {
    y < top + HEADER_HEIGHT || y < HEADER_HEIGHT + rows as f32 * ROW_HEIGHT
}

/// Draw the alignment; the column that was clicked in this frame, if one was.
pub fn show(ui: &mut Ui, v: &View) -> Option<usize> {
    let font = FontId::monospace(FONT_SIZE);
    let cell = letter_width(ui, &font);
    let (n, width) = (v.rows.len(), v.columns.len());
    let content = vec2(
        NAME_WIDTH + width as f32 * cell + 24.0,
        HEADER_HEIGHT + n as f32 * ROW_HEIGHT + 12.0,
    );
    let selected = v.selected.filter(|&k| k < width);
    let mut area = ScrollArea::both().auto_shrink([false, false]);
    if let Some(at) = v.scroll {
        area = area.scroll_offset(at);
    } else if let Some(k) = selected.filter(|_| v.reveal) {
        // sideways only: the rows stay where they were
        let x = offset_to_show(k, cell, ui.available_width(), content.x);
        area = area.scroll_offset(vec2(x, scrolled(ui.ctx()).y));
    }
    let mut clicked = None;
    let shown = area.show_viewport(ui, |ui, view| {
        let (rect, response) = ui.allocate_exact_size(content, Sense::click());
        let origin = rect.min;
        if let Some(pointer) = response
            .interact_pointer_pos()
            .filter(|_| response.clicked())
        {
            let at = pointer - origin;
            if over_columns(at.y, view.min.y, n) {
                clicked = column_at(at.x, view.min.x, cell, width);
            }
        }
        let p = ui.painter_at(rect);
        let look = v.look;
        let at = |x: f32, y: f32| origin + vec2(x, y);
        let cols = visible(view.min.x, view.max.x, 0.0, cell, width);
        let rows = visible(view.min.y, view.max.y, HEADER_HEIGHT, ROW_HEIGHT, n);

        // the letters: one run of colour and one piece of text for each visible stretch of a row
        for r in rows.clone() {
            let y = HEADER_HEIGHT + r as f32 * ROW_HEIGHT;
            let bases = &v.rows[r].1;
            let end = cols.end.min(bases.len());
            let start = cols.start.min(end);
            let fills: Vec<Option<Color32>> = (start..end).map(|k| v.fill(bases[k], k)).collect();
            let mut run: Option<(usize, Color32)> = None;
            for k in start..=end {
                // k == end is one past the last cell: it closes the run that is still open
                let fill = fills.get(k - start).copied().flatten();
                if let Some((from, colour)) = run {
                    if fill != Some(colour) {
                        let rect = Rect::from_min_max(
                            at(NAME_WIDTH + from as f32 * cell, y),
                            at(NAME_WIDTH + k as f32 * cell, y + ROW_HEIGHT),
                        );
                        p.rect_filled(rect, 0.0, colour);
                        run = None;
                    }
                }
                if run.is_none() {
                    run = fill.map(|colour| (k, colour));
                }
            }
            // the letters are in the text colour, as they always were; the colour-blind-safe cells (a light yellow on a dark page) and
            // high contrast need black or white by the cell. A stretch of one colour is one piece of text.
            let inks: Vec<Color32> = fills
                .iter()
                .map(|fill| design::letter_ink(look, *fill))
                .collect();
            let mut from = start;
            while from < end {
                let ink = inks[from - start];
                let to = (from..end).find(|&k| inks[k - start] != ink).unwrap_or(end);
                p.text(
                    at(NAME_WIDTH + from as f32 * cell, y + ROW_HEIGHT / 2.0),
                    Align2::LEFT_CENTER,
                    String::from_utf8_lossy(&bases[from..to]).into_owned(),
                    font.clone(),
                    ink,
                );
                from = to;
            }
        }

        // the selected column, outlined over all the rows; the names and the top are painted over it where it is under them
        if let Some(k) = selected {
            let x = NAME_WIDTH + k as f32 * cell;
            p.rect_stroke(
                Rect::from_min_max(
                    at(x, HEADER_HEIGHT),
                    at(x + cell, HEADER_HEIGHT + n as f32 * ROW_HEIGHT),
                ),
                0.0,
                Stroke::new(SELECTED_LINE, look.accent),
                StrokeKind::Inside,
            );
        }

        // the stretches marked for oligos, outlined over all the rows
        for m in v.marks.iter().filter(|m| m.last < width) {
            let x1 = NAME_WIDTH + m.first as f32 * cell;
            let x2 = NAME_WIDTH + (m.last + 1) as f32 * cell;
            p.rect_stroke(
                Rect::from_min_max(
                    at(x1, HEADER_HEIGHT),
                    at(x2, HEADER_HEIGHT + n as f32 * ROW_HEIGHT),
                ),
                0.0,
                Stroke::new(SELECTED_LINE, m.role.colour(look)),
                StrokeKind::Inside,
            );
        }

        // the names stay at the left
        let left = view.min.x;
        let names_area =
            Rect::from_min_max(at(left, HEADER_HEIGHT), at(left + NAME_WIDTH, content.y));
        p.rect_filled(names_area, 0.0, look.surface);
        p.line_segment(
            [
                at(left + NAME_WIDTH, HEADER_HEIGHT),
                at(left + NAME_WIDTH, content.y),
            ],
            Stroke::new(1.0, look.border),
        );
        for r in rows {
            let y = HEADER_HEIGHT + r as f32 * ROW_HEIGHT + ROW_HEIGHT / 2.0;
            let name = shorten(&v.rows[r].0, 22);
            p.text(
                at(left + 8.0, y),
                Align2::LEFT_CENTER,
                name,
                FontId::proportional(12.5),
                look.text,
            );
        }

        // the ruler, the conservation bars and the consensus stay at the top
        let top = view.min.y;
        p.rect_filled(
            Rect::from_min_max(at(left, top), at(view.max.x, top + HEADER_HEIGHT)),
            0.0,
            look.surface,
        );
        p.line_segment(
            [
                at(left, top + HEADER_HEIGHT),
                at(view.max.x, top + HEADER_HEIGHT),
            ],
            Stroke::new(1.0, look.border),
        );
        for k in cols {
            let x = NAME_WIDTH + k as f32 * cell;
            if selected == Some(k) {
                // its number, bar and consensus letter stand on a tint, so that the column is found from the top too
                p.rect_filled(
                    Rect::from_min_max(at(x, top), at(x + cell, top + HEADER_HEIGHT)),
                    0.0,
                    look.accent_soft,
                );
            }
            if (k + 1) % RULER_STEP == 0 {
                p.line_segment(
                    [
                        at(x + cell, top + RULER_HEIGHT - 5.0),
                        at(x + cell, top + RULER_HEIGHT - 1.0),
                    ],
                    Stroke::new(1.0, look.dim),
                );
                p.text(
                    at(x + cell, top + RULER_HEIGHT / 2.0 - 3.0),
                    Align2::CENTER_CENTER,
                    (k + 1).to_string(),
                    FontId::proportional(10.0),
                    look.dim,
                );
            }
            let column = v.columns[k];
            let base = top + RULER_HEIGHT + BARS_HEIGHT - 2.0;
            match v.bars.get(k).copied().flatten() {
                Some(score) => {
                    let bar = score * (BARS_HEIGHT - 4.0);
                    let tone = if score >= v.highlight {
                        look.accent
                    } else {
                        look.accent.gamma_multiply(0.45)
                    };
                    p.rect_filled(
                        Rect::from_min_max(at(x + 1.0, base - bar), at(x + cell - 1.0, base)),
                        1.0,
                        tone,
                    );
                }
                // not scored (too many gaps): a thin mark, so that it is not mistaken for a score of nothing
                None => {
                    p.rect_filled(
                        Rect::from_min_max(at(x + 3.0, base - 1.5), at(x + cell - 3.0, base)),
                        0.0,
                        look.dim.gamma_multiply(0.5),
                    );
                }
            }
            p.text(
                pos2(
                    origin.x + x + cell / 2.0,
                    origin.y + top + RULER_HEIGHT + BARS_HEIGHT + CONSENSUS_HEIGHT / 2.0,
                ),
                Align2::CENTER_CENTER,
                char::from(column.consensus).to_string(),
                font.clone(),
                look.dim,
            );
        }
        // a bar for each marked stretch in the spare room at the foot of the header
        for m in v.marks.iter().filter(|m| m.last < width) {
            let x1 = NAME_WIDTH + m.first as f32 * cell;
            let x2 = NAME_WIDTH + (m.last + 1) as f32 * cell;
            p.rect_filled(
                Rect::from_min_max(
                    at(x1, top + HEADER_HEIGHT - 5.0),
                    at(x2, top + HEADER_HEIGHT - 1.0),
                ),
                1.0,
                m.role.colour(look),
            );
        }
        // the corner above the names
        p.rect_filled(
            Rect::from_min_max(at(left, top), at(left + NAME_WIDTH, top + HEADER_HEIGHT)),
            0.0,
            look.surface,
        );
        p.text(
            at(left + 8.0, top + RULER_HEIGHT + BARS_HEIGHT * 0.5),
            Align2::LEFT_CENTER,
            "Conservation",
            FontId::proportional(11.5),
            look.dim,
        );
        p.text(
            at(
                left + 8.0,
                top + RULER_HEIGHT + BARS_HEIGHT + CONSENSUS_HEIGHT / 2.0,
            ),
            Align2::LEFT_CENTER,
            "Consensus",
            FontId::proportional(11.5),
            look.dim,
        );
        p.line_segment(
            [
                at(left + NAME_WIDTH, top),
                at(left + NAME_WIDTH, top + HEADER_HEIGHT),
            ],
            Stroke::new(1.0, look.border),
        );
    });
    if v.scroll.is_none() {
        ui.ctx()
            .data_mut(|d| d.insert_temp(scroll_id(), shown.state.offset));
    }
    clicked
}

/// `name` cut to `limit` characters, with an ellipsis when something was cut.
fn shorten(name: &str, limit: usize) -> String {
    if name.chars().count() <= limit {
        return name.to_string();
    }
    let kept: String = name.chars().take(limit - 1).collect();
    format!("{kept}…")
}

#[cfg(test)]
mod tests {
    use super::*;
    use eframe::egui;

    #[test]
    fn only_the_visible_items_are_painted_with_one_spare() {
        // 10-point items from 0: points 25 to 61 show items 2 to 6, and item 7 is the spare
        assert_eq!(visible(25.0, 61.0, 0.0, 10.0, 100), 2..8);
        // the viewport can end beyond the last item, and start before the first
        assert_eq!(visible(-30.0, 25.0, 0.0, 10.0, 4), 0..4);
        assert_eq!(visible(-30.0, 25.0, 0.0, 10.0, 100), 0..4);
        assert_eq!(visible(500.0, 900.0, 0.0, 10.0, 4), 4..4);
        // an offset (a header) shifts them
        assert_eq!(visible(100.0, 130.0, 80.0, 10.0, 100), 2..6);
        // above the first item only the spare one is painted
        assert_eq!(visible(0.0, 50.0, 80.0, 10.0, 100), 0..1);
    }

    /// How far a run of `n` letters, laid out as one piece of text, reaches in a window that is zoomed by `zoom` on a screen of `native`
    /// pixels to the point, and how far `n` cells of the viewer's width reach.
    fn run_and_cells(native: f32, zoom: f32, n: usize) -> (f32, f32) {
        let ctx = egui::Context::default();
        crate::ui::fonts::install(&ctx);
        ctx.set_zoom_factor(zoom);
        let font = FontId::monospace(FONT_SIZE);
        let mut out = (0.0, 0.0);
        for _ in 0..3 {
            let input = egui::RawInput {
                viewports: std::iter::once((
                    egui::ViewportId::ROOT,
                    egui::ViewportInfo {
                        native_pixels_per_point: Some(native),
                        ..Default::default()
                    },
                ))
                .collect(),
                screen_rect: Some(Rect::from_min_size(pos2(0.0, 0.0), vec2(800.0, 600.0))),
                ..Default::default()
            };
            let mut output = ctx.run_ui(input, |ui| {
                let run = ui
                    .painter()
                    .layout_no_wrap("M".repeat(n), font.clone(), Color32::WHITE)
                    .size()
                    .x;
                out = (run, letter_width(ui, &font) * n as f32);
            });
            output.textures_delta.clear();
        }
        out
    }

    #[test]
    fn the_cells_are_as_wide_as_the_letters_at_every_size_the_window_can_be_zoomed_to() {
        // a letter is a whole number of pixels wide, which at some sizes is a little less than the points it was measured in: a cell
        // as wide as the rounded-up width of one letter would then drift away from the letters of a long row
        for native in [1.0, 1.25, 1.5, 2.0] {
            for zoom in [0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 1.75, 2.0] {
                let (run, cells) = run_and_cells(native, zoom, 200);
                assert!(
                    (run - cells).abs() < 0.5,
                    "{native} pixels to the point, zoom {zoom}: 200 letters reach {run}, 200 cells reach {cells}"
                );
            }
        }
    }

    #[test]
    fn the_pointer_is_over_a_column_only_where_its_letters_are() {
        // columns 8 points wide after the names
        let x = |k: f32| NAME_WIDTH + k * 8.0;
        assert_eq!(column_at(x(0.0), 0.0, 8.0, 40), Some(0));
        assert_eq!(column_at(x(0.99), 0.0, 8.0, 40), Some(0));
        assert_eq!(column_at(x(1.0), 0.0, 8.0, 40), Some(1));
        assert_eq!(column_at(x(39.5), 0.0, 8.0, 40), Some(39));
        // past the last column, and over the names
        assert_eq!(column_at(x(40.0), 0.0, 8.0, 40), None);
        assert_eq!(column_at(NAME_WIDTH - 0.5, 0.0, 8.0, 40), None);
        assert_eq!(column_at(3.0, 0.0, 8.0, 40), None);
        // scrolled 100 points to the right, the names lie over the first columns: a click there is on a name
        assert_eq!(column_at(x(5.0), 100.0, 8.0, 40), None);
        assert_eq!(column_at(100.0 + NAME_WIDTH - 0.5, 100.0, 8.0, 40), None);
        assert_eq!(column_at(100.0 + NAME_WIDTH, 100.0, 8.0, 40), Some(12));
        assert_eq!(column_at(x(0.0), 0.0, 8.0, 0), None, "no columns");
    }

    #[test]
    fn the_pointer_is_over_the_columns_on_the_rows_and_on_the_top_that_stays_but_not_under_the_last_row(
    ) {
        let rows_end = HEADER_HEIGHT + 3.0 * ROW_HEIGHT;
        assert!(over_columns(0.0, 0.0, 3));
        assert!(over_columns(HEADER_HEIGHT + 1.0, 0.0, 3));
        assert!(over_columns(rows_end - 0.5, 0.0, 3));
        assert!(!over_columns(rows_end, 0.0, 3));
        assert!(!over_columns(rows_end + 11.0, 0.0, 3));
        // scrolled down, the top lies over the last rows and the room under them: it is the top of the columns still
        assert!(over_columns(rows_end + 5.0, rows_end - 10.0, 3));
        assert!(!over_columns(
            rows_end + 5.0 + HEADER_HEIGHT,
            rows_end - 10.0,
            3
        ));
        // no rows: only the top
        assert!(over_columns(HEADER_HEIGHT - 1.0, 0.0, 0));
        assert!(!over_columns(HEADER_HEIGHT, 0.0, 0));
    }

    #[test]
    fn a_column_is_scrolled_to_the_middle_of_the_letters_and_never_past_an_end() {
        // 400 columns of 8 points in a view 970 wide, of which the names take 170: the letters have 800, a hundred columns
        let content = NAME_WIDTH + 400.0 * 8.0 + 24.0;
        let shown = |k: usize| offset_to_show(k, 8.0, 970.0, content);
        // column 200 in the middle: fifty columns on each side of it
        assert_eq!(shown(200), 200.5 * 8.0 - 400.0);
        // near the start there is nothing to scroll back to, and near the end the content ends
        assert_eq!(shown(0), 0.0);
        assert_eq!(shown(30), 0.0);
        assert_eq!(shown(399), content - 970.0);
        // an alignment that fits is not scrolled at all
        assert_eq!(offset_to_show(5, 8.0, 970.0, 500.0), 0.0);
    }

    #[test]
    fn long_names_are_cut_with_an_ellipsis() {
        assert_eq!(shorten("Isolate 1", 22), "Isolate 1");
        assert_eq!(
            shorten("a_very_long_sequence_name_indeed", 10),
            "a_very_lo…"
        );
    }
}
