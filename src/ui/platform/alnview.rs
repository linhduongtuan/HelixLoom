//! The alignment viewer: names kept at the left and a ruler, conservation bars and the consensus kept at the
//! top while the letters scroll, with only the visible part painted so a long alignment stays quick.

use crate::seq::alnstats::Column;
use crate::ui::design::{self, Look};
use eframe::egui::{pos2, vec2, Align2, Color32, FontId, Rect, ScrollArea, Sense, Stroke, Ui};

const FONT_SIZE: f32 = 13.0;
const ROW_HEIGHT: f32 = 20.0;
const NAME_WIDTH: f32 = 170.0;
const RULER_HEIGHT: f32 = 18.0;
const BARS_HEIGHT: f32 = 30.0;
const CONSENSUS_HEIGHT: f32 = 22.0;
const HEADER_HEIGHT: f32 = RULER_HEIGHT + BARS_HEIGHT + CONSENSUS_HEIGHT + 6.0;
/// A tick and a number on the ruler for every this many columns.
const RULER_STEP: usize = 10;

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
    pub protein: bool,
    pub colouring: Colouring,
    pub look: &'a Look,
}

impl View<'_> {
    /// The background of the cell holding `c` in column `k`, if it has one.
    fn fill(&self, c: u8, k: usize) -> Option<Color32> {
        match self.colouring {
            Colouring::Plain => None,
            Colouring::ByType if self.protein => design::residue_fill(c, self.look.dark),
            Colouring::ByType => design::base_fill(c, self.look.dark),
            Colouring::Differences => {
                let same = c.eq_ignore_ascii_case(&self.columns[k].consensus);
                (!same && !matches!(c, b'-' | b'.')).then_some(self.look.accent_soft)
            }
        }
    }
}

/// How wide a letter is in the viewer's monospaced font.
fn letter_width(ui: &Ui, font: &FontId) -> f32 {
    ui.painter()
        .layout_no_wrap("M".to_string(), font.clone(), Color32::WHITE)
        .size()
        .x
        .ceil()
}

/// The range of items of `size` points whose positions (counted from `offset`) fall in `lo..hi`.
fn visible(lo: f32, hi: f32, offset: f32, size: f32, count: usize) -> std::ops::Range<usize> {
    let first = ((lo - offset) / size).floor().max(0.0) as usize;
    let last = (((hi - offset) / size).ceil().max(0.0) as usize + 1).min(count);
    first.min(last)..last
}

pub fn show(ui: &mut Ui, v: &View) {
    let font = FontId::monospace(FONT_SIZE);
    let cell = letter_width(ui, &font);
    let (n, width) = (v.rows.len(), v.columns.len());
    let content = vec2(
        NAME_WIDTH + width as f32 * cell + 24.0,
        HEADER_HEIGHT + n as f32 * ROW_HEIGHT + 12.0,
    );
    ScrollArea::both()
        .auto_shrink([false, false])
        .show_viewport(ui, |ui, view| {
            let (rect, _) = ui.allocate_exact_size(content, Sense::hover());
            let origin = rect.min;
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
                let mut run: Option<(usize, Color32)> = None;
                for k in start..=end {
                    // k == end is one past the last cell: it closes the run that is still open
                    let fill = bases
                        .get(k)
                        .filter(|_| k < end)
                        .and_then(|&base| v.fill(base, k));
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
                let text = String::from_utf8_lossy(&bases[start..end]).into_owned();
                p.text(
                    at(NAME_WIDTH + start as f32 * cell, y + ROW_HEIGHT / 2.0),
                    Align2::LEFT_CENTER,
                    text,
                    font.clone(),
                    look.text,
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
                let bar = column.conservation * (BARS_HEIGHT - 4.0);
                let base = top + RULER_HEIGHT + BARS_HEIGHT - 2.0;
                let tone = if column.conservation >= 0.999 {
                    look.accent
                } else {
                    look.accent.gamma_multiply(0.45)
                };
                p.rect_filled(
                    Rect::from_min_max(at(x + 1.0, base - bar), at(x + cell - 1.0, base)),
                    1.0,
                    tone,
                );
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

    #[test]
    fn long_names_are_cut_with_an_ellipsis() {
        assert_eq!(shorten("Isolate 1", 22), "Isolate 1");
        assert_eq!(
            shorten("a_very_long_sequence_name_indeed", 10),
            "a_very_lo…"
        );
    }
}
