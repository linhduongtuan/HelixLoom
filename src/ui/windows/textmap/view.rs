//! The text of a Text Map in its window: a line to a row, with the colours of the features behind the characters they cover, and
//! the lines of DNA answering a click.

use super::MapRun;
use crate::seq::features::Rgb;
use crate::seq::textmap::{TextMap, MARGIN};
use crate::ui::design::{self, Look};
use eframe::egui::{
    self, text::LayoutJob, Color32, FontId, RichText, ScrollArea, Sense, TextFormat, TextStyle,
    TextWrapMode, Ui,
};

/// A coloured stretch of a line: from, up to but not including, and the colour behind it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Piece {
    pub from: usize,
    pub to: usize,
    pub color: Rgb,
}

/// The lines of a page as characters, and where each has colour behind it.
pub(super) struct Painted {
    lines: Vec<Vec<char>>,
    pieces: Vec<Vec<Piece>>,
}

impl Painted {
    /// The colours: each tag with a background, the one drawn over the others last, put on the characters it covers.
    pub(super) fn of(page: &TextMap) -> Painted {
        let lines: Vec<Vec<char>> = page.text.split('\n').map(|l| l.chars().collect()).collect();
        let mut colors: Vec<Vec<Option<Rgb>>> = lines.iter().map(|l| vec![None; l.len()]).collect();
        for tag in page.tags.iter().filter(|t| t.background.is_some()) {
            for &(from, to) in &tag.ranges {
                for line in from.0..=to.0.min(lines.len()) {
                    let Some(row) = colors.get_mut(line - 1) else {
                        continue;
                    };
                    let start = if line == from.0 { from.1 } else { 0 };
                    let end = if line == to.0 { to.1 } else { row.len() };
                    for cell in row.iter_mut().take(end).skip(start) {
                        *cell = tag.background;
                    }
                }
            }
        }
        let pieces = colors.iter().map(|row| runs(row)).collect();
        Painted { lines, pieces }
    }

    /// The coloured stretches of each line (for the tests).
    #[cfg(test)]
    pub(super) fn pieces_in_rows(&self) -> &[Vec<Piece>] {
        &self.pieces
    }

    /// The lines of the text and the characters of the longest.
    pub(super) fn shape(&self) -> (usize, usize) {
        let longest = self.lines.iter().map(Vec::len).max();
        (self.lines.len(), longest.unwrap_or(0))
    }

    /// A line as text with its colours. An empty line still takes a row.
    fn job(&self, n: usize, font: &FontId, look: &Look) -> LayoutJob {
        let chars = &self.lines[n];
        let mut job = LayoutJob::default();
        let mut put = |text: String, color: Option<Rgb>| {
            let format = TextFormat {
                font_id: font.clone(),
                color: color.map_or(look.text, |_| Color32::BLACK),
                background: color.map_or(Color32::TRANSPARENT, |c| {
                    Color32::from_rgb(c[0], c[1], c[2])
                }),
                ..TextFormat::default()
            };
            job.append(&text, 0.0, format);
        };
        let mut at = 0;
        for piece in &self.pieces[n] {
            if piece.from > at {
                put(chars[at..piece.from].iter().collect(), None);
            }
            put(
                chars[piece.from..piece.to].iter().collect(),
                Some(piece.color),
            );
            at = piece.to;
        }
        let rest: String = chars[at..].iter().collect();
        put(
            if rest.is_empty() && at == 0 {
                " ".into()
            } else {
                rest
            },
            None,
        );
        job
    }
}

/// The stretches of a row of cells that have a colour: the neighbours of the same colour are one.
fn runs(row: &[Option<Rgb>]) -> Vec<Piece> {
    let mut out: Vec<Piece> = Vec::new();
    for (at, cell) in row.iter().enumerate() {
        let Some(color) = *cell else {
            continue;
        };
        match out.last_mut() {
            Some(last) if last.to == at && last.color == color => last.to = at + 1,
            _ => out.push(Piece {
                from: at,
                to: at + 1,
                color,
            }),
        }
    }
    out
}

/// A click on the text: the line (from 1), the character on it (from 0), and whether Shift was held.
pub(super) struct Click {
    pub line: usize,
    pub col: usize,
    pub shift: bool,
}

/// The buttons and the note that the sequence has changed; (make again, copy as text).
pub(super) fn head(ui: &mut Ui, look: &Look, live: bool) -> (bool, bool) {
    let mut pressed = (false, false);
    ui.horizontal(|ui| {
        pressed.0 = design::primary_button(ui, look, "Make again").clicked();
        pressed.1 = design::ghost_button(ui, look, "Copy as text").clicked();
    });
    if !live {
        let text =
            "The sequence has changed since this page was made; make it again to refresh it.";
        ui.label(RichText::new(text).small().color(look.warn));
    }
    ui.separator();
    pressed
}

/// The rows, scrolling both ways: the lines are not wrapped, as the original's are not.
pub(super) fn rows(ui: &mut Ui, look: &Look, run: &MapRun) -> Option<Click> {
    let font = TextStyle::Monospace.resolve(ui.style());
    let advance = ui.ctx().fonts_mut(|f| f.glyph_width(&font, 'M'));
    let row = ui.text_style_height(&TextStyle::Monospace);
    let mut click = None;
    ScrollArea::both().auto_shrink([false, false]).show_rows(
        ui,
        row,
        run.painted.lines.len(),
        |ui, range| {
            for n in range {
                let job = run.painted.job(n, &font, look);
                let label = egui::Label::new(job)
                    .wrap_mode(TextWrapMode::Extend)
                    .selectable(false)
                    .sense(Sense::click());
                let response = ui.add(label);
                let at = response
                    .interact_pointer_pos()
                    .filter(|_| response.clicked());
                if let Some(pos) = at {
                    let col = ((pos.x - response.rect.min.x) / advance).floor().max(0.0) as usize;
                    click = Some(Click {
                        line: n + 1,
                        col,
                        shift: ui.input(|i| i.modifiers.shift),
                    });
                }
            }
        },
    );
    click
}

/// The base of the sequence (from 0) that a click on a line of DNA is on; none for a click anywhere else.
pub(super) fn base_at(run: &MapRun, line: usize, col: usize) -> Option<usize> {
    let row = run.page.dna_rows.iter().find(|r| r.line == line)?;
    if col < MARGIN || col >= MARGIN + row.bases {
        return None;
    }
    let base = row.number - 1 + (col - MARGIN);
    Some(base.min(run.seq_len.saturating_sub(1)))
}
