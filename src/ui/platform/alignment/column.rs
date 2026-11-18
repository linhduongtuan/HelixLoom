//! The selected column of the alignment: what the page says of it, and the line it says it on (`docs/ROADMAP.md` 1.4). A
//! click on a letter, or on the ruler above it, selects the column; a click on it again lets it go.

use super::AlignmentPage;
use crate::seq::alnstats;
use crate::ui::design::{self, Look};
use eframe::egui::Ui;

/// What the line says while no column is selected.
pub const HINT: &str = "Click a column to see what it holds.";
/// The link that lets the column go.
pub const CLEAR: &str = "Deselect";
/// So many of a column's letters are named; the rest are counted.
const MOST_SYMBOLS: usize = 6;

/// What the page says of one column.
#[derive(Clone, Debug, PartialEq)]
pub struct Facts {
    /// The column, counted from 0, of `of` columns.
    pub column: usize,
    pub of: usize,
    /// The commonest letter, or a gap where gaps are half of the rows or more.
    pub consensus: u8,
    /// The letters of the column with the rows that carry each, the commonest first.
    pub symbols: Vec<(u8, usize)>,
    pub gaps: usize,
    /// The score that the bars show, by its name, and the column's (none: it is not scored).
    pub score_name: &'static str,
    pub score: Option<f32>,
    /// Whether the bar of the column is drawn in full colour: it is among the conserved.
    pub conserved: bool,
    /// The conserved block the column is in, first and last column (counted from 0).
    pub block: Option<(usize, usize)>,
}

impl Facts {
    /// The facts in a line, as the page writes them.
    pub fn words(&self) -> String {
        let consensus = match self.consensus {
            b'-' => "a gap".to_string(),
            letter => char::from(letter).to_string(),
        };
        let named: Vec<String> = self
            .symbols
            .iter()
            .take(MOST_SYMBOLS)
            .map(|&(letter, rows)| format!("{} {rows}", char::from(letter)))
            .collect();
        let more = match self.symbols.len().saturating_sub(MOST_SYMBOLS) {
            0 => String::new(),
            more => format!(" and {more} more"),
        };
        let gaps = match self.gaps {
            0 => "no gaps".to_string(),
            1 => "1 gap".to_string(),
            gaps => format!("{gaps} gaps"),
        };
        let score = match self.score {
            None => "not scored (too many gaps)".to_string(),
            Some(score) if self.conserved => {
                format!("{} {score:.2}, among the conserved", self.score_name)
            }
            Some(score) => format!("{} {score:.2}", self.score_name),
        };
        let mut parts = vec![
            format!("Column {} of {}", self.column + 1, self.of),
            format!("consensus {consensus}"),
        ];
        // a column that is a gap in every row has no letters to name
        if !named.is_empty() {
            parts.push(format!("{}{more}", named.join(", ")));
        }
        parts.extend([gaps, score]);
        if let Some((from, to)) = self.block {
            parts.push(format!("in the conserved block {}-{}", from + 1, to + 1));
        }
        parts.join(" · ")
    }
}

impl AlignmentPage {
    /// The column that is selected, counted from 0.
    pub fn selected_column(&self) -> Option<usize> {
        self.selected
    }

    /// Select a column (none, or one that the alignment does not have: no column). `reveal` asks the viewer to scroll to it
    /// when it is next drawn: for a column that was selected somewhere else than in the viewer.
    pub fn select_column(&mut self, column: Option<usize>, reveal: bool) {
        let width = self.stats.as_ref().map_or(0, |s| s.columns.len());
        self.selected = column.filter(|&k| k < width);
        self.reveal = reveal && self.selected.is_some();
    }

    /// What the page says of column `k`; none without an alignment, or past its last column.
    pub(in crate::ui) fn column_facts(&self, k: usize) -> Option<Facts> {
        let (result, stats) = (self.result.as_ref()?, self.stats.as_ref()?);
        let column = stats.columns.get(k)?;
        let score = stats.bars.values.get(k).copied().flatten();
        Some(Facts {
            column: k,
            of: stats.columns.len(),
            consensus: column.consensus,
            // read where they are: this is asked in every frame while a column is selected
            symbols: alnstats::symbols_of(result.rows.iter().map(|row| &row.1), k),
            gaps: column.gaps,
            score_name: stats.bars.score.label(),
            score,
            conserved: score.is_some_and(|s| s >= stats.bars.highlight),
            block: stats
                .bars
                .blocks
                .iter()
                .find(|b| (b.from..=b.to).contains(&k))
                .map(|b| (b.from, b.to)),
        })
    }
}

/// What was clicked on the line of the selected column.
#[derive(Default)]
pub struct Asks {
    pub clear: bool,
    /// Go to the Structure page, where the residue of the column is picked.
    pub to_structure: bool,
}

/// The words of the link to the residues of the structure that are in the column (`residues`, each as the structure page
/// names it); none when no residue is.
pub fn link_words(residues: &[String]) -> Option<String> {
    let first = residues.first()?;
    let more = match residues.len() - 1 {
        0 => String::new(),
        1 => " (and 1 more chain)".to_string(),
        more => format!(" (and {more} more chains)"),
    };
    Some(format!("Show {first} in 3D{more}"))
}

/// The line above the viewer: what the selected column holds, or how to select one. `residues` are the residues of the
/// structure on the Structure page that were lined up with the column.
pub fn line(ui: &mut Ui, look: &Look, facts: Option<&Facts>, residues: &[String]) -> Asks {
    let words = facts.map(Facts::words);
    let to_structure = link_words(residues);
    let mut links: Vec<(&str, &str)> = Vec::new();
    if let Some(words) = &to_structure {
        links.push((
            words,
            "Go to the Structure page, where the residue is picked",
        ));
    }
    links.push((CLEAR, "Or click the column again"));
    let clicked = design::selection_line(ui, look, HINT, words.as_deref(), &links);
    Asks {
        clear: clicked == Some(links.len() - 1),
        to_structure: to_structure.is_some() && clicked == Some(0),
    }
}

#[cfg(test)]
mod tests;
