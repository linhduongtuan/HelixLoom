//! One selection, shown on every page that holds it (`docs/ROADMAP.md` 1.4): a column of the alignment and the residue of
//! the structure that was lined up with it. The two are linked from the moment the structure is coloured by the alignment
//! (*Colour by conservation* lines each chain up with a row) until the alignment or its score changes, which takes the
//! colours and the columns off again.

use super::{Platform, Workspace};
use crate::pipeline::steps::Aligned;
use crate::ui::app::App;

/// The sequence of the example structure (ubiquitin, PDB 1UBQ).
const UBIQUITIN: &[u8] =
    b"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";
/// How many made-up variants stand beside it in the picture of the link, and the column that is selected there (from 0).
const VARIANTS: usize = 11;
const PICTURED_COLUMN: usize = 43;

/// Ubiquitin and variants of it that are made up for the picture (each has a few residues changed, by a rule and not from
/// any organism): the rows say so by their names.
fn pictured_family() -> Aligned {
    let variant = |v: usize| {
        let mut row = UBIQUITIN.to_vec();
        // the residues 5v + 2, 5v + 13, … of the variant, so that the columns differ in how many rows changed them
        for at in ((5 * v + 2)..row.len()).step_by(11 + v) {
            row[at] = if row[at] == b'A' { b'S' } else { b'A' };
        }
        (format!("made-up variant {}", v + 1), row)
    };
    let rows = std::iter::once(("ubiquitin (1UBQ)".to_string(), UBIQUITIN.to_vec()))
        .chain((0..VARIANTS).map(variant))
        .collect();
    Aligned {
        rows,
        engine: "made up for this picture".into(),
        seconds: 0.0,
    }
}

impl App {
    /// For a picture of the link: the example structure lined up with ubiquitin and made-up variants of it, one column
    /// selected, on the Structure page or the Alignment page.
    pub(in crate::ui) fn show_linked_example(&mut self, on_structure: bool) {
        self.platform.structure.load_example();
        self.platform
            .alignment
            .show_aligned(pictured_family(), "An alignment made up for this picture.");
        self.colour_structure_by_alignment();
        self.platform.select_column(Some(PICTURED_COLUMN));
        self.workspace = if on_structure {
            Workspace::Structure
        } else {
            Workspace::Alignment
        };
    }
}

impl Platform {
    /// Select a column of the alignment (none: let it go), as a click on the Alignment page does. A structure that is lined
    /// up with the alignment picks the residue of that column; one that is not is left as it is.
    pub fn select_column(&mut self, column: Option<usize>) {
        self.alignment.select_column(column, false);
        self.structure.pick_column(self.alignment.selected_column());
    }

    /// The residue that is picked on the Structure page has changed (another, or none). When the structure is lined up with
    /// the alignment, the column of that residue is selected there, and the Alignment page scrolls to it when it is next
    /// shown; a residue that no row has selects no column. A structure that is not lined up leaves the column alone.
    pub fn residue_picked(&mut self) {
        if let Some(column) = self.structure.column_of_pick() {
            self.alignment.select_column(column, true);
        }
    }
}

#[cfg(test)]
mod tests;
