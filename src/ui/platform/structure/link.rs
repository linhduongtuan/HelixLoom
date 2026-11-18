//! The structure page's side of the one selection (`ui/platform/linked.rs`): the picked residue, and the column of the
//! alignment that it was lined up with.

use super::StructurePage;
use crate::seq::fold::Fold;

impl StructurePage {
    /// The residue that is picked, as its place among the residues of the structure.
    pub fn picked(&self) -> Option<usize> {
        self.style.picked
    }

    /// Pick a residue (none: let it go), as a click on it does.
    #[cfg(test)]
    pub(in crate::ui) fn pick(&mut self, residue: Option<usize>) {
        self.style.picked = residue;
    }

    /// The structure, if it is lined up with the alignment: some residue knows its column (`Fold::apply_conservation`).
    fn lined_up(&self) -> Option<&Fold> {
        let fold = &self.loaded.as_ref()?.fold;
        fold.residues
            .iter()
            .any(|r| r.column.is_some())
            .then_some(fold)
    }

    /// Pick the residue that was lined up with `column` of the alignment (the first, when the copies of a chain all were);
    /// no column, or a column that no residue is in, picks none. A structure that is not lined up is left as it is.
    pub fn pick_column(&mut self, column: Option<usize>) {
        let Some(fold) = self.lined_up() else {
            return;
        };
        self.style.picked = column.and_then(|k| fold.in_column(k).first().copied());
    }

    /// The column of the alignment that the picked residue was lined up with (none: nothing is picked, or no row has that
    /// residue); nothing at all for a structure that is not lined up with the alignment.
    pub fn column_of_pick(&self) -> Option<Option<usize>> {
        let fold = self.lined_up()?;
        let picked = self.style.picked.and_then(|k| fold.residues.get(k));
        Some(picked.and_then(|r| r.column))
    }

    /// The residues that were lined up with column `k` of the alignment, each as the page names it: the one that is picked
    /// first (the copies of a chain are all in the column, and the link of the Alignment page names the first).
    pub fn residues_in_column(&self, k: usize) -> Vec<String> {
        let Some(loaded) = &self.loaded else {
            return Vec::new();
        };
        let fold = &loaded.fold;
        let mut places = fold.in_column(k);
        // a stable sort: the others keep the order of the file
        places.sort_by_key(|&place| Some(place) != self.style.picked);
        places
            .into_iter()
            .map(|place| {
                let r = &fold.residues[place];
                format!("{} {} of chain {}", char::from(r.letter), r.number, r.chain)
            })
            .collect()
    }
}
