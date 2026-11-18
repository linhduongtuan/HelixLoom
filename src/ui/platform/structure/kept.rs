//! The Structure page as the session keeps it: the structure that is shown, as the PDB text it was read from, under its name, with
//! the record of the prediction that made it when there is one. How it is drawn and turned is not kept.

use super::StructurePage;
use crate::kept::{fields, line};
use crate::provenance::Provenance;

impl StructurePage {
    /// A number that is another one whenever what `kept` gives has changed; none when no structure is shown. Asked at every frame.
    pub(in crate::ui) fn stamp(&self) -> Option<u64> {
        self.loaded
            .as_ref()
            .map(|_| self.shown * 2 + u64::from(self.made.is_some()))
    }

    /// The structure as text that `put_back` reads; none when no structure is shown.
    pub(in crate::ui) fn kept(&self) -> Option<String> {
        let loaded = self.loaded.as_ref()?;
        let mut text = line(&["name", &loaded.name]) + &line(&["pdb", &loaded.text]);
        if let Some(made) = &self.made {
            text += &line(&["made", &made.kept()]);
        }
        Some(text)
    }

    /// Show the structure that `text` (of `kept`) holds. A text with no structure that can be read changes nothing: why is said.
    pub(in crate::ui) fn put_back(&mut self, text: &str) -> Result<(), String> {
        let (mut name, mut pdb, mut made) = (String::new(), None, None);
        for kept in text.lines() {
            let kept = fields(kept);
            let kept: Vec<&str> = kept.iter().map(String::as_str).collect();
            match kept.as_slice() {
                ["name", held] => name = held.to_string(),
                ["pdb", held] => pdb = Some(held.to_string()),
                ["made", record] => made = Provenance::from_kept(record),
                _ => {}
            }
        }
        let pdb = pdb.ok_or("Nothing of the structure that was kept can be read.")?;
        self.show_pdb(&name, &pdb)?;
        self.made = made;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::{Scheme, EXAMPLE_PDB};
    use super::*;
    use crate::provenance::{Asked, Input};

    fn record() -> Provenance {
        Asked {
            analysis: "structure prediction",
            inputs: vec![Input::of("protein", b"MQIFVKTLTG")],
            parameters: vec![("service".to_string(), "ESMFold\tof a lab".to_string())],
            ran_on: "the workstation",
        }
        .made(vec!["ESMFold v1".to_string()], 12.5)
    }

    fn page() -> StructurePage {
        let mut page = StructurePage::default();
        page.show_pdb("a name\twith a tab", EXAMPLE_PDB)
            .expect("the example");
        page.made = Some(record());
        page
    }

    #[test]
    fn the_structure_comes_back_as_the_text_it_was_read_from_under_its_name_with_its_record() {
        let text = page().kept().expect("a structure");
        let mut back = StructurePage::default();
        back.put_back(&text).expect("put back");
        let loaded = back.loaded.as_ref().expect("a structure");
        assert_eq!(loaded.name, "a name\twith a tab");
        assert_eq!(loaded.text, EXAMPLE_PDB);
        assert_eq!(loaded.fold.residues.len(), 76);
        assert_eq!(back.made, Some(record()));
        assert_eq!(back.kept().as_deref(), Some(text.as_str()));
        // one that was opened, which nothing here made, comes back with no record
        let mut opened = page();
        opened.made = None;
        back.put_back(&opened.kept().expect("a structure"))
            .expect("put back");
        assert!(back.loaded.is_some() && back.made.is_none());
    }

    #[test]
    fn the_stamp_is_another_one_for_every_structure_shown_and_for_its_record_and_not_for_how_it_is_drawn(
    ) {
        let mut shown = StructurePage::default();
        assert_eq!((shown.stamp(), shown.kept()), (None, None));
        shown.show_pdb("one", EXAMPLE_PDB).expect("the example");
        let first = shown.stamp().expect("a stamp");
        assert_eq!(shown.stamp(), Some(first));
        shown.made = Some(record());
        let with_record = shown.stamp().expect("a stamp");
        // the same text shown again is still another showing: the record went with the one before
        shown.show_pdb("one", EXAMPLE_PDB).expect("the example");
        let again = shown.stamp().expect("a stamp");
        assert!(first != with_record && again != first && again != with_record);
        (shown.style.scheme, shown.spin) = (Scheme::Chain, true);
        shown.message = "said".to_string();
        assert_eq!(shown.stamp(), Some(again));
    }

    #[test]
    fn a_text_with_no_structure_that_can_be_read_changes_nothing() {
        let mut shown = page();
        let before = (shown.stamp(), shown.kept());
        for (text, why) in [
            ("", "Nothing of the structure"),
            ("name\tx\nmade\tnot a record\n", "Nothing of the structure"),
            ("name\tx\npdb\tnot a structure\n", ""),
            ("name\tx\npdb\tHETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n", ""),
        ] {
            let said = shown.put_back(text).expect_err(text);
            assert!(!said.is_empty() && said.contains(why), "{said}");
            assert_eq!((shown.stamp(), shown.kept()), before);
        }
        // a line from a later version is passed over, and a record that cannot be read is no record
        let odd = format!(
            "{}later\tfields\nmade\tnot a record\n",
            page().kept().expect("a structure")
        );
        shown.put_back(&odd).expect("put back");
        assert!(shown.loaded.is_some() && shown.made.is_none());
    }
}
