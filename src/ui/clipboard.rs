//! Copy and paste the way the original does it: features travel with the bases
//! (`clip_cut_copy`, `clip_paste`) and letters that are not bases are asked about first
//! (`atgcfilter`).

use super::app::App;
use super::state::DocView;
use crate::seq::copyspecial::{bankit, convert, fasta, Kind};
use crate::seq::dna::revcom;
use crate::seq::features::Feature;
use crate::seq::io::to_genbank;
use eframe::egui;

/// Letters a paste lets through without asking: the original's allowed keys, with `*` already
/// turned into N.
const PASTE_BASES: &[u8] = b"ACGTNacgtn";

/// Letters that mark a pasted FASTA record (the original's check for "open this as a file").
const RECORD_LETTERS: &[u8] = b"ACGTNBDHKMRSUVWYacgtnbdhkmrsuvwy";

/// The disallowed-letters message shows this many letters, then "...".
const SHOWN_DISALLOWED: usize = 100;

/// What the last copy put on the clipboard, with the features that sat on those bases.
#[derive(Clone)]
struct Copied {
    text: Vec<u8>,
    features: Vec<Feature>,
}

fn memory_id() -> egui::Id {
    egui::Id::new("helixloom_copied_features")
}

/// Pasted text after the original's clean-up, waiting to be inserted.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PasteText {
    /// Spaces, digits and line breaks removed, `*` turned into N.
    pub text: Vec<u8>,
    /// Characters that are not bases, each listed once.
    pub disallowed: Vec<u8>,
}

impl PasteText {
    pub fn new(text: &str) -> PasteText {
        let text: Vec<u8> = text
            .bytes()
            .filter(|c| !(c.is_ascii_whitespace() || *c == 0x0B || c.is_ascii_digit()))
            .map(|c| if c == b'*' { b'N' } else { c })
            .collect();
        let mut disallowed: Vec<u8> = text
            .iter()
            .copied()
            .filter(|c| !PASTE_BASES.contains(c))
            .collect();
        disallowed.sort_unstable();
        disallowed.dedup();
        PasteText { text, disallowed }
    }

    pub fn without_disallowed(&self) -> Vec<u8> {
        let keep = |c: &&u8| PASTE_BASES.contains(c);
        self.text.iter().filter(keep).copied().collect()
    }

    pub fn disallowed_as_n(&self) -> Vec<u8> {
        let fix = |&c: &u8| if PASTE_BASES.contains(&c) { c } else { b'N' };
        self.text.iter().map(fix).collect()
    }

    fn shown_disallowed(&self) -> String {
        let letters = String::from_utf8_lossy(&self.disallowed);
        if letters.chars().count() > SHOWN_DISALLOWED {
            let head: String = letters.chars().take(SHOWN_DISALLOWED + 1).collect();
            format!("{head}...")
        } else {
            letters.into_owned()
        }
    }
}

/// What a paste leaves for the application to finish.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PasteJob {
    /// The text is a GenBank or FASTA record: open it as a document instead of inserting it.
    Record(String),
    /// The text holds letters that are not bases: ask whether to remove them, make them N, or stop.
    Ask { text: PasteText, reverse: bool },
}

/// GenBank (LOCUS … bp … ORIGIN … //) or FASTA (">" then a line of bases) text.
pub fn looks_like_record(text: &str) -> bool {
    let genbank = text.find("LOCUS").is_some_and(|i| {
        let rest = &text[i..];
        let after = |s: &str, from: usize| rest[from..].find(s).map(|p| from + p);
        after(" bp", 0)
            .and_then(|bp| after("ORIGIN", bp))
            .and_then(|origin| after("//", origin))
            .is_some()
    });
    let fasta = text.find('>').is_some_and(|i| {
        text.as_bytes()[i..]
            .windows(2)
            .any(|w| w[0] == b'\n' && RECORD_LETTERS.contains(&w[1]))
    });
    genbank || fasta
}

/// The items of Edit → Copy Special.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Special {
    /// The whole document as a GenBank record.
    Genbank,
    Convert {
        kind: Kind,
        reverse_complement: bool,
    },
    Fasta,
    /// The features as an NCBI Bankit table.
    Bankit,
}

/// The text Copy Special puts on the clipboard: the selection, or everything when nothing is selected.
pub fn copy_special(dv: &DocView, what: Special, table: &str) -> String {
    match what {
        Special::Genbank => to_genbank(&dv.doc),
        Special::Bankit => bankit(&dv.doc),
        Special::Convert {
            kind,
            reverse_complement,
        } => convert(kind, dv.selected_or_all(), reverse_complement, table),
        Special::Fasta => {
            let (a, b) = dv.selection().unwrap_or((0, dv.doc.len()));
            let name = std::path::Path::new(&dv.doc.name)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or(&dv.doc.name);
            fasta(name, a, b, dv.doc.slice(a, b))
        }
    }
}

pub fn copy(ctx: &egui::Context, dv: &DocView, reverse: bool) {
    let Some((a, b)) = dv.selection() else { return };
    let bases = dv.doc.slice(a, b);
    let text = if reverse {
        revcom(bases)
    } else {
        bases.to_vec()
    };
    let features = dv.doc.features_in(a, b, reverse);
    let copied = Copied {
        text: text.clone(),
        features,
    };
    ctx.data_mut(|d| d.insert_temp(memory_id(), copied));
    ctx.copy_text(String::from_utf8_lossy(&text).into_owned());
}

/// Handle pasted text: a record is handed to the application, odd letters are asked about,
/// anything else goes straight in.
pub fn paste(ctx: &egui::Context, dv: &mut DocView, text: &str, reverse: bool) {
    let text = text.replace("\r\n", "\n").replace('\r', "\n");
    if looks_like_record(&text) {
        dv.paste_job = Some(PasteJob::Record(text));
        return;
    }
    let pasted = PasteText::new(&text);
    if pasted.disallowed.is_empty() {
        insert(ctx, dv, pasted.text, reverse);
    } else {
        dv.paste_job = Some(PasteJob::Ask {
            text: pasted,
            reverse,
        });
    }
}

/// Put clean bases over the selection (or at the caret) and select them. Features recorded by the
/// last copy come along when the pasted bases are exactly the ones that were copied.
pub fn insert(ctx: &egui::Context, dv: &mut DocView, bases: Vec<u8>, reverse: bool) {
    if bases.is_empty() {
        return;
    }
    let bases = if reverse { revcom(&bases) } else { bases };
    let features = ctx
        .data(|d| d.get_temp::<Copied>(memory_id()))
        .filter(|c| {
            let source = if reverse {
                revcom(&bases)
            } else {
                bases.clone()
            };
            c.text == source
        })
        .map(|c| {
            let mut features = c.features;
            if reverse {
                features.iter_mut().for_each(|f| f.mirror(0, bases.len()));
            }
            features
        })
        .unwrap_or_default();
    let (a, b) = dv.selection().unwrap_or((dv.caret, dv.caret));
    dv.edit(|d| d.paste(a, b, &bases, features));
    dv.select(a, a + bases.len());
    dv.scroll_to = Some(a + bases.len());
}

impl App {
    /// Finish what the last paste into the current document left open.
    pub(super) fn paste_follow_up(&mut self, ctx: &egui::Context) {
        let Some(i) = self.current_doc() else { return };
        match &self.docs[i].paste_job {
            Some(PasteJob::Record(_)) => {
                if let Some(PasteJob::Record(text)) = self.docs[i].paste_job.take() {
                    self.open_text(&text, "Pasted");
                }
            }
            Some(PasteJob::Ask { .. }) => self.ask_about_pasted_letters(ctx, i),
            None => {}
        }
    }

    fn ask_about_pasted_letters(&mut self, ctx: &egui::Context, i: usize) {
        let Some(PasteJob::Ask { text, reverse }) = self.docs[i].paste_job.clone() else {
            return;
        };
        let mut answer: Option<Option<Vec<u8>>> = None;
        crate::ui::design::window("Disallowed characters being inserted...")
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                ui.label(format!(
                    "Remove disallowed characters: {} ?",
                    text.shown_disallowed()
                ));
                ui.label("You may remove them, convert them to N, or cancel the insertion.");
                ui.horizontal(|ui| {
                    if ui.button("Remove them").clicked() {
                        answer = Some(Some(text.without_disallowed()));
                    }
                    if ui.button("Convert to N").clicked() {
                        answer = Some(Some(text.disallowed_as_n()));
                    }
                    if ui.button("Cancel").clicked() {
                        answer = Some(None);
                    }
                });
            });
        if let Some(choice) = answer {
            self.docs[i].paste_job = None;
            if let Some(bases) = choice {
                insert(ctx, &mut self.docs[i], bases, reverse);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::doc::Doc;

    fn view_with_gene() -> DocView {
        let mut doc = Doc::new("t", b"AAAAAGGGGGCCCCCTTTTT".to_vec());
        doc.features.push(Feature::new("g", "gene", 5, 10, false));
        DocView::new(doc)
    }

    #[test]
    fn copy_special_works_on_the_selection_or_on_everything() {
        const CODE: &str = crate::seq::dna::STANDARD_CODE;
        let mut dv = DocView::new(Doc::new("pDemo.gb", b"ATGgctaGCtaa".to_vec()));
        let upper = Special::Convert {
            kind: Kind::Uppercase,
            reverse_complement: false,
        };
        assert_eq!(copy_special(&dv, upper, CODE), "ATGGC");
        assert_eq!(
            copy_special(&dv, Special::Fasta, CODE),
            ">pDemo 0-12\nATGgctaGCtaa"
        );
        dv.select(3, 9);
        assert_eq!(copy_special(&dv, upper, CODE), "GC");
        assert_eq!(
            copy_special(&dv, Special::Fasta, CODE),
            ">pDemo 3-9\ngctaGC"
        );
        let translated = Special::Convert {
            kind: Kind::Translated,
            reverse_complement: true,
        };
        // gctaGC reverse-complemented is GCtagc, which reads GCT AGC: Ala, Ser
        assert_eq!(copy_special(&dv, translated, CODE), "AS");
        assert!(copy_special(&dv, Special::Genbank, CODE).starts_with("LOCUS"));
        assert_eq!(copy_special(&dv, Special::Bankit, CODE), "");
    }

    #[test]
    fn pasting_drops_spaces_digits_and_line_breaks_and_turns_stars_into_n() {
        let p = PasteText::new("gg\n cc 12*");
        assert_eq!(p.text, b"ggccN");
        assert!(p.disallowed.is_empty());
    }

    #[test]
    fn letters_that_are_not_bases_are_listed_once_and_can_be_removed_or_made_n() {
        let p = PasteText::new("ACGTRYXRX");
        assert_eq!(p.disallowed, b"RXY");
        assert_eq!(p.without_disallowed(), b"ACGT");
        assert_eq!(p.disallowed_as_n(), b"ACGTNNNNN");
    }

    #[test]
    fn a_long_list_of_odd_letters_is_cut_without_splitting_a_character() {
        let many = PasteText {
            text: Vec::new(),
            disallowed: (128u8..=255).collect(),
        };
        let shown = many.shown_disallowed();
        assert!(shown.ends_with("..."));
        assert_eq!(shown.chars().count(), SHOWN_DISALLOWED + 1 + 3);
    }

    #[test]
    fn plain_bases_go_straight_in_and_end_up_selected() {
        let ctx = egui::Context::default();
        let mut dv = DocView::new(Doc::new("t", b"AAAA".to_vec()));
        dv.place_caret(2);
        paste(&ctx, &mut dv, "gg\n cc 12", false);
        assert_eq!(dv.doc.seq, b"AAggccAA");
        assert_eq!(dv.selection(), Some((2, 6)));
        let mut dv = DocView::new(Doc::new("t", b"AAAA".to_vec()));
        dv.place_caret(4);
        paste(&ctx, &mut dv, "GAATTCC", true);
        assert_eq!(dv.doc.seq, b"AAAAGGAATTC");
    }

    #[test]
    fn odd_letters_wait_for_an_answer_and_records_open_as_documents() {
        let ctx = egui::Context::default();
        let mut dv = DocView::new(Doc::new("t", b"AAAA".to_vec()));
        paste(&ctx, &mut dv, "ACGRT", false);
        assert_eq!(dv.doc.seq, b"AAAA");
        assert!(matches!(dv.paste_job, Some(PasteJob::Ask { .. })));
        dv.paste_job = None;
        paste(&ctx, &mut dv, ">seq1\nACGTACGT\n", false);
        assert_eq!(dv.doc.seq, b"AAAA");
        assert!(matches!(dv.paste_job, Some(PasteJob::Record(_))));
    }

    #[test]
    fn records_are_recognised_by_their_markers() {
        assert!(looks_like_record(
            "LOCUS       pX   20 bp    DNA\nORIGIN\n  1 acgt\n//\n"
        ));
        assert!(looks_like_record("x\n>name\nACGT"));
        assert!(!looks_like_record("ACGT GGCC 1234"));
        assert!(!looks_like_record("LOCUS but nothing else"));
    }

    #[test]
    fn copied_features_come_back_only_with_the_same_bases() {
        let ctx = egui::Context::default();
        let mut dv = view_with_gene();
        dv.select(5, 10);
        copy(&ctx, &dv, false);
        dv.place_caret(0);
        paste(&ctx, &mut dv, "GGGGG", false);
        assert_eq!(dv.doc.features.len(), 2);
        assert_eq!(dv.doc.features[1].spans, vec![(0, 5)]);
        // different bases: no features
        dv.place_caret(0);
        paste(&ctx, &mut dv, "TTTTT", false);
        assert_eq!(dv.doc.features.len(), 2);
    }

    #[test]
    fn a_reverse_complement_paste_mirrors_the_copied_features() {
        let ctx = egui::Context::default();
        let mut dv = view_with_gene();
        dv.select(3, 12); // AAGGGGGCC
        copy(&ctx, &dv, false);
        dv.place_caret(20);
        paste(&ctx, &mut dv, "AAGGGGGCC", true);
        let pasted = dv.doc.features.last().unwrap();
        assert!(pasted.reverse);
        assert_eq!(pasted.spans, vec![(22, 27)]);
        assert_eq!(dv.doc.feature_bases(pasted), b"GGGGG");
    }
}
