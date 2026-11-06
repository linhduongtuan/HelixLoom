//! One open sequence: bases, features and the edits that keep them consistent.

use super::dna::revcom;
use super::features::{color_from_name, Feature, DEFAULT_FWD};
use std::path::PathBuf;

/// Type and colour (both strands) of the feature the original's `insert_enzyme_site` adds.
const SITE_FEATURE_TYPE: &str = "protein_bind";
const SITE_COLOR: &str = "#C90000";

#[derive(Clone, Debug, PartialEq)]
pub struct Doc {
    pub name: String,
    pub seq: Vec<u8>,
    pub circular: bool,
    pub features: Vec<Feature>,
    pub definition: String,
    /// The file comment, one line after the other: the GenBank COMMENT lines other than ApE's own
    /// `ApEinfo:` ones (a PCR product keeps its primers and warnings here).
    pub comment: String,
    /// Dam/Dcm methylated DNA (plasmid) as opposed to PCR product.
    pub methylated: bool,
    pub path: Option<PathBuf>,
    pub dirty: bool,
}

impl Doc {
    pub fn new(name: impl Into<String>, seq: Vec<u8>) -> Doc {
        Doc {
            name: name.into(),
            seq,
            circular: false,
            features: Vec::new(),
            definition: String::new(),
            comment: String::new(),
            methylated: true,
            path: None,
            dirty: false,
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Insert bases before position `at`.
    pub fn insert(&mut self, at: usize, bases: &[u8]) {
        let at = at.min(self.seq.len());
        self.seq.splice(at..at, bases.iter().copied());
        self.features
            .iter_mut()
            .for_each(|f| f.after_insert(at, bases.len()));
        self.dirty = true;
    }

    /// Delete bases `[a, b)`.
    pub fn delete(&mut self, a: usize, b: usize) {
        let (a, b) = (a.min(self.seq.len()), b.min(self.seq.len()));
        if b <= a {
            return;
        }
        self.seq.drain(a..b);
        self.features.retain_mut(|f| f.after_delete(a, b));
        self.dirty = true;
    }

    pub fn replace(&mut self, a: usize, b: usize, bases: &[u8]) {
        self.delete(a, b);
        self.insert(a.min(self.seq.len()), bases);
    }

    pub fn slice(&self, a: usize, b: usize) -> &[u8] {
        &self.seq[a.min(self.seq.len())..b.min(self.seq.len())]
    }

    /// Put `bases` over the bases from `at` on, one for one; features stay where they are
    /// (the original's change-case takes the text out and back in but restores every tag).
    pub fn overwrite(&mut self, at: usize, bases: &[u8]) {
        let end = (at + bases.len()).min(self.seq.len());
        if at >= end {
            return;
        }
        self.seq[at..end].copy_from_slice(&bases[..end - at]);
        self.dirty = true;
    }

    /// Reverse-complement `[a, b)` in place. Features inside it are mirrored and change direction;
    /// one that crosses the edge is split: the outside part stays, the inside part becomes a new
    /// mirrored feature (how the original's `rev_com_window` restores tags through `put_tags`).
    pub fn reverse_complement_range(&mut self, a: usize, b: usize) {
        let (a, b) = (a.min(self.seq.len()), b.min(self.seq.len()));
        if b <= a {
            return;
        }
        let rc = revcom(&self.seq[a..b]);
        self.seq[a..b].copy_from_slice(&rc);
        let mut split = Vec::new();
        for f in &mut self.features {
            let Some(mut inner) = f.inside(a, b) else {
                continue;
            };
            inner.mirror(a, b);
            if f.remove_range(a, b) {
                split.push(inner);
            } else {
                *f = inner;
            }
        }
        self.features.extend(split);
        self.dirty = true;
    }

    /// Features on `[a, b)`, cut to that range and counted from `a`; mirrored when `reverse`.
    /// This is what a copy keeps beside the bases (`clip_cut_copy`).
    pub fn features_in(&self, a: usize, b: usize, reverse: bool) -> Vec<Feature> {
        self.features
            .iter()
            .filter_map(|f| f.inside(a, b))
            .map(|mut f| {
                for span in &mut f.spans {
                    *span = (span.0 - a, span.1 - a);
                }
                if reverse {
                    f.mirror(0, b - a);
                }
                f
            })
            .collect()
    }

    /// A new document made of `range`, or of all of this one (`duplicate_window`, AppMain.tcl L7750 @0x5c69b). A stretch
    /// becomes a linear document with the features inside it; the whole keeps the circle, every feature and the methylation.
    /// The name is the original's `New_DNA`, and the comment starts with where it came from, `name from a to b`, as the
    /// original writes it (the file's name without its extension, or the name of the window), then the old comment.
    pub fn duplicate(&self, range: Option<(usize, usize)>) -> Doc {
        let (a, b) = match range {
            Some((a, b)) => (a.min(b.min(self.len())), b.min(self.len())),
            None => (0, self.len()),
        };
        let source = match self.path.as_ref().and_then(|p| p.file_stem()) {
            Some(stem) => stem.to_string_lossy().into_owned(),
            None => self.name.clone(),
        };
        let mut copy = Doc::new("New_DNA", self.slice(a, b).to_vec());
        copy.circular = range.is_none() && self.circular;
        copy.features = match range {
            Some(_) => self.features_in(a, b, false),
            None => self.features.clone(),
        };
        copy.definition = self.definition.clone();
        copy.methylated = self.methylated;
        // the original's comment box ends without a newline when the old comment is empty; here no empty line is kept
        let from = format!("{source} from {} to {b}", a + 1);
        copy.comment = [from.as_str(), self.comment.as_str()]
            .join("\n")
            .trim_end_matches('\n')
            .to_string();
        copy
    }

    /// Replace `[a, b)` with `bases` and lay the copied `features` over the new bases (`clip_paste`).
    pub fn paste(&mut self, a: usize, b: usize, bases: &[u8], features: Vec<Feature>) {
        self.replace(a, b, bases);
        self.features
            .extend(features.into_iter().map(|f| f.shifted(a)));
    }

    /// Insert an enzyme's recognition sequence and mark it with a `protein_bind` feature named
    /// after the enzyme (`insert_enzyme_site`); a feature the site lands inside grows with it.
    pub fn insert_site(&mut self, at: usize, enzyme: &str, bases: &[u8], reverse: bool) {
        self.insert(at, bases);
        let mut f = Feature::new(enzyme, SITE_FEATURE_TYPE, at, at + bases.len(), reverse);
        let color = color_from_name(SITE_COLOR).unwrap_or(DEFAULT_FWD);
        (f.fwd_color, f.rev_color) = (color, color);
        self.features.push(f);
    }

    /// Rotate a circular sequence so base `at` becomes the first base.
    pub fn set_origin(&mut self, at: usize) {
        let n = self.seq.len();
        if n == 0 || at == 0 || at >= n {
            return;
        }
        self.seq.rotate_left(at);
        self.features.iter_mut().for_each(|f| f.rotate(at, n));
        self.dirty = true;
    }

    pub fn reverse_complement(&mut self) {
        let n = self.seq.len();
        self.seq = revcom(&self.seq);
        self.features.iter_mut().for_each(|f| f.flip(n));
        self.dirty = true;
    }

    /// The feature's bases in reading order (joined spans, reverse-complemented for the minus strand).
    pub fn feature_bases(&self, f: &Feature) -> Vec<u8> {
        let mut spans = f.spans.clone();
        spans.sort_unstable();
        // a feature crossing the origin reads to the end of the sequence, then continues at base 1
        if spans.len() > 1 && spans[0].0 == 0 && spans.last().is_some_and(|s| s.1 == self.seq.len())
        {
            spans.rotate_right(1);
        }
        let joined: Vec<u8> = spans
            .iter()
            .flat_map(|&(s, e)| self.slice(s, e).iter().copied())
            .collect();
        if f.reverse {
            revcom(&joined)
        } else {
            joined
        }
    }

    pub fn title(&self) -> String {
        if self.dirty {
            format!("{} •", self.name)
        } else {
            self.name.clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc() -> Doc {
        let mut d = Doc::new("t", b"AAAAAGGGGGCCCCCTTTTT".to_vec());
        d.features.push(Feature::new("g", "gene", 5, 10, false));
        d.features.push(Feature::new("c", "gene", 10, 15, true));
        d
    }

    #[test]
    fn inserting_before_a_feature_moves_it() {
        let mut d = doc();
        d.insert(2, b"TT");
        assert_eq!(d.features[0].spans, vec![(7, 12)]);
        assert_eq!(d.slice(7, 12), b"GGGGG");
        assert!(d.dirty);
    }

    #[test]
    fn deleting_inside_shrinks_and_deleting_all_removes() {
        let mut d = doc();
        d.delete(6, 8);
        assert_eq!(d.features[0].spans, vec![(5, 8)]);
        assert_eq!(d.features[1].spans, vec![(8, 13)]);
        d.delete(8, 13);
        assert_eq!(d.features.len(), 1);
        assert_eq!(d.features[0].name, "g");
    }

    #[test]
    fn rotating_the_origin_keeps_features_on_their_bases() {
        let mut d = doc();
        d.set_origin(5);
        assert_eq!(&d.seq[..5], b"GGGGG");
        assert_eq!(d.features[0].spans, vec![(0, 5)]);
        assert_eq!(d.slice(0, 5), b"GGGGG");
        assert_eq!(d.features[1].spans, vec![(5, 10)]);
    }

    #[test]
    fn reverse_complement_flips_features_with_the_bases() {
        let mut d = doc();
        let before = d.feature_bases(&d.features[1]);
        d.reverse_complement();
        let g = d.features.iter().find(|f| f.name == "c").unwrap();
        assert!(!g.reverse);
        assert_eq!(d.feature_bases(g), before);
        assert_eq!(d.seq, revcom(b"AAAAAGGGGGCCCCCTTTTT"));
    }

    #[test]
    fn features_crossing_the_origin_read_through_it() {
        let mut d = Doc::new("w", b"TTTAAAAAAAAAAGGG".to_vec());
        let mut f = Feature::new("wrap", "gene", 13, 16, false);
        f.spans = vec![(0, 3), (13, 16)];
        assert_eq!(d.feature_bases(&f), b"GGGTTT");
        f.reverse = true;
        assert_eq!(d.feature_bases(&f), b"AAACCC");
        d.circular = true;
    }

    #[test]
    fn overwriting_keeps_every_feature_in_place() {
        let mut d = doc();
        d.overwrite(5, b"gggggccccc");
        assert_eq!(d.seq, b"AAAAAgggggcccccTTTTT");
        assert_eq!(d.features.len(), 2);
        assert_eq!(d.features[0].spans, vec![(5, 10)]);
        assert_eq!(d.features[1].spans, vec![(10, 15)]);
        d.overwrite(18, b"acgtacgt"); // runs past the end: only what fits is written
        assert_eq!(&d.seq[18..], b"ac");
    }

    #[test]
    fn reverse_complementing_a_range_mirrors_the_features_inside_it() {
        let mut d = doc();
        let before: Vec<Vec<u8>> = d.features.iter().map(|f| d.feature_bases(f)).collect();
        d.reverse_complement_range(5, 15);
        assert_eq!(d.features.len(), 2);
        let (g, c) = (&d.features[0], &d.features[1]);
        assert_eq!((g.spans.clone(), g.reverse), (vec![(10, 15)], true));
        assert_eq!((c.spans.clone(), c.reverse), (vec![(5, 10)], false));
        // each feature still reads the bases it read before
        assert_eq!(d.feature_bases(g), before[0]);
        assert_eq!(d.feature_bases(c), before[1]);
    }

    #[test]
    fn a_feature_crossing_the_edge_of_the_range_is_split_in_two() {
        let mut d = doc();
        d.reverse_complement_range(2, 8); // AAAGGG -> CCCTTT, the gene at 5..10 is cut at 8
        assert_eq!(d.seq, b"AACCCTTTGGCCCCCTTTTT");
        assert_eq!(d.features.len(), 3);
        assert_eq!(d.features[0].spans, vec![(8, 10)]);
        assert!(!d.features[0].reverse);
        let inner = &d.features[2];
        assert_eq!(inner.name, "g");
        assert_eq!((inner.spans.clone(), inner.reverse), (vec![(2, 5)], true));
        assert_eq!(d.feature_bases(inner), b"GGG");
    }

    #[test]
    fn copied_features_come_back_with_the_pasted_bases() {
        let mut d = doc();
        let copied = d.features_in(5, 15, false);
        assert_eq!(copied[0].spans, vec![(0, 5)]);
        assert_eq!(
            (copied[1].spans.clone(), copied[1].reverse),
            (vec![(5, 10)], true)
        );
        let mirrored = d.features_in(5, 15, true);
        assert_eq!(
            (mirrored[0].spans.clone(), mirrored[0].reverse),
            (vec![(5, 10)], true)
        );
        assert_eq!(
            (mirrored[1].spans.clone(), mirrored[1].reverse),
            (vec![(0, 5)], false)
        );
        d.paste(0, 0, b"GGGGGCCCCC", copied);
        assert_eq!(d.features.len(), 4);
        assert_eq!(d.features[0].spans, vec![(15, 20)]); // the old ones moved along
        assert_eq!(d.features[2].spans, vec![(0, 5)]); // the new ones sit on the pasted bases
        assert_eq!(d.slice(0, 5), b"GGGGG");
    }

    #[test]
    fn an_inserted_site_gets_a_protein_bind_feature_and_grows_the_feature_around_it() {
        let mut d = doc();
        d.insert_site(7, "EcoRI", b"GAATTC", true);
        assert_eq!(d.features[0].spans, vec![(5, 16)]); // the gene grew by the site
        let site = d.features.last().unwrap();
        assert_eq!(
            (site.name.as_str(), site.kind.as_str()),
            ("EcoRI", "protein_bind")
        );
        assert_eq!(site.spans, vec![(7, 13)]);
        assert!(site.reverse);
        assert_eq!(
            (site.fwd_color, site.rev_color),
            ([0xC9, 0, 0], [0xC9, 0, 0])
        );
        assert_eq!(d.slice(7, 13), b"GAATTC");
    }

    #[test]
    fn a_duplicated_stretch_is_linear_and_keeps_the_features_inside_it() {
        let mut d = doc();
        d.circular = true;
        d.comment = "made by hand".into();
        let copy = d.duplicate(Some((3, 12)));
        assert_eq!(copy.name, "New_DNA");
        assert_eq!(copy.seq, b"AAGGGGGCC".to_vec());
        assert!(!copy.circular, "a stretch is never a circle");
        // the gene 5..10 sits at 2..7 in the copy; the second feature, 10..15, is cut to 10..12 and so is 7..9
        assert_eq!(copy.features.len(), 2);
        assert_eq!(copy.features[0].spans, vec![(2, 7)]);
        assert_eq!(copy.features[1].spans, vec![(7, 9)]);
        assert_eq!(copy.comment, "t from 4 to 12\nmade by hand");
    }

    #[test]
    fn a_duplicated_whole_keeps_the_circle_the_features_and_the_methylation() {
        let mut d = doc();
        d.circular = true;
        d.methylated = false;
        d.definition = "a plasmid".into();
        let copy = d.duplicate(None);
        assert_eq!(copy.seq, d.seq);
        assert!(copy.circular && !copy.methylated);
        assert_eq!(copy.features, d.features);
        assert_eq!(copy.definition, "a plasmid");
        assert_eq!(copy.comment, "t from 1 to 20", "no empty line after it");
    }

    #[test]
    fn a_duplicate_names_the_file_without_its_extension_and_survives_a_bad_range() {
        let mut d = doc();
        d.path = Some(PathBuf::from("/data/pUC19.gb"));
        assert_eq!(d.duplicate(None).comment, "pUC19 from 1 to 20");
        // the range is cut to the sequence, and an upside-down one is empty
        assert_eq!(d.duplicate(Some((15, 99))).seq, b"TTTTT".to_vec());
        assert!(d.duplicate(Some((30, 10))).seq.is_empty());
    }

    #[test]
    fn feature_bases_read_the_minus_strand() {
        let d = doc();
        assert_eq!(d.feature_bases(&d.features[1]), b"GGGGG");
        assert_eq!(d.feature_bases(&d.features[0]), b"GGGGG");
    }
}
