//! The colours of the features on the DNA lines of a Text Map ("Copy Highlighting"): `copy_tags` (AppMain.tcl L12999). Every
//! tag of the sequence but the selection and the `fn` tags is put on the characters of the line that lie over the same bases.

use super::features::{MapFeature, SeqTag};
use crate::seq::features::Rgb;

/// The letters `copy_tags` counts as bases when it moves a tag from the sequence over to the line.
const BASES: &str = "abcdghkmnrstuvwyABCDGHKMNRSTUVWY";

fn is_base(c: Option<&char>) -> bool {
    c.is_some_and(|c| BASES.contains(*c))
}

/// Where the line is cut: every place from `first + 1` to `last` where any tag of the sequence begins or ends. The original
/// asks the text widget for all of them, the tags it does not copy too (the tags that span a feature of several parts, and the
/// selection, if there is one).
fn cuts(
    tags: &[SeqTag],
    selection: Option<(usize, usize)>,
    first: usize,
    last: usize,
) -> Vec<usize> {
    let ends = |(a, b): (usize, usize)| [a, b];
    let mut found: Vec<usize> = tags
        .iter()
        .filter(|t| t.overlaps(first, last))
        .flat_map(|t| t.turns_between(first, last))
        .chain(
            selection
                .into_iter()
                .flat_map(ends)
                .filter(|&p| p > first && p <= last),
        )
        .collect();
    found.sort_unstable();
    found.dedup();
    found
}

/// The tags to put on a line of DNA made from the bases `first` to `last` (inclusive): each as its name, the colour behind it,
/// and the characters of the line it covers (from, up to but not including). The line is cut where a tag begins or ends, and
/// each piece gets the tags that cover its last base. The last piece of a line that ends where the sequence ends gets none
/// when `none_at_end` says so: if the length of the sequence is a whole number of the rows of the sequence window, the place
/// after its last base is the end of the text itself, and the original asks which tags are on the character before that, the
/// line break every text widget has there, and no tag is on that.
pub(super) fn copy_tags(
    tags: &[SeqTag],
    selection: Option<(usize, usize)>,
    (first, last): (usize, usize),
    line: &str,
    features: &[MapFeature],
    seq_length: usize,
    none_at_end: bool,
) -> Vec<(String, Option<Rgb>, usize, usize)> {
    let copied: Vec<&SeqTag> = tags.iter().filter(|t| !t.is_fn).collect();
    let mut cuts = cuts(tags, selection, first, last);
    cuts.push(last + 1);
    let line: Vec<char> = line.chars().collect();
    let (mut to, mut counted) = (0usize, 0usize);
    let mut marks = Vec::new();
    for cut in cuts {
        while to < line.len() && !is_base(line.get(to)) {
            to += 1;
        }
        let from = to;
        while to < line.len() && counted < cut - first {
            to += 1;
            if is_base(line.get(to)) {
                counted += 1;
            }
        }
        let at_end = none_at_end && cut == seq_length;
        for tag in copied.iter().filter(|t| !at_end && t.covers(cut - 1)) {
            marks.push((tag.name.clone(), features[tag.feature].color(), from, to));
        }
    }
    marks
}
