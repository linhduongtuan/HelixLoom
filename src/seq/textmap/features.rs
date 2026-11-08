//! The features in a Text Map. The original keeps each feature as tags of the text widget that holds the sequence: `f<n>#` on
//! the stretches of the feature, and, for a feature of several stretches, `fn<n>#` from its first base to its last, which
//! `do_text_map` reads to tell the gaps from the parts. `get_tags_exact` (AppMain.tcl L12530) says which tags lie in a row
//! and where, and the lines are made from that (L19794-19845 for the translation of a CDS, L20000-20133 for the rows).

use super::genes::translate;
use super::out::Out;
use super::tcl::{chars, divide, modulo, range, repeat, replace, spaces, text};
use super::{Settings, MARGIN};
use crate::seq::features::Rgb;
use std::collections::BTreeMap;

/// A feature as the page needs it.
#[derive(Clone, Debug, PartialEq)]
pub struct MapFeature {
    pub name: String,
    pub kind: String,
    /// 0-based, half-open stretches, in the order the feature lists them.
    pub spans: Vec<(usize, usize)>,
    pub reverse: bool,
    /// A hidden feature has no colour on the sequence; the page leaves it out unless it is asked to show them.
    pub hidden: bool,
    pub fwd: Rgb,
    pub rev: Rgb,
}

impl MapFeature {
    /// The colour behind the feature on the sequence: none for a hidden one.
    pub fn color(&self) -> Option<Rgb> {
        (!self.hidden).then(|| self.direction_color())
    }

    fn direction_color(&self) -> Rgb {
        if self.reverse {
            self.rev
        } else {
            self.fwd
        }
    }
}

/// A tag of the sequence's text widget.
pub(super) struct SeqTag {
    pub name: String,
    pub feature: usize,
    /// Whether this is the tag that spans a feature of several stretches (`fn<n>#`).
    pub is_fn: bool,
    /// What `tag ranges` says: in order, and a range that touches the next is one with it.
    pub ranges: Vec<(usize, usize)>,
}

impl SeqTag {
    pub(super) fn covers(&self, place: usize) -> bool {
        self.ranges.iter().any(|&(s, e)| s <= place && place < e)
    }

    /// Whether the tag has anything in the bases from `first` to `last`.
    pub(super) fn overlaps(&self, first: usize, last: usize) -> bool {
        match (self.ranges.first(), self.ranges.last()) {
            (Some(a), Some(b)) => a.0 <= last && b.1 > first,
            _ => false,
        }
    }

    /// The places from `first + 1` to `last` where the tag begins or ends. Its ranges do not touch, so these are all of them.
    pub(super) fn turns_between(
        &self,
        first: usize,
        last: usize,
    ) -> impl Iterator<Item = usize> + '_ {
        self.ranges
            .iter()
            .flat_map(|&(s, e)| [s, e])
            .filter(move |&p| p > first && p <= last)
    }
}

fn merge(list: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut sorted: Vec<(usize, usize)> = list.iter().copied().filter(|r| r.0 < r.1).collect();
    sorted.sort_unstable();
    let mut out: Vec<(usize, usize)> = Vec::new();
    for (s, e) in sorted {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

/// The tags of the sequence for a list of features, in the order a file's tags are put on (`put_tags`, L12847): the tags of
/// all the features first, then the `fn` tags of those with several stretches.
pub(super) fn sequence_tags(features: &[MapFeature], len: usize, circular: bool) -> Vec<SeqTag> {
    let mut parts = Vec::new();
    let mut joins = Vec::new();
    for (i, f) in features.iter().enumerate() {
        let mut spans: Vec<(usize, usize)> = f
            .spans
            .iter()
            .map(|&(s, e)| (s, e.min(len)))
            .filter(|&(s, e)| s < e)
            .collect();
        if spans.is_empty() {
            continue;
        }
        if !circular {
            spans.sort_unstable();
        }
        parts.push(SeqTag {
            name: format!("f{i}#"),
            feature: i,
            is_fn: false,
            ranges: merge(&spans),
        });
        if let Some(ranges) = whole_ranges(&spans, len, circular) {
            joins.push(SeqTag {
                name: format!("fn{i}#"),
                feature: i,
                is_fn: true,
                ranges,
            });
        }
    }
    parts.extend(joins);
    parts
}

/// The `fn` tag of a feature of several stretches (`feature_table_to_taglist`, L41228-41245): from the first base to the last,
/// or, when the stretches run through the origin of a circle, from the first base to the end and from the start to the last.
/// Two stretches that are one location through the origin (`a..len` and `0..b`) have none.
fn whole_ranges(
    spans: &[(usize, usize)],
    len: usize,
    circular: bool,
) -> Option<Vec<(usize, usize)>> {
    let (first, last) = (spans.first()?, spans.last()?);
    if spans.len() < 2 {
        return None;
    }
    if circular && first.0 > last.1 {
        if spans.len() == 2 && first.1 == len && last.0 == 0 {
            return None;
        }
        return Some(merge(&[(first.0, len), (0, last.1)]));
    }
    Some(merge(&[(first.0, last.1)]))
}

/// A tag in a stretch of the sequence, as `get_tags_exact` gives it: the places where it begins and ends, counted from the
/// first base of the stretch.
pub(super) struct Win {
    pub tag: usize,
    pub offsets: Vec<usize>,
}

/// The index the original's `bp2ix` makes of base `p` in a text window with rows of `row` bases: line and column, as the
/// numbers they are written with.
fn index_of(p: usize, row: usize) -> (usize, usize) {
    (p / row + 1, p % row)
}

/// Whether two bases have the same index to Tcl's `!=`, which compares two things that look like numbers as numbers:
/// "1.2" and "1.20" (the second base of a row and the twentieth) are the same number.
fn same_number(a: usize, b: usize, row: usize) -> bool {
    let ((line_a, col_a), (line_b, col_b)) = (index_of(a, row), index_of(b, row));
    let digits = |col: usize| col.to_string().trim_end_matches('0').to_string();
    line_a == line_b && digits(col_a) == digits(col_b)
}

/// `get_tags_exact` for the bases from `first` to `last`: every tag on `first`, every place from `first + 1` to `last` where a
/// tag begins or ends, and `last` for a tag that goes on past it. A tag that begins on `last` is not in the stretch, and
/// neither is a tag that begins in a place with the same number as `last` (the original compares the two indexes with `!=`).
pub(super) fn window(tags: &[SeqTag], first: usize, last: usize, row: usize) -> Vec<Win> {
    let mut found = Vec::new();
    for (index, tag) in tags
        .iter()
        .enumerate()
        .filter(|(_, t)| t.overlaps(first, last))
    {
        let mut places = Vec::new();
        if tag.covers(first) {
            places.push(first);
        }
        places.extend(tag.turns_between(first, last));
        if tag.covers(last) {
            if places.first().is_some_and(|&p| same_number(p, last, row)) {
                continue;
            }
            places.push(last);
        }
        if !places.is_empty() {
            found.push(Win {
                tag: index,
                offsets: places.iter().map(|p| p - first).collect(),
            });
        }
    }
    found
}

/// The translation of a CDS, kept while the rows that show it are made: `f_store` of the original.
pub(super) struct CdsStore {
    /// The amino acids, one for each three bases, in the order of the bases of the sequence (turned round to begin at the
    /// origin for a feature through it), and two spaces after the last.
    pub translation: Vec<char>,
    /// How far into it the rows have got.
    pub start: i64,
    /// How many bases of the feature lie before the origin.
    pub split: usize,
    pub reverse: bool,
}

/// `get_cds` (L14593): the bases of a feature, joined, from its first base to its last.
fn cds_text(tag: &SeqTag, whole: Option<&SeqTag>, seq: &[u8]) -> Vec<u8> {
    let piece = |&(s, e): &(usize, usize)| seq[s.min(seq.len())..e.min(seq.len())].to_vec();
    let mut ranges = tag.ranges.clone();
    match whole {
        Some(w) => {
            if w.ranges.len() == 2 {
                ranges = rotated(&ranges, w.ranges[1].0);
            }
        }
        None if ranges.len() == 2 => ranges = vec![ranges[1], ranges[0]],
        None if ranges.len() != 1 => return Vec::new(),
        None => {}
    }
    ranges.iter().flat_map(piece).collect()
}

/// The ranges from the one that begins at `pivot`, then the ones before it.
fn rotated(ranges: &[(usize, usize)], pivot: usize) -> Vec<(usize, usize)> {
    match ranges.iter().position(|r| r.0 == pivot) {
        Some(at) => ranges[at..].iter().chain(&ranges[..at]).copied().collect(),
        None => ranges.to_vec(),
    }
}

/// The translations of the CDS features in `region`, for the rows to take from (`do_text_map`, L19794-19845).
pub(super) fn cds_stores(
    tags: &[SeqTag],
    features: &[MapFeature],
    seq: &[u8],
    region: (usize, usize),
    show_hidden: bool,
    code: &[u8],
    row: usize,
) -> BTreeMap<usize, CdsStore> {
    let mut stores = BTreeMap::new();
    for w in window(tags, region.0, region.1, row) {
        let tag = &tags[w.tag];
        let feature = &features[tag.feature];
        if tag.is_fn
            || (feature.hidden && !show_hidden)
            || !feature.kind.eq_ignore_ascii_case("cds")
        {
            continue;
        }
        let whole = tags.iter().find(|t| t.is_fn && t.feature == tag.feature);
        let bases = cds_text(tag, whole, seq);
        let mut translation = chars(&translate(&bases, code, false, !feature.reverse, 2));
        translation.extend([' ', ' ']);
        let split = split_at_origin(tag, whole);
        if split > 0 {
            let mut turned = range(&translation, split as i64, translation.len() as i64);
            turned.extend(range(&translation, 0, split as i64 - 1));
            translation = turned;
        }
        stores.insert(
            w.tag,
            CdsStore {
                translation,
                start: bases_before(tag, region.0),
                split,
                reverse: feature.reverse,
            },
        );
    }
    stores
}

/// How many bases of a feature lie before the origin: nothing for one that does not run through it. The original tells a
/// feature through the origin by its tag `fn<n>#` having two ranges, or, when it has none, by the tag of the feature having two.
fn split_at_origin(tag: &SeqTag, whole: Option<&SeqTag>) -> usize {
    let length = |r: &(usize, usize)| r.1 - r.0;
    match whole.map_or(0, |w| w.ranges.len()) {
        0 if tag.ranges.len() == 2 => length(&tag.ranges[1]),
        2 => {
            let pivot = whole.map_or(0, |w| w.ranges[1].0);
            let from = tag.ranges.iter().position(|r| r.0 == pivot).unwrap_or(0);
            tag.ranges[from..].iter().map(length).sum()
        }
        _ => 0,
    }
}

/// How many bases of the feature come before the base `start`: the ones the rows that begin there have gone past.
fn bases_before(tag: &SeqTag, start: usize) -> i64 {
    let mut before = 0i64;
    for &(f0, f1) in &tag.ranges {
        if start < f0 {
            break;
        }
        before += if start < f1 { start - f0 } else { f1 - f0 } as i64;
    }
    before
}

/// A line of a row of features: the marks (`>` and `<`, and `.` where a feature of several parts has a gap), the name under
/// them, and, for a CDS, its amino acids above with their numbers.
#[derive(Default)]
struct Line {
    marks: Vec<char>,
    name: Vec<char>,
    /// The tags to put on the marks: the tag, and from where to where.
    tags: Vec<(usize, usize, usize)>,
    amino: Vec<char>,
    numbers: Vec<char>,
}

/// The lines of the features in one row of DNA.
pub(super) struct Rows<'a> {
    tags: &'a [SeqTag],
    features: &'a [MapFeature],
    show_hidden: bool,
    crop: bool,
    width: usize,
    /// How many bases a row of the sequence window has.
    row: usize,
    lines: Vec<Line>,
}

impl<'a> Rows<'a> {
    pub(super) fn new(
        tags: &'a [SeqTag],
        features: &'a [MapFeature],
        s: &Settings,
        row: usize,
    ) -> Rows<'a> {
        Rows {
            tags,
            features,
            show_hidden: s.show_hidden_features,
            crop: s.crop_features,
            width: s.width,
            row,
            lines: Vec::new(),
        }
    }

    /// Work out the lines for the `linelength` bases from `first`: the tags are taken in the order of where they begin in the
    /// row (and, for those that begin together, in the order the sequence has them).
    pub(super) fn collect(
        &mut self,
        first: usize,
        linelength: usize,
        cds: &mut BTreeMap<usize, CdsStore>,
    ) {
        let (tags, features) = (self.tags, self.features);
        let list = window(tags, first, first + linelength, self.row);
        let mut sorted: Vec<&Win> = list.iter().collect();
        sorted.sort_by_key(|w| w.offsets[0]);
        for w in sorted {
            let tag = &tags[w.tag];
            let whole = list
                .iter()
                .find(|o| tags[o.tag].is_fn && tags[o.tag].feature == tag.feature);
            if !tag.is_fn {
                self.feature_line(w, whole, cds);
            } else if !list
                .iter()
                .any(|o| !tags[o.tag].is_fn && tags[o.tag].feature == tag.feature)
            {
                // the row lies inside the gap of a feature of several parts
                let f = &features[tag.feature];
                self.lines.push(Line {
                    marks: repeat(".", linelength as i64),
                    name: chars(&f.name),
                    ..Line::default()
                });
            }
        }
    }

    /// The first line, counting from the top, that the next feature fits on: one whose marks and name end before it begins.
    fn free_line(&mut self, begins: i64) -> usize {
        let taken = |l: &Line| l.marks.len() as i64 > begins || l.name.len() as i64 > begins;
        let at = self
            .lines
            .iter()
            .position(|l| !taken(l))
            .unwrap_or(self.lines.len());
        if at == self.lines.len() {
            self.lines.push(Line::default());
        }
        at
    }

    fn feature_line(&mut self, w: &Win, whole: Option<&Win>, cds: &mut BTreeMap<usize, CdsStore>) {
        let features = self.features;
        let f = &features[self.tags[w.tag].feature];
        if f.hidden && !self.show_hidden {
            return;
        }
        let (start, end) = match whole {
            Some(o) => (o.offsets[0], o.offsets[1]),
            None => (w.offsets[0], w.offsets[1]),
        };
        let (begins, ends) = (start as i64, end as i64);
        let at = self.free_line(begins);
        let mark = if f.reverse { "<" } else { ">" };
        let name = chars(&f.name);
        let line = &mut self.lines[at];
        line.marks.extend(spaces(begins - line.marks.len() as i64));
        line.name.extend(spaces(begins - line.name.len() as i64));
        if whole.is_some() {
            line.marks.extend(repeat(".", ends - begins));
            for pair in w.offsets.chunks(2) {
                if let &[f0, f1] = pair {
                    let marks = repeat(mark, f1 as i64 - f0 as i64);
                    line.marks = replace(&line.marks, f0 as i64, f1 as i64 - 1, &marks);
                    line.tags.push((w.tag, f0, f1));
                }
            }
        } else {
            line.marks.extend(repeat(mark, ends - begins));
            line.tags.push((w.tag, start, end));
        }
        if let Some(store) = cds.get_mut(&w.tag) {
            add_amino_acids(line, store, w, whole.is_some(), (begins, ends));
        }
        line.name.extend(name);
    }

    /// Put the lines down: for each, the amino acids and their numbers above if it has any, the marks, the name.
    pub(super) fn write(&self, out: &mut Out) {
        let margin = " ".repeat(MARGIN);
        for line in &self.lines {
            if !line.amino.is_empty() {
                out.insert(&format!("{margin}{}\n", text(&line.amino)), &[]);
                out.insert(&format!("{margin}{}\n", text(&line.numbers)), &[]);
            }
            let (row, _) = out.line_start();
            out.insert(&format!("{margin}{}\n", text(&line.marks)), &[]);
            for &(t, from, to) in &line.tags {
                let tag = &self.tags[t];
                if !out.has_ranges(&tag.name) {
                    let f = &self.features[tag.feature];
                    out.set_background(&tag.name, Some(f.direction_color()));
                }
                let (a, b) = (
                    out.plus_chars((row, 0), from + MARGIN),
                    out.plus_chars((row, 0), to + MARGIN),
                );
                out.tag_add(&tag.name, a, b);
            }
            let name = if self.crop {
                range(&line.name, 0, self.width as i64 - 1)
            } else {
                line.name.clone()
            };
            out.insert(&format!("{margin}{}\n", text(&name)), &[]);
        }
    }
}

/// The amino acids of a CDS under its marks, and the number of the first (counted along the feature, from the other end for
/// a feature on the other strand).
fn add_amino_acids(line: &mut Line, store: &mut CdsStore, w: &Win, joined: bool, span: (i64, i64)) {
    let (begins, ends) = span;
    line.amino.extend(spaces(begins - line.amino.len() as i64));
    line.numbers
        .extend(spaces(begins - line.numbers.len() as i64));
    let before = range(&store.translation, 0, store.start - 1);
    let counted = before.iter().filter(|c| !c.is_whitespace()).count() as i64;
    let length = store.translation.len() as i64;
    let codons = divide(3.max(length), 3);
    let mut number = modulo(counted + divide(store.split as i64 + 2, 3), codons) + 1;
    if store.reverse {
        number = divide(length, 3) - number + 1;
    }
    if joined {
        line.amino.extend(spaces(ends - begins));
        line.numbers
            .extend(spaces(w.offsets[0] as i64 - 2 - line.numbers.len() as i64));
        for pair in w.offsets.chunks(2) {
            if let &[f0, f1] = pair {
                let size = f1 as i64 - f0 as i64;
                let piece = range(&store.translation, store.start, store.start + size - 1);
                store.start += size;
                line.amino = replace(&line.amino, f0 as i64, f1 as i64 - 1, &piece);
            }
        }
    } else {
        let size = ends - begins;
        line.amino.extend(range(
            &store.translation,
            store.start,
            store.start + size - 1,
        ));
        store.start += size;
    }
    line.numbers.extend(chars(&number.to_string()));
}
