//! What the original writes a Text Map into: a Tk text widget that only grows at its end, with tags on stretches of it. The tag
//! ranges are kept the way `tag ranges` gives them (sorted, with touching ranges made one), so that what this program says
//! can be compared with what the original's widget says, and so that the window can paint the same stretches.
//!
//! One oddity of Tk is kept: a range that runs to `end` (the original's feature lines ask for one when a mark lies beyond the
//! end of its line) goes on over the text that is put after it unless that text is put with a list of tags of its own.

use crate::seq::features::Rgb;
use std::collections::{BTreeMap, BTreeSet};

/// A place in the text, as Tk numbers it: the line (from 1) and the character on it (from 0).
pub type Pos = (usize, usize);

/// The index `end`, which a range can reach but which has no line of its own yet: the text still grows, so where it is
/// is not known until the text is done.
const END: Pos = (usize::MAX, 0);

/// One tag of the finished text: where it lies and the colour behind it (none when the tag has no background).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TagInfo {
    pub name: String,
    pub ranges: Vec<(Pos, Pos)>,
    pub background: Option<Rgb>,
}

/// The ranges of one tag: the closed ones, and where the one that runs to `end` begins, if there is one.
#[derive(Default)]
struct Ranges {
    closed: Vec<(Pos, Pos)>,
    open: Option<Pos>,
}

#[derive(Default)]
pub(super) struct Out {
    text: String,
    /// How many characters each finished line has, without its line break.
    lengths: Vec<usize>,
    line: usize,
    col: usize,
    ranges: BTreeMap<String, Ranges>,
    background: BTreeMap<String, Option<Rgb>>,
    /// The tags in the order they were made, which is their priority: a later one is drawn over an earlier one.
    made: Vec<String>,
    known: BTreeSet<String>,
}

impl Out {
    pub(super) fn new() -> Out {
        Out {
            line: 1,
            ..Out::default()
        }
    }

    /// `insert end text` with the tags the text gets. With none, the text is put as `insert end text` does with no list of
    /// tags: it gets the tags that are on the characters on both sides of it.
    pub(super) fn insert(&mut self, s: &str, tags: &[&str]) {
        self.put(s, tags, !tags.is_empty());
    }

    /// `insert end text {}`: the text gets no tag at all.
    pub(super) fn insert_untagged(&mut self, s: &str) {
        self.put(s, &[], true);
    }

    fn put(&mut self, s: &str, tags: &[&str], own_tags: bool) {
        let from = (self.line, self.col);
        for c in s.chars() {
            self.text.push(c);
            if c == '\n' {
                self.lengths.push(self.col);
                self.line += 1;
                self.col = 0;
            } else {
                self.col += 1;
            }
        }
        let to = (self.line, self.col);
        for ranges in self.ranges.values_mut() {
            let Some(begins) = ranges.open else {
                continue;
            };
            // the range that runs to the end covers the character after this text (the line break every text has at its end)
            // and, unless the text has tags of its own, the one before it if that is in the range
            if own_tags && begins < from {
                ranges.closed.push((begins, from));
            }
            if own_tags || begins >= from {
                ranges.open = Some(to);
            }
        }
        for tag in tags {
            self.tag_add(tag, from, to);
        }
    }

    /// Where `index "end-1char linestart"` is: the start of the line the next text goes on.
    pub(super) fn line_start(&self) -> Pos {
        (self.line, 0)
    }

    /// Where `index "<base> + <n>chars"` is: the characters are counted over the ends of lines, and an index that falls after
    /// the last character is `end`.
    pub(super) fn plus_chars(&self, base: Pos, n: usize) -> Pos {
        let (mut line, mut col, mut left) = (base.0, base.1, n);
        loop {
            let length = self.lengths.get(line - 1).copied().unwrap_or(self.col);
            let to_end = length.saturating_sub(col) + 1;
            if left < to_end {
                return (line, col + left);
            }
            left -= to_end;
            if line >= self.line {
                return END;
            }
            line += 1;
            col = 0;
        }
    }

    pub(super) fn tag_add(&mut self, tag: &str, from: Pos, to: Pos) {
        if from >= to {
            return;
        }
        self.make(tag);
        if let Some(ranges) = self.ranges.get_mut(tag) {
            if to == END {
                ranges.open = Some(ranges.open.map_or(from, |old| old.min(from)));
            } else {
                ranges.closed.push((from, to));
            }
        }
    }

    fn make(&mut self, tag: &str) {
        if self.known.insert(tag.to_string()) {
            self.made.push(tag.to_string());
            self.ranges.insert(tag.to_string(), Ranges::default());
        }
    }

    /// Whether the tag has any text yet (`tag ranges` is not empty).
    pub(super) fn has_ranges(&self, tag: &str) -> bool {
        self.ranges
            .get(tag)
            .is_some_and(|r| !r.closed.is_empty() || r.open.is_some())
    }

    /// `tag configure -background`.
    pub(super) fn set_background(&mut self, tag: &str, color: Option<Rgb>) {
        self.make(tag);
        self.background.insert(tag.to_string(), color);
    }

    /// `tag raise`: the tag is drawn over the others from now on.
    pub(super) fn raise(&mut self, tag: &str) {
        if let Some(at) = self.made.iter().position(|t| t == tag) {
            let name = self.made.remove(at);
            self.made.push(name);
        }
    }

    pub(super) fn has_tag(&self, tag: &str) -> bool {
        self.known.contains(tag)
    }

    pub(super) fn finish(self) -> (String, Vec<TagInfo>) {
        let end = (self.line + 1, 0);
        let mut tags = Vec::new();
        for name in &self.made {
            let mut list = Vec::new();
            if let Some(r) = self.ranges.get(name) {
                list.extend(r.closed.iter().copied());
                list.extend(r.open.map(|from| (from, end)));
            }
            let background = self.background.get(name).copied().flatten();
            tags.push(TagInfo {
                name: name.clone(),
                ranges: merged(&list),
                background,
            });
        }
        (self.text, tags)
    }
}

/// The ranges in order, with ranges that overlap or touch made one: how a text widget reports a tag's ranges.
fn merged(list: &[(Pos, Pos)]) -> Vec<(Pos, Pos)> {
    let mut sorted = list.to_vec();
    sorted.sort_unstable();
    let mut out: Vec<(Pos, Pos)> = Vec::new();
    for (from, to) in sorted {
        match out.last_mut() {
            Some(last) if from <= last.1 => last.1 = last.1.max(to),
            _ => out.push((from, to)),
        }
    }
    out
}
