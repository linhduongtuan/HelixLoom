//! The consensus of a template and the reads that support it. NECAT uses the consensus of FALCON, which its paper names without describing;
//! this one is **the program's own**: a vote at each position of the template, with the template's base as one vote.
//!
//! * A position that fewer than `min_support` supporting reads cover is **not corrected**.
//! * Where more of the reads delete the position than have a base there (the template counted), it is left out.
//! * Otherwise the base with the most evidence from the qualities (`pileup::weight`, the template's base counted at the default quality) is
//!   written; two that tie keep the template's.
//! * After it, an insertion is written when over half of the reads that have a base there (the template counted, which has none) insert the
//!   same bases.
//!
//! In **every round but the last** (the first of the two) a position that is not corrected keeps the template's base, so the whole read comes
//! out; and each alignment counts only over its longest run of blocks of 500 bases whose error is at most half (the paper aligns blocks
//! outwards from a seed and stops at a bad one). In the **last round** the whole alignments count and only the corrected stretches come out,
//! each a read of its own.

use super::select::Aligned;
use super::Params;
use crate::seq::variants::pileup::{base_index, pile, weight, Filter, Pileup};
use crate::seq::variants::sam::{uses_read, uses_target, Record};

/// The part of a record that lies on `[from, to)` of the target, with its bases; `None` when nothing is aligned there.
pub fn clip(r: &Record, from: usize, to: usize) -> Option<Record> {
    let (mut at, mut read_at) = (r.pos, 0);
    let (mut cigar, mut seq, mut qual, mut start) = (Vec::new(), Vec::new(), Vec::new(), None);
    for &(op, n) in &r.cigar {
        let take = match op {
            'M' | '=' | 'X' | 'D' => {
                let (a, b) = (at.max(from), (at + n).min(to));
                (a < b).then(|| (a - at, b - a))
            }
            // an insertion belongs between two bases of the target: it is kept when both lie inside
            'I' => (at > from && at < to).then_some((0, n)),
            _ => None,
        };
        if let Some((skip, length)) = take {
            if uses_target(op) {
                start.get_or_insert(at + skip);
            }
            if uses_read(op) {
                let first = read_at + skip;
                seq.extend_from_slice(r.seq.get(first..first + length)?);
                if !r.qual.is_empty() {
                    qual.extend_from_slice(r.qual.get(first..first + length)?);
                }
            }
            cigar.push((op, length));
        }
        if uses_read(op) {
            read_at += n;
        }
        if uses_target(op) {
            at += n;
        }
    }
    Some(Record {
        pos: start?,
        cigar,
        seq,
        qual,
        ..r.clone()
    })
}

/// The longest run of blocks of the template (`block` bases each) in which the record's error is at most `max_error`, as positions of the
/// template; `None` when no block is.
pub fn good_run(
    r: &Record,
    template: &[u8],
    block: usize,
    max_error: f64,
) -> Option<(usize, usize)> {
    let blocks = template.len().div_ceil(block.max(1));
    let (mut columns, mut errors) = (vec![0usize; blocks], vec![0usize; blocks]);
    let (mut at, mut read_at) = (r.pos, 0);
    for &(op, n) in &r.cigar {
        for i in 0..n {
            let t = at + if uses_target(op) { i } else { 0 };
            let b = (t / block.max(1)).min(blocks.saturating_sub(1));
            let wrong = match op {
                'M' | '=' | 'X' => {
                    let read = r.seq.get(read_at + i);
                    Some(
                        read.zip(template.get(t))
                            .is_none_or(|(x, y)| !x.eq_ignore_ascii_case(y)),
                    )
                }
                'I' | 'D' => Some(true),
                _ => None,
            };
            if let Some(wrong) = wrong {
                columns[b] += 1;
                errors[b] += usize::from(wrong);
            }
        }
        if uses_read(op) {
            read_at += n;
        }
        if uses_target(op) {
            at += n;
        }
    }
    let mut best: Option<(usize, usize)> = None;
    let mut start = None;
    for b in 0..=blocks {
        let good =
            b < blocks && columns[b] > 0 && errors[b] as f64 <= max_error * columns[b] as f64;
        match (good, start) {
            (true, None) => start = Some(b),
            (false, Some(s)) => {
                if best.is_none_or(|(x, y)| b - s > y - x) {
                    best = Some((s, b));
                }
                start = None;
            }
            _ => {}
        }
    }
    best.map(|(a, b)| (a * block, (b * block).min(template.len())))
}

/// The template as it is voted: its corrected bases, and for each position whether it was corrected.
#[derive(Debug, PartialEq)]
pub struct Voted {
    /// For each position of the template: the bases that stand for it (none where it is deleted, several where an insertion follows), and
    /// whether enough reads covered it.
    pub positions: Vec<(Vec<u8>, bool)>,
}

fn vote_at(
    pile: &Pileup,
    p: usize,
    template: u8,
    default_quality: u8,
    min_support: usize,
) -> (Vec<u8>, bool) {
    if pile.depth(p) < min_support {
        return (vec![template], false);
    }
    let mut out = Vec::new();
    let cover = usize::from(pile.cover[p]) + 1;
    if usize::from(pile.deleted[p]) <= cover {
        let own = base_index(template);
        let mut scores = pile.scores[p];
        if let Some(b) = own {
            scores[b] += weight(default_quality);
        }
        // the template's base first, so that it stays on a tie
        let first = own.unwrap_or(0);
        let best = (0..4).fold(first, |b, a| if scores[a] > scores[b] { a } else { b });
        out.push(if scores[best] > 0.0 {
            b"ACGT"[best]
        } else {
            template
        });
    }
    if let Some(found) = pile.inserted.get(&p) {
        if let Some((bases, &n)) = found.iter().max_by_key(|(_, &n)| n) {
            if 2 * usize::from(n) > cover {
                out.extend_from_slice(bases);
            }
        }
    }
    (out, true)
}

/// Votes the template over the alignments of its supports; with `blocks`, each counts only over its longest run of good blocks.
pub fn vote(template: &[u8], supports: &[Aligned], blocks: bool, p: &Params) -> Voted {
    let records: Vec<Record> = supports
        .iter()
        .filter_map(|a| {
            if blocks {
                let (from, to) = good_run(&a.record, template, p.block, p.max_block_error)?;
                clip(&a.record, from, to)
            } else {
                Some(a.record.clone())
            }
        })
        .collect();
    let filter = Filter {
        min_base_quality: 0,
        max_depth: usize::from(u16::MAX),
        ..Filter::default()
    };
    // the records name the template as the index of it does
    let target = supports
        .first()
        .and_then(|a| a.record.target.clone())
        .unwrap_or_default();
    let pile = pile(&[(target, template.len())], &records, &filter).remove(0);
    Voted {
        positions: template
            .iter()
            .enumerate()
            .map(|(i, &b)| vote_at(&pile, i, b, filter.default_quality, p.min_support))
            .collect(),
    }
}

impl Voted {
    /// The whole template with its corrections.
    pub fn whole(&self) -> Vec<u8> {
        self.positions
            .iter()
            .flat_map(|(b, _)| b.iter().copied())
            .collect()
    }

    /// The stretches that were corrected, those of at least `shortest` bases.
    pub fn pieces(&self, shortest: usize) -> Vec<Vec<u8>> {
        self.positions
            .split(|(_, corrected)| !corrected)
            .map(|run| {
                run.iter()
                    .flat_map(|(b, _)| b.iter().copied())
                    .collect::<Vec<u8>>()
            })
            .filter(|piece| !piece.is_empty() && piece.len() >= shortest)
            .collect()
    }
}

#[cfg(test)]
mod tests;
