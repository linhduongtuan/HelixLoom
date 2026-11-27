//! The pileup: for each position of a target, what the aligned reads hold there. Nothing is kept per read: a position keeps, for each of the
//! four bases, how many reads have it and the sum of the evidence that they give it (see `weight`), and for the reads without a base there, how
//! many have a deletion, and what the reads that have an insertion after it inserted.
//!
//! The choices that are this program's: a base counts only when its quality is at least `min_base_quality` (13, as `bcftools mpileup` does by
//! default), its quality is held at `max_base_quality` (30: a base is never trusted more than that, because the errors of one read are not
//! independent of the others'), a read counts only when its mapping quality is at least `min_mapq` (0, as there), and at most `max_depth` reads
//! count at a position (250, as there: the first that come). Base alignment quality (BAQ) is not done.

use super::sam::{uses_read, uses_target, Record};
use std::collections::HashMap;

/// What counts.
#[derive(Clone, Debug, PartialEq)]
pub struct Filter {
    pub min_mapq: u8,
    pub min_base_quality: u8,
    pub max_base_quality: u8,
    pub max_depth: usize,
    /// The quality of every base of a read that has none (a read from FASTA, or a SAM line with `*`).
    pub default_quality: u8,
}

impl Default for Filter {
    fn default() -> Self {
        Filter {
            min_mapq: 0,
            min_base_quality: 13,
            max_base_quality: 30,
            max_depth: 250,
            default_quality: 15,
        }
    }
}

/// The index of a base among A, C, G, T.
pub fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

/// The evidence that a read base of quality `q` gives for its base over each of the three others: the log of the ratio of the chance of
/// seeing it if it is right (`1 - e`) to the chance of seeing it if the base is another one (`e / 3`, the error being any of the three),
/// where `e = 10^(-q/10)`. A column's score for a base is the sum of these over the reads that have it; the reads are taken to err
/// independently of each other. A base is never evidence against its own base: the weight is not below 0 (a base of quality 0 or 1, which is
/// more likely wrong than right, gives nothing).
pub fn weight(q: u8) -> f32 {
    let e = 10f64.powf(-f64::from(q) / 10.0);
    let w = (3.0 * (1.0 - e) / e).ln();
    if w.is_finite() {
        w.max(0.0) as f32
    } else {
        0.0
    }
}

/// The weights by quality, 0 to 93.
fn weights(f: &Filter) -> Vec<f32> {
    (0..=93u8)
        .map(|q| weight(q.min(f.max_base_quality)))
        .collect()
}

#[derive(Clone, Debug, PartialEq)]
pub struct Pileup {
    pub length: usize,
    /// Reads with each of A, C, G, T at the position (with a base of enough quality).
    pub counts: Vec<[u16; 4]>,
    /// The sum of the weights of those reads, for each base.
    pub scores: Vec<[f32; 4]>,
    /// Reads that have a base at the position, whatever its quality, and reads that have a deletion there.
    pub cover: Vec<u16>,
    pub deleted: Vec<u16>,
    /// After a position: the insertions that reads have there, and how many reads have each.
    pub inserted: HashMap<usize, HashMap<Vec<u8>, u16>>,
    /// Reads that were used.
    pub reads: usize,
}

impl Pileup {
    pub fn new(length: usize) -> Pileup {
        Pileup {
            length,
            counts: vec![[0; 4]; length],
            scores: vec![[0.0; 4]; length],
            cover: vec![0; length],
            deleted: vec![0; length],
            inserted: HashMap::new(),
            reads: 0,
        }
    }

    /// The reads that lie over the position, with a base or with a deletion.
    pub fn depth(&self, p: usize) -> usize {
        usize::from(self.cover[p]) + usize::from(self.deleted[p])
    }

    /// The reads that have a base of enough quality at the position.
    pub fn bases(&self, p: usize) -> usize {
        self.counts[p].iter().map(|&n| usize::from(n)).sum()
    }

    fn full(&self, p: usize, f: &Filter) -> bool {
        self.depth(p) >= f.max_depth
    }

    /// Adds the read of a record (one that `usable` says counts and that lies on this target; anything else is passed over).
    fn add(&mut self, r: &Record, f: &Filter, table: &[f32]) -> bool {
        if !r.usable() || r.mapq < f.min_mapq {
            return false;
        }
        let (mut read_at, mut at) = (0usize, r.pos);
        let mut after_base = false;
        for &(op, n) in &r.cigar {
            match op {
                'M' | '=' | 'X' => {
                    for i in 0..n {
                        let (p, i_read) = (at + i, read_at + i);
                        after_base = false;
                        if p >= self.length || i_read >= r.seq.len() || self.full(p, f) {
                            continue;
                        }
                        after_base = true;
                        self.cover[p] = self.cover[p].saturating_add(1);
                        let q = r.qual.get(i_read).copied().unwrap_or(f.default_quality);
                        if let Some(b) = base_index(r.seq[i_read]) {
                            if q >= f.min_base_quality {
                                self.counts[p][b] = self.counts[p][b].saturating_add(1);
                                self.scores[p][b] += table[usize::from(q.min(93))];
                            }
                        }
                    }
                }
                'I' => {
                    if after_base && at > 0 && at - 1 < self.length {
                        if let Some(bases) = r.seq.get(read_at..read_at + n) {
                            let slot = self
                                .inserted
                                .entry(at - 1)
                                .or_default()
                                .entry(bases.to_vec())
                                .or_insert(0);
                            *slot = slot.saturating_add(1);
                        }
                    }
                }
                'D' => {
                    for p in at..(at + n).min(self.length) {
                        if !self.full(p, f) {
                            self.deleted[p] = self.deleted[p].saturating_add(1);
                        }
                    }
                    after_base = false;
                }
                'N' => after_base = false,
                _ => {}
            }
            if uses_read(op) {
                read_at += n;
            }
            if uses_target(op) {
                at += n;
            }
        }
        self.reads += 1;
        true
    }
}

/// A pileup for each target, from the records that lie on it. `targets` are the names and lengths, in the order that the files give.
pub fn pile(targets: &[(String, usize)], records: &[Record], f: &Filter) -> Vec<Pileup> {
    let table = weights(f);
    let mut piles: Vec<Pileup> = targets.iter().map(|(_, n)| Pileup::new(*n)).collect();
    let by_name: HashMap<&str, usize> = targets
        .iter()
        .enumerate()
        .map(|(i, (n, _))| (n.as_str(), i))
        .collect();
    for r in records {
        if let Some(&t) = r.target.as_deref().and_then(|n| by_name.get(n)) {
            piles[t].add(r, f, &table);
        }
    }
    piles
}

#[cfg(test)]
mod tests;
