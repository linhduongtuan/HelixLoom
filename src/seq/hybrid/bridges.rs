//! Where long reads join two single-copy segments. Each read is mapped to the segments that are single copy (the program's mapper,
//! `seq::minimap`); two alignments that follow each other along a read, the first reaching the end of its segment and the second starting at
//! the start of its, are a join: the read leaves the first, crosses a gap, and enters the second. The reads that make the same join are
//! one candidate.

use crate::seq::debruijn::graph::{handle, rc, Handle};
use crate::seq::dna::revcom;
use crate::seq::minimap::{map_read, Hit, Index, Options as MapOptions};
use crate::seq::par::par_map;
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Options {
    /// The shortest alignment, in bases of the read, that counts.
    pub min_align: usize,
    /// How far from the end of its segment an alignment may stop (and from the start, one begin), in bases.
    pub end_slack: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            min_align: 300,
            end_slack: 200,
        }
    }
}

/// One read's evidence for a join.
#[derive(Clone, Debug, PartialEq)]
pub struct ReadGap {
    /// The bases between the end of the first segment and the start of the second, by this read (negative: they overlap).
    pub gap: i64,
    /// The read's own bases in the gap, read from the first segment to the second.
    pub bases: Vec<u8>,
    /// The mean of the shares of the two alignments that the minimizers cover.
    pub covered: f64,
}

/// A join from the end of `from` to the start of `to` (as handles: the segment read in the direction it is left or entered), with the reads
/// that show it. The join and its reverse complement are one candidate, kept as the smaller of the two.
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    pub from: Handle,
    pub to: Handle,
    pub reads: Vec<ReadGap>,
}

/// Bases between where the alignment stops and the end of the segment that the read leaves by.
fn tail(index: &Index, h: &Hit) -> usize {
    if h.rev {
        h.tstart
    } else {
        index.length(h.target) - h.tend
    }
}

/// Bases between the start of the segment that the read enters and where the alignment begins.
fn head(index: &Index, h: &Hit) -> usize {
    if h.rev {
        index.length(h.target) - h.tend
    } else {
        h.tstart
    }
}

fn joins_of(
    index: &Index,
    read: &[u8],
    count: &[u32],
    o: &Options,
    map: &MapOptions,
) -> Vec<((Handle, Handle), ReadGap)> {
    let mut hits: Vec<Hit> = map_read(index, read, map, &|t| count[t] != 1)
        .into_iter()
        .filter(|h| h.primary && h.qend - h.qstart >= o.min_align)
        .collect();
    hits.sort_by_key(|h| (h.qstart, h.qend));
    let covered = |h: &Hit| h.covered as f64 / (h.qend - h.qstart).max(1) as f64;
    let mut out = Vec::new();
    for pair in hits.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        let (leave, enter) = (tail(index, a), head(index, b));
        if leave > o.end_slack || enter > o.end_slack {
            continue;
        }
        let gap = b.qstart as i64 - a.qend as i64 - leave as i64 - enter as i64;
        let (from, to) = (handle(a.target, a.rev), handle(b.target, b.rev));
        let (lo, hi) = (a.qend + leave, b.qstart.saturating_sub(enter));
        let bases = if gap > 0 && lo <= hi && hi <= read.len() {
            read[lo..hi].to_vec()
        } else {
            Vec::new()
        };
        let evidence = ReadGap {
            gap,
            bases,
            covered: (covered(a) + covered(b)) / 2.0,
        };
        // the join read from the other side is the same join: keep the smaller form
        if (rc(to), rc(from)) < (from, to) {
            out.push((
                (rc(to), rc(from)),
                ReadGap {
                    bases: revcom(&evidence.bases),
                    ..evidence
                },
            ));
        } else {
            out.push(((from, to), evidence));
        }
    }
    out
}

/// The joins that the reads show, each with its reads, in order of the handles (so that the answer does not depend on the order of the reads).
pub fn candidates(
    index: &Index,
    reads: &[(String, Vec<u8>)],
    count: &[u32],
    o: &Options,
    map: &MapOptions,
) -> Vec<Candidate> {
    let per_read = par_map(reads, |r| joins_of(index, &r.1, count, o, map));
    let mut by: BTreeMap<(Handle, Handle), Vec<ReadGap>> = BTreeMap::new();
    for joins in per_read {
        for (key, evidence) in joins {
            by.entry(key).or_default().push(evidence);
        }
    }
    by.into_iter()
        .map(|((from, to), reads)| Candidate { from, to, reads })
        .collect()
}

#[cfg(test)]
mod tests;
