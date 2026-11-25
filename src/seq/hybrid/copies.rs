//! How many copies of the genome each segment stands for, from its depth and from the graph (the paper's *multiplicity*): a segment of one
//! copy has the depth of the genome as a whole; a repeat has a multiple of it, and more ways in and out.
//!
//! The paper's steps are kept: the depth that a single copy has (the median of the graph's), seeds that are single copy (that depth, and
//! at most one way on at either end), a greedy spreading through the links where the numbers add up, and, when that stops, the largest
//! segment that could be single copy set to one. The paper gives no tolerance for "near" the depth; the ones here are this program's.
//! The depth of a plasmid that has another number of copies than the chromosome does not fit (the paper says so too): its segments come
//! out as repeats unless the depth of a single copy is given.

use crate::seq::debruijn::graph::{handle, seg, Graph};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Options {
    /// The depth of a single copy; none: the median depth of the bases of the graph.
    pub single_depth: Option<f64>,
    /// How far from one a seed's depth (as a share of a single copy's) may be.
    pub tolerance: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            single_depth: None,
            tolerance: 0.25,
        }
    }
}

/// How far from a whole number of copies a depth may be, as a share of that number, for the links to be believed.
const SPREAD_TOLERANCE: f64 = 0.4;
/// How far from a single copy a segment may be, as a share, to be taken as one when nothing else decides.
const STALL_TOLERANCE: f64 = 0.5;

#[derive(Clone, Debug, PartialEq)]
pub struct Copies {
    pub single_depth: f64,
    /// The copies of each segment, by its number; 0 for one that is not there.
    pub count: Vec<u32>,
    /// How many were set by their depth alone, nothing in the links agreeing.
    pub by_depth: usize,
}

/// The depth that half of the bases of the graph are at or below.
fn median_depth(g: &Graph) -> f64 {
    let mut by: Vec<(f64, usize)> = g
        .ids()
        .into_iter()
        .map(|i| (g.coverage(i), g.length(i)))
        .collect();
    by.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: usize = by.iter().map(|b| b.1).sum();
    let mut seen = 0;
    by.iter()
        .find(|b| {
            seen += b.1;
            seen * 2 >= total
        })
        .map_or(1.0, |b| b.0)
}

pub fn estimate(g: &Graph, o: &Options) -> Copies {
    let slots = g.handles() / 2;
    let single = o
        .single_depth
        .unwrap_or_else(|| median_depth(g))
        .max(f64::MIN_POSITIVE);
    let ratio = |i: usize| g.coverage(i) / single;
    let mut count = vec![0u32; slots];
    // seeds: the depth of one copy, and no more than one way on at either end
    for i in g.ids() {
        let (f, r) = (handle(i, false), handle(i, true));
        if (ratio(i) - 1.0).abs() <= o.tolerance && g.outdeg(f) <= 1 && g.outdeg(r) <= 1 {
            count[i] = 1;
        }
    }
    loop {
        let mut changed = false;
        for u in g.ids() {
            for h in [handle(u, false), handle(u, true)] {
                // the ways on from this end, when each of them is reached from here only: the copies that leave are the copies that arrive
                let ways = g.succ(h);
                if ways.is_empty() || ways.iter().any(|&x| g.indeg(x) != 1) {
                    continue;
                }
                let unknown: Vec<usize> = ways
                    .iter()
                    .map(|&x| seg(x))
                    .filter(|&s| count[s] == 0)
                    .collect();
                let known: u32 = ways.iter().map(|&x| count[seg(x)]).sum();
                let (target, value) = match (count[u], unknown.len()) {
                    (0, 0) => (u, known),
                    (m, 1) if m > known => (unknown[0], m - known),
                    _ => continue,
                };
                if value >= 1
                    && (ratio(target) - f64::from(value)).abs()
                        <= SPREAD_TOLERANCE * f64::from(value)
                    && count[target] == 0
                {
                    count[target] = value;
                    changed = true;
                }
            }
        }
        if changed {
            continue;
        }
        // nothing more follows from the links: the largest segment that could be single copy is
        let largest = g
            .ids()
            .into_iter()
            .filter(|&i| count[i] == 0 && (ratio(i) - 1.0).abs() <= STALL_TOLERANCE)
            .max_by_key(|&i| (g.length(i), std::cmp::Reverse(i)));
        match largest {
            Some(i) => count[i] = 1,
            None => break,
        }
    }
    let mut by_depth = 0;
    for i in g.ids() {
        if count[i] == 0 {
            count[i] = ratio(i).round().max(1.0) as u32;
            by_depth += 1;
        }
    }
    Copies {
        single_depth: single,
        count,
        by_depth,
    }
}

#[cfg(test)]
mod tests;
