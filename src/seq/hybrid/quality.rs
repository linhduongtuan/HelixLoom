//! The score of a bridge, from 0 to 100. The paper says what a bridge's quality depends on and not how: for a bridge made from long reads, the
//! number of reads, how well the reads' bases agree with the graph's path, the lengths of the two segments, the lengths and identities of the
//! reads' alignments to them, and how alike their depths are; a bridge made from the long reads' own bases (no path in the graph) is
//! rated lower. **The formula here is this program's own**: the product of one factor for each, each between 0 and 1, times 100.
//! The thresholds of the modes (25, 10, 1) are the paper's numbers, on this other scale.

/// A segment of this many bases or more counts in full; a shorter one in proportion.
pub const FULL_LENGTH: f64 = 1000.0;
/// The share of an alignment that the minimizers cover at which its quality counts in full.
pub const FULL_COVER: f64 = 0.4;
/// What a bridge of the reads' own bases, not of a path through the graph, is multiplied by.
pub const READS_ONLY: f64 = 0.5;

#[derive(Clone, Debug, PartialEq)]
pub struct Facts {
    /// Reads that support the bridge.
    pub reads: usize,
    /// How far the reads' gaps and the path's distance differ, in bases, and how much that may be (`route::tolerance`).
    pub deviation: f64,
    pub tolerance: f64,
    /// The bases of the two segments.
    pub lengths: (usize, usize),
    /// The mean of the shares of each alignment that the minimizers cover, 0 to 1.
    pub covered: f64,
    /// The depths of the two segments.
    pub depths: (f64, f64),
    /// The bridge is a path through the graph, and not the reads' bases.
    pub through_graph: bool,
}

pub fn quality(f: &Facts) -> f64 {
    let n = f.reads as f64;
    let support = n / (n + 1.0);
    let agreement = 1.0 - 0.5 * (f.deviation / f.tolerance.max(1.0)).clamp(0.0, 1.0);
    let length =
        (f.lengths.0 as f64 / FULL_LENGTH).min(1.0) * (f.lengths.1 as f64 / FULL_LENGTH).min(1.0);
    let alignment = (f.covered / FULL_COVER).clamp(0.0, 1.0);
    let (lo, hi) = (f.depths.0.min(f.depths.1), f.depths.0.max(f.depths.1));
    let depth = if hi > 0.0 { lo / hi } else { 1.0 };
    let origin = if f.through_graph { 1.0 } else { READS_ONLY };
    100.0 * support * agreement * length * alignment * depth * origin
}

#[cfg(test)]
mod tests;
