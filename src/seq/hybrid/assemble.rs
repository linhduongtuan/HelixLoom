//! The steps in order: the copies of each segment, the joins that long reads show, the bridges that they make and how good each is, the
//! best put into the graph, and what is then unbranched made one segment.

use super::apply::{apply, Bridge, Mode, Via};
use super::bridges::{candidates, Candidate, Options as BridgeOptions};
use super::copies::{estimate, Options as CopyOptions};
use super::quality::{quality, Facts};
use super::route::{routes, tolerance};
use crate::seq::debruijn::graph::{seg, Graph};
use crate::seq::minimap::{Index, Options as MapOptions};

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub mode: Mode,
    /// The lowest quality of a bridge; none: the mode's.
    pub min_quality: Option<f64>,
    /// The fewest reads that make a bridge.
    pub min_reads: usize,
    pub copies: CopyOptions,
    pub bridges: BridgeOptions,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            mode: Mode::Normal,
            min_quality: None,
            min_reads: 1,
            copies: CopyOptions::default(),
            bridges: BridgeOptions::default(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    pub single_depth: f64,
    pub segments: usize,
    pub single_copy: usize,
    pub repeats: usize,
    /// Segments whose number of copies came from their depth alone.
    pub by_depth: usize,
    pub joins_seen: usize,
    pub bridges_found: usize,
    pub bridges_good_enough: usize,
    pub bridges_made: usize,
    pub segments_after: usize,
    pub circular: usize,
}

pub struct Outcome {
    pub graph: Graph,
    pub report: Report,
}

fn median(mut v: Vec<i64>) -> i64 {
    v.sort_unstable();
    v[v.len() / 2]
}

/// The bridge that a candidate makes: through a route of repeats if one fits the reads' gap, else from the reads' own bases.
fn resolve(g: &Graph, count: &[u32], c: &Candidate) -> Option<Bridge> {
    let gaps: Vec<i64> = c.reads.iter().map(|r| r.gap).collect();
    let gap = median(gaps.clone());
    let (a, b) = (seg(c.from), seg(c.to));
    let covered = c.reads.iter().map(|r| r.covered).sum::<f64>() / c.reads.len() as f64;
    let mut facts = Facts {
        reads: c.reads.len(),
        deviation: 0.0,
        tolerance: tolerance(gap) as f64,
        lengths: (g.length(a), g.length(b)),
        covered,
        depths: (g.coverage(a), g.coverage(b)),
        through_graph: true,
    };
    let route = routes(g, count, c.from, c.to, gap).into_iter().next();
    let via = match route {
        Some(r) => {
            facts.deviation = (r.distance - gap).abs() as f64;
            if r.via.is_empty() {
                Via::Link
            } else {
                Via::Route(r.via)
            }
        }
        None => {
            // no route: the bases of the read whose gap is the median one
            let best = c
                .reads
                .iter()
                .min_by_key(|r| ((r.gap - gap).abs(), r.gap))?;
            if gap <= 0 || best.bases.is_empty() {
                return None;
            }
            let (lo, hi) = (gaps.iter().min()?, gaps.iter().max()?);
            facts.deviation = (hi - lo) as f64 / 2.0;
            facts.through_graph = false;
            Via::Reads(best.bases.clone())
        }
    };
    Some(Bridge {
        from: c.from,
        to: c.to,
        via,
        quality: quality(&facts),
        reads: c.reads.len(),
    })
}

/// Join the single-copy segments of `graph` with the long `reads`.
pub fn run(graph: Graph, reads: &[(String, Vec<u8>)], o: &Options) -> Result<Outcome, String> {
    let copies = estimate(&graph, &o.copies);
    let targets: Vec<(String, Vec<u8>)> = (0..graph.handles() / 2)
        .map(|i| (i.to_string(), graph.segment(i).seq.clone()))
        .collect();
    let map = MapOptions::default();
    let index = Index::build(targets, &map)?;
    let found = candidates(&index, reads, &copies.count, &o.bridges, &map);
    let min_quality = o.min_quality.unwrap_or_else(|| o.mode.min_quality());
    let mut bridges: Vec<Bridge> = found
        .iter()
        .filter(|c| c.reads.len() >= o.min_reads)
        .filter_map(|c| resolve(&graph, &copies.count, c))
        .collect();
    let bridges_found = bridges.len();
    bridges.retain(|b| b.quality >= min_quality);
    bridges.sort_by(|x, y| {
        y.quality
            .total_cmp(&x.quality)
            .then((x.from, x.to).cmp(&(y.from, y.to)))
    });
    let single_copy = graph
        .ids()
        .into_iter()
        .filter(|&i| copies.count[i] == 1)
        .count();
    let segments = graph.count();
    let applied = apply(graph, &copies, &bridges, o.mode);
    let graph = applied.graph;
    Ok(Outcome {
        report: Report {
            single_depth: copies.single_depth,
            segments,
            single_copy,
            repeats: segments - single_copy,
            by_depth: copies.by_depth,
            joins_seen: found.len(),
            bridges_found,
            bridges_good_enough: bridges.len(),
            bridges_made: applied.made,
            segments_after: graph.count(),
            circular: graph
                .ids()
                .into_iter()
                .filter(|&i| graph.circular(i))
                .count(),
        },
        graph,
    })
}

#[cfg(test)]
mod tests;
