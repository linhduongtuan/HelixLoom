//! Mapping one read to an index: its minimizers are looked up, the places they share with each target are chained, and each chain is a
//! mapping. Mappings that overlap on the read are one place with a best mapping (primary) and weaker ones (secondary); the quality of the
//! primary says how much better it is than the best secondary (the formula of minimap2's paper).

use super::base::{along, Base};
use super::chain::{chain, Anchor, Params};
use super::index::Index;
use super::sketch::sketch;
use super::Options;
use crate::seq::dna::revcom;
use std::collections::HashMap;

#[derive(Clone, Debug, PartialEq)]
pub struct Hit {
    pub target: usize,
    /// The read lies on the other strand of the target.
    pub rev: bool,
    /// Where the mapping is on the read (as the read is written) and on the target, the end not included.
    pub qstart: usize,
    pub qend: usize,
    pub tstart: usize,
    pub tend: usize,
    pub score: f64,
    /// The anchors of the chain, on the read as it lies against the target.
    pub anchors: Vec<Anchor>,
    pub primary: bool,
    /// For a primary: the score of its best secondary; for a secondary: 0.
    pub best_secondary: f64,
    pub mapq: u8,
    /// The divergence that the minimizers of the chain suggest, `ln(n / m) / k` (the paper's estimate): none when it cannot be told.
    pub divergence: Option<f64>,
    /// The bases that the anchors cover on the read.
    pub covered: usize,
    pub base: Option<Base>,
}

/// minimap2's mapping quality, `40 (1 - f2/f1) min(1, m/10) ln f1`, from the score `f1` of the best chain, `f2` of the second best, and the
/// number `m` of anchors of the best, held to 0 to 60.
pub fn mapq(f1: f64, f2: f64, m: usize) -> u8 {
    if f1 <= 1.0 {
        return 0;
    }
    let q = 40.0 * (1.0 - (f2 / f1).min(1.0)) * (m as f64 / 10.0).min(1.0) * f1.ln();
    q.clamp(0.0, 60.0).round() as u8
}

/// The bases of the read that a set of anchors covers: the union of their k-mers.
fn covered(anchors: &[Anchor], k: usize) -> usize {
    let mut ys: Vec<usize> = anchors.iter().map(|a| a.y).collect();
    ys.sort_unstable();
    let (mut total, mut end) = (0, 0);
    for y in ys {
        total += (y + k).saturating_sub(end.max(y));
        end = end.max(y + k);
    }
    total
}

fn overlap(a: (usize, usize), b: (usize, usize)) -> usize {
    a.1.min(b.1).saturating_sub(a.0.max(b.0))
}

/// The mappings of the read `q` to `index`, the primary of each place first and then its secondaries. A target for which `skip` is true is
/// not looked at. With `o.overlaps` every chain is reported and none is labelled secondary.
pub fn map_read(index: &Index, q: &[u8], o: &Options, skip: &dyn Fn(usize) -> bool) -> Vec<Hit> {
    let k = index.k;
    let qlen = q.len();
    let minimizers = sketch(q, k, index.w);
    // the anchors by target and strand, the read counted along the target
    let mut groups: HashMap<(usize, bool), Vec<Anchor>> = HashMap::new();
    for m in &minimizers {
        for s in index
            .hits(m.hash)
            .iter()
            .filter(|s| !skip(s.target as usize))
        {
            let reverse = m.rev != s.rev;
            let y = if reverse { qlen - k - m.pos } else { m.pos };
            groups
                .entry((s.target as usize, reverse))
                .or_default()
                .push(Anchor {
                    x: s.pos as usize,
                    y,
                });
        }
    }
    let params = Params {
        k,
        max_gap: o.max_gap,
        bandwidth: o.bandwidth,
        lookback: o.lookback,
        min_anchors: o.min_anchors,
        min_score: o.min_score,
    };
    let mut hits: Vec<Hit> = Vec::new();
    let mut keys: Vec<&(usize, bool)> = groups.keys().collect();
    keys.sort();
    for &(target, rev) in keys {
        for c in chain(&groups[&(target, rev)], &params) {
            // the anchors of a chain rise along the target and along the read
            let (first, last) = (c.anchors[0], c.anchors[c.anchors.len() - 1]);
            let (ys, ye) = (first.y, last.y + k);
            let (qstart, qend) = if rev {
                (qlen - ye, qlen - ys)
            } else {
                (ys, ye)
            };
            // the minimizers of the read between the ends of the chain, against those that are in it
            let seeds = minimizers
                .iter()
                .filter(|m| {
                    let y = if rev { qlen - k - m.pos } else { m.pos };
                    y >= ys && y + k <= ye
                })
                .count();
            let divergence = (seeds >= c.anchors.len() && !c.anchors.is_empty())
                .then(|| (seeds as f64 / c.anchors.len() as f64).ln() / k as f64);
            hits.push(Hit {
                target,
                rev,
                qstart,
                qend,
                tstart: first.x,
                tend: last.x + k,
                score: c.score,
                covered: covered(&c.anchors, k),
                anchors: c.anchors,
                primary: true,
                best_secondary: 0.0,
                mapq: 255,
                divergence,
                base: None,
            });
        }
    }
    if !o.overlaps {
        hits = label(hits, o);
    }
    if o.cigar {
        for h in &mut hits {
            let oriented = if h.rev {
                revcom(q)
            } else {
                q.to_ascii_uppercase()
            };
            h.base = along(&h.anchors, &oriented, index.sequence(h.target), k);
        }
    }
    hits
}

/// Primary and secondary mappings: from the best score down, a mapping that overlaps a primary on the read by at least `mask_level` of the
/// shorter of the two is its secondary; the quality of a primary comes from its best secondary; secondaries below the ratio or over the
/// number asked for are left out.
fn label(mut hits: Vec<Hit>, o: &Options) -> Vec<Hit> {
    hits.sort_by(|a, b| {
        b.score
            .total_cmp(&a.score)
            .then((a.target, a.qstart).cmp(&(b.target, b.qstart)))
    });
    let mut primaries: Vec<usize> = Vec::new();
    let mut secondaries: Vec<Vec<usize>> = Vec::new();
    for (i, h) in hits.iter().enumerate() {
        let span = (h.qstart, h.qend);
        let parent = primaries.iter().position(|&p| {
            let other = (hits[p].qstart, hits[p].qend);
            let shorter = (span.1 - span.0).min(other.1 - other.0);
            overlap(span, other) as f64 >= o.mask_level * shorter as f64
        });
        match parent {
            Some(p) => secondaries[p].push(i),
            None => {
                primaries.push(i);
                secondaries.push(Vec::new());
            }
        }
    }
    let mut out = Vec::with_capacity(hits.len());
    for (p, &i) in primaries.iter().enumerate() {
        let mut primary = hits[i].clone();
        let second = secondaries[p]
            .iter()
            .map(|&s| hits[s].score)
            .fold(0.0, f64::max);
        primary.best_secondary = second;
        primary.mapq = mapq(primary.score, second, primary.anchors.len());
        let keep: Vec<Hit> = secondaries[p]
            .iter()
            .filter(|&&s| hits[s].score >= o.secondary_ratio * primary.score)
            .take(o.max_secondary)
            .map(|&s| Hit {
                primary: false,
                mapq: 0,
                best_secondary: 0.0,
                ..hits[s].clone()
            })
            .collect();
        out.push(primary);
        out.extend(keep);
    }
    out
}

#[cfg(test)]
mod tests;
