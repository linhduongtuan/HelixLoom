//! A circular replicon turned so that it starts at a start gene. The paper searches the finished replicons for *dnaA* or *repA* with
//! TBLASTN and turns and flips the sequence so that the gene is first, on the forward strand. **TBLASTN is not done here**: the genes are
//! nucleotide sequences that the person gives (no gene is bundled), and are found with the program's mapper, which finds a gene of the same
//! species or genus, not a distant one. The thresholds are this program's.

use crate::seq::dna::revcom;
use crate::seq::minimap::{map_read, Index, Options as MapOptions};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Options {
    /// The share of a gene that its alignment must cover, by the minimizers, to be taken as found.
    pub min_cover: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options { min_cover: 0.3 }
    }
}

/// A mapper's options for a gene that is not quite the one in the replicon.
fn gene_options() -> MapOptions {
    MapOptions {
        k: 11,
        w: 5,
        occ_fraction: 0.0,
        max_gap: 2000,
        bandwidth: 200,
        min_anchors: 3,
        min_score: 20.0,
        ..MapOptions::default()
    }
}

/// `circle` turned to start at the best-matching gene of `genes`, on the other strand if the gene lies on it; none when no gene is found.
pub fn rotate(circle: &[u8], genes: &[(String, Vec<u8>)], o: &Options) -> Option<Vec<u8>> {
    let n = circle.len();
    if n == 0 {
        return None;
    }
    // written twice, so that a gene across the join is whole
    let twice: Vec<u8> = circle.iter().chain(circle).copied().collect();
    let map = gene_options();
    let index = Index::build(vec![("circle".into(), twice)], &map).ok()?;
    let mut best: Option<(f64, usize, bool)> = None;
    for (_, gene) in genes {
        for h in map_read(&index, gene, &map, &|_| false) {
            let share = h.covered as f64 / gene.len().max(1) as f64;
            if h.primary && share >= o.min_cover && best.is_none_or(|b| h.score > b.0) {
                // where the gene begins on the circle: the part of it before the alignment lies before it, or after it on the other strand
                let at = if h.rev {
                    h.tend + h.qstart
                } else {
                    h.tstart + n - h.qstart % n
                };
                best = Some((h.score, at % n, h.rev));
            }
        }
    }
    let (_, at, reverse) = best?;
    let (seq, start) = if reverse {
        (revcom(circle), (n - at) % n)
    } else {
        (circle.to_vec(), at)
    };
    let mut turned = seq[start..].to_vec();
    turned.extend_from_slice(&seq[..start]);
    Some(turned)
}

#[cfg(test)]
mod tests;
