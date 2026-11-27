//! Calling the differences between the reads and a reference, and the consensus, for a sample with one copy of the sequence (a plasmid, a
//! bacterial chromosome, a haploid genome). **The model is this program's own**: the paper that the approach follows (Li, *Bioinformatics*
//! 2011, 27:2987: a likelihood of the bases at a site from their qualities, errors taken as independent, and a prior over the alleles) gave its
//! formulas as images that could not be read here, so what is written below is the simplest model of that kind, stated in full.
//!
//! **A substitution.** At a position the reads show bases with Phred qualities. If the sample has the base `a`, a read base `b` of error
//! chance `e` is seen with chance `1 - e` when `b = a`, and `e / 3` otherwise; the reads are independent. So the log-likelihood of `a` is a
//! constant plus the sum, over the reads that show `a`, of `ln(3 (1 - e) / e)` (`pileup::weight`). The prior of the reference base is
//! `1 - π` and each of the three others has `π / 3` (`π`, `snp_prior`, 0.001). The posterior of the reference base is the normalised product;
//! a site is a variant when another base has the best posterior, and its quality is `-10 log10` of the posterior of the reference base,
//! the Phred scale of the chance that the site is not a variant, as a VCF's QUAL is (held at 999).
//!
//! **An insertion or a deletion.** The reads do not give a quality for these. A deletion at a position is seen by the reads that have a gap
//! there, an insertion after a position by the reads that have the same inserted bases there. With `n` of `N` reads showing it, and a chance
//! `ε` (`indel_error`, 0.1: nanopore reads err this often in a place) that a read shows an indel that is not there, or misses one that is:
//! the log-likelihood ratio of "the sample has it" to "the sample has not" is `(2n - N) ln((1 - ε) / ε)`, and the prior of an indel is `ρ`
//! (`indel_prior`, 0.0001). A run of positions that are deleted is one deletion; its quality, depth and counts are those of its weakest
//! position. Of the insertions after a position the commonest is taken (the reads that insert something else are counted against it): an
//! insertion that reads write in different ways (a nanopore read's errors inside it) is split among them and is found less often.
//!
//! **The consensus** is the reference with the calls made (quality at least `min_quality`, 20) applied: a base that is called, a deletion
//! that is called leaving out its positions, an insertion that is called put after its position, and `N` where fewer than `min_depth` (5)
//! reads have a base of enough quality and no deletion was called. Where the reference has no known base (N), the base the reads give the
//! most evidence for is written when it is ahead of the next by `min_quality` (as a likelihood ratio), and `N` otherwise; no VCF line. **Not done:** diploid or mixed samples, the dependence of errors between
//! reads at one site (which makes these qualities too confident at high depth: the cap on base quality eases it), base alignment quality,
//! the bases at the ends of a read beyond the first and last anchor of the mapping, and phasing neighbours into one event.

use super::pileup::{base_index, Pileup};
use std::collections::HashMap;

/// The highest quality that is written.
pub const MAX_QUALITY: f64 = 999.0;

#[derive(Clone, Debug, PartialEq)]
pub struct Params {
    /// Fewest reads with a base to call a site or to write a base of the consensus; fewest reads over a position to call an indel there.
    pub min_depth: usize,
    /// Lowest quality of a call that is made.
    pub min_quality: f64,
    pub snp_prior: f64,
    pub indel_prior: f64,
    pub indel_error: f64,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            min_depth: 5,
            min_quality: 20.0,
            snp_prior: 0.001,
            indel_prior: 0.0001,
            indel_error: 0.1,
        }
    }
}

impl Params {
    pub fn check(&self) -> Result<(), String> {
        let open = |x: f64| x > 0.0 && x < 1.0;
        if self.min_depth == 0 {
            return Err("The least depth must be at least 1.".into());
        }
        if !self.min_quality.is_finite() || self.min_quality < 0.0 {
            return Err("The least quality must be a number that is not negative.".into());
        }
        if !open(self.snp_prior) || !open(self.indel_prior) {
            return Err("The priors must be between 0 and 1, neither included.".into());
        }
        if !(self.indel_error > 0.0 && self.indel_error < 0.5) {
            return Err("The chance of a wrong indel must be above 0 and under 0.5.".into());
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Kind {
    Snp,
    Insertion,
    Deletion,
}

/// One difference. `pos` is where `reference` begins on the target, from 0 (a file counts from 1). For an insertion or a deletion, the first
/// base of both `reference` and `alt` is the base before it, as VCF writes them.
#[derive(Clone, Debug, PartialEq)]
pub struct Variant {
    pub kind: Kind,
    pub pos: usize,
    pub reference: Vec<u8>,
    pub alt: Vec<u8>,
    pub quality: f64,
    /// The reads at the site, those that have the reference, and those that have the alternate.
    pub depth: usize,
    pub ref_reads: usize,
    pub alt_reads: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Calls {
    pub variants: Vec<Variant>,
    pub consensus: Vec<u8>,
    /// The reads over each position of the reference.
    pub depth: Vec<u32>,
    /// Positions of the reference that the consensus writes as `N`.
    pub unknown: usize,
    pub mean_depth: f64,
}

/// `-10 log10` of the chance that a site is the reference when the two hypotheses have these log-weights: the quality of "not the reference".
fn phred(log_alt: f64, log_ref: f64) -> f64 {
    let top = log_alt.max(log_ref);
    let total = top + ((log_alt - top).exp() + (log_ref - top).exp()).ln();
    ((total - log_ref) * 10.0 / std::f64::consts::LN_10).clamp(0.0, MAX_QUALITY)
}

/// The best of the four bases at a position and the quality of its being not the reference; `None` when it is the reference.
fn substitution(scores: &[f32; 4], reference: usize, prior: f64) -> Option<(usize, f64)> {
    let log_weight = |a: usize| {
        let p = if a == reference {
            1.0 - prior
        } else {
            prior / 3.0
        };
        p.ln() + f64::from(scores[a])
    };
    let logs: Vec<f64> = (0..4).map(log_weight).collect();
    let best = (0..4).fold(reference, |b, a| if logs[a] > logs[b] { a } else { b });
    if best == reference {
        return None;
    }
    let top = logs[best];
    let total: f64 = logs.iter().map(|l| (l - top).exp()).sum();
    let p_reference = (logs[reference] - top).exp() / total;
    Some((best, (-10.0 * p_reference.log10()).clamp(0.0, MAX_QUALITY)))
}

/// The quality of an indel that `alt` of `total` reads show, when it is called.
fn indel(alt: usize, total: usize, p: &Params) -> Option<f64> {
    if alt == 0 || total < p.min_depth {
        return None;
    }
    let step = ((1.0 - p.indel_error) / p.indel_error).ln();
    let ratio = (2.0 * alt as f64 - total as f64) * step;
    let (log_alt, log_ref) = (p.indel_prior.ln() + ratio, (1.0 - p.indel_prior).ln());
    if log_alt <= log_ref {
        return None;
    }
    let q = phred(log_alt, log_ref);
    (q >= p.min_quality).then_some(q)
}

/// A run of positions that are deleted, with the numbers of its weakest.
struct Run {
    first: usize,
    last: usize,
    quality: f64,
    depth: usize,
    alt: usize,
}

fn deletions(pile: &Pileup, n: usize, p: &Params) -> Vec<Run> {
    let mut runs: Vec<Run> = Vec::new();
    for at in 0..n {
        let (alt, depth) = (usize::from(pile.deleted[at]), pile.depth(at));
        let Some(quality) = indel(alt, depth, p) else {
            continue;
        };
        match runs.last_mut() {
            Some(run) if run.last + 1 == at => {
                run.last = at;
                if quality < run.quality {
                    (run.quality, run.depth, run.alt) = (quality, depth, alt);
                }
            }
            _ => runs.push(Run {
                first: at,
                last: at,
                quality,
                depth,
                alt,
            }),
        }
    }
    runs
}

/// An insertion that is called: the position that it follows, what it inserts, its quality, the reads over the position and the reads
/// that have it.
struct Inserted {
    at: usize,
    bases: Vec<u8>,
    quality: f64,
    depth: usize,
    count: usize,
}

/// The insertions that are called, by the position that they follow.
fn insertions(pile: &Pileup, deleted: &[bool], p: &Params) -> Vec<Inserted> {
    let mut places: Vec<&usize> = pile.inserted.keys().collect();
    places.sort();
    let mut out = Vec::new();
    for &at in places {
        if at >= deleted.len() || deleted[at] {
            continue;
        }
        // an insertion that is called has more than half the reads, so two that tie are never called, and which is taken does not matter
        let best = pile.inserted[&at]
            .iter()
            .filter(|(bases, _)| bases.iter().all(|&b| base_index(b).is_some()))
            .max_by_key(|(_, &count)| count);
        if let Some((bases, &count)) = best {
            let depth = usize::from(pile.cover[at]);
            if let Some(q) = indel(usize::from(count), depth, p) {
                out.push(Inserted {
                    at,
                    bases: bases.clone(),
                    quality: q,
                    depth,
                    count: usize::from(count),
                });
            }
        }
    }
    out
}

/// The calls and the consensus for one target: `reference` is its sequence, `pile` the pileup of the reads on it. The two must be of the same
/// length.
pub fn call_target(reference: &[u8], pile: &Pileup, p: &Params) -> Calls {
    let reference: Vec<u8> = reference.to_ascii_uppercase();
    let n = reference.len().min(pile.length);
    let mut variants = Vec::new();
    let mut deleted = vec![false; n];
    for run in deletions(pile, n, p) {
        let (pos, end) = if run.first > 0 {
            (run.first - 1, run.last + 1)
        } else if run.last + 1 < n {
            (0, run.last + 2)
        } else {
            continue;
        };
        for flag in &mut deleted[run.first..=run.last] {
            *flag = true;
        }
        variants.push(Variant {
            kind: Kind::Deletion,
            pos,
            reference: reference[pos..end].to_vec(),
            alt: vec![reference[if run.first > 0 { pos } else { end - 1 }]],
            quality: run.quality,
            depth: run.depth,
            ref_reads: run.depth - run.alt,
            alt_reads: run.alt,
        });
    }
    let mut inserted: HashMap<usize, Vec<u8>> = HashMap::new();
    for ins in insertions(pile, &deleted, p) {
        let mut alt = vec![reference[ins.at]];
        alt.extend_from_slice(&ins.bases);
        variants.push(Variant {
            kind: Kind::Insertion,
            pos: ins.at,
            reference: vec![reference[ins.at]],
            alt,
            quality: ins.quality,
            depth: ins.depth,
            ref_reads: ins.depth.saturating_sub(ins.count),
            alt_reads: ins.count,
        });
        inserted.insert(ins.at, ins.bases);
    }
    let mut consensus = Vec::with_capacity(n);
    let mut unknown = 0;
    for at in 0..n {
        if deleted[at] {
            continue;
        }
        if pile.bases(at) < p.min_depth {
            consensus.push(b'N');
            unknown += 1;
        } else if let Some(r) = base_index(reference[at]) {
            match substitution(&pile.scores[at], r, p.snp_prior) {
                Some((best, quality)) if quality >= p.min_quality => {
                    consensus.push(b"ACGT"[best]);
                    variants.push(Variant {
                        kind: Kind::Snp,
                        pos: at,
                        reference: vec![reference[at]],
                        alt: vec![b"ACGT"[best]],
                        quality,
                        depth: pile.bases(at),
                        ref_reads: usize::from(pile.counts[at][r]),
                        alt_reads: usize::from(pile.counts[at][best]),
                    });
                }
                _ => consensus.push(reference[at]),
            }
        } else {
            // the reference has no known base here: the base that the reads give the most evidence for, when it is clearly ahead of the
            // next (by the lowest quality of a call), and otherwise N
            let mut order = [0usize, 1, 2, 3];
            order.sort_by(|&a, &b| pile.scores[at][b].total_cmp(&pile.scores[at][a]));
            let lead = f64::from(pile.scores[at][order[0]] - pile.scores[at][order[1]]);
            if lead * 10.0 / std::f64::consts::LN_10 >= p.min_quality {
                consensus.push(b"ACGT"[order[0]]);
            } else {
                consensus.push(b'N');
                unknown += 1;
            }
        }
        if let Some(bases) = inserted.get(&at) {
            consensus.extend_from_slice(bases);
        }
    }
    variants.sort_by_key(|v| (v.pos, v.kind));
    let depth: Vec<u32> = (0..n).map(|at| pile.depth(at) as u32).collect();
    let mean_depth = if n == 0 {
        0.0
    } else {
        depth.iter().map(|&d| f64::from(d)).sum::<f64>() / n as f64
    };
    Calls {
        variants,
        consensus,
        depth,
        unknown,
        mean_depth,
    }
}

#[cfg(test)]
mod tests;
