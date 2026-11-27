//! Which reads support the correction of one read (the template), after NECAT's method (Chen et al., *Nature Communications* 2021): the reads
//! that overlap it best first, an error threshold that is the template's own when its best overlaps are good, and no more reads over a stretch
//! that enough reads cover already.
//!
//! * **Candidates**: every other read that the program's mapper chains to the template (the paper chains k-mer pairs it scores by their
//!   distance; this uses minimizer chains and their scores in its place), the best chain of each, from the best score down.
//! * **The threshold**: of the 50 best candidates, those whose alignment covers over 60 % of the template or of the candidate are measured, and
//!   the lowest identity (matching columns over all columns) that an alignment may have is their mean less five standard deviations, and never
//!   under 0.5 (an error of half). The paper's formula for the deviation could not be read, and it does not say what is measured: an identity is
//!   the reading under which its two thresholds behave as it describes.
//! * **Coverage**: a candidate is passed over when every position of its span on the template has more than 12 reads already, and only an
//!   alignment that meets the threshold adds to the count; the 200 best candidates are looked at.

use super::Params;
use crate::seq::minimap::{map_read, Hit, Index, Options};
use crate::seq::variants::reads::Read;
use crate::seq::variants::sam::Record;
use std::collections::HashSet;

/// A candidate read: which one, and where its chain lies on the template.
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    pub read: usize,
    pub score: f64,
    pub from: usize,
    pub to: usize,
}

/// The other reads that chain to the template, the best chain of each, from the best score down (reads that tie in the order they came).
pub fn candidates(
    all: &Index,
    template: usize,
    sequence: &[u8],
    overlap: &Options,
) -> Vec<Candidate> {
    let hits: Vec<Hit> = map_read(all, sequence, overlap, &|t| t == template);
    let mut best: Vec<Candidate> = hits
        .iter()
        .map(|h| Candidate {
            read: h.target,
            score: h.score,
            from: h.qstart,
            to: h.qend,
        })
        .collect();
    best.sort_by(|a, b| b.score.total_cmp(&a.score).then(a.read.cmp(&b.read)));
    let mut seen = HashSet::new();
    best.retain(|c| seen.insert(c.read));
    best
}

/// A candidate's alignment on the template.
#[derive(Clone, Debug)]
pub struct Aligned {
    pub record: Record,
    pub identity: f64,
    /// Where it lies on the template.
    pub from: usize,
    pub to: usize,
    /// Does it cover over the share that makes it one of those the template's threshold is measured on?
    pub long: bool,
}

/// The best alignment of `read` on the template (an index that holds the template alone), if it chains to it. (The mapper gives the best
/// place first.)
pub fn align(template: &Index, read: &Read, mapping: &Options, p: &Params) -> Option<Aligned> {
    let hit = map_read(template, &read.seq, mapping, &|_| false)
        .into_iter()
        .next()?;
    let base = hit.base.as_ref()?;
    let identity = base.matches as f64 / base.block.max(1) as f64;
    let long = (hit.tend - hit.tstart) as f64 > p.long_share * template.length(0) as f64
        || (hit.qend - hit.qstart) as f64 > p.long_share * read.seq.len() as f64;
    Some(Aligned {
        record: Record::from_hit(read, template, &hit)?,
        identity,
        from: hit.tstart,
        to: hit.tend,
        long,
    })
}

/// The lowest identity that an alignment to this template may have: the mean of `measured` less `sigmas` standard deviations, and never
/// below `floor`. With fewer than two to measure, the floor.
pub fn threshold(measured: &[f64], sigmas: f64, floor: f64) -> f64 {
    if measured.len() < 2 {
        return floor;
    }
    let n = measured.len() as f64;
    let mean = measured.iter().sum::<f64>() / n;
    let deviation = (measured.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n).sqrt();
    (mean - sigmas * deviation).max(floor)
}

/// The alignments that support the correction of the template, in the order of their candidates.
pub fn supports(
    template: &Index,
    candidates: &[Candidate],
    reads: &[Read],
    mapping: &Options,
    p: &Params,
) -> Vec<Aligned> {
    // the first of the candidates are aligned once, for the template's threshold, and kept for the choice
    let mut first: Vec<Option<Aligned>> = candidates
        .iter()
        .take(p.measured)
        .map(|c| align(template, &reads[c.read], mapping, p))
        .collect();
    let measured: Vec<f64> = first
        .iter()
        .flatten()
        .filter(|a| a.long)
        .map(|a| a.identity)
        .collect();
    let lowest = threshold(&measured, p.sigmas, p.min_identity);
    let mut covered = vec![0u32; template.length(0)];
    let mut chosen = Vec::new();
    for (i, c) in candidates.iter().enumerate().take(p.max_candidates) {
        let span = c.from.min(covered.len())..c.to.min(covered.len());
        if !span.is_empty() && covered[span].iter().all(|&n| n > p.coverage) {
            continue;
        }
        let aligned = match first.get_mut(i) {
            Some(kept) => kept.take(),
            None => align(template, &reads[c.read], mapping, p),
        };
        let Some(a) = aligned else { continue };
        if a.identity >= lowest {
            for n in &mut covered[a.from..a.to] {
                *n += 1;
            }
            chosen.push(a);
        }
    }
    chosen
}

#[cfg(test)]
mod tests;
