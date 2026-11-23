//! Where the alleles of a locus lie in the contigs. The program's mapper (`seq::minimap`) finds the places where an allele shares
//! minimizers with a contig, on either strand; the alleles that land on one place are the candidates for it; the best few are aligned
//! over the whole of the place (the program's aligner, with the overhanging ends free, and a margin of contig on each side), and the
//! stretch of contig that the alignment covers is what is found there. An end of the allele that has no contig under it (the contig
//! ends) is left out of the stretch; so may an end base that does not match, which the free overhang makes possible (an allele with a changed
//! base at its very end can then come out as a partial match, not a novel one).
//!
//! The found stretch is looked up among the alleles: when it is one of them, the call there is exact. Otherwise the closest allele and how
//! much of it was found and how alike it is say whether it is novel (all of the allele, alike enough), partial (some of it) or nothing.
//! A base `N` in a contig is a mismatch, as any other base that is not the allele's.

use super::scheme::Locus;
use crate::seq::align::{align, Scoring};
use crate::seq::dna::revcom;
use crate::seq::minimap::{map_read, Hit, Index, Options};
use std::collections::HashMap;

/// Contig bases to take on each side of a place, so that the alignment of the whole allele has room.
const MARGIN: usize = 40;
/// Alleles aligned at a place: the best few by how many minimizers they share with it (alleles of a locus are alike, and all of them
/// are not worth a table each).
const CANDIDATES: usize = 8;

/// The contigs of a genome, with the index that the search needs.
pub struct Contigs {
    index: Index,
}

/// The options of the mapper for short, alike sequences: 11-mers in windows of 3, a chain of two anchors will do, and every chain counts.
fn options() -> Options {
    Options {
        k: 11,
        w: 3,
        occ_fraction: 0.0,
        max_gap: 1000,
        bandwidth: 100,
        min_anchors: 2,
        min_score: 15.0,
        overlaps: true,
        ..Options::default()
    }
}

impl Contigs {
    pub fn new(contigs: Vec<(String, Vec<u8>)>) -> Result<Contigs, String> {
        Ok(Contigs {
            index: Index::build(contigs, &options())?,
        })
    }
}

/// What was found at one place for a locus.
#[derive(Clone, Debug, PartialEq)]
pub struct Found {
    /// The allele closest to what is there.
    pub allele: String,
    /// When the stretch found is exactly an allele: its id (the closest one is then the same).
    pub exact: Option<String>,
    /// The bases of the allele that were aligned, of all of them (`allele_len`).
    pub covered: usize,
    pub allele_len: usize,
    /// Columns that are the same, of the columns between the first and the last aligned base.
    pub matches: usize,
    pub columns: usize,
    /// The contig bases that the alignment covers, written as the allele is (the other strand when the allele lies on the reverse strand).
    pub sequence: Vec<u8>,
}

impl Found {
    pub fn identity(&self) -> f64 {
        if self.columns == 0 {
            0.0
        } else {
            self.matches as f64 / self.columns as f64
        }
    }

    pub fn coverage(&self) -> f64 {
        self.covered as f64 / self.allele_len.max(1) as f64
    }
}

/// A place on a contig, on one strand, that some alleles are alike at.
struct Site {
    contig: usize,
    rev: bool,
    start: usize,
    end: usize,
    candidates: Vec<(usize, f64)>,
}

/// The places on the contigs that the alleles of `locus` land on.
fn sites(index: &Index, locus: &Locus) -> Vec<Site> {
    let mut hits: Vec<(usize, Hit)> = Vec::new();
    for (a, allele) in locus.alleles.iter().enumerate() {
        hits.extend(
            map_read(index, &allele.seq, &options(), &|_| false)
                .into_iter()
                .map(|h| (a, h)),
        );
    }
    hits.sort_by_key(|(_, h)| (h.target, h.rev, h.tstart));
    let mut sites: Vec<Site> = Vec::new();
    for (a, h) in hits {
        match sites.last_mut() {
            Some(s) if s.contig == h.target && s.rev == h.rev && h.tstart < s.end => {
                s.end = s.end.max(h.tend);
                s.candidates.push((a, h.score));
            }
            _ => sites.push(Site {
                contig: h.target,
                rev: h.rev,
                start: h.tstart,
                end: h.tend,
                candidates: vec![(a, h.score)],
            }),
        }
    }
    sites
}

/// An allele aligned over a window of contig, the ends that overhang being free: how much of the allele and how many columns are in the aligned stretch, and
/// the contig bases of it. The stretch is from the first column where both have a base to the last; an end of the allele that has no
/// contig under it (the contig ends) is not in it, and does not count as covered.
fn aligned_over(allele: &[u8], window: &[u8]) -> Option<(usize, usize, usize, Vec<u8>)> {
    let al = align(allele, window, &Scoring::default(), true).ok()?;
    let both = |i: usize| al.a[i] != b'-' && al.b[i] != b'-';
    let first = (0..al.a.len()).find(|&i| both(i))?;
    let last = (0..al.a.len()).rev().find(|&i| both(i))?;
    let (a, b) = (&al.a[first..=last], &al.b[first..=last]);
    let covered = a.iter().filter(|&&c| c != b'-').count();
    let matches = a
        .iter()
        .zip(b)
        .filter(|(x, y)| **x != b'-' && x.eq_ignore_ascii_case(y))
        .count();
    Some((
        covered,
        matches,
        a.len(),
        b.iter().copied().filter(|&c| c != b'-').collect(),
    ))
}

/// What is found at each place where an allele of `locus` lies in the contigs.
pub fn find(contigs: &Contigs, locus: &Locus) -> Vec<Found> {
    let known: HashMap<&[u8], &str> = locus
        .alleles
        .iter()
        .map(|a| (a.seq.as_slice(), a.id.as_str()))
        .collect();
    let mut found = Vec::new();
    for site in sites(&contigs.index, locus) {
        let contig = contigs.index.sequence(site.contig);
        let window =
            &contig[site.start.saturating_sub(MARGIN)..(site.end + MARGIN).min(contig.len())];
        let mut ranked = site.candidates.clone();
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
        ranked.dedup_by_key(|c| c.0);
        let mut best: Option<Found> = None;
        for &(a, _) in ranked.iter().take(CANDIDATES) {
            let allele = &locus.alleles[a];
            // the allele as it lies against the forward strand of the contig
            let oriented = if site.rev {
                revcom(&allele.seq)
            } else {
                allele.seq.clone()
            };
            let Some((covered, matches, columns, bases)) = aligned_over(&oriented, window) else {
                continue;
            };
            let sequence = if site.rev { revcom(&bases) } else { bases };
            let exact = known.get(sequence.as_slice()).map(|id| id.to_string());
            let candidate = Found {
                allele: exact.clone().unwrap_or_else(|| allele.id.clone()),
                exact,
                covered,
                allele_len: allele.seq.len(),
                matches,
                columns,
                sequence,
            };
            let better = best.as_ref().is_none_or(|b| {
                (candidate.exact.is_some(), candidate.matches) > (b.exact.is_some(), b.matches)
            });
            if better {
                best = Some(candidate);
            }
        }
        found.extend(best);
    }
    found
}

#[cfg(test)]
mod tests;
