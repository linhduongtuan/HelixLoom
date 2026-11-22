//! PAF, the format minimap2 writes (the Pairwise mApping Format): one line to a mapping, twelve columns separated by tabs, then tags.
//!
//! 1 read name, 2 read length, 3 start on the read (from 0), 4 end, 5 strand (`+` or `-`), 6 target name, 7 target length, 8 start, 9 end,
//! 10 matching bases, 11 bases of the mapping with its gaps, 12 mapping quality (255 where there is none). With the base-level alignment, 10 and
//! 11 are counted from it; without it, 10 is the bases that the anchors of the chain cover and 11 the longer of the two spans.
//!
//! Tags: `tp` (`P` primary, `S` secondary), `cm` anchors in the chain, `s1` its score, `s2` the score of its best secondary, `dv` the
//! divergence that the chain suggests, and with the alignment `NM` (mismatches and bases of gaps) and `cg` (the CIGAR).

use super::index::Index;
use super::map::Hit;

pub fn line(name: &str, qlen: usize, index: &Index, h: &Hit, overlaps: bool) -> String {
    let (matches, block) = match &h.base {
        Some(b) => (b.matches, b.block),
        None => (h.covered, (h.qend - h.qstart).max(h.tend - h.tstart)),
    };
    let mapq = if overlaps { 255 } else { h.mapq };
    let mut out = format!(
        "{name}\t{qlen}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{matches}\t{block}\t{mapq}",
        h.qstart,
        h.qend,
        if h.rev { '-' } else { '+' },
        index.name(h.target),
        index.length(h.target),
        h.tstart,
        h.tend,
    );
    out.push_str(&format!(
        "\ttp:A:{}\tcm:i:{}\ts1:i:{}",
        if h.primary { 'P' } else { 'S' },
        h.anchors.len(),
        h.score.round() as i64
    ));
    if h.primary && h.best_secondary > 0.0 && !overlaps {
        out.push_str(&format!("\ts2:i:{}", h.best_secondary.round() as i64));
    }
    if let Some(d) = h.divergence {
        out.push_str(&format!("\tdv:f:{d:.4}"));
    }
    if let Some(b) = &h.base {
        out.push_str(&format!("\tNM:i:{}\tcg:Z:{}", b.differences, b.cigar));
    }
    out
}

#[cfg(test)]
mod tests;
