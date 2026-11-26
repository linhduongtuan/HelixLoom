//! Two checks on the contigs of a long-read assembly of a plasmid, in the manner of the NanoSeq route (its README, which gives a 30 base
//! overlap and 95 % identity for the ends, was read by someone else on 19 November 2026 and could not be found again to be read for this;
//! `docs/ROUTES.md` says so): a **concatemer** is a contig that is several head-to-tail copies of one sequence, which assemblers of
//! reads that go round a small circle more than once make; a **circle** has its end overlapping its start, and is made one by
//! cutting the overlap off.
//!
//! The way each is found is this program's own: the ends are compared base to base (no gaps) for every length from the shortest to the
//! longest overlap allowed; a concatemer is looked for by mapping each of the `n` equal parts of the contig to the first, with the
//! program's mapper, for the largest `n` that fits.

use crate::seq::minimap::{map_read, Index, Options as MapOptions};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Overlap {
    /// The shortest overlap, in bases, and how alike the end and the start must be over it (a share, 0 to 1).
    pub min_length: usize,
    pub min_identity: f64,
    /// The longest overlap that is looked for (and never more than half the contig).
    pub max_length: usize,
}

impl Default for Overlap {
    fn default() -> Self {
        Overlap {
            min_length: 30,
            min_identity: 0.95,
            max_length: 5000,
        }
    }
}

/// The length of the overlap between the end of `seq` and its start: the length (from `min_length` up) over which the last bases are most
/// alike the first, the longer when they tie, provided that they are alike enough; none when there is no such length.
pub fn end_overlap(seq: &[u8], o: &Overlap) -> Option<usize> {
    let n = seq.len();
    let most = o.max_length.min(n / 2);
    let mut best: Option<(f64, usize)> = None;
    for l in o.min_length.max(1)..=most {
        let same = (0..l)
            .filter(|&i| seq[n - l + i].eq_ignore_ascii_case(&seq[i]))
            .count();
        let identity = same as f64 / l as f64;
        // a hair of slack, so that 95 % of 20 bases (19 of 20) is not lost to the rounding of 0.95
        if identity + 1e-9 >= o.min_identity
            && best.is_none_or(|(b, bl)| identity > b || (identity == b && l > bl))
        {
            best = Some((identity, l));
        }
    }
    best.map(|b| b.1)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Concatemer {
    /// The shortest unit that is looked for.
    pub min_unit: usize,
    /// The most copies that are looked for.
    pub max_copies: usize,
    /// How alike each copy must be to the first, and how much of it the alignment must cover (shares, 0 to 1).
    pub min_identity: f64,
    pub min_covered: f64,
}

impl Default for Concatemer {
    fn default() -> Self {
        Concatemer {
            min_unit: 200,
            max_copies: 10,
            min_identity: 0.95,
            min_covered: 0.9,
        }
    }
}

/// A contig that is `copies` head-to-tail copies of a unit of `unit` bases.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Found {
    pub unit: usize,
    pub copies: usize,
}

/// What an alignment by minimizers may fall short of the ends of what it aligns, at both ends together (the first and the last minimizer are
/// not at the very ends), so that a unit that is a whole copy still counts as covered. Little, so that a part that is only partly the unit
/// (a sequence that repeats at another period) does not.
const END_SLOP: usize = 30;

/// Whether `part` maps to the unit that `index` holds, over most of its length and alike enough.
fn alike(index: &Index, part: &[u8], o: &Concatemer, map: &MapOptions) -> bool {
    map_read(index, part, map, &|_| false)
        .into_iter()
        .filter(|h| h.primary)
        .any(|h| {
            (h.qend - h.qstart + END_SLOP) as f64 >= o.min_covered * part.len() as f64
                && 1.0 - h.divergence.unwrap_or(0.0) + 1e-9 >= o.min_identity
        })
}

/// Whether `seq` is two or more copies of one unit, and the most copies for which it is: the contig is cut into that many parts of one
/// length (the last takes what is left), and each of the others must be alike the first.
pub fn find_concatemer(seq: &[u8], o: &Concatemer) -> Option<Found> {
    let map = MapOptions::default();
    for copies in (2..=o.max_copies).rev() {
        let unit = seq.len() / copies;
        if unit < o.min_unit {
            continue;
        }
        let Ok(index) = Index::build(vec![("unit".to_string(), seq[..unit].to_vec())], &map) else {
            continue;
        };
        let every_part_alike = (1..copies).all(|i| {
            let end = if i + 1 == copies {
                seq.len()
            } else {
                (i + 1) * unit
            };
            alike(&index, &seq[i * unit..end], o, &map)
        });
        if every_part_alike {
            return Some(Found { unit, copies });
        }
    }
    None
}

#[cfg(test)]
mod tests;
