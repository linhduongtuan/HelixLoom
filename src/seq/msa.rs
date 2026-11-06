//! Multiple sequence alignment by the centre-star method: the sequence that shares the most
//! 6-mers with the others becomes the centre, every other sequence is aligned to it with the
//! pairwise aligner (free ends, the original's alignment scores), and the pairwise alignments are
//! merged so that a gap opened for one sequence stays open for all ("once a gap, always a gap").
//! A small stand-in for MAFFT or MUSCLE; it suits sequences that are clearly related.

use super::align::{align, Scoring};
use super::par::par_map;
use std::collections::HashSet;

const KMER: usize = 6;
const GAP: u8 = b'-';

#[derive(Clone, Debug, PartialEq)]
pub struct Msa {
    pub names: Vec<String>,
    /// Aligned rows of equal length, gaps written `-`.
    pub rows: Vec<Vec<u8>>,
}

impl Msa {
    pub fn columns(&self) -> usize {
        self.rows.first().map_or(0, Vec::len)
    }
}

/// How one sequence sits against the centre: what is inserted before each centre base, and what
/// stands opposite each centre base.
struct Placed {
    before: Vec<Vec<u8>>,
    opposite: Vec<u8>,
}

fn place(centre_aligned: &[u8], other_aligned: &[u8], centre_len: usize) -> Placed {
    let mut before = vec![Vec::new(); centre_len + 1];
    let mut opposite = Vec::with_capacity(centre_len);
    for (&c, &o) in centre_aligned.iter().zip(other_aligned) {
        if c == GAP {
            before[opposite.len()].push(o);
        } else {
            opposite.push(o);
        }
    }
    Placed { before, opposite }
}

/// The sequence sharing the most 6-mers with all the others (the first one on a tie).
fn centre(seqs: &[Vec<u8>]) -> usize {
    let upper: Vec<Vec<u8>> = seqs.iter().map(|s| s.to_ascii_uppercase()).collect();
    let sets: Vec<HashSet<&[u8]>> = upper.iter().map(|s| s.windows(KMER).collect()).collect();
    let shared = |i: usize| -> usize {
        (0..sets.len())
            .filter(|&j| j != i)
            .map(|j| sets[i].intersection(&sets[j]).count())
            .sum()
    };
    let scores: Vec<usize> = (0..seqs.len()).map(shared).collect();
    let best = scores.iter().max().copied().unwrap_or(0);
    scores.iter().position(|&s| s == best).unwrap_or(0)
}

/// Align `seqs` (name, bases). `Err` when there is nothing to align or a pair is too long.
pub fn center_star(seqs: &[(String, Vec<u8>)]) -> Result<Msa, String> {
    let (names, bases): (Vec<String>, Vec<Vec<u8>>) = seqs.iter().cloned().unzip();
    if bases.is_empty() || bases.iter().any(Vec::is_empty) {
        return Err("Every sequence needs at least one base.".into());
    }
    let c = centre(&bases);
    let centre_len = bases[c].len();
    let others: Vec<usize> = (0..bases.len()).filter(|&i| i != c).collect();
    let aligned = par_map(&others, |&i| {
        align(&bases[c], &bases[i], &Scoring::default(), true)
            .map(|al| place(&al.a, &al.b, centre_len))
            .map_err(|e| format!("{} and {}: {e}", names[c], names[i]))
    });
    let mut placed: Vec<Option<Placed>> = (0..bases.len()).map(|_| None).collect();
    for (&i, p) in others.iter().zip(aligned) {
        placed[i] = Some(p?);
    }
    // gap columns needed before each centre base: the widest insertion any sequence has there
    let width: Vec<usize> = (0..=centre_len)
        .map(|p| {
            placed
                .iter()
                .flatten()
                .map(|pl| pl.before[p].len())
                .max()
                .unwrap_or(0)
        })
        .collect();
    let rows = (0..bases.len())
        .map(|i| {
            let mut row = Vec::new();
            for p in 0..=centre_len {
                let inserted: &[u8] = placed[i].as_ref().map_or(&[], |pl| &pl.before[p]);
                row.extend_from_slice(inserted);
                row.extend(std::iter::repeat_n(GAP, width[p] - inserted.len()));
                if p < centre_len {
                    row.push(placed[i].as_ref().map_or(bases[c][p], |pl| pl.opposite[p]));
                }
            }
            row
        })
        .collect();
    Ok(Msa { names, rows })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seqs(list: &[(&str, &str)]) -> Vec<(String, Vec<u8>)> {
        list.iter()
            .map(|(n, s)| (n.to_string(), s.as_bytes().to_vec()))
            .collect()
    }

    fn degapped(row: &[u8]) -> Vec<u8> {
        row.iter().copied().filter(|&c| c != GAP).collect()
    }

    #[test]
    fn identical_sequences_align_without_gaps() {
        let s = "ACGTACGGTCAGTCAGGTACCAGT";
        let msa = center_star(&seqs(&[("a", s), ("b", s), ("c", s)])).unwrap();
        assert!(msa.rows.iter().all(|r| r == s.as_bytes()));
        assert_eq!(msa.columns(), s.len());
    }

    #[test]
    fn a_deletion_in_one_sequence_becomes_a_gap_and_every_row_keeps_its_bases() {
        let full = "ACGTTGCAAGGCTTAACCGGTTACGTACGGCATCAG";
        let cut = "ACGTTGCAAGGCTTACGGTTACGTACGGCATCAG"; // "AC" missing after TTA
        let other = "ACGTTGCAAGGCTTAACCGGTTACGTACGGCTTCAG"; // one substitution
        let input = seqs(&[("full", full), ("cut", cut), ("other", other)]);
        let msa = center_star(&input).unwrap();
        let width = msa.columns();
        assert!(msa.rows.iter().all(|r| r.len() == width));
        for (row, (_, original)) in msa.rows.iter().zip(&input) {
            assert_eq!(&degapped(row), original);
        }
        let cut_row = &msa.rows[1];
        assert_eq!(cut_row.iter().filter(|&&c| c == GAP).count(), 2);
    }

    #[test]
    fn an_insertion_in_one_sequence_opens_a_gap_in_all_the_others() {
        let base = "TTGACCAGGTCAGTTACGGATCCAGTAGCATGCAAC";
        let longer = "TTGACCAGGTCAGTTAAACGGATCCAGTAGCATGCAAC"; // "AA" added
        let msa = center_star(&seqs(&[("a", base), ("b", longer), ("c", base)])).unwrap();
        let width = msa.columns();
        assert_eq!(width, longer.len());
        assert_eq!(msa.rows[0].iter().filter(|&&c| c == GAP).count(), 2);
        assert_eq!(msa.rows[1].iter().filter(|&&c| c == GAP).count(), 0);
    }

    #[test]
    fn empty_input_is_an_error_and_one_sequence_is_itself() {
        assert!(center_star(&[]).is_err());
        assert!(center_star(&seqs(&[("a", "")])).is_err());
        let one = center_star(&seqs(&[("a", "ACGT")])).unwrap();
        assert_eq!(one.rows, vec![b"ACGT".to_vec()]);
    }
}
