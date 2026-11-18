//! Multiple sequence alignment by the centre-star method: the sequence that shares the most
//! 6-mers with the others becomes the centre, every other sequence is aligned to it with the
//! pairwise aligner (free ends, the original's alignment scores), and the pairwise alignments are
//! merged so that a gap opened for one sequence stays open for all ("once a gap, always a gap").
//! A small stand-in for MAFFT or MUSCLE; it suits sequences that are clearly related.

use super::align::{align_many, Scoring};
use std::collections::{HashMap, HashSet};

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
    first_best(&shared_scores(&upper))
}

/// For each sequence (upper case), how many different 6-mers it shares with each of the others, added up. A 6-mer that `df` sequences
/// hold is shared by a sequence that has it with the `df - 1` others that have it, so the sum is, over the different 6-mers of the
/// sequence, `df - 1`: one pass over all the 6-mers, where comparing every pair of sets takes time as the square of the number of
/// sequences (it took half a minute for 600 sequences of 1,500 bases).
pub(crate) fn shared_scores(upper: &[Vec<u8>]) -> Vec<u64> {
    let sets: Vec<HashSet<&[u8]>> = upper.iter().map(|s| s.windows(KMER).collect()).collect();
    let mut holders: HashMap<&[u8], u64> = HashMap::new();
    for set in &sets {
        for &word in set {
            *holders.entry(word).or_insert(0) += 1;
        }
    }
    sets.iter()
        .map(|set| set.iter().map(|word| holders[word] - 1).sum())
        .collect()
}

/// The position of the largest score (the first on a tie).
fn first_best(scores: &[u64]) -> usize {
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
    // every other sequence against the centre, eight at a time and on all cores (`align_many`)
    let others_bases: Vec<&[u8]> = others.iter().map(|&i| bases[i].as_slice()).collect();
    let aligned = align_many(&bases[c], &others_bases, &Scoring::default(), true);
    let mut placed: Vec<Option<Placed>> = (0..bases.len()).map(|_| None).collect();
    for (&i, al) in others.iter().zip(aligned) {
        let al = al.map_err(|e| format!("{} and {}: {e}", names[c], names[i]))?;
        placed[i] = Some(place(&al.a, &al.b, centre_len));
    }
    Ok(merge(names, &bases, c, &placed))
}

/// The pairwise alignments against the centre `c` (`None` for the centre itself) merged into one alignment of all the
/// sequences.
fn merge(names: Vec<String>, bases: &[Vec<u8>], c: usize, placed: &[Option<Placed>]) -> Msa {
    let centre_len = bases[c].len();
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
    Msa { names, rows }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::align::align;

    /// The score as the first version worked it out, by intersecting the sets of every pair: the oracle for `shared_scores`.
    fn pairwise_scores(upper: &[Vec<u8>]) -> Vec<u64> {
        let sets: Vec<HashSet<&[u8]>> = upper.iter().map(|s| s.windows(KMER).collect()).collect();
        (0..sets.len())
            .map(|i| {
                (0..sets.len())
                    .filter(|&j| j != i)
                    .map(|j| sets[i].intersection(&sets[j]).count() as u64)
                    .sum()
            })
            .collect()
    }

    #[test]
    fn the_one_pass_score_is_the_sum_over_every_pair_of_the_shared_six_mers() {
        let mut rng = crate::seq::tree::Rng::new(41);
        for (count, length) in [(1usize, 30usize), (2, 4), (7, 60), (40, 200), (25, 500)] {
            let seqs: Vec<Vec<u8>> = (0..count)
                .map(|_| (0..length).map(|_| b"ACGT"[rng.below(4)]).collect())
                .collect();
            assert_eq!(
                shared_scores(&seqs),
                pairwise_scores(&seqs),
                "{count} x {length}"
            );
        }
        // related sequences (copies with changes), one shorter than a 6-mer, one repeated, one empty
        let base: Vec<u8> = (0..300).map(|_| b"ACGT"[rng.below(4)]).collect();
        let mut seqs: Vec<Vec<u8>> = (0..12)
            .map(|_| {
                base.iter()
                    .map(|&c| {
                        if rng.below(10) == 0 {
                            b"ACGT"[rng.below(4)]
                        } else {
                            c
                        }
                    })
                    .collect()
            })
            .collect();
        seqs.push(b"ACG".to_vec());
        seqs.push(seqs[3].clone());
        seqs.push(Vec::new());
        assert_eq!(shared_scores(&seqs), pairwise_scores(&seqs));
        // a protein-sized alphabet, the same
        let letters = b"ACDEFGHIKLMNPQRSTVWY";
        let proteins: Vec<Vec<u8>> = (0..15)
            .map(|_| (0..120).map(|_| letters[rng.below(20)]).collect())
            .collect();
        assert_eq!(shared_scores(&proteins), pairwise_scores(&proteins));
    }

    #[test]
    fn sequences_that_share_nothing_score_nothing_and_the_centre_is_the_first_then() {
        let apart = vec![
            b"AAAAAAAA".to_vec(),
            b"CCCCCCCC".to_vec(),
            b"GGGGGGGG".to_vec(),
        ];
        assert_eq!(shared_scores(&apart), vec![0, 0, 0]);
        assert_eq!(centre(&apart), 0);
    }

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

    /// `count` related sequences of about `length` bases: a random one, then copies with one base in 30 changed, and short
    /// insertions and deletions.
    fn related(count: usize, length: usize, seed: u64) -> Vec<(String, Vec<u8>)> {
        let mut rng = crate::seq::tree::Rng::new(seed);
        let base: Vec<u8> = (0..length).map(|_| b"ACGT"[rng.below(4)]).collect();
        (0..count)
            .map(|k| {
                let mut s = Vec::with_capacity(length + 16);
                for &c in &base {
                    match rng.below(120) {
                        0..=3 => s.push(b"ACGT"[rng.below(4)]),
                        4 => {}
                        5 => {
                            s.push(c);
                            s.extend((0..1 + rng.below(4)).map(|_| b"ACGT"[rng.below(4)]));
                        }
                        _ => s.push(c),
                    }
                }
                (format!("s{k}"), s)
            })
            .collect()
    }

    #[test]
    #[ignore = "timing: cargo test [--release] timing_of_center_star -- --ignored --nocapture"]
    fn timing_of_center_star() {
        let input = related(200, 1500, 11);
        let bases: Vec<Vec<u8>> = input.iter().map(|(_, s)| s.clone()).collect();
        let best = |f: &dyn Fn()| {
            (0..3)
                .map(|_| {
                    let start = std::time::Instant::now();
                    f();
                    start.elapsed().as_secs_f64()
                })
                .fold(f64::MAX, f64::min)
        };
        let whole = best(&|| {
            std::hint::black_box(center_star(&input).unwrap());
        });
        let choose = best(&|| {
            std::hint::black_box(centre(&bases));
        });
        let c = centre(&bases);
        let others: Vec<&[u8]> = (0..bases.len())
            .filter(|&i| i != c)
            .map(|i| bases[i].as_slice())
            .collect();
        let pairs = best(&|| {
            std::hint::black_box(align_many(&bases[c], &others, &Scoring::default(), true));
        });
        println!(
            "center_star of 200 x 1500: {whole:.2} s (choosing the centre {choose:.2} s, the 199 alignments {pairs:.2} s)"
        );
    }

    #[test]
    fn the_alignments_eight_at_a_time_give_the_same_msa_as_one_pair_at_a_time() {
        // sets of one group and a part, of two groups and a part, with insertions and deletions so the lengths differ
        for (count, length, seed) in [
            (2usize, 60usize, 1u64),
            (9, 150, 2),
            (20, 200, 3),
            (33, 90, 4),
        ] {
            let input = related(count, length, seed);
            let msa = center_star(&input).unwrap();
            let bases: Vec<Vec<u8>> = input.iter().map(|(_, s)| s.clone()).collect();
            let c = centre(&bases);
            let placed: Vec<Option<Placed>> = (0..count)
                .map(|i| {
                    (i != c).then(|| {
                        let al = align(&bases[c], &bases[i], &Scoring::default(), true).unwrap();
                        place(&al.a, &al.b, bases[c].len())
                    })
                })
                .collect();
            let names = input.iter().map(|(n, _)| n.clone()).collect();
            assert_eq!(msa, merge(names, &bases, c, &placed), "{count} x {length}");
            let width = msa.columns();
            for (row, original) in msa.rows.iter().zip(&bases) {
                assert_eq!(row.len(), width);
                assert_eq!(&degapped(row), original);
            }
        }
    }

    #[test]
    fn empty_input_is_an_error_and_one_sequence_is_itself() {
        assert!(center_star(&[]).is_err());
        assert!(center_star(&seqs(&[("a", "")])).is_err());
        let one = center_star(&seqs(&[("a", "ACGT")])).unwrap();
        assert_eq!(one.rows, vec![b"ACGT".to_vec()]);
    }
}
