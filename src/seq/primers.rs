//! Primer tools: candidate search with the original's default limits, the original's
//! self-complementarity score, and PCR product prediction.

use super::dna::{gc_percent, melting_temp, melting_temp_raw, revcom};

#[derive(Clone, Debug)]
pub struct PrimerParams {
    pub len_min: usize,
    pub len_max: usize,
    pub tm_min: f64,
    pub tm_max: f64,
    pub gc_min: f64,
    pub gc_max: f64,
    /// Consecutive G/C bases at the primer's 3' end.
    pub clamp_min: usize,
    pub clamp_max: usize,
    pub self_total_max: i32,
    pub self_3prime_max: i32,
    pub salt_mm: f64,
    pub dna_nm: f64,
}

impl Default for PrimerParams {
    /// The original's defaults (ApE_Defaults.txt).
    fn default() -> Self {
        PrimerParams {
            len_min: 20,
            len_max: 25,
            tm_min: 55.0,
            tm_max: 60.0,
            gc_min: 45.0,
            gc_max: 60.0,
            clamp_min: 1,
            clamp_max: 2,
            self_total_max: 10,
            self_3prime_max: 3,
            salt_mm: 50.0,
            dna_nm: 250.0,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Primer {
    /// 5' to 3' as it would be ordered.
    pub seq: Vec<u8>,
    /// Template coordinates of the bases it matches (0-based, end exclusive).
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
    pub tm: f64,
    pub gc: f64,
    pub clamp: usize,
    pub self_3prime: i32,
    pub self_total: i32,
}

/// Self-complementarity as `(3' score, best total)`, ported from the original's `align`:
/// the primer is aligned against its own reverse complement (match +1, mismatch -1, gap -2).
pub fn self_scores(primer: &[u8]) -> (i32, i32) {
    let p = primer.to_ascii_uppercase();
    let t = revcom(&p);
    let mut prev = vec![0i32; p.len()];
    let (mut best_3, mut best_total) = (-1, -1);
    for &tc in &t {
        let mut col = Vec::with_capacity(p.len());
        let (mut diag, mut above) = (0i32, 0i32);
        for (r, &pc) in p.iter().enumerate() {
            let step = diag + if tc == pc { 1 } else { -1 };
            diag = prev[r];
            above = step.max(above - 2).max(prev[r] - 2);
            col.push(above);
        }
        if above > best_total {
            best_total = above;
        }
        if p.last() == Some(&tc) && above > best_3 {
            best_3 = above;
        }
        prev = col;
    }
    (best_3, best_total)
}

fn gc_clamp(three_prime_first: impl Iterator<Item = u8>) -> usize {
    three_prime_first
        .take_while(|c| matches!(c, b'G' | b'C'))
        .count()
}

/// Candidate primers whose bases lie inside `[from, to)`. Forward primers match the top strand
/// and end at their 3' end; reverse primers are the reverse complement of the window.
pub fn find_primers(
    template: &[u8],
    from: usize,
    to: usize,
    reverse: bool,
    p: &PrimerParams,
) -> Vec<Primer> {
    let to = to.min(template.len());
    let mut out = Vec::new();
    for len in p.len_min..=p.len_max {
        if to < from + len {
            break;
        }
        for start in from..=to - len {
            let window = template[start..start + len].to_ascii_uppercase();
            if window.iter().any(|c| !b"ACGT".contains(c)) {
                continue;
            }
            let Some(tm_raw) = melting_temp_raw(&window, p.dna_nm, p.salt_mm) else {
                continue;
            };
            // the original accepts Tm within [min - 0.5, max + 0.5)
            if tm_raw <= p.tm_min - 0.5 || tm_raw > p.tm_max + 0.4999999999 {
                continue;
            }
            let gc = gc_percent(&window);
            let clamp = if reverse {
                gc_clamp(window.iter().copied())
            } else {
                gc_clamp(window.iter().rev().copied())
            };
            if gc < p.gc_min || gc > p.gc_max || clamp < p.clamp_min || clamp > p.clamp_max {
                continue;
            }
            let seq = if reverse { revcom(&window) } else { window };
            let (self_3prime, self_total) = self_scores(&seq);
            if self_total > p.self_total_max || self_3prime > p.self_3prime_max {
                continue;
            }
            out.push(Primer {
                tm: melting_temp(&seq, p.dna_nm, p.salt_mm).unwrap_or(0.0),
                seq,
                start,
                end: start + len,
                reverse,
                gc,
                clamp,
                self_3prime,
                self_total,
            });
        }
    }
    let mid = (p.tm_min + p.tm_max) / 2.0;
    out.sort_by(|a, b| {
        ((a.tm - mid).abs(), a.start)
            .partial_cmp(&((b.tm - mid).abs(), b.start))
            .unwrap()
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pseudo_dna(len: usize, mut seed: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[((seed >> 33) % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn candidates_respect_every_limit() {
        let t = pseudo_dna(600, 42);
        let p = PrimerParams::default();
        for reverse in [false, true] {
            let found = find_primers(&t, 0, 300, reverse, &p);
            assert!(!found.is_empty(), "no candidates (reverse={reverse})");
            for c in &found {
                assert!((20..=25).contains(&c.seq.len()));
                // the original accepts an unrounded Tm in (min - 0.5, max + 0.5)
                assert!(c.tm >= 54.4 && c.tm <= 60.6, "Tm {}", c.tm);
                assert!(c.gc >= 45.0 && c.gc <= 60.0);
                assert!((1..=2).contains(&c.clamp));
                assert!(c.self_total <= 10 && c.self_3prime <= 3);
                let window = &t[c.start..c.end];
                let expect = if reverse {
                    revcom(window)
                } else {
                    window.to_vec()
                };
                assert_eq!(c.seq, expect);
                let three_prime: Vec<u8> = c.seq.iter().rev().copied().collect();
                assert_eq!(gc_clamp(three_prime.into_iter()), c.clamp);
            }
        }
    }

    #[test]
    fn a_palindromic_primer_scores_high_and_a_random_one_low() {
        // (3' score, total) as printed by the original's own align routine
        let reference: [(&[u8], (i32, i32)); 8] = [
            (b"GATCGATCGATCGATC", (16, 16)),
            (b"ATGGCGTTTGAACGGCTCAA", (7, 7)),
            (b"AAAAAAAAAAAAAAAAAAAA", (-1, -1)),
            (b"TAATACGACTCACTATAGGG", (6, 6)),
            (b"CAGGAAACAGCTATGAC", (2, 2)),
            (b"GCGCGCGCGCGCGCGCGCGC", (20, 20)),
            (b"GTAAAACGACGGCCAGT", (4, 4)),
            (b"ACGTACGTACGTACGTACGT", (20, 20)),
        ];
        for (primer, expected) in reference {
            assert_eq!(
                self_scores(primer),
                expected,
                "{}",
                String::from_utf8_lossy(primer)
            );
        }
    }
}
