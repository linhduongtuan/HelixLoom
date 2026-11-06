//! Primer pair design after primer3's generic task (github.com/primer3-org/primer3).
//!
//! Defaults come from `pr_set_default_global_args_1/_2` in `libprimer3.cc`, the melting
//! temperature from `oligotm()` in `oligotm.c` (SantaLucia 1998 table, end and symmetry
//! corrections, divalent salt turned into monovalent), the penalties from `p_obj_fn` and `obj_fn`.
//! Written from reading those rules and checked against a build of primer3 itself; no code was
//! copied. Not ported: the thermodynamic dimer and hairpin checks (the older +1/-1/-2
//! complementarity scores are used instead, as with primer3's own non-thermodynamic mode),
//! template mispriming, internal oligos, excluded regions and several product-size ranges.

use super::dna::{gc_percent, nearest_neighbour, revcom};
use std::cmp::{Ordering, Reverse};
use std::collections::BinaryHeap;

/// primer3's defaults, named after its tags.
#[derive(Clone, Debug, PartialEq)]
pub struct Settings {
    /// PRIMER_OPT_SIZE, PRIMER_MIN_SIZE, PRIMER_MAX_SIZE
    pub opt_size: usize,
    pub min_size: usize,
    pub max_size: usize,
    /// PRIMER_OPT_TM, PRIMER_MIN_TM, PRIMER_MAX_TM
    pub opt_tm: f64,
    pub min_tm: f64,
    pub max_tm: f64,
    /// PRIMER_MIN_GC, PRIMER_MAX_GC (percent)
    pub min_gc: f64,
    pub max_gc: f64,
    /// PRIMER_SALT_MONOVALENT, PRIMER_SALT_DIVALENT, PRIMER_DNTP_CONC (mM)
    pub salt_mm: f64,
    pub divalent_mm: f64,
    pub dntp_mm: f64,
    /// PRIMER_DNA_CONC (nM)
    pub dna_nm: f64,
    /// PRIMER_MAX_POLY_X, PRIMER_GC_CLAMP, PRIMER_MAX_END_GC
    pub max_poly_x: usize,
    pub gc_clamp: usize,
    pub max_end_gc: usize,
    /// PRIMER_MAX_SELF_ANY, PRIMER_MAX_SELF_END
    pub max_self_any: i32,
    pub max_self_end: i32,
    /// PRIMER_PAIR_MAX_COMPL_ANY, PRIMER_PAIR_MAX_COMPL_END
    pub pair_compl_any: i32,
    pub pair_compl_end: i32,
    /// PRIMER_PAIR_MAX_DIFF_TM
    pub max_diff_tm: f64,
    /// PRIMER_PRODUCT_SIZE_RANGE
    pub product_min: usize,
    pub product_max: usize,
    /// PRIMER_NUM_RETURN
    pub num_return: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            opt_size: 20,
            min_size: 18,
            max_size: 27,
            opt_tm: 60.0,
            min_tm: 57.0,
            max_tm: 63.0,
            min_gc: 20.0,
            max_gc: 80.0,
            salt_mm: 50.0,
            divalent_mm: 1.5,
            dntp_mm: 0.6,
            dna_nm: 50.0,
            max_poly_x: 5,
            gc_clamp: 0,
            max_end_gc: 5,
            max_self_any: 8,
            max_self_end: 3,
            pair_compl_any: 8,
            pair_compl_end: 3,
            max_diff_tm: 100.0,
            product_min: 100,
            product_max: 300,
            num_return: 5,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Oligo {
    /// 5' to 3' as it would be ordered (the reverse complement of the template for a right primer).
    pub seq: Vec<u8>,
    /// Template bases it covers, 0-based, end exclusive.
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
    pub tm: f64,
    pub gc: f64,
    pub penalty: f64,
    pub self_any: i32,
    pub self_end: i32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Pair {
    pub left: Oligo,
    pub right: Oligo,
    pub penalty: f64,
    pub product_size: usize,
    pub compl_any: i32,
    pub compl_end: i32,
}

/// primer3 counts 120 monovalent mM for every sqrt(mM) of divalent salt left over after dNTPs.
fn divalent_as_monovalent(divalent: f64, dntp: f64) -> Option<f64> {
    let dntp = if divalent == 0.0 { 0.0 } else { dntp };
    if divalent < 0.0 || dntp < 0.0 {
        return None;
    }
    Some(120.0 * (divalent.max(dntp) - dntp).sqrt())
}

/// Melting temperature the way primer3 computes it; `None` unless `seq` is at least two bases
/// of A, C, G and T.
pub fn tm(seq: &[u8], s: &Settings) -> Option<f64> {
    let seq = seq.to_ascii_uppercase();
    if seq.len() < 2 || !seq.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
        return None;
    }
    let symmetric = revcom(&seq) == seq;
    let (mut dh, mut ds) = (0i32, if symmetric { 14 } else { 0 });
    for end in [seq[0], seq[seq.len() - 1]] {
        if matches!(end, b'A' | b'T') {
            (dh, ds) = (dh - 23, ds - 41);
        } else {
            (dh, ds) = (dh - 1, ds + 28);
        }
    }
    for w in seq.windows(2) {
        let (h, e) = nearest_neighbour(w[0].to_ascii_lowercase(), w[1].to_ascii_lowercase());
        (dh, ds) = (dh + h, ds + e);
    }
    let mono = s.salt_mm + divalent_as_monovalent(s.divalent_mm, s.dntp_mm)?;
    let delta_h = f64::from(dh) * -100.0;
    let delta_s = f64::from(ds) * -0.1 + 0.368 * (seq.len() - 1) as f64 * (mono / 1000.0).ln();
    let strands = if symmetric { 1e9 } else { 4e9 };
    Some(delta_h / (delta_s + 1.987 * (s.dna_nm / strands).ln()) - 273.15)
}

/// primer3's own dimer score (its `dpal` with the default nucleotide scores): `a` against the
/// reverse complement of `b`, match +1, mismatch -1, a gap of one base -2 and no longer gaps.
/// Returns `(any, end)`: the best local alignment, and the best alignment that starts at the
/// start of either sequence and finishes on `a`'s last base, its 3' end. Neither is below 0.
fn dimer_scores(a: &[u8], b: &[u8]) -> (i32, i32) {
    let (x, y) = (a, revcom(b));
    let (m, n) = (x.len(), y.len());
    if m == 0 || n == 0 {
        return (0, 0);
    }
    let run = |local: bool| -> i32 {
        let mut s = vec![vec![0i32; n]; m];
        let mut best = i32::MIN;
        for i in 0..m {
            for j in 0..n {
                let mut before = 0;
                if i > 0 && j > 0 {
                    before = s[i - 1][j - 1];
                    if i > 1 {
                        before = before.max(s[i - 2][j - 1] - 2);
                    }
                    if j > 1 {
                        before = before.max(s[i - 1][j - 2] - 2);
                    }
                }
                let hit = if x[i] == y[j] { 1 } else { -1 };
                s[i][j] = if local {
                    (before + hit).max(0)
                } else {
                    before + hit
                };
                if local || i == m - 1 {
                    best = best.max(s[i][j]);
                }
            }
        }
        best.max(0)
    };
    (run(true), run(false))
}

fn longest_run(seq: &[u8]) -> usize {
    let (mut best, mut run, mut prev) = (0, 0, 0u8);
    for &c in seq {
        run = if c == prev { run + 1 } else { 1 };
        (prev, best) = (c, best.max(run));
    }
    best
}

/// One candidate oligo, or `None` when it breaks a limit. Its penalty is the distance from the
/// optimum Tm plus the distance from the optimum size (primer3's default weights are all 1).
fn evaluate(seq: Vec<u8>, span: (usize, usize), reverse: bool, s: &Settings) -> Option<Oligo> {
    if !seq.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
        return None;
    }
    let gc = gc_percent(&seq);
    if gc < s.min_gc || gc > s.max_gc || longest_run(&seq) > s.max_poly_x {
        return None;
    }
    let tm = tm(&seq, s)?;
    if tm < s.min_tm || tm > s.max_tm {
        return None;
    }
    let is_gc = |c: &u8| matches!(c, b'G' | b'C');
    let end_gc = seq[seq.len().saturating_sub(5)..]
        .iter()
        .filter(|c| is_gc(c))
        .count();
    let clamp = &seq[seq.len().saturating_sub(s.gc_clamp)..];
    if end_gc > s.max_end_gc || !clamp.iter().all(is_gc) {
        return None;
    }
    let (self_any, self_end) = dimer_scores(&seq, &seq);
    if self_any > s.max_self_any || self_end > s.max_self_end {
        return None;
    }
    let penalty = (tm - s.opt_tm).abs() + (seq.len() as f64 - s.opt_size as f64).abs();
    Some(Oligo {
        seq,
        start: span.0,
        end: span.1,
        reverse,
        tm,
        gc,
        penalty,
        self_any,
        self_end,
    })
}

/// Candidates for one side. Left primers are listed by start, right primers by end, shortest first.
fn candidates(
    template: &[u8],
    s: &Settings,
    reverse: bool,
    keep: impl Fn(usize, usize) -> bool,
) -> Vec<Oligo> {
    let n = template.len();
    let mut out = Vec::new();
    for edge in 0..=n {
        for len in s.min_size..=s.max_size {
            let span = if reverse {
                edge.checked_sub(len).map(|a| (a, edge))
            } else {
                (edge + len <= n).then_some((edge, edge + len))
            };
            let Some((a, b)) = span.filter(|&(a, b)| keep(a, b)) else {
                continue;
            };
            let seq = if reverse {
                revcom(&template[a..b])
            } else {
                template[a..b].to_vec()
            };
            out.extend(evaluate(seq, (a, b), reverse, s));
        }
    }
    out
}

#[derive(PartialEq)]
struct Key(f64);
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Key {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

/// The best pairs, lowest penalty first. A pair's penalty is the sum of its two oligo penalties;
/// pairs are taken in that order until `num_return` also pass the dimer limits. `target` is
/// `(start, length)` of bases the product must contain; primers may not overlap it.
pub fn design(template: &[u8], target: Option<(usize, usize)>, s: &Settings) -> Vec<Pair> {
    let template = template.to_ascii_uppercase();
    let (t_start, t_end) = target.map_or((template.len(), 0), |(a, l)| (a, a + l));
    let left_ok =
        |a: usize, b: usize| target.is_none() || (b <= t_start && a + s.product_max >= t_end);
    let right_ok =
        |a: usize, b: usize| target.is_none() || (a >= t_end && b <= t_start + s.product_max);
    let lefts = candidates(&template, s, false, left_ok);
    let rights = candidates(&template, s, true, right_ok);
    // rights a left primer can pair with, by product size, as indexes sorted by penalty
    let window = |l: &Oligo| -> Vec<usize> {
        let lo = rights.partition_point(|r| r.end < l.start + s.product_min);
        let hi = rights.partition_point(|r| r.end <= l.start + s.product_max);
        let fits =
            |j: &usize| rights[*j].start >= l.end && (l.tm - rights[*j].tm).abs() <= s.max_diff_tm;
        let mut ids: Vec<usize> = (lo..hi).filter(fits).collect();
        ids.sort_by(|&a, &b| {
            rights[a]
                .penalty
                .total_cmp(&rights[b].penalty)
                .then(a.cmp(&b))
        });
        ids
    };
    let mut heap = BinaryHeap::new();
    for (i, l) in lefts.iter().enumerate() {
        if let Some(&j) = window(l).first() {
            heap.push(Reverse((Key(l.penalty + rights[j].penalty), i, 0)));
        }
    }
    let mut windows: Vec<Option<Vec<usize>>> = vec![None; lefts.len()];
    let mut found = Vec::new();
    while let Some(Reverse((Key(penalty), i, k))) = heap.pop() {
        let list = windows[i].get_or_insert_with(|| window(&lefts[i]));
        let j = list[k];
        if let Some(&next) = list.get(k + 1) {
            heap.push(Reverse((
                Key(lefts[i].penalty + rights[next].penalty),
                i,
                k + 1,
            )));
        }
        found.extend(compatible(&lefts[i], &rights[j], penalty, s));
        if found.len() >= s.num_return {
            break;
        }
    }
    found
}

/// The pair, if the two primers don't stick to each other too much. The 3' end score counts both
/// ways round: either primer's 3' end on the other.
fn compatible(left: &Oligo, right: &Oligo, penalty: f64, s: &Settings) -> Option<Pair> {
    let (compl_any, end_left) = dimer_scores(&left.seq, &right.seq);
    let (_, end_right) = dimer_scores(&right.seq, &left.seq);
    let compl_end = end_left.max(end_right);
    (compl_any <= s.pair_compl_any && compl_end <= s.pair_compl_end).then(|| Pair {
        left: left.clone(),
        right: right.clone(),
        penalty,
        product_size: right.end - left.start,
        compl_any,
        compl_end,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 1,000 bases from a small deterministic generator (the same string was given to primer3).
    fn template() -> Vec<u8> {
        let mut x = 777u32;
        (0..1000)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                b"ACGT"[((x >> 16) & 3) as usize]
            })
            .collect()
    }

    // Reference values come from `oligotm()` of primer3 itself, compiled from its source, with
    // 50 nM DNA, 50 mM monovalent, 1.5 mM divalent and 0.6 mM dNTP salt.
    #[test]
    fn melting_temperature_matches_primer3() {
        let reference = [
            ("ATGGCGTTTGAACGGCTCAA", 61.1696613315),
            ("GATCAAGGAA", 21.5295375494),
            ("GCGCGCGCGCGCGCGCGCGC", 83.4052962664),
            ("GAATTCGAATTC", 30.9483296473),
            ("TAATACGACTCACTATAGGG", 50.3240382171),
            ("CAGGAAACAGCTATGAC", 49.0953229079),
            ("GTAAAACGACGGCCAGT", 54.6954946911),
            ("AAAAAAAAAAAAAAAAAAAA", 41.5085435030),
            ("GGGGGGGGGGGGGGGGGGGG", 77.8277859500),
            ("ACGTACGTACGTACGTACGT", 59.5974899721),
        ];
        for (seq, want) in reference {
            let got = tm(seq.as_bytes(), &Settings::default()).unwrap();
            assert!((got - want).abs() < 1e-8, "{seq}: {got} vs {want}");
        }
        assert_eq!(tm(b"ACGN", &Settings::default()), None);
    }

    /// (left start, left length, right end, right length, penalty, product size) as primer3 printed them.
    fn summary(pairs: &[Pair]) -> Vec<(usize, usize, usize, usize, f64, usize)> {
        pairs
            .iter()
            .map(|p| {
                let (l, r) = (&p.left, &p.right);
                (
                    l.start,
                    l.seq.len(),
                    r.end,
                    r.seq.len(),
                    p.penalty,
                    p.product_size,
                )
            })
            .collect()
    }

    fn assert_same(got: &[Pair], want: &[(usize, usize, usize, usize, f64, usize)]) {
        let got = summary(got);
        assert_eq!(got.len(), want.len(), "{got:?}");
        for (g, w) in got.iter().zip(want) {
            assert_eq!(
                (g.0, g.1, g.2, g.3, g.5),
                (w.0, w.1, w.2, w.3, w.5),
                "{got:?}"
            );
            assert!((g.4 - w.4).abs() < 1e-5, "penalty {} vs {}", g.4, w.4);
        }
    }

    #[test]
    fn pairs_around_a_target_match_primer3() {
        let pairs = design(&template(), Some((480, 40)), &Settings::default());
        assert_same(
            &pairs,
            &[
                (456, 20, 656, 20, 0.146610, 200),
                (373, 20, 656, 20, 0.219692, 283),
                (354, 20, 639, 20, 0.278366, 285),
                (456, 20, 655, 20, 0.281953, 199),
                (456, 20, 675, 20, 0.285056, 219),
            ],
        );
        let first = &pairs[0];
        assert_eq!(first.left.seq, b"TATGTCCAGTGCCGCATACC");
        assert_eq!(first.right.seq, b"TCTCCTCCTTAACCGCTGGA");
        assert!((first.left.tm - 59.894).abs() < 5e-4 && (first.right.tm - 59.959).abs() < 5e-4);
        assert_eq!((first.left.self_any, first.left.self_end), (5, 0));
        assert_eq!((first.right.self_any, first.right.self_end), (4, 0));
        assert_eq!(dimers(&pairs), [(5, 1), (5, 3), (5, 2), (6, 0), (4, 3)]);
    }

    /// (COMPL_ANY, COMPL_END) of each pair as primer3 printed them.
    fn dimers(pairs: &[Pair]) -> Vec<(i32, i32)> {
        pairs.iter().map(|p| (p.compl_any, p.compl_end)).collect()
    }

    #[test]
    fn the_best_pairs_anywhere_in_a_template_match_primer3() {
        let pairs = design(&template(), None, &Settings::default());
        assert_same(
            &pairs,
            &[
                (148, 20, 334, 20, 0.056221, 186),
                (42, 20, 168, 20, 0.065822, 126),
                (42, 20, 334, 20, 0.066348, 292),
                (687, 20, 905, 20, 0.070760, 218),
                (751, 20, 905, 20, 0.131065, 154),
            ],
        );
        assert_eq!(dimers(&pairs), [(5, 3), (3, 1), (4, 0), (3, 0), (4, 2)]);
    }

    #[test]
    fn limits_are_obeyed() {
        let s = Settings {
            product_min: 150,
            product_max: 200,
            num_return: 3,
            ..Settings::default()
        };
        let pairs = design(&template(), Some((480, 40)), &s);
        assert!(!pairs.is_empty() && pairs.len() <= 3);
        for p in &pairs {
            assert!((150..=200).contains(&p.product_size));
            assert!(p.left.end <= 480 && p.right.start >= 520);
            assert!((s.min_tm..=s.max_tm).contains(&p.left.tm));
            assert!(longest_run(&p.left.seq) <= s.max_poly_x);
        }
        assert!(pairs.windows(2).all(|w| w[0].penalty <= w[1].penalty));
        assert!(design(b"ACGT", None, &Settings::default()).is_empty());
    }
}
