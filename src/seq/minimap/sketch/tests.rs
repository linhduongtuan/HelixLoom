use super::*;
use crate::seq::dna::revcom;

fn random(n: usize, seed: u32) -> Vec<u8> {
    let mut x = seed.wrapping_mul(2_654_435_761).wrapping_add(11);
    (0..n)
        .map(|_| {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((x >> 16) & 3) as usize]
        })
        .collect()
}

/// The same thing written out plainly: every window of `w` k-mers of a run, every k-mer of the least hash in it.
fn plain(seq: &[u8], k: usize, w: usize) -> Vec<Minimizer> {
    let mut found: Vec<Minimizer> = Vec::new();
    let mut start = 0;
    while start < seq.len() {
        if code(seq[start]).is_none() {
            start += 1;
            continue;
        }
        let end = seq[start..]
            .iter()
            .position(|&b| code(b).is_none())
            .map_or(seq.len(), |n| start + n);
        if end - start >= k {
            let all = hashed(&seq[start..end], k);
            let windows: Vec<std::ops::Range<usize>> = if all.len() < w {
                std::iter::once(0..all.len()).collect()
            } else {
                (0..=all.len() - w).map(|i| i..i + w).collect()
            };
            for win in windows {
                let best = all[win.clone()].iter().flatten().map(|e| e.0).min();
                for i in win {
                    if let (Some(b), Some((h, rev))) = (best, all[i]) {
                        if h == b {
                            found.push(Minimizer {
                                hash: h,
                                pos: start + i,
                                rev,
                            });
                        }
                    }
                }
            }
        }
        start = end;
    }
    found.sort_by_key(|m| m.pos);
    found.dedup();
    found
}

#[test]
fn the_sketch_is_the_plain_one_for_every_k_and_w_with_breaks_and_lower_case() {
    for seed in 0..40 {
        let mut seq = random(300 + 7 * seed as usize, seed);
        // a break, a lower-case stretch, an RNA stretch, a short run and an ambiguous base
        seq[100] = b'N';
        seq[104] = b'N';
        seq[150..180].make_ascii_lowercase();
        seq[200] = b'U';
        for (k, w) in [
            (1, 1),
            (3, 1),
            (5, 4),
            (11, 7),
            (15, 10),
            (21, 5),
            (31, 20),
            (32, 3),
            (8, 400),
        ] {
            assert_eq!(
                sketch(&seq, k, w),
                plain(&seq, k, w),
                "seed {seed} k {k} w {w}"
            );
        }
    }
}

#[test]
fn every_window_has_a_minimizer_and_the_minimizers_lie_in_order() {
    let seq = random(2000, 5);
    for (k, w) in [(15, 10), (11, 5), (21, 30)] {
        let m = sketch(&seq, k, w);
        assert!(m.windows(2).all(|p| p[0].pos < p[1].pos));
        // each window of w k-mers has one of them inside
        for start in 0..=(seq.len() - k + 1).saturating_sub(w) {
            assert!(
                m.iter().any(|x| x.pos >= start && x.pos < start + w),
                "k {k} w {w} window {start}"
            );
        }
        // and they are about 2 / (w + 1) of the k-mers
        let share = m.len() as f64 / (seq.len() - k + 1) as f64;
        assert!(
            (share - 2.0 / (w as f64 + 1.0)).abs() < 0.5 * 2.0 / (w as f64 + 1.0),
            "{share}"
        );
    }
}

#[test]
fn the_reverse_complement_has_the_same_minimizers_mirrored_with_the_other_strand() {
    for seed in 0..30 {
        let mut seq = random(400, seed);
        seq[50] = b'N';
        let back = revcom(&seq);
        for (k, w) in [(15, 10), (9, 4), (21, 15), (12, 700)] {
            let (a, b) = (sketch(&seq, k, w), sketch(&back, k, w));
            let mirrored: Vec<(u64, usize, bool)> = {
                let mut m: Vec<_> = b
                    .iter()
                    .map(|x| (x.hash, seq.len() - k - x.pos, !x.rev))
                    .collect();
                m.sort_by_key(|e| e.1);
                m
            };
            let direct: Vec<(u64, usize, bool)> =
                a.iter().map(|x| (x.hash, x.pos, x.rev)).collect();
            assert_eq!(direct, mirrored, "seed {seed} k {k} w {w}");
        }
    }
}

#[test]
fn a_k_mer_that_is_its_own_reverse_complement_is_not_taken_and_short_sequences_give_nothing() {
    // ACGT is its own reverse complement: with k = 4 there is nothing to take
    assert!(sketch(b"ACGT", 4, 1).is_empty());
    assert!(sketch(b"ACGTACGT", 4, 1).iter().all(|m| m.pos != 0));
    assert!(sketch(b"ACG", 4, 1).is_empty());
    assert!(sketch(b"", 15, 10).is_empty());
    assert!(sketch(b"NNNNNNNNNNNNNNNNNNNN", 5, 3).is_empty());
    // a run shorter than a window gives its smallest k-mer
    assert_eq!(
        sketch(&random(20, 1), 15, 10).len(),
        plain(&random(20, 1), 15, 10).len()
    );
    assert!(!sketch(&random(20, 1), 15, 10).is_empty());
}

#[test]
fn the_hash_is_one_to_one_on_a_sample_and_the_strand_follows_the_smaller_number() {
    let mut seen = std::collections::HashSet::new();
    assert!((0..100_000u64).all(|x| seen.insert(mix(x))));
    // AAAAC against its reverse complement GTTTT: the smaller number is AAAAC, so the forward k-mer is canonical
    let fwd = sketch(b"AAAAC", 5, 1);
    let rev = sketch(b"GTTTT", 5, 1);
    assert_eq!((fwd[0].rev, rev[0].rev), (false, true));
    assert_eq!(fwd[0].hash, rev[0].hash);
}
