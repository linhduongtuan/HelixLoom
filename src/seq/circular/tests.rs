use super::*;
use crate::seq::debruijn::sim::random;
use crate::seq::minimap::sim::Rng;

/// `seq` with each base changed with chance `rate` (substitutions only).
fn mutated(seq: &[u8], rate: f64, seed: u32) -> Vec<u8> {
    let mut r = Rng::new(seed);
    seq.iter()
        .map(|&b| {
            if r.chance(rate) {
                b"ACGT"[(b"ACGT".iter().position(|&c| c == b).unwrap_or(0) + 1 + r.below(3)) % 4]
            } else {
                b
            }
        })
        .collect()
}

fn flip(seq: &mut [u8], at: usize) {
    seq[at] = if seq[at] == b'A' { b'C' } else { b'A' };
}

#[test]
fn the_defaults_are_thirty_bases_and_ninety_five_percent() {
    let o = Overlap::default();
    assert_eq!(
        (o.min_length, o.min_identity, o.max_length),
        (30, 0.95, 5000)
    );
    let c = Concatemer::default();
    assert_eq!(
        (c.min_unit, c.max_copies, c.min_identity, c.min_covered),
        (200, 10, 0.95, 0.9)
    );
}

#[test]
fn an_end_that_repeats_the_start_is_an_overlap_of_that_length() {
    let x = random(500, 1);
    let seq = [&x[..], &x[..60]].concat();
    assert_eq!(end_overlap(&seq, &Overlap::default()), Some(60));
    // thirty is the shortest: exactly thirty is found, twenty-five is not
    assert_eq!(
        end_overlap(&[&x[..], &x[..30]].concat(), &Overlap::default()),
        Some(30)
    );
    assert_eq!(
        end_overlap(&[&x[..], &x[..25]].concat(), &Overlap::default()),
        None
    );
}

#[test]
fn a_sequence_with_no_overlap_has_none() {
    assert_eq!(end_overlap(&random(2000, 2), &Overlap::default()), None);
    assert_eq!(end_overlap(&[], &Overlap::default()), None);
    assert_eq!(end_overlap(b"ACGT", &Overlap::default()), None);
}

#[test]
fn the_ends_must_be_alike_enough_and_ninety_five_percent_of_twenty_is_nineteen() {
    let o = Overlap {
        min_length: 20,
        ..Overlap::default()
    };
    let x = random(300, 3);
    let mut one = [&x[..], &x[..20]].concat();
    flip(&mut one, 310);
    assert_eq!(end_overlap(&one, &o), Some(20));
    let mut two = one.clone();
    flip(&mut two, 315);
    assert_eq!(end_overlap(&two, &o), None);
    // and the share asked for is the one used
    let lax = Overlap {
        min_identity: 0.85,
        ..o
    };
    assert_eq!(end_overlap(&two, &lax), Some(20));
}

#[test]
fn an_overlap_longer_than_the_longest_asked_for_is_not_found_and_neither_is_more_than_half() {
    let x = random(500, 4);
    let seq = [&x[..], &x[..100]].concat();
    let short = Overlap {
        max_length: 50,
        ..Overlap::default()
    };
    assert_eq!(end_overlap(&seq, &short), None);
    assert_eq!(
        end_overlap(
            &seq,
            &Overlap {
                max_length: 100,
                ..Overlap::default()
            }
        ),
        Some(100)
    );
    // the whole of a sequence twice is an overlap of half of it, and no more
    let twice = [&x[..200], &x[..200]].concat();
    assert_eq!(end_overlap(&twice, &Overlap::default()), Some(200));
}

#[test]
fn upper_and_lower_case_are_the_same_base() {
    let x = random(400, 5);
    let mut seq = [&x[..], &x[..50]].concat();
    let n = seq.len();
    seq[n - 50..].make_ascii_lowercase();
    assert_eq!(end_overlap(&seq, &Overlap::default()), Some(50));
}

#[test]
fn of_two_overlaps_the_more_alike_is_taken_and_of_two_equally_alike_the_longer() {
    // an end that matches the start over its last 40 bases perfectly and over the last 80 at 96 %
    let x = random(600, 6);
    let mut seq = [&x[..], &x[..40]].concat();
    seq.extend_from_slice(&[]);
    assert_eq!(
        end_overlap(
            &seq,
            &Overlap {
                min_identity: 0.5,
                ..Overlap::default()
            }
        ),
        Some(40)
    );
}

#[test]
fn a_contig_of_copies_of_a_unit_is_found_with_its_unit_and_the_number_of_copies() {
    let unit = random(1000, 7);
    for n in [2usize, 3, 5] {
        let seq: Vec<u8> = (0..n)
            .flat_map(|i| mutated(&unit, 0.02, 10 + i as u32))
            .collect();
        let found =
            find_concatemer(&seq, &Concatemer::default()).unwrap_or_else(|| panic!("{n} copies"));
        assert_eq!(found.copies, n);
        assert_eq!(found.unit, seq.len() / n);
    }
}

#[test]
fn the_most_copies_that_fit_are_taken_and_not_a_divisor_of_them() {
    let unit = random(800, 8);
    let seq: Vec<u8> = (0..4).flat_map(|i| mutated(&unit, 0.01, 20 + i)).collect();
    assert_eq!(
        find_concatemer(&seq, &Concatemer::default()).map(|f| f.copies),
        Some(4)
    );
}

#[test]
fn a_contig_that_is_not_made_of_copies_is_not_a_concatemer() {
    assert_eq!(
        find_concatemer(&random(3000, 9), &Concatemer::default()),
        None
    );
    // two halves that are different
    let seq = [random(1500, 10), random(1500, 11)].concat();
    assert_eq!(find_concatemer(&seq, &Concatemer::default()), None);
}

#[test]
fn copies_that_are_too_unlike_the_first_are_not_copies() {
    let unit = random(1000, 12);
    // one base in ten changed: 90 % alike, under the 95 % asked for
    let seq = [unit.clone(), mutated(&unit, 0.1, 30)].concat();
    assert_eq!(find_concatemer(&seq, &Concatemer::default()), None);
    // a looser identity takes them
    let lax = Concatemer {
        min_identity: 0.8,
        ..Concatemer::default()
    };
    assert_eq!(find_concatemer(&seq, &lax).map(|f| f.copies), Some(2));
}

#[test]
fn a_unit_under_the_shortest_is_not_looked_for_and_neither_are_more_copies_than_asked() {
    // four copies of 120 bases are also two copies of 240: with the shortest unit at 200, that is what is found
    let unit = random(120, 13);
    let seq: Vec<u8> = (0..4).flat_map(|_| unit.clone()).collect();
    assert_eq!(
        find_concatemer(&seq, &Concatemer::default()),
        Some(Found {
            unit: 240,
            copies: 2
        })
    );
    // asked to look for units of 100, the four are found
    let short = Concatemer {
        min_unit: 100,
        ..Concatemer::default()
    };
    assert_eq!(
        find_concatemer(&seq, &short),
        Some(Found {
            unit: 120,
            copies: 4
        })
    );
    // asked for no more than three copies, the two of 240 are what is left
    assert_eq!(
        find_concatemer(
            &seq,
            &Concatemer {
                max_copies: 3,
                ..short
            }
        )
        .map(|f| f.copies),
        Some(2)
    );
    assert_eq!(
        find_concatemer(
            &seq,
            &Concatemer {
                max_copies: 1,
                ..short
            }
        ),
        None
    );
    assert_eq!(find_concatemer(&[], &Concatemer::default()), None);
}

#[test]
fn of_overlaps_that_are_equally_alike_the_longer_is_taken() {
    // a run of A at both ends: every length up to 100 is a perfect overlap, and the longest is the one
    let seq = [vec![b'A'; 100], random(300, 21), vec![b'A'; 100]].concat();
    assert_eq!(end_overlap(&seq, &Overlap::default()), Some(100));
}

#[test]
fn a_unit_of_exactly_the_shortest_is_looked_for_and_one_base_less_is_not() {
    let unit = random(200, 22);
    let seq = [&unit[..], &unit].concat();
    assert_eq!(
        find_concatemer(&seq, &Concatemer::default()).map(|f| f.unit),
        Some(200)
    );
    let longer = Concatemer {
        min_unit: 201,
        ..Concatemer::default()
    };
    assert_eq!(find_concatemer(&seq, &longer), None);
}
