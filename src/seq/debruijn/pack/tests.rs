use super::*;
use crate::seq::debruijn::sim::random;
use crate::seq::dna::revcom;

#[test]
fn k_must_be_odd_and_between_three_and_sixty_three() {
    assert!(Pack::new(21).is_ok() && Pack::new(3).is_ok() && Pack::new(63).is_ok());
    for k in [0, 1, 2, 4, 20, 64, 65, 127] {
        assert!(Pack::new(k).is_err(), "k = {k}");
    }
    assert!(Pack::new(22).unwrap_err().contains("odd"));
}

#[test]
fn a_k_mer_packs_and_unpacks_and_a_base_that_is_not_known_does_not_pack() {
    let p = Pack::new(7).unwrap();
    let kmer = p.pack(b"ACGTACG").unwrap();
    assert_eq!(p.unpack(kmer), b"ACGTACG");
    // A is 0, C 1, G 2, T 3: two bits a base, the first base highest
    assert_eq!(p.pack(b"AAAAAAC"), Some(1));
    assert_eq!(p.pack(b"CAAAAAA"), Some(1 << 12));
    assert_eq!(p.pack(b"acgtacg"), Some(kmer));
    assert_eq!(p.pack(b"ACGUACG"), p.pack(b"ACGTACG"));
    assert_eq!(p.pack(b"ACGNACG"), None);
    assert_eq!(p.pack(b"ACGTAC"), None);
    // only the first k bases are packed
    assert_eq!(p.pack(b"ACGTACGTTT"), Some(kmer));
}

#[test]
fn the_reverse_complement_is_what_revcom_gives_for_every_length() {
    for k in [3, 5, 21, 31, 33, 55, 63] {
        let p = Pack::new(k).unwrap();
        let seq = random(k + 40, k as u32);
        for window in seq.windows(k) {
            let kmer = p.pack(window).unwrap();
            assert_eq!(p.unpack(p.rc(kmer)), revcom(window), "k = {k}");
            assert_eq!(p.rc(p.rc(kmer)), kmer);
        }
    }
}

#[test]
fn a_k_mer_and_its_reverse_complement_have_one_canonical_form_the_smaller() {
    let p = Pack::new(9).unwrap();
    let a = p.pack(b"AACGTTCCA").unwrap();
    let b = p.rc(a);
    assert_ne!(a, b);
    assert_eq!(p.canonical(a), p.canonical(b));
    assert_eq!(p.canonical(a), a.min(b));
}

#[test]
fn adding_a_base_at_the_end_and_at_the_start_move_the_window() {
    let p = Pack::new(5).unwrap();
    let x = p.pack(b"ACGTA").unwrap();
    assert_eq!(p.unpack(p.append(x, 2)), b"CGTAG");
    assert_eq!(p.unpack(p.prepend(x, 3)), b"TACGT");
    // they undo each other when the base that fell off is put back
    assert_eq!(p.prepend(p.append(x, 2), 0), x);
    assert_eq!(p.append(p.prepend(x, 3), 0), x);
}

#[test]
fn the_windows_of_a_sequence_skip_what_has_an_unknown_base_and_come_with_their_reverse_complements()
{
    let p = Pack::new(5).unwrap();
    let seq = b"ACGTACGNACGTTTGCA";
    let mut seen = Vec::new();
    p.for_each(seq, |f, r| seen.push((p.unpack(f), p.unpack(r))));
    let expected: Vec<Vec<u8>> = seq
        .windows(5)
        .filter(|w| !w.contains(&b'N'))
        .map(|w| w.to_vec())
        .collect();
    assert_eq!(
        seen.iter().map(|s| s.0.clone()).collect::<Vec<_>>(),
        expected
    );
    assert!(seen.iter().all(|(f, r)| revcom(f) == *r));
    assert_eq!(seen.len(), 3 + 5);
    // too short: none
    let mut n = 0;
    p.for_each(b"ACGT", |_, _| n += 1);
    assert_eq!(n, 0);
}

#[test]
fn the_count_gives_both_strands_one_k_mer_and_equals_a_plain_count() {
    let p = Pack::new(11).unwrap();
    let genome = random(300, 5);
    let pieces = vec![
        genome[..200].to_vec(),
        revcom(&genome[100..]),
        genome[50..250].to_vec(),
    ];
    let counted = count(&pieces, p);
    let mut plain: std::collections::HashMap<u128, f32> = std::collections::HashMap::new();
    for piece in &pieces {
        for w in piece.windows(11) {
            let k = p.canonical(p.pack(w).unwrap());
            *plain.entry(k).or_insert(0.0) += 1.0;
        }
    }
    assert_eq!(counted.len(), plain.len());
    assert!(plain.iter().all(|(k, n)| counted[k] == *n));
    // a stretch that is in all three pieces is seen three times, once on the other strand
    let shared = p.canonical(p.pack(&genome[120..131]).unwrap());
    assert_eq!(counted[&shared], 3.0);
    assert!(count(&[], p).is_empty());
}

#[test]
fn many_pieces_are_counted_the_same_whatever_the_cores_share_out() {
    let p = Pack::new(15).unwrap();
    let genome = random(2000, 6);
    let pieces: Vec<Vec<u8>> = (0..300)
        .map(|i| genome[i * 3..i * 3 + 100].to_vec())
        .collect();
    let counted = count(&pieces, p);
    let mut plain: std::collections::HashMap<u128, f32> = std::collections::HashMap::new();
    for piece in &pieces {
        p.for_each(piece, |f, r| *plain.entry(f.min(r)).or_insert(0.0) += 1.0);
    }
    assert_eq!(counted.len(), plain.len());
    assert!(plain.iter().all(|(k, n)| counted[k] == *n));
}
