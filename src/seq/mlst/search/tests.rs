use super::*;
use crate::seq::mlst::tests::{allele, contigs, put, scheme, LENGTH};

fn window(seq: &[u8], before: usize, after: usize) -> Vec<u8> {
    let mut w = crate::seq::minimap::sim::genome(before, 71);
    w.extend_from_slice(seq);
    w.extend(crate::seq::minimap::sim::genome(after, 72));
    w
}

#[test]
fn an_allele_in_a_window_is_aligned_whole_with_the_margins_left_out() {
    let a = allele(0, 3);
    let (covered, matches, columns, bases) = aligned_over(&a, &window(&a, 40, 40)).unwrap();
    assert_eq!((covered, matches, columns), (LENGTH, LENGTH, LENGTH));
    assert_eq!(bases, a);
    // a base changed in the middle is a column that does not match, in a stretch that is still whole
    let mut changed = a.clone();
    changed[200] = if changed[200] == b'A' { b'C' } else { b'A' };
    let (covered, matches, columns, bases) = aligned_over(&a, &window(&changed, 40, 40)).unwrap();
    assert_eq!(
        (covered, matches, columns, bases),
        (LENGTH, LENGTH - 1, LENGTH, changed)
    );
    // three more bases in the contig: a gap, and columns that are more than the allele's
    let mut more = a[..250].to_vec();
    more.extend_from_slice(b"TTT");
    more.extend_from_slice(&a[250..]);
    let (covered, matches, columns, _) = aligned_over(&a, &window(&more, 40, 40)).unwrap();
    assert_eq!((covered, matches, columns), (LENGTH, LENGTH, LENGTH + 3));
}

#[test]
fn what_does_not_match_at_an_end_and_what_the_window_does_not_reach_is_not_covered() {
    let a = allele(0, 3);
    // the last 100 bases are not in the window (the contig ends)
    let (covered, matches, columns, bases) = aligned_over(&a, &window(&a[..350], 40, 0)).unwrap();
    assert_eq!((covered, matches, columns), (350, 350, 350));
    assert_eq!(bases, a[..350]);
    // a changed first base: the aligner may keep it as a column that does not match or leave it out as an overhang, which it does not
    // say in advance; either way every other base is a match
    let mut changed = a.clone();
    changed[0] = if changed[0] == b'A' { b'C' } else { b'A' };
    let mut w = window(&changed, 40, 40);
    w[39] = *b"ACGT"
        .iter()
        .find(|&&b| b != a[0] && b != changed[0])
        .unwrap();
    let (covered, matches, _, bases) = aligned_over(&a, &w).unwrap();
    assert_eq!(matches, LENGTH - 1);
    assert!(covered == LENGTH || covered == LENGTH - 1, "{covered}");
    assert!(
        bases == changed || bases == changed[1..],
        "the stretch is the contig's bases that the alignment covers"
    );
}

#[test]
fn an_exact_allele_is_found_on_either_strand_and_each_place_is_found_once() {
    let s = scheme("synth", 3, 6);
    let l = &s.loci[1];
    let a = allele(1, 4);
    for rev in [false, true] {
        let found = find(
            &Contigs::new(contigs(1, &[put(a.clone(), rev, 1)])).unwrap(),
            l,
        );
        assert_eq!(found.len(), 1, "{found:?}");
        let f = &found[0];
        assert_eq!(
            (
                f.exact.as_deref(),
                f.allele.as_str(),
                f.covered,
                f.allele_len
            ),
            (Some("4"), "4", LENGTH, LENGTH)
        );
        assert_eq!((f.identity(), f.coverage()), (1.0, 1.0));
        assert_eq!(
            f.sequence, a,
            "the stretch is written as the allele is, whichever strand"
        );
    }
    // two copies in two places: two finds
    let two = find(
        &Contigs::new(contigs(
            2,
            &[put(a.clone(), false, 0), put(a.clone(), true, 2)],
        ))
        .unwrap(),
        l,
    );
    assert_eq!(two.len(), 2);
    assert!(two.iter().all(|f| f.exact.as_deref() == Some("4")));
    // a genome with nothing of the locus has no place
    assert!(find(
        &Contigs::new(contigs(3, &[put(allele(2, 2), false, 0)])).unwrap(),
        l
    )
    .is_empty());
}

#[test]
fn the_closest_allele_of_a_novel_one_is_named_and_the_numbers_say_how_close() {
    let s = scheme("synth", 1, 6);
    let mut changed = allele(0, 5);
    for at in [60, 160, 260, 360] {
        changed[at] = if changed[at] == b'T' { b'G' } else { b'T' };
    }
    let found = find(
        &Contigs::new(contigs(4, &[put(changed.clone(), false, 0)])).unwrap(),
        &s.loci[0],
    );
    assert_eq!(found.len(), 1);
    let f = &found[0];
    assert_eq!(
        (f.exact.clone(), f.covered, f.columns, f.matches),
        (None, LENGTH, LENGTH, LENGTH - 4)
    );
    assert!(f.identity() > 0.99 && f.identity() < 1.0);
    assert_eq!(f.sequence, changed);
}

#[test]
fn identity_and_coverage_are_shares_and_an_empty_find_has_none() {
    let f = Found {
        allele: "1".into(),
        exact: None,
        covered: 90,
        allele_len: 100,
        matches: 95,
        columns: 100,
        sequence: Vec::new(),
    };
    assert_eq!((f.identity(), f.coverage()), (0.95, 0.9));
    let none = Found {
        columns: 0,
        allele_len: 0,
        covered: 0,
        matches: 0,
        ..f
    };
    assert_eq!((none.identity(), none.coverage()), (0.0, 0.0));
}
