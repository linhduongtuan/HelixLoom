use super::*;
use crate::seq::debruijn::graph::{rc, seg};
use crate::seq::debruijn::pack::count;
use crate::seq::debruijn::sim::random;
use crate::seq::dna::revcom;

fn graph_of(pieces: &[Vec<u8>], k: usize) -> (Graph, Pack) {
    let pack = Pack::new(k).unwrap();
    (unitigs(&count(pieces, pack), pack), pack)
}

/// Every k-mer of every segment, canonical.
fn kmers_of(g: &Graph, pack: Pack) -> Vec<Kmer> {
    let mut out = Vec::new();
    for id in g.ids() {
        pack.for_each(&g.segment(id).seq, |f, r| out.push(f.min(r)));
    }
    out.sort_unstable();
    out
}

#[test]
fn a_sequence_without_repeats_is_one_segment_whichever_strand_it_is_read_from() {
    let seq = random(300, 1);
    let (g, _) = graph_of(std::slice::from_ref(&seq), 21);
    assert_eq!(g.count(), 1);
    let s = &g.segment(0).seq;
    assert!(*s == seq || *s == revcom(&seq));
    assert_eq!(g.weight_total(), 280.0);
    let (g2, _) = graph_of(&[revcom(&seq)], 21);
    assert_eq!(g2.count(), 1);
    assert!(g2.segment(0).seq == seq || g2.segment(0).seq == revcom(&seq));
}

#[test]
fn reads_of_one_sequence_on_both_strands_add_their_counts_to_one_segment() {
    let seq = random(200, 2);
    let (g, _) = graph_of(&[seq.clone(), revcom(&seq), seq[50..150].to_vec()], 15);
    assert_eq!(g.count(), 1);
    // the k-mers of the middle are seen three times, the ends twice
    let cov = g.coverage(0);
    assert!(cov > 2.0 && cov < 3.0, "{cov}");
}

#[test]
fn a_repeat_makes_branches_and_every_k_mer_is_in_a_segment_once() {
    let (x, r, y, z) = (
        random(100, 3),
        random(14, 4),
        random(100, 5),
        random(100, 6),
    );
    let seq: Vec<u8> = [&x[..], &r, &y, &r, &z].concat();
    let (g, pack) = graph_of(std::slice::from_ref(&seq), 7);
    assert!(g.count() > 1);
    let mut want: Vec<Kmer> = Vec::new();
    pack.for_each(&seq, |f, rv| want.push(f.min(rv)));
    want.sort_unstable();
    want.dedup();
    assert_eq!(kmers_of(&g, pack), want);
    // every link joins segments whose ends share k-1 bases, in the direction of the link
    assert!(!g.links().is_empty());
    for (a, b) in g.links() {
        let (oa, ob) = (g.oriented(a), g.oriented(b));
        assert_eq!(oa[oa.len() - 6..], ob[..6], "link {a} {b}");
        assert!(g.succ(rc(b)).contains(&rc(a)));
    }
    // no two segments can be joined without a branch: the graph is compact
    assert_eq!(g.compact().count(), g.count());
}

#[test]
fn a_circle_is_one_segment_that_is_k_minus_one_bases_longer_than_its_k_mers_and_leads_to_itself() {
    let circle = random(120, 7);
    let k = 15;
    let round: Vec<u8> = circle.iter().chain(&circle[..k - 1]).copied().collect();
    let (g, _) = graph_of(&[round], k);
    assert_eq!(g.count(), 1);
    assert!(g.circular(0));
    assert_eq!(g.length(0), 120 + k - 1);
    assert_eq!(g.kmers(0), 120);
    // its end and its start are the same k-1 bases
    let s = &g.segment(0).seq;
    assert_eq!(s[..k - 1], s[s.len() - (k - 1)..]);
}

#[test]
fn sequences_that_share_nothing_are_separate_segments_with_no_links() {
    let (g, _) = graph_of(&[random(150, 8), random(150, 9)], 21);
    assert_eq!(g.count(), 2);
    assert!(g.links().is_empty());
}

#[test]
fn a_difference_of_one_base_makes_a_bubble_of_two_paths_between_the_same_ends() {
    let a = random(120, 10);
    let mut b = a.clone();
    b[60] = if b[60] == b'A' { b'C' } else { b'A' };
    let (g, _) = graph_of(&[a, b], 15);
    // before, the two paths, after: 4 segments
    assert_eq!(g.count(), 4);
    let mut lengths: Vec<usize> = g.ids().into_iter().map(|i| g.length(i)).collect();
    lengths.sort_unstable();
    // each path is the k k-mers that hold the changed base: 2k-1 bases
    assert_eq!(lengths[0..2], [29, 29]);
    let branching = g
        .ids()
        .into_iter()
        .flat_map(|i| [handle(i, false), handle(i, true)])
        .filter(|&h| g.outdeg(h) == 2)
        .count();
    assert_eq!(branching, 2);
}

#[test]
fn the_segment_of_a_handle_is_found_from_the_handle() {
    let (g, _) = graph_of(&[random(100, 11)], 21);
    assert_eq!(seg(handle(0, true)), 0);
    assert_eq!(g.count(), 1);
}

#[test]
fn a_hairpin_is_not_walked_round_for_ever() {
    // a sequence followed by its own reverse complement: its k-mers are in the graph twice over, as one
    let a = random(60, 12);
    let hairpin: Vec<u8> = a.iter().copied().chain(revcom(&a)).collect();
    let (g, pack) = graph_of(&[hairpin], 11);
    assert!(g.count() >= 1);
    let kmers = kmers_of(&g, pack);
    let mut distinct = kmers.clone();
    distinct.dedup();
    assert_eq!(kmers, distinct);
}
