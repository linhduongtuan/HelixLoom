use super::*;

const K: usize = 4;

/// A graph of segments given by their bases, with the links `(a, b)` between handles.
fn graph(seqs: &[&str], links: &[(Handle, Handle)]) -> Graph {
    let mut g = Graph::new(K);
    for s in seqs {
        g.add(s.as_bytes().to_vec(), (s.len() + 1 - K) as f64 * 10.0);
    }
    for &(a, b) in links {
        g.add_link(a, b);
    }
    g
}

fn f(id: usize) -> Handle {
    handle(id, false)
}

fn r(id: usize) -> Handle {
    handle(id, true)
}

#[test]
fn a_handle_is_a_segment_in_one_direction_and_its_reverse_is_the_other() {
    assert_eq!((f(3), r(3)), (6, 7));
    assert_eq!((rc(6), rc(7)), (7, 6));
    assert_eq!((seg(6), seg(7)), (3, 3));
}

#[test]
fn a_link_is_also_the_link_between_the_reverse_complements() {
    let g = graph(&["ACGTAC", "TACGGA"], &[(f(0), f(1))]);
    assert_eq!(g.succ(f(0)), [f(1)]);
    assert_eq!(g.succ(r(1)), [r(0)]);
    assert_eq!(g.pred(f(1)), [f(0)]);
    assert_eq!(g.pred(r(0)), [r(1)]);
    assert_eq!((g.outdeg(f(0)), g.indeg(f(1)), g.outdeg(f(1))), (1, 1, 0));
    assert_eq!(g.links(), [(f(0), f(1))]);
}

#[test]
fn a_link_added_twice_or_as_its_complement_is_one() {
    let mut g = graph(&["ACGTAC", "TACGGA"], &[]);
    g.add_link(f(0), f(1));
    g.add_link(f(0), f(1));
    g.add_link(r(1), r(0));
    assert_eq!(g.links().len(), 1);
    assert_eq!(g.succ(f(0)).len(), 1);
}

#[test]
fn a_link_from_a_segment_to_its_own_reverse_is_one_link() {
    let mut g = graph(&["ACGTAC"], &[]);
    g.add_link(f(0), r(0));
    assert_eq!(g.succ(f(0)), [r(0)]);
    assert_eq!(g.links(), [(f(0), r(0))]);
}

#[test]
fn removing_a_segment_takes_its_links_with_it() {
    let mut g = graph(
        &["ACGTAC", "TACGGA", "GGATTT"],
        &[(f(0), f(1)), (f(1), f(2)), (f(0), f(2))],
    );
    g.remove(1);
    assert_eq!(g.ids(), [0, 2]);
    assert_eq!(g.count(), 2);
    assert!(!g.alive(1) && g.alive(0));
    assert_eq!(g.succ(f(0)), [f(2)]);
    assert_eq!(g.succ(r(2)), [r(0)]);
    assert_eq!(g.links(), [(f(0), f(2))]);
    assert!(g.succ(f(1)).is_empty() && g.succ(r(1)).is_empty());
}

#[test]
fn coverage_is_the_weight_over_the_k_mers_and_the_bases_read_either_way() {
    let mut g = graph(&["ACGTAC"], &[]);
    assert_eq!((g.length(0), g.kmers(0)), (6, 3));
    assert_eq!(g.coverage(0), 10.0);
    g.add_weight(0, 30.0);
    assert_eq!(g.coverage(0), 20.0);
    assert_eq!(g.oriented(f(0)), b"ACGTAC");
    assert_eq!(g.oriented(r(0)), b"GTACGT");
    assert_eq!(g.total_length(), 6);
}

#[test]
fn a_run_with_no_branch_becomes_one_segment_with_the_bases_joined_and_the_weights_added() {
    // ACGTAC > TACGGA (share TAC) > GGATTT (share GGA)
    let g = graph(
        &["ACGTAC", "TACGGA", "GGATTT"],
        &[(f(0), f(1)), (f(1), f(2))],
    );
    let c = g.compact();
    assert_eq!(c.count(), 1);
    assert_eq!(c.segment(0).seq, b"ACGTACGGATTT");
    assert_eq!(c.segment(0).weight, 30.0 + 30.0 + 30.0);
    assert!(c.links().is_empty());
}

#[test]
fn a_run_is_joined_whichever_way_its_segments_are_written() {
    // the middle segment is written on the other strand
    let g = graph(
        &["ACGTAC", "TCCGTA", "GGATTT"],
        &[(f(0), r(1)), (r(1), f(2))],
    );
    let c = g.compact();
    assert_eq!(c.count(), 1);
    let seq = c.oriented(f(0));
    assert!(seq == b"ACGTACGGATTT" || seq == revcom(b"ACGTACGGATTT"));
}

#[test]
fn a_branch_stops_the_joining_and_the_links_are_kept() {
    // 0 leads to 1 and to 2; 3 leads on from 1
    let g = graph(
        &["ACGTAC", "TACGGA", "TACAAA", "GGATTT"],
        &[(f(0), f(1)), (f(0), f(2)), (f(1), f(3))],
    );
    let c = g.compact();
    // 1 and 3 are one run; 0 and 2 stand alone
    assert_eq!(c.count(), 3);
    assert_eq!(c.links().len(), 2);
    let lengths: Vec<usize> = c.ids().into_iter().map(|i| c.length(i)).collect();
    assert!(lengths.contains(&9));
}

#[test]
fn two_ways_into_a_segment_stop_the_joining_too() {
    let g = graph(
        &["ACGTAC", "GGGTAC", "TACGGA"],
        &[(f(0), f(2)), (f(1), f(2))],
    );
    assert_eq!(g.compact().count(), 3);
}

#[test]
fn a_ring_of_segments_becomes_one_segment_that_leads_to_itself() {
    // four segments, each leading to the next and the last to the first: one circle
    let g = graph(
        &["ACGTC", "GTCAA", "CAAAC", "AACGT"],
        &[(f(0), f(1)), (f(1), f(2)), (f(2), f(3)), (f(3), f(0))],
    );
    let c = g.compact();
    assert_eq!(c.count(), 1);
    assert!(c.circular(0));
    assert_eq!(c.succ(f(0)), [f(0)]);
    assert_eq!(c.length(0), 5 + 2 + 2 + 2);
    assert_eq!(c.segment(0).weight, 4.0 * 20.0);
}

#[test]
fn a_segment_that_leads_to_itself_alone_is_a_circle_and_one_with_other_links_is_not() {
    let mut g = graph(&["ACGTACG", "GTACGT"], &[]);
    g.add_link(f(0), f(0));
    assert!(g.circular(0) && !g.circular(1));
    g.add_link(f(0), f(1));
    assert!(!g.circular(0));
    // and it is not merged with what follows it
    assert_eq!(g.compact().count(), 2);
}

#[test]
fn compacting_a_graph_that_is_compact_changes_nothing() {
    let g = graph(
        &["ACGTAC", "TACGGA", "TACAAA"],
        &[(f(0), f(1)), (f(0), f(2))],
    );
    let c = g.compact();
    assert_eq!(c.count(), 3);
    assert_eq!(c.compact().links(), c.links());
    assert_eq!(c.compact().count(), 3);
}

#[test]
fn removed_segments_are_left_out_when_compacting() {
    let mut g = graph(
        &["ACGTAC", "TACGGA", "TACAAA", "GGATTT"],
        &[(f(0), f(1)), (f(0), f(2)), (f(1), f(3))],
    );
    g.remove(2);
    let c = g.compact();
    // what is left is one run
    assert_eq!(c.count(), 1);
    assert_eq!(c.segment(0).seq, b"ACGTACGGATTT");
}
