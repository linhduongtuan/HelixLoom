use super::*;
use crate::seq::debruijn::graph::handle;
use crate::seq::debruijn::sim::random;
use crate::seq::hybrid::copies::{estimate, Options as CopyOptions};

const K: usize = 5;

fn f(id: usize) -> Handle {
    handle(id, false)
}

/// The line `X R Y R Z` of random bases as a graph with k = 5: segments `X+`, `Y`, `Z` that are single copy and the repeat `R` (ids 0, 1, 2,
/// 3), the links of a de Bruijn graph that cannot tell which way through the repeat is right.
fn chain() -> (Graph, Vec<u8>) {
    let (x, y, z, r) = (random(30, 1), random(30, 2), random(30, 3), random(16, 4));
    let genome = [&x[..], &r, &y, &r, &z].concat();
    let mut g = Graph::new(K);
    let seg = |g: &mut Graph, s: Vec<u8>, depth: f64| {
        let w = depth * (s.len() + 1 - K) as f64;
        g.add(s, w)
    };
    seg(&mut g, [&x[..], &r[..4]].concat(), 10.0);
    seg(&mut g, [&r[12..], &y[..], &r[..4]].concat(), 10.0);
    seg(&mut g, [&r[12..], &z[..]].concat(), 10.0);
    seg(&mut g, r.clone(), 20.0);
    for (a, b) in [(0, 3), (1, 3), (3, 1), (3, 2)] {
        g.add_link(f(a), f(b));
    }
    (g, genome)
}

fn through_r(from: usize, to: usize, quality: f64) -> Bridge {
    Bridge {
        from: f(from),
        to: f(to),
        via: Via::Route(vec![f(3)]),
        quality,
        reads: 5,
    }
}

fn copies(g: &Graph) -> Copies {
    estimate(g, &CopyOptions::default())
}

#[test]
fn the_copies_of_the_chain_are_one_one_one_and_two() {
    let (g, _) = chain();
    assert_eq!(copies(&g).count, [1, 1, 1, 2]);
}

#[test]
fn two_bridges_through_the_repeat_make_the_line_one_segment_and_use_the_repeat_up() {
    let (g, genome) = chain();
    let c = copies(&g);
    let a = apply(
        g,
        &c,
        &[through_r(0, 1, 80.0), through_r(1, 2, 70.0)],
        Mode::Normal,
    );
    assert_eq!(a.made, 2);
    assert_eq!(a.graph.count(), 1);
    let s = a.graph.segment(a.graph.ids()[0]);
    assert!(s.seq == genome || s.seq == crate::seq::dna::revcom(&genome));
    assert!(!a.graph.circular(a.graph.ids()[0]));
}

#[test]
fn the_depth_of_a_bridge_is_the_mean_of_the_two_segments_it_joins() {
    let (mut g, _) = chain();
    // z is 20 deep, y 10: the bridge between them is 15
    g.add_weight(2, 10.0 * 30.0);
    // the repeat is three copies deep, so that after one use two are left and it is not joined to its neighbours
    let c = Copies {
        single_depth: 10.0,
        count: vec![1, 1, 1, 3],
        by_depth: 0,
    };
    let a = apply(g, &c, &[through_r(1, 2, 70.0)], Mode::Normal);
    assert_eq!(a.made, 1);
    // y, the bridge and z are one: 34 + 12 + 30 k-mers of 10, 15 and 20
    let joined = a
        .graph
        .ids()
        .into_iter()
        .find(|&i| a.graph.length(i) == 80)
        .expect("y, the bridge and z are one segment");
    assert!((a.graph.coverage(joined) - (340.0 + 180.0 + 600.0) / 76.0).abs() < 1e-9);
}

#[test]
fn a_bridge_that_would_take_an_end_that_is_used_is_skipped_and_the_better_one_wins() {
    let (g, _) = chain();
    let c = copies(&g);
    // both leave the end of x; the better is first
    let good = through_r(0, 1, 80.0);
    let bad = through_r(0, 2, 20.0);
    let a = apply(g.clone(), &c, &[good.clone(), bad.clone()], Mode::Normal);
    assert_eq!(a.made, 1);
    // and the other way round, the order given is the order applied: the caller sorts
    let b = apply(g, &c, &[bad, good], Mode::Normal);
    assert_eq!(b.made, 1);
}

#[test]
fn a_bridge_between_the_two_ends_of_one_segment_is_not_made_when_they_are_one_end() {
    let (g, _) = chain();
    let c = copies(&g);
    // from the end of x to the reverse of x: the two ends are the same two ends only if from == reverse of to
    let silly = Bridge {
        from: f(0),
        to: handle(0, true),
        via: Via::Route(vec![f(3)]),
        quality: 90.0,
        reads: 5,
    };
    let a = apply(g, &c, &[silly], Mode::Normal);
    // the ends are x's end and x's end: refused
    assert_eq!(a.made, 0);
}

#[test]
fn a_bridge_of_the_reads_own_bases_is_put_between_two_segments_that_have_no_link() {
    let (x, z, gap) = (random(40, 5), random(40, 6), random(20, 7));
    let genome = [&x[..], &gap, &z].concat();
    let mut g = Graph::new(K);
    g.add(x, 36.0 * 10.0);
    g.add(z, 36.0 * 10.0);
    let c = Copies {
        single_depth: 10.0,
        count: vec![1, 1],
        by_depth: 0,
    };
    let b = Bridge {
        from: f(0),
        to: f(1),
        via: Via::Reads(gap),
        quality: 50.0,
        reads: 3,
    };
    let a = apply(g, &c, &[b], Mode::Normal);
    assert_eq!(a.made, 1);
    assert_eq!(a.graph.count(), 1);
    assert_eq!(a.graph.segment(a.graph.ids()[0]).seq, genome);
}

#[test]
fn a_bridge_that_is_a_link_takes_the_other_links_at_its_ends_away() {
    let base = random(60, 8);
    let mut g = Graph::new(K);
    let w = |len: usize| (len + 1 - K) as f64 * 10.0;
    g.add(base[0..30].to_vec(), w(30)); // x
    g.add(base[26..56].to_vec(), w(30)); // y
    g.add([&base[26..30], &random(26, 9)[..]].concat(), w(30)); // w, the other way on from x
    g.add_link(f(0), f(1));
    g.add_link(f(0), f(2));
    let c = Copies {
        single_depth: 10.0,
        count: vec![1, 1, 1],
        by_depth: 0,
    };
    let link = Bridge {
        from: f(0),
        to: f(1),
        via: Via::Link,
        quality: 60.0,
        reads: 4,
    };
    let a = apply(g.clone(), &c, std::slice::from_ref(&link), Mode::Normal);
    assert_eq!(a.made, 1);
    // x and y are one; w stands alone
    let mut lengths: Vec<usize> = a
        .graph
        .ids()
        .into_iter()
        .map(|i| a.graph.length(i))
        .collect();
    lengths.sort_unstable();
    assert_eq!(lengths, [30, 56]);
    assert!(a.graph.links().is_empty());
    // a link that is not there is not a bridge
    let mut unlinked = g;
    unlinked.remove_link(f(0), f(1));
    let b = apply(unlinked, &c, &[link], Mode::Normal);
    assert_eq!(b.made, 0);
}

/// `A R B` with the repeat in the middle of a straight line, `R` two copies deep.
fn straight() -> (Graph, Copies) {
    let g = random(44, 10);
    let mut graph = Graph::new(K);
    let w = |len: usize, depth: f64| (len + 1 - K) as f64 * depth;
    graph.add(g[0..16].to_vec(), w(16, 10.0));
    graph.add(g[12..32].to_vec(), w(20, 20.0));
    graph.add(g[28..44].to_vec(), w(16, 10.0));
    graph.add_link(f(0), f(1));
    graph.add_link(f(1), f(2));
    let c = Copies {
        single_depth: 10.0,
        count: vec![1, 2, 1],
        by_depth: 0,
    };
    (graph, c)
}

#[test]
fn bold_joins_everything_unbranched_and_the_other_modes_only_what_is_single_copy() {
    let (g, c) = straight();
    let count = |mode| apply(g.clone(), &c, &[], mode).graph.count();
    assert_eq!(count(Mode::Bold), 1);
    assert_eq!(count(Mode::Normal), 3);
    assert_eq!(count(Mode::Conservative), 3);
}

#[test]
fn single_copy_segments_that_are_unbranched_are_joined_in_every_mode() {
    let (g, mut c) = straight();
    c.count = vec![1, 1, 1];
    for mode in [Mode::Conservative, Mode::Normal, Mode::Bold] {
        assert_eq!(apply(g.clone(), &c, &[], mode).graph.count(), 1, "{mode:?}");
    }
}

#[test]
fn the_modes_ask_for_a_quality_of_25_10_and_1() {
    assert_eq!(Mode::Conservative.min_quality(), 25.0);
    assert_eq!(Mode::Normal.min_quality(), 10.0);
    assert_eq!(Mode::Bold.min_quality(), 1.0);
}

#[test]
fn one_bridge_through_a_repeat_of_two_copies_leaves_one_copy_that_the_rest_is_joined_through() {
    let (g, genome) = chain();
    let c = copies(&g);
    let a = apply(g, &c, &[through_r(0, 1, 80.0)], Mode::Normal);
    assert_eq!(a.made, 1);
    // x and y are joined by the bridge, and the repeat that is left (one copy, between y and z) is single copy now: it is the whole line
    assert_eq!(a.graph.count(), 1);
    assert_eq!(a.graph.segment(a.graph.ids()[0]).seq, genome);
}

#[test]
fn a_repeat_as_short_as_k_makes_a_bridge_of_k_bases() {
    // R is 5 bases, k = 5: the bridge is those 5 bases, one k-mer
    let (x, y, r) = (random(30, 21), random(30, 22), random(5, 23));
    let genome = [&x[..], &r, &y].concat();
    let mut g = Graph::new(K);
    let w = |len: usize, depth: f64| (len + 1 - K) as f64 * depth;
    let xs = [&x[..], &r[..4]].concat();
    let ys = [&r[1..], &y[..]].concat();
    g.add(xs.clone(), w(xs.len(), 10.0));
    g.add(ys.clone(), w(ys.len(), 10.0));
    g.add(r.clone(), w(5, 20.0));
    g.add_link(f(0), f(2));
    g.add_link(f(2), f(1));
    let c = Copies {
        single_depth: 10.0,
        count: vec![1, 1, 2],
        by_depth: 0,
    };
    let a = apply(
        g,
        &c,
        &[Bridge {
            from: f(0),
            to: f(1),
            via: Via::Route(vec![f(2)]),
            quality: 50.0,
            reads: 3,
        }],
        Mode::Normal,
    );
    assert_eq!(a.made, 1);
    // x, the bridge and y are one segment (the repeat, which has a copy left, stands alone)
    assert!(a
        .graph
        .ids()
        .into_iter()
        .any(|i| a.graph.segment(i).seq == genome));
}
