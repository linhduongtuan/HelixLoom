use super::*;
use crate::seq::debruijn::graph::rc;
use crate::seq::debruijn::sim::random;

const K: usize = 5;

fn settings() -> Settings {
    Settings {
        tip_length: 100,
        tip_ratio: 2.5,
        bulge_length: 150,
        bulge_diff: 0.1,
    }
}

/// A segment of `len` bases (made from `seed`) seen `cov` times on average.
fn add(g: &mut Graph, len: usize, cov: f64, seed: u32) -> usize {
    g.add(random(len, seed), cov * (len + 1 - K) as f64)
}

fn f(id: usize) -> Handle {
    handle(id, false)
}

fn r(id: usize) -> Handle {
    handle(id, true)
}

/// A junction `a`, two ways on: `b`, which goes on to a long `c`, and `tip`, which goes nowhere.
fn forked(tip_len: usize, tip_cov: f64, b_cov: f64) -> (Graph, usize, usize, usize) {
    let mut g = Graph::new(K);
    let a = add(&mut g, 60, 40.0, 1);
    let b = add(&mut g, 60, b_cov, 2);
    let tip = add(&mut g, tip_len, tip_cov, 3);
    let c = add(&mut g, 200, 40.0, 4);
    g.add_link(f(a), f(b));
    g.add_link(f(a), f(tip));
    g.add_link(f(b), f(c));
    (g, a, b, tip)
}

#[test]
fn a_short_dead_end_with_little_coverage_beside_a_way_on_is_a_tip_and_is_removed() {
    let (mut g, a, b, tip) = forked(12, 2.0, 40.0);
    assert_eq!(remove_tips(&mut g, &settings()), 1);
    assert!(!g.alive(tip) && g.alive(a) && g.alive(b));
    assert_eq!(g.succ(f(a)), [f(b)]);
}

#[test]
fn a_tip_with_more_than_the_ratio_times_the_coverage_beside_it_is_kept() {
    // 2.5 times 40 is 100: a tip seen 100 times is kept, one seen 101 times too, one seen 99 times goes
    let (mut g, _, _, tip) = forked(12, 101.0, 40.0);
    assert_eq!(remove_tips(&mut g, &settings()), 0);
    assert!(g.alive(tip));
    let (mut g, _, _, tip) = forked(12, 99.0, 40.0);
    assert_eq!(remove_tips(&mut g, &settings()), 1);
    assert!(!g.alive(tip));
    // the limit itself goes: not more than the ratio
    let (mut g, _, _, tip) = forked(12, 100.0, 40.0);
    assert_eq!(remove_tips(&mut g, &settings()), 1);
    assert!(!g.alive(tip));
}

#[test]
fn a_tip_is_judged_against_the_best_covered_way_on() {
    let mut g = Graph::new(K);
    let a = add(&mut g, 60, 40.0, 1);
    let weak = add(&mut g, 60, 1.0, 2);
    let strong = add(&mut g, 60, 40.0, 3);
    let tip = add(&mut g, 12, 90.0, 4);
    let end = add(&mut g, 200, 40.0, 5);
    for x in [weak, strong, tip] {
        g.add_link(f(a), f(x));
    }
    for x in [weak, strong] {
        g.add_link(f(x), f(end));
    }
    // 90 is under 2.5 times 40, though over 2.5 times 1
    assert_eq!(remove_tips(&mut g, &settings()), 1);
    assert!(!g.alive(tip) && g.alive(weak) && g.alive(strong));
}

#[test]
fn a_dead_end_that_is_too_long_is_not_a_tip_and_the_length_is_in_k_mers() {
    // of K = 5: 105 bases are 101 k-mers, 104 are 100, which is the longest tip
    let (mut g, _, _, tip) = forked(105, 2.0, 40.0);
    assert_eq!(remove_tips(&mut g, &settings()), 0);
    assert!(g.alive(tip));
    let (mut g, _, _, tip) = forked(104, 2.0, 40.0);
    assert_eq!(remove_tips(&mut g, &settings()), 1);
    assert!(!g.alive(tip));
}

#[test]
fn the_start_of_a_path_that_leads_into_a_junction_is_a_tip_too() {
    // the tip leads into `b`, as `a` does: it has nothing before it
    let mut g = Graph::new(K);
    let a = add(&mut g, 60, 40.0, 1);
    let b = add(&mut g, 60, 40.0, 2);
    let tip = add(&mut g, 12, 2.0, 3);
    g.add_link(f(a), f(b));
    g.add_link(f(tip), f(b));
    assert_eq!(remove_tips(&mut g, &settings()), 1);
    assert!(!g.alive(tip));
    assert_eq!(g.pred(f(b)), [f(a)]);
}

#[test]
fn a_segment_with_nothing_on_either_side_or_a_single_way_is_not_a_tip() {
    let mut g = Graph::new(K);
    let lone = add(&mut g, 12, 2.0, 1);
    let a = add(&mut g, 60, 40.0, 2);
    let end = add(&mut g, 12, 2.0, 3);
    g.add_link(f(a), f(end));
    // `a` has one way on only: removing `end` would leave `a` dead
    assert_eq!(remove_tips(&mut g, &settings()), 0);
    assert!(g.alive(lone) && g.alive(end));
}

#[test]
fn tips_are_removed_shortest_first_and_each_is_judged_with_what_is_left() {
    // two dead ends beside one way on: both go
    let mut g = Graph::new(K);
    let a = add(&mut g, 60, 40.0, 1);
    let b = add(&mut g, 60, 40.0, 2);
    let c = add(&mut g, 200, 40.0, 5);
    let t1 = add(&mut g, 10, 2.0, 3);
    let t2 = add(&mut g, 20, 2.0, 4);
    for x in [b, t1, t2] {
        g.add_link(f(a), f(x));
    }
    g.add_link(f(b), f(c));
    assert_eq!(remove_tips(&mut g, &settings()), 2);
    assert_eq!(g.ids(), [a, b, c]);
    // only dead ends beside each other: the first is removed (the other is a way on), the last is then alone and stays
    let mut g = Graph::new(K);
    let a = add(&mut g, 60, 40.0, 1);
    let t1 = add(&mut g, 10, 2.0, 3);
    let t2 = add(&mut g, 20, 2.0, 4);
    g.add_link(f(a), f(t1));
    g.add_link(f(a), f(t2));
    assert_eq!(remove_tips(&mut g, &settings()), 1);
    assert!(!g.alive(t1) && g.alive(t2));
}

#[test]
fn an_infinite_ratio_removes_every_short_tip() {
    let (mut g, _, _, tip) = forked(12, 1000.0, 1.0);
    let s = Settings {
        tip_ratio: f64::INFINITY,
        ..settings()
    };
    assert_eq!(remove_tips(&mut g, &s), 1);
    assert!(!g.alive(tip));
}

/// `a`, then two paths of the given bases to `b`.
fn bubble(x: &[u8], y: &[u8], cx: f64, cy: f64) -> (Graph, usize, usize, usize, usize) {
    let mut g = Graph::new(K);
    let a = add(&mut g, 60, 40.0, 1);
    let b = add(&mut g, 60, 40.0, 2);
    let px = g.add(x.to_vec(), cx * (x.len() + 1 - K) as f64);
    let py = g.add(y.to_vec(), cy * (y.len() + 1 - K) as f64);
    for p in [px, py] {
        g.add_link(f(a), f(p));
        g.add_link(f(p), f(b));
    }
    (g, a, b, px, py)
}

fn changed(seq: &[u8], at: usize) -> Vec<u8> {
    let mut s = seq.to_vec();
    s[at] = if s[at] == b'A' { b'C' } else { b'A' };
    s
}

#[test]
fn of_two_paths_that_spell_nearly_the_same_the_less_covered_goes_and_the_other_gains_its_coverage()
{
    let x = random(40, 5);
    let (mut g, _, _, px, py) = bubble(&x, &changed(&x, 20), 30.0, 2.0);
    assert_eq!(pop_bulges(&mut g, &settings()), 1);
    assert!(g.alive(px) && !g.alive(py));
    assert_eq!(g.coverage(px), 32.0);
    // either way round
    let (mut g, _, _, px, py) = bubble(&x, &changed(&x, 20), 2.0, 30.0);
    assert_eq!(pop_bulges(&mut g, &settings()), 1);
    assert!(!g.alive(px) && g.alive(py));
    assert_eq!(g.coverage(py), 32.0);
}

#[test]
fn a_bulge_goes_when_it_is_the_shorter_of_paths_that_differ_by_a_base_put_in() {
    let x = random(40, 6);
    let mut y = x.clone();
    y.insert(20, b'G');
    let (mut g, _, _, px, py) = bubble(&x, &y, 3.0, 30.0);
    assert_eq!(pop_bulges(&mut g, &settings()), 1);
    assert!(!g.alive(px) && g.alive(py));
}

#[test]
fn paths_that_differ_by_more_than_the_share_allowed_stay() {
    let x = random(40, 7);
    // 4 of 40 bases is the tenth; 5 is over it
    let mut four = x.clone();
    for at in [5, 15, 25, 35] {
        four = changed(&four, at);
    }
    let (mut g, _, _, _, py) = bubble(&x, &four, 30.0, 2.0);
    assert_eq!(pop_bulges(&mut g, &settings()), 1);
    assert!(!g.alive(py));
    let mut five = four.clone();
    five = changed(&five, 30);
    let (mut g, _, _, px, py) = bubble(&x, &five, 30.0, 2.0);
    assert_eq!(pop_bulges(&mut g, &settings()), 0);
    assert!(g.alive(px) && g.alive(py));
}

#[test]
fn paths_that_are_not_alike_or_too_long_are_not_bulges() {
    let (mut g, ..) = bubble(&random(40, 8), &random(40, 9), 30.0, 2.0);
    assert_eq!(pop_bulges(&mut g, &settings()), 0);
    let x = random(155, 10);
    let (mut g, _, _, px, py) = bubble(&x, &changed(&x, 80), 30.0, 2.0);
    assert_eq!(pop_bulges(&mut g, &settings()), 0);
    assert!(g.alive(px) && g.alive(py));
    // the longest allowed is 150 k-mers (154 bases of K = 5), itself too
    let x = random(154, 10);
    let (mut g, _, _, _, py) = bubble(&x, &changed(&x, 80), 30.0, 2.0);
    assert_eq!(pop_bulges(&mut g, &settings()), 1);
    assert!(!g.alive(py));
}

#[test]
fn of_three_paths_the_ones_like_the_best_covered_go() {
    let x = random(40, 11);
    let mut g = Graph::new(K);
    let a = add(&mut g, 60, 40.0, 1);
    let b = add(&mut g, 60, 40.0, 2);
    let best = g.add(x.clone(), 20.0 * 36.0);
    let near = g.add(changed(&x, 10), 3.0 * 36.0);
    let other = g.add(random(40, 12), 25.0 * 36.0);
    for p in [best, near, other] {
        g.add_link(f(a), f(p));
        g.add_link(f(p), f(b));
    }
    assert_eq!(pop_bulges(&mut g, &settings()), 1);
    // `other` is covered most and unlike the others; `near` is like `best`, which has more coverage than it: it goes, not `best`
    assert!(!g.alive(near) && g.alive(best) && g.alive(other));
    assert_eq!(g.coverage(best), 23.0);
}

#[test]
fn a_path_that_other_places_lead_into_is_not_part_of_a_bulge() {
    let x = random(40, 13);
    let (mut g, _, _, px, py) = bubble(&x, &changed(&x, 20), 30.0, 2.0);
    let c = add(&mut g, 60, 40.0, 14);
    g.add_link(f(c), f(py));
    assert_eq!(pop_bulges(&mut g, &settings()), 0);
    assert!(g.alive(px) && g.alive(py));
}

#[test]
fn the_edit_distance_counts_bases_replaced_put_in_and_taken_out() {
    assert_eq!(edit_distance(b"kitten", b"sitting"), 3);
    assert_eq!(edit_distance(b"ACGT", b"ACGT"), 0);
    assert_eq!(edit_distance(b"", b"ACG"), 3);
    assert_eq!(edit_distance(b"ACG", b""), 3);
    assert_eq!(edit_distance(b"ACGT", b"AGT"), 1);
    assert_eq!(edit_distance(b"ACGT", b"ACCT"), 1);
    assert_eq!(edit_distance(b"AAAA", b"TTTT"), 4);
}

#[test]
fn what_is_seen_fewer_times_than_the_cutoff_is_removed_and_what_is_seen_so_often_is_not() {
    let mut g = Graph::new(K);
    let low = add(&mut g, 60, 2.9, 1);
    let edge = add(&mut g, 60, 3.0, 2);
    let high = add(&mut g, 60, 30.0, 3);
    g.add_link(f(low), f(high));
    assert_eq!(remove_below(&mut g, 3.0), 1);
    assert!(!g.alive(low) && g.alive(edge) && g.alive(high));
    assert!(g.succ(f(high)).is_empty() && g.pred(f(high)).is_empty());
    assert_eq!(remove_below(&mut g, 0.0), 0);
}

fn long_segments(covs: &[(usize, f64)]) -> Graph {
    let mut g = Graph::new(K);
    for (i, &(len, cov)) in covs.iter().enumerate() {
        add(&mut g, len, cov, i as u32 + 1);
    }
    g
}

#[test]
fn the_automatic_cutoff_is_a_tenth_of_the_coverage_most_bases_are_at_but_not_under_one_and_a_half()
{
    assert!((auto_cutoff(&long_segments(&[(500, 40.0)])) - 4.0).abs() < 1e-9);
    assert!((auto_cutoff(&long_segments(&[(500, 100.0)])) - 10.0).abs() < 1e-9);
    // a tenth of 10 is 1: not under 1.5
    assert!((auto_cutoff(&long_segments(&[(500, 10.0)])) - 1.5).abs() < 1e-9);
    // and not over half the coverage: a tenth of 2 is under 1.5, half of 2 is 1
    assert!((auto_cutoff(&long_segments(&[(500, 2.0)])) - 1.0).abs() < 1e-9);
    assert_eq!(auto_cutoff(&Graph::new(K)), 0.0);
}

#[test]
fn the_automatic_cutoff_goes_by_the_bases_and_not_by_the_number_of_segments() {
    // two long segments at 40 and three short pieces at 400: half of the bases are at 40 or less
    let g = long_segments(&[
        (1000, 40.0),
        (1000, 42.0),
        (8, 400.0),
        (8, 400.0),
        (8, 400.0),
    ]);
    // none of the short pieces is two k long (10): they are not counted while there are long ones
    assert!((auto_cutoff(&g) - 4.0).abs() < 1e-9);
    let only_short = long_segments(&[(8, 1000.0), (8, 20.0), (8, 20.0)]);
    // 2k is 10, nothing is so long: all count, the median (by bases) is 20
    assert!((auto_cutoff(&only_short) - 2.0).abs() < 1e-9);
}

#[test]
fn removing_leaves_the_reverse_links_of_what_is_left_alone() {
    let (mut g, a, b, _) = forked(12, 2.0, 40.0);
    remove_tips(&mut g, &settings());
    assert_eq!(g.succ(r(b)), [r(a)]);
    assert_eq!(g.succ(rc(f(b))), [r(a)]);
}

#[test]
fn a_short_segment_that_is_joined_to_nothing_is_removed_and_one_with_a_link_or_that_is_long_is_not()
{
    let mut g = Graph::new(K);
    let short = add(&mut g, 100, 2.0, 1);
    let edge = add(&mut g, 200, 2.0, 2);
    let long = add(&mut g, 300, 2.0, 3);
    let (p, q) = (add(&mut g, 100, 2.0, 4), add(&mut g, 100, 2.0, 5));
    g.add_link(f(p), f(q));
    // the circle of one segment is joined to itself
    let circle = add(&mut g, 100, 2.0, 6);
    g.add_link(f(circle), f(circle));
    assert_eq!(remove_isolated(&mut g, 200), 1);
    assert!(!g.alive(short));
    // 200 bases is not shorter than 200
    assert!(g.alive(edge) && g.alive(long) && g.alive(p) && g.alive(q) && g.alive(circle));
    assert_eq!(remove_isolated(&mut g, 0), 0);
    assert_eq!(remove_isolated(&mut g, 201), 1);
    assert!(!g.alive(edge));
}

#[test]
fn a_dead_end_that_a_segment_with_no_other_way_on_leads_into_is_not_a_tip() {
    // `tip` is reached from the junction `a` and from `q` (too long to be a tip itself), which has no other way on: removing `tip` would
    // leave `q` dead
    let mut g = Graph::new(K);
    let a = add(&mut g, 60, 40.0, 1);
    let b = add(&mut g, 60, 40.0, 2);
    let c = add(&mut g, 200, 40.0, 3);
    let tip = add(&mut g, 12, 2.0, 4);
    let q = add(&mut g, 300, 40.0, 5);
    g.add_link(f(a), f(b));
    g.add_link(f(b), f(c));
    g.add_link(f(a), f(tip));
    g.add_link(f(q), f(tip));
    assert_eq!(remove_tips(&mut g, &settings()), 0);
    assert!(g.alive(tip) && g.alive(q));
}

#[test]
fn a_segment_that_leads_to_its_own_reverse_is_not_a_tip() {
    // t leads to its reverse complement (a hairpin) and to x, which has another way in: nothing about it is a dead end to clip
    let mut g = Graph::new(K);
    let t = add(&mut g, 12, 2.0, 1);
    let x = add(&mut g, 60, 40.0, 2);
    let w = add(&mut g, 60, 40.0, 3);
    g.add_link(f(t), r(t));
    g.add_link(f(t), f(x));
    g.add_link(f(w), f(x));
    assert_eq!(remove_tips(&mut g, &settings()), 0);
    assert!(g.alive(t));
}

#[test]
fn the_share_of_bases_two_paths_may_differ_by_is_rounded_down() {
    let x = random(45, 20);
    let (mut four, mut five) = (x.clone(), x.clone());
    for at in [5, 15, 25, 35] {
        four = changed(&four, at);
        five = changed(&five, at);
    }
    five = changed(&five, 42);
    // a tenth of 45 is 4.5: four changes are allowed and five are not
    assert!(similar(&x, &four, 0.1));
    assert!(!similar(&x, &five, 0.1));
    // a path with a base put in is one change
    let mut longer = x.clone();
    longer.insert(10, b'T');
    assert!(similar(&x, &longer, 0.1));
}

#[test]
fn a_bulge_in_a_circle_whose_two_paths_both_lead_back_to_the_segment_they_left_is_removed() {
    // a, then x or y, then a again: a circle with one place where two copies differ
    let x = random(40, 21);
    let mut g = Graph::new(K);
    let a = add(&mut g, 100, 40.0, 1);
    let px = g.add(x.clone(), 30.0 * 36.0);
    let py = g.add(changed(&x, 20), 2.0 * 36.0);
    for p in [px, py] {
        g.add_link(f(a), f(p));
        g.add_link(f(p), f(a));
    }
    assert_eq!(pop_bulges(&mut g, &settings()), 1);
    assert!(g.alive(px) && !g.alive(py));
    let c = g.compact();
    assert_eq!(c.count(), 1);
    assert!(c.circular(0));
}
