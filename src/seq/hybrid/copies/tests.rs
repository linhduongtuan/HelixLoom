use super::*;

const K: usize = 5;

fn f(id: usize) -> usize {
    handle(id, false)
}

/// A segment of `len` bases seen `depth` times on average.
fn add(g: &mut Graph, len: usize, depth: f64) -> usize {
    g.add(vec![b'A'; len], depth * (len + 1 - K) as f64)
}

/// The circle `X R Y R Z` as a graph: three single-copy segments of depth 10 and a repeat of depth 20.
fn ring() -> (Graph, [usize; 4]) {
    let mut g = Graph::new(K);
    let x = add(&mut g, 300, 10.0);
    let y = add(&mut g, 300, 10.0);
    let z = add(&mut g, 300, 10.0);
    let r = add(&mut g, 100, 20.0);
    for (a, b) in [(x, r), (y, r), (r, y), (r, z), (z, x)] {
        g.add_link(f(a), f(b));
    }
    (g, [x, y, z, r])
}

#[test]
fn single_copy_segments_have_one_copy_and_a_repeat_has_the_sum_of_what_leads_into_it() {
    let (g, [x, y, z, r]) = ring();
    let c = estimate(&g, &Options::default());
    assert_eq!(
        (c.count[x], c.count[y], c.count[z], c.count[r]),
        (1, 1, 1, 2)
    );
    assert_eq!(c.by_depth, 0);
    assert!((c.single_depth - 10.0).abs() < 1e-9);
}

#[test]
fn the_depth_of_a_single_copy_can_be_given() {
    let (g, [x, _, _, r]) = ring();
    // told that a single copy is 5 deep, the segments are 2 and 4 copies by depth, and the links do not agree
    let c = estimate(
        &g,
        &Options {
            single_depth: Some(5.0),
            ..Options::default()
        },
    );
    assert_eq!(c.single_depth, 5.0);
    assert_eq!((c.count[x], c.count[r]), (2, 4));
    assert_eq!(c.by_depth, 4);
}

#[test]
fn a_segment_with_nothing_to_go_by_but_its_depth_gets_the_nearest_whole_number() {
    let mut g = Graph::new(K);
    let a = add(&mut g, 500, 10.0);
    let lone = add(&mut g, 200, 31.0);
    let c = estimate(&g, &Options::default());
    assert_eq!((c.count[a], c.count[lone]), (1, 3));
    assert_eq!(c.by_depth, 1);
    // under half a copy is still one
    let mut g = Graph::new(K);
    let big = add(&mut g, 500, 10.0);
    let small = add(&mut g, 50, 2.0);
    let c = estimate(&g, &Options::default());
    assert_eq!((c.count[big], c.count[small]), (1, 1));
}

#[test]
fn a_seed_must_be_near_the_single_depth_and_have_one_way_on_at_each_end() {
    let o = Options::default();
    // 1.25 times the median is a seed; a hair over is not, and is then set by the stall rule (it is under 1.5)
    let mut g = Graph::new(K);
    add(&mut g, 1000, 10.0);
    let near = add(&mut g, 100, 12.5);
    let over = add(&mut g, 100, 12.6);
    let c = estimate(&g, &o);
    assert_eq!(c.count[near], 1);
    assert_eq!(c.count[over], 1);
    assert_eq!(c.by_depth, 0);
    // two ways on at an end: no seed, whatever the depth
    let mut g = Graph::new(K);
    let a = add(&mut g, 500, 10.0);
    let b = add(&mut g, 100, 10.0);
    let c2 = add(&mut g, 100, 10.0);
    g.add_link(f(a), f(b));
    g.add_link(f(a), f(c2));
    let cc = estimate(&g, &o);
    // the stall rule then sets the largest, and the others follow by it
    assert_eq!((cc.count[a], cc.count[b], cc.count[c2]), (1, 1, 1));
}

#[test]
fn the_largest_segment_that_could_be_single_copy_is_set_when_the_links_decide_nothing() {
    // p leads to q and r, which are half as deep: nothing seeds, p is the largest that could be one copy, and q, r then are too
    let mut g = Graph::new(K);
    let p = add(&mut g, 300, 10.0);
    let q = add(&mut g, 100, 5.0);
    let r = add(&mut g, 100, 5.0);
    g.add_link(f(p), f(q));
    g.add_link(f(p), f(r));
    let c = estimate(&g, &Options::default());
    assert_eq!((c.count[p], c.count[q], c.count[r]), (1, 1, 1));
    assert_eq!(c.by_depth, 0);
}

#[test]
fn a_depth_that_the_links_do_not_agree_with_is_not_taken_from_them() {
    // y leads only into x and x is single copy by seed, but y is three times as deep: the links would say 1
    let mut g = Graph::new(K);
    let x = add(&mut g, 400, 10.0);
    let y = add(&mut g, 400, 30.0);
    g.add_link(f(y), f(x));
    let c = estimate(&g, &Options::default());
    assert_eq!(c.count[x], 1);
    assert_eq!(c.count[y], 3);
    assert_eq!(c.by_depth, 1);
}

#[test]
fn the_median_depth_is_that_of_half_the_bases_and_a_removed_segment_is_not_counted() {
    let mut g = Graph::new(K);
    add(&mut g, 100, 50.0);
    add(&mut g, 100, 10.0);
    let gone = add(&mut g, 1000, 3.0);
    g.remove(gone);
    // 100 bases at 10 and 100 at 50: half the bases are at 10 or less
    assert!((median_depth(&g) - 10.0).abs() < 1e-9);
    // an empty graph has a depth of 1, and no copy count is made
    let c = estimate(&Graph::new(K), &Options::default());
    assert_eq!((c.single_depth, c.count.len()), (1.0, 0));
}

#[test]
fn the_defaults_are_a_quarter_for_the_seeds_and_nothing_given_for_the_depth() {
    assert_eq!(
        Options::default(),
        Options {
            single_depth: None,
            tolerance: 0.25
        }
    );
}

#[test]
fn what_leads_into_a_segment_with_another_way_in_is_not_taken_for_all_of_it() {
    // x and y both lead only into r, which is 1.4 times as deep as they are: r is the sum of what leads in (2), not what x leads on to
    let mut g = Graph::new(K);
    let x = add(&mut g, 300, 10.0);
    let y = add(&mut g, 300, 10.0);
    let r = add(&mut g, 300, 14.0);
    g.add_link(f(x), f(r));
    g.add_link(f(y), f(r));
    let c = estimate(&g, &Options::default());
    assert_eq!((c.count[x], c.count[y], c.count[r]), (1, 1, 2));
}

#[test]
fn a_depth_that_rounds_up_is_rounded_up() {
    let mut g = Graph::new(K);
    let a = add(&mut g, 500, 10.0);
    let lone = add(&mut g, 200, 35.6);
    let c = estimate(&g, &Options::default());
    assert_eq!((c.count[a], c.count[lone]), (1, 4));
}
