use super::*;
use crate::seq::debruijn::graph::handle;

const K: usize = 5;

fn f(id: usize) -> Handle {
    handle(id, false)
}

fn add(g: &mut Graph, len: usize) -> usize {
    g.add(vec![b'A'; len], 10.0)
}

/// `a`, `b`, and between them any number of repeats; the first two are single copy (count 1), the rest are not.
fn graph(repeats: &[usize]) -> (Graph, Vec<u32>) {
    let mut g = Graph::new(K);
    add(&mut g, 100);
    add(&mut g, 100);
    let mut count = vec![1, 1];
    for &len in repeats {
        add(&mut g, len);
        count.push(2);
    }
    (g, count)
}

#[test]
fn the_tolerance_is_a_hundred_bases_and_a_fifth_of_the_gap() {
    assert_eq!(tolerance(0), 100);
    assert_eq!(tolerance(1000), 300);
    assert_eq!(tolerance(-1000), 300);
    assert_eq!(tolerance(-4), 100);
}

#[test]
fn a_direct_link_is_a_route_with_nothing_between_and_a_distance_of_minus_the_overlap() {
    let (mut g, count) = graph(&[]);
    g.add_link(f(0), f(1));
    let r = routes(&g, &count, f(0), f(1), -4);
    assert_eq!(
        r,
        [Route {
            via: vec![],
            distance: -4
        }]
    );
    // the reads put them far apart: no
    assert!(routes(&g, &count, f(0), f(1), 500).is_empty());
    // not linked: no
    let (g2, count2) = graph(&[]);
    assert!(routes(&g2, &count2, f(0), f(1), -4).is_empty());
}

#[test]
fn a_repeat_between_is_passed_and_adds_its_length_less_two_overlaps() {
    let (mut g, count) = graph(&[300]);
    g.add_link(f(0), f(2));
    g.add_link(f(2), f(1));
    // 300 - 2 x 4
    let r = routes(&g, &count, f(0), f(1), 292);
    assert_eq!(
        r,
        [Route {
            via: vec![f(2)],
            distance: 292
        }]
    );
    // the reads' gap may be off by 100 and a fifth of itself: 490 is just within (198 of 198), 491 is not; 160 is (132 of 132), 159 is not
    assert_eq!(routes(&g, &count, f(0), f(1), 490).len(), 1);
    assert!(routes(&g, &count, f(0), f(1), 491).is_empty());
    assert_eq!(routes(&g, &count, f(0), f(1), 160).len(), 1);
    assert!(routes(&g, &count, f(0), f(1), 159).is_empty());
}

#[test]
fn only_segments_that_are_not_single_copy_are_passed() {
    let (mut g, mut count) = graph(&[300]);
    g.add_link(f(0), f(2));
    g.add_link(f(2), f(1));
    count[2] = 1;
    assert!(routes(&g, &count, f(0), f(1), 292).is_empty());
    // a segment that has not been given a number of copies is not one either
    count[2] = 0;
    assert_eq!(routes(&g, &count, f(0), f(1), 292).len(), 1);
}

#[test]
fn of_two_routes_the_one_nearer_the_gap_comes_first_and_a_shorter_one_wins_a_tie() {
    let (mut g, count) = graph(&[300, 340]);
    for r in [2, 3] {
        g.add_link(f(0), f(r));
        g.add_link(f(r), f(1));
    }
    // 292 and 332
    let near_short = routes(&g, &count, f(0), f(1), 300);
    assert_eq!(
        near_short.iter().map(|r| r.distance).collect::<Vec<_>>(),
        [292, 332]
    );
    let near_long = routes(&g, &count, f(0), f(1), 330);
    assert_eq!(
        near_long.iter().map(|r| r.distance).collect::<Vec<_>>(),
        [332, 292]
    );
}

#[test]
fn a_route_through_several_repeats_adds_them_all() {
    let (mut g, count) = graph(&[100, 200]);
    g.add_link(f(0), f(2));
    g.add_link(f(2), f(3));
    g.add_link(f(3), f(1));
    // 100 + 200 - 3 x 4
    let r = routes(&g, &count, f(0), f(1), 288);
    assert_eq!(
        r,
        [Route {
            via: vec![f(2), f(3)],
            distance: 288
        }]
    );
}

#[test]
fn a_route_may_pass_at_most_eight_repeats() {
    let lens = [100usize; 9];
    let (mut g, count) = graph(&lens);
    let ids: Vec<usize> = (2..11).collect();
    g.add_link(f(0), f(ids[0]));
    for w in ids.windows(2) {
        g.add_link(f(w[0]), f(w[1]));
    }
    g.add_link(f(ids[8]), f(1));
    // nine: too many; eight: fine
    let nine = 9 * 100 - 10 * 4;
    assert!(routes(&g, &count, f(0), f(1), nine as i64).is_empty());
    let (mut g8, count8) = graph(&lens[..8]);
    let ids8: Vec<usize> = (2..10).collect();
    g8.add_link(f(0), f(ids8[0]));
    for w in ids8.windows(2) {
        g8.add_link(f(w[0]), f(w[1]));
    }
    g8.add_link(f(ids8[7]), f(1));
    let eight = 8 * 100 - 9 * 4;
    assert_eq!(routes(&g8, &count8, f(0), f(1), eight as i64).len(), 1);
}

#[test]
fn a_repeat_that_leads_to_itself_does_not_make_the_search_go_on_for_ever() {
    let (mut g, count) = graph(&[100]);
    g.add_link(f(0), f(2));
    g.add_link(f(2), f(2));
    g.add_link(f(2), f(1));
    // through it once, twice... each a longer route: the ones within the tolerance are found and the search ends
    let r = routes(&g, &count, f(0), f(1), 92);
    assert!(!r.is_empty());
    assert_eq!(r[0].via, [f(2)]);
    assert!(r.len() < 10);
}
