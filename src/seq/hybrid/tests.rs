//! The hybrid assembler on a made-up circle with a repeat that short reads cannot span, and long reads that can.

use super::apply::Mode;
use super::assemble::{run, Options};
use super::sim::{is_circle_of, long_reads, ring, short_graph};
use crate::seq::debruijn::graph::Graph;

fn options(mode: Mode) -> Options {
    Options {
        mode,
        ..Options::default()
    }
}

fn lengths(g: &Graph) -> Vec<usize> {
    let mut l: Vec<usize> = g.ids().into_iter().map(|i| g.length(i)).collect();
    l.sort_unstable();
    l
}

#[test]
fn the_graph_of_short_reads_has_the_repeat_once_and_two_ways_through_it() {
    let r = ring(3000, 1200, 1);
    let g = short_graph(&r.genome, true, 2);
    // the three flanks join into two segments (the one that closes the circle is with its neighbour), and the repeat
    assert_eq!(g.count(), 3);
    assert_eq!(lengths(&g), [1201, 3107, 6107]);
    assert_eq!(g.links().len(), 4);
}

#[test]
fn long_reads_that_span_the_repeat_make_the_circle_one_circular_segment_in_every_mode() {
    let r = ring(3000, 1200, 1);
    let g = short_graph(&r.genome, true, 2);
    let reads = long_reads(&r.genome, true, 5000, 40.0, 0.02, 3);
    for mode in [Mode::Conservative, Mode::Normal, Mode::Bold] {
        let out = run(g.clone(), &reads, &options(mode)).unwrap();
        assert_eq!(out.graph.count(), 1, "{mode:?}");
        let id = out.graph.ids()[0];
        assert!(out.graph.circular(id), "{mode:?}");
        // once round: without the k-1 bases that the end shares with the start
        let k1 = out.graph.k() - 1;
        let once = &out.graph.segment(id).seq[..out.graph.length(id) - k1];
        assert!(is_circle_of(once, &r.genome), "{mode:?}");
        let rep = &out.report;
        assert_eq!((rep.single_copy, rep.repeats, rep.bridges_made), (2, 1, 2));
        assert_eq!((rep.segments_after, rep.circular), (1, 1));
    }
}

#[test]
fn a_quality_that_no_bridge_reaches_leaves_the_graph_as_it_was() {
    let r = ring(3000, 1200, 1);
    let g = short_graph(&r.genome, true, 2);
    let reads = long_reads(&r.genome, true, 5000, 40.0, 0.02, 3);
    let out = run(
        g.clone(),
        &reads,
        &Options {
            min_quality: Some(100.0),
            ..Options::default()
        },
    )
    .unwrap();
    assert_eq!(out.report.bridges_made, 0);
    assert_eq!(out.report.bridges_found, 2);
    assert_eq!(lengths(&out.graph), lengths(&g));
}

#[test]
fn reads_that_are_not_from_the_genome_change_nothing() {
    let r = ring(3000, 1200, 1);
    let g = short_graph(&r.genome, true, 2);
    let other = ring(3000, 1200, 50);
    let reads = long_reads(&other.genome, true, 5000, 20.0, 0.02, 4);
    let out = run(g.clone(), &reads, &Options::default()).unwrap();
    assert_eq!(out.report.joins_seen, 0);
    assert_eq!(lengths(&out.graph), lengths(&g));
}

#[test]
fn one_read_is_enough_but_a_higher_minimum_asks_for_more() {
    let r = ring(3000, 1200, 1);
    let g = short_graph(&r.genome, true, 2);
    // few reads: about one of 5 kb in every 2 kb
    let reads = long_reads(&r.genome, true, 5000, 4.0, 0.02, 5);
    let low = run(g.clone(), &reads, &Options::default()).unwrap();
    let strict = run(
        g,
        &reads,
        &Options {
            min_reads: 1000,
            ..Options::default()
        },
    )
    .unwrap();
    assert!(low.report.bridges_made >= 1);
    assert_eq!(strict.report.bridges_made, 0);
}
