use super::*;
use crate::seq::debruijn::graph::Graph;
use crate::seq::debruijn::sim::random;
use crate::seq::hybrid::sim::long_reads;

const K: usize = 55;

fn weight(len: usize, depth: f64) -> f64 {
    depth * (len + 1 - K) as f64
}

#[test]
fn the_defaults_are_those_of_the_documentation() {
    let o = Options::default();
    assert_eq!(
        (o.mode, o.min_quality, o.min_reads),
        (Mode::Normal, None, 1)
    );
    assert_eq!(o.copies, CopyOptions::default());
    assert_eq!(o.bridges, BridgeOptions::default());
}

#[test]
fn two_segments_with_a_gap_between_them_that_no_graph_has_are_joined_by_the_reads_bases() {
    let (x, z, gap) = (random(4000, 1), random(4000, 2), random(300, 3));
    let genome = [&x[..], &gap, &z].concat();
    let mut g = Graph::new(K);
    g.add(x.clone(), weight(4000, 20.0));
    g.add(z.clone(), weight(4000, 20.0));
    let reads = long_reads(&genome, false, 5000, 30.0, 0.02, 4);
    for mode in [Mode::Conservative, Mode::Normal, Mode::Bold] {
        let out = run(
            g.clone(),
            &reads,
            &Options {
                mode,
                ..Options::default()
            },
        )
        .unwrap();
        assert_eq!(out.report.bridges_made, 1, "{mode:?}");
        assert_eq!(out.graph.count(), 1);
        let s = &out.graph.segment(out.graph.ids()[0]).seq;
        // the two ends are the graph's own bases; the gap is a read's, with a read's errors
        assert_eq!(&s[..4000], &x[..]);
        assert_eq!(&s[s.len() - 4000..], &z[..]);
        assert!(
            (s.len() as i64 - genome.len() as i64).abs() <= 150,
            "{}",
            s.len()
        );
    }
}

#[test]
fn a_segment_whose_two_ends_the_reads_run_between_is_closed_into_a_circle() {
    let circle = random(5000, 5);
    // the assembly lacks 200 bases of the circle, where the ends should meet
    let mut g = Graph::new(K);
    g.add(circle[..4800].to_vec(), weight(4800, 20.0));
    let reads = long_reads(&circle, true, 3000, 30.0, 0.02, 6);
    let out = run(g, &reads, &Options::default()).unwrap();
    assert_eq!(out.report.bridges_made, 1);
    assert_eq!((out.graph.count(), out.report.circular), (1, 1));
    let id = out.graph.ids()[0];
    let once = out.graph.length(id) - (K - 1);
    assert!((once as i64 - 5000).abs() <= 150, "{once}");
    // the bases that the graph had are all there, in order, on one strand or the other
    let seq = &out.graph.segment(id).seq;
    let has = |needle: &[u8]| seq.windows(needle.len()).any(|w| w == needle);
    assert!(has(&circle[..4800]) || has(&crate::seq::dna::revcom(&circle[..4800])));
}

#[test]
fn a_gap_of_zero_or_less_with_no_link_and_no_route_is_not_a_bridge() {
    // two segments that follow each other in the genome with no gap and no link: the reads say 0, there is nothing to put in
    let (x, z) = (random(4000, 7), random(4000, 8));
    let genome = [&x[..], &z].concat();
    let mut g = Graph::new(K);
    g.add(x, weight(4000, 20.0));
    g.add(z, weight(4000, 20.0));
    let reads = long_reads(&genome, false, 5000, 30.0, 0.02, 9);
    let out = run(g, &reads, &Options::default()).unwrap();
    assert_eq!(out.report.joins_seen, 1);
    assert_eq!(out.report.bridges_found, 0);
    assert_eq!(out.graph.count(), 2);
}

#[test]
fn no_reads_leave_the_graph_as_it_is_and_the_report_counts_what_was_there() {
    let (x, y) = (random(3000, 10), random(1500, 11));
    let mut g = Graph::new(K);
    g.add(x, weight(3000, 20.0));
    g.add(y, weight(1500, 60.0));
    let out = run(g, &[], &Options::default()).unwrap();
    let r = &out.report;
    assert_eq!((r.segments, r.single_copy, r.repeats), (2, 1, 1));
    assert_eq!((r.joins_seen, r.bridges_found, r.bridges_made), (0, 0, 0));
    assert_eq!((r.segments_after, r.circular), (2, 0));
    assert!((r.single_depth - 20.0).abs() < 1e-9);
}

#[test]
fn the_median_of_the_gaps_is_taken() {
    assert_eq!(median(vec![5]), 5);
    assert_eq!(median(vec![9, 1, 5]), 5);
    assert_eq!(median(vec![1, 2, 3, 100]), 3);
}

/// The bases that a read from the end of `x` to the start of `z` has, with `gap` between them, exactly.
fn exact_read(x: &[u8], gap: &[u8], z: &[u8]) -> Vec<u8> {
    [&x[x.len() - 1500..], gap, &z[..1500]].concat()
}

#[test]
fn the_fewest_reads_is_a_number_of_reads_that_is_enough() {
    let (x, z, gap) = (random(3000, 21), random(3000, 22), random(300, 23));
    let mut g = Graph::new(K);
    g.add(x.clone(), weight(3000, 20.0));
    g.add(z.clone(), weight(3000, 20.0));
    let reads: Vec<(String, Vec<u8>)> = (0..3)
        .map(|i| (format!("r{i}"), exact_read(&x, &gap, &z)))
        .collect();
    let made = |min_reads| {
        run(
            g.clone(),
            &reads,
            &Options {
                min_reads,
                ..Options::default()
            },
        )
        .unwrap()
        .report
        .bridges_made
    };
    assert_eq!(made(3), 1);
    assert_eq!(made(4), 0);
}

#[test]
fn a_bridge_of_the_reads_own_bases_is_rated_by_how_much_the_reads_gaps_differ_and_by_being_the_reads(
) {
    let (x, z) = (random(3000, 24), random(3000, 25));
    let mut g = Graph::new(K);
    g.add(x.clone(), weight(3000, 20.0));
    g.add(z.clone(), weight(3000, 20.0));
    let f = |id: usize| crate::seq::debruijn::graph::handle(id, false);
    let gaps = [
        (100, random(100, 26)),
        (140, random(140, 27)),
        (120, random(120, 28)),
    ];
    let candidate = Candidate {
        from: f(0),
        to: f(1),
        reads: gaps
            .iter()
            .map(|(gap, bases)| crate::seq::hybrid::bridges::ReadGap {
                gap: *gap,
                bases: bases.clone(),
                covered: 0.5,
            })
            .collect(),
    };
    let b = resolve(&g, &[1, 1], &candidate).unwrap();
    // the read whose gap is the median one gives the bases
    assert_eq!(b.via, Via::Reads(gaps[2].1.clone()));
    // the gaps differ by 40: a deviation of 20, against a tolerance of 124, and it is not a route through the graph
    let expected = quality(&Facts {
        reads: 3,
        deviation: 20.0,
        tolerance: tolerance(120) as f64,
        lengths: (3000, 3000),
        covered: 0.5,
        depths: (20.0, 20.0),
        through_graph: false,
    });
    assert!(
        (b.quality - expected).abs() < 1e-9,
        "{} {expected}",
        b.quality
    );
    // no gap at all is not a bridge
    let none = Candidate {
        reads: vec![crate::seq::hybrid::bridges::ReadGap {
            gap: 0,
            bases: Vec::new(),
            covered: 0.5,
        }],
        ..candidate
    };
    assert!(resolve(&g, &[1, 1], &none).is_none());
}

#[test]
fn of_two_bridges_for_one_end_the_better_supported_is_made() {
    let (x, y, z) = (random(3000, 31), random(3000, 32), random(3000, 33));
    let (g1, g2) = (random(300, 34), random(300, 35));
    let mut g = Graph::new(K);
    for s in [&x, &y, &z] {
        g.add(s.to_vec(), weight(3000, 20.0));
    }
    // many reads from x on to y, and a few from x on to z
    let mut reads: Vec<(String, Vec<u8>)> = (0..12)
        .map(|i| (format!("y{i}"), exact_read(&x, &g1, &y)))
        .collect();
    reads.extend((0..2).map(|i| (format!("z{i}"), exact_read(&x, &g2, &z))));
    let out = run(g, &reads, &Options::default()).unwrap();
    assert_eq!(out.report.bridges_found, 2);
    assert_eq!(out.report.bridges_made, 1);
    // x is joined to y, whose last bases end the joined segment; z stands alone
    let joined = out
        .graph
        .ids()
        .into_iter()
        .find(|&i| out.graph.length(i) > 3500)
        .unwrap();
    let seq = out.graph.segment(joined).seq.clone();
    let has = |needle: &[u8]| seq.windows(needle.len()).any(|w| w == needle);
    assert!(
        has(&y[1000..2000]) && !has(&z[1000..2000]),
        "x must be joined to y"
    );
}
