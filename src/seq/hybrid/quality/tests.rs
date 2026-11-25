use super::*;

fn perfect() -> Facts {
    Facts {
        reads: 9,
        deviation: 0.0,
        tolerance: 100.0,
        lengths: (5000, 5000),
        covered: 0.5,
        depths: (10.0, 10.0),
        through_graph: true,
    }
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

#[test]
fn a_bridge_with_nothing_against_it_scores_the_share_that_its_reads_give() {
    // 9 reads: 9/10
    assert!(near(quality(&perfect()), 90.0));
    let one = Facts {
        reads: 1,
        ..perfect()
    };
    assert!(near(quality(&one), 50.0));
    let none = Facts {
        reads: 0,
        ..perfect()
    };
    assert!(near(quality(&none), 0.0));
}

#[test]
fn more_reads_score_higher_and_never_over_a_hundred() {
    let q = |n| {
        quality(&Facts {
            reads: n,
            ..perfect()
        })
    };
    assert!(q(1) < q(2) && q(2) < q(10) && q(10) < q(1000));
    assert!(q(1_000_000) < 100.0);
}

#[test]
fn a_gap_that_the_route_misses_scores_less_down_to_half_at_the_tolerance_and_no_less() {
    let q = |d| {
        quality(&Facts {
            deviation: d,
            ..perfect()
        })
    };
    assert!(near(q(0.0), 90.0));
    assert!(near(q(50.0), 90.0 * 0.75));
    assert!(near(q(100.0), 45.0));
    assert!(near(q(500.0), 45.0));
}

#[test]
fn short_segments_score_in_proportion_to_their_length_each_up_to_a_thousand_bases() {
    let q = |a, b| {
        quality(&Facts {
            lengths: (a, b),
            ..perfect()
        })
    };
    assert!(near(q(1000, 1000), 90.0));
    assert!(near(q(500, 5000), 45.0));
    assert!(near(q(500, 500), 22.5));
    assert!(near(q(10_000, 1000), 90.0));
}

#[test]
fn alignments_that_the_minimizers_cover_little_of_score_less_up_to_four_tenths() {
    let q = |c| {
        quality(&Facts {
            covered: c,
            ..perfect()
        })
    };
    assert!(near(q(0.4), 90.0));
    assert!(near(q(0.9), 90.0));
    assert!(near(q(0.2), 45.0));
    assert!(near(q(0.0), 0.0));
}

#[test]
fn depths_that_differ_score_as_the_lower_over_the_higher_whichever_is_which() {
    let q = |a, b| {
        quality(&Facts {
            depths: (a, b),
            ..perfect()
        })
    };
    assert!(near(q(10.0, 5.0), 45.0));
    assert!(near(q(5.0, 10.0), 45.0));
    assert!(near(q(0.0, 0.0), 90.0));
}

#[test]
fn a_bridge_of_the_reads_own_bases_scores_half_of_one_through_the_graph() {
    let q = |g| {
        quality(&Facts {
            through_graph: g,
            ..perfect()
        })
    };
    assert!(near(q(false), q(true) * 0.5));
}
