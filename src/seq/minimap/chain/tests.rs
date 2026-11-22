use super::*;

fn params() -> Params {
    Params {
        k: 15,
        max_gap: 10_000,
        bandwidth: 500,
        lookback: 50,
        min_anchors: 3,
        min_score: 40.0,
    }
}

fn line(from: (usize, usize), step: usize, n: usize) -> Vec<Anchor> {
    (0..n)
        .map(|i| Anchor {
            x: from.0 + i * step,
            y: from.1 + i * step,
        })
        .collect()
}

#[test]
fn the_gap_cost_is_the_papers() {
    assert_eq!(gap_cost(0, 15), 0.0);
    assert!((gap_cost(10, 15) - (0.01 * 15.0 * 10.0 + 0.5 * 10f64.log2())).abs() < 1e-12);
    assert!((gap_cost(1, 15) - 0.15).abs() < 1e-12);
    // it grows with the difference, and with k
    assert!(gap_cost(100, 15) > gap_cost(10, 15) && gap_cost(10, 21) > gap_cost(10, 15));
}

#[test]
fn anchors_on_one_line_are_one_chain_that_scores_the_bases_they_add() {
    let chains = chain(&line((1000, 50), 10, 20), &params());
    assert_eq!(chains.len(), 1);
    // the first anchor gives its k-mer, each of the 19 after it the 10 bases that it adds
    assert_eq!(chains[0].score, 15.0 + 19.0 * 10.0);
    assert_eq!(chains[0].anchors.len(), 20);
    assert!(chains[0]
        .anchors
        .windows(2)
        .all(|w| w[0].x < w[1].x && w[0].y < w[1].y));
    // anchors a whole k-mer or more apart gain a whole k-mer each
    let sparse = chain(&line((0, 0), 40, 10), &params());
    assert_eq!(sparse[0].score, 15.0 * 10.0);
}

#[test]
fn a_difference_between_the_gaps_costs_what_the_paper_says_and_a_great_one_breaks_the_chain() {
    let mut a = line((0, 0), 10, 10);
    // the read has 7 extra bases after the tenth anchor: the next anchors lie 7 further on the read
    a.extend(line((100, 107), 10, 10));
    let c = chain(&a, &params());
    assert_eq!(c.len(), 1);
    assert!(
        (c[0].score - (15.0 + 19.0 * 10.0 - gap_cost(7, 15))).abs() < 1e-9,
        "{}",
        c[0].score
    );
    // 600 extra bases are past the bandwidth of 500: two chains
    let mut far = line((0, 0), 10, 10);
    far.extend(line((100, 700), 10, 10));
    let two = chain(&far, &params());
    assert_eq!(two.len(), 2);
    assert!(two.iter().all(|c| c.anchors.len() == 10));
    // and a wider bandwidth joins them
    let wide = chain(
        &far,
        &Params {
            bandwidth: 1000,
            ..params()
        },
    );
    assert_eq!(wide.len(), 1);
}

#[test]
fn two_places_for_the_same_read_are_two_chains_and_no_anchor_is_in_both() {
    let mut a = line((1000, 0), 12, 15);
    a.extend(line((9000, 0), 12, 12));
    let chains = chain(&a, &params());
    assert_eq!(chains.len(), 2);
    assert!(chains[0].score >= chains[1].score);
    let mut all: Vec<Anchor> = chains.iter().flat_map(|c| c.anchors.clone()).collect();
    let n = all.len();
    all.sort();
    all.dedup();
    assert_eq!(all.len(), n);
    assert_eq!(n, 27);
}

#[test]
fn an_anchor_that_is_off_the_line_is_left_out_and_a_cross_over_does_not_chain() {
    let mut a = line((0, 0), 10, 12);
    a.push(Anchor { x: 55, y: 300 });
    let c = chain(&a, &params());
    assert_eq!(c[0].anchors.len(), 12);
    // the read runs backwards along the target: no two can follow each other
    let backwards: Vec<Anchor> = (0..10)
        .map(|i| Anchor {
            x: i * 20,
            y: 200 - i * 20,
        })
        .collect();
    assert!(chain(&backwards, &params()).is_empty());
    // anchors further apart than the most are not joined
    let apart: Vec<Anchor> = (0..6)
        .map(|i| Anchor {
            x: i * 20_000,
            y: i * 20_000,
        })
        .collect();
    assert!(chain(&apart, &params()).is_empty());
}

#[test]
fn a_chain_with_too_few_anchors_or_too_low_a_score_is_left_out_and_the_input_may_be_in_any_order() {
    let two = line((0, 0), 10, 2);
    assert!(chain(&two, &params()).is_empty());
    assert_eq!(
        chain(
            &two,
            &Params {
                min_anchors: 2,
                min_score: 10.0,
                ..params()
            }
        )
        .len(),
        1
    );
    assert!(chain(
        &line((0, 0), 10, 5),
        &Params {
            min_score: 100.0,
            ..params()
        }
    )
    .is_empty());
    let mut shuffled = line((300, 20), 11, 25);
    shuffled.reverse();
    shuffled.push(shuffled[3]);
    let c = chain(&shuffled, &params());
    assert_eq!((c.len(), c[0].anchors.len()), (1, 25));
}

#[test]
fn looking_back_further_does_not_change_chains_that_are_regular() {
    let mut a = line((0, 0), 9, 80);
    a.extend(line((2000, 3000), 13, 60));
    let few = chain(
        &a,
        &Params {
            lookback: 5,
            ..params()
        },
    );
    let many = chain(
        &a,
        &Params {
            lookback: 500,
            ..params()
        },
    );
    assert_eq!(
        few.iter().map(|c| c.anchors.len()).collect::<Vec<_>>(),
        many.iter().map(|c| c.anchors.len()).collect::<Vec<_>>()
    );
    assert_eq!(few.len(), 2);
}

#[test]
fn the_limits_of_a_gap_and_of_the_bandwidth_are_themselves_allowed() {
    let few = Params {
        min_anchors: 2,
        min_score: 10.0,
        max_gap: 100,
        ..params()
    };
    // two anchors exactly the longest gap apart on both sequences are joined, one base further are not (each alone is no chain of two)
    assert_eq!(
        chain(&[Anchor { x: 0, y: 0 }, Anchor { x: 100, y: 100 }], &few).len(),
        1
    );
    assert_eq!(
        chain(&[Anchor { x: 0, y: 0 }, Anchor { x: 101, y: 101 }], &few).len(),
        0
    );
    // and the gaps may differ by exactly the bandwidth, but not by one more (a gap of 30 against 20 differs by 10: the two still gain
    // more than the difference costs)
    let band = Params {
        bandwidth: 10,
        min_anchors: 2,
        min_score: 10.0,
        ..params()
    };
    assert_eq!(
        chain(&[Anchor { x: 0, y: 0 }, Anchor { x: 30, y: 20 }], &band).len(),
        1
    );
    assert_eq!(
        chain(&[Anchor { x: 0, y: 0 }, Anchor { x: 31, y: 20 }], &band).len(),
        0
    );
}

#[test]
fn a_chain_that_runs_into_one_taken_before_it_counts_only_what_it_adds() {
    // twenty-five anchors on a line, and five more beside its end, a few bases further along the target and 400 further along the read: each of
    // the five joins the line (a hop of 400 costs 64), so a chain of them can only run back into the line, which has been taken. What they add
    // is their own bases less that cost, which is no chain; counted with the line's score they would be one.
    let mut a = line((0, 0), 10, 25);
    a.extend((0..5).map(|j| Anchor {
        x: 205 + 10 * j,
        y: 205 + 10 * j + 400,
    }));
    let lenient = Params {
        min_anchors: 1,
        min_score: 20.0,
        ..params()
    };
    let chains = chain(&a, &lenient);
    assert_eq!(
        chains.len(),
        1,
        "{:?}",
        chains
            .iter()
            .map(|c| (c.score, c.anchors.len()))
            .collect::<Vec<_>>()
    );
    assert_eq!(chains[0].anchors.len(), 25);
    // and a chain's score is never more than k for each anchor that it holds
    for c in chain(
        &a,
        &Params {
            min_anchors: 1,
            min_score: -1000.0,
            ..params()
        },
    ) {
        assert!(
            c.score <= 15.0 * c.anchors.len() as f64 + 1e-9,
            "{} for {}",
            c.score,
            c.anchors.len()
        );
    }
}
