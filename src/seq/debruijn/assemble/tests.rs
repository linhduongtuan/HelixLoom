use super::*;
use crate::seq::debruijn::sim::{random, reads, tiled};

fn options(ks: &[usize]) -> Options {
    Options {
        ks: ks.to_vec(),
        ..Options::default()
    }
}

#[test]
fn the_k_that_are_used_when_none_is_asked_for_are_those_that_fit_in_the_longest_read() {
    assert_eq!(default_ks(150), [21, 33, 55]);
    assert_eq!(default_ks(55), [21, 33, 55]);
    assert_eq!(default_ks(54), [21, 33]);
    assert_eq!(default_ks(33), [21, 33]);
    assert_eq!(default_ks(32), [21]);
    assert_eq!(default_ks(21), [21]);
    assert!(default_ks(20).is_empty());
}

#[test]
fn the_defaults_are_those_of_the_documentation() {
    let o = Options::default();
    assert_eq!(
        (
            o.min_count,
            o.tip_ratio,
            o.bulge_diff,
            &o.cutoff,
            &o.tip_length,
            &o.bulge_length
        ),
        (1, 2.5, 0.1, &Cutoff::Auto, &None, &None)
    );
    assert_eq!(o.isolated_length, 200);
    assert!(o.ks.is_empty());
}

#[test]
fn k_that_are_not_odd_in_range_and_rising_are_refused() {
    let pieces = vec![random(100, 1)];
    for ks in [[22].as_slice(), &[2], &[65], &[33, 21], &[21, 21]] {
        assert!(assemble(&pieces, &options(ks)).is_err(), "{ks:?}");
    }
    assert!(assemble(&pieces, &options(&[21, 33])).is_ok());
    assert!(assemble(&pieces, &options(&[21, 21]))
        .unwrap_err()
        .contains("rising order"));
}

#[test]
fn reads_too_short_for_the_smallest_default_k_say_so() {
    let e = assemble(&[random(20, 1)], &options(&[])).unwrap_err();
    assert!(
        e.contains("too short") && e.contains("21") && e.contains("20 bases"),
        "{e}"
    );
    assert!(assemble(&[], &options(&[])).is_err());
    // with a k that fits it goes
    assert!(assemble(&[random(20, 1)], &options(&[15])).is_ok());
}

#[test]
fn nothing_is_made_when_no_k_mer_is_long_enough_for_a_chosen_k() {
    let e = assemble(&[random(10, 1)], &options(&[21])).unwrap_err();
    assert!(e.contains("no k-mer of 21"), "{e}");
}

#[test]
fn one_step_is_reported_for_each_k_with_what_was_done() {
    let genome = random(2000, 2);
    let pieces = reads(&genome, false, 100, 30.0, 0.0, 3);
    let a = assemble(&pieces, &options(&[21, 33])).unwrap();
    assert_eq!(a.steps.iter().map(|s| s.k).collect::<Vec<_>>(), [21, 33]);
    assert_eq!(a.pieces, pieces.len());
    assert_eq!(a.bases, pieces.len() * 100);
    for s in &a.steps {
        assert_eq!((s.segments, s.bases), (1, 2000));
        assert_eq!(s.kmers, 2000 - s.k + 1);
        assert_eq!((s.tips, s.bulges, s.low_coverage), (0, 0, 0));
    }
    assert_eq!(a.graph.k(), 33);
}

#[test]
fn a_k_mer_seen_fewer_times_than_the_minimum_is_dropped() {
    let genome = random(1500, 4);
    let pieces = reads(&genome, false, 100, 30.0, 0.01, 5);
    let all = assemble(&pieces, &options(&[21])).unwrap();
    let min = assemble(
        &pieces,
        &Options {
            min_count: 3,
            ..options(&[21])
        },
    )
    .unwrap();
    assert!(min.steps[0].kmers < all.steps[0].kmers);
    // dropping the k-mers seen once or twice leaves the graph without the errors to clean
    assert!(min.steps[0].segments_made < all.steps[0].segments_made);
}

#[test]
fn the_cutoff_that_is_asked_for_is_the_one_used() {
    let genome = random(1500, 6);
    let pieces = reads(&genome, false, 100, 30.0, 0.0, 7);
    let used = |cutoff| {
        assemble(
            &pieces,
            &Options {
                cutoff,
                ..options(&[21])
            },
        )
        .unwrap()
        .steps[0]
            .cutoff
    };
    assert_eq!(used(Cutoff::Off), 0.0);
    assert_eq!(used(Cutoff::Value(7.5)), 7.5);
    let auto = used(Cutoff::Auto);
    assert!(auto > 1.5 && auto < 6.0, "{auto}");
}

#[test]
fn the_k_mers_of_the_graph_before_are_put_in_with_their_coverage_where_the_reads_have_none() {
    let pack = Pack::new(5).unwrap();
    let mut before = Graph::new(3);
    // 12 bases: 10 k-mers of 3, seen 4 times each on average
    let seq = random(12, 31);
    before.add(seq.clone(), 40.0);
    let known = pack.canonical(pack.pack(&seq[..5]).unwrap());
    let mut counts: Table<f32> = Table::default();
    counts.insert(known, 9.0);
    carry(&mut counts, &before, pack);
    // the distinct k-mers of 5 in the 12 bases
    let mut distinct = Vec::new();
    pack.for_each(&seq, |f, r| distinct.push(f.min(r)));
    distinct.sort_unstable();
    distinct.dedup();
    assert_eq!(counts.len(), distinct.len());
    assert_eq!(counts[&known], 9.0);
    assert!(counts
        .iter()
        .filter(|(k, _)| **k != known)
        .all(|(_, n)| *n == 4.0));
}

#[test]
fn a_segment_shorter_than_the_new_k_adds_no_k_mer() {
    let pack = Pack::new(9).unwrap();
    let mut before = Graph::new(3);
    before.add(b"ACGTTG".to_vec(), 16.0);
    let mut counts: Table<f32> = Table::default();
    carry(&mut counts, &before, pack);
    assert!(counts.is_empty());
}

#[test]
fn a_short_stretch_that_no_read_joins_to_the_rest_is_removed_at_the_last_k_unless_asked_to_stay() {
    let genome = random(3000, 61);
    let island = random(150, 62);
    let mut pieces = reads(&genome, false, 100, 30.0, 0.0, 63);
    pieces.extend(reads(&island, false, 100, 30.0, 0.0, 64));
    let a = assemble(&pieces, &options(&[21, 33])).unwrap();
    assert_eq!(a.steps[0].isolated, 0);
    assert_eq!((a.steps[1].isolated, a.steps[1].segments), (1, 1));
    let kept = assemble(
        &pieces,
        &Options {
            isolated_length: 0,
            ..options(&[21, 33])
        },
    )
    .unwrap();
    assert_eq!((kept.steps[1].isolated, kept.steps[1].segments), (0, 2));
}

#[test]
fn a_k_mer_seen_as_often_as_the_minimum_is_kept_and_one_seen_fewer_times_is_not() {
    // reads of 100 every 20 bases: each k-mer of 21 that is not near an end is in exactly 4 reads
    let genome = random(1000, 71);
    let pieces = tiled(&genome, 100, 20);
    let kept = assemble(
        &pieces,
        &Options {
            min_count: 4,
            ..options(&[21])
        },
    )
    .unwrap();
    assert!(kept.steps[0].kmers > 800);
    let none = assemble(
        &pieces,
        &Options {
            min_count: 5,
            ..options(&[21])
        },
    );
    assert!(none.unwrap_err().contains("no k-mer of 21"));
}
