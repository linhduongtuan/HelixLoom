use super::*;
use crate::seq::debruijn::sim::random;

fn index(segments: &[&[u8]]) -> Index {
    let targets = segments
        .iter()
        .enumerate()
        .map(|(i, s)| (i.to_string(), s.to_vec()))
        .collect();
    Index::build(targets, &MapOptions::default()).unwrap()
}

fn named(reads: Vec<Vec<u8>>) -> Vec<(String, Vec<u8>)> {
    reads
        .into_iter()
        .enumerate()
        .map(|(i, r)| (format!("r{i}"), r))
        .collect()
}

fn f(id: usize) -> Handle {
    handle(id, false)
}

/// Two segments of 3000 bases, and a read that goes from the last 1000 of the first, across 100 bases of its own, into the first 1000 of the second.
fn setup() -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    let (a, b, gap) = (random(3000, 1), random(3000, 2), random(100, 3));
    let read = [&a[2000..], &gap[..], &b[..1000]].concat();
    (a, b, gap, read)
}

#[test]
fn a_read_that_leaves_one_segment_and_enters_another_is_a_join_with_the_bases_between() {
    let (a, b, gap, read) = setup();
    let idx = index(&[&a, &b]);
    let c = candidates(
        &idx,
        &named(vec![read]),
        &[1, 1],
        &Options::default(),
        &MapOptions::default(),
    );
    assert_eq!(c.len(), 1);
    assert_eq!((c[0].from, c[0].to), (f(0), f(1)));
    let r = &c[0].reads[0];
    // the alignments stop a few bases short of the ends, which the gap counts back in
    assert!((r.gap - 100).abs() <= 5, "{}", r.gap);
    assert_eq!(r.bases, gap);
    assert!(r.covered > 0.3 && r.covered <= 1.0, "{}", r.covered);
}

#[test]
fn the_same_read_from_the_other_strand_is_the_same_join_with_the_same_bases() {
    let (a, b, gap, read) = setup();
    let idx = index(&[&a, &b]);
    let flipped = revcom(&read);
    let c = candidates(
        &idx,
        &named(vec![read, flipped]),
        &[1, 1],
        &Options::default(),
        &MapOptions::default(),
    );
    assert_eq!(c.len(), 1, "{c:?}");
    assert_eq!((c[0].from, c[0].to), (f(0), f(1)));
    assert_eq!(c[0].reads.len(), 2);
    assert!(c[0].reads.iter().all(|r| r.bases == gap));
}

#[test]
fn a_join_to_the_reverse_of_a_segment_is_kept_with_its_orientation() {
    let (a, b, gap, _) = setup();
    // the read goes on into the reverse complement of b
    let read = [&a[2000..], &gap[..], &revcom(&b[2000..])[..]].concat();
    let idx = index(&[&a, &b]);
    let c = candidates(
        &idx,
        &named(vec![read]),
        &[1, 1],
        &Options::default(),
        &MapOptions::default(),
    );
    assert_eq!(c.len(), 1);
    assert_eq!((c[0].from, c[0].to), (f(0), handle(1, true)));
}

#[test]
fn a_segment_that_is_not_single_copy_is_not_a_target() {
    let (a, b, _, read) = setup();
    let idx = index(&[&a, &b]);
    let c = candidates(
        &idx,
        &named(vec![read]),
        &[1, 2],
        &Options::default(),
        &MapOptions::default(),
    );
    assert!(c.is_empty());
}

#[test]
fn an_alignment_that_stops_far_from_the_end_does_not_join() {
    let (a, b, gap, _) = setup();
    // the read leaves a 700 bases before its end
    let read = [&a[1000..2300], &gap[..], &b[..1000]].concat();
    let idx = index(&[&a, &b]);
    let strict = Options {
        end_slack: 200,
        ..Options::default()
    };
    assert!(candidates(
        &idx,
        &named(vec![read.clone()]),
        &[1, 1],
        &strict,
        &MapOptions::default()
    )
    .is_empty());
    let lax = Options {
        end_slack: 800,
        ..Options::default()
    };
    assert_eq!(
        candidates(
            &idx,
            &named(vec![read]),
            &[1, 1],
            &lax,
            &MapOptions::default()
        )
        .len(),
        1
    );
}

#[test]
fn an_alignment_that_is_too_short_does_not_count() {
    let (a, b, gap, _) = setup();
    // 250 bases of b only
    let read = [&a[2000..], &gap[..], &b[..250]].concat();
    let idx = index(&[&a, &b]);
    let c = |min| {
        candidates(
            &idx,
            &named(vec![read.clone()]),
            &[1, 1],
            &Options {
                min_align: min,
                ..Options::default()
            },
            &MapOptions::default(),
        )
    };
    assert!(c(300).is_empty());
    assert_eq!(c(200).len(), 1);
}

#[test]
fn reads_that_make_the_same_join_are_one_candidate_in_the_order_of_the_handles() {
    let (a, b, gap, read) = setup();
    let c2 = random(3000, 9);
    let read2 = [&b[2000..], &random(80, 5)[..], &c2[..1000]].concat();
    let idx = index(&[&a, &b, &c2]);
    let reads = named(vec![read2.clone(), read.clone(), read]);
    let c = candidates(
        &idx,
        &reads,
        &[1, 1, 1],
        &Options::default(),
        &MapOptions::default(),
    );
    assert_eq!(c.len(), 2);
    assert_eq!((c[0].from, c[0].to, c[0].reads.len()), (f(0), f(1), 2));
    assert_eq!((c[1].from, c[1].to, c[1].reads.len()), (f(1), f(2), 1));
    assert!(c[0].reads.iter().all(|r| r.bases == gap));
}

#[test]
fn a_read_in_one_piece_that_lies_in_one_segment_or_is_not_there_makes_no_join() {
    let (a, b, _, _) = setup();
    let idx = index(&[&a, &b]);
    let inside = a[500..2500].to_vec();
    let alien = random(3000, 77);
    let c = candidates(
        &idx,
        &named(vec![inside, alien]),
        &[1, 1],
        &Options::default(),
        &MapOptions::default(),
    );
    assert!(c.is_empty());
}

#[test]
fn a_secondary_alignment_to_a_segment_that_shares_the_bases_is_not_a_join() {
    // b holds all of what the read lies on in a, and more: the read maps to both, the second as a secondary alignment
    let a = random(3000, 11);
    let b = [&a[1000..], &random(1000, 12)[..]].concat();
    let read = a[1000..].to_vec();
    let idx = index(&[&a, &b]);
    let c = candidates(
        &idx,
        &named(vec![read]),
        &[1, 1],
        &Options::default(),
        &MapOptions::default(),
    );
    assert!(c.is_empty(), "{c:?}");
}
