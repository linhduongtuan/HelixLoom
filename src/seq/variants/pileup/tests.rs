use super::*;
use crate::seq::variants::sam::{DUPLICATE, SECONDARY, UNMAPPED};
use crate::seq::variants::sim::record;

fn one(length: usize, records: &[Record], f: &Filter) -> Pileup {
    pile(&[("t".to_string(), length)], records, f).remove(0)
}

#[test]
fn the_weight_of_a_base_is_the_log_of_the_odds_that_it_is_right_over_each_of_the_others() {
    // e = 0.1: ln(3 * 0.9 / 0.1) = ln 27
    assert!((f64::from(weight(10)) - 27f64.ln()).abs() < 1e-4);
    assert!((f64::from(weight(20)) - 297f64.ln()).abs() < 1e-4);
    assert!((f64::from(weight(30)) - 2997f64.ln()).abs() < 1e-3);
    assert!(weight(13) < weight(20) && weight(20) < weight(30));
    // a base that is more likely wrong than right gives nothing, and not less than nothing
    assert_eq!((weight(0), weight(1)), (0.0, 0.0));
    assert!(weight(5) > 0.0);
}

#[test]
fn bases_are_counted_by_position_with_the_evidence_of_their_qualities() {
    let reads = [
        record("a", 0, "4M", "ACGT", 30),
        record("b", 0, "4M", "ACGT", 20),
        record("c", 1, "3M", "TGT", 30),
    ];
    let p = one(6, &reads, &Filter::default());
    assert_eq!(p.reads, 3);
    assert_eq!(p.counts[0], [2, 0, 0, 0]);
    assert_eq!(p.counts[1], [0, 2, 0, 1]);
    assert_eq!(p.counts[3], [0, 0, 0, 3]);
    assert_eq!(p.counts[4], [0, 0, 0, 0]);
    assert_eq!(p.cover, [2, 3, 3, 3, 0, 0]);
    let expect = weight(30) + weight(20);
    assert!((p.scores[0][0] - expect).abs() < 1e-4);
    assert_eq!(p.depth(2), 3);
    assert_eq!(p.bases(1), 3);
}

#[test]
fn a_base_below_the_lowest_quality_covers_the_position_but_gives_no_base() {
    let reads = [
        record("low", 0, "2M", "AC", 12),
        record("edge", 0, "2M", "AC", 13),
    ];
    let p = one(2, &reads, &Filter::default());
    assert_eq!(p.cover, [2, 2]);
    assert_eq!(p.counts[0], [1, 0, 0, 0]);
    // the lowest quality is a choice that can be moved
    let f = Filter {
        min_base_quality: 20,
        ..Filter::default()
    };
    assert_eq!(one(2, &reads, &f).counts[0], [0; 4]);
}

#[test]
fn a_quality_above_the_highest_is_taken_as_the_highest() {
    let reads = [record("a", 0, "1M", "A", 60)];
    let p = one(1, &reads, &Filter::default());
    assert!((p.scores[0][0] - weight(30)).abs() < 1e-6);
    let f = Filter {
        max_base_quality: 40,
        ..Filter::default()
    };
    assert!((one(1, &reads, &f).scores[0][0] - weight(40)).abs() < 1e-6);
}

#[test]
fn a_read_with_no_qualities_has_the_default_quality_on_every_base() {
    let mut r = record("a", 0, "2M", "AC", 0);
    r.qual.clear();
    let p = one(2, &[r.clone()], &Filter::default());
    assert!((p.scores[0][0] - weight(15)).abs() < 1e-6);
    let f = Filter {
        default_quality: 10,
        ..Filter::default()
    };
    // 10 is under the lowest quality of 13: no base
    assert_eq!(one(2, &[r], &f).counts[0], [0; 4]);
}

#[test]
fn a_read_counts_only_when_it_is_usable_and_its_mapping_quality_is_enough() {
    let good = record("g", 0, "2M", "AC", 30);
    let mut low = record("l", 0, "2M", "AC", 30);
    low.mapq = 5;
    let mut second = record("s", 0, "2M", "AC", 30);
    second.flag = SECONDARY;
    let mut dup = record("d", 0, "2M", "AC", 30);
    dup.flag = DUPLICATE;
    let mut lost = record("u", 0, "2M", "AC", 30);
    lost.flag = UNMAPPED;
    let all = [good, low, second, dup, lost];
    assert_eq!(one(2, &all, &Filter::default()).reads, 2);
    let f = Filter {
        min_mapq: 6,
        ..Filter::default()
    };
    let p = one(2, &all, &f);
    assert_eq!((p.reads, p.cover[0]), (1, 1));
    // exactly the lowest counts
    let f = Filter {
        min_mapq: 5,
        ..Filter::default()
    };
    assert_eq!(one(2, &all, &f).reads, 2);
}

#[test]
fn at_most_the_most_reads_count_at_a_position() {
    let reads: Vec<Record> = (0..5)
        .map(|i| record(&format!("r{i}"), 0, "2M", "AC", 30))
        .collect();
    let f = Filter {
        max_depth: 3,
        ..Filter::default()
    };
    let p = one(2, &reads, &f);
    assert_eq!((p.cover[0], p.counts[0][0]), (3, 3));
    // a read that arrives at a full position is not counted at it, and its insertion is not counted either
    let mut with_insertion: Vec<Record> = (0..3)
        .map(|i| record(&format!("r{i}"), 0, "2M", "AC", 30))
        .collect();
    with_insertion.push(record("late", 0, "1M2I1M", "AGGC", 30));
    let p = one(2, &with_insertion, &f);
    assert!(p.inserted.is_empty());
}

#[test]
fn a_deletion_is_counted_on_each_position_it_covers_and_the_bases_after_it_are_the_reads() {
    let p = one(
        10,
        &[record("a", 2, "2M3D3M", "ACGTA", 30)],
        &Filter::default(),
    );
    assert_eq!(p.deleted, [0, 0, 0, 0, 1, 1, 1, 0, 0, 0]);
    assert_eq!(p.cover, [0, 0, 1, 1, 0, 0, 0, 1, 1, 1]);
    assert_eq!(p.depth(5), 1);
    assert_eq!(p.bases(5), 0);
    assert_eq!(p.counts[7], [0, 0, 1, 0]);
}

#[test]
fn an_insertion_is_kept_after_the_position_before_it_and_counted_by_what_it_inserts() {
    let reads = [
        record("a", 10, "3M2I3M", "ACGTTACG", 30),
        record("b", 10, "3M2I3M", "ACGTTACG", 30),
        record("c", 10, "3M2I3M", "ACGAAACG", 30),
    ];
    let p = one(30, &reads, &Filter::default());
    let at = &p.inserted[&12];
    assert_eq!(at[&b"TT".to_vec()], 2);
    assert_eq!(at[&b"AA".to_vec()], 1);
    assert_eq!(p.inserted.len(), 1);
    // the bases of an insertion are not bases of the target
    assert_eq!(p.cover[13], 3);
}

#[test]
fn an_insertion_with_no_base_before_it_is_not_counted() {
    // the read opens with it, or it follows a deletion, or a soft clip that opens the read
    let reads = [
        record("a", 5, "2I3M", "TTACG", 30),
        record("b", 5, "2M1D2I2M", "ACTTAC", 30),
        record("c", 0, "1M2I1M", "AGGC", 30),
    ];
    let p = one(20, &reads, &Filter::default());
    assert_eq!(p.inserted.len(), 1);
    assert!(p.inserted.contains_key(&0));
    // soft clips use read bases and no target
    let clipped = one(
        10,
        &[record("s", 3, "2S3M", "GGACG", 30)],
        &Filter::default(),
    );
    assert_eq!(clipped.cover[3], 1);
    assert_eq!(clipped.counts[3], [1, 0, 0, 0]);
    assert_eq!(clipped.cover[2], 0);
}

#[test]
fn a_skip_moves_along_the_target_and_an_n_in_the_read_covers_but_gives_no_base() {
    let p = one(
        12,
        &[record("a", 0, "2M4N2M", "ACNT", 30)],
        &Filter::default(),
    );
    assert_eq!(p.cover, [1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]);
    assert_eq!(p.bases(6), 0);
    assert_eq!(p.bases(7), 1);
}

#[test]
fn what_runs_past_the_end_of_the_target_is_dropped() {
    let p = one(4, &[record("a", 2, "6M", "ACGTAC", 30)], &Filter::default());
    assert_eq!(p.cover, [0, 0, 1, 1]);
    let d = one(4, &[record("a", 1, "1M5D", "A", 30)], &Filter::default());
    assert_eq!(d.deleted, [0, 0, 1, 1]);
}

#[test]
fn records_go_to_the_target_they_name_and_one_that_is_not_in_the_list_is_passed_over() {
    let mut on_u = record("u", 0, "2M", "GG", 30);
    on_u.target = Some("u".to_string());
    let mut elsewhere = record("z", 0, "2M", "TT", 30);
    elsewhere.target = Some("zzz".to_string());
    let piles = pile(
        &[("t".to_string(), 2), ("u".to_string(), 2)],
        &[record("t", 0, "2M", "AC", 30), on_u, elsewhere],
        &Filter::default(),
    );
    assert_eq!(piles[0].counts[0], [1, 0, 0, 0]);
    assert_eq!(piles[1].counts[0], [0, 0, 1, 0]);
    assert_eq!((piles[0].reads, piles[1].reads), (1, 1));
}

#[test]
fn the_filter_has_the_defaults_of_the_manual_and_the_bases_are_four() {
    let f = Filter::default();
    assert_eq!(
        (
            f.min_mapq,
            f.min_base_quality,
            f.max_depth,
            f.max_base_quality,
            f.default_quality
        ),
        (0, 13, 250, 30, 15)
    );
    assert_eq!(
        (*b"AcGtN-").map(base_index),
        [Some(0), Some(1), Some(2), Some(3), None, None]
    );
}

#[test]
fn the_deletions_too_count_only_up_to_the_most_reads() {
    let reads: Vec<Record> = (0..5)
        .map(|i| record(&format!("r{i}"), 0, "2M2D2M", "ACGT", 30))
        .collect();
    let f = Filter {
        max_depth: 3,
        ..Filter::default()
    };
    let p = one(6, &reads, &f);
    assert_eq!((p.deleted[2], p.deleted[3]), (3, 3));
    assert_eq!(p.cover[0], 3);
}

#[test]
fn an_insertion_right_after_a_skip_follows_no_base_and_is_not_counted() {
    let p = one(
        10,
        &[record("a", 0, "2M2N1I2M", "ACGTT", 30)],
        &Filter::default(),
    );
    assert!(p.inserted.is_empty(), "{:?}", p.inserted);
    assert_eq!(p.cover, [1, 1, 0, 0, 1, 1, 0, 0, 0, 0]);
}

#[test]
fn u_in_a_sam_line_is_a_t() {
    assert_eq!(base_index(b'U'), Some(3));
    assert_eq!(base_index(b'u'), Some(3));
}
