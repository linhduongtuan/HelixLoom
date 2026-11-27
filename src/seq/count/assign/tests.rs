use super::*;
use crate::seq::count::annotation;
use crate::seq::variants::sam::{parse_cigar, REVERSE, UNMAPPED};

/// g1 has the exons 101-200 and 301-400 on the plus strand, g2 has 351-450 on the minus, g3 has 1001-1100 with no strand, g4 is on chr u.
const SAF: &str = "GeneID\tChr\tStart\tEnd\tStrand\n\
g1\tt\t101\t200\t+\n\
g1\tt\t301\t400\t+\n\
g2\tt\t351\t450\t-\n\
g3\tt\t1001\t1100\t.\n\
g4\tu\t1\t50\t+\n";

fn ann(text: &str) -> annotation::Annotation {
    annotation::parse(text, &["exon".to_string()], "gene_id").unwrap()
}

fn read(pos: usize, cigar: &str) -> Record {
    let cigar = parse_cigar(cigar).unwrap();
    let length: usize = cigar
        .iter()
        .filter(|c| matches!(c.0, 'M' | 'I' | 'S' | '=' | 'X'))
        .map(|c| c.1)
        .sum();
    Record {
        name: "r".into(),
        flag: 0,
        target: Some("t".into()),
        pos,
        mapq: 60,
        cigar,
        seq: vec![b'A'; length],
        qual: Vec::new(),
        tags: Vec::new(),
    }
}

fn assigned(r: &Record, rules: &Rules) -> (Status, Vec<String>) {
    let a = ann(SAF);
    let index = Index::new(&a, rules.feature_level);
    let (status, targets) = index.assign(r, rules);
    let names = targets
        .iter()
        .map(|&t| {
            if rules.feature_level {
                format!(
                    "{}:{}",
                    a.names[a.features[t].meta],
                    a.features[t].start + 1
                )
            } else {
                a.names[t].clone()
            }
        })
        .collect();
    (status, names)
}

fn status(r: &Record) -> Status {
    assigned(r, &Rules::default()).0
}

#[test]
fn one_base_in_a_feature_is_an_overlap_and_none_is_not() {
    // g1's first exon holds the bases 100 to 199, counted from 0
    assert_eq!(status(&read(199, "10M")), Status::Assigned);
    assert_eq!(status(&read(200, "10M")), Status::NoFeatures);
    assert_eq!(status(&read(90, "10M")), Status::NoFeatures);
    assert_eq!(status(&read(91, "10M")), Status::Assigned);
    assert_eq!(assigned(&read(199, "10M"), &Rules::default()).1, ["g1"]);
}

#[test]
fn only_aligned_bases_overlap() {
    // skipping over an exon, deleting it, inserting at it, or clipping bases that would lie in it is not overlapping it
    assert_eq!(status(&read(50, "10M150N10M")), Status::NoFeatures);
    assert_eq!(status(&read(50, "10M150D10M")), Status::NoFeatures);
    assert_eq!(status(&read(80, "10M50I")), Status::NoFeatures);
    assert_eq!(status(&read(60, "40S")), Status::NoFeatures);
    // a clip before the aligned part does not move it
    assert_eq!(status(&read(50, "100S10M")), Status::NoFeatures);
    assert_eq!(status(&read(195, "50S10M")), Status::Assigned);
}

#[test]
fn a_read_over_two_exons_of_one_gene_counts_once_for_it() {
    let r = read(150, "50M100N50M");
    assert_eq!(
        assigned(&r, &Rules::default()),
        (Status::Assigned, vec!["g1".to_string()])
    );
    // at feature level it overlaps two exons: ambiguous, or both with -O
    let features = Rules {
        feature_level: true,
        ..Rules::default()
    };
    assert_eq!(assigned(&r, &features).0, Status::Ambiguity);
    let both = Rules {
        multi_overlap: true,
        ..features
    };
    assert_eq!(assigned(&r, &both).1, ["g1:101", "g1:301"]);
}

#[test]
fn a_base_of_a_gene_in_two_of_its_features_counts_once() {
    let a = ann("GeneID\tChr\tStart\tEnd\tStrand\nx\tt\t1\t100\t+\nx\tt\t51\t150\t+\n");
    let index = Index::new(&a, false);
    assert_eq!(index.lengths, [150]);
    let rules = |n: usize| Rules {
        min_overlap: n,
        ..Rules::default()
    };
    let r = read(0, "120M");
    assert_eq!(index.assign(&r, &rules(120)).0, Status::Assigned);
    assert_eq!(index.assign(&r, &rules(121)).0, Status::OverlappingLength);
}

#[test]
fn a_gene_on_two_sequences_has_the_length_of_both() {
    let a = ann("GeneID\tChr\tStart\tEnd\tStrand\nx\ta\t1\t10\t+\nx\tb\t1\t10\t+\n");
    assert_eq!(Index::new(&a, false).lengths, [20]);
    assert_eq!(Index::new(&a, true).lengths, [10, 10]);
}

#[test]
fn a_long_feature_that_starts_before_a_short_one_is_found_past_it() {
    let a = ann("GeneID\tChr\tStart\tEnd\tStrand\nlong\tt\t1\t1000\t+\nshort\tt\t11\t20\t+\n");
    let index = Index::new(&a, false);
    let (status, targets) = index.assign(&read(500, "10M"), &Rules::default());
    assert_eq!((status, targets), (Status::Assigned, vec![0]));
}

#[test]
fn a_read_on_two_genes_is_ambiguous_unless_each_is_allowed() {
    let r = read(360, "30M");
    assert_eq!(assigned(&r, &Rules::default()).0, Status::Ambiguity);
    let each = Rules {
        multi_overlap: true,
        ..Rules::default()
    };
    assert_eq!(assigned(&r, &each).1, ["g1", "g2"]);
}

#[test]
fn the_largest_overlap_takes_the_read_and_two_equal_ones_stay_ambiguous() {
    // 340..380: g1 (to 400) holds 40 bases, g2 (from 350) holds 30
    let largest = Rules {
        largest_overlap: true,
        ..Rules::default()
    };
    assert_eq!(assigned(&read(340, "40M"), &largest).1, ["g1"]);
    // 390..420: g1 holds 10, g2 holds 30
    assert_eq!(assigned(&read(390, "30M"), &largest).1, ["g2"]);
    // 360..390: both hold 30
    assert_eq!(assigned(&read(360, "30M"), &largest).0, Status::Ambiguity);
    let both = Rules {
        multi_overlap: true,
        ..largest
    };
    assert_eq!(assigned(&read(360, "30M"), &both).1, ["g1", "g2"]);
}

#[test]
fn the_strand_of_the_read_picks_the_gene_when_the_library_is_stranded() {
    let r = read(360, "30M");
    let mut reverse = r.clone();
    reverse.flag = REVERSE;
    let s = |n: u8| Rules {
        strandedness: n,
        ..Rules::default()
    };
    assert_eq!(assigned(&r, &s(1)).1, ["g1"]);
    assert_eq!(assigned(&reverse, &s(1)).1, ["g2"]);
    assert_eq!(assigned(&r, &s(2)).1, ["g2"]);
    assert_eq!(assigned(&reverse, &s(2)).1, ["g1"]);
    assert_eq!(assigned(&r, &s(0)).0, Status::Ambiguity);
    // a gene with no strand fits either
    let on_g3 = read(1010, "20M");
    let mut on_g3_reverse = on_g3.clone();
    on_g3_reverse.flag = REVERSE;
    for rules in [s(1), s(2)] {
        assert_eq!(assigned(&on_g3, &rules).1, ["g3"]);
        assert_eq!(assigned(&on_g3_reverse, &rules).1, ["g3"]);
    }
    // a read on the wrong strand of the only gene there overlaps nothing
    assert_eq!(assigned(&read(120, "10M"), &s(2)).0, Status::NoFeatures);
}

#[test]
fn the_filters_come_in_the_order_of_the_guide() {
    let mut unmapped = read(120, "10M");
    unmapped.flag = UNMAPPED;
    unmapped.mapq = 0;
    let q = Rules {
        min_mapq: 10,
        ..Rules::default()
    };
    assert_eq!(assigned(&unmapped, &q).0, Status::Unmapped);
    let mut low = read(120, "10M");
    low.mapq = 9;
    low.tags = vec!["NH:i:3".into()];
    low.flag = DUPLICATE;
    let all = Rules {
        min_mapq: 10,
        ignore_duplicates: true,
        ..Rules::default()
    };
    assert_eq!(assigned(&low, &all).0, Status::MappingQuality);
    low.mapq = 10;
    assert_eq!(assigned(&low, &all).0, Status::Duplicate);
    low.flag = SECONDARY;
    assert_eq!(assigned(&low, &all).0, Status::MultiMapping);
    let primary = Rules {
        primary_only: true,
        ..all
    };
    // with --primary the multi-mapping filter is not looked at, and the secondary is left out
    assert_eq!(assigned(&low, &primary).0, Status::Secondary);
    low.flag = 0;
    assert_eq!(assigned(&low, &primary).0, Status::Assigned);
}

#[test]
fn duplicates_are_left_out_only_when_asked() {
    let mut d = read(120, "10M");
    d.flag = DUPLICATE;
    assert_eq!(status(&d), Status::Assigned);
}

#[test]
fn multi_mapping_reads_are_left_out_unless_counted_and_one_line_with_nh_one_is_not_multi() {
    let mut r = read(120, "10M");
    r.tags = vec!["NH:i:1".into()];
    assert_eq!(status(&r), Status::Assigned);
    r.tags = vec!["NM:i:0".into(), "NH:i:2".into()];
    assert_eq!(status(&r), Status::MultiMapping);
    let m = Rules {
        multi_mapping: true,
        ..Rules::default()
    };
    assert_eq!(assigned(&r, &m).0, Status::Assigned);
    // without --primary a secondary alignment is counted like any other
    r.flag = SECONDARY;
    assert_eq!(assigned(&r, &m).0, Status::Assigned);
}

#[test]
fn alignments_with_a_junction_can_be_kept_or_left() {
    let spliced = read(150, "50M100N50M");
    let whole = read(120, "20M");
    let split = Rules {
        splice: Splice::Split,
        ..Rules::default()
    };
    let only_whole = Rules {
        splice: Splice::Whole,
        ..Rules::default()
    };
    assert_eq!(assigned(&spliced, &split).0, Status::Assigned);
    assert_eq!(assigned(&whole, &split).0, Status::NonSplit);
    assert_eq!(assigned(&spliced, &only_whole).0, Status::Split);
    assert_eq!(assigned(&whole, &only_whole).0, Status::Assigned);
}

#[test]
fn an_overlap_must_be_enough_in_every_way_that_is_asked() {
    // 190..210: 10 bases in g1, whose length is 200, of a read of 30 (10 of them clipped)
    let r = read(190, "10S20M");
    let with = |f: &dyn Fn(&mut Rules)| {
        let mut rules = Rules::default();
        f(&mut rules);
        assigned(&r, &rules).0
    };
    assert_eq!(with(&|x| x.min_overlap = 10), Status::Assigned);
    assert_eq!(with(&|x| x.min_overlap = 11), Status::OverlappingLength);
    // a share of the read counts the clipped bases in its length: 10 of 30
    assert_eq!(with(&|x| x.frac_overlap = 10.0 / 30.0), Status::Assigned);
    assert_eq!(with(&|x| x.frac_overlap = 0.34), Status::OverlappingLength);
    assert_eq!(with(&|x| x.frac_overlap_feature = 0.05), Status::Assigned);
    assert_eq!(
        with(&|x| x.frac_overlap_feature = 0.051),
        Status::OverlappingLength
    );
    assert_eq!(with(&|x| x.non_overlap = Some(20)), Status::Assigned);
    assert_eq!(
        with(&|x| x.non_overlap = Some(19)),
        Status::OverlappingLength
    );
    assert_eq!(
        with(&|x| x.non_overlap_feature = Some(190)),
        Status::Assigned
    );
    assert_eq!(
        with(&|x| x.non_overlap_feature = Some(189)),
        Status::OverlappingLength
    );
    // both must be met
    assert_eq!(
        with(&|x| {
            x.min_overlap = 5;
            x.frac_overlap = 0.9;
        }),
        Status::OverlappingLength
    );
}

#[test]
fn a_share_that_is_exactly_met_counts_whatever_the_rounding() {
    // 7 of 100 bases and a share of 0.07, whose product in binary is 7.000000000000001, a little over 7
    let a = ann("GeneID\tChr\tStart\tEnd\tStrand\nx\tt\t94\t100\t+\n");
    let index = Index::new(&a, false);
    let rules = Rules {
        frac_overlap: 0.07,
        ..Rules::default()
    };
    assert_eq!(index.assign(&read(0, "100M"), &rules).0, Status::Assigned);
}

#[test]
fn a_feature_that_ends_where_a_read_begins_does_not_overlap_it() {
    // "end" stops where the read starts; "other" is long, on the strand that is not counted, so the search looks back past it
    let a = ann("GeneID\tChr\tStart\tEnd\tStrand\nother\tt\t1\t1000\t-\nend\tt\t11\t200\t+\n");
    let index = Index::new(&a, false);
    let stranded = Rules {
        strandedness: 1,
        ..Rules::default()
    };
    assert_eq!(
        index.assign(&read(200, "10M"), &stranded).0,
        Status::NoFeatures
    );
    assert_eq!(
        index.assign(&read(199, "10M"), &stranded).0,
        Status::Assigned
    );
}

fn tally(records: &[Record], rules: &Rules) -> Tally {
    let a = ann(SAF);
    count(&Index::new(&a, rules.feature_level), records, rules).unwrap()
}

#[test]
fn each_assigned_alignment_adds_one_to_each_of_its_targets_and_the_statuses_are_counted() {
    let t = tally(
        &[
            read(120, "10M"),
            read(130, "10M"),
            read(1010, "10M"),
            read(600, "10M"),
            read(360, "30M"),
        ],
        &Rules::default(),
    );
    assert_eq!(t.counts, [2.0, 0.0, 1.0, 0.0]);
    assert_eq!(t.of(Status::Assigned), 3);
    assert_eq!(t.of(Status::NoFeatures), 1);
    assert_eq!(t.of(Status::Ambiguity), 1);
    assert_eq!(t.statuses.iter().sum::<u64>(), 5);
}

#[test]
fn a_fraction_shares_a_count_between_the_places_and_the_genes() {
    let mut multi = read(120, "10M");
    multi.tags = vec!["NH:i:4".into()];
    let both = read(360, "30M");
    let rules = Rules {
        multi_mapping: true,
        multi_overlap: true,
        fraction: true,
        ..Rules::default()
    };
    let t = tally(&[multi.clone(), both.clone()], &rules);
    assert_eq!(t.counts, [0.25 + 0.5, 0.5, 0.0, 0.0]);
    // a read that is both: 1 / (4 * 2)
    let mut mo = both.clone();
    mo.tags = vec!["NH:i:4".into()];
    let t = tally(&[mo], &rules);
    assert_eq!(t.counts, [0.125, 0.125, 0.0, 0.0]);
    // without the fraction each gets one
    let whole = Rules {
        fraction: false,
        ..rules.clone()
    };
    assert_eq!(
        tally(&[multi.clone(), both.clone()], &whole).counts,
        [2.0, 1.0, 0.0, 0.0]
    );
    // a fraction of the overlaps only leaves the multi-mapping count whole
    let overlaps = Rules {
        multi_mapping: false,
        ..rules.clone()
    };
    assert_eq!(
        tally(std::slice::from_ref(&both), &overlaps).counts,
        [0.5, 0.5, 0.0, 0.0]
    );
    // with --primary the multi-mapping share is not taken
    let primary = Rules {
        primary_only: true,
        multi_overlap: false,
        ..rules
    };
    assert_eq!(tally(&[multi], &primary).counts, [1.0, 0.0, 0.0, 0.0]);
}

#[test]
fn reads_of_pairs_are_refused_unless_each_is_asked_to_be_counted() {
    let mut r = read(120, "10M");
    r.flag = PAIRED;
    let a = ann(SAF);
    let index = Index::new(&a, false);
    let e = count(&index, &[r.clone()], &Rules::default()).unwrap_err();
    assert!(e.contains("pairs") && e.contains("--each-read"), "{e}");
    let each = Rules {
        each_read: true,
        ..Rules::default()
    };
    assert_eq!(count(&index, &[r], &each).unwrap().counts[0], 1.0);
}

#[test]
fn a_read_on_a_sequence_that_the_annotation_does_not_name_overlaps_nothing() {
    let mut r = read(120, "10M");
    r.target = Some("chrZ".into());
    assert_eq!(status(&r), Status::NoFeatures);
    r.target = Some("u".into());
    r.pos = 0;
    assert_eq!(assigned(&r, &Rules::default()).1, ["g4"]);
}

#[test]
fn the_statuses_are_in_order_with_the_names_that_featurecounts_writes() {
    let names: Vec<&str> = Status::ALL.iter().map(|s| s.name()).collect();
    assert_eq!(
        names,
        [
            "Assigned",
            "Unassigned_Unmapped",
            "Unassigned_MappingQuality",
            "Unassigned_Duplicate",
            "Unassigned_MultiMapping",
            "Unassigned_Secondary",
            "Unassigned_Split",
            "Unassigned_NonSplit",
            "Unassigned_NoFeatures",
            "Unassigned_Overlapping_Length",
            "Unassigned_Ambiguity"
        ]
    );
    // the number of each is its place
    assert!(Status::ALL
        .iter()
        .enumerate()
        .all(|(i, &s)| s as usize == i));
}

#[test]
fn the_defaults_are_those_of_the_guide() {
    let r = Rules::default();
    assert_eq!((r.min_overlap, r.min_mapq, r.strandedness), (1, 0, 0));
    assert!(
        !r.multi_overlap && !r.multi_mapping && !r.fraction && !r.feature_level && !r.primary_only
    );
    assert_eq!(r.splice, Splice::Any);
    assert_eq!(r.non_overlap, None);
    assert_eq!(r.non_overlap_feature, None);
}
