use super::*;
use crate::seq::minimap::sim::{genome, read, Rng};
use crate::seq::variants::sam::{REVERSE, SECONDARY, SUPPLEMENTARY, UNMAPPED};

fn options() -> Options {
    Options {
        cigar: true,
        ..Options::default()
    }
}

#[test]
fn every_read_gives_its_lines_in_the_order_it_came_and_a_read_that_does_not_map_gives_one_line_of_its_own(
) {
    let g = genome(8000, 31);
    let mut r = Rng::new(9);
    let junk = genome(1500, 777);
    let reads = vec![
        Read {
            name: "a".into(),
            seq: read(&g, 500, 2000, false, 0.03, &mut r),
            qual: Vec::new(),
        },
        Read {
            name: "junk".into(),
            seq: junk,
            qual: Vec::new(),
        },
        Read {
            name: "c".into(),
            seq: read(&g, 4000, 2000, true, 0.03, &mut r),
            qual: Vec::new(),
        },
    ];
    let sam = map_to_sam(&[("ref".to_string(), g)], &reads, &options()).unwrap();
    assert_eq!(sam.targets, [("ref".to_string(), 8000)]);
    let names: Vec<&str> = sam.records.iter().map(|r| r.name.as_str()).collect();
    let pos = |n: &str| names.iter().position(|x| *x == n).unwrap();
    assert!(pos("a") < pos("junk") && pos("junk") < pos("c"));
    let junk = &sam.records[pos("junk")];
    assert!(!junk.mapped());
    assert_eq!(junk.flag, UNMAPPED);
    // the others lie on the reference, the second on the other strand
    assert!(sam.records[pos("a")].usable());
    let c = &sam.records[pos("c")];
    assert!(c.usable() && c.flag & REVERSE != 0);
}

#[test]
fn sam_needs_the_alignment_of_the_bases_and_does_not_do_overlaps() {
    let g = genome(2000, 32);
    let reads = vec![Read {
        name: "a".into(),
        seq: g[100..900].to_vec(),
        qual: Vec::new(),
    }];
    let t = [("ref".to_string(), g)];
    let e = map_to_sam(&t, &reads, &Options::default()).unwrap_err();
    assert!(e.contains("base-level alignment"), "{e}");
    let e = map_to_sam(
        &t,
        &reads,
        &Options {
            overlaps: true,
            ..options()
        },
    )
    .unwrap_err();
    assert!(e.contains("Overlaps"), "{e}");
    assert!(map_to_sam(&[], &reads, &options()).is_err());
}

#[test]
fn a_read_that_lies_in_two_places_has_one_primary_line_and_a_supplementary_one() {
    let g = genome(20_000, 33);
    // two thousand bases from one place and then two thousand from another, as a chimera, or a read across the start of a circle, would be
    let seq = [&g[2000..4000], &g[14_000..16_000]].concat();
    let reads = vec![Read {
        name: "split".into(),
        seq,
        qual: Vec::new(),
    }];
    let sam = map_to_sam(&[("ref".to_string(), g)], &reads, &options()).unwrap();
    assert_eq!(sam.records.len(), 2, "{:?}", sam.records);
    let flags: Vec<u16> = sam.records.iter().map(|r| r.flag).collect();
    assert_eq!(flags, [0, SUPPLEMENTARY]);
    // both have the bases and are counted, and each place has one line
    assert!(sam.records.iter().all(|r| r.usable() && !r.seq.is_empty()));
    assert!(sam
        .records
        .iter()
        .all(|r| r.tags.iter().any(|t| t == "NH:i:1")));
}

#[test]
fn nh_counts_the_lines_of_a_place() {
    // two copies of one stretch: the read's place has a secondary
    let unit = genome(3000, 34);
    let mut g = genome(2000, 35);
    g.extend_from_slice(&unit);
    g.extend(genome(2000, 36));
    g.extend_from_slice(&unit);
    g.extend(genome(2000, 37));
    let reads = vec![Read {
        name: "twice".into(),
        seq: unit[200..2800].to_vec(),
        qual: Vec::new(),
    }];
    let sam = map_to_sam(&[("ref".to_string(), g)], &reads, &options()).unwrap();
    assert_eq!(sam.records.len(), 2, "{:?}", sam.records);
    assert_eq!(sam.records[0].flag, 0);
    assert_eq!(sam.records[1].flag, SECONDARY);
    assert!(sam
        .records
        .iter()
        .all(|r| r.tags.iter().any(|t| t == "NH:i:2")));
}
