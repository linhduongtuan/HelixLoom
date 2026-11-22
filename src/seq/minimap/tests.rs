use super::*;
use crate::seq::dna::revcom;
use crate::seq::minimap::sim::{genome, read, Rng};

#[test]
fn the_presets_set_what_the_manual_says_and_the_defaults_are_those_of_map_ont() {
    let d = Options::default();
    assert_eq!(
        (
            d.k,
            d.w,
            d.max_gap,
            d.bandwidth,
            d.min_anchors,
            d.max_secondary,
            d.lookback
        ),
        (15, 10, 10_000, 500, 3, 5, 50)
    );
    assert_eq!(
        (d.occ_fraction, d.min_score, d.secondary_ratio, d.mask_level),
        (0.0002, 40.0, 0.8, 0.5)
    );
    assert!(!d.cigar && !d.overlaps);
    assert_eq!(Options::preset("map-ont"), Some(d.clone()));
    let ava = Options::preset("ava-ont").unwrap();
    assert_eq!(
        (ava.k, ava.w, ava.min_score, ava.bandwidth, ava.overlaps),
        (15, 5, 100.0, 2000, true)
    );
    assert_eq!(ava.max_gap, d.max_gap);
    assert_eq!(Options::preset("splice"), None);
}

fn columns(paf: &str) -> Vec<Vec<&str>> {
    paf.lines().map(|l| l.split('\t').collect()).collect()
}

#[test]
fn reads_are_mapped_in_the_order_they_came_and_each_line_is_paf() {
    let g = genome(30_000, 41);
    let targets = vec![("ref".to_string(), g.clone())];
    let mut r = Rng::new(6);
    let reads: Vec<(String, Vec<u8>)> = (0..12)
        .map(|n| {
            (
                format!("read{n}"),
                read(&g, 1000 + 2000 * n, 2500, n % 2 == 0, 0.03, &mut r),
            )
        })
        .collect();
    let paf = map_text(&targets, &reads, &Options::default()).unwrap();
    let rows = columns(&paf);
    assert!(rows.len() >= 12);
    assert!(rows
        .iter()
        .all(|f| f.len() >= 12 && f[5] == "ref" && f[6] == "30000"));
    // the first line of each read is its primary, in the order of the reads
    let first: Vec<&str> = rows
        .iter()
        .filter(|f| f.contains(&"tp:A:P"))
        .map(|f| f[0])
        .collect();
    assert_eq!(
        first,
        reads.iter().map(|r| r.0.as_str()).collect::<Vec<_>>()
    );
    for (n, f) in rows.iter().filter(|f| f.contains(&"tp:A:P")).enumerate() {
        let start: usize = f[7].parse().unwrap();
        assert!(start.abs_diff(1000 + 2000 * n) < 300, "read{n} at {start}");
        assert_eq!(f[4], if n % 2 == 0 { "-" } else { "+" });
    }
    assert!(map_text(&[], &reads, &Options::default())
        .unwrap_err()
        .contains("no target"));
    assert_eq!(map_text(&targets, &[], &Options::default()).unwrap(), "");
}

#[test]
fn overlaps_are_found_once_for_each_pair_with_the_earlier_read_as_the_query_and_never_with_itself()
{
    let g = genome(20_000, 51);
    let mut r = Rng::new(13);
    // 0 and 1 share 4,000 to 6,000; 1 and 2 share 7,000 to 9,000 (2 is on the other strand); 3 is alone
    let spans = [(1000, 5000), (4000, 5000), (7000, 5000), (15_000, 4000)];
    let reads: Vec<(String, Vec<u8>)> = spans
        .iter()
        .enumerate()
        .map(|(n, &(a, len))| (format!("r{n}"), read(&g, a, len, n == 2, 0.03, &mut r)))
        .collect();
    let paf = overlap_text(&reads, &Options::preset("ava-ont").unwrap()).unwrap();
    let rows = columns(&paf);
    assert!(rows.iter().all(|f| f[0] != f[5]), "a read with itself");
    let pair = |a: &str, b: &str| {
        rows.iter()
            .filter(|f| f[0] == a && f[5] == b)
            .collect::<Vec<_>>()
    };
    // each pair is there once, the earlier read being the query
    assert_eq!(pair("r0", "r1").len(), 1, "{paf}");
    assert_eq!(pair("r1", "r2").len(), 1, "{paf}");
    assert!(pair("r1", "r0").is_empty() && pair("r2", "r1").is_empty());
    assert!(
        rows.iter().all(|f| f[0] != "r3" && f[5] != "r3"),
        "r3 overlaps nobody"
    );
    assert!(pair("r0", "r2").is_empty());
    // where: r0 to r1 is the last 2,000 of r0 and the first 2,000 of r1, both on one strand
    let f = &pair("r0", "r1")[0];
    let num = |i: usize| f[i].parse::<usize>().unwrap();
    assert_eq!(f[4], "+");
    assert!(
        num(2).abs_diff(3000) < 250 && num(3).abs_diff(5000) < 250,
        "{f:?}"
    );
    assert!(num(7) < 200 && num(8).abs_diff(2000) < 250, "{f:?}");
    // r1 to r2: the other strand
    assert_eq!(pair("r1", "r2")[0][4], "-");
    // no quality of mapping: there is none for an overlap
    assert!(rows.iter().all(|f| f[11] == "255"));
    // one read has no overlap with anyone
    assert_eq!(overlap_text(&reads[..1], &Options::default()).unwrap(), "");
}

#[test]
fn a_reverse_complement_read_overlaps_as_the_forward_one_does_with_the_strand_the_other_way() {
    let g = genome(15_000, 61);
    let a = g[0..5000].to_vec();
    let b = g[3000..8000].to_vec();
    let o = Options::preset("ava-ont").unwrap();
    let forward = overlap_text(&[("a".into(), a.clone()), ("b".into(), b.clone())], &o).unwrap();
    let flipped = overlap_text(&[("a".into(), a), ("b".into(), revcom(&b))], &o).unwrap();
    let (f, r) = (columns(&forward), columns(&flipped));
    assert_eq!((f.len(), r.len()), (1, 1));
    assert_eq!((f[0][4], r[0][4]), ("+", "-"));
    assert_eq!(f[0][2], r[0][2]);
    assert_eq!(f[0][3], r[0][3]);
}

/// How well and how fast it maps 1,500 reads of nanopore size to a 1 Mb genome, with the numbers printed. Left out of the usual run (it takes
/// seconds): `cargo test --release minimap_accuracy -- --ignored --nocapture`.
#[test]
#[ignore = "takes a few seconds and prints numbers"]
fn minimap_accuracy_and_speed_on_a_million_bases() {
    let g = genome(1_000_000, 81);
    let o = Options::default();
    let mut r = Rng::new(5);
    let reads: Vec<(usize, usize, bool, Vec<u8>)> = (0..1500)
        .map(|n| {
            let len = 1000 + r.below(9000);
            let start = r.below(g.len() - len);
            (
                start,
                len,
                n % 2 == 0,
                read(&g, start, len, n % 2 == 0, 0.03, &mut r),
            )
        })
        .collect();
    let t = std::time::Instant::now();
    let index = Index::build(vec![("g".to_string(), g.clone())], &o).unwrap();
    let built = t.elapsed();
    let t = std::time::Instant::now();
    let hits: Vec<Vec<map::Hit>> =
        crate::seq::par::par_map(&reads, |(_, _, _, q)| map_read(&index, q, &o, &|_| false));
    let mapped = t.elapsed();
    let (mut right, mut none, mut wrong_strong) = (0, 0, 0);
    for ((start, len, rev, _), h) in reads.iter().zip(&hits) {
        match h.first() {
            None => none += 1,
            Some(b) => {
                let overlap = b.tend.min(start + len).saturating_sub(b.tstart.max(*start));
                if b.rev == *rev && overlap as f64 >= 0.8 * *len as f64 {
                    right += 1;
                } else if b.mapq >= 30 {
                    wrong_strong += 1;
                }
            }
        }
    }
    let bases: usize = reads.iter().map(|r| r.3.len()).sum();
    println!(
        "{} reads, {} Mb of reads: index {:.2} s ({} minimizers), mapping {:.2} s ({:.1} Mb/s); right {right}, unmapped {none}, wrong with quality 30 or more {wrong_strong}",
        reads.len(),
        bases / 1_000_000,
        built.as_secs_f64(),
        index.minimizers(),
        mapped.as_secs_f64(),
        bases as f64 / 1e6 / mapped.as_secs_f64(),
    );
    assert!(right >= reads.len() * 95 / 100);
    assert!(wrong_strong <= 2);
}
