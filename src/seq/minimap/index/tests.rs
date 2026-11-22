use super::*;
use crate::seq::minimap::sketch::sketch;

fn random(n: usize, seed: u32) -> Vec<u8> {
    let mut x = seed.wrapping_mul(2_654_435_761).wrapping_add(11);
    (0..n)
        .map(|_| {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((x >> 16) & 3) as usize]
        })
        .collect()
}

#[test]
fn every_minimizer_of_every_target_is_found_by_its_hash_with_its_place() {
    let o = Options::default();
    let targets = vec![
        ("one".to_string(), random(3000, 1)),
        ("two".to_string(), random(2000, 2)),
    ];
    let index = Index::build(targets.clone(), &o).unwrap();
    assert_eq!(
        (
            index.targets(),
            index.name(1),
            index.length(0),
            index.k,
            index.w
        ),
        (2, "two", 3000, 15, 10)
    );
    let mut total = 0;
    for (t, (_, seq)) in targets.iter().enumerate() {
        for m in sketch(seq, o.k, o.w) {
            total += 1;
            let found = index.hits(m.hash);
            assert!(
                found.contains(&Seed {
                    hash: m.hash,
                    target: t as u32,
                    pos: m.pos as u32,
                    rev: m.rev
                }),
                "{m:?}"
            );
        }
    }
    assert_eq!(index.minimizers(), total);
    // a hash that nobody has
    assert!(index.hits(12345).is_empty());
}

#[test]
fn what_the_two_targets_share_is_found_under_both() {
    let shared = random(400, 9);
    let a = [random(300, 3), shared.clone()].concat();
    let b = [shared, random(300, 4)].concat();
    let index = Index::build(vec![("a".into(), a), ("b".into(), b)], &Options::default()).unwrap();
    let both = sketch(&random(400, 9), 15, 10)
        .into_iter()
        .filter(|m| {
            index.hits(m.hash).iter().any(|s| s.target == 0)
                && index.hits(m.hash).iter().any(|s| s.target == 1)
        })
        .count();
    assert!(both > 50, "{both}");
}

#[test]
fn a_minimizer_with_too_many_holders_is_ignored_and_the_limit_follows_the_share_asked() {
    // five thousand bases of their own, and a unit of 50 bases that is there a hundred times
    let unit = random(50, 7);
    let target = [random(5000, 5), unit.repeat(100)].concat();
    let mut o = Options::default();
    let plain = Index::build(vec![("t".into(), target.clone())], &o).unwrap();
    // 0.02 % of about 5,000 is no minimizer: the hundred-fold ones are indexed as any
    let repeat_hash = sketch(&unit.repeat(3), 15, 10)[0].hash;
    assert!(
        plain.hits(repeat_hash).len() >= 90,
        "{}",
        plain.hits(repeat_hash).len()
    );
    o.occ_fraction = 0.05;
    let filtered = Index::build(vec![("t".into(), target)], &o).unwrap();
    assert_eq!(filtered.limit(), 10);
    assert!(filtered.hits(repeat_hash).is_empty());
    // a hash of one holder is still there
    let own = sketch(&random(5000, 5), 15, 10)[3].hash;
    assert!(!filtered.hits(own).is_empty());
}

#[test]
fn targets_are_kept_in_upper_case_and_a_choice_that_is_none_is_refused() {
    let index = Index::build(
        vec![("t".into(), b"acgtacgtnnacgtacgtacgtacgt".to_vec())],
        &Options::default(),
    )
    .unwrap();
    assert_eq!(index.sequence(0), b"ACGTACGTNNACGTACGTACGTACGT");
    let one = vec![("t".to_string(), random(100, 1))];
    for (k, w, why) in [
        (0, 10, "k must be"),
        (33, 10, "k must be"),
        (15, 0, "w must"),
    ] {
        let o = Options {
            k,
            w,
            ..Options::default()
        };
        assert!(
            Index::build(one.clone(), &o).err().unwrap().contains(why),
            "{k} {w}"
        );
    }
    assert!(Index::build(Vec::new(), &Options::default())
        .err()
        .unwrap()
        .contains("no target"));
}

#[test]
fn the_limit_is_the_count_at_the_foot_of_the_top_share_between_the_least_and_the_most() {
    // a hundred minimizers: ninety of one holder, nine of twenty and one of five hundred
    let counts: Vec<usize> = std::iter::repeat_n(1, 90)
        .chain(std::iter::repeat_n(20, 9))
        .chain([500])
        .collect();
    // the top 1.5 % of 100 is the place 98.5: the count there (counting from 0 and rounding down) is the 99th, which is 20
    assert_eq!(occurrence_limit(counts.clone(), 0.015), 20);
    assert_eq!(occurrence_limit(counts.clone(), 0.01), 500);
    assert_eq!(occurrence_limit(counts.clone(), 0.0), 500);
    // under the least that the manual allows, the least; over the most, the most
    assert_eq!(occurrence_limit(counts.clone(), 0.5), 10);
    assert_eq!(occurrence_limit(vec![5_000_000], 0.0), 1_000_000);
    // a share over 1 or under 0 is held to them, and nothing indexed is the least
    assert_eq!(occurrence_limit(counts, 7.0), 10);
    assert_eq!(occurrence_limit(Vec::new(), 0.0002), 10);
}
