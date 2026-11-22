use super::*;
use crate::seq::dna::revcom;
use crate::seq::minimap::sim::{genome, read, Rng};

fn index(g: &[u8], o: &Options) -> Index {
    Index::build(vec![("chr".to_string(), g.to_vec())], o).unwrap()
}

fn none(_: usize) -> bool {
    false
}

/// A mapping in a line, for the message of a test that fails.
fn brief(h: &Hit) -> String {
    format!(
        "target {} {} query {}..{} target {}..{} score {:.0} anchors {} {} mapq {} second {:.0}",
        h.target,
        if h.rev { "-" } else { "+" },
        h.qstart,
        h.qend,
        h.tstart,
        h.tend,
        h.score,
        h.anchors.len(),
        if h.primary { "primary" } else { "secondary" },
        h.mapq,
        h.best_secondary
    )
}

#[test]
fn the_quality_of_a_mapping_is_the_papers_formula_held_to_zero_and_sixty() {
    // 40 (1 - f2/f1) min(1, m/10) ln f1
    assert_eq!(mapq(100.0, 0.0, 20), 60);
    assert_eq!(mapq(100.0, 100.0, 20), 0);
    assert_eq!(mapq(100.0, 120.0, 20), 0);
    assert_eq!(
        mapq(100.0, 50.0, 5),
        (40.0f64 * 0.5 * 0.5 * 100f64.ln()).round() as u8
    );
    assert_eq!(
        mapq(30.0, 25.0, 10),
        (40.0f64 * (1.0 / 6.0) * 30f64.ln()).round() as u8
    );
    // below the ceiling it is what the formula gives, and at it the ceiling
    assert_eq!(mapq(30.0, 10.0, 10), 60);
    assert_eq!(mapq(1.0, 0.0, 50), 0);
    assert_eq!(mapq(0.0, 0.0, 50), 0);
    // fewer than ten anchors lower it in proportion
    assert!(mapq(40.0, 0.0, 3) < mapq(40.0, 0.0, 10));
}

#[test]
fn reads_from_a_known_place_are_mapped_there_on_the_strand_they_came_from() {
    let g = genome(30_000, 1);
    let o = Options::default();
    let idx = index(&g, &o);
    let mut r = Rng::new(4);
    let (mut right, mut strong, total) = (0, 0, 80);
    for n in 0..total {
        let len = 1500 + r.below(3500);
        let start = r.below(g.len() - len);
        let rev = n % 2 == 1;
        let q = read(&g, start, len, rev, 0.03, &mut r);
        let hits = map_read(&idx, &q, &o, &none);
        let Some(best) = hits.first().filter(|h| h.primary) else {
            continue;
        };
        let overlap = best
            .tend
            .min(start + len)
            .saturating_sub(best.tstart.max(start));
        if best.rev == rev && overlap as f64 >= 0.8 * len as f64 {
            right += 1;
        }
        if best.mapq >= 30 {
            strong += 1;
        }
    }
    assert!(right >= total * 95 / 100, "{right} of {total}");
    assert!(strong >= total * 90 / 100, "{strong} of {total}");
}

#[test]
fn an_exact_read_is_placed_to_within_a_window_on_either_strand_with_the_read_counted_as_written() {
    let g = genome(20_000, 2);
    let o = Options::default();
    let idx = index(&g, &o);
    let slack = o.k + o.w;
    let forward = &g[4000..6000];
    let h = &map_read(&idx, forward, &o, &none)[0];
    assert!(!h.rev && h.primary);
    assert!(
        h.tstart >= 4000 && h.tstart <= 4000 + slack && h.tend <= 6000 && h.tend >= 6000 - slack,
        "{}",
        brief(h)
    );
    assert_eq!((h.qstart, h.qend), (h.tstart - 4000, h.tend - 4000));
    // the reverse complement: the same place, the other strand, and the read's own coordinates run the other way
    let back = revcom(forward);
    let h = &map_read(&idx, &back, &o, &none)[0];
    assert!(h.rev);
    assert!(
        h.tstart >= 4000 && h.tstart <= 4000 + slack && h.tend <= 6000 && h.tend >= 6000 - slack,
        "{}",
        brief(h)
    );
    assert_eq!((h.qstart, h.qend), (6000 - h.tend, 6000 - h.tstart));
    assert!(h.mapq == 60 && h.score > 1500.0);
}

#[test]
fn a_read_from_a_repeat_has_a_second_place_and_a_low_quality_while_one_across_its_edge_is_certain()
{
    // 8,000 bases, a repeat of 3,000, 8,000 of their own, the repeat again with two bases changed, and 8,000 more
    let (a, repeat, b, c) = (
        genome(8000, 11),
        genome(3000, 12),
        genome(8000, 13),
        genome(8000, 14),
    );
    let mut again = repeat.clone();
    for i in (0..again.len()).step_by(1500) {
        again[i] = if again[i] == b'A' { b'C' } else { b'A' };
    }
    let g = [a.clone(), repeat.clone(), b.clone(), again, c].concat();
    let o = Options::default();
    let idx = index(&g, &o);
    let mut r = Rng::new(3);
    let inside = read(&g, 8500, 2000, false, 0.02, &mut r);
    let hits = map_read(&idx, &inside, &o, &none);
    let primaries: Vec<&Hit> = hits.iter().filter(|h| h.primary).collect();
    assert_eq!(
        primaries.len(),
        1,
        "{:?}",
        hits.iter().map(brief).collect::<Vec<_>>()
    );
    assert!(primaries[0].best_secondary > 0.0, "the other copy scores");
    assert!(primaries[0].mapq <= 10, "{}", primaries[0].mapq);
    assert!(hits.iter().any(|h| !h.primary));
    // a read that runs out of the repeat into the 8,000 after the first copy is the first copy's, and certain
    let across = read(&g, 10_000, 3000, false, 0.02, &mut r);
    let best = &map_read(&idx, &across, &o, &none)[0];
    assert!(best.tstart < 10_500 && best.mapq >= 30, "{}", brief(best));
}

#[test]
fn a_read_made_of_two_places_has_two_primaries_that_do_not_overlap_on_the_read() {
    let g = genome(40_000, 3);
    let o = Options::default();
    let idx = index(&g, &o);
    let mut r = Rng::new(8);
    let mut q = read(&g, 3000, 2500, false, 0.03, &mut r);
    let first = q.len();
    q.extend(read(&g, 25_000, 2500, true, 0.03, &mut r));
    let hits = map_read(&idx, &q, &o, &none);
    let primaries: Vec<&Hit> = hits.iter().filter(|h| h.primary).collect();
    assert_eq!(
        primaries.len(),
        2,
        "{:?}",
        hits.iter().map(brief).collect::<Vec<_>>()
    );
    let (x, y) = (primaries[0], primaries[1]);
    assert!(
        x.qend <= y.qstart + 100 || y.qend <= x.qstart + 100,
        "{} / {}",
        brief(x),
        brief(y)
    );
    let (one, two) = if x.tstart < y.tstart { (x, y) } else { (y, x) };
    assert!(!one.rev && two.rev);
    assert!(one.tstart.abs_diff(3000) < 300 && two.tstart.abs_diff(25_000) < 300);
    assert!((one.qend as f64 - first as f64).abs() < 300.0);
}

#[test]
fn what_is_not_in_the_target_a_short_read_and_a_read_of_n_map_nowhere() {
    let g = genome(20_000, 5);
    let o = Options::default();
    let idx = index(&g, &o);
    assert!(map_read(&idx, &genome(3000, 99), &o, &none).is_empty());
    assert!(map_read(&idx, &g[100..110], &o, &none).is_empty());
    assert!(map_read(&idx, b"", &o, &none).is_empty());
    assert!(map_read(&idx, &vec![b'N'; 3000], &o, &none).is_empty());
}

#[test]
fn a_target_that_is_skipped_is_not_mapped_to_and_the_other_is() {
    let a = genome(10_000, 21);
    let b = genome(10_000, 22);
    let o = Options::default();
    let idx = Index::build(vec![("a".into(), a.clone()), ("b".into(), b)], &o).unwrap();
    let q = &a[2000..4000];
    assert_eq!(map_read(&idx, q, &o, &none)[0].target, 0);
    assert!(map_read(&idx, q, &o, &|t| t == 0).is_empty());
}

fn hit(target: usize, q: (usize, usize), score: f64) -> Hit {
    Hit {
        target,
        rev: false,
        qstart: q.0,
        qend: q.1,
        tstart: 0,
        tend: q.1 - q.0,
        score,
        anchors: vec![Anchor { x: 0, y: 0 }; 20],
        primary: true,
        best_secondary: 0.0,
        mapq: 255,
        divergence: None,
        covered: 0,
        base: None,
    }
}

#[test]
fn a_mapping_that_overlaps_a_better_one_by_half_of_the_shorter_is_its_secondary() {
    let o = Options::default();
    let hits = vec![
        hit(0, (0, 1000), 800.0),
        // the same stretch of the read elsewhere, 90 % of the best: secondary, kept
        hit(1, (0, 1000), 720.0),
        // overlaps the shorter of the two (800) by 300, which is under half: another primary
        hit(2, (700, 1500), 300.0),
        // overlaps 600 of its 1,000 against the first: secondary, but under the ratio
        hit(3, (400, 1400), 200.0),
        // far from all: a primary of its own
        hit(4, (2000, 2600), 150.0),
    ];
    let out = label(hits, &o);
    let kinds: Vec<(usize, bool)> = out.iter().map(|h| (h.target, h.primary)).collect();
    // the order: the primary, its secondaries, the next primary ...
    assert_eq!(kinds, [(0, true), (1, false), (2, true), (4, true)]);
    let best = &out[0];
    assert_eq!(best.best_secondary, 720.0);
    assert_eq!(best.mapq, mapq(800.0, 720.0, 20));
    // a secondary has none of its own
    assert_eq!((out[1].mapq, out[1].best_secondary), (0, 0.0));
    // an overlap of exactly the share (500 of the shorter's 1,000) is enough
    let half = label(
        vec![hit(0, (0, 1000), 800.0), hit(1, (500, 1500), 790.0)],
        &o,
    );
    assert_eq!(
        half.iter().map(|h| h.primary).collect::<Vec<_>>(),
        [true, false]
    );
    let just_under = label(
        vec![hit(0, (0, 1000), 800.0), hit(1, (501, 1501), 790.0)],
        &o,
    );
    assert_eq!(
        just_under.iter().map(|h| h.primary).collect::<Vec<_>>(),
        [true, true]
    );
    // few are kept: at most N
    let many: Vec<Hit> = std::iter::once(hit(0, (0, 1000), 800.0))
        .chain((1..9).map(|t| hit(t, (0, 1000), 790.0 - t as f64)))
        .collect();
    let o2 = Options {
        max_secondary: 3,
        ..Options::default()
    };
    assert_eq!(label(many, &o2).len(), 4);
    // and a lower ratio keeps what the higher left out
    let low = Options {
        secondary_ratio: 0.2,
        mask_level: 0.5,
        ..Options::default()
    };
    let kept = label(
        vec![hit(0, (0, 1000), 800.0), hit(3, (400, 1400), 200.0)],
        &low,
    );
    assert_eq!(kept.iter().filter(|h| !h.primary).count(), 1);
}

#[test]
fn with_the_alignment_asked_for_the_cigar_walks_both_sequences_and_agrees_with_the_counts() {
    let g = genome(20_000, 7);
    let o = Options {
        cigar: true,
        ..Options::default()
    };
    let idx = index(&g, &o);
    let mut r = Rng::new(12);
    for rev in [false, true] {
        let q = read(&g, 6000, 3000, rev, 0.03, &mut r);
        let hits = map_read(&idx, &q, &o, &none);
        let h = hits.iter().find(|h| h.primary).expect("a mapping");
        let base = h.base.as_ref().expect("an alignment");
        let (mut qi, mut ti, mut n, mut matches) = (h.qstart, h.tstart, 0usize, 0usize);
        let oriented = if h.rev { revcom(&q) } else { q.clone() };
        let (qo, mut qpos) = (oriented, if h.rev { q.len() - h.qend } else { h.qstart });
        for c in base.cigar.chars() {
            match c {
                '0'..='9' => n = n * 10 + c.to_digit(10).unwrap() as usize,
                'M' => {
                    matches += (0..n).filter(|&i| qo[qpos + i] == g[ti + i]).count();
                    (qpos, ti, qi, n) = (qpos + n, ti + n, qi + n, 0);
                }
                'I' => (qpos, qi, n) = (qpos + n, qi + n, 0),
                _ => (ti, n) = (ti + n, 0),
            }
        }
        assert_eq!(matches, base.matches);
        assert_eq!(qi, h.qend);
        assert_eq!(ti, h.tend);
        assert!(base.matches as f64 / base.block as f64 > 0.8);
    }
}

#[test]
fn the_bases_the_anchors_cover_are_the_union_of_their_k_mers() {
    let a = |ys: &[usize]| -> Vec<Anchor> { ys.iter().map(|&y| Anchor { x: y, y }).collect() };
    // 0 to 20 (two overlapping k-mers) and 40 to 55
    assert_eq!(covered(&a(&[0, 5, 40]), 15), 35);
    assert_eq!(covered(&a(&[40, 0, 5]), 15), 35);
    assert_eq!(covered(&a(&[10]), 15), 15);
    assert_eq!(covered(&[], 15), 0);
    assert_eq!(covered(&a(&[0, 15, 30]), 15), 45);
}

#[test]
fn the_divergence_of_a_chain_is_small_for_an_exact_read_and_grows_with_the_errors() {
    let g = genome(30_000, 91);
    let o = Options::default();
    let idx = index(&g, &o);
    let exact = map_read(&idx, &g[5000..9000], &o, &none)
        .remove(0)
        .divergence
        .expect("an estimate");
    assert!(exact.abs() < 0.005, "{exact}");
    let mut r = Rng::new(7);
    let mut found = Vec::new();
    for rate in [0.02, 0.05] {
        let q = read(&g, 12_000, 4000, false, rate, &mut r);
        found.push(
            map_read(&idx, &q, &o, &none)
                .remove(0)
                .divergence
                .expect("an estimate"),
        );
    }
    // the more errors, the greater it is, and it is of the size of the errors (6 % and 15 % of the bases are wrong)
    assert!(exact < found[0] && found[0] < found[1], "{exact} {found:?}");
    assert!((0.03..0.2).contains(&found[0]), "{found:?}");
    assert!((0.08..0.4).contains(&found[1]), "{found:?}");
}
