use super::*;
use crate::seq::minimap::sim::{genome, read, Rng};

fn plain(name: &str, seq: Vec<u8>) -> Read {
    Read {
        name: name.into(),
        seq,
        qual: Vec::new(),
    }
}

fn index_of(seqs: &[&[u8]], o: &Options) -> Index {
    Index::build(
        seqs.iter()
            .enumerate()
            .map(|(i, s)| (format!("r{i}"), s.to_vec()))
            .collect(),
        o,
    )
    .unwrap()
}

#[test]
fn the_threshold_is_the_mean_less_five_deviations_and_never_under_the_floor() {
    // 0.8 and 0.9: a mean of 0.85 and a deviation of 0.05
    assert!((threshold(&[0.8, 0.9], 5.0, 0.5) - 0.6).abs() < 1e-9);
    assert!((threshold(&[0.8, 0.9], 1.0, 0.5) - 0.8).abs() < 1e-9);
    // when all are alike there is no deviation
    assert!((threshold(&[0.9, 0.9, 0.9], 5.0, 0.5) - 0.9).abs() < 1e-9);
    // a loose set goes no lower than the floor
    assert_eq!(threshold(&[0.55, 0.65], 5.0, 0.5), 0.5);
    // fewer than two: the floor
    assert_eq!(threshold(&[0.99], 5.0, 0.5), 0.5);
    assert_eq!(threshold(&[], 5.0, 0.42), 0.42);
}

#[test]
fn the_candidates_are_the_other_reads_that_chain_best_first_one_to_a_read() {
    let g = genome(12_000, 11);
    let template = g[2000..6000].to_vec();
    let seqs: [&[u8]; 4] = [&template, &g[5000..7000], &g[2500..5500], &genome(3000, 99)];
    let p = Params::default();
    let all = index_of(&seqs, &p.overlap());
    let found = candidates(&all, 0, &template, &p.overlap());
    let ids: Vec<usize> = found.iter().map(|c| c.read).collect();
    assert_eq!(ids, [2, 1], "{found:?}");
    assert!(found[0].score > found[1].score);
    // a read that chains to the template in two places is one candidate
    let swapped = [&g[4500..5500], &g[2500..3500]].concat();
    let two: [&[u8]; 2] = [&template, &swapped];
    let all_two = index_of(&two, &p.overlap());
    let found_two = candidates(&all_two, 0, &template, &p.overlap());
    assert_eq!(found_two.len(), 1, "{found_two:?}");
    // the second begins about 3000 bases along the template and runs to its end
    assert!(
        found[1].from >= 2900 && found[1].to > 3900,
        "{:?}",
        found[1]
    );
}

#[test]
fn an_alignment_says_where_it_lies_how_alike_it_is_and_whether_it_is_long() {
    let g = genome(10_000, 12);
    let template = g[0..3000].to_vec();
    let p = Params::default();
    let own = index_of(&[&template], &p.mapping());
    let part = plain("a", g[500..2500].to_vec());
    let a = align(&own, &part, &p.mapping(), &p).unwrap();
    assert!((a.identity - 1.0).abs() < 1e-9, "{}", a.identity);
    assert!(a.from >= 500 && a.from < 600 && a.to > 2400 && a.to <= 2500);
    assert!(a.long);
    assert_eq!(a.record.target.as_deref(), Some("r0"));
    // a read of which a fifth lies on a third of the template is not long
    let mut tail = g[0..1000].to_vec();
    tail.extend(genome(4000, 13));
    let b = align(&own, &plain("b", tail), &p.mapping(), &p).unwrap();
    assert!(!b.long, "{b:?}");
    // a short read that lies wholly on the template is long: it covers itself
    let c = align(&own, &plain("c", g[1000..1600].to_vec()), &p.mapping(), &p).unwrap();
    assert!(c.long);
    // a read that covers most of the template but less than half of itself is long for the template
    let mut over = g[0..2500].to_vec();
    over.extend(genome(3500, 16));
    let e = align(&own, &plain("e", over), &p.mapping(), &p).unwrap();
    assert!(e.long, "{e:?}");
    // one that does not chain gives nothing
    assert!(align(&own, &plain("d", genome(2000, 14)), &p.mapping(), &p).is_none());
}

/// A template, and `n` reads that are copies of it with errors at `rate`.
fn copies(n: usize, rate: f64) -> (Vec<u8>, Vec<Read>) {
    let g = genome(3000, 15);
    let mut r = Rng::new(8);
    let reads = (0..n)
        .map(|i| plain(&format!("c{i}"), read(&g, 0, 3000, false, rate, &mut r)))
        .collect();
    (g, reads)
}

fn chosen(template: &[u8], reads: &[Read], p: &Params) -> Vec<Aligned> {
    let mut seqs: Vec<&[u8]> = vec![template];
    seqs.extend(reads.iter().map(|r| r.seq.as_slice()));
    let all = index_of(&seqs, &p.overlap());
    let found = candidates(&all, 0, template, &p.overlap());
    let own = index_of(&[template], &p.mapping());
    // the reads are numbered from 1 in the index of all, the template being 0
    let mut with_template = vec![plain("t", template.to_vec())];
    with_template.extend_from_slice(reads);
    supports(&own, &found, &with_template, &p.mapping(), p)
}

#[test]
fn a_stretch_with_more_reads_than_the_coverage_takes_no_more() {
    let (g, reads) = copies(30, 0.0);
    let p = Params::default();
    // every read covers the whole template: after thirteen it has more than twelve everywhere
    let got = chosen(&g, &reads, &p);
    assert_eq!(got.len(), 13);
    let five = Params {
        coverage: 5,
        ..Params::default()
    };
    assert_eq!(chosen(&g, &reads, &five).len(), 6);
}

#[test]
fn only_the_best_candidates_are_looked_at() {
    let (g, reads) = copies(30, 0.0);
    let p = Params {
        max_candidates: 4,
        ..Params::default()
    };
    assert_eq!(chosen(&g, &reads, &p).len(), 4);
}

#[test]
fn an_alignment_under_the_templates_threshold_is_not_a_support() {
    // twenty exact copies and one with a quarter of its bases wrong
    let (g, mut reads) = copies(20, 0.0);
    let mut r = Rng::new(21);
    reads.push(plain("bad", read(&g, 0, 3000, false, 0.08, &mut r)));
    let tight = Params {
        sigmas: 1.0,
        coverage: 100,
        ..Params::default()
    };
    let got = chosen(&g, &reads, &tight);
    assert_eq!(
        got.len(),
        20,
        "{:?}",
        got.iter().map(|a| a.identity).collect::<Vec<_>>()
    );
    // a floor over every identity leaves none
    let none = Params {
        min_identity: 1.01,
        ..Params::default()
    };
    assert!(chosen(&g, &reads, &none).is_empty());
    // and with five deviations the bad one, which the mean is not far from, comes in
    let loose = Params {
        coverage: 100,
        ..Params::default()
    };
    assert_eq!(chosen(&g, &reads, &loose).len(), 21);
}

#[test]
fn only_long_alignments_are_measured_for_the_threshold() {
    // twenty exact copies, and five reads with errors that cover a third of the template and a third of themselves
    let (g, mut reads) = copies(20, 0.0);
    let mut r = Rng::new(22);
    for i in 0..5 {
        let mut seq = genome(2000, 300 + i);
        seq.extend(read(&g, 0, 1000, false, 0.05, &mut r));
        reads.push(plain(&format!("short{i}"), seq));
    }
    let p = Params {
        coverage: 100,
        ..Params::default()
    };
    // the long ones are all exact: nothing with an error is alike enough
    let got = chosen(&g, &reads, &p);
    assert_eq!(
        got.len(),
        20,
        "{:?}",
        got.iter().map(|a| (a.identity, a.long)).collect::<Vec<_>>()
    );
}

#[test]
fn a_candidate_is_passed_over_only_when_all_of_its_span_is_covered() {
    let g = genome(3000, 17);
    // ten reads over the first 2500 bases, then three over the last 1000: the first six of the ten are taken (coverage 5), and the three
    // reach positions that no read covers yet
    let mut reads: Vec<Read> = (0..10)
        .map(|i| plain(&format!("a{i}"), g[0..2500].to_vec()))
        .collect();
    reads.extend((0..3).map(|i| plain(&format!("b{i}"), g[2000..3000].to_vec())));
    let p = Params {
        coverage: 5,
        ..Params::default()
    };
    let got = chosen(&g, &reads, &p);
    assert_eq!(
        got.len(),
        9,
        "{:?}",
        got.iter().map(|a| (a.from, a.to)).collect::<Vec<_>>()
    );
}
