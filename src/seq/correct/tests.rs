use super::*;
use crate::seq::dna::revcom;
use crate::seq::minimap::map_read;
use crate::seq::minimap::sim::{genome, read, Rng};

/// Reads of `genome`, 1500 to 2000 bases, `cover` times it, with `rate` of each kind of error; and where each came from.
fn simulated(
    g: &[u8],
    cover: usize,
    rate: f64,
    seed: u32,
) -> (Vec<Read>, Vec<(usize, usize, bool)>) {
    let mut r = Rng::new(seed);
    let n = g.len() * cover / 1750;
    let mut origins = Vec::new();
    let reads = (0..n)
        .map(|i| {
            let len = 1500 + r.below(500);
            let start = r.below(g.len() - len);
            origins.push((start, len, i % 2 == 0));
            let seq = read(g, start, len, i % 2 == 0, rate, &mut r);
            Read {
                name: format!("r{i}"),
                qual: vec![15; seq.len()],
                seq,
            }
        })
        .collect();
    (reads, origins)
}

/// The share of the aligned columns of each read against the bases it came from that match, over all reads.
fn identity(g: &[u8], set: &[Read], origins: &[(usize, usize, bool)]) -> f64 {
    let o = Options {
        cigar: true,
        ..Options::default()
    };
    let (mut matches, mut block) = (0usize, 0usize);
    for r in set {
        let i: usize = r.name[1..].split('_').next().unwrap().parse().unwrap();
        let (start, len, reverse) = origins[i];
        let truth = &g[start..start + len];
        let truth = if reverse {
            revcom(truth)
        } else {
            truth.to_vec()
        };
        let index = Index::build(vec![("t".to_string(), truth)], &o).unwrap();
        if let Some(h) = map_read(&index, &r.seq, &o, &|_| false)
            .into_iter()
            .find(|h| h.primary)
        {
            let base = h.base.unwrap();
            matches += base.matches;
            block += base.block;
        }
    }
    matches as f64 / block.max(1) as f64
}

#[test]
fn the_defaults_are_the_papers_numbers_and_this_programs_choices() {
    let p = Params::default();
    assert_eq!(
        (
            p.rounds,
            p.measured,
            p.sigmas,
            p.min_identity,
            p.max_candidates,
            p.coverage
        ),
        (2, 50, 5.0, 0.5, 200, 12)
    );
    assert_eq!((p.long_share, p.block, p.max_block_error), (0.6, 500, 0.5));
    assert_eq!((p.min_support, p.min_piece, p.k, p.w), (4, 500, 15, 5));
    let overlap = p.overlap();
    assert!(overlap.overlaps && !overlap.cigar);
    assert_eq!((overlap.k, overlap.w), (15, 5));
    let mapping = p.mapping();
    assert!(mapping.cigar && !mapping.overlaps);
}

#[test]
fn reads_of_a_genome_with_errors_come_out_with_far_fewer() {
    let g = genome(8000, 101);
    let (reads, origins) = simulated(&g, 30, 0.04, 3);
    let before = identity(&g, &reads, &origins);
    let out = correct(&reads, &Params::default()).unwrap();
    let after = identity(&g, &out, &origins);
    assert!(before < 0.92, "{before}");
    assert!(after > 0.99, "{before} to {after}");
    // most of the bases come out, and a read whose pieces were all kept has the name it had
    let bases: usize = out.iter().map(|r| r.seq.len()).sum();
    let had: usize = reads.iter().map(|r| r.seq.len()).sum();
    assert!(bases * 10 > had * 8, "{bases} of {had}");
    assert!(out.iter().all(|r| r.qual.is_empty() && r.seq.len() >= 500));
}

#[test]
fn one_round_keeps_every_read_and_its_name() {
    let g = genome(6000, 102);
    let (reads, origins) = simulated(&g, 25, 0.04, 4);
    let one = Params {
        rounds: 1,
        ..Params::default()
    };
    let out = correct(&reads, &one).unwrap();
    assert_eq!(out.len(), reads.len());
    assert!(out.iter().zip(&reads).all(|(a, b)| a.name == b.name));
    assert!(identity(&g, &out, &origins) > identity(&g, &reads, &origins) + 0.05);
}

#[test]
fn a_read_that_overlaps_no_other_is_kept_by_one_round_and_dropped_by_two() {
    let g = genome(6000, 103);
    let (mut reads, _) = simulated(&g, 20, 0.03, 5);
    let alone = Read {
        name: "alone".into(),
        seq: genome(2000, 555),
        qual: Vec::new(),
    };
    reads.push(alone.clone());
    let one = Params {
        rounds: 1,
        ..Params::default()
    };
    let out = correct(&reads, &one).unwrap();
    let kept = out.iter().find(|r| r.name == "alone").unwrap();
    assert_eq!(kept.seq, alone.seq);
    let out = correct(&reads, &Params::default()).unwrap();
    assert!(out.iter().all(|r| !r.name.starts_with("alone")));
}

#[test]
fn a_read_with_an_unsupported_middle_comes_out_in_two_pieces_with_numbers() {
    // two stretches of a genome that other reads cover, joined in one read by bases of nothing
    let g = genome(12_000, 104);
    let (mut reads, _) = simulated(&g[..5000], 20, 0.02, 6);
    let (right, _) = simulated(&g[7000..], 20, 0.02, 7);
    reads.extend(right.into_iter().map(|r| Read {
        name: format!("x{}", r.name),
        ..r
    }));
    let mut joined = g[1000..2500].to_vec();
    joined.extend(genome(800, 777));
    joined.extend_from_slice(&g[8000..9500]);
    reads.push(Read {
        name: "joined".into(),
        seq: joined,
        qual: Vec::new(),
    });
    let out = correct(&reads, &Params::default()).unwrap();
    let names: Vec<&str> = out
        .iter()
        .map(|r| r.name.as_str())
        .filter(|n| n.starts_with("joined"))
        .collect();
    assert_eq!(names, ["joined_1", "joined_2"]);
}

#[test]
fn no_reads_or_no_rounds_is_an_error() {
    assert!(correct(&[], &Params::default())
        .unwrap_err()
        .contains("no reads"));
    let r = Read {
        name: "a".into(),
        seq: genome(100, 1),
        qual: Vec::new(),
    };
    let none = Params {
        rounds: 0,
        ..Params::default()
    };
    assert!(correct(std::slice::from_ref(&r), &none)
        .unwrap_err()
        .contains("one round"));
    let bad = Params {
        k: 40,
        ..Params::default()
    };
    assert!(correct(&[r], &bad).unwrap_err().contains("k must be"));
}
