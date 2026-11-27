//! Samples simulated for the tests: a reference, a sample that differs from it in known places, and reads from the sample with the errors of
//! nanopore reads. The way to know what a caller should find without another program's answer.

use super::reads::Read;
use crate::seq::dna::revcom;
use crate::seq::minimap::sim::{noisy, Rng};

/// The other base that a substitution makes of `b`.
pub fn other(b: u8) -> u8 {
    match b {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

/// `reference` with the substitutions at `snps`, the deletions `(start, length)`, and the insertions `(after, bases)`, each by its place on
/// the reference.
pub fn plant(
    reference: &[u8],
    snps: &[usize],
    deletions: &[(usize, usize)],
    insertions: &[(usize, &[u8])],
) -> Vec<u8> {
    let mut out = Vec::new();
    let mut at = 0;
    while at < reference.len() {
        if let Some(&(_, length)) = deletions.iter().find(|d| d.0 == at) {
            at += length;
            continue;
        }
        let b = reference[at];
        out.push(if snps.contains(&at) { other(b) } else { b });
        if let Some((_, bases)) = insertions.iter().find(|i| i.0 == at) {
            out.extend_from_slice(bases);
        }
        at += 1;
    }
    out
}

/// `n` reads of `length` bases from `sample` (a circle is not simulated: reads lie inside it), either strand, with `rate` of each of the
/// three kinds of error, every base of quality `quality`.
pub fn reads(
    sample: &[u8],
    n: usize,
    length: usize,
    rate: f64,
    quality: u8,
    seed: u32,
) -> Vec<Read> {
    let mut r = Rng::new(seed);
    (0..n)
        .map(|i| {
            let start = r.below(sample.len() - length + 1);
            let seq = noisy(&sample[start..start + length], rate, rate, rate, &mut r);
            let seq = if r.chance(0.5) { revcom(&seq) } else { seq };
            Read {
                name: format!("read{i}"),
                qual: vec![quality; seq.len()],
                seq,
            }
        })
        .collect()
}

/// A record by hand: on the target `t`, from `pos` (0-based), with every base of quality `q`, mapping quality 60.
pub fn record(name: &str, pos: usize, cigar: &str, seq: &str, q: u8) -> super::sam::Record {
    super::sam::Record {
        name: name.to_string(),
        flag: 0,
        target: Some("t".to_string()),
        pos,
        mapq: 60,
        cigar: super::sam::parse_cigar(cigar).unwrap(),
        seq: seq.bytes().collect(),
        qual: vec![q; seq.len()],
        tags: Vec::new(),
    }
}

/// The first place at or after `from` where a deletion of two bases, or an insertion of the base that `other_than` is not, cannot be moved
/// along the reference: its neighbours differ, so that the alignment has one way to put it.
pub fn plain(reference: &[u8], from: usize) -> usize {
    (from.max(3)..reference.len() - 6)
        .find(|&i| {
            reference[i - 1] != reference[i + 1]
                && reference[i] != reference[i + 2]
                && reference[i + 1] != reference[i + 3]
                && reference[i - 2] != reference[i]
        })
        .expect("a plain place")
}
