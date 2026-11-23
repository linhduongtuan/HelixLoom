//! Reads simulated from a known sequence, for the tests.

use crate::seq::dna::revcom;
use crate::seq::minimap::sim::{genome, Rng};

pub fn random(n: usize, seed: u32) -> Vec<u8> {
    genome(n, seed)
}

/// Reads of `len` bases from either strand of `genome` at about `coverage` times, each base wrong with chance `error`; a circle is read
/// across its join. The first and last `len` bases of a line are read at least once, which random starts would leave to chance.
pub fn reads(
    genome: &[u8],
    circular: bool,
    len: usize,
    coverage: f64,
    error: f64,
    seed: u32,
) -> Vec<Vec<u8>> {
    let mut r = Rng::new(seed);
    let source: Vec<u8> = if circular {
        genome.iter().chain(&genome[..len - 1]).copied().collect()
    } else {
        genome.to_vec()
    };
    let starts = genome.len() - if circular { 0 } else { len - 1 };
    let n = (coverage * genome.len() as f64 / len as f64).ceil() as usize;
    let mut out = Vec::new();
    let mut starts_at: Vec<usize> = (0..n).map(|_| r.below(starts)).collect();
    if !circular {
        starts_at.extend([0, starts - 1]);
    }
    for at in starts_at {
        let mut read = source[at..at + len].to_vec();
        for b in read.iter_mut() {
            if r.chance(error) {
                let i = b"ACGT".iter().position(|c| c == b).unwrap_or(0);
                *b = b"ACGT"[(i + 1 + r.below(3)) % 4];
            }
        }
        out.push(if r.chance(0.5) { revcom(&read) } else { read });
    }
    out
}

/// Reads of `len` bases every `step` bases along `genome` (the last one at the end), every other one on the reverse strand, with no errors:
/// the coverage is the same all along.
pub fn tiled(genome: &[u8], len: usize, step: usize) -> Vec<Vec<u8>> {
    let mut starts: Vec<usize> = (0..=genome.len() - len).step_by(step).collect();
    if starts.last() != Some(&(genome.len() - len)) {
        starts.push(genome.len() - len);
    }
    starts
        .into_iter()
        .enumerate()
        .map(|(i, at)| {
            let read = genome[at..at + len].to_vec();
            if i % 2 == 1 {
                revcom(&read)
            } else {
                read
            }
        })
        .collect()
}
