//! A genome with repeats, the graph that short reads give it, and long reads from it, for the tests.

use crate::seq::debruijn::assemble::{assemble, Cutoff, Options};
use crate::seq::debruijn::graph::Graph;
use crate::seq::debruijn::sim::{random, reads};
use crate::seq::minimap::sim::{read, Rng};

/// A circle of `X R Y R Z`, `R` a repeat of `repeat` bases that is twice in it, longer than the reads of a short-read assembler can span.
pub struct Ring {
    pub genome: Vec<u8>,
}

pub fn ring(unit: usize, repeat: usize, seed: u32) -> Ring {
    let x = random(unit, seed);
    let y = random(unit, seed + 1);
    let z = random(unit, seed + 2);
    let r = random(repeat, seed + 3);
    Ring {
        genome: [&x[..], &r, &y, &r, &z].concat(),
    }
}

/// The graph that error-free short reads of `genome` give a de Bruijn assembler, at k = 21, 33 and 55.
pub fn short_graph(genome: &[u8], circular: bool, seed: u32) -> Graph {
    let pieces = reads(genome, circular, 100, 40.0, 0.0, seed);
    let o = Options {
        ks: vec![21, 33, 55],
        cutoff: Cutoff::Off,
        isolated_length: 0,
        ..Options::default()
    };
    assemble(&pieces, &o).expect("assembly").graph
}

/// Long reads from `genome` (a circle when `circular`), of `len` bases each, so many that `coverage` is reached, on either strand, each kind
/// of error at `rate`.
pub fn long_reads(
    genome: &[u8],
    circular: bool,
    len: usize,
    coverage: f64,
    rate: f64,
    seed: u32,
) -> Vec<(String, Vec<u8>)> {
    let mut r = Rng::new(seed);
    let source: Vec<u8> = if circular {
        genome.iter().chain(&genome[..len]).copied().collect()
    } else {
        genome.to_vec()
    };
    let starts = genome.len() - if circular { 0 } else { len };
    let n = (coverage * genome.len() as f64 / len as f64).ceil() as usize;
    (0..n)
        .map(|i| {
            let at = r.below(starts + 1);
            (
                format!("long{i}"),
                read(&source, at, len, r.chance(0.5), rate, &mut r),
            )
        })
        .collect()
}

/// Is `seq` the genome as a circle turned and perhaps flipped?
pub fn is_circle_of(seq: &[u8], genome: &[u8]) -> bool {
    use crate::seq::dna::revcom;
    let twice: Vec<u8> = [genome, genome].concat();
    let find = |hay: &[u8]| hay.windows(seq.len()).any(|w| w == seq);
    seq.len() == genome.len() && (find(&twice) || find(&revcom(&twice)))
}
