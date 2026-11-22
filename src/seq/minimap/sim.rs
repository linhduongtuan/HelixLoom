//! Reads simulated from a known sequence, for the tests: the way to know where a read belongs without another program's answer.

use crate::seq::dna::revcom;

/// A small deterministic generator (SplitMix64), so that the tests do not change from one run to the next.
pub struct Rng(u64);

impl Rng {
    /// Seeds that are near each other give streams that have nothing in common (a stream that began at `seed` times the step would be the
    /// stream of the next seed one place on, and a "different" genome would hold the other's).
    pub fn new(seed: u32) -> Rng {
        let mut z = u64::from(seed) ^ 0xD1B5_4A32_D192_ED03;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        Rng(z ^ (z >> 31))
    }

    pub fn next(&mut self) -> u32 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 32) as u32
    }

    pub fn below(&mut self, n: usize) -> usize {
        self.next() as usize % n.max(1)
    }

    pub fn chance(&mut self, p: f64) -> bool {
        f64::from(self.next() % 1_000_000) < p * 1_000_000.0
    }
}

pub fn genome(n: usize, seed: u32) -> Vec<u8> {
    let mut r = Rng::new(seed);
    (0..n).map(|_| b"ACGT"[r.below(4)]).collect()
}

/// `segment` as a read would have it: each base replaced with `sub`, dropped with `del`, and followed by an extra base with `ins`.
pub fn noisy(segment: &[u8], sub: f64, ins: f64, del: f64, r: &mut Rng) -> Vec<u8> {
    let mut out = Vec::with_capacity(segment.len());
    for &b in segment {
        if r.chance(del) {
            continue;
        }
        out.push(if r.chance(sub) {
            *b"ACGT"
                .iter()
                .find(|&&c| c != b && r.chance(0.5))
                .unwrap_or(&b'A')
        } else {
            b
        });
        if r.chance(ins) {
            out.push(b"ACGT"[r.below(4)]);
        }
    }
    out
}

/// A read from `genome[start..start + len]`, with the errors of nanopore reads of about 10 % (`rate` is the share of each kind of error),
/// on the reverse strand when `rev`.
pub fn read(genome: &[u8], start: usize, len: usize, rev: bool, rate: f64, r: &mut Rng) -> Vec<u8> {
    let part = noisy(&genome[start..start + len], rate, rate, rate, r);
    if rev {
        revcom(&part)
    } else {
        part
    }
}
