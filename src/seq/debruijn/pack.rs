//! K-mers packed two bits a base into a `u128` (A 0, C 1, G 2, T 3), so that k can be up to 63, and the table that counts them.

use crate::seq::par::par_map;
use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};

/// The longest k: two bits a base in 128, and an odd length (an odd k-mer is never its own reverse complement).
pub const MAX_K: usize = 63;
pub const MIN_K: usize = 3;

pub type Kmer = u128;

/// A fast hasher for packed k-mers: the standard one is made to resist an attacker who picks the keys, which a table of k-mers does not need.
#[derive(Default)]
pub struct MixHasher(u64);

impl Hasher for MixHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ u64::from(b)).wrapping_mul(0x517c_c1b7_2722_0a95);
        }
    }

    fn write_u128(&mut self, x: u128) {
        let mut z = (x as u64) ^ ((x >> 64) as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        self.0 = z ^ (z >> 31);
    }
}

pub type Table<V> = HashMap<Kmer, V, BuildHasherDefault<MixHasher>>;
pub type Set = HashSet<Kmer, BuildHasherDefault<MixHasher>>;

/// The arithmetic of k-mers of one length.
#[derive(Clone, Copy, Debug)]
pub struct Pack {
    k: usize,
    mask: u128,
}

impl Pack {
    pub fn new(k: usize) -> Result<Pack, String> {
        if !(MIN_K..=MAX_K).contains(&k) || k.is_multiple_of(2) {
            return Err(format!("k must be odd, from {MIN_K} to {MAX_K}, not {k}."));
        }
        Ok(Pack {
            k,
            mask: (1u128 << (2 * k)) - 1,
        })
    }

    pub fn k(&self) -> usize {
        self.k
    }

    /// The number of a base, none for any other letter.
    pub fn code(b: u8) -> Option<u8> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' | b'U' | b'u' => Some(3),
            _ => None,
        }
    }

    /// The k-mer that follows `kmer` when `code` is added at its end.
    pub fn append(&self, kmer: Kmer, code: u8) -> Kmer {
        ((kmer << 2) | u128::from(code)) & self.mask
    }

    /// The k-mer that comes before `kmer` when `code` is put at its start.
    pub fn prepend(&self, kmer: Kmer, code: u8) -> Kmer {
        (kmer >> 2) | (u128::from(code) << (2 * (self.k - 1)))
    }

    /// The reverse complement: the order of the two-bit groups of the whole word is reversed, the k that matter are shifted down, and each
    /// is complemented (T is 3 − A, G is 3 − C: the same as flipping both bits).
    pub fn rc(&self, kmer: Kmer) -> Kmer {
        let mut x = kmer.swap_bytes();
        x = ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F)
            | ((x & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F) << 4);
        x = ((x >> 2) & 0x3333_3333_3333_3333_3333_3333_3333_3333)
            | ((x & 0x3333_3333_3333_3333_3333_3333_3333_3333) << 2);
        (x >> (128 - 2 * self.k)) ^ self.mask
    }

    /// The smaller, as a number, of a k-mer and its reverse complement: the two are one k-mer of the graph.
    pub fn canonical(&self, kmer: Kmer) -> Kmer {
        kmer.min(self.rc(kmer))
    }

    /// The first k bases of `seq` packed; none when there are fewer or one is not A, C, G or T.
    pub fn pack(&self, seq: &[u8]) -> Option<Kmer> {
        let head = seq.get(..self.k)?;
        head.iter().try_fold(0u128, |acc, &b| {
            Pack::code(b).map(|c| (acc << 2) | u128::from(c))
        })
    }

    pub fn unpack(&self, kmer: Kmer) -> Vec<u8> {
        (0..self.k)
            .map(|i| b"ACGT"[((kmer >> (2 * (self.k - 1 - i))) & 3) as usize])
            .collect()
    }

    /// Every window of k bases of `seq` that holds only A, C, G and T, as the k-mer and its reverse complement, in order.
    pub fn for_each(&self, seq: &[u8], mut f: impl FnMut(Kmer, Kmer)) {
        let (mut fwd, mut rev, mut run) = (0u128, 0u128, 0usize);
        for &b in seq {
            match Pack::code(b) {
                Some(c) => {
                    fwd = self.append(fwd, c);
                    rev = self.prepend(rev, 3 - c);
                    run += 1;
                    if run >= self.k {
                        f(fwd, rev);
                    }
                }
                None => run = 0,
            }
        }
    }
}

/// How many times each canonical k-mer occurs in the pieces (reads, on either strand), the work shared among the cores.
pub fn count(pieces: &[Vec<u8>], pack: Pack) -> Table<f32> {
    let workers = std::thread::available_parallelism().map_or(1, |w| w.get());
    let parts: Vec<&[Vec<u8>]> = pieces
        .chunks(pieces.len().div_ceil(workers).max(1))
        .collect();
    let tables = par_map(&parts, |part| {
        let mut table: Table<f32> = Table::default();
        for piece in *part {
            pack.for_each(piece, |f, r| {
                *table.entry(f.min(r)).or_insert(0.0) += 1.0;
            });
        }
        table
    });
    let mut tables = tables.into_iter();
    let mut all = tables.next().unwrap_or_default();
    for table in tables {
        for (kmer, n) in table {
            *all.entry(kmer).or_insert(0.0) += n;
        }
    }
    all
}

#[cfg(test)]
mod tests;
