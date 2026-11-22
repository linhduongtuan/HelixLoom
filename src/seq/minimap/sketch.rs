//! (w,k)-minimizers: of every window of `w` consecutive k-mers the one with the smallest hash is kept (Roberts et al. 2004, as minimap2's
//! paper uses them). A k-mer is taken as the smaller, as a number, of itself and its reverse complement, so that a sequence and its reverse
//! complement give the same minimizers, and the strand each had is kept. A k-mer that is its own reverse complement has no strand and is
//! not taken. Letters other than A, C, G, T (U counts as T) break the sequence into runs; a run is sketched on its own, and a run with fewer
//! than `w` k-mers gives the one smallest.
//!
//! The hash is the finalizer of SplitMix64: any one-to-one mixing of 64 bits would do, and the paper does not name its own. When several
//! k-mers of a window share the smallest hash, all of them are kept (so the sketch of the reverse complement is the mirror image).

use std::collections::VecDeque;

/// The longest k-mer: two bits a base in 64.
pub const MAX_K: usize = 32;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Minimizer {
    pub hash: u64,
    /// Where the k-mer starts, counted from 0.
    pub pos: usize,
    /// The k-mer as it is written in the sequence is the reverse complement of the one that was hashed.
    pub rev: bool,
}

fn mix(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// The hash and strand of each k-mer of one run of A, C, G, T; none for a k-mer that is its own reverse complement.
fn hashed(run: &[u8], k: usize) -> Vec<Option<(u64, bool)>> {
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let (mut fwd, mut rev) = (0u64, 0u64);
    let mut out = Vec::with_capacity(run.len().saturating_sub(k - 1));
    for (i, &b) in run.iter().enumerate() {
        let c = code(b).unwrap_or(0);
        fwd = ((fwd << 2) | c) & mask;
        rev = (rev >> 2) | ((3 - c) << (2 * (k - 1)));
        if i + 1 >= k {
            out.push((fwd != rev).then(|| (mix(fwd.min(rev)), rev < fwd)));
        }
    }
    out
}

/// The minimizers of `seq` for k-mers of `k` bases and windows of `w` k-mers, by position.
pub fn sketch(seq: &[u8], k: usize, w: usize) -> Vec<Minimizer> {
    let (k, w) = (k.clamp(1, MAX_K), w.max(1));
    let mut out = Vec::new();
    let mut start = 0;
    while start < seq.len() {
        if code(seq[start]).is_none() {
            start += 1;
            continue;
        }
        let end = seq[start..]
            .iter()
            .position(|&b| code(b).is_none())
            .map_or(seq.len(), |n| start + n);
        if end - start >= k {
            run_minimizers(&hashed(&seq[start..end], k), start, w, &mut out);
        }
        start = end;
    }
    out
}

/// The minimizers of one run, whose k-mers are `kmers`, the first at `offset` in the sequence.
fn run_minimizers(
    kmers: &[Option<(u64, bool)>],
    offset: usize,
    w: usize,
    out: &mut Vec<Minimizer>,
) {
    let n = kmers.len();
    // the smallest hashes of the window, in order, with the equal ones kept: the front is the smallest
    let mut window: VecDeque<(u64, usize, bool)> = VecDeque::new();
    let mut last: Option<usize> = None;
    let mut emit = |window: &VecDeque<(u64, usize, bool)>, last: &mut Option<usize>| {
        let Some(&(best, ..)) = window.front() else {
            return;
        };
        for &(h, i, rev) in window.iter().take_while(|e| e.0 == best) {
            if last.is_none_or(|l| i > l) {
                out.push(Minimizer {
                    hash: h,
                    pos: offset + i,
                    rev,
                });
            }
        }
        if let Some(m) = window.iter().take_while(|e| e.0 == best).map(|e| e.1).max() {
            *last = Some(last.map_or(m, |l| l.max(m)));
        }
    };
    for (i, kmer) in kmers.iter().enumerate() {
        if let Some((h, rev)) = *kmer {
            while window.back().is_some_and(|b| b.0 > h) {
                window.pop_back();
            }
            window.push_back((h, i, rev));
        }
        while window.front().is_some_and(|f| f.1 + w <= i) {
            window.pop_front();
        }
        if i + 1 >= w {
            emit(&window, &mut last);
        }
    }
    // a run with fewer k-mers than a window has one window: the whole of it
    if n < w {
        emit(&window, &mut last);
    }
}

#[cfg(test)]
mod tests;
