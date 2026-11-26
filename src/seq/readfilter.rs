//! The selection of long reads that the multiplex-MinION paper (Wick et al., *Microb Genom* 2017, 3:e000132, Methods, "Base-calling and read
//! preparation") describes: each set is brought down to a target number of bases by dropping the reads under a length first, and then,
//! if there are still too many bases, the reads that have the worst stretch of quality, judged by the mean Phred quality over a window of 100
//! bases, until the target is met. The paper's script is not available to this program (and was not looked for): the paper's numbers (2 kb,
//! 500 Mbp, a 100 base window) are the defaults, and *how* a read's worst stretch is scored (the lowest mean over any window) and how a tie
//! is broken (the later read goes first) are this program's own.

use crate::seq::fastq::FastqRead;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Filter {
    /// Reads shorter than this are dropped, whatever the bases.
    pub min_length: usize,
    /// The most bases to keep; no more reads are dropped once the kept are this many or fewer.
    pub target_bases: u64,
    /// The length of the window whose mean quality scores a read.
    pub window: usize,
}

impl Default for Filter {
    fn default() -> Self {
        Filter {
            min_length: 2000,
            target_bases: 500_000_000,
            window: 100,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Filtered {
    pub reads: Vec<FastqRead>,
    pub too_short: usize,
    pub low_quality: usize,
    pub bases_in: u64,
    pub bases_out: u64,
}

/// The lowest mean quality over any window of `window` bases of the read (the mean of the whole read when it is shorter than the window).
pub fn worst_window(qual: &[u8], window: usize) -> f64 {
    let w = window.max(1).min(qual.len().max(1));
    let mut sum: u64 = qual.iter().take(w).map(|&q| u64::from(q)).sum();
    let mut lowest = sum;
    for i in w..qual.len() {
        sum = sum + u64::from(qual[i]) - u64::from(qual[i - w]);
        lowest = lowest.min(sum);
    }
    lowest as f64 / w as f64
}

pub fn filter_long(reads: &[FastqRead], f: &Filter) -> Filtered {
    let bases_in: u64 = reads.iter().map(|r| r.seq.len() as u64).sum();
    let long: Vec<usize> = (0..reads.len())
        .filter(|&i| reads[i].seq.len() >= f.min_length)
        .collect();
    let too_short = reads.len() - long.len();
    let mut total: u64 = long.iter().map(|&i| reads[i].seq.len() as u64).sum();
    let mut keep = vec![false; reads.len()];
    long.iter().for_each(|&i| keep[i] = true);
    let mut low_quality = 0;
    if total > f.target_bases {
        // the worst first; of two as bad, the later read
        let mut order = long.clone();
        order.sort_by(|&a, &b| {
            worst_window(&reads[a].qual, f.window)
                .total_cmp(&worst_window(&reads[b].qual, f.window))
                .then(b.cmp(&a))
        });
        for i in order {
            if total <= f.target_bases {
                break;
            }
            keep[i] = false;
            total -= reads[i].seq.len() as u64;
            low_quality += 1;
        }
    }
    let kept: Vec<FastqRead> = (0..reads.len())
        .filter(|&i| keep[i])
        .map(|i| reads[i].clone())
        .collect();
    Filtered {
        reads: kept,
        too_short,
        low_quality,
        bases_in,
        bases_out: total,
    }
}

#[cfg(test)]
mod tests;
