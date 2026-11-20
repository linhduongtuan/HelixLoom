//! K-mer analysis in the manner of KAT, the K-mer Analysis Toolkit (Mapleson, Garcia Accinelli, Kettleborough, Wright and Clavijo,
//! Bioinformatics 2017, 33:574): count the k-mers of reads or of an assembly, and read what the counts say.
//!
//! Written from the paper and the manual (`kat.readthedocs.io`); KAT is GPL licensed, and none of its source was read or copied. What is
//! here:
//!
//! * **counting**: every window of `k` bases of A, C, G, T (`U` counts as T, any other letter breaks the window: a k-mer holding it is
//!   not counted), canonical, that is, a k-mer and its reverse complement are one, whichever is the smaller as a number;
//! * **`hist`**: the spectrum, how many distinct k-mers occur once, twice, and so on, with a last bin that holds all the rest;
//! * **`gcp`**: the same split by the number of G and C in the k-mer, a matrix;
//! * **`comp`**: the multiplicity in one set against the multiplicity in another, a matrix with a row and a column for 0 (what one has
//!   and the other has not);
//! * **`sect`**: the multiplicity of every k-mer along each sequence of a file, with the median, mean and share of k-mers that occur;
//! * **a reading of the spectrum**: the valley after the errors, the peak of the genome, and the size that follows from them.
//!
//! What the paper describes in words (the error peak at low multiplicity, the heterozygous peak at half the homozygous one, the repeats
//! above it) is read here by the plain rules in `summarize`; the fit of distributions that KAT's `kat_distanalysis.py` makes is not
//! described in what was read, and is not done. Counting is in memory: bacterial genomes and plasmids, and small eukaryotes, not a human
//! genome.

use super::par::par_map;
use std::collections::HashMap;

pub mod cli;
#[cfg(test)]
mod tests;

/// The longest k-mer: two bits a base in 64.
pub const MAX_K: usize = 32;

/// The number of a base, none for any other letter.
fn code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' | b'U' | b'u' => Some(3),
        _ => None,
    }
}

/// Why `k` is not a length of k-mer that can be counted.
pub fn check_k(k: usize) -> Result<(), String> {
    if (1..=MAX_K).contains(&k) {
        Ok(())
    } else {
        Err(format!("k must be from 1 to {MAX_K}, not {k}."))
    }
}

/// One entry for each window of `k` bases of `seq`, in order: the canonical number of its k-mer, or none for a window that holds a
/// letter other than A, C, G, T and U.
pub fn windows(seq: &[u8], k: usize) -> Vec<Option<u64>> {
    if check_k(k).is_err() || seq.len() < k {
        return Vec::new();
    }
    let mut out = vec![None; seq.len() - k + 1];
    let mask = if k == MAX_K {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let (mut fwd, mut rev, mut run) = (0u64, 0u64, 0usize);
    for (j, &b) in seq.iter().enumerate() {
        match code(b) {
            Some(c) => {
                fwd = ((fwd << 2) | c) & mask;
                rev = (rev >> 2) | ((3 - c) << (2 * (k - 1)));
                run += 1;
                if run >= k {
                    out[j + 1 - k] = Some(fwd.min(rev));
                }
            }
            None => (fwd, rev, run) = (0, 0, 0),
        }
    }
    out
}

/// The number of G and C in a k-mer of `k` bases.
fn gc_of(kmer: u64, k: usize) -> usize {
    (0..k)
        .filter(|i| matches!((kmer >> (2 * i)) & 3, 1 | 2))
        .count()
}

/// How often each canonical k-mer occurs in a set of sequences.
#[derive(Clone, Debug, Default)]
pub struct Counts {
    pub k: usize,
    map: HashMap<u64, u32>,
}

impl Counts {
    /// Count the k-mers of `seqs` (on all the cores, a share of the sequences each).
    pub fn count(seqs: &[&[u8]], k: usize) -> Result<Counts, String> {
        check_k(k)?;
        let groups: Vec<&[&[u8]]> = seqs.chunks(seqs.len().div_ceil(16).max(1)).collect();
        let parts = par_map(&groups, |group| {
            let mut local: HashMap<u64, u32> = HashMap::new();
            for seq in group.iter() {
                for kmer in windows(seq, k).into_iter().flatten() {
                    let n = local.entry(kmer).or_insert(0);
                    *n = n.saturating_add(1);
                }
            }
            local
        });
        let mut map: HashMap<u64, u32> = HashMap::new();
        for part in parts {
            for (kmer, n) in part {
                let slot = map.entry(kmer).or_insert(0);
                *slot = slot.saturating_add(n);
            }
        }
        Ok(Counts { k, map })
    }

    /// How often the k-mer (given as its canonical number) occurs.
    pub fn get(&self, kmer: u64) -> u32 {
        self.map.get(&kmer).copied().unwrap_or(0)
    }

    pub fn distinct(&self) -> usize {
        self.map.len()
    }

    /// All the k-mers counted, with repeats.
    pub fn total(&self) -> u64 {
        self.map.values().map(|&n| u64::from(n)).sum()
    }
}

/// The spectrum: entry `m` is how many distinct k-mers occur `m` times, for `low <= m <= high`; those that occur more often are all in
/// the entry `high`, and those that occur less often are left out. Entry 0 is always 0.
pub fn spectrum(c: &Counts, low: u32, high: u32) -> Vec<u64> {
    let mut bins = vec![0u64; high as usize + 1];
    for &n in c.map.values() {
        if n >= low.max(1) {
            bins[n.min(high) as usize] += 1;
        }
    }
    bins
}

/// The matrix of `gcp`: row `g` is the k-mers with `g` bases G or C (0 to k), column `m` how many of them occur `m` times (1 to `high`,
/// the last holding those that occur more often); column 0 is always 0.
pub fn gc_matrix(c: &Counts, high: u32) -> Vec<Vec<u64>> {
    let mut rows = vec![vec![0u64; high as usize + 1]; c.k + 1];
    for (&kmer, &n) in &c.map {
        rows[gc_of(kmer, c.k)][n.min(high) as usize] += 1;
    }
    rows
}

/// The matrix of `comp`: row `i`, column `j` is how many distinct k-mers occur `i` times in `a` and `j` times in `b` (row 0 the k-mers
/// only `b` has, column 0 those only `a` has; the last row and column hold all the larger). Both must be of the same k.
pub fn comparison(
    a: &Counts,
    b: &Counts,
    high_a: u32,
    high_b: u32,
) -> Result<Vec<Vec<u64>>, String> {
    if a.k != b.k {
        return Err(format!(
            "The two sets were counted with different k ({} and {}).",
            a.k, b.k
        ));
    }
    let mut rows = vec![vec![0u64; high_b as usize + 1]; high_a as usize + 1];
    for (&kmer, &n) in &a.map {
        rows[n.min(high_a) as usize][b.get(kmer).min(high_b) as usize] += 1;
    }
    for (&kmer, &n) in &b.map {
        if !a.map.contains_key(&kmer) {
            rows[0][n.min(high_b) as usize] += 1;
        }
    }
    Ok(rows)
}

/// What `sect` says of one sequence: how often each of its k-mers occurs in the counts (none for a k-mer with a letter that is no base),
/// and the numbers that sum it up.
#[derive(Clone, Debug, PartialEq)]
pub struct Section {
    pub coverage: Vec<Option<u32>>,
    pub median: f64,
    pub mean: f64,
    /// Percent of the bases that are G or C.
    pub gc: f64,
    pub invalid: usize,
    pub non_zero: usize,
}

pub fn section(c: &Counts, seq: &[u8]) -> Section {
    let coverage: Vec<Option<u32>> = windows(seq, c.k)
        .into_iter()
        .map(|w| w.map(|kmer| c.get(kmer)))
        .collect();
    let mut values: Vec<u32> = coverage.iter().flatten().copied().collect();
    values.sort_unstable();
    let median = match values.len() {
        0 => 0.0,
        n if n % 2 == 1 => f64::from(values[n / 2]),
        n => (f64::from(values[n / 2 - 1]) + f64::from(values[n / 2])) / 2.0,
    };
    let mean = if values.is_empty() {
        0.0
    } else {
        values.iter().map(|&v| f64::from(v)).sum::<f64>() / values.len() as f64
    };
    let bases = seq.iter().filter(|&&b| code(b).is_some()).count();
    let gc = seq
        .iter()
        .filter(|&&b| matches!(code(b), Some(1 | 2)))
        .count();
    Section {
        invalid: coverage.iter().filter(|v| v.is_none()).count(),
        non_zero: values.iter().filter(|&&v| v > 0).count(),
        median,
        mean,
        gc: if bases == 0 {
            0.0
        } else {
            100.0 * gc as f64 / bases as f64
        },
        coverage,
    }
}

/// What the spectrum says, by the plain rules described under `summarize`.
#[derive(Clone, Debug, PartialEq)]
pub struct Summary {
    pub total: u64,
    pub distinct: u64,
    /// The multiplicity at the bottom of the valley after the errors, and at the top of the peak after it.
    pub valley: Option<usize>,
    pub peak: Option<usize>,
    /// The k-mers that occur from the valley on, over the peak: the size of the genome, when the reads cover it once on each copy.
    pub genome_size: Option<u64>,
    /// The share of the distinct k-mers that occur fewer times than the valley says (most of them are errors).
    pub error_share: f64,
    /// A second peak at about half of the first: a genome with two different copies of much of it.
    pub half_peak: Option<usize>,
}

/// Read a spectrum (`spectrum`'s entries). The counts are smoothed over three entries; the **valley** is where they stop falling after
/// the errors (the first rise, from multiplicity 2); the **peak** is the highest smoothed entry from the valley on; the **genome size**
/// is the sum of multiplicity times count from the valley on, over the peak (the textbook estimate: how many k-mers of the genome there
/// are in the reads, over how often each is read; k-mers in the last bin are counted at its multiplicity, which underestimates a genome
/// with long repeats); and a **half peak** is a local maximum of the smoothed counts within one of half the peak that is at least a tenth
/// as high as it. A spectrum that never rises after its start (no peak: too few reads) has no valley and no size.
pub fn summarize(spec: &[u64]) -> Summary {
    let total: u64 = spec.iter().enumerate().map(|(m, &n)| m as u64 * n).sum();
    let distinct: u64 = spec.iter().sum();
    let smooth: Vec<f64> = (0..spec.len())
        .map(|m| {
            let (a, b) = (m.saturating_sub(1).max(1), (m + 1).min(spec.len() - 1));
            if a > b || m == 0 {
                return 0.0;
            }
            spec[a..=b].iter().map(|&n| n as f64).sum::<f64>() / (b - a + 1) as f64
        })
        .collect();
    let none = Summary {
        total,
        distinct,
        valley: None,
        peak: None,
        genome_size: None,
        error_share: 0.0,
        half_peak: None,
    };
    let top = spec.len().saturating_sub(1);
    let mut valley = 2;
    while valley < top && smooth[valley + 1] <= smooth[valley] {
        valley += 1;
    }
    if valley >= top {
        return none;
    }
    let peak = (valley + 1..=top).fold(
        valley + 1,
        |best, m| if smooth[m] > smooth[best] { m } else { best },
    );
    let genome_kmers: u64 = (valley..=top).map(|m| m as u64 * spec[m]).sum();
    let below: u64 = (1..valley).map(|m| spec[m]).sum();
    let half_target = peak / 2;
    let half_peak = (half_target.saturating_sub(1).max(valley + 1)..=(half_target + 1).min(top))
        .filter(|&m| {
            m > valley
                && m < peak
                && m < top
                && smooth[m] >= smooth[m - 1]
                && smooth[m] >= smooth[m + 1]
        })
        .filter(|&m| smooth[m] >= 0.1 * smooth[peak])
        .max_by(|&a, &b| smooth[a].total_cmp(&smooth[b]));
    Summary {
        valley: Some(valley),
        peak: Some(peak),
        genome_size: Some(genome_kmers / peak as u64),
        error_share: if distinct == 0 {
            0.0
        } else {
            below as f64 / distinct as f64
        },
        half_peak,
        ..none
    }
}
