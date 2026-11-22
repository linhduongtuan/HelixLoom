//! The index of the target sequences: every minimizer of every target, found again by its hash. Minimizers that occur too often (the
//! repeats) are ignored when a read is looked up, as minimap2's `-f` does.

use super::sketch::{sketch, MAX_K};
use super::Options;
use crate::seq::par::par_map;
use std::collections::HashMap;

/// One minimizer of a target.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Seed {
    pub hash: u64,
    pub target: u32,
    /// Where its k-mer starts in the target.
    pub pos: u32,
    pub rev: bool,
}

/// The least and the most that the occurrence limit may be: a minimizer of a few holders is never too frequent, and one with a million is
/// always (the man page's `-U` bounds, 10 and 1,000,000).
const MOST_OCCURRENCES: (usize, usize) = (10, 1_000_000);

/// How many holders a minimizer may have: the number at the foot of the top `fraction` of the minimizers when they are put in order of how often
/// they occur, but never under the least or over the most that the man page's `-U` allows.
fn occurrence_limit(mut counts: Vec<usize>, fraction: f64) -> usize {
    counts.sort_unstable();
    let at = ((counts.len() as f64) * (1.0 - fraction.clamp(0.0, 1.0))).floor() as usize;
    let quantile = counts
        .get(at.min(counts.len().saturating_sub(1)))
        .copied()
        .unwrap_or(0);
    quantile.min(MOST_OCCURRENCES.1).max(MOST_OCCURRENCES.0)
}

pub struct Index {
    pub k: usize,
    pub w: usize,
    targets: Vec<(String, Vec<u8>)>,
    seeds: Vec<Seed>,
    /// Where the seeds of each hash are in `seeds`: first and one past the last.
    table: HashMap<u64, (usize, usize)>,
    limit: usize,
}

impl Index {
    /// Index `targets` (name and bases). `Err` when there is none, a name or a length cannot be used, or `k` is not a length of k-mer.
    pub fn build(targets: Vec<(String, Vec<u8>)>, o: &Options) -> Result<Index, String> {
        if !(1..=MAX_K).contains(&o.k) {
            return Err(format!("k must be from 1 to {MAX_K}, not {}.", o.k));
        }
        if o.w == 0 {
            return Err("w must be 1 or more.".into());
        }
        if targets.is_empty() {
            return Err("There is no target sequence to map to.".into());
        }
        if let Some((name, _)) = targets.iter().find(|t| t.1.len() >= u32::MAX as usize) {
            return Err(format!(
                "{name} is longer than 4 billion bases, which the index cannot hold."
            ));
        }
        let targets: Vec<(String, Vec<u8>)> = targets
            .into_iter()
            .map(|(n, mut s)| {
                s.make_ascii_uppercase();
                (n, s)
            })
            .collect();
        let ids: Vec<usize> = (0..targets.len()).collect();
        let sketches = par_map(&ids, |&t| sketch(&targets[t].1, o.k, o.w));
        let mut seeds: Vec<Seed> = sketches
            .into_iter()
            .enumerate()
            .flat_map(|(t, s)| {
                s.into_iter().map(move |m| Seed {
                    hash: m.hash,
                    target: t as u32,
                    pos: m.pos as u32,
                    rev: m.rev,
                })
            })
            .collect();
        seeds.sort_unstable_by_key(|s| (s.hash, s.target, s.pos));
        let mut table: HashMap<u64, (usize, usize)> = HashMap::new();
        let mut counts: Vec<usize> = Vec::new();
        let mut from = 0;
        while from < seeds.len() {
            let hash = seeds[from].hash;
            let to = from + seeds[from..].iter().take_while(|s| s.hash == hash).count();
            table.insert(hash, (from, to));
            counts.push(to - from);
            from = to;
        }
        let limit = occurrence_limit(counts, o.occ_fraction);
        Ok(Index {
            k: o.k,
            w: o.w,
            targets,
            seeds,
            table,
            limit,
        })
    }

    #[cfg(test)]
    pub fn targets(&self) -> usize {
        self.targets.len()
    }

    pub fn name(&self, target: usize) -> &str {
        &self.targets[target].0
    }

    pub fn length(&self, target: usize) -> usize {
        self.targets[target].1.len()
    }

    pub fn sequence(&self, target: usize) -> &[u8] {
        &self.targets[target].1
    }

    /// How many minimizers there are in all.
    #[cfg(test)]
    pub fn minimizers(&self) -> usize {
        self.seeds.len()
    }

    /// A minimizer with more holders than this is ignored.
    #[cfg(test)]
    pub fn limit(&self) -> usize {
        self.limit
    }

    /// The seeds of one hash: none when it has no holder or has too many.
    pub fn hits(&self, hash: u64) -> &[Seed] {
        match self.table.get(&hash) {
            Some(&(a, b)) if b - a <= self.limit => &self.seeds[a..b],
            _ => &[],
        }
    }
}

#[cfg(test)]
mod tests;
