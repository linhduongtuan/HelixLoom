//! Mapping reads to a reference and finding overlaps between reads, in the manner of minimap2 (Li, *Bioinformatics* 2018, 34:3094), written
//! from its paper, its manual and its README (minimap2 is MIT licensed, and none of its source was read or copied). What is here:
//!
//! * **minimizers** (`sketch`): the k-mer of the smallest hash in each window of `w` k-mers, canonical, with the strand kept;
//! * **an index** (`index`) of every minimizer of every target, with the frequent ones (repeats) ignored (`-f`);
//! * **chaining** (`chain`): the paper's dynamic programme over the places where a read and a target share a minimizer, with its gap cost,
//!   50 anchors of look-back, and each anchor in one chain;
//! * **mapping** (`map`): primary and secondary chains (a chain that overlaps a better one on the read by half is its secondary) and the
//!   mapping quality of the paper;
//! * **base-level alignment** (`base`, with `-c`): the gaps between anchors aligned by the program's own aligner, written as a CIGAR;
//! * **PAF** (`paf`), and **overlaps** between reads (`-x ava-ont`: each pair once, no read with itself).
//!
//! Not done: the SIMD kernel and the extension of the ends by a drop-off rule (the alignment runs from the first anchor to the last), splice
//! mode, short reads in pairs, homopolymer-compressed minimizers (`-H`, for PacBio CLR), the other presets, and the second occurrence
//! filter that samples frequent minimizers (`-e`). The output is **not** compared with minimap2's, which was not run (that needs a
//! download): the algorithms are the paper's, the numbers have not been set against the program's.

mod base;
mod chain;
pub mod cli;
mod index;
mod map;
mod paf;
#[cfg(test)]
pub(crate) mod sim;
mod sketch;

use index::Index;
use map::map_read;
use sketch::MAX_K;

use crate::seq::par::par_map;

/// How to map. The defaults are minimap2's for noisy long reads (`-x map-ont`).
#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    /// The length of a k-mer and the number of k-mers in a window (`-k`, `-w`).
    pub k: usize,
    pub w: usize,
    /// The share of the most frequent minimizers that are ignored (`-f`).
    pub occ_fraction: f64,
    /// The most that two anchors of a chain may be apart (`-g`), and the most that the gaps of two may differ (`-r`).
    pub max_gap: usize,
    pub bandwidth: usize,
    /// The fewest anchors (`-n`) and the lowest score (`-m`) of a chain.
    pub min_anchors: usize,
    pub min_score: f64,
    /// A secondary mapping must score this share of its primary (`-p`); at most this many are kept (`-N`); and a mapping is secondary to
    /// one that it overlaps on the read by this share of the shorter (`-M`).
    pub secondary_ratio: f64,
    pub max_secondary: usize,
    pub mask_level: f64,
    /// Anchors looked at before one is given up on (the paper's `h`).
    pub lookback: usize,
    /// Align the bases along each chain and write the CIGAR (`-c`).
    pub cigar: bool,
    /// Overlaps between reads (`-X`): every chain is written, none is secondary, each pair once.
    pub overlaps: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            k: 15,
            w: 10,
            occ_fraction: 0.0002,
            max_gap: 10_000,
            bandwidth: 500,
            min_anchors: 3,
            min_score: 40.0,
            secondary_ratio: 0.8,
            max_secondary: 5,
            mask_level: 0.5,
            lookback: 50,
            cigar: false,
            overlaps: false,
        }
    }
}

impl Options {
    /// The options of a preset (`-x`): `map-ont` (the defaults) and `ava-ont` (`-k15 -Xw5 -e0 -m100 -r2k` in the manual).
    pub fn preset(name: &str) -> Option<Options> {
        match name {
            "map-ont" => Some(Options::default()),
            "ava-ont" => Some(Options {
                w: 5,
                min_score: 100.0,
                bandwidth: 2000,
                overlaps: true,
                ..Options::default()
            }),
            _ => None,
        }
    }
}

pub type Named = (String, Vec<u8>);

/// The mappings of every read in `queries` to `targets`, as PAF text. `Err` when there is no target or an option is not one.
pub fn map_text(targets: &[Named], queries: &[Named], o: &Options) -> Result<String, String> {
    let index = Index::build(targets.to_vec(), o)?;
    Ok(run(&index, queries, o, &|_, _| false))
}

/// The overlaps between the reads: each pair of reads that share a chain once (the one that comes first is the query), none for a read with
/// itself, as PAF text.
pub fn overlap_text(reads: &[Named], o: &Options) -> Result<String, String> {
    let index = Index::build(
        reads.to_vec(),
        &Options {
            overlaps: true,
            ..o.clone()
        },
    )?;
    Ok(run(
        &index,
        reads,
        &Options {
            overlaps: true,
            ..o.clone()
        },
        &|query, target| target <= query,
    ))
}

fn run(
    index: &Index,
    queries: &[Named],
    o: &Options,
    skip: &(dyn Fn(usize, usize) -> bool + Sync),
) -> String {
    let ids: Vec<usize> = (0..queries.len()).collect();
    let lines = par_map(&ids, |&i| {
        let (name, bases) = &queries[i];
        map_read(index, bases, o, &|t| skip(i, t))
            .iter()
            .map(|h| paf::line(name, bases.len(), index, h, o.overlaps))
            .collect::<Vec<_>>()
    });
    lines.into_iter().flatten().map(|l| l + "\n").collect()
}

#[cfg(test)]
mod tests;
