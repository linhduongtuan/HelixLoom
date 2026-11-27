//! Correction of long noisy reads by the reads that overlap them, after NECAT (Chen et al., *Nature Communications* 2021, 12:60). `helixloom
//! correct`.
//!
//! * **choosing the supporting reads** (`select`): overlaps by the program's mapper, a threshold of identity for each read, a limit of
//!   coverage;
//! * **the consensus** (`consensus`): a vote at each position;
//! * **two rounds** (`--rounds`): every round but the last writes each whole read with the positions that enough reads cover corrected, from
//!   the parts of alignments whose blocks have an error of at most half; the last corrects those reads again with the whole alignments and
//!   writes **only the corrected stretches** (a read can give several; one with none is dropped). With one round, it is the first kind.
//!
//! NECAT states no licence, so **nothing of its code or text was used**: this follows the method as the paper describes it, and the parts the
//! paper leaves out (the consensus, the support a position needs, the shortest piece kept, what "difference" its threshold measures) are this
//! program's own, said where they are made. NECAT was not run, so nothing says the result is the same, and the paper's accuracy is not claimed.
//! Made for plasmids and small genomes: every read is held in memory.

pub mod cli;
pub mod consensus;
pub mod select;

use crate::seq::minimap::{Index, Options};
use crate::seq::par::par_map;
use crate::seq::variants::reads::Read;

#[derive(Clone, Debug, PartialEq)]
pub struct Params {
    /// How many rounds; the last writes only the corrected stretches.
    pub rounds: usize,
    /// The number of the best candidates that the threshold is measured on (50).
    pub measured: usize,
    /// The share of a read that an alignment must cover to be measured (0.6).
    pub long_share: f64,
    /// Standard deviations under the mean (5).
    pub sigmas: f64,
    /// The identity that no threshold goes under (0.5).
    pub min_identity: f64,
    /// The candidates that are looked at (200).
    pub max_candidates: usize,
    /// A candidate is passed over when its whole span has more reads than this (12).
    pub coverage: u32,
    /// The size of a block of the first round (500) and the most error of a block that counts (0.5).
    pub block: usize,
    pub max_block_error: f64,
    /// The fewest supporting reads over a position to correct it (this program's choice, 4).
    pub min_support: usize,
    /// The shortest corrected stretch that the last round writes (this program's choice, 500).
    pub min_piece: usize,
    /// The length of the minimizers' k-mers and the window, for the overlaps and the alignments.
    pub k: usize,
    pub w: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            rounds: 2,
            measured: 50,
            long_share: 0.6,
            sigmas: 5.0,
            min_identity: 0.5,
            max_candidates: 200,
            coverage: 12,
            block: 500,
            max_block_error: 0.5,
            min_support: 4,
            min_piece: 500,
            k: 15,
            w: 5,
        }
    }
}

impl Params {
    /// The options of the overlap search: every chain is reported.
    fn overlap(&self) -> Options {
        Options {
            k: self.k,
            w: self.w,
            overlaps: true,
            ..Options::default()
        }
    }

    /// The options of the alignment of a candidate on its template.
    fn mapping(&self) -> Options {
        Options {
            k: self.k,
            w: self.w,
            cigar: true,
            ..Options::default()
        }
    }
}

/// One round: each read corrected by the others. With `whole`, every read comes out whole and the alignments count over their good blocks;
/// otherwise only the corrected stretches come out.
fn round(reads: &[Read], whole: bool, p: &Params) -> Result<Vec<Read>, String> {
    let all = Index::build(
        reads
            .iter()
            .map(|r| (r.name.clone(), r.seq.clone()))
            .collect(),
        &p.overlap(),
    )?;
    let ids: Vec<usize> = (0..reads.len()).collect();
    let out = par_map(&ids, |&i| {
        let read = &reads[i];
        let candidates = select::candidates(&all, i, &read.seq, &p.overlap());
        let own = Index::build(
            vec![("template".to_string(), read.seq.clone())],
            &p.mapping(),
        );
        let supports = match own {
            Ok(own) => select::supports(&own, &candidates, reads, &p.mapping(), p),
            Err(_) => Vec::new(),
        };
        let voted = consensus::vote(&read.seq, &supports, whole, p);
        let named = |name: String, seq: Vec<u8>| Read {
            name,
            seq,
            qual: Vec::new(),
        };
        if whole {
            vec![named(read.name.clone(), voted.whole())]
        } else {
            let pieces = voted.pieces(p.min_piece);
            let several = pieces.len() > 1;
            pieces
                .into_iter()
                .enumerate()
                .map(|(k, seq)| {
                    let name = if several {
                        format!("{}_{}", read.name, k + 1)
                    } else {
                        read.name.clone()
                    };
                    named(name, seq)
                })
                .collect()
        }
    });
    Ok(out
        .into_iter()
        .flatten()
        .filter(|r| !r.seq.is_empty())
        .collect())
}

/// The reads corrected in `p.rounds` rounds. `Err` when there is no read, or an option is not one that can be used.
pub fn correct(reads: &[Read], p: &Params) -> Result<Vec<Read>, String> {
    if reads.is_empty() {
        return Err("There are no reads to correct.".into());
    }
    if p.rounds == 0 {
        return Err("At least one round is needed.".into());
    }
    let mut current = reads.to_vec();
    for r in 0..p.rounds {
        // the last round writes the corrected stretches; a single round keeps the reads whole
        let whole = r + 1 < p.rounds || p.rounds == 1;
        current = round(&current, whole, p)?;
        if current.is_empty() {
            break;
        }
    }
    Ok(current)
}

#[cfg(test)]
mod tests;
