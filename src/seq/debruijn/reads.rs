//! Reads for the assembler: FASTA or FASTQ files (not compressed; several joined with commas), cut into the pieces that have only bases
//! that are known and good enough, since a k-mer with a base that is wrong or not known would be an error of the graph.

use super::pack::Pack;
use crate::pipeline::steps::read_fasta;
use crate::seq::fastq::parse_fastq;

#[derive(Clone, Debug, PartialEq)]
pub struct Read {
    pub seq: Vec<u8>,
    /// Phred quality of each base; none for a read from FASTA.
    pub qual: Vec<u8>,
}

pub fn load(spec: &str) -> Result<Vec<Read>, String> {
    let mut out = Vec::new();
    for path in spec.split(',').filter(|p| !p.is_empty()) {
        let text =
            std::fs::read_to_string(path).map_err(|e| format!("Could not read {path}: {e}"))?;
        match text.trim_start().chars().next() {
            Some('>') => out.extend(read_fasta(&text).into_iter().map(|(_, seq)| Read {
                seq,
                qual: Vec::new(),
            })),
            Some('@') => out.extend(
                parse_fastq(&text)
                    .map_err(|e| format!("{path}: {e}"))?
                    .into_iter()
                    .map(|r| Read {
                        seq: r.seq,
                        qual: r.qual,
                    }),
            ),
            _ => {
                return Err(format!(
                    "{path} is neither FASTA nor FASTQ: it should begin with > or @."
                ))
            }
        }
    }
    if out.is_empty() {
        return Err(format!("{spec} holds no reads."));
    }
    Ok(out)
}

/// The runs of each read that hold only A, C, G and T (U as T) with a quality of at least `min_qual`, those of at least `min_length` bases.
pub fn pieces(reads: &[Read], min_qual: u8, min_length: usize) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    for read in reads {
        let mut run: Vec<u8> = Vec::new();
        for (i, &b) in read.seq.iter().enumerate() {
            let good = Pack::code(b).is_some() && read.qual.get(i).is_none_or(|&q| q >= min_qual);
            if good {
                run.push(b.to_ascii_uppercase());
            } else if !run.is_empty() {
                if run.len() >= min_length {
                    out.push(std::mem::take(&mut run));
                }
                run.clear();
            }
        }
        if run.len() >= min_length {
            out.push(run);
        }
    }
    out
}

#[cfg(test)]
mod tests;
