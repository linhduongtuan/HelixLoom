//! Reads from FASTA or FASTQ files, each with the Phred qualities that a FASTQ gives it (none for FASTA: the pileup then uses a quality that
//! is set once). A read is named by the first word of its header, as SAM names it. Bases come in capitals, and **U is read as T**, so that
//! the reads of RNA are voted and turned round like any others.

use crate::pipeline::steps::read_fasta;
use crate::seq::fastq::parse_fastq;

#[derive(Clone, Debug, PartialEq)]
pub struct Read {
    pub name: String,
    pub seq: Vec<u8>,
    /// Phred scores, one to a base; empty for a read from FASTA.
    pub qual: Vec<u8>,
}

fn first_word(header: &str) -> String {
    header
        .split_whitespace()
        .next()
        .unwrap_or("unnamed")
        .to_string()
}

/// The reads of the files of `spec` (FASTA or FASTQ, not compressed, joined with commas), bases in capitals.
pub fn load(spec: &str) -> Result<Vec<Read>, String> {
    let mut out = Vec::new();
    for path in spec.split(',').filter(|p| !p.is_empty()) {
        let text =
            std::fs::read_to_string(path).map_err(|e| format!("Could not read {path}: {e}"))?;
        out.extend(from_text(&text).map_err(|e| format!("{path}: {e}"))?);
    }
    if out.is_empty() {
        return Err(format!("{spec} holds no sequences."));
    }
    Ok(out)
}

/// A sequence in capitals with every U a T.
fn dna(seq: Vec<u8>) -> Vec<u8> {
    seq.into_iter()
        .map(|b| match b.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        })
        .collect()
}

/// The reads of one text: FASTA when it begins with `>`, FASTQ when it begins with `@`.
pub fn from_text(text: &str) -> Result<Vec<Read>, String> {
    match text.trim_start().chars().next() {
        Some('>') => Ok(read_fasta(text)
            .into_iter()
            .map(|(header, seq)| Read {
                name: first_word(&header),
                seq: dna(seq),
                qual: Vec::new(),
            })
            .collect()),
        Some('@') => Ok(parse_fastq(text)?
            .into_iter()
            .map(|r| Read {
                name: first_word(&r.name),
                seq: dna(r.seq),
                qual: r.qual,
            })
            .collect()),
        _ => Err("it is neither FASTA nor FASTQ: it should begin with > or @.".into()),
    }
}

#[cfg(test)]
mod tests;
