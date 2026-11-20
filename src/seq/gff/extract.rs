//! The sequences that a transcript cuts out of a genome: spliced exons, the coding sequence (with the phase of its first segment
//! taken off) and the protein it translates to. Everything is cut in capitals.

use super::{Strand, Transcript};
use crate::seq::dna::{revcom, translate, STANDARD_CODE};
use std::collections::HashMap;

/// The sequences of a genome, found by the first word of the name of each.
#[derive(Clone, Debug, Default)]
pub struct Genome {
    sequences: HashMap<String, Vec<u8>>,
}

impl Genome {
    pub fn from_fasta(text: &str) -> Genome {
        let mut sequences: HashMap<String, Vec<u8>> = HashMap::new();
        let mut current: Option<String> = None;
        for line in text.lines() {
            if let Some(header) = line.strip_prefix('>') {
                let name = header.split_whitespace().next().unwrap_or("").to_string();
                sequences.entry(name.clone()).or_default();
                current = Some(name);
            } else if let Some(name) = &current {
                if let Some(bases) = sequences.get_mut(name) {
                    bases.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
                }
            }
        }
        Genome { sequences }
    }

    pub fn sequence(&self, name: &str) -> Option<&[u8]> {
        self.sequences.get(name).map(Vec::as_slice)
    }
}

/// Bases `start..=end` (counted from 1) of sequence `name`.
pub(super) fn cut<'a>(
    genome: &'a Genome,
    name: &str,
    start: usize,
    end: usize,
) -> Result<&'a [u8], String> {
    let bases = genome
        .sequence(name)
        .ok_or_else(|| format!("The sequence {name} is not in the genome FASTA."))?;
    if end > bases.len() {
        return Err(format!(
            "A feature reaches base {end}, beyond the end of {name} ({} bases).",
            bases.len()
        ));
    }
    Ok(&bases[start - 1..end])
}

fn on_strand(mut bases: Vec<u8>, strand: Strand) -> Vec<u8> {
    if strand == Strand::Minus {
        bases = revcom(&bases);
    }
    bases
}

/// The exons joined from the 5' end, on the strand of the transcript.
pub fn spliced(t: &Transcript, genome: &Genome) -> Result<Vec<u8>, String> {
    let mut out = Vec::new();
    for &(a, b) in &t.exons {
        out.extend_from_slice(cut(genome, &t.seqname, a, b)?);
    }
    Ok(on_strand(out, t.strand))
}

/// The coding sequence: the CDS segments joined from the 5' end, less the bases that the phase of the first of them says; none for a
/// transcript without CDS.
pub fn cds_sequence(t: &Transcript, genome: &Genome) -> Result<Option<Vec<u8>>, String> {
    if t.cds.is_empty() {
        return Ok(None);
    }
    let mut out = Vec::new();
    for c in &t.cds {
        out.extend_from_slice(cut(genome, &t.seqname, c.start, c.end)?);
    }
    let out = on_strand(out, t.strand);
    // the 5' end of a minus strand transcript is its right end
    let first = if t.strand == Strand::Minus {
        t.cds.last()
    } else {
        t.cds.first()
    };
    let skip = first.map_or(0, |c| usize::from(c.phase)).min(out.len());
    Ok(Some(out[skip..].to_vec()))
}

/// The protein of the coding sequence by the standard code. A stop codon at the end is not in it (as a GFF3 CDS includes it); one
/// inside is written as `stop`. A last codon that is not whole is left out.
pub fn protein(t: &Transcript, genome: &Genome, stop: u8) -> Result<Option<Vec<u8>>, String> {
    let Some(cds) = cds_sequence(t, genome)? else {
        return Ok(None);
    };
    let mut aa = translate(&cds, STANDARD_CODE);
    if aa.last() == Some(&b'*') {
        aa.pop();
    }
    Ok(Some(
        aa.into_iter()
            .map(|c| if c == b'*' { stop } else { c })
            .collect(),
    ))
}
