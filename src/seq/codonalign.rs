//! Aligning DNA by its translation: each coding sequence is translated, the proteins are aligned by whatever aligns them, and the
//! alignment is carried back to the codons, so that every gap is a whole codon and the codons of a column stand for one amino acid.
//! This is what people do by hand, or with a tool such as pal2nal or TranslatorX; this module is written from that description and does
//! no aligning itself.
//!
//! A sequence is read from one reading frame, the same for every sequence: the first `frame` bases are a head that is kept as it is, the
//! codons that follow are translated, and the 0 to 2 bases that are left at the end are a tail. Heads are all as long as the frame, so
//! they line up; tails are put in a block at the end and padded with gaps. A stop codon is written `X` in the protein, so that an
//! aligner that does not know `*` can take it (the codon itself is kept, and comes back); so is a codon with an ambiguous base.

use super::dna::{translate, GENETIC_CODES};

/// How the DNA is read: the bases to leave before the first codon (0 to 2, that is reading frames 1 to 3) and the genetic code (a place
/// in `dna::GENETIC_CODES`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Codons {
    pub frame: usize,
    pub code: usize,
}

impl Codons {
    /// In words, as a record of how an alignment was made says it: `reading frame 2, genetic code 1 (Standard)`.
    pub fn describe(&self) -> String {
        let name = GENETIC_CODES.get(self.code).map_or("unknown", |c| c.name);
        let id = GENETIC_CODES.get(self.code).map_or(0, |c| c.id);
        format!(
            "reading frame {}, genetic code {id} ({name})",
            self.frame + 1
        )
    }

    /// What `describe` wrote, read back; none for anything else.
    pub fn parse(text: &str) -> Option<Codons> {
        let rest = text.strip_prefix("reading frame ")?;
        let (frame, rest) = rest.split_once(", genetic code ")?;
        let (id, _) = rest.split_once(' ')?;
        let frame: usize = frame.parse::<usize>().ok()?.checked_sub(1)?;
        let id: u8 = id.parse().ok()?;
        let code = GENETIC_CODES.iter().position(|c| c.id == id)?;
        (frame <= 2).then_some(Codons { frame, code })
    }
}

/// One sequence cut into its head, its codons and its tail, as the letters it was written with (upper case, no gaps).
#[derive(Clone, Debug, PartialEq)]
struct Cut {
    head: Vec<u8>,
    codons: Vec<[u8; 3]>,
    tail: Vec<u8>,
}

/// The sequences with their proteins, and what is needed to put the codons back.
#[derive(Clone, Debug, PartialEq)]
pub struct Translated {
    /// Name and protein of each sequence, to be aligned.
    pub proteins: Vec<(String, Vec<u8>)>,
    cuts: Vec<Cut>,
}

fn is_gap(b: u8) -> bool {
    matches!(b, b'-' | b'.') || b.is_ascii_whitespace()
}

/// Whether most of the letters are those of DNA or RNA (with ambiguity codes); a protein is not.
fn looks_like_dna(bases: &[u8]) -> bool {
    let dna = bases
        .iter()
        .filter(|b| b"ACGTUNRYKMSWBDHV".contains(b))
        .count();
    !bases.is_empty() && dna * 10 >= bases.len() * 9
}

/// Translate every sequence in the reading frame `c`. `Err` when the choice is not one, a sequence is not DNA, or one has no whole codon.
pub fn translate_all(seqs: &[(String, Vec<u8>)], c: &Codons) -> Result<Translated, String> {
    let table = GENETIC_CODES
        .get(c.code)
        .ok_or("That genetic code is not one this program has.")?
        .table;
    if c.frame > 2 {
        return Err("The reading frame is 1, 2 or 3.".into());
    }
    let mut cuts = Vec::with_capacity(seqs.len());
    let mut proteins = Vec::with_capacity(seqs.len());
    for (name, bases) in seqs {
        let clean: Vec<u8> = bases
            .iter()
            .filter(|&&b| !is_gap(b))
            .map(u8::to_ascii_uppercase)
            .collect();
        if !looks_like_dna(&clean) {
            return Err(format!(
                "{name} does not look like DNA: align it as a protein, or choose DNA sequences."
            ));
        }
        let skip = c.frame.min(clean.len());
        let (head, rest) = clean.split_at(skip);
        let (chunks, tail) = rest.as_chunks::<3>();
        if chunks.is_empty() {
            return Err(format!(
                "{name} has no whole codon in reading frame {}.",
                c.frame + 1
            ));
        }
        let flat: Vec<u8> = chunks.iter().flatten().copied().collect();
        let protein = translate(&flat, table)
            .into_iter()
            .map(|a| if a == b'*' { b'X' } else { a })
            .collect();
        proteins.push((name.clone(), protein));
        cuts.push(Cut {
            head: head.to_vec(),
            codons: chunks.to_vec(),
            tail: tail.to_vec(),
        });
    }
    Ok(Translated { proteins, cuts })
}

/// The DNA alignment that an alignment of the proteins stands for: every residue becomes its codon, every gap three gaps, the heads come
/// first and the tails, padded to one length, last. `Err` when the alignment is not of these proteins (an aligner that changed a
/// sequence, or left one out) or its rows are not one length.
pub fn back_translate(
    t: &Translated,
    aligned: &[(String, Vec<u8>)],
) -> Result<Vec<(String, Vec<u8>)>, String> {
    if aligned.len() != t.cuts.len() {
        return Err(format!(
            "The aligner returned {} sequences for {}.",
            aligned.len(),
            t.cuts.len()
        ));
    }
    let width = aligned.first().map_or(0, |r| r.1.len());
    if aligned.iter().any(|r| r.1.len() != width) {
        return Err("The rows of the protein alignment are not all of one length.".into());
    }
    let tail_width = t.cuts.iter().map(|c| c.tail.len()).max().unwrap_or(0);
    let mut rows = Vec::with_capacity(aligned.len());
    for (((name, protein), cut), (_, row)) in t.proteins.iter().zip(&t.cuts).zip(aligned) {
        let mut out = cut.head.clone();
        let mut next = 0;
        for &residue in row {
            if is_gap(residue) {
                out.extend_from_slice(b"---");
                continue;
            }
            let codon = cut.codons.get(next);
            let same = protein
                .get(next)
                .is_some_and(|p| p.eq_ignore_ascii_case(&residue));
            let Some(codon) = codon.filter(|_| same) else {
                return Err(format!(
                    "The aligner changed the sequence {name}: its residue {} is not what was sent.",
                    next + 1
                ));
            };
            out.extend_from_slice(codon);
            next += 1;
        }
        if next != cut.codons.len() {
            return Err(format!("The aligner left out part of the sequence {name}."));
        }
        out.extend_from_slice(&cut.tail);
        out.resize(cut.head.len() + 3 * width + tail_width, b'-');
        rows.push((name.clone(), out));
    }
    Ok(rows)
}

#[cfg(test)]
mod tests;
