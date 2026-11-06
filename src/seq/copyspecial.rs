//! The original's Edit → Copy Special (`clip_copy_special`, AppMain.tcl line 12354, file offset 0x90351):
//! the conversions it puts on the clipboard, and the FASTA and NCBI Bankit table forms
//! (`bankit_text`, line 12136). The conversions are tested against the original's own code.

use super::dna::{complement, revcom, translate};
use super::doc::Doc;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    /// Only the upper-case letters (lower case marks text the person set apart).
    Uppercase,
    Translated,
    /// Translated after the lower-case letters are dropped.
    TranslatedUppercase,
    Reverse,
    Complement,
    /// T becomes U, t becomes u.
    Rna,
}

/// What Copy Special puts on the clipboard for `text`, reverse-complemented first when asked.
pub fn convert(kind: Kind, text: &[u8], reverse_complement: bool, table: &str) -> String {
    let text = if reverse_complement {
        revcom(text)
    } else {
        text.to_vec()
    };
    let converted: Vec<u8> = match kind {
        Kind::Uppercase => text.into_iter().filter(u8::is_ascii_uppercase).collect(),
        Kind::Translated => translate(&text, table),
        Kind::TranslatedUppercase => {
            let upper: Vec<u8> = text
                .into_iter()
                .filter(|c| !c.is_ascii_lowercase())
                .collect();
            translate(&upper, table)
        }
        Kind::Reverse => text.into_iter().rev().collect(),
        Kind::Complement => text.iter().map(|&c| complement(c)).collect(),
        Kind::Rna => text
            .iter()
            .map(|&c| match c {
                b't' => b'u',
                b'T' => b'U',
                other => other,
            })
            .collect(),
    };
    String::from_utf8_lossy(&converted).into_owned()
}

/// Copy as FASTA: `>name first-last`, then the bases on one line. The numbers are what the original
/// prints: the count of bases before the selection, then the count up to its end.
pub fn fasta(name: &str, before: usize, upto: usize, text: &[u8]) -> String {
    format!(">{name} {before}-{upto}\n{}", String::from_utf8_lossy(text))
}

/// Copy Features as NCBI Bankit table: for each feature a line with its first location and type, a line
/// for every further location, and a `label` line with its name. A feature on the minus strand lists its
/// locations in reverse with each pair written end first; one across the origin starts with the part at the
/// end of the sequence. Features come in the order of the feature list.
pub fn bankit(doc: &Doc) -> String {
    let mut out = String::new();
    for f in &doc.features {
        let mut spans = f.spans.clone();
        spans.sort_unstable();
        if spans.len() > 1 && spans[0].0 == 0 && spans.last().is_some_and(|s| s.1 == doc.seq.len())
        {
            spans.rotate_right(1);
        }
        let mut places: Vec<(usize, usize)> = spans.iter().map(|&(a, b)| (a + 1, b)).collect();
        if f.reverse {
            places.reverse();
        }
        let mut places = places
            .into_iter()
            .map(|(a, b)| if f.reverse { (b, a) } else { (a, b) });
        let Some((a, b)) = places.next() else {
            continue;
        };
        out.push_str(&format!("{a}\t{b}\t{}\n", f.kind));
        for (a, b) in places {
            out.push_str(&format!("{a}\t{b}\n"));
        }
        out.push_str(&format!("\t\t\tlabel\t{}\n", f.name));
    }
    out
}

#[cfg(test)]
mod golden;
#[cfg(test)]
mod tests;
