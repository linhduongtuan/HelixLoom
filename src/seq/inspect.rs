//! What the inspector panel says about the thing that is selected: a feature, or else a range of bases.
//!
//! No UI code here. The panel asks `inspect` once per frame and draws the answer, so the numbers can be checked without a window.

use super::dna::{gc_percent, melting_temp, orf_summary, translate};
use super::doc::Doc;
use super::features::Feature;

/// Residues of a translation the panel shows before it says "…" (a plasmid's biggest gene is a few hundred).
pub const MOST_PROTEIN: usize = 60;

/// What is being looked at.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Subject {
    /// The feature with this index in `Doc::features`.
    Feature(usize),
    /// The bases from `a` up to, not including, `b` (0-based).
    Range(usize, usize),
}

/// The start of a feature's translation, on the strand the feature is on.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Protein {
    /// The first `MOST_PROTEIN` residues, one letter each (`*` is a stop).
    pub shown: String,
    /// How many residues there are in all.
    pub residues: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Inspection {
    pub subject: Subject,
    /// The feature's name, or "Selection".
    pub title: String,
    /// The feature's type (`CDS`, `promoter`…), or "Range" for bases without a feature.
    pub kind: String,
    /// 1-based, as GenBank writes it: `10..20`, `complement(10..20)`, `join(…)`.
    pub location: String,
    pub length: usize,
    /// "forward" or "reverse" for a feature, empty for a range.
    pub strand: &'static str,
    pub gc: f64,
    /// `None` when the original's formula has no answer: nothing, or more than 100 bases.
    pub tm: Option<f64>,
    /// The original's selection summary (`M>85*/<85`), for the bases as they read on the feature's strand.
    pub orf: String,
    /// Only for a CDS.
    pub protein: Option<Protein>,
    pub qualifiers: Vec<(String, String)>,
    /// Indices of the other features that share a base with this one, in order of where they start.
    pub overlaps: Vec<usize>,
}

impl Inspection {
    /// The label and value of each line of the panel's table, in order. A line with nothing to say is left out.
    pub fn rows(&self) -> Vec<(&'static str, String)> {
        let length = match &self.protein {
            Some(p) => format!("{} bp · {} codons", self.length, p.residues),
            None => format!("{} bp", self.length),
        };
        let mut rows = vec![
            ("Type", self.kind.clone()),
            ("Location", self.location.clone()),
            ("Length", length),
        ];
        if !self.strand.is_empty() {
            rows.push(("Strand", self.strand.to_string()));
        }
        rows.push(("GC", format!("{:.0}%", self.gc)));
        if let Some(tm) = self.tm {
            rows.push(("Tm", format!("{tm:.1}° C")));
        }
        if !self.orf.is_empty() {
            rows.push(("ORF", self.orf.clone()));
        }
        rows
    }
}

/// The inspector's answer for the selected feature if there is one, else for the selected range, else `None`.
/// `table` is the genetic code in use (64 letters, `dna::GENETIC_CODES`).
pub fn inspect(
    doc: &Doc,
    feature: Option<usize>,
    range: Option<(usize, usize)>,
    dna_nm: f64,
    salt_mm: f64,
    table: &str,
) -> Option<Inspection> {
    if let Some(fi) = feature.filter(|&fi| fi < doc.features.len()) {
        return Some(of_feature(doc, fi, dna_nm, salt_mm, table));
    }
    let (a, b) = range?;
    let b = b.min(doc.len());
    (a < b).then(|| of_range(doc, a, b, dna_nm, salt_mm))
}

fn of_feature(doc: &Doc, fi: usize, dna_nm: f64, salt_mm: f64, table: &str) -> Inspection {
    let f = &doc.features[fi];
    // everything below reads the bases the way the feature does: its pieces in reading order (a feature over the origin of a
    // circle starts with its tail), and the other strand reverse-complemented
    let bases = doc.feature_bases(f);
    let protein = (f.kind == "CDS").then(|| {
        let aa = translate(&bases, table);
        Protein {
            shown: String::from_utf8_lossy(&aa[..aa.len().min(MOST_PROTEIN)]).into_owned(),
            residues: aa.len(),
        }
    });
    let mut overlaps: Vec<usize> = (0..doc.features.len())
        .filter(|&oi| oi != fi && shares_a_base(f, &doc.features[oi]))
        .collect();
    overlaps.sort_by_key(|&oi| (doc.features[oi].start(), oi));
    Inspection {
        subject: Subject::Feature(fi),
        title: f.name.clone(),
        kind: f.kind.clone(),
        location: f.location_text(),
        length: bases.len(),
        strand: if f.reverse { "reverse" } else { "forward" },
        gc: gc_percent(&bases),
        tm: melting_temp(&bases, dna_nm, salt_mm),
        orf: orf_summary(&bases),
        protein,
        qualifiers: f.qualifiers.clone(),
        overlaps,
    }
}

fn of_range(doc: &Doc, a: usize, b: usize, dna_nm: f64, salt_mm: f64) -> Inspection {
    let bases = doc.slice(a, b);
    let mut overlaps: Vec<usize> = (0..doc.features.len())
        .filter(|&oi| doc.features[oi].spans.iter().any(|&(s, e)| s < b && a < e))
        .collect();
    overlaps.sort_by_key(|&oi| (doc.features[oi].start(), oi));
    Inspection {
        subject: Subject::Range(a, b),
        title: "Selection".to_string(),
        kind: "Range".to_string(),
        location: format!("{}..{}", a + 1, b),
        length: bases.len(),
        strand: "",
        gc: gc_percent(bases),
        tm: melting_temp(bases, dna_nm, salt_mm),
        orf: orf_summary(bases),
        protein: None,
        qualifiers: Vec::new(),
        overlaps,
    }
}

fn shares_a_base(a: &Feature, b: &Feature) -> bool {
    a.spans
        .iter()
        .any(|&(s1, e1)| b.spans.iter().any(|&(s2, e2)| s1 < e2 && s2 < e1))
}

#[cfg(test)]
mod tests;
