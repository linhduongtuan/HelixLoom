//! The features that reads are counted over: lines of a GTF or GFF3 file of the types that are asked for, grouped by an attribute into
//! meta-features (exons into genes), or the rows of a SAF table (`GeneID Chr Start End Strand`, positions from 1, both ends included).

use crate::seq::gff::{self, Strand};
use std::collections::HashMap;

/// One interval of a sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct Feature {
    pub seqname: String,
    /// From 0, the end not included.
    pub start: usize,
    pub end: usize,
    pub strand: Strand,
    /// The meta-feature that it belongs to (its index in `Annotation::names`).
    pub meta: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Annotation {
    pub features: Vec<Feature>,
    /// The names of the meta-features, in the order each was first met.
    pub names: Vec<String>,
}

#[derive(Default)]
struct Builder {
    features: Vec<Feature>,
    names: Vec<String>,
    index: HashMap<String, usize>,
}

impl Builder {
    fn add(&mut self, name: &str, seqname: &str, start: usize, end: usize, strand: Strand) {
        let meta = match self.index.get(name) {
            Some(&m) => m,
            None => {
                self.names.push(name.to_string());
                self.index.insert(name.to_string(), self.names.len() - 1);
                self.names.len() - 1
            }
        };
        self.features.push(Feature {
            seqname: seqname.to_string(),
            start: start - 1,
            end,
            strand,
            meta,
        });
    }

    fn done(self) -> Annotation {
        Annotation {
            features: self.features,
            names: self.names,
        }
    }
}

/// A SAF table begins with its header, `GeneID` and a tab.
fn is_saf(text: &str) -> bool {
    text.lines()
        .find(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .is_some_and(|l| l.starts_with("GeneID\t"))
}

fn saf(text: &str) -> Result<Annotation, String> {
    let mut b = Builder::default();
    let mut header = true;
    for (k, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        if header {
            header = false;
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            return Err(format!(
                "Line {} of the SAF table has {} columns, and it needs GeneID, Chr, Start, End and Strand.",
                k + 1,
                f.len()
            ));
        }
        let position = |s: &str| {
            s.trim()
                .parse::<usize>()
                .ok()
                .filter(|v| *v >= 1)
                .ok_or_else(|| format!("Line {} of the SAF table has {s} for a position.", k + 1))
        };
        let (start, end) = (position(f[2])?, position(f[3])?);
        if start > end {
            return Err(format!(
                "Line {} of the SAF table starts at {start}, after its end {end}.",
                k + 1
            ));
        }
        let strand = match f[4].trim() {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => Strand::Unknown,
        };
        b.add(f[0], f[1], start, end, strand);
    }
    if b.features.is_empty() {
        return Err("The SAF table has a header and no rows.".into());
    }
    Ok(b.done())
}

/// The features of a GTF, GFF3 or SAF text. `types` are the third column's values that count (`exon`), `attribute` names the meta-feature
/// (`gene_id`); a SAF table has both already, so neither is looked at there.
pub fn parse(text: &str, types: &[String], attribute: &str) -> Result<Annotation, String> {
    if is_saf(text) {
        return saf(text);
    }
    let lines = gff::lines(text)?;
    let mut b = Builder::default();
    for line in lines.iter().filter(|l| types.contains(&l.kind)) {
        let Some((_, name)) = line.attributes.iter().find(|(k, _)| k == attribute) else {
            return Err(format!(
                "The {} at {}:{}-{} has no attribute {attribute}; name the one that says which gene it belongs to with -g.",
                line.kind, line.seqname, line.start, line.end
            ));
        };
        b.add(name, &line.seqname, line.start, line.end, line.strand);
    }
    if b.features.is_empty() {
        let mut found: Vec<&str> = lines.iter().map(|l| l.kind.as_str()).collect();
        found.sort_unstable();
        found.dedup();
        return Err(format!(
            "No line of the annotation is of the type {}; it has {}. Name the type with -t.",
            types.join(" or "),
            found.join(", ")
        ));
    }
    Ok(b.done())
}

#[cfg(test)]
mod tests;
