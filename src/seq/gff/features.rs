//! An annotation as features of a sequence that is open: each coding transcript a `CDS` (and an `mRNA` when its exons say more than
//! its CDS does), each other transcript a feature of its own kind, on the strand it is on, with the gene's name where there is one.

use super::{Annotation, Strand, Transcript};
use crate::seq::features::Feature;

/// What was made, how many were left out because they reach beyond the sequence, and the name the annotation has for the sequence.
#[derive(Debug, PartialEq)]
pub struct Imported {
    pub features: Vec<Feature>,
    pub skipped: usize,
    pub seqname: String,
}

fn spans(ranges: impl Iterator<Item = (usize, usize)>) -> Vec<(usize, usize)> {
    ranges.map(|(a, b)| (a - 1, b)).collect()
}

fn feature(t: &Transcript, kind: &str, spans: Vec<(usize, usize)>) -> Feature {
    let name = t
        .gene_name
        .clone()
        .or_else(|| {
            t.attributes
                .iter()
                .find(|(k, _)| k == "Name")
                .map(|(_, v)| v.clone())
        })
        .unwrap_or_else(|| t.id.clone());
    let mut f = Feature::new(
        name,
        kind,
        spans[0].0,
        spans[0].1,
        t.strand == Strand::Minus,
    );
    f.spans = spans;
    if let Some(g) = &t.gene_name {
        f.qualifiers.push(("gene".into(), g.clone()));
    }
    if let Some((_, product)) = t.attributes.iter().find(|(k, _)| k == "product") {
        f.qualifiers.push(("product".into(), product.clone()));
    }
    f.qualifiers
        .push(("note".into(), format!("From the annotation: {}", t.id)));
    f
}

/// The kind of feature that a transcript of this kind is when it has no CDS.
fn rna_kind(kind: &str) -> &'static str {
    match kind {
        "ncRNA" => "ncRNA",
        "tRNA" => "tRNA",
        "rRNA" => "rRNA",
        "mRNA" | "transcript" => "mRNA",
        _ => "misc_RNA",
    }
}

/// The features of the transcripts on the sequence called `doc_name` (the annotation's name for it may be the first word of it); when
/// none has that name and the annotation names only one sequence, that one. `length` is the length of the sequence.
pub fn features_for(a: &Annotation, doc_name: &str, length: usize) -> Result<Imported, String> {
    let first_word = doc_name.split_whitespace().next().unwrap_or(doc_name);
    let mut names: Vec<&str> = Vec::new();
    for t in &a.transcripts {
        if !names.contains(&t.seqname.as_str()) {
            names.push(&t.seqname);
        }
    }
    let seqname = match names.iter().find(|n| **n == doc_name || **n == first_word) {
        Some(n) => (*n).to_string(),
        None if names.len() == 1 => names[0].to_string(),
        None => {
            return Err(format!(
                "None of the sequences of the annotation ({}) is called {doc_name}.",
                names.join(", ")
            ))
        }
    };
    let mut features = Vec::new();
    let mut skipped = 0;
    for t in a.transcripts.iter().filter(|t| t.seqname == seqname) {
        let made: Vec<Feature> = if t.is_coding() {
            let cds = spans(t.cds.iter().map(|c| (c.start, c.end)));
            let exons = spans(t.exons.iter().copied());
            let mut v = vec![feature(t, "CDS", cds.clone())];
            if exons != cds {
                v.push(feature(t, "mRNA", exons));
            }
            v
        } else {
            vec![feature(
                t,
                rna_kind(&t.kind),
                spans(t.exons.iter().copied()),
            )]
        };
        if made
            .iter()
            .any(|f| f.spans.iter().any(|(_, end)| *end > length))
        {
            skipped += 1;
        } else {
            features.extend(made);
        }
    }
    features.sort_by_key(|f| (f.start(), f.name.clone()));
    Ok(Imported {
        features,
        skipped,
        seqname,
    })
}
