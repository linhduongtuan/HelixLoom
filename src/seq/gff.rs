//! Gene annotation in GFF3 and GTF: genes, transcripts, exons and coding sequences read from a file, the sequences they cut out of a
//! genome, the same written as GFF3 or GTF, and filters over them (`helixloom gff`, after `gffread`'s documented behaviour).
//!
//! Written from the GFF3 specification (Sequence Ontology, version 1.26), the GTF format as the Ensembl and BioPerl pages describe
//! it, and the rules that `gffread`'s usage text and the paper give (`gffread` is MIT licensed; none of its source was copied):
//!
//! * A **transcript** is what exons and CDS segments name as their `Parent` (GFF3) or `transcript_id` (GTF). A CDS with no parent,
//!   as bacterial annotations write them, is a transcript of its own; a record of a transcript kind (`mRNA`, `ncRNA`, `tRNA` …)
//!   with no children is one single exon.
//! * When exons are given, UTR lines add nothing; when CDS is given, start and stop codon lines add nothing; with CDS and UTR only,
//!   the exons are made from them (segments that touch are one), and with CDS only, from the CDS.
//! * Coordinates are 1-based and inclusive. Exons and CDS are kept from the left whichever the strand, and the 5' end of a minus
//!   strand transcript is its right end.
//!
//! Not done, from `gffread`: merging and clustering of transcripts into loci (`-M`, `-K`, `-Q`, `-Y`), BED and TLF, streaming,
//! Ensembl-style GTF conversion (`-L`), and ranges of other sequences than one.

use std::collections::HashMap;

mod extract;
mod features;
mod filter;
mod write;

pub mod cli;
pub use extract::{cds_sequence, protein, spliced, Genome};
pub use features::features_for;
pub use filter::{Filters, Range};
pub use write::{write_gff3, write_gtf};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl Strand {
    fn of(text: &str) -> Strand {
        match text {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => Strand::Unknown,
        }
    }

    pub fn symbol(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
            Strand::Unknown => '.',
        }
    }
}

/// A stretch of coding sequence, with the number of bases to skip from its 5' end to reach a whole codon (GFF3's phase, GTF's frame).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Cds {
    pub start: usize,
    pub end: usize,
    pub phase: u8,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Transcript {
    pub id: String,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub seqname: String,
    pub source: String,
    /// `mRNA`, `transcript`, `ncRNA` … as the file says.
    pub kind: String,
    pub strand: Strand,
    pub start: usize,
    pub end: usize,
    /// From the left, whichever the strand.
    pub exons: Vec<(usize, usize)>,
    pub cds: Vec<Cds>,
    /// The other attributes of the transcript's own line, decoded.
    pub attributes: Vec<(String, String)>,
}

impl Transcript {
    pub fn is_coding(&self) -> bool {
        !self.cds.is_empty()
    }

    /// The lengths of the introns, from the left.
    pub fn introns(&self) -> Vec<(usize, usize)> {
        self.exons
            .windows(2)
            .map(|w| (w[0].1 + 1, w[1].0 - 1))
            .filter(|(a, b)| a <= b)
            .collect()
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Annotation {
    pub transcripts: Vec<Transcript>,
}

impl Annotation {
    /// By sequence in the order each was first met, then by start and end.
    pub fn sort_by_location(&mut self) {
        let mut met: HashMap<String, usize> = HashMap::new();
        for t in &self.transcripts {
            let next = met.len();
            met.entry(t.seqname.clone()).or_insert(next);
        }
        self.transcripts
            .sort_by_key(|t| (met[&t.seqname], t.start, t.end));
    }
}

/// One line of the file.
struct Rec {
    seqname: String,
    source: String,
    feature: String,
    start: usize,
    end: usize,
    strand: Strand,
    phase: Option<u8>,
    attrs: String,
}

fn records(text: &str) -> Result<Vec<Rec>, String> {
    let mut out = Vec::new();
    for (k, line) in text.lines().enumerate() {
        let n = k + 1;
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            return Err(format!(
                "Line {n}: expected nine columns with a tab between them, found {}.",
                f.len()
            ));
        }
        let number = |what: &str, s: &str| -> Result<usize, String> {
            s.trim()
                .parse::<usize>()
                .ok()
                .filter(|v| *v >= 1)
                .ok_or_else(|| {
                    format!("Line {n}: the {what} {s:?} is not a position (they count from 1).")
                })
        };
        let (start, end) = (number("start", f[3])?, number("end", f[4])?);
        if start > end {
            return Err(format!(
                "Line {n}: the start {start} is after the end {end}."
            ));
        }
        out.push(Rec {
            seqname: f[0].to_string(),
            source: f[1].to_string(),
            feature: f[2].to_string(),
            start,
            end,
            strand: Strand::of(f[6]),
            phase: f[7].trim().parse::<u8>().ok().filter(|p| *p <= 2),
            attrs: f[8].to_string(),
        });
    }
    Ok(out)
}

/// What the file is: the directive says it, or the attributes do.
#[derive(PartialEq)]
enum Format {
    Gff3,
    Gtf,
}

fn format_of(text: &str, recs: &[Rec]) -> Result<Format, String> {
    if text.lines().any(|l| l.starts_with("##gff-version 3")) {
        return Ok(Format::Gff3);
    }
    if recs.iter().any(|r| r.attrs.contains("transcript_id")) {
        return Ok(Format::Gtf);
    }
    if recs.iter().any(|r| r.attrs.contains('=')) {
        return Ok(Format::Gff3);
    }
    Err("Could not tell whether this is GFF3 or GTF: no attribute of either kind was found.".into())
}

fn decode(text: &str) -> String {
    let bytes = text.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        let hex = |b: u8| (b as char).to_digit(16);
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            if let (Some(h), Some(l)) = (hex(bytes[i + 1]), hex(bytes[i + 2])) {
                out.push((h * 16 + l) as u8);
                i += 3;
                continue;
            }
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

/// The attributes of a GFF3 line as written (values still escaped).
fn gff3_attributes(attrs: &str) -> Vec<(String, String)> {
    attrs
        .split(';')
        .filter_map(|part| part.split_once('='))
        .map(|(k, v)| (k.trim().to_string(), v.trim().to_string()))
        .collect()
}

/// The attributes of a GTF line: `key "value";` or `key value;`, a `;` inside quotes belonging to the value.
fn gtf_attributes(attrs: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut rest = attrs.trim();
    while !rest.is_empty() {
        let Some((key, after)) = rest.split_once(char::is_whitespace) else {
            break;
        };
        let after = after.trim_start();
        let (value, tail) = if let Some(quoted) = after.strip_prefix('"') {
            match quoted.split_once('"') {
                Some((v, t)) => (v.to_string(), t),
                None => (quoted.to_string(), ""),
            }
        } else {
            match after.split_once(';') {
                Some((v, t)) => (v.trim().to_string(), t),
                None => (after.trim().to_string(), ""),
            }
        };
        out.push((key.to_string(), value));
        rest = tail.trim_start_matches([';', ' ', '\t']);
    }
    out
}

fn value_of(attrs: &[(String, String)], key: &str) -> Option<String> {
    attrs.iter().find(|(k, _)| k == key).map(|(_, v)| decode(v))
}

fn list_of(attrs: &[(String, String)], key: &str) -> Vec<String> {
    attrs
        .iter()
        .find(|(k, _)| k == key)
        .map(|(_, v)| v.split(',').map(decode).collect())
        .unwrap_or_default()
}

const TRANSCRIPT_KINDS: [&str; 12] = [
    "mRNA",
    "transcript",
    "ncRNA",
    "lnc_RNA",
    "tRNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "miRNA",
    "primary_transcript",
    "pseudogenic_transcript",
    "antisense_RNA",
];

fn is_part(feature: &str) -> bool {
    matches!(
        feature,
        "exon"
            | "CDS"
            | "UTR"
            | "five_prime_UTR"
            | "three_prime_UTR"
            | "5UTR"
            | "3UTR"
            | "start_codon"
            | "stop_codon"
    )
}

/// The parts of one transcript, as the records that name it.
struct Parts<'a> {
    own: Option<&'a Rec>,
    parts: Vec<&'a Rec>,
}

/// The transcript made of `parts`, with the rules in the module's text.
fn assemble(
    id: &str,
    p: &Parts,
    gene_id: Option<String>,
    gene_name: Option<String>,
    attributes: Vec<(String, String)>,
) -> Option<Transcript> {
    let first = p.own.or_else(|| p.parts.first().copied())?;
    let mut exons: Vec<(usize, usize)> = Vec::new();
    let mut cds: Vec<Cds> = Vec::new();
    let mut utrs: Vec<(usize, usize)> = Vec::new();
    for r in &p.parts {
        match r.feature.as_str() {
            "exon" => exons.push((r.start, r.end)),
            "CDS" => cds.push(Cds {
                start: r.start,
                end: r.end,
                phase: r.phase.unwrap_or(0),
            }),
            f if f.contains("UTR") => utrs.push((r.start, r.end)),
            _ => {}
        }
    }
    cds.sort_by_key(|c| (c.start, c.end));
    cds.dedup();
    if exons.is_empty() {
        exons = cds.iter().map(|c| (c.start, c.end)).chain(utrs).collect();
    }
    if exons.is_empty() {
        // a record of a transcript kind with nothing under it is one exon
        let own = p.own?;
        if !TRANSCRIPT_KINDS.contains(&own.feature.as_str()) {
            return None;
        }
        exons.push((own.start, own.end));
    }
    exons.sort_unstable();
    exons.dedup();
    // segments that touch or overlap are one
    let mut joined: Vec<(usize, usize)> = Vec::new();
    for (a, b) in exons {
        match joined.last_mut() {
            Some(last) if a <= last.1 + 1 => last.1 = last.1.max(b),
            _ => joined.push((a, b)),
        }
    }
    let kind = match p.own {
        // a gene with CDS under it stands for its transcript
        Some(own) if own.feature == "gene" && cds.is_empty() => "transcript".to_string(),
        Some(own) if own.feature == "gene" => "mRNA".to_string(),
        Some(own) => own.feature.clone(),
        None if cds.is_empty() => "transcript".to_string(),
        None => "mRNA".to_string(),
    };
    Some(Transcript {
        id: id.to_string(),
        gene_id,
        gene_name,
        seqname: first.seqname.clone(),
        source: first.source.clone(),
        kind,
        strand: first.strand,
        start: joined[0].0,
        end: joined[joined.len() - 1].1,
        exons: joined,
        cds,
        attributes,
    })
}

fn from_gff3(recs: &[Rec]) -> Vec<Transcript> {
    let attrs: Vec<Vec<(String, String)>> =
        recs.iter().map(|r| gff3_attributes(&r.attrs)).collect();
    let mut by_id: HashMap<String, usize> = HashMap::new();
    for (i, a) in attrs.iter().enumerate() {
        if let Some(id) = value_of(a, "ID") {
            by_id.entry(id).or_insert(i);
        }
    }
    // who the parts belong to; a CDS with no parent is its own transcript
    let mut kids: HashMap<String, Vec<usize>> = HashMap::new();
    let mut loose: Vec<usize> = Vec::new();
    let mut dangling: Vec<String> = Vec::new();
    for (i, r) in recs.iter().enumerate() {
        if !is_part(&r.feature) {
            continue;
        }
        let parents = list_of(&attrs[i], "Parent");
        if parents.is_empty() {
            if r.feature == "CDS" {
                loose.push(i);
            }
            continue;
        }
        for p in parents {
            if !by_id.contains_key(&p) && !kids.contains_key(&p) {
                dangling.push(p.clone());
            }
            kids.entry(p).or_default().push(i);
        }
    }
    let mut out = Vec::new();
    let mut done: std::collections::HashSet<String> = std::collections::HashSet::new();
    for (i, r) in recs.iter().enumerate() {
        if is_part(&r.feature) {
            if loose.contains(&i) {
                let id =
                    value_of(&attrs[i], "ID").unwrap_or_else(|| format!("cds{}", out.len() + 1));
                let parts = Parts {
                    own: None,
                    parts: vec![r],
                };
                let others = other_attributes(&attrs[i]);
                if let Some(t) = assemble(&id, &parts, None, None, others) {
                    out.push(t);
                }
            }
            continue;
        }
        let Some(id) = value_of(&attrs[i], "ID") else {
            continue;
        };
        if by_id.get(&id) != Some(&i) || done.contains(&id) {
            continue;
        }
        let parts = Parts {
            own: Some(r),
            parts: kids
                .get(&id)
                .map(|k| k.iter().map(|&j| &recs[j]).collect())
                .unwrap_or_default(),
        };
        // the gene it belongs to, and its name
        let (mut gene_id, mut gene_name) = (
            value_of(&attrs[i], "gene_id").or_else(|| value_of(&attrs[i], "geneID")),
            None,
        );
        for parent in list_of(&attrs[i], "Parent") {
            if let Some(&g) = by_id.get(&parent) {
                gene_id = Some(parent.clone());
                gene_name =
                    value_of(&attrs[g], "Name").or_else(|| value_of(&attrs[g], "gene_name"));
            }
        }
        gene_name = gene_name.or_else(|| value_of(&attrs[i], "gene_name"));
        // a gene that has CDS directly under it (bacterial annotations) is the gene of its transcript
        if r.feature == "gene" {
            gene_id = gene_id.or_else(|| Some(id.clone()));
            gene_name = gene_name.or_else(|| value_of(&attrs[i], "Name"));
        }
        // what the CDS lines say of the product and the like, where the gene's own line does not
        let mut attributes = other_attributes(&attrs[i]);
        for part in parts.parts.iter().filter(|p| p.feature == "CDS") {
            for (k, v) in other_attributes(&gff3_attributes(&part.attrs)) {
                if !attributes.iter().any(|(key, _)| *key == k) {
                    attributes.push((k, v));
                }
            }
        }
        if let Some(t) = assemble(&id, &parts, gene_id, gene_name, attributes) {
            done.insert(id);
            out.push(t);
        }
    }
    // parts whose parent has no line of its own
    for id in dangling {
        if done.insert(id.clone()) {
            let parts = Parts {
                own: None,
                parts: kids[&id].iter().map(|&j| &recs[j]).collect(),
            };
            out.extend(assemble(&id, &parts, None, None, Vec::new()));
        }
    }
    out
}

/// The attributes of a GFF3 line that the transcript does not carry as fields.
fn other_attributes(attrs: &[(String, String)]) -> Vec<(String, String)> {
    attrs
        .iter()
        .filter(|(k, _)| {
            !matches!(
                k.as_str(),
                "ID" | "Parent" | "geneID" | "gene_id" | "gene_name"
            )
        })
        .map(|(k, v)| (k.clone(), decode(v)))
        .collect()
}

fn from_gtf(recs: &[Rec]) -> Vec<Transcript> {
    let mut order: Vec<String> = Vec::new();
    let mut groups: HashMap<String, Vec<&Rec>> = HashMap::new();
    let mut gene: HashMap<String, (Option<String>, Option<String>)> = HashMap::new();
    for r in recs {
        let a = gtf_attributes(&r.attrs);
        let Some(id) = a
            .iter()
            .find(|(k, _)| k == "transcript_id")
            .map(|(_, v)| v.clone())
        else {
            continue;
        };
        if !is_part(&r.feature) {
            continue;
        }
        if !groups.contains_key(&id) {
            order.push(id.clone());
        }
        let g = gene.entry(id.clone()).or_default();
        g.0 = g.0.take().or_else(|| {
            a.iter()
                .find(|(k, _)| k == "gene_id")
                .map(|(_, v)| v.clone())
        });
        g.1 = g.1.take().or_else(|| {
            a.iter()
                .find(|(k, _)| k == "gene_name")
                .map(|(_, v)| v.clone())
        });
        groups.entry(id).or_default().push(r);
    }
    order
        .into_iter()
        .filter_map(|id| {
            let parts = Parts {
                own: None,
                parts: groups.remove(&id)?,
            };
            let (gene_id, gene_name) = gene.remove(&id).unwrap_or_default();
            assemble(&id, &parts, gene_id, gene_name, Vec::new())
        })
        .collect()
}

/// The transcripts of a GFF3 or GTF text; why not, when a line cannot be read or the file holds none.
pub fn parse(text: &str) -> Result<Annotation, String> {
    let recs = records(text)?;
    if recs.is_empty() {
        return Err("There is no annotation in this file: no line with nine columns.".into());
    }
    let transcripts = match format_of(text, &recs)? {
        Format::Gff3 => from_gff3(&recs),
        Format::Gtf => from_gtf(&recs),
    };
    if transcripts.is_empty() {
        return Err("This file holds no transcripts: no exon or CDS names a parent.".into());
    }
    Ok(Annotation { transcripts })
}

#[cfg(test)]
mod tests;
