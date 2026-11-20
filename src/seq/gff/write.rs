//! The annotation written as GFF3 or GTF: a line for each transcript, then its exons and CDS. Like `gffread`'s, the output holds the
//! transcripts only (no gene lines).

use super::{Annotation, Transcript};

/// A GFF3 value with what the format reserves written as `%` and two hexadecimal digits.
fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '%' | ';' | '=' | '&' | ',' | '\t' | '\n' | '\r' => {
                out.push_str(&format!("%{:02X}", c as u32))
            }
            c => out.push(c),
        }
    }
    out
}

fn source(t: &Transcript) -> &str {
    if t.source.is_empty() {
        "."
    } else {
        &t.source
    }
}

/// GFF3: the transcript with its `ID`, `geneID` and `gene_name` (and, with `keep_attributes`, the rest of its attributes), then an
/// `exon` line and a `CDS` line for each, naming it as their parent.
pub fn write_gff3(a: &Annotation, keep_attributes: bool) -> String {
    let mut out = String::from("##gff-version 3\n");
    for t in &a.transcripts {
        let mut attrs = vec![format!("ID={}", escape(&t.id))];
        if let Some(g) = &t.gene_id {
            attrs.push(format!("geneID={}", escape(g)));
        }
        if let Some(n) = &t.gene_name {
            attrs.push(format!("gene_name={}", escape(n)));
        }
        if keep_attributes {
            attrs.extend(
                t.attributes
                    .iter()
                    .map(|(k, v)| format!("{k}={}", escape(v))),
            );
        }
        let (src, strand) = (source(t), t.strand.symbol());
        out += &format!(
            "{}\t{src}\t{}\t{}\t{}\t.\t{strand}\t.\t{}\n",
            t.seqname,
            t.kind,
            t.start,
            t.end,
            attrs.join(";")
        );
        for (s, e) in &t.exons {
            out += &format!(
                "{}\t{src}\texon\t{s}\t{e}\t.\t{strand}\t.\tParent={}\n",
                t.seqname,
                escape(&t.id)
            );
        }
        for c in &t.cds {
            out += &format!(
                "{}\t{src}\tCDS\t{}\t{}\t.\t{strand}\t{}\tParent={}\n",
                t.seqname,
                c.start,
                c.end,
                c.phase,
                escape(&t.id)
            );
        }
    }
    out
}

/// GTF: a `transcript` line, then `exon` and `CDS` lines, each with `gene_id` and `transcript_id` (the transcript's own for a gene it
/// has none of) and `gene_name` when there is one; the frame of a CDS is its phase.
pub fn write_gtf(a: &Annotation, keep_attributes: bool) -> String {
    let quote = |v: &str| v.replace('"', "'");
    let mut out = String::new();
    for t in &a.transcripts {
        let gene = t.gene_id.as_deref().unwrap_or(&t.id);
        let mut attrs = format!(
            "gene_id \"{}\"; transcript_id \"{}\";",
            quote(gene),
            quote(&t.id)
        );
        if let Some(n) = &t.gene_name {
            attrs += &format!(" gene_name \"{}\";", quote(n));
        }
        if keep_attributes {
            for (k, v) in &t.attributes {
                attrs += &format!(" {k} \"{}\";", quote(v));
            }
        }
        let (src, strand) = (source(t), t.strand.symbol());
        out += &format!(
            "{}\t{src}\ttranscript\t{}\t{}\t.\t{strand}\t.\t{attrs}\n",
            t.seqname, t.start, t.end
        );
        for (s, e) in &t.exons {
            out += &format!(
                "{}\t{src}\texon\t{s}\t{e}\t.\t{strand}\t.\t{attrs}\n",
                t.seqname
            );
        }
        for c in &t.cds {
            out += &format!(
                "{}\t{src}\tCDS\t{}\t{}\t.\t{strand}\t{}\t{attrs}\n",
                t.seqname, c.start, c.end, c.phase
            );
        }
    }
    out
}
