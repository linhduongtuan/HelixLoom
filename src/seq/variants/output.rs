//! What the variant caller writes: a VCF (version 4.3) of the differences with one haploid sample, the consensus as FASTA, the depth at each
//! position as a table, and a summary line for each target.

use super::call::{Calls, Kind, Params};

/// How many bases to a line of FASTA.
const FASTA_WIDTH: usize = 60;

/// One target with its calls.
pub struct Target<'a> {
    pub name: &'a str,
    pub length: usize,
    pub calls: &'a Calls,
}

/// The VCF: `##` header lines, the columns, and a line for each difference (the position of a line counts from 1; the sample's genotype is
/// `1`, the one copy it has holds the alternate). QUAL has one decimal; FILTER is `PASS` (a call under the lowest quality is not made).
pub fn vcf(results: &[Target<'_>], sample: &str, p: &Params) -> String {
    let mut out = String::from("##fileformat=VCFv4.3\n");
    out.push_str(&format!(
        "##source=helixloom {}\n",
        env!("CARGO_PKG_VERSION")
    ));
    for r in results {
        out.push_str(&format!("##contig=<ID={},length={}>\n", r.name, r.length));
    }
    out.push_str(&format!(
        "##variantCaller=haploid; lowest quality {}; lowest depth {}\n",
        p.min_quality, p.min_depth
    ));
    out.push_str("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Reads at the site\">\n");
    out.push_str("##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"The difference is an insertion or a deletion\">\n");
    out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    out.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Reads at the site\">\n");
    out.push_str("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Reads with the reference and with the alternate\">\n");
    out.push_str(&format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    ));
    for r in results {
        for v in &r.calls.variants {
            let info = if v.kind == Kind::Snp {
                format!("DP={}", v.depth)
            } else {
                format!("DP={};INDEL", v.depth)
            };
            out.push_str(&format!(
                "{}\t{}\t.\t{}\t{}\t{:.1}\tPASS\t{info}\tGT:DP:AD\t1:{}:{},{}\n",
                r.name,
                v.pos + 1,
                String::from_utf8_lossy(&v.reference),
                String::from_utf8_lossy(&v.alt),
                v.quality,
                v.depth,
                v.ref_reads,
                v.alt_reads,
            ));
        }
    }
    out
}

/// The consensus of each target, named `NAME consensus`, in lines of 60.
pub fn consensus(results: &[Target<'_>]) -> String {
    let mut out = String::new();
    for r in results {
        out.push_str(&format!(">{} consensus\n", r.name));
        for line in r.calls.consensus.chunks(FASTA_WIDTH) {
            out.push_str(&String::from_utf8_lossy(line));
            out.push('\n');
        }
    }
    out
}

/// Target, position (from 1), the reads over it.
pub fn depth_table(results: &[Target<'_>]) -> String {
    let mut out = String::from("target\tposition\tdepth\n");
    for r in results {
        for (i, d) in r.calls.depth.iter().enumerate() {
            out.push_str(&format!("{}\t{}\t{d}\n", r.name, i + 1));
        }
    }
    out
}

/// One line to a target: its length, the mean depth, the positions written as N, and the differences by kind.
pub fn summary(results: &[Target<'_>]) -> String {
    let mut out = String::new();
    for r in results {
        let count = |k: Kind| r.calls.variants.iter().filter(|v| v.kind == k).count();
        out.push_str(&format!(
            "{}: {} bases, mean depth {:.1}, {} written as N, {} substitutions, {} insertions, {} deletions\n",
            r.name,
            r.length,
            r.calls.mean_depth,
            r.calls.unknown,
            count(Kind::Snp),
            count(Kind::Insertion),
            count(Kind::Deletion),
        ));
    }
    out
}

#[cfg(test)]
mod tests;
