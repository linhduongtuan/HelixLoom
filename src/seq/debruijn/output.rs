//! The result as files: the contigs as FASTA, named as the SPAdes manual names them (`NODE_<number>_length_<bases>_cov_<coverage>`, the
//! coverage being the average k-mer coverage at the last k), and the graph as GFA (segments, and links that overlap by k-1 bases).

use super::assemble::Assembly;
use super::graph::{rc, seg, Graph};

/// A segment of the graph, numbered.
#[derive(Clone, Debug, PartialEq)]
pub struct Contig {
    /// From 1, the longest first.
    pub number: usize,
    /// The segment's place in the graph.
    pub id: usize,
    /// The bases; of a circle, the bases once round (the segment holds the k-1 that repeat as well).
    pub seq: Vec<u8>,
    /// The average number of times a k-mer of the segment was seen.
    pub coverage: f64,
    pub circular: bool,
}

/// Every segment of the graph numbered, the longest first (those of one length by their bases, so that the numbers depend on nothing else).
pub fn number(g: &Graph) -> Vec<Contig> {
    let mut ids = g.ids();
    ids.sort_by(|&a, &b| {
        g.length(b)
            .cmp(&g.length(a))
            .then_with(|| g.segment(a).seq.cmp(&g.segment(b).seq))
    });
    ids.into_iter()
        .enumerate()
        .map(|(i, id)| {
            let circular = g.circular(id);
            let mut seq = g.segment(id).seq.clone();
            if circular {
                seq.truncate(seq.len() - (g.k() - 1));
            }
            Contig {
                number: i + 1,
                id,
                seq,
                coverage: g.coverage(id),
                circular,
            }
        })
        .collect()
}

const WIDTH: usize = 60;

/// The contigs of at least `min_length` bases as FASTA.
pub fn fasta(contigs: &[Contig], min_length: usize) -> String {
    let mut out = String::new();
    for c in contigs.iter().filter(|c| c.seq.len() >= min_length) {
        out += &format!(
            ">NODE_{}_length_{}_cov_{:.6}{}\n",
            c.number,
            c.seq.len(),
            c.coverage,
            if c.circular { " circular" } else { "" }
        );
        for line in c.seq.chunks(WIDTH) {
            out += &String::from_utf8_lossy(line);
            out.push('\n');
        }
    }
    out
}

/// The whole graph as GFA 1.0: a segment for each contig (named by its number, with its length `LN`, its k-mer count `KC` and its depth `DP`)
/// and a link for each pair of ends that meet, overlapping by k-1 bases.
pub fn gfa(g: &Graph, contigs: &[Contig]) -> String {
    let mut number = vec![0usize; g.handles() / 2];
    for c in contigs {
        number[c.id] = c.number;
    }
    let mut out = String::from("H\tVN:Z:1.0\n");
    for c in contigs {
        let s = g.segment(c.id);
        out += &format!(
            "S\t{}\t{}\tLN:i:{}\tKC:i:{}\tDP:f:{:.6}\n",
            c.number,
            String::from_utf8_lossy(&s.seq),
            s.seq.len(),
            s.weight.round() as u64,
            c.coverage
        );
    }
    let orient = |h: usize| if h & 1 == 0 { '+' } else { '-' };
    let written = |a: usize, b: usize| (number[seg(a)], orient(a), number[seg(b)], orient(b));
    // a link is the same link read from the other end: written as the form that starts forward, the lower number first
    let mut links: Vec<(usize, char, usize, char)> = g
        .links()
        .into_iter()
        .map(|(a, b)| {
            let (one, other) = (written(a, b), written(rc(b), rc(a)));
            if (one.1 == '-', one) <= (other.1 == '-', other) {
                one
            } else {
                other
            }
        })
        .collect();
    links.sort_unstable();
    for (a, ao, b, bo) in links {
        out += &format!("L\t{a}\t{ao}\t{b}\t{bo}\t{}M\n", g.k() - 1);
    }
    out
}

/// The shortest of the lengths whose sum with those above it is at least half the total: N50.
pub fn n50(lengths: &[usize]) -> usize {
    let mut sorted = lengths.to_vec();
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let total: usize = sorted.iter().sum();
    let mut sum = 0;
    for l in sorted {
        sum += l;
        if sum * 2 >= total {
            return l;
        }
    }
    0
}

/// What the assembly did, in lines for a person: each k, then the contigs of at least `min_length` bases.
pub fn summary(a: &Assembly, contigs: &[Contig], min_length: usize) -> String {
    let mut out = format!("{} read pieces, {} bases\n", a.pieces, a.bases);
    for s in &a.steps {
        out += &format!(
            "k = {}: {} k-mers, {} segments; {} tips and {} bulges removed; {} under the depth cutoff {:.2} removed; {} short isolated ones removed; {} segments, {} bases\n",
            s.k, s.kmers, s.segments_made, s.tips, s.bulges, s.low_coverage, s.cutoff, s.isolated, s.segments, s.bases
        );
    }
    let kept: Vec<usize> = contigs
        .iter()
        .map(|c| c.seq.len())
        .filter(|&l| l >= min_length)
        .collect();
    out += &format!(
        "{} contigs of {min_length} bases or more: {} bases, the longest {}, N50 {}; {} circular\n",
        kept.len(),
        kept.iter().sum::<usize>(),
        kept.iter().max().copied().unwrap_or(0),
        n50(&kept),
        contigs
            .iter()
            .filter(|c| c.circular && c.seq.len() >= min_length)
            .count()
    );
    out
}

#[cfg(test)]
mod tests;
