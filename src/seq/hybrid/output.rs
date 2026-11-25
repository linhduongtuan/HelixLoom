//! The result as files: the replicons as FASTA, named by number with their length, depth and whether they are circles, and the graph as GFA
//! (the same writer as `helixloom dbg`'s). The form of the FASTA name is this program's choice (the README that shows Unicycler's could not be
//! read).

use super::rotate::{rotate, Options as RotateOptions};
use crate::seq::debruijn::graph::Graph;
use crate::seq::debruijn::output::{number, Contig};

#[derive(Clone, Debug, PartialEq)]
pub struct Replicon {
    pub number: usize,
    pub id: usize,
    /// The bases; a circle once round, and turned to its start gene when it was found.
    pub seq: Vec<u8>,
    pub depth: f64,
    pub circular: bool,
    pub rotated: bool,
}

/// Every segment of the graph as a replicon, the longest first; circles are turned to start at a gene of `genes`.
pub fn replicons(g: &Graph, genes: &[(String, Vec<u8>)], o: &RotateOptions) -> Vec<Replicon> {
    number(g)
        .into_iter()
        .map(|c: Contig| {
            let turned = if c.circular && !genes.is_empty() {
                rotate(&c.seq, genes, o)
            } else {
                None
            };
            Replicon {
                number: c.number,
                id: c.id,
                rotated: turned.is_some(),
                seq: turned.unwrap_or(c.seq),
                depth: c.coverage,
                circular: c.circular,
            }
        })
        .collect()
}

const WIDTH: usize = 70;

pub fn fasta(replicons: &[Replicon]) -> String {
    let mut out = String::new();
    for r in replicons {
        out += &format!(
            ">{} length={} depth={:.2}x{}\n",
            r.number,
            r.seq.len(),
            r.depth,
            if r.circular { " circular=true" } else { "" }
        );
        for line in r.seq.chunks(WIDTH) {
            out += &String::from_utf8_lossy(line);
            out.push('\n');
        }
    }
    out
}

#[cfg(test)]
mod tests;
