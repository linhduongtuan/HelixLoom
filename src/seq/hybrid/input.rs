//! The assembly graph as GFA: the segments with their bases and depth, and the links, which must all overlap by the same k-1 bases (a graph
//! from `helixloom dbg` or from SPAdes does). The depth is the `DP` tag, or the k-mer count `KC` over the k-mers of the segment.

use crate::seq::debruijn::graph::{handle, Graph};
use std::collections::HashMap;

pub struct Input {
    pub graph: Graph,
    /// Every segment said how deep it was.
    pub has_depth: bool,
}

struct Seg {
    name: String,
    seq: Vec<u8>,
    depth: Option<f64>,
    kmer_count: Option<f64>,
}

fn tag<'a>(fields: &'a [&str], key: &str) -> Option<&'a str> {
    fields
        .iter()
        .find_map(|f| f.strip_prefix(key).and_then(|v| v.strip_prefix(':')))
        .and_then(|v| v.split_once(':').map(|x| x.1))
}

fn orientation(sign: &str, line: usize) -> Result<bool, String> {
    match sign {
        "+" => Ok(false),
        "-" => Ok(true),
        other => Err(format!(
            "Line {line}: an orientation is + or -, not {other:?}."
        )),
    }
}

/// The overlap of a link as a number of bases: `54M` is 54.
fn overlap(cigar: &str, line: usize) -> Result<usize, String> {
    cigar
        .strip_suffix('M')
        .and_then(|n| n.parse().ok())
        .ok_or_else(|| format!("Line {line}: a link must overlap by a number of matching bases such as 54M, not {cigar:?}."))
}

/// The graph in `text`. `overlap_if_none` is the k-1 to use when the file has no link to say it.
pub fn read(text: &str, overlap_if_none: usize) -> Result<Input, String> {
    let mut segs: Vec<Seg> = Vec::new();
    let mut by_name: HashMap<String, usize> = HashMap::new();
    let mut links: Vec<(usize, usize, bool, usize, bool, usize)> = Vec::new();
    for (n, line) in text.lines().enumerate().map(|(n, l)| (n + 1, l.trim_end())) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "S" => {
                if f.len() < 3 {
                    return Err(format!("Line {n}: a segment needs a name and its bases."));
                }
                if f[2] == "*" {
                    return Err(format!(
                        "Line {n}: segment {} has no bases (*): the bases are needed to join it to others.",
                        f[1]
                    ));
                }
                if by_name.insert(f[1].to_string(), segs.len()).is_some() {
                    return Err(format!("Line {n}: segment {} is there twice.", f[1]));
                }
                segs.push(Seg {
                    name: f[1].to_string(),
                    seq: f[2].to_ascii_uppercase().into_bytes(),
                    depth: tag(&f[3..], "DP").and_then(|v| v.parse().ok()),
                    kmer_count: tag(&f[3..], "KC").and_then(|v| v.parse().ok()),
                });
            }
            "L" => {
                if f.len() < 6 {
                    return Err(format!(
                        "Line {n}: a link needs two segments with their orientations and an overlap."
                    ));
                }
                let (a, b) = (by_name.get(f[1]), by_name.get(f[3]));
                let (Some(&a), Some(&b)) = (a, b) else {
                    return Err(format!(
                        "Line {n}: the link {} to {} names a segment that is not (yet) in the file; segments must come before their links.",
                        f[1], f[3]
                    ));
                };
                links.push((
                    a,
                    b,
                    orientation(f[2], n)?,
                    overlap(f[5], n)?,
                    orientation(f[4], n)?,
                    n,
                ));
            }
            _ => {}
        }
    }
    if segs.is_empty() {
        return Err("The graph has no segment.".into());
    }
    let ov = links.first().map_or(overlap_if_none, |l| l.3);
    if let Some(bad) = links.iter().find(|l| l.3 != ov) {
        return Err(format!(
            "Line {}: this link overlaps by {} bases but the first overlaps by {ov}: the links of a de Bruijn graph all overlap by k-1.",
            bad.5, bad.3
        ));
    }
    let k = ov + 1;
    let mut graph = Graph::new(k);
    for s in &segs {
        if s.seq.len() < k {
            return Err(format!(
                "Segment {} has {} bases, fewer than the {k} that the overlap of {ov} implies.",
                s.name,
                s.seq.len()
            ));
        }
    }
    let mut has_depth = true;
    for s in &segs {
        let kmers = (s.seq.len() + 1 - k) as f64;
        let depth = s
            .depth
            .or_else(|| s.kmer_count.map(|c| c / kmers))
            .unwrap_or_else(|| {
                has_depth = false;
                1.0
            });
        graph.add(s.seq.clone(), depth * kmers);
    }
    for &(a, b, ra, _, rb, n) in &links {
        let (ha, hb) = (handle(a, ra), handle(b, rb));
        let (oa, ob) = (graph.oriented(ha), graph.oriented(hb));
        if oa[oa.len() - ov..] != ob[..ov] {
            return Err(format!(
                "Line {n}: the last {ov} bases of {} do not match the first {ov} of {}.",
                segs[a].name, segs[b].name
            ));
        }
        graph.add_link(ha, hb);
    }
    Ok(Input { graph, has_depth })
}

#[cfg(test)]
mod tests;
