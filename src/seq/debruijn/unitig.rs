//! The condensed graph from the k-mers: every maximal unbranched path of k-mers is one segment. A path goes on from a k-mer to the one
//! that follows it only when that is the *only* k-mer that follows, and the one that follows has no other k-mer before it.

use super::graph::{handle, Graph};
use super::pack::{Kmer, Pack, Set, Table};

/// What a walk to the right from a k-mer collected: the bases added, and the sum of the counts of the k-mers passed.
struct Walk {
    bases: Vec<u8>,
    weight: f64,
}

fn present(counts: &Table<f32>, pack: Pack, kmer: Kmer) -> bool {
    counts.contains_key(&pack.canonical(kmer))
}

/// Walk from the k-mer `start` (as read in one direction) while the way on is one k-mer and no other k-mer leads into it, and it has not
/// been passed before.
fn walk(counts: &Table<f32>, pack: Pack, seen: &mut Set, start: Kmer) -> Walk {
    let mut w = Walk {
        bases: Vec::new(),
        weight: 0.0,
    };
    let mut cur = start;
    loop {
        let mut ahead = (0..4u8).filter(|&c| present(counts, pack, pack.append(cur, c)));
        let (Some(c), None) = (ahead.next(), ahead.next()) else {
            break;
        };
        let next = pack.append(cur, c);
        if (0..4u8)
            .filter(|&d| present(counts, pack, pack.prepend(next, d)))
            .count()
            != 1
        {
            break;
        }
        // a k-mer that has been passed ends the walk: on a circle, the first k-mer is where it comes round to
        let canon = pack.canonical(next);
        if !seen.insert(canon) {
            break;
        }
        w.bases.push(c);
        w.weight += f64::from(counts[&canon]);
        cur = next;
    }
    w
}

/// The graph of the k-mers in `counts` (canonical k-mers with the number of times each was seen).
pub fn unitigs(counts: &Table<f32>, pack: Pack) -> Graph {
    let mut graph = Graph::new(pack.k());
    let mut keys: Vec<Kmer> = counts.keys().copied().collect();
    keys.sort_unstable();
    let mut seen = Set::default();
    for first in keys {
        if !seen.insert(first) {
            continue;
        }
        let right = walk(counts, pack, &mut seen, first);
        // on a circle the walk to the right has been all the way round, and the one to the left finds nothing new
        let left = walk(counts, pack, &mut seen, pack.rc(first));
        let weight = f64::from(counts[&first]) + right.weight + left.weight;
        // what was added on the way left, read from the other strand, is the start of the segment
        let mut seq: Vec<u8> = left
            .bases
            .iter()
            .rev()
            .map(|&c| b"ACGT"[usize::from(3 - c)])
            .collect();
        seq.extend(pack.unpack(first));
        seq.extend(right.bases.iter().map(|&c| b"ACGT"[usize::from(c)]));
        graph.add(seq, weight);
    }
    link_ends(&mut graph, pack);
    graph
}

/// The links: the k-mer after the last k-mer of a segment, if it is there, is the first k-mer of the segment that follows.
fn link_ends(graph: &mut Graph, pack: Pack) {
    let k = pack.k();
    let mut starts: Table<usize> = Table::default();
    let mut ends: Vec<(usize, Kmer)> = Vec::new();
    for id in graph.ids() {
        let seq = &graph.segment(id).seq;
        let (head, tail) = (pack.pack(seq), pack.pack(&seq[seq.len() - k..]));
        if let (Some(head), Some(tail)) = (head, tail) {
            starts.insert(head, handle(id, false));
            starts.insert(pack.rc(tail), handle(id, true));
            ends.push((handle(id, false), tail));
            ends.push((handle(id, true), pack.rc(head)));
        }
    }
    for (h, last) in ends {
        for c in 0..4u8 {
            if let Some(&next) = starts.get(&pack.append(last, c)) {
                graph.add_link(h, next);
            }
        }
    }
}

#[cfg(test)]
mod tests;
