//! A bridge made real: the best first, each put into the graph as a segment of its own between the two single-copy segments it joins (the
//! bases of the route through the repeats, or the long reads' own bases), and every other link at the two ends taken away, since a single
//! copy has one way on at each end. Each end takes one bridge. The repeats that bridges used up are removed; what is then unbranched is made
//! one segment.

use super::copies::Copies;
use crate::seq::debruijn::graph::{handle, rc, seg, Graph, Handle};
use std::collections::HashSet;

/// What a bridge is made of.
#[derive(Clone, Debug, PartialEq)]
pub enum Via {
    /// The two segments are linked already.
    Link,
    /// A route through repeat segments (as handles).
    Route(Vec<Handle>),
    /// The long reads' own bases between the two segments.
    Reads(Vec<u8>),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Bridge {
    pub from: Handle,
    pub to: Handle,
    pub via: Via,
    pub quality: f64,
    pub reads: usize,
}

/// How to join: the modes of the paper.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// Only bridges of high quality, and only single-copy segments and their bridges are made one.
    Conservative,
    Normal,
    /// Bridges of any quality, and everything that is unbranched is made one, whatever it is.
    Bold,
}

impl Mode {
    /// The lowest quality a bridge needs (the paper's numbers, 25, 10 and 1, on this program's scale).
    pub fn min_quality(self) -> f64 {
        match self {
            Mode::Conservative => 25.0,
            Mode::Normal => 10.0,
            Mode::Bold => 1.0,
        }
    }
}

pub struct Applied {
    pub graph: Graph,
    /// How many bridges were put in.
    pub made: usize,
}

/// The bases of the segment that a bridge makes: the last k-1 of the first segment, the bases between, and the first k-1 of the second. A
/// route's last segment already ends with those k-1 (the links overlap by that much), so only the reads' own bases need them put on.
fn bridge_bases(g: &Graph, b: &Bridge) -> Option<Vec<u8>> {
    let o = g.k() - 1;
    let mut seq = {
        let a = g.oriented(b.from);
        a[a.len() - o..].to_vec()
    };
    match &b.via {
        Via::Link => return None,
        Via::Route(via) => via
            .iter()
            .for_each(|&h| seq.extend_from_slice(&g.oriented(h)[o..])),
        Via::Reads(bases) => {
            seq.extend_from_slice(bases);
            seq.extend_from_slice(&g.oriented(b.to)[..o]);
        }
    }
    Some(seq)
}

/// Put the bridges (best first) into the graph. `copies` is how many copies each segment stands for; a repeat that bridges have used as
/// many times as it has copies is removed. Returns the graph after joining what is unbranched, as the mode allows.
pub fn apply(mut g: Graph, copies: &Copies, bridges: &[Bridge], mode: Mode) -> Applied {
    let mut left = copies.count.clone();
    let mut is_bridge = vec![false; left.len()];
    let mut used: HashSet<Handle> = HashSet::new();
    let mut made = 0;
    for b in bridges {
        // the end of the first segment and the start of the second are two ends that leave: `from` and the reverse of `to`
        let ends = [b.from, rc(b.to)];
        if ends[0] == ends[1] || ends.iter().any(|e| used.contains(e)) {
            continue;
        }
        let between = match &b.via {
            Via::Link => {
                if !g.succ(b.from).contains(&b.to) {
                    continue;
                }
                None
            }
            _ => bridge_bases(&g, b).filter(|s| s.len() >= g.k()),
        };
        let joined = match between {
            Some(seq) => {
                let depth = (g.coverage(seg(b.from)) + g.coverage(seg(b.to))) / 2.0;
                let weight = depth * (seq.len() + 1 - g.k()) as f64;
                let id = g.add(seq, weight);
                left.push(1);
                is_bridge.push(true);
                let h = handle(id, false);
                g.add_link(b.from, h);
                g.add_link(h, b.to);
                if let Via::Route(via) = &b.via {
                    for &v in via {
                        left[seg(v)] = left[seg(v)].saturating_sub(1);
                    }
                }
                h
            }
            None if b.via == Via::Link => b.to,
            None => continue,
        };
        for x in g.succ(b.from).to_vec() {
            if x != joined {
                g.remove_link(b.from, x);
            }
        }
        for p in g.pred(b.to) {
            if p != if joined == b.to { b.from } else { joined } {
                g.remove_link(p, b.to);
            }
        }
        used.extend(ends);
        made += 1;
    }
    // a repeat that every one of its copies has been used in a bridge is not needed
    for id in g.ids() {
        if id < copies.count.len() && copies.count[id] > 1 && left[id] == 0 {
            g.remove(id);
        }
    }
    let single = |id: usize| left[id] == 1 || is_bridge[id];
    let graph = match mode {
        Mode::Bold => g.compact(),
        _ => g.compact_where(&|a, b| single(a) && single(b)),
    };
    Applied { graph, made }
}

#[cfg(test)]
mod tests;
