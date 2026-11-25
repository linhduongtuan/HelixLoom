//! The ways through the repeats of the graph from the end of one single-copy segment to the start of another: a route is the repeat
//! segments passed, and the distance it makes between the two (the bases between the end of the first and the start of the second, so −(k−1)
//! when they are linked directly and overlap).

use crate::seq::debruijn::graph::{seg, Graph, Handle};

/// The most segments a route may pass, and the most steps of the search: a graph of repeats can branch without end.
const MAX_VIA: usize = 8;
const MAX_STEPS: usize = 20_000;

#[derive(Clone, Debug, PartialEq)]
pub struct Route {
    /// The repeat segments passed, as handles, in order (none for a direct link).
    pub via: Vec<Handle>,
    pub distance: i64,
}

/// How far the distance of a route may be from what the reads say, in bases: a fixed part and a share of the gap (long reads have indels).
pub fn tolerance(gap: i64) -> i64 {
    100 + gap.abs() / 5
}

/// The routes from `from` (leaving through its end) to `to` (entered at its start) that pass only segments that are not single copy and whose
/// distance is within `tolerance(gap)` of `gap`, the closest first (a shorter route first among equals).
pub fn routes(g: &Graph, count: &[u32], from: Handle, to: Handle, gap: i64) -> Vec<Route> {
    let overlap = g.k() as i64 - 1;
    let tol = tolerance(gap);
    let mut found: Vec<Route> = Vec::new();
    let mut steps = 0;
    let mut stack: Vec<(Handle, Vec<Handle>, i64)> = vec![(from, Vec::new(), 0)];
    while let Some((at, via, distance)) = stack.pop() {
        steps += 1;
        if steps > MAX_STEPS {
            break;
        }
        for &x in g.succ(at) {
            if x == to {
                let total = distance - overlap;
                if (total - gap).abs() <= tol {
                    found.push(Route {
                        via: via.clone(),
                        distance: total,
                    });
                }
                continue;
            }
            let repeat = count[seg(x)] != 1 && seg(x) != seg(from) && seg(x) != seg(to);
            let next = distance + g.length(seg(x)) as i64 - overlap;
            // what is left to add can only take bases off (the last overlap), so a route already too long is given up
            if repeat && via.len() < MAX_VIA && next - overlap <= gap + tol {
                let mut longer = via.clone();
                longer.push(x);
                stack.push((x, longer, next));
            }
        }
    }
    found.sort_by_key(|r| ((r.distance - gap).abs(), r.via.len(), r.via.clone()));
    found.dedup();
    found
}

#[cfg(test)]
mod tests;
