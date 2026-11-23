//! The removal of what is wrong with the graph of reads: tips, bulges, and what is seen too few times.
//!
//! From the papers: a **tip** is a short stray run that leaves a junction and goes nowhere (an error near the end of a read makes one);
//! tips are taken from the shortest, and a tip is kept only when its coverage is far above that of what it is alongside. A **bulge** is two
//! short paths between the same two places that spell nearly the same bases (an error in the middle of a read makes one, so does a base
//! at which two copies differ); the path with less coverage is removed, and what it was seen is added to the other. **Low coverage**: what
//! is seen fewer times than a cutoff goes. What the papers leave open (the numbers of the bulge test, the automatic cutoff) is this
//! program's choice, and `docs/DEBRUIJN.md` says which. Lengths of tips and bulges are counted in k-mers, as the length of an edge is in
//! the papers: a segment of n k-mers has n + k - 1 bases, so that an error makes the same tip at every k.

use super::graph::{handle, seg, Graph, Handle};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Settings {
    /// The longest a tip may be, in k-mers.
    pub tip_length: usize,
    /// A tip is removed unless its coverage is more than this many times the coverage of the path that it is alongside.
    pub tip_ratio: f64,
    /// The longest a path may be to be removed as a bulge, in k-mers.
    pub bulge_length: usize,
    /// The share of its bases (the longer path's) by which the two paths of a bulge may differ (edits).
    pub bulge_diff: f64,
}

/// Is the segment a tip: dead at one end, joined at the other to places that have a way on without it.
fn is_tip(g: &Graph, id: usize, s: &Settings) -> bool {
    [handle(id, false), handle(id, true)].into_iter().any(|h| {
        let preds = g.pred(h);
        if !g.succ(h).is_empty() || preds.is_empty() {
            return false;
        }
        if preds.iter().any(|&p| seg(p) == id || g.outdeg(p) < 2) {
            return false;
        }
        let alongside = preds
            .iter()
            .flat_map(|&p| g.succ(p).iter().copied())
            .filter(|&x| seg(x) != id)
            .map(|x| g.coverage(seg(x)))
            .fold(0.0, f64::max);
        g.coverage(id) <= s.tip_ratio * alongside
    })
}

/// Remove the tips, the shortest first; how many were removed.
pub fn remove_tips(g: &mut Graph, s: &Settings) -> usize {
    let mut ids = g.ids();
    ids.sort_by_key(|&id| (g.kmers(id), id));
    let mut removed = 0;
    for id in ids {
        if g.alive(id) && g.kmers(id) <= s.tip_length && is_tip(g, id, s) {
            g.remove(id);
            removed += 1;
        }
    }
    removed
}

/// The fewest changes (a base replaced, put in or taken out) that turn `a` into `b`.
fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let mut row: Vec<usize> = (0..=b.len()).collect();
    for (i, &x) in a.iter().enumerate() {
        let mut diagonal = row[0];
        row[0] = i + 1;
        for (j, &y) in b.iter().enumerate() {
            let above = row[j + 1];
            row[j + 1] = (diagonal + usize::from(x != y))
                .min(row[j] + 1)
                .min(above + 1);
            diagonal = above;
        }
    }
    row[b.len()]
}

fn similar(a: &[u8], b: &[u8], diff: f64) -> bool {
    let allowed = (diff * a.len().max(b.len()) as f64).floor() as usize;
    a.len().abs_diff(b.len()) <= allowed && edit_distance(a, b) <= allowed
}

/// The paths of one step that lead from the end of `a` to the same handle, grouped by where they lead.
fn arms(g: &Graph, a: Handle) -> Vec<(Handle, Vec<Handle>)> {
    let mut groups: Vec<(Handle, Vec<Handle>)> = Vec::new();
    for &x in g.succ(a) {
        if seg(x) == seg(a) || g.indeg(x) != 1 || g.outdeg(x) != 1 {
            continue;
        }
        let b = g.succ(x)[0];
        // a path that leads back to the segment it left is a path of a circle: its other path may be the same bulge
        if seg(b) == seg(x) {
            continue;
        }
        match groups.iter_mut().find(|(to, _)| *to == b) {
            Some((_, arms)) => arms.push(x),
            None => groups.push((b, vec![x])),
        }
    }
    groups
}

/// Remove the bulges: of the paths between two places, a path that is short and spells nearly what a path with more coverage spells goes,
/// and that path gains what it was seen. How many paths were removed.
pub fn pop_bulges(g: &mut Graph, s: &Settings) -> usize {
    let mut removed = 0;
    for a in 0..g.handles() {
        if !g.alive(seg(a)) || g.outdeg(a) < 2 {
            continue;
        }
        for (_, mut group) in arms(g, a) {
            group.sort_by(|&x, &y| {
                g.coverage(seg(y))
                    .total_cmp(&g.coverage(seg(x)))
                    .then(seg(x).cmp(&seg(y)))
            });
            for (i, &x) in group.iter().enumerate().skip(1) {
                if !g.alive(seg(x)) || g.kmers(seg(x)) > s.bulge_length {
                    continue;
                }
                let like = group[..i].iter().copied().find(|&y| {
                    g.alive(seg(y))
                        && seg(y) != seg(x)
                        && similar(&g.oriented(y), &g.oriented(x), s.bulge_diff)
                });
                if let Some(y) = like {
                    let gain = g.coverage(seg(x)) * g.kmers(seg(y)) as f64;
                    g.add_weight(seg(y), gain);
                    g.remove(seg(x));
                    removed += 1;
                }
            }
        }
    }
    removed
}

/// Remove the segments whose coverage is under `cutoff`; how many.
pub fn remove_below(g: &mut Graph, cutoff: f64) -> usize {
    let low: Vec<usize> = g
        .ids()
        .into_iter()
        .filter(|&id| g.coverage(id) < cutoff)
        .collect();
    for &id in &low {
        g.remove(id);
    }
    low.len()
}

/// Remove the segments that are joined to nothing and are shorter than `max_length` bases: what is left of an error that had no company, or
/// of a stretch that no read joins to the rest. (The paper removes such short isolated paths last of all.) How many.
pub fn remove_isolated(g: &mut Graph, max_length: usize) -> usize {
    let alone: Vec<usize> = g
        .ids()
        .into_iter()
        .filter(|&id| {
            g.length(id) < max_length
                && g.succ(handle(id, false)).is_empty()
                && g.succ(handle(id, true)).is_empty()
        })
        .collect();
    for &id in &alone {
        g.remove(id);
    }
    alone.len()
}

/// The cutoff for "automatic": a tenth of the coverage that half the bases are at or above (the segments of two k or more, or all of
/// them if none is so long), but at least 1.5, which takes away what was seen once or twice, and never more than half of that coverage.
pub fn auto_cutoff(g: &Graph) -> f64 {
    let long = |min: usize| -> Vec<(f64, usize)> {
        g.ids()
            .into_iter()
            .filter(|&id| g.length(id) >= min)
            .map(|id| (g.coverage(id), g.length(id)))
            .collect()
    };
    let mut by_cov = long(2 * g.k());
    if by_cov.is_empty() {
        by_cov = long(0);
    }
    by_cov.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: usize = by_cov.iter().map(|c| c.1).sum();
    let mut seen = 0;
    let median = by_cov
        .iter()
        .find(|c| {
            seen += c.1;
            seen * 2 >= total
        })
        .map_or(0.0, |c| c.0);
    (0.1 * median).max(1.5).min(0.5 * median)
}

#[cfg(test)]
mod tests;
