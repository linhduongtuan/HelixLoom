//! The search for the species tree. ASTRAL's dynamic programme finds, among the species trees whose every cluster is in a
//! given set, the one with the most quartets in common with the gene trees: the best tree on a cluster is the best way to
//! split it in two clusters of the set, the best trees of the two parts, and the quartets that the split itself agrees on
//! (`weight::tripartition`). A last polish swaps neighbouring groups (nearest-neighbour interchanges) while that helps.

use super::bits::Bits;
use super::complete::{self, Space};
use super::genes::Genes;
use super::species::Sp;
use super::weight::tripartition;
use crate::seq::par::par_map;

/// The most taxa the exact search takes: it looks at every group of taxa, which doubles with each taxon.
pub const EXACT_LIMIT: usize = 12;
/// Most sweeps of the polish.
const POLISH_SWEEPS: usize = 30;

pub struct Options {
    /// Look at every group of taxa, not only the ones in the enriched set of the gene trees' splits.
    pub exact: bool,
    /// 0: the gene trees' (completed) splits alone; 1: with the ways to resolve polytomies.
    pub extra: u8,
    pub polylimit: usize,
    /// Swap neighbouring groups after the programme as long as the score rises.
    pub polish: bool,
    /// The taxon the tree hangs from; it only decides how the tree is written.
    pub root: usize,
}

pub struct Found {
    pub sp: Sp,
    pub score: u64,
    /// How many clusters the programme could choose from.
    pub clusters: usize,
}

/// Every non-empty group of the taxa other than the root one.
fn exact_space(n: usize, root: usize) -> Result<Space, String> {
    if n > EXACT_LIMIT {
        return Err(format!(
            "The exact search takes at most {EXACT_LIMIT} taxa; there are {n}."
        ));
    }
    let others: Vec<usize> = (0..n).filter(|&t| t != root).collect();
    let mut space = Space::new();
    for mask in 1u32..1 << others.len() {
        let set = others
            .iter()
            .enumerate()
            .filter(|(k, _)| mask >> k & 1 == 1)
            .map(|(_, &t)| t);
        space.add(Bits::of(n, set));
    }
    Ok(space)
}

pub fn search(genes: &Genes, o: &Options) -> Result<Found, String> {
    let n = genes.taxa.len();
    let space = if o.exact {
        exact_space(n, o.root)?
    } else {
        complete::build(genes, o.root, o.extra, o.polylimit)
    };
    let sp = best_tree(genes, &space, o.root)?;
    let sp = if o.polish { polish(sp, genes) } else { sp };
    let score = sp.score(genes);
    Ok(Found {
        sp,
        score,
        clusters: space.sets.len(),
    })
}

/// Every way to split a cluster of the space in two clusters of the space: (cluster, one part, the other).
fn splits(space: &Space) -> Vec<(usize, usize, usize)> {
    let sizes: Vec<u32> = space.sets.iter().map(Bits::count).collect();
    let mut order: Vec<usize> = (0..sizes.len()).collect();
    order.sort_by_key(|&i| sizes[i]);
    let mut out = Vec::new();
    for &v in &order {
        if sizes[v] < 2 {
            continue;
        }
        for &a in order.iter().take_while(|&&a| 2 * sizes[a] <= sizes[v]) {
            if !space.sets[a].is_subset_of(&space.sets[v]) {
                continue;
            }
            let rest = space.sets[v].minus(&space.sets[a]);
            let Some(b) = space.find(&rest) else { continue };
            // the same pair is met twice when the two parts are as big as each other
            if sizes[a] == sizes[b] && a > b {
                continue;
            }
            out.push((v, a, b));
        }
    }
    out
}

fn best_tree(genes: &Genes, space: &Space, root: usize) -> Result<Sp, String> {
    let n = genes.taxa.len();
    let triples = splits(space);
    let weights = par_map(&triples, |&(_, a, b)| {
        tripartition(&genes.genes, &space.sets[a], &space.sets[b])
    });
    // best[v]: the highest sum of weights below the cluster, and the split that gives it
    let mut best: Vec<Option<(u64, usize, usize)>> = space
        .sets
        .iter()
        .map(|s| (s.count() == 1).then_some((0, 0, 0)))
        .collect();
    for (&(v, a, b), &w) in triples.iter().zip(&weights) {
        let (Some((x, _, _)), Some((y, _, _))) = (best[a], best[b]) else {
            continue;
        };
        let total = x + y + w;
        if best[v].is_none_or(|(have, _, _)| total > have) {
            best[v] = Some((total, a, b));
        }
    }
    let whole = Bits::of(n, (0..n).filter(|&t| t != root));
    let top = space
        .find(&whole)
        .filter(|&v| best[v].is_some())
        .ok_or("No species tree can be built from the clusters the gene trees give.")?;
    Ok(rebuild(space, &best, top, n, root))
}

/// The tree that the programme's choices describe, from the cluster `top` down.
fn rebuild(
    space: &Space,
    best: &[Option<(u64, usize, usize)>],
    top: usize,
    n: usize,
    root: usize,
) -> Sp {
    let mut kids: Vec<Option<(usize, usize)>> = vec![None];
    let mut taxon: Vec<Option<usize>> = vec![None];
    let mut work = vec![(top, 0usize)];
    while let Some((cluster, here)) = work.pop() {
        match best[cluster] {
            Some((_, a, b)) if space.sets[cluster].count() > 1 => {
                let (ia, ib) = (kids.len(), kids.len() + 1);
                kids.extend([None, None]);
                taxon.extend([None, None]);
                kids[here] = Some((ia, ib));
                work.push((a, ia));
                work.push((b, ib));
            }
            _ => taxon[here] = space.sets[cluster].first(),
        }
    }
    Sp::from_shape(n, root, &kids, &taxon, 0)
}

/// Swap neighbouring groups while the score rises: at a branch with the groups `l r | s o`, the other two ways to pair them
/// change only the weights of the two nodes at its ends.
pub fn polish(mut sp: Sp, genes: &Genes) -> Sp {
    let weight = |a: &Bits, b: &Bits| tripartition(&genes.genes, a, b);
    for _ in 0..POLISH_SWEEPS {
        let mut changed = false;
        for v in 0..sp.nodes.len() {
            if let Some(better) = better_swap(&sp, v, &weight) {
                sp = better;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    sp
}

/// The tree after the swap at the branch above `v` that raises the score most, if any does.
fn better_swap(sp: &Sp, v: usize, weight: &impl Fn(&Bits, &Bits) -> u64) -> Option<Sp> {
    let (l, r) = sp.nodes[v].kids?;
    let u = sp.nodes[v].parent?;
    let (x, y) = sp.nodes[u].kids?;
    let s = if x == v { y } else { x };
    let set = |i: usize| &sp.nodes[i].set;
    let now = weight(set(l), set(r)) + weight(set(v), set(s));
    // moving `a` out of `v` and `s` in, or `b` out
    let try_swap = |a: usize, b: usize| {
        let joined = set(b).union(set(s));
        (weight(set(b), set(s)) + weight(&joined, set(a)), a, b)
    };
    let (gain, out, stay) = [try_swap(l, r), try_swap(r, l)]
        .into_iter()
        .max_by_key(|t| t.0)?;
    if gain <= now {
        return None;
    }
    let mut kids: Vec<Option<(usize, usize)>> = sp.nodes.iter().map(|n| n.kids).collect();
    kids[v] = Some((stay, s));
    kids[u] = Some((v, out));
    let taxon: Vec<Option<usize>> = sp.nodes.iter().map(|n| n.taxon).collect();
    Some(Sp::from_shape(sp.n, sp.root_taxon, &kids, &taxon, sp.root))
}
