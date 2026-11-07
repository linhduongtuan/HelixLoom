//! The set of clusters the species-tree search may build from. ASTRAL looks only at species trees whose every branch splits
//! the taxa the way some tree in a set of trees does; the set is the gene trees' splits, enriched (ASTRAL-II) in three ways:
//! gene trees with missing taxa are completed by putting the missing ones where a guide tree puts them, a polytomy in the
//! greedy consensus of the completed trees gets all the ways to join two of its groups, and so do the polytomies of the gene
//! trees. The guide tree joins the taxa by their average distance in the gene trees (neighbour joining).
//!
//! A cluster is a group of taxa that never has the root taxon: the tree is rooted at that leaf, so every split is read as
//! its side without it.

use super::bits::Bits;
use super::genes::{GTree, Genes};
use crate::seq::par::par_map;
use crate::seq::tree::{neighbor_joining, splits};
use std::collections::HashMap;

/// Most children of a polytomy of which the pairs are joined.
pub const POLYLIMIT: usize = 8;
/// More taxa than this and the guide tree is left out (neighbour joining takes cubic time).
pub const GUIDE_LIMIT: usize = 1_500;

#[derive(Debug)]
pub struct Space {
    pub sets: Vec<Bits>,
    index: HashMap<Bits, usize>,
}

impl Space {
    pub fn new() -> Space {
        Space {
            sets: Vec::new(),
            index: HashMap::new(),
        }
    }

    /// Add a cluster; false when it was there already (or empty).
    pub fn add(&mut self, set: Bits) -> bool {
        if set.is_empty() || self.index.contains_key(&set) {
            return false;
        }
        self.index.insert(set.clone(), self.sets.len());
        self.sets.push(set);
        true
    }

    pub fn find(&self, set: &Bits) -> Option<usize> {
        self.index.get(set).copied()
    }
}

/// The side of a split that does not have the root taxon.
fn side(set: Bits, root: usize, n: usize) -> Bits {
    if set.contains(root) {
        set.complement(n)
    } else {
        set
    }
}

/// The most clusters the search may choose from. It looks at pairs of them, so a million is hours of work and memory; more
/// happens when the gene trees share so few species that what each lacks is put back in a way of its own (above all with more
/// species than `GUIDE_LIMIT`, when there is no guide tree to put them back by).
pub const MAX_CLUSTERS: usize = 1_000_000;

/// The clusters for a search rooted at taxon `root`, and whether the gene trees were completed (they are not without a guide
/// tree, which more than `GUIDE_LIMIT` species leave out: what each gene lacks would then be put back in a way of its own, and
/// the groups to choose from would be as many as the genes times the species). With `extra` 0 only the (completed) gene
/// trees' splits, the guide tree's, the single taxa and the whole are in; with 1 the polytomy resolutions are added. More than
/// `MAX_CLUSTERS` of them is refused, as soon as it is known.
pub fn build(
    genes: &Genes,
    root: usize,
    extra: u8,
    polylimit: usize,
) -> Result<(Space, bool), String> {
    build_within(genes, root, extra, polylimit, MAX_CLUSTERS)
}

pub(super) fn build_within(
    genes: &Genes,
    root: usize,
    extra: u8,
    polylimit: usize,
    most: usize,
) -> Result<(Space, bool), String> {
    let n = genes.taxa.len();
    let mut space = Space::new();
    for t in (0..n).filter(|&t| t != root) {
        space.add(Bits::of(n, [t]));
    }
    space.add(Bits::of(n, (0..n).filter(|&t| t != root)));
    let guide = guide_clusters(genes, root);
    for c in &guide {
        space.add(c.clone());
    }
    let completed = !guide.is_empty();
    let trees = if completed {
        completed_trees(genes, &guide, root)
    } else {
        // so that some tree can always be built from the groups, however few the gene trees give: the taxa added one by one
        let mut so_far = Bits::new(n);
        for t in (0..n).filter(|&t| t != root) {
            so_far.insert(t);
            space.add(so_far.clone());
        }
        genes.trees.clone()
    };
    let limit = if extra >= 1 { polylimit } else { 0 };
    for t in &trees {
        add_tree(&mut space, t, n, root, limit);
        too_many(&space, most)?;
    }
    if extra >= 1 {
        consensus_resolutions(&mut space, &trees, n, root, polylimit);
        too_many(&space, most)?;
    }
    Ok((space, completed))
}

fn too_many(space: &Space, most: usize) -> Result<(), String> {
    if space.sets.len() > most {
        return Err(format!(
            "The gene trees give more than {most} groups of species to choose a species tree from, which is more than the search can look through. This happens when the gene trees share few species with each other: use genes that have more species in common, or fewer species."
        ));
    }
    Ok(())
}

/// The clusters of the guide tree: neighbour joining on the average number of branches between two taxa over the gene
/// trees that have both (a pair no tree has both of gets the overall average).
pub(super) fn guide_clusters(genes: &Genes, root: usize) -> Vec<Bits> {
    let n = genes.taxa.len();
    if n > GUIDE_LIMIT {
        return Vec::new();
    }
    let mut sum = vec![vec![0.0; n]; n];
    let mut seen = vec![vec![0.0; n]; n];
    for tree in &genes.trees {
        add_distances(tree, &mut sum, &mut seen);
    }
    let (total, pairs) = (
        sum.iter().flatten().sum::<f64>(),
        seen.iter().flatten().sum::<f64>(),
    );
    let fill = if pairs > 0.0 { total / pairs } else { 1.0 };
    let dist: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| match (i == j, seen[i][j] > 0.0) {
                    (true, _) => 0.0,
                    (false, true) => sum[i][j] / seen[i][j],
                    (false, false) => fill,
                })
                .collect()
        })
        .collect();
    let Some(joined) = neighbor_joining(&genes.taxa, &dist) else {
        return Vec::new();
    };
    // a hash set has no order, and the order the clusters are added in decides which of several equally good trees is
    // built and how its branches are listed: sorted, the same gene trees give the same tree every time
    let mut found: Vec<Vec<usize>> = splits(&joined, &genes.taxa).into_iter().collect();
    found.sort();
    found
        .into_iter()
        .map(|s| side(Bits::of(n, s), root, n))
        .collect()
}

/// Add the number of branches between every two leaves of a tree to the sums (kept in the upper triangle, mirrored).
fn add_distances(tree: &GTree, sum: &mut [Vec<f64>], seen: &mut [Vec<f64>]) {
    let mut depth = vec![0u32; tree.nodes.len()];
    let mut stack = vec![tree.root];
    while let Some(i) = stack.pop() {
        for &c in &tree.nodes[i].children {
            depth[c] = depth[i] + 1;
            stack.push(c);
        }
    }
    let leaves: Vec<(usize, usize)> = tree
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(i, n)| n.taxon.map(|t| (i, t)))
        .collect();
    for (k, &(a, ta)) in leaves.iter().enumerate() {
        for &(b, tb) in &leaves[k + 1..] {
            let d = f64::from(path_length(tree, &depth, a, b));
            for (x, y) in [(ta, tb), (tb, ta)] {
                sum[x][y] += d;
                seen[x][y] += 1.0;
            }
        }
    }
}

fn path_length(tree: &GTree, depth: &[u32], mut a: usize, mut b: usize) -> u32 {
    let mut steps = 0;
    while a != b {
        let up = if depth[a] >= depth[b] { &mut a } else { &mut b };
        match tree.nodes[*up].parent {
            Some(p) => *up = p,
            None => break,
        }
        steps += 1;
    }
    steps
}

/// Every gene tree with the taxa it lacks put back, each on the branch above the smallest group of the gene tree that holds
/// its nearest relatives in the guide tree.
pub(super) fn completed_trees(genes: &Genes, guide: &[Bits], root: usize) -> Vec<GTree> {
    let n = genes.taxa.len();
    let mut by_size: Vec<&Bits> = guide.iter().collect();
    by_size.sort_by_key(|c| c.count());
    // every tree is completed on its own, so they are shared out among the cores, and come back in their order
    par_map(&genes.trees, |t| complete_one(t, &by_size, root, n))
}

/// One gene tree with the taxa it lacks put back, the guide tree's clusters being given smallest first.
fn complete_one(gene: &GTree, by_size: &[&Bits], root: usize, n: usize) -> GTree {
    let mut tree = gene.clone();
    // the taxa below each node, kept up to date as taxa are put back instead of worked out again each time (which made the
    // completion cost the cube of the species for genes that share few of them)
    let mut clusters = tree.clusters(n);
    let mut present = clusters[tree.root].clone();
    while let Some((missing, near)) = next_missing(&present, by_size, root, n) {
        let at = (0..tree.nodes.len())
            .filter(|&i| near.is_subset_of(&clusters[i]))
            .min_by_key(|&i| clusters[i].count())
            .unwrap_or(tree.root);
        tree.insert_above(at, missing);
        // the new leaf and the join above `at` (which `insert_above` numbers in this order), and every node above them has
        // the taxon too
        let mut join = clusters[at].clone();
        join.insert(missing);
        clusters.push(Bits::of(n, [missing]));
        clusters.push(join);
        let mut above = tree.nodes[clusters.len() - 1].parent;
        while let Some(p) = above {
            clusters[p].insert(missing);
            above = tree.nodes[p].parent;
        }
        present.insert(missing);
    }
    tree
}

/// A taxon the tree lacks (`present` has the taxa it has), with the taxa of the tree that are its nearest relatives in the
/// guide tree.
fn next_missing(present: &Bits, guide: &[&Bits], root: usize, n: usize) -> Option<(usize, Bits)> {
    let missing = (0..n).find(|&t| !present.contains(t))?;
    let wanted = guide
        .iter()
        .find(|c| c.contains(missing) && c.and_count(present) > 0)
        .map(|c| c.minus(&c.minus(present)));
    // with no guide tree, or when the nearest relative is the root taxon itself
    let near = wanted.unwrap_or_else(|| {
        if present.contains(root) && present.count() == 1 {
            present.clone()
        } else {
            present.minus(&Bits::of(n, [root]))
        }
    });
    Some((
        missing,
        if near.is_empty() {
            present.clone()
        } else {
            near
        },
    ))
}

/// Add the splits of a tree, and, with a `limit` above 0, every way to join two branches of a node that has more than
/// three (a polytomy), as the unrooted tree sees its branches: its children and the rest of the tree.
fn add_tree(space: &mut Space, tree: &GTree, n: usize, root: usize, limit: usize) {
    let clusters = tree.clusters(n);
    let all = clusters[tree.root].clone();
    for (i, node) in tree.nodes.iter().enumerate() {
        space.add(side(clusters[i].clone(), root, n));
        if limit > 0 && node.children.len() > 2 {
            let mut branches: Vec<Bits> =
                node.children.iter().map(|&c| clusters[c].clone()).collect();
            branches.push(all.minus(&clusters[i]));
            branches.retain(|b| !b.is_empty());
            join_pairs(space, &branches, n, root, limit);
        }
    }
}

/// Every union of two of the branches of a node, for a node with more than three branches and at most `limit`.
fn join_pairs(space: &mut Space, branches: &[Bits], n: usize, root: usize, limit: usize) {
    if branches.len() <= 3 || branches.len() > limit {
        return;
    }
    for (k, a) in branches.iter().enumerate() {
        for b in &branches[k + 1..] {
            space.add(side(a.union(b), root, n));
        }
    }
}

/// The greedy consensus of the trees' splits (the most frequent first, each kept when it fits those kept), and the ways
/// to join two groups of each polytomy it has left.
fn consensus_resolutions(space: &mut Space, trees: &[GTree], n: usize, root: usize, limit: usize) {
    let mut count: HashMap<Bits, u32> = HashMap::new();
    for t in trees {
        for c in t.clusters(n) {
            *count.entry(side(c, root, n)).or_insert(0) += 1;
        }
    }
    let mut ranked: Vec<(Bits, u32)> = count.into_iter().collect();
    ranked.sort_by(|a, b| {
        b.1.cmp(&a.1)
            .then(b.0.count().cmp(&a.0.count()))
            .then(a.0.cmp(&b.0))
    });
    let mut kept: Vec<Bits> = Vec::new();
    for (c, _) in ranked {
        let fits = kept
            .iter()
            .all(|k| k.is_disjoint(&c) || c.is_subset_of(k) || k.is_subset_of(&c));
        if fits {
            kept.push(c);
        }
    }
    kept.push(Bits::of(n, (0..n).filter(|&t| t != root)));
    for k in &kept {
        let inside: Vec<&Bits> = kept
            .iter()
            .filter(|c| c.is_subset_of(k) && *c != k)
            .collect();
        let children: Vec<&Bits> = inside
            .iter()
            .copied()
            .filter(|c| !inside.iter().any(|d| *d != *c && c.is_subset_of(d)))
            .collect();
        let covered = children.iter().fold(Bits::new(n), |all, c| all.union(c));
        let mut branches: Vec<Bits> = children.iter().map(|c| (*c).clone()).collect();
        branches.extend(k.minus(&covered).iter().map(|t| Bits::of(n, [t])));
        branches.push(k.complement(n));
        join_pairs(space, &branches, n, root, limit);
    }
}
