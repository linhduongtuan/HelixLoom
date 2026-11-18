//! Gene family trees, in the way of ASTRAL-Pro (Zhang, Scornavacca, Molloy and Mirarab 2020, *ASTRAL-Pro: quartet-based
//! species-tree inference despite paralogy*, `docs/SPECIES.md`): trees in which a species can occur several times, because
//! genes were duplicated or lost. The leaf labels are the species names, repeated once for each copy.
//!
//! * A family tree is rooted. A tree whose top joins two branches is taken as it is written; one with three (an unrooted
//!   tree) is rooted where the cost is least, counting one for every duplication and one for every loss that the species of
//!   the two branches of a node ask for (`rooting_cost`).
//! * A node is a **duplication** when the species below its two branches overlap, and a **speciation** otherwise.
//! * What the score reads of a family is, for every speciation node `w`, three sets of species: those below its two branches
//!   (`M1`, `M2`) and those of the leaves outside `w` whose last common ancestor with `w` is a speciation (`M3`): see
//!   `pscore.rs`.

use super::bits::Bits;
use super::genes::{check_size, GNode, GTree, Genes};
use crate::seq::newick::NwTree;
use std::collections::{BTreeSet, HashMap};

/// How many single-copy trees are sampled from each family for the groups of species the search chooses from: twice the
/// paper's default `C` of 4.
const SAMPLES: usize = 8;

/// One family tree, rooted and binary.
#[derive(Clone, Debug)]
pub struct Rooted {
    pub kids: Vec<Option<(usize, usize)>>,
    pub taxon: Vec<Option<usize>>,
    pub root: usize,
}

/// What the score reads at one speciation node of a family tree: the species below its two branches and the species of the
/// leaves outside it that meet it at a speciation.
#[derive(Clone, Debug)]
pub struct ProNode {
    pub m: [Bits; 3],
}

/// The family trees, as the score and the search read them.
#[derive(Debug)]
pub struct Families {
    pub taxa: Vec<String>,
    pub nodes: Vec<ProNode>,
    /// Single-copy trees sampled from the families, for the groups of species to choose from.
    pub sample: Genes,
    /// The family trees that were kept, rooted (for the tests that count their quartets one by one).
    #[cfg(test)]
    pub rooted: Vec<Rooted>,
    pub families: usize,
    /// Family trees left out for having too few species.
    pub skipped: usize,
    pub duplications: usize,
    /// Trees that had no root and were rooted by the cost of duplications and losses.
    pub rooted_by_cost: usize,
    /// The quartet classes of all the families, whatever the species tree (see `pscore.rs`).
    pub classes: u128,
}

fn c2(n: u128) -> u128 {
    n * n.saturating_sub(1) / 2
}

/// The undirected tree of a parsed family tree, with the root of a rooted tree taken out, the species of its leaves, and the
/// one branch that a rooted tree is rooted on (none for an unrooted tree).
struct Unrooted {
    adj: Vec<Vec<usize>>,
    taxon: Vec<Option<usize>>,
    fixed: Option<(usize, usize)>,
}

fn unrooted(t: &NwTree, index: &HashMap<&str, usize>, k: usize) -> Result<Unrooted, String> {
    let mut adj: Vec<Vec<usize>> = t
        .nodes
        .iter()
        .map(|n| n.children.iter().copied().chain(n.parent).collect())
        .collect();
    let taxon: Vec<Option<usize>> = t
        .nodes
        .iter()
        .map(|n| {
            n.children
                .is_empty()
                .then(|| index.get(n.name.as_str()).copied())
                .flatten()
        })
        .collect();
    let mut fixed = None;
    let top = &t.nodes[t.root];
    if top.children.len() == 2 {
        let (x, y) = (top.children[0], top.children[1]);
        for (from, to) in [(x, y), (y, x)] {
            if let Some(slot) = adj[from].iter_mut().find(|e| **e == t.root) {
                *slot = to;
            }
        }
        adj[t.root].clear();
        fixed = Some((x.min(y), x.max(y)));
    }
    for (v, edges) in adj.iter().enumerate() {
        let leaf = t.nodes[v].children.is_empty();
        let removed = v == t.root && fixed.is_some();
        if !removed && !leaf && edges.len() != 3 {
            return Err(format!(
                "Gene family tree {k} is not binary: a node joins {} branches. Resolve the polytomies first (an unrooted tree may have three branches at its top).",
                edges.len()
            ));
        }
    }
    Ok(Unrooted { adj, taxon, fixed })
}

/// The tree rooted on the branch between `a` and `b`.
fn rooted_on(u: &Unrooted, a: usize, b: usize) -> Rooted {
    let mut kids: Vec<Option<(usize, usize)>> = vec![None];
    let mut taxon: Vec<Option<usize>> = vec![None];
    // (node of the undirected tree, the node it was reached from, its place in the new tree)
    let mut work = Vec::new();
    for (node, from) in [(a, b), (b, a)] {
        let here = kids.len();
        kids.push(None);
        taxon.push(u.taxon[node]);
        work.push((node, from, here));
    }
    kids[0] = Some((1, 2));
    while let Some((v, from, here)) = work.pop() {
        let next: Vec<usize> = u.adj[v].iter().copied().filter(|&x| x != from).collect();
        if let [x, y] = next[..] {
            let (ix, iy) = (kids.len(), kids.len() + 1);
            kids.extend([None, None]);
            taxon.extend([u.taxon[x], u.taxon[y]]);
            kids[here] = Some((ix, iy));
            work.push((x, v, ix));
            work.push((y, v, iy));
        }
    }
    Rooted {
        kids,
        taxon,
        root: 0,
    }
}

impl Rooted {
    /// Nodes children before parents.
    pub(super) fn postorder(&self) -> Vec<usize> {
        let mut order = Vec::with_capacity(self.kids.len());
        let mut stack = vec![self.root];
        while let Some(v) = stack.pop() {
            order.push(v);
            if let Some((a, b)) = self.kids[v] {
                stack.extend([a, b]);
            }
        }
        order.reverse();
        order
    }

    /// The species below every node.
    pub(super) fn species(&self, n: usize) -> Vec<Bits> {
        let mut sets = vec![Bits::new(n); self.kids.len()];
        for v in self.postorder() {
            sets[v] = match (self.kids[v], self.taxon[v]) {
                (Some((a, b)), _) => sets[a].union(&sets[b]),
                (None, Some(t)) => Bits::of(n, [t]),
                (None, None) => Bits::new(n),
            };
        }
        sets
    }

    /// One for every duplication, and for every loss that the species of the two branches of a duplication ask for: none when
    /// they are the same, one when those of one branch are among those of the other, two otherwise (the way the paper roots a
    /// tree: speciations cost nothing, and so does incomplete lineage sorting).
    fn rooting_cost(&self, n: usize) -> usize {
        let sets = self.species(n);
        self.kids
            .iter()
            .flatten()
            .map(|&(a, b)| {
                let (x, y) = (&sets[a], &sets[b]);
                if x.is_disjoint(y) {
                    0
                } else if x == y {
                    1
                } else if x.is_subset_of(y) || y.is_subset_of(x) {
                    2
                } else {
                    3
                }
            })
            .sum()
    }

    /// The nodes that are speciations, by node number (a leaf is none), and how many nodes are duplications.
    pub(super) fn speciations(&self, sets: &[Bits]) -> (Vec<bool>, usize) {
        let mut dups = 0;
        let tags = self
            .kids
            .iter()
            .map(|k| match k {
                Some((a, b)) if sets[*a].is_disjoint(&sets[*b]) => true,
                Some(_) => {
                    dups += 1;
                    false
                }
                None => false,
            })
            .collect();
        (tags, dups)
    }

    /// What the score reads at each speciation, and how many duplications there are.
    fn pro_nodes(&self, n: usize) -> (Vec<ProNode>, usize) {
        let sets = self.species(n);
        let (speciation, dups) = self.speciations(&sets);
        // the species outside each node that meet it at a speciation, from the top down
        let mut outside = vec![Bits::new(n); self.kids.len()];
        let mut found = Vec::new();
        let mut stack = vec![self.root];
        while let Some(w) = stack.pop() {
            let Some((a, b)) = self.kids[w] else { continue };
            if speciation[w] {
                found.push(ProNode {
                    m: [sets[a].clone(), sets[b].clone(), outside[w].clone()],
                });
            }
            for (child, sibling) in [(a, b), (b, a)] {
                outside[child] = if speciation[w] {
                    outside[w].union(&sets[sibling])
                } else {
                    outside[w].clone()
                };
                stack.push(child);
            }
        }
        (found, dups)
    }

    /// The quartet classes this tree has, whatever the species tree: for each speciation, the pairs of its two branches
    /// against each other, and a pair in one branch with one species of the other and one from outside.
    fn classes(nodes: &[ProNode]) -> u128 {
        nodes
            .iter()
            .map(|p| {
                let [l, r, o] = [0, 1, 2].map(|i| u128::from(p.m[i].count()));
                c2(l) * c2(r) + (c2(l) * r + l * c2(r)) * o
            })
            .sum()
    }

    /// A tree with one copy of each species it keeps: at every duplication one of the two branches is dropped, chosen by
    /// `pick`. The tree of the leaves that are left, rooted.
    fn single_copy(&self, n: usize, mut pick: impl FnMut() -> bool) -> Option<GTree> {
        let sets = self.species(n);
        let (speciation, _) = self.speciations(&sets);
        // the nodes that are kept, parents before children: below a duplication only one of the two branches is
        let mut kept = Vec::new();
        let mut stack = vec![self.root];
        while let Some(v) = stack.pop() {
            kept.push(v);
            match self.kids[v] {
                Some((a, b)) if speciation[v] => stack.extend([a, b]),
                Some((a, b)) => stack.push(if pick() { a } else { b }),
                None => {}
            }
        }
        // the tree of what is kept, a node with one kept branch being the branch itself
        let mut nodes: Vec<GNode> = Vec::new();
        let mut made: Vec<Option<usize>> = vec![None; self.kids.len()];
        for &v in kept.iter().rev() {
            made[v] = match self.kids[v] {
                None => {
                    nodes.push(GNode {
                        parent: None,
                        children: Vec::new(),
                        taxon: self.taxon[v],
                    });
                    Some(nodes.len() - 1)
                }
                Some((a, b)) => match (made[a], made[b]) {
                    (Some(x), Some(y)) => {
                        nodes.push(GNode {
                            parent: None,
                            children: vec![x, y],
                            taxon: None,
                        });
                        let me = nodes.len() - 1;
                        nodes[x].parent = Some(me);
                        nodes[y].parent = Some(me);
                        Some(me)
                    }
                    (one, other) => one.or(other),
                },
            };
        }
        let root = made[self.root]?;
        (nodes.iter().filter(|g| g.children.is_empty()).count() >= 4)
            .then_some(GTree { nodes, root })
    }
}

/// A source of choices that is the same on every run.
struct Choices(u64);

impl Choices {
    fn next(&mut self) -> bool {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 40) & 1 == 1
    }
}

/// The most memory the sets of species of all the family trees may take, as for gene trees.
const MAX_SET_BYTES: u128 = 8 << 30;

/// Read family trees (one `;`-ended Newick tree after the other, the leaf labels being species names that may repeat).
/// A tree with fewer than `min_leaves` species (and always fewer than four) is left out.
pub fn read(text: &str, min_leaves: usize) -> Result<Families, String> {
    let parsed = NwTree::parse_all(text)?;
    let least = min_leaves.max(super::genes::MIN_LEAVES);
    let mut names = BTreeSet::new();
    let mut kept: Vec<(usize, &NwTree)> = Vec::new();
    let mut skipped = 0;
    for (k, tree) in parsed.iter().enumerate() {
        let leaves: Vec<&str> = tree
            .nodes
            .iter()
            .filter(|n| n.children.is_empty())
            .map(|n| n.name.as_str())
            .collect();
        if leaves.iter().any(|n| n.is_empty()) {
            return Err(format!("A leaf of gene family tree {} has no name.", k + 1));
        }
        let species: BTreeSet<&str> = leaves.iter().copied().collect();
        if species.len() < least {
            skipped += 1;
            continue;
        }
        names.extend(species.iter().map(|s| s.to_string()));
        kept.push((k + 1, tree));
    }
    if names.len() < super::genes::MIN_LEAVES {
        return Err("A species tree needs at least four species in the gene family trees.".into());
    }
    let taxa: Vec<String> = names.into_iter().collect();
    let n = taxa.len();
    let leaves: u128 = kept.iter().map(|(_, t)| t.nodes.len() as u128).sum();
    check_size(kept.len(), leaves * (1 + SAMPLES as u128), n, MAX_SET_BYTES)?;
    let index: HashMap<&str, usize> = taxa
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let (mut nodes, mut duplications, mut rooted_by_cost, mut classes) = (Vec::new(), 0, 0, 0u128);
    let (mut trees, mut sampled) = (Vec::new(), Vec::new());
    for (position, (k, t)) in kept.iter().enumerate() {
        let u = unrooted(t, &index, *k)?;
        let rooted = match u.fixed {
            Some((a, b)) => rooted_on(&u, a, b),
            None => {
                rooted_by_cost += 1;
                best_rooting(&u, n)
            }
        };
        let (pro, dups) = rooted.pro_nodes(n);
        duplications += dups;
        classes += Rooted::classes(&pro);
        nodes.extend(pro);
        let mut pick = Choices(position as u64 + 1);
        for _ in 0..SAMPLES {
            if let Some(tree) = rooted.single_copy(n, || pick.next()) {
                sampled.push(tree);
            }
        }
        trees.push(rooted);
    }
    if nodes.is_empty() {
        return Err("No quartet of four different species can be read from these gene family trees: no node of them joins branches with no species in common.".into());
    }
    if !super::genes::counts_fit(n, nodes.len()) {
        return Err(format!(
            "There are too many species ({n}) and speciation nodes ({}) for the search to count their quartets. Use fewer of them.",
            nodes.len()
        ));
    }
    let genes = sampled.iter().map(|t| t.gene(n)).collect();
    Ok(Families {
        taxa: taxa.clone(),
        nodes,
        sample: Genes {
            taxa,
            trees: sampled,
            genes,
            wtrees: Vec::new(),
            notes: Vec::new(),
            skipped: 0,
            input_trees: 0,
            origin: Vec::new(),
        },
        families: trees.len(),
        #[cfg(test)]
        rooted: trees,
        skipped,
        duplications,
        rooted_by_cost,
        classes,
    })
}

/// The rooting of an unrooted tree that costs least; the first of those that cost the same.
fn best_rooting(u: &Unrooted, n: usize) -> Rooted {
    let mut best: Option<(usize, Rooted)> = None;
    for a in 0..u.adj.len() {
        for &b in u.adj[a].iter().filter(|&&b| a < b) {
            let candidate = rooted_on(u, a, b);
            let cost = candidate.rooting_cost(n);
            if best.as_ref().is_none_or(|(have, _)| cost < *have) {
                best = Some((cost, candidate));
            }
        }
    }
    best.map(|(_, r)| r).unwrap_or(Rooted {
        kids: vec![None],
        taxon: vec![None],
        root: 0,
    })
}
