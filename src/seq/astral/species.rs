//! The species tree: a binary tree on the taxa that is scored against the gene trees and annotated branch by branch.
//!
//! The tree is held rooted at one leaf, the root taxon: the node under that leaf holds every other taxon, and every node
//! knows the taxa below it. Which leaf the tree hangs from is a matter of display only; the score and the support of a
//! branch do not depend on it.

use super::bits::Bits;
use super::genes::Genes;
use super::special::{coalescent_length, posteriors};
use super::weight::{around_branch, tripartition};
use crate::seq::newick::NwTree;
use crate::seq::par::par_map;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct SpNode {
    /// The taxa below the node (never the root taxon).
    pub set: Bits,
    pub kids: Option<(usize, usize)>,
    pub taxon: Option<usize>,
    pub parent: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct Sp {
    pub n: usize,
    pub root_taxon: usize,
    pub nodes: Vec<SpNode>,
    /// The node under the root taxon's leaf.
    pub root: usize,
}

/// What the gene trees say about one inner branch of the species tree.
#[derive(Clone, Debug)]
pub struct Branch {
    /// The node under the branch; its taxa are one side of it.
    pub node: usize,
    /// The share of the genes' quartets around the branch for each topology: the branch's own, then the two others.
    pub f: [f64; 3],
    /// The same as fractions of `en`.
    pub q: [f64; 3],
    /// The local posterior probability of each topology.
    pub pp: [f64; 3],
    /// The effective number of genes.
    pub en: f64,
    /// The number of quartets the four groups around the branch make.
    pub quartets: u64,
    /// The length in coalescent units.
    pub length: f64,
}

impl Sp {
    /// Build from the shape of a tree: `kids[i]` are the two children of node `i` (none for a leaf, whose taxon is
    /// `taxon[i]`), and `root` is the node under the root taxon's leaf.
    pub fn from_shape(
        n: usize,
        root_taxon: usize,
        kids: &[Option<(usize, usize)>],
        taxon: &[Option<usize>],
        root: usize,
    ) -> Sp {
        let mut nodes: Vec<SpNode> = kids
            .iter()
            .zip(taxon)
            .map(|(&k, &t)| SpNode {
                set: t.map_or_else(|| Bits::new(n), |t| Bits::of(n, [t])),
                kids: k,
                taxon: t,
                parent: None,
            })
            .collect();
        // children before parents, without recursion
        let mut order = Vec::with_capacity(nodes.len());
        let mut stack = vec![root];
        while let Some(i) = stack.pop() {
            order.push(i);
            if let Some((a, b)) = nodes[i].kids {
                stack.extend([a, b]);
            }
        }
        for &i in order.iter().rev() {
            if let Some((a, b)) = nodes[i].kids {
                nodes[i].set = nodes[a].set.union(&nodes[b].set);
                nodes[a].parent = Some(i);
                nodes[b].parent = Some(i);
            }
        }
        Sp {
            n,
            root_taxon,
            nodes,
            root,
        }
    }

    /// Read a binary species tree (rooted or not) over exactly the gene trees' taxa and hang it from the first leaf.
    pub fn from_newick(text: &str, taxa: &[String]) -> Result<Sp, String> {
        let tree = NwTree::parse_all(text)?
            .into_iter()
            .next()
            .ok_or("No species tree found.")?;
        let index: HashMap<&str, usize> = taxa
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        let mut adj = adjacency(&tree);
        contract(&mut adj, &tree)?;
        let leaves: Vec<usize> = (0..tree.nodes.len()).filter(|&i| tree.is_leaf(i)).collect();
        let mut taxon_of = vec![None; tree.nodes.len()];
        for &l in &leaves {
            let name = tree.nodes[l].name.as_str();
            let t = *index.get(name).ok_or_else(|| {
                format!("The species tree has the taxon {name}, which no gene tree has.")
            })?;
            if taxon_of.iter().flatten().any(|&x| x == t) {
                return Err(format!("The species tree has the taxon {name} twice."));
            }
            taxon_of[l] = Some(t);
        }
        if let Some(missing) = (0..taxa.len()).find(|t| !taxon_of.iter().flatten().any(|x| x == t))
        {
            return Err(format!(
                "The species tree lacks the taxon {}, which the gene trees have.",
                taxa[missing]
            ));
        }
        let first = *leaves.first().ok_or("The species tree has no leaves.")?;
        Ok(hang_from(&adj, &taxon_of, first, taxa.len()))
    }

    /// Half the sum, over the inner nodes, of the quartets of the gene trees that agree with the node: the number of
    /// quartets of the species tree that the gene trees share (each counted once per gene tree).
    pub fn score(&self, genes: &Genes) -> u64 {
        let inner: Vec<(usize, usize)> = self.nodes.iter().filter_map(|n| n.kids).collect();
        let weights = par_map(&inner, |&(a, b)| {
            tripartition(&genes.genes, &self.nodes[a].set, &self.nodes[b].set)
        });
        weights.iter().sum::<u64>() / 2
    }

    /// What the gene trees say about every inner branch (all but the one under the root taxon's leaf).
    pub fn annotate(&self, genes: &Genes, lambda: f64) -> Vec<Branch> {
        let all = Bits::of(self.n, 0..self.n);
        let ends: Vec<(usize, Bits, Bits, Bits, Bits)> = self
            .nodes
            .iter()
            .enumerate()
            .filter(|&(i, n)| n.kids.is_some() && i != self.root)
            .filter_map(|(i, n)| {
                let (l, r) = n.kids?;
                let parent = &self.nodes[n.parent?];
                let sibling = match parent.kids? {
                    (a, _) if a != i => a,
                    (_, b) => b,
                };
                let (s, l, r) = (&self.nodes[sibling], &self.nodes[l], &self.nodes[r]);
                let o = all.minus(&n.set.union(&s.set));
                Some((i, l.set.clone(), r.set.clone(), s.set.clone(), o))
            })
            .collect();
        par_map(&ends, |(i, l, r, s, o)| {
            let a = around_branch(&genes.genes, l, r, s, o);
            let q = a.f.map(|x| if a.en > 0.0 { x / a.en } else { 0.0 });
            Branch {
                node: *i,
                f: a.f,
                q,
                pp: posteriors(a.f),
                en: a.en,
                quartets: a.quartets,
                length: coalescent_length(a.f[0], a.en, lambda),
            }
        })
    }
}

/// The edges of a parsed tree as lists of neighbours (a node's children first, then its parent).
fn adjacency(t: &NwTree) -> Vec<Vec<usize>> {
    t.nodes
        .iter()
        .map(|n| n.children.iter().copied().chain(n.parent).collect())
        .collect()
}

/// Join the two neighbours of every node that has just two, so that the root of a rooted tree and a node with one child
/// disappear; then every inner node must join three branches.
fn contract(adj: &mut [Vec<usize>], t: &NwTree) -> Result<(), String> {
    for v in 0..adj.len() {
        if adj[v].len() != 2 {
            continue;
        }
        let (x, y) = (adj[v][0], adj[v][1]);
        for (from, to) in [(x, y), (y, x)] {
            if let Some(slot) = adj[from].iter_mut().find(|e| **e == v) {
                *slot = to;
            }
        }
        adj[v].clear();
    }
    for (v, edges) in adj.iter().enumerate() {
        let leaf = t.is_leaf(v);
        if !edges.is_empty() && !leaf && edges.len() != 3 {
            return Err("The species tree must be binary: every inner node has to join three branches (no polytomies).".into());
        }
    }
    Ok(())
}

/// The tree seen from the leaf `first`, as a tree rooted at that leaf.
fn hang_from(adj: &[Vec<usize>], taxon_of: &[Option<usize>], first: usize, n: usize) -> Sp {
    let root_taxon = taxon_of[first].unwrap_or(0);
    let mut kids: Vec<Option<(usize, usize)>> = Vec::new();
    let mut taxon: Vec<Option<usize>> = Vec::new();
    // (node of the parsed tree, the node it came from, its place in the new tree)
    let top = adj[first][0];
    let mut work = vec![(top, first, 0usize)];
    kids.push(None);
    taxon.push(taxon_of[top]);
    while let Some((v, from, here)) = work.pop() {
        let next: Vec<usize> = adj[v].iter().copied().filter(|&x| x != from).collect();
        if let [a, b] = next[..] {
            let (ia, ib) = (kids.len(), kids.len() + 1);
            kids.extend([None, None]);
            taxon.extend([taxon_of[a], taxon_of[b]]);
            kids[here] = Some((ia, ib));
            work.push((a, v, ia));
            work.push((b, v, ib));
        }
    }
    Sp::from_shape(n, root_taxon, &kids, &taxon, 0)
}
