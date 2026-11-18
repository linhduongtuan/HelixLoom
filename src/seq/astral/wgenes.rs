//! Gene trees with what weighted ASTRAL reads of them: the support and the length of every branch (Zhang and Mirarab 2022,
//! *Weighting by gene tree uncertainty improves accuracy of quartet-based species trees*, `docs/SPECIES.md`).
//!
//! The weight of a quartet `ab|cd` of a gene tree is, by the paper, the probability that the tree is right about it, times
//! how little the tree has to change to join the species of each pair:
//!
//! * by support: `1 − Π (1 − s(e))` over the branches `e` between the two nodes where the pairs meet;
//! * by length: `exp(−(l(a,b) + l(c,d)))`, `l` being the whole length of the path between two leaves;
//! * both: the product of the two.
//!
//! A quartet that the tree does not resolve (a polytomy joins its four species) has no weight, as in ASTRAL.

use crate::seq::newick::NwTree;
use std::collections::HashMap;

/// How the quartets of the gene trees are weighted.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Weights {
    /// Every quartet counts one: ASTRAL.
    #[default]
    None,
    /// By the support of the branches between the pairs (wASTRAL-s in the paper).
    Support,
    /// By the lengths of the paths of the pairs (wASTRAL-bl).
    Length,
    /// By both (wASTRAL-h, the paper's general choice).
    Hybrid,
}

impl Weights {
    pub const ALL: [Weights; 4] = [
        Weights::None,
        Weights::Support,
        Weights::Length,
        Weights::Hybrid,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Weights::None => "none",
            Weights::Support => "support",
            Weights::Length => "length",
            Weights::Hybrid => "hybrid",
        }
    }

    pub fn parse(word: &str) -> Option<Weights> {
        let word = word.trim().to_ascii_lowercase();
        Weights::ALL.into_iter().find(|w| w.name() == word)
    }

    pub fn uses_support(self) -> bool {
        matches!(self, Weights::Support | Weights::Hybrid)
    }

    pub fn uses_length(self) -> bool {
        matches!(self, Weights::Length | Weights::Hybrid)
    }

    /// What the log says of a run with these weights.
    pub fn describe(self) -> &'static str {
        match self {
            Weights::None => "every quartet counts one",
            Weights::Support => {
                "quartets weighted by the support of the branches between their pairs"
            }
            Weights::Length => "quartets weighted by the lengths of the paths of their pairs",
            Weights::Hybrid => "quartets weighted by both the support and the lengths",
        }
    }
}

/// The branch from a node to its parent, as the weights use it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WEdge {
    /// `exp(−length)` when lengths weigh, else 1.
    pub decay: f64,
    /// `1 − support` of an inner branch when supports weigh, else 1 (and always 1 for a branch to a leaf, which no pair's path
    /// between anchors ever holds).
    pub miss: f64,
}

#[derive(Clone, Debug)]
pub struct WNode {
    pub taxon: Option<usize>,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    /// The branch to the parent (meaningless for the root).
    pub edge: WEdge,
}

/// A gene tree held unrooted, as a tree rooted at an inner node: the nodes are numbered with parents before children (the
/// root is 0), and every node but the root has a parent edge. The root of a rooted tree and every node with one child are
/// taken out, their two branches made one.
#[derive(Clone, Debug)]
pub struct WTree {
    pub nodes: Vec<WNode>,
}

/// What the trees did not say, and the weights took as they are when something is not said: a branch without a support is
/// taken as fully supported, a branch without a length as of length 0.
#[derive(Clone, Copy, Debug, Default)]
pub struct Missing {
    pub supports: usize,
    pub lengths: usize,
}

impl Missing {
    pub fn notes(&self) -> Vec<String> {
        let one = |n: usize, what: &str, taken: &str| {
            let (es, was) = if n == 1 { ("", "was") } else { ("es", "were") };
            format!("{n} {what}{es} of the gene trees {was} not given and taken as {taken}")
        };
        let mut notes = Vec::new();
        if self.supports > 0 {
            notes.push(one(self.supports, "inner branch support", "full support"));
        }
        if self.lengths > 0 {
            notes.push(one(self.lengths, "branch length", "0"));
        }
        notes
    }
}

/// A branch on the way to being a `WEdge`.
#[derive(Clone, Copy)]
struct Link {
    to: usize,
    length: Option<f64>,
    support: Option<f64>,
}

/// The one value of two branches that are made one, if either has it.
fn merged(a: Option<f64>, b: Option<f64>, join: fn(f64, f64) -> f64) -> Option<f64> {
    match (a, b) {
        (Some(x), Some(y)) => Some(join(x, y)),
        (x, None) | (None, x) => x,
    }
}

/// The support a label gives: none for no label, the number otherwise.
fn support_of(label: &str, tree: usize) -> Result<Option<f64>, String> {
    let label = label.trim();
    if label.is_empty() {
        return Ok(None);
    }
    if label.contains('/') {
        return Err(format!(
            "Gene tree {tree} has the branch label {label}, which holds more than one number: keep one support value in each label."
        ));
    }
    match label.parse::<f64>() {
        Ok(s) if s.is_finite() && s >= 0.0 => Ok(Some(s)),
        _ => Err(format!(
            "Gene tree {tree} has the branch label {label}, which is not a support value (a number of at least 0)."
        )),
    }
}

/// The largest support in any of the trees, each with its number (counted from 1) among all the gene trees: values above 1
/// are percentages.
pub fn largest_support(trees: &[(usize, &NwTree)]) -> Result<f64, String> {
    let mut most: f64 = 0.0;
    for &(k, t) in trees {
        for (v, n) in t.nodes.iter().enumerate() {
            if v != t.root && !n.children.is_empty() {
                if let Some(s) = support_of(&n.name, k)? {
                    most = most.max(s);
                }
            }
        }
    }
    Ok(most)
}

impl WTree {
    /// The weights' view of gene tree number `k` (counted from 1), read with `NwTree::parse_all_raw` (a length that is not
    /// written is NaN). `scale` turns a support into a probability (1, or 0.01 for percentages); what the tree did not say is
    /// added to `missing`.
    pub fn build(
        t: &NwTree,
        index: &HashMap<&str, usize>,
        w: Weights,
        scale: f64,
        k: usize,
        missing: &mut Missing,
    ) -> Result<WTree, String> {
        let m = t.nodes.len();
        let leaf: Vec<bool> = t.nodes.iter().map(|n| n.children.is_empty()).collect();
        let mut adj: Vec<Vec<Link>> = vec![Vec::new(); m];
        for (v, n) in t.nodes.iter().enumerate() {
            let Some(p) = n.parent else { continue };
            if n.dist.is_infinite() || n.dist < 0.0 {
                return Err(format!(
                    "Gene tree {k} has a branch of length {}: lengths are numbers of at least 0.",
                    n.dist
                ));
            }
            let length = (!n.dist.is_nan()).then_some(n.dist);
            let support = if leaf[v] || !w.uses_support() {
                None
            } else {
                support_of(&n.name, k)?.map(|s| s * scale)
            };
            if support.is_some_and(|s| s > 1.0 + 1e-9) {
                return Err(format!(
                    "Gene tree {k} has the support {}: supports are probabilities (0 to 1) or percentages (0 to 100).",
                    n.name.trim()
                ));
            }
            for (from, to) in [(v, p), (p, v)] {
                adj[from].push(Link {
                    to,
                    length,
                    support,
                });
            }
        }
        // a node that is not a leaf and joins one branch only (the root of `((A,B,C,D))`) is cut away, then a node that joins
        // two (the root of a rooted tree, or a node with one child) is bridged: its two branches become one
        let mut work: Vec<usize> = (0..m).filter(|&v| !leaf[v] && adj[v].len() == 1).collect();
        while let Some(v) = work.pop() {
            if let Some(link) = adj[v].pop() {
                adj[link.to].retain(|l| l.to != v);
                if !leaf[link.to] && adj[link.to].len() == 1 {
                    work.push(link.to);
                }
            }
        }
        for v in 0..m {
            if adj[v].len() != 2 || leaf[v] {
                continue;
            }
            let (x, y) = (adj[v][0], adj[v][1]);
            let (length, support) = (
                merged(x.length, y.length, |a, b| a + b),
                merged(x.support, y.support, f64::max),
            );
            for (from, to) in [(x.to, y.to), (y.to, x.to)] {
                if let Some(l) = adj[from].iter_mut().find(|l| l.to == v) {
                    *l = Link {
                        to,
                        length,
                        support,
                    };
                }
            }
            adj[v].clear();
        }
        let first = (0..m)
            .find(|&v| !leaf[v] && adj[v].len() >= 3)
            .ok_or_else(|| format!("Gene tree {k} has no inner node that joins three branches."))?;
        // the tree from `first`, numbered with parents before children
        let mut up: Vec<Option<(usize, Link)>> = vec![None; m];
        let mut seen = vec![false; m];
        let mut order = vec![first];
        seen[first] = true;
        let mut at = 0;
        while at < order.len() {
            let v = order[at];
            at += 1;
            for &l in &adj[v] {
                if !seen[l.to] {
                    seen[l.to] = true;
                    up[l.to] = Some((v, l));
                    order.push(l.to);
                }
            }
        }
        let mut place = vec![usize::MAX; m];
        for (i, &v) in order.iter().enumerate() {
            place[v] = i;
        }
        let mut nodes: Vec<WNode> = order
            .iter()
            .map(|&v| WNode {
                taxon: leaf[v]
                    .then(|| index.get(t.nodes[v].name.as_str()).copied())
                    .flatten(),
                parent: up[v].map(|(p, _)| place[p]),
                children: Vec::new(),
                edge: WEdge {
                    decay: 1.0,
                    miss: 1.0,
                },
            })
            .collect();
        let (mut supports, mut inner, mut lengths, mut edges) = (0usize, 0usize, 0usize, 0usize);
        let mut here = Missing::default();
        for (i, &v) in order.iter().enumerate() {
            let Some((p, link)) = up[v] else { continue };
            nodes[place[p]].children.push(i);
            edges += 1;
            let length = match link.length {
                Some(l) => {
                    lengths += 1;
                    l
                }
                None => {
                    here.lengths += 1;
                    0.0
                }
            };
            let miss = if !w.uses_support() || leaf[v] {
                1.0
            } else {
                inner += 1;
                match link.support {
                    Some(s) => {
                        supports += 1;
                        1.0 - s
                    }
                    None => {
                        here.supports += 1;
                        0.0
                    }
                }
            };
            let decay = if w.uses_length() {
                (-length).exp()
            } else {
                1.0
            };
            nodes[i].edge = WEdge { decay, miss };
        }
        if w.uses_support() && inner > 0 && supports == 0 {
            return Err(format!(
                "Gene tree {k} has no branch support: the labels of its inner branches must hold the supports (bootstrap, aLRT, aBayes and so on) for the weights \"support\" and \"hybrid\"."
            ));
        }
        if w.uses_length() && lengths == 0 && edges > 0 {
            return Err(format!(
                "Gene tree {k} has no branch lengths, which the weights \"length\" and \"hybrid\" use."
            ));
        }
        missing.supports += here.supports;
        // lengths only count when they weigh
        if w.uses_length() {
            missing.lengths += here.lengths;
        }
        Ok(WTree { nodes })
    }
}
