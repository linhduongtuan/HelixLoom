//! The weighted quartet score of one node of a species tree (weighted ASTRAL, `wgenes.rs`): the clusters `A`, `B` and `C`
//! around it, and for every quartet of a gene tree with two species in one cluster and one in each of the other two, resolved
//! as the pair against the two others, its weight.
//!
//! A quartet `x1 x2 | y z` has two nodes where its pairs meet: `u` for `x1 x2` and `m` for `y z`. Its weight is
//!
//! ```text
//! (1 − Π (1 − s(e)) over the branches e from u to m)  ·  exp(−d(x1,u) − d(x2,u) − d(y,m) − d(z,m))
//! ```
//!
//! (either factor left out, as the weights ask). The factors fall apart at the two nodes, so the sum over all quartets is
//! found without listing them: at every node `u` of the gene tree and every branch `k` leaving it, the pairs `x1 x2` on two
//! other branches (a sum of products of what each branch carries) times what the tree behind `k` gives for `y z`: the sum
//! over the nodes `m` behind `k` of the pairs that meet at `m`, with the product of `1 − s` from `u` to `m`. What lies behind
//! each branch, in each direction, is found for all branches in two passes over the tree.

use super::bits::Bits;
use super::wgenes::{WTree, Weights};

/// What lies behind one branch, seen from the node it leaves. Cluster `x` of the three (`A`, `B`, `C`) is `x` below; the
/// pair of clusters that joins at a node is named by the cluster that is left out of it, so that `role` is also the cluster
/// whose two species are the pair at the other end.
#[derive(Clone, Copy, Default)]
struct Behind {
    /// Species of the cluster behind the branch, each counted `exp(−distance)` from the node it is seen from.
    e: [f64; 3],
    /// The pairs `y z` (`y` of the first, `z` of the second cluster other than `role`) that meet behind the branch, each
    /// counted with `exp(−d(y,m) − d(z,m))`.
    pairs: [f64; 3],
    /// The same, each also multiplied by `Π (1 − s)` from the node the branch leaves to where the pair meets.
    missed: [f64; 3],
}

/// The two clusters other than `role`.
const OTHERS: [(usize, usize); 3] = [(1, 2), (0, 2), (0, 1)];

/// What a node holds of the branches around it that are summed, leaving one out: the clusters, the products of the two
/// clusters of a role on one branch, and the pairs.
#[derive(Clone, Copy, Default)]
struct Totals {
    e: [f64; 3],
    both: [f64; 3],
    pairs: [f64; 3],
    missed: [f64; 3],
}

impl Totals {
    fn add(&mut self, b: &Behind) {
        for (r, &(y, z)) in OTHERS.iter().enumerate() {
            self.e[r] += b.e[r];
            self.both[r] += b.e[y] * b.e[z];
            self.pairs[r] += b.pairs[r];
            self.missed[r] += b.missed[r];
        }
    }

    /// What the branch to the parent of a node shows, when it leaves out `skip` (one of the branches at the node, or none)
    /// and `miss`, `decay` are the branch's own factors.
    fn behind(&self, skip: &Behind, miss: f64, decay: f64) -> Behind {
        let mut out = Behind::default();
        let rest: [f64; 3] = std::array::from_fn(|x| self.e[x] - skip.e[x]);
        for (r, &(y, z)) in OTHERS.iter().enumerate() {
            out.e[r] = decay * rest[r];
            let meeting = rest[y] * rest[z] - (self.both[r] - skip.e[y] * skip.e[z]);
            out.pairs[r] = meeting + (self.pairs[r] - skip.pairs[r]);
            out.missed[r] = miss * (meeting + (self.missed[r] - skip.missed[r]));
        }
        out
    }
}

/// Which of the three clusters a species is in: `A`, `B`, or everything else.
fn cluster(taxon: usize, a: &Bits, b: &Bits) -> usize {
    if a.contains(taxon) {
        0
    } else if b.contains(taxon) {
        1
    } else {
        2
    }
}

/// The weighted quartets of one gene tree for the node with clusters `a`, `b` and the rest.
fn gene(
    t: &WTree,
    a: &Bits,
    b: &Bits,
    support: bool,
    down: &mut [Behind],
    up: &mut [Behind],
) -> f64 {
    // up the tree: what lies below each node, seen from its parent
    for v in (1..t.nodes.len()).rev() {
        let node = &t.nodes[v];
        down[v] = match node.taxon {
            Some(x) => {
                let mut leaf = Behind::default();
                leaf.e[cluster(x, a, b)] = node.edge.decay;
                leaf
            }
            None => {
                let mut sum = Totals::default();
                for &c in &node.children {
                    sum.add(&down[c]);
                }
                sum.behind(&Behind::default(), node.edge.miss, node.edge.decay)
            }
        };
    }
    // down the tree: what lies above each node, seen from it. At a node, the branches are its children and its own `up`.
    let mut total = 0.0;
    for p in 0..t.nodes.len() {
        let node = &t.nodes[p];
        if node.taxon.is_some() {
            continue;
        }
        let mut sum = Totals::default();
        for &c in &node.children {
            sum.add(&down[c]);
        }
        if node.parent.is_some() {
            sum.add(&up[p]);
        }
        for &c in &node.children {
            let edge = t.nodes[c].edge;
            up[c] = sum.behind(&down[c], edge.miss, edge.decay);
        }
        total += at_node(
            &sum,
            node.children
                .iter()
                .map(|&c| &down[c])
                .chain(node.parent.map(|_| &up[p])),
            support,
        );
    }
    total
}

/// The quartets whose pair `x1 x2` meets at one node, given the branches around it: for each branch `k`, the pairs of species
/// of a cluster on two other branches times the weight of the pairs behind `k`.
fn at_node<'a>(
    sum: &Totals,
    branches: impl Iterator<Item = &'a Behind> + Clone,
    support: bool,
) -> f64 {
    let squares: [f64; 3] =
        std::array::from_fn(|r| branches.clone().map(|k| k.e[r] * k.e[r]).sum());
    let mut total = 0.0;
    for k in branches {
        for (r, square) in squares.iter().enumerate() {
            let others = sum.e[r] - k.e[r];
            let pairs = (others * others - (square - k.e[r] * k.e[r])) / 2.0;
            let behind = k.pairs[r] - if support { k.missed[r] } else { 0.0 };
            total += pairs * behind;
        }
    }
    total
}

/// The weighted quartets of all the gene trees for the node with clusters `a`, `b` and the rest: each quartet is seen from the
/// two nodes where its pairs meet and counted from both, as in `weight::tripartition`.
pub fn tripartition(trees: &[WTree], a: &Bits, b: &Bits, weights: Weights) -> f64 {
    let most = trees.iter().map(|t| t.nodes.len()).max().unwrap_or(0);
    let mut down = vec![Behind::default(); most];
    let mut up = vec![Behind::default(); most];
    trees
        .iter()
        .map(|t| gene(t, a, b, weights.uses_support(), &mut down, &mut up))
        .sum()
}
