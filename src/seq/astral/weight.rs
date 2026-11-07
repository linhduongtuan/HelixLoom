//! Counting quartets of gene trees against the nodes and branches of a species tree.
//!
//! A quartet `ab|cd` of a tree has two nodes that join a pair: one where `a` and `b` meet, with `c` and `d` together on the
//! third branch, and one where `c` and `d` meet. So every count below looks at the node that joins the pair it names and
//! asks only that the other two taxa are together on one branch: each quartet is then counted once.

use super::bits::Bits;
use super::genes::{Branch, Gene};

/// The quartets of the gene trees that agree with the node of the species tree that has the clusters `a` and `b` below it
/// and everything else (`c`) above. Counted are the quartets with two taxa in one of the three clusters and one in each of
/// the others, resolved by the gene tree as that pair against the other two. A quartet of the species tree is seen from two
/// of its nodes, so the quartet score of a tree is half the sum of this over its inner nodes.
pub fn tripartition(genes: &[Gene], a: &Bits, b: &Bits) -> u64 {
    let mut scratch = Vec::with_capacity(8);
    genes
        .iter()
        .flat_map(|g| g.nodes.iter())
        .map(|branches| node_weight(branches, a, b, &mut scratch))
        .sum()
}

/// Taxa of each of the three clusters on one branch.
type Counts = (u64, u64, u64);

fn node_weight(branches: &[Branch], a: &Bits, b: &Bits, scratch: &mut Vec<Counts>) -> u64 {
    scratch.clear();
    let (mut ta, mut tb, mut tc) = (0u64, 0u64, 0u64);
    let (mut qa, mut qb, mut qc) = (0u64, 0u64, 0u64);
    for br in branches {
        let x = u64::from(a.and_count(&br.set));
        let y = u64::from(b.and_count(&br.set));
        let z = br.size - x - y;
        scratch.push((x, y, z));
        (ta, tb, tc) = (ta + x, tb + y, tc + z);
        (qa, qb, qc) = (qa + x * x, qb + y * y, qc + z * z);
    }
    // two taxa of a cluster on two different branches other than k: the sum of s_i s_j over those pairs
    let across =
        |total: u64, squares: u64, k: u64| ((total - k) * (total - k) - (squares - k * k)) / 2;
    scratch
        .iter()
        .map(|&(x, y, z)| {
            across(ta, qa, x) * y * z + across(tb, qb, y) * x * z + across(tc, qc, z) * x * y
        })
        .sum()
}

/// The quartets around one inner branch of the species tree, with the four clusters on its ends: `l` and `r` on one end,
/// `s` and `o` on the other. Of the three ways to pair them `lr|so` is the branch's own; `rs|lo` and `ro|ls` are the others.
pub struct Around {
    /// A gene tree contributes the share of its quartets around the branch that agree with each pairing; the three add up to
    /// 1 for a tree that resolves them all, less when a polytomy leaves some unresolved, and 0 for one with none around.
    pub f: [f64; 3],
    /// The sum of the three: the effective number of genes.
    pub en: f64,
    /// How many quartets the four clusters make, `|l| |r| |s| |o|`.
    pub quartets: u64,
}

pub fn around_branch(genes: &[Gene], l: &Bits, r: &Bits, s: &Bits, o: &Bits) -> Around {
    let mut f = [0.0; 3];
    let mut scratch: Vec<[u64; 4]> = Vec::with_capacity(8);
    for g in genes {
        let m = [l, r, s, o].map(|c| u64::from(c.and_count(&g.present)));
        let possible = m[0] * m[1] * m[2] * m[3];
        if possible == 0 {
            continue;
        }
        let n = gene_counts(g, [l, r, s, o], &mut scratch);
        for (total, found) in f.iter_mut().zip(n) {
            *total += found as f64 / possible as f64;
        }
    }
    Around {
        en: f.iter().sum(),
        f,
        quartets: [l, r, s, o].iter().map(|c| u64::from(c.count())).product(),
    }
}

/// For one gene tree, the quartets (one taxon from each of `l r s o`) it resolves as `lr|so`, `rs|lo` and `ro|ls`.
fn gene_counts(g: &Gene, sets: [&Bits; 4], scratch: &mut Vec<[u64; 4]>) -> [u64; 3] {
    let mut n = [0u64; 3];
    for branches in &g.nodes {
        scratch.clear();
        let mut total = [0u64; 4];
        for br in branches {
            let c = sets.map(|s| u64::from(s.and_count(&br.set)));
            for (t, x) in total.iter_mut().zip(c) {
                *t += x;
            }
            scratch.push(c);
        }
        let paired = |i: usize, j: usize| scratch.iter().map(|c| c[i] * c[j]).sum::<u64>();
        // ordered pairs of taxa of two clusters on different branches that are not `k`
        let across = |i: usize, j: usize, k: &[u64; 4]| {
            (total[i] - k[i]) * (total[j] - k[j]) - (paired(i, j) - k[i] * k[j])
        };
        for k in scratch.iter() {
            n[0] += k[2] * k[3] * across(0, 1, k);
            n[1] += k[0] * k[3] * across(1, 2, k);
            n[2] += k[0] * k[2] * across(1, 3, k);
        }
    }
    n
}
