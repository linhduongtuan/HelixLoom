//! The ASTRAL-Pro score of one node of a species tree against gene family trees (`family.rs`): the quartet classes of the
//! family trees that its three clusters `P1 P2 P3` agree with.
//!
//! A quartet of a family tree is *speciation-driven* when its four leaves are of four different species and every three of
//! them have a speciation for last common ancestor; two such quartets are one *class* when they have the same four species
//! and the same node where the pairs of the quartet meet (the last common ancestor of the quartet's two anchors), and the
//! score of a species tree is the number of classes whose topology it has. At a speciation `w` with species `M1`, `M2` below
//! its two branches and `M3` outside (disjoint sets), the classes that meet at `w` are of two kinds:
//!
//! * two species of `M1` against two of `M2` (`ab|cd`, `a b` in `M1`);
//! * two species of one branch, one of the other, and one of `M3`: (`aa'|bz`).
//!
//! Each class is seen from two nodes of the species tree (where its two pairs are), and is counted at one of them: where
//! the pair of the two species of the same branch of `w` has its species together in a cluster and the other two are in the
//! other two. With `I[m][j]` the species of `M(m)` in cluster `Pj`, that is
//!
//! ```text
//! Σ_i C(I[1][i],2) · I[2][j] · I[2][k]                      (j < k: the other two clusters)
//! + Σ over (i,j,k) a permutation of 1 2 3 of  C(I[1][i],2) · I[2][j] · I[3][k]  +  C(I[2][j],2) · I[1][i] · I[3][k]
//! ```
//!
//! and the score of a species tree is the sum of this over its nodes and the speciations of the family trees, with no half.

use super::bits::Bits;
use super::family::ProNode;

const PERMUTATIONS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

fn c2(n: u64) -> u64 {
    n * n.saturating_sub(1) / 2
}

/// The classes of one speciation for the clusters `a`, `b` and the rest.
fn node(p: &ProNode, a: &Bits, b: &Bits) -> u64 {
    let counts: [[u64; 3]; 3] = std::array::from_fn(|m| {
        let (x, y) = (
            u64::from(p.m[m].and_count(a)),
            u64::from(p.m[m].and_count(b)),
        );
        [x, y, u64::from(p.m[m].count()) - x - y]
    });
    let [m1, m2, m3] = counts;
    // two species of the first branch against two of the second
    let mut total: u64 = (0..3)
        .map(|i| {
            let (j, k) = match i {
                0 => (1, 2),
                1 => (0, 2),
                _ => (0, 1),
            };
            c2(m1[i]) * m2[j] * m2[k]
        })
        .sum();
    for [i, j, k] in PERMUTATIONS {
        total += c2(m1[i]) * m2[j] * m3[k] + c2(m2[j]) * m1[i] * m3[k];
    }
    total
}

/// The classes of all the family trees for the node with clusters `a`, `b` and the rest.
pub fn tripartition(nodes: &[ProNode], a: &Bits, b: &Bits) -> u64 {
    nodes.iter().map(|p| node(p, a, b)).sum()
}
