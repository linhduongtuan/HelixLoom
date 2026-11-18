//! The greedy consensus of the species trees of a bootstrap (what ASTRAL writes after the replicate trees): the splits are taken in
//! the order of how many of the trees have them, the most first, and each is kept if it fits in one tree with all that were kept
//! before. It has every split that more than half of the trees have, and as many others as can be added; a branch is labelled with
//! the share of the trees (in percent) that have it. A split is the set of taxa below a branch, on the side without the root taxon.

use super::bits::Bits;
use super::report::{name, support_text};
use std::collections::HashMap;

/// Two splits fit in one tree when one is inside the other or they share no taxon.
fn compatible(a: &Bits, b: &Bits) -> bool {
    a.is_disjoint(b) || a.is_subset_of(b) || b.is_subset_of(a)
}

/// The splits to keep, with the number of trees that have each, from the splits of each tree. Ties are settled by the splits'
/// taxa in order, so the answer does not depend on the order of the trees.
pub fn greedy(trees: &[Vec<Bits>]) -> Vec<(Bits, usize)> {
    let mut count: HashMap<&Bits, usize> = HashMap::new();
    for split in trees.iter().flatten() {
        *count.entry(split).or_insert(0) += 1;
    }
    let mut all: Vec<(Vec<usize>, &Bits, usize)> = count
        .into_iter()
        .map(|(split, c)| (split.iter().collect(), split, c))
        .collect();
    all.sort_by(|a, b| b.2.cmp(&a.2).then_with(|| a.0.cmp(&b.0)));
    let mut kept: Vec<(Bits, usize)> = Vec::new();
    for (_, split, c) in all {
        if kept.iter().all(|(k, _)| compatible(split, k)) {
            kept.push((split.clone(), c));
        }
    }
    kept
}

/// What a node of the consensus has below it.
enum Child {
    Node(usize),
    Leaf(usize),
}

/// The consensus as Newick, hung from the taxon `root`: the kept splits as nodes, each labelled with the share of the `total` trees
/// that have it, and the rest as polytomies. The children of a node come in the order of their lowest taxon.
pub fn newick(chosen: &[(Bits, usize)], total: usize, taxa: &[String], root: usize) -> String {
    let n = taxa.len();
    // node 0 is everything but the root taxon; the kept splits follow, the biggest first, so that the nearest one before a
    // node that holds it is the smallest that does
    let mut order: Vec<&(Bits, usize)> = chosen.iter().collect();
    order.sort_by_key(|(b, _)| std::cmp::Reverse(b.count()));
    let mut sets: Vec<Bits> = vec![Bits::of(n, (0..n).filter(|&t| t != root))];
    sets.extend(order.iter().map(|(b, _)| b.clone()));
    let mut share: Vec<usize> = vec![0];
    share.extend(order.iter().map(|(_, c)| *c));
    let mut kids: Vec<Vec<(usize, Child)>> = (0..sets.len()).map(|_| Vec::new()).collect();
    for i in 1..sets.len() {
        let parent = (0..i)
            .rev()
            .find(|&j| sets[i].is_subset_of(&sets[j]))
            .unwrap_or(0);
        let lowest = sets[i].first().unwrap_or(0);
        kids[parent].push((lowest, Child::Node(i)));
    }
    for t in (0..n).filter(|&t| t != root) {
        let parent = (0..sets.len())
            .rev()
            .find(|&j| sets[j].contains(t))
            .unwrap_or(0);
        kids[parent].push((t, Child::Leaf(t)));
    }
    for list in &mut kids {
        list.sort_by_key(|(lowest, _)| *lowest);
    }
    // written without recursion, as the species tree is: a tree can be as deep as it is wide
    enum Step {
        Node(usize),
        Comma,
        Close(usize),
        Leaf(usize),
    }
    let mut out = format!("({},", name(&taxa[root]));
    let mut stack = vec![Step::Node(0)];
    while let Some(step) = stack.pop() {
        match step {
            Step::Comma => out.push(','),
            Step::Leaf(t) => out.push_str(&name(&taxa[t])),
            Step::Close(i) => {
                out.push(')');
                if i > 0 && total > 0 {
                    out.push_str(&support_text(100.0 * share[i] as f64 / total as f64));
                }
            }
            Step::Node(i) => {
                out.push('(');
                stack.push(Step::Close(i));
                let last = kids[i].len().saturating_sub(1);
                for (k, (_, child)) in kids[i].iter().enumerate().rev() {
                    stack.push(match child {
                        Child::Node(j) => Step::Node(*j),
                        Child::Leaf(t) => Step::Leaf(*t),
                    });
                    if k != 0 && k <= last {
                        stack.push(Step::Comma);
                    }
                }
            }
        }
    }
    out.push_str(");");
    out
}
