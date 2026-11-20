//! Editing a tree as FigTree and TreeGraph 2 do, written from their documentation (no source of either was read): putting the root on a
//! branch, rotating and ordering the branches, collapsing the weakly supported ones, cutting out a clade, renaming tips, and carrying the
//! support values of other trees over to a reference tree.
//!
//! Annotations move with the branch they describe. A branch is stored on the node below it (its length, and for an inner node its label
//! and its `[&...]` annotations), so a reversed branch hands them to the node that is now below it.

use super::newick::{number, NwNode, NwTree};
use std::collections::HashMap;

/// A tree rebuilt from the root down, so that nodes that were left out are gone and every parent comes before its children.
fn compact(t: &NwTree, root: usize) -> NwTree {
    let mut nodes: Vec<NwNode> = Vec::new();
    let mut stack = vec![(root, None::<usize>)];
    while let Some((i, parent)) = stack.pop() {
        let me = nodes.len();
        let mut n = t.nodes[i].clone();
        (n.parent, n.children) = (parent, Vec::new());
        nodes.push(n);
        if let Some(p) = parent {
            nodes[p].children.push(me);
        }
        stack.extend(t.nodes[i].children.iter().rev().map(|&c| (c, Some(me))));
    }
    NwTree { nodes, root: 0 }
}

fn length_of(t: &NwTree, i: usize) -> f64 {
    let d = t.nodes[i].dist;
    if d.is_finite() {
        d
    } else {
        0.0
    }
}

/// The tree rooted on the branch above `node`, `fraction` of the way along it from `node` (a half is the middle), as FigTree's and
/// TreeGraph's "root here" does. The branch's length is shared between the two new branches that meet at the root, and its label and
/// annotations are on both. Every branch between the old root and `node` is reversed, with what is written on it; an old root that is left
/// with one child is dissolved.
pub fn reroot_above(t: &NwTree, node: usize, fraction: f64) -> Result<NwTree, String> {
    if node >= t.nodes.len() {
        return Err("There is no such node.".into());
    }
    if node == t.root {
        return Err("The root has no branch above it to put the root on.".into());
    }
    let f = if fraction.is_finite() {
        fraction.clamp(0.0, 1.0)
    } else {
        0.5
    };
    let mut path = vec![node];
    while let Some(p) = t.nodes[*path.last().unwrap_or(&node)].parent {
        path.push(p);
    }
    path.reverse();
    let k = path.len() - 1;
    let mut nodes = t.nodes.clone();
    let new_root = nodes.len();
    nodes.push(NwNode {
        name: String::new(),
        dist: 0.0,
        parent: None,
        children: vec![node, path[k - 1]],
        attrs: Vec::new(),
    });
    // the branch above `node`: both halves carry its label (not a tip's name) and annotations
    let whole = length_of(t, node);
    let above = &mut nodes[path[k - 1]];
    above.name = if t.is_leaf(node) {
        String::new()
    } else {
        t.nodes[node].name.clone()
    };
    above.attrs = t.nodes[node].attrs.clone();
    (above.dist, nodes[node].dist) = ((1.0 - f) * whole, f * whole);
    // each branch on the way down is now above the node that was above it
    for i in 1..k {
        let from = &t.nodes[path[i]];
        let to = &mut nodes[path[i - 1]];
        (to.name, to.attrs, to.dist) = (from.name.clone(), from.attrs.clone(), from.dist);
    }
    for i in 0..k {
        let next = path[i + 1];
        let mut kids: Vec<usize> = t.nodes[path[i]]
            .children
            .iter()
            .copied()
            .filter(|&c| c != next)
            .collect();
        if i > 0 {
            kids.push(path[i - 1]);
        }
        nodes[path[i]].children = kids;
    }
    let mut tree = NwTree {
        nodes,
        root: new_root,
    };
    dissolve_old_root(&mut tree, t, path[0], k, f, whole);
    Ok(compact(&tree, tree.root))
}

/// The old root had two children, one of them on the path: it is now a node with one child, which is no node at all. Its child takes its
/// place, and the lengths are made one.
fn dissolve_old_root(
    tree: &mut NwTree,
    old: &NwTree,
    was_root: usize,
    k: usize,
    f: f64,
    whole: f64,
) {
    if tree.nodes[was_root].children.len() != 1 {
        return;
    }
    let only = tree.nodes[was_root].children[0];
    if k == 1 {
        // the two branches of the old root were one branch: its whole length is shared out again at the new root
        let total = whole + length_of(old, only);
        let top = tree.root;
        let kids: Vec<usize> = tree.nodes[top]
            .children
            .iter()
            .map(|&c| if c == was_root { only } else { c })
            .collect();
        tree.nodes[top].children = kids;
        let node = tree.nodes[top].children[0];
        tree.nodes[node].dist = f * total;
        tree.nodes[only].dist = (1.0 - f) * total;
        return;
    }
    let carried = tree.nodes[was_root].dist;
    let below = old_parent_after_reversal(tree, was_root);
    if let Some(p) = below {
        let kids: Vec<usize> = tree.nodes[p]
            .children
            .iter()
            .map(|&c| if c == was_root { only } else { c })
            .collect();
        tree.nodes[p].children = kids;
    }
    tree.nodes[only].dist = length_of(tree, only) + if carried.is_finite() { carried } else { 0.0 };
}

/// The node that has `child` among its children.
fn old_parent_after_reversal(t: &NwTree, child: usize) -> Option<usize> {
    (0..t.nodes.len()).find(|&i| t.nodes[i].children.contains(&child))
}

/// The branches of `node` in the opposite order (FigTree's "rotate").
pub fn rotate(t: &NwTree, node: usize) -> NwTree {
    let mut out = t.clone();
    if let Some(n) = out.nodes.get_mut(node) {
        n.children.reverse();
    }
    out
}

/// Every node's branches ordered by the number of tips they lead to, the smaller first when `increasing`, else the larger first
/// (FigTree's "order nodes"); branches with as many tips keep their order.
pub fn order(t: &NwTree, increasing: bool) -> NwTree {
    let mut out = t.clone();
    let mut tips = vec![0usize; out.nodes.len()];
    for i in out.postorder() {
        tips[i] = if out.is_leaf(i) {
            1
        } else {
            out.nodes[i].children.iter().map(|&c| tips[c]).sum()
        };
    }
    for n in &mut out.nodes {
        if increasing {
            n.children.sort_by_key(|&c| tips[c]);
        } else {
            n.children.sort_by_key(|&c| std::cmp::Reverse(tips[c]));
        }
    }
    out
}

/// The value that decides whether an inner node's branch is supported: the annotation `key`, or the label when `key` is empty.
pub fn support_value(t: &NwTree, i: usize, key: &str) -> Option<f64> {
    if key.is_empty() {
        t.nodes[i]
            .name
            .trim()
            .parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
    } else {
        t.attr_number(i, key)
    }
}

/// The tree with each inner branch whose value (`support_value`) is below `least` taken out, its branches going to its parent as TreeGraph
/// 2's "collapse nodes by support" does, with its length added to theirs; and how many were. The root and the tips stay, and a branch with no
/// value stays.
pub fn collapse_below(t: &NwTree, key: &str, least: f64) -> (NwTree, usize) {
    let mut out = t.clone();
    let mut gone = 0;
    for i in t.preorder() {
        let Some(parent) = out.nodes[i].parent.filter(|_| !out.is_leaf(i)) else {
            continue;
        };
        if support_value(t, i, key).is_none_or(|v| v >= least) {
            continue;
        }
        let carried = length_of(&out, i);
        let kids = std::mem::take(&mut out.nodes[i].children);
        for &c in &kids {
            out.nodes[c].parent = Some(parent);
            out.nodes[c].dist = length_of(&out, c) + carried;
        }
        let at = out.nodes[parent]
            .children
            .iter()
            .position(|&c| c == i)
            .unwrap_or(0);
        out.nodes[parent].children.splice(at..=at, kids);
        gone += 1;
    }
    (compact(&out, out.root), gone)
}

/// The clade under `node` as a tree of its own.
pub fn extract_clade(t: &NwTree, node: usize) -> Result<NwTree, String> {
    if node >= t.nodes.len() {
        return Err("There is no such node.".into());
    }
    let mut out = compact(t, node);
    out.nodes[0].dist = 0.0;
    Ok(out)
}

/// The tips renamed by a table of `(old, new)` names; how many were. A name not in the table stays.
pub fn rename_tips(t: &NwTree, table: &[(String, String)]) -> (NwTree, usize) {
    let by: HashMap<&str, &str> = table.iter().map(|(a, b)| (a.trim(), b.trim())).collect();
    let mut out = t.clone();
    let mut done = 0;
    for n in out.nodes.iter_mut().filter(|n| n.children.is_empty()) {
        if let Some(new) = by.get(n.name.trim()).filter(|new| !new.is_empty()) {
            n.name = (*new).to_string();
            done += 1;
        }
    }
    (out, done)
}

/// A table of names from text: one pair to a line, the two separated by a tab, a comma or a semicolon; blank lines and lines of one
/// name are left out.
pub fn parse_name_table(text: &str) -> Vec<(String, String)> {
    text.lines()
        .filter_map(|l| {
            let (a, b) = l.split_once(['\t', ',', ';'])?;
            let (a, b) = (a.trim().trim_matches('"'), b.trim().trim_matches('"'));
            (!a.is_empty() && !b.is_empty()).then(|| (a.to_string(), b.to_string()))
        })
        .collect()
}

/// A split with a bit set for each of `total` tips.
fn full_mask(total: usize) -> Vec<u64> {
    let words = total.div_ceil(64).max(1);
    let mut full = vec![u64::MAX; words];
    if !total.is_multiple_of(64) {
        full[words - 1] = (1u64 << (total % 64)) - 1;
    }
    full
}

/// A branch of a tree as the set of tips on one side of it: the side that does not hold the first common tip, as a row of bits.
type Split = Vec<u64>;

/// For every node of `t`, its split over the tips that `index` knows (by name); none for a node that is trivial there (fewer than two
/// tips on a side). Nodes of a rooted tree's two root branches have the same split.
fn splits(t: &NwTree, index: &HashMap<&str, usize>) -> Vec<Option<Split>> {
    let words = index.len().div_ceil(64).max(1);
    let mut below: Vec<Split> = vec![vec![0; words]; t.nodes.len()];
    for i in t.postorder() {
        if t.is_leaf(i) {
            if let Some(&b) = index.get(t.nodes[i].name.as_str()) {
                below[i][b / 64] |= 1 << (b % 64);
            }
        }
        if let Some(p) = t.nodes[i].parent {
            let kid = below[i].clone();
            below[p].iter_mut().zip(kid).for_each(|(a, b)| *a |= b);
        }
    }
    let total = index.len();
    (0..t.nodes.len())
        .map(|i| {
            let mut s = below[i].clone();
            if s[0] & 1 == 1 {
                // the side with the first tip is never the one that is kept
                let full = full_mask(total);
                s.iter_mut().zip(full).for_each(|(a, b)| *a = !*a & b);
            }
            let ones: u32 = s.iter().map(|w| w.count_ones()).sum();
            (t.nodes[i].parent.is_some() && ones >= 2 && total - ones as usize >= 2).then_some(s)
        })
        .collect()
}

fn names_of(t: &NwTree) -> Vec<&str> {
    t.preorder()
        .into_iter()
        .filter(|&i| t.is_leaf(i))
        .map(|i| t.nodes[i].name.as_str())
        .collect()
}

/// The tips that every one of the trees has, in the reference's order.
fn shared_tips<'a>(
    reference: &'a NwTree,
    others: &[NwTree],
) -> Result<HashMap<&'a str, usize>, String> {
    let mut common: Vec<&str> = names_of(reference);
    for o in others {
        let theirs = names_of(o);
        common.retain(|n| theirs.contains(n));
    }
    if common.len() < 4 {
        return Err(format!(
            "The trees share {} tips: at least four are needed to compare their branches.",
            common.len()
        ));
    }
    Ok(common
        .into_iter()
        .enumerate()
        .map(|(i, n)| (n, i))
        .collect())
}

/// Whether two splits can be branches of one tree: one of the four sets they cut out is empty.
fn compatible(a: &Split, b: &Split, total: usize) -> bool {
    let words = a.len();
    let full = full_mask(total);
    let empty = |f: &dyn Fn(u64, u64, u64) -> u64| (0..words).all(|w| f(a[w], b[w], full[w]) == 0);
    empty(&|x, y, _| x & y)
        || empty(&|x, y, _| x & !y)
        || empty(&|x, y, _| !x & y)
        || empty(&|x, y, f| !x & !y & f)
}

/// What a set of trees says of each branch of a reference tree.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Support {
    /// The share (0 to 1) of the trees that have the branch.
    pub frequency: f64,
    /// The share of the trees that have the branch that conflicts with it most often: one that cannot be in a tree with it.
    pub conflict: f64,
}

/// The most distinct branches that are searched for a conflict with a branch.
const MOST_RIVALS: usize = 5_000;

/// How often each branch of `reference` (matched by the tips on its two sides, the root's place ignored) is among the branches of `others`,
/// say bootstrap replicates, and how often its strongest rival is; none for a branch that splits off fewer than two of the tips the trees
/// share. Tips that not every tree has are left out of the comparison.
pub fn frequencies(reference: &NwTree, others: &[NwTree]) -> Result<Vec<Option<Support>>, String> {
    if others.is_empty() {
        return Err("There are no other trees to take the support from.".into());
    }
    let index = shared_tips(reference, others)?;
    let total = index.len();
    let mut count: HashMap<Split, usize> = HashMap::new();
    for o in others {
        let mut seen: Vec<Split> = splits(o, &index).into_iter().flatten().collect();
        seen.sort();
        seen.dedup();
        for s in seen {
            *count.entry(s).or_insert(0) += 1;
        }
    }
    let n = others.len() as f64;
    let mut rivals: Vec<(&Split, usize)> = count.iter().map(|(s, &c)| (s, c)).collect();
    rivals.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
    rivals.truncate(MOST_RIVALS);
    Ok(splits(reference, &index)
        .into_iter()
        .map(|split| {
            let split = split?;
            let frequency = count.get(&split).copied().unwrap_or(0) as f64 / n;
            let conflict = rivals
                .iter()
                .find(|(r, _)| !compatible(&split, r, total))
                .map_or(0.0, |&(_, c)| c as f64 / n);
            Some(Support {
                frequency,
                conflict,
            })
        })
        .collect())
}

/// For each branch of `reference`, the value written on the branch of `other` that cuts off the same tips (its annotation `key`, or its
/// label when `key` is empty): the way TreeGraph 2 imports the posterior probabilities of another analysis. None where `other` has no such
/// branch or no value on it.
pub fn values_from(
    reference: &NwTree,
    other: &NwTree,
    key: &str,
) -> Result<Vec<Option<f64>>, String> {
    let index = shared_tips(reference, std::slice::from_ref(other))?;
    let theirs: HashMap<Split, usize> = splits(other, &index)
        .into_iter()
        .enumerate()
        .filter_map(|(i, s)| Some((s?, i)))
        .collect();
    Ok(splits(reference, &index)
        .into_iter()
        .map(|s| support_value(other, *theirs.get(&s?)?, key))
        .collect())
}

/// The tree with `values` written on its nodes as the annotation `key` (a node with none gets nothing). A key that is there is replaced.
pub fn with_values(t: &NwTree, key: &str, values: &[Option<f64>]) -> NwTree {
    let mut out = t.clone();
    for (i, v) in values.iter().enumerate().take(out.nodes.len()) {
        out.nodes[i].attrs.retain(|(k, _)| k != key || v.is_none());
        if let Some(v) = v {
            out.nodes[i].attrs.retain(|(k, _)| k != key);
            out.nodes[i].attrs.push((key.to_string(), number(*v)));
        }
    }
    out
}

#[cfg(test)]
mod tests;
