//! Tree tools for Newick trees, the jobs ETE and DendroPy are most used for: midpoint rooting,
//! ladderizing, patristic distances and the Robinson–Foulds distance between two trees.
//! Every function returns a new tree or a table and leaves its input alone.

use super::newick::{number, NwNode, NwTree};
use std::collections::{BTreeMap, BTreeSet};

const NEEDS_TWO_TIPS: &str = "The tree needs at least two tips.";

/// Counts and sizes of one tree.
#[derive(Clone, Debug, PartialEq)]
pub struct Summary {
    pub tips: usize,
    pub inner: usize,
    /// Sum of all branch lengths.
    pub length: f64,
    /// Longest distance from the root to a tip.
    pub height: f64,
    /// Every inner node has exactly two children.
    pub binary: bool,
}

pub fn summary(t: &NwTree) -> Summary {
    let order = t.levelorder();
    let from_root = root_distances(t);
    let tips = order.iter().filter(|&&i| t.is_leaf(i)).count();
    Summary {
        tips,
        inner: order.len() - tips,
        length: order
            .iter()
            .filter(|&&i| i != t.root)
            .map(|&i| t.nodes[i].dist)
            .sum(),
        height: order
            .iter()
            .filter(|&&i| t.is_leaf(i))
            .map(|&i| from_root[i])
            .fold(0.0, f64::max),
        binary: order
            .iter()
            .all(|&i| t.is_leaf(i) || t.nodes[i].children.len() == 2),
    }
}

/// Distance from the root to every node under it (0 elsewhere).
fn root_distances(t: &NwTree) -> Vec<f64> {
    let mut d = vec![0.0; t.nodes.len()];
    for i in t.levelorder() {
        if let Some(p) = t.nodes[i].parent {
            d[i] = d[p] + t.nodes[i].dist;
        }
    }
    d
}

/// The tree without a root edge: a root with a single child is replaced by that child.
fn strip_root_edge(t: &NwTree) -> NwTree {
    let mut out = t.clone();
    while let [child] = out.nodes[out.root].children[..] {
        out.nodes[child].parent = None;
        out.root = child;
    }
    out
}

/// Tips left to right as the tree is drawn.
fn tips_in_order(t: &NwTree) -> Vec<usize> {
    let mut out = Vec::new();
    let mut stack = vec![t.root];
    while let Some(i) = stack.pop() {
        if t.is_leaf(i) {
            out.push(i);
        } else {
            stack.extend(t.nodes[i].children.iter().rev());
        }
    }
    out
}

/// Tip names in tree order and the length of the path between every two tips.
// ponytail: climbs from each tip to the common ancestor, O(tips² × depth); fine for the hundreds
// of tips the windows take, a lowest-common-ancestor table if thousands are ever needed.
pub fn patristic(t: &NwTree) -> (Vec<String>, Vec<Vec<f64>>) {
    let from_root = root_distances(t);
    let mut level = vec![0usize; t.nodes.len()];
    for i in t.levelorder() {
        if let Some(p) = t.nodes[i].parent {
            level[i] = level[p] + 1;
        }
    }
    let common_ancestor = |mut a: usize, mut b: usize| {
        while a != b {
            let lower = if level[a] >= level[b] { &mut a } else { &mut b };
            match t.nodes[*lower].parent {
                Some(p) => *lower = p,
                None => break,
            }
        }
        a
    };
    let tips = tips_in_order(t);
    let mut table = vec![vec![0.0; tips.len()]; tips.len()];
    for (i, &a) in tips.iter().enumerate() {
        for (j, &b) in tips.iter().enumerate().skip(i + 1) {
            let d = from_root[a] + from_root[b] - 2.0 * from_root[common_ancestor(a, b)];
            table[i][j] = d;
            table[j][i] = d;
        }
    }
    (
        tips.iter().map(|&i| t.nodes[i].name.clone()).collect(),
        table,
    )
}

/// The patristic distances as tab-separated text, ready to paste into a spreadsheet.
pub fn distance_table(t: &NwTree) -> String {
    let (names, table) = patristic(t);
    let mut out = format!("\t{}\n", names.join("\t"));
    for (name, row) in names.iter().zip(&table) {
        let cells: Vec<String> = row.iter().map(|&d| number(d)).collect();
        out.push_str(&format!("{name}\t{}\n", cells.join("\t")));
    }
    out
}

/// Smaller clades first at every node (ETE's `ladderize`); clades of equal size keep their order.
pub fn ladderize(t: &NwTree) -> NwTree {
    let mut out = t.clone();
    let mut tips = vec![0usize; out.nodes.len()];
    for i in out.postorder() {
        tips[i] = if out.is_leaf(i) {
            1
        } else {
            out.nodes[i].children.iter().map(|&c| tips[c]).sum()
        };
    }
    for node in &mut out.nodes {
        node.children.sort_by_key(|&c| tips[c]);
    }
    out
}

struct Edge {
    a: usize,
    b: usize,
    len: f64,
    /// What the inner node below this edge was called (often a support value).
    label: String,
}

/// A tree without a root: nodes joined by edges that have lengths.
struct Unrooted {
    names: Vec<String>,
    tip: Vec<bool>,
    edges: Vec<Edge>,
    /// Edges at each node.
    at: Vec<Vec<usize>>,
}

impl Unrooted {
    /// A root with two children is not a node of the unrooted tree: its two branches become one edge.
    fn of(t: &NwTree) -> Unrooted {
        let order = t.levelorder();
        let drop_root = t.nodes[t.root].children.len() == 2;
        let mut index = vec![usize::MAX; t.nodes.len()];
        let (mut names, mut tip) = (Vec::new(), Vec::new());
        for &i in order.iter().filter(|&&i| !(drop_root && i == t.root)) {
            index[i] = names.len();
            names.push(t.nodes[i].name.clone());
            tip.push(t.is_leaf(i));
        }
        let label = |i: usize| {
            if t.is_leaf(i) {
                String::new()
            } else {
                t.nodes[i].name.clone()
            }
        };
        let mut edges: Vec<Edge> = order
            .iter()
            .filter_map(|&i| {
                let p = t.nodes[i].parent.filter(|&p| !(drop_root && p == t.root))?;
                Some(Edge {
                    a: index[i],
                    b: index[p],
                    len: t.nodes[i].dist,
                    label: label(i),
                })
            })
            .collect();
        if let [x, y] = t.nodes[t.root].children[..] {
            let named = if label(x).is_empty() { y } else { x };
            edges.push(Edge {
                a: index[x],
                b: index[y],
                len: t.nodes[x].dist + t.nodes[y].dist,
                label: label(named),
            });
        }
        let mut at = vec![Vec::new(); names.len()];
        for (e, edge) in edges.iter().enumerate() {
            at[edge.a].push(e);
            at[edge.b].push(e);
        }
        Unrooted {
            names,
            tip,
            edges,
            at,
        }
    }

    fn other_end(&self, e: usize, from: usize) -> usize {
        let edge = &self.edges[e];
        if edge.a == from {
            edge.b
        } else {
            edge.a
        }
    }

    /// Distance from `start` to every node, and the edge each node was reached by.
    fn walk(&self, start: usize) -> (Vec<f64>, Vec<Option<usize>>) {
        let mut dist = vec![f64::NAN; self.names.len()];
        let mut via = vec![None; self.names.len()];
        dist[start] = 0.0;
        let mut stack = vec![start];
        while let Some(g) = stack.pop() {
            for &e in &self.at[g] {
                let o = self.other_end(e, g);
                if dist[o].is_nan() {
                    dist[o] = dist[g] + self.edges[e].len;
                    via[o] = Some(e);
                    stack.push(o);
                }
            }
        }
        (dist, via)
    }

    /// The tip farthest away (the first of any that tie).
    fn farthest_tip(&self, dist: &[f64]) -> Option<usize> {
        (0..self.names.len())
            .filter(|&g| self.tip[g] && !dist[g].is_nan())
            .fold(None, |best, g| match best {
                Some(b) if dist[b] >= dist[g] => Some(b),
                _ => Some(g),
            })
    }

    /// Root the tree on edge `e`, `at` from its first end and the rest from its second.
    fn rooted_on(&self, e: usize, at: f64) -> NwTree {
        let mut nodes = vec![NwNode {
            name: String::new(),
            dist: 0.0,
            parent: None,
            children: Vec::new(),
        }];
        let edge = &self.edges[e];
        // (node, edge it hangs from, its parent in the new tree, length of that edge)
        let mut stack = vec![(edge.b, e, 0, edge.len - at), (edge.a, e, 0, at)];
        while let Some((g, via, parent, dist)) = stack.pop() {
            let me = nodes.len();
            nodes.push(NwNode {
                name: if self.tip[g] {
                    self.names[g].clone()
                } else {
                    self.edges[via].label.clone()
                },
                dist,
                parent: Some(parent),
                children: Vec::new(),
            });
            nodes[parent].children.push(me);
            for &next in self.at[g].iter().rev().filter(|&&n| n != via) {
                stack.push((self.other_end(next, g), next, me, self.edges[next].len));
            }
        }
        NwTree { nodes, root: 0 }
    }
}

/// The two tips farthest apart and how far.
#[derive(Clone, Debug, PartialEq)]
pub struct Midpoint {
    pub tip_a: String,
    pub tip_b: String,
    pub diameter: f64,
}

/// Root the tree halfway along the longest path between two tips. Branch labels stay with their
/// branches; the old root disappears if it only joined two branches.
pub fn midpoint_root(t: &NwTree) -> Result<(NwTree, Midpoint), String> {
    let g = Unrooted::of(&strip_root_edge(t));
    let first = (0..g.names.len())
        .find(|&i| g.tip[i])
        .ok_or(NEEDS_TWO_TIPS)?;
    // the farthest tip from any tip is one end of the longest path; the farthest from that, the other
    let a = g.farthest_tip(&g.walk(first).0).ok_or(NEEDS_TWO_TIPS)?;
    let (from_a, via) = g.walk(a);
    let b = g.farthest_tip(&from_a).ok_or(NEEDS_TWO_TIPS)?;
    let diameter = from_a[b];
    if a == b {
        return Err(NEEDS_TWO_TIPS.into());
    }
    if diameter <= 0.0 {
        return Err("Every branch has length zero, so there is no midpoint.".into());
    }
    let half = diameter / 2.0;
    // walk from b back towards a to the edge that holds the midpoint
    let mut node = b;
    while let Some(e) = via[node] {
        let toward_a = g.other_end(e, node);
        if from_a[toward_a] <= half {
            let from_there = (half - from_a[toward_a]).clamp(0.0, g.edges[e].len);
            let at = if g.edges[e].a == toward_a {
                from_there
            } else {
                g.edges[e].len - from_there
            };
            let found = Midpoint {
                tip_a: g.names[a].clone(),
                tip_b: g.names[b].clone(),
                diameter,
            };
            return Ok((g.rooted_on(e, at), found));
        }
        node = toward_a;
    }
    Err("Could not place the midpoint.".into())
}

/// How far apart two trees on the same tips are.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RobinsonFoulds {
    /// Splits found in only one of the two trees.
    pub distance: usize,
    /// The most the distance could be: every split of both trees.
    pub max: usize,
    /// Splits found in both.
    pub shared: usize,
}

fn tip_names(t: &NwTree) -> Result<Vec<String>, String> {
    let mut names: Vec<String> = tips_in_order(t)
        .iter()
        .map(|&i| t.nodes[i].name.clone())
        .collect();
    names.sort();
    if names.first().is_some_and(String::is_empty) {
        return Err("A tip has no name.".into());
    }
    if let Some(pair) = names.windows(2).find(|w| w[0] == w[1]) {
        return Err(format!("The tip name {} is used twice.", pair[0]));
    }
    Ok(names)
}

/// The splits of the tips a branch makes, each as the set of tips on the side without the first
/// name, so a tree's two branches at a root (one split) count once. Single-tip splits are left out.
fn splits(t: &NwTree, names: &[String]) -> BTreeSet<Vec<u64>> {
    let n = names.len();
    let words = n.div_ceil(64);
    let slot: BTreeMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let mut below: Vec<Vec<u64>> = vec![Vec::new(); t.nodes.len()];
    let mut found = BTreeSet::new();
    for i in t.postorder() {
        let mut set = vec![0u64; words];
        if t.is_leaf(i) {
            if let Some(&k) = slot.get(t.nodes[i].name.as_str()) {
                set[k / 64] |= 1 << (k % 64);
            }
        } else {
            for &c in &t.nodes[i].children {
                set.iter_mut().zip(&below[c]).for_each(|(w, x)| *w |= x);
            }
        }
        let size: usize = set.iter().map(|w| w.count_ones() as usize).sum();
        if i != t.root && size >= 2 && size + 2 <= n {
            let mut key = set.clone();
            if key[0] & 1 == 1 {
                key.iter_mut().for_each(|w| *w = !*w);
                if !n.is_multiple_of(64) {
                    key[words - 1] &= (1u64 << (n % 64)) - 1;
                }
            }
            found.insert(key);
        }
        below[i] = set;
    }
    found
}

/// The Robinson–Foulds distance between two trees on the same tips, ignoring where each is rooted.
pub fn robinson_foulds(a: &NwTree, b: &NwTree) -> Result<RobinsonFoulds, String> {
    let (names, other) = (tip_names(a)?, tip_names(b)?);
    if names != other {
        let only = |from: &[String], not_in: &[String]| -> Vec<String> {
            from.iter()
                .filter(|n| !not_in.contains(n))
                .take(3)
                .cloned()
                .collect()
        };
        return Err(format!(
            "The two trees need the same tips. Only in the first: {:?}. Only in the second: {:?}.",
            only(&names, &other),
            only(&other, &names)
        ));
    }
    let (sa, sb) = (splits(a, &names), splits(b, &names));
    let shared = sa.intersection(&sb).count();
    Ok(RobinsonFoulds {
        distance: sa.len() + sb.len() - 2 * shared,
        max: sa.len() + sb.len(),
        shared,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Answers from ETE3 (ladderize, Robinson–Foulds, distances) and DendroPy (midpoint rooting) for five
    /// random ten-tip trees: the tree, its ladderized form, its midpoint-rooted form, and the distances A–B,
    /// C–J, D–H and E–F.
    type Case = (&'static str, &'static str, &'static str, [f64; 4]);
    const CASES: [Case; 5] = [
        (
            "((((B:1.277,A:1.568):1.702,J:0.051):0.831,(I:0.459,H:1.825):1.701):0.187,((G:0.966,F:1.962):0.804,(E:0.825,(D:0.576,C:0.22):0.192):1.918):1.544);",
            "(((I:0.459,H:1.825):1.701,(J:0.051,(B:1.277,A:1.568):1.702):0.831):0.187,((G:0.966,F:1.962):0.804,(E:0.825,(D:0.576,C:0.22):0.192):1.918):1.544);",
            "(((G:0.966,F:1.962):0.804,(E:0.825,(D:0.576,C:0.22):0.192):1.918):1.5330000000000006,(((B:1.277,A:1.568):1.702,J:0.051):0.831,(I:0.459,H:1.825):1.701):0.19799999999999945);",
            [2.8449999999999998, 4.943, 7.943, 5.509],
        ),
        (
            "(((F:1.13,(E:0.263,D:1.049):1.39):1.24,(C:0.807,(B:0.548,A:1.493):0.991):1.465):0.656,(J:0.359,(I:0.203,(H:0.814,G:0.871):0.142):1.277):0.495);",
            "((J:0.359,(I:0.203,(H:0.814,G:0.871):0.142):1.277):0.495,((F:1.13,(E:0.263,D:1.049):1.39):1.24,(C:0.807,(B:0.548,A:1.493):0.991):1.465):0.656);",
            "((C:0.807,(B:0.548,A:1.493):0.991):1.3299999999999996,((F:1.13,(E:0.263,D:1.049):1.39):1.24,(J:0.359,(I:0.203,(H:0.814,G:0.871):0.142):1.277):1.151):0.13500000000000045);",
            [2.0410000000000004, 3.7820000000000005, 7.063, 2.7829999999999995],
        ),
        (
            "(((H:1.592,G:1.734):1.386,(F:1.167,(E:1.297,(D:1.766,C:1.47):1.987):1.269):1.807):1.834,((B:0.796,A:1.186):1.749,(J:1.237,I:0.206):0.863):0.736);",
            "(((B:0.796,A:1.186):1.749,(J:1.237,I:0.206):0.863):0.736,((H:1.592,G:1.734):1.386,(F:1.167,(E:1.297,(D:1.766,C:1.47):1.987):1.269):1.807):1.834);",
            "((F:1.167,(E:1.297,(D:1.766,C:1.47):1.987):1.269):1.1450000000000007,((H:1.592,G:1.734):1.386,((B:0.796,A:1.186):1.749,(J:1.237,I:0.206):0.863):2.5700000000000003):0.6619999999999993);",
            [1.982, 11.203, 9.807, 3.7329999999999997],
        ),
        (
            "(((H:0.665,(G:1.432,(F:0.653,E:1.593):1.163):1.686):0.467,(D:1.965,C:0.932):0.857):0.929,((B:0.196,A:0.111):1.27,(J:1.752,I:0.131):1.013):0.587);",
            "(((B:0.196,A:0.111):1.27,(J:1.752,I:0.131):1.013):0.587,((D:1.965,C:0.932):0.857,(H:0.665,(G:1.432,(F:0.653,E:1.593):1.163):1.686):0.467):0.929);",
            "((H:0.665,(G:1.432,(F:0.653,E:1.593):1.163):1.686):0.1529999999999998,((D:1.965,C:0.932):0.857,((B:0.196,A:0.111):1.27,(J:1.752,I:0.131):1.013):1.516):0.3140000000000002);",
            [0.307, 6.07, 3.9540000000000006, 2.246],
        ),
        (
            "((I:1.966,H:0.243):0.719,((G:1.852,(F:0.221,E:0.654):0.291):1.739,((D:1.955,(C:0.561,(B:0.081,A:1.034):1.451):1.868):0.522,J:0.811):1.607):0.687);",
            "((I:1.966,H:0.243):0.719,((G:1.852,(F:0.221,E:0.654):0.291):1.739,(J:0.811,(D:1.955,(C:0.561,(B:0.081,A:1.034):1.451):1.868):0.522):1.607):0.687);",
            "(((D:1.955,(C:0.561,(B:0.081,A:1.034):1.451):1.868):0.522,J:0.811):0.1615000000000002,((G:1.852,(F:0.221,E:0.654):0.291):1.739,(I:1.966,H:0.243):1.4060000000000001):1.4454999999999998);",
            [1.115, 3.7620000000000005, 5.7330000000000005, 0.875],
        ),
    ];
    /// ETE3's unrooted Robinson–Foulds distance and its maximum for case i against case i + 1.
    const RF: [(usize, usize); 4] = [(12, 14), (10, 14), (6, 14), (8, 14)];

    fn tree(text: &str) -> NwTree {
        NwTree::parse_all(text).unwrap().remove(0)
    }

    /// Smallest tip name below each node, to put children in a fixed order.
    fn smallest_tips(t: &NwTree) -> Vec<String> {
        let mut m = vec![String::new(); t.nodes.len()];
        for i in t.postorder() {
            m[i] = if t.is_leaf(i) {
                t.nodes[i].name.clone()
            } else {
                t.nodes[i]
                    .children
                    .iter()
                    .map(|&c| m[c].clone())
                    .min()
                    .unwrap()
            };
        }
        m
    }

    /// Whether two trees have the same shape, names and lengths (to `1e-9`); with `any_order`, the
    /// order of the children does not matter.
    fn same(a: &NwTree, b: &NwTree, any_order: bool) -> bool {
        struct Side<'a> {
            t: &'a NwTree,
            smallest: Vec<String>,
        }
        fn kids(s: &Side, i: usize, any_order: bool) -> Vec<usize> {
            let mut k = s.t.nodes[i].children.clone();
            if any_order {
                k.sort_by(|&x, &y| s.smallest[x].cmp(&s.smallest[y]));
            }
            k
        }
        fn go(a: &Side, i: usize, b: &Side, j: usize, any_order: bool) -> bool {
            let (x, y) = (&a.t.nodes[i], &b.t.nodes[j]);
            if x.children.len() != y.children.len() {
                return false;
            }
            let root = i == a.t.root;
            if !root && (x.dist - y.dist).abs() > 1e-9 {
                return false;
            }
            if x.children.is_empty() {
                return x.name == y.name;
            }
            kids(a, i, any_order)
                .into_iter()
                .zip(kids(b, j, any_order))
                .all(|(p, q)| go(a, p, b, q, any_order))
        }
        let (sa, sb) = (
            Side {
                t: a,
                smallest: smallest_tips(a),
            },
            Side {
                t: b,
                smallest: smallest_tips(b),
            },
        );
        go(&sa, a.root, &sb, b.root, any_order)
    }

    #[test]
    fn ladderizing_matches_ete() {
        for (text, ladder, _, _) in CASES {
            assert!(
                same(&ladderize(&tree(text)), &tree(ladder), false),
                "{text}"
            );
        }
    }

    #[test]
    fn midpoint_rooting_matches_dendropy() {
        for (text, _, mid, _) in CASES {
            let (rooted, found) = midpoint_root(&tree(text)).unwrap();
            assert!(
                same(&rooted, &tree(mid), true),
                "{text}\n got {}\nwant {mid}",
                rooted.newick()
            );
            // the two ends of the longest path are the same distance from the new root
            let (names, table) = patristic(&rooted);
            let at = |n: &str| names.iter().position(|x| x == n).unwrap();
            let longest = table[at(&found.tip_a)][at(&found.tip_b)];
            assert!((longest - found.diameter).abs() < 1e-9);
            let root_to = |n: &str| {
                let mut i = rooted.nodes.iter().position(|x| x.name == n).unwrap();
                let mut d = 0.0;
                while let Some(p) = rooted.nodes[i].parent {
                    d += rooted.nodes[i].dist;
                    i = p;
                }
                d
            };
            assert!((root_to(&found.tip_a) - root_to(&found.tip_b)).abs() < 1e-9);
        }
    }

    #[test]
    fn patristic_distances_match_ete() {
        for (text, _, _, wanted) in CASES {
            let (names, table) = patristic(&tree(text));
            let at = |n: &str| names.iter().position(|x| x == n).unwrap();
            for ((a, b), want) in [("A", "B"), ("C", "J"), ("D", "H"), ("E", "F")]
                .into_iter()
                .zip(wanted)
            {
                let got = table[at(a)][at(b)];
                assert!((got - want).abs() < 1e-9, "{a}-{b}: {got} vs {want}");
                assert_eq!(got, table[at(b)][at(a)]);
            }
        }
    }

    #[test]
    fn robinson_foulds_matches_ete() {
        for (i, (distance, max)) in RF.into_iter().enumerate() {
            let got = robinson_foulds(&tree(CASES[i].0), &tree(CASES[i + 1].0)).unwrap();
            assert_eq!((got.distance, got.max), (distance, max), "case {i}");
            assert_eq!(got.max - got.distance, 2 * got.shared);
        }
        // a tree is at distance 0 from itself, from its ladderized form and from its midpoint rooting
        let t = tree(CASES[0].0);
        for other in [t.clone(), ladderize(&t), midpoint_root(&t).unwrap().0] {
            assert_eq!(robinson_foulds(&t, &other).unwrap().distance, 0);
        }
    }

    #[test]
    fn robinson_foulds_refuses_different_tips() {
        let e = robinson_foulds(
            &tree("((A:1,B:1):1,(C:1,D:1):1);"),
            &tree("((A:1,B:1):1,(C:1,E:1):1);"),
        );
        assert!(e.unwrap_err().contains("same tips"));
        let twice = robinson_foulds(
            &tree("((A:1,A:1):1,(C:1,D:1):1);"),
            &tree("((A:1,B:1):1,(C:1,D:1):1);"),
        );
        assert!(twice.unwrap_err().contains("twice"));
    }

    /// Each label with the tips on its side of the split that leaves the first tip name out.
    fn label_splits(t: &NwTree) -> BTreeMap<String, Vec<String>> {
        let all = tip_names(t).unwrap();
        let mut out = BTreeMap::new();
        for i in t.levelorder() {
            if t.is_leaf(i) || i == t.root || t.nodes[i].name.is_empty() {
                continue;
            }
            let mut below: Vec<String> = Vec::new();
            let mut stack = vec![i];
            while let Some(j) = stack.pop() {
                if t.is_leaf(j) {
                    below.push(t.nodes[j].name.clone());
                } else {
                    stack.extend(&t.nodes[j].children);
                }
            }
            below.sort();
            if below.contains(&all[0]) {
                below = all.iter().filter(|n| !below.contains(n)).cloned().collect();
            }
            out.insert(t.nodes[i].name.clone(), below);
        }
        out
    }

    #[test]
    fn labels_stay_on_their_own_splits_when_the_root_moves() {
        // F's long branch holds the midpoint, so the root moves into the middle of the tree; the old
        // root's two labels (95 and 60) describe one split, which keeps the first of them
        let t = tree("((A:1,B:1)95:1,((C:1,D:1)80:1,(E:1,F:6)70:1)60:1);");
        let (rooted, found) = midpoint_root(&t).unwrap();
        assert_eq!((found.tip_a.as_str(), found.diameter), ("F", 10.0));
        let before = label_splits(&t);
        let after = label_splits(&rooted);
        assert_eq!(after.len(), 3, "{}", rooted.newick());
        for label in ["95", "80", "70"] {
            assert_eq!(after[label], before[label], "label {label}");
        }
    }

    #[test]
    fn small_and_degenerate_trees() {
        let (two, m) = midpoint_root(&tree("(A:1,B:3);")).unwrap();
        assert_eq!((m.diameter, two.leaf_count()), (4.0, 2));
        let lengths: Vec<f64> = two.nodes[two.root]
            .children
            .iter()
            .map(|&c| two.nodes[c].dist)
            .collect();
        assert_eq!(lengths, [2.0, 2.0]);
        assert!(midpoint_root(&tree("(A:0,B:0);")).is_err());
        assert!(midpoint_root(&tree("A;")).is_err());
        // a root edge and a three-way root are both fine
        assert!(midpoint_root(&tree("((A:1,B:2):1);")).is_ok());
        let (star, _) = midpoint_root(&tree("(A:1,B:2,C:5);")).unwrap();
        assert_eq!(star.leaf_count(), 3);
    }

    #[test]
    fn summary_counts_and_the_distance_table_is_square() {
        let t = tree("((A:1,B:2):3,C:4);");
        let s = summary(&t);
        assert_eq!((s.tips, s.inner, s.binary), (3, 2, true));
        assert_eq!((s.length, s.height), (10.0, 5.0));
        assert!(!summary(&tree("(A:1,B:1,C:1);")).binary);
        let table = distance_table(&t);
        let rows: Vec<&str> = table.lines().collect();
        assert_eq!(rows[0], "\tA\tB\tC");
        assert_eq!(rows[1], "A\t0\t3\t8");
        assert_eq!(rows[3], "C\t8\t9\t0");
    }
}
