//! Distance trees: distances from a multiple alignment, then neighbour joining (Saitou & Nei 1987)
//! or UPGMA, with bootstrap support. The first MEGA-style pieces of the platform. Written from the
//! published formulas; no code from MEGA, IQ-TREE or any other tree program is used.

use super::msa::Msa;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DistanceModel {
    /// Fraction of sites that differ.
    P,
    /// Jukes-Cantor 1969: d = -3/4 ln(1 - 4p/3).
    JukesCantor,
    /// Kimura 2-parameter: d = -1/2 ln(1 - 2P - Q) - 1/4 ln(1 - 2Q), P transitions, Q transversions.
    Kimura2,
}

impl DistanceModel {
    pub const ALL: [DistanceModel; 3] = [
        DistanceModel::P,
        DistanceModel::JukesCantor,
        DistanceModel::Kimura2,
    ];

    pub fn label(self) -> &'static str {
        match self {
            DistanceModel::P => "p-distance",
            DistanceModel::JukesCantor => "Jukes-Cantor",
            DistanceModel::Kimura2 => "Kimura 2-parameter",
        }
    }
}

fn is_purine(c: u8) -> bool {
    c == b'A' || c == b'G'
}

/// Distance between two aligned sequences of equal length. Columns with a gap or an ambiguity
/// code in either sequence are skipped (pairwise deletion). `None` when no column can be
/// compared or the sequences differ too much for the model (the logarithm is undefined).
pub fn aligned_distance(a: &[u8], b: &[u8], model: DistanceModel) -> Option<f64> {
    let (mut sites, mut transitions, mut transversions) = (0u32, 0u32, 0u32);
    for (&x, &y) in a.iter().zip(b) {
        let (x, y) = (x.to_ascii_uppercase(), y.to_ascii_uppercase());
        if !(b"ACGT".contains(&x) && b"ACGT".contains(&y)) {
            continue;
        }
        sites += 1;
        if x != y {
            if is_purine(x) == is_purine(y) {
                transitions += 1;
            } else {
                transversions += 1;
            }
        }
    }
    distance_from_counts(sites, transitions, transversions, model)
}

/// The distance between two sequences from how many sites were compared and how many of them differ by a transition and by a
/// transversion. `None` when no site was compared or the sequences differ too much for the model (the logarithm is
/// undefined). Counts are whole numbers, so whoever counts them (the processor, a GPU) gets the same distance.
pub fn distance_from_counts(
    sites: u32,
    transitions: u32,
    transversions: u32,
    model: DistanceModel,
) -> Option<f64> {
    if sites == 0 {
        return None;
    }
    let n = f64::from(sites);
    let (ts, tv) = (f64::from(transitions) / n, f64::from(transversions) / n);
    let p = ts + tv;
    let d = match model {
        DistanceModel::P => Some(p),
        DistanceModel::JukesCantor => {
            let arg = 1.0 - 4.0 * p / 3.0;
            (arg > 0.0).then(|| -0.75 * arg.ln())
        }
        DistanceModel::Kimura2 => {
            let (first, second) = (1.0 - 2.0 * ts - tv, 1.0 - 2.0 * tv);
            (first > 0.0 && second > 0.0).then(|| -0.5 * first.ln() - 0.25 * second.ln())
        }
    };
    d.map(f64::abs) // turns -0.0 into 0.0
}

/// The distances of every pair from the GPUs, or the first pair that cannot be compared; none if the cards failed.
type Distances = Option<Result<Vec<Vec<f64>>, (usize, usize)>>;

/// What is said of a pair of rows that can't be compared.
fn not_comparable(msa: &Msa, i: usize, j: usize, model: DistanceModel) -> String {
    format!(
        "{} and {} are too different to give a {} distance",
        msa.names[i],
        msa.names[j],
        model.label()
    )
}

/// Distances between every pair of rows of an alignment, on the GPUs when there are some and there are enough pairs and columns
/// to be worth it (the same numbers to the last digit: `accel::distance`), else on the processor. `Err` names the first pair that
/// can't be compared.
pub fn msa_distances(msa: &Msa, model: DistanceModel) -> Result<Vec<Vec<f64>>, String> {
    if let Ok(Some(cards)) = Cards::for_work(msa, 1) {
        if let Some(d) = cards.distances(None, model) {
            return d.map_err(|(i, j)| not_comparable(msa, i, j, model));
        }
    }
    cpu_distances(msa, model)
}

/// The same on the processor, whatever the size (the reference that the GPUs' answers are checked against).
pub(crate) fn cpu_distances(msa: &Msa, model: DistanceModel) -> Result<Vec<Vec<f64>>, String> {
    let n = msa.rows.len();
    let mut d = vec![vec![0.0; n]; n];
    for (i, j) in (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))) {
        let v = aligned_distance(&msa.rows[i], &msa.rows[j], model)
            .ok_or_else(|| not_comparable(msa, i, j, model))?;
        (d[i][j], d[j][i]) = (v, v);
    }
    Ok(d)
}

/// The columns of an alignment on the GPUs, for any number of weightings of them (the resamples of a bootstrap): nothing where
/// there is no GPU, it is not worth it, or the alignment does not fit.
#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
pub struct Cards {
    engine: crate::accel::Distances,
    /// Why a card failed, the first time one did: the processor does everything after that.
    failure: std::sync::Mutex<Option<String>>,
}

#[cfg(not(all(feature = "gpu", not(target_arch = "wasm32"))))]
pub struct Cards;

#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
impl Cards {
    fn new(engine: crate::accel::Distances) -> Cards {
        Cards {
            engine,
            failure: std::sync::Mutex::new(None),
        }
    }

    /// The columns of `msa` on the GPUs, if `times` sets of distances of this size are worth it (see `accel::threshold`):
    /// `Ok(None)` where they are not (or there is no GPU), and `Err` with the reason where they were meant to be used and could
    /// not be (said on the standard error too).
    fn for_work(msa: &Msa, times: usize) -> Result<Option<Cards>, String> {
        use crate::accel::{note_failure, threshold, worth_it, Distances};
        let n = msa.rows.len() as f64;
        let cells = times as f64 * n * n * msa.columns() as f64;
        if !worth_it(cells, threshold::distance_cells()) {
            return Ok(None);
        }
        match Distances::open(&msa.rows) {
            Ok(engine) => Ok(engine.map(Cards::new)),
            Err(why) => {
                note_failure(
                    "the GPUs could not take the distances",
                    &format!("{why}; the processor does them"),
                );
                Err(why)
            }
        }
    }

    /// Why a card failed, if one has.
    fn failure(&self) -> Option<String> {
        self.failure
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone()
    }

    /// The distances with each column taken `weights[column]` times: `Err` is the first pair that cannot be compared; `None` if a
    /// card failed (this time or earlier: a card that failed is not asked again), and the processor must do it.
    fn distances(&self, weights: Option<&[u32]>, model: DistanceModel) -> Distances {
        if self.failure().is_some() {
            return None;
        }
        match self.engine.counts(weights) {
            Ok(counts) => Some(counts.distances(model)),
            Err(why) => {
                crate::accel::note_failure(
                    "the GPUs failed on distances",
                    &format!("{why}; the processor does them"),
                );
                self.failure
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner)
                    .get_or_insert(why);
                None
            }
        }
    }

    /// Where the distances of `resamples` bootstrap resamples were computed, `done` of them on the cards.
    fn report(&self, done: usize, resamples: usize) -> String {
        let place = self.engine.describe();
        match self.failure() {
            None => place,
            Some(why) if done == 0 => format!("the processor (a GPU failed: {why})"),
            Some(why) => format!(
                "{place} for {done} of {resamples} resamples, then a GPU failed ({why}) and the processor did the rest"
            ),
        }
    }
}

#[cfg(not(all(feature = "gpu", not(target_arch = "wasm32"))))]
impl Cards {
    fn for_work(_: &Msa, _: usize) -> Result<Option<Cards>, String> {
        Ok(None)
    }

    fn distances(&self, _: Option<&[u32]>, _: DistanceModel) -> Distances {
        None
    }

    fn report(&self, _: usize, _: usize) -> String {
        String::new()
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Node {
    /// Sequence name for a leaf, empty inside the tree.
    pub name: String,
    /// Length of the branch leading to this node.
    pub length: f64,
    /// Bootstrap or posterior support in percent, on inner branches only.
    pub support: Option<f64>,
    pub children: Vec<Node>,
}

impl Node {
    fn leaf(name: &str) -> Node {
        Node {
            name: name.to_string(),
            length: 0.0,
            support: None,
            children: Vec::new(),
        }
    }

    pub fn inner(children: Vec<Node>) -> Node {
        Node {
            name: String::new(),
            length: 0.0,
            support: None,
            children,
        }
    }

    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    pub fn leaf_count(&self) -> usize {
        if self.is_leaf() {
            1
        } else {
            self.children.iter().map(Node::leaf_count).sum()
        }
    }

    /// Longest path from this node down to a leaf, including its own branch.
    pub fn depth(&self) -> f64 {
        self.length + self.children.iter().map(Node::depth).fold(0.0, f64::max)
    }

    /// Newick text with a closing `;`; the root has no branch length. Support values follow the
    /// closing bracket of their branch, as most tree programs write them.
    pub fn newick(&self) -> String {
        format!("{};", self.newick_inner(false))
    }

    /// Newick text that reads back (with `NwTree::parse_all` and `to_node`) as exactly this tree. `newick` makes names
    /// plain, lengths six decimals and supports whole numbers, which is right for a file to read; this quotes a name that
    /// holds anything Newick uses and writes lengths and supports in full, so that a tree sent to another computer and
    /// back is the same tree.
    pub fn newick_exact(&self) -> String {
        format!("{};", self.exact_inner(false))
    }

    fn exact_inner(&self, with_length: bool) -> String {
        let mut text = if self.is_leaf() {
            quoted_name(&self.name)
        } else {
            let kids: Vec<String> = self.children.iter().map(|c| c.exact_inner(true)).collect();
            format!("({})", kids.join(","))
        };
        if let (Some(s), false) = (self.support, self.is_leaf()) {
            text.push_str(&s.to_string());
        }
        if with_length {
            text.push_str(&format!(":{}", self.length));
        }
        text
    }

    fn newick_inner(&self, with_length: bool) -> String {
        let label: String = self
            .name
            .chars()
            .map(|c| {
                if c.is_whitespace() || "()[]:;,'".contains(c) {
                    '_'
                } else {
                    c
                }
            })
            .collect();
        let mut text = if self.is_leaf() {
            label
        } else {
            let kids: Vec<String> = self.children.iter().map(|c| c.newick_inner(true)).collect();
            format!("({})", kids.join(","))
        };
        if let (Some(s), false) = (self.support, self.is_leaf()) {
            text.push_str(&format!("{s:.0}"));
        }
        if with_length {
            text.push_str(&format!(":{:.6}", self.length));
        }
        text
    }
}

/// A name as Newick writes it: as it is when it holds nothing that Newick uses, else in single quotes with each quote in it
/// doubled.
fn quoted_name(name: &str) -> String {
    let plain = |c: char| !(c.is_whitespace() || "()[]:;,'".contains(c));
    if name.chars().all(plain) {
        name.to_string()
    } else {
        format!("'{}'", name.replace('\'', "''"))
    }
}

/// One join of neighbour-joining: the positions, in the current list of clusters, of the two that are joined (`i` before `j`), and the
/// lengths of their branches to the new cluster, which goes to the end of the list.
#[derive(Clone, Debug, PartialEq)]
pub struct Join {
    pub i: usize,
    pub j: usize,
    pub li: f64,
    pub lj: f64,
}

/// Neighbour-joining tree from a symmetric distance matrix. The tree is unrooted, drawn from a
/// root with three children (two for two sequences). Negative branch lengths are set to zero.
pub fn neighbor_joining(names: &[String], dist: &[Vec<f64>]) -> Option<Node> {
    let n = names.len();
    if n < 2 || dist.len() != n || dist.iter().any(|row| row.len() != n) {
        return None;
    }
    // many sequences: the GPUs make the joins (the same ones, but for ties that rounding decides)
    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    if let Some((joins, last)) = crate::accel::neighbor_joins(dist) {
        return Some(join_tree(names, &joins, &last));
    }
    let (joins, last) = cpu_joins(dist);
    Some(join_tree(names, &joins, &last))
}

// The joins of neighbour-joining on the processor, and the distances among the (two or three) clusters that are left: in
// `tree/joins.rs`, which makes the same joins to the last bit as the plain way and shares the work among the cores.
mod joins;
pub(crate) use joins::cpu_joins;

/// The tree that `joins` make of the sequences `names`, the last two or three clusters being `last` apart.
pub(crate) fn join_tree(names: &[String], joins: &[Join], last: &[Vec<f64>]) -> Node {
    let mut nodes: Vec<Node> = names.iter().map(|s| Node::leaf(s)).collect();
    for join in joins {
        let (mut right, mut left) = (nodes.remove(join.j), nodes.remove(join.i));
        (left.length, right.length) = (join.li, join.lj);
        nodes.push(Node::inner(vec![left, right]));
    }
    let lengths: Vec<f64> = if nodes.len() == 2 {
        vec![last[0][1] / 2.0; 2]
    } else {
        let (ab, ac, bc) = (last[0][1], last[0][2], last[1][2]);
        vec![
            (ab + ac - bc) / 2.0,
            (ab + bc - ac) / 2.0,
            (ac + bc - ab) / 2.0,
        ]
    };
    for (node, length) in nodes.iter_mut().zip(lengths) {
        node.length = length.max(0.0);
    }
    Node::inner(nodes)
}

/// Drop rows and columns `i < j` of `d` and append the merged cluster with distances `to_new`.
fn join_rows(d: &mut Vec<Vec<f64>>, i: usize, j: usize, to_new: Vec<f64>) {
    d.remove(j);
    d.remove(i);
    for (row, &x) in d.iter_mut().zip(&to_new) {
        row.remove(j);
        row.remove(i);
        row.push(x);
    }
    let mut last = to_new;
    last.push(0.0);
    d.push(last);
}

/// UPGMA: clusters are joined by average distance, giving a rooted tree whose leaves all end at
/// the same height (a molecular clock).
pub fn upgma(names: &[String], dist: &[Vec<f64>]) -> Option<Node> {
    let n = names.len();
    if n < 2 || dist.len() != n || dist.iter().any(|row| row.len() != n) {
        return None;
    }
    // each cluster: its tree, height above the leaves, and number of leaves
    let mut clusters: Vec<(Node, f64, f64)> =
        names.iter().map(|s| (Node::leaf(s), 0.0, 1.0)).collect();
    let mut d: Vec<Vec<f64>> = dist.to_vec();
    while clusters.len() > 1 {
        let k = clusters.len();
        let (dij, i, j) = (0..k)
            .flat_map(|i| (i + 1..k).map(move |j| (i, j)))
            .map(|(i, j)| (d[i][j], i, j))
            .fold(
                (f64::INFINITY, 0, 1),
                |best, c| if c.0 < best.0 { c } else { best },
            );
        let height = dij / 2.0;
        let to_new: Vec<f64> = (0..clusters.len())
            .filter(|&k| k != i && k != j)
            .map(|k| {
                let (si, sj) = (clusters[i].2, clusters[j].2);
                (si * d[i][k] + sj * d[j][k]) / (si + sj)
            })
            .collect();
        let ((mut right, hr, sr), (mut left, hl, sl)) = (clusters.remove(j), clusters.remove(i));
        (left.length, right.length) = ((height - hl).max(0.0), (height - hr).max(0.0));
        clusters.push((Node::inner(vec![left, right]), height, sl + sr));
        join_rows(&mut d, i, j, to_new);
    }
    clusters.pop().map(|(root, _, _)| root)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    NeighborJoining,
    Upgma,
}

impl Method {
    pub fn tree(self, names: &[String], dist: &[Vec<f64>]) -> Option<Node> {
        match self {
            Method::NeighborJoining => neighbor_joining(names, dist),
            Method::Upgma => upgma(names, dist),
        }
    }
}

/// A small deterministic random generator (xorshift64*), so a bootstrap can be repeated exactly.
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed.max(1))
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n.max(1) as u64) as usize
    }

    /// A number in [0, 1).
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// The leaves below each inner branch (all but the root), as sorted indexes into `names`, taken
/// from the side without leaf 0 so the same split reads the same however a tree is rooted.
pub fn splits(tree: &Node, names: &[String]) -> HashSet<Vec<usize>> {
    let index: HashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let mut out = HashSet::new();
    for child in &tree.children {
        collect_splits(child, &index, names.len(), &mut out);
    }
    out
}

fn collect_splits(
    node: &Node,
    index: &HashMap<&str, usize>,
    total: usize,
    out: &mut HashSet<Vec<usize>>,
) -> Vec<usize> {
    if node.is_leaf() {
        return index.get(node.name.as_str()).copied().into_iter().collect();
    }
    let mut below: Vec<usize> = node
        .children
        .iter()
        .flat_map(|c| collect_splits(c, index, total, out))
        .collect();
    below.sort_unstable();
    let key = if below.first() == Some(&0) {
        (0..total)
            .filter(|i| below.binary_search(i).is_err())
            .collect()
    } else {
        below.clone()
    };
    out.insert(key);
    below
}

fn set_support(
    node: &mut Node,
    index: &HashMap<&str, usize>,
    total: usize,
    counts: &HashMap<Vec<usize>, f64>,
) {
    fn leaves(node: &Node, index: &HashMap<&str, usize>, out: &mut Vec<usize>) {
        if node.is_leaf() {
            out.extend(index.get(node.name.as_str()));
        }
        node.children.iter().for_each(|c| leaves(c, index, out));
    }
    for child in &mut node.children {
        if !child.is_leaf() {
            let mut below = Vec::new();
            leaves(child, index, &mut below);
            below.sort_unstable();
            let key: Vec<usize> = if below.first() == Some(&0) {
                (0..total)
                    .filter(|i| below.binary_search(i).is_err())
                    .collect()
            } else {
                below
            };
            child.support = counts.get(&key).copied();
            set_support(child, index, total, counts);
        }
    }
}

/// Bootstrap support (percent) for the inner branches of `tree`: the alignment's columns are
/// resampled with replacement `replicates` times, a tree is built from each resample with the
/// same method, and a branch scores the share of resampled trees that contain its split.
pub fn add_bootstrap(
    tree: &mut Node,
    msa: &Msa,
    model: DistanceModel,
    method: Method,
    replicates: usize,
    seed: u64,
) -> Option<String> {
    match Cards::for_work(msa, replicates) {
        Ok(cards) => add_bootstrap_with(tree, msa, model, method, replicates, seed, cards),
        Err(why) => {
            add_bootstrap_with(tree, msa, model, method, replicates, seed, None);
            Some(format!("the processor (the GPUs could not be used: {why})"))
        }
    }
}

/// The same with the columns on the GPUs (or not) as `cards` says: every resample is then a weight for each column, and its
/// distances come from the cards without the columns being sent again.
fn add_bootstrap_with(
    tree: &mut Node,
    msa: &Msa,
    model: DistanceModel,
    method: Method,
    replicates: usize,
    seed: u64,
    cards: Option<Cards>,
) -> Option<String> {
    let (width, mut rng) = (msa.columns(), Rng::new(seed));
    let mut counts: HashMap<Vec<usize>, f64> = HashMap::new();
    let mut usable = 0.0;
    let mut on_cards_count = 0;
    for _ in 0..replicates {
        let picks: Vec<usize> = (0..width).map(|_| rng.below(width)).collect();
        let on_cards = cards.as_ref().and_then(|c| {
            let mut weights = vec![0u32; width];
            picks.iter().for_each(|&p| weights[p] += 1);
            c.distances(Some(&weights), model)
        });
        on_cards_count += usize::from(on_cards.is_some());
        let d = match on_cards {
            Some(Ok(d)) => d,
            // two sequences that cannot be compared in this resample: it is skipped, as on the processor
            Some(Err(_)) => continue,
            None => {
                let rows = msa
                    .rows
                    .iter()
                    .map(|r| picks.iter().map(|&c| r[c]).collect())
                    .collect();
                let resampled = Msa {
                    names: msa.names.clone(),
                    rows,
                };
                let Ok(d) = cpu_distances(&resampled, model) else {
                    continue;
                };
                d
            }
        };
        if let Some(t) = method.tree(&msa.names, &d) {
            usable += 1.0;
            for key in splits(&t, &msa.names) {
                *counts.entry(key).or_insert(0.0) += 100.0;
            }
        }
    }
    if usable > 0.0 {
        counts.values_mut().for_each(|c| *c /= usable);
        let index: HashMap<&str, usize> = msa
            .names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i))
            .collect();
        set_support(tree, &index, msa.names.len(), &counts);
    }
    cards.as_ref().map(|c| c.report(on_cards_count, replicates))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::msa::center_star;

    /// The branch supports of a tree in the order of a walk down it.
    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    fn supports(node: &Node, out: &mut Vec<Option<f64>>) {
        out.push(node.support);
        node.children.iter().for_each(|c| supports(c, out));
    }

    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    #[test]
    fn a_bootstrap_on_the_cards_gives_the_support_of_a_bootstrap_on_the_processor_for_the_same_seed(
    ) {
        use crate::accel::distance::Engine;
        use crate::accel::pool::Member;
        use crate::accel::Pool;
        // three groups of four sequences, each group a copy of its own ancestor with a few changes: a tree with something to support
        let mut rng = Rng::new(77);
        let letters = b"ACGT";
        let ancestor: Vec<u8> = (0..300).map(|_| letters[rng.below(4)]).collect();
        let mut rows: Vec<Vec<u8>> = Vec::new();
        for _ in 0..3 {
            let group: Vec<u8> = ancestor
                .iter()
                .map(|&c| {
                    if rng.below(100) < 18 {
                        letters[rng.below(4)]
                    } else {
                        c
                    }
                })
                .collect();
            for _ in 0..4 {
                rows.push(
                    group
                        .iter()
                        .map(|&c| match rng.below(60) {
                            0 => b'-',
                            1..=6 => letters[rng.below(4)],
                            _ => c,
                        })
                        .collect(),
                );
            }
        }
        let msa = Msa {
            names: (0..rows.len()).map(|k| format!("s{k}")).collect(),
            rows,
        };
        let model = DistanceModel::JukesCantor;
        let start = neighbor_joining(&msa.names, &cpu_distances(&msa, model).unwrap()).unwrap();
        let mut pools = vec![
            Pool::new(vec![Member::new(
                candle_core::Device::Cpu,
                "CPU".into(),
                1.0,
            )]),
            Pool::new(vec![
                Member::new(candle_core::Device::Cpu, "CPU".into(), 2.0),
                Member::new(candle_core::Device::Cpu, "CPU".into(), 1.0),
            ]),
        ];
        pools.extend(Pool::global().cloned());
        for method in [Method::NeighborJoining, Method::Upgma] {
            let mut on_cpu = start.clone();
            add_bootstrap_with(&mut on_cpu, &msa, model, method, 40, 5, None);
            let mut want = Vec::new();
            supports(&on_cpu, &mut want);
            assert!(
                want.iter().flatten().any(|&s| s > 0.0),
                "there is something to support"
            );
            for pool in &pools {
                let cards = Cards::new(crate::accel::Distances::Candle(
                    Engine::new(pool, &msa.rows).unwrap(),
                ));
                let mut on_cards = start.clone();
                add_bootstrap_with(&mut on_cards, &msa, model, method, 40, 5, Some(cards));
                let mut got = Vec::new();
                supports(&on_cards, &mut got);
                assert_eq!(got, want, "{} cards", pool.len());
            }
        }
    }

    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    #[test]
    fn a_card_that_fails_in_a_bootstrap_hands_every_resample_to_the_processor_and_the_report_says_so(
    ) {
        use crate::accel::distance::Engine;
        use crate::accel::pool::Member;
        use crate::accel::Pool;
        let mut rng = Rng::new(3);
        let letters = b"ACGT";
        let mut rows = |columns: usize| -> Vec<Vec<u8>> {
            (0..6)
                .map(|_| (0..columns).map(|_| letters[rng.below(4)]).collect())
                .collect()
        };
        let names: Vec<String> = (0..6).map(|k| format!("s{k}")).collect();
        let msa = Msa {
            names: names.clone(),
            rows: rows(120),
        };
        // columns for a different alignment: a resample of this one is the wrong size for them, which the cards refuse
        let others = rows(60);
        let pool = Pool::new(vec![Member::new(
            candle_core::Device::Cpu,
            "CPU".into(),
            1.0,
        )]);
        let model = DistanceModel::P;
        let start = neighbor_joining(&names, &cpu_distances(&msa, model).unwrap()).unwrap();
        let mut on_cpu = start.clone();
        add_bootstrap_with(
            &mut on_cpu,
            &msa,
            model,
            Method::NeighborJoining,
            12,
            5,
            None,
        );
        let mut with_failure = start.clone();
        let cards = Cards::new(crate::accel::Distances::Candle(
            Engine::new(&pool, &others).unwrap(),
        ));
        let said = add_bootstrap_with(
            &mut with_failure,
            &msa,
            model,
            Method::NeighborJoining,
            12,
            5,
            Some(cards),
        )
        .expect("a report");
        assert!(said.contains("a GPU failed"), "{said}");
        assert!(said.contains("one weight for each column"), "{said}");
        // the processor did every resample, so the support is what it gives alone
        let (mut got, mut want) = (Vec::new(), Vec::new());
        supports(&with_failure, &mut got);
        supports(&on_cpu, &mut want);
        assert_eq!(got, want);
    }

    fn names(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    /// (leaf name, branch lengths from the root down to it, each with a node number)
    fn paths(
        node: &Node,
        trail: &mut Vec<(usize, f64)>,
        next: &mut usize,
        out: &mut Vec<(String, Vec<(usize, f64)>)>,
    ) {
        trail.push((*next, node.length));
        *next += 1;
        if node.is_leaf() {
            out.push((node.name.clone(), trail.clone()));
        }
        for child in &node.children {
            paths(child, trail, next, out);
        }
        trail.pop();
    }

    /// Length of the path between two leaves: the branches not shared by their root paths.
    fn between(tree: &Node, x: &str, y: &str) -> f64 {
        let mut all = Vec::new();
        paths(tree, &mut Vec::new(), &mut 0, &mut all);
        let find = |name: &str| all.iter().find(|(n, _)| n == name).unwrap().1.clone();
        let (a, b) = (find(x), find(y));
        let only = |from: &[(usize, f64)], other: &[(usize, f64)]| -> f64 {
            from.iter()
                .filter(|(id, _)| !other.iter().any(|(o, _)| o == id))
                .map(|e| e.1)
                .sum()
        };
        only(&a, &b) + only(&b, &a)
    }

    #[test]
    fn distance_models_follow_their_formulas() {
        let a = b"AAAAAAAAAAAAAAAAAAAA";
        let b = b"GGCAAAAAAAAAAAAAAAAA"; // 2 transitions, 1 transversion in 20 sites
        let d = |m| aligned_distance(a, b, m).unwrap();
        assert!((d(DistanceModel::P) - 0.15).abs() < 1e-12);
        assert!((d(DistanceModel::JukesCantor) - 0.167_357_66).abs() < 1e-7);
        assert!((d(DistanceModel::Kimura2) - 0.170_181_17).abs() < 1e-7);
        assert_eq!(
            aligned_distance(a, a, DistanceModel::JukesCantor),
            Some(0.0)
        );
    }

    #[test]
    fn gaps_and_ambiguity_codes_are_skipped_and_saturation_gives_none() {
        let d = aligned_distance(b"ACGT--NA", b"ACGTAAGA", DistanceModel::P);
        assert_eq!(d, Some(0.0)); // five comparable columns, all identical
        assert_eq!(aligned_distance(b"---", b"ACG", DistanceModel::P), None);
        // p = 0.75 or more has no Jukes-Cantor distance
        assert_eq!(
            aligned_distance(b"AAAA", b"CCCG", DistanceModel::JukesCantor),
            None
        );
    }

    #[test]
    fn neighbour_joining_reproduces_an_additive_matrix() {
        let n = names(&["a", "b", "c", "d", "e"]);
        let d = vec![
            vec![0.0, 5.0, 9.0, 9.0, 8.0],
            vec![5.0, 0.0, 10.0, 10.0, 9.0],
            vec![9.0, 10.0, 0.0, 8.0, 7.0],
            vec![9.0, 10.0, 8.0, 0.0, 3.0],
            vec![8.0, 9.0, 7.0, 3.0, 0.0],
        ];
        let tree = neighbor_joining(&n, &d).unwrap();
        assert_eq!(tree.leaf_count(), 5);
        for i in 0..5 {
            for j in i + 1..5 {
                let got = between(&tree, &n[i], &n[j]);
                assert!((got - d[i][j]).abs() < 1e-9, "{} {}: {got}", n[i], n[j]);
            }
        }
    }

    #[test]
    fn three_sequences_make_a_star_and_two_make_a_split_branch() {
        let star = neighbor_joining(
            &names(&["x", "y", "z"]),
            &[
                vec![0.0, 3.0, 4.0],
                vec![3.0, 0.0, 5.0],
                vec![4.0, 5.0, 0.0],
            ],
        )
        .unwrap();
        assert_eq!(star.newick(), "(x:1.000000,y:2.000000,z:3.000000);");
        let two = neighbor_joining(&names(&["p", "q"]), &[vec![0.0, 4.0], vec![4.0, 0.0]]).unwrap();
        assert_eq!(two.newick(), "(p:2.000000,q:2.000000);");
        assert!(neighbor_joining(&names(&["only"]), &[vec![0.0]]).is_none());
        assert!((star.depth() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn names_are_made_safe_for_newick() {
        let t = neighbor_joining(
            &names(&["a b", "c:d", "e(f)"]),
            &[
                vec![0.0, 1.0, 1.0],
                vec![1.0, 0.0, 1.0],
                vec![1.0, 1.0, 0.0],
            ],
        )
        .unwrap();
        assert_eq!(t.newick(), "(a_b:0.500000,c_d:0.500000,e_f_:0.500000);");
    }

    #[test]
    fn upgma_joins_the_closest_pairs_and_ends_every_leaf_at_the_same_height() {
        let d = vec![
            vec![0.0, 4.0, 8.0, 8.0],
            vec![4.0, 0.0, 8.0, 8.0],
            vec![8.0, 8.0, 0.0, 4.0],
            vec![8.0, 8.0, 4.0, 0.0],
        ];
        let tree = upgma(&names(&["a", "b", "c", "d"]), &d).unwrap();
        assert_eq!(
            tree.newick(),
            "((a:2.000000,b:2.000000):2.000000,(c:2.000000,d:2.000000):2.000000);"
        );
        assert!((tree.depth() - 4.0).abs() < 1e-12);
    }

    /// 60 bases from a small deterministic generator.
    fn base() -> Vec<u8> {
        let mut x = 12345u32;
        (0..60)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                b"ACGT"[((x >> 16) & 3) as usize]
            })
            .collect()
    }

    fn mutated(positions: &[usize]) -> Vec<u8> {
        let mut s = base();
        for &p in positions {
            s[p] = match s[p] {
                b'A' => b'G',
                b'G' => b'A',
                b'C' => b'T',
                _ => b'C',
            };
        }
        s
    }

    fn four_related() -> Msa {
        let seqs = vec![
            ("A".to_string(), mutated(&[])),
            ("B".to_string(), mutated(&[3, 20])),
            (
                "C".to_string(),
                mutated(&[5, 9, 14, 22, 30, 41, 50, 55, 57, 58, 59, 1]),
            ),
            (
                "D".to_string(),
                mutated(&[5, 9, 14, 22, 30, 41, 50, 55, 57, 58, 59, 1, 33, 44]),
            ),
        ];
        center_star(&seqs).unwrap()
    }

    #[test]
    fn related_sequences_group_with_their_closest_relatives() {
        let msa = four_related();
        let d = msa_distances(&msa, DistanceModel::JukesCantor).unwrap();
        assert!(d[0][1] < d[0][2] && d[2][3] < d[0][3]);
        let tree = neighbor_joining(&msa.names, &d).unwrap();
        let found = splits(&tree, &msa.names);
        // A and B against C and D, written as the side without leaf 0
        assert!(found.contains(&vec![2, 3]), "{found:?}");
    }

    #[test]
    fn a_pair_that_cannot_be_compared_is_named_in_the_error() {
        let msa = Msa {
            names: names(&["one", "two"]),
            rows: vec![b"AAAA".to_vec(), b"NNNN".to_vec()],
        };
        let err = msa_distances(&msa, DistanceModel::P).unwrap_err();
        assert!(err.contains("one") && err.contains("two"), "{err}");
    }

    #[test]
    fn bootstrap_supports_a_clear_split_and_repeats_with_the_same_seed() {
        let msa = four_related();
        let d = msa_distances(&msa, DistanceModel::JukesCantor).unwrap();
        let mut tree = neighbor_joining(&msa.names, &d).unwrap();
        add_bootstrap(
            &mut tree,
            &msa,
            DistanceModel::JukesCantor,
            Method::NeighborJoining,
            40,
            7,
        );
        let support: Vec<f64> = tree.children.iter().filter_map(|c| c.support).collect();
        assert_eq!(support.len(), 1, "one inner branch in a four-leaf tree");
        assert!(support[0] > 60.0, "{support:?}");
        let mut again = neighbor_joining(&msa.names, &d).unwrap();
        add_bootstrap(
            &mut again,
            &msa,
            DistanceModel::JukesCantor,
            Method::NeighborJoining,
            40,
            7,
        );
        assert_eq!(tree, again);
        assert!(tree.newick().contains(')'));
    }
}
