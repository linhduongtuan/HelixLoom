//! Maximum-likelihood trees in the manner of IQ-TREE, much smaller: site patterns from an
//! alignment, Felsenstein's pruning algorithm for the DNA models JC69, K80, F81 and HKY85, branch
//! lengths and the transition/transversion ratio by golden-section search, the model chosen by
//! BIC, and a nearest-neighbour-interchange search started from the neighbour-joining tree.
//! Not included: rate variation between sites, partitions, ultrafast bootstrap. Written from the
//! published methods; no IQ-TREE code was used.

use super::msa::Msa;
use super::par::par_map;
use super::tree::{msa_distances, neighbor_joining, DistanceModel, Node};
use std::collections::{BTreeMap, HashMap};

pub mod bayes;
#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
pub mod check;
#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
mod engine;
#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
pub mod gpu;
#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
mod plan;
#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
pub mod torch;

type Matrix = [[f64; 4]; 4];

const MIN_BRANCH: f64 = 1e-6;
const MAX_BRANCH: f64 = 5.0;
/// Golden-section tolerance on the natural log of a branch length or of kappa.
const LOG_TOLERANCE: f64 = 1e-3;
const KAPPA_RANGE: (f64, f64) = (0.1, 100.0);
const TAYLOR_TERMS: u32 = 12;
/// Most topology changes accepted in one search, and the smallest log-likelihood gain that counts.
const NNI_MOVES: usize = 20;
const NNI_MIN_GAIN: f64 = 1e-3;
/// The share of the log-likelihood that a change must exceed on the GPUs (see `Sites::smallest_gain`).
#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
const GPU_NOISE_SHARE: f64 = 1e-6;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Jc69,
    K80,
    F81,
    Hky85,
}

impl Kind {
    pub const ALL: [Kind; 4] = [Kind::Jc69, Kind::K80, Kind::F81, Kind::Hky85];

    pub fn label(self) -> &'static str {
        match self {
            Kind::Jc69 => "JC69",
            Kind::K80 => "K80",
            Kind::F81 => "F81",
            Kind::Hky85 => "HKY85",
        }
    }

    fn free_kappa(self) -> bool {
        matches!(self, Kind::K80 | Kind::Hky85)
    }

    fn free_freqs(self) -> bool {
        matches!(self, Kind::F81 | Kind::Hky85)
    }
}

/// The alignment as distinct columns (patterns) with how often each occurs. States: 0..=3 for
/// A, C, G, T; 4 for a gap or ambiguity code, which fits any base.
struct Sites {
    tips: usize,
    patterns: Vec<Vec<u8>>,
    counts: Vec<f64>,
    total: f64,
    /// Base frequencies of the alignment (every base counted once more, so none is zero).
    freqs: [f64; 4],
    /// The columns on the GPUs, when there are some and the alignment is big enough to be worth it (`accelerate`).
    #[cfg_attr(
        not(all(feature = "gpu", not(target_arch = "wasm32"))),
        allow(dead_code)
    )]
    accel: Accel,
}

#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
type Accel = Option<engine::Engine>;
#[cfg(not(all(feature = "gpu", not(target_arch = "wasm32"))))]
type Accel = ();

impl Sites {
    /// Put the columns on the GPUs if there are some and the likelihood is big enough to be worth it (patterns times tips: see
    /// `accel::threshold`); the processor does it otherwise. Says where it will be computed, for the report of the result: where
    /// the GPUs were meant to be used and could not be, that too (and on the standard error, which a job's log keeps).
    fn accelerate(&mut self) -> Option<String> {
        #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
        {
            use crate::accel::{note_failure, threshold, worth_it};
            let cells = self.patterns.len() as f64 * self.tips as f64;
            if worth_it(cells, threshold::likelihood_cells()) {
                match engine::Engine::open(self) {
                    Ok(Some(engine)) => {
                        let place = engine.describe();
                        self.accel = Some(engine);
                        return Some(place);
                    }
                    // `worth_it` says there is a device
                    Ok(None) => {}
                    Err(why) => {
                        note_failure(
                            "the GPUs could not take the likelihood",
                            &format!("{why}; the processor does it"),
                        );
                        return Some(format!("the processor (the GPUs could not be used: {why})"));
                    }
                }
            }
        }
        None
    }

    /// Whether the likelihood is worked out on the processor: there are no GPUs for it, or a card failed and the processor has taken over.
    fn on_processor(&self) -> bool {
        #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
        let on_the_cards = self.accel.as_ref().is_some_and(|e| !e.broken());
        #[cfg(not(all(feature = "gpu", not(target_arch = "wasm32"))))]
        let on_the_cards = false;
        !on_the_cards
    }

    /// Where the likelihood was computed, as `accelerate` said, or that a card failed and the processor finished it.
    fn accelerator(&self, said: Option<String>) -> Option<String> {
        #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
        if let Some(why) = self.accel.as_ref().and_then(engine::Engine::failure) {
            return Some(format!(
                "a GPU failed ({why}), and the processor finished the likelihood"
            ));
        }
        #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
        if let Some(why) = self.accel.as_ref().and_then(engine::Engine::retried) {
            return Some(format!(
                "{}; a GPU failed ({why}) and another GPU did its share",
                said.as_deref().unwrap_or("the GPUs")
            ));
        }
        said
    }

    /// The least a change of topology must raise the log-likelihood by to count as an improvement. On the processor, any rise
    /// above rounding does. On GPUs the likelihood is summed from 32-bit logarithms, whose noise (estimated, not measured on
    /// a CUDA card) grows with the size of the alignment and is a fixed fraction of the likelihood itself; a rise smaller than
    /// that fraction could be noise taken for a better tree. 0.1 log-likelihood units is also the stopping tolerance that
    /// IQ-TREE uses, so little real signal is lost at that size.
    fn smallest_gain(&self, ln_l: f64) -> f64 {
        #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
        if self.accel.as_ref().is_some_and(|e| !e.broken()) {
            return NNI_MIN_GAIN.max(GPU_NOISE_SHARE * ln_l.abs());
        }
        let _ = ln_l;
        NNI_MIN_GAIN
    }
}

fn state(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4,
    }
}

fn sites(msa: &Msa) -> Option<Sites> {
    let mut seen: BTreeMap<Vec<u8>, f64> = BTreeMap::new();
    let mut bases = [1.0f64; 4];
    for col in 0..msa.columns() {
        let pattern: Vec<u8> = msa.rows.iter().map(|r| state(r[col])).collect();
        if pattern.iter().all(|&s| s == 4) {
            continue;
        }
        pattern
            .iter()
            .filter(|&&s| s < 4)
            .for_each(|&s| bases[usize::from(s)] += 1.0);
        *seen.entry(pattern).or_insert(0.0) += 1.0;
    }
    if seen.is_empty() {
        return None;
    }
    let total = seen.values().sum();
    let sum: f64 = bases.iter().sum();
    let (patterns, counts) = seen.into_iter().unzip();
    Some(Sites {
        tips: msa.rows.len(),
        patterns,
        counts,
        total,
        freqs: bases.map(|b| b / sum),
        accel: Accel::default(),
    })
}

/// HKY85 rates scaled to one expected substitution per unit of branch length; kappa is the
/// transition/transversion rate ratio (1 gives F81, equal frequencies as well give JC69).
fn rate_matrix(kappa: f64, freqs: &[f64; 4]) -> Matrix {
    let mut q = [[0.0; 4]; 4];
    for (i, row) in q.iter_mut().enumerate() {
        for (j, rate) in row.iter_mut().enumerate() {
            // A=0, C=1, G=2, T=3: A-G and C-T are transitions (even index sum)
            if i != j {
                *rate = freqs[j] * if (i + j) % 2 == 0 { kappa } else { 1.0 };
            }
        }
        row[i] = -row.iter().sum::<f64>();
    }
    let scale: f64 = -(0..4).map(|i| freqs[i] * q[i][i]).sum::<f64>();
    q.map(|row| row.map(|x| x / scale))
}

fn mat_mul(a: &Matrix, b: &Matrix) -> Matrix {
    let mut out = [[0.0; 4]; 4];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            *cell = (0..4).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    out
}

/// exp(q t) by scaling and squaring.
fn expm(q: &Matrix, t: f64) -> Matrix {
    let norm = q
        .iter()
        .map(|row| row.iter().map(|x| x.abs()).sum::<f64>())
        .fold(0.0, f64::max)
        * t;
    let squarings = if norm > 0.5 {
        (norm / 0.5).log2().ceil() as u32
    } else {
        0
    };
    let step = t / f64::from(1u32 << squarings);
    let a = q.map(|row| row.map(|x| x * step));
    let identity: Matrix =
        std::array::from_fn(|i| std::array::from_fn(|j| f64::from(u8::from(i == j))));
    let (mut term, mut sum) = (identity, identity);
    for k in 1..=TAYLOR_TERMS {
        term = mat_mul(&term, &a).map(|row| row.map(|x| x / f64::from(k)));
        for (s_row, t_row) in sum.iter_mut().zip(&term) {
            s_row.iter_mut().zip(t_row).for_each(|(s, t)| *s += t);
        }
    }
    for _ in 0..squarings {
        sum = mat_mul(&sum, &sum);
    }
    sum
}

/// An unrooted tree: tips are nodes `0..tips` (in alignment order), inner nodes have three
/// neighbours, every edge is stored at both ends with its length.
#[derive(Clone)]
struct Tree {
    adj: Vec<Vec<(usize, f64)>>,
    tips: usize,
}

impl Tree {
    fn from_node(root: &Node, names: &[String]) -> Option<Tree> {
        let index: HashMap<&str, usize> = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i))
            .collect();
        let mut tree = Tree {
            adj: vec![Vec::new(); names.len()],
            tips: names.len(),
        };
        tree.attach(root, &index)?;
        Some(tree)
    }

    fn attach(&mut self, node: &Node, index: &HashMap<&str, usize>) -> Option<usize> {
        if node.is_leaf() {
            return index.get(node.name.as_str()).copied();
        }
        let id = self.adj.len();
        self.adj.push(Vec::new());
        for child in &node.children {
            let c = self.attach(child, index)?;
            self.adj[id].push((c, child.length));
            self.adj[c].push((id, child.length));
        }
        Some(id)
    }

    /// The tree drawn from the first inner node, which has three children.
    fn to_node(&self, names: &[String]) -> Node {
        self.node_at(self.tips, usize::MAX, 0.0, names)
    }

    fn node_at(&self, id: usize, from: usize, length: f64, names: &[String]) -> Node {
        let children: Vec<Node> = if id < self.tips {
            Vec::new()
        } else {
            self.adj[id]
                .iter()
                .filter(|&&(n, _)| n != from)
                .map(|&(n, l)| self.node_at(n, id, l, names))
                .collect()
        };
        Node {
            name: names
                .get(id)
                .filter(|_| id < self.tips)
                .cloned()
                .unwrap_or_default(),
            length,
            support: None,
            children,
        }
    }

    fn edges(&self) -> Vec<(usize, usize)> {
        (0..self.adj.len())
            .flat_map(|a| {
                self.adj[a]
                    .iter()
                    .filter(move |e| a < e.0)
                    .map(move |e| (a, e.0))
            })
            .collect()
    }

    /// Every edge that touches `u` or `v`.
    fn around(&self, u: usize, v: usize) -> Vec<(usize, usize)> {
        let mut out: Vec<(usize, usize)> = [u, v]
            .iter()
            .flat_map(|&a| self.adj[a].iter().map(move |e| (a.min(e.0), a.max(e.0))))
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }

    fn set_length(&mut self, a: usize, b: usize, t: f64) {
        for (from, to) in [(a, b), (b, a)] {
            if let Some(e) = self.adj[from].iter_mut().find(|e| e.0 == to) {
                e.1 = t;
            }
        }
    }

    fn length(&self, a: usize, b: usize) -> f64 {
        self.adj[a].iter().find(|e| e.0 == b).map_or(0.0, |e| e.1)
    }

    /// Swap one subtree on each side of the inner edge `u`-`v`: `which` picks which of v's two
    /// other neighbours trades places with the first other neighbour of u.
    fn swap_across(&mut self, u: usize, v: usize, which: usize) -> bool {
        let x = self.adj[u].iter().map(|e| e.0).find(|&n| n != v);
        let y = self.adj[v]
            .iter()
            .map(|e| e.0)
            .filter(|&n| n != u)
            .nth(which);
        let (Some(x), Some(y)) = (x, y) else {
            return false;
        };
        let (lx, ly) = (self.length(u, x), self.length(v, y));
        self.adj[u].retain(|e| e.0 != x);
        self.adj[x].retain(|e| e.0 != u);
        self.adj[v].retain(|e| e.0 != y);
        self.adj[y].retain(|e| e.0 != v);
        for (a, b, l) in [(u, y, ly), (v, x, lx)] {
            self.adj[a].push((b, l));
            self.adj[b].push((a, l));
        }
        true
    }

    /// Log-likelihood of the alignment, rooted at tip 0 (any root gives the same value).
    fn ln_likelihood(&self, sites: &Sites, q: &Matrix, freqs: &[f64; 4]) -> f64 {
        #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
        if let Some(value) = sites
            .accel
            .as_ref()
            .and_then(|engine| engine.ln_likelihood(self, q, freqs))
        {
            return value;
        }
        let (first, length) = self.adj[0][0];
        let below = self.partial(first, 0, sites, q);
        let p = expm(q, length);
        sites
            .patterns
            .iter()
            .zip(&below)
            .zip(&sites.counts)
            .map(|((pattern, sub), count)| {
                let l: f64 = (0..4)
                    .filter(|&s| pattern[0] == 4 || usize::from(pattern[0]) == s)
                    .map(|s| freqs[s] * (0..4).map(|x| p[s][x] * sub[x]).sum::<f64>())
                    .sum();
                count * l.ln()
            })
            .sum()
    }

    /// Likelihood of the subtree at `node` (away from `from`) for every pattern, given each
    /// state at `node`.
    fn partial(&self, node: usize, from: usize, sites: &Sites, q: &Matrix) -> Vec<[f64; 4]> {
        if node < self.tips {
            return sites
                .patterns
                .iter()
                .map(|p| {
                    std::array::from_fn(|s| {
                        f64::from(u8::from(p[node] == 4 || usize::from(p[node]) == s))
                    })
                })
                .collect();
        }
        let mut out = vec![[1.0; 4]; sites.patterns.len()];
        for &(child, length) in self.adj[node].iter().filter(|e| e.0 != from) {
            let sub = self.partial(child, node, sites, q);
            let p = expm(q, length);
            for (o, s) in out.iter_mut().zip(&sub) {
                for (i, slot) in o.iter_mut().enumerate() {
                    *slot *= (0..4).map(|x| p[i][x] * s[x]).sum::<f64>();
                }
            }
        }
        out
    }
}

/// The `x` in `[lo, hi]` where `f` is largest, for a function with one peak.
fn golden(mut f: impl FnMut(f64) -> f64, mut lo: f64, mut hi: f64, tolerance: f64) -> f64 {
    const R: f64 = 0.618_033_988_749_895;
    let (mut c, mut d) = (hi - R * (hi - lo), lo + R * (hi - lo));
    let (mut fc, mut fd) = (f(c), f(d));
    while hi - lo > tolerance {
        if fc > fd {
            (hi, d, fd) = (d, c, fc);
            c = hi - R * (hi - lo);
            fc = f(c);
        } else {
            (lo, c, fc) = (c, d, fd);
            d = lo + R * (hi - lo);
            fd = f(d);
        }
    }
    (lo + hi) / 2.0
}

/// The log-likelihood as a function of the length of one branch, everything else of the tree held fixed. What lies beyond each end of
/// the branch is worked out once, as a likelihood for each state at that end and each pattern, and with the stationary frequency of the
/// state at one end multiplied in, the 16 products of a pattern are all that a length tried costs. The pruning of the whole tree that
/// `Tree::ln_likelihood` does for every length (about twenty of them for a branch) is done twice for the branch, as one pruning. The same
/// likelihood to rounding (the sums are made in another order), whatever the branch.
struct EdgeScan {
    /// For each pattern, `freqs[s]` times the likelihood beyond one end given state `s` there, times the likelihood beyond the other end
    /// given state `t` there, at `4 s + t`.
    weights: Vec<[f64; 16]>,
}

impl EdgeScan {
    fn new(
        tree: &Tree,
        a: usize,
        b: usize,
        sites: &Sites,
        q: &Matrix,
        freqs: &[f64; 4],
    ) -> EdgeScan {
        let beyond_a = tree.partial(a, b, sites, q);
        let beyond_b = tree.partial(b, a, sites, q);
        let weights = beyond_a
            .iter()
            .zip(&beyond_b)
            .map(|(x, y)| std::array::from_fn(|k| freqs[k / 4] * x[k / 4] * y[k % 4]))
            .collect();
        EdgeScan { weights }
    }

    /// The log-likelihood of the alignment with the branch `length` long.
    fn ln_likelihood(&self, sites: &Sites, q: &Matrix, length: f64) -> f64 {
        let p = expm(q, length);
        self.weights
            .iter()
            .zip(&sites.counts)
            .map(|(w, count)| {
                let l: f64 = w.iter().zip(p.iter().flatten()).map(|(w, p)| w * p).sum();
                count * l.ln()
            })
            .sum()
    }
}

/// The best length of the branch `a`-`b`, tried on the processor through its `EdgeScan`.
fn best_length_on_processor(
    tree: &Tree,
    (a, b): (usize, usize),
    sites: &Sites,
    q: &Matrix,
    freqs: &[f64; 4],
) -> f64 {
    let scan = EdgeScan::new(tree, a, b, sites, q, freqs);
    golden(
        |x| scan.ln_likelihood(sites, q, x.exp()),
        MIN_BRANCH.ln(),
        MAX_BRANCH.ln(),
        LOG_TOLERANCE,
    )
}

/// The same with the whole tree's likelihood for every length (what the GPUs do where they do not keep the partial likelihoods).
fn best_length_by_whole_tree(
    tree: &mut Tree,
    (a, b): (usize, usize),
    sites: &Sites,
    q: &Matrix,
    freqs: &[f64; 4],
) -> f64 {
    golden(
        |x| {
            tree.set_length(a, b, x.exp());
            tree.ln_likelihood(sites, q, freqs)
        },
        MIN_BRANCH.ln(),
        MAX_BRANCH.ln(),
        LOG_TOLERANCE,
    )
}

/// The same on the GPUs, where they keep what lies beyond both ends of the branch (PyTorch): `None` where they do not, or a card failed
/// on the way (and the processor has taken over).
#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
fn best_length_on_cards(
    tree: &Tree,
    (a, b): (usize, usize),
    sites: &Sites,
    q: &Matrix,
    freqs: &[f64; 4],
) -> Option<f64> {
    let scan = sites.accel.as_ref()?.edge_scan(tree, a, b, q, freqs)?;
    let mut failed = false;
    let best = golden(
        |x| {
            scan.ln_likelihood(q, x.exp()).unwrap_or_else(|| {
                failed = true;
                f64::NEG_INFINITY
            })
        },
        MIN_BRANCH.ln(),
        MAX_BRANCH.ln(),
        LOG_TOLERANCE,
    );
    (!failed).then_some(best)
}

#[cfg(not(all(feature = "gpu", not(target_arch = "wasm32"))))]
fn best_length_on_cards(
    _: &Tree,
    _: (usize, usize),
    _: &Sites,
    _: &Matrix,
    _: &[f64; 4],
) -> Option<f64> {
    None
}

/// Set each listed branch to its best length, `sweeps` times round.
fn optimize_edges(
    tree: &mut Tree,
    edges: &[(usize, usize)],
    sweeps: usize,
    sites: &Sites,
    q: &Matrix,
    freqs: &[f64; 4],
) {
    for _ in 0..sweeps {
        for &edge in edges {
            let on_cards = if sites.on_processor() {
                None
            } else {
                best_length_on_cards(tree, edge, sites, q, freqs)
            };
            let best = match on_cards {
                Some(best) => best,
                // the processor, or a card failed on the way and the processor has taken over
                None if sites.on_processor() => {
                    best_length_on_processor(tree, edge, sites, q, freqs)
                }
                None => best_length_by_whole_tree(tree, edge, sites, q, freqs),
            };
            tree.set_length(edge.0, edge.1, best.exp());
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Fit {
    pub kind: Kind,
    pub ln_l: f64,
    pub bic: f64,
    pub kappa: f64,
    pub freqs: [f64; 4],
}

impl Fit {
    pub fn describe(&self) -> String {
        format!(
            "{}: ln L {:.2}, BIC {:.2}{}",
            self.kind.label(),
            self.ln_l,
            self.bic,
            if self.kind.free_kappa() {
                format!(", kappa {:.2}", self.kappa)
            } else {
                String::new()
            }
        )
    }
}

/// Fit one model: branch lengths (and kappa) for the topology of `tree`.
fn fit(kind: Kind, tree: &mut Tree, sites: &Sites) -> Fit {
    let freqs = if kind.free_freqs() {
        sites.freqs
    } else {
        [0.25; 4]
    };
    let q_for = |kappa: f64| rate_matrix(kappa, &freqs);
    let edges = tree.edges();
    let mut kappa = 1.0;
    optimize_edges(tree, &edges, 2, sites, &q_for(kappa), &freqs);
    if kind.free_kappa() {
        for _ in 0..2 {
            let x = golden(
                |x| tree.ln_likelihood(sites, &q_for(x.exp()), &freqs),
                KAPPA_RANGE.0.ln(),
                KAPPA_RANGE.1.ln(),
                LOG_TOLERANCE,
            );
            kappa = x.exp();
            optimize_edges(tree, &edges, 1, sites, &q_for(kappa), &freqs);
        }
    }
    let ln_l = tree.ln_likelihood(sites, &q_for(kappa), &freqs);
    let free = (2 * sites.tips - 3)
        + usize::from(kind.free_kappa())
        + if kind.free_freqs() { 3 } else { 0 };
    Fit {
        kind,
        ln_l,
        bic: free as f64 * sites.total.ln() - 2.0 * ln_l,
        kappa,
        freqs,
    }
}

/// Hill-climb over nearest-neighbour interchanges: take any swap that raises the likelihood
/// (after re-fitting the five branches around it), until none does.
fn improve_topology(tree: &mut Tree, sites: &Sites, q: &Matrix, freqs: &[f64; 4]) {
    let mut best = tree.ln_likelihood(sites, q, freqs);
    for _ in 0..NNI_MOVES {
        let mut moved = false;
        'edges: for (u, v) in tree.edges() {
            if u < tree.tips || v < tree.tips {
                continue;
            }
            for which in 0..2 {
                let mut candidate = tree.clone();
                if !candidate.swap_across(u, v, which) {
                    continue;
                }
                let local = candidate.around(u, v);
                optimize_edges(&mut candidate, &local, 1, sites, q, freqs);
                let l = candidate.ln_likelihood(sites, q, freqs);
                if l > best + sites.smallest_gain(best) {
                    (*tree, best, moved) = (candidate, l, true);
                    break 'edges;
                }
            }
        }
        if !moved {
            break;
        }
    }
}

/// Fit every model to a tree given from outside (another program's, say) the way the search fits its starting
/// tree: branch lengths and kappa. The tips are matched to the alignment by name; the tree needs three
/// branches at its first inner node, as unrooted trees are written.
// ponytail: only the tests compare other programs' trees this way.
#[cfg(test)]
pub fn fit_models(msa: &Msa, topology: &Node) -> Result<Vec<Fit>, String> {
    let sites = sites(msa).ok_or("No column of the alignment holds a base.")?;
    let tree = Tree::from_node(topology, &msa.names)
        .ok_or("The tree's tips do not match the alignment's names.")?;
    Ok(par_map(&Kind::ALL, |&kind| {
        let mut copy = tree.clone();
        fit(kind, &mut copy, &sites)
    }))
}

pub struct MlResult {
    pub tree: Node,
    pub best: Fit,
    /// Every model tried on the starting tree, in the order JC69, K80, F81, HKY85.
    pub fits: Vec<Fit>,
    /// Where the likelihood was computed when it was not on the processor: the GPUs, with the share of each.
    pub accelerator: Option<String>,
}

/// The maximum-likelihood tree of an alignment of at least three sequences.
pub fn ml_tree(msa: &Msa) -> Result<MlResult, String> {
    ml_tree_using(msa, Sites::accelerate)
}

/// The same, with the putting of the columns on the GPUs (or not) done by `accelerate`: what the tests use to try a pool of their own.
fn ml_tree_using(
    msa: &Msa,
    accelerate: impl FnOnce(&mut Sites) -> Option<String>,
) -> Result<MlResult, String> {
    if msa.rows.len() < 3 {
        return Err("A likelihood tree needs at least three sequences.".into());
    }
    let mut sites = sites(msa).ok_or("No column of the alignment holds a base.")?;
    let start = [
        DistanceModel::Kimura2,
        DistanceModel::JukesCantor,
        DistanceModel::P,
    ]
    .iter()
    .find_map(|&m| msa_distances(msa, m).ok())
    .ok_or("The sequences are too different to start a tree from.")?;
    let start_tree = neighbor_joining(&msa.names, &start)
        .and_then(|t| Tree::from_node(&t, &msa.names))
        .ok_or("The starting tree could not be built.")?;
    // the columns go to the GPUs after the starting tree's distances, which have had theirs and let them go: the two are not on a
    // card together
    let placed = accelerate(&mut sites);
    // the four models are fitted side by side, each on its own copy of the starting tree
    let fitted = par_map(&Kind::ALL, |&kind| {
        let mut tree = start_tree.clone();
        (fit(kind, &mut tree, &sites), tree)
    });
    let fits: Vec<Fit> = fitted.iter().map(|(f, _)| f.clone()).collect();
    let (chosen, mut tree) = fitted
        .into_iter()
        .min_by(|a, b| a.0.bic.total_cmp(&b.0.bic))
        .ok_or("No model could be fitted.")?;
    let q = rate_matrix(chosen.kappa, &chosen.freqs);
    improve_topology(&mut tree, &sites, &q, &chosen.freqs);
    let edges = tree.edges();
    optimize_edges(&mut tree, &edges, 2, &sites, &q, &chosen.freqs);
    let ln_l = tree.ln_likelihood(&sites, &q, &chosen.freqs);
    let free = (2 * sites.tips - 3)
        + usize::from(chosen.kind.free_kappa())
        + if chosen.kind.free_freqs() { 3 } else { 0 };
    let best = Fit {
        ln_l,
        bic: free as f64 * sites.total.ln() - 2.0 * ln_l,
        ..chosen
    };
    Ok(MlResult {
        tree: tree.to_node(&msa.names),
        best,
        fits,
        accelerator: sites.accelerator(placed),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::tree::{splits, Rng};

    /// A tree from its edges; an unrooted tree has one node more than it has edges.
    fn tree_of(tips: usize, edges: &[(usize, usize, f64)]) -> Tree {
        let mut adj = vec![Vec::new(); edges.len() + 1];
        for &(a, b, l) in edges {
            adj[a].push((b, l));
            adj[b].push((a, l));
        }
        Tree { adj, tips }
    }

    fn jc() -> (Matrix, [f64; 4]) {
        (rate_matrix(1.0, &[0.25; 4]), [0.25; 4])
    }

    #[test]
    fn transition_probabilities_follow_jukes_cantor() {
        let (q, _) = jc();
        let p = expm(&q, 0.3);
        assert!((p[0][0] - 0.752_740_034_5).abs() < 1e-9);
        assert!((p[0][1] - 0.082_419_988_5).abs() < 1e-9);
        assert!(p
            .iter()
            .all(|row| (row.iter().sum::<f64>() - 1.0).abs() < 1e-12));
        let long = expm(&q, 8.0);
        assert!(
            (long[2][1] - 0.25).abs() < 1e-4,
            "long branches forget their start"
        );
    }

    #[test]
    fn hky_rates_favour_transitions_and_keep_the_stationary_frequencies() {
        let freqs = [0.1, 0.2, 0.3, 0.4];
        let p = expm(&rate_matrix(4.0, &freqs), 0.2);
        assert!(p[0][2] > p[0][1], "A to G (a transition) beats A to C");
        // pi P(t) = pi
        for j in 0..4 {
            let sum: f64 = (0..4).map(|i| freqs[i] * p[i][j]).sum();
            assert!((sum - freqs[j]).abs() < 1e-10);
        }
    }

    // Reference value from a brute-force sum over the inner states, computed separately in Python.
    #[test]
    fn the_likelihood_matches_a_brute_force_sum() {
        // tips 0..=3, inner nodes 4 (with tips 0, 1) and 5 (with tips 2, 3)
        let tree = tree_of(
            4,
            &[
                (0, 4, 0.1),
                (1, 4, 0.2),
                (4, 5, 0.05),
                (2, 5, 0.3),
                (3, 5, 0.15),
            ],
        );
        let sites = Sites {
            tips: 4,
            patterns: vec![vec![0, 0, 0, 0], vec![0, 0, 1, 1], vec![0, 2, 1, 3]],
            counts: vec![10.0, 3.0, 1.0],
            total: 14.0,
            freqs: [0.25; 4],
            accel: Default::default(),
        };
        let (q, freqs) = jc();
        let got = tree.ln_likelihood(&sites, &q, &freqs);
        assert!((got - -47.852_488_656_3).abs() < 1e-8, "{got}");
        // the same tree with the tips numbered the other way round is rooted at another tip
        let reversed = tree_of(
            4,
            &[
                (3, 4, 0.1),
                (2, 4, 0.2),
                (4, 5, 0.05),
                (1, 5, 0.3),
                (0, 5, 0.15),
            ],
        );
        let flipped = Sites {
            patterns: sites
                .patterns
                .iter()
                .map(|p| p.iter().rev().copied().collect())
                .collect(),
            counts: sites.counts.clone(),
            ..sites
        };
        assert!((reversed.ln_likelihood(&flipped, &q, &freqs) - got).abs() < 1e-9);
    }

    #[test]
    fn a_gap_fits_any_base() {
        let tree = tree_of(3, &[(0, 3, 0.1), (1, 3, 0.1), (2, 3, 0.1)]);
        let (q, freqs) = jc();
        let mk = |last: u8| Sites {
            tips: 3,
            patterns: vec![vec![0, 0, last]],
            counts: vec![1.0],
            total: 1.0,
            freqs: [0.25; 4],
            accel: Default::default(),
        };
        let with_gap = tree.ln_likelihood(&mk(4), &q, &freqs);
        let two_only = Tree {
            adj: vec![vec![(1, 0.2)], vec![(0, 0.2)]],
            tips: 2,
        }
        .ln_likelihood(
            &Sites {
                tips: 2,
                patterns: vec![vec![0, 0]],
                counts: vec![1.0],
                total: 1.0,
                freqs: [0.25; 4],
                accel: Default::default(),
            },
            &q,
            &freqs,
        );
        assert!(
            (with_gap - two_only).abs() < 1e-9,
            "{with_gap} vs {two_only}"
        );
    }

    #[test]
    fn golden_section_finds_the_peak() {
        let x = golden(|x| -(x - 1.7).powi(2), -5.0, 5.0, 1e-6);
        assert!((x - 1.7).abs() < 1e-4);
    }

    /// Four sequences evolved under JC69 along ((A,B),(C,D)): terminal branches 0.05, inner 0.2.
    pub(super) fn simulated() -> Msa {
        let mut rng = Rng::new(2024);
        let root: Vec<u8> = (0..600).map(|_| b"ACGT"[rng.below(4)]).collect();
        let mut evolve = |from: &[u8], t: f64| -> Vec<u8> {
            let change = 0.75 * (1.0 - (-4.0 * t / 3.0).exp());
            from.iter()
                .map(|&c| {
                    if rng.unit() < change {
                        let others: Vec<u8> = b"ACGT".iter().copied().filter(|&b| b != c).collect();
                        others[rng.below(3)]
                    } else {
                        c
                    }
                })
                .collect()
        };
        let (x, y) = (evolve(&root, 0.1), evolve(&root, 0.1));
        let rows = vec![
            evolve(&x, 0.05),
            evolve(&x, 0.05),
            evolve(&y, 0.05),
            evolve(&y, 0.05),
        ];
        Msa {
            names: ["A", "B", "C", "D"].map(String::from).to_vec(),
            rows,
        }
    }

    #[test]
    fn a_swap_across_an_inner_edge_changes_the_topology_and_the_search_undoes_a_wrong_one() {
        let msa = simulated();
        let sites = sites(&msa).unwrap();
        // wrong start: (A,C) with (B,D)
        let mut tree = tree_of(
            4,
            &[
                (0, 4, 0.1),
                (2, 4, 0.1),
                (4, 5, 0.1),
                (1, 5, 0.1),
                (3, 5, 0.1),
            ],
        );
        let (q, freqs) = jc();
        let edges = tree.edges();
        optimize_edges(&mut tree, &edges, 2, &sites, &q, &freqs);
        let wrong = tree.ln_likelihood(&sites, &q, &freqs);
        improve_topology(&mut tree, &sites, &q, &freqs);
        assert!(tree.ln_likelihood(&sites, &q, &freqs) > wrong + 1.0);
        let names = msa.names.clone();
        let found = splits(&tree.to_node(&names), &names);
        assert!(
            found.contains(&vec![2, 3]),
            "the A,B | C,D split: {found:?}"
        );
    }

    #[test]
    fn the_full_search_returns_the_right_tree_and_a_model_by_bic() {
        let msa = simulated();
        let result = ml_tree(&msa).unwrap();
        let found = splits(&result.tree, &msa.names);
        assert!(found.contains(&vec![2, 3]), "{found:?}");
        assert_eq!(result.fits.len(), 4);
        let lowest = result
            .fits
            .iter()
            .map(|f| f.bic)
            .fold(f64::INFINITY, f64::min);
        assert!(
            result.best.bic <= lowest + 5.0,
            "re-fitting after the search only improves the fit"
        );
        // JC69 is nested in HKY85, so the richer model never fits worse
        assert!(result.fits[3].ln_l >= result.fits[0].ln_l - 0.05);
        // the inner branch is the longest, as simulated
        let inner = result.tree.children.iter().find(|c| !c.is_leaf());
        assert!(
            inner.is_some_and(|c| c.length > 0.1 && c.length < 0.35),
            "{:?}",
            inner.map(|c| c.length)
        );
        assert!(ml_tree(&Msa {
            names: vec!["a".into(), "b".into()],
            rows: vec![b"AC".to_vec(); 2]
        })
        .is_err());
    }

    // ---- the likelihood of one branch's length alone, against the whole tree's pruning for every length ----

    /// Setting each listed branch to its best length by trying the whole tree for every length: the way it was done before `EdgeScan`.
    fn optimize_edges_reference(
        tree: &mut Tree,
        edges: &[(usize, usize)],
        sweeps: usize,
        sites: &Sites,
        q: &Matrix,
        freqs: &[f64; 4],
    ) {
        for _ in 0..sweeps {
            for &(a, b) in edges {
                let best = golden(
                    |x| {
                        tree.set_length(a, b, x.exp());
                        tree.ln_likelihood(sites, q, freqs)
                    },
                    MIN_BRANCH.ln(),
                    MAX_BRANCH.ln(),
                    LOG_TOLERANCE,
                );
                tree.set_length(a, b, best.exp());
            }
        }
    }

    /// The search of topologies as it was, with the reference's way of setting branches.
    fn improve_topology_reference(tree: &mut Tree, sites: &Sites, q: &Matrix, freqs: &[f64; 4]) {
        let mut best = tree.ln_likelihood(sites, q, freqs);
        for _ in 0..NNI_MOVES {
            let mut moved = false;
            'edges: for (u, v) in tree.edges() {
                if u < tree.tips || v < tree.tips {
                    continue;
                }
                for which in 0..2 {
                    let mut candidate = tree.clone();
                    if !candidate.swap_across(u, v, which) {
                        continue;
                    }
                    let local = candidate.around(u, v);
                    optimize_edges_reference(&mut candidate, &local, 1, sites, q, freqs);
                    let l = candidate.ln_likelihood(sites, q, freqs);
                    if l > best + sites.smallest_gain(best) {
                        (*tree, best, moved) = (candidate, l, true);
                        break 'edges;
                    }
                }
            }
            if !moved {
                break;
            }
        }
    }

    /// A random unrooted tree over `tips` tips (each tip added in the middle of a random branch) with random branch lengths, and an
    /// alignment of `columns` columns evolved along it under JC69, with a few gaps and ambiguity codes in it.
    fn random_case(tips: usize, columns: usize, seed: u64) -> (Tree, Msa) {
        let mut rng = Rng::new(seed);
        let mut edges: Vec<(usize, usize, f64)> = vec![
            (0, tips, 0.05 + 0.2 * rng.unit()),
            (1, tips, 0.05 + 0.2 * rng.unit()),
            (2, tips, 0.05 + 0.2 * rng.unit()),
        ];
        for k in 3..tips {
            let at = rng.below(edges.len());
            let (a, b, l) = edges[at];
            let middle = tips + k - 2;
            edges[at] = (a, middle, l / 2.0);
            edges.push((middle, b, l / 2.0));
            edges.push((k, middle, 0.02 + 0.2 * rng.unit()));
        }
        let tree = tree_of(tips, &edges);
        let mut rows: Vec<Vec<u8>> = vec![Vec::new(); tips];
        let root: Vec<u8> = (0..columns).map(|_| b"ACGT"[rng.below(4)]).collect();
        // walk out from tip 0, evolving the sequence along each branch
        let mut stack = vec![(0usize, usize::MAX, root)];
        while let Some((node, from, seq)) = stack.pop() {
            if node < tips {
                rows[node] = seq.clone();
            }
            for &(next, length) in tree.adj[node].iter().filter(|e| e.0 != from) {
                let change = 0.75 * (1.0 - (-4.0 * length / 3.0).exp());
                let evolved = seq
                    .iter()
                    .map(|&c| {
                        if rng.unit() < change {
                            let others: Vec<u8> =
                                b"ACGT".iter().copied().filter(|&b| b != c).collect();
                            others[rng.below(3)]
                        } else {
                            c
                        }
                    })
                    .collect();
                stack.push((next, node, evolved));
            }
        }
        for row in rows.iter_mut().skip(1) {
            for c in row.iter_mut() {
                match rng.below(60) {
                    0 => *c = b'N',
                    1 => *c = b'-',
                    2 => *c = b'R',
                    _ => {}
                }
            }
        }
        let names = (0..tips).map(|k| format!("t{k}")).collect();
        (tree, Msa { names, rows })
    }

    /// A tree whose inner nodes form a chain, the tips hung on it in order, every branch 0.1: a start that the search must change.
    fn caterpillar(tips: usize) -> Tree {
        let mut edges = vec![(0, tips, 0.1), (1, tips, 0.1)];
        for k in 2..tips - 1 {
            edges.push((tips + k - 2, tips + k - 1, 0.1));
            edges.push((k, tips + k - 1, 0.1));
        }
        edges.push((tips - 2, 2 * tips - 3, 0.1));
        edges.push((tips - 1, 2 * tips - 3, 0.1));
        tree_of(tips, &edges)
    }

    #[test]
    fn the_scan_of_one_branch_gives_the_likelihood_of_the_whole_tree_for_every_length() {
        for (tips, columns, seed) in [
            (3usize, 120usize, 1u64),
            (4, 200, 2),
            (7, 300, 3),
            (12, 150, 4),
        ] {
            let (tree, msa) = random_case(tips, columns, seed);
            let sites = sites(&msa).unwrap();
            let q = rate_matrix(2.3, &sites.freqs);
            for (a, b) in tree.edges() {
                let scan = EdgeScan::new(&tree, a, b, &sites, &q, &sites.freqs);
                for length in [MIN_BRANCH, 0.01, 0.1, 0.7, MAX_BRANCH] {
                    let mut set = tree.clone();
                    set.set_length(a, b, length);
                    let whole = set.ln_likelihood(&sites, &q, &sites.freqs);
                    let one = scan.ln_likelihood(&sites, &q, length);
                    assert!(
                        (one - whole).abs() <= 1e-9 * whole.abs().max(1.0),
                        "{tips} tips, branch {a}-{b} at {length}: {one} against {whole}"
                    );
                }
            }
        }
    }

    #[test]
    fn setting_the_branches_by_the_scan_reaches_the_likelihood_that_trying_the_whole_tree_does() {
        for (tips, columns, seed) in [(4usize, 300usize, 5u64), (8, 400, 6), (13, 250, 7)] {
            let (truth, msa) = random_case(tips, columns, seed);
            let sites = sites(&msa).unwrap();
            let q = rate_matrix(2.0, &sites.freqs);
            // from branches that are all wrong, as the search's candidates are
            let flat = tree_of(
                tips,
                &truth
                    .edges()
                    .iter()
                    .map(|&(a, b)| (a, b, 0.1))
                    .collect::<Vec<_>>(),
            );
            let edges = flat.edges();
            let (mut by_scan, mut whole) = (flat.clone(), flat);
            optimize_edges(&mut by_scan, &edges, 2, &sites, &q, &sites.freqs);
            optimize_edges_reference(&mut whole, &edges, 2, &sites, &q, &sites.freqs);
            let (x, y) = (
                by_scan.ln_likelihood(&sites, &q, &sites.freqs),
                whole.ln_likelihood(&sites, &q, &sites.freqs),
            );
            assert!((x - y).abs() < 1e-6, "{tips} tips: {x} against {y}");
        }
    }

    #[test]
    fn from_a_tree_one_swap_off_the_search_with_the_scan_puts_it_right_as_the_search_without_it_does(
    ) {
        // (A greedy search takes the first change that helps, so two ways of working out a likelihood that differ by 1e-4 can take
        // different paths from a bad start, and end on trees that are equally good by luck: over 40 caterpillar starts the final
        // likelihoods were the same in 21, higher with the scan in 11 and lower in 8, on average 12.8 higher. From one swap off the
        // truth there is one way back.)
        for (tips, columns, seed) in [
            (7usize, 600usize, 21u64),
            (8, 700, 22),
            (9, 700, 23),
            (10, 800, 24),
        ] {
            let (truth, msa) = random_case(tips, columns, seed);
            let sites = sites(&msa).unwrap();
            let q = rate_matrix(1.5, &sites.freqs);
            let names = msa.names.clone();
            let edges = truth.edges();
            let mut right = truth.clone();
            optimize_edges_reference(&mut right, &edges, 6, &sites, &q, &sites.freqs);
            // one swap across the first inner branch that has one
            let mut start = truth.clone();
            let swapped = edges
                .iter()
                .filter(|&&(u, v)| u >= tips && v >= tips)
                .any(|&(u, v)| start.swap_across(u, v, 0));
            assert!(swapped);
            let (mut fast, mut slow) = (start.clone(), start);
            optimize_edges(&mut fast, &edges, 2, &sites, &q, &sites.freqs);
            optimize_edges_reference(&mut slow, &edges, 2, &sites, &q, &sites.freqs);
            improve_topology(&mut fast, &sites, &q, &sites.freqs);
            improve_topology_reference(&mut slow, &sites, &q, &sites.freqs);
            assert_eq!(
                splits(&fast.to_node(&names), &names),
                splits(&slow.to_node(&names), &names),
                "{tips} tips"
            );
            let edges = fast.edges();
            optimize_edges_reference(&mut fast, &edges, 8, &sites, &q, &sites.freqs);
            optimize_edges_reference(&mut slow, &edges, 8, &sites, &q, &sites.freqs);
            let (x, y) = (
                fast.ln_likelihood(&sites, &q, &sites.freqs),
                slow.ln_likelihood(&sites, &q, &sites.freqs),
            );
            assert!(
                (x - y).abs() < 1e-6,
                "{tips} tips, converged: {x} against {y}"
            );
            // and the search got somewhere: at least as high as the truth, which the data need not favour
            assert!(
                x >= right.ln_likelihood(&sites, &q, &sites.freqs) - 0.5,
                "{tips} tips"
            );
        }
    }

    #[test]
    #[ignore = "timing: cargo test --release timing_of_the_topology_search -- --ignored --nocapture"]
    fn timing_of_the_topology_search() {
        for (tips, columns) in [(10usize, 600usize), (20, 1000), (40, 1500)] {
            let (_, msa) = random_case(tips, columns, 99);
            let sites = sites(&msa).unwrap();
            let q = rate_matrix(2.0, &sites.freqs);
            let start = caterpillar(tips);
            let time = |reference: bool| {
                let mut tree = start.clone();
                let begin = std::time::Instant::now();
                if reference {
                    improve_topology_reference(&mut tree, &sites, &q, &sites.freqs);
                } else {
                    improve_topology(&mut tree, &sites, &q, &sites.freqs);
                }
                (
                    begin.elapsed().as_secs_f64(),
                    tree.ln_likelihood(&sites, &q, &sites.freqs),
                )
            };
            let (slow, slow_l) = time(true);
            let (fast, fast_l) = time(false);
            println!(
                "{tips} tips x {columns} columns: {slow:.2} s before, {fast:.2} s now ({:.1} times), ln L {slow_l:.3} and {fast_l:.3}",
                slow / fast
            );
        }
    }
}
