//! Bayesian trees in the manner of MrBayes, much smaller: a Metropolis-Hastings sampler over the
//! tree topology, every branch length and the transition/transversion ratio under HKY85 with the
//! alignment's base frequencies. Priors: all topologies equally likely, an exponential prior (mean
//! 0.1, MrBayes' default exp(10)) on each branch length, uniform in log on kappa. Independent runs
//! go on separate cores and their agreement is reported as the average standard deviation of
//! split frequencies. Not included: Metropolis coupling, relaxed clocks, tip dating, the
//! coalescent and the other BEAST2 and MrBayes models. Written from the published method; no
//! MrBayes or BEAST2 code was used.

use super::{rate_matrix, sites, Sites, Tree, KAPPA_RANGE, MAX_BRANCH, MIN_BRANCH};
use crate::seq::msa::Msa;
use crate::seq::par::par_map;
use crate::seq::tree::{msa_distances, neighbor_joining, DistanceModel, Node, Rng};
use std::collections::{HashMap, HashSet};

/// Mean of the exponential prior on every branch length.
const BRANCH_PRIOR_MEAN: f64 = 0.1;
/// Width of the multiplier windows: a branch length (or kappa) changes by a factor within e^±(step/2).
const BRANCH_STEP: f64 = 2.0;
const KAPPA_STEP: f64 = 1.0;
/// Share of moves that change a branch length, then that swap subtrees; the rest change kappa.
const SHARE_BRANCH: f64 = 0.6;
const SHARE_TOPOLOGY: f64 = 0.9;
/// Random swaps that scatter the starting tree of every run after the first, per run.
const SCATTER_SWAPS: usize = 12;
/// Splits seen in at least this share of a run's samples count towards the convergence measure.
const ASDSF_MIN_FREQUENCY: f64 = 0.1;
/// Splits with a posterior probability above this make the consensus tree.
const CONSENSUS_ABOVE: f64 = 0.5;
/// Below this average standard deviation of split frequencies the runs agree (MrBayes' guideline).
pub const CONVERGED_BELOW: f64 = 0.01;
/// Leaves are bits of a u128, with one bit spare.
const MAX_TAXA: usize = 127;

#[derive(Clone, Debug)]
pub struct Settings {
    pub generations: usize,
    /// A tree is sampled every this many generations after the burn-in.
    pub sample_every: usize,
    /// Share of the generations thrown away at the start.
    pub burn_in: f64,
    /// Independent runs, one per core.
    pub runs: usize,
    pub seed: u64,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            generations: 10_000,
            sample_every: 10,
            burn_in: 0.25,
            runs: 2,
            seed: 2024,
        }
    }
}

pub struct BayesResult {
    /// The majority-rule consensus tree: support is the posterior probability in percent, branch
    /// lengths are means over the sampled trees that have the branch.
    pub tree: Node,
    /// Average standard deviation of split frequencies between runs; `None` for a single run.
    pub asdsf: Option<f64>,
    pub samples: usize,
    pub mean_ln_l: f64,
    pub mean_kappa: f64,
}

struct Model<'a> {
    sites: &'a Sites,
}

impl Model<'_> {
    fn ln_l(&self, tree: &Tree, kappa: f64) -> f64 {
        let q = rate_matrix(kappa, &self.sites.freqs);
        tree.ln_likelihood(self.sites, &q, &self.sites.freqs)
    }
}

fn inner_edges(tree: &Tree) -> Vec<(usize, usize)> {
    tree.edges()
        .into_iter()
        .filter(|&(a, b)| a >= tree.tips && b >= tree.tips)
        .collect()
}

struct Chain {
    tree: Tree,
    kappa: f64,
    ln_l: f64,
    rng: Rng,
}

impl Chain {
    fn new(start: &Tree, model: &Model, seed: u64, scatter: usize) -> Chain {
        let (mut rng, mut tree) = (Rng::new(seed), start.clone());
        for _ in 0..scatter {
            let edges = inner_edges(&tree);
            if edges.is_empty() {
                break;
            }
            let (u, v) = edges[rng.below(edges.len())];
            tree.swap_across(u, v, rng.below(2));
        }
        let kappa = 1.0;
        let ln_l = model.ln_l(&tree, kappa);
        Chain {
            tree,
            kappa,
            ln_l,
            rng,
        }
    }

    /// One Metropolis-Hastings move.
    fn generation(&mut self, model: &Model) {
        let pick = self.rng.unit();
        if pick < SHARE_BRANCH {
            self.change_length(model);
        } else if pick < SHARE_TOPOLOGY {
            self.swap_subtrees(model);
        } else {
            self.change_kappa(model);
        }
    }

    fn accept(&mut self, ln_ratio: f64) -> bool {
        self.rng.unit().ln() < ln_ratio
    }

    fn change_length(&mut self, model: &Model) {
        let edges = self.tree.edges();
        let (a, b) = edges[self.rng.below(edges.len())];
        let old = self.tree.length(a, b);
        let new = old * (BRANCH_STEP * (self.rng.unit() - 0.5)).exp();
        if !(MIN_BRANCH..=MAX_BRANCH).contains(&new) {
            return;
        }
        self.tree.set_length(a, b, new);
        let ln_l = model.ln_l(&self.tree, self.kappa);
        // likelihood ratio, exponential prior ratio, and the multiplier move's proposal ratio
        let ln_ratio = ln_l - self.ln_l - (new - old) / BRANCH_PRIOR_MEAN + (new / old).ln();
        if self.accept(ln_ratio) {
            self.ln_l = ln_l;
        } else {
            self.tree.set_length(a, b, old);
        }
    }

    fn swap_subtrees(&mut self, model: &Model) {
        let edges = inner_edges(&self.tree);
        if edges.is_empty() {
            return; // three sequences have one topology
        }
        let (u, v) = edges[self.rng.below(edges.len())];
        let before = self.tree.clone();
        if !self.tree.swap_across(u, v, self.rng.below(2)) {
            return;
        }
        let ln_l = model.ln_l(&self.tree, self.kappa);
        // every topology has n - 3 inner edges with two swaps each, so the proposal is symmetric
        if self.accept(ln_l - self.ln_l) {
            self.ln_l = ln_l;
        } else {
            self.tree = before;
        }
    }

    fn change_kappa(&mut self, model: &Model) {
        let new = self.kappa * (KAPPA_STEP * (self.rng.unit() - 0.5)).exp();
        if !(KAPPA_RANGE.0..=KAPPA_RANGE.1).contains(&new) {
            return;
        }
        let ln_l = model.ln_l(&self.tree, new);
        // a uniform prior on ln(kappa) and a multiplier move cancel, leaving the likelihood ratio
        if self.accept(ln_l - self.ln_l) {
            (self.kappa, self.ln_l) = (new, ln_l);
        }
    }
}

/// What a run saw in the trees it sampled.
#[derive(Default)]
struct Tally {
    samples: f64,
    /// Inner splits (the leaves on the side without leaf 0, one bit each): times seen and summed
    /// branch length.
    splits: HashMap<u128, (f64, f64)>,
    /// Summed length of each tip's own branch.
    tips: Vec<f64>,
    kappa: f64,
    ln_l: f64,
}

/// The leaves beyond `node`, looking away from `from`, one bit each.
fn leaves_beyond(tree: &Tree, node: usize, from: usize) -> u128 {
    if node < tree.tips {
        return 1u128 << node;
    }
    tree.adj[node]
        .iter()
        .filter(|e| e.0 != from)
        .fold(0, |mask, &(next, _)| mask | leaves_beyond(tree, next, node))
}

impl Tally {
    fn record(&mut self, tree: &Tree, kappa: f64, ln_l: f64) {
        let everything = (1u128 << tree.tips) - 1;
        self.samples += 1.0;
        self.kappa += kappa;
        self.ln_l += ln_l;
        for (a, b) in tree.edges() {
            let length = tree.length(a, b);
            if a < tree.tips {
                self.tips[a] += length;
                continue;
            }
            let beyond = leaves_beyond(tree, b, a);
            let key = if beyond & 1 == 1 {
                everything ^ beyond
            } else {
                beyond
            };
            let seen = self.splits.entry(key).or_insert((0.0, 0.0));
            (seen.0, seen.1) = (seen.0 + 1.0, seen.1 + length);
        }
    }
}

fn pooled(runs: &[Tally]) -> Tally {
    let mut all = Tally {
        tips: vec![0.0; runs.first().map_or(0, |t| t.tips.len())],
        ..Tally::default()
    };
    for t in runs {
        all.samples += t.samples;
        all.kappa += t.kappa;
        all.ln_l += t.ln_l;
        all.tips.iter_mut().zip(&t.tips).for_each(|(a, b)| *a += b);
        for (key, &(count, sum)) in &t.splits {
            let seen = all.splits.entry(*key).or_insert((0.0, 0.0));
            (seen.0, seen.1) = (seen.0 + count, seen.1 + sum);
        }
    }
    all
}

/// Average standard deviation of split frequencies between runs, over the splits that reach 10 %
/// in some run; `None` with a single run or no such split.
fn asdsf(runs: &[Tally]) -> Option<f64> {
    if runs.len() < 2 {
        return None;
    }
    let keys: HashSet<u128> = runs
        .iter()
        .flat_map(|t| {
            t.splits
                .iter()
                .filter(|(_, v)| v.0 / t.samples >= ASDSF_MIN_FREQUENCY)
                .map(|(k, _)| *k)
        })
        .collect();
    if keys.is_empty() {
        return None;
    }
    let deviation = |key: &u128| -> f64 {
        let f: Vec<f64> = runs
            .iter()
            .map(|t| t.splits.get(key).map_or(0.0, |v| v.0) / t.samples)
            .collect();
        let mean = f.iter().sum::<f64>() / f.len() as f64;
        let spread: f64 = f.iter().map(|x| (x - mean).powi(2)).sum();
        (spread / (f.len() - 1) as f64).sqrt()
    };
    Some(keys.iter().map(deviation).sum::<f64>() / keys.len() as f64)
}

/// The majority-rule consensus: leaf 0 and every clade above 50 % hang from one root.
fn consensus(names: &[String], all: &Tally) -> Node {
    let leaf = |i: usize| Node {
        name: names[i].clone(),
        length: all.tips[i] / all.samples,
        support: None,
        children: Vec::new(),
    };
    let mut clades: Vec<(u128, f64, f64)> = all
        .splits
        .iter()
        .map(|(&mask, &(count, sum))| (mask, count / all.samples, sum / count))
        .filter(|c| c.1 > CONSENSUS_ABOVE)
        .collect();
    // smallest clades first, so each one gathers the finished pieces inside it
    clades.sort_by_key(|c| (c.0.count_ones(), c.0));
    let mut open: Vec<(u128, Node)> = (1..names.len()).map(|i| (1u128 << i, leaf(i))).collect();
    for (mask, probability, length) in clades {
        let (inside, outside): (Vec<_>, Vec<_>) =
            open.into_iter().partition(|(m, _)| m & !mask == 0);
        let mut node = Node::inner(inside.into_iter().map(|(_, n)| n).collect());
        (node.length, node.support) = (length, Some(100.0 * probability));
        open = outside;
        open.push((mask, node));
    }
    Node::inner(
        std::iter::once(leaf(0))
            .chain(open.into_iter().map(|(_, n)| n))
            .collect(),
    )
}

fn run_chain(model: &Model, start: &Tree, s: &Settings, run: usize) -> Tally {
    let seed = s
        .seed
        .wrapping_add(run as u64 + 1)
        .wrapping_mul(0x9E37_79B9_7F4A_7C15);
    let mut chain = Chain::new(start, model, seed, SCATTER_SWAPS * run);
    let burn = (s.generations as f64 * s.burn_in) as usize;
    let every = s.sample_every.max(1);
    let mut tally = Tally {
        tips: vec![0.0; start.tips],
        ..Tally::default()
    };
    for g in 1..=s.generations {
        chain.generation(model);
        if g > burn && g % every == 0 {
            tally.record(&chain.tree, chain.kappa, chain.ln_l);
        }
    }
    tally
}

/// Sample trees from the posterior of an alignment of 3 to 127 sequences.
pub fn bayes_tree(msa: &Msa, settings: &Settings) -> Result<BayesResult, String> {
    let n = msa.rows.len();
    if !(3..=MAX_TAXA).contains(&n) {
        return Err(format!(
            "A Bayesian tree needs between 3 and {MAX_TAXA} sequences."
        ));
    }
    let data = sites(msa).ok_or("No column of the alignment holds a base.")?;
    let distances = [
        DistanceModel::Kimura2,
        DistanceModel::JukesCantor,
        DistanceModel::P,
    ]
    .iter()
    .find_map(|&m| msa_distances(msa, m).ok())
    .ok_or("The sequences are too different to start a tree from.")?;
    let start = neighbor_joining(&msa.names, &distances)
        .and_then(|t| Tree::from_node(&t, &msa.names))
        .ok_or("The starting tree could not be built.")?;
    let model = Model { sites: &data };
    let runs: Vec<usize> = (0..settings.runs.max(1)).collect();
    let tallies = par_map(&runs, |&r| run_chain(&model, &start, settings, r));
    let all = pooled(&tallies);
    if all.samples == 0.0 {
        return Err("No tree was sampled: use more generations or a shorter burn-in.".into());
    }
    Ok(BayesResult {
        tree: consensus(&msa.names, &all),
        asdsf: asdsf(&tallies),
        samples: all.samples as usize,
        mean_ln_l: all.ln_l / all.samples,
        mean_kappa: all.kappa / all.samples,
    })
}

#[cfg(test)]
mod tests {
    use super::super::tests::simulated;
    use super::*;
    use crate::seq::tree::splits;

    fn four_tips() -> Tree {
        let edges = [
            (0, 4, 0.1),
            (1, 4, 0.1),
            (4, 5, 0.1),
            (2, 5, 0.1),
            (3, 5, 0.1),
        ];
        let mut adj = vec![Vec::new(); 6];
        for (a, b, l) in edges {
            adj[a].push((b, l));
            adj[b].push((a, l));
        }
        Tree { adj, tips: 4 }
    }

    // A sampler is right when, given no data at all, it hands back its priors.
    #[test]
    fn with_no_data_the_sampler_returns_its_priors() {
        let data = Sites {
            tips: 4,
            patterns: vec![vec![4, 4, 4, 4]],
            counts: vec![1.0],
            total: 1.0,
            freqs: [0.25; 4],
        };
        let model = Model { sites: &data };
        let mut chain = Chain::new(&four_tips(), &model, 7, 0);
        let mut tally = Tally {
            tips: vec![0.0; 4],
            ..Tally::default()
        };
        for g in 1..=60_000 {
            chain.generation(&model);
            if g % 5 == 0 {
                tally.record(&chain.tree, chain.kappa, chain.ln_l);
            }
        }
        // the three topologies of four sequences, each a third of the time
        for mask in [0b1100u128, 0b1010, 0b0110] {
            let share = tally.splits.get(&mask).map_or(0.0, |v| v.0) / tally.samples;
            assert!((share - 1.0 / 3.0).abs() < 0.05, "{mask:b}: {share}");
        }
        // every branch length is exponential with mean 0.1
        let mean = tally.tips[1] / tally.samples;
        assert!((mean - BRANCH_PRIOR_MEAN).abs() < 0.015, "{mean}");
    }

    #[test]
    fn the_sampler_finds_the_true_split_and_the_runs_agree() {
        let msa = simulated();
        let settings = Settings {
            generations: 4_000,
            sample_every: 4,
            ..Settings::default()
        };
        let result = bayes_tree(&msa, &settings).unwrap();
        assert!(splits(&result.tree, &msa.names).contains(&vec![2, 3]));
        let inner = result.tree.children.iter().find(|c| !c.is_leaf()).unwrap();
        assert!(inner.support.unwrap() > 95.0, "{:?}", inner.support);
        assert!(
            inner.length > 0.1 && inner.length < 0.35,
            "{}",
            inner.length
        );
        assert!(
            result.asdsf.is_some_and(|a| a < CONVERGED_BELOW),
            "{:?}",
            result.asdsf
        );
        assert!(result.mean_kappa > 0.1 && result.mean_ln_l < 0.0);
        assert_eq!(result.samples, 2 * 750);
        // the same seed gives the same tree
        assert_eq!(bayes_tree(&msa, &settings).unwrap().tree, result.tree);
        let two = Msa {
            names: vec!["a".into(), "b".into()],
            rows: vec![b"AC".to_vec(); 2],
        };
        assert!(bayes_tree(&two, &settings).is_err());
    }
}
