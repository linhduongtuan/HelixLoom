//! What `helixloom gpu-check` and the tests need to try the likelihood on a pool of cards: an alignment evolved along a random
//! tree, a tree for it, and the comparison with the processor's likelihood, for correctness and for time.

use super::engine::Engine;
use super::{rate_matrix, sites, Kind, Matrix, Sites, Tree};
use crate::accel::Target;
use crate::seq::msa::Msa;
use crate::seq::par::par_map;
use crate::seq::tree::{msa_distances, neighbor_joining, DistanceModel, Rng};
use std::time::Instant;

/// `n` sequences related by a random tree: each new one is a copy of an earlier one with about `rate` of its columns changed, and
/// gaps and ambiguity codes here and there (a tip says "any base" there).
pub(super) fn simulate(n: usize, columns: usize, rate: f64, seed: u64) -> Msa {
    let mut rng = Rng::new(seed);
    let letters = b"ACGT";
    let mut all: Vec<Vec<u8>> = vec![(0..columns).map(|_| letters[rng.below(4)]).collect()];
    while all.len() < n {
        let parent = all[rng.below(all.len())].clone();
        let child = parent
            .iter()
            .map(|&c| {
                if (rng.below(1000) as f64) < rate * 1000.0 {
                    letters[rng.below(4)]
                } else {
                    c
                }
            })
            .collect();
        all.push(child);
    }
    let rows = all
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|c| match rng.below(50) {
                    0 => b'-',
                    1 => b'N',
                    _ => c,
                })
                .collect()
        })
        .collect();
    Msa {
        names: (0..n).map(|k| format!("s{k}")).collect(),
        rows,
    }
}

/// The columns of an alignment and a tree for it (the neighbour-joining tree, with every branch given a length that is usable).
pub(super) fn setup(msa: &Msa) -> Result<(Sites, Tree), String> {
    let sites = sites(msa).ok_or("No column holds a base.")?;
    let d = msa_distances(msa, DistanceModel::P)?;
    let node = neighbor_joining(&msa.names, &d).ok_or("No tree could be made.")?;
    let mut tree = Tree::from_node(&node, &msa.names).ok_or("No tree could be made.")?;
    for edges in &mut tree.adj {
        for e in edges.iter_mut() {
            e.1 = e.1.clamp(1e-3, 2.0);
        }
    }
    // both ends of a branch carry its length: make them agree after the clamping
    for a in 0..tree.adj.len() {
        for k in 0..tree.adj[a].len() {
            let (b, l) = tree.adj[a][k];
            tree.set_length(a, b, l);
        }
    }
    Ok((sites, tree))
}

/// The rate matrix and base frequencies a model fits with (kappa 2.7 where it has one).
pub(super) fn parameters(kind: Kind, sites: &Sites) -> (Matrix, [f64; 4]) {
    let freqs = if kind.free_freqs() {
        sites.freqs
    } else {
        [0.25; 4]
    };
    let kappa = if kind.free_kappa() { 2.7 } else { 1.0 };
    (rate_matrix(kappa, &freqs), freqs)
}

/// How well the likelihood on the cards agrees with the processor's.
#[derive(Clone, Copy, Debug)]
pub struct Agreement {
    /// The largest difference over the models, as a share of the log-likelihood (or of 1 if that is smaller).
    pub worst_relative: f64,
    /// The same difference in log-likelihood units: what a search that compares two trees by their likelihoods sees.
    pub worst_absolute: f64,
    pub models: usize,
    pub patterns: usize,
}

/// The reason the cards did not give a likelihood, in words.
fn missing(engine: &Engine) -> String {
    format!(
        "The cards did not give a likelihood ({}).",
        engine
            .failure()
            .unwrap_or_else(|| "no reason was kept".into())
    )
}

/// Evolve an alignment of `tips` sequences and `columns` columns, and compare the likelihood of its tree on `target` with the
/// processor's, under each of the four models: one after the other from this thread (the cards share the patterns), and again from
/// the threads of `par_map`, each of which is given one card for the whole tree, as the model fits and the sampler's chains are.
pub fn agree(target: &Target, tips: usize, columns: usize, seed: u64) -> Result<Agreement, String> {
    let msa = simulate(tips, columns, 0.12, seed);
    let (sites, tree) = setup(&msa)?;
    let engine = Engine::on(target, &sites)?;
    let (mut relative, mut absolute) = (0.0f64, 0.0f64);
    let mut compare = |got: f64, want: f64| {
        relative = relative.max((got - want).abs() / want.abs().max(1.0));
        absolute = absolute.max((got - want).abs());
    };
    for kind in Kind::ALL {
        let (q, freqs) = parameters(kind, &sites);
        let want = tree.ln_likelihood(&sites, &q, &freqs);
        let got = engine
            .ln_likelihood(&tree, &q, &freqs)
            .ok_or_else(|| missing(&engine))?;
        compare(got, want);
    }
    let from_workers = par_map(&Kind::ALL, |&kind| {
        let (q, freqs) = parameters(kind, &sites);
        let want = tree.ln_likelihood(&sites, &q, &freqs);
        engine
            .ln_likelihood(&tree, &q, &freqs)
            .map(|got| (got, want))
    });
    for pair in from_workers {
        let (got, want) = pair.ok_or_else(|| missing(&engine))?;
        compare(got, want);
    }
    Ok(Agreement {
        worst_relative: relative,
        worst_absolute: absolute,
        models: Kind::ALL.len(),
        patterns: sites.patterns.len(),
    })
}

/// How long one evaluation of the likelihood of a tree takes on the processor and on `target` (the average of a few, after one to
/// warm the cards up), for an alignment of this size.
#[derive(Clone, Copy, Debug)]
pub struct Timing {
    pub patterns: usize,
    pub tips: usize,
    pub cpu_seconds: f64,
    pub cards_seconds: f64,
}

pub fn time(target: &Target, tips: usize, columns: usize) -> Result<Timing, String> {
    let msa = simulate(tips, columns, 0.12, 17);
    let (sites, tree) = setup(&msa)?;
    let engine = Engine::on(target, &sites)?;
    let (q, freqs) = parameters(Kind::Hky85, &sites);
    let start = Instant::now();
    let want = tree.ln_likelihood(&sites, &q, &freqs);
    let cpu_seconds = start.elapsed().as_secs_f64();
    let run = || {
        engine
            .ln_likelihood(&tree, &q, &freqs)
            .ok_or_else(|| missing(&engine))
    };
    run()?;
    const REPEATS: usize = 3;
    let start = Instant::now();
    let mut got = 0.0;
    for _ in 0..REPEATS {
        got = run()?;
    }
    let cards_seconds = start.elapsed().as_secs_f64() / REPEATS as f64;
    if (got - want).abs() > 1e-4 * want.abs().max(1.0) {
        return Err(format!("The cards gave {got}, the processor {want}."));
    }
    Ok(Timing {
        patterns: sites.patterns.len(),
        tips,
        cpu_seconds,
        cards_seconds,
    })
}
