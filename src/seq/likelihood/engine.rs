//! The engine that puts the likelihood on the GPUs, whichever library the person chose for the work (`accel::backend`): candle, which is
//! built into the program, or PyTorch, in Python workers.

use super::{gpu, torch, Matrix, Sites, Tree};
use crate::accel::torch::TorchPool;
use crate::accel::{backend, Backend, Pool, Target};

pub enum Engine {
    Candle(gpu::Engine),
    Torch(torch::Engine),
}

impl Engine {
    /// Put the columns of `sites` on the devices of the backend in use: `Ok(None)` where it has none, `Err` with the reason where it
    /// has some and cannot take them.
    pub(super) fn open(sites: &Sites) -> Result<Option<Engine>, String> {
        match backend() {
            Backend::Candle => Pool::global()
                .map(|pool| gpu::Engine::new(pool, sites).map(Engine::Candle))
                .transpose(),
            Backend::Torch => TorchPool::global()
                .map(|pool| torch::Engine::new(pool, sites).map(Engine::Torch))
                .transpose(),
        }
    }

    /// The same on the devices of `target` (all of them), whichever backend it is of.
    pub(super) fn on(target: &Target, sites: &Sites) -> Result<Engine, String> {
        match target {
            Target::Candle(pool) => gpu::Engine::new(pool, sites).map(Engine::Candle),
            Target::Torch(pool) => torch::Engine::new(pool, sites).map(Engine::Torch),
        }
    }

    /// Which devices, and what share each does, in words.
    pub fn describe(&self) -> String {
        match self {
            Engine::Candle(e) => e.describe(),
            Engine::Torch(e) => e.describe(),
        }
    }

    pub fn broken(&self) -> bool {
        match self {
            Engine::Candle(e) => e.broken(),
            Engine::Torch(e) => e.broken(),
        }
    }

    /// What went wrong on a card whose share another card then did, if that happened (PyTorch only): the work went on on the cards.
    pub fn retried(&self) -> Option<String> {
        match self {
            Engine::Candle(_) => None,
            Engine::Torch(e) => e.retried(),
        }
    }

    /// What went wrong with the device that failed first, if one did.
    pub fn failure(&self) -> Option<String> {
        match self {
            Engine::Candle(e) => e.failure(),
            Engine::Torch(e) => e.failure(),
        }
    }

    /// The likelihood as a function of the length of the branch `a`-`b` alone, kept on the devices so that each length costs one product for
    /// each pattern and not the pruning of the whole tree: where PyTorch does the work (candle does not keep partial likelihoods yet), and
    /// `None` otherwise, or where a device has failed, when the caller tries the whole tree for every length.
    pub(super) fn edge_scan(
        &self,
        tree: &Tree,
        a: usize,
        b: usize,
        q: &Matrix,
        freqs: &[f64; 4],
    ) -> Option<torch::EdgeScan<'_>> {
        match self {
            Engine::Candle(_) => None,
            Engine::Torch(e) => e.edge_scan(tree, a, b, q, freqs),
        }
    }

    /// The log-likelihood of `tree`, or `None` if a device failed (and from then on).
    pub(super) fn ln_likelihood(&self, tree: &Tree, q: &Matrix, freqs: &[f64; 4]) -> Option<f64> {
        match self {
            Engine::Candle(e) => e.ln_likelihood(tree, q, freqs),
            Engine::Torch(e) => e.ln_likelihood(tree, q, freqs),
        }
    }
}
