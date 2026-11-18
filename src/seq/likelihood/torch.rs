//! The likelihood of a tree on PyTorch workers: the same pruning as the candle engine (`gpu`), with the same host-side plan
//! (`plan`), done by `src/accel/torch/worker.py`. Every worker holds the columns of the alignment (a number for each tip at each
//! pattern) on its card. A call from a thread of `par_map` (a model being fitted, a chain of the sampler) goes to the worker that
//! thread was given; any other call is cut in ranges of patterns, one for each worker by its speed, which are done at once. Each
//! worker sums the 32-bit logarithms of its patterns, weighted by how often each occurs, in 64 bits and sends the one number.

use super::plan::{Beyond, Plan};
use super::{expm, Matrix, Sites, Tree};
use crate::accel::torch::wire::{Array, Message};
use crate::accel::torch::{next_set, release, TorchPool};
use crate::seq::par::worker;
use std::ops::Range;
use std::sync::{Arc, Mutex, PoisonError};

/// What the tips, written out, may take on one card, in bytes (a 32-bit number each).
const RESIDENT_LIMIT: usize = 4 << 30;

pub struct Engine {
    pool: Arc<TorchPool>,
    id: u32,
    patterns: usize,
    tips: usize,
    /// Why a card failed, the first time it did: the processor does the rest from then on.
    failure: Mutex<Option<String>>,
    /// Why a card failed on a share that another card did, the first time: what went wrong, for the report, though the work went on.
    retried: Mutex<Option<String>>,
}

impl Engine {
    /// Put the columns of `sites` on every worker's card.
    pub(super) fn new(pool: &Arc<TorchPool>, sites: &Sites) -> Result<Engine, String> {
        if sites.tips < 3 || sites.patterns.is_empty() {
            return Err("There are fewer than three tips or no column.".into());
        }
        let patterns = sites.patterns.len();
        let bytes = sites.tips * patterns * 4;
        if bytes > RESIDENT_LIMIT {
            return Err(format!(
                "The columns written out need {} MB on each card, more than the {} MB that are used.",
                bytes >> 20,
                RESIDENT_LIMIT >> 20
            ));
        }
        // a byte for each tip at each pattern, tip after tip
        let states: Vec<u8> = (0..sites.tips)
            .flat_map(|t| sites.patterns.iter().map(move |p| p[t]))
            .collect();
        let id = next_set();
        let load = Message::new("load_likelihood")
            .with_number("id", id as usize)
            .with_number("tips", sites.tips)
            .with_number("patterns", patterns)
            .with("states", Array::U8(states))
            .with("counts", Array::F64(sites.counts.clone()));
        let loaded = pool.run(|_, w| w.call(&load).map(|_| ()));
        // built first, so that the workers that did take the columns are told to let go if one did not
        let engine = Engine {
            pool: Arc::clone(pool),
            id,
            patterns,
            tips: sites.tips,
            failure: Mutex::new(None),
            retried: Mutex::new(None),
        };
        loaded.into_iter().collect::<Result<Vec<()>, String>>()?;
        Ok(engine)
    }

    /// Which devices, and what share each does, in words.
    pub fn describe(&self) -> String {
        format!(
            "PyTorch: {} ({} tips, {} patterns on each)",
            self.pool.likelihood_description(),
            self.tips,
            self.patterns
        )
    }

    /// Whether a card has failed, so that the processor now does the work.
    pub fn broken(&self) -> bool {
        self.failure().is_some()
    }

    /// What went wrong with the card that failed first, if one did.
    pub fn failure(&self) -> Option<String> {
        self.failure
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .clone()
    }

    /// The log-likelihood of `tree`, or `None` if a card has failed (and from then on): the caller goes on with the processor's.
    pub(super) fn ln_likelihood(&self, tree: &Tree, q: &Matrix, freqs: &[f64; 4]) -> Option<f64> {
        if self.broken() {
            return None;
        }
        let why = match self.evaluate(tree, q, freqs) {
            Ok(v) if v.is_finite() => return Some(v),
            Ok(v) => format!("a likelihood of {v}"),
            Err(why) => why,
        };
        crate::accel::note_failure(
            "the GPUs failed on a likelihood",
            &format!("{why}; the processor finishes it"),
        );
        let mut failure = self.failure.lock().unwrap_or_else(PoisonError::into_inner);
        failure.get_or_insert(why);
        None
    }

    fn evaluate(&self, tree: &Tree, q: &Matrix, freqs: &[f64; 4]) -> Result<f64, String> {
        let plan = Plan::new(tree, q, freqs)?;
        let message = |range: Range<usize>| {
            Message::new("likelihood")
                .with_number("id", self.id as usize)
                .with_number("first", range.start)
                .with_number("length", range.len())
                .with_number("first_node", plan.first)
                .with("steps", Array::U32(plan.steps_table()))
                .with("inner_p", Array::F32(plan.inner_p.clone()))
                .with("tip_p", Array::F32(plan.tip_p.clone()))
                .with("root_p", Array::F32(plan.root_p.to_vec()))
                .with("root_table", Array::F32(plan.root_table.to_vec()))
        };
        let on = |k: usize, range: Range<usize>| -> Result<f64, String> {
            let worker = &self.pool.workers()[k];
            let reply = worker.call(&message(range))?;
            reply
                .f64s("sum")
                .and_then(|s| s.first().copied())
                .ok_or_else(|| format!("{} gave no sum.", worker.name))
        };
        // a thread of `par_map` has a worker of its own; any other call is shared among the workers
        let partial: Vec<Result<f64, String>> = match worker() {
            Some(w) => vec![self.retrying(w % self.pool.len(), 0..self.patterns, &on)],
            None if self.pool.len() == 1 => vec![on(0, 0..self.patterns)],
            None => {
                let shares = self.pool.likelihood_shares(self.patterns);
                self.pool.run(|k, _| {
                    if shares[k].is_empty() {
                        Ok(0.0)
                    } else {
                        self.retrying(k, shares[k].clone(), &on)
                    }
                })
            }
        };
        partial.into_iter().sum()
    }

    /// `on(first, range)`, and where that worker fails, the same patterns on each of the other workers in turn: every worker holds all the
    /// columns, so any can do another's share, and a card that failed costs a share's time and not the cards' work. `Err` (the first
    /// worker's reason) only where none of them can. What went wrong is kept for the report of the result.
    fn retrying(
        &self,
        first: usize,
        range: Range<usize>,
        on: &impl Fn(usize, Range<usize>) -> Result<f64, String>,
    ) -> Result<f64, String> {
        let why = match on(first, range.clone()) {
            Ok(v) => return Ok(v),
            Err(why) => why,
        };
        for k in (0..self.pool.len()).filter(|&k| k != first) {
            if let Ok(v) = on(k, range.clone()) {
                crate::accel::note_failure(
                    "a GPU failed on a share of a likelihood",
                    &format!("{why}; another card did it"),
                );
                self.retried
                    .lock()
                    .unwrap_or_else(PoisonError::into_inner)
                    .get_or_insert(why);
                return Ok(v);
            }
        }
        Err(why)
    }

    /// Why a card failed on a share that another card then did, the first time it did: the likelihood went on on the cards.
    pub fn retried(&self) -> Option<String> {
        self.retried
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .clone()
    }
}

/// The likelihood of the alignment as a function of the length of one branch, kept on the workers (`edge_begin` in `worker.py`): what
/// lies beyond both ends of the branch is pruned once, and each length asked for then costs one product for each pattern, not the pruning
/// of the whole tree. Where a thread of `par_map` made it, one worker holds it all; otherwise every worker holds its range of patterns.
pub(super) struct EdgeScan<'a> {
    engine: &'a Engine,
    id: u32,
    /// The workers that hold a part, and the patterns of each.
    parts: Vec<(usize, Range<usize>)>,
}

impl Engine {
    /// Keep what lies beyond both ends of the branch `a`-`b` of `tree` on the workers. `None` where a card has failed or the tree could not
    /// be planned: the caller goes on as before, with the whole tree for every length.
    pub(super) fn edge_scan(
        &self,
        tree: &Tree,
        a: usize,
        b: usize,
        q: &Matrix,
        freqs: &[f64; 4],
    ) -> Option<EdgeScan<'_>> {
        if self.broken() {
            return None;
        }
        let (Ok(end_a), Ok(end_b)) = (Beyond::new(tree, a, b, q), Beyond::new(tree, b, a, q))
        else {
            return None;
        };
        let id = next_set();
        let message = |range: Range<usize>| {
            let mut m = Message::new("edge_begin")
                .with_number("id", self.id as usize)
                .with_number("scan", id as usize)
                .with_number("first", range.start)
                .with_number("length", range.len())
                .with(
                    "freqs",
                    Array::F32(freqs.iter().map(|&f| f as f32).collect()),
                );
            for (name, end) in [("a", &end_a), ("b", &end_b)] {
                m = m
                    .with_number(&format!("{name}_tip"), usize::from(end.is_tip))
                    .with_number(&format!("{name}_node"), end.node)
                    .with(&format!("{name}_steps"), Array::U32(end.steps_table()))
                    .with(&format!("{name}_inner_p"), Array::F32(end.inner_p.clone()))
                    .with(&format!("{name}_tip_p"), Array::F32(end.tip_p.clone()));
            }
            m
        };
        let parts: Vec<(usize, Range<usize>)> = match worker() {
            Some(w) => vec![(w % self.pool.len(), 0..self.patterns)],
            None if self.pool.len() == 1 => vec![(0, 0..self.patterns)],
            None => self
                .pool
                .likelihood_shares(self.patterns)
                .into_iter()
                .enumerate()
                .filter(|(_, range)| !range.is_empty())
                .collect(),
        };
        let scan = EdgeScan {
            engine: self,
            id,
            parts,
        };
        let begun = scan.on_parts(|k, range| {
            self.pool.workers()[k]
                .call(&message(range.clone()))
                .map(|_| ())
        });
        // `scan` lets go of what it holds on the workers if it is dropped here. A worker that cannot keep its part is not a failure of the
        // cards yet: the caller goes on with the whole tree for each length, whose shares are retried on the other workers, and which says
        // so if the cards cannot do it.
        begun
            .into_iter()
            .collect::<Result<Vec<()>, String>>()
            .ok()
            .map(|_| scan)
    }
}

impl EdgeScan<'_> {
    /// `f` for each worker that holds a part (with its patterns), all at once where there are several.
    fn on_parts<T: Send>(&self, f: impl Fn(usize, &Range<usize>) -> T + Sync) -> Vec<T> {
        if let [(k, range)] = self.parts.as_slice() {
            return vec![f(*k, range)];
        }
        let results = self.engine.pool.run(|k, _| {
            self.parts
                .iter()
                .find(|(part, _)| *part == k)
                .map(|(_, range)| f(k, range))
        });
        results.into_iter().flatten().collect()
    }

    /// The log-likelihood of the alignment with the branch `length` long, or `None` if a card has failed (and from then on).
    pub(super) fn ln_likelihood(&self, q: &Matrix, length: f64) -> Option<f64> {
        if self.engine.broken() {
            return None;
        }
        let p: Vec<f32> = expm(q, length)
            .iter()
            .flatten()
            .map(|&x| x as f32)
            .collect();
        let ask = Message::new("edge_eval")
            .with_number("scan", self.id as usize)
            .with("p", Array::F32(p));
        let sums = self.on_parts(|k, _| {
            let worker = &self.engine.pool.workers()[k];
            worker
                .call(&ask)?
                .f64s("sum")
                .and_then(|s| s.first().copied())
                .ok_or_else(|| format!("{} gave no sum.", worker.name))
        });
        // a worker that fails, or a number that is not one, is left to the whole tree's likelihood to find out about and report
        sums.into_iter()
            .sum::<Result<f64, String>>()
            .ok()
            .filter(|v| v.is_finite())
    }
}

impl Drop for EdgeScan<'_> {
    fn drop(&mut self) {
        let end = Message::new("edge_end").with_number("scan", self.id as usize);
        // a worker that has failed has nothing to forget
        self.on_parts(|k, _| {
            let _ = self.engine.pool.workers()[k].call(&end);
        });
    }
}

impl Drop for Engine {
    fn drop(&mut self) {
        release(&self.pool, self.id);
    }
}

#[cfg(test)]
mod tests;
