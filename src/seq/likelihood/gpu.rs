//! The likelihood of a tree over the columns of an alignment on the GPUs, which is what maximum-likelihood search, model fitting
//! and Bayesian sampling spend their time on (`Tree::ln_likelihood`). The same pruning as the processor's, written as tensor
//! operations over all the distinct columns (patterns) at once:
//!
//! * a tip is one number for each pattern, the place of its base (A, C, G, T; a gap or an ambiguity code is a fifth place). What
//!   a tip says to its parent for a base is a column of the transition matrix `P`, so it is looked up in a table of five rows
//!   (the fifth all ones: "any base"): no matrix product;
//! * an inner node multiplies what its two children say to it; what a child that is itself an inner node says is `P` times its
//!   vector, written as an elementwise product and a sum rather than a matrix product, which a card may do with rounded inputs
//!   (CUDA's TF32) and which would cost the likelihood its digits;
//! * every inner node's vector is divided by its largest element, and the logarithm of that carried along, so that a deep tree
//!   does not underflow in 32 bits.
//!
//! The columns live on every card (a copy each: the tips are one 32-bit number a pattern each). A call from one of the threads
//! of `par_map` (a model being fitted, a chain of the sampler) is computed by one card, the one that thread was given, so that
//! the cards work on different models or chains side by side; any other call (a search, one big tree) is cut in ranges of
//! patterns, one for each card in proportion to its weight, which are computed at once and added. The sum over the patterns is
//! made in 64 bits on the host from the 32-bit logarithms. All of this agrees with the processor to about seven digits of each
//! pattern's likelihood, not to the last: the tests say how closely.

use super::plan::Plan;
use super::{Matrix, Sites, Tree};
use crate::accel::pool::Pool;
use crate::seq::par::worker;
use candle_core::{Device, Tensor};
use std::ops::Range;
use std::sync::{Mutex, PoisonError};

/// What the tips, written out, may take on one card, in bytes.
const RESIDENT_LIMIT: usize = 4 << 30;

fn bad(e: candle_core::Error) -> String {
    format!("GPU: {e}")
}

struct Card {
    device: Device,
    /// For each tip, the place of its base at every pattern (0 to 3, or 4 for a gap or an ambiguity code).
    tips: Vec<Tensor>,
}

pub struct Engine {
    pool: Pool,
    cards: Vec<Card>,
    /// How many times each pattern occurs, on the host.
    counts: Vec<f64>,
    patterns: usize,
    tips: usize,
    /// Why a card failed, the first time it did: the processor does the rest from then on.
    failure: Mutex<Option<String>>,
}

impl Engine {
    /// Put the columns of `sites` on every card of `pool`.
    pub(super) fn new(pool: &Pool, sites: &Sites) -> Result<Engine, String> {
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
        let states: Vec<Vec<u32>> = (0..sites.tips)
            .map(|t| sites.patterns.iter().map(|p| u32::from(p[t])).collect())
            .collect();
        let cards = pool
            .members()
            .iter()
            .map(|m| {
                let _card = m.exclusive();
                let tips = states
                    .iter()
                    .map(|s| Tensor::from_vec(s.clone(), patterns, &m.device).map_err(bad))
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(Card {
                    device: m.device.clone(),
                    tips,
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        Ok(Engine {
            pool: pool.clone(),
            cards,
            counts: sites.counts.clone(),
            patterns,
            tips: sites.tips,
            failure: Mutex::new(None),
        })
    }

    /// Which cards, and what share each does, in words.
    pub fn describe(&self) -> String {
        let shares: Vec<String> = self
            .pool
            .members()
            .iter()
            .map(|m| format!("{} {:.0} %", m.name, 100.0 * m.weight))
            .collect();
        format!(
            "{} ({} tips, {} patterns on each)",
            shares.join(", "),
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

    /// The log-likelihood of `tree` (the sum over the patterns of their counts times the log of their likelihood), or `None` if
    /// a card has failed (and from then on): the caller goes on with the processor's.
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
        // a thread of `par_map` has a card of its own; any other call is shared among the cards
        let partial: Vec<Result<f64, String>> = match worker() {
            Some(w) => vec![self.on(w % self.cards.len(), 0..self.patterns, &plan)],
            None if self.cards.len() == 1 => vec![self.on(0, 0..self.patterns, &plan)],
            None => {
                let shares = self.pool.shares(self.patterns);
                self.pool.run(|k, _| {
                    if shares[k].is_empty() {
                        Ok(0.0)
                    } else {
                        self.on(k, shares[k].clone(), &plan)
                    }
                })
            }
        };
        partial.into_iter().sum()
    }

    /// The sum over `range` of the patterns' counts times the logarithms of their likelihoods, worked out on card `k`.
    fn on(&self, k: usize, range: Range<usize>, plan: &Plan) -> Result<f64, String> {
        // one thread at a time works on a card (see `Member::exclusive`)
        let _card = self.pool.members()[k].exclusive();
        let log = self.run(&self.cards[k], range.clone(), plan).map_err(bad)?;
        Ok(log
            .iter()
            .zip(&self.counts[range])
            .map(|(&l, &c)| c * f64::from(l))
            .sum())
    }

    /// The logarithm of the likelihood of each pattern of `range`.
    fn run(&self, card: &Card, range: Range<usize>, plan: &Plan) -> candle_core::Result<Vec<f32>> {
        let (device, len) = (&card.device, range.len());
        let inner_p =
            Tensor::from_vec(plan.inner_p.clone(), (plan.inner_p.len() / 16, 16), device)?;
        let tip_p = Tensor::from_vec(plan.tip_p.clone(), (plan.tip_p.len() / 20, 20), device)?;
        let missing = || candle_core::Error::Msg("a node was not done before its parent".into());
        // what each inner node holds: its vector for every pattern, and the logarithm of what it was divided by
        let mut done: Vec<Option<(Tensor, Tensor)>> = (0..plan.nodes).map(|_| None).collect();
        for step in &plan.steps {
            let (mut product, mut scale): (Option<Tensor>, Option<Tensor>) = (None, None);
            for child in &step.children {
                let (said, scaled) = match child.tip {
                    Some(tip) => {
                        let table = tip_p.narrow(0, child.slot, 1)?.reshape((5, 4))?;
                        let at = card.tips[tip].narrow(0, range.start, len)?;
                        (table.index_select(&at, 0)?, None)
                    }
                    None => {
                        let (vector, scale) = done[child.node].take().ok_or_else(missing)?;
                        let p = inner_p.narrow(0, child.slot, 1)?.reshape((1, 4, 4))?;
                        (vector.unsqueeze(1)?.broadcast_mul(&p)?.sum(2)?, Some(scale))
                    }
                };
                product = Some(match product {
                    None => said,
                    Some(other) => (other * said)?,
                });
                scale = match (scale, scaled) {
                    (Some(a), Some(b)) => Some((a + b)?),
                    (a, b) => a.or(b),
                };
            }
            let vector = product.ok_or_else(missing)?;
            let top = vector.max_keepdim(1)?;
            let vector = vector.broadcast_div(&top)?;
            let logs = top.log()?;
            let scale = match scale {
                Some(s) => (s + logs)?,
                None => logs,
            };
            done[step.node] = Some((vector, scale));
        }
        let (vector, scale) = done[plan.first].take().ok_or_else(missing)?;
        let root = Tensor::from_vec(plan.root_p.to_vec(), (1, 4, 4), device)?;
        let said = vector.unsqueeze(1)?.broadcast_mul(&root)?.sum(2)?;
        let table = Tensor::from_vec(plan.root_table.to_vec(), (5, 4), device)?;
        let first_tip = card.tips[0].narrow(0, range.start, len)?;
        let weights = table.index_select(&first_tip, 0)?;
        let likelihood = (weights * said)?.sum(1)?;
        (likelihood.log()? + scale.squeeze(1)?)?.to_vec1::<f32>()
    }
}

#[cfg(test)]
mod tests;
