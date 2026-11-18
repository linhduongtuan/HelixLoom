//! Distances between all the sequences of an alignment on PyTorch workers: the same three products as `accel::distance` (counted
//! as whole numbers in 32 bits, so the distances equal the processor's to the last digit), with the columns written out as 0 and 1
//! on each card by its worker, from the letters, which is all that is sent.

use super::wire::{Array, Message};
use super::{next_set, release, TorchPool};
use crate::accel::distance::{check_rows, check_weights, Counts};
use std::sync::Arc;

/// The alignment on every worker's card, ready for any number of weightings of its columns.
pub struct Engine {
    pool: Arc<TorchPool>,
    id: u32,
    n: usize,
    columns: usize,
}

impl std::fmt::Debug for Engine {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("torch::distance::Engine")
            .field("sequences", &self.n)
            .field("columns", &self.columns)
            .field("workers", &self.pool.len())
            .finish()
    }
}

impl Engine {
    /// Send the rows of an alignment (equal in length, `-` for a gap) to every worker.
    pub fn new(pool: &Arc<TorchPool>, rows: &[Vec<u8>]) -> Result<Engine, String> {
        let (n, columns) = check_rows(rows)?;
        let id = next_set();
        let flat: Vec<u8> = rows.iter().flatten().copied().collect();
        let load = Message::new("load_distances")
            .with_number("id", id as usize)
            .with_number("n", n)
            .with_number("columns", columns)
            .with("rows", Array::U8(flat));
        let loaded = pool.run(|_, worker| worker.call(&load).map(|_| ()));
        let engine = Engine {
            pool: Arc::clone(pool),
            id,
            n,
            columns,
        };
        // a worker that failed to take them: the ones that did are told to let go (by `drop`) and the reason is given
        loaded.into_iter().collect::<Result<Vec<()>, String>>()?;
        Ok(engine)
    }

    /// Which devices, and what share of the rows each does, in words.
    pub fn describe(&self) -> String {
        format!("PyTorch: {}", self.pool.distance_description())
    }

    /// The counts for every pair, each column taken `weights[column]` times (none: once). The rows of the answer are shared out
    /// among the workers by their speed at matrix products.
    pub fn counts(&self, weights: Option<&[u32]>) -> Result<Counts, String> {
        if let Some(w) = weights {
            check_weights(w, self.columns)?;
        }
        let n = self.n;
        let shares = self.pool.distance_shares(n);
        let parts = self.pool.run(|k, worker| -> Result<[Vec<f32>; 3], String> {
            let rows = shares[k].clone();
            if rows.is_empty() {
                return Ok([Vec::new(), Vec::new(), Vec::new()]);
            }
            let mut ask = Message::new("counts")
                .with_number("id", self.id as usize)
                .with_number("first", rows.start)
                .with_number("length", rows.len());
            if let Some(w) = weights {
                ask = ask.with("weights", Array::U32(w.to_vec()));
            }
            let reply = worker.call(&ask)?;
            let take = |name: &str| -> Result<Vec<f32>, String> {
                reply
                    .f32s(name)
                    .filter(|v| v.len() == rows.len() * n)
                    .map(<[f32]>::to_vec)
                    .ok_or_else(|| {
                        format!(
                            "{} gave no \"{name}\" counts of the right size.",
                            worker.name
                        )
                    })
            };
            Ok([take("valid")?, take("identical")?, take("same_class")?])
        });
        let (mut sites, mut identical, mut same_class) = (Vec::new(), Vec::new(), Vec::new());
        for part in parts {
            let [v, e, c] = part?;
            sites.extend(v);
            identical.extend(e);
            same_class.extend(c);
        }
        Ok(Counts::from_products(n, &sites, &identical, &same_class))
    }
}

impl Drop for Engine {
    fn drop(&mut self) {
        release(&self.pool, self.id);
    }
}
