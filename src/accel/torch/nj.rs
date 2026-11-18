//! Neighbour-joining on a PyTorch worker: the joins of `seq::tree::cpu_joins`, made over the whole distance matrix on a card in 64 bits.
//! The processor needs about n cubed steps; the card does each join in a few passes over the matrix.

use super::wire::{Array, Message};
use super::TorchPool;
use crate::seq::tree::Join;

/// The most sequences one card is asked to join: the host holds several copies of the n by n matrix of 64-bit numbers on the way (about
/// 19 GB at 20,000), and a request that takes more than the half hour a worker is given is killed.
const MOST_SEQUENCES: usize = 15_000;
/// The largest distance taken: the scores are (n - 2) times a distance, which must not overflow.
const LARGEST_DISTANCE: f64 = 1e100;

/// Whether `n` sequences are within what is asked of a card.
fn fits(n: usize) -> bool {
    n <= MOST_SEQUENCES
}

/// The joins that neighbour-joining makes of `dist`, and the distances among the three clusters left (in the order of the list), made on
/// the pool's fastest worker. `Err` says what failed (a device with no 64-bit floats, such as a Mac's GPU, is one).
pub fn joins(pool: &TorchPool, dist: &[Vec<f64>]) -> Result<(Vec<Join>, Vec<Vec<f64>>), String> {
    let n = dist.len();
    if n <= 3 {
        return Err("There are no joins to make for fewer than four sequences.".into());
    }
    // an infinite or missing distance has no neighbour-joining that the two sides would agree on
    if dist
        .iter()
        .flatten()
        .any(|d| !d.is_finite() || d.abs() > LARGEST_DISTANCE)
        || dist.iter().any(|row| row.len() != n)
    {
        return Err("The distances are not an n by n matrix of finite numbers.".into());
    }
    if !fits(n) {
        return Err(format!(
            "{n} sequences are more than the {MOST_SEQUENCES} that a card is asked to join."
        ));
    }
    let worker = &pool.workers()[pool.fastest()];
    let reply = worker.call(
        &Message::new("neighbor_joining")
            .with_number("n", n)
            .with("dist", Array::F64(dist.iter().flatten().copied().collect())),
    )?;
    let (first, second) = (reply.u32s("first"), reply.u32s("second"));
    let (li, lj, last) = (reply.f64s("li"), reply.f64s("lj"), reply.f64s("last"));
    let (Some(first), Some(second), Some(li), Some(lj), Some(last)) = (first, second, li, lj, last)
    else {
        return Err(format!("{} gave no joins.", worker.name));
    };
    let steps = n - 3;
    if [first.len(), second.len(), li.len(), lj.len()] != [steps; 4] || last.len() != 9 {
        return Err(format!("{} gave joins of the wrong size.", worker.name));
    }
    let joins = (0..steps)
        .map(|k| Join {
            i: first[k] as usize,
            j: second[k] as usize,
            li: li[k],
            lj: lj[k],
        })
        .collect();
    Ok((joins, last.chunks(3).map(<[f64]>::to_vec).collect()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_card_is_asked_to_join_only_so_many_sequences() {
        assert!(fits(4) && fits(MOST_SEQUENCES));
        assert!(!fits(MOST_SEQUENCES + 1));
    }
}
