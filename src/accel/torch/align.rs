//! Pairwise alignments on PyTorch workers: the tables of `seq::align` filled along their anti-diagonals for many pairs at once
//! (`align_batch` in `worker.py`), in 32 bits with the processor's additions, comparisons and ties, so that every alignment comes
//! out the same to the last bit. **Opt-in, and slower than the processor wherever it has been measured** (a Mac's GPU, a worker on the
//! CPU device: the dynamic programming is a chain of small operations for every diagonal and the processor does it eight pairs to a
//! vector instruction), so it is used only where the person asks (`HELIXLOOM_ACCEL=force`, or a size in
//! `HELIXLOOM_ACCEL_ALIGN_CELLS`), to be tried on cards that nobody here has: `helixloom gpu-check --bench` says where it wins, if it does.

use super::wire::{Array, Message};
use super::TorchPool;
use crate::accel::pool::shares;
use crate::seq::align::{alignment_from_device, Alignment, Scoring};

/// The most bytes of trace that one request asks a worker for and gets back (a byte for each cell of each pair): the pairs of a worker
/// are sent in as many requests as that takes.
const MOST_TRACE_BYTES: usize = 256 << 20;

/// The alignments of `a` with each of `bs` (all of them of at least one base, `a` too, and no pair too large for `align`), made by the
/// pool's workers, each doing a run of the pairs (by length, so that the pairs of a request are alike) in proportion to its speed, and
/// returned in the order of `bs`. `Err` says what failed: the caller leaves the pairs to the processor.
pub fn alignments(
    pool: &TorchPool,
    a: &[u8],
    bs: &[&[u8]],
    sc: &Scoring,
    free_ends: bool,
) -> Result<Vec<Alignment>, String> {
    if a.is_empty() || bs.iter().any(|b| b.is_empty()) {
        return Err("A sequence with no base has no table to fill on a device.".into());
    }
    let mut order: Vec<usize> = (0..bs.len()).collect();
    order.sort_by_key(|&i| bs[i].len());
    let cuts = shares(&pool.likelihood_weights, order.len());
    let done = pool.run(|k, worker| -> Result<Vec<(usize, Alignment)>, String> {
        let mine = &order[cuts[k].clone()];
        let mut out = Vec::with_capacity(mine.len());
        let mut rest = mine;
        while !rest.is_empty() {
            // as many of the pairs as fit, the longest last (the pairs are sorted by length)
            let mut count = 1;
            while count < rest.len()
                && (count + 1) * (a.len() + 1) * (bs[rest[count]].len() + 1) <= MOST_TRACE_BYTES
            {
                count += 1;
            }
            let (now, later) = rest.split_at(count);
            rest = later;
            let sent: Vec<&[u8]> = now.iter().map(|&i| bs[i]).collect();
            let width = sent.iter().map(|b| b.len()).max().unwrap_or(1) + 1;
            let reply = worker.call(
                &Message::new("align_batch")
                    .with("centre", Array::U8(a.to_vec()))
                    .with(
                        "lens",
                        Array::U32(sent.iter().map(|b| b.len() as u32).collect()),
                    )
                    .with("bases", Array::U8(sent.concat()))
                    .with_number("free_ends", usize::from(free_ends))
                    .with(
                        "scoring",
                        Array::F32(vec![
                            crate::seq::align::MATCH,
                            sc.mismatch,
                            sc.gap_open,
                            sc.gap_extend,
                        ]),
                    ),
            )?;
            let (Some(trace), Some(last_col), Some(last_row)) = (
                reply.u8s("trace"),
                reply.f32s("last_col"),
                reply.f32s("last_row"),
            ) else {
                return Err(format!("{} gave no tables.", worker.name));
            };
            let (table, column, row) = ((a.len() + 1) * width, (a.len() + 1) * 3, width * 3);
            if trace.len() != count * table
                || last_col.len() != count * column
                || last_row.len() != count * row
            {
                return Err(format!("{} gave tables of the wrong size.", worker.name));
            }
            for (p, (&i, b)) in now.iter().zip(&sent).enumerate() {
                out.push((
                    i,
                    alignment_from_device(
                        a,
                        b,
                        width,
                        &trace[p * table..(p + 1) * table],
                        &last_col[p * column..(p + 1) * column],
                        &last_row[p * row..(p + 1) * row],
                        free_ends,
                    ),
                ));
            }
        }
        Ok(out)
    });
    let mut placed: Vec<Option<Alignment>> = (0..bs.len()).map(|_| None).collect();
    for part in done {
        for (i, alignment) in part? {
            placed[i] = Some(alignment);
        }
    }
    placed
        .into_iter()
        .collect::<Option<Vec<Alignment>>>()
        .ok_or_else(|| "A device left a pair out.".to_string())
}
