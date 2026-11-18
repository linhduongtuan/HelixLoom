//! The joins of neighbour-joining on the processor.
//!
//! Every join looks for the pair of clusters with the smallest `q = (m - 2) d(i, j) - r(i) - r(j)`, where `m` is the number of
//! clusters left and `r(i)` the sum of row `i`, so it is two passes over the whole table of distances: the sums, then the
//! search. What this does differently from the plain way of writing it, without changing one bit of the answer:
//! * the table is one flat array that is never shifted: a join does not cut two rows and columns out of every row, it
//!   takes two clusters out of the list that says which rows are alive, in which order, and puts the new cluster in the
//!   place of one of them;
//! * the sums are made four rows at a time (each row's own sum is added in the same order as ever, but the four chains of
//!   additions overlap, where one chain waits for each addition before it can start the next);
//! * the sums and the search are shared among the cores when there are many clusters.
//!
//! The results are the same to the last bit because each row's sum is added in the order of the list, each `q` is worked out
//! with the same operations, and the pair taken is the first of the smallest in the order of the list, as the plain way
//! takes it (checked against a copy of the plain way, below, on every kind of matrix: ties, zeros of either sign, matrices that are not
//! symmetric, infinities and not-a-numbers).

use super::Join;
use crate::seq::par::{par_queue, worker};
use std::ops::Range;

/// From how many clusters the sums and the search are shared among the cores. Each join starts the threads twice (the sums, then
/// the search), about 0.4 ms for ten threads on the machine that this was measured on, and below this a join takes less than twice that.
const SHARE_FROM: usize = 800;

/// The joins of neighbour-joining on the processor, and the distances among the (two or three) clusters that are left.
/// `dist` is square.
pub(crate) fn cpu_joins(dist: &[Vec<f64>]) -> (Vec<Join>, Vec<Vec<f64>>) {
    run(dist, SHARE_FROM, cores())
}

/// How many cores to share the work among: one when this is already a thread of `par_map` (a bootstrap, say, that has a
/// tree for every core), or when there is only one.
fn cores() -> usize {
    if worker().is_some() {
        return 1;
    }
    std::thread::available_parallelism().map_or(1, |c| c.get())
}

/// How many runs of rows for each thread: the threads take the runs one by one, so a core that is slower does fewer of them.
const RUNS_FOR_EACH: usize = 4;

/// `cpu_joins`, sharing the work among `parts` threads from `share_from` clusters on (for the tests to choose both).
fn run(dist: &[Vec<f64>], share_from: usize, parts: usize) -> (Vec<Join>, Vec<Vec<f64>>) {
    let mut n = dist.len();
    assert!(dist.iter().all(|row| row.len() == n), "a square matrix");
    // row `s` of the table, `n` to a row, is the cluster in slot `s`; `order` holds the slots of the clusters that are left, in the
    // order of the list
    let mut d: Vec<f64> = dist.iter().flatten().copied().collect();
    let mut order: Vec<usize> = (0..n).collect();
    let mut joins = Vec::new();
    while order.len() > 3 {
        let size = order.len();
        let m = size as f64;
        let parts = if size >= share_from { parts } else { 1 };
        let r = row_sums(&d, n, &order, parts);
        let (i, j) = smallest_q(&d, n, &order, &r, m - 2.0, parts);
        let (si, sj) = (order[i], order[j]);
        let dij = d[si * n + sj];
        let li = (dij / 2.0 + (r[i] - r[j]) / (2.0 * (m - 2.0))).max(0.0);
        let lj = (dij - li).max(0.0);
        joins.push(Join { i, j, li, lj });
        // distance from the new node to every node that is left, in their current order; the new node goes to the end of the
        // list and into the slot of `j`
        let others: Vec<usize> = order
            .iter()
            .enumerate()
            .filter(|&(k, _)| k != i && k != j)
            .map(|(_, &s)| s)
            .collect();
        let to_new: Vec<f64> = others
            .iter()
            .map(|&sk| (d[si * n + sk] + d[sj * n + sk] - dij) / 2.0)
            .collect();
        for (&sk, &x) in others.iter().zip(&to_new) {
            d[sk * n + sj] = x;
            d[sj * n + sk] = x;
        }
        d[sj * n + sj] = 0.0;
        order = others;
        order.push(sj);
        // when a quarter of the slots are empty the clusters that are left lie scattered, a few in each cache line that a pass
        // reads: pack them into a smaller table, in the order of the list
        if order.len() * 4 <= n * 3 {
            let packed: Vec<f64> = order
                .iter()
                .flat_map(|&s| order.iter().map(move |&c| (s, c)))
                .map(|(s, c)| d[s * n + c])
                .collect();
            n = order.len();
            d = packed;
            order = (0..n).collect();
        }
    }
    let last = order
        .iter()
        .map(|&s| order.iter().map(|&c| d[s * n + c]).collect())
        .collect();
    (joins, last)
}

/// The sum of each row of the clusters in the list, each added in the order of the list.
fn row_sums(d: &[f64], n: usize, order: &[usize], parts: usize) -> Vec<f64> {
    let size = order.len();
    if parts <= 1 {
        return sums_of(d, n, order, 0..size);
    }
    let ranges = balanced(size, parts * RUNS_FOR_EACH, |_| 1);
    par_queue(&ranges, |rows| sums_of(d, n, order, rows.clone()))
        .into_iter()
        .flatten()
        .collect()
}

/// The sums of the rows at the positions `rows` of the list, four at a time.
fn sums_of(d: &[f64], n: usize, order: &[usize], rows: Range<usize>) -> Vec<f64> {
    // what `Iterator::sum` starts from, whichever it is in this version of Rust (it is -0.0 in the newer ones)
    let zero: f64 = std::iter::empty::<f64>().sum();
    let row = |p: usize| &d[order[p] * n..(order[p] + 1) * n];
    let mut sums = Vec::with_capacity(rows.len());
    let mut p = rows.start;
    while p + 4 <= rows.end {
        let (r0, r1, r2, r3) = (row(p), row(p + 1), row(p + 2), row(p + 3));
        let mut acc = [zero; 4];
        for &c in order {
            acc[0] += r0[c];
            acc[1] += r1[c];
            acc[2] += r2[c];
            acc[3] += r3[c];
        }
        sums.extend(acc);
        p += 4;
    }
    for p in p..rows.end {
        let r = row(p);
        sums.push(order.iter().map(|&c| r[c]).sum());
    }
    sums
}

/// The positions in the list of the pair with the smallest `q`: the first of them in the order of the list (row by row, `i` before
/// `j`) when several are equal, and the first two when none is below infinity.
fn smallest_q(
    d: &[f64],
    n: usize,
    order: &[usize],
    r: &[f64],
    m2: f64,
    parts: usize,
) -> (usize, usize) {
    let size = order.len();
    let ranges = if parts <= 1 {
        std::iter::once(0..size).collect()
    } else {
        // row `i` has `size - 1 - i` pairs
        balanced(size, parts * RUNS_FOR_EACH, |i| size - 1 - i)
    };
    let scan = |rows: &Range<usize>| {
        let mut best = (f64::INFINITY, 0, 1);
        for i in rows.clone() {
            let row = &d[order[i] * n..(order[i] + 1) * n];
            let ri = r[i];
            for j in i + 1..size {
                let q = m2 * row[order[j]] - ri - r[j];
                if q < best.0 {
                    best = (q, i, j);
                }
            }
        }
        best
    };
    let found = if ranges.len() == 1 {
        vec![scan(&ranges[0])]
    } else {
        par_queue(&ranges, scan)
    };
    // the ranges are in the order of the list, so the first that is smaller wins, as in one pass
    let mut best = (f64::INFINITY, 0, 1);
    for found in found {
        if found.0 < best.0 {
            best = found;
        }
    }
    (best.1, best.2)
}

/// Cut the rows `0..size` into at most `parts` runs that each take about the same work, a row `i` taking `weight(i)`.
fn balanced(size: usize, parts: usize, weight: impl Fn(usize) -> usize) -> Vec<Range<usize>> {
    let total: usize = (0..size).map(&weight).sum();
    let mut ranges = Vec::with_capacity(parts);
    let (mut start, mut done) = (0, 0);
    for row in 0..size {
        done += weight(row);
        if ranges.len() + 1 < parts && done * parts >= total * (ranges.len() + 1) {
            ranges.push(start..row + 1);
            start = row + 1;
        }
    }
    if start < size {
        ranges.push(start..size);
    }
    ranges
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::tree::Rng;

    // ---- the oracle: the joins as they were made before, copied unchanged ----

    fn join_rows_reference(d: &mut Vec<Vec<f64>>, i: usize, j: usize, to_new: Vec<f64>) {
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

    fn cpu_joins_reference(dist: &[Vec<f64>]) -> (Vec<Join>, Vec<Vec<f64>>) {
        let mut d: Vec<Vec<f64>> = dist.to_vec();
        let mut joins = Vec::new();
        while d.len() > 3 {
            let size = d.len();
            let m = size as f64;
            let r: Vec<f64> = d.iter().map(|row| row.iter().sum()).collect();
            let mut best = (f64::INFINITY, 0, 1);
            for i in 0..size {
                for j in i + 1..size {
                    let q = (m - 2.0) * d[i][j] - r[i] - r[j];
                    if q < best.0 {
                        best = (q, i, j);
                    }
                }
            }
            let (_, i, j) = best;
            let dij = d[i][j];
            let li = (dij / 2.0 + (r[i] - r[j]) / (2.0 * (m - 2.0))).max(0.0);
            let lj = (dij - li).max(0.0);
            let to_new: Vec<f64> = (0..size)
                .filter(|&k| k != i && k != j)
                .map(|k| (d[i][k] + d[j][k] - dij) / 2.0)
                .collect();
            joins.push(Join { i, j, li, lj });
            join_rows_reference(&mut d, i, j, to_new);
        }
        (joins, d)
    }

    /// Everything the joins hold, the numbers as their bits: equal keys are equal results to the last bit.
    type Key = (Vec<(usize, usize, u64, u64)>, Vec<Vec<u64>>);

    fn key((joins, last): &(Vec<Join>, Vec<Vec<f64>>)) -> Key {
        (
            joins
                .iter()
                .map(|j| (j.i, j.j, j.li.to_bits(), j.lj.to_bits()))
                .collect(),
            last.iter()
                .map(|row| row.iter().map(|x| x.to_bits()).collect())
                .collect(),
        )
    }

    /// A symmetric matrix from `f(i, j)` for `i < j`, with `diagonal` on the diagonal.
    fn symmetric(n: usize, diagonal: f64, mut f: impl FnMut(usize, usize) -> f64) -> Vec<Vec<f64>> {
        let upper: Vec<Vec<f64>> = (0..n)
            .map(|i| (i + 1..n).map(|j| f(i, j)).collect())
            .collect();
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| match i.cmp(&j) {
                        std::cmp::Ordering::Less => upper[i][j - i - 1],
                        std::cmp::Ordering::Greater => upper[j][i - j - 1],
                        std::cmp::Ordering::Equal => diagonal,
                    })
                    .collect()
            })
            .collect()
    }

    /// Matrices of every kind the joins meet, in every size from none to 130 clusters.
    fn matrices() -> Vec<(String, Vec<Vec<f64>>)> {
        let mut rng = Rng::new(31_415);
        let mut out = Vec::new();
        for n in [0usize, 1, 2, 3, 4, 5, 6, 7, 8, 10, 13, 20, 33, 50, 80, 130] {
            let points: Vec<(f64, f64)> = (0..n).map(|_| (rng.unit(), rng.unit())).collect();
            let line: Vec<f64> = (0..n).map(|_| rng.below(50) as f64).collect();
            let mut add = |kind: &str, d: Vec<Vec<f64>>| out.push((format!("{kind}, {n}"), d));
            add(
                "plane",
                symmetric(n, 0.0, |i, j| {
                    ((points[i].0 - points[j].0).powi(2) + (points[i].1 - points[j].1).powi(2))
                        .sqrt()
                }),
            );
            // whole numbers: a tie for almost every pair
            add(
                "whole numbers",
                symmetric(n, 0.0, |_, _| rng.below(5) as f64),
            );
            add("all equal", symmetric(n, 0.0, |_, _| 1.0));
            add("all zero", symmetric(n, 0.0, |_, _| 0.0));
            // identical sequences under Jukes-Cantor give a distance of -0.0
            add("negative zeros", symmetric(n, -0.0, |_, _| -0.0));
            add(
                "negative zeros and ones",
                symmetric(n, -0.0, |_, _| if rng.below(2) == 0 { -0.0 } else { 1.0 }),
            );
            // points on a line: additive, with many ties and some branches that come out negative
            add(
                "a line",
                symmetric(n, 0.0, |i, j| (line[i] - line[j]).abs()),
            );
            add("random", symmetric(n, 0.0, |_, _| 10.0 * rng.unit()));
            // not symmetric, with something on the diagonal: not a distance matrix, but the joins do not look
            add(
                "not symmetric",
                (0..n)
                    .map(|_| (0..n).map(|_| rng.unit() - 0.2).collect())
                    .collect(),
            );
            add(
                "infinities",
                symmetric(n, 0.0, |_, _| {
                    if rng.below(4) == 0 {
                        f64::INFINITY
                    } else {
                        rng.unit()
                    }
                }),
            );
            add(
                "not a number",
                symmetric(n, 0.0, |_, _| {
                    if rng.below(9) == 0 {
                        f64::NAN
                    } else {
                        rng.unit()
                    }
                }),
            );
            add(
                "tiny differences",
                symmetric(n, 0.0, |_, _| 1.0 + f64::EPSILON * rng.below(3) as f64),
            );
        }
        out
    }

    #[test]
    fn the_joins_equal_the_plain_way_to_the_last_bit_on_every_kind_of_matrix() {
        let all = matrices();
        assert!(all.len() > 150);
        for (what, dist) in &all {
            let want = key(&cpu_joins_reference(dist));
            // as it is called, then sharing the work among 2, 3, 4 and 7 threads from the smallest size, so the merging and
            // the cutting of the rows are tried on everything
            assert_eq!(key(&cpu_joins(dist)), want, "{what}");
            for parts in [1usize, 2, 3, 4, 7] {
                assert_eq!(key(&run(dist, 0, parts)), want, "{what}, {parts} parts");
            }
        }
    }

    #[test]
    fn many_clusters_are_joined_as_the_plain_way_does_when_the_work_is_shared() {
        let mut rng = Rng::new(7);
        // from 300 clusters on, then from the size that is used (a join of more than that many is shared)
        for (n, from) in [(450usize, 300usize), (SHARE_FROM + 1, SHARE_FROM)] {
            let points: Vec<(f64, f64)> = (0..n).map(|_| (rng.unit(), rng.unit())).collect();
            let dist = symmetric(n, 0.0, |i, j| {
                // rounded, so that there are ties among thousands of pairs
                (((points[i].0 - points[j].0).powi(2) + (points[i].1 - points[j].1).powi(2)).sqrt()
                    * 40.0)
                    .round()
            });
            assert_eq!(
                key(&run(&dist, from, 4)),
                key(&cpu_joins_reference(&dist)),
                "{n} clusters"
            );
        }
    }

    #[test]
    fn the_rows_are_cut_into_runs_that_cover_them_all_once_and_share_the_work() {
        for size in [0usize, 1, 2, 3, 4, 5, 9, 10, 11, 100, 1001] {
            for parts in [1usize, 2, 3, 4, 7, 10, 64] {
                for triangle in [false, true] {
                    let weight = |i: usize| if triangle { size - 1 - i } else { 1 };
                    let ranges = balanced(size, parts, weight);
                    assert!(ranges.len() <= parts.max(1), "{size} rows, {parts} parts");
                    let mut next = 0;
                    for range in &ranges {
                        assert_eq!(range.start, next);
                        assert!(range.end > range.start);
                        next = range.end;
                    }
                    assert_eq!(next, size);
                    // no run takes much more than its share (a run can overshoot by one row's weight)
                    let total: usize = (0..size).map(weight).sum();
                    let widest = (0..size).map(weight).max().unwrap_or(0);
                    for range in &ranges {
                        let work: usize = range.clone().map(weight).sum();
                        assert!(work <= total / parts.min(size.max(1)) + widest + 1);
                    }
                }
            }
        }
    }

    // ---- timing, not part of the suite: cargo test [--release] timing_of_the_joins -- --ignored --nocapture ----

    fn plane(n: usize, seed: u64) -> Vec<Vec<f64>> {
        let mut rng = Rng::new(seed);
        let points: Vec<(f64, f64)> = (0..n).map(|_| (rng.unit(), rng.unit())).collect();
        symmetric(n, 0.0, |i, j| {
            ((points[i].0 - points[j].0).powi(2) + (points[i].1 - points[j].1).powi(2)).sqrt()
        })
    }

    #[test]
    #[ignore = "timing"]
    fn timing_of_the_joins() {
        let time = |f: &dyn Fn()| {
            let start = std::time::Instant::now();
            f();
            start.elapsed().as_secs_f64()
        };
        for n in [300usize, 600, 1200, 2400, 4000] {
            let dist = plane(n, 9);
            // the work: about n cubed over 3 additions and n cubed over 6 values of q
            let steps = (n as f64).powi(3) / 3.0;
            let old = time(&|| {
                std::hint::black_box(cpu_joins_reference(&dist));
            });
            let one = time(&|| {
                std::hint::black_box(run(&dist, usize::MAX, 1));
            });
            let all = time(&|| {
                std::hint::black_box(cpu_joins(&dist));
            });
            println!(
                "{n:>5} clusters: before {old:7.3} s | one thread {one:7.3} s ({:.2} ns/step) | all cores {all:7.3} s",
                one * 1e9 / steps
            );
        }
    }
}
