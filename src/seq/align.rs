//! Pairwise alignment with affine gaps (Gotoh). The original's Align tool scores match +1, mismatch -0.1, first gap position
//! -1 and each further gap position -0.2; here the same scores are ten times as large, in whole numbers (match +10, mismatch
//! -1, first gap -10, further gaps -2), which a 32-bit float holds exactly: every sum, on the processor and on the cards, is
//! exact, so the same alignment is found whatever the order of the additions, and a tie is a tie. The alignments are those of
//! the original's scores except where rounding of 0.1 and 0.2 in a float broke a tie (roadmap question 10).

use super::par::par_map;

/// What a match scores (the original's +1, ten times as large). It is not a setting: the others are in units of it.
pub const MATCH: f32 = 10.0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Scoring {
    pub mismatch: f32,
    pub gap_open: f32,
    pub gap_extend: f32,
}

impl Default for Scoring {
    fn default() -> Self {
        Scoring {
            mismatch: -1.0,
            gap_open: -10.0,
            gap_extend: -2.0,
        }
    }
}

/// Largest dynamic-programming table we are willing to fill.
pub const MAX_CELLS: usize = 40_000_000;

#[derive(Clone, Debug)]
pub struct Alignment {
    pub a: Vec<u8>,
    pub b: Vec<u8>,
    pub score: f32,
    pub matches: usize,
    /// Columns that count toward identity (end overhangs are excluded when they were free).
    pub columns: usize,
}

impl Alignment {
    pub fn identity(&self) -> f64 {
        if self.columns == 0 {
            0.0
        } else {
            self.matches as f64 / self.columns as f64
        }
    }
}

const NEG: f32 = f32::NEG_INFINITY;
const GAP: u8 = b'-';

fn same(x: u8, y: u8) -> bool {
    x.eq_ignore_ascii_case(&y)
}

/// Best predecessor among (M, X, Y) scores: value and index, preferring M then X then Y. M against X first, the winner
/// against Y after: the same answer as testing each against both others (also with a NaN in the way), in two selects the
/// compiler turns into conditional moves, with nothing for a processor to mispredict, and a shorter chain along a row.
#[inline(always)]
fn best3(m: f32, x: f32, y: f32) -> (f32, u8) {
    let (v, t) = if m >= x { (m, 0) } else { (x, 1) };
    if v >= y {
        (v, t)
    } else {
        (y, 2)
    }
}

fn check_size(n: usize, m: usize) -> Result<(), String> {
    if (n + 1).saturating_mul(m + 1) > MAX_CELLS {
        return Err(format!(
            "{n} x {m} bases is too large to align in one piece (limit about 6,000 x 6,000)."
        ));
    }
    Ok(())
}

/// Align `a` with `b`. With `free_ends`, overhanging ends cost nothing (overlap alignment).
pub fn align(a: &[u8], b: &[u8], sc: &Scoring, free_ends: bool) -> Result<Alignment, String> {
    check_size(a.len(), b.len())?;
    Ok(align_lanes::<1>(a, &[b], sc, free_ends).swap_remove(0))
}

/// How many pairs `align_many` works on in lock step: eight scores of 32 bits are two 128-bit vector registers.
const LANES: usize = 8;

/// The most bytes of traceback table (one byte per cell and pair) a group of `LANES` pairs may fill: 64 MB, a little more
/// than the 40 MB one pair may. A group over that is aligned one pair at a time.
const LANES_MAX_TRACE: usize = 64 << 20;
const _: () = assert!(LANES_MAX_TRACE / LANES <= MAX_CELLS);

/// Align `a` with every one of `bs`: the same as `align(a, b, sc, free_ends)` for each `b` (the same alignment to the last
/// bit, the same error for a pair too large), and several times faster for many pairs, because eight pairs share one pass
/// over the sequence `a`, a lane each, and the processor does the same step for the eight in one instruction. All cores
/// share the groups of eight; pairs of about the same length are grouped, since a group's table is as wide as its longest.
pub fn align_many(
    a: &[u8],
    bs: &[&[u8]],
    sc: &Scoring,
    free_ends: bool,
) -> Vec<Result<Alignment, String>> {
    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    if let Some(done) = on_devices(a, bs, sc, free_ends, |chosen| {
        crate::accel::align_pairs(a, chosen, sc, free_ends)
    }) {
        return done;
    }
    align_many_on_processor(a, bs, sc, free_ends)
}

/// The pairs on the PyTorch workers, where the person has asked for that (see `accel::align_pairs`: not unless they have, as the
/// processor is faster everywhere it has been measured): the same alignments to the last bit. `None` where the processor is to do them, and
/// when the workers fail (said on the standard error). The pairs the workers do not take (an empty sequence, a pair too large) are left
/// to the processor, which refuses or aligns them in its own words.
#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
fn on_devices(
    a: &[u8],
    bs: &[&[u8]],
    sc: &Scoring,
    free_ends: bool,
    devices: impl FnOnce(&[&[u8]]) -> Option<Vec<Alignment>>,
) -> Option<Vec<Result<Alignment, String>>> {
    let taken: Vec<usize> = (0..bs.len())
        .filter(|&i| !a.is_empty() && !bs[i].is_empty() && check_size(a.len(), bs[i].len()).is_ok())
        .collect();
    let chosen: Vec<&[u8]> = taken.iter().map(|&i| bs[i]).collect();
    let done = devices(&chosen)?;
    let mut out: Vec<Option<Result<Alignment, String>>> = (0..bs.len()).map(|_| None).collect();
    for (&i, alignment) in taken.iter().zip(done) {
        out[i] = Some(Ok(alignment));
    }
    let left: Vec<usize> = (0..bs.len()).filter(|&i| out[i].is_none()).collect();
    let left_bs: Vec<&[u8]> = left.iter().map(|&i| bs[i]).collect();
    for (&i, result) in left
        .iter()
        .zip(align_many_on_processor(a, &left_bs, sc, free_ends))
    {
        out[i] = Some(result);
    }
    Some(out.into_iter().flatten().collect())
}

/// `align_many` on the processor: groups of eight in lock step, on all cores (also what `gpu-check --bench` times the cards against).
#[cfg_attr(
    not(all(feature = "gpu", not(target_arch = "wasm32"))),
    allow(dead_code)
)]
pub(crate) fn align_many_on_processor(
    a: &[u8],
    bs: &[&[u8]],
    sc: &Scoring,
    free_ends: bool,
) -> Vec<Result<Alignment, String>> {
    let mut order: Vec<usize> = (0..bs.len()).collect();
    order.sort_by_key(|&i| bs[i].len());
    let groups: Vec<&[usize]> = order.chunks(LANES).collect();
    let done = par_map(&groups, |g| {
        align_group(a, bs, g, sc, free_ends, LANES_MAX_TRACE)
    });
    let mut out: Vec<Option<Result<Alignment, String>>> = (0..bs.len()).map(|_| None).collect();
    for (group, results) in groups.iter().zip(done) {
        for (&i, r) in group.iter().zip(results) {
            out[i] = Some(r);
        }
    }
    out.into_iter().flatten().collect()
}

/// The pairs `a` with `bs[i]` for `i` in `idx` (at most `LANES`), in that order; in lock step unless that takes more than
/// `cap` bytes of traceback or there is only one pair.
fn align_group(
    a: &[u8],
    bs: &[&[u8]],
    idx: &[usize],
    sc: &Scoring,
    free_ends: bool,
    cap: usize,
) -> Vec<Result<Alignment, String>> {
    let group: Vec<&[u8]> = idx.iter().map(|&i| bs[i]).collect();
    let widest = group.iter().map(|b| b.len()).max().unwrap_or(0) + 1;
    let bytes = (a.len() + 1).saturating_mul(widest).saturating_mul(LANES);
    if group.len() > 1 && bytes <= cap {
        align_lanes::<LANES>(a, &group, sc, free_ends)
            .into_iter()
            .map(Ok)
            .collect()
    } else {
        group.iter().map(|b| align(a, b, sc, free_ends)).collect()
    }
}

/// The scores of the three states at one cell of the table, for each of the `L` pairs: the best alignment ending with two
/// bases aligned (`m`), with a base of `a` against a gap (`x`), with a base of `b` against a gap (`y`).
#[derive(Clone, Copy)]
struct Cell<const L: usize> {
    m: [f32; L],
    x: [f32; L],
    y: [f32; L],
}

/// Where an alignment ends: its score, the row and column of its last cell, and the state there.
type End = (f32, usize, usize, u8);

/// Align `a` with each of the (at most `L`) sequences `bs` at once, lane `k` for `bs[k]`. Every lane does the sums and
/// the choices that `align` always did, in the same order, so what comes out is the same to the last bit; lanes whose
/// sequence is shorter than the longest have columns on the right that nothing reads. With `L = 1` this is the
/// aligner for one pair. The caller has checked that the tables are not too large.
fn align_lanes<const L: usize>(
    a: &[u8],
    bs: &[&[u8]],
    sc: &Scoring,
    free_ends: bool,
) -> Vec<Alignment> {
    let n = a.len();
    let lanes = bs.len();
    let lens: [usize; L] = std::array::from_fn(|k| bs.get(k).map_or(0, |b| b.len()));
    let w = lens.iter().fold(0, |most, &m| most.max(m)) + 1;
    let gap_run = |len: usize| {
        if free_ends {
            0.0
        } else {
            sc.gap_open + (len as f32 - 1.0) * sc.gap_extend
        }
    };
    // the bases compared without regard to case, one column of the table holding a base of each lane's `b`
    let lower: Vec<u8> = a.iter().map(u8::to_ascii_lowercase).collect();
    let mut cols = vec![[0u8; L]; w];
    for (k, b) in bs.iter().enumerate() {
        for (col, c) in cols.iter_mut().skip(1).zip(b.iter()) {
            col[k] = c.to_ascii_lowercase();
        }
    }
    let mut trace = vec![[0u8; L]; (n + 1) * w];
    // row 0, then the row last done: each row is written over the one before it
    let mut cells = vec![
        Cell {
            m: [NEG; L],
            x: [NEG; L],
            y: [NEG; L],
        };
        w
    ];
    cells[0].m = [0.0; L];
    for (j, cell) in cells.iter_mut().enumerate().skip(1) {
        cell.y = [gap_run(j); L];
    }
    let mut ends: [End; L] = std::array::from_fn(|k| (NEG, n, lens[k], 0u8));
    let consider = |ends: &mut [End; L], k: usize, i: usize, j: usize, c: &Cell<L>| {
        let (s, st) = best3(c.m[k], c.x[k], c.y[k]);
        if s > ends[k].0 {
            ends[k] = (s, i, j, st);
        }
    };
    if n == 0 {
        for k in 0..lanes {
            consider(&mut ends, k, 0, lens[k], &cells[lens[k]]);
        }
    }
    for i in 1..=n {
        let row = &mut trace[i * w..(i + 1) * w];
        fill_row(lower[i - 1], gap_run(i), &cols, &mut cells, row, sc);
        for k in 0..lanes {
            if free_ends {
                consider(&mut ends, k, i, lens[k], &cells[lens[k]]);
                if i == n {
                    (0..=lens[k]).for_each(|j| consider(&mut ends, k, n, j, &cells[j]));
                }
            } else if i == n {
                consider(&mut ends, k, n, lens[k], &cells[lens[k]]);
            }
        }
    }
    (0..lanes)
        .map(|k| trace_back(a, bs[k], &trace, w, k, ends[k], free_ends))
        .collect()
}

/// Row of the table for `a[i - 1] == ai` (lower case), for all lanes: `cells` holds the row above on the way in and this
/// row on the way out, `first_x` is the score of the cell at its left end, `trace` receives the choices made.
///
/// A function of its own and never inlined, because measured otherwise: inside the larger function the compiler leaves
/// the loop over the lanes as scalar code and the row takes three times as long as here, where it makes vector code.
#[inline(never)]
fn fill_row<const L: usize>(
    ai: u8,
    first_x: f32,
    cols: &[[u8; L]],
    cells: &mut [Cell<L>],
    trace: &mut [[u8; L]],
    sc: &Scoring,
) {
    let (open, extend, mismatch) = (sc.gap_open, sc.gap_extend, sc.mismatch);
    let mut diag = cells[0];
    let mut left = Cell {
        m: [NEG; L],
        x: [first_x; L],
        y: [NEG; L],
    };
    cells[0] = left;
    for ((cell, t), col) in cells[1..].iter_mut().zip(&mut trace[1..]).zip(&cols[1..]) {
        let up = *cell;
        let mut new = up;
        for k in 0..L {
            let s = if ai == col[k] { MATCH } else { mismatch };
            let (bm, tm) = best3(diag.m[k], diag.x[k], diag.y[k]);
            let (bx, tx) = best3(up.m[k] + open, up.x[k] + extend, up.y[k] + open);
            // same state order as the others (0 = M, 1 = X, 2 = Y) so the trace decodes right
            let (by, ty) = best3(left.m[k] + open, left.x[k] + open, left.y[k] + extend);
            new.m[k] = bm + s;
            new.x[k] = bx;
            new.y[k] = by;
            t[k] = tm | (tx << 2) | (ty << 4);
        }
        *cell = new;
        diag = up;
        left = new;
    }
}

/// The alignment of `a` with `b` from what a device made of their table (`align_batch` in `accel/torch/worker.py`): `trace`, the
/// bytes of the cells (`width` to a row, as `trace_back` reads them), `last_col`, the three scores at the last column of each row
/// (three for each row, from row 0, which is not used), and `last_row`, the three at each column of the last row. The end cell is
/// chosen in the order that `align_lanes` considers them (the last column row by row, then the last row column by column, the first
/// that is greater than all before), so the alignment is the processor's. Both sequences have a base.
#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
pub(crate) fn alignment_from_device(
    a: &[u8],
    b: &[u8],
    width: usize,
    trace: &[u8],
    last_col: &[f32],
    last_row: &[f32],
    free_ends: bool,
) -> Alignment {
    let (n, m) = (a.len(), b.len());
    let at = |cells: &[f32], k: usize| (cells[3 * k], cells[3 * k + 1], cells[3 * k + 2]);
    let mut end: End = (NEG, n, m, 0);
    let mut candidates: Vec<(usize, usize, (f32, f32, f32))> = Vec::new();
    if free_ends {
        candidates.extend((1..=n).map(|i| (i, m, at(last_col, i))));
        candidates.extend((0..=m).map(|j| (n, j, at(last_row, j))));
    } else {
        candidates.push((n, m, at(last_row, m)));
    }
    for (i, j, (cm, cx, cy)) in candidates {
        let (score, state) = best3(cm, cx, cy);
        if score > end.0 {
            end = (score, i, j, state);
        }
    }
    let table: Vec<[u8; 1]> = trace.iter().map(|&t| [t]).collect();
    trace_back(a, b, &table, width, 0, end, free_ends)
}

/// Read the alignment of lane `lane` back from the traceback table, `w` cells to a row, one byte for each lane in a cell.
fn trace_back<const L: usize>(
    a: &[u8],
    b: &[u8],
    trace: &[[u8; L]],
    w: usize,
    lane: usize,
    end: End,
    free_ends: bool,
) -> Alignment {
    let (score, mut i, mut j, mut state) = end;
    let (n, m) = (a.len(), b.len());
    let (mut ra, mut rb): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    // trailing overhang beyond the chosen end cell
    for k in (i..n).rev() {
        ra.push(a[k]);
        rb.push(GAP);
    }
    for k in (j..m).rev() {
        ra.push(GAP);
        rb.push(b[k]);
    }
    let overhang_end = ra.len();
    while i > 0 && j > 0 {
        let t = trace[i * w + j][lane];
        match state {
            0 => {
                ra.push(a[i - 1]);
                rb.push(b[j - 1]);
                state = t & 3;
                i -= 1;
                j -= 1;
            }
            1 => {
                ra.push(a[i - 1]);
                rb.push(GAP);
                state = (t >> 2) & 3;
                i -= 1;
            }
            _ => {
                ra.push(GAP);
                rb.push(b[j - 1]);
                state = (t >> 4) & 3;
                j -= 1;
            }
        }
    }
    let core_end = ra.len();
    while i > 0 {
        ra.push(a[i - 1]);
        rb.push(GAP);
        i -= 1;
    }
    while j > 0 {
        ra.push(GAP);
        rb.push(b[j - 1]);
        j -= 1;
    }
    let leading = ra.len() - core_end;
    ra.reverse();
    rb.reverse();
    let total = ra.len();
    let matches = ra
        .iter()
        .zip(&rb)
        .filter(|(x, y)| **x != GAP && **y != GAP && same(**x, **y))
        .count();
    let columns = if free_ends {
        total - overhang_end - leading
    } else {
        total
    };
    Alignment {
        a: ra,
        b: rb,
        score,
        matches,
        columns,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::tree::Rng;

    fn run(a: &str, b: &str, free: bool) -> Alignment {
        align(a.as_bytes(), b.as_bytes(), &Scoring::default(), free).unwrap()
    }

    #[test]
    fn identical_sequences_align_without_gaps() {
        let r = run("ACGTACGT", "ACGTACGT", false);
        assert_eq!((r.matches, r.columns), (8, 8));
        assert_eq!(r.score, 80.0);
        assert_eq!(r.a, r.b);
    }

    #[test]
    fn a_mismatch_costs_a_tenth_of_a_match() {
        let r = run("ACGTACGT", "ACGAACGT", false);
        assert_eq!(r.matches, 7);
        // 7 matches of 10 and a mismatch of -1: whole numbers, so exact
        assert_eq!(r.score, 69.0);
    }

    #[test]
    fn a_single_gap_costs_as_much_as_a_match_earns_in_one_column_over_the_gap() {
        let r = run("ACGTTACGT", "ACGTACGT", false);
        assert_eq!(r.matches, 8);
        // 8 matches of 10 and a gap of -10
        assert_eq!(r.score, 70.0);
        assert_eq!(r.a.len(), 9);
        assert_eq!(r.b.iter().filter(|&&c| c == b'-').count(), 1);
    }

    #[test]
    fn longer_gaps_use_the_cheaper_extension() {
        let r = run("ACGTAAAACGT", "ACGTCGT", false);
        // 7 matches of 10, one gap of 4 = -10 + 3 * -2 = -16
        assert_eq!(r.score, 70.0 - 16.0, "score {}", r.score);
    }

    #[test]
    fn a_gap_in_either_sequence_gives_an_alignment_that_reads_back_as_the_inputs() {
        let (short, long) = ("TTGACCAGGTCAGTTACGGATCC", "TTGACCAGGTCAGTTAAACGGATCC");
        for (first, second) in [(short, long), (long, short)] {
            for free in [false, true] {
                let r = run(first, second, free);
                let plain = |row: &[u8]| -> Vec<u8> {
                    row.iter().copied().filter(|&c| c != b'-').collect()
                };
                assert_eq!(plain(&r.a), first.as_bytes());
                assert_eq!(plain(&r.b), second.as_bytes());
                assert_eq!(r.a.len(), long.len());
                assert_eq!(r.a.iter().chain(&r.b).filter(|&&c| c == b'-').count(), 2);
            }
        }
    }

    #[test]
    fn free_ends_ignore_overhangs() {
        let r = run("TTTTTACGTACGTTTTT", "ACGTACGT", true);
        assert_eq!(r.matches, 8);
        assert!(
            (r.identity() - 1.0).abs() < 1e-9,
            "identity {}",
            r.identity()
        );
        let global = run("TTTTTACGTACGTTTTT", "ACGTACGT", false);
        assert!(global.score < r.score);
        assert_eq!(r.a.len(), r.b.len());
    }

    #[test]
    fn oversized_inputs_are_refused_with_a_message() {
        let big = vec![b'A'; 7000];
        assert!(align(&big, &big, &Scoring::default(), false).is_err());
    }

    #[test]
    fn empty_inputs_align_to_gaps() {
        let r = run("", "ACG", false);
        assert_eq!((r.a.as_slice(), r.b.as_slice()), (&b"---"[..], &b"ACG"[..]));
        let r = run("ACG", "", false);
        assert_eq!((r.a.as_slice(), r.b.as_slice()), (&b"ACG"[..], &b"---"[..]));
    }

    // ---- the oracle: the aligner as it was before the faster versions, copied unchanged ----

    fn best3_reference(m: f32, x: f32, y: f32) -> (f32, u8) {
        if m >= x && m >= y {
            (m, 0)
        } else if x >= y {
            (x, 1)
        } else {
            (y, 2)
        }
    }

    fn align_reference(
        a: &[u8],
        b: &[u8],
        sc: &Scoring,
        free_ends: bool,
    ) -> Result<Alignment, String> {
        let (n, m) = (a.len(), b.len());
        if (n + 1).saturating_mul(m + 1) > MAX_CELLS {
            return Err(format!(
                "{n} x {m} bases is too large to align in one piece (limit about 6,000 x 6,000)."
            ));
        }
        let w = m + 1;
        let gap_run = |len: usize| {
            if free_ends {
                0.0
            } else {
                sc.gap_open + (len as f32 - 1.0) * sc.gap_extend
            }
        };
        let mut trace = vec![0u8; (n + 1) * w];
        let (mut pm, mut px, mut py) = (vec![NEG; w], vec![NEG; w], vec![NEG; w]);
        pm[0] = 0.0;
        for (j, y) in py.iter_mut().enumerate().skip(1) {
            *y = gap_run(j);
        }
        let (mut cm, mut cx, mut cy) = (vec![NEG; w], vec![NEG; w], vec![NEG; w]);
        let mut end = (NEG, n, m, 0u8);
        let mut consider = |i: usize, j: usize, m: f32, x: f32, y: f32| {
            let (s, st) = best3_reference(m, x, y);
            if s > end.0 {
                end = (s, i, j, st);
            }
        };
        if n == 0 {
            consider(0, m, pm[m], px[m], py[m]);
        }
        for i in 1..=n {
            cm[0] = NEG;
            cy[0] = NEG;
            cx[0] = gap_run(i);
            for j in 1..=m {
                let s = if same(a[i - 1], b[j - 1]) {
                    MATCH
                } else {
                    sc.mismatch
                };
                let (bm, tm) = best3_reference(pm[j - 1], px[j - 1], py[j - 1]);
                cm[j] = bm + s;
                let (bx, tx) = best3_reference(
                    pm[j] + sc.gap_open,
                    px[j] + sc.gap_extend,
                    py[j] + sc.gap_open,
                );
                cx[j] = bx;
                let (by, ty) = best3_reference(
                    cm[j - 1] + sc.gap_open,
                    cx[j - 1] + sc.gap_open,
                    cy[j - 1] + sc.gap_extend,
                );
                cy[j] = by;
                trace[i * w + j] = tm | (tx << 2) | (ty << 4);
            }
            if free_ends {
                consider(i, m, cm[m], cx[m], cy[m]);
                if i == n {
                    (0..=m).for_each(|j| consider(n, j, cm[j], cx[j], cy[j]));
                }
            } else if i == n {
                consider(n, m, cm[m], cx[m], cy[m]);
            }
            std::mem::swap(&mut pm, &mut cm);
            std::mem::swap(&mut px, &mut cx);
            std::mem::swap(&mut py, &mut cy);
        }
        Ok(trace_back_reference(a, b, &trace, w, end, free_ends))
    }

    fn trace_back_reference(
        a: &[u8],
        b: &[u8],
        trace: &[u8],
        w: usize,
        end: (f32, usize, usize, u8),
        free_ends: bool,
    ) -> Alignment {
        let (score, mut i, mut j, mut state) = end;
        let (n, m) = (a.len(), b.len());
        let (mut ra, mut rb): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
        for k in (i..n).rev() {
            ra.push(a[k]);
            rb.push(GAP);
        }
        for k in (j..m).rev() {
            ra.push(GAP);
            rb.push(b[k]);
        }
        let overhang_end = ra.len();
        while i > 0 && j > 0 {
            let t = trace[i * w + j];
            match state {
                0 => {
                    ra.push(a[i - 1]);
                    rb.push(b[j - 1]);
                    state = t & 3;
                    i -= 1;
                    j -= 1;
                }
                1 => {
                    ra.push(a[i - 1]);
                    rb.push(GAP);
                    state = (t >> 2) & 3;
                    i -= 1;
                }
                _ => {
                    ra.push(GAP);
                    rb.push(b[j - 1]);
                    state = (t >> 4) & 3;
                    j -= 1;
                }
            }
        }
        let core_end = ra.len();
        while i > 0 {
            ra.push(a[i - 1]);
            rb.push(GAP);
            i -= 1;
        }
        while j > 0 {
            ra.push(GAP);
            rb.push(b[j - 1]);
            j -= 1;
        }
        let leading = ra.len() - core_end;
        ra.reverse();
        rb.reverse();
        let total = ra.len();
        let matches = ra
            .iter()
            .zip(&rb)
            .filter(|(x, y)| **x != GAP && **y != GAP && same(**x, **y))
            .count();
        let columns = if free_ends {
            total - overhang_end - leading
        } else {
            total
        };
        Alignment {
            a: ra,
            b: rb,
            score,
            matches,
            columns,
        }
    }

    /// Everything an `Alignment` holds, the score as its bits: equal keys are equal results to the last bit.
    fn key(r: &Alignment) -> (Vec<u8>, Vec<u8>, u32, usize, usize) {
        (
            r.a.clone(),
            r.b.clone(),
            r.score.to_bits(),
            r.matches,
            r.columns,
        )
    }

    fn random_bases(rng: &mut Rng, len: usize, letters: &[u8]) -> Vec<u8> {
        (0..len)
            .map(|_| letters[rng.below(letters.len())])
            .collect()
    }

    /// `base` with about one base in `rate` changed, and some short insertions and deletions.
    fn mutated(rng: &mut Rng, base: &[u8], rate: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(base.len() + 8);
        for &c in base {
            match rng.below(rate * 4) {
                0 => out.push(b"ACGT"[rng.below(4)]),
                1 => {} // a deletion
                2 => {
                    // an insertion of one to four bases
                    out.push(c);
                    let run = 1 + rng.below(4);
                    out.extend(random_bases(rng, run, b"ACGT"));
                }
                _ => out.push(c),
            }
        }
        out
    }

    /// A few hundred pairs of every kind the aligner meets: tiny, empty, related, repeats, poly-A, mixed case.
    fn oracle_pairs() -> Vec<(Vec<u8>, Vec<u8>)> {
        let mut rng = Rng::new(20_240_607);
        let mut pairs: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        let mut add = |a: &[u8], b: &[u8]| pairs.push((a.to_vec(), b.to_vec()));
        // lengths 0, 1, 2 and a spread of others, on DNA, on two letters (ties everywhere), on a protein alphabet
        let lengths = [
            0usize, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 33, 64, 100, 129,
        ];
        for (alphabet, rounds) in [
            (&b"ACGT"[..], 2),
            (&b"AC"[..], 2),
            (&b"ACGTacgtN"[..], 1),
            (&b"ACDEFGHIKLMNPQRSTVWY"[..], 1),
        ] {
            for _ in 0..rounds {
                for &la in &lengths {
                    for &lb in &[0usize, 1, 2, 6, 17, 40, 128] {
                        add(
                            &random_bases(&mut rng, la, alphabet),
                            &random_bases(&mut rng, lb, alphabet),
                        );
                    }
                }
            }
        }
        // related pairs: a copy with changes, insertions and deletions, in several lengths and rates
        for (len, rate) in [(30usize, 5usize), (80, 10), (200, 20), (200, 3), (400, 15)] {
            for _ in 0..6 {
                let base = random_bases(&mut rng, len, b"ACGT");
                let copy = mutated(&mut rng, &base, rate);
                add(&base, &copy);
                add(&copy, &base);
                add(&base, &base.to_ascii_lowercase());
                // one of them cut at an end: an overlap, as the free ends are for
                add(&base[len / 3..], &copy[..copy.len() * 2 / 3]);
            }
        }
        // ties: repeats, poly-A, a run against a longer run, a repeat against its own shift
        let acs: Vec<u8> = b"AC".iter().copied().cycle().take(60).collect();
        let tail: Vec<u8> = b"A".repeat(7);
        for (a, b) in [
            (b"AAAAAAAAAA".to_vec(), b"AAAAAAAAAAAAAAA".to_vec()),
            (b"A".repeat(50), b"A".repeat(33)),
            (
                b"A".repeat(40),
                [b"A".repeat(20), b"C".to_vec(), b"A".repeat(20)].concat(),
            ),
            (acs.clone(), acs[2..].to_vec()),
            (acs.clone(), acs[1..].to_vec()),
            (
                [b"GATTACA".repeat(9), tail.clone()].concat(),
                b"GATTACA".repeat(8),
            ),
            (b"ACGTACGTACGT".to_vec(), b"acgtACGTacgt".to_vec()),
            (b"TTTTT".to_vec(), b"AAAAA".to_vec()),
            (b"NNNNNNNN".to_vec(), b"NNNNNNNN".to_vec()),
        ] {
            add(&a, &b);
            add(&b, &a);
        }
        pairs
    }

    #[test]
    fn the_aligner_equals_the_reference_on_every_pair_whole() {
        let sc = Scoring::default();
        let pairs = oracle_pairs();
        assert!(pairs.len() >= 500, "{} pairs", pairs.len());
        for (a, b) in &pairs {
            for free in [false, true] {
                let want = align_reference(a, b, &sc, free).unwrap();
                let got = align(a, b, &sc, free).unwrap();
                assert_eq!(
                    key(&got),
                    key(&want),
                    "{} vs {} (free ends: {free})",
                    String::from_utf8_lossy(a),
                    String::from_utf8_lossy(b)
                );
            }
        }
    }

    #[test]
    fn eight_at_a_time_equals_the_reference_on_every_pair_whole() {
        let sc = Scoring::default();
        let mut rng = Rng::new(99);
        let pairs = oracle_pairs();
        // sets of sequences against one `a`: every `a` of the oracle with the `b`s that follow it, in several sizes of set
        // (one pair, a partial group, a full one, more than one group), lengths unequal, and some empty
        for (set, &size) in [1usize, 2, 3, 7, 8, 9, 16, 23]
            .iter()
            .cycle()
            .enumerate()
            .take(60)
        {
            let a = &pairs[rng.below(pairs.len())].0;
            let bs: Vec<&[u8]> = (0..size)
                .map(|_| pairs[rng.below(pairs.len())].1.as_slice())
                .collect();
            for free in [false, true] {
                let got = align_many(a, &bs, &sc, free);
                assert_eq!(got.len(), bs.len());
                for (b, got) in bs.iter().zip(got) {
                    let want = align_reference(a, b, &sc, free).unwrap();
                    assert_eq!(
                        key(&got.unwrap()),
                        key(&want),
                        "set {set}: {} vs {} (free ends: {free})",
                        String::from_utf8_lossy(a),
                        String::from_utf8_lossy(b)
                    );
                }
            }
        }
    }

    #[test]
    fn a_group_too_big_for_the_lanes_is_aligned_one_pair_at_a_time_with_the_same_results() {
        let sc = Scoring::default();
        let mut rng = Rng::new(5);
        let a = random_bases(&mut rng, 90, b"ACGT");
        let bs: Vec<Vec<u8>> = (0..5).map(|k| mutated(&mut rng, &a[k * 5..], 10)).collect();
        let refs: Vec<&[u8]> = bs.iter().map(Vec::as_slice).collect();
        let idx: Vec<usize> = (0..refs.len()).collect();
        for cap in [0, usize::MAX] {
            let got = align_group(&a, &refs, &idx, &sc, true, cap);
            for (b, got) in refs.iter().zip(got) {
                let want = align_reference(&a, b, &sc, true).unwrap();
                assert_eq!(key(&got.unwrap()), key(&want), "cap {cap}");
            }
        }
    }

    #[test]
    fn a_pair_too_large_is_refused_in_a_set_as_it_is_alone_and_the_rest_are_aligned() {
        let sc = Scoring::default();
        let a = vec![b'A'; 7000];
        let big = vec![b'A'; 7000];
        let ok = vec![b'A'; 30];
        let got = align_many(
            &a,
            &[ok.as_slice(), big.as_slice(), ok.as_slice()],
            &sc,
            true,
        );
        assert_eq!(got.len(), 3);
        assert_eq!(
            got[1].as_ref().err(),
            align(&a, &big, &sc, true).as_ref().err()
        );
        assert!(got[0].is_ok() && got[2].is_ok());
        assert!(align_many(&a, &[], &sc, true).is_empty());
    }

    /// Sets of sequences for the tests of the devices: short ones, each with a base, against short `a`s.
    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    fn device_sets(sets: usize, seed: u64) -> Vec<(Vec<u8>, Vec<Vec<u8>>)> {
        let mut rng = Rng::new(seed);
        let pairs: Vec<(Vec<u8>, Vec<u8>)> = oracle_pairs()
            .into_iter()
            .filter(|(a, b)| (1..=70).contains(&a.len()) && (1..=70).contains(&b.len()))
            .collect();
        [1usize, 2, 5, 9]
            .iter()
            .cycle()
            .take(sets)
            .map(|&size| {
                let a = pairs[rng.below(pairs.len())].0.clone();
                let bs = (0..size)
                    .map(|_| pairs[rng.below(pairs.len())].1.clone())
                    .collect();
                (a, bs)
            })
            .collect()
    }

    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    #[test]
    fn the_pytorch_workers_align_as_the_reference_to_the_last_bit() {
        use crate::accel::torch::{align::alignments, tests::pools};
        let sc = Scoring::default();
        for pool in pools() {
            for (a, bs) in device_sets(12, 404) {
                let refs: Vec<&[u8]> = bs.iter().map(Vec::as_slice).collect();
                for free in [false, true] {
                    let got = alignments(&pool, &a, &refs, &sc, free).expect("alignments");
                    assert_eq!(got.len(), refs.len());
                    for (b, got) in refs.iter().zip(&got) {
                        let want = align_reference(&a, b, &sc, free).unwrap();
                        assert_eq!(
                            key(got),
                            key(&want),
                            "{} workers: {} vs {} (free ends: {free})",
                            pool.len(),
                            String::from_utf8_lossy(&a),
                            String::from_utf8_lossy(b)
                        );
                    }
                }
            }
        }
    }

    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    #[test]
    fn a_dead_worker_makes_the_pytorch_alignments_fail_cleanly_and_the_processor_is_left_to_do_them(
    ) {
        use crate::accel::torch::{configured_python, tests::kill_worker, TorchPool};
        let Some(python) = configured_python() else {
            return;
        };
        let pool = std::sync::Arc::new(
            TorchPool::start_weighted(
                &python,
                Some(vec!["cpu".into(), "cpu".into()]),
                Some(vec![1.0, 1.0]),
            )
            .expect("workers")
            .expect("a pool"),
        );
        let (a, bs) = device_sets(1, 7).remove(0);
        let mut bs = bs;
        bs.push(b"ACGTAC".to_vec());
        bs.push(b"ACGTTC".to_vec());
        let refs: Vec<&[u8]> = bs.iter().map(Vec::as_slice).collect();
        kill_worker(&pool, 1);
        let sc = Scoring::default();
        assert!(crate::accel::torch::align::alignments(&pool, &a, &refs, &sc, true).is_err());
        // empty sequences are refused before anything is sent
        assert!(crate::accel::torch::align::alignments(&pool, b"", &refs, &sc, true).is_err());
    }

    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    #[test]
    fn pairs_the_workers_do_not_take_are_left_to_the_processor_and_everything_comes_back_in_order()
    {
        use crate::accel::torch::{align::alignments, tests::pools};
        let Some(pool) = pools().into_iter().last() else {
            return;
        };
        let sc = Scoring::default();
        let (a, mut bs) = device_sets(1, 3).remove(0);
        // an empty sequence and a pair too large for any table, among the ones that the workers take
        bs.insert(1, Vec::new());
        bs.push(vec![b'A'; 7000]);
        let big_a = vec![b'A'; 7000];
        let refs: Vec<&[u8]> = bs.iter().map(Vec::as_slice).collect();
        for centre in [&a, &big_a] {
            let calls = std::cell::Cell::new(0);
            let got = on_devices(centre, &refs, &sc, true, |chosen| {
                calls.set(calls.get() + 1);
                // only what the processor could align itself, none empty, is offered
                assert!(chosen.iter().all(|b| !b.is_empty()));
                assert!(chosen
                    .iter()
                    .all(|b| check_size(centre.len(), b.len()).is_ok()));
                Some(alignments(&pool, centre, chosen, &sc, true).expect("alignments"))
            });
            if centre == &big_a {
                // every pair of it is empty or too large except the shortest ones: the rest are still in order
                assert!(got.is_some() || calls.get() == 0);
            }
            let Some(got) = got else {
                continue;
            };
            assert_eq!(got.len(), refs.len());
            for (b, got) in refs.iter().zip(got) {
                match align(centre, b, &sc, true) {
                    Ok(want) => assert_eq!(key(&got.expect("aligned")), key(&want)),
                    Err(why) => assert_eq!(got.err(), Some(why)),
                }
            }
        }
    }

    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    #[test]
    #[ignore = "timing: HELIXLOOM_TORCH_PYTHON=... cargo test timing_of_the_pytorch_aligner -- --ignored --nocapture"]
    fn timing_of_the_pytorch_aligner() {
        use crate::accel::torch::{align::alignments, configured_python, TorchPool};
        let Some(python) = configured_python() else {
            return;
        };
        let sc = Scoring::default();
        let mut rng = Rng::new(3);
        for device in ["cpu", "mps"] {
            let Ok(Some(pool)) =
                TorchPool::start_weighted(&python, Some(vec![device.to_string()]), Some(vec![1.0]))
            else {
                println!("no {device} worker");
                continue;
            };
            for (count, len) in [(32usize, 300usize), (64, 500)] {
                let a = random_bases(&mut rng, len, b"ACGT");
                let bs: Vec<Vec<u8>> = (0..count)
                    .map(|_| random_bases(&mut rng, len, b"ACGT"))
                    .collect();
                let refs: Vec<&[u8]> = bs.iter().map(Vec::as_slice).collect();
                let cells = (count * (len + 1) * (len + 1)) as f64;
                // once to warm up the device (its kernels are made the first time), then timed
                let _ = alignments(&pool, &a, &refs, &sc, true).expect("alignments");
                let on_device = best_of(2, || {
                    std::hint::black_box(
                        alignments(&pool, &a, &refs, &sc, true).expect("alignments"),
                    );
                });
                let on_processor = best_of(2, || {
                    std::hint::black_box(align_many_on_processor(&a, &refs, &sc, true));
                });
                println!(
                    "{device}: {count} pairs of {len} x {len}: {:.1} ns a cell and pair on the worker, {:.2} on the processor's cores ({:.0} times slower)",
                    on_device * 1e9 / cells,
                    on_processor * 1e9 / cells,
                    on_device / on_processor
                );
            }
        }
    }

    #[test]
    fn the_devices_do_not_get_the_alignments_unless_asked() {
        // opt-in: the size from which they go is infinite until a person sets it
        #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
        assert!(crate::accel::threshold::ALIGN_CELLS.is_infinite());
    }

    // ---- timing, not part of the suite: cargo test [--release] timing_of -- --ignored --nocapture ----

    fn best_of(runs: usize, mut f: impl FnMut()) -> f64 {
        (0..runs)
            .map(|_| {
                let start = std::time::Instant::now();
                f();
                start.elapsed().as_secs_f64()
            })
            .fold(f64::MAX, f64::min)
    }

    #[test]
    #[ignore = "timing"]
    fn timing_of_the_pairwise_aligner() {
        let sc = Scoring::default();
        let mut rng = Rng::new(7);
        for len in [500usize, 1000, 2000, 4000] {
            let a = random_bases(&mut rng, len, b"ACGT");
            let unrelated = random_bases(&mut rng, len, b"ACGT");
            let related = mutated(&mut rng, &a, 20);
            let cells = ((len + 1) * (len + 1)) as f64;
            let runs = if len >= 2000 { 2 } else { 4 };
            for (what, b) in [("unrelated", &unrelated), ("related", &related)] {
                let cells = if what == "related" {
                    ((len + 1) * (b.len() + 1)) as f64
                } else {
                    cells
                };
                let t = best_of(runs, || {
                    std::hint::black_box(align(&a, b, &sc, true).unwrap());
                });
                println!(
                    "align        {len:>5} {what:<9} {:7.2} ns/cell ({:8.2} ms)",
                    t * 1e9 / cells,
                    t * 1e3
                );
            }
        }
    }

    /// One group of `LANES` unrelated pairs, so one thread, per cell and pair: `align_many` as it is, and the lock step
    /// forced (it is not used above `LANES_MAX_TRACE`).
    #[test]
    #[ignore = "timing"]
    fn timing_of_the_lock_step_aligner() {
        let sc = Scoring::default();
        let mut rng = Rng::new(7);
        for len in [500usize, 1000, 2000, 4000] {
            let a = random_bases(&mut rng, len, b"ACGT");
            let bs: Vec<Vec<u8>> = (0..LANES)
                .map(|_| random_bases(&mut rng, len, b"ACGT"))
                .collect();
            let refs: Vec<&[u8]> = bs.iter().map(Vec::as_slice).collect();
            let idx: Vec<usize> = (0..LANES).collect();
            let cells = ((len + 1) * (len + 1) * LANES) as f64;
            let runs = if len >= 2000 { 2 } else { 4 };
            let many = best_of(runs, || {
                std::hint::black_box(align_many(&a, &refs, &sc, true));
            });
            let forced = best_of(runs, || {
                std::hint::black_box(align_group(&a, &refs, &idx, &sc, true, usize::MAX));
            });
            println!(
                "align_many x{LANES} {len:>5}  {:7.2} ns/cell ({:8.2} ms for {LANES} pairs); lock step forced {:7.2} ns/cell",
                many * 1e9 / cells,
                many * 1e3,
                forced * 1e9 / cells
            );
        }
    }
}
