//! Pairwise alignment with affine gaps (Gotoh), scored like the original's Align tool:
//! match +1, mismatch -0.1, first gap position -1, each further gap position -0.2.

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Scoring {
    pub mismatch: f32,
    pub gap_open: f32,
    pub gap_extend: f32,
}

impl Default for Scoring {
    fn default() -> Self {
        Scoring {
            mismatch: -0.1,
            gap_open: -1.0,
            gap_extend: -0.2,
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

/// Best predecessor among (M, X, Y) scores: value and index, preferring M then X then Y.
fn best3(m: f32, x: f32, y: f32) -> (f32, u8) {
    if m >= x && m >= y {
        (m, 0)
    } else if x >= y {
        (x, 1)
    } else {
        (y, 2)
    }
}

/// Align `a` with `b`. With `free_ends`, overhanging ends cost nothing (overlap alignment).
pub fn align(a: &[u8], b: &[u8], sc: &Scoring, free_ends: bool) -> Result<Alignment, String> {
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
    // best end cell: (score, i, j, state)
    let mut end = (NEG, n, m, 0u8);
    let mut consider = |i: usize, j: usize, m: f32, x: f32, y: f32| {
        let (s, st) = best3(m, x, y);
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
                1.0
            } else {
                sc.mismatch
            };
            let (bm, tm) = best3(pm[j - 1], px[j - 1], py[j - 1]);
            cm[j] = bm + s;
            let (bx, tx) = best3(
                pm[j] + sc.gap_open,
                px[j] + sc.gap_extend,
                py[j] + sc.gap_open,
            );
            cx[j] = bx;
            // same state order as the others (0 = M, 1 = X, 2 = Y) so the trace decodes right
            let (by, ty) = best3(
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
    Ok(trace_back(a, b, &trace, w, end, free_ends))
}

fn trace_back(
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

#[cfg(test)]
mod tests {
    use super::*;

    fn run(a: &str, b: &str, free: bool) -> Alignment {
        align(a.as_bytes(), b.as_bytes(), &Scoring::default(), free).unwrap()
    }

    #[test]
    fn identical_sequences_align_without_gaps() {
        let r = run("ACGTACGT", "ACGTACGT", false);
        assert_eq!((r.matches, r.columns), (8, 8));
        assert!((r.score - 8.0).abs() < 1e-5);
        assert_eq!(r.a, r.b);
    }

    #[test]
    fn a_mismatch_costs_a_tenth() {
        let r = run("ACGTACGT", "ACGAACGT", false);
        assert_eq!(r.matches, 7);
        assert!((r.score - 6.9).abs() < 1e-5);
    }

    #[test]
    fn a_single_gap_costs_one() {
        let r = run("ACGTTACGT", "ACGTACGT", false);
        assert_eq!(r.matches, 8);
        assert!((r.score - 7.0).abs() < 1e-5);
        assert_eq!(r.a.len(), 9);
        assert_eq!(r.b.iter().filter(|&&c| c == b'-').count(), 1);
    }

    #[test]
    fn longer_gaps_use_the_cheaper_extension() {
        let r = run("ACGTAAAACGT", "ACGTCGT", false);
        // 7 matches, one gap of 4 = -1 + 3 * -0.2 = -1.6
        assert!((r.score - (7.0 - 1.6)).abs() < 1e-4, "score {}", r.score);
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
}
