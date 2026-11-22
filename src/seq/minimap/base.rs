//! The alignment of the bases along a chain: the anchors are exact matches, and what lies between two of them is aligned with the program's
//! own aligner (`seq::align`, affine gaps), then written as a CIGAR. The alignment begins at the first anchor and ends at the last: nothing
//! is added at the ends of the chain (minimap2 extends them with a drop-off rule; this does not).

use super::chain::Anchor;
use crate::seq::align::{align, Scoring, MAX_CELLS};

#[derive(Clone, Debug, PartialEq)]
pub struct Base {
    pub cigar: String,
    /// Columns where the two have the same base.
    pub matches: usize,
    /// Columns in all: matches, mismatches and the bases of the gaps.
    pub block: usize,
    /// Mismatches and the bases of the gaps.
    pub differences: usize,
}

/// The operations so far, in runs: `M` for a base against a base, `I` for a base of the read with none of the target, `D` the other way.
#[derive(Default)]
struct Ops {
    runs: Vec<(char, usize)>,
    matches: usize,
    block: usize,
}

impl Ops {
    fn push(&mut self, op: char, n: usize) {
        if n == 0 {
            return;
        }
        self.block += n;
        match self.runs.last_mut() {
            Some((last, count)) if *last == op => *count += n,
            _ => self.runs.push((op, n)),
        }
    }

    /// `n` columns of the read against the target, counting the equal ones.
    fn pairs(&mut self, q: &[u8], t: &[u8]) {
        self.matches += q
            .iter()
            .zip(t)
            .filter(|(a, b)| a.eq_ignore_ascii_case(b))
            .count();
        self.push('M', q.len());
    }

    /// The part of the read and of the target between two anchors.
    fn gap(&mut self, q: &[u8], t: &[u8]) {
        match (q.is_empty(), t.is_empty()) {
            (true, true) => {}
            (true, false) => self.push('D', t.len()),
            (false, true) => self.push('I', q.len()),
            _ if (q.len() + 1).saturating_mul(t.len() + 1) <= MAX_CELLS => {
                match align(q, t, &Scoring::default(), false) {
                    Ok(al) => self.columns(&al.a, &al.b),
                    Err(_) => self.rough(q, t),
                }
            }
            _ => self.rough(q, t),
        }
    }

    fn columns(&mut self, a: &[u8], b: &[u8]) {
        for (&x, &y) in a.iter().zip(b) {
            match (x == b'-', y == b'-') {
                (false, false) => {
                    self.matches += usize::from(x.eq_ignore_ascii_case(&y));
                    self.push('M', 1);
                }
                (true, _) => self.push('D', 1),
                (_, true) => self.push('I', 1),
            }
        }
    }

    /// A gap too large to align in one piece: the shorter side against the first of the longer, the rest a gap.
    fn rough(&mut self, q: &[u8], t: &[u8]) {
        let n = q.len().min(t.len());
        self.pairs(&q[..n], &t[..n]);
        self.push(
            if q.len() > n { 'I' } else { 'D' },
            q.len().max(t.len()) - n,
        );
    }
}

/// The alignment along `anchors` (in the order of the target), of the read `q` as it lies against the target `t` (the reverse complement
/// of the read for the other strand); anchors that do not fit the ones before them are passed over.
pub fn along(anchors: &[Anchor], q: &[u8], t: &[u8], k: usize) -> Option<Base> {
    let first = anchors.first()?;
    let mut ops = Ops::default();
    let (mut qe, mut te) = (first.y, first.x);
    for a in anchors {
        let (qs, ts) = (a.y, a.x);
        let (qend, tend) = ((qs + k).min(q.len()), (ts + k).min(t.len()));
        if qs >= qe && ts >= te {
            ops.gap(&q[qe..qs], &t[te..ts]);
            ops.pairs(&q[qs..qend], &t[ts..tend]);
            (qe, te) = (qend, tend);
        } else if qs + te == ts + qe && qend > qe && tend > te {
            // on the line of the block before it, and reaching past its end: the block grows
            ops.pairs(&q[qe..qend], &t[te..tend]);
            (qe, te) = (qend, tend);
        }
    }
    let differences = ops.block - ops.matches;
    let cigar = ops.runs.iter().map(|(op, n)| format!("{n}{op}")).collect();
    Some(Base {
        cigar,
        matches: ops.matches,
        block: ops.block,
        differences,
    })
}

#[cfg(test)]
mod tests;
