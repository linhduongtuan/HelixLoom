//! Chaining: from the places where a read and a target share a minimizer (anchors) to runs of them that lie on one line, as minimap2's
//! paper describes (Li 2018, section 2.1): a dynamic programme over the anchors in the order of their place on the target, in which an
//! anchor joins the best earlier one, gaining the bases it adds (at most the length of the k-mer) and paying for the difference between the
//! gap on the read and the gap on the target, `0.01 * k * |l| + 0.5 * log2 |l|`. Only the 50 anchors before it that did not improve the
//! score are looked at (the paper's `h`). Chains are then taken from the best end backwards, each anchor going to one chain.

/// A shared k-mer: where it starts on the target and on the read (the read, for the reverse strand, counted on its reverse complement).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Anchor {
    pub x: usize,
    pub y: usize,
}

#[derive(Clone, Copy, Debug)]
pub struct Params {
    pub k: usize,
    /// The most that two anchors of a chain may be apart on either sequence.
    pub max_gap: usize,
    /// The most that the gap on one sequence may differ from the gap on the other.
    pub bandwidth: usize,
    /// Anchors looked at before one is given up on without a better score.
    pub lookback: usize,
    pub min_anchors: usize,
    pub min_score: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Chain {
    pub score: f64,
    /// In the order of the target.
    pub anchors: Vec<Anchor>,
}

/// What a difference of `l` bases between the gaps costs.
pub fn gap_cost(l: usize, k: usize) -> f64 {
    if l == 0 {
        0.0
    } else {
        0.01 * k as f64 * l as f64 + 0.5 * (l as f64).log2()
    }
}

/// The chains of `anchors` (of one target and one strand), best first; those with fewer anchors or a lower score than the limits are left out.
pub fn chain(anchors: &[Anchor], p: &Params) -> Vec<Chain> {
    let mut a: Vec<Anchor> = anchors.to_vec();
    a.sort_unstable();
    a.dedup();
    let n = a.len();
    let (mut f, mut from): (Vec<f64>, Vec<Option<usize>>) = (vec![0.0; n], vec![None; n]);
    for (i, ai) in a.iter().enumerate() {
        let (mut best, mut pred, mut idle) = (p.k as f64, None, 0);
        for j in (0..i).rev() {
            let dx = ai.x - a[j].x;
            if dx > p.max_gap {
                break;
            }
            idle += 1;
            let dy =
                ai.y.checked_sub(a[j].y)
                    .filter(|&d| d > 0 && d <= p.max_gap);
            if let Some(dy) = dy.filter(|&dy| dx.abs_diff(dy) <= p.bandwidth) {
                let gained = dx.min(dy).min(p.k) as f64;
                let score = f[j] + gained - gap_cost(dx.abs_diff(dy), p.k);
                if score > best {
                    (best, pred, idle) = (score, Some(j), 0);
                }
            }
            if idle > p.lookback {
                break;
            }
        }
        (f[i], from[i]) = (best, pred);
    }
    // from the best end back: an anchor that is used goes to the chain that used it first
    let mut ends: Vec<usize> = (0..n).collect();
    ends.sort_by(|&i, &j| f[j].total_cmp(&f[i]).then(j.cmp(&i)));
    let mut used = vec![false; n];
    let mut chains = Vec::new();
    for end in ends {
        if used[end] {
            continue;
        }
        let mut path = Vec::new();
        let mut at = Some(end);
        while let Some(t) = at.filter(|&t| !used[t]) {
            path.push(t);
            at = from[t];
        }
        // a chain that runs into one taken before it counts only what it adds
        let score = f[end] - at.map_or(0.0, |t| f[t]);
        path.iter().for_each(|&t| used[t] = true);
        if path.len() >= p.min_anchors && score >= p.min_score {
            path.reverse();
            chains.push(Chain {
                score,
                anchors: path.into_iter().map(|t| a[t]).collect(),
            });
        }
    }
    chains.sort_by(|x, y| y.score.total_cmp(&x.score));
    chains
}

#[cfg(test)]
mod tests;
