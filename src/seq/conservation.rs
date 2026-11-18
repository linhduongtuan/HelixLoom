//! How conserved each column of an alignment is, by the scores that are used to find the residues that matter for a
//! protein's function: the Jensen-Shannon divergence of Capra and Singh (2007) and the Shannon entropy. The column's
//! residues are weighted so that a crowd of near copies does not outvote the rest (Henikoff and Henikoff 1994), the
//! score is lowered by the share of gaps, and a window of the neighbouring columns is mixed in, because conserved
//! residues lie close together. Written from those papers; see `CREDITS.md`.

use super::align::{align, Scoring};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// One minus the entropy of the column over its largest possible value.
    Entropy,
    /// The Jensen-Shannon divergence of the column from the background distribution of the alphabet.
    Jsd,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Settings {
    pub method: Method,
    /// Columns on each side that are mixed in (0: none).
    pub window: usize,
    /// How much of the score of the window goes into the score of a column, 0 to 1.
    pub lambda: f64,
    /// Lower the score of a column by the share of its rows that have a gap there.
    pub gap_penalty: bool,
    /// Weight the rows so that near copies count less.
    pub weights: bool,
    /// A column in which more than this share of the rows have a gap is not scored.
    pub max_gaps: f64,
}

impl Default for Settings {
    /// The choices of Capra and Singh: window of 3 and lambda 0.5, gap penalty on, and no score for a column that is
    /// more than 30 % gaps.
    fn default() -> Settings {
        Settings {
            method: Method::Jsd,
            window: 3,
            lambda: 0.5,
            gap_penalty: true,
            weights: true,
            max_gaps: 0.3,
        }
    }
}

// in the order of the background table below (BLOSUM's own: ... T W Y V)
const PROTEIN: &[u8] = b"ARNDCQEGHILKMFPSTWYV";
const DNA: &[u8] = b"ACGT";

/// The background frequencies of the amino acids in the BLOSUM62 matrix, in the order of `PROTEIN`, as Capra and Singh use
/// them (they are made to sum to 1 below).
const BLOSUM62_BACKGROUND: [f64; 20] = [
    0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024,
    0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072,
];

/// Added to every count so that no frequency is zero (as Capra and Singh do).
const PSEUDOCOUNT: f64 = 1e-7;

fn alphabet(protein: bool) -> &'static [u8] {
    if protein {
        PROTEIN
    } else {
        DNA
    }
}

/// The place of each byte in `alphabet` (upper or lower case), or -1; with `u_is_t`, U is T.
const fn index_table(alphabet: &[u8], u_is_t: bool) -> [i8; 256] {
    let mut table = [-1i8; 256];
    let mut i = 0;
    while i < alphabet.len() {
        table[alphabet[i] as usize] = i as i8;
        table[(alphabet[i] + 32) as usize] = i as i8;
        i += 1;
    }
    if u_is_t {
        table[b'U' as usize] = table[b'T' as usize];
        table[b'u' as usize] = table[b'T' as usize];
    }
    table
}

static PROTEIN_INDEX: [i8; 256] = index_table(PROTEIN, false);
static DNA_INDEX: [i8; 256] = index_table(DNA, true);

/// The place of a letter in the alphabet; `None` for a gap and for any letter that is not a residue of it (X, B, Z, N...),
/// which count as gaps. U is T.
fn symbol(c: u8, protein: bool) -> Option<usize> {
    let table = if protein { &PROTEIN_INDEX } else { &DNA_INDEX };
    usize::try_from(table[usize::from(c)]).ok()
}

fn background(protein: bool) -> Vec<f64> {
    if protein {
        let total: f64 = BLOSUM62_BACKGROUND.iter().sum();
        BLOSUM62_BACKGROUND.iter().map(|b| b / total).collect()
    } else {
        vec![0.25; DNA.len()]
    }
}

/// The weight of each row, summing to 1 (position-based weights: in each column a row gets 1 / (the number of different
/// residues there times the number of rows that have its residue), and the sums are made to total 1).
pub fn weights(rows: &[Vec<u8>], protein: bool) -> Vec<f64> {
    let n = rows.len();
    let width = rows.iter().map(Vec::len).max().unwrap_or(0);
    let mut w = vec![0.0; n];
    for k in 0..width {
        let syms: Vec<Option<usize>> = rows
            .iter()
            .map(|r| r.get(k).and_then(|&c| symbol(c, protein)))
            .collect();
        let mut counts = vec![0usize; alphabet(protein).len()];
        for &s in syms.iter().flatten() {
            counts[s] += 1;
        }
        let distinct = counts.iter().filter(|&&c| c > 0).count();
        for (wi, s) in w.iter_mut().zip(&syms) {
            if let Some(s) = s {
                *wi += 1.0 / (distinct * counts[*s]) as f64;
            }
        }
    }
    let total: f64 = w.iter().sum();
    if total <= 0.0 {
        vec![1.0 / n.max(1) as f64; n]
    } else {
        w.iter().map(|x| x / total).collect()
    }
}

/// The score of one column, before the window: `syms` are its rows' residues, `w` their weights.
fn column_score(syms: &[Option<usize>], w: &[f64], q: &[f64], settings: &Settings) -> Option<f64> {
    let k = q.len();
    let gap: f64 = syms
        .iter()
        .zip(w)
        .filter(|(s, _)| s.is_none())
        .map(|(_, w)| w)
        .sum();
    // the cutoff is the share of the rows that have a gap; only the penalty counts them by weight
    let gap_rows = syms.iter().filter(|s| s.is_none()).count();
    if gap_rows as f64 / syms.len().max(1) as f64 > settings.max_gaps + 1e-12 {
        return None;
    }
    let mut counts = vec![0.0; k];
    for (s, w) in syms.iter().zip(w) {
        if let Some(s) = s {
            counts[*s] += w;
        }
    }
    let total: f64 = counts.iter().sum();
    if total <= 0.0 {
        return None;
    }
    let p: Vec<f64> = counts
        .iter()
        .map(|c| (c + PSEUDOCOUNT) / (total + k as f64 * PSEUDOCOUNT))
        .collect();
    let score = match settings.method {
        Method::Entropy => {
            let entropy: f64 = -p.iter().map(|&x| x * x.log2()).sum::<f64>();
            1.0 - entropy / (k as f64).log2()
        }
        Method::Jsd => {
            // the mean of the two divergences from the mixture of the column and the background
            let divergence = |a: &[f64]| -> f64 {
                a.iter()
                    .zip(p.iter().zip(q))
                    .map(|(&x, (&pi, &qi))| {
                        let r = (pi + qi) / 2.0;
                        if x > 0.0 {
                            x * (x / r).log2()
                        } else {
                            0.0
                        }
                    })
                    .sum()
            };
            0.5 * divergence(&p) + 0.5 * divergence(q)
        }
    };
    Some(if settings.gap_penalty {
        score * (1.0 - gap)
    } else {
        score
    })
}

/// Mix the scores of the columns within `window` on each side into each score, `lambda` of the weight to them (their
/// mean) and the rest to the column itself. Columns without a score stay without and are not counted among the neighbours;
/// a column with no neighbour that has a score keeps its own. (Differences from the authors' program, which as far as is
/// known leaves the first and last `window` columns as they are, and weighs `lambda` the other way round: nothing at the
/// 0.5 used here, and for the ends see `CREDITS.md`.)
fn smooth(scores: &[Option<f64>], window: usize, lambda: f64) -> Vec<Option<f64>> {
    scores
        .iter()
        .enumerate()
        .map(|(i, s)| {
            let own = (*s)?;
            let lo = i.saturating_sub(window);
            let hi = (i + window + 1).min(scores.len());
            let near: Vec<f64> = (lo..hi)
                .filter(|&j| j != i)
                .filter_map(|j| scores[j])
                .collect();
            Some(if near.is_empty() {
                own
            } else {
                (1.0 - lambda) * own + lambda * near.iter().sum::<f64>() / near.len() as f64
            })
        })
        .collect()
}

/// The score of every column (0 to 1, higher is more conserved), or `None` for a column that is not scored. Fewer than two
/// rows have nothing to compare, so no column is scored.
pub fn scores(rows: &[Vec<u8>], protein: bool, settings: &Settings) -> Vec<Option<f32>> {
    let width = rows.iter().map(Vec::len).max().unwrap_or(0);
    if rows.len() < 2 {
        return vec![None; width];
    }
    let w = if settings.weights {
        weights(rows, protein)
    } else {
        vec![1.0 / rows.len() as f64; rows.len()]
    };
    let q = background(protein);
    let raw: Vec<Option<f64>> = (0..width)
        .map(|k| {
            let syms: Vec<Option<usize>> = rows
                .iter()
                .map(|r| r.get(k).and_then(|&c| symbol(c, protein)))
                .collect();
            column_score(&syms, &w, &q, settings)
        })
        .collect();
    let mixed = if settings.window > 0 {
        smooth(&raw, settings.window, settings.lambda)
    } else {
        raw
    };
    mixed
        .into_iter()
        .map(|s| s.map(|v| v.clamp(0.0, 1.0) as f32))
        .collect()
}

/// A run of conserved columns, first to last (both counted), and the mean score in it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Block {
    pub from: usize,
    pub to: usize,
    pub mean: f32,
}

/// The runs of at least `min_len` columns in a row whose scores are all at least `threshold` (a column with no score ends a run).
pub fn blocks(scores: &[Option<f32>], threshold: f32, min_len: usize) -> Vec<Block> {
    let mut found = Vec::new();
    let mut start = 0;
    let mut sum = 0.0;
    // one step past the end closes a run that reaches it
    for k in 0..=scores.len() {
        match scores.get(k).copied().flatten().filter(|&s| s >= threshold) {
            Some(s) => {
                if k == 0 || scores[k - 1].is_none_or(|p| p < threshold) {
                    (start, sum) = (k, 0.0);
                }
                sum += s;
            }
            // a run ends here, and is long enough
            None if k > start
                && k > 0
                && scores[k - 1].is_some_and(|p| p >= threshold)
                && k - start >= min_len.max(1) =>
            {
                found.push(Block {
                    from: start,
                    to: k - 1,
                    mean: sum / (k - start) as f32,
                });
            }
            None => {}
        }
    }
    found
}

/// The score that a quarter of the scored columns reach or pass: what makes a column "conserved" for this alignment.
pub fn top_quarter(scores: &[Option<f32>]) -> Option<f32> {
    let mut scored: Vec<f32> = scores.iter().flatten().copied().collect();
    if scored.is_empty() {
        return None;
    }
    scored.sort_by(|a, b| b.total_cmp(a));
    let quarter = scored.len().div_ceil(4).max(1);
    Some(scored[quarter - 1])
}

/// The column of each residue of `row` (its letters that are not gaps), in order.
pub fn columns_of(row: &[u8]) -> Vec<usize> {
    row.iter()
        .enumerate()
        .filter(|(_, &c)| !super::alnstats::is_gap(c))
        .map(|(k, _)| k)
        .collect()
}

/// How the residues of a structure's chain line up with those of one row of an alignment.
#[derive(Clone, Debug, PartialEq)]
pub struct Mapping {
    /// For each residue of the chain, the residue of the row (counted without gaps) that it matches, if any.
    pub to_row: Vec<Option<usize>>,
    /// The share of the matched residues that are the same letter, 0 to 1.
    pub identity: f64,
    /// How many residues of the chain are matched with one of the row, and how many of those are the same letter: a row that
    /// covers more of the chain is the better one even if a short piece of another is more alike.
    pub matched: usize,
    pub same: usize,
}

/// Line a chain up with a row (`row` without its gaps), or `None` when they are not the same protein: less than 70 %
/// identical, or matched over less than half of the shorter.
pub fn map_to_chain(row: &[u8], chain: &[u8]) -> Option<Mapping> {
    if row.is_empty() || chain.is_empty() {
        return None;
    }
    // overhangs cost nothing: a tag on the chain, or a stretch of the row that the experiment did not see
    let aligned = align(chain, row, &Scoring::default(), true).ok()?;
    let mut to_row = vec![None; chain.len()];
    let (mut c, mut r) = (0, 0);
    let (mut matched, mut same) = (0usize, 0usize);
    for (&x, &y) in aligned.a.iter().zip(&aligned.b) {
        let (has_chain, has_row) = (x != b'-', y != b'-');
        if has_chain && has_row {
            to_row[c] = Some(r);
            matched += 1;
            same += usize::from(x.eq_ignore_ascii_case(&y));
        }
        c += usize::from(has_chain);
        r += usize::from(has_row);
    }
    let identity = same as f64 / matched.max(1) as f64;
    (matched * 2 >= row.len().min(chain.len()) && identity >= 0.7).then_some(Mapping {
        to_row,
        identity,
        matched,
        same,
    })
}

#[cfg(test)]
mod tests;
