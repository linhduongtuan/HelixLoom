//! Distances between all the sequences of an alignment on the GPUs. For every pair the processor counts the columns where both
//! have a base, and how many of those differ by a transition and by a transversion (`seq::tree::aligned_distance`); here the
//! same three numbers for every pair come from three matrix products of the columns written as 0 and 1:
//!
//! * `valid · validᵀ` is the number of columns where both have a base (a gap or an ambiguity code is 0),
//! * `onehot · onehotᵀ` (A, C, G and T each a place) is the number of those where they have the same base,
//! * `class · classᵀ` (purine or pyrimidine) is the number where they are of the same class, the same base included.
//!
//! so a pair differs by a transition in `class - onehot` columns and by a transversion in `valid - class`. A product of 0 and 1
//! is a whole number that a 32-bit float holds exactly (up to 2^24 columns), so the counts are exactly the processor's, and
//! the distance is worked out from them by the same function (`distance_from_counts`): the two agree to the last digit.
//!
//! A bootstrap resample is the same columns with each taken as many times as it was drawn, which is a weight for each column
//! on the left side of the products: the replicates of a bootstrap all reuse the columns that are on the cards.

use super::pool::Pool;
use crate::seq::tree::{distance_from_counts, DistanceModel};
use candle_core::{DType, Device, Tensor};

/// A 32-bit float holds a whole number exactly up to here, and so the number of columns, with their weights, may not pass it.
pub(super) const MOST_COLUMNS: usize = 1 << 24;
/// The heaviest column of a resample: a product with a weight must stay exact on a card that rounds its inputs to a short
/// mantissa for speed (TF32 keeps whole numbers to 2048).
pub(super) const MOST_WEIGHT: u32 = 1024;
/// What the columns, written out, may take on one card, in bytes (the three products' left and right sides).
pub(super) const RESIDENT_LIMIT: usize = 4 << 30;
/// About this many bytes of one-hot columns are made and sent at a time.
const CHUNK_BYTES: usize = 256 << 20;
pub(super) const FLOATS_PER_COLUMN: usize = 7;

fn bad(e: candle_core::Error) -> String {
    format!("GPU: {e}")
}

/// The base of a letter: A, C, G, T (either case) as 0 to 3; anything else (a gap, N, ...) is not compared.
fn base(c: u8) -> Option<usize> {
    match c.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// For each pair of sequences, how many columns were compared and how many of those differ by a transition and by a
/// transversion. `n` by `n`, row by row.
#[derive(Debug)]
pub struct Counts {
    pub n: usize,
    pub sites: Vec<u32>,
    pub transitions: Vec<u32>,
    pub transversions: Vec<u32>,
}

impl Counts {
    /// The counts from the three products, each `n` by `n` row by row and holding whole numbers: how many columns both sequences
    /// have a base in, in how many they have the same base, and in how many the same class of base (the same base included).
    pub fn from_products(n: usize, valid: &[f32], identical: &[f32], same_class: &[f32]) -> Counts {
        let whole = |x: &f32| x.round() as u32;
        let (sites, identical, same_class): (Vec<u32>, Vec<u32>, Vec<u32>) = (
            valid.iter().map(whole).collect(),
            identical.iter().map(whole).collect(),
            same_class.iter().map(whole).collect(),
        );
        Counts {
            n,
            transitions: same_class
                .iter()
                .zip(&identical)
                .map(|(&c, &e)| c.saturating_sub(e))
                .collect(),
            transversions: sites
                .iter()
                .zip(&same_class)
                .map(|(&v, &c)| v.saturating_sub(c))
                .collect(),
            sites,
        }
    }

    /// The distances of every pair under `model`, or the first pair (the smaller number first) that cannot be compared: no column
    /// has a base in both, or they differ too much for the model.
    pub fn distances(&self, model: DistanceModel) -> Result<Vec<Vec<f64>>, (usize, usize)> {
        let n = self.n;
        let mut d = vec![vec![0.0; n]; n];
        for (i, j) in (0..n).flat_map(|i| (i + 1..n).map(move |j| (i, j))) {
            let at = i * n + j;
            let v = distance_from_counts(
                self.sites[at],
                self.transitions[at],
                self.transversions[at],
                model,
            )
            .ok_or((i, j))?;
            (d[i][j], d[j][i]) = (v, v);
        }
        Ok(d)
    }
}

/// The number of sequences and of columns of `rows`, or why the cards cannot count them: too few, not all the same length, more
/// columns than a 32-bit float counts exactly, or written out as 0 and 1 more than a card is given.
pub(super) fn check_rows(rows: &[Vec<u8>]) -> Result<(usize, usize), String> {
    let n = rows.len();
    let columns = rows.first().map_or(0, Vec::len);
    if n < 2 || columns == 0 {
        return Err("There are fewer than two sequences or no column to compare.".into());
    }
    if rows.iter().any(|r| r.len() != columns) {
        return Err("The rows of the alignment are not all the same length.".into());
    }
    if columns > MOST_COLUMNS {
        return Err(format!(
            "{columns} columns are more than the GPUs count exactly ({MOST_COLUMNS})."
        ));
    }
    let bytes = n * columns * FLOATS_PER_COLUMN * 4;
    if bytes > RESIDENT_LIMIT {
        return Err(format!(
            "The alignment written out needs {} MB on each card, more than the {} MB that are used.",
            bytes >> 20,
            RESIDENT_LIMIT >> 20
        ));
    }
    Ok((n, columns))
}

/// Whether `weights` (how many times each column is taken) can be counted exactly.
pub(super) fn check_weights(weights: &[u32], columns: usize) -> Result<(), String> {
    if weights.len() != columns {
        return Err("There is not one weight for each column.".into());
    }
    if weights.iter().any(|&x| x > MOST_WEIGHT) {
        return Err(format!("A column is taken more than {MOST_WEIGHT} times."));
    }
    if weights.iter().map(|&x| x as usize).sum::<usize>() > MOST_COLUMNS {
        return Err("The weights add up to more than the GPUs count exactly.".into());
    }
    Ok(())
}

/// One chunk of the columns, on one card.
struct Chunk {
    /// The columns of the whole alignment that this chunk holds.
    columns: std::ops::Range<usize>,
    /// (sequences, 4 places for each column), (sequences, 2 for each), (sequences, 1 for each).
    onehot: Tensor,
    class: Tensor,
    valid: Tensor,
}

/// The alignment, written out as 0 and 1 on every card of a pool, ready for any number of weightings of its columns.
pub struct Engine {
    pool: Pool,
    n: usize,
    columns: usize,
    /// For each member of the pool, its copy of the chunks.
    resident: Vec<Vec<Chunk>>,
}

impl std::fmt::Debug for Engine {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Engine")
            .field("sequences", &self.n)
            .field("columns", &self.columns)
            .field("cards", &self.pool.len())
            .finish()
    }
}

impl Engine {
    /// Put the rows of an alignment (equal in length, `-` for a gap) on every card of `pool`.
    pub fn new(pool: &Pool, rows: &[Vec<u8>]) -> Result<Engine, String> {
        Engine::with_chunk(pool, rows, None)
    }

    /// The same with the number of columns made and sent at a time given (for the tests, which want many small chunks).
    pub fn with_chunk(
        pool: &Pool,
        rows: &[Vec<u8>],
        chunk_columns: Option<usize>,
    ) -> Result<Engine, String> {
        let (n, columns) = check_rows(rows)?;
        let per_chunk = chunk_columns
            .unwrap_or_else(|| (CHUNK_BYTES / (n * FLOATS_PER_COLUMN * 4)).clamp(64, 8192))
            .max(1);
        let mut resident: Vec<Vec<Chunk>> = pool.members().iter().map(|_| Vec::new()).collect();
        let mut start = 0;
        while start < columns {
            let end = (start + per_chunk).min(columns);
            let (onehot, class, valid) = write_out(rows, start..end);
            let width = end - start;
            for (member, chunks) in pool.members().iter().zip(&mut resident) {
                let device: &Device = &member.device;
                let _card = member.exclusive();
                chunks.push(Chunk {
                    columns: start..end,
                    onehot: Tensor::from_vec(onehot.clone(), (n, 4 * width), device)
                        .map_err(bad)?,
                    class: Tensor::from_vec(class.clone(), (n, 2 * width), device).map_err(bad)?,
                    valid: Tensor::from_vec(valid.clone(), (n, width), device).map_err(bad)?,
                });
            }
            start = end;
        }
        Ok(Engine {
            pool: pool.clone(),
            n,
            columns,
            resident,
        })
    }

    /// Which cards, and what share of the rows each does, in words.
    pub fn describe(&self) -> String {
        self.pool
            .members()
            .iter()
            .map(|m| format!("{} {:.0} %", m.name, 100.0 * m.weight))
            .collect::<Vec<_>>()
            .join(", ")
    }

    /// The counts for every pair, each column taken `weights[column]` times (none: once). The rows of the answer are shared
    /// out among the cards by their weights, and each card goes through all the chunks of columns.
    pub fn counts(&self, weights: Option<&[u32]>) -> Result<Counts, String> {
        if let Some(w) = weights {
            check_weights(w, self.columns)?;
        }
        let n = self.n;
        let shares = self.pool.shares(n);
        let parts = self.pool.run(|k, member| -> Result<[Vec<f32>; 3], String> {
            let rows = shares[k].clone();
            if rows.is_empty() {
                return Ok([Vec::new(), Vec::new(), Vec::new()]);
            }
            let device: &Device = &member.device;
            let _card = member.exclusive();
            let zeros = || Tensor::zeros((rows.len(), n), DType::F32, device).map_err(bad);
            let (mut identical, mut same_class, mut valid_both) = (zeros()?, zeros()?, zeros()?);
            for chunk in &self.resident[k] {
                let w = weights.map(|w| &w[chunk.columns.clone()]);
                let product = |whole: &Tensor, per_column: usize| -> Result<Tensor, String> {
                    let mut left = whole.narrow(0, rows.start, rows.len()).map_err(bad)?;
                    if let Some(w) = w {
                        // each column's weight sits on all the places that the column has
                        let repeated: Vec<f32> = w
                            .iter()
                            .flat_map(|&x| std::iter::repeat_n(x as f32, per_column))
                            .collect();
                        let width = repeated.len();
                        let on_card =
                            Tensor::from_vec(repeated, (1, width), device).map_err(bad)?;
                        left = left.broadcast_mul(&on_card).map_err(bad)?;
                    }
                    left.matmul(&whole.t().map_err(bad)?).map_err(bad)
                };
                identical = (identical + product(&chunk.onehot, 4)?).map_err(bad)?;
                same_class = (same_class + product(&chunk.class, 2)?).map_err(bad)?;
                valid_both = (valid_both + product(&chunk.valid, 1)?).map_err(bad)?;
            }
            let flat = |t: &Tensor| -> Result<Vec<f32>, String> {
                t.flatten_all().map_err(bad)?.to_vec1::<f32>().map_err(bad)
            };
            Ok([flat(&valid_both)?, flat(&identical)?, flat(&same_class)?])
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

/// The columns `range` of every row written as 0 and 1 for the three products: (n by 4 for each column), (n by 2), (n by 1).
fn write_out(rows: &[Vec<u8>], range: std::ops::Range<usize>) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let (n, width) = (rows.len(), range.len());
    let mut onehot = vec![0.0f32; n * 4 * width];
    let mut class = vec![0.0f32; n * 2 * width];
    let mut valid = vec![0.0f32; n * width];
    for (i, row) in rows.iter().enumerate() {
        for (c, &letter) in row[range.clone()].iter().enumerate() {
            if let Some(b) = base(letter) {
                onehot[i * 4 * width + c * 4 + b] = 1.0;
                // A and G (0 and 2) are purines, C and T (1 and 3) pyrimidines
                class[i * 2 * width + c * 2 + b % 2] = 1.0;
                valid[i * width + c] = 1.0;
            }
        }
    }
    (onehot, class, valid)
}

#[cfg(test)]
pub(crate) mod tests;
