//! Progressive multiple sequence alignment in the manner of ClustalW (Thompson, Higgins and Gibson 1994), written from the published
//! method and not from ClustalW's code, and so not ClustalW: it follows the same steps and its answers will differ from ClustalW's.
//!
//! 1. **Distances.** Every pair of sequences is aligned (a global alignment with affine gaps) and the distance is worked out from the
//!    aligned columns: Kimura's two-parameter distance for nucleotides, Kimura's correction of the fraction of different residues
//!    (`-ln(1 - D - 0.2 D^2)`) for proteins. For many sequences, where all the pairs would take too long, the fraction of the 3-mers
//!    (6-mers for nucleotides) that two sequences share gives the distance instead (Edgar 2004), as ClustalW's own "fast" mode does
//!    with words.
//! 2. **Guide tree.** Neighbour-joining (Saitou and Nei 1987) rooted at the middle of its longest path, or UPGMA.
//! 3. **Sequence weights** from the tree (ClustalW's): the weight of a sequence is the sum, along its path to the root, of each
//!    branch's length divided by the number of sequences below it, so that a crowd of near copies counts for about as much as one
//!    sequence that stands alone.
//! 4. **Progressive alignment.** Going up the tree, the two groups below each join are aligned as profiles: a column of one against a
//!    column of the other scores the weighted average of the substitution scores of the residues in them (BLOSUM62 for proteins, +5 and
//!    -4 for nucleotides; gaps count nothing), and gaps cost an opening and an extension, with ClustalW's position-specific changes:
//!    opening a gap where sequences already have gaps costs less, in proportion to how many of them (by weight) do not, and opening one
//!    within eight columns of a gap costs from about two and a quarter to three and three quarters times as much. A gap is never taken out of a group once it is in.
//! 5. **Refinement** (optional, in the manner of MUSCLE): each branch of the tree divides the sequences in two; the two halves are
//!    aligned again as profiles, and the new alignment is kept when its weighted sum-of-pairs score is higher. Repeated until nothing
//!    improves or the passes are used up.
//!
//! Not done, from ClustalW: the choice among several matrices by how divergent the sequences are (BLOSUM62 only), the residue-specific
//! gap penalties and the hydrophobic-stretch rule, and delaying the divergent sequences.

use super::align::MAX_CELLS;
use super::alnstats::is_protein;
use super::msa::Msa;
use super::par::par_map;
use super::tree::{distance_from_counts, neighbor_joining, upgma, DistanceModel, Node};

const GAP: u8 = b'-';
const AMINO: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";
const NUCLEO: &[u8; 4] = b"ACGT";

/// BLOSUM62 (Henikoff and Henikoff 1992), in the order of `AMINO`. Typed from the published matrix; the tests check that it is
/// symmetric and that its diagonal and its best known entries are the published ones.
#[rustfmt::skip]
const BLOSUM62: [[i8; 20]; 20] = [
    [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0],
    [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3],
    [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3],
    [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2],
    [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2],
    [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3],
    [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1],
    [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2],
    [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2],
    [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1],
    [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4],
];

/// What a gap costs, in the units of the matrices: the first position of a gap the first number, each further position the second;
/// and what a mismatch of nucleotides scores. Chosen on sequences that were evolved by a program (the test `tuning_the_gap_costs...`,
/// with seeds that the other tests do not use): of the costs tried these aligned the most pairs as the program had made them. Those
/// sequences have more short insertions and deletions than most real alignments do, so real sequences may do better with dearer gaps.
const DNA_COSTS: (f32, f32, f32) = (10.0, 2.0, -4.0);
const PROTEIN_COSTS: (f32, f32) = (8.0, 0.5);
/// A gap opened within this many columns of a gap that is there costs more (ClustalW's gap separation distance).
const GAP_DISTANCE: usize = 8;
/// More sequences than this are not aligned: the distances are all the pairs.
pub const MAX_SEQUENCES: usize = 2000;
/// Up to this many cells of dynamic programming (the pairs times their lengths squared) the distances come from aligning every pair;
/// beyond it from the 3-mers or 6-mers the sequences share.
const FULL_BUDGET: f64 = 4.0e8;
/// Refinement is left out above this many sequences (it adds up the scores of all pairs).
const REFINE_MOST: usize = 60;
/// The most pairs times columns that the score of an alignment may add up, for refinement to be done.
const REFINE_CELLS: usize = 40_000_000;

/// The most passes of refinement a request may ask for.
pub const MOST_REFINE: usize = 10;

/// Which sort of residues the sequences are.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Dna,
    Protein,
}

/// The tree that tells in which order the sequences are joined.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Guide {
    NeighborJoining,
    Upgma,
}

impl Guide {
    /// The word a message and a kept text give it.
    pub fn name(self) -> &'static str {
        match self {
            Guide::NeighborJoining => "nj",
            Guide::Upgma => "upgma",
        }
    }

    pub fn from_name(name: &str) -> Option<Guide> {
        [Guide::NeighborJoining, Guide::Upgma]
            .into_iter()
            .find(|g| g.name() == name)
    }

    fn label(self) -> &'static str {
        match self {
            Guide::NeighborJoining => "neighbour-joining",
            Guide::Upgma => "UPGMA",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Options {
    pub guide: Guide,
    /// Gaps at the ends of a sequence cost nothing in the progressive alignment (for sequences of different lengths, such as
    /// fragments). The distances that make the guide tree always come from alignments in which they cost, as in ClustalW's pairwise
    /// stage.
    pub free_end_gaps: bool,
    /// Passes of refinement; 0 for none.
    pub refine: usize,
    /// The distances come from the words that two sequences share and not from aligning every pair (ClustalW's "fast" pairwise
    /// alignment): quicker, a little less exact. It is so anyway when all the pairs would take too long.
    pub fast: bool,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            guide: Guide::NeighborJoining,
            free_end_gaps: true,
            refine: 0,
            fast: false,
        }
    }
}

/// The scores of a kind of residue: the matrix, and the place of each byte in it. The last place is every other letter or sign.
struct Scores {
    kind: Kind,
    size: usize,
    matrix: Vec<f32>,
    index: [u8; 256],
    open: f32,
    extend: f32,
}

impl Scores {
    fn new(kind: Kind) -> Scores {
        match kind {
            Kind::Dna => Scores::with(kind, DNA_COSTS.0, DNA_COSTS.1, DNA_COSTS.2),
            Kind::Protein => Scores::with(kind, PROTEIN_COSTS.0, PROTEIN_COSTS.1, 0.0),
        }
    }

    /// The scores with the given gap opening, gap extension and (for nucleotides) mismatch.
    fn with(kind: Kind, open: f32, extend: f32, mismatch: f32) -> Scores {
        let letters: &[u8] = match kind {
            Kind::Dna => NUCLEO,
            Kind::Protein => AMINO,
        };
        let size = letters.len() + 1;
        let mut index = [(size - 1) as u8; 256];
        for (k, &l) in letters.iter().enumerate() {
            index[l as usize] = k as u8;
            index[l.to_ascii_lowercase() as usize] = k as u8;
        }
        if kind == Kind::Dna {
            index[b'U' as usize] = 3;
            index[b'u' as usize] = 3;
        }
        let mut matrix = vec![-1.0f32; size * size];
        for a in 0..size - 1 {
            for b in 0..size - 1 {
                matrix[a * size + b] = match kind {
                    Kind::Dna if a == b => 5.0,
                    Kind::Dna => mismatch,
                    Kind::Protein => f32::from(BLOSUM62[a][b]),
                };
            }
        }
        Scores {
            kind,
            size,
            matrix,
            index,
            open,
            extend,
        }
    }

    fn score(&self, a: usize, b: usize) -> f32 {
        self.matrix[a * self.size + b]
    }

    fn place(&self, byte: u8) -> usize {
        usize::from(self.index[byte as usize])
    }
}

/// Aligned sequences of one group: which sequences (by their place in the input) and their rows, of one length.
#[derive(Clone)]
struct Group {
    ids: Vec<usize>,
    rows: Vec<Vec<u8>>,
}

impl Group {
    fn single(id: usize, row: Vec<u8>) -> Group {
        Group {
            ids: vec![id],
            rows: vec![row],
        }
    }

    fn columns(&self) -> usize {
        self.rows.first().map_or(0, Vec::len)
    }
}

/// A group as columns of weighted residue counts, and what it costs to open a gap beside each column.
struct Profile {
    freq: Vec<f32>,
    total: f32,
    opening: Vec<f32>,
}

fn profile(group: &Group, weights: &[f32], sc: &Scores) -> Profile {
    let (len, size) = (group.columns(), sc.size);
    let mut freq = vec![0.0f32; len * size];
    let mut gaps = vec![0.0f32; len];
    let mut total = 0.0f32;
    for (&id, row) in group.ids.iter().zip(&group.rows) {
        let w = weights[id];
        total += w;
        for (c, &byte) in row.iter().enumerate() {
            if byte == GAP {
                gaps[c] += w;
            } else {
                freq[c * size + sc.place(byte)] += w;
            }
        }
    }
    // how far each column is from the nearest column that has a gap (0 for one that has)
    let mut near = vec![usize::MAX; len];
    let mut last: Option<usize> = None;
    for c in 0..len {
        if gaps[c] > 0.0 {
            last = Some(c);
        }
        near[c] = last.map_or(usize::MAX, |l| c - l);
    }
    last = None;
    for c in (0..len).rev() {
        if gaps[c] > 0.0 {
            last = Some(c);
        }
        if let Some(l) = last {
            near[c] = near[c].min(l - c);
        }
    }
    let opening = (0..len)
        .map(|c| {
            if gaps[c] > 0.0 {
                // where sequences already have gaps a new one is cheaper, by the share that have none
                (1.0 - gaps[c] / total).max(0.05)
            } else if near[c] < GAP_DISTANCE {
                // and near a gap it is dearer
                2.0 + (GAP_DISTANCE - near[c]) as f32 * 2.0 / GAP_DISTANCE as f32
            } else {
                1.0
            }
        })
        .collect();
    Profile {
        freq,
        total,
        opening,
    }
}

const NEG: f32 = f32::NEG_INFINITY;

/// The best of three scores and which it was (0, 1 or 2), the first on a tie.
fn best3(a: f32, b: f32, c: f32) -> (f32, u8) {
    let (v, t) = if a >= b { (a, 0) } else { (b, 1) };
    if v >= c {
        (v, t)
    } else {
        (c, 2)
    }
}

/// What a gap costs, as (opening, extension), at the `at`th place of a profile of `opening.len()` columns (`at` columns lie before
/// it): nothing at either end when the ends are free; elsewhere the opening cost of the scores times the factor of the column the
/// gap lies before (of the last column, for a gap after it).
fn gap_cost(sc: &Scores, opening: &[f32], at: usize, free_ends: bool) -> (f32, f32) {
    let columns = opening.len();
    if free_ends && (at == 0 || at == columns) {
        (0.0, 0.0)
    } else {
        (sc.open * opening[at.min(columns - 1)], sc.extend)
    }
}

/// Align two groups as profiles (affine gaps, three states: a column of each together, a column of `p` against a gap, a column of
/// `q` against a gap), and join their rows. With `free_ends`, gaps before and after everything cost nothing.
fn merge(
    p: &Group,
    q: &Group,
    weights: &[f32],
    sc: &Scores,
    free_ends: bool,
) -> Result<Group, String> {
    let (m, n) = (p.columns(), q.columns());
    if (m + 1).saturating_mul(n + 1) > MAX_CELLS {
        return Err(format!(
            "{m} x {n} columns is too large to align in one piece (limit about 6,000 x 6,000)."
        ));
    }
    let (pp, qq) = (profile(p, weights, sc), profile(q, weights, sc));
    let size = sc.size;
    // what each column of p would score against each residue, and the residues each column of q holds
    let scale = 1.0 / (pp.total * qq.total);
    let mut against = vec![0.0f32; m * size];
    for i in 0..m {
        for a in 0..size {
            let f = pp.freq[i * size + a];
            if f != 0.0 {
                for b in 0..size {
                    against[i * size + b] += f * sc.score(a, b) * scale;
                }
            }
        }
    }
    let held: Vec<Vec<(usize, f32)>> = (0..n)
        .map(|j| {
            (0..size)
                .filter_map(|b| {
                    let f = qq.freq[j * size + b];
                    (f > 0.0).then_some((b, f))
                })
                .collect()
        })
        .collect();
    // a gap in q opposite a column of p, with j columns of q behind it; a gap in p, with i columns of p behind it
    let gap_in_q = |j: usize| gap_cost(sc, &qq.opening, j, free_ends);
    let gap_in_p = |i: usize| gap_cost(sc, &pp.opening, i, free_ends);
    let w = n + 1;
    // for each cell, where each of its three states came from: two bits each (0 the diagonal state, 1 gap in q, 2 gap in p)
    let mut trace = vec![0u8; (m + 1) * w];
    let (mut pm, mut px, mut py) = (vec![NEG; w], vec![NEG; w], vec![NEG; w]);
    let (mut cm, mut cx, mut cy) = (vec![NEG; w], vec![NEG; w], vec![NEG; w]);
    pm[0] = 0.0;
    for j in 1..=n {
        let (open, ext) = gap_in_p(0);
        let (v, t) = if j == 1 {
            (pm[0] - open, 0)
        } else {
            (py[j - 1] - ext, 2)
        };
        py[j] = v;
        trace[j] = t << 4;
    }
    for i in 1..=m {
        cm[0] = NEG;
        cy[0] = NEG;
        let (open, ext) = gap_in_q(0);
        let (v, t) = if i == 1 {
            (pm[0] - open, 0)
        } else {
            (px[0] - ext, 1)
        };
        cx[0] = v;
        trace[i * w] = t << 2;
        for j in 1..=n {
            let s: f32 = held[j - 1]
                .iter()
                .map(|&(b, f)| against[(i - 1) * size + b] * f)
                .sum();
            let (from, tm) = best3(pm[j - 1], px[j - 1], py[j - 1]);
            cm[j] = from + s;
            let (open_x, ext_x) = gap_in_q(j);
            let (vx, tx) = best3(pm[j] - open_x, px[j] - ext_x, py[j] - open_x);
            cx[j] = vx;
            let (open_y, ext_y) = gap_in_p(i);
            let (vy, ty) = best3(cm[j - 1] - open_y, cx[j - 1] - open_y, cy[j - 1] - ext_y);
            // best3 of (M, X, Y) answers 0, 1 or 2 in that order, which is the order of the arguments above
            cy[j] = vy;
            trace[i * w + j] = tm | (tx << 2) | (ty << 4);
        }
        std::mem::swap(&mut pm, &mut cm);
        std::mem::swap(&mut px, &mut cx);
        std::mem::swap(&mut py, &mut cy);
    }
    let (_, mut state) = best3(pm[n], px[n], py[n]);
    // walk back: 0 a column of each, 1 a column of p alone, 2 a column of q alone
    let mut ops = Vec::with_capacity(m + n);
    let (mut i, mut j) = (m, n);
    while i > 0 || j > 0 {
        let t = trace[i * w + j];
        ops.push(state);
        match state {
            0 => {
                state = t & 3;
                i -= 1;
                j -= 1;
            }
            1 => {
                state = (t >> 2) & 3;
                i -= 1;
            }
            _ => {
                state = (t >> 4) & 3;
                j -= 1;
            }
        }
    }
    ops.reverse();
    let mut rows: Vec<Vec<u8>> = vec![Vec::with_capacity(ops.len()); p.rows.len() + q.rows.len()];
    let (mut ip, mut iq) = (0, 0);
    for &op in &ops {
        let (take_p, take_q) = (op != 2, op != 1);
        for (k, row) in p.rows.iter().enumerate() {
            rows[k].push(if take_p { row[ip] } else { GAP });
        }
        for (k, row) in q.rows.iter().enumerate() {
            rows[p.rows.len() + k].push(if take_q { row[iq] } else { GAP });
        }
        ip += usize::from(take_p);
        iq += usize::from(take_q);
    }
    Ok(Group {
        ids: p.ids.iter().chain(&q.ids).copied().collect(),
        rows,
    })
}

/// What a distance is when two sequences differ too much for the model to say how much: nucleotides' Kimura distance is capped here,
/// and a protein's when the share of different residues is too large to take the logarithm of.
const DNA_SATURATED: f64 = 2.0;
const PROTEIN_SATURATED: f64 = 5.0;

/// A nucleotide as the distance reads it: in capitals, with U as T.
fn as_dna(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    }
}

/// The distance between two sequences from their global alignment, in which the end gaps cost: Kimura's for nucleotides (U read as
/// T), Kimura's correction of the share of different residues for proteins. The residues of the shorter sequence that the alignment
/// leaves opposite gaps count as differences (the largest there is, for nucleotides), so that two sequences that have only a letter
/// or two lined up, because a gap costs less than the letters that would be set against each other, are far apart, not alike.
fn pair_distance(a: &[u8], b: &[u8], sc: &Scores) -> Result<f64, String> {
    let unit = [1.0f32, 1.0];
    let both = merge(
        &Group::single(0, a.to_vec()),
        &Group::single(1, b.to_vec()),
        &unit,
        sc,
        false,
    )?;
    let (x, y) = (&both.rows[0], &both.rows[1]);
    let shorter = a.len().min(b.len());
    let (mut compared, mut same) = (0usize, 0usize);
    let (mut sites, mut transitions, mut transversions) = (0u32, 0u32, 0u32);
    for (&u, &v) in x.iter().zip(y) {
        if u == GAP || v == GAP {
            continue;
        }
        compared += 1;
        same += usize::from(u.eq_ignore_ascii_case(&v));
        let (u, v) = (as_dna(u), as_dna(v));
        if b"ACGT".contains(&u) && b"ACGT".contains(&v) {
            sites += 1;
            if u != v {
                if matches!(u, b'A' | b'G') == matches!(v, b'A' | b'G') {
                    transitions += 1;
                } else {
                    transversions += 1;
                }
            }
        }
    }
    let unmatched = shorter.saturating_sub(compared);
    Ok(match sc.kind {
        Kind::Dna => {
            let unmatched = u32::try_from(unmatched).unwrap_or(u32::MAX);
            distance_from_counts(
                sites.saturating_add(unmatched),
                transitions,
                transversions.saturating_add(unmatched),
                DistanceModel::Kimura2,
            )
            .unwrap_or(DNA_SATURATED)
        }
        Kind::Protein => {
            let d = 1.0 - same as f64 / shorter.max(1) as f64;
            let inside = 1.0 - d - 0.2 * d * d;
            if inside < 1e-4 {
                PROTEIN_SATURATED
            } else {
                -inside.ln()
            }
        }
    })
}

/// The words of `k` residues of a sequence, as numbers, in order.
fn words(seq: &[u8], k: usize, sc: &Scores) -> Vec<u32> {
    let base = sc.size as u32;
    let mut out: Vec<u32> = seq
        .windows(k)
        .map(|w| w.iter().fold(0u32, |c, &b| c * base + sc.place(b) as u32))
        .collect();
    out.sort_unstable();
    out
}

/// One minus the share of the words of the shorter sequence that the two have in common, each word counted as often as both have
/// it (Edgar 2004).
fn word_distance(a: &[u32], b: &[u32], k: usize, lengths: (usize, usize)) -> f64 {
    let shorter = lengths.0.min(lengths.1);
    if shorter < k {
        return 1.0;
    }
    let (mut i, mut j, mut common) = (0, 0, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                common += 1;
                i += 1;
                j += 1;
            }
        }
    }
    (1.0 - common as f64 / (shorter - k + 1) as f64).clamp(0.0, 1.0)
}

/// The distances between all the sequences, and whether they came from the words the sequences share.
fn distances(seqs: &[Vec<u8>], sc: &Scores, fast: bool) -> Result<(Vec<Vec<f64>>, bool), String> {
    let n = seqs.len();
    let mean = seqs.iter().map(Vec::len).sum::<usize>() as f64 / n as f64;
    let pairs = (n * (n - 1) / 2) as f64;
    let by_words = fast || pairs * mean * mean > FULL_BUDGET;
    let mut d = vec![vec![0.0f64; n]; n];
    if by_words {
        let k = if sc.kind == Kind::Dna { 6 } else { 3 };
        let set: Vec<Vec<u32>> = seqs.iter().map(|s| words(s, k, sc)).collect();
        for i in 0..n {
            for j in i + 1..n {
                let x = word_distance(&set[i], &set[j], k, (seqs[i].len(), seqs[j].len()));
                (d[i][j], d[j][i]) = (x, x);
            }
        }
        return Ok((d, true));
    }
    let list: Vec<(usize, usize)> = (0..n)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .collect();
    let found = par_map(&list, |&(i, j)| pair_distance(&seqs[i], &seqs[j], sc));
    for (&(i, j), x) in list.iter().zip(found) {
        let x = x?;
        (d[i][j], d[j][i]) = (x, x);
    }
    Ok((d, false))
}

/// A tree with the sequences at its leaves and every join below the root.
struct Tree {
    kids: Vec<Vec<usize>>,
    leaf: Vec<Option<usize>>,
    len: Vec<f64>,
    root: usize,
}

impl Tree {
    fn push(&mut self, leaf: Option<usize>, len: f64) -> usize {
        self.kids.push(Vec::new());
        self.leaf.push(leaf);
        self.len.push(len);
        self.kids.len() - 1
    }

    /// The nodes with every node after the node it hangs from.
    fn preorder(&self) -> Vec<usize> {
        let mut order = Vec::with_capacity(self.kids.len());
        let mut stack = vec![self.root];
        while let Some(v) = stack.pop() {
            order.push(v);
            stack.extend(self.kids[v].iter().rev());
        }
        order
    }

    /// How many sequences hang below each node (a leaf has one).
    fn leaves_below(&self) -> Vec<usize> {
        let mut count = vec![0usize; self.kids.len()];
        for &v in self.preorder().iter().rev() {
            count[v] = if self.leaf[v].is_some() {
                1
            } else {
                self.kids[v].iter().map(|&k| count[k]).sum()
            };
        }
        count
    }

    /// The sequences below `v`.
    fn sequences_below(&self, v: usize) -> Vec<usize> {
        let mut out = Vec::new();
        let mut stack = vec![v];
        while let Some(u) = stack.pop() {
            match self.leaf[u] {
                Some(s) => out.push(s),
                None => stack.extend(&self.kids[u]),
            }
        }
        out
    }
}

/// A rooted tree from a tree of `Node`s whose leaves are named by the number of the sequence.
fn rooted(node: &Node, tree: &mut Tree) -> usize {
    let leaf = node
        .is_leaf()
        .then(|| node.name.parse::<usize>().unwrap_or(0));
    let v = tree.push(leaf, node.length);
    let kids: Vec<usize> = node.children.iter().map(|k| rooted(k, tree)).collect();
    tree.kids[v] = kids;
    v
}

/// The tree without a root, as edges.
struct Graph {
    edges: Vec<Vec<(usize, f64)>>,
    leaf: Vec<Option<usize>>,
}

fn unrooted(node: &Node, parent: Option<(usize, f64)>, g: &mut Graph) -> usize {
    let v = g.edges.len();
    g.edges.push(Vec::new());
    g.leaf.push(
        node.is_leaf()
            .then(|| node.name.parse::<usize>().unwrap_or(0)),
    );
    if let Some((p, len)) = parent {
        g.edges[v].push((p, len));
        g.edges[p].push((v, len));
    }
    for kid in &node.children {
        unrooted(kid, Some((v, kid.length)), g);
    }
    v
}

/// How far each node is from `from`, and the node before it on the way back.
fn walk(g: &Graph, from: usize) -> (Vec<f64>, Vec<usize>) {
    let mut dist = vec![f64::NAN; g.edges.len()];
    let mut back = vec![usize::MAX; g.edges.len()];
    dist[from] = 0.0;
    let mut stack = vec![from];
    while let Some(v) = stack.pop() {
        for &(u, len) in &g.edges[v] {
            if dist[u].is_nan() {
                dist[u] = dist[v] + len;
                back[u] = v;
                stack.push(u);
            }
        }
    }
    (dist, back)
}

/// The leaf that is farthest from the distances `dist` start at (the first on a tie).
fn farthest(g: &Graph, dist: &[f64]) -> usize {
    let mut best = usize::MAX;
    for (v, leaf) in g.leaf.iter().enumerate() {
        if leaf.is_some() && (best == usize::MAX || dist[v] > dist[best]) {
            best = v;
        }
    }
    best
}

/// Hang the part of the tree beyond `v` (away from `from`) below a new node of `tree`.
fn hang(g: &Graph, v: usize, from: usize, len: f64, tree: &mut Tree) -> usize {
    let me = tree.push(g.leaf[v], len);
    let kids: Vec<usize> = g.edges[v]
        .iter()
        .filter(|&&(u, _)| u != from)
        .map(|&(u, l)| hang(g, u, v, l, tree))
        .collect();
    tree.kids[me] = kids;
    me
}

/// The tree rooted halfway along the longest path between two leaves.
fn midpoint_rooted(g: &Graph) -> Tree {
    let first = g.leaf.iter().position(Option::is_some).unwrap_or(0);
    let (from_first, _) = walk(g, first);
    let a = farthest(g, &from_first);
    let (from_a, back) = walk(g, a);
    let mut b = farthest(g, &from_a);
    if b == a {
        // every leaf is where the first is: any other leaf will do for the other end
        b = (0..g.leaf.len())
            .rev()
            .find(|&v| g.leaf[v].is_some() && v != a)
            .unwrap_or(b);
    }
    let half = from_a[b] / 2.0;
    // on the way from b back to a: the first node that is no farther from a than the middle
    let mut u = b;
    while from_a[back[u]] > half && back[u] != a {
        u = back[u];
    }
    let v = back[u];
    let edge = from_a[u] - from_a[v];
    let up = (half - from_a[v]).clamp(0.0, edge);
    let mut tree = Tree {
        kids: Vec::new(),
        leaf: Vec::new(),
        len: Vec::new(),
        root: 0,
    };
    let root = tree.push(None, 0.0);
    let below = hang(g, u, v, edge - up, &mut tree);
    let above = hang(g, v, u, up, &mut tree);
    tree.kids[root] = vec![below, above];
    tree
}

/// The guide tree of the sequences from their distances.
fn guide_tree(guide: Guide, dist: &[Vec<f64>]) -> Result<Tree, String> {
    let n = dist.len();
    let names: Vec<String> = (0..n).map(|i| i.to_string()).collect();
    let tree = Tree {
        kids: Vec::new(),
        leaf: Vec::new(),
        len: Vec::new(),
        root: 0,
    };
    match guide {
        Guide::Upgma => {
            let node = upgma(&names, dist).ok_or("The guide tree could not be made.")?;
            let mut tree = tree;
            tree.root = rooted(&node, &mut tree);
            Ok(tree)
        }
        Guide::NeighborJoining => {
            let node = neighbor_joining(&names, dist).ok_or("The guide tree could not be made.")?;
            let mut g = Graph {
                edges: Vec::new(),
                leaf: Vec::new(),
            };
            unrooted(&node, None, &mut g);
            Ok(midpoint_rooted(&g))
        }
    }
}

/// The weight of each sequence, from the tree (ClustalW's): along the path from the sequence to the root, the length of each branch
/// shared out among the sequences below it. They add up to the number of sequences (a little more when one would weigh next to
/// nothing, which is raised to a floor of a thousandth); all one when the tree has no length.
fn weights_from(tree: &Tree, n: usize) -> Vec<f32> {
    let below = tree.leaves_below();
    let mut weight = vec![0.0f64; tree.kids.len()];
    for &v in &tree.preorder() {
        let carried = if v == tree.root { 0.0 } else { weight[v] };
        for &k in &tree.kids[v] {
            weight[k] = carried + tree.len[k] / below[k] as f64;
        }
    }
    let mut w = vec![0.0f64; n];
    for (v, leaf) in tree.leaf.iter().enumerate() {
        if let Some(s) = *leaf {
            w[s] = weight[v];
        }
    }
    let total: f64 = w.iter().sum();
    if total <= 0.0 || !total.is_finite() {
        return vec![1.0; n];
    }
    // no sequence counts for nothing: a group of such would have no weight to score with
    w.iter()
        .map(|x| (x / total * n as f64).max(1e-3) as f32)
        .collect()
}

/// Rows without the columns that are all gaps.
fn without_empty_columns(rows: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let len = rows.first().map_or(0, Vec::len);
    let keep: Vec<usize> = (0..len)
        .filter(|&c| rows.iter().any(|r| r[c] != GAP))
        .collect();
    rows.iter()
        .map(|r| keep.iter().map(|&c| r[c]).collect())
        .collect()
}

/// The weighted sum of the scores of all the pairs of rows: residues against residues by the matrix, each run of gaps in one against
/// residues of the other an opening and an extension for each position after the first (nothing at the ends when `free_ends`),
/// columns that are gaps in both left out.
fn sum_of_pairs(
    rows: &[Vec<u8>],
    by: &[usize],
    weights: &[f32],
    sc: &Scores,
    free_ends: bool,
) -> f64 {
    let n = rows.len();
    let len = rows.first().map_or(0, Vec::len);
    let span: Vec<(usize, usize)> = rows
        .iter()
        .map(|r| {
            let first = r.iter().position(|&b| b != GAP).unwrap_or(len);
            let last = r.iter().rposition(|&b| b != GAP).map_or(0, |p| p + 1);
            (first, last)
        })
        .collect();
    let mut total = 0.0f64;
    for i in 0..n {
        for j in i + 1..n {
            let mut score = 0.0f64;
            // which row the gap that runs is in (0 none, 1 row i, 2 row j)
            let mut run = 0u8;
            for (c, (&x, &y)) in rows[i].iter().zip(&rows[j]).enumerate() {
                if x == GAP && y == GAP {
                    continue;
                }
                if x != GAP && y != GAP {
                    score += f64::from(sc.score(sc.place(x), sc.place(y)));
                    run = 0;
                    continue;
                }
                let (who, own) = if x == GAP { (1u8, i) } else { (2u8, j) };
                let at_end = c < span[own].0 || c >= span[own].1;
                if !(free_ends && at_end) {
                    score -= f64::from(if run == who { sc.extend } else { sc.open });
                }
                run = who;
            }
            total += f64::from(weights[by[i]]) * f64::from(weights[by[j]]) * score;
        }
    }
    total
}

/// Put the rows of a group in the order of the input.
fn in_input_order(group: Group) -> Group {
    let mut order: Vec<usize> = (0..group.ids.len()).collect();
    order.sort_by_key(|&k| group.ids[k]);
    Group {
        ids: order.iter().map(|&k| group.ids[k]).collect(),
        rows: order.iter().map(|&k| group.rows[k].clone()).collect(),
    }
}

/// Refine: for each branch of the tree, align its two sides again as profiles, and keep the result when the sum-of-pairs score is
/// higher. At most `passes` passes over the branches, fewer when a pass changes nothing.
fn refine(
    mut group: Group,
    tree: &Tree,
    weights: &[f32],
    sc: &Scores,
    free_ends: bool,
    passes: usize,
) -> Result<Group, String> {
    let n = group.ids.len();
    let cells = n * n / 2 * group.columns();
    if !(3..=REFINE_MOST).contains(&n) || cells > REFINE_CELLS {
        return Ok(group);
    }
    let score = |g: &Group| sum_of_pairs(&g.rows, &g.ids, weights, sc, free_ends);
    let mut best = score(&group);
    for _ in 0..passes {
        let mut improved = false;
        // the two branches under the root cut the sequences in the same two: one of them is enough
        let twin = tree.kids[tree.root].get(1).copied();
        for v in tree.preorder() {
            if v == tree.root || Some(v) == twin {
                continue;
            }
            let side: std::collections::HashSet<usize> =
                tree.sequences_below(v).into_iter().collect();
            if side.is_empty() || side.len() == n {
                continue;
            }
            let split = |inside: bool| -> Group {
                let keep: Vec<usize> = (0..n)
                    .filter(|&k| side.contains(&group.ids[k]) == inside)
                    .collect();
                let rows: Vec<Vec<u8>> = keep.iter().map(|&k| group.rows[k].clone()).collect();
                Group {
                    ids: keep.iter().map(|&k| group.ids[k]).collect(),
                    rows: without_empty_columns(&rows),
                }
            };
            let (one, two) = (split(true), split(false));
            // a halving that is too large to align in one piece is left as it is: refining is an extra, and the alignment is whole
            let Ok(merged) = merge(&one, &two, weights, sc, free_ends) else {
                continue;
            };
            let again = in_input_order(merged);
            let now = score(&again);
            if now > best + 1e-6 {
                (group, best, improved) = (again, now, true);
            }
        }
        if !improved {
            break;
        }
    }
    Ok(group)
}

/// What an alignment says of itself: how the distances were found, and the settings that were used.
#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    pub msa: Msa,
    /// The distances came from shared words: asked for, or because aligning every pair would take too long.
    pub by_words: bool,
    /// Refinement was asked for and left out (too many sequences or columns).
    pub refine_skipped: bool,
    /// The weighted sum-of-pairs score of the alignment, and of the alignment before refinement; none when there are too many
    /// sequences or columns to add up the pairs.
    pub score: Option<f64>,
    pub unrefined_score: Option<f64>,
}

/// What the engine line of a progressive alignment begins with (so that a result can be told from another method's).
pub const NAME: &str = "progressive alignment after ClustalW";

/// What an alignment says of how it was made, in a line (the name of the engine in a log and in the record of a result).
pub fn describe(o: &Options, r: &Report) -> String {
    let mut parts = vec![format!("{} guide tree", o.guide.label())];
    parts.push(if r.by_words {
        "distances from shared words".to_string()
    } else {
        "distances from aligned pairs".to_string()
    });
    parts.push(if o.free_end_gaps {
        "free end gaps".to_string()
    } else {
        "end gaps cost".to_string()
    });
    match (o.refine, r.refine_skipped) {
        (0, _) => {}
        (_, true) => parts.push(format!(
            "refinement left out (it needs 3 to {REFINE_MOST} sequences and not too many columns)"
        )),
        (n, false) => {
            let n = n.min(MOST_REFINE);
            parts.push(format!(
                "refined, up to {n} pass{}",
                if n == 1 { "" } else { "es" }
            ));
        }
    }
    format!("{NAME} (built in; {})", parts.join(", "))
}

/// Align `seqs` (name, residues; a gap sign in a sequence is ignored). `Err` when there is nothing to align, there are too many
/// sequences, or two groups are too long to align in one piece.
pub fn align(seqs: &[(String, Vec<u8>)], o: &Options) -> Result<Report, String> {
    align_scored(seqs, o, Scores::new)
}

/// `align`, with the scores of each kind of residue made by `scores`.
fn align_scored(
    seqs: &[(String, Vec<u8>)],
    o: &Options,
    scores: impl Fn(Kind) -> Scores,
) -> Result<Report, String> {
    if seqs.is_empty() {
        return Err("There is nothing to align.".into());
    }
    if seqs.len() > MAX_SEQUENCES {
        return Err(format!(
            "At most {MAX_SEQUENCES} sequences can be aligned this way; use the centre-star method or MAFFT for more."
        ));
    }
    let clean: Vec<Vec<u8>> = seqs
        .iter()
        .map(|(_, s)| {
            s.iter()
                .copied()
                .filter(|&b| b != GAP && b != b'.')
                .collect()
        })
        .collect();
    if clean.iter().any(Vec::is_empty) {
        return Err("Every sequence needs at least one residue.".into());
    }
    let names: Vec<String> = seqs.iter().map(|(n, _)| n.clone()).collect();
    let n = clean.len();
    let kind = if is_protein(&clean) {
        Kind::Protein
    } else {
        Kind::Dna
    };
    let sc = scores(kind);
    if n == 1 {
        return Ok(Report {
            msa: Msa { names, rows: clean },
            by_words: false,
            refine_skipped: o.refine > 0,
            score: Some(0.0),
            unrefined_score: Some(0.0),
        });
    }
    let (dist, by_words) = distances(&clean, &sc, o.fast)?;
    let (tree, weights) = if n == 2 {
        // two sequences have no tree to speak of
        let mut t = Tree {
            kids: Vec::new(),
            leaf: Vec::new(),
            len: Vec::new(),
            root: 0,
        };
        let root = t.push(None, 0.0);
        let a = t.push(Some(0), dist[0][1] / 2.0);
        let b = t.push(Some(1), dist[0][1] / 2.0);
        t.kids[root] = vec![a, b];
        let w = weights_from(&t, 2);
        (t, w)
    } else {
        let t = guide_tree(o.guide, &dist)?;
        let w = weights_from(&t, n);
        (t, w)
    };
    // up the tree: every leaf is a group of one, every join a merge of the groups below it
    let mut groups: Vec<Option<Group>> = vec![None; tree.kids.len()];
    for &v in tree.preorder().iter().rev() {
        let made = match tree.leaf[v] {
            Some(s) => Group::single(s, clean[s].clone()),
            None => {
                let mut kids = tree.kids[v].iter().filter_map(|&k| groups[k].take());
                let first = kids
                    .next()
                    .ok_or("The guide tree has a join with nothing below it.")?;
                kids.try_fold(first, |joined, next| {
                    merge(&joined, &next, &weights, &sc, o.free_end_gaps)
                })?
            }
        };
        groups[v] = Some(made);
    }
    let mut group = in_input_order(groups[tree.root].take().ok_or("Nothing was aligned.")?);
    // the pairs can be added up when there are not too many of them
    let countable = n <= REFINE_MOST && n * n / 2 * group.columns() <= REFINE_CELLS;
    let score_of = |g: &Group| {
        countable.then(|| sum_of_pairs(&g.rows, &g.ids, &weights, &sc, o.free_end_gaps))
    };
    let unrefined_score = score_of(&group);
    let mut refine_skipped = false;
    if o.refine > 0 {
        refine_skipped = n < 3 || !countable;
        group = in_input_order(refine(
            group,
            &tree,
            &weights,
            &sc,
            o.free_end_gaps,
            o.refine.min(MOST_REFINE),
        )?);
    }
    let score = score_of(&group);
    Ok(Report {
        msa: Msa {
            names,
            rows: without_empty_columns(&group.rows),
        },
        by_words,
        refine_skipped,
        score,
        unrefined_score,
    })
}

pub mod aln;

#[cfg(test)]
mod tests;
