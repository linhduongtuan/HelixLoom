//! The assembly: for each k in turn, count the k-mers of the reads, make the graph, clean it, and keep its segments for the next k.
//!
//! The iteration is that of the paper in outline: the graph for the next, larger, k is made from the reads and from the segments of the
//! graph before (the segments are put in as reads that were seen as often as the segment was), so that a place where the reads are too few
//! for the larger k keeps what the smaller one made of it, and a repeat that is shorter than the larger k is resolved by it. Errors in the
//! reads come back at each k and are cleaned again.

use super::clean::{auto_cutoff, pop_bulges, remove_below, remove_isolated, remove_tips, Settings};
use super::graph::Graph;
use super::pack::{count, Pack, Table};
use super::unitig::unitigs;

/// The k that are used when none is asked for, as many of them as the reads are long enough for.
pub const DEFAULT_KS: [usize; 3] = [21, 33, 55];
/// The most rounds of tips and bulges at one k (each round can make new ones of what it leaves).
const MAX_ROUNDS: usize = 20;
/// The longest a bulge may be, when not given, in k-mers: a single changed base makes two paths of k k-mers, and k is at most 63.
const BULGE_LENGTH: usize = 150;

#[derive(Clone, Debug, PartialEq)]
pub enum Cutoff {
    Off,
    Auto,
    Value(f64),
}

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    /// The k, odd and rising; none: chosen from the length of the reads.
    pub ks: Vec<usize>,
    /// A k-mer that is seen fewer times is dropped before the graph is made (1: none is).
    pub min_count: u32,
    /// The longest tip, in k-mers; none: the length of the longest read.
    pub tip_length: Option<usize>,
    pub tip_ratio: f64,
    /// The longest bulge, in k-mers; none: 150.
    pub bulge_length: Option<usize>,
    pub bulge_diff: f64,
    pub cutoff: Cutoff,
    /// At the last k, a segment that is joined to nothing and is shorter than this many bases is removed (0: none is).
    pub isolated_length: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            ks: Vec::new(),
            min_count: 1,
            tip_length: None,
            tip_ratio: 2.5,
            bulge_length: None,
            bulge_diff: 0.1,
            cutoff: Cutoff::Auto,
            isolated_length: 200,
        }
    }
}

/// What one k did.
#[derive(Clone, Debug, PartialEq)]
pub struct Step {
    pub k: usize,
    /// Distinct k-mers (canonical), and segments before any cleaning.
    pub kmers: usize,
    pub segments_made: usize,
    pub tips: usize,
    pub bulges: usize,
    /// The cutoff of the depth filter (0 when it was off) and what it took away.
    pub cutoff: f64,
    pub low_coverage: usize,
    /// Short segments joined to nothing that were removed (at the last k only).
    pub isolated: usize,
    pub segments: usize,
    pub bases: usize,
}

#[derive(Clone, Debug)]
pub struct Assembly {
    pub graph: Graph,
    pub steps: Vec<Step>,
    /// The pieces of reads that were used, and their bases.
    pub pieces: usize,
    pub bases: usize,
}

/// The k to use for reads of at most `longest` bases: those of `DEFAULT_KS` that fit in a read.
pub fn default_ks(longest: usize) -> Vec<usize> {
    DEFAULT_KS.into_iter().filter(|&k| k <= longest).collect()
}

fn check_ks(ks: &[usize]) -> Result<(), String> {
    for &k in ks {
        Pack::new(k)?;
    }
    if ks.windows(2).any(|w| w[0] >= w[1]) {
        return Err("The k must be given in rising order, none twice.".into());
    }
    Ok(())
}

/// Tips and bulges, round after round, until a round finds none; the graph is made compact after each.
fn simplify(g: &mut Graph, s: &Settings, step: &mut Step) {
    for _ in 0..MAX_ROUNDS {
        let tips = remove_tips(g, s);
        let bulges = pop_bulges(g, s);
        step.tips += tips;
        step.bulges += bulges;
        if tips + bulges == 0 {
            break;
        }
        *g = g.compact();
    }
}

/// The segments of the graph before, put in with the coverage they had, where the k-mers of the reads have no count.
fn carry(counts: &mut Table<f32>, before: &Graph, pack: Pack) {
    for id in before.ids() {
        let coverage = before.coverage(id) as f32;
        pack.for_each(&before.segment(id).seq, |f, r| {
            counts.entry(f.min(r)).or_insert(coverage);
        });
    }
}

pub fn assemble(pieces: &[Vec<u8>], o: &Options) -> Result<Assembly, String> {
    let longest = pieces.iter().map(Vec::len).max().unwrap_or(0);
    let ks = if o.ks.is_empty() {
        default_ks(longest)
    } else {
        o.ks.clone()
    };
    if ks.is_empty() {
        return Err(format!(
            "The reads are too short for k = {} (the longest piece is {longest} bases): give a smaller -k, or longer reads.",
            DEFAULT_KS[0]
        ));
    }
    check_ks(&ks)?;
    let mut before: Option<Graph> = None;
    let mut steps = Vec::new();
    for &k in &ks {
        let pack = Pack::new(k)?;
        let mut counts = count(pieces, pack);
        if o.min_count > 1 {
            counts.retain(|_, n| *n >= o.min_count as f32);
        }
        if let Some(g) = &before {
            carry(&mut counts, g, pack);
        }
        if counts.is_empty() {
            return Err(format!(
                "There is no k-mer of {k} bases in the reads (are they shorter than {k}, or all below the quality limit?)."
            ));
        }
        let mut g = unitigs(&counts, pack);
        let mut step = Step {
            k,
            kmers: counts.len(),
            segments_made: g.count(),
            tips: 0,
            bulges: 0,
            cutoff: 0.0,
            low_coverage: 0,
            isolated: 0,
            segments: 0,
            bases: 0,
        };
        drop(counts);
        let settings = Settings {
            tip_length: o.tip_length.unwrap_or(longest),
            tip_ratio: o.tip_ratio,
            bulge_length: o.bulge_length.unwrap_or(BULGE_LENGTH),
            bulge_diff: o.bulge_diff,
        };
        simplify(&mut g, &settings, &mut step);
        step.cutoff = match o.cutoff {
            Cutoff::Off => 0.0,
            Cutoff::Auto => auto_cutoff(&g),
            Cutoff::Value(v) => v,
        };
        step.low_coverage = remove_below(&mut g, step.cutoff);
        if step.low_coverage > 0 {
            g = g.compact();
            simplify(&mut g, &settings, &mut step);
        }
        if k == ks[ks.len() - 1] {
            step.isolated = remove_isolated(&mut g, o.isolated_length);
        }
        step.segments = g.count();
        step.bases = g.total_length();
        steps.push(step);
        before = Some(g);
    }
    let graph = before.ok_or("There is nothing to assemble.")?;
    Ok(Assembly {
        graph,
        steps,
        pieces: pieces.len(),
        bases: pieces.iter().map(Vec::len).sum(),
    })
}

#[cfg(test)]
mod tests;
