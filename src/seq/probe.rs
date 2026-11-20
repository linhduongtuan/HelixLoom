//! TaqMan (hydrolysis) probe design between the two primers of a pair.
//!
//! The rules are the published ones for real-time PCR probes (Applied Biosystems, *Primer Express
//! Version 1.5 and TaqMan MGB Probes for Allelic Discrimination*, and its TaqMan gene-expression
//! design notes): no G at the 5' end of the probe (a G next to the reporter dye quenches it even
//! after cleavage), no runs of identical bases and none of four or more G (the guidance says only
//! to avoid runs; the limit of four for the other bases is this program's choice), G+C between 30
//! and 80 %, a length of 30 or less, a melting temperature about 10 °C above the primers', and primers that sit
//! close to the probe without overlapping it. Probes of both strands are tried.
//!
//! The melting temperature is primer3's (`primer3::tm`, SantaLucia), the same one the primers get,
//! so the probe is placed against the primers by a difference, not against a number that belongs to
//! Primer Express's own calculator. The weights of the penalty are this program's own choice, and
//! minor-groove-binder probes are not modelled: their Tm cannot be had from nearest-neighbour
//! values.

use super::dna::{gc_percent, revcom};
use super::primer3::{self, Oligo, Pair, Settings as Primers};

/// How a probe is chosen. The defaults follow the guidelines named above.
#[derive(Clone, Debug, PartialEq)]
pub struct Settings {
    pub min_size: usize,
    pub opt_size: usize,
    pub max_size: usize,
    /// How far above the higher of the two primers' Tm the probe melts, in °C: least, best, most.
    pub min_above: f64,
    pub opt_above: f64,
    pub max_above: f64,
    /// G+C of the probe, percent.
    pub min_gc: f64,
    pub max_gc: f64,
    /// The longest run of one base allowed, and the longest run of G.
    pub max_run: usize,
    pub max_g_run: usize,
    /// The probe's first base (as it is ordered, 5' to 3') may not be G.
    pub no_five_prime_g: bool,
    /// The most the probe may stick to itself (primer3's scores, as for primers).
    pub max_self_any: i32,
    pub max_self_end: i32,
    /// Bases between the probe and the nearer primer, at the least.
    pub min_gap: usize,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            min_size: 18,
            opt_size: 24,
            max_size: 30,
            min_above: 6.0,
            opt_above: 10.0,
            max_above: 14.0,
            min_gc: 30.0,
            max_gc: 80.0,
            max_run: 4,
            max_g_run: 3,
            no_five_prime_g: true,
            max_self_any: 12,
            max_self_end: 12,
            min_gap: 1,
        }
    }
}

/// A pair of primers with the probe that goes with them.
#[derive(Clone, Debug, PartialEq)]
pub struct Assay {
    pub pair: Pair,
    pub probe: Oligo,
    /// The pair's penalty and the probe's.
    pub penalty: f64,
}

/// The weights of the penalty. The Tm term has weight 1, as primer3's does; the others are this program's own choice.
const SIZE_WEIGHT: f64 = 0.25;
const G_OVER_C_WEIGHT: f64 = 0.5;
const DISTANCE_WEIGHT: f64 = 0.05;
/// Pairs asked of primer3 for each assay wanted (a pair with no probe is dropped), and the fewest asked for.
const SPARE_PAIRS: usize = 8;
const LEAST_PAIRS: usize = 40;

/// The longest run of `base` in `seq`.
fn longest_run_of(seq: &[u8], base: u8) -> usize {
    seq.split(|&c| c != base)
        .map(<[u8]>::len)
        .max()
        .unwrap_or(0)
}

/// One candidate, judged by what does not depend on the primers; `None` when it breaks a limit. Its penalty so far is the distance
/// from the best size and the surplus of G over C (the strand with more C than G is the better one to order).
fn evaluate(
    seq: Vec<u8>,
    span: (usize, usize),
    reverse: bool,
    primers: &Primers,
    s: &Settings,
) -> Option<Oligo> {
    if seq.len() < 2 || !seq.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
        return None;
    }
    let gc = gc_percent(&seq);
    if (s.no_five_prime_g && seq[0] == b'G')
        || gc < s.min_gc
        || gc > s.max_gc
        || primer3::longest_run(&seq) > s.max_run
        || longest_run_of(&seq, b'G') > s.max_g_run
    {
        return None;
    }
    let (self_any, self_end) = primer3::dimer_scores(&seq, &seq);
    if self_any > s.max_self_any || self_end > s.max_self_end {
        return None;
    }
    let tm = primer3::tm(&seq, primers)?;
    let count = |base: u8| seq.iter().filter(|&&c| c == base).count();
    let surplus_g = count(b'G').saturating_sub(count(b'C'));
    let penalty = SIZE_WEIGHT * (seq.len() as f64 - s.opt_size as f64).abs()
        + G_OVER_C_WEIGHT * surplus_g as f64;
    Some(Oligo {
        seq,
        start: span.0,
        end: span.1,
        reverse,
        tm,
        gc,
        penalty,
        self_any,
        self_end,
    })
}

/// Every probe of either strand that lies in `from..to` of the template and breaks no limit of its own.
fn candidates(
    template: &[u8],
    from: usize,
    to: usize,
    primers: &Primers,
    s: &Settings,
) -> Vec<Oligo> {
    let to = to.min(template.len());
    let mut out = Vec::new();
    for len in s.min_size.max(2)..=s.max_size {
        for start in from..=to.saturating_sub(len) {
            if start + len > to {
                continue;
            }
            let plus = &template[start..start + len];
            for (seq, reverse) in [(plus.to_vec(), false), (revcom(plus), true)] {
                out.extend(evaluate(seq, (start, start + len), reverse, primers, s));
            }
        }
    }
    out
}

/// The probe's penalty with `pair`, or `None` when it does not belong between these primers: it must lie between them with the
/// gap asked for, melt the asked number of degrees above the higher primer, and not stick to either primer (the limits are the
/// primers' own, for a pair).
fn fits_pair(probe: &Oligo, pair: &Pair, primers: &Primers, s: &Settings) -> Option<f64> {
    let (left, right) = (&pair.left, &pair.right);
    if probe.start < left.end + s.min_gap || probe.end + s.min_gap > right.start {
        return None;
    }
    let above = probe.tm - left.tm.max(right.tm);
    if above < s.min_above || above > s.max_above {
        return None;
    }
    for primer in [left, right] {
        let (any, on_probe) = primer3::dimer_scores(&primer.seq, &probe.seq);
        let (_, on_primer) = primer3::dimer_scores(&probe.seq, &primer.seq);
        if any > primers.pair_compl_any || on_probe.max(on_primer) > primers.pair_compl_end {
            return None;
        }
    }
    let gap = (probe.start - left.end).min(right.start - probe.end);
    Some(probe.penalty + (above - s.opt_above).abs() + DISTANCE_WEIGHT * gap as f64)
}

/// The candidates that fit `pair`, with their penalties as the pair has them, best first.
fn best_first(candidates: &[Oligo], pair: &Pair, primers: &Primers, s: &Settings) -> Vec<Oligo> {
    let mut fit: Vec<Oligo> = candidates
        .iter()
        .filter_map(|c| {
            fits_pair(c, pair, primers, s).map(|penalty| Oligo {
                penalty,
                ..c.clone()
            })
        })
        .collect();
    fit.sort_by(|a, b| {
        a.penalty
            .total_cmp(&b.penalty)
            .then((a.start, a.reverse).cmp(&(b.start, b.reverse)))
    });
    fit
}

/// The probes that fit between the primers of `pair`, best first: every one, which `design_assays` takes the first of. (Used by the
/// tests, which compare it with a plain statement of the rules.)
#[cfg(test)]
fn probes_for(template: &[u8], pair: &Pair, primers: &Primers, s: &Settings) -> Vec<Oligo> {
    let template = template.to_ascii_uppercase();
    let found = candidates(&template, pair.left.end, pair.right.start, primers, s);
    best_first(&found, pair, primers, s)
}

/// Stretches of the template, given in any order, as the fewest that cover the same bases: those that overlap or touch are one.
fn joined(mut stretches: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    stretches.sort_unstable();
    let mut out: Vec<(usize, usize)> = Vec::new();
    for (a, b) in stretches {
        match out.last_mut() {
            Some(last) if a <= last.1 => last.1 = last.1.max(b),
            _ => out.push((a, b)),
        }
    }
    out
}

/// The best assays, lowest penalty first: pairs as `primer3::design` picks them (`primers.num_return` of them are wanted), each
/// with its best probe; a pair that has none is left out, and the others are asked for in spare.
pub fn design_assays(
    template: &[u8],
    target: Option<(usize, usize)>,
    primers: &Primers,
    s: &Settings,
) -> Vec<Assay> {
    let template = template.to_ascii_uppercase();
    let asked = Primers {
        num_return: (primers.num_return * SPARE_PAIRS).max(LEAST_PAIRS),
        ..primers.clone()
    };
    let pairs = primer3::design(&template, target, &asked);
    // the stretches between primers, joined where they meet, so that each base is looked at once
    let gaps: Vec<(usize, usize)> = pairs
        .iter()
        .map(|p| (p.left.end, p.right.start))
        .filter(|(a, b)| a < b)
        .collect();
    let found: Vec<Oligo> = joined(gaps)
        .into_iter()
        .flat_map(|(a, b)| candidates(&template, a, b, primers, s))
        .collect();
    let mut assays: Vec<Assay> = pairs
        .into_iter()
        .filter_map(|pair| {
            let probe = best_first(&found, &pair, primers, s).into_iter().next()?;
            Some(Assay {
                penalty: pair.penalty + probe.penalty,
                pair,
                probe,
            })
        })
        .collect();
    assays.sort_by(|a, b| a.penalty.total_cmp(&b.penalty));
    assays.truncate(primers.num_return);
    assays
}

#[cfg(test)]
mod tests;
