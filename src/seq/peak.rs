//! What a trace says of the bases that are picked in it (`docs/ROADMAP.md` 1.4): of one base its call, how sure the machine
//! was of it and the signal of the four channels at its peak; of several, how good they are together.

use super::abif::Trace;

/// Below this quality a call is counted as poor: a chance of one in a hundred or more that it is wrong.
pub const POOR: u8 = 20;

/// One base of a trace.
#[derive(Clone, Debug, PartialEq)]
pub struct Peak {
    /// The base, counted from 0 in the order of the file.
    pub at: usize,
    pub call: u8,
    /// How sure the machine was of the call; none from a file that gives no qualities.
    pub quality: Option<u8>,
    /// The scan the peak lies at; none when the file names no peak for the base.
    pub scan: Option<u32>,
    /// The signal of each channel at that scan with the base the channel stands for: the strongest first, and those as
    /// strong in the order of the file. Empty without a peak.
    pub signal: Vec<(u8, i16)>,
}

/// Several bases of a trace, `from` up to `to` (which is not one of them).
#[derive(Clone, Debug, PartialEq)]
pub struct Stretch {
    pub from: usize,
    pub to: usize,
    /// The mean of their qualities; none from a file that gives no qualities, as the two below.
    pub mean: Option<f64>,
    /// The lowest quality, and the first base that has it.
    pub lowest: Option<(u8, usize)>,
    /// How many of the calls are poor (below `POOR`).
    pub poor: Option<usize>,
    /// How many of the calls are not A, C, G or T.
    pub unsure: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Picked {
    One(Peak),
    Several(Stretch),
}

/// How many calls of this quality hold one wrong one: the quality is ten times the logarithm of it (10 is one in ten, 20
/// one in a hundred).
pub fn one_in(quality: u8) -> u64 {
    10f64.powf(f64::from(quality) / 10.0).round() as u64
}

/// What the trace says of its bases `from` up to `to` (not counted); none when it has none of them.
pub fn picked(t: &Trace, from: usize, to: usize) -> Option<Picked> {
    let to = to.min(t.bases.len());
    if from >= to {
        return None;
    }
    let scored = t.has_quality();
    // a base past the end of the qualities (a file that is cut short) has none
    let quality = |at: usize| t.quality.get(at).copied().filter(|_| scored);
    if to - from == 1 {
        let scan = t.peaks.get(from).copied();
        let mut signal: Vec<(u8, i16)> = scan.map_or_else(Vec::new, |scan| {
            (0..4)
                .map(|k| {
                    let here = t.channels[k].get(scan as usize).copied().unwrap_or(0);
                    (t.order[k], here)
                })
                .collect()
        });
        // a stable sort: those as strong keep the order of the file
        signal.sort_by_key(|&(_, s)| std::cmp::Reverse(s));
        return Some(Picked::One(Peak {
            at: from,
            call: t.bases[from],
            quality: quality(from),
            scan,
            signal,
        }));
    }
    let qualities: Vec<(u8, usize)> = (from..to)
        .filter_map(|at| quality(at).map(|q| (q, at)))
        .collect();
    let any = !qualities.is_empty();
    Some(Picked::Several(Stretch {
        from,
        to,
        mean: any.then(|| {
            qualities.iter().map(|&(q, _)| f64::from(q)).sum::<f64>() / qualities.len() as f64
        }),
        // the first of the lowest: `min_by_key` keeps the first of equals
        lowest: qualities.iter().copied().min_by_key(|&(q, _)| q),
        poor: any.then(|| qualities.iter().filter(|&&(q, _)| q < POOR).count()),
        unsure: t.bases[from..to]
            .iter()
            .filter(|b| !b"ACGT".contains(&b.to_ascii_uppercase()))
            .count(),
    }))
}

#[cfg(test)]
mod tests;
