//! Assigning each alignment to a feature or a meta-feature, with the rules that featureCounts' users' guide gives.
//!
//! * An alignment **overlaps** a feature when at least one of its aligned bases (the `M`, `=` and `X` of its CIGAR: not what it deletes, skips,
//!   inserts or clips) lies in it. For a meta-feature, the overlap is the number of its aligned bases that lie in any of its features, each base
//!   counted once.
//! * The **filters** are applied in the guide's order, and an alignment that none stops is assigned: unmapped, mapping quality, duplicate,
//!   multi-mapping (`NH` above 1, unless `-M` or `--primary`), secondary (`--primary`), split or not split, no feature, overlap too short,
//!   ambiguous (more than one target and no `-O`). The paired filters of the guide (read type, singleton, chimera, fragment length) are not here:
//!   pairs are not counted as fragments yet.
//! * A **count** is 1 for each target of an assigned alignment; with `--fraction`, it is divided by `NH` (with `-M`) and by the number of
//!   targets (with `-O`).

use super::annotation::Annotation;
use crate::seq::gff::Strand;
use crate::seq::variants::sam::{Record, DUPLICATE, REVERSE, SECONDARY};
use std::collections::HashMap;

/// The flag of a read that is one of a pair.
pub const PAIRED: u16 = 0x1;

/// Which alignments are counted by whether their CIGAR has a skip (`N`, an exon junction).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Splice {
    Any,
    /// `--splitOnly`
    Split,
    /// `--nonSplitOnly`
    Whole,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Rules {
    /// Count features, not meta-features (`-f`).
    pub feature_level: bool,
    pub multi_overlap: bool,
    pub multi_mapping: bool,
    pub fraction: bool,
    pub primary_only: bool,
    pub min_mapq: u8,
    pub ignore_duplicates: bool,
    pub splice: Splice,
    pub min_overlap: usize,
    pub frac_overlap: f64,
    pub frac_overlap_feature: f64,
    pub non_overlap: Option<usize>,
    pub non_overlap_feature: Option<usize>,
    pub largest_overlap: bool,
    /// 0: either strand; 1: the read's strand is the feature's; 2: it is the other.
    pub strandedness: u8,
    /// Count each read of a pair on its own.
    pub each_read: bool,
}

impl Default for Rules {
    fn default() -> Self {
        Rules {
            feature_level: false,
            multi_overlap: false,
            multi_mapping: false,
            fraction: false,
            primary_only: false,
            min_mapq: 0,
            ignore_duplicates: false,
            splice: Splice::Any,
            min_overlap: 1,
            frac_overlap: 0.0,
            frac_overlap_feature: 0.0,
            non_overlap: None,
            non_overlap_feature: None,
            largest_overlap: false,
            strandedness: 0,
            each_read: false,
        }
    }
}

/// What became of an alignment, in the order of the filters.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Assigned,
    Unmapped,
    MappingQuality,
    Duplicate,
    MultiMapping,
    Secondary,
    Split,
    NonSplit,
    NoFeatures,
    OverlappingLength,
    Ambiguity,
}

impl Status {
    pub const ALL: [Status; 11] = [
        Status::Assigned,
        Status::Unmapped,
        Status::MappingQuality,
        Status::Duplicate,
        Status::MultiMapping,
        Status::Secondary,
        Status::Split,
        Status::NonSplit,
        Status::NoFeatures,
        Status::OverlappingLength,
        Status::Ambiguity,
    ];

    /// The name in the summary, as featureCounts writes it.
    pub fn name(self) -> &'static str {
        match self {
            Status::Assigned => "Assigned",
            Status::Unmapped => "Unassigned_Unmapped",
            Status::MappingQuality => "Unassigned_MappingQuality",
            Status::Duplicate => "Unassigned_Duplicate",
            Status::MultiMapping => "Unassigned_MultiMapping",
            Status::Secondary => "Unassigned_Secondary",
            Status::Split => "Unassigned_Split",
            Status::NonSplit => "Unassigned_NonSplit",
            Status::NoFeatures => "Unassigned_NoFeatures",
            Status::OverlappingLength => "Unassigned_Overlapping_Length",
            Status::Ambiguity => "Unassigned_Ambiguity",
        }
    }
}

/// The total length of the union of intervals.
fn union_length(intervals: &mut [(usize, usize)]) -> usize {
    intervals.sort_unstable();
    let (mut total, mut end) = (0, 0);
    for &mut (a, b) in intervals {
        total += b.saturating_sub(a.max(end));
        end = end.max(b);
    }
    total
}

/// The features of one sequence by start, with the greatest end among each and those before it.
#[derive(Default)]
struct Sorted {
    features: Vec<usize>,
    reach: Vec<usize>,
}

pub struct Index<'a> {
    pub annotation: &'a Annotation,
    sequences: HashMap<&'a str, Sorted>,
    /// The length of each target: a feature's, or the bases of a meta-feature's features counted once.
    pub lengths: Vec<usize>,
    feature_level: bool,
}

impl<'a> Index<'a> {
    pub fn new(annotation: &'a Annotation, feature_level: bool) -> Index<'a> {
        let mut sequences: HashMap<&str, Sorted> = HashMap::new();
        for (i, f) in annotation.features.iter().enumerate() {
            sequences.entry(&f.seqname).or_default().features.push(i);
        }
        for s in sequences.values_mut() {
            s.features
                .sort_by_key(|&i| (annotation.features[i].start, i));
            let mut reach = 0;
            s.reach = s
                .features
                .iter()
                .map(|&i| {
                    reach = reach.max(annotation.features[i].end);
                    reach
                })
                .collect();
        }
        let lengths = if feature_level {
            annotation
                .features
                .iter()
                .map(|f| f.end - f.start)
                .collect()
        } else {
            // by meta-feature and sequence: features on different sequences share no base
            let mut parts: Vec<HashMap<&str, Vec<(usize, usize)>>> =
                vec![HashMap::new(); annotation.names.len()];
            for f in &annotation.features {
                parts[f.meta]
                    .entry(&f.seqname)
                    .or_default()
                    .push((f.start, f.end));
            }
            parts
                .into_iter()
                .map(|by_seq| by_seq.into_values().map(|mut v| union_length(&mut v)).sum())
                .collect()
        };
        Index {
            annotation,
            sequences,
            lengths,
            feature_level,
        }
    }

    /// Does the annotation have a feature on this sequence?
    pub fn names(&self, seqname: &str) -> bool {
        self.sequences.contains_key(seqname)
    }

    /// The features of `seqname` that share a base with `[from, to)`.
    fn overlapping(&self, seqname: &str, from: usize, to: usize) -> Vec<usize> {
        let Some(s) = self.sequences.get(seqname) else {
            return Vec::new();
        };
        let features = &self.annotation.features;
        let before = s.features.partition_point(|&i| features[i].start < to);
        (0..before)
            .rev()
            .take_while(|&k| s.reach[k] > from)
            .map(|k| s.features[k])
            .filter(|&i| features[i].end > from)
            .collect()
    }

    /// Does the strand of the read fit the feature's? A feature with no strand fits either.
    fn strand_fits(strand: Strand, reverse: bool, rules: &Rules) -> bool {
        match (rules.strandedness, strand) {
            (0, _) | (_, Strand::Unknown) => true,
            (1, s) => (s == Strand::Minus) == reverse,
            (_, s) => (s == Strand::Minus) != reverse,
        }
    }

    /// What becomes of one alignment, and the targets it is assigned to (when it is).
    pub fn assign(&self, r: &Record, rules: &Rules) -> (Status, Vec<usize>) {
        let stop = |status| (status, Vec::new());
        if !r.mapped() {
            return stop(Status::Unmapped);
        }
        if r.mapq < rules.min_mapq {
            return stop(Status::MappingQuality);
        }
        if rules.ignore_duplicates && r.flag & DUPLICATE != 0 {
            return stop(Status::Duplicate);
        }
        if !rules.multi_mapping && !rules.primary_only && hits(r) > 1 {
            return stop(Status::MultiMapping);
        }
        if rules.primary_only && r.flag & SECONDARY != 0 {
            return stop(Status::Secondary);
        }
        let split = r.cigar.iter().any(|&(op, _)| op == 'N');
        match rules.splice {
            Splice::Whole if split => return stop(Status::Split),
            Splice::Split if !split => return stop(Status::NonSplit),
            _ => {}
        }
        let name = r.target.as_deref().unwrap_or_default();
        let reverse = r.flag & REVERSE != 0;
        let mut found: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        let (mut at, mut read_length) = (r.pos, 0);
        for &(op, n) in &r.cigar {
            if matches!(op, 'M' | '=' | 'X') {
                for i in self.overlapping(name, at, at + n) {
                    let f = &self.annotation.features[i];
                    if Index::strand_fits(f.strand, reverse, rules) {
                        let target = if self.feature_level { i } else { f.meta };
                        found
                            .entry(target)
                            .or_default()
                            .push((at.max(f.start), (at + n).min(f.end)));
                    }
                }
            }
            if matches!(op, 'M' | 'I' | 'S' | '=' | 'X') {
                read_length += n;
            }
            if matches!(op, 'M' | 'D' | 'N' | '=' | 'X') {
                at += n;
            }
        }
        if found.is_empty() {
            return stop(Status::NoFeatures);
        }
        let mut passing: Vec<(usize, usize)> = found
            .into_iter()
            .filter_map(|(target, mut parts)| {
                let overlap = union_length(&mut parts);
                self.enough(overlap, read_length, self.lengths[target], rules)
                    .then_some((target, overlap))
            })
            .collect();
        if passing.is_empty() {
            return stop(Status::OverlappingLength);
        }
        if rules.largest_overlap {
            let most = passing.iter().map(|p| p.1).max().unwrap_or(0);
            passing.retain(|p| p.1 == most);
        }
        if passing.len() > 1 && !rules.multi_overlap {
            return stop(Status::Ambiguity);
        }
        let mut targets: Vec<usize> = passing.into_iter().map(|p| p.0).collect();
        targets.sort_unstable();
        (Status::Assigned, targets)
    }

    /// Is the overlap enough, in all the ways that were asked for?
    fn enough(&self, overlap: usize, read: usize, target: usize, rules: &Rules) -> bool {
        // a share that is exactly met counts, whatever the rounding of the product
        let share = |part: f64, of: usize| overlap as f64 + 1e-9 >= part * of as f64;
        overlap >= rules.min_overlap
            && share(rules.frac_overlap, read)
            && share(rules.frac_overlap_feature, target)
            && rules
                .non_overlap
                .is_none_or(|n| read.saturating_sub(overlap) <= n)
            && rules
                .non_overlap_feature
                .is_none_or(|n| target.saturating_sub(overlap) <= n)
    }
}

/// `NH`, the number of places that the read was reported at: 1 when the line does not say.
fn hits(r: &Record) -> usize {
    r.tags
        .iter()
        .find_map(|t| t.strip_prefix("NH:i:"))
        .and_then(|n| n.parse().ok())
        .unwrap_or(1)
}

/// The counts of one file: one for each target, and how many alignments had each status.
#[derive(Clone, Debug, PartialEq)]
pub struct Tally {
    pub counts: Vec<f64>,
    pub statuses: [u64; 11],
}

impl Tally {
    pub fn of(&self, status: Status) -> u64 {
        self.statuses[status as usize]
    }
}

/// Counts the alignments of one file. `Err` when they are of pairs and each read was not asked to be counted on its own: a pair would be
/// counted twice.
pub fn count(index: &Index<'_>, records: &[Record], rules: &Rules) -> Result<Tally, String> {
    if !rules.each_read && records.iter().any(|r| r.flag & PAIRED != 0) {
        return Err("The reads are of pairs, and counting a pair as one fragment is not built yet: --each-read counts each read of a pair on its own.".into());
    }
    let mut tally = Tally {
        counts: vec![0.0; index.lengths.len()],
        statuses: [0; 11],
    };
    for r in records {
        let (status, targets) = index.assign(r, rules);
        tally.statuses[status as usize] += 1;
        let mut weight = 1.0;
        if rules.fraction {
            if rules.multi_mapping && !rules.primary_only {
                weight /= hits(r).max(1) as f64;
            }
            if rules.multi_overlap {
                weight /= targets.len().max(1) as f64;
            }
        }
        for t in targets {
            tally.counts[t] += weight;
        }
    }
    Ok(tally)
}

#[cfg(test)]
mod tests;
