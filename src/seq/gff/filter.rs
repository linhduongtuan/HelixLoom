//! Which transcripts to keep: `gffread`'s filters `-C --nc -U -i -l -r -R -N -V -J --ids --nids`, as its usage text describes them.

use super::extract::{cds_sequence, cut, Genome};
use super::{Strand, Transcript};
use crate::seq::dna::{revcom, translate, STANDARD_CODE};
use std::collections::HashSet;

/// `[strand]seqname:start-end`.
#[derive(Clone, Debug, PartialEq)]
pub struct Range {
    pub seqname: String,
    pub start: usize,
    pub end: usize,
    pub strand: Option<Strand>,
}

impl Range {
    pub fn parse(text: &str) -> Result<Range, String> {
        let bad = || {
            format!("{text:?} is not a range: write [strand]sequence:start-end, for example chr1:1000-2000.")
        };
        let (strand, rest) = match text.as_bytes().first() {
            Some(b'+') => (Some(Strand::Plus), &text[1..]),
            Some(b'-') => (Some(Strand::Minus), &text[1..]),
            _ => (None, text),
        };
        let (seqname, span) = rest.rsplit_once(':').ok_or_else(bad)?;
        let (start, end) = span.split_once('-').ok_or_else(bad)?;
        let (start, end): (usize, usize) = (
            start.replace(',', "").parse().map_err(|_| bad())?,
            end.replace(',', "").parse().map_err(|_| bad())?,
        );
        if seqname.is_empty() || start == 0 || start > end {
            return Err(bad());
        }
        Ok(Range {
            seqname: seqname.to_string(),
            start,
            end,
            strand,
        })
    }
}

#[derive(Clone, Debug, Default)]
pub struct Filters {
    /// `-C`: only transcripts that have CDS.
    pub coding_only: bool,
    /// `--nc`: only transcripts that have none.
    pub noncoding_only: bool,
    /// `-U`: no transcript of one exon.
    pub no_single_exon: bool,
    /// `-l`: no transcript with an intron shorter than this; `-i`: none with one longer.
    pub min_intron: Option<usize>,
    pub max_intron: Option<usize>,
    /// `-r`: transcripts that overlap the range; with `true` (`-R`), only those that lie wholly inside it.
    pub range: Option<(Range, bool)>,
    pub only_ids: Option<HashSet<String>>,
    pub not_ids: Option<HashSet<String>>,
    /// `-N`: no transcript with an intron that does not begin and end GT-AG, GC-AG or AT-AC (needs the genome).
    pub canonical_splice: bool,
    /// `-V`: no transcript whose CDS has a stop inside it (needs the genome).
    pub no_inframe_stop: bool,
    /// `-J`: only transcripts whose CDS begins with ATG, ends with a stop codon, and has none inside (needs the genome).
    pub complete_cds: bool,
}

const SPLICE_SITES: [(&[u8; 2], &[u8; 2]); 3] = [(b"GT", b"AG"), (b"GC", b"AG"), (b"AT", b"AC")];

impl Filters {
    fn needs_genome(&self) -> bool {
        self.canonical_splice || self.no_inframe_stop || self.complete_cds
    }

    /// Whether the transcript is kept; why not, when a filter needs the genome and none was given, or the genome lacks what is cut.
    pub fn keeps(&self, t: &Transcript, genome: Option<&Genome>) -> Result<bool, String> {
        if self.needs_genome() && genome.is_none() {
            return Err(
                "The filters -N, -V and -J read the sequence: give the genome with -g.".into(),
            );
        }
        if (self.coding_only && !t.is_coding()) || (self.noncoding_only && t.is_coding()) {
            return Ok(false);
        }
        if self.no_single_exon && t.exons.len() < 2 {
            return Ok(false);
        }
        let introns = t.introns();
        let long = |n: usize| introns.iter().any(|(a, b)| b - a + 1 > n);
        let short = |n: usize| introns.iter().any(|(a, b)| b - a + 1 < n);
        if self.max_intron.is_some_and(long) || self.min_intron.is_some_and(short) {
            return Ok(false);
        }
        if let Some((r, inside)) = &self.range {
            let strand_ok = r.strand.is_none_or(|s| s == t.strand);
            let place_ok = if *inside {
                r.start <= t.start && t.end <= r.end
            } else {
                t.start <= r.end && t.end >= r.start
            };
            if t.seqname != r.seqname || !strand_ok || !place_ok {
                return Ok(false);
            }
        }
        if self
            .only_ids
            .as_ref()
            .is_some_and(|ids| !ids.contains(&t.id))
            || self.not_ids.as_ref().is_some_and(|ids| ids.contains(&t.id))
        {
            return Ok(false);
        }
        let Some(genome) = genome else {
            return Ok(true);
        };
        if self.canonical_splice && !canonical(t, genome)? {
            return Ok(false);
        }
        if self.no_inframe_stop || self.complete_cds {
            let cds = cds_sequence(t, genome)?;
            if self.complete_cds && !complete(cds.as_deref()) {
                return Ok(false);
            }
            if self.no_inframe_stop && !cds.as_deref().is_none_or(no_stop_inside) {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

/// Whether every intron begins and ends as an intron usually does, read on the strand of the transcript.
fn canonical(t: &Transcript, genome: &Genome) -> Result<bool, String> {
    for (a, b) in t.introns() {
        let mut intron = cut(genome, &t.seqname, a, b)?.to_vec();
        if t.strand == Strand::Minus {
            intron = revcom(&intron);
        }
        if intron.len() < 4 {
            return Ok(false);
        }
        let (donor, acceptor) = (&intron[..2], &intron[intron.len() - 2..]);
        if !SPLICE_SITES
            .iter()
            .any(|(d, a)| donor == &d[..] && acceptor == &a[..])
        {
            return Ok(false);
        }
    }
    Ok(true)
}

/// No stop codon before the last codon.
fn no_stop_inside(cds: &[u8]) -> bool {
    let aa = translate(cds, STANDARD_CODE);
    let inside = aa.len().saturating_sub(1);
    !aa[..inside].contains(&b'*')
}

/// A CDS that is whole codons, begins with a start codon and ends with a stop, with none between.
fn complete(cds: Option<&[u8]>) -> bool {
    let Some(cds) = cds else {
        return false;
    };
    let aa = translate(cds, STANDARD_CODE);
    cds.len() % 3 == 0
        && aa.first() == Some(&b'M')
        && cds.len() >= 3
        && cds[..3].eq_ignore_ascii_case(b"ATG")
        && aa.last() == Some(&b'*')
        && !aa[..aa.len() - 1].contains(&b'*')
}
