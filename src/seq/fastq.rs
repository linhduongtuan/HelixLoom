//! FASTQ reading, quality statistics of the kind FastQC reports, and read cleaning in the manner
//! of fastp. The defaults, the adapter match and the sliding-window cuts were read from fastp's
//! source (github.com/OpenGene/fastp: options.h, filter.cpp, adaptertrimmer.cpp) and re-written;
//! no fastp code was copied. Not ported: paired reads, adapter detection (an adapter is given),
//! fastp's one-gap adapter fallback, poly-G/X trimming, the low-complexity filter, deduplication.

/// Phred+33, the encoding of current Illumina data.
const PHRED_OFFSET: u8 = 33;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastqRead {
    pub name: String,
    pub seq: Vec<u8>,
    /// Phred quality of each base.
    pub qual: Vec<u8>,
}

/// Records of four lines: `@name`, bases, `+`, qualities.
pub fn parse_fastq(text: &str) -> Result<Vec<FastqRead>, String> {
    let lines: Vec<&str> = text.lines().map(|l| l.trim_end_matches('\r')).collect();
    let mut reads = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        if lines[i].is_empty() {
            i += 1;
            continue;
        }
        let number = reads.len() + 1;
        let record = lines
            .get(i..i + 4)
            .ok_or_else(|| format!("Record {number} is cut short."))?;
        let (head, seq, plus, qual) = (record[0], record[1], record[2], record[3]);
        let name = head
            .strip_prefix('@')
            .ok_or_else(|| format!("Record {number} doesn't start with @."))?;
        if !plus.starts_with('+') {
            return Err(format!("Record {number} has no + line."));
        }
        if seq.len() != qual.len() {
            return Err(format!(
                "Record {number}: bases and qualities differ in length."
            ));
        }
        if qual.bytes().any(|q| q < PHRED_OFFSET) {
            return Err(format!(
                "Record {number}: quality characters below '!' (not Phred+33)."
            ));
        }
        reads.push(FastqRead {
            name: name.to_string(),
            seq: seq.bytes().map(|b| b.to_ascii_uppercase()).collect(),
            qual: qual.bytes().map(|q| q - PHRED_OFFSET).collect(),
        });
        i += 4;
    }
    if reads.is_empty() {
        return Err("No FASTQ records found.".into());
    }
    Ok(reads)
}

pub fn to_fastq(reads: &[FastqRead]) -> String {
    reads
        .iter()
        .map(|r| {
            let qual: String = r
                .qual
                .iter()
                .map(|&q| char::from(q + PHRED_OFFSET))
                .collect();
            format!(
                "@{}\n{}\n+\n{qual}\n",
                r.name,
                String::from_utf8_lossy(&r.seq)
            )
        })
        .collect()
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Stats {
    pub reads: usize,
    pub bases: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub gc_percent: f64,
    pub mean_quality: f64,
    /// Share of bases with Phred quality of at least 20 and 30.
    pub q20_percent: f64,
    pub q30_percent: f64,
    /// Mean quality at each position, over the reads that reach it.
    pub per_position: Vec<f64>,
}

pub fn stats(reads: &[FastqRead]) -> Stats {
    let bases: usize = reads.iter().map(|r| r.seq.len()).sum();
    if bases == 0 {
        return Stats {
            reads: reads.len(),
            ..Stats::default()
        };
    }
    let share = |count: usize| 100.0 * count as f64 / bases as f64;
    let all = || {
        reads
            .iter()
            .flat_map(|r| r.qual.iter().map(|&q| usize::from(q)))
    };
    let longest = reads.iter().map(|r| r.seq.len()).max().unwrap_or(0);
    let (mut sum, mut count) = (vec![0usize; longest], vec![0usize; longest]);
    for r in reads {
        for (p, &q) in r.qual.iter().enumerate() {
            sum[p] += usize::from(q);
            count[p] += 1;
        }
    }
    Stats {
        reads: reads.len(),
        bases,
        min_length: reads.iter().map(|r| r.seq.len()).min().unwrap_or(0),
        max_length: longest,
        gc_percent: share(
            reads
                .iter()
                .flat_map(|r| &r.seq)
                .filter(|b| matches!(b, b'G' | b'C'))
                .count(),
        ),
        mean_quality: all().sum::<usize>() as f64 / bases as f64,
        q20_percent: share(all().filter(|&q| q >= 20).count()),
        q30_percent: share(all().filter(|&q| q >= 30).count()),
        per_position: sum
            .iter()
            .zip(&count)
            .map(|(&s, &c)| s as f64 / c as f64)
            .collect(),
    }
}

/// fastp's defaults (`options.h`) with the sliding-window cuts off, as there. fastp has no default
/// adapter (it detects one); the Illumina TruSeq read-1 adapter is offered here.
#[derive(Clone, Debug, PartialEq)]
pub struct Cleaning {
    /// Adapter to cut from the 3' end; empty for none.
    pub adapter: Vec<u8>,
    pub cut_front: bool,
    pub cut_tail: bool,
    pub cut_right: bool,
    /// `cut_window_size` and `cut_mean_quality`.
    pub window: usize,
    pub window_quality: usize,
    /// `qualified_quality_phred`: a base below this is low quality.
    pub qualified_quality: u8,
    /// `unqualified_percent_limit`: reads with more low-quality bases than this share are dropped.
    pub unqualified_percent: usize,
    /// `n_base_limit`: reads with more N than this are dropped.
    pub n_base_limit: usize,
    /// `average_qual`: reads with a lower mean quality are dropped; 0 for off.
    pub average_quality: usize,
    /// `length_required`: shorter reads are dropped.
    pub min_length: usize,
}

impl Default for Cleaning {
    fn default() -> Self {
        Cleaning {
            adapter: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_vec(),
            cut_front: false,
            cut_tail: false,
            cut_right: false,
            window: 4,
            window_quality: 20,
            qualified_quality: 15,
            unqualified_percent: 40,
            n_base_limit: 5,
            average_quality: 0,
            min_length: 15,
        }
    }
}

/// Where an adapter starts in a read (fastp's `trimBySequence`): the adapter slides along the read,
/// beginning a few bases before it for a long adapter, needs four bases of overlap and may differ
/// at one base in eight. A negative start means the read begins inside the adapter.
fn adapter_start(read: &[u8], adapter: &[u8]) -> Option<isize> {
    const OVERLAP_REQUIRED: isize = 4;
    const ONE_MISMATCH_FOR_EACH: isize = 8;
    let (rlen, alen) = (read.len() as isize, adapter.len() as isize);
    if alen < OVERLAP_REQUIRED {
        return None;
    }
    let first = match alen {
        16.. => -4,
        12.. => -3,
        8.. => -2,
        _ => 0,
    };
    (first..rlen - OVERLAP_REQUIRED).find(|&pos| {
        let compared = (rlen - pos).min(alen);
        let skip = (-pos).max(0);
        let differ = (skip..compared)
            .filter(|&i| adapter[i as usize] != read[(i + pos) as usize])
            .count() as isize;
        differ <= compared / ONE_MISMATCH_FOR_EACH
    })
}

/// The sliding-window cuts of fastp's `trimAndCut`: `(first base kept, bases kept)`, or `None`
/// when nothing is left. A window is a run of `window` bases; its mean quality must reach
/// `quality`.
fn quality_cut(read: &FastqRead, c: &Cleaning) -> Option<(usize, usize)> {
    let (l, w, q) = (read.seq.len(), c.window.max(1), c.window_quality);
    let qual = |i: usize| usize::from(read.qual[i]);
    let seq = &read.seq;
    let (mut front, mut keep) = (0usize, l);
    if l <= w && (c.cut_front || c.cut_tail || c.cut_right) {
        return None;
    }
    if c.cut_front {
        let mut s = 0;
        while s + w < l {
            if (s..s + w).map(qual).sum::<usize>() >= w * q {
                break;
            }
            s += 1;
        }
        if s > 0 {
            s += w - 1;
        }
        while s < l && seq[s] == b'N' {
            s += 1;
        }
        front = s;
        keep = l.saturating_sub(front);
    }
    if c.cut_right {
        let mut s = front;
        let low_window = loop {
            if s + w >= l {
                break false;
            }
            if (s..s + w).map(qual).sum::<usize>() < w * q {
                break true;
            }
            s += 1;
        };
        if low_window {
            while s + 1 < l && qual(s) >= q {
                s += 1;
            }
            keep = s.saturating_sub(front);
        }
    } else if c.cut_tail {
        // from the 3' end backwards: the first window, going left, that reaches the quality
        let mut t = l - 1;
        while t >= front + w {
            if (t + 1 - w..=t).map(qual).sum::<usize>() >= w * q {
                break;
            }
            t -= 1;
        }
        if t < l - 1 {
            t = (t + 1).saturating_sub(w);
        }
        let mut end = t as isize;
        while end >= 0 && seq[end as usize] == b'N' {
            end -= 1;
        }
        keep = (end + 1 - front as isize).max(0) as usize;
    }
    (keep > 0 && front + 1 < l).then_some((front, keep))
}

/// What happened to the reads.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Report {
    pub too_many_low_quality: usize,
    pub too_many_n: usize,
    pub too_short: usize,
    pub adapter_trimmed: usize,
}

/// Cut, trim and filter reads: sliding-window cuts, the adapter, then fastp's filters in its
/// order (low-quality share, mean quality, N count, length). Returns the reads that pass.
pub fn clean(reads: &[FastqRead], c: &Cleaning) -> (Vec<FastqRead>, Report) {
    let mut report = Report::default();
    let mut kept = Vec::new();
    let adapter = c.adapter.to_ascii_uppercase();
    for read in reads {
        let Some((from, len)) = quality_cut(read, c) else {
            report.too_short += 1;
            continue;
        };
        let (mut seq, mut qual) = (
            read.seq[from..from + len].to_vec(),
            read.qual[from..from + len].to_vec(),
        );
        if let Some(pos) = adapter_start(&seq, &adapter) {
            report.adapter_trimmed += 1;
            let keep = usize::try_from(pos).unwrap_or(0);
            seq.truncate(keep);
            qual.truncate(keep);
        }
        if seq.is_empty() {
            report.too_short += 1;
            continue;
        }
        let n = seq.len();
        let low = qual.iter().filter(|&&q| q < c.qualified_quality).count();
        let total: usize = qual.iter().map(|&q| usize::from(q)).sum();
        let ns = seq.iter().filter(|&&b| b == b'N').count();
        if low as f64 > (c.unqualified_percent * n) as f64 / 100.0
            || (c.average_quality > 0 && total / n < c.average_quality)
        {
            report.too_many_low_quality += 1;
        } else if ns > c.n_base_limit {
            report.too_many_n += 1;
        } else if n < c.min_length {
            report.too_short += 1;
        } else {
            kept.push(FastqRead {
                name: read.name.clone(),
                seq,
                qual,
            });
        }
    }
    (kept, report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(seq: &str, qual: &[u8]) -> FastqRead {
        FastqRead {
            name: "r".into(),
            seq: seq.bytes().collect(),
            qual: qual.to_vec(),
        }
    }

    fn flat(seq: &str, q: u8) -> FastqRead {
        read(seq, &vec![q; seq.len()])
    }

    #[test]
    fn records_parse_with_an_at_sign_in_the_quality_line_and_round_trip() {
        let text = "@r1 first\nACGTN\n+\n@III#\n@r2\nacgt\n+r2\nIIII\n\n";
        let reads = parse_fastq(text).unwrap();
        assert_eq!(reads.len(), 2);
        assert_eq!(
            (reads[0].name.as_str(), reads[0].seq.as_slice()),
            ("r1 first", &b"ACGTN"[..])
        );
        assert_eq!(reads[0].qual, vec![31, 40, 40, 40, 2]);
        assert_eq!(reads[1].seq, b"ACGT");
        assert_eq!(to_fastq(&reads[..1]), "@r1 first\nACGTN\n+\n@III#\n");
        assert!(parse_fastq("@r\nACGT\n+\nIII\n").is_err());
        assert!(parse_fastq("r\nACGT\n+\nIIII\n").is_err());
        assert!(parse_fastq("@r\nACGT\n").is_err());
        assert!(parse_fastq("").is_err());
    }

    #[test]
    fn statistics_count_bases_quality_and_gc() {
        let reads = [read("GGCC", &[40, 40, 30, 10]), read("AT", &[20, 20])];
        let s = stats(&reads);
        assert_eq!((s.reads, s.bases, s.min_length, s.max_length), (2, 6, 2, 4));
        assert!((s.gc_percent - 400.0 / 6.0).abs() < 1e-9);
        assert!((s.mean_quality - 160.0 / 6.0).abs() < 1e-9);
        assert!((s.q30_percent - 300.0 / 6.0).abs() < 1e-9);
        assert!((s.q20_percent - 500.0 / 6.0).abs() < 1e-9);
        assert_eq!(s.per_position, vec![30.0, 30.0, 30.0, 10.0]);
        assert_eq!(stats(&[]).reads, 0);
    }

    // fastp's own self-test for its adapter trimmer
    #[test]
    fn an_adapter_is_cut_where_it_starts_even_with_a_mismatch() {
        let read = b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG";
        let adapter = b"TTTTCCACGGGGATACTACTG";
        assert_eq!(adapter_start(read, adapter), Some(44));
        assert_eq!(adapter_start(b"ACGTACGTACGTACGTACGTACGT", adapter), None);
        // a read that begins inside the adapter is dropped whole (negative start)
        assert!(adapter_start(
            b"GATCGGAAGAGCACACGTCTGAACTCCAGTCAAAA",
            &Cleaning::default().adapter
        )
        .is_some_and(|p| p <= 0));
    }

    #[test]
    fn the_default_filters_are_fastps() {
        let c = Cleaning::default();
        assert_eq!(
            (
                c.qualified_quality,
                c.unqualified_percent,
                c.n_base_limit,
                c.min_length
            ),
            (15, 40, 5, 15)
        );
        assert_eq!((c.window, c.window_quality, c.average_quality), (4, 20, 0));
        let good = flat("ACGTACGTACGTACGTACGT", 35);
        let short = flat("ACGTACGT", 35);
        let ns = flat("ACGTNNNNNNACGTACGTAC", 35);
        let low = read("ACGTACGTACGTACGTACGT", &[10; 20]);
        let mostly_ok = read(
            "ACGTACGTACGTACGTACGT",
            &[
                10, 10, 10, 10, 10, 10, 10, 10, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
            ],
        );
        let (kept, report) = clean(&[good.clone(), short, ns, low, mostly_ok.clone()], &c);
        assert_eq!(kept, vec![good, mostly_ok]); // 8 of 20 low bases is 40 %, not more than 40 %
        assert_eq!(
            report,
            Report {
                too_many_low_quality: 1,
                too_many_n: 1,
                too_short: 1,
                adapter_trimmed: 0
            }
        );
    }

    #[test]
    fn adapters_are_trimmed_then_the_short_remainder_is_dropped() {
        let adapter = Cleaning::default().adapter;
        let insert = b"ACGTTGCAAGGCTTAACCGGTTACGTAC"; // 28 bases
        let mut seq = insert.to_vec();
        seq.extend_from_slice(&adapter[..20]);
        let long = FastqRead {
            name: "long".into(),
            qual: vec![35; seq.len()],
            seq,
        };
        let mut tiny = b"ACGTACGTAC".to_vec(); // 10 bases of insert: too short once trimmed
        tiny.extend_from_slice(&adapter[..20]);
        let tiny = FastqRead {
            name: "tiny".into(),
            qual: vec![35; tiny.len()],
            seq: tiny,
        };
        let (kept, report) = clean(&[long, tiny], &Cleaning::default());
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].seq, insert);
        assert_eq!((report.adapter_trimmed, report.too_short), (2, 1));
    }

    #[test]
    fn the_sliding_window_cuts_follow_fastp() {
        // ten good bases, then ten bad ones
        let mut qual = vec![30u8; 10];
        qual.extend(vec![5u8; 10]);
        let r = read("ACGTACGTACGTACGTACGT", &qual);
        let right = Cleaning {
            cut_right: true,
            min_length: 5,
            ..Cleaning::default()
        };
        assert_eq!(quality_cut(&r, &right), Some((0, 10)));
        let tail = Cleaning {
            cut_tail: true,
            ..right.clone()
        };
        let (kept, _) = clean(std::slice::from_ref(&r), &tail);
        assert_eq!(kept.len(), 1);
        assert!(
            kept[0].seq.len() >= 6 && kept[0].seq.len() <= 10,
            "{}",
            kept[0].seq.len()
        );
        // good qualities are left alone, and with no cut enabled nothing changes
        assert_eq!(
            quality_cut(&flat("ACGTACGTACGT", 35), &right),
            Some((0, 12))
        );
        assert_eq!(quality_cut(&r, &Cleaning::default()), Some((0, 20)));
    }
}
