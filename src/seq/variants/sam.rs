//! SAM, the text format of alignments (the specification of the SAM/BAM working group, "Sequence Alignment/Map Format Specification", written
//! from it: nothing of samtools was read). One line to an alignment, eleven columns separated by tabs, then tags: 1 read name, 2 flag, 3 target
//! name, 4 position (from 1), 5 mapping quality, 6 CIGAR, 7 to 9 the mate (`*`, 0, 0 here: no pairs), 10 bases, 11 qualities (Phred + 33).
//!
//! What is here: the lines that the mapper's hits make (`Record::from_hit`, the read as it lies against the target with the ends that were not
//! aligned as soft clips), the reading of a SAM file (`parse`), and the header. Not done: BAM and CRAM (a binary format in blocks of deflate),
//! the mate of a pair, the `SA` tag of a split read, hard clips on output, the `MD` tag.

use crate::seq::dna::revcom;
use crate::seq::minimap::{Hit, Index};
use crate::seq::variants::reads::Read;

pub const UNMAPPED: u16 = 0x4;
pub const REVERSE: u16 = 0x10;
pub const SECONDARY: u16 = 0x100;
pub const QC_FAIL: u16 = 0x200;
pub const DUPLICATE: u16 = 0x400;
pub const SUPPLEMENTARY: u16 = 0x800;

/// A CIGAR as runs: the letter of the operation and how many.
pub type Cigar = Vec<(char, usize)>;

/// Does the operation use up bases of the read (`M I S = X`)?
pub fn uses_read(op: char) -> bool {
    matches!(op, 'M' | 'I' | 'S' | '=' | 'X')
}

/// Does the operation use up bases of the target (`M D N = X`)?
pub fn uses_target(op: char) -> bool {
    matches!(op, 'M' | 'D' | 'N' | '=' | 'X')
}

/// `*` for no CIGAR.
pub fn cigar_text(cigar: &Cigar) -> String {
    if cigar.is_empty() {
        return "*".into();
    }
    cigar.iter().map(|(op, n)| format!("{n}{op}")).collect()
}

pub fn parse_cigar(text: &str) -> Result<Cigar, String> {
    if text == "*" {
        return Ok(Vec::new());
    }
    let (mut out, mut number): (Cigar, Option<usize>) = (Vec::new(), None);
    for c in text.chars() {
        if let Some(digit) = c.to_digit(10) {
            number = Some(
                number
                    .unwrap_or(0)
                    .checked_mul(10)
                    .and_then(|n| n.checked_add(digit as usize))
                    .ok_or_else(|| format!("The CIGAR {text} has a number that is too large."))?,
            );
        } else if "MIDNSHP=X".contains(c) {
            let n = number
                .take()
                .ok_or_else(|| format!("The CIGAR {text} has {c} with no number before it."))?;
            if n > 0 {
                out.push((c, n));
            }
        } else {
            return Err(format!(
                "The CIGAR {text} has {c}, which is not an operation (M I D N S H P = X)."
            ));
        }
    }
    if number.is_some() || out.is_empty() {
        return Err(format!(
            "The CIGAR {text} is not runs of a number and an operation."
        ));
    }
    Ok(out)
}

/// One alignment, or one read that did not align.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub name: String,
    pub flag: u16,
    /// `None` for `*`.
    pub target: Option<String>,
    /// From 0 (the file counts from 1).
    pub pos: usize,
    pub mapq: u8,
    pub cigar: Cigar,
    /// Empty for `*`.
    pub seq: Vec<u8>,
    /// Phred scores; empty for `*`.
    pub qual: Vec<u8>,
    /// The tags as they are written (`NM:i:3`).
    pub tags: Vec<String>,
}

impl Record {
    /// A read that did not map: flag 4, no target, no CIGAR.
    pub fn unmapped(read: &Read) -> Record {
        Record {
            name: read.name.clone(),
            flag: UNMAPPED,
            target: None,
            pos: 0,
            mapq: 0,
            cigar: Vec::new(),
            seq: read.seq.clone(),
            qual: read.qual.clone(),
            tags: Vec::new(),
        }
    }

    /// The line of a mapping that has its base-level alignment (`-c`): the read as it lies against the target (the other strand's
    /// reverse complement, qualities turned round), the bases before the first anchor and after the last as soft clips. A secondary mapping has
    /// flag 256 and no bases (`*`), as the specification lets it. `None` for a hit that has no alignment.
    pub fn from_hit(read: &Read, index: &Index, h: &Hit) -> Option<Record> {
        let base = h.base.as_ref()?;
        let mut cigar = parse_cigar(&base.cigar).ok()?;
        let used: usize = cigar
            .iter()
            .filter(|(op, _)| uses_read(*op))
            .map(|c| c.1)
            .sum();
        let qlen = read.seq.len();
        let left = if h.rev { qlen - h.qend } else { h.qstart };
        let right = qlen.saturating_sub(left + used);
        if left > 0 {
            cigar.insert(0, ('S', left));
        }
        if right > 0 {
            cigar.push(('S', right));
        }
        let (seq, qual) = if h.rev {
            (revcom(&read.seq), read.qual.iter().rev().copied().collect())
        } else {
            (read.seq.clone(), read.qual.clone())
        };
        let mut flag = if h.rev { REVERSE } else { 0 };
        if !h.primary {
            flag |= SECONDARY;
        }
        let mut tags = vec![
            format!("NM:i:{}", base.differences),
            format!("tp:A:{}", if h.primary { 'P' } else { 'S' }),
            format!("cm:i:{}", h.anchors.len()),
            format!("s1:i:{}", h.score.round() as i64),
        ];
        if let Some(d) = h.divergence {
            tags.push(format!("dv:f:{d:.4}"));
        }
        let (seq, qual) = if h.primary {
            (seq, qual)
        } else {
            (Vec::new(), Vec::new())
        };
        Some(Record {
            name: read.name.clone(),
            flag,
            target: Some(index.name(h.target).to_string()),
            pos: h.tstart,
            mapq: h.mapq,
            cigar,
            seq,
            qual,
            tags,
        })
    }

    /// Does it lie on a target?
    pub fn mapped(&self) -> bool {
        self.flag & UNMAPPED == 0 && self.target.is_some() && !self.cigar.is_empty()
    }

    /// Does it count for a pileup: mapped, and not a secondary, a duplicate or a read that failed the sequencer's checks? (A supplementary,
    /// the other part of a split read, counts, as `bcftools mpileup` counts it by default.)
    pub fn usable(&self) -> bool {
        self.mapped() && self.flag & (SECONDARY | DUPLICATE | QC_FAIL) == 0 && !self.seq.is_empty()
    }

    pub fn line(&self) -> String {
        let seq = if self.seq.is_empty() {
            "*".to_string()
        } else {
            String::from_utf8_lossy(&self.seq).into_owned()
        };
        let qual: String = if self.qual.is_empty() {
            "*".into()
        } else {
            self.qual
                .iter()
                .map(|&q| char::from(q.min(93) + 33))
                .collect()
        };
        let mut out = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t*\t0\t0\t{seq}\t{qual}",
            self.name,
            self.flag,
            self.target.as_deref().unwrap_or("*"),
            if self.target.is_some() {
                self.pos + 1
            } else {
                0
            },
            self.mapq,
            cigar_text(&self.cigar),
        );
        for tag in &self.tags {
            out.push('\t');
            out.push_str(tag);
        }
        out
    }
}

fn parse_line(line: &str, number: usize) -> Result<Record, String> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 11 {
        return Err(format!(
            "Line {number} of the SAM has {} columns, and a line has at least 11.",
            f.len()
        ));
    }
    let bad =
        |what: &str, value: &str| format!("Line {number} of the SAM has {value} for the {what}.");
    let flag: u16 = f[1].parse().map_err(|_| bad("flag", f[1]))?;
    let pos: usize = f[3].parse().map_err(|_| bad("position", f[3]))?;
    let mapq: u8 = f[4].parse().map_err(|_| bad("mapping quality", f[4]))?;
    let cigar = parse_cigar(f[5]).map_err(|e| format!("Line {number}: {e}"))?;
    let seq: Vec<u8> = if f[9] == "*" {
        Vec::new()
    } else {
        f[9].bytes().map(|b| b.to_ascii_uppercase()).collect()
    };
    let qual: Vec<u8> = if f[10] == "*" || seq.is_empty() {
        Vec::new()
    } else if f[10].len() != seq.len() {
        return Err(format!(
            "Line {number} of the SAM has {} bases and {} qualities.",
            seq.len(),
            f[10].len()
        ));
    } else if f[10].bytes().any(|q| q < 33) {
        return Err(format!("Line {number} of the SAM has a quality below '!'."));
    } else {
        f[10].bytes().map(|q| q - 33).collect()
    };
    let in_read: usize = cigar
        .iter()
        .filter(|(op, _)| uses_read(*op))
        .map(|c| c.1)
        .sum();
    if !seq.is_empty() && !cigar.is_empty() && in_read != seq.len() {
        return Err(format!(
            "Line {number} of the SAM has a CIGAR of {in_read} bases and a read of {}.",
            seq.len()
        ));
    }
    // a position of 0 is no coordinate: such a line lies on no target
    let target = (f[2] != "*" && pos > 0).then(|| f[2].to_string());
    Ok(Record {
        name: f[0].to_string(),
        flag,
        pos: pos.saturating_sub(1),
        target,
        mapq,
        cigar,
        seq,
        qual,
        tags: f[11..].iter().map(|t| t.to_string()).collect(),
    })
}

/// A SAM file: the targets of its `@SQ` lines, and its alignments.
#[derive(Clone, Debug, PartialEq)]
pub struct Sam {
    pub targets: Vec<(String, usize)>,
    pub records: Vec<Record>,
}

fn tag_value<'a>(fields: &[&'a str], tag: &str) -> Option<&'a str> {
    fields.iter().find_map(|f| f.strip_prefix(tag))
}

/// Reads SAM text. An alignment that names a target that the `@SQ` lines (when there are any) do not hold is an error.
pub fn parse(text: &str) -> Result<Sam, String> {
    let (mut targets, mut records) = (Vec::new(), Vec::new());
    for (i, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('@') {
            if rest.starts_with("SQ") {
                let f: Vec<&str> = rest.split('\t').collect();
                let name = tag_value(&f, "SN:")
                    .ok_or_else(|| format!("Line {} (@SQ) has no SN.", i + 1))?;
                let length: usize = tag_value(&f, "LN:")
                    .and_then(|l| l.parse().ok())
                    .ok_or_else(|| format!("Line {} (@SQ) has no LN that is a number.", i + 1))?;
                targets.push((name.to_string(), length));
            }
            continue;
        }
        records.push(parse_line(line, i + 1)?);
    }
    if !targets.is_empty() {
        if let Some(r) = records.iter().find(|r| {
            r.target
                .as_ref()
                .is_some_and(|t| !targets.iter().any(|(n, _)| n == t))
        }) {
            return Err(format!(
                "The read {} names the target {}, which no @SQ line holds.",
                r.name,
                r.target.as_deref().unwrap_or("*")
            ));
        }
    }
    Ok(Sam { targets, records })
}

/// The header: version, no order, the targets, and this program.
pub fn header(targets: &[(String, usize)]) -> String {
    let mut out = String::from("@HD\tVN:1.6\tSO:unsorted\n");
    for (name, length) in targets {
        out.push_str(&format!("@SQ\tSN:{name}\tLN:{length}\n"));
    }
    out.push_str(&format!(
        "@PG\tID:helixloom\tPN:helixloom\tVN:{}\n",
        env!("CARGO_PKG_VERSION")
    ));
    out
}

pub fn write(sam: &Sam) -> String {
    let mut out = header(&sam.targets);
    for r in &sam.records {
        out.push_str(&r.line());
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests;
