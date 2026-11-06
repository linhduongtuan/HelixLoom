//! The original's ORFs → Find Next / Find Previous (`next_ORF`, AppMain.tcl line 15025, file offset 0xa9b56)
//! and the list of ORFs it walks (`orf_list`, line 15100, offset 0xaa4a3), quirks included. Tested against
//! the original's own code (`nextorf/golden.txt`).

use super::dna::{revcom, StartRule, Strand};

/// The original reads the position of a stop codon before it has found any (a sequence with no stop
/// codon at all) and stops with a Tcl error.
#[derive(Debug, PartialEq, Eq)]
pub struct NoStopCodon;

/// One entry of the original's list: the first and last base of the ORF on the top strand (counting from 0,
/// both included; an ORF with no stop in frame ends at `len`) and its strand.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Entry {
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
}

/// The ORF menu's settings: topology, strand, shortest ORF and what an ORF starts with.
#[derive(Clone, Copy, Debug)]
pub struct Settings {
    pub circular: bool,
    pub strand: Strand,
    pub min: usize,
    pub rule: StartRule,
}

/// Where the next ORF is and what the original selects for it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Found {
    pub entry: Entry,
    /// First base and one past the last base to select (may run one past the end of the sequence).
    pub select: (usize, usize),
}

/// Start positions of every match of `pattern` (upper case) in `text`, left to right without overlap.
fn matches(text: &[u8], pattern: &[&[u8]]) -> Vec<i64> {
    let mut found = Vec::new();
    let mut p = 0;
    while p + 3 <= text.len() {
        let here = text[p..p + 3].to_ascii_uppercase();
        if pattern.contains(&here.as_slice()) {
            found.push(p as i64);
            p += 3;
        } else {
            p += 1;
        }
    }
    found
}

const STOPS: [&[u8]; 3] = [b"TGA", b"TAA", b"TAG"];

/// The original's two lists for one strand: where ORFs may begin and where stop codons are. With the
/// "Stop+1" rule an ORF begins right after any stop codon.
fn starts_and_stops(text: &[u8], rule: StartRule) -> (Vec<i64>, Vec<i64>) {
    let stops = matches(text, &STOPS);
    let starts = match rule {
        StartRule::Met => matches(text, &[b"ATG"]),
        StartRule::AfterStop => stops.iter().map(|s| s + 3).collect(),
    };
    (starts, stops)
}

/// The original's `orf_list`. The stop position `c` it tests keeps its value from one start to the next, and
/// from the top strand to the bottom one, so a start with no stop in its frame is judged by whichever stop it
/// looked at last: with the top strand it then ends at the end of the sequence, with the bottom strand it
/// is listed as if that stop were its own.
pub fn orf_list(text: &[u8], settings: &Settings) -> Result<Vec<Entry>, NoStopCodon> {
    let Settings {
        circular,
        strand,
        min,
        rule,
    } = *settings;
    let end = text.len() as i64;
    if end == 0 {
        return Ok(Vec::new());
    }
    let min = min as i64;
    let both_copies = if circular {
        [text, text].concat()
    } else {
        text.to_vec()
    };
    let mut result = Vec::new();
    let mut c: Option<i64> = None;
    let mut strand_pass = |text: &[u8], reverse: bool, result: &mut Vec<Entry>| {
        let (starts, stops) = starts_and_stops(text, rule);
        for &b in starts.iter().take_while(|&&b| b <= end) {
            let frame = b % 3;
            for &stop in &stops {
                c = Some(stop);
                if stop >= b && stop % 3 == frame {
                    break;
                }
            }
            let Some(c) = c else {
                return Err(NoStopCodon);
            };
            if c + 2 - b <= min {
                continue;
            }
            result.push(if reverse {
                Entry {
                    start: (end - c - 3).rem_euclid(end) as usize,
                    end: (end - b - 1).rem_euclid(end) as usize,
                    reverse: true,
                }
            } else if c % 3 == frame {
                Entry {
                    start: (b % end) as usize,
                    end: ((c + 2) % end) as usize,
                    reverse: false,
                }
            } else {
                Entry {
                    start: (b % end) as usize,
                    end: end as usize,
                    reverse: false,
                }
            });
        }
        Ok(())
    };
    if strand != Strand::Bottom {
        strand_pass(&both_copies, false, &mut result)?;
    }
    if strand != Strand::Top {
        strand_pass(&revcom(&both_copies), true, &mut result)?;
    }
    Ok(result)
}

/// What Find Next ORF (`forward`) or Find Previous ORF selects, starting from the selection's first base or,
/// with nothing selected, from the caret. The list is taken in order of start; going forward from a selection
/// looks from one base after it, going back from one base before it.
pub fn find_next(
    text: &[u8],
    settings: &Settings,
    caret: usize,
    selection_start: Option<usize>,
    forward: bool,
) -> Result<Option<Found>, NoStopCodon> {
    let len = text.len() as i64;
    if len == 0 {
        return Ok(None);
    }
    let from = match selection_start {
        Some(s) if forward => (s as i64 + 1).rem_euclid(len),
        Some(s) => (s as i64 - 1).rem_euclid(len),
        None => caret as i64,
    };
    let mut list = orf_list(text, settings)?;
    list.sort_by_key(|e| e.start);
    let ahead = list.iter().take_while(|e| from > e.start as i64).count() as i64;
    let at = if forward { ahead } else { ahead - 1 };
    let Some(entry) = list.get(at.rem_euclid(list.len().max(1) as i64) as usize) else {
        return Ok(None);
    };
    // an ORF across the origin is only selected up to the end of the sequence
    let last = if entry.start > entry.end {
        len as usize
    } else {
        entry.end
    };
    Ok(Some(Found {
        entry: *entry,
        select: (entry.start, last + 1),
    }))
}

#[cfg(test)]
mod tests;
