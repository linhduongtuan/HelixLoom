//! SCF files (`.scf`), the second kind of trace file the original opens, read and written by the rules of its own
//! procedures: `read_abi_info` (AppMain.tcl:30572, 0x1684b9) reads versions 2 and 3, and `abi_save_scf`
//! (AppMain.tcl:31933, 0x177560) writes version 3.00. All numbers are big-endian.
//!
//! A version 3 file keeps the signal of each channel, one channel after another, as the change of the change from one
//! scan to the next; a version 2 file keeps the four channels together, scan by scan, as they are. The channels are
//! always in the order A, C, G, T.

use super::abif::{Stamp, Trace};

const MAGIC: &[u8; 4] = b".scf";
const ORDER: [u8; 4] = *b"ACGT";
/// What the original writes: a header of 128 bytes, two-byte samples and code set 2.
const WRITTEN_HEADER: usize = 128;
const WRITTEN_SAMPLE_BYTES: usize = 2;
const WRITTEN_CODE_SET: usize = 2;
/// In version 3, for each base: where it is (4 bytes), how sure the file is that it is each of A, C, G and T (one byte
/// each), the base itself and three spare bytes. A version 2 record is the same twelve bytes in a different order.
const BASE_BYTES: usize = 12;
const CUT: &str = "The trace file is cut short: its sections run past the end.";

fn be_u32(bytes: &[u8], at: usize) -> Option<u32> {
    let end = at.checked_add(4)?;
    Some(u32::from_be_bytes(bytes.get(at..end)?.try_into().ok()?))
}

/// The header fields the original reads, in the order it scans them (`IIIIIIIIa4IIII` after the magic).
struct Header {
    samples: usize,
    samples_offset: usize,
    bases: usize,
    bases_offset: usize,
    comments_size: usize,
    comments_offset: usize,
    version: f64,
    /// 1 or 2: the original reads anything but 1 as two-byte samples.
    sample_bytes: usize,
}

/// The version as the original compares it (`$version >= 3.00`): a number, and for text that is not one, version 3.
fn version(text: &[u8]) -> f64 {
    std::str::from_utf8(text)
        .ok()
        .and_then(|t| t.trim().parse().ok())
        .unwrap_or(3.0)
}

fn header(bytes: &[u8]) -> Result<Header, String> {
    // after the magic: samples, where they start, bases, two clip points, where the bases start, the size and place of
    // the comments, the version (four characters), the sample size, the code set, the size and place of private data
    let field = |k: usize| {
        be_u32(bytes, 4 + 4 * k)
            .and_then(|v| usize::try_from(v).ok())
            .ok_or(CUT)
    };
    Ok(Header {
        samples: field(0)?,
        samples_offset: field(1)?,
        bases: field(2)?,
        bases_offset: field(5)?,
        comments_size: field(6)?,
        comments_offset: field(7)?,
        version: version(bytes.get(36..40).ok_or(CUT)?),
        sample_bytes: if field(9)? == 1 { 1 } else { 2 },
    })
}

/// `count` items of `size` bytes from `at`, or the error for a file that is cut short.
fn section(bytes: &[u8], at: usize, count: usize, size: usize) -> Result<&[u8], String> {
    let end = count
        .checked_mul(size)
        .and_then(|n| n.checked_add(at))
        .ok_or(CUT)?;
    bytes.get(at..end).ok_or_else(|| CUT.to_string())
}

/// Sample `i` of `raw`: a signed byte or a signed two-byte number.
fn sample(raw: &[u8], i: usize, size: usize) -> i16 {
    if size == 1 {
        i16::from(raw[i] as i8)
    } else {
        i16::from_be_bytes([raw[2 * i], raw[2 * i + 1]])
    }
}

/// Version 3: each channel in turn, every sample the change from the one before of the change before that, so two
/// running sums give the signal (wrapping, as the sixteen-bit numbers they were made from).
fn channels_v3(bytes: &[u8], h: &Header) -> Result<[Vec<i16>; 4], String> {
    let size = h.sample_bytes;
    let raw = section(
        bytes,
        h.samples_offset,
        h.samples.checked_mul(4).ok_or(CUT)?,
        size,
    )?;
    Ok(std::array::from_fn(|k| {
        let channel = &raw[k * h.samples * size..];
        let (mut change, mut level) = (0i16, 0i16);
        (0..h.samples)
            .map(|i| {
                change = change.wrapping_add(sample(channel, i, size));
                level = level.wrapping_add(change);
                level
            })
            .collect()
    }))
}

/// Version 2: for each scan the four channels, as they are.
fn channels_v2(bytes: &[u8], h: &Header) -> Result<[Vec<i16>; 4], String> {
    let size = h.sample_bytes;
    let raw = section(
        bytes,
        h.samples_offset,
        h.samples.checked_mul(4).ok_or(CUT)?,
        size,
    )?;
    Ok(std::array::from_fn(|k| {
        (0..h.samples)
            .map(|i| sample(raw, i * 4 + k, size))
            .collect()
    }))
}

/// The base calls of a file: where each is, the base, and how sure the file is of it.
struct Calls {
    peaks: Vec<u32>,
    bases: Vec<u8>,
    quality: Vec<u8>,
}

/// How sure the file is of a base: the chance it gives for that letter, and none for a letter that is not A, C, G or T.
fn sureness(letter: u8, chance: impl Fn(usize) -> u8) -> u8 {
    ORDER
        .iter()
        .position(|&b| b == letter.to_ascii_uppercase())
        .map_or(0, chance)
}

/// Version 3: where the bases are, the four arrays of chances, then the bases.
fn calls_v3(bytes: &[u8], h: &Header) -> Result<Calls, String> {
    let n = h.bases;
    let places = section(bytes, h.bases_offset, n, 4)?;
    let after_places = h
        .bases_offset
        .checked_add(n.checked_mul(4).ok_or(CUT)?)
        .ok_or(CUT)?;
    let chances = section(bytes, after_places, n, 4)?;
    let letters = section(
        bytes,
        after_places
            .checked_add(n.checked_mul(4).ok_or(CUT)?)
            .ok_or(CUT)?,
        n,
        1,
    )?;
    Ok(Calls {
        peaks: places
            .as_chunks::<4>()
            .0
            .iter()
            .map(|p| u32::from_be_bytes(*p))
            .collect(),
        bases: letters.to_vec(),
        quality: letters
            .iter()
            .enumerate()
            .map(|(i, &b)| sureness(b, |k| chances[k * n + i]))
            .collect(),
    })
}

/// Version 2: a record of twelve bytes for each base, right after the samples.
fn calls_v2(bytes: &[u8], h: &Header) -> Result<Calls, String> {
    let start = h
        .samples
        .checked_mul(4 * h.sample_bytes)
        .and_then(|n| n.checked_add(h.samples_offset))
        .ok_or(CUT)?;
    let records = section(bytes, start, h.bases, BASE_BYTES)?;
    let mut calls = Calls {
        peaks: Vec::new(),
        bases: Vec::new(),
        quality: Vec::new(),
    };
    for r in records.as_chunks::<BASE_BYTES>().0 {
        calls
            .peaks
            .push(u32::from_be_bytes([r[0], r[1], r[2], r[3]]));
        calls.bases.push(r[8]);
        calls.quality.push(sureness(r[8], |k| r[4 + k]));
    }
    Ok(calls)
}

/// What the comments say that a trace keeps: the instrument, when the run began and ended, and the signal strengths.
#[derive(Default)]
struct Notes {
    instrument: String,
    started: Option<Stamp>,
    ended: Option<Stamp>,
    strength: [i16; 4],
}

/// The comments are lines of `KEY=value`. The original writes a run's date and its time on two lines (a slip of its
/// list handling), so a time may sit on the line after the date as well as after it on the same line.
fn notes(comments: &str) -> Notes {
    let lines: Vec<&str> = comments
        .lines()
        .map(|l| l.trim_matches(|c: char| c == '\0' || c.is_whitespace()))
        .collect();
    let mut notes = Notes::default();
    for (i, line) in lines.iter().enumerate() {
        let when = |rest: &str| match rest.split_once(' ') {
            Some((date, time)) => Stamp::from_texts(date, time.trim()),
            None => Stamp::from_texts(rest, lines.get(i + 1).copied().unwrap_or("")),
        };
        if let Some(rest) = line.strip_prefix("STRT=") {
            notes.started = when(rest);
        } else if let Some(rest) = line.strip_prefix("STOP=") {
            notes.ended = when(rest);
        } else if let Some(model) = line.strip_prefix("MODL=") {
            notes.instrument = model.trim().to_string();
        } else if let Some(list) = line.strip_prefix("SIGN=") {
            for part in list.split(',') {
                let Some((letter, value)) = part.split_once('=') else {
                    continue;
                };
                let slot = ORDER.iter().position(|&b| letter.as_bytes() == [b]);
                if let (Some(k), Ok(v)) = (slot, value.trim().parse()) {
                    notes.strength[k] = v;
                }
            }
        }
    }
    notes
}

/// Read a trace from the bytes of an `.scf` file.
pub fn parse(bytes: &[u8]) -> Result<Trace, String> {
    if bytes.get(..4) != Some(MAGIC.as_slice()) {
        return Err("This is not an SCF file: it does not start with \".scf\".".into());
    }
    let h = header(bytes)?;
    let (channels, calls) = if h.version >= 3.0 {
        (channels_v3(bytes, &h)?, calls_v3(bytes, &h)?)
    } else {
        (channels_v2(bytes, &h)?, calls_v2(bytes, &h)?)
    };
    let comments = bytes
        .get(h.comments_offset..h.comments_offset.saturating_add(h.comments_size))
        .map(|c| String::from_utf8_lossy(c).into_owned())
        .unwrap_or_default();
    let notes = notes(&comments);
    Ok(Trace {
        sample: String::new(),
        instrument: notes.instrument,
        order: ORDER,
        channels,
        bases: calls.bases,
        quality: calls.quality,
        peaks: calls.peaks,
        strength: notes.strength,
        started: notes.started,
        ended: notes.ended,
    })
}

/// A run's date and time as the original puts them in the comments, "01/01/00" and "00:00:00" when it has none.
fn stamp_texts(stamp: Option<Stamp>) -> (String, String) {
    stamp.map_or(("01/01/00".into(), "00:00:00".into()), |s| {
        (s.date_text(), s.time_text())
    })
}

/// The comments `abi_save_scf` writes: the run's start and end, each with its time on a line of its own, the strength of
/// each channel in the trace's own order, and the instrument.
fn comments_text(t: &Trace) -> String {
    let ((start_date, start_time), (end_date, end_time)) =
        (stamp_texts(t.started), stamp_texts(t.ended));
    let strengths: Vec<String> = t
        .order
        .iter()
        .zip(t.strength)
        .map(|(&b, s)| format!("{}={s}", char::from(b)))
        .collect();
    let mut lines = vec![
        format!("STRT={start_date}"),
        start_time,
        format!("STOP={end_date}"),
        end_time,
        format!("SIGN={}", strengths.join(",")),
    ];
    if !t.instrument.is_empty() {
        lines.push(format!("MODL={}", t.instrument));
    }
    lines.join("\n")
}

fn put(out: &mut Vec<u8>, value: usize) {
    out.extend((value as u32).to_be_bytes());
}

/// The trace as the bytes of an SCF version 3.00 file, byte for byte as the original's "Export as SCF" writes it:
/// the channels in the order A, C, G, T, each as the change of the change, the base calls with the chance of the called
/// base under its own letter and nothing under the others, and the comments.
pub fn to_scf(t: &Trace) -> Vec<u8> {
    let n = t.bases.len().min(t.peaks.len());
    let samples = t.channels[0].len();
    let comments = comments_text(t);
    let bases_offset = WRITTEN_HEADER + samples * 4 * WRITTEN_SAMPLE_BYTES;
    let comments_offset = bases_offset + n * BASE_BYTES;
    let mut out = MAGIC.to_vec();
    for value in [
        samples,
        WRITTEN_HEADER,
        n,
        0,
        0,
        bases_offset,
        comments.len(),
        comments_offset,
    ] {
        put(&mut out, value);
    }
    out.extend(b"3.00");
    for value in [
        WRITTEN_SAMPLE_BYTES,
        WRITTEN_CODE_SET,
        0,
        comments_offset + comments.len(),
    ] {
        put(&mut out, value);
    }
    out.resize(WRITTEN_HEADER, 0);
    for letter in ORDER {
        let channel = t
            .order
            .iter()
            .position(|&b| b == letter)
            .map_or(&[][..], |k| &t.channels[k][..]);
        let (mut before, mut change_before) = (0i32, 0i32);
        for i in 0..samples {
            let value = i32::from(channel.get(i).copied().unwrap_or(0));
            let change = value - before;
            // sixteen bits: what does not fit wraps, as the original's `binary format S` does
            out.extend(((change - change_before) as i16).to_be_bytes());
            (before, change_before) = (value, change);
        }
    }
    for &peak in &t.peaks[..n] {
        out.extend(peak.to_be_bytes());
    }
    for letter in ORDER {
        out.extend((0..n).map(|i| {
            if t.bases[i].to_ascii_uppercase() == letter {
                t.quality.get(i).copied().unwrap_or(0)
            } else {
                0
            }
        }));
    }
    out.extend(&t.bases[..n]);
    out.resize(out.len() + n * 3, 0);
    out.extend(comments.as_bytes());
    out
}

#[cfg(test)]
mod tests;
