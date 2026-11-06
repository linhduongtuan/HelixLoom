use super::*;
use std::collections::HashMap;

/// A made-up trace in text, and what the original's own procedures made of it (cut out of AppMain.tcl and run with
/// tclsh, `abi_save_scf` and `read_abi_info`, unchanged): its SCF export, and its reading of that export. The two
/// version 2 files were written by a small script; `golden_v2_*.read` is what `read_abi_info` made of each.
const TRACE_TEXT: &str = include_str!("golden_trace.txt");
const WRITTEN: &[u8] = include_bytes!("golden.scf");
const WRITTEN_READ: &str = include_str!("golden.read");
const V2_BYTES: &[u8] = include_bytes!("golden_v2_bytes.scf");
const V2_BYTES_READ: &str = include_str!("golden_v2_bytes.read");
const V2_SHORTS: &[u8] = include_bytes!("golden_v2_shorts.scf");
const V2_SHORTS_READ: &str = include_str!("golden_v2_shorts.read");

fn fields(text: &str) -> HashMap<&str, &str> {
    text.lines().filter_map(|l| l.split_once('=')).collect()
}

fn numbers(text: &str) -> Vec<i64> {
    text.split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().expect("a number"))
        .collect()
}

/// The trace the golden text describes: four channels in the order G, A, T, C, forty bases (one an N), peaks fourteen
/// scans apart, qualities, signal strengths and run times.
fn golden_trace() -> Trace {
    let f = fields(TRACE_TEXT);
    let channel = |k: usize| -> Vec<i16> {
        numbers(f[format!("data{k}").as_str()])
            .iter()
            .map(|&v| v as i16)
            .collect()
    };
    Trace {
        sample: String::new(),
        instrument: f["modl"].to_string(),
        order: f["wheel"].as_bytes().try_into().expect("four letters"),
        channels: std::array::from_fn(channel),
        bases: f["bases"].as_bytes().to_vec(),
        quality: numbers(f["prob"]).iter().map(|&q| q as u8).collect(),
        peaks: numbers(f["loc"]).iter().map(|&p| p as u32).collect(),
        strength: numbers(f["scale"])
            .iter()
            .map(|&s| s as i16)
            .collect::<Vec<_>>()
            .try_into()
            .expect("four strengths"),
        started: Stamp::from_texts(f["rund1"], f["runt1"]),
        ended: Stamp::from_texts(f["rund2"], f["runt2"]),
    }
}

/// What a trace must be to equal the original's reading of a file (`dump_read.tcl`'s fields).
fn assert_reads_as_the_original_does(t: &Trace, original: &str, what: &str) {
    let f = fields(original);
    assert_eq!(f["status"], "ok", "{what}");
    assert_eq!(
        String::from_utf8_lossy(&t.order),
        f["wheel"],
        "{what}: order"
    );
    assert_eq!(
        String::from_utf8_lossy(&t.bases),
        f["bases"],
        "{what}: bases"
    );
    let as_i64 = |v: &[u32]| v.iter().map(|&x| i64::from(x)).collect::<Vec<_>>();
    assert_eq!(as_i64(&t.peaks), numbers(f["loc"]), "{what}: places");
    let quality: Vec<i64> = t.quality.iter().map(|&q| i64::from(q)).collect();
    assert_eq!(quality, numbers(f["prob"]), "{what}: how sure");
    for (k, channel) in t.channels.iter().enumerate() {
        let values: Vec<i64> = channel.iter().map(|&v| i64::from(v)).collect();
        assert_eq!(
            values,
            numbers(f[format!("data{k}").as_str()]),
            "{what}: channel {k}"
        );
    }
}

fn first_difference(a: &[u8], b: &[u8]) -> Option<usize> {
    (0..a.len().max(b.len())).find(|&i| a.get(i) != b.get(i))
}

#[test]
fn the_export_is_byte_for_byte_what_the_originals_writer_makes() {
    let mine = to_scf(&golden_trace());
    assert_eq!(
        first_difference(&mine, WRITTEN),
        None,
        "the lengths are {} and {}",
        mine.len(),
        WRITTEN.len()
    );
}

#[test]
fn the_originals_own_export_reads_as_its_own_reader_reads_it() {
    let t = parse(WRITTEN).expect("a trace");
    assert_reads_as_the_original_does(&t, WRITTEN_READ, "version 3");
}

#[test]
fn version_two_files_with_one_and_two_byte_samples_read_as_the_original_reads_them() {
    for (file, original, what) in [
        (V2_BYTES, V2_BYTES_READ, "one-byte samples"),
        (V2_SHORTS, V2_SHORTS_READ, "two-byte samples"),
    ] {
        let t = parse(file).expect("a trace");
        assert_reads_as_the_original_does(&t, original, what);
        // the comments say what instrument it was
        assert_eq!(t.instrument, "310", "{what}");
    }
}

#[test]
fn a_trace_written_and_read_again_is_the_same_trace_in_the_order_of_the_file() {
    let t = golden_trace();
    let back = parse(&to_scf(&t)).expect("a trace");
    assert_eq!(back.order, *b"ACGT");
    for (k, &letter) in back.order.iter().enumerate() {
        let from = t
            .order
            .iter()
            .position(|&b| b == letter)
            .expect("a channel");
        assert_eq!(
            back.channels[k],
            t.channels[from],
            "channel of {}",
            char::from(letter)
        );
        assert_eq!(
            back.strength[k],
            t.strength[from],
            "strength of {}",
            char::from(letter)
        );
    }
    assert_eq!(back.bases, t.bases);
    assert_eq!(back.peaks, t.peaks);
    // the file keeps the chance of the called base only for A, C, G and T: an N is left with none
    let expected: Vec<u8> = t
        .bases
        .iter()
        .zip(&t.quality)
        .map(|(&b, &q)| if b == b'N' { 0 } else { q })
        .collect();
    assert_eq!(back.quality, expected);
    assert_eq!((back.started, back.ended), (t.started, t.ended));
    assert_eq!(back.instrument, "3730");
    // and writing it again changes nothing
    assert_eq!(
        to_scf(&back),
        to_scf(&parse(&to_scf(&back)).expect("a trace"))
    );
}

#[test]
fn the_comments_put_a_runs_time_on_a_line_of_its_own_as_the_original_does() {
    // the header's eighth number after the magic is the comments' place
    let at = be_u32(WRITTEN, 32).expect("the comments' place") as usize;
    assert_eq!(
        String::from_utf8_lossy(&WRITTEN[at..]),
        "STRT=2/19/2009\n1:19:30\nSTOP=2/19/2009\n4:04:15\nSIGN=G=93,A=134,T=134,C=52\nMODL=3730"
    );
    // a run with no dates gets the original's defaults
    let bare = Trace::default();
    let text = comments_text(&bare);
    assert!(
        text.starts_with("STRT=01/01/00\n00:00:00\nSTOP=01/01/00\n00:00:00\n"),
        "{text}"
    );
    assert!(!text.contains("MODL"));
}

#[test]
fn the_comments_are_read_whether_the_time_follows_on_the_line_or_on_the_next() {
    let own_line = notes("STRT=2/19/2009\n1:19:30\nSTOP=3/1/2010\n23:05:09\nSIGN=A=1,C=2,G=3,T=4\nMODL=ABI PRISM 310");
    let same_line = notes("STRT=2/19/2009 1:19:30\nSTOP=3/1/2010 23:05:09\nSIGN=A=1,C=2,G=3,T=4\nMODL=ABI PRISM 310\0");
    for n in [own_line, same_line] {
        let begin = n.started.expect("a start");
        assert_eq!((begin.year, begin.month, begin.day), (2009, 2, 19));
        assert_eq!((begin.hour, begin.minute, begin.second), (1, 19, 30));
        assert_eq!(n.ended.expect("an end").hour, 23);
        assert_eq!(n.strength, [1, 2, 3, 4]);
        assert_eq!(n.instrument, "ABI PRISM 310");
    }
    // a date alone, or a time that is not one, leaves the moment unknown
    assert!(notes("STRT=2/19/2009").started.is_none());
    assert!(notes("STRT=2/19/2009\nnot a time").started.is_none());
}

#[test]
fn the_version_is_compared_as_the_original_does() {
    assert_eq!(version(b"3.00"), 3.0);
    assert_eq!(version(b"3.10"), 3.1);
    assert_eq!(version(b"2.90"), 2.9);
    assert_eq!(version(b"2.00"), 2.0);
    // text that is not a number reads as version 3
    assert_eq!(version(b"v3?!"), 3.0);
}

#[test]
fn signals_that_jump_further_than_sixteen_bits_still_come_back_exactly() {
    let t = Trace {
        order: *b"ACGT",
        channels: [
            vec![0, 30000, -30000, 0, 32767, -32768, 5],
            vec![1; 7],
            vec![-1; 7],
            vec![0; 7],
        ],
        bases: b"AC".to_vec(),
        quality: vec![40, 50],
        peaks: vec![2, 5],
        ..Trace::default()
    };
    let back = parse(&to_scf(&t)).expect("a trace");
    assert_eq!(back.channels, t.channels);
}

#[test]
fn an_empty_trace_is_written_and_read_without_trouble() {
    let t = Trace::default();
    let back = parse(&to_scf(&t)).expect("a trace");
    assert!(back.bases.is_empty() && back.peaks.is_empty() && back.channels[0].is_empty());
}

#[test]
fn damaged_files_make_an_error_or_a_trace_and_never_a_failure() {
    // cut at every length
    for len in (0..WRITTEN.len()).step_by(11) {
        let _ = parse(&WRITTEN[..len]);
    }
    // every byte of the header scrambled in turn
    for k in 0..WRITTEN_HEADER {
        for value in [0u8, 0x7f, 0xff] {
            let mut bad = WRITTEN.to_vec();
            bad[k] = value;
            let _ = parse(&bad);
        }
    }
    // counts and places far beyond the file
    for field in [4usize, 8, 12, 24, 28, 32] {
        let mut bad = WRITTEN.to_vec();
        bad[field..field + 4].copy_from_slice(&u32::MAX.to_be_bytes());
        let _ = parse(&bad);
    }
    assert!(parse(b"").is_err());
    assert!(parse(b".scf").is_err());
    assert!(
        parse(b"ABIF and not a trace at all, long enough for a header, so it cannot be one")
            .unwrap_err()
            .contains("not an SCF")
    );
}

#[test]
fn either_kind_of_file_is_told_by_its_first_four_bytes() {
    use crate::seq::abif::read;
    assert!(read(WRITTEN).is_ok());
    assert!(read(V2_BYTES).is_ok());
    assert!(read(b"FASTA text")
        .unwrap_err()
        .contains("neither as an ABIF"));
    assert!(read(b"").is_err());
}
