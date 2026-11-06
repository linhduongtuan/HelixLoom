use super::*;

/// One block of data as a file stores it: its name and number, its element type and size, and its bytes.
struct Block {
    name: [u8; 4],
    number: i32,
    kind: i16,
    element_size: i16,
    count: i32,
    data: Vec<u8>,
}

fn shorts(name: &[u8; 4], number: i32, values: &[i16]) -> Block {
    Block {
        name: *name,
        number,
        kind: SHORT,
        element_size: 2,
        count: values.len() as i32,
        data: values.iter().flat_map(|v| v.to_be_bytes()).collect(),
    }
}

fn bytes(name: &[u8; 4], number: i32, values: &[u8]) -> Block {
    Block {
        name: *name,
        number,
        kind: BYTE,
        element_size: 1,
        count: values.len() as i32,
        data: values.to_vec(),
    }
}

fn chars(name: &[u8; 4], number: i32, text: &str) -> Block {
    Block {
        name: *name,
        number,
        kind: CHAR,
        element_size: 1,
        count: text.len() as i32,
        data: text.as_bytes().to_vec(),
    }
}

fn pstring(name: &[u8; 4], number: i32, text: &str) -> Block {
    let mut data = vec![text.len() as u8];
    data.extend(text.as_bytes());
    Block {
        name: *name,
        number,
        kind: PSTRING,
        element_size: 1,
        count: data.len() as i32,
        data,
    }
}

/// An ABIF file holding `blocks`, laid out the way the sequencers write them.
fn file(blocks: &[Block]) -> Vec<u8> {
    let mut out = b"ABIF".to_vec();
    out.extend(101u16.to_be_bytes());
    out.extend([0u8; ENTRY]);
    let mut entries = Vec::new();
    for b in blocks {
        let mut e = Vec::new();
        e.extend(b.name);
        e.extend(b.number.to_be_bytes());
        e.extend(b.kind.to_be_bytes());
        e.extend(b.element_size.to_be_bytes());
        e.extend(b.count.to_be_bytes());
        e.extend((b.data.len() as i32).to_be_bytes());
        if b.data.len() <= INLINE {
            let mut raw = [0u8; 4];
            raw[..b.data.len()].copy_from_slice(&b.data);
            e.extend(raw);
        } else {
            e.extend((out.len() as i32).to_be_bytes());
            out.extend(&b.data);
        }
        e.extend(0i32.to_be_bytes());
        entries.push(e);
    }
    let at = out.len();
    out.extend(entries.concat());
    let mut head = b"tdir".to_vec();
    head.extend(1i32.to_be_bytes());
    head.extend(1023i16.to_be_bytes());
    head.extend((ENTRY as i16).to_be_bytes());
    head.extend((blocks.len() as i32).to_be_bytes());
    head.extend(((blocks.len() * ENTRY) as i32).to_be_bytes());
    head.extend((at as i32).to_be_bytes());
    head.extend(0i32.to_be_bytes());
    out[HEADER..HEADER + ENTRY].copy_from_slice(&head);
    out
}

/// A small run: six bases, 30 scans of four channels with a peak under each base.
fn run(extra: Vec<Block>) -> Vec<u8> {
    let channel = |k: usize| -> Vec<i16> { (0..30).map(|s| (s * 10 + k as i32) as i16).collect() };
    let mut blocks = vec![
        chars(b"FWO_", 1, "GATC"),
        pstring(b"SMPL", 1, "sample 1"),
        chars(b"MODL", 1, "3730"),
        chars(b"PBAS", 1, "GATTAC"),
        bytes(b"PCON", 1, &[10, 20, 30, 40, 50, 60]),
        shorts(b"PLOC", 1, &[3, 8, 14, 19, 24, 28]),
    ];
    blocks.extend((0..4).map(|k| shorts(b"DATA", 9 + k as i32, &channel(k))));
    blocks.extend(extra);
    file(&blocks)
}

#[test]
fn a_trace_is_read_with_its_bases_qualities_peaks_and_channels() {
    let t = parse(&run(Vec::new())).expect("a trace");
    assert_eq!(t.sample, "sample 1");
    assert_eq!(t.instrument, "3730");
    assert_eq!(t.order, *b"GATC");
    assert_eq!(t.bases, b"GATTAC");
    assert_eq!(t.quality, vec![10, 20, 30, 40, 50, 60]);
    assert_eq!(t.peaks, vec![3, 8, 14, 19, 24, 28]);
    assert_eq!(t.scans(), 30);
    assert_eq!(t.channels[2][3], 32);
    assert_eq!(t.tallest(), 293);
    assert_eq!(t.channel_of(b'a'), Some(1));
    assert_eq!(t.channel_of(b'N'), None);
}

#[test]
fn the_calls_of_set_one_are_the_ones_shown_as_in_the_original() {
    // the original shows PBAS, PLOC and PCON number 1 (it calls them "edited") and keeps number 2 apart
    let second = vec![
        chars(b"PBAS", 2, "GATNAC"),
        bytes(b"PCON", 2, &[1, 2, 3, 4, 5, 6]),
        shorts(b"PLOC", 2, &[4, 9, 15, 20, 25, 29]),
    ];
    let t = parse(&run(second)).expect("a trace");
    assert_eq!(t.bases, b"GATTAC");
    assert_eq!(t.quality, vec![10, 20, 30, 40, 50, 60]);
    assert_eq!(t.peaks, vec![3, 8, 14, 19, 24, 28]);
}

#[test]
fn set_two_is_used_only_when_a_file_has_no_set_one() {
    let blocks = vec![
        chars(b"PBAS", 2, "ACG"),
        bytes(b"PCON", 2, &[7, 8, 9]),
        shorts(b"PLOC", 2, &[1, 5, 9]),
        shorts(b"DATA", 9, &[1, 2, 3]),
        shorts(b"DATA", 10, &[1, 2, 3]),
        shorts(b"DATA", 11, &[1, 2, 3]),
        shorts(b"DATA", 12, &[1, 2, 3]),
    ];
    let t = parse(&file(&blocks)).expect("a trace");
    assert_eq!(t.bases, b"ACG");
    assert_eq!(t.quality, vec![7, 8, 9]);
    assert_eq!(t.peaks, vec![1, 5, 9]);
}

/// A date block (`RUND`: the year in two bytes, the month and the day) or a time block (`RUNT`: hour, minute, second
/// and hundredths), which sit in the entry itself.
fn stamp_block(name: &[u8; 4], number: i32, data: [u8; 4]) -> Block {
    Block {
        name: *name,
        number,
        kind: if name == b"RUND" { 10 } else { 11 },
        element_size: 4,
        count: 1,
        data: data.to_vec(),
    }
}

#[test]
fn run_dates_times_and_signal_strengths_are_read_as_the_original_reads_them() {
    let year = 2009u16.to_be_bytes();
    let blocks = vec![
        stamp_block(b"RUND", 1, [year[0], year[1], 2, 19]),
        stamp_block(b"RUNT", 1, [1, 19, 30, 0]),
        stamp_block(b"RUND", 2, [year[0], year[1], 2, 19]),
        stamp_block(b"RUNT", 2, [4, 4, 15, 0]),
        shorts(b"S/N%", 1, &[93, 134, 134, 52]),
    ];
    let t = parse(&run(blocks)).expect("a trace");
    let (begin, end) = (t.started.expect("a start"), t.ended.expect("an end"));
    assert_eq!(
        (begin.date_text(), begin.time_text()),
        ("2/19/2009".into(), "1:19:30".to_string())
    );
    assert_eq!(
        (end.date_text(), end.time_text()),
        ("2/19/2009".into(), "4:04:15".to_string())
    );
    assert_eq!(t.strength, [93, 134, 134, 52]);
    // a file with none of them has none, and a run with a date only keeps the original's time of day
    let plain = parse(&run(Vec::new())).expect("a trace");
    assert_eq!(
        (plain.started, plain.ended, plain.strength),
        (None, None, [0; 4])
    );
    let date_only = parse(&run(vec![stamp_block(
        b"RUND",
        1,
        [year[0], year[1], 12, 31],
    )]))
    .expect("a trace");
    assert_eq!(date_only.started.expect("a start").time_text(), "0:00:00");
}

#[test]
fn a_moment_is_written_in_words_as_the_originals_file_info_does() {
    let stamp = Stamp::from_texts("2/19/2009", "1:19:30").expect("a moment");
    // Thursday, as `clock format` gives it
    assert_eq!(stamp.long_text(), "1:19AM Thu, Feb 19, 2009");
    let late = Stamp::from_texts("12/31/1999", "23:05:09").expect("a moment");
    assert_eq!(late.long_text(), "11:05PM Fri, Dec 31, 1999");
    let noon = Stamp::from_texts("1/1/2000", "12:00:00").expect("a moment");
    assert_eq!(noon.long_text(), "12:00PM Sat, Jan 01, 2000");
    let midnight = Stamp::from_texts("01/01/00", "00:00:00").expect("the default");
    assert_eq!(midnight, Stamp::default());
    assert_eq!(midnight.long_text(), "12:00AM Sat, Jan 01, 2000");
    // things that are not moments
    assert!(Stamp::from_texts("13/1/2000", "1:00:00").is_none());
    assert!(Stamp::from_texts("1/1/2000", "1:00").is_none());
    assert!(Stamp::from_texts("a/b/c", "1:00:00").is_none());
}

#[test]
fn data_of_four_bytes_or_less_is_read_from_the_entry_itself() {
    // three bases fit in the entry, as do the order and nothing else
    let t = parse(&file(&[
        chars(b"FWO_", 1, "GATC"),
        chars(b"PBAS", 1, "AC"),
        bytes(b"PCON", 1, &[30, 40]),
        shorts(b"PLOC", 1, &[1, 2]),
        shorts(b"DATA", 9, &[1, 2, 3]),
        shorts(b"DATA", 10, &[1, 2, 3]),
        shorts(b"DATA", 11, &[1, 2, 3]),
        shorts(b"DATA", 12, &[1, 2, 3]),
    ]))
    .expect("a trace");
    assert_eq!(t.bases, b"AC");
    assert_eq!(t.quality, vec![30, 40]);
    assert_eq!(t.order, *b"GATC");
}

#[test]
fn the_raw_signal_is_used_when_there_is_no_analysed_one() {
    let raw: Vec<Block> = (0..4).map(|k| shorts(b"DATA", 1 + k, &[5, 6, 7])).collect();
    let mut blocks = vec![
        chars(b"PBAS", 1, "ACG"),
        bytes(b"PCON", 1, &[9, 9, 9]),
        shorts(b"PLOC", 1, &[0, 1, 2]),
    ];
    blocks.extend(raw);
    let t = parse(&file(&blocks)).expect("a trace");
    assert_eq!(t.channels[0], vec![5, 6, 7]);
    // with no order in the file the original's default is assumed
    assert_eq!(t.order, *b"ACGT");
}

#[test]
fn a_file_without_quality_values_reads_as_all_zero() {
    let blocks = vec![
        chars(b"PBAS", 1, "ACGTT"),
        shorts(b"PLOC", 1, &[1, 2, 3, 4, 5]),
        shorts(b"DATA", 9, &[1]),
        shorts(b"DATA", 10, &[1]),
        shorts(b"DATA", 11, &[1]),
        shorts(b"DATA", 12, &[1]),
    ];
    let t = parse(&file(&blocks)).expect("a trace");
    assert_eq!(t.quality, vec![0; 5]);
    assert!(!t.has_quality());
    assert_eq!(t.good_stretch(0.05), None);
}

#[test]
fn files_that_are_not_traces_are_refused_in_plain_words() {
    assert!(parse(b"").unwrap_err().contains("not an ABIF"));
    assert!(
        parse(b"FASTA text, not a trace at all, long enough to hold a header....")
            .unwrap_err()
            .contains("not an ABIF")
    );
    let good = run(Vec::new());
    // cut short in the middle of the directory
    assert!(parse(&good[..good.len() - 40])
        .unwrap_err()
        .contains("cut short"));
    // a directory with nothing of a sequencing run in it
    let bare = file(&[chars(b"SMPL", 1, "nothing here")]);
    assert!(parse(&bare).unwrap_err().contains("no base calls"));
    let no_signal = file(&[chars(b"PBAS", 1, "ACGT"), bytes(b"PCON", 1, &[1, 2, 3, 4])]);
    assert!(parse(&no_signal).unwrap_err().contains("no signal"));
}

#[test]
fn damaged_pointers_do_not_make_the_reader_fail() {
    // every byte of the directory scrambled in turn: an error or a trace, never a panic
    let good = run(Vec::new());
    let at = be_i32(&good[HEADER + 20..]) as usize;
    for k in at..good.len() {
        for value in [0u8, 0xff] {
            let mut bad = good.clone();
            bad[k] = value;
            let _ = parse(&bad);
        }
    }
    // a pointer far past the end of the file
    let mut bad = good.clone();
    let first = at + 20;
    bad[first..first + 4].copy_from_slice(&i32::MAX.to_be_bytes());
    assert!(parse(&bad).is_err() || parse(&bad).is_ok());
}

#[test]
fn the_good_stretch_is_the_run_of_bases_most_likely_right() {
    let with = |quality: &[u8]| Trace {
        order: *b"GATC",
        bases: vec![b'A'; quality.len()],
        quality: quality.to_vec(),
        peaks: (0..quality.len() as u32).collect(),
        ..Trace::default()
    };
    // a ragged start and end around a solid middle
    let t = with(&[3, 5, 2, 30, 40, 40, 38, 35, 40, 12, 4, 2, 3]);
    assert_eq!(t.good_stretch(0.05), Some((3, 9)));
    // everything good keeps everything; everything poor keeps nothing
    assert_eq!(with(&[40; 10]).good_stretch(0.05), Some((0, 10)));
    assert_eq!(with(&[2; 10]).good_stretch(0.05), None);
    assert_eq!(with(&[]).good_stretch(0.05), None);
    // a mediocre base (Q15 is still more likely right than the 0.05 cutoff) does not split a good stretch ...
    assert_eq!(
        with(&[40, 40, 40, 15, 40, 40, 40]).good_stretch(0.05),
        Some((0, 7))
    );
    // ... but a base that is probably wrong (Q5) does: the first of two equal runs wins
    assert_eq!(
        with(&[40, 40, 40, 5, 40, 40, 40]).good_stretch(0.05),
        Some((0, 3))
    );
}

#[test]
fn the_share_of_good_bases_is_counted_in_percent() {
    let t = parse(&run(Vec::new())).expect("a trace");
    assert_eq!(t.percent_at_least(30), 100.0 * 4.0 / 6.0);
    assert_eq!(t.percent_at_least(0), 100.0);
    assert!(t.has_quality());
}

/// The real files: with `HELIXLOOM_AB1_DIR` set to a folder holding `310`, `3100` and `3730` as `.ab1` files, each
/// with the `.ref` file that Biopython's reader made from it, every value must equal Biopython's (set 1 of the calls;
/// in these three files it equals set 2).
#[test]
#[ignore = "needs real trace files; see the doc comment"]
fn real_files_give_what_biopython_gives() {
    let dir = std::env::var("HELIXLOOM_AB1_DIR").expect("HELIXLOOM_AB1_DIR names the folder");
    for name in ["310", "3100", "3730"] {
        let data = std::fs::read(format!("{dir}/{name}.ab1")).expect("the trace file");
        let reference =
            std::fs::read_to_string(format!("{dir}/{name}.ref")).expect("the reference file");
        let value = |key: &str| -> String {
            reference
                .lines()
                .find_map(|l| l.strip_prefix(&format!("{key}=")))
                .unwrap_or_default()
                .to_string()
        };
        let numbers = |key: &str| -> Vec<i64> {
            value(key)
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().expect("a number"))
                .collect()
        };
        let t = parse(&data).expect("a trace");
        assert_eq!(t.sample, value("sample").trim(), "{name}: sample");
        assert_eq!(
            String::from_utf8_lossy(&t.order),
            value("order"),
            "{name}: order"
        );
        assert_eq!(
            String::from_utf8_lossy(&t.bases),
            value("bases1"),
            "{name}: bases"
        );
        let quality: Vec<i64> = t.quality.iter().map(|&q| i64::from(q)).collect();
        assert_eq!(quality, numbers("quality1"), "{name}: quality");
        let peaks: Vec<i64> = t.peaks.iter().map(|&p| i64::from(p)).collect();
        assert_eq!(peaks, numbers("peaks1"), "{name}: peaks");
        for (k, channel) in t.channels.iter().enumerate() {
            let values: Vec<i64> = channel.iter().map(|&v| i64::from(v)).collect();
            assert_eq!(
                values,
                numbers(&format!("data{}", 9 + k)),
                "{name}: channel {k}"
            );
        }
        assert!(
            t.instrument.starts_with(value("MODL1").trim()),
            "{name}: instrument {}",
            t.instrument
        );
    }
}

/// The real files against the original's own reader: with `HELIXLOOM_AB1_DIR` holding `310.orig`, `3100.orig` and
/// `3730.orig`, which `tclsh dump_read.tcl <file>` printed after running `read_abi_info` (cut out of AppMain.tcl,
/// unchanged) on each `.ab1`, every value the original reads must equal this reader's.
#[test]
#[ignore = "needs real trace files and the original's dumps; see the doc comment"]
fn real_files_read_as_the_original_reads_them() {
    let dir = std::env::var("HELIXLOOM_AB1_DIR").expect("HELIXLOOM_AB1_DIR names the folder");
    for name in ["310", "3100", "3730"] {
        let data = std::fs::read(format!("{dir}/{name}.ab1")).expect("the trace file");
        let dump =
            std::fs::read_to_string(format!("{dir}/{name}.orig")).expect("the original's dump");
        let f: std::collections::HashMap<&str, &str> =
            dump.lines().filter_map(|l| l.split_once('=')).collect();
        let numbers = |key: &str| -> Vec<i64> {
            f[key]
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| s.parse().expect("a number"))
                .collect()
        };
        let t = parse(&data).expect("a trace");
        assert_eq!(f["status"], "ok", "{name}");
        assert_eq!(t.sample, f["sample"], "{name}: sample");
        assert_eq!(
            String::from_utf8_lossy(&t.order),
            f["wheel"],
            "{name}: order"
        );
        assert_eq!(
            String::from_utf8_lossy(&t.bases),
            f["bases"],
            "{name}: bases"
        );
        let peaks: Vec<i64> = t.peaks.iter().map(|&p| i64::from(p)).collect();
        assert_eq!(peaks, numbers("loc"), "{name}: places");
        let quality: Vec<i64> = t.quality.iter().map(|&q| i64::from(q)).collect();
        assert_eq!(quality, numbers("prob"), "{name}: quality");
        for (k, channel) in t.channels.iter().enumerate() {
            let values: Vec<i64> = channel.iter().map(|&v| i64::from(v)).collect();
            assert_eq!(values, numbers(&format!("data{k}")), "{name}: channel {k}");
        }
        let strength: Vec<i64> = t.strength.iter().map(|&s| i64::from(s)).collect();
        assert_eq!(strength, numbers("scale"), "{name}: signal strength");
        let (begin, end) = (t.started.expect("a start"), t.ended.expect("an end"));
        assert_eq!(
            (begin.date_text(), begin.time_text()),
            (f["rund1"].to_string(), f["runt1"].to_string()),
            "{name}: start"
        );
        assert_eq!(
            (end.date_text(), end.time_text()),
            (f["rund2"].to_string(), f["runt2"].to_string()),
            "{name}: end"
        );
    }
}
