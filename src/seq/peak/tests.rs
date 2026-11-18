use super::*;

/// Six bases, GATTNC, with a peak every ten scans from scan 5; the channels stand for G, A, T, C in that order, and each
/// base's own channel is 100 times its place (from 1) at its peak, the others less, as written below.
fn trace() -> Trace {
    let mut channels: [Vec<i16>; 4] = std::array::from_fn(|_| vec![0; 70]);
    // at each peak: G, A, T, C
    let at_peaks: [[i16; 4]; 6] = [
        [100, 7, 3, 1],   // G
        [9, 200, 9, 4],   // A: G and T are as strong as each other
        [0, 0, 300, 0],   // T
        [0, 0, 400, 20],  // T
        [50, 60, 55, 58], // N: no channel stands out
        [0, 700, 0, 600], // called C, though A is stronger
    ];
    for (base, signal) in at_peaks.iter().enumerate() {
        for (k, &s) in signal.iter().enumerate() {
            channels[k][5 + 10 * base] = s;
        }
    }
    Trace {
        order: *b"GATC",
        channels,
        bases: b"GATTNC".to_vec(),
        quality: vec![42, 20, 19, 60, 3, 3],
        peaks: vec![5, 15, 25, 35, 45, 55],
        ..Trace::default()
    }
}

fn one(t: &Trace, at: usize) -> Peak {
    match picked(t, at, at + 1) {
        Some(Picked::One(peak)) => peak,
        other => panic!("base {at}: {other:?}"),
    }
}

fn several(t: &Trace, from: usize, to: usize) -> Stretch {
    match picked(t, from, to) {
        Some(Picked::Several(stretch)) => stretch,
        other => panic!("bases {from} to {to}: {other:?}"),
    }
}

#[test]
fn one_base_is_its_call_its_quality_its_scan_and_the_four_signals_the_strongest_first() {
    let t = trace();
    let g = one(&t, 0);
    assert_eq!(
        (g.at, g.call, g.quality, g.scan),
        (0, b'G', Some(42), Some(5))
    );
    assert_eq!(g.signal, [(b'G', 100), (b'A', 7), (b'T', 3), (b'C', 1)]);
    // signals as strong as each other stay in the order of the file: G before T
    let a = one(&t, 1);
    assert_eq!(a.signal, [(b'A', 200), (b'G', 9), (b'T', 9), (b'C', 4)]);
    // the signal is read at the base's own peak, not at a neighbour's
    assert_eq!(one(&t, 3).signal[0], (b'T', 400));
    assert_eq!(one(&t, 3).scan, Some(35));
    // the strongest is listed first whatever was called
    let c = one(&t, 5);
    assert_eq!(
        (c.call, c.signal[0], c.signal[1]),
        (b'C', (b'A', 700), (b'C', 600))
    );
}

#[test]
fn several_bases_are_their_mean_their_lowest_and_how_many_are_poor_or_unsure() {
    let t = trace();
    let all = several(&t, 0, 6);
    assert_eq!((all.from, all.to), (0, 6));
    let mean = all.mean.expect("a mean");
    assert!(
        (mean - (42.0 + 20.0 + 19.0 + 60.0 + 3.0 + 3.0) / 6.0).abs() < 1e-9,
        "{mean}"
    );
    // the first of the two bases of quality 3
    assert_eq!(all.lowest, Some((3, 4)));
    // 19, 3 and 3 are below 20; 20 itself is not
    assert_eq!(all.poor, Some(3));
    assert_eq!(all.unsure, 1);
    // a part of them: only what is in it counts
    let part = several(&t, 0, 3);
    assert_eq!(
        (part.lowest, part.poor, part.unsure),
        (Some((19, 2)), Some(1), 0)
    );
    let two = several(&t, 1, 3);
    assert!((two.mean.expect("a mean") - 19.5).abs() < 1e-9);
    // a lower-case call that is a base is not unsure
    let mut lower = trace();
    lower.bases = b"gattnc".to_vec();
    assert_eq!(several(&lower, 0, 6).unsure, 1);
}

#[test]
fn a_file_with_no_qualities_gives_none_and_a_base_without_a_peak_has_no_signal() {
    let mut t = trace();
    t.quality = vec![0; 6];
    assert_eq!(one(&t, 0).quality, None);
    let all = several(&t, 0, 6);
    assert_eq!((all.mean, all.lowest, all.poor), (None, None, None));
    assert_eq!(all.unsure, 1, "the calls are counted all the same");
    // a quality of nothing among others is a quality: the file gives them
    let mut some = trace();
    some.quality[1] = 0;
    assert_eq!(one(&some, 1).quality, Some(0));
    assert_eq!(several(&some, 0, 3).lowest, Some((0, 1)));
    // a file whose peaks or qualities stop short of its bases
    let mut short = trace();
    short.peaks.truncate(4);
    short.quality.truncate(5);
    let last = one(&short, 5);
    assert_eq!((last.scan, last.quality), (None, None));
    assert!(last.signal.is_empty());
    assert_eq!(several(&short, 3, 6).lowest, Some((3, 4)));
    // the mean is of the bases that have a quality (60 and 3), not of all three that were picked
    assert_eq!(several(&short, 3, 6).mean, Some(31.5));
    assert_eq!(several(&short, 3, 6).poor, Some(1));
    // a peak past the end of the channels reads as no signal, not as a fault
    let mut beyond = trace();
    beyond.peaks[5] = 500;
    assert!(one(&beyond, 5).signal.iter().all(|&(_, s)| s == 0));
}

#[test]
fn bases_that_are_not_there_are_not_picked() {
    let t = trace();
    assert_eq!(picked(&t, 3, 3), None);
    assert_eq!(picked(&t, 4, 2), None);
    assert_eq!(picked(&t, 6, 9), None);
    // what reaches past the end is cut to the end
    assert_eq!(several(&t, 4, 60).to, 6);
    assert!(matches!(picked(&t, 5, 60), Some(Picked::One(_))));
    assert_eq!(picked(&Trace::default(), 0, 1), None);
}

#[test]
fn a_quality_is_one_wrong_call_in_so_many() {
    assert_eq!(one_in(0), 1);
    assert_eq!(one_in(3), 2);
    assert_eq!(one_in(10), 10);
    assert_eq!(one_in(20), 100);
    assert_eq!(one_in(42), 15_849);
    assert_eq!(one_in(60), 1_000_000);
}
