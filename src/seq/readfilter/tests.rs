use super::*;

fn read(len: usize, q: u8, name: &str) -> FastqRead {
    FastqRead {
        name: name.to_string(),
        seq: vec![b'A'; len],
        qual: vec![q; len],
    }
}

#[test]
fn the_defaults_are_the_papers_two_kilobases_five_hundred_megabases_and_a_hundred_base_window() {
    assert_eq!(
        Filter::default(),
        Filter {
            min_length: 2000,
            target_bases: 500_000_000,
            window: 100
        }
    );
}

#[test]
fn the_worst_window_is_the_lowest_mean_over_any_run_of_the_window_s_length() {
    let mut q = vec![30u8; 300];
    assert_eq!(worst_window(&q, 100), 30.0);
    q[100..200].iter_mut().for_each(|x| *x = 10);
    assert_eq!(worst_window(&q, 100), 10.0);
    // half of a window low: the mean of the worst window is between
    let mut half = vec![30u8; 300];
    half[100..150].iter_mut().for_each(|x| *x = 10);
    assert_eq!(worst_window(&half, 100), 20.0);
    // shorter than the window: the mean of the read
    assert_eq!(worst_window(&[10, 20, 30], 100), 20.0);
    assert_eq!(worst_window(&[], 100), 0.0);
    // a window of nothing is of one base
    assert_eq!(worst_window(&[5, 9], 0), 5.0);
}

#[test]
fn reads_under_the_length_are_dropped_and_one_of_exactly_the_length_stays() {
    let reads = vec![
        read(1999, 30, "a"),
        read(2000, 30, "b"),
        read(5000, 30, "c"),
    ];
    let f = filter_long(&reads, &Filter::default());
    assert_eq!(
        f.reads.iter().map(|r| r.name.as_str()).collect::<Vec<_>>(),
        ["b", "c"]
    );
    assert_eq!((f.too_short, f.low_quality), (1, 0));
    assert_eq!((f.bases_in, f.bases_out), (8999, 7000));
}

#[test]
fn over_the_target_the_reads_with_the_worst_stretch_go_first_until_it_is_met() {
    let reads = vec![
        read(3000, 30, "good"),
        read(3000, 12, "bad"),
        read(3000, 20, "middle"),
        read(3000, 25, "fair"),
    ];
    let f = Filter {
        target_bases: 9000,
        ..Filter::default()
    };
    let kept = filter_long(&reads, &f);
    assert_eq!(
        kept.reads
            .iter()
            .map(|r| r.name.as_str())
            .collect::<Vec<_>>(),
        ["good", "middle", "fair"]
    );
    assert_eq!((kept.low_quality, kept.bases_out), (1, 9000));
    // two must go
    let f = Filter {
        target_bases: 8999,
        ..f
    };
    let kept = filter_long(&reads, &f);
    assert_eq!(
        kept.reads
            .iter()
            .map(|r| r.name.as_str())
            .collect::<Vec<_>>(),
        ["good", "fair"]
    );
    assert_eq!(kept.low_quality, 2);
}

#[test]
fn a_read_is_scored_by_its_worst_window_not_its_mean() {
    // a has a mean of 27 and a stretch of 10; b is 22 all through
    let mut a = read(3000, 30, "a");
    a.qual[1000..1100].iter_mut().for_each(|q| *q = 10);
    let b = read(3000, 22, "b");
    let f = filter_long(
        &[a, b],
        &Filter {
            target_bases: 3000,
            ..Filter::default()
        },
    );
    assert_eq!(
        f.reads.iter().map(|r| r.name.as_str()).collect::<Vec<_>>(),
        ["b"]
    );
}

#[test]
fn of_two_reads_as_bad_the_later_goes_first_and_the_order_that_is_left_is_the_original() {
    let reads = vec![
        read(3000, 20, "first"),
        read(3000, 20, "second"),
        read(3000, 35, "third"),
    ];
    let f = filter_long(
        &reads,
        &Filter {
            target_bases: 6000,
            ..Filter::default()
        },
    );
    assert_eq!(
        f.reads.iter().map(|r| r.name.as_str()).collect::<Vec<_>>(),
        ["first", "third"]
    );
}

#[test]
fn under_the_target_nothing_is_dropped_for_quality_and_nothing_is_a_nothing() {
    let reads = vec![read(3000, 5, "poor"), read(3000, 6, "poor2")];
    let f = filter_long(&reads, &Filter::default());
    assert_eq!((f.reads.len(), f.low_quality), (2, 0));
    let none = filter_long(&[], &Filter::default());
    assert_eq!((none.reads.len(), none.bases_in, none.bases_out), (0, 0, 0));
}

#[test]
fn a_short_read_does_not_count_toward_the_target() {
    // 1000 bases of short reads and 6000 of long ones, a target of 6000: nothing is dropped for quality
    let reads = vec![
        read(1000, 40, "short"),
        read(3000, 10, "a"),
        read(3000, 10, "b"),
    ];
    let f = filter_long(
        &reads,
        &Filter {
            target_bases: 6000,
            ..Filter::default()
        },
    );
    assert_eq!((f.too_short, f.low_quality, f.reads.len()), (1, 0, 2));
}
