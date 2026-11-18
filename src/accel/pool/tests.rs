use super::*;

fn cpu(weight: f64) -> Member {
    Member::new(Device::Cpu, "CPU".into(), weight)
}

#[test]
fn the_weights_of_a_pool_add_up_to_one_and_a_bad_one_counts_as_the_mean() {
    let pool = Pool::new(vec![cpu(3.0), cpu(1.0)]);
    assert!((pool.members()[0].weight - 0.75).abs() < 1e-12);
    assert!((pool.members()[1].weight - 0.25).abs() < 1e-12);
    // zero, negative and not-a-number are mistakes: the member gets the mean of the good ones
    for bad in [0.0, -2.0, f64::NAN, f64::INFINITY] {
        let pool = Pool::new(vec![cpu(2.0), cpu(bad), cpu(4.0)]);
        let w: Vec<f64> = pool.members().iter().map(|m| m.weight).collect();
        assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12, "{w:?}");
        assert!((w[1] - 3.0 / 9.0).abs() < 1e-12, "{bad}: {w:?}");
    }
    // nothing good at all: equal shares
    let equal = Pool::new(vec![cpu(f64::NAN), cpu(0.0)]);
    assert!(equal
        .members()
        .iter()
        .all(|m| (m.weight - 0.5).abs() < 1e-12));
    assert!(Pool::new(Vec::new()).is_empty());
}

#[test]
fn the_work_is_cut_in_proportion_to_the_weights_and_covers_it_exactly() {
    // a card twice as fast as the other gets twice the work
    let pool = Pool::new(vec![cpu(2.0), cpu(1.0)]);
    assert_eq!(pool.shares(300), vec![0..200, 200..300]);
    assert_eq!(pool.shares(0), vec![0..0, 0..0]);
    assert_eq!(pool.shares(1), vec![0..1, 1..1]);
    // what does not divide goes to the largest fractions, and nothing is lost or counted twice
    let three = Pool::new(vec![cpu(1.0), cpu(1.0), cpu(1.0)]);
    for total in [0, 1, 2, 3, 4, 5, 100, 101, 102, 999_999] {
        let shares = three.shares(total);
        assert_eq!(shares.iter().map(|r| r.len()).sum::<usize>(), total);
        assert_eq!(shares.first().map(|r| r.start), Some(0));
        assert!(
            shares.windows(2).all(|w| w[0].end == w[1].start),
            "{shares:?}"
        );
        let (min, max) = (
            shares.iter().map(|r| r.len()).min().unwrap_or(0),
            shares.iter().map(|r| r.len()).max().unwrap_or(0),
        );
        assert!(max - min <= 1, "{total}: {shares:?}");
    }
    // the shares of an RTX 3090 beside an RTX 5060 Ti at 5:3 (what the bandwidths of the two make likely) for a million columns
    let cards = Pool::new(vec![cpu(5.0), cpu(3.0)]);
    assert_eq!(
        cards.shares(1_000_000),
        vec![0..625_000, 625_000..1_000_000]
    );
    assert_eq!(Pool::new(vec![cpu(1.0)]).shares(7), vec![0..7]);
}

#[test]
fn every_member_runs_on_a_thread_of_its_own_and_the_answers_keep_the_order_of_the_members() {
    let pool = Pool::new(vec![cpu(1.0), cpu(2.0), cpu(3.0)]);
    let answers = pool.run(|k, member| (k, member.weight, std::thread::current().id()));
    assert_eq!(
        answers.iter().map(|a| a.0).collect::<Vec<_>>(),
        vec![0, 1, 2]
    );
    assert!((answers[2].1 - 0.5).abs() < 1e-12);
    let mut threads: Vec<_> = answers.iter().map(|a| a.2).collect();
    threads.dedup();
    assert_eq!(threads.len(), 3);
    assert!(!threads.contains(&std::thread::current().id()));
}

#[test]
fn what_the_person_writes_for_the_mode_and_the_weights_is_read_strictly() {
    for (text, mode) in [
        ("", Mode::Auto),
        ("auto", Mode::Auto),
        ("nonsense", Mode::Auto),
        ("off", Mode::Off),
        (" OFF ", Mode::Off),
        ("cpu", Mode::Off),
        ("0", Mode::Off),
        ("force", Mode::Force),
        ("Always", Mode::Force),
    ] {
        assert_eq!(parse_mode(text), mode, "{text:?}");
    }
    assert_eq!(parse_weights("0.6,0.4", 2), Some(vec![0.6, 0.4]));
    assert_eq!(parse_weights(" 5 , 3 ", 2), Some(vec![5.0, 3.0]));
    // the wrong number of weights, a weight that is not a number, zero, negative: all ignored, never half used
    for bad in [
        "0.6",
        "0.6,0.4,0.2",
        "a,b",
        "0.6,",
        "0,1",
        "-1,2",
        "nan,1",
        "inf,1",
        "",
    ] {
        assert_eq!(parse_weights(bad, 2), None, "{bad:?}");
    }
}

#[test]
fn a_gpu_that_can_run_the_operations_has_a_speed() {
    // the CPU is too slow to time like this in a test; the first GPU there is (the Apple GPU on a Mac, a CUDA card on a
    // workstation), if there is one, answers with a positive speed
    if let Some(device) = devices().into_iter().next() {
        let speed = speed(&device).expect("a speed");
        assert!(speed.is_finite() && speed > 0.0);
    }
}

#[test]
fn a_card_that_cannot_run_the_speed_test_is_left_out_and_the_others_keep_their_weights() {
    let devices = || vec![Device::Cpu, Device::Cpu, Device::Cpu];
    // the second card gets no number: it is dropped, and the numbers of the others are theirs
    let calls = std::cell::Cell::new(0);
    let measure = |_: &Device| {
        calls.set(calls.get() + 1);
        (calls.get() != 2).then_some(calls.get() as f64)
    };
    let members = members_from(devices(), None, measure);
    assert_eq!(members.len(), 2);
    assert_eq!((members[0].weight, members[1].weight), (1.0, 3.0));
    // weights the person gave are taken as they are, and nothing is measured
    let given = [0.6, 0.4, 0.2];
    let never = |_: &Device| -> Option<f64> { panic!("measured") };
    let members = members_from(devices(), Some(&given), never);
    assert_eq!(members.len(), 3);
    assert_eq!(members[1].weight, 0.4);
    // no card can run it: no pool
    assert!(members_from(devices(), None, |_| None).is_empty());
}
