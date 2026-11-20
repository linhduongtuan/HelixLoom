use super::*;
use crate::seq::dna::{gc_percent, revcom};

/// The same 1,000 bases that primer3's own tests use.
fn template() -> Vec<u8> {
    let mut x = 777u32;
    (0..1000)
        .map(|_| {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((x >> 16) & 3) as usize]
        })
        .collect()
}

/// Products of the length that real-time PCR uses.
fn primers() -> Primers {
    Primers {
        product_min: 70,
        product_max: 150,
        ..Primers::default()
    }
}

/// Every limit on what one oligo is, and none on how it sticks to itself.
fn lenient() -> Settings {
    Settings {
        max_self_any: 1000,
        max_self_end: 1000,
        ..Settings::default()
    }
}

fn is_ok(seq: &[u8], s: &Settings) -> bool {
    evaluate(seq.to_vec(), (0, seq.len()), false, &primers(), s).is_some()
}

#[test]
fn every_probe_found_follows_the_rules_and_lies_between_its_primers() {
    let (t, p, s) = (template(), primers(), Settings::default());
    let found = design_assays(&t, None, &p, &s);
    assert!(!found.is_empty() && found.len() <= p.num_return);
    for a in &found {
        let (l, r, o) = (&a.pair.left, &a.pair.right, &a.probe);
        assert!(
            o.start >= l.end + s.min_gap && o.end + s.min_gap <= r.start,
            "{a:?}"
        );
        // the probe is the template there, or its other strand, as it would be ordered
        let slice = &t[o.start..o.end];
        let ordered = if o.reverse {
            revcom(slice)
        } else {
            slice.to_vec()
        };
        assert_eq!(o.seq, ordered);
        assert!((s.min_size..=s.max_size).contains(&o.seq.len()));
        assert_ne!(o.seq[0], b'G', "a G at the 5' end quenches the reporter");
        assert!((s.min_gc..=s.max_gc).contains(&o.gc));
        assert_eq!(o.gc, gc_percent(&o.seq));
        assert!(primer3::longest_run(&o.seq) <= s.max_run);
        assert!(!o
            .seq
            .windows(s.max_g_run + 1)
            .any(|w| w.iter().all(|&c| c == b'G')));
        // it melts the stated number of degrees above the higher of the two primers
        assert_eq!(primer3::tm(&o.seq, &p), Some(o.tm));
        let above = o.tm - l.tm.max(r.tm);
        assert!((s.min_above..=s.max_above).contains(&above), "{above}");
        assert!(o.self_any <= s.max_self_any && o.self_end <= s.max_self_end);
        assert!((a.penalty - (a.pair.penalty + o.penalty)).abs() < 1e-9);
        // it is the best of all those that fit its pair
        assert_eq!(probes_for(&t, &a.pair, &p, &s).first(), Some(o));
    }
    assert!(found.windows(2).all(|w| w[0].penalty <= w[1].penalty));
    // asked twice it says the same
    assert_eq!(design_assays(&t, None, &p, &s), found);
}

#[test]
fn each_rule_on_its_own_refuses_what_it_names_and_nothing_else() {
    let s = lenient();
    let good = b"CATCAGTCAGCATGACTAGCATCA";
    assert!(is_ok(good, &s), "the baseline passes");
    // no G at the 5' end
    let g_first = b"GATCAGTCAGCATGACTAGCATCA";
    assert!(!is_ok(g_first, &s));
    assert!(is_ok(
        g_first,
        &Settings {
            no_five_prime_g: false,
            ..s.clone()
        }
    ));
    // runs of G: three are allowed, four are not
    assert!(is_ok(b"CATCAGTCAGGGATGACTAGCATC", &s));
    let four_g = b"CATCAGTCAGGGGATGACTAGCAT";
    assert!(!is_ok(four_g, &s));
    assert!(is_ok(
        four_g,
        &Settings {
            max_g_run: 4,
            ..s.clone()
        }
    ));
    // runs of any base: four are allowed, five are not
    assert!(is_ok(b"CATCAGTCAAAATGACTAGCATCA", &s));
    assert!(!is_ok(b"CATCAGTCAAAAATGACTAGCATC", &s));
    // G+C between the limits
    let few = b"ATTATATTCATAATATTAGATTTA";
    assert!(
        !is_ok(few, &s)
            && is_ok(
                few,
                &Settings {
                    min_gc: 0.0,
                    ..s.clone()
                }
            )
    );
    let many = b"CCGCGCGCCGCCGCGCCGCGCGCC";
    assert!(!is_ok(many, &s));
    assert!(is_ok(
        many,
        &Settings {
            max_gc: 100.0,
            max_run: 9,
            max_g_run: 9,
            ..s.clone()
        }
    ));
    // letters that are not bases
    assert!(!is_ok(b"CATCAGTCAGCATNACTAGCATCA", &s));
}

#[test]
fn a_probe_that_sticks_to_itself_is_refused_unless_the_limit_is_lifted() {
    let palindrome = b"ACGTACGTACGTACGTACGTACGT";
    assert!(!is_ok(palindrome, &Settings::default()));
    assert!(is_ok(palindrome, &lenient()));
    // the limit on the 3' end alone: the one on the whole is out of the way
    let end_only = Settings {
        max_self_any: 1000,
        max_self_end: 5,
        ..Settings::default()
    };
    assert!(!is_ok(palindrome, &end_only));
}

#[test]
fn the_strand_that_does_not_begin_with_g_is_the_one_that_is_chosen() {
    let p = primers();
    let s = Settings {
        min_size: 24,
        max_size: 24,
        ..lenient()
    };
    // the plus strand begins with G, the minus strand does not
    let plus_g = b"GCATCAGTCAGCATGACTAGCATA".to_vec();
    let found = candidates(&plus_g, 0, plus_g.len(), &p, &s);
    assert_eq!(found.len(), 1);
    assert!(found[0].reverse && found[0].seq == revcom(&plus_g));
    // and the other way round
    let minus_g = revcom(&plus_g);
    let found = candidates(&minus_g, 0, minus_g.len(), &p, &s);
    assert_eq!(found.len(), 1);
    assert!(!found[0].reverse && found[0].seq == minus_g);
    // with the rule off, both strands are candidates
    let both = Settings {
        no_five_prime_g: false,
        ..s
    };
    assert_eq!(candidates(&plus_g, 0, plus_g.len(), &p, &both).len(), 2);
}

/// Every placement, strand and length, judged by a plain statement of the rules; the penalty as the module's text gives it.
fn brute_force(
    t: &[u8],
    pair: &Pair,
    p: &Primers,
    s: &Settings,
) -> Vec<((usize, usize, bool), f64)> {
    let (from, to) = (pair.left.end, pair.right.start);
    let mut out = Vec::new();
    for len in s.min_size..=s.max_size {
        for start in from..=to.saturating_sub(len) {
            if start + len > to {
                continue;
            }
            for reverse in [false, true] {
                let seq = if reverse {
                    revcom(&t[start..start + len])
                } else {
                    t[start..start + len].to_vec()
                };
                let count = |b: u8| seq.iter().filter(|&&c| c == b).count();
                let gc = 100.0 * (count(b'G') + count(b'C')) as f64 / len as f64;
                let run_of = |b: Option<u8>| {
                    let (mut best, mut now) = (0, 0);
                    for (i, &c) in seq.iter().enumerate() {
                        now = if i > 0 && c == seq[i - 1] { now + 1 } else { 1 };
                        if b.is_none_or(|b| b == c) {
                            best = best.max(now);
                        }
                    }
                    best
                };
                let tm = primer3::tm(&seq, p).unwrap();
                let above = tm - pair.left.tm.max(pair.right.tm);
                let (any, end) = primer3::dimer_scores(&seq, &seq);
                let sticks = [&pair.left, &pair.right].iter().any(|o| {
                    let (a, e1) = primer3::dimer_scores(&o.seq, &seq);
                    let (_, e2) = primer3::dimer_scores(&seq, &o.seq);
                    a > p.pair_compl_any || e1.max(e2) > p.pair_compl_end
                });
                let gap = (start - pair.left.end).min(pair.right.start - (start + len));
                let rules = (!s.no_five_prime_g || seq[0] != b'G')
                    && (s.min_gc..=s.max_gc).contains(&gc)
                    && run_of(None) <= s.max_run
                    && run_of(Some(b'G')) <= s.max_g_run
                    && (s.min_above..=s.max_above).contains(&above)
                    && any <= s.max_self_any
                    && end <= s.max_self_end
                    && !sticks
                    && start >= pair.left.end + s.min_gap
                    && start + len + s.min_gap <= pair.right.start;
                if rules {
                    let surplus = count(b'G').saturating_sub(count(b'C')) as f64;
                    out.push((
                        (start, start + len, reverse),
                        0.25 * (len as f64 - s.opt_size as f64).abs()
                            + 0.5 * surplus
                            + (above - s.opt_above).abs()
                            + 0.05 * gap as f64,
                    ));
                }
            }
        }
    }
    out.sort_by_key(|a| a.0);
    out
}

#[test]
fn the_probes_for_a_pair_are_all_of_those_that_pass_and_the_best_comes_first() {
    let (t, p) = (template(), primers());
    for s in [
        Settings::default(),
        Settings {
            min_above: 2.0,
            max_above: 20.0,
            min_gap: 4,
            ..Settings::default()
        },
        Settings {
            no_five_prime_g: false,
            min_size: 20,
            max_size: 22,
            ..Settings::default()
        },
        // a probe may touch the primer next to it
        Settings {
            min_gap: 0,
            ..Settings::default()
        },
    ] {
        let mut passing = 0;
        let pairs = primer3::design(
            &t,
            None,
            &Primers {
                num_return: 8,
                ..p.clone()
            },
        );
        assert!(!pairs.is_empty());
        for pair in &pairs {
            let found = probes_for(&t, pair, &p, &s);
            // best first
            assert!(found.windows(2).all(|w| w[0].penalty <= w[1].penalty));
            let mut got: Vec<_> = found
                .iter()
                .map(|o| ((o.start, o.end, o.reverse), o.penalty))
                .collect();
            got.sort_by_key(|a| a.0);
            let want = brute_force(&t, pair, &p, &s);
            passing += want.len();
            let (a, b): (Vec<_>, Vec<_>) = (
                got.iter().map(|x| x.0).collect(),
                want.iter().map(|x| x.0).collect(),
            );
            assert_eq!(a, b, "{pair:?}");
            assert!(got.iter().zip(&want).all(|(a, b)| (a.1 - b.1).abs() < 1e-9));
        }
        // the comparison is not of nothing with nothing
        assert!(passing > 0, "no probe passes for {s:?}");
    }
}

#[test]
fn the_g_plus_c_limits_are_inclusive_at_thirty_and_eighty_percent() {
    let s = lenient();
    // twenty bases: 6 G or C is 30 %, 5 is 25 %; 16 is 80 %, 17 is 85 %
    assert!(is_ok(b"ATCATGATCAGTCATAGATA", &s));
    assert!(!is_ok(b"ATCATGATCAGTAATAGATA", &s));
    assert!(is_ok(b"CCGCGCATGCGCCGATCGCG", &s));
    assert!(!is_ok(b"CCGCGCATGCGCCGCTCGCG", &s));
}

#[test]
fn the_assays_are_the_best_pairs_each_with_its_best_probe_for_any_case_of_the_template() {
    let (t, p) = (template(), primers());
    for s in [
        Settings::default(),
        Settings {
            min_gap: 0,
            ..Settings::default()
        },
        Settings {
            min_above: 2.0,
            max_above: 20.0,
            no_five_prime_g: false,
            ..Settings::default()
        },
    ] {
        // the pairs the module asks primer3 for (eight for each assay wanted, forty at least), each with the best probe that a plain
        // search finds, put in order by their totals
        let asked = Primers {
            num_return: (p.num_return * 8).max(40),
            ..p.clone()
        };
        let pairs = primer3::design(&t, None, &asked);
        let mut want: Vec<(f64, &Pair)> = pairs
            .iter()
            .filter_map(|pair| {
                let best = brute_force(&t, pair, &p, &s)
                    .into_iter()
                    .map(|x| x.1)
                    .min_by(f64::total_cmp)?;
                Some((pair.penalty + best, pair))
            })
            .collect();
        want.sort_by(|a, b| a.0.total_cmp(&b.0));
        want.truncate(p.num_return);
        assert!(!want.is_empty(), "no assay for {s:?}");
        let lower = t.to_ascii_lowercase();
        for template in [&t, &lower] {
            let got = design_assays(template, None, &p, &s);
            assert_eq!(got.len(), want.len(), "{s:?}");
            for (g, (penalty, pair)) in got.iter().zip(&want) {
                assert!((g.penalty - penalty).abs() < 1e-9, "{s:?}");
                assert_eq!(&g.pair, *pair, "{s:?}");
            }
        }
    }
}

#[test]
fn a_probe_that_would_stick_to_a_primer_is_refused() {
    let (t, p, s) = (template(), primers(), Settings::default());
    let a = &design_assays(&t, None, &p, &s)[0];
    // the probe as it was before the pair gave it its penalty
    let own = candidates(&t, a.pair.left.end, a.pair.right.start, &p, &s)
        .into_iter()
        .find(|c| (c.start, c.end, c.reverse) == (a.probe.start, a.probe.end, a.probe.reverse))
        .expect("the candidate");
    assert_eq!(fits_pair(&own, &a.pair, &p, &s), Some(a.probe.penalty));
    // a probe that is the complement of a primer ends up on it
    for primer in [&a.pair.left, &a.pair.right] {
        let stuck = Oligo {
            seq: revcom(&primer.seq),
            ..own.clone()
        };
        assert_eq!(fits_pair(&stuck, &a.pair, &p, &s), None);
    }
}

#[test]
fn a_probe_that_overlaps_a_primer_or_is_too_close_or_melts_too_little_or_too_much_is_refused() {
    let (t, p, s) = (template(), primers(), Settings::default());
    let a = &design_assays(&t, None, &p, &s)[0];
    let moved = |start: usize, end: usize| Oligo {
        start,
        end,
        ..a.probe.clone()
    };
    let len = a.probe.end - a.probe.start;
    let (l, r) = (&a.pair.left, &a.pair.right);
    assert_eq!(
        fits_pair(&moved(l.end - 1, l.end - 1 + len), &a.pair, &p, &s),
        None
    );
    assert_eq!(
        fits_pair(&moved(r.start - len + 1, r.start + 1), &a.pair, &p, &s),
        None
    );
    // a gap of the limit is allowed, one base less is not
    let at_limit = moved(l.end + s.min_gap, l.end + s.min_gap + len);
    assert!(fits_pair(&at_limit, &a.pair, &p, &s).is_some());
    assert_eq!(
        fits_pair(
            &moved(l.end + s.min_gap - 1, l.end + s.min_gap - 1 + len),
            &a.pair,
            &p,
            &s
        ),
        None
    );
    let hot = Oligo {
        tm: l.tm.max(r.tm) + s.max_above + 0.01,
        ..a.probe.clone()
    };
    let cold = Oligo {
        tm: l.tm.max(r.tm) + s.min_above - 0.01,
        ..a.probe.clone()
    };
    assert_eq!(fits_pair(&hot, &a.pair, &p, &s), None);
    assert_eq!(fits_pair(&cold, &a.pair, &p, &s), None);
}

#[test]
fn a_gap_of_no_probe_is_no_assay_and_the_limits_that_cannot_be_met_leave_none() {
    let t = template();
    // products too short to hold a probe between two primers
    let short = Primers {
        product_min: 40,
        product_max: 50,
        ..Primers::default()
    };
    assert!(
        !primer3::design(&t, None, &short).is_empty(),
        "pairs exist, which is the point"
    );
    assert!(design_assays(&t, None, &short, &Settings::default()).is_empty());
    // a probe that must melt 40 to 50 degrees above the primers is none
    let hot = Settings {
        min_above: 40.0,
        max_above: 50.0,
        ..Settings::default()
    };
    assert!(design_assays(&t, None, &primers(), &hot).is_empty());
    // as many as asked for, and no more
    let few = Primers {
        num_return: 2,
        ..primers()
    };
    assert!(design_assays(&t, None, &few, &Settings::default()).len() <= 2);
    // nothing in, nothing out
    assert!(design_assays(b"", None, &primers(), &Settings::default()).is_empty());
}

#[test]
fn a_target_is_between_the_primers_and_the_probe_may_cover_it() {
    let (t, p, s) = (template(), primers(), Settings::default());
    let found = design_assays(&t, Some((500, 20)), &p, &s);
    assert!(!found.is_empty());
    for a in &found {
        // primers keep off the target (as primer3's do) and the probe lies between them
        assert!(a.pair.left.end <= 500 && a.pair.right.start >= 520);
        assert!(a.probe.start >= a.pair.left.end && a.probe.end <= a.pair.right.start);
    }
}

#[test]
fn the_defaults_are_the_published_guidelines_and_the_choices_the_module_names() {
    let s = Settings::default();
    // Applied Biosystems: no G at the 5' end, no run of four G, G+C 30 to 80 %, at most 30 bases, about 10 degrees above the primers
    assert!(s.no_five_prime_g);
    assert_eq!(s.max_g_run, 3);
    assert_eq!((s.min_gc, s.max_gc), (30.0, 80.0));
    assert_eq!(s.max_size, 30);
    assert!(s.min_size <= s.opt_size && s.opt_size <= s.max_size);
    // this program's own choices: the window around the 10 degrees, runs of other bases, the gap to the primers
    assert_eq!((s.min_above, s.opt_above, s.max_above), (6.0, 10.0, 14.0));
    assert_eq!((s.max_run, s.min_gap), (4, 1));
}

#[test]
fn stretches_that_overlap_touch_or_lie_inside_one_another_are_joined_and_the_rest_are_kept_apart() {
    assert_eq!(joined(vec![]), vec![]);
    // one inside another: the outer stands whichever comes first
    assert_eq!(joined(vec![(100, 200), (120, 150)]), vec![(100, 200)]);
    assert_eq!(joined(vec![(120, 150), (100, 200)]), vec![(100, 200)]);
    // overlapping, and touching
    assert_eq!(joined(vec![(100, 150), (120, 200)]), vec![(100, 200)]);
    assert_eq!(joined(vec![(150, 200), (100, 150)]), vec![(100, 200)]);
    // apart, in any order, and a chain of three
    assert_eq!(
        joined(vec![(300, 350), (100, 150)]),
        vec![(100, 150), (300, 350)]
    );
    assert_eq!(
        joined(vec![(100, 150), (140, 180), (170, 220), (400, 410)]),
        vec![(100, 220), (400, 410)]
    );
}
