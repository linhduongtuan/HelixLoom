use super::*;
use crate::seq::dna::revcom;

/// The same 1,000 bases that primer3's and the probe's own tests use.
fn template_bases() -> Vec<u8> {
    let mut x = 777u32;
    (0..1000)
        .map(|_| {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((x >> 16) & 3) as usize]
        })
        .collect()
}

fn rows_of(n: usize) -> Vec<(String, Vec<u8>)> {
    (0..n)
        .map(|i| (format!("s{}", i + 1), template_bases()))
        .collect()
}

/// A base that is not the one there.
fn other(b: u8) -> u8 {
    if b == b'A' {
        b'C'
    } else {
        b'A'
    }
}

fn change(rows: &mut [(String, Vec<u8>)], row: usize, column: usize) {
    rows[row].1[column] = other(rows[row].1[column]);
}

fn letters(rows: &[(String, Vec<u8>)]) -> Vec<Vec<u8>> {
    rows.iter().map(|r| r.1.clone()).collect()
}

fn sites(d: &Design) -> Vec<&Site> {
    [Some(&d.left), Some(&d.right), d.probe.as_ref()]
        .into_iter()
        .flatten()
        .collect()
}

#[test]
fn rows_that_agree_everywhere_give_the_pairs_primer3_gives_for_the_sequence() {
    let t = template_bases();
    let report = design(&rows_of(3), &Settings::default()).unwrap();
    let plain = primer3::design(&t, None, &Primers::default());
    assert_eq!(report.designs.len(), plain.len());
    for (d, p) in report.designs.iter().zip(&plain) {
        assert_eq!((&d.left.oligo, &d.right.oligo), (&p.left, &p.right));
        assert_eq!((d.exact, d.tolerated), (3, 3));
        assert!(d.uncovered.is_empty());
        assert!((d.penalty - p.penalty).abs() < 1e-12);
    }
    assert_eq!((report.rows, report.columns), (3, 1000));
    assert_eq!((report.usable, report.longest_run), (1000, 1000));
    assert_eq!(report.consensus, t);
}

#[test]
fn a_column_where_the_rows_differ_is_in_no_oligo_when_all_must_agree() {
    let mut rows = rows_of(3);
    change(&mut rows, 1, 160);
    let report = design(&rows, &Settings::default()).unwrap();
    assert_eq!((report.usable, report.longest_run), (999, 839));
    assert!(!report.designs.is_empty());
    for d in &report.designs {
        for s in sites(d) {
            assert!(!(s.oligo.start..s.oligo.end).contains(&160), "{s:?}");
        }
        assert_eq!((d.exact, d.tolerated), (3, 3));
    }
}

#[test]
fn with_a_lower_share_an_oligo_may_lie_on_a_variable_column_and_the_counts_say_how_many_rows_it_fits(
) {
    let mut rows = rows_of(5);
    change(&mut rows, 3, 160);
    let s = Settings {
        min_identity: 0.6,
        ..Settings::default()
    };
    let report = design(&rows, &s).unwrap();
    // the best pair of the plain sequence has the left primer on columns 148 to 167: one row differs inside it, away from its 3' end
    let best = &report.designs[0];
    assert_eq!((best.left.oligo.start, best.left.oligo.end), (148, 168));
    assert_eq!((best.left.exact, best.left.tolerated), (4, 5));
    assert_eq!((best.exact, best.tolerated), (4, 5));
    assert!(best.uncovered.is_empty());
}

#[test]
fn a_mismatch_in_the_last_five_bases_of_a_primer_is_not_tolerated_but_one_in_a_probe_is() {
    let t = template_bases();
    let pair = &primer3::design(&t, None, &Primers::default())[0];
    let (left, right) = (&pair.left, &pair.right);
    assert_eq!(
        (left.start, left.end, right.start, right.end),
        (148, 168, 314, 334)
    );
    let s = Settings::default();
    // the left primer's 3' end is its right end (columns 163 to 167); the right primer's is its left end (314 to 318)
    let cases = [
        (left, 167, true, (false, false)),
        (left, 163, true, (false, false)),
        (left, 162, true, (false, true)),
        (left, 148, true, (false, true)),
        (right, 314, true, (false, false)),
        (right, 318, true, (false, false)),
        (right, 319, true, (false, true)),
        (right, 333, true, (false, true)),
        (left, 167, false, (false, true)),
    ];
    for (oligo, column, primer, want) in cases {
        let mut rows = rows_of(2);
        change(&mut rows, 1, column);
        let got = fits(oligo, &t, &letters(&rows), &s, primer);
        assert_eq!(got, [(true, true), want], "column {column} primer {primer}");
    }
}

#[test]
fn a_gap_is_a_mismatch_and_the_most_allowed_counts_them() {
    let t = template_bases();
    let left = primer3::design(&t, None, &Primers::default())[0]
        .left
        .clone();
    let mut rows = rows_of(3);
    rows[1].1[150] = b'-';
    rows[2].1[150] = b'-';
    rows[2].1[155] = b'-';
    let l = letters(&rows);
    let one = Settings::default();
    let none = Settings {
        max_mismatches: 0,
        ..Settings::default()
    };
    assert_eq!(
        fits(&left, &t, &l, &one, true),
        [(true, true), (false, true), (false, false)]
    );
    assert_eq!(
        fits(&left, &t, &l, &none, true),
        [(true, true), (false, false), (false, false)]
    );
}

#[test]
fn rna_and_lower_case_are_read_as_dna() {
    let mut rows = rows_of(3);
    rows[1].1 = rows[1]
        .1
        .iter()
        .map(|&b| {
            if b == b'T' {
                b'u'
            } else {
                b.to_ascii_lowercase()
            }
        })
        .collect();
    let report = design(&rows, &Settings::default()).unwrap();
    assert_eq!((report.usable, report.designs[0].exact), (1000, 3));
}

#[test]
fn a_column_that_is_mostly_gaps_is_left_out_of_the_consensus() {
    let mut rows = rows_of(4);
    for row in rows.iter_mut().take(3) {
        row.1[500] = b'-';
    }
    let report = design(&rows, &Settings::default()).unwrap();
    assert_eq!(report.consensus[500], b'-');
    assert_eq!(report.usable, 999);
    let (bases, at) = ungapped(&report.consensus);
    assert_eq!(bases.len(), 999);
    assert_eq!((at[499], at[500], at[501]), (Some(499), None, Some(500)));
    // exactly half the rows with a gap is not a majority: the column stays, though one row differs there
    let mut half = rows_of(4);
    half[0].1[500] = b'-';
    half[1].1[500] = b'-';
    let kept = design(&half, &Settings::default()).unwrap();
    assert_ne!(kept.consensus[500], b'-');
}

#[test]
fn the_product_is_counted_in_the_bases_of_each_row() {
    let mut rows = rows_of(4);
    // three bases of one row are missing inside the amplicon (columns 200 to 202)
    (200..203).for_each(|k| rows[2].1[k] = b'-');
    let s = Settings {
        min_identity: 0.5,
        ..Settings::default()
    };
    let best = &design(&rows, &s).unwrap().designs[0];
    assert_eq!((best.left.oligo.start, best.right.oligo.end), (148, 334));
    assert_eq!(best.tolerated, 4);
    assert_eq!(best.product, Some((183, 186)));
}

#[test]
fn positions_in_the_reference_row_count_the_bases_of_that_row() {
    assert_eq!(in_row(b"--AC-GT", 0, 2), None);
    assert_eq!(in_row(b"--AC-GT", 1, 4), Some((1, 2)));
    assert_eq!(in_row(b"--AC-GT", 4, 7), Some((3, 4)));
    assert_eq!(in_row(b"ACGT", 0, 4), Some((1, 4)));
    let mut rows = rows_of(4);
    (0..10).for_each(|k| rows[0].1[k] = b'-');
    let s = Settings {
        min_identity: 0.5,
        ..Settings::default()
    };
    let best = &design(&rows, &s).unwrap().designs[0];
    // the first ten columns are not in the reference row: 148 columns in, 138 of its own bases come before the primer
    assert_eq!(best.left.in_reference, Some((139, 158)));
    assert_eq!(best.right.in_reference, Some((305, 324)));
}

#[test]
fn the_sequence_of_each_oligo_is_the_consensus_or_its_reverse_complement_as_ordered() {
    let t = template_bases();
    let report = design(&rows_of(2), &Settings::default()).unwrap();
    for d in &report.designs {
        let (l, r) = (&d.left.oligo, &d.right.oligo);
        assert!(!l.reverse && r.reverse);
        assert_eq!(l.seq, t[l.start..l.end]);
        assert_eq!(r.seq, revcom(&t[r.start..r.end]));
    }
}

#[test]
fn a_probe_lies_between_the_primers_and_on_columns_all_the_rows_share() {
    let mut rows = rows_of(3);
    change(&mut rows, 2, 100);
    let s = Settings {
        with_probe: true,
        primers: Primers {
            product_min: 70,
            product_max: 150,
            ..Primers::default()
        },
        ..Settings::default()
    };
    let report = design(&rows, &s).unwrap();
    assert!(!report.designs.is_empty());
    for d in &report.designs {
        let p = d.probe.as_ref().expect("a probe");
        assert!(p.oligo.start > d.left.oligo.end && p.oligo.end < d.right.oligo.start);
        assert_eq!((p.exact, d.exact), (3, 3));
        for s in sites(d) {
            assert!(!(s.oligo.start..s.oligo.end).contains(&100));
        }
    }
    // a pair with no probe is not a design when a probe is asked for
    let none = Settings {
        probe: ProbeSettings {
            min_gc: 99.0,
            max_gc: 100.0,
            ..ProbeSettings::default()
        },
        ..s
    };
    assert!(design(&rows, &none).unwrap().designs.is_empty());
}

#[test]
fn the_search_can_be_kept_to_a_region() {
    let s = Settings {
        window: Some((400, 700)),
        ..Settings::default()
    };
    let report = design(&rows_of(2), &s).unwrap();
    assert_eq!((report.usable, report.longest_run), (300, 300));
    assert!(!report.designs.is_empty());
    for d in &report.designs {
        for s in sites(d) {
            assert!(s.oligo.start >= 400 && s.oligo.end <= 700, "{s:?}");
        }
    }
}

#[test]
fn rows_that_share_nothing_give_no_design_and_a_report_that_says_so() {
    let rows: Vec<(String, Vec<u8>)> = ["A", "C", "G", "T"]
        .iter()
        .map(|b| (format!("row{b}"), b.repeat(200).into_bytes()))
        .collect();
    let report = design(&rows, &Settings::default()).unwrap();
    assert_eq!((report.usable, report.longest_run), (0, 0));
    assert!(report.designs.is_empty());
    // a quarter of the rows agree, which is enough when the share asked for is that
    let low = Settings {
        min_identity: 0.25,
        ..Settings::default()
    };
    assert_eq!(design(&rows, &low).unwrap().usable, 200);
}

fn dummy(tolerated: usize, penalty: f64) -> Design {
    let oligo = Oligo {
        seq: b"ACGT".to_vec(),
        start: 0,
        end: 4,
        reverse: false,
        tm: 60.0,
        gc: 50.0,
        penalty: 0.0,
        self_any: 0,
        self_end: 0,
    };
    let site = Site {
        oligo,
        exact: tolerated,
        tolerated,
        in_reference: None,
    };
    Design {
        left: site.clone(),
        right: site,
        probe: None,
        penalty,
        exact: tolerated,
        tolerated,
        uncovered: Vec::new(),
        product: None,
    }
}

#[test]
fn designs_that_cover_more_rows_come_before_those_with_less_penalty() {
    let ranked = rank(
        vec![dummy(3, 0.1), dummy(5, 0.9), dummy(5, 0.4), dummy(4, 0.0)],
        3,
    );
    let order: Vec<(usize, f64)> = ranked.iter().map(|d| (d.tolerated, d.penalty)).collect();
    assert_eq!(order, [(5, 0.4), (5, 0.9), (4, 0.0)]);
}

#[test]
fn what_is_not_an_alignment_of_dna_is_refused_with_a_reason() {
    let one = &rows_of(1);
    let why = |rows: &[(String, Vec<u8>)], s: &Settings| design(rows, s).unwrap_err();
    assert!(why(one, &Settings::default()).contains("at least two"));
    let mut uneven = rows_of(2);
    uneven[1].1.pop();
    assert!(why(&uneven, &Settings::default()).contains("one length"));
    let protein: Vec<(String, Vec<u8>)> = ["MKVLAAGIVGLLLAQP", "MKVLAAGIVGLLLAQP"]
        .iter()
        .map(|p| ("p".to_string(), p.as_bytes().to_vec()))
        .collect();
    assert!(why(&protein, &Settings::default()).contains("protein"));
    for bad in [0.0, 1.5, f64::NAN] {
        let s = Settings {
            min_identity: bad,
            ..Settings::default()
        };
        assert!(why(&rows_of(2), &s).contains("at most 1"), "{bad}");
    }
    for window in [(10, 10), (600, 400), (0, 1001)] {
        let s = Settings {
            window: Some(window),
            ..Settings::default()
        };
        assert!(why(&rows_of(2), &s).contains("not inside"), "{window:?}");
    }
    let long: Vec<(String, Vec<u8>)> = (0..2)
        .map(|i| (format!("l{i}"), b"ACGT".repeat(MAX_COLUMNS / 4 + 1)))
        .collect();
    assert!(why(&long, &Settings::default()).contains("too many"));
}

#[test]
fn the_table_has_a_line_for_each_oligo_with_the_rows_that_are_not_covered() {
    let mut rows = rows_of(3);
    change(&mut rows, 2, 160);
    let s = Settings {
        min_identity: 0.6,
        with_probe: true,
        primers: Primers {
            product_min: 70,
            product_max: 150,
            ..Primers::default()
        },
        ..Settings::default()
    };
    let report = design(&rows, &s).unwrap();
    let names: Vec<String> = rows.iter().map(|r| r.0.clone()).collect();
    let text = table_text(&report, &names);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 1 + 3 * report.designs.len());
    assert!(lines[0].starts_with("design\trole\tsequence_5to3\tcolumns\tstrand"));
    let first: Vec<&str> = lines[1].split('\t').collect();
    let d = &report.designs[0];
    assert_eq!(first.len(), 14);
    assert_eq!(first[0..2], ["1", "left primer"]);
    assert_eq!(first[2].as_bytes(), d.left.oligo.seq);
    assert_eq!(
        first[3],
        format!("{}-{}", d.left.oligo.start + 1, d.left.oligo.end)
    );
    assert_eq!(first[4], "+");
    assert_eq!(lines[2].split('\t').nth(4), Some("-"));
    assert_eq!(lines[3].split('\t').nth(1), Some("probe"));
    // the numbers are the design's
    let (low, high) = d.product.unwrap();
    assert_eq!(
        (first[10], first[11], first[12]),
        (
            low.to_string().as_str(),
            high.to_string().as_str(),
            d.tolerated.to_string().as_str()
        )
    );
    assert_eq!(first[13], "");
}

#[test]
fn a_dot_is_a_gap_too_and_a_tie_between_bases_is_broken_the_same_way_every_time() {
    let mut rows = rows_of(4);
    for row in rows.iter_mut().take(3) {
        row.1[500] = b'.';
    }
    let report = design(&rows, &Settings::default()).unwrap();
    assert_eq!((report.consensus[500], report.usable), (b'-', 999));
    // two rows that differ everywhere: half and half, and the first of A, C, G, T wins
    let tie: Vec<(String, Vec<u8>)> = vec![
        ("a".into(), b"ACGTACGT".to_vec()),
        ("b".into(), b"CATGCATG".to_vec()),
    ];
    let low = Settings {
        min_identity: 0.5,
        ..Settings::default()
    };
    assert_eq!(design(&tie, &low).unwrap().consensus, b"AAGGAAGG");
}

#[test]
fn the_design_that_covers_every_row_is_chosen_over_one_with_less_penalty_that_misses_one() {
    let mut rows = rows_of(5);
    // the best pair by penalty has its left primer on columns 148 to 167: one row differs in the last of them, a primer's 3' end
    change(&mut rows, 4, 167);
    let s = Settings {
        min_identity: 0.5,
        primers: Primers {
            num_return: 1,
            ..Primers::default()
        },
        ..Settings::default()
    };
    let report = design(&rows, &s).unwrap();
    assert_eq!(report.designs.len(), 1);
    let best = &report.designs[0];
    assert_eq!((best.tolerated, best.uncovered.len()), (5, 0));
    let plain = &primer3::design(&template_bases(), None, &Primers::default())[0];
    assert_eq!(plain.left.end, 168);
    assert_ne!(
        (best.left.oligo.start, best.right.oligo.end),
        (plain.left.start, plain.right.end)
    );
}
