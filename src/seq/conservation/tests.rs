use super::*;

fn rows(text: &[&str]) -> Vec<Vec<u8>> {
    text.iter().map(|r| r.as_bytes().to_vec()).collect()
}

/// Scores as the reference gives them: the same columns, each within a small tolerance; `None` where it says so.
fn assert_scores(got: &[Option<f32>], want: &[f64]) {
    assert_eq!(got.len(), want.len(), "{got:?}");
    for (k, (g, w)) in got.iter().zip(want).enumerate() {
        let g = g.unwrap_or_else(|| panic!("column {k} has no score")) as f64;
        assert!((g - w).abs() < 1e-4, "column {k}: {g} against {w}");
    }
}

/// Eight proteins of 22 columns, with gaps; the numbers below come from a separate, plain implementation of the same
/// published formulas (`reference.py`, next to this file: run it to see them), written for the purpose and not part of the
/// program.
fn protein() -> Vec<Vec<u8>> {
    rows(&[
        "MKTAYIAKQRQISFVKSHFSRQ",
        "MKTAYIAKQRQISFVKSHFSRQ",
        "MKSAYIAKQRQVSFVKSHFSRQ",
        "MKTAYLAKQRQISFVRSHFARQ",
        "MRTAYIAKQK-ISFVKSHFSRQ",
        "MKTGYIAKQRQISFVKTHF-RQ",
        "MKTAYIVKQRQISFAKSHFSRE",
        "MQTAYIAKQRQ-SFVKSHFSRQ",
    ])
}

fn dna() -> Vec<Vec<u8>> {
    rows(&[
        "ACGTACGTAC",
        "ACGTACGAAC",
        "ACGAACGTAC",
        "ACGTTCGTAC",
        "ACG-ACGTAC",
        "ACGTACGTGC",
    ])
}

fn plain(method: Method) -> Settings {
    Settings {
        method,
        window: 0,
        ..Settings::default()
    }
}

#[test]
fn sequence_weights_follow_henikoff_and_henikoff_by_hand() {
    // columns: A in all three (each row 1/3); D in two and E in one (1/4, 1/4 and 1/2): 7/12, 7/12, 10/12, over 2
    let w = weights(&rows(&["AD", "AD", "AE"]), true);
    for (got, want) in w.iter().zip([7.0 / 24.0, 7.0 / 24.0, 10.0 / 24.0]) {
        assert!((got - want).abs() < 1e-12, "{w:?}");
    }
    assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    // rows that are all alike weigh the same; with no letters at all every row weighs the same
    let alike = weights(&rows(&["ACD", "ACD", "ACD", "ACD"]), true);
    assert!(alike.iter().all(|w| (w - 0.25).abs() < 1e-12));
    let none = weights(&rows(&["--", "--"]), true);
    assert_eq!(none, vec![0.5, 0.5]);
    // a gap adds nothing to its row
    let gap = weights(&rows(&["A-", "AD", "AD"]), true);
    assert!(gap[0] < gap[1]);
}

#[test]
fn the_weights_of_the_protein_example_match_the_reference() {
    let w = weights(&protein(), true);
    let want = [
        0.096952, 0.096952, 0.135372, 0.154852, 0.122565, 0.132125, 0.155393, 0.10579,
    ];
    for (got, want) in w.iter().zip(want) {
        assert!((got - want).abs() < 1e-5, "{w:?}");
    }
}

#[test]
fn entropy_is_worked_out_by_hand_for_nucleotides() {
    // a column of one letter is fully conserved, an even mix has no conservation, half and half is one bit of two
    let s = Settings {
        weights: false,
        max_gaps: 1.0,
        ..plain(Method::Entropy)
    };
    let c = |a: &str, b: &str, c: &str, d: &str| scores(&rows(&[a, b, c, d]), false, &s)[0];
    for (col, want) in [
        (["A", "A", "A", "A"], 1.0),
        (["A", "C", "G", "T"], 0.0),
        (["A", "A", "C", "C"], 0.5),
    ] {
        let got = c(col[0], col[1], col[2], col[3]).expect("a score") as f64;
        assert!((got - want).abs() < 1e-4, "{col:?}: {got}");
    }
    // U is T, and case does not matter
    let u = scores(&rows(&["t", "u", "T", "U"]), false, &s)[0].expect("a score");
    assert!(u > 0.999, "{u}");
}

#[test]
fn the_jensen_shannon_divergence_of_the_protein_example_matches_the_reference() {
    assert_scores(
        &scores(&protein(), true, &plain(Method::Jsd)),
        &[
            0.918022, 0.720812, 0.764030, 0.702613, 0.892278, 0.719602, 0.705155, 0.842276,
            0.892278, 0.777667, 0.782914, 0.657878, 0.835983, 0.868614, 0.711136, 0.765831,
            0.760168, 0.915333, 0.868614, 0.634036, 0.853011, 0.805458,
        ],
    );
}

#[test]
fn the_window_of_neighbouring_columns_is_mixed_in_as_the_paper_says() {
    assert_scores(
        &scores(&protein(), true, &Settings::default()),
        &[
            0.823587, 0.770024, 0.777348, 0.744631, 0.817347, 0.759687, 0.754804, 0.818629,
            0.819930, 0.781874, 0.797682, 0.734655, 0.798328, 0.810466, 0.755885, 0.796236,
            0.777048, 0.840400, 0.828794, 0.737277, 0.829436, 0.795339,
        ],
    );
}

#[test]
fn the_entropy_of_the_protein_example_matches_the_reference() {
    assert_scores(
        &scores(&protein(), true, &plain(Method::Entropy)),
        &[
            0.999989, 0.768011, 0.867644, 0.869669, 0.999989, 0.856107, 0.855801, 0.999989,
            0.999989, 0.875812, 0.877424, 0.767308, 0.999989, 0.999989, 0.855801, 0.856107,
            0.869669, 0.999989, 0.999989, 0.731994, 0.999989, 0.855801,
        ],
    );
}

#[test]
fn without_weights_and_gap_penalty_the_scores_are_those_of_the_plain_column() {
    let settings = Settings {
        weights: false,
        gap_penalty: false,
        max_gaps: 1.0,
        ..plain(Method::Jsd)
    };
    assert_scores(
        &scores(&protein(), true, &settings),
        &[
            0.918022, 0.717735, 0.766365, 0.704730, 0.892278, 0.728996, 0.712351, 0.842276,
            0.892278, 0.777112, 0.892276, 0.737725, 0.835983, 0.868614, 0.719084, 0.771175,
            0.761660, 0.915333, 0.868614, 0.738989, 0.853011, 0.812318,
        ],
    );
}

#[test]
fn nucleotides_are_scored_against_an_even_background() {
    assert_scores(
        &scores(&dna(), false, &plain(Method::Jsd)),
        &[
            0.548792, 0.548792, 0.548792, 0.312026, 0.368437, 0.548792, 0.548792, 0.368437,
            0.368437, 0.548792,
        ],
    );
    assert_scores(
        &scores(&dna(), false, &Settings::default()),
        &[
            0.509331, 0.496652, 0.507080, 0.415379, 0.423854, 0.483972, 0.483972, 0.422543,
            0.436070, 0.488673,
        ],
    );
    assert_scores(
        &scores(&dna(), false, &plain(Method::Entropy)),
        &[
            0.999996, 0.999996, 0.999996, 0.550051, 0.657240, 0.999996, 0.999996, 0.657240,
            0.657240, 0.999996,
        ],
    );
}

#[test]
fn a_column_with_too_many_gaps_has_no_score_and_does_not_lend_one_to_its_neighbours() {
    // the middle column is a gap in three of eight rows
    let a = rows(&["MKT", "MKT", "MKT", "MKT", "MKT", "M-T", "M-T", "M-T"]);
    let s = scores(&a, true, &Settings::default());
    assert!(s[0].is_some() && s[1].is_none() && s[2].is_some(), "{s:?}");
    // the cutoff counts rows, not weights: three of eight is more than 30 % whatever the rows weigh
    let by_row = Settings {
        weights: false,
        ..Settings::default()
    };
    assert!(scores(&a, true, &by_row)[1].is_none());
    // it can be asked for: and then it is lowered by its gaps, and by nothing else
    let lenient = |gap_penalty| Settings {
        max_gaps: 1.0,
        gap_penalty,
        ..plain(Method::Jsd)
    };
    let lowered = scores(&a, true, &lenient(true));
    let whole = scores(&a, true, &lenient(false));
    assert!(
        lowered[1].expect("a score") < whole[1].expect("a score"),
        "{lowered:?} {whole:?}"
    );
    assert_eq!(lowered[0], whole[0], "a column with no gap is not lowered");
    // all gaps, one row, and no rows: nothing to score
    assert_eq!(
        scores(&rows(&["-", "-"]), true, &Settings::default()),
        vec![None]
    );
    assert_eq!(
        scores(&rows(&["MKT"]), true, &Settings::default()),
        vec![None; 3]
    );
    assert!(scores(&[], true, &Settings::default()).is_empty());
}

#[test]
fn letters_that_are_not_residues_count_as_gaps() {
    // X, B, Z and the like have no place in the alphabet, and the stop sign is no residue
    let a = rows(&[
        "MKT", "MKT", "MXT", "MKT", "MKT", "MKT", "MKT", "MKT", "MKT", "MKT",
    ]);
    let b = rows(&[
        "MKT", "MKT", "M-T", "MKT", "MKT", "MKT", "MKT", "MKT", "MKT", "MKT",
    ]);
    assert_eq!(
        scores(&a, true, &Settings::default()),
        scores(&b, true, &Settings::default())
    );
}

#[test]
fn a_conserved_column_scores_higher_than_a_variable_one_whatever_the_method() {
    // twenty rows: column 0 all L, column 1 all twenty different residues
    let letters = "ARNDCQEGHILKMFPSTVWY";
    let a: Vec<Vec<u8>> = letters.bytes().map(|c| vec![b'L', c]).collect();
    for method in [Method::Jsd, Method::Entropy] {
        let s = scores(&a, true, &plain(method));
        let (conserved, variable) = (s[0].expect("a score"), s[1].expect("a score"));
        // (a conserved column of a common residue such as L scores about 0.77 by JSD: a rare residue scores more)
        assert!(conserved > 0.75 && variable < 0.2, "{method:?}: {s:?}");
    }
}

#[test]
fn blocks_are_runs_of_conserved_columns_that_a_column_without_a_score_ends() {
    let s = [
        Some(0.9),
        Some(0.8),
        Some(0.85),
        Some(0.2),
        Some(0.95),
        Some(0.9),
        None,
        Some(0.9),
        Some(0.9),
        Some(0.9),
        Some(0.9),
    ];
    let found = blocks(&s, 0.7, 3);
    assert_eq!(found.len(), 2, "{found:?}");
    assert_eq!((found[0].from, found[0].to), (0, 2));
    assert!((found[0].mean - 0.85).abs() < 1e-6);
    assert_eq!((found[1].from, found[1].to), (7, 10));
    // shorter runs are left out, and nothing conserved gives none
    assert_eq!(blocks(&s, 0.7, 4).len(), 1);
    assert!(blocks(&s, 0.99, 1).is_empty());
    assert!(blocks(&[], 0.5, 1).is_empty());
}

#[test]
fn a_quarter_of_the_scored_columns_reach_the_top_quarter_score() {
    let s: Vec<Option<f32>> = (1..=8)
        .map(|k| Some(k as f32 / 10.0))
        .chain([None, None])
        .collect();
    let top = top_quarter(&s).expect("a score");
    let reaching = s.iter().flatten().filter(|&&v| v >= top).count();
    assert_eq!(reaching, 2, "{top}");
    assert_eq!(top_quarter(&[None, None]), None);
    assert_eq!(top_quarter(&[]), None);
}

#[test]
fn the_residues_of_a_row_are_in_its_columns_that_are_not_gaps() {
    // the gaps, written either way, are passed over
    assert_eq!(columns_of(b"A-CDE"), [0, 2, 3, 4]);
    assert_eq!(columns_of(b".a--c."), [1, 4]);
    assert_eq!(columns_of(b"-----"), Vec::<usize>::new());
    assert_eq!(columns_of(b""), Vec::<usize>::new());
}

#[test]
fn a_structure_chain_is_lined_up_with_its_row_even_where_residues_are_missing() {
    let row = b"MKTAYIAKQRQISFVKSHFSRQ";
    // the same protein: each residue to its own
    let same = map_to_chain(row, row).expect("a mapping");
    assert_eq!(same.to_row, (0..row.len()).map(Some).collect::<Vec<_>>());
    assert!((same.identity - 1.0).abs() < 1e-12);
    // a chain that starts with a three-residue tag the alignment does not have, and whose residue I (the sixth of the row)
    // was not seen in the experiment
    let chain = b"GSHMKTAYAKQRQISFVKSHFSRQ";
    let m = map_to_chain(row, chain).expect("a mapping");
    assert_eq!(m.to_row[..3], [None, None, None]);
    assert_eq!(m.to_row[3], Some(0));
    // the chain goes from Y (the fifth of the row) straight to A: that A is the seventh of the row
    assert_eq!(m.to_row[3 + 4], Some(4));
    assert_eq!(m.to_row[3 + 5], Some(6), "{:?}", m.to_row);
    assert!(m.identity > 0.99);
    // a different protein is not mapped
    assert!(map_to_chain(row, b"GGGGGGGGGGGGGGGGGGGGGG").is_none());
    assert!(map_to_chain(row, b"").is_none());
    assert!(map_to_chain(b"", row).is_none());
}

#[test]
fn a_column_of_one_residue_scores_what_its_background_frequency_says_whatever_order_the_table_is_kept_in(
) {
    // Each residue's background frequency is written here BY NAME (the BLOSUM62 background as the authors' table gives it,
    // before it is made to sum to 1), and not taken from the program's own table: a table and an alphabet in two orders
    // (V, W and Y once swapped places) would pass every number that the reference makes, because both took the same
    // constants. For a column of one residue with background q the score is worked out by hand:
    //   JSD = 1/2 (log2(2 / (1 + q)) + q log2(2 q / (1 + q)) + (1 - q))
    let by_name = [
        ('A', 0.078),
        ('C', 0.024),
        ('K', 0.056),
        ('L', 0.092),
        ('M', 0.024),
        ('V', 0.072),
        ('W', 0.014),
        ('Y', 0.034),
    ];
    let settings = Settings {
        weights: false,
        ..plain(Method::Jsd)
    };
    for (aa, written) in by_name {
        let q: f64 = written / 1.002;
        let want = 0.5 * ((2.0 / (1.0 + q)).log2() + q * (2.0 * q / (1.0 + q)).log2() + (1.0 - q));
        let column: Vec<Vec<u8>> = (0..20).map(|_| vec![aa as u8]).collect();
        let got = scores(&column, true, &settings)[0].expect("a score") as f64;
        assert!((got - want).abs() < 2e-4, "{aa}: {got} against {want}");
    }
    // the consequence that a reader sees: the rarer the residue, the more a column of it says
    let one = |aa: char| {
        scores(
            &(0..20).map(|_| vec![aa as u8]).collect::<Vec<_>>(),
            true,
            &settings,
        )[0]
        .expect("a score")
    };
    assert!(one('W') > one('Y') && one('Y') > one('V') && one('V') > one('L'));
}

#[test]
fn a_score_exactly_at_the_threshold_is_in_the_block() {
    let s = [Some(0.7), Some(0.7), Some(0.7), Some(0.6999)];
    let found = blocks(&s, 0.7, 3);
    assert_eq!(found.len(), 1);
    assert_eq!((found[0].from, found[0].to), (0, 2));
}

#[test]
fn lower_case_letters_and_u_are_residues_like_their_capitals() {
    let upper = rows(&["ACGU", "ACGT", "ACGT", "ACGT", "ACGT", "ACGT"]);
    let lower = rows(&["acgu", "acgt", "ACGt", "acgT", "ACGT", "acgt"]);
    assert_eq!(
        scores(&upper, false, &Settings::default()),
        scores(&lower, false, &Settings::default())
    );
    let protein_upper = protein();
    let protein_lower: Vec<Vec<u8>> = protein_upper
        .iter()
        .map(|r| r.to_ascii_lowercase())
        .collect();
    assert_eq!(
        scores(&protein_upper, true, &Settings::default()),
        scores(&protein_lower, true, &Settings::default())
    );
}

#[test]
fn a_chain_is_matched_with_the_row_that_covers_most_of_it_not_the_one_that_is_most_alike_in_a_piece(
) {
    let chain = b"MKTAYIAKQRQISFVKSHFSRQMKTAYIAKQRQISFVKSHFSRQ";
    // a piece of 12 residues that is the chain exactly, and the whole chain with two residues changed
    let piece = &chain[..12];
    let mut whole = chain.to_vec();
    whole[5] = b'W';
    whole[30] = b'W';
    let (fragment, full) = (
        map_to_chain(piece, chain).expect("a mapping"),
        map_to_chain(&whole, chain).expect("a mapping"),
    );
    assert!(
        fragment.identity > full.identity,
        "the piece is the more alike"
    );
    assert!(
        full.same > fragment.same,
        "{} against {}",
        full.same,
        fragment.same
    );
    assert_eq!((full.matched, fragment.matched), (chain.len(), piece.len()));
}
