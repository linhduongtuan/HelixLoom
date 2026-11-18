use super::*;
use crate::seq::tree::Rng;

fn random_bases(rng: &mut Rng, length: usize) -> Vec<u8> {
    (0..length).map(|_| BASES[rng.below(4)]).collect()
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// `count` random sequences of `length` bases with `word` put in each at a random place, each base changed with chance `error`,
/// on a random strand. The places and strands, so that the search can be checked.
fn planted(
    word: &[u8],
    count: usize,
    length: usize,
    error: f64,
    seed: u64,
) -> (Vec<Vec<u8>>, Vec<(usize, bool)>) {
    let mut rng = Rng::new(seed);
    let mut seqs = Vec::new();
    let mut places = Vec::new();
    for _ in 0..count {
        let mut seq = random_bases(&mut rng, length);
        let mut copy: Vec<u8> = word
            .iter()
            .map(|&b| {
                if rng.unit() < error {
                    BASES[rng.below(4)]
                } else {
                    b
                }
            })
            .collect();
        let reverse = rng.below(2) == 1;
        if reverse {
            copy = reverse_complement(&copy);
        }
        let at = rng.below(length - word.len() + 1);
        seq[at..at + word.len()].copy_from_slice(&copy);
        seqs.push(seq);
        places.push((at, reverse));
    }
    (seqs, places)
}

fn mismatches(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

const EQUAL: [f64; 4] = [0.25; 4];

#[test]
fn a_motif_of_sites_has_the_consensus_of_its_sites_columns_that_add_up_to_one_and_the_information_it_says(
) {
    let sites: [&[u8]; 4] = [b"ACGT", b"ACGT", b"ACGA", b"ACGT"];
    let motif = Motif::from_sites("m", &sites).expect("a motif");
    assert_eq!(motif.consensus(), "ACGT");
    assert_eq!(motif.width(), 4);
    for column in &motif.probs {
        assert!((column.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(column.iter().all(|&p| p > 0.0));
    }
    // three columns that never change say nearly two bits each; the last, three to one, says less
    let bits = motif.information();
    assert!(bits > 5.0 && bits < 8.0, "{bits}");
    // columns of equal chances say nothing
    let flat = Motif {
        name: "flat".into(),
        probs: vec![[0.25; 4]; 5],
    };
    assert!(flat.information().abs() < 1e-12);
}

#[test]
fn sites_that_do_not_make_a_motif_are_refused_with_the_reason() {
    assert!(Motif::from_sites("m", &[]).is_err());
    assert!(Motif::from_sites("m", &[b"ACGT", b"ACG"]).is_err());
    let error = Motif::from_sites("m", &[b"ACNT"]).expect_err("not a base");
    assert!(error.contains('N'), "{error}");
}

#[test]
fn the_motif_of_the_other_strand_reads_as_the_reverse_complement_and_back_again() {
    let motif = Motif::from_sites("m", &[b"AACGT", b"AACGT", b"AACCT"]).expect("a motif");
    let other = motif.reverse_complement();
    assert_eq!(
        other.consensus().as_bytes(),
        reverse_complement(motif.consensus().as_bytes())
    );
    assert_eq!(other.reverse_complement(), motif);
}

#[test]
fn a_scan_finds_a_site_on_each_strand_where_it_is_and_skips_windows_with_unknown_letters() {
    let motif = Motif::from_sites("m", &[&b"GATTACA"[..]; 5]).expect("a motif");
    let site = b"GATTACA";
    let mut seq = b"CCCCCCCCCC".to_vec();
    seq.extend(site); // forward at 10
    seq.extend(b"TTTTTTTTTT");
    seq.extend(reverse_complement(site)); // reverse at 27
    seq.extend(b"CCCCCCCCCC");
    let hits = scan(&seq, &motif, &EQUAL, 8.0, true);
    assert_eq!(hits.len(), 2, "{hits:?}");
    assert_eq!((hits[0].position, hits[0].reverse), (10, false));
    assert_eq!((hits[1].position, hits[1].reverse), (27, true));
    assert!(hits.iter().all(|h| h.score > 8.0));
    // one strand only: only the forward site
    let forward = scan(&seq, &motif, &EQUAL, 8.0, false);
    assert_eq!(forward.len(), 1);
    // a letter that is not a base inside the window: that window is not scored
    let mut broken = seq.clone();
    broken[13] = b'N';
    assert_eq!(scan(&broken, &motif, &EQUAL, 8.0, true).len(), 1);
    // a sequence shorter than the motif has nothing
    assert!(scan(b"GATT", &motif, &EQUAL, -100.0, true).is_empty());
}

#[test]
fn the_p_value_of_a_score_is_the_share_of_all_words_that_score_so_much_and_the_threshold_for_a_p_value_is_the_lowest_that_fits(
) {
    // a motif of four columns, so that all 256 words can be scored one by one
    let motif =
        Motif::from_sites("m", &[b"ACGT", b"ACGA", b"ACCT", b"TCGT", b"ACGT"]).expect("a motif");
    let bg = [0.3, 0.2, 0.2, 0.3];
    let odds = motif.log_odds(&bg);
    let hundredths: Vec<[i64; 4]> = odds
        .iter()
        .map(|c| std::array::from_fn(|b| (c[b] * 100.0).round() as i64))
        .collect();
    let mut everything: Vec<(i64, f64)> = Vec::new();
    for word in 0..256usize {
        let letters: Vec<usize> = (0..4).map(|k| (word >> (2 * k)) & 3).collect();
        let score: i64 = letters.iter().zip(&hundredths).map(|(&b, c)| c[b]).sum();
        let chance: f64 = letters.iter().map(|&b| bg[b]).product();
        everything.push((score, chance));
    }
    let distribution = ScoreDistribution::new(&motif, &bg);
    for &(score, _) in &everything {
        let want: f64 = everything
            .iter()
            .filter(|(s, _)| *s >= score)
            .map(|(_, c)| c)
            .sum();
        let got = distribution.p_value(score as f64 / 100.0);
        assert!(
            (got - want).abs() < 1e-12,
            "score {score}: {got} against {want}"
        );
    }
    // everything scores at least the lowest, and nothing scores more than the highest
    assert!((distribution.p_value(-1000.0) - 1.0).abs() < 1e-12);
    assert_eq!(distribution.p_value(1000.0), 0.0);
    for p in [0.5, 0.1, 0.01, 0.001] {
        let threshold = distribution.score_for(p);
        assert!(distribution.p_value(threshold) <= p + 1e-12, "{p}");
        assert!(distribution.p_value(threshold - 0.01) > p, "{p}");
    }
}

#[test]
fn a_motif_planted_in_random_sequences_on_either_strand_is_found_where_it_was_put() {
    let word = b"ACCGTTAGCA";
    let (seqs, places) = planted(word, 40, 120, 0.1, 7);
    let found = discover(
        &seqs,
        &Settings {
            width: 10,
            ..Settings::default()
        },
    )
    .expect("a search");
    assert_eq!(found.len(), 1);
    let first = &found[0];
    // the consensus is the word, or its other strand's, give or take a base
    let (a, b) = (
        mismatches(first.motif.consensus().as_bytes(), word),
        mismatches(
            first.motif.consensus().as_bytes(),
            &reverse_complement(word),
        ),
    );
    assert!(
        a.min(b) <= 1,
        "{} against {}",
        first.motif.consensus(),
        String::from_utf8_lossy(word)
    );
    // nearly every planted site is found, at its place and on its strand
    let right = first
        .sites
        .iter()
        .filter(|s| {
            places[s.sequence] == (s.position, s.reverse) || {
                // a motif found as the other strand reads each site the other way round
                let (at, reverse) = places[s.sequence];
                at == s.position && reverse != s.reverse
            }
        })
        .count();
    assert!(right >= 36, "{right} of 40 sites right");
    assert!(first.gain > 20.0, "{}", first.gain);
    assert!(first.motif.information() > 10.0);
}

#[test]
fn two_motifs_in_the_same_sequences_are_found_one_after_the_other() {
    let (mut one, _) = planted(b"ACCGTTAGCA", 40, 100, 0.08, 11);
    let (two, _) = planted(b"TTGGCCAATG", 40, 100, 0.08, 12);
    // each sequence carries both words: the second word is put over a stretch of the first's sequence away from the first word
    for (a, b) in one.iter_mut().zip(&two) {
        let word = &b[..100];
        // take a random-looking piece of the second sequence that holds its word by using the whole sequence beside the first
        a.extend(word);
    }
    let found = discover(
        &one,
        &Settings {
            width: 10,
            motifs: 2,
            ..Settings::default()
        },
    )
    .expect("a search");
    assert_eq!(found.len(), 2);
    let near = |motif: &Motif, word: &[u8]| {
        mismatches(motif.consensus().as_bytes(), word).min(mismatches(
            motif.consensus().as_bytes(),
            &reverse_complement(word),
        ))
    };
    let (a, b) = (&found[0].motif, &found[1].motif);
    let both = (near(a, b"ACCGTTAGCA") <= 2 && near(b, b"TTGGCCAATG") <= 2)
        || (near(a, b"TTGGCCAATG") <= 2 && near(b, b"ACCGTTAGCA") <= 2);
    assert!(both, "{} and {}", a.consensus(), b.consensus());
}

#[test]
fn a_search_that_cannot_be_made_says_why() {
    let short = vec![b"ACGT".to_vec(), b"AC".to_vec()];
    assert!(discover(
        &short,
        &Settings {
            width: 8,
            ..Settings::default()
        }
    )
    .is_err());
    assert!(discover(
        &short,
        &Settings {
            width: 2,
            ..Settings::default()
        }
    )
    .is_err());
    assert!(discover(
        &short,
        &Settings {
            width: 41,
            ..Settings::default()
        }
    )
    .is_err());
    // unknown letters only: nothing to start from, so no motif, and no error either
    let unknown = vec![b"NNNNNNNNNNNN".to_vec(); 3];
    assert!(discover(&unknown, &Settings::default())
        .expect("a search")
        .is_empty());
}

#[test]
fn motifs_are_compared_by_correlation_on_either_strand_and_at_any_offset() {
    let a = Motif::from_sites("a", &[&b"GGATTACAGG"[..]; 6]).expect("a motif");
    assert!((compare(&a, &a, 6) - 1.0).abs() < 1e-9);
    // the other strand of the same motif
    assert!((compare(&a, &a.reverse_complement(), 6) - 1.0).abs() < 1e-9);
    // the same motif with flanks cut off, found at its offset
    let inner = Motif::from_sites("b", &[&b"ATTACA"[..]; 6]).expect("a motif");
    assert!(compare(&a, &inner, 6) > 0.99);
    // unrelated
    let other = Motif::from_sites("c", &[&b"CCCCCCCCCC"[..]; 6]).expect("a motif");
    assert!(compare(&a, &other, 6) < 0.4);
}

#[test]
fn a_motif_is_written_in_the_text_format_of_meme_and_reads_back_from_it() {
    let motif = Motif::from_sites("my motif", &[&b"GATTACA"[..]; 4]).expect("a motif");
    let text = motif.to_meme(&EQUAL, 4);
    assert!(text.starts_with("MEME version 4\n"));
    assert!(text.contains("MOTIF my_motif\n"), "{text}");
    assert!(text.contains("alength= 4 w= 7 nsites= 4"), "{text}");
    let rows: Vec<Vec<f64>> = text
        .lines()
        .skip_while(|l| !l.starts_with("letter-probability"))
        .skip(1)
        .map(|l| {
            l.split_whitespace()
                .map(|x| x.parse().expect("a number"))
                .collect()
        })
        .collect();
    assert_eq!(rows.len(), 7);
    for (row, column) in rows.iter().zip(&motif.probs) {
        assert!((row.iter().sum::<f64>() - 1.0).abs() < 1e-5);
        for (a, b) in row.iter().zip(column) {
            assert!((a - b).abs() < 1e-6);
        }
    }
}

#[test]
fn a_small_set_is_searched_in_full_and_a_big_one_less_so_that_the_work_stays_bounded() {
    let settings = Settings::default();
    // 12 sequences of 90 bases: every start, every seed, every iteration
    let small = plan(2.0 * 12.0 * 80.0 * 10.0, &settings);
    assert_eq!(small, (3000, settings.seeds, settings.iterations));
    // the work grows with the set, and the plan keeps the number of passes about what the budget allows
    let mut last = usize::MAX;
    for windows in [1e4, 1e5, 1e6, 1e7, 1e8] {
        let (candidates, seeds, iterations) = plan(windows * 10.0, &settings);
        let passes = (candidates + seeds * iterations) as f64;
        assert!(
            passes <= 1.0 + WORK / (windows * 10.0) + 20.0 + 5.0,
            "{windows}: {passes} passes"
        );
        assert!(candidates <= last);
        last = candidates;
        assert!(seeds >= 1 && iterations >= 5 && candidates >= 10);
    }
}

#[test]
fn a_set_of_hundreds_of_thousands_of_bases_is_searched_in_seconds_and_still_finds_a_strong_motif() {
    // 20 sequences of 5,000 bases with a motif in each: 200,000 windows on both strands
    let word = b"ACCGTTAGCA";
    let (seqs, _) = planted(word, 20, 5000, 0.05, 31);
    let start = std::time::Instant::now();
    let found = discover(
        &seqs,
        &Settings {
            width: 10,
            ..Settings::default()
        },
    )
    .expect("a search");
    assert!(
        start.elapsed() < std::time::Duration::from_secs(60),
        "{:?}",
        start.elapsed()
    );
    let consensus = found[0].motif.consensus();
    let off = |w: &[u8]| {
        mismatches(consensus.as_bytes(), w)
            .min(mismatches(consensus.as_bytes(), &reverse_complement(w)))
    };
    assert!(off(word) <= 2, "{consensus}");
}

#[test]
fn more_effort_buys_more_passes_over_a_big_set_and_a_small_set_is_searched_in_full_either_way() {
    let base = Settings::default();
    let (more, less) = (
        Settings {
            effort: 10.0,
            ..base
        },
        Settings {
            effort: 0.1,
            ..base
        },
    );
    let passes = |cost: f64, s: &Settings| {
        let (candidates, seeds, iterations) = plan(cost, s);
        candidates + seeds * iterations
    };
    // a big set: 10 million windows of 10 columns
    let big = 1e8;
    assert_eq!(plan(big, &base), (10, 1, 5), "the budget of the window");
    assert!(passes(big, &more) >= 3 * passes(big, &base));
    assert!(passes(big, &less) <= passes(big, &base));
    // 10 times the effort over a set 10 times as big is the same plan as the default over the small one
    assert_eq!(plan(big, &more), plan(big / 10.0, &base));
    // a small set gets every start, seed and iteration whatever the effort
    let small = 2.0 * 12.0 * 80.0 * 10.0;
    assert_eq!(plan(small, &base), plan(small, &more));
}

#[test]
fn an_effort_that_is_no_number_or_out_of_range_is_refused() {
    let (seqs, _) = planted(b"ACCGTTAGCA", 6, 60, 0.0, 3);
    for effort in [0.0, -1.0, 0.001, 5000.0, f64::NAN, f64::INFINITY] {
        let said = discover(
            &seqs,
            &Settings {
                effort,
                ..Settings::default()
            },
        )
        .err()
        .unwrap_or_else(|| panic!("{effort} should be refused"));
        assert!(said.contains("effort"), "{said}");
    }
    for effort in [0.01, 1.0, 1000.0] {
        assert!(discover(
            &seqs,
            &Settings {
                effort,
                width: 10,
                ..Settings::default()
            }
        )
        .is_ok());
    }
}

#[test]
fn several_motifs_are_one_file_in_the_text_format_of_meme_with_one_header() {
    let a = Motif::from_sites("first", &[&b"GATTACA"[..]; 3]).expect("a motif");
    let b = Motif::from_sites("second one", &[&b"TTGACC"[..]; 3]).expect("a motif");
    let site = |sequence| Site {
        sequence,
        position: 0,
        reverse: false,
        score: 5.0,
        probability: 0.9,
    };
    let found = [
        Found {
            motif: a.clone(),
            sites: vec![site(0), site(1), site(2)],
            gain: 12.0,
        },
        Found {
            motif: b,
            sites: vec![site(0), site(2)],
            gain: 7.0,
        },
    ];
    let text = meme_text(&found, &EQUAL);
    assert_eq!(text.matches("MEME version 4").count(), 1, "{text}");
    assert_eq!(text.matches("ALPHABET= ACGT").count(), 1);
    assert_eq!(text.matches("Background letter frequencies").count(), 1);
    assert_eq!(text.matches("MOTIF ").count(), 2);
    assert!(text.contains("MOTIF first\n") && text.contains("MOTIF second_one\n"));
    // each motif says how many sites it was made from
    assert!(text.contains("w= 7 nsites= 3"), "{text}");
    assert!(text.contains("w= 6 nsites= 2"), "{text}");
    // one motif alone is what `to_meme` writes
    assert_eq!(meme_text(&found[..1], &EQUAL), a.to_meme(&EQUAL, 3));
    // no motif is the header alone
    assert_eq!(
        meme_text(&[], &EQUAL),
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.2500 C 0.2500 G 0.2500 T 0.2500\n"
    );
}

#[test]
fn sequences_are_read_from_fasta_or_one_to_a_line_and_two_are_the_least() {
    let (names, seqs) = read_sequences(">one a note\nACGT\nAC\n>two\nGGGG\n").expect("fasta");
    assert_eq!(names, ["one", "two"]);
    assert_eq!(seqs, [b"ACGTAC".to_vec(), b"GGGG".to_vec()]);
    let (names, seqs) = read_sequences("ac gt\n\n123 TTT\n").expect("lines");
    assert_eq!(names, ["sequence 1", "sequence 2"]);
    assert_eq!(seqs, [b"acgt".to_vec(), b"TTT".to_vec()]);
    for one in ["", ">only\nACGT\n", "ACGT\n"] {
        let said = read_sequences(one).unwrap_err();
        assert!(said.contains("at least two sequences"), "{said}");
    }
}

#[test]
fn a_search_whose_windows_would_not_fit_in_memory_is_refused_before_any_is_made() {
    // 16 million bases at 8 columns on both strands: 32 million windows of 128 bytes, about 3.8 GB, which is taken
    assert!(check_size(16_000_000, 8, true).is_ok());
    // a few more, or a wider motif, is not; one strand is half
    assert!(check_size(17_000_000, 8, true).is_err());
    assert!(check_size(16_000_000, 40, true).is_err());
    assert!(check_size(33_000_000, 8, false).is_ok());
    assert!(check_size(34_000_000, 8, false).is_err());
    let said = check_size(60_000_000, 40, true).unwrap_err();
    assert!(
        said.contains("60000000 bases") && said.contains("at most 4 GB"),
        "{said}"
    );
    // and the search itself says so without holding them: 3 sequences of 30 million bases at 40 columns
    let big = vec![vec![b'A'; 30_000_000]; 3];
    let said = discover(
        &big,
        &Settings {
            width: 40,
            ..Settings::default()
        },
    )
    .unwrap_err();
    assert!(said.contains("use fewer bases"), "{said}");
}
