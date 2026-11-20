use super::*;

fn named(rows: &[(&str, &str)]) -> Vec<(String, Vec<u8>)> {
    rows.iter()
        .map(|(n, s)| (n.to_string(), s.as_bytes().to_vec()))
        .collect()
}

fn dna(t: &Translated, rows: &[(&str, &str)]) -> Vec<(String, String)> {
    let aligned = named(rows);
    back_translate(t, &aligned)
        .unwrap()
        .into_iter()
        .map(|(n, r)| (n, String::from_utf8(r).unwrap()))
        .collect()
}

#[test]
fn codons_become_amino_acids_and_a_stop_or_an_ambiguous_base_becomes_x() {
    let t = translate_all(
        &named(&[("a", "ATGGCCAAATAA"), ("b", "atgnccaaatga")]),
        &Codons::default(),
    )
    .unwrap();
    assert_eq!(t.proteins[0].1, b"MAKX");
    // the ambiguous codon is X, a stop is X; lower case is the same as upper case
    assert_eq!(t.proteins[1].1, b"MXKX");
    assert_eq!(t.proteins[0].0, "a");
}

#[test]
fn the_reading_frame_leaves_a_head_and_a_tail_and_the_genetic_code_changes_the_protein() {
    let seqs = named(&[("a", "GATGGCCAA")]);
    let t = translate_all(&seqs, &Codons { frame: 1, code: 0 }).unwrap();
    assert_eq!(t.proteins[0].1, b"MA");
    assert_eq!(t.cuts[0].head, b"G");
    // GATGGCCAA: the G, the codons ATG and GCC, and the two bases AA
    assert_eq!(t.cuts[0].tail, b"AA");
    // TGA is a stop in the standard code and tryptophan in the vertebrate mitochondrial code (the second in the list)
    let tga = named(&[("a", "ATGTGAAAA")]);
    assert_eq!(
        translate_all(&tga, &Codons { frame: 0, code: 0 })
            .unwrap()
            .proteins[0]
            .1,
        b"MXK"
    );
    assert_eq!(
        translate_all(&tga, &Codons { frame: 0, code: 1 })
            .unwrap()
            .proteins[0]
            .1,
        b"MWK"
    );
    // gaps and spaces in what was given are not part of the sequence, and U is T
    let spaced = named(&[("a", "AUG-GCC AAA")]);
    assert_eq!(
        translate_all(&spaced, &Codons::default()).unwrap().proteins[0].1,
        b"MAK"
    );
}

#[test]
fn the_alignment_of_the_proteins_is_carried_back_as_whole_codons() {
    let t = translate_all(
        &named(&[("a", "ATGGCCAAA"), ("b", "ATGAAA")]),
        &Codons::default(),
    )
    .unwrap();
    assert_eq!(
        dna(&t, &[("a", "MAK"), ("b", "M-K")]),
        [
            ("a".into(), "ATGGCCAAA".into()),
            ("b".into(), "ATG---AAA".into())
        ]
    );
    // gaps at the ends are codons too
    assert_eq!(
        dna(&t, &[("a", "-MAK"), ("b", "MK--")])[1].1,
        "ATGAAA------"
    );
}

#[test]
fn heads_line_up_and_tails_are_padded_to_one_length_at_the_end() {
    let t = translate_all(
        &named(&[("a", "GATGGCCAA"), ("b", "GATGGCCA"), ("c", "GATGGCC")]),
        &Codons { frame: 1, code: 0 },
    )
    .unwrap();
    let rows = dna(&t, &[("a", "MA"), ("b", "MA"), ("c", "MA")]);
    assert_eq!(
        rows.iter().map(|r| r.1.as_str()).collect::<Vec<_>>(),
        ["GATGGCCAA", "GATGGCCA-", "GATGGCC--"]
    );
}

#[test]
fn what_the_aligner_changes_or_leaves_out_is_refused_with_the_sequence_named() {
    let t = translate_all(
        &named(&[("a", "ATGGCCAAA"), ("b", "ATGAAA")]),
        &Codons::default(),
    )
    .unwrap();
    let why = |rows: &[(&str, &str)]| back_translate(&t, &named(rows)).unwrap_err();
    assert!(why(&[("a", "MAK")]).contains("returned 1 sequences for 2"));
    assert!(
        why(&[("a", "MAK"), ("b", "M-")]).contains("not all of one length")
            || why(&[("a", "MAK"), ("b", "M-")]).contains("left out part of the sequence b")
    );
    assert!(why(&[("a", "MAK"), ("b", "M-KK")]).contains("one length"));
    assert!(why(&[("a", "MAK"), ("b", "MA-")]).contains("changed the sequence b"));
    assert!(why(&[("a", "MAKK"), ("b", "M-K-")]).contains("changed the sequence a"));
    assert!(why(&[("a", "MA-"), ("b", "M-K")]).contains("left out part of the sequence a"));
    // an aligner that writes lower case is not changing it
    assert_eq!(dna(&t, &[("a", "mak"), ("b", "m-k")])[1].1, "ATG---AAA");
}

#[test]
fn a_choice_that_is_none_a_protein_and_a_sequence_without_a_codon_are_refused() {
    let one = named(&[("a", "ATGGCC")]);
    assert!(translate_all(&one, &Codons { frame: 3, code: 0 })
        .unwrap_err()
        .contains("1, 2 or 3"));
    assert!(translate_all(&one, &Codons { frame: 0, code: 99 })
        .unwrap_err()
        .contains("genetic code"));
    assert!(
        translate_all(&named(&[("p", "MKVLAAGIVGLLL")]), &Codons::default())
            .unwrap_err()
            .contains("does not look like DNA")
    );
    assert!(translate_all(
        &named(&[("a", "ATGGCC"), ("short", "AT")]),
        &Codons::default()
    )
    .unwrap_err()
    .contains("short has no whole codon in reading frame 1"));
    assert!(
        translate_all(&named(&[("a", "ATGG")]), &Codons { frame: 2, code: 0 })
            .unwrap_err()
            .contains("reading frame 3")
    );
    assert!(translate_all(&named(&[("empty", "")]), &Codons::default()).is_err());
}

#[test]
fn the_choice_is_written_in_words_and_read_back() {
    for frame in 0..3 {
        for code in 0..GENETIC_CODES.len() {
            let c = Codons { frame, code };
            assert_eq!(Codons::parse(&c.describe()), Some(c), "{}", c.describe());
        }
    }
    assert_eq!(
        Codons { frame: 1, code: 0 }.describe(),
        "reading frame 2, genetic code 1 (Standard)"
    );
    for bad in [
        "",
        "reading frame 4, genetic code 1 (Standard)",
        "reading frame 0, genetic code 1 (Standard)",
        "reading frame 1, genetic code 99 (x)",
        "frame 1",
    ] {
        assert_eq!(Codons::parse(bad), None, "{bad}");
    }
}

/// Random coding sequences, and a stand-in aligner that puts gaps into the proteins here and there.
fn random_case(seed: u32) -> (Vec<(String, Vec<u8>)>, Codons) {
    let mut x = seed.wrapping_mul(2_654_435_761).wrapping_add(1);
    let mut next = move |n: u32| {
        x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
        (x >> 16) % n
    };
    let codons = Codons {
        frame: next(3) as usize,
        code: next(2) as usize,
    };
    let rows = (0..2 + next(5))
        .map(|k| {
            let n = 4 + next(30) as usize;
            let mut s: Vec<u8> = (0..codons.frame)
                .map(|_| b"ACGT"[next(4) as usize])
                .collect();
            s.extend(
                (0..3 * n + next(3) as usize)
                    .map(|_| b"ACGTN"[next(5).min(3) as usize + usize::from(next(40) == 0)]),
            );
            (format!("s{k}"), s)
        })
        .collect();
    (rows, codons)
}

fn fake_alignment(proteins: &[(String, Vec<u8>)], seed: u32) -> Vec<(String, Vec<u8>)> {
    let mut x = seed.wrapping_add(7);
    let mut next = move |n: usize| {
        x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
        ((x >> 16) as usize) % n
    };
    let width = proteins.iter().map(|p| p.1.len()).max().unwrap_or(0) + 4;
    proteins
        .iter()
        .map(|(name, protein)| {
            // the residues go to increasing places among `width` columns
            let mut places: Vec<usize> = (0..width).collect();
            while places.len() > protein.len() {
                places.remove(next(places.len()));
            }
            let mut row = vec![b'-'; width];
            for (&p, &r) in places.iter().zip(protein) {
                row[p] = r;
            }
            (name.clone(), row)
        })
        .collect()
}

#[test]
fn for_random_sequences_no_base_is_lost_or_moved_and_every_gap_is_a_whole_codon() {
    for seed in 0..300 {
        let (rows, codons) = random_case(seed);
        let t = translate_all(&rows, &codons).unwrap();
        let aligned = fake_alignment(&t.proteins, seed);
        let width = aligned[0].1.len();
        let out = back_translate(&t, &aligned).unwrap();
        for ((name, original), (_, row)) in rows.iter().zip(&out) {
            // without its gaps the row is the sequence that was given
            let ungapped: Vec<u8> = row.iter().copied().filter(|&b| b != b'-').collect();
            assert_eq!(&ungapped, original, "seed {seed} {name}");
            // the codons of the alignment sit at the same places in every row, and a gap is three
            let head = codons.frame;
            for col in 0..width {
                let cell = &row[head + 3 * col..head + 3 * col + 3];
                let gaps = cell.iter().filter(|&&b| b == b'-').count();
                assert!(
                    gaps == 0 || gaps == 3,
                    "seed {seed} {name} column {col}: {cell:?}"
                );
            }
        }
        // every row has one length
        assert!(
            out.iter().all(|r| r.1.len() == out[0].1.len()),
            "seed {seed}"
        );
    }
}
