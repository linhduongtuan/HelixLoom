use super::*;

fn rows(of: &[(&str, &str)]) -> Vec<(String, Vec<u8>)> {
    of.iter()
        .map(|(n, s)| (n.to_string(), s.as_bytes().to_vec()))
        .collect()
}

fn marks(of: &[(&str, &str)]) -> String {
    conservation(&rows(of))
}

/// Residues that no nucleotide code is, so that rows of a few letters are read as protein.
const FILLER: &str = "EEEEEEEEEE";

#[test]
fn nucleotides_get_stars_where_every_row_agrees_and_blanks_otherwise() {
    // a column with a gap is blank even when the others agree, and so is one where the bases differ
    assert_eq!(
        marks(&[("a", "ACGT-AC-"), ("b", "ACGA-A-C"), ("c", "ACGTTAC-")]),
        "***  *  "
    );
    // one row: every column is conserved
    assert_eq!(marks(&[("a", "ACGT")]), "****");
    // upper and lower case agree
    assert_eq!(marks(&[("a", "acgt"), ("b", "ACGT")]), "****");
    // a strong group means nothing for nucleotides
    assert_eq!(marks(&[("a", "ACGT"), ("b", "TGCA")]), "    ");
}

#[test]
fn proteins_get_clustalws_strong_and_weak_groups_and_a_gap_blanks_the_column() {
    // columns: identical | STA | NEQ | FYW | CSA | FVL | H Y and a gap | identical | D W | K R A | all gaps
    let got = marks(&[
        ("a", "MSNFCFHADK-"),
        ("b", "MTEYSVYAWR-"),
        ("c", "MAQWAL.ADA-"),
    ]);
    assert_eq!(got, "*:::.. *   ");
}

#[test]
fn a_residue_in_both_a_strong_and_a_weak_group_gets_the_strong_mark_and_a_group_must_hold_every_residue(
) {
    // S and A are in STA (strong) and in CSA (weak)
    let pad = |s: &str| format!("{s}{FILLER}");
    let star = "*".repeat(FILLER.len());
    assert_eq!(
        marks(&[("a", &pad("SW")), ("b", &pad("AW"))]),
        format!(":*{star}")
    );
    // most of a column in a group is not all of it
    assert_eq!(
        marks(&[("a", &pad("SSS")), ("b", &pad("SSA")), ("c", &pad("SSW"))]),
        format!("** {star}")
    );
    // weak only: C S A is in CSA, and in no strong group
    assert_eq!(
        marks(&[("a", &pad("C")), ("b", &pad("S")), ("c", &pad("A"))]),
        format!(".{star}")
    );
}

#[test]
fn a_text_is_written_in_blocks_of_sixty_with_counts_and_marks_under_each() {
    let a: String = "ACGT".repeat(32); // 128 columns
    let mut b = a.clone();
    b.replace_range(10..12, "--"); // two gaps in the first block
    b.replace_range(100..101, "T"); // a difference in the third
    let text = write(&rows(&[("first one", &a), ("second", &b)])).expect("written");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "CLUSTAL format alignment by HelixLoom");
    assert!(lines[1].is_empty() && lines[2].is_empty());
    // three blocks, each a blank line, two rows and a line of marks
    assert_eq!(lines.len(), 2 + 3 * 4);
    let name_width = "first_one".len().max(NAME_WIDTH) + 1;
    for (k, (count_a, count_b)) in [(60, 58), (120, 118), (128, 126)].into_iter().enumerate() {
        let at = 3 + k * 4;
        assert!(lines[at - 1].is_empty(), "{k}");
        let (first, second, under) = (lines[at], lines[at + 1], lines[at + 2]);
        assert!(first.starts_with("first_one "), "{first}");
        assert!(second.starts_with("second "), "{second}");
        assert!(first.ends_with(&format!(" {count_a}")), "{first}");
        assert!(second.ends_with(&format!(" {count_b}")), "{second}");
        // the marks lie under the letters
        assert!(under.starts_with(&" ".repeat(name_width)));
        assert!(under[name_width..].chars().count() <= BLOCK);
    }
    // a gap of two columns in the second row is blank in the marks of the first block
    let first_marks = &lines[5][name_width..];
    assert_eq!(&first_marks[10..12], "  ");
    assert_eq!(&first_marks[0..10], "**********");
    // a name with white space has underscores in it, and a name with nothing in it has a name
    let odd = write(&rows(&[(" a \t b ", "AC"), ("", "AC")])).expect("written");
    assert!(odd.lines().any(|l| l.starts_with("a_b ")), "{odd}");
    assert!(odd.lines().any(|l| l.starts_with("sequence ")), "{odd}");
}

#[test]
fn what_is_written_is_read_back_whole() {
    let a: String = "ACGTTGCAAC".repeat(14); // 140 columns, three blocks
    let mut b = a.clone();
    b.replace_range(0..3, "---");
    b.replace_range(139..140, "-");
    let rows = rows(&[("one", &a), ("two", &b), ("three", &a)]);
    let text = write(&rows).expect("written");
    assert!(is_clustal(&text));
    assert_eq!(read(&text), Ok(rows.clone()));
    // a protein as well
    let protein = self::rows(&[("p", "MKVLAAGIVGLLLAQP"), ("q", "MKV-AAGIVGLL-AQP")]);
    assert_eq!(read(&write(&protein).expect("written")), Ok(protein));
}

#[test]
fn a_file_of_clustalw_with_its_version_marks_counts_and_blank_lines_is_read() {
    let text = [
        "CLUSTAL W (1.83) multiple sequence alignment",
        "",
        "",
        "seq_a           ACGT-ACGTA     9",
        "seq_b           ACGTTACGTA    10",
        "                **** *****",
        "",
        "seq_a           CCGG   13",
        "seq_b           CCGA   14",
        "                *** ",
        "",
    ]
    .join("\n");
    let text = text.as_str();
    let found = read(text).expect("read");
    assert_eq!(
        found,
        rows(&[("seq_a", "ACGT-ACGTACCGG"), ("seq_b", "ACGTTACGTACCGA")])
    );
}

#[test]
fn files_of_other_programs_are_read_with_no_counts_with_no_marks_and_with_their_own_headers() {
    for header in [
        "CLUSTAL format alignment by MAFFT FFT-NS-2 (v7.526)",
        "CLUSTAL O(1.2.4) multiple sequence alignment",
        "MUSCLE (3.8) multiple sequence alignment",
        "clustal w multiple sequence alignment",
    ] {
        let text = format!("{header}\n\n\nx  AC-GT\ny  ACAGT\n\nx  TT\ny  TA\n");
        assert_eq!(
            read(&text),
            Ok(rows(&[("x", "AC-GTTT"), ("y", "ACAGTTA")])),
            "{header}"
        );
    }
    // windows line ends
    let text = "CLUSTAL W\r\n\r\nx  ACGT\r\ny  ACGA\r\n\r\n";
    assert_eq!(read(text), Ok(rows(&[("x", "ACGT"), ("y", "ACGA")])));
}

#[test]
fn a_text_that_is_not_a_clustal_alignment_is_refused_in_words() {
    assert!(!is_clustal("") && !is_clustal(">a\nACGT\n") && !is_clustal("x ACGT\n"));
    assert!(is_clustal("\n\n  CLUSTAL W\n"));
    // other programs' headers, one that stands on its own over rows, and a byte order mark in front
    assert!(is_clustal(
        "Kalign (3.4) alignment in ClustalW format\n\nx ACGT\n"
    ));
    assert!(is_clustal("PROBCONS multiple sequence alignment\n"));
    assert!(is_clustal(
        "\u{feff}CLUSTAL W (1.83) multiple sequence alignment\n"
    ));
    assert!(is_clustal("CLUSTAL\n\nx  ACGT\ny  ACGA\n"));
    // the letters of a protein are not a header, however they spell
    assert!(
        !is_clustal("CLUSTAL\nKKKK\n")
            && !is_clustal("MUSCLEMUSCLE\n")
            && !is_clustal("CLUSTALW KKKK")
    );
    for (text, why) in [
        ("", "first line should say CLUSTAL"),
        (">a\nACGT\n", "first line should say CLUSTAL"),
        ("CLUSTAL W\n\n", "no sequences"),
        // rows that end up of different lengths
        ("CLUSTAL W\n\nx ACGT\ny ACG\n", "not of one length"),
        // a later block with the names in another order or with another name
        ("CLUSTAL W\n\nx AC\ny AC\n\ny GT\nx GT\n", "Line 6"),
        ("CLUSTAL W\n\nx AC\ny AC\n\nx GT\nz GT\n", "Line 7"),
        // a row in a later block that the first did not have
        (
            "CLUSTAL W\n\nx AC\ny AC\n\nx GT\ny GT\nz GT\n",
            "not one of the sequences of the first block",
        ),
        // not residues
        (
            "CLUSTAL W\n\nx AC#T\ny ACGT\n",
            "neither a residue nor a gap",
        ),
        (
            "CLUSTAL W\n\nx\ny ACGT\n",
            "expected a name and its residues",
        ),
        ("CLUSTAL W\n\nx A C G\n", "expected a name and its residues"),
        // a third word that is not a count
        (
            "CLUSTAL W\n\nx ACGT foo\ny ACGA\n",
            "expected a name and its residues",
        ),
        // an indented line that is not marks
        (
            "CLUSTAL W\n\nx ACGT\n   some text\n",
            "marks of conservation",
        ),
    ] {
        let said = read(text).expect_err(text);
        assert!(said.contains(why), "{text:?}: {said}");
    }
    // writing is refused for rows that are no alignment
    assert!(write(&[]).is_err());
    assert!(write(&rows(&[("a", "ACGT"), ("b", "ACG")]))
        .unwrap_err()
        .contains("not of one length"));
    assert!(write(&rows(&[("a", "")])).is_err());
}

#[test]
fn every_group_that_clustalw_documents_gives_its_mark_and_a_letter_outside_it_takes_it_away() {
    // the lists as the documentation gives them (written here again, apart from the module's)
    let strong = [
        "STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW",
    ];
    let weak = [
        "CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY",
    ];
    let in_a_strong_group = |letters: &str| {
        strong
            .iter()
            .any(|g| letters.chars().all(|c| g.contains(c)))
    };
    // a column with one row for each letter; the rest of each row is residues that make it a protein
    let column = |letters: &str| -> char {
        let rows: Vec<(String, String)> = letters
            .chars()
            .enumerate()
            .map(|(k, c)| (format!("r{k}"), format!("{c}{FILLER}")))
            .collect();
        let refs: Vec<(&str, &str)> = rows.iter().map(|(n, s)| (n.as_str(), s.as_str())).collect();
        marks(&refs).chars().next().expect("a mark")
    };
    for group in strong {
        assert_eq!(column(group), ':', "{group}");
        // one residue more, from outside it, and no group has them all (W is in none of the others with these)
        assert_ne!(column(&format!("{group}P")), ':', "{group} and P");
    }
    for group in weak {
        let want = if in_a_strong_group(group) { ':' } else { '.' };
        assert_eq!(column(group), want, "{group}");
    }
}

#[test]
fn a_line_of_marks_that_lost_its_indent_is_still_marks_and_a_byte_order_mark_is_left_out() {
    let text = "\u{feff}CLUSTAL W (1.83) multiple sequence alignment\n\nx  ACGT\ny  ACGA\n*** \n\nx  GG\ny  GG\n**\n";
    assert_eq!(read(text), Ok(rows(&[("x", "ACGTGG"), ("y", "ACGAGG")])));
    // a bare header over rows
    let bare = "CLUSTAL\n\nx  ACGT\ny  ACGA\n";
    assert_eq!(read(bare), Ok(rows(&[("x", "ACGT"), ("y", "ACGA")])));
}

#[test]
fn the_counts_at_the_ends_of_lines_catch_a_block_that_is_missing_and_a_count_of_columns_is_taken_too(
) {
    // the second block of this file is lost: the count says 30 residues, there are 20
    let lost =
        "CLUSTAL W\n\nx  AAAAAAAAAA 10\ny  AAAAAAAAAA 10\n\nx  CCCCCCCCCC 30\ny  CCCCCCCCCC 30\n";
    let said = read(lost).expect_err("a block is missing");
    assert!(
        said.contains("Line 6") && said.contains("30") && said.contains("missing"),
        "{said}"
    );
    assert!(read("CLUSTAL W\n\nx ACGT 99\ny ACGA 4\n").is_err());
    // counts of residues, and of columns (a gap counts in these), are both right
    assert!(read("CLUSTAL W\n\nx AC-GT 4\ny ACAGT 5\n").is_ok());
    assert!(read("CLUSTAL W\n\nx AC-GT 5\ny ACAGT 5\n").is_ok());
    // what is written is what is read, with the counts it wrote
    let a: String = "ACGTACGTAC".repeat(13);
    let b = a.replace("GT", "G-");
    assert!(read(&write(&rows(&[("a", &a), ("b", &b)])).expect("written")).is_ok());
}
