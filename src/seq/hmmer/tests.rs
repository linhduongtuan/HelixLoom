use super::*;

/// A table made to the format the HMMER 3.4 User's Guide describes (22 fields and a description), by hand: it is not output of the
/// program. Sequence `q1` has two fibronectin domains and a kinase (given out of order), `q2` a globin.
const TABLE: &str = "\
# target name        accession   tlen query name  accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- ----- ---------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
fn3                  PF00041.26    85 q1         -            700   1.2e-40  130.5   0.3   2   3   4.0e-14   6.5e-11   39.8   0.2     3    85   130   215   128   216 0.97 Fibronectin type III domain
Pkinase              PF00069.30   264 q1         -            700   1.2e-40  130.5   0.3   1   1   1e-86   2.2e-86  290.1   0.9     1   264   420   690   418   692 0.98 Protein kinase domain

fn3                  PF00041.26    85 q1         -            700   1.2e-40  130.5   0.3   1   3   2.1e-14   3.5e-11   41.2   0.1     4    84    30   112    28   114 0.95 Fibronectin type III domain
Globin               PF00042.25   109 q2         -            150   3.0e-20   70.0   0.0   1   1   1e-17     2e-17     60.0   0.0     5   105    10   120     8   122 0.90 Globin
Mystery              -            120 q2         -            150   1.0e-3    12.0   0.0   1   1   0.5       0.9       10.0   0.0     1    30   125   148   124   150 0.80 
# Program:         hmmscan
# [ok]
";

fn parse() -> Vec<DomainHit> {
    parse_domtblout(TABLE, &[("q1", 700), ("q2", 150)]).expect("a table")
}

#[test]
fn the_domains_of_a_table_are_read_with_their_places_scores_and_descriptions() {
    let hits = parse();
    assert_eq!(
        hits.len(),
        5,
        "comments, the rule line and blank lines are not domains"
    );
    let kinase = &hits[1];
    assert_eq!(
        (
            kinase.sequence,
            kinase.family.as_str(),
            kinase.accession.as_str()
        ),
        (0, "Pkinase", "PF00069.30")
    );
    assert_eq!(
        (
            kinase.ali_from,
            kinase.ali_to,
            kinase.env_from,
            kinase.env_to
        ),
        (420, 690, 418, 692)
    );
    assert_eq!(
        (kinase.hmm_from, kinase.hmm_to, kinase.family_length),
        (1, 264, 264)
    );
    assert_eq!(
        (kinase.score, kinase.i_evalue, kinase.c_evalue),
        (290.1, 2.2e-86, 1e-86)
    );
    assert_eq!(kinase.accuracy, 0.98);
    assert_eq!(kinase.description, "Protein kinase domain");
    // the description is the rest of the line, spaces and all; none is none; a family with no accession has an empty one
    assert_eq!(hits[0].description, "Fibronectin type III domain");
    let mystery = &hits[4];
    assert_eq!(
        (
            mystery.sequence,
            mystery.accession.as_str(),
            mystery.description.as_str()
        ),
        (1, "", "")
    );
    assert_eq!(hits[3].sequence, 1);
}

#[test]
fn the_domains_of_a_sequence_come_in_the_order_they_lie_along_it() {
    let hits = parse();
    let along_q1: Vec<(&str, usize)> = along(&hits, 0)
        .iter()
        .map(|h| (h.family.as_str(), h.ali_from))
        .collect();
    assert_eq!(along_q1, [("fn3", 30), ("fn3", 130), ("Pkinase", 420)]);
    assert_eq!(along(&hits, 1).len(), 2);
    assert!(along(&hits, 7).is_empty());
}

#[test]
fn envelopes_that_share_a_residue_overlap_only_in_the_same_sequence() {
    let hits = parse();
    // the globin envelope 8-122 and the other family's 124-150 do not; the two fibronectin domains do not
    assert!(!hits[3].overlaps(&hits[4]));
    assert!(!hits[2].overlaps(&hits[0]));
    let mut moved = hits[4].clone();
    moved.env_from = 122;
    assert!(hits[3].overlaps(&moved) && moved.overlaps(&hits[3]));
    moved.sequence = 0;
    assert!(!hits[3].overlaps(&moved));
}

#[test]
fn a_table_that_does_not_fit_the_sequences_that_were_searched_is_refused_with_the_line() {
    let one = "fn3 PF00041.26 85 q1 - 700 1e-40 130.5 0.3 1 1 2e-14 3e-11 41.2 0.1 4 84 30 112 28 114 0.95 Fibronectin";
    let ok = |text: &str| parse_domtblout(text, &[("q1", 700)]);
    assert!(ok(one).is_ok());
    let refused = |text: &str, why: &str| {
        let said = ok(text).expect_err(why);
        assert!(said.contains(why), "{why}: {said}");
    };
    // a name that was not searched, and a sequence of another length
    refused(&one.replace(" q1 ", " q9 "), "not searched");
    refused(&one.replace(" 700 ", " 701 "), "not 700");
    // too few fields, and fields that are no numbers
    refused("fn3 PF00041.26 85 q1 - 700", "6 fields");
    refused(&one.replace(" 85 ", " eighty "), "length of the profile");
    refused(&one.replace(" 41.2 ", " high "), "domain score");
    // coordinates that are not inside the sequence or the profile, or are backwards
    refused(&one.replace(" 30 112 ", " 30 701 "), "sequence of 700");
    refused(&one.replace(" 28 114 ", " 0 114 "), "sequence of 700");
    refused(&one.replace(" 30 112 ", " 112 30 "), "sequence of 700");
    refused(&one.replace(" 4 84 ", " 4 86 "), "profile of 85");
    refused(&one.replace(" 4 84 ", " 0 84 "), "profile of 85");
    // a number that cannot be a score, an E-value or a probability
    refused(&one.replace(" 3e-11 ", " -1 "), "out of range");
    refused(&one.replace(" 3e-11 ", " NaN "), "out of range");
    refused(&one.replace(" 0.95 ", " 1.5 "), "out of range");
    // the line numbers count the comments
    let said = ok(&format!("# a comment\n\n{}", one.replace(" q1 ", " q9 "))).unwrap_err();
    assert!(said.starts_with("Line 3 "), "{said}");
}

#[test]
fn the_domains_make_a_table_for_a_spreadsheet() {
    let hits = parse();
    let names = ["first\tprotein".to_string(), "globin".to_string()];
    let text = to_table(&hits[3..4], &names);
    let mut lines = text.lines();
    assert!(lines
        .next()
        .is_some_and(|h| h.starts_with("sequence\tfamily\taccession\tfrom\tto\tscore\ti_evalue")));
    assert_eq!(
        lines.next(),
        Some("globin\tGlobin\tPF00042.25\t10\t120\t60\t2e-17\t1e-17\t5\t105\t109\tGlobin")
    );
    assert_eq!(lines.next(), None);
    let two = to_table(&hits[..1], &names);
    assert!(two.contains("first protein\tfn3\t"), "{two}");
}

#[test]
fn protein_sequences_are_read_from_fasta_or_one_to_a_line_in_capitals() {
    let found = read_proteins(">sp|P1|NAME_HUMAN a protein\nmkv lle*\nGA\n>second\nAC12D\n")
        .expect("fasta");
    assert_eq!(
        found,
        [
            ("sp|P1|NAME_HUMAN".to_string(), b"MKVLLEGA".to_vec()),
            ("second".to_string(), b"ACD".to_vec())
        ]
    );
    let lines = read_proteins("mkvlle\n\n 123 \nGA\n").expect("lines");
    assert_eq!(
        lines,
        [
            ("sequence 1".to_string(), b"MKVLLE".to_vec()),
            ("sequence 2".to_string(), b"GA".to_vec())
        ]
    );
    for nothing in ["", "   \n", ">only a header\n"] {
        assert!(read_proteins(nothing).unwrap_err().contains("at least one"));
    }
    let long = format!(">long\n{}\n", "A".repeat(MAX_LENGTH + 1));
    assert!(read_proteins(&long).unwrap_err().contains("longest"));
    let many = "MKV\n".repeat(MAX_SEQUENCES + 1);
    assert!(read_proteins(&many).unwrap_err().contains("more than"));
}
