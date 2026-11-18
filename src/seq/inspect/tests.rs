use super::*;
use crate::seq::dna::STANDARD_CODE;

const NM: f64 = 250.0;
const MM: f64 = 50.0;

fn look(doc: &Doc, f: Option<usize>, r: Option<(usize, usize)>) -> Option<Inspection> {
    inspect(doc, f, r, NM, MM, STANDARD_CODE)
}

/// `CC` + ATGAAATAA (a three-residue gene, M K stop) + `GG`.
fn gene_doc() -> Doc {
    let mut d = Doc::new("g", b"CCATGAAATAAGG".to_vec());
    d.features.push(Feature::new("lacZ", "CDS", 2, 11, false));
    d
}

#[test]
fn nothing_selected_says_nothing() {
    assert_eq!(look(&gene_doc(), None, None), None);
}

#[test]
fn a_feature_is_described_by_its_own_bases() {
    let d = gene_doc();
    let i = look(&d, Some(0), None).unwrap();
    assert_eq!(i.subject, Subject::Feature(0));
    assert_eq!(i.title, "lacZ");
    assert_eq!(i.kind, "CDS");
    assert_eq!(i.location, "3..11");
    assert_eq!(i.length, 9);
    assert_eq!(i.strand, "forward");
    // ATGAAATAA: 2 of 9 are G or C
    assert!((i.gc - 100.0 * 1.0 / 9.0).abs() < 1e-9, "{}", i.gc);
    assert_eq!(i.orf, orf_summary(b"ATGAAATAA"));
    assert_eq!(
        i.protein,
        Some(Protein {
            shown: "MK*".to_string(),
            residues: 3
        })
    );
}

#[test]
fn a_reverse_feature_is_read_on_its_own_strand() {
    // ATGAAATAA written on the other strand
    let mut d = Doc::new("g", b"CCTTATTTCATGG".to_vec());
    d.features.push(Feature::new("rev", "CDS", 2, 11, true));
    let i = look(&d, Some(0), None).unwrap();
    assert_eq!(i.strand, "reverse");
    assert_eq!(i.location, "complement(3..11)");
    assert_eq!(i.protein.unwrap().shown, "MK*");
    assert_eq!(i.orf, orf_summary(b"ATGAAATAA"));
}

#[test]
fn a_feature_over_the_origin_joins_its_pieces() {
    // the gene starts at base 10 and wraps: ATG | AAATAA across the end of a 13-base circle. The pieces are in the order the
    // GenBank reader keeps them (sorted), which is not the order they are read in.
    let mut d = Doc::new("c", b"AAATAACCCCATG".to_vec());
    d.circular = true;
    let mut f = Feature::new("wrap", "CDS", 0, 6, false);
    f.spans.push((10, 13));
    d.features.push(f);
    let i = look(&d, Some(0), None).unwrap();
    assert_eq!(i.length, 9);
    assert_eq!(i.location, "join(1..6,11..13)");
    assert_eq!(i.protein.unwrap().shown, "MK*");
    assert_eq!(i.orf, orf_summary(b"ATGAAATAA"));
}

#[test]
fn a_reverse_feature_over_the_origin_is_read_from_its_tail_on_the_other_strand() {
    // ATGAAATAA on the other strand, over the origin: TTA at the end of the circle, TTTCAT at its start
    let mut d = Doc::new("c", b"TTTCATCCCCTTA".to_vec());
    d.circular = true;
    let mut f = Feature::new("wrap", "CDS", 0, 6, true);
    f.spans.push((10, 13));
    d.features.push(f);
    let i = look(&d, Some(0), None).unwrap();
    assert_eq!(i.protein.unwrap().shown, "MK*");
    assert_eq!(i.orf, orf_summary(b"ATGAAATAA"));
}

#[test]
fn a_feature_read_from_a_genbank_file_is_described_the_same_way() {
    // the reader's own route: location text in, spans sorted
    let text = "LOCUS       t  13 bp  DNA  circular\nFEATURES             Location/Qualifiers\n     CDS             join(11..13,1..6)\n                     /gene=\"x\"\nORIGIN\n        1 aaataacccc atg\n//\n";
    let d = crate::seq::io::parse_any(text, "t").unwrap().remove(0);
    let i = look(&d, Some(0), None).unwrap();
    assert_eq!(i.protein.unwrap().shown, "MK*");
}

#[test]
fn only_a_cds_has_a_protein() {
    let mut d = gene_doc();
    d.features[0].kind = "promoter".to_string();
    assert_eq!(look(&d, Some(0), None).unwrap().protein, None);
}

#[test]
fn a_long_protein_is_cut_where_the_panel_can_hold_it() {
    let mut seq = b"ATG".to_vec();
    seq.extend(b"GCT".repeat(200));
    let mut d = Doc::new("big", seq);
    d.features.push(Feature::new("big", "CDS", 0, 603, false));
    let p = look(&d, Some(0), None).unwrap().protein.unwrap();
    assert_eq!(p.residues, 201);
    assert_eq!(p.shown.len(), MOST_PROTEIN);
    assert!(p.shown.starts_with("MAAA"));
}

#[test]
fn the_genetic_code_in_use_changes_the_protein() {
    use crate::seq::dna::GENETIC_CODES;
    let mut d = Doc::new("m", b"TGA".to_vec());
    d.features.push(Feature::new("x", "CDS", 0, 3, false));
    let std = inspect(&d, Some(0), None, NM, MM, GENETIC_CODES[0].table).unwrap();
    let mito = inspect(&d, Some(0), None, NM, MM, GENETIC_CODES[1].table).unwrap();
    assert_eq!(std.protein.unwrap().shown, "*");
    assert_eq!(mito.protein.unwrap().shown, "W");
}

#[test]
fn a_range_without_a_feature_is_described_as_a_range() {
    let d = gene_doc();
    let i = look(&d, None, Some((0, 4))).unwrap();
    assert_eq!(i.subject, Subject::Range(0, 4));
    assert_eq!(i.title, "Selection");
    assert_eq!(i.kind, "Range");
    assert_eq!(i.location, "1..4");
    assert_eq!(i.length, 4);
    assert_eq!(i.strand, "");
    assert_eq!(i.protein, None);
    assert!((i.gc - 50.0).abs() < 1e-9); // CCAT: two of four are G or C
    assert_eq!(i.tm, melting_temp(b"CCAT", NM, MM));
}

#[test]
fn the_selected_feature_wins_over_the_range() {
    let d = gene_doc();
    let i = look(&d, Some(0), Some((0, 4))).unwrap();
    assert_eq!(i.subject, Subject::Feature(0));
}

#[test]
fn a_stale_feature_index_falls_back_to_the_range() {
    let d = gene_doc();
    assert_eq!(look(&d, Some(7), None), None);
    assert_eq!(
        look(&d, Some(7), Some((0, 3))).unwrap().subject,
        Subject::Range(0, 3)
    );
}

#[test]
fn an_empty_or_backwards_or_too_long_range() {
    let d = gene_doc();
    assert_eq!(look(&d, None, Some((5, 5))), None);
    assert_eq!(look(&d, None, Some((9, 3))), None);
    assert_eq!(look(&d, None, Some((99, 120))), None);
    // a range that runs past the end is cut at the end
    assert_eq!(look(&d, None, Some((10, 99))).unwrap().length, 3);
}

#[test]
fn the_melting_temperature_has_no_answer_past_a_hundred_bases() {
    let d = Doc::new("long", b"ACGT".repeat(40));
    assert_eq!(look(&d, None, Some((0, 160))).unwrap().tm, None);
    assert!(look(&d, None, Some((0, 20))).unwrap().tm.is_some());
}

#[test]
fn qualifiers_come_along() {
    let mut d = gene_doc();
    d.features[0].qualifiers = vec![("gene".into(), "lacZ".into()), ("note".into(), "x".into())];
    let i = look(&d, Some(0), None).unwrap();
    assert_eq!(i.qualifiers.len(), 2);
    assert_eq!(i.qualifiers[1], ("note".to_string(), "x".to_string()));
}

#[test]
fn overlapping_features_are_listed_by_where_they_start() {
    let mut d = Doc::new("o", vec![b'A'; 100]);
    d.features.push(Feature::new("main", "CDS", 20, 60, false)); // 0
    d.features
        .push(Feature::new("late", "promoter", 50, 70, false)); // 1
    d.features
        .push(Feature::new("early", "terminator", 10, 25, true)); // 2
    d.features
        .push(Feature::new("apart", "rep_origin", 80, 90, false)); // 3
    d.features
        .push(Feature::new("touch", "misc", 60, 65, false)); // 4: starts where main ends, shares nothing
    let i = look(&d, Some(0), None).unwrap();
    assert_eq!(i.overlaps, vec![2, 1]);
}

#[test]
fn a_range_lists_the_features_it_touches() {
    let mut d = Doc::new("o", vec![b'A'; 100]);
    d.features.push(Feature::new("a", "CDS", 20, 30, false));
    d.features.push(Feature::new("b", "CDS", 30, 40, false));
    let i = look(&d, None, Some((29, 31))).unwrap();
    assert_eq!(i.overlaps, vec![0, 1]);
    let i = look(&d, None, Some((30, 31))).unwrap();
    assert_eq!(i.overlaps, vec![1]);
}

#[test]
fn a_feature_never_overlaps_itself_even_with_two_pieces() {
    let mut d = Doc::new("c", vec![b'A'; 20]);
    let mut f = Feature::new("wrap", "CDS", 15, 20, false);
    f.spans.push((0, 4));
    d.features.push(f);
    assert!(look(&d, Some(0), None).unwrap().overlaps.is_empty());
}

fn labels(i: &Inspection) -> Vec<&'static str> {
    i.rows().into_iter().map(|(k, _)| k).collect()
}

#[test]
fn a_feature_gives_every_line() {
    let i = look(&gene_doc(), Some(0), None).unwrap();
    let rows = i.rows();
    assert_eq!(
        labels(&i),
        ["Type", "Location", "Length", "Strand", "GC", "Tm", "ORF"]
    );
    assert_eq!(rows[0].1, "CDS");
    assert_eq!(rows[1].1, "3..11");
    assert_eq!(rows[2].1, "9 bp · 3 codons");
    assert_eq!(rows[3].1, "forward");
    assert_eq!(rows[4].1, "11%");
}

#[test]
fn a_range_has_no_strand_line_and_a_long_one_no_tm() {
    let d = Doc::new("long", b"ACGT".repeat(40));
    let i = look(&d, None, Some((0, 160))).unwrap();
    let l = labels(&i);
    assert!(!l.contains(&"Strand"));
    assert!(!l.contains(&"Tm"));
    assert_eq!(i.rows()[2].1, "160 bp");
}
