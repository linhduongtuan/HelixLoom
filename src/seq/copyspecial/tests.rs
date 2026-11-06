use super::golden::CASES;
use super::*;
use crate::seq::dna::STANDARD_CODE;
use crate::seq::features::Feature;

#[test]
fn every_conversion_equals_what_the_originals_own_code_copies() {
    for (n, c) in CASES.iter().enumerate() {
        let kind = match c.function {
            "upper" => Kind::Uppercase,
            "translated" => Kind::Translated,
            "translated_uc" => Kind::TranslatedUppercase,
            "reverse" => Kind::Reverse,
            "complement" => Kind::Complement,
            "rna" => Kind::Rna,
            other => panic!("unknown function {other}"),
        };
        let got = convert(kind, c.text.as_bytes(), c.rev_com, STANDARD_CODE);
        assert_eq!(
            got,
            c.copied,
            "case {n}: {} of {:?}{}",
            c.function,
            c.text,
            if c.rev_com {
                " reverse-complemented"
            } else {
                ""
            }
        );
    }
}

#[test]
fn fasta_has_the_originals_header_and_one_line_of_bases() {
    assert_eq!(fasta("pDemo", 10, 20, b"ACGT"), ">pDemo 10-20\nACGT");
}

fn doc_with(features: Vec<Feature>) -> Doc {
    let mut doc = Doc::new("t", vec![b'A'; 100]);
    doc.features = features;
    doc
}

#[test]
fn the_bankit_table_lists_locations_and_a_label_per_feature() {
    let forward = Feature::new("T7 promoter", "promoter", 150, 170, false);
    let reverse = Feature::new("rev gene", "CDS", 200, 230, true);
    let mut split = Feature::new("exons", "CDS", 10, 20, false);
    split.spans.push((40, 55));
    let mut split_reverse = split.clone();
    split_reverse.reverse = true;
    split_reverse.name = "exons rev".into();
    let doc = doc_with(vec![forward, reverse, split, split_reverse]);
    assert_eq!(
        bankit(&doc),
        "151\t170\tpromoter\n\t\t\tlabel\tT7 promoter\n\
         230\t201\tCDS\n\t\t\tlabel\trev gene\n\
         11\t20\tCDS\n41\t55\n\t\t\tlabel\texons\n\
         55\t41\tCDS\n20\t11\n\t\t\tlabel\texons rev\n"
    );
    assert_eq!(bankit(&doc_with(Vec::new())), "");
}

#[test]
fn a_feature_across_the_origin_starts_with_the_part_at_the_end() {
    let mut f = Feature::new("wrap", "CDS", 0, 10, false);
    f.spans.push((90, 100));
    assert_eq!(
        bankit(&doc_with(vec![f])),
        "91\t100\tCDS\n1\t10\n\t\t\tlabel\twrap\n"
    );
}
