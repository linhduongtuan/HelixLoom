use super::golden::CASES;
use super::*;
use crate::seq::dna::STANDARD_CODE;

fn options(c: &golden::Case) -> Options {
    Options {
        letters: if c.one_letter {
            Letters::One
        } else {
            Letters::Three
        },
        codon_spacing: c.spacing,
        reverse_complement: c.revcom,
        width: c.width,
        numbers: match c.numbers {
            "left" => Numbers::Left,
            "right" => Numbers::Right,
            "both" => Numbers::Both,
            _ => Numbers::None,
        },
        dna_line: match c.dna_line {
            "above" => DnaLine::Above,
            "below" => DnaLine::Below,
            _ => DnaLine::None,
        },
    }
}

use super::golden;

/// The header and the lines as one text, every line ending in a newline (what the original's window holds).
fn text(out: &Output) -> String {
    out.header
        .iter()
        .chain(&out.lines)
        .map(|l| format!("{l}\n"))
        .collect()
}

#[test]
fn the_text_equals_what_the_originals_own_code_builds() {
    for (n, c) in CASES.iter().enumerate() {
        let source = if c.cds.is_empty() {
            Source::Range { first: c.first }
        } else {
            Source::Feature { name: c.cds }
        };
        let got = text(&build(c.dna.as_bytes(), STANDARD_CODE, source, &options(c)));
        assert_eq!(
            got,
            c.text,
            "case {n}: {} bases from {}",
            c.dna.len(),
            c.first
        );
    }
}

#[test]
fn the_defaults_are_the_originals() {
    let d = Options::default();
    assert_eq!(
        (d.letters, d.codon_spacing, d.reverse_complement, d.width),
        (Letters::Three, false, false, 20)
    );
    assert_eq!((d.numbers, d.dna_line), (Numbers::None, DnaLine::None));
}

#[test]
fn molecular_weight_adds_residue_masses_and_water() {
    // G + A + S + K + M = 474.58, plus water
    assert_eq!(molecular_weight(b"MASKG"), 492.6);
    assert_eq!(molecular_weight(b""), 18.02);
    // X and the stop count nothing; lower case counts like upper case
    assert_eq!(molecular_weight(b"MX*"), molecular_weight(b"m"));
    assert_eq!(tcl_number(492.6), "492.6");
    assert_eq!(tcl_number(100.0), "100.0");
}

#[test]
fn a_feature_ignores_reverse_complement_and_keeps_a_trailing_partial_codon_in_its_length() {
    let o = Options {
        reverse_complement: true,
        ..Options::default()
    };
    let out = build(b"ATGGCTA", STANDARD_CODE, Source::Feature { name: "x" }, &o);
    assert_eq!(out.header[0], "x ");
    assert_eq!(out.lines, ["MetAla"]);
}
