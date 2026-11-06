use super::*;
use crate::seq::enzymes::{default_enzymes, find_sites};
use crate::seq::sitesearch::golden;

/// 30 cases made by running the original's `do_diagnostic_site` on random sequences: the sequence, the stretch selected,
/// the choice of enzymes and the number of mismatches, the sites the original had already found, and the text of its window.
const GOLDEN: &str = include_str!("golden.txt");

fn enzyme(name: &str) -> Enzyme {
    default_enzymes()
        .into_iter()
        .find(|e| e.name == name)
        .expect("a default enzyme")
}

#[test]
fn reports_equal_the_originals_on_every_case() {
    let all = default_enzymes();
    let cases = golden::parse(GOLDEN);
    assert_eq!(cases.len(), 30);
    for (n, case) in cases.iter().enumerate() {
        let chosen: Vec<&Enzyme> = if case.mode == "absent" {
            all.iter()
                .filter(|e| find_sites(e, &case.seq, case.circular).is_empty())
                .collect()
        } else {
            case.chosen
                .iter()
                .filter_map(|name| all.iter().find(|e| &e.name == name))
                .collect()
        };
        let (from, to) = case.region();
        let dna = &case.seq[from..to];
        let found = find(dna, &chosen, case.mismatches);
        let text = report(("DATE", "TITLE"), (from, to - 1), from, dna, &found.found);
        assert_eq!(text, case.text, "case {n}");
    }
}

#[test]
fn site_lists_equal_the_originals() {
    golden::check_sites(&golden::parse(GOLDEN));
}

#[test]
fn a_site_one_or_two_bases_away_is_found_on_the_strand_it_is_nearer() {
    let bsi = enzyme("BsiWI");
    let two = find(b"TTCGAATGTT", &[&bsi], 2);
    assert_eq!(two.found.len(), 1);
    let f = &two.found[0];
    assert_eq!((f.index, f.mismatches), (2, 2));
    assert_eq!(f.pattern, b"CGTACG".to_vec());
    assert!(find(b"TTCGAATGTT", &[&bsi], 1).found.is_empty());
}

#[test]
fn an_exact_site_is_not_a_diagnostic_site() {
    let eco = enzyme("EcoRI");
    assert!(find(b"AAGAATTCAA", &[&eco], 3).found.is_empty());
}

#[test]
fn a_palindrome_is_tried_once_and_other_sites_on_both_strands() {
    let eco = enzyme("EcoRI");
    let bsa = enzyme("BsaI");
    let (a, b) = (find(b"GAATTG", &[&eco], 1), find(b"GAGACT", &[&bsa], 1));
    assert_eq!(a.found.len(), 1);
    // GGTCTC read the other way round is GAGACC: the DNA GAGACT is one base off it
    assert_eq!(b.found.len(), 1);
    assert_eq!(b.found[0].pattern, b"GAGACC".to_vec());
}

#[test]
fn enzymes_that_need_methylated_dna_are_left_out() {
    let dpn = enzyme("DpnI");
    assert!(SKIPPED.contains(&"DpnI"));
    assert!(find(b"GATA", &[&dpn], 1).found.is_empty());
}

#[test]
fn results_come_in_order_of_position_and_equal_positions_in_order_of_enzyme() {
    let (eco, bam) = (enzyme("EcoRI"), enzyme("BamHI"));
    let found = find(b"GAATTGGGATCAAA", &[&eco, &bam], 1).found;
    let order: Vec<(usize, &str)> = found.iter().map(|f| (f.index, f.enzyme.as_str())).collect();
    assert_eq!(order, vec![(0, "EcoRI"), (6, "BamHI")]);
}

#[test]
fn ambiguity_codes_in_the_site_accept_the_bases_they_cover() {
    let sgr = enzyme("SgrAI");
    let dna = b"GACAGGCG";
    let found = find(dna, &[&sgr], 2).found;
    assert_eq!(found.len(), 1);
    let shown: String = marks(&found[0], dna)
        .iter()
        .map(|&(c, changed)| {
            if changed {
                format!("<{c}>")
            } else {
                c.to_string()
            }
        })
        .collect();
    assert_eq!(shown, "<C>AC<C>GGCG");
}

#[test]
fn a_short_dna_has_no_window() {
    let sgr = enzyme("SgrAI");
    assert!(find(b"GAC", &[&sgr], 3).found.is_empty());
    assert!(find(b"", &[&sgr], 3).found.is_empty());
}

#[test]
fn a_search_that_finds_too_many_places_stops_at_the_limit_and_keeps_them_in_order() {
    let all = default_enzymes();
    let enzymes: Vec<&Enzyme> = all.iter().collect();
    let dna: Vec<u8> = (0..30_000).map(|i| b"ACGT"[(i * 7 + i / 3) % 4]).collect();
    let search = find(&dna, &enzymes, 3);
    assert!(search.truncated);
    assert_eq!(search.found.len(), MAX_FOUND);
    assert!(search.found.windows(2).all(|w| w[0].index <= w[1].index));
}
