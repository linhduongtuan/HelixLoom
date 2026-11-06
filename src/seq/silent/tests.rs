use super::*;
use crate::seq::dna::STANDARD_CODE;
use crate::seq::enzymes::{default_enzymes, find_sites};
use crate::seq::sitesearch::golden;

/// 80 cases made by running the original's `silent_sites` on random sequences: the sequence, the stretch selected, the
/// enzymes chosen, the sites the original had already found, and the text of its window.
const GOLDEN: &str = include_str!("golden.txt");

fn enzyme(name: &str) -> Enzyme {
    default_enzymes()
        .into_iter()
        .find(|e| e.name == name)
        .expect("a default enzyme")
}

fn run(seq: &[u8], region: (usize, usize), enzymes: &[&Enzyme], existing: &[Vec<usize>]) -> Search {
    search(&Query {
        seq,
        region,
        table: STANDARD_CODE,
        enzymes,
        existing,
    })
}

#[test]
fn reports_equal_the_originals_on_every_case() {
    let all = default_enzymes();
    let cases = golden::parse(GOLDEN);
    assert_eq!(cases.len(), 80);
    for (n, case) in cases.iter().enumerate() {
        let chosen: Vec<&Enzyme> = if case.chosen.is_empty() {
            all.iter().collect()
        } else {
            case.chosen
                .iter()
                .filter_map(|name| all.iter().find(|e| &e.name == name))
                .collect()
        };
        let existing: Vec<Vec<usize>> = chosen
            .iter()
            .map(|e| {
                let sites = find_sites(e, &case.seq, case.circular);
                let mut starts: Vec<usize> = sites.iter().map(|s| s.start + 1).collect();
                starts.sort_unstable();
                starts
            })
            .collect();
        let region = case.region();
        let found = run(&case.seq, region, &chosen, &existing);
        let text = report("DATE", "TITLE", region.0 + 1, region.1, &found.hits);
        assert_eq!(text, case.text, "case {n}");
    }
}

#[test]
fn site_lists_equal_the_originals() {
    golden::check_sites(&golden::parse(GOLDEN));
}

#[test]
fn a_site_is_read_in_three_frames_with_a_set_of_residues_per_codon() {
    let frames = degenerate_translation(b"GAATTC");
    // frame 0 starts the site at a codon: GAA TTC is Glu Phe
    let glu = Residues::Of(residue_bit(b'E'));
    let phe = Residues::Of(residue_bit(b'F'));
    assert_eq!(frames[2], vec![glu, phe]);
    // the other two frames have a codon with a base before the site or after it
    assert_eq!(frames[0].len(), 3);
    assert_eq!(frames[1].len(), 3);
    assert_eq!(Residues::of_triplet(b"NNN"), Residues::Any);
    assert!(Residues::of_triplet(b"TAR").allows(b'*'));
    assert!(!Residues::of_triplet(b"TAR").allows(b'X'));
    assert!(Residues::Any.allows(b'X'));
}

#[test]
fn matches_do_not_overlap_like_regexp_all() {
    let any = [Residues::Any, Residues::Any];
    assert_eq!(scan(b"ABCDE", &any), vec![0, 2]);
    assert_eq!(scan(b"A", &any), Vec::<usize>::new());
    assert_eq!(scan(b"ABC", &[]), Vec::<usize>::new());
}

#[test]
fn the_fewest_changes_that_keep_the_residues_are_made() {
    // Ser Stop -> AluI AGCT: only the first codon changes, to another Ser codon
    assert_eq!(
        coordinate(b"TCATAA", b"AGCTNN", STANDARD_CODE),
        vec![b"AGCTAA".to_vec()]
    );
    // an ambiguity code the DNA already fits stays as the DNA has it
    assert_eq!(
        coordinate(b"GGGGCC", b"GGNGCC", STANDARD_CODE),
        vec![b"GGGGCC".to_vec()]
    );
    // changes that cannot keep the residue give the original's marker
    assert_eq!(fit_codon(b"ATG", b"TGA", STANDARD_CODE), b"xxx".to_vec());
}

#[test]
fn a_star_or_question_mark_in_the_dna_matches_anything_as_in_the_originals_glob() {
    // `string match "*$orig_char*" $possible`: with `*` or `?` as the DNA's base the answer is yes whatever the site says
    assert_eq!(substitute(b"*CG", b"RYS"), b"*CG".to_vec());
    assert_eq!(substitute(b"?TA", b"NNN"), b"?TA".to_vec());
    // an ordinary base the site's letter does not allow is replaced by the letter
    assert_eq!(substitute(b"ACG", b"RYS"), b"ACG".to_vec());
    assert_eq!(substitute(b"CAT", b"RYS"), b"RYS".to_vec());
}

#[test]
fn ambiguous_results_are_spelled_out_with_the_first_letter_changing_fastest() {
    assert_eq!(
        expand_all(b"RYA"),
        vec![
            b"ACA".to_vec(),
            b"GCA".to_vec(),
            b"ATA".to_vec(),
            b"GTA".to_vec()
        ]
    );
    assert_eq!(expand_all(b"ACG"), vec![b"ACG".to_vec()]);
    assert_eq!(expand_all(b"S"), vec![b"G".to_vec(), b"C".to_vec()]);
}

#[test]
fn a_place_that_already_is_a_site_is_not_offered() {
    let eco = enzyme("EcoRI");
    let seq = b"GAATTC";
    let with_site = run(seq, (0, 6), &[&eco], &[vec![1]]);
    assert!(with_site.hits.is_empty());
    let without = run(seq, (0, 6), &[&eco], &[Vec::new()]);
    assert_eq!(without.hits.len(), 1);
    let hit = &without.hits[0];
    assert_eq!((hit.position, hit.site_start, hit.frame), (1, 1, 0));
    assert_eq!(hit.variants, vec![b"GAATTC".to_vec()]);
    // the original pads the site with three spaces before it cuts it into threes
    assert_eq!(hit.layout, "GAA TTC    ");
    assert_eq!(hit.residues, "E   F");
    assert_eq!(hit.codons_at, (0, 6));
}

#[test]
fn a_stretch_is_searched_alone_and_numbered_along_the_whole_sequence() {
    let eco = enzyme("EcoRI");
    let seq = b"AAAAAAGAATTCAAA";
    let hits = run(seq, (6, 12), &[&eco], &[Vec::new()]).hits;
    assert_eq!(hits.len(), 1);
    assert_eq!((hits[0].position, hits[0].codons_at), (7, (6, 12)));
}

#[test]
fn a_stretch_outside_the_sequence_finds_nothing() {
    let eco = enzyme("EcoRI");
    assert_eq!(run(b"GAATTC", (3, 9), &[&eco], &[]), Search::default());
    assert_eq!(run(b"GAATTC", (4, 2), &[&eco], &[]), Search::default());
}

#[test]
fn positions_wrap_like_seq_len_mod() {
    assert_eq!(wrap(0, 5), 5);
    assert_eq!(wrap(5, 5), 5);
    assert_eq!(wrap(7, 5), 2);
    assert_eq!(wrap(3, 5), 3);
}

#[test]
fn changed_bases_are_marked_in_the_report_text() {
    let m = marks(b"AGCTAA", b"tcataa");
    assert_eq!(
        m.iter().map(|&(_, changed)| changed).collect::<Vec<_>>(),
        vec![true, true, true, false, false, false]
    );
    assert_eq!(spaced(&m), "<A><G><C> TAA");
    let none = run(b"AAAAAA", (0, 6), &[&enzyme("EcoRI")], &[]);
    assert_eq!(
        report("D", "T", 1, 6, &none.hits),
        "D\nT \nFrom 1 to 6. \nSilent Sites\n\nNo sites Found"
    );
}

#[test]
fn a_search_that_finds_too_many_places_stops_at_the_limit_and_keeps_them_in_order() {
    let all = default_enzymes();
    let enzymes: Vec<&Enzyme> = all.iter().collect();
    let existing = vec![Vec::new(); enzymes.len()];
    let mut state = 7u32;
    let seq: Vec<u8> = (0..150_000)
        .map(|_| {
            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            b"ACGT"[(state >> 16) as usize % 4]
        })
        .collect();
    let found = run(&seq, (0, seq.len()), &enzymes, &existing);
    assert!(found.truncated);
    assert_eq!(found.hits.len(), MAX_PLACES);
    assert!(found
        .hits
        .windows(2)
        .all(|w| w[0].site_start <= w[1].site_start));
}
