use super::*;
use crate::seq::pdb::{parse, Structure};

const UBIQUITIN: &str = include_str!("../../../assets/structure/1ubq.pdb");

/// What PyDSSP 0.9.1 gives for ubiquitin (PDB 1UBQ): H helix, E strand, - coil.
const UBIQUITIN_PYDSSP: &str =
    "-EEEEEE----EEEEE-----EHHHHHHHHHHHH---HHHEEEEE--EE-----E-HHH------EEEEEE-----";

fn states(s: &Structure) -> (String, Vec<usize>) {
    let backbone = s.backbone();
    let residues: Vec<(Backbone, bool)> = backbone
        .iter()
        .map(|(_, b, linked)| (*b, *linked))
        .collect();
    let text = assign(&residues).iter().map(|s| s.letter()).collect();
    (text, backbone.iter().map(|(i, _, _)| *i).collect())
}

fn agreement(a: &str, b: &str) -> f64 {
    let same = a.chars().zip(b.chars()).filter(|(x, y)| x == y).count();
    100.0 * same as f64 / a.len().max(1) as f64
}

#[test]
fn ubiquitin_comes_out_as_pydssp_has_it() {
    let s = parse(UBIQUITIN).expect("a structure");
    let (mine, _) = states(&s);
    assert_eq!(mine.len(), UBIQUITIN_PYDSSP.len());
    let different: Vec<usize> = mine
        .chars()
        .zip(UBIQUITIN_PYDSSP.chars())
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, _)| i + 1)
        .collect();
    assert!(
        different.len() <= 2,
        "residues that differ: {different:?}\nmine   {mine}\npydssp {UBIQUITIN_PYDSSP}"
    );
}

#[test]
fn ubiquitin_matches_the_helices_and_strands_the_file_lists() {
    let s = parse(UBIQUITIN).expect("a structure");
    let (mine, residues) = states(&s);
    let state_of = |number: i32| -> char {
        let k = residues
            .iter()
            .position(|&r| s.residues[r].number == number)
            .expect("a residue");
        mine.chars().nth(k).expect("a state")
    };
    let share = |spans: &[(i32, i32)], want: char| -> f64 {
        let all: Vec<i32> = spans.iter().flat_map(|&(a, b)| a..=b).collect();
        100.0 * all.iter().filter(|&&n| state_of(n) == want).count() as f64 / all.len() as f64
    };
    let helices: Vec<(i32, i32)> = s.helices.iter().map(|h| (h.from, h.to)).collect();
    let strands: Vec<(i32, i32)> = s.strands.iter().map(|h| (h.from, h.to)).collect();
    // the file's 3-10 helix (56-59) is shorter than a minimal alpha helix; the long one must be found
    assert_eq!(share(&[helices[0]], 'H'), 100.0);
    // the file's strands are drawn a little wider than the bonds that make them (PyDSSP agrees: 73 %)
    assert!(share(&strands, 'E') >= 70.0, "{}", share(&strands, 'E'));
    for &(a, b) in &strands {
        assert!(share(&[(a, b)], 'E') >= 40.0, "strand {a}-{b}");
    }
}

#[test]
fn nothing_or_too_little_to_judge_is_coil() {
    assert!(assign(&[]).is_empty());
    let one = [(
        Backbone {
            n: [0.0; 3],
            ca: [1.0, 0.0, 0.0],
            c: [2.0, 0.0, 0.0],
            o: [2.0, 1.0, 0.0],
            proline: false,
        },
        false,
    )];
    assert_eq!(assign(&one), vec![Ss::Coil]);
    // a few residues in a line bond to nothing
    let line: Vec<(Backbone, bool)> = (0..8)
        .map(|i| {
            let x = 3.8 * i as f32;
            (
                Backbone {
                    n: [x, 0.0, 0.0],
                    ca: [x + 1.4, 0.0, 0.0],
                    c: [x + 2.5, 0.0, 0.0],
                    o: [x + 2.5, 1.2, 0.0],
                    proline: false,
                },
                i > 0,
            )
        })
        .collect();
    assert!(assign(&line).iter().all(|&s| s == Ss::Coil));
    assert_eq!(
        (Ss::Helix.letter(), Ss::Strand.letter(), Ss::Coil.letter()),
        ('H', 'E', '-')
    );
}

/// With `HELIXLOOM_PDB_DIR` naming a folder of `<id>.pdb` files with PyDSSP's answers as `<id>.dssp`, the answers must
/// agree on at least 90 of 100 residues of every structure. Measured: 1UBQ 100, 1LYZ 95, 3GB1 100, 1TEN 100,
/// 2LZM 100, 1CRN 91 (46 residues, so a boundary residue is 2 %).
#[test]
#[ignore = "needs real structures; see the doc comment"]
fn real_structures_agree_with_pydssp() {
    let dir = std::env::var("HELIXLOOM_PDB_DIR").expect("HELIXLOOM_PDB_DIR names the folder");
    for id in ["1UBQ", "1LYZ", "3GB1", "1TEN", "2LZM", "1CRN"] {
        let text = std::fs::read_to_string(format!("{dir}/{id}.pdb")).expect("the structure");
        let reference = std::fs::read_to_string(format!("{dir}/{id}.dssp")).expect("the reference");
        let (mine, _) = states(&parse(&text).expect("a structure"));
        let reference = reference.trim();
        assert_eq!(mine.len(), reference.len(), "{id}: residues");
        let share = agreement(&mine, reference);
        eprintln!("{id}: {share:.1} % of {} residues agree", mine.len());
        assert!(
            share >= 90.0,
            "{id}: {share:.1} %\nmine   {mine}\npydssp {reference}"
        );
    }
}
