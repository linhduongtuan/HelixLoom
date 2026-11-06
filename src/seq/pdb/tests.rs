use super::*;

/// Ubiquitin, PDB entry 1UBQ (the example structure shipped with the program).
const UBIQUITIN: &str = include_str!("../../../assets/structure/1ubq.pdb");

const UBIQUITIN_SEQUENCE: &str =
    "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

/// One ATOM line in the fixed columns.
#[allow(clippy::too_many_arguments)] // one argument for each field of the line, as the format has them
fn atom(
    record: &str,
    serial: u32,
    name: &str,
    alt: char,
    residue: &str,
    chain: char,
    number: i32,
    xyz: [f32; 3],
    b: f32,
) -> String {
    format!(
        "{record:<6}{serial:>5} {name:<4}{alt}{residue:>3} {chain}{number:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
        xyz[0],
        xyz[1],
        xyz[2],
        1.0,
        b,
        name.chars().next().unwrap_or('C')
    )
}

#[test]
fn ubiquitin_is_read_with_its_residues_chain_sequence_and_annotated_elements() {
    let s = parse(UBIQUITIN).expect("a structure");
    assert_eq!(
        s.title,
        "STRUCTURE OF UBIQUITIN REFINED AT 1.8 ANGSTROMS RESOLUTION"
    );
    assert!(s.experimental && !s.looks_predicted());
    assert_eq!(s.atoms.len(), 602, "602 atoms; the 58 waters are left out");
    assert_eq!(s.chains(), vec!['A']);
    assert_eq!(
        String::from_utf8_lossy(&s.sequence('A')),
        UBIQUITIN_SEQUENCE
    );
    assert_eq!(s.amino_acids().len(), 76);
    assert_eq!(
        s.helices,
        vec![
            Span {
                chain: 'A',
                from: 23,
                to: 34
            },
            Span {
                chain: 'A',
                from: 56,
                to: 59
            }
        ]
    );
    let strands: Vec<(i32, i32)> = s.strands.iter().map(|x| (x.from, x.to)).collect();
    assert_eq!(
        strands,
        vec![(10, 17), (1, 7), (64, 72), (40, 45), (48, 50)]
    );
    // the first atom of the file
    let first = &s.atoms[0];
    assert_eq!((first.name.as_str(), first.element.as_str()), ("N", "N"));
    assert_eq!(s.residues[0].name, "MET");
}

#[test]
fn the_backbone_has_four_atoms_per_residue_and_every_residue_follows_the_one_before() {
    let s = parse(UBIQUITIN).expect("a structure");
    let backbone = s.backbone();
    assert_eq!(backbone.len(), 76);
    assert!(!backbone[0].2, "the first residue follows nothing");
    assert!(backbone[1..].iter().all(|b| b.2), "one unbroken chain");
    // a peptide bond: C(i) to N(i+1) is about 1.33 angstroms; CA to CA about 3.8
    let d = |a: [f32; 3], b: [f32; 3]| distance(a, b);
    for w in backbone.windows(2) {
        assert!((d(w[0].1.c, w[1].1.n) - 1.33).abs() < 0.2);
        assert!((d(w[0].1.ca, w[1].1.ca) - 3.8).abs() < 0.3);
    }
    assert_eq!(backbone.iter().filter(|b| b.1.proline).count(), 3);
    let (centre, radius) = s.centre_and_radius();
    assert!((10.0..40.0).contains(&radius), "radius {radius}");
    assert!(s
        .atoms
        .iter()
        .all(|a| distance(a.pos, centre) <= radius + 1e-3));
}

#[test]
fn only_the_first_model_and_the_first_conformer_are_read_and_water_is_left_out() {
    let text = [
        atom("ATOM", 1, "N", ' ', "ALA", 'A', 1, [0.0, 0.0, 0.0], 10.0),
        atom("ATOM", 2, "CA", 'A', "ALA", 'A', 1, [1.5, 0.0, 0.0], 10.0),
        atom("ATOM", 3, "CA", 'B', "ALA", 'A', 1, [9.9, 9.9, 9.9], 10.0),
        atom(
            "HETATM",
            4,
            "O",
            ' ',
            "HOH",
            'A',
            101,
            [5.0, 5.0, 5.0],
            30.0,
        ),
        "ENDMDL".to_string(),
        atom("ATOM", 5, "N", ' ', "GLY", 'A', 2, [3.0, 0.0, 0.0], 10.0),
    ]
    .join("\n");
    let s = parse(&text).expect("a structure");
    assert_eq!(s.atoms.len(), 2);
    assert_eq!(s.atoms[1].pos, [1.5, 0.0, 0.0]);
    assert_eq!(s.residues.len(), 1);
}

#[test]
fn a_ligand_is_not_part_of_the_chain_but_selenomethionine_is() {
    let text = [
        atom("ATOM", 1, "CA", ' ', "ALA", 'A', 1, [0.0, 0.0, 0.0], 10.0),
        atom("HETATM", 2, "CA", ' ', "MSE", 'A', 2, [3.8, 0.0, 0.0], 10.0),
        atom(
            "HETATM",
            3,
            "FE",
            ' ',
            "HEM",
            'A',
            300,
            [9.0, 0.0, 0.0],
            10.0,
        ),
    ]
    .join("\n");
    let s = parse(&text).expect("a structure");
    assert_eq!(s.residues.len(), 3);
    assert_eq!(String::from_utf8_lossy(&s.sequence('A')), "AM");
    assert!(s.residues[2].hetero && !s.residues[1].hetero);
}

#[test]
fn text_that_is_not_a_structure_is_refused_and_short_lines_do_not_fail() {
    assert!(parse("").unwrap_err().contains("No atoms"));
    assert!(parse(">sequence\nACGT\n").unwrap_err().contains("PDB"));
    let ragged = "ATOM      1  CA\nHELIX   1\nSHEET\nATOM      1  N   ALA A   1       1.000   2.000   x.xxx\n";
    assert!(parse(ragged).is_err());
}

#[test]
fn a_predictors_confidence_is_read_from_the_temperature_factor_in_either_scale() {
    let make = |b: f32| {
        let text = [
            atom("ATOM", 1, "CA", ' ', "ALA", 'A', 1, [0.0, 0.0, 0.0], b),
            atom(
                "ATOM",
                2,
                "CA",
                ' ',
                "GLY",
                'A',
                2,
                [3.8, 0.0, 0.0],
                b / 2.0,
            ),
        ]
        .join("\n");
        parse(&text).expect("a structure")
    };
    // 0 to 100, as AlphaFold and the ESM Atlas write it
    let percent = make(90.0);
    assert_eq!(percent.confidence(0), Some(90.0));
    assert!(percent.looks_predicted());
    // 0 to 1 is scaled up
    let fraction = make(0.9);
    assert!((fraction.confidence(0).unwrap() - 90.0).abs() < 1e-3);
    assert_eq!(fraction.confidence(9), None);
    // an experiment's temperature factors can be anything
    assert!(!make(250.0).looks_predicted());
}

#[test]
fn letters_cover_the_twenty_amino_acids_and_modified_ones() {
    let twenty = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    ];
    let letters: String = twenty.iter().map(|n| char::from(letter(n))).collect();
    assert_eq!(letters, "ARNDCQEGHILKMFPSTWYV");
    assert_eq!(letter("MSE"), b'M');
    assert_eq!(letter("HEM"), b'X');
}
