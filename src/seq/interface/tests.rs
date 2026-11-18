use super::*;
use crate::seq::pdb::{self, Atom, Residue};

const UBIQUITIN: &str = include_str!("../../../assets/structure/1ubq.pdb");

fn ubiquitin() -> Structure {
    pdb::parse(UBIQUITIN).expect("the example structure")
}

/// Ubiquitin cut into two chains between residues 38 and 39: the first part is chain A, the rest chain B.
fn cut_ubiquitin() -> Structure {
    let mut s = ubiquitin();
    for r in &mut s.residues {
        r.chain = if r.number <= 38 { 'A' } else { 'B' };
    }
    s
}

fn atom(name: &str, element: &str, pos: [f32; 3], residue: usize) -> Atom {
    Atom {
        name: name.into(),
        element: element.into(),
        pos,
        bfactor: 0.0,
        residue,
    }
}

fn residue(name: &str, chain: char, number: i32, atoms: Vec<usize>, hetero: bool) -> Residue {
    Residue {
        name: name.into(),
        chain,
        number,
        insertion: ' ',
        atoms,
        hetero,
    }
}

/// Two residues, one atom each, `d` apart: the first in chain A, the second in chain B.
fn pair(first: (&str, &str, &str), second: (&str, &str, &str), d: f32) -> Structure {
    let mut s = Structure::default();
    s.atoms.push(atom(first.2, first.1, [0.0, 0.0, 0.0], 0));
    s.atoms.push(atom(second.2, second.1, [d, 0.0, 0.0], 1));
    s.residues.push(residue(first.0, 'A', 1, vec![0], false));
    s.residues.push(residue(second.0, 'B', 1, vec![1], false));
    s
}

fn close(a: f64, b: f64, relative: f64) -> bool {
    (a - b).abs() <= b.abs() * relative
}

// --- the surface ------------------------------------------------------------------------------------------------------------------

#[test]
fn the_points_are_unit_vectors_spread_evenly_over_the_sphere() {
    let p = sphere_points(960);
    assert_eq!(p.len(), 960);
    for u in &p {
        assert!((u[0] * u[0] + u[1] * u[1] + u[2] * u[2] - 1.0).abs() < 1e-12);
    }
    let mean = |k: usize| p.iter().map(|u| u[k]).sum::<f64>() / p.len() as f64;
    assert!(mean(0).abs() < 2e-3 && mean(1).abs() < 2e-3 && mean(2).abs() < 2e-3);
    // an even spread has a third of the squared height on average
    let z2 = p.iter().map(|u| u[2] * u[2]).sum::<f64>() / p.len() as f64;
    assert!((z2 - 1.0 / 3.0).abs() < 2e-3, "{z2}");
}

#[test]
fn a_lone_atom_has_the_area_of_the_sphere_of_its_radius_and_the_probe() {
    let a = atom_areas(&[[1.0, 2.0, 3.0]], &[1.7], 1.4, 960);
    assert!(close(a[0], 4.0 * PI * 3.1 * 3.1, 1e-12), "{}", a[0]);
    assert!(atom_areas(&[], &[], 1.4, 960).is_empty());
}

#[test]
fn two_atoms_that_meet_bury_the_cap_that_geometry_says() {
    // sphere 1 of radius R1, sphere 2 of R2, centres d apart: the points of sphere 1 that sphere 2 holds are the cap where
    // cos(theta) > (R1^2 + d^2 - R2^2) / (2 R1 d), whose area is 2 pi R1^2 (1 - cos(theta))
    let hidden = |r1: f64, r2: f64, d: f64| {
        2.0 * PI * r1 * r1 * (1.0 - (r1 * r1 + d * d - r2 * r2) / (2.0 * r1 * d))
    };
    for (ra, rb, d) in [(1.7, 1.7, 3.5), (1.7, 1.55, 3.0), (1.52, 1.8, 2.4)] {
        let (r1, r2) = (ra + PROBE, rb + PROBE);
        let areas = atom_areas(&[[0.0; 3], [d, 0.0, 0.0]], &[ra, rb], PROBE, 6000);
        let want1 = 4.0 * PI * r1 * r1 - hidden(r1, r2, d);
        let want2 = 4.0 * PI * r2 * r2 - hidden(r2, r1, d);
        assert!(
            close(areas[0], want1, 0.01),
            "{ra} {rb} {d}: {} against {want1}",
            areas[0]
        );
        assert!(
            close(areas[1], want2, 0.01),
            "{ra} {rb} {d}: {} against {want2}",
            areas[1]
        );
    }
}

#[test]
fn atoms_too_far_apart_for_their_spheres_to_meet_do_not_touch_each_others_area() {
    let (r1, r2) = (1.7 + PROBE, 1.55 + PROBE);
    let d = r1 + r2 + 0.01;
    let a = atom_areas(&[[0.0; 3], [d, 0.0, 0.0]], &[1.7, 1.55], PROBE, 960);
    assert!(close(a[0], 4.0 * PI * r1 * r1, 1e-12) && close(a[1], 4.0 * PI * r2 * r2, 1e-12));
}

/// Every atom against every other, with no grid: the answer the grid must give.
fn brute(coords: &[[f64; 3]], radii: &[f64]) -> Vec<f64> {
    let spots = sphere_points(POINTS);
    (0..coords.len())
        .map(|i| {
            let ri = radii[i] + PROBE;
            let free = spots
                .iter()
                .filter(|u| {
                    let p = [
                        coords[i][0] + ri * u[0],
                        coords[i][1] + ri * u[1],
                        coords[i][2] + ri * u[2],
                    ];
                    (0..coords.len())
                        .filter(|&j| j != i)
                        .all(|j| dist2(&p, &coords[j]) >= (radii[j] + PROBE).powi(2))
                })
                .count();
            4.0 * PI * ri * ri * free as f64 / spots.len() as f64
        })
        .collect()
}

fn heavy(s: &Structure) -> (Vec<[f64; 3]>, Vec<f64>) {
    s.atoms
        .iter()
        .filter_map(|a| radius(&a.element).map(|r| (a.pos.map(f64::from), r)))
        .unzip()
}

#[test]
fn the_grid_gives_what_comparing_every_pair_of_atoms_gives() {
    let (coords, radii) = heavy(&ubiquitin());
    let fast = atom_areas(&coords, &radii, PROBE, POINTS);
    let slow = brute(&coords, &radii);
    for (k, (f, s)) in fast.iter().zip(&slow).enumerate() {
        assert!((f - s).abs() < 1e-9, "atom {k}: {f} against {s}");
    }
}

#[test]
fn the_total_surface_of_ubiquitin_is_what_an_independent_program_gives() {
    // worked out by a separate implementation (own parser, numpy, the same points and radii): 602 atoms, 4872.4617
    let (coords, radii) = heavy(&ubiquitin());
    assert_eq!(coords.len(), 602);
    let total: f64 = atom_areas(&coords, &radii, PROBE, POINTS).iter().sum();
    assert!(close(total, 4872.4617, 5e-4), "{total}");
}

#[test]
fn the_surface_does_not_depend_on_how_the_molecule_is_turned_or_moved() {
    let (coords, radii) = heavy(&ubiquitin());
    let (ca, sa, cb, sb) = (0.6f64.cos(), 0.6f64.sin(), 1.1f64.cos(), 1.1f64.sin());
    let turned: Vec<[f64; 3]> = coords
        .iter()
        .map(|p| {
            let (x, y) = (p[0] * ca - p[1] * sa, p[0] * sa + p[1] * ca);
            let (x, z) = (x * cb + p[2] * sb, -x * sb + p[2] * cb);
            [x + 12.0, y - 7.5, z + 3.25]
        })
        .collect();
    let total = |c: &[[f64; 3]]| atom_areas(c, &radii, PROBE, POINTS).iter().sum::<f64>();
    assert!(
        close(total(&turned), total(&coords), 3e-3),
        "{} against {}",
        total(&turned),
        total(&coords)
    );
}

#[test]
fn water_loving_residues_are_more_exposed_than_the_ones_in_the_core() {
    let s = ubiquitin();
    let (coords, radii) = heavy(&s);
    let areas = atom_areas(&coords, &radii, PROBE, POINTS);
    // atoms in the order of the structure's, hydrogens (none here) left out: the same index as `heavy`
    let mean_of = |names: &[&str]| {
        let mut sum = 0.0;
        let mut n = 0;
        for r in s
            .residues
            .iter()
            .filter(|r| names.contains(&r.name.as_str()))
        {
            sum += r.atoms.iter().map(|&a| areas[a]).sum::<f64>();
            n += 1;
        }
        sum / n as f64
    };
    let core = mean_of(&["LEU", "ILE", "VAL", "PHE"]);
    let outside = mean_of(&["LYS", "GLU", "ASP", "ARG"]);
    assert!(core < outside * 0.7, "core {core}, outside {outside}");
}

#[test]
fn hydrogens_are_left_out_and_every_other_element_has_a_radius() {
    assert_eq!(radius("H"), None);
    assert_eq!(radius("d"), None);
    assert_eq!(radius(" N "), Some(1.55));
    assert_eq!(radius("SE"), Some(1.90));
    assert_eq!(
        radius("ZN"),
        Some(1.70),
        "an element not listed has the radius of carbon"
    );
    let mut s = cut_ubiquitin();
    let with = analyse(&s, &['A'], &['B']).unwrap();
    let extra = s.atoms.len();
    s.atoms.push(atom(
        "HA",
        "H",
        [
            s.atoms[0].pos[0] + 0.5,
            s.atoms[0].pos[1],
            s.atoms[0].pos[2],
        ],
        0,
    ));
    s.residues[0].atoms.push(extra);
    let again = analyse(&s, &['A'], &['B']).unwrap();
    assert_eq!(again.area_bound, with.area_bound);
}

// --- the interface ---------------------------------------------------------------------------------------------------------------

#[test]
fn what_two_atoms_bury_is_the_cap_each_loses_and_the_interface_is_half_of_it() {
    let s = pair(("ALA", "C", "CB"), ("ALA", "C", "CB"), 3.5);
    let i = analyse(&s, &['A'], &['B']).unwrap();
    let r = 1.7 + PROBE;
    let cap = 2.0 * PI * r * r * (1.0 - 3.5 / (2.0 * r));
    assert!(
        close(i.a[0].buried(), cap, 0.02),
        "{} against {cap}",
        i.a[0].buried()
    );
    assert!(close(i.b[0].buried(), cap, 0.02));
    assert!(close(i.buried_total(), 2.0 * cap, 0.02));
    assert!(close(i.interface_area(), cap, 0.02));
    assert!(close(i.area_a, 4.0 * PI * r * r, 1e-12) && close(i.area_b, 4.0 * PI * r * r, 1e-12));
    assert_eq!(i.interface_residues(&i.a).len(), 1);
}

#[test]
fn residues_with_atoms_within_four_and_a_half_angstroms_are_in_contact() {
    let near = analyse(
        &pair(("ALA", "C", "CB"), ("ALA", "C", "CB"), 4.4),
        &['A'],
        &['B'],
    )
    .unwrap();
    assert_eq!(
        near.contacts,
        [Contact {
            a: 0,
            b: 1,
            atom_pairs: 1,
            closest: f64::from(4.4f32)
        }]
    );
    let far = analyse(
        &pair(("ALA", "C", "CB"), ("ALA", "C", "CB"), 4.6),
        &['A'],
        &['B'],
    )
    .unwrap();
    assert!(far.contacts.is_empty());
    // out of contact, the spheres of the atoms and the probe still meet: some surface is buried without a contact
    assert!(far.buried_total() > 0.0);
}

#[test]
fn a_donor_and_an_acceptor_close_enough_are_a_candidate_hydrogen_bond() {
    let bond = |a: (&str, &str, &str), b: (&str, &str, &str), d: f32| {
        analyse(&pair(a, b, d), &['A'], &['B'])
            .unwrap()
            .hbonds
            .len()
    };
    // a backbone nitrogen gives a hydrogen, a backbone oxygen takes it
    assert_eq!(bond(("ALA", "N", "N"), ("ALA", "O", "O"), 2.9), 1);
    assert_eq!(
        bond(("ALA", "O", "O"), ("ALA", "N", "N"), 2.9),
        1,
        "either side may give"
    );
    assert_eq!(
        bond(("ALA", "N", "N"), ("ALA", "O", "O"), 3.5),
        1,
        "3.5 is still near enough"
    );
    assert_eq!(bond(("ALA", "N", "N"), ("ALA", "O", "O"), 3.6), 0);
    assert_eq!(
        bond(("ALA", "N", "N"), ("ALA", "N", "N"), 2.9),
        0,
        "two donors"
    );
    assert_eq!(
        bond(("ALA", "O", "O"), ("ALA", "O", "O"), 2.9),
        0,
        "two acceptors"
    );
    assert_eq!(
        bond(("PRO", "N", "N"), ("ALA", "O", "O"), 2.9),
        0,
        "a proline nitrogen has no hydrogen"
    );
    assert_eq!(
        bond(("ALA", "C", "CB"), ("ALA", "O", "O"), 2.9),
        0,
        "a carbon neither gives nor takes"
    );
    // side chains by name
    assert_eq!(bond(("LYS", "N", "NZ"), ("ASP", "O", "OD1"), 2.8), 1);
    assert_eq!(
        bond(("SER", "O", "OG"), ("THR", "O", "OG1"), 2.8),
        1,
        "both give and take"
    );
    assert_eq!(
        bond(("ASP", "O", "OD1"), ("GLU", "O", "OE2"), 2.8),
        0,
        "two carboxylates"
    );
    // a ligand's nitrogens and oxygens are taken for either
    let mut s = pair(("ALA", "N", "N"), ("NAG", "O", "O6"), 2.8);
    s.residues[1].hetero = true;
    assert_eq!(analyse(&s, &['A'], &['B']).unwrap().hbonds.len(), 1);
}

#[test]
fn the_hydrogen_bonds_are_listed_closest_first_with_the_atoms_of_each_side() {
    let mut s = Structure::default();
    for (k, d) in [3.3f32, 2.7, 3.0].into_iter().enumerate() {
        let y = k as f32 * 20.0;
        s.atoms.push(atom("N", "N", [0.0, y, 0.0], 2 * k));
        s.atoms.push(atom("O", "O", [d, y, 0.0], 2 * k + 1));
        s.residues
            .push(residue("ALA", 'A', k as i32 + 1, vec![2 * k], false));
        s.residues
            .push(residue("ALA", 'B', k as i32 + 1, vec![2 * k + 1], false));
    }
    let i = analyse(&s, &['A'], &['B']).unwrap();
    let got: Vec<(usize, usize)> = i.hbonds.iter().map(|h| (h.a_atom, h.b_atom)).collect();
    assert_eq!(got, [(2, 3), (4, 5), (0, 1)]);
    assert!(i.hbonds.windows(2).all(|w| w[0].distance <= w[1].distance));
}

#[test]
fn ubiquitin_cut_in_two_gives_the_interface_an_independent_program_finds() {
    // the same cut (chain A is residues 1 to 38) worked out by the separate implementation
    let s = cut_ubiquitin();
    let i = analyse(&s, &['A'], &['B']).unwrap();
    assert!(close(i.area_a, 3220.3827, 5e-4), "{}", i.area_a);
    assert!(close(i.area_b, 3673.2258, 5e-4), "{}", i.area_b);
    assert!(close(i.area_bound, 4872.4617, 5e-4), "{}", i.area_bound);
    assert!(
        close(i.buried_total(), 2021.1468, 2e-3),
        "{}",
        i.buried_total()
    );
    let buried_of = |n: i32| {
        i.a.iter()
            .chain(&i.b)
            .find(|b| s.residues[b.residue].number == n)
            .map(Buried::buried)
            .unwrap()
    };
    assert!((buried_of(23) - 118.8708).abs() < 0.5, "{}", buried_of(23));
    assert!((buried_of(38) - 55.0363).abs() < 0.5, "{}", buried_of(38));
    assert!((buried_of(39) - 48.0579).abs() < 0.5, "{}", buried_of(39));
    assert!((buried_of(40) - 39.4794).abs() < 0.5, "{}", buried_of(40));
    assert_eq!(
        buried_of(33),
        0.0,
        "a residue far from the cut buries nothing"
    );
    assert_eq!(buried_of(45), 0.0);
    assert_eq!(i.contacts.len(), 65);
    assert_eq!(i.contacts.iter().map(|c| c.atom_pairs).sum::<usize>(), 377);
    // 21 candidates: a rule table written afresh in another program, on the same coordinates, finds the same number, and the
    // three closest pairs. The closest, 2.2 angstroms, is residue 38's oxygen and residue 39's nitrogen: a peptide bond, which
    // the made-up cut passes through (in a real complex the two sides are not joined)
    assert_eq!(i.hbonds.len(), 21);
    let pair_of = |h: &HBond| {
        let (a, b) = (&s.atoms[h.a_atom], &s.atoms[h.b_atom]);
        (
            s.residues[a.residue].number,
            a.name.clone(),
            s.residues[b.residue].number,
            b.name.clone(),
        )
    };
    let wanted = [
        (2.227, (38, "O", 39, "N")),
        (2.783, (2, "O", 64, "N")),
        (2.819, (4, "O", 67, "N")),
    ];
    for (h, (distance, (an, aname, bn, bname))) in i.hbonds.iter().zip(wanted) {
        assert!((h.distance - distance).abs() < 0.002, "{}", h.distance);
        assert_eq!(pair_of(h), (an, aname.to_string(), bn, bname.to_string()));
    }
}

#[test]
fn what_the_two_sides_bury_adds_up_to_what_the_complex_hides_and_does_not_depend_on_which_is_first()
{
    let s = cut_ubiquitin();
    let ab = analyse(&s, &['A'], &['B']).unwrap();
    let ba = analyse(&s, &['B'], &['A']).unwrap();
    let sum: f64 = ab.a.iter().chain(&ab.b).map(Buried::buried).sum();
    assert!(
        close(sum, ab.buried_total(), 1e-9),
        "{sum} against {}",
        ab.buried_total()
    );
    assert!(close(ba.buried_total(), ab.buried_total(), 1e-12));
    assert_eq!(ba.contacts.len(), ab.contacts.len());
    assert_eq!(ba.hbonds.len(), ab.hbonds.len());
    assert_eq!(ba.a.len(), ab.b.len());
    // each hydrogen bond has one atom in each group
    for h in &ab.hbonds {
        assert!(s.residues[s.atoms[h.a_atom].residue].number <= 38);
        assert!(s.residues[s.atoms[h.b_atom].residue].number > 38);
    }
}

#[test]
fn the_residues_at_the_interface_are_those_that_bury_something_and_the_most_buried_come_first() {
    let s = cut_ubiquitin();
    let i = analyse(&s, &['A'], &['B']).unwrap();
    let at = i.interface_residues(&i.a);
    assert!(at.len() > 5 && at.len() < i.a.len());
    assert!(at.iter().all(|r| r.buried() >= MIN_BURIED));
    assert!(at.windows(2).all(|w| w[0].buried() >= w[1].buried()));
    assert_eq!(s.residues[at[0].residue].number, 23, "the most buried");
}

#[test]
fn a_third_chain_is_not_looked_at_and_asking_for_it_is_the_way_to_include_it() {
    let mut s = cut_ubiquitin();
    for r in &mut s.residues {
        if r.number > 60 {
            r.chain = 'C';
        }
    }
    let two = analyse(&s, &['A'], &['B']).unwrap();
    assert!(two.b.iter().all(|r| s.residues[r.residue].chain == 'B'));
    let with = analyse(&s, &['A'], &['B', 'C']).unwrap();
    assert!(with.b.iter().any(|r| s.residues[r.residue].chain == 'C'));
    assert!(with.area_bound > two.area_bound);
    assert_eq!(with.b_chains, ['B', 'C']);
}

#[test]
fn a_side_with_no_chain_a_chain_on_both_sides_and_a_chain_that_is_not_there_are_said_so() {
    let s = cut_ubiquitin();
    assert!(analyse(&s, &[], &['B'])
        .unwrap_err()
        .contains("at least one chain for each side"));
    assert!(analyse(&s, &['A'], &[]).is_err());
    assert!(analyse(&s, &['A', 'B'], &['B'])
        .unwrap_err()
        .contains("Chain B is on both sides"));
    assert!(analyse(&s, &['A'], &['Q'])
        .unwrap_err()
        .contains("no chain Q"));
}

#[test]
fn the_table_lists_the_interface_residues_of_each_group_with_their_contacts_and_hydrogen_bonds() {
    let s = cut_ubiquitin();
    let i = analyse(&s, &['A'], &['B']).unwrap();
    let rows = i.rows(&s);
    let first: Vec<&Row> = rows.iter().filter(|r| r.side == 1).collect();
    let second: Vec<&Row> = rows.iter().filter(|r| r.side == 2).collect();
    assert_eq!(first.len(), i.interface_residues(&i.a).len());
    assert_eq!(second.len(), i.interface_residues(&i.b).len());
    assert!(first.iter().all(|r| r.chain == 'A') && second.iter().all(|r| r.chain == 'B'));
    assert!(first.windows(2).all(|w| w[0].buried >= w[1].buried));
    assert_eq!((first[0].number, first[0].name.as_str()), (23, "ILE"));
    // each contact has a residue on each side; the table lists the residues that bury something, and a residue can touch the
    // other group and bury less than MIN_BURIED, so a side's contacts add up to those of its listed residues
    let listed_a: Vec<usize> = first.iter().map(|r| r.residue).collect();
    let listed_b: Vec<usize> = second.iter().map(|r| r.residue).collect();
    assert_eq!(
        first.iter().map(|r| r.contacts).sum::<usize>(),
        i.contacts
            .iter()
            .filter(|c| listed_a.contains(&c.a))
            .count()
    );
    assert_eq!(
        second.iter().map(|r| r.contacts).sum::<usize>(),
        i.contacts
            .iter()
            .filter(|c| listed_b.contains(&c.b))
            .count()
    );
    assert!(i.contacts.len() >= first.iter().map(|r| r.contacts).sum::<usize>());
    assert_eq!(
        first.iter().map(|r| r.hbonds).sum::<usize>(),
        i.hbonds
            .iter()
            .filter(|h| listed_a.contains(&s.atoms[h.a_atom].residue))
            .count()
    );
    assert_eq!(
        second.iter().map(|r| r.hbonds).sum::<usize>(),
        i.hbonds
            .iter()
            .filter(|h| listed_b.contains(&s.atoms[h.b_atom].residue))
            .count()
    );
    assert!(rows.iter().all(|r| r.buried >= MIN_BURIED));
    // every row is the residue it says it is
    for r in &rows {
        let res = &s.residues[r.residue];
        assert_eq!(
            (res.chain, res.name.as_str(), res.number),
            (r.chain, r.name.as_str(), r.number)
        );
    }
}

#[test]
fn the_table_goes_to_a_spreadsheet_as_text() {
    let rows = vec![
        Row {
            side: 1,
            residue: 0,
            chain: 'H',
            name: "TYR".into(),
            number: 33,
            insertion: ' ',
            buried: 91.234,
            contacts: 4,
            hbonds: 2,
        },
        Row {
            side: 2,
            residue: 1,
            chain: 'A',
            name: "LYS".into(),
            number: 97,
            insertion: 'A',
            buried: 12.0,
            contacts: 1,
            hbonds: 0,
        },
    ];
    assert_eq!(
        tsv(&rows),
        "side\tchain\tresidue\tnumber\tburied_A2\tcontacts\thbond_candidates\n1\tH\tTYR\t33\t91.2\t4\t2\n2\tA\tLYS\t97A\t12.0\t1\t0\n"
    );
    assert_eq!(tsv(&[]).lines().count(), 1, "a head line alone");
}

#[test]
fn an_insertion_code_is_kept_in_the_table_and_the_text() {
    let mut s = pair(("ALA", "C", "CB"), ("ALA", "C", "CB"), 3.8);
    s.residues[0].insertion = 'A';
    let i = analyse(&s, &['A'], &['B']).unwrap();
    let rows = i.rows(&s);
    assert_eq!(rows[0].insertion, 'A');
    assert_eq!(rows[0].number_text(), "1A");
    assert_eq!(rows[1].number_text(), "1");
    assert!(tsv(&rows).contains("\tALA\t1A\t"));
    assert!(tsv(&rows).contains("\tALA\t1\t"));
}

#[test]
fn an_ion_named_like_nitrogen_is_not_a_hydrogen_bond_partner_but_a_ligands_oxygen_is() {
    let mut s = pair(("ALA", "O", "O"), ("NA", "NA", "NA"), 2.5);
    s.residues[1].hetero = true;
    assert!(
        analyse(&s, &['A'], &['B']).unwrap().hbonds.is_empty(),
        "a sodium ion by its element"
    );
    let mut s = pair(("ALA", "N", "N"), ("NAG", "O", "O6"), 2.8);
    s.residues[1].hetero = true;
    assert_eq!(analyse(&s, &['A'], &['B']).unwrap().hbonds.len(), 1);
}

#[test]
fn an_atom_with_no_coordinates_is_refused_and_named_not_measured_as_if_it_were_not_there() {
    for bad in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
        let mut s = pair(("ALA", "C", "CB"), ("GLY", "C", "CA"), 3.8);
        s.atoms[1].pos[1] = bad;
        let why = analyse(&s, &['A'], &['B']).unwrap_err();
        assert!(
            why.contains("GLY 1")
                && why.contains("chain B")
                && why.contains("no usable coordinates"),
            "{why}"
        );
        // a chain that is left out may hold such an atom: it is not looked at
        let mut s = pair(("ALA", "C", "CB"), ("GLY", "C", "CA"), 3.8);
        s.atoms.push(atom("CA", "C", [f32::NAN, 0.0, 0.0], 2));
        s.residues.push(residue("ALA", 'C', 1, vec![2], false));
        assert!(analyse(&s, &['A'], &['B']).is_ok());
    }
}

#[test]
fn the_grid_is_right_for_atoms_of_very_different_sizes_too() {
    // radii from a hydrogen-bonding oxygen to a large ion, scattered in a box that makes them overlap in all sorts of ways
    let sizes = [1.2, 1.55, 1.7, 1.9, 2.75, 3.4];
    let mut seed = 12345u64;
    let mut next = || {
        seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (seed >> 33) as f64 / (1u64 << 31) as f64
    };
    let coords: Vec<[f64; 3]> = (0..220)
        .map(|_| [next() * 18.0, next() * 18.0, next() * 18.0])
        .collect();
    let radii: Vec<f64> = (0..220).map(|k| sizes[k % sizes.len()]).collect();
    let fast = atom_areas(&coords, &radii, PROBE, POINTS);
    let slow = brute(&coords, &radii);
    for (k, (f, s)) in fast.iter().zip(&slow).enumerate() {
        assert!(
            (f - s).abs() < 1e-9,
            "atom {k} (radius {}): {f} against {s}",
            radii[k]
        );
    }
    // some atoms are buried and some are free, so the test has something to tell apart
    assert!(fast.iter().any(|&a| a < 1.0) && fast.iter().any(|&a| a > 100.0));
}

// --- a real complex ---------------------------------------------------------------------------------------------------------------

/// The antibody D1.3 Fv (chains A and B, light and heavy) bound to hen egg-white lysozyme (chain C): PDB 1VFB (Bhat et al. 1994),
/// as the PDB gives it, kept in `assets/structure/1vfb.pdb`.
const ANTIBODY: &str = include_str!("../../../assets/structure/1vfb.pdb");

#[test]
fn an_antibody_and_its_antigen_from_the_pdb() {
    let s = pdb::parse(ANTIBODY).unwrap();
    assert_eq!(s.chains(), ['A', 'B', 'C']);
    assert_eq!(s.atoms.len(), 3408);
    let i = analyse(&s, &['A', 'B'], &['C']).unwrap();
    // as measured on the first run, so that a change in the method shows: 1372 A2 buried in all, 39 contacts, 18 candidates
    assert!(
        (i.buried_total() - 1372.0).abs() < 8.0,
        "{}",
        i.buried_total()
    );
    assert_eq!(i.contacts.len(), 39);
    assert_eq!(i.hbonds.len(), 18);
    let rows = i.rows(&s);
    // what is known of this complex (D1.3 and lysozyme, Amit et al. 1986 and the structures after it), and not worked out by this
    // program: the paratope is in the CDRs (Kabat: light 24-34, 50-56, 89-97; heavy 31-35, 50-65, 95-102), taken with a residue to
    // spare each side; the epitope is two segments of lysozyme, 18-27 and 116-129; the interface is some hundreds of square
    // angstroms on each side
    let in_cdr = |chain: char, n: i32| {
        let windows: &[(i32, i32)] = if chain == 'A' {
            &[(23, 35), (49, 57), (88, 98)]
        } else {
            &[(30, 36), (49, 66), (94, 103)]
        };
        windows.iter().any(|&(from, to)| (from..=to).contains(&n))
    };
    let paratope: Vec<&Row> = rows.iter().filter(|r| r.side == 1).collect();
    let in_cdrs = paratope
        .iter()
        .filter(|r| in_cdr(r.chain, r.number))
        .count();
    println!(
        "{in_cdrs} of {} residues of the paratope are in a CDR",
        paratope.len()
    );
    assert!(
        in_cdrs * 100 >= paratope.len() * 85,
        "{in_cdrs} of {}",
        paratope.len()
    );
    let epitope: Vec<&Row> = rows.iter().filter(|r| r.side == 2).collect();
    let in_segments = epitope
        .iter()
        .filter(|r| (18..=27).contains(&r.number) || (116..=129).contains(&r.number))
        .count();
    assert!(
        in_segments * 100 >= epitope.len() * 85,
        "{in_segments} of {}",
        epitope.len()
    );
    let buries = |chain: char, name: &str, n: i32| {
        rows.iter()
            .find(|r| (r.chain, r.name.as_str(), r.number) == (chain, name, n))
            .map(|r| r.buried)
    };
    for (chain, name, n) in [
        ('A', "TYR", 32),
        ('A', "TYR", 50),
        ('A', "TRP", 92),
        ('B', "TRP", 52),
        ('B', "ASP", 100),
        ('B', "TYR", 101),
        ('C', "GLN", 121),
        ('C', "ASP", 119),
        ('C', "LYS", 116),
        ('C', "GLY", 117),
    ] {
        assert!(
            buries(chain, name, n).is_some_and(|b| b > 20.0),
            "{chain} {name} {n}: {:?}",
            buries(chain, name, n)
        );
    }
    let deepest = epitope
        .iter()
        .max_by(|x, y| x.buried.total_cmp(&y.buried))
        .unwrap();
    assert_eq!(
        (deepest.name.as_str(), deepest.number),
        ("GLN", 121),
        "the glutamine that goes into the paratope"
    );
    assert!(
        (500.0..900.0).contains(&i.interface_area()),
        "{}",
        i.interface_area()
    );
}
