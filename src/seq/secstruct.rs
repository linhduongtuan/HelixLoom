//! Secondary structure from the backbone atoms, in three states (helix, strand, coil), by the rules of Kabsch &
//! Sander's DSSP (Biopolymers 1983): a hydrogen bond is found by an electrostatic energy between the C=O of one
//! residue and the N-H of another; runs of bonds between residues three, four or five apart make helices, and
//! pairs of bonds between two stretches make the bridges of a sheet. Written from the published definition; checked
//! against PyDSSP's answers on real structures.

use super::pdb::{distance, Backbone};
use std::collections::HashSet;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ss {
    Helix,
    Strand,
    Coil,
}

impl Ss {
    /// DSSP's letter for the state, as the checks against PyDSSP write it.
    #[cfg(test)]
    pub fn letter(self) -> char {
        match self {
            Ss::Helix => 'H',
            Ss::Strand => 'E',
            Ss::Coil => '-',
        }
    }
}

/// The energy (kcal/mol) below which DSSP counts a hydrogen bond.
const BOND_BELOW: f32 = -0.5;
/// q1 * q2 * f of DSSP's electrostatic model: 0.42 e and 0.20 e, with 332 to turn it into kcal/mol.
const COUPLING: f32 = 0.084 * 332.0;
/// How far the N-H hydrogen sits from the nitrogen, in angstroms.
const NH_LENGTH: f32 = 1.0;
/// Pairs of residues whose alpha carbons are further apart than this cannot bond.
const REACH: f32 = 9.0;
/// Residues this close along the chain are not tested for bonds (they are bonded covalently, or in a tight turn).
const MIN_SEPARATION: usize = 3;
/// Atoms closer than this make the energy meaningless (it is cut at DSSP's floor of -9.9).
const MIN_DISTANCE: f32 = 0.5;

fn sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Where the amide hydrogen of the residue after `previous` sits: on the line from the previous residue's O to
/// its C, one bond length from the nitrogen (DSSP's way).
fn hydrogen(previous: &Backbone, n: [f32; 3]) -> [f32; 3] {
    let v = sub(previous.c, previous.o);
    let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
        .sqrt()
        .max(f32::EPSILON);
    [
        n[0] + NH_LENGTH * v[0] / length,
        n[1] + NH_LENGTH * v[1] / length,
        n[2] + NH_LENGTH * v[2] / length,
    ]
}

/// Every hydrogen bond as (acceptor, donor): the C=O of the first to the N-H of the second.
fn bonds(res: &[(Backbone, bool)]) -> HashSet<(usize, usize)> {
    let mut found = HashSet::new();
    for d in 1..res.len() {
        let (donor, linked) = (&res[d].0, res[d].1);
        if !linked || donor.proline {
            continue;
        }
        let h = hydrogen(&res[d - 1].0, donor.n);
        for (a, (acceptor, _)) in res.iter().enumerate() {
            let gap = a.abs_diff(d);
            if gap < MIN_SEPARATION || distance(acceptor.ca, donor.ca) > REACH {
                continue;
            }
            let near = |x: f32| x.max(MIN_DISTANCE);
            let (on, ch) = (
                near(distance(acceptor.o, donor.n)),
                near(distance(acceptor.c, h)),
            );
            let (oh, cn) = (
                near(distance(acceptor.o, h)),
                near(distance(acceptor.c, donor.n)),
            );
            let energy = (COUPLING * (1.0 / on + 1.0 / ch - 1.0 / oh - 1.0 / cn)).max(-9.9);
            if energy < BOND_BELOW {
                found.insert((a, d));
            }
        }
    }
    found
}

/// The state of every residue of `res` (the backbone atoms of each, and whether it follows the one before it).
pub fn assign(res: &[(Backbone, bool)]) -> Vec<Ss> {
    let n = res.len();
    let bonded = bonds(res);
    let has = |a: isize, d: isize| a >= 0 && d >= 0 && bonded.contains(&(a as usize, d as usize));
    // an n-turn at i: the C=O of i bonded to the N-H of i+n
    let turn = |width: isize, i: isize| has(i, i + width);
    // two turns in a row make a helix of that kind over the residues the turns hold (4: alpha, 3: 3-10, 5: pi)
    let helix_start =
        |width: isize, i: usize| i > 0 && turn(width, i as isize - 1) && turn(width, i as isize);
    let mut helix = vec![false; n];
    for (width, span) in [(4isize, 4usize), (3, 3), (5, 5)] {
        // a 3-10 helix is only counted where no alpha helix is already, and a pi helix where neither is
        let already = helix.clone();
        for i in (0..n).filter(|&i| !already[i] && helix_start(width, i)) {
            helix[i..(i + span).min(n)]
                .iter_mut()
                .for_each(|h| *h = true);
        }
    }
    // a bridge between i and j: two bonds in the pattern of a parallel or an antiparallel pair
    let mut strand = vec![false; n];
    for i in 1..n.saturating_sub(1) {
        for j in (i + MIN_SEPARATION)..n.saturating_sub(1) {
            let (a, b) = (i as isize, j as isize);
            let parallel = (has(a - 1, b) && has(b, a + 1)) || (has(b - 1, a) && has(a, b + 1));
            let anti = (has(a, b) && has(b, a)) || (has(a - 1, b + 1) && has(b - 1, a + 1));
            if parallel || anti {
                strand[i] = true;
                strand[j] = true;
            }
        }
    }
    (0..n)
        .map(|i| match (helix[i], strand[i]) {
            (true, _) => Ss::Helix,
            (false, true) => Ss::Strand,
            _ => Ss::Coil,
        })
        .collect()
}

#[cfg(test)]
mod tests;
