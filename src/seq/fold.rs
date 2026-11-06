//! A structure prepared for viewing: its protein residues in unbroken runs, the secondary structure of each, how
//! confident a predictor was about it, and the sphere that holds it all.

use super::pdb::{distance, Backbone, Structure};
use super::secstruct::{assign, Ss};

/// A bond between two alpha carbons is about 3.8 angstroms; further apart than this the chain is broken there.
const CA_CA_LIMIT: f32 = 4.6;

#[derive(Clone, Debug, PartialEq)]
pub struct Res {
    /// Index into `Structure::residues`.
    pub index: usize,
    pub letter: u8,
    pub number: i32,
    pub chain: char,
    /// Position of the chain among the structure's chains, for colouring.
    pub chain_index: usize,
    pub ca: [f32; 3],
    pub ss: Ss,
    /// 0 to 100 when the structure looks like a prediction, else none.
    pub confidence: Option<f32>,
}

#[derive(Clone, Debug)]
pub struct Fold {
    pub structure: Structure,
    pub residues: Vec<Res>,
    /// The backbone atoms of each of `residues`.
    pub backbone: Vec<Backbone>,
    /// Runs of consecutive `residues` that belong to one unbroken chain.
    pub runs: Vec<Vec<usize>>,
    pub centre: [f32; 3],
    /// The sphere about `centre` that holds every atom.
    pub radius: f32,
    /// The smaller sphere about `centre` that holds every alpha carbon: what a ribbon is fitted to.
    pub backbone_radius: f32,
    pub predicted: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Summary {
    pub residues: usize,
    pub chains: usize,
    pub atoms: usize,
    pub helix: f64,
    pub strand: f64,
    pub confidence: Option<f64>,
}

impl Fold {
    pub fn new(structure: Structure) -> Fold {
        let backbone = structure.backbone();
        let atoms: Vec<_> = backbone
            .iter()
            .map(|(_, b, linked)| (*b, *linked))
            .collect();
        let states = assign(&atoms);
        let chains = structure.chains();
        let predicted = structure.looks_predicted();
        let residues: Vec<Res> = backbone
            .iter()
            .zip(&states)
            .map(|((index, b, _), &ss)| {
                let r = &structure.residues[*index];
                Res {
                    index: *index,
                    letter: super::pdb::letter(&r.name),
                    number: r.number,
                    chain: r.chain,
                    chain_index: chains.iter().position(|&c| c == r.chain).unwrap_or(0),
                    ca: b.ca,
                    ss,
                    confidence: predicted.then(|| structure.confidence(*index)).flatten(),
                }
            })
            .collect();
        let mut runs: Vec<Vec<usize>> = Vec::new();
        for (k, r) in residues.iter().enumerate() {
            let continues = k > 0
                && residues[k - 1].chain == r.chain
                && distance(residues[k - 1].ca, r.ca) <= CA_CA_LIMIT;
            match runs.last_mut() {
                Some(run) if continues => run.push(k),
                _ => runs.push(vec![k]),
            }
        }
        let (centre, radius) = structure.centre_and_radius();
        let backbone_radius = residues
            .iter()
            .map(|r| distance(r.ca, centre))
            .fold(1.0, f32::max);
        let backbone = backbone.iter().map(|(_, b, _)| *b).collect();
        Fold {
            structure,
            residues,
            backbone,
            runs,
            centre,
            radius,
            backbone_radius,
            predicted,
        }
    }

    pub fn summary(&self) -> Summary {
        let n = self.residues.len().max(1) as f64;
        let share = |ss: Ss| 100.0 * self.residues.iter().filter(|r| r.ss == ss).count() as f64 / n;
        let known: Vec<f64> = self
            .residues
            .iter()
            .filter_map(|r| r.confidence.map(f64::from))
            .collect();
        Summary {
            residues: self.residues.len(),
            chains: self.structure.chains().len(),
            atoms: self.structure.atoms.len(),
            helix: share(Ss::Helix),
            strand: share(Ss::Strand),
            confidence: (!known.is_empty()).then(|| known.iter().sum::<f64>() / known.len() as f64),
        }
    }

    /// The sequence of each chain, one letter per residue, in file order.
    pub fn chain_sequences(&self) -> Vec<(char, Vec<u8>)> {
        let mut out: Vec<(char, Vec<u8>)> = Vec::new();
        for r in &self.residues {
            match out.last_mut() {
                Some((chain, letters)) if *chain == r.chain => letters.push(r.letter),
                _ => out.push((r.chain, vec![r.letter])),
            }
        }
        out
    }

    /// Pairs of atoms that are bonded: closer than a covalent bond is long (1.9 angstroms, 2.2 when one is
    /// sulphur). Found with a grid, so a large structure is quick.
    pub fn bonds(&self) -> Vec<(usize, usize)> {
        let atoms = &self.structure.atoms;
        let cell = |p: [f32; 3]| {
            [
                (p[0] / BOND_CELL).floor() as i32,
                (p[1] / BOND_CELL).floor() as i32,
                (p[2] / BOND_CELL).floor() as i32,
            ]
        };
        let mut grid: std::collections::HashMap<[i32; 3], Vec<usize>> =
            std::collections::HashMap::new();
        for (i, a) in atoms.iter().enumerate() {
            grid.entry(cell(a.pos)).or_default().push(i);
        }
        let mut bonds = Vec::with_capacity(atoms.len());
        for (i, a) in atoms.iter().enumerate() {
            let c = cell(a.pos);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let Some(near) = grid.get(&[c[0] + dx, c[1] + dy, c[2] + dz]) else {
                            continue;
                        };
                        for &j in near.iter().filter(|&&j| j > i) {
                            let sulphur = a.element == "S" || atoms[j].element == "S";
                            let limit = if sulphur { BOND_SULPHUR } else { BOND_LIMIT };
                            if distance(a.pos, atoms[j].pos) < limit {
                                bonds.push((i, j));
                            }
                        }
                    }
                }
            }
        }
        bonds
    }
}

/// Longest covalent bond between two heavy atoms, and between one of them and sulphur.
const BOND_LIMIT: f32 = 1.9;
const BOND_SULPHUR: f32 = 2.2;
/// The grid's cell is as wide as the longest bond, so bonded atoms are in the same or a neighbouring cell.
const BOND_CELL: f32 = BOND_SULPHUR;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::pdb::parse;

    const UBIQUITIN: &str = include_str!("../../assets/structure/1ubq.pdb");

    #[test]
    fn ubiquitin_is_one_run_of_seventy_six_residues_with_helix_and_strands() {
        let f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        assert_eq!(f.residues.len(), 76);
        assert_eq!(f.runs.len(), 1);
        assert_eq!(f.runs[0].len(), 76);
        assert!(!f.predicted);
        assert!(f.residues.iter().all(|r| r.confidence.is_none()));
        let s = f.summary();
        assert_eq!((s.residues, s.chains, s.atoms), (76, 1, 602));
        // PyDSSP: 18 helix and 26 strand residues of 76
        assert!((s.helix - 100.0 * 18.0 / 76.0).abs() < 3.0, "{}", s.helix);
        assert!((s.strand - 100.0 * 26.0 / 76.0).abs() < 3.0, "{}", s.strand);
        assert_eq!(s.confidence, None);
        let chains = f.chain_sequences();
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].0, 'A');
        assert!(chains[0].1.starts_with(
            b"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
        ));
        assert_eq!(f.residues[0].letter, b'M');
        assert_eq!(f.residues[75].number, 76);
        // the alpha carbons sit inside the sphere of all the atoms, and not far inside it
        assert!(f.backbone_radius < f.radius);
        assert!(
            f.backbone_radius > 0.6 * f.radius,
            "{} {}",
            f.backbone_radius,
            f.radius
        );
    }

    #[test]
    fn bonds_join_neighbouring_atoms_once_and_each_atom_has_a_sensible_number() {
        let f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        let bonds = f.bonds();
        let mut degree = vec![0usize; f.structure.atoms.len()];
        for &(a, b) in &bonds {
            assert!(a < b);
            degree[a] += 1;
            degree[b] += 1;
        }
        let mut unique = bonds.clone();
        unique.sort_unstable();
        unique.dedup();
        assert_eq!(unique.len(), bonds.len(), "each bond is listed once");
        // heavy atoms have one to four bonded neighbours (a few may touch a fifth in a crowded spot)
        assert!(
            degree.iter().all(|&d| (1..=5).contains(&d)),
            "{:?}",
            degree
                .iter()
                .enumerate()
                .filter(|(_, &d)| !(1..=5).contains(&d))
                .collect::<Vec<_>>()
        );
        // a chain of 602 atoms in 76 residues: 602 minus the number of separate molecules, plus 76 for rings and the like
        assert!((600..=700).contains(&bonds.len()), "{}", bonds.len());
        // the backbone's N to CA bond of the first residue is among them
        let first = &f.structure.residues[0];
        let (n, ca) = (first.atoms[0], first.atoms[1]);
        assert!(bonds.contains(&(n.min(ca), n.max(ca))));
        assert_eq!(f.backbone.len(), f.residues.len());
    }

    #[test]
    fn a_gap_in_the_chain_starts_a_new_run() {
        let line = |serial: u32, number: i32, x: f32, b: f32| {
            let atom = |name: &str, dx: f32| {
                format!(
                    "ATOM  {serial:>5} {name:<4} ALA A{number:>4}    {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}           C",
                    x + dx,
                    0.0,
                    0.0,
                    1.0,
                    b
                )
            };
            [
                atom("N", -1.2),
                atom("CA", 0.0),
                atom("C", 1.2),
                atom("O", 1.2),
            ]
            .join("\n")
        };
        // residues 1 and 2 are 3.8 apart, 3 is far off, 4 follows 3
        let text = [
            line(1, 1, 0.0, 90.0),
            line(5, 2, 3.8, 80.0),
            line(9, 3, 30.0, 70.0),
            line(13, 4, 33.8, 60.0),
        ]
        .join("\n");
        let f = Fold::new(parse(&text).expect("a structure"));
        assert_eq!(f.runs, vec![vec![0, 1], vec![2, 3]]);
        // temperature factors of 0 to 100 on a file with no EXPDTA read as a predictor's confidence
        assert!(f.predicted);
        assert_eq!(f.residues[0].confidence, Some(90.0));
        assert_eq!(f.summary().confidence, Some(75.0));
    }
}
