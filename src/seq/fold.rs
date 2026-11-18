//! A structure prepared for viewing: its protein residues in unbroken runs, the secondary structure of each, how
//! confident a predictor was about it, and the sphere that holds it all.

use super::conservation;
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
    /// How conserved the residue is in an alignment of its family, 0 (the most variable of the alignment) to 1 (the most
    /// conserved); none until `apply_conservation` has set it, and for a residue that no row of the alignment matches.
    pub conservation: Option<f32>,
    /// The column of that alignment that the residue was lined up with (counted from 0), whether the column has a score or
    /// not: what links a residue to its column. It is set and taken off with `conservation`.
    pub column: Option<usize>,
    /// What the residue buries where its group of chains meets another, in square angstroms (`apply_interface`); none until the
    /// interface has been analysed, and for a residue that is in neither group.
    pub buried: Option<f32>,
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

/// The most rows of an alignment that a chain is lined up with: lining up is the slow part (a table of residues of the chain
/// times residues of the row), and a window must not stop while a hundred rows are tried.
const MOST_CANDIDATES: usize = 4;
/// A row shares this many residues in a run with a chain to be taken for a likely match.
const RUN: usize = 4;

/// The rows (their places in `plain`, which holds each row's residues in upper case, with no gaps) worth lining `letters` up with:
/// those with the most runs of four residues in common (a close relative has many: 70 % identity leaves about a quarter of the
/// runs of four whole), and, for a chain so short or rows so unlike that none shares one, the rows of about the chain's length.
fn candidate_rows(plain: &[Vec<u8>], letters: &[u8]) -> Vec<usize> {
    use std::collections::HashSet;
    let runs: HashSet<&[u8]> = letters.windows(RUN).collect();
    let mut shared: Vec<(usize, usize)> = plain
        .iter()
        .enumerate()
        .map(|(k, row)| (row.windows(RUN).filter(|w| runs.contains(w)).count(), k))
        .filter(|&(n, _)| n > 0)
        .collect();
    if !shared.is_empty() {
        // the most in common first; on a tie the earlier row
        shared.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        return shared
            .into_iter()
            .take(MOST_CANDIDATES)
            .map(|(_, k)| k)
            .collect();
    }
    (0..plain.len())
        .filter(|&k| plain[k].len() * 2 >= letters.len() && plain[k].len() <= letters.len() * 2)
        .take(16)
        .collect()
}

/// A chain of the structure and the row of an alignment that it was lined up with.
#[derive(Clone, Debug, PartialEq)]
pub struct Matched {
    pub chain: char,
    pub row: String,
    /// The share of the matched residues that are the same letter, 0 to 1.
    pub identity: f64,
    /// Residues of the chain that got a conservation, of `residues`.
    pub coloured: usize,
    pub residues: usize,
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
                    conservation: None,
                    column: None,
                    buried: None,
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

    /// Give each residue the conservation of its column in an alignment: each chain is lined up with the row of `rows` that it
    /// matches best (the same protein or a close relative: at least 70 % identical), and a residue gets the score
    /// (`scores`, one for each column, as `conservation::scores` gives them) of the column of the residue it matches, made
    /// to run from 0 at the lowest score of the alignment to 1 at the highest. A chain that no row matches gets nothing.
    pub fn apply_conservation(
        &mut self,
        rows: &[(String, Vec<u8>)],
        scores: &[Option<f32>],
    ) -> Vec<Matched> {
        // so that the colours use the whole range of this alignment, whatever the score and the alignment
        let scored: Vec<f32> = scores.iter().flatten().copied().collect();
        let low = scored.iter().copied().fold(f32::MAX, f32::min);
        let high = scored.iter().copied().fold(f32::MIN, f32::max);
        let scale = |score: f32| {
            if high - low < 1e-6 {
                0.5
            } else {
                (score - low) / (high - low)
            }
        };
        self.clear_conservation();
        let mut chains: Vec<char> = Vec::new();
        for r in &self.residues {
            if !chains.contains(&r.chain) {
                chains.push(r.chain);
            }
        }
        // ponytail: this runs on the window's own thread: at most a few lining-ups (`MOST_CANDIDATES`) for each distinct chain, but
        // a chain of thousands of residues makes each a table of millions of cells (seconds); make it a `Job` if that is felt
        //
        // each row's residues once, in upper case, whatever the case of the alignment (the lining up ignores case, so the quick
        // look at shared runs must too)
        let plain: Vec<Vec<u8>> = rows
            .iter()
            .map(|(_, row)| {
                row.iter()
                    .filter(|&&c| c != b'-' && c != b'.')
                    .map(u8::to_ascii_uppercase)
                    .collect()
            })
            .collect();
        // chains with one sequence (the copies of a homo-oligomer) are lined up once
        let mut done: std::collections::HashMap<Vec<u8>, Option<(usize, conservation::Mapping)>> =
            std::collections::HashMap::new();
        let mut matched = Vec::new();
        for chain in chains {
            let places: Vec<usize> = (0..self.residues.len())
                .filter(|&k| self.residues[k].chain == chain)
                .collect();
            let letters: Vec<u8> = places.iter().map(|&k| self.residues[k].letter).collect();
            let best = done
                .entry(letters.clone())
                .or_insert_with(|| {
                    // the row that covers most of this chain with the same residues (on a tie the more alike, then the first)
                    candidate_rows(&plain, &letters)
                        .into_iter()
                        .filter_map(|k| {
                            conservation::map_to_chain(&plain[k], &letters).map(|m| (k, m))
                        })
                        .reduce(|best, next| {
                            let more =
                                (next.1.same, next.1.identity) > (best.1.same, best.1.identity);
                            if more {
                                next
                            } else {
                                best
                            }
                        })
                })
                .clone();
            let Some((k, mapping)) = best else {
                continue;
            };
            let (name, row) = &rows[k];
            let columns = conservation::columns_of(row);
            let mut coloured = 0;
            for (j, &place) in places.iter().enumerate() {
                let column = mapping.to_row[j].and_then(|i| columns.get(i).copied());
                let score = column.and_then(|k| scores.get(k).copied().flatten());
                self.residues[place].column = column;
                self.residues[place].conservation = score.map(scale);
                coloured += usize::from(score.is_some());
            }
            matched.push(Matched {
                chain,
                row: name.clone(),
                identity: mapping.identity,
                coloured,
                residues: places.len(),
            });
        }
        matched
    }

    /// Set what every residue of the two groups of an interface buries.
    pub fn apply_interface(&mut self, i: &super::interface::Interface) {
        self.clear_interface();
        let buried: std::collections::HashMap<usize, f32> =
            i.a.iter()
                .chain(&i.b)
                .map(|b| (b.residue, b.buried() as f32))
                .collect();
        for r in &mut self.residues {
            r.buried = buried.get(&r.index).copied();
        }
    }

    pub fn clear_interface(&mut self) {
        for r in &mut self.residues {
            r.buried = None;
        }
    }

    /// Take the conservation off every residue, and the column it came from.
    pub fn clear_conservation(&mut self) {
        for r in &mut self.residues {
            (r.conservation, r.column) = (None, None);
        }
    }

    /// The residues (their places in `residues`) that were lined up with column `k` of the alignment: one for each chain
    /// that a row matches, so the copies of a chain are all there.
    pub fn in_column(&self, k: usize) -> Vec<usize> {
        (0..self.residues.len())
            .filter(|&place| self.residues[place].column == Some(k))
            .collect()
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
    const ANTIBODY: &str = include_str!("../../assets/structure/1vfb.pdb");

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

    /// Ubiquitin with `gaps` dashes put in after the residues at those places: a row of an alignment.
    fn ubiquitin_row(gaps_after: &[usize]) -> Vec<u8> {
        let f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        let mut row = Vec::new();
        for (k, r) in f.residues.iter().enumerate() {
            row.push(r.letter);
            // one gap for every time the residue is named
            row.extend(std::iter::repeat_n(
                b'-',
                gaps_after.iter().filter(|&&g| g == k).count(),
            ));
        }
        row
    }

    /// The score of column `k`: a different one for each of ten, so that a residue given the wrong column shows.
    fn scores_for(width: usize, none_at: &[usize]) -> Vec<Option<f32>> {
        (0..width)
            .map(|k| (!none_at.contains(&k)).then(|| (k % 10) as f32 / 10.0))
            .collect()
    }

    #[test]
    fn a_residue_gets_the_score_of_its_own_column_through_the_gaps_of_the_row() {
        let mut f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        assert!(f.residues.iter().all(|r| r.conservation.is_none()));
        // two gaps after residue 9 and one after 40 (counting from 0), and a column with no score
        let row = ubiquitin_row(&[9, 9, 40]);
        let rows = vec![("ubiquitin".to_string(), row)];
        let scores = scores_for(79, &[12]);
        let matched = f.apply_conservation(&rows, &scores);
        assert_eq!(matched.len(), 1, "{matched:?}");
        assert_eq!(
            (matched[0].chain, matched[0].row.as_str()),
            ('A', "ubiquitin")
        );
        assert!((matched[0].identity - 1.0).abs() < 1e-9);
        // the column of residue j: j, and 2 more after residue 9, and 1 more after residue 40 (the row has the gaps there)
        let column = |j: usize| j + 2 * usize::from(j > 9) + usize::from(j > 40);
        for (j, r) in f.residues.iter().enumerate() {
            // the residue knows its column, the one with no score too
            assert_eq!(r.column, Some(column(j)), "residue {j}");
            let want = scores[column(j)].map(|s| s / 0.9);
            match (r.conservation, want) {
                (Some(got), Some(want)) => assert!(
                    (got - want).abs() < 1e-5,
                    "residue {j}: {got} against {want}"
                ),
                (None, None) => {}
                other => panic!("residue {j}: {other:?}"),
            }
        }
        // the lowest score of the alignment is 0 and the highest 1
        let given: Vec<f32> = f.residues.iter().filter_map(|r| r.conservation).collect();
        assert!(given.iter().cloned().fold(f32::MAX, f32::min) < 1e-6);
        assert!((given.iter().cloned().fold(f32::MIN, f32::max) - 1.0).abs() < 1e-6);
        assert_eq!(matched[0].coloured, given.len());
        assert_eq!(matched[0].residues, 76);
        assert!(
            matched[0].coloured < 76,
            "the residue in the column with no score gets none"
        );
        // and a column is asked for its residues: one here, none for a gap of the row or a column past the end
        assert_eq!(f.in_column(column(12)), [12]);
        assert_eq!(f.in_column(column(41)), [41]);
        assert!(f.in_column(10).is_empty(), "a gap of the row");
        assert!(f.in_column(79).is_empty());
    }

    #[test]
    fn the_copies_of_a_chain_are_all_in_the_column() {
        // ubiquitin twice: chain A, and a copy of it called chain B (the chain is the twenty-second character of a line)
        let atoms: Vec<&str> = UBIQUITIN
            .lines()
            .filter(|l| l.starts_with("ATOM"))
            .collect();
        let copy = atoms.iter().map(|l| format!("{}B{}", &l[..21], &l[22..]));
        let text = atoms
            .iter()
            .map(|l| l.to_string())
            .chain(copy)
            .collect::<Vec<_>>()
            .join("\n");
        let mut f = Fold::new(parse(&text).expect("a structure"));
        assert_eq!(f.residues.len(), 152);
        // a row with one gap after residue 9
        let rows = [("ubiquitin".to_string(), ubiquitin_row(&[9]))];
        let matched = f.apply_conservation(&rows, &scores_for(77, &[]));
        assert_eq!(matched.len(), 2, "{matched:?}");
        // residue 20 is in column 21, in both chains, and residue 5 in column 5; the gap holds neither
        assert_eq!(f.in_column(21), [20, 96]);
        assert_eq!(f.in_column(5), [5, 81]);
        assert!(f.in_column(10).is_empty());
    }

    #[test]
    fn a_chain_that_no_row_matches_is_in_no_column_and_another_alignment_takes_the_old_columns_away(
    ) {
        // the antibody with its antigen: the alignment holds the antigen only (chain C)
        let mut f = Fold::new(parse(ANTIBODY).expect("a structure"));
        let antigen: Vec<u8> = f
            .residues
            .iter()
            .filter(|r| r.chain == 'C')
            .map(|r| r.letter)
            .collect();
        let width = antigen.len();
        let matched = f.apply_conservation(
            &[("lysozyme".to_string(), antigen)],
            &scores_for(width, &[]),
        );
        assert_eq!(matched.len(), 1, "{matched:?}");
        let first_of_c = f
            .residues
            .iter()
            .position(|r| r.chain == 'C')
            .expect("chain C");
        assert_eq!(f.in_column(0), [first_of_c]);
        assert!(f
            .residues
            .iter()
            .all(|r| (r.chain == 'C') == r.column.is_some()));
        // another alignment takes the old columns away before it gives its own
        let other = f.apply_conservation(&[("x".to_string(), b"WWWWWWWWWW".to_vec())], &[None; 10]);
        assert!(other.is_empty());
        assert!(f.residues.iter().all(|r| r.column.is_none()));
    }

    #[test]
    fn a_close_relative_is_matched_and_the_row_that_fits_best_is_the_one_used() {
        let mut f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        // a relative with three residues changed, and an unrelated row that comes first
        let mut relative = ubiquitin_row(&[]);
        for k in [5, 30, 60] {
            relative[k] = b'W';
        }
        let rows = vec![
            ("unrelated".to_string(), vec![b'G'; 76]),
            ("relative".to_string(), relative),
            ("the same".to_string(), ubiquitin_row(&[])),
        ];
        let matched = f.apply_conservation(&rows, &scores_for(76, &[]));
        assert_eq!(matched.len(), 1);
        assert_eq!(matched[0].row, "the same", "{matched:?}");
        // the relative alone is matched, a little less than fully
        let mut g = Fold::new(parse(UBIQUITIN).expect("a structure"));
        let alone = g.apply_conservation(&rows[1..2], &scores_for(76, &[]));
        assert_eq!(alone[0].row, "relative");
        assert!(
            alone[0].identity < 1.0 && alone[0].identity > 0.9,
            "{}",
            alone[0].identity
        );
    }

    #[test]
    fn an_alignment_of_another_protein_colours_nothing_and_clears_what_was_there() {
        let mut f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        f.apply_conservation(
            &[("ubiquitin".to_string(), ubiquitin_row(&[]))],
            &scores_for(76, &[]),
        );
        assert!(f.residues.iter().any(|r| r.conservation.is_some()));
        let other = vec![(
            "other".to_string(),
            b"MSDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEY".to_vec(),
        )];
        assert!(f
            .apply_conservation(&other, &scores_for(50, &[]))
            .is_empty());
        assert!(f.residues.iter().all(|r| r.conservation.is_none()));
        assert!(f.apply_conservation(&[], &[]).is_empty());
    }

    /// 400 rows that each have their own pattern, so that none shares a run of four residues with ubiquitin, and then ubiquitin.
    fn hundreds_of_rows_and_the_match_last() -> Vec<(String, Vec<u8>)> {
        let letters = b"ACDEFGHIKLMNPQRSTVWY";
        let mut rows: Vec<(String, Vec<u8>)> = (0..400)
            .map(|k| {
                let row: Vec<u8> = (0..76)
                    .map(|j| letters[(k * 7 + j * (k % 13 + 2)) % 20])
                    .collect();
                (format!("other {k}"), row)
            })
            .collect();
        rows.push(("ubiquitin".to_string(), ubiquitin_row(&[])));
        rows
    }

    #[test]
    fn only_a_few_likely_rows_are_lined_up_with_a_chain_and_the_right_one_is_among_them() {
        let rows = hundreds_of_rows_and_the_match_last();
        let chain = Fold::new(parse(UBIQUITIN).expect("a structure")).chain_sequences()[0]
            .1
            .clone();
        let plain: Vec<Vec<u8>> = rows.iter().map(|(_, r)| r.clone()).collect();
        let candidates = candidate_rows(&plain, &chain);
        assert!(candidates.len() <= MOST_CANDIDATES, "{}", candidates.len());
        assert!(candidates.contains(&400), "{candidates:?}");
        // and so the whole is found, and coloured
        let mut f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        let matched = f.apply_conservation(&rows, &scores_for(76, &[]));
        assert_eq!(matched.len(), 1, "{matched:?}");
        assert_eq!(
            (matched[0].row.as_str(), matched[0].coloured),
            ("ubiquitin", 76)
        );
    }

    #[test]
    fn an_alignment_in_lower_case_or_with_dots_for_gaps_is_matched_like_any_other() {
        // twenty distractors besides the match: more than the rows tried when nothing shares a run, so the match must be
        // found by the runs it shares, in whatever case it is written
        let mut rows = hundreds_of_rows_and_the_match_last();
        rows.drain(20..400);
        let (name, row) = rows.pop().expect("the match");
        let lower: Vec<u8> = row.to_ascii_lowercase();
        rows.push((name, lower));
        let mut f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        let matched = f.apply_conservation(&rows, &scores_for(76, &[]));
        assert_eq!(matched.len(), 1, "{matched:?}");
        assert_eq!(matched[0].coloured, 76);
        // dots and dashes are both gaps, and a residue gets the score of its column either way
        let with_dashes = ubiquitin_row(&[9, 9, 40]);
        let with_dots: Vec<u8> = with_dashes
            .iter()
            .map(|&c| if c == b'-' { b'.' } else { c })
            .collect();
        let scores = scores_for(79, &[]);
        let (mut a, mut b) = (
            Fold::new(parse(UBIQUITIN).expect("a structure")),
            Fold::new(parse(UBIQUITIN).expect("a structure")),
        );
        a.apply_conservation(&[("x".to_string(), with_dashes)], &scores);
        b.apply_conservation(&[("x".to_string(), with_dots)], &scores);
        assert!(a
            .residues
            .iter()
            .zip(&b.residues)
            .all(|(x, y)| x.conservation == y.conservation));
        assert!(a.residues.iter().all(|r| r.conservation.is_some()));
    }

    #[test]
    fn the_row_that_covers_the_chain_beats_a_short_piece_that_is_more_alike() {
        let mut f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        // the first 40 residues exactly (and gaps for the rest), and the whole chain with two residues changed
        let exact = ubiquitin_row(&[]);
        let mut piece = exact.clone();
        piece[40..].fill(b'-');
        let mut whole = exact;
        whole[10] = b'W';
        whole[60] = b'W';
        let rows = vec![("piece".to_string(), piece), ("whole".to_string(), whole)];
        let matched = f.apply_conservation(&rows, &scores_for(76, &[]));
        assert_eq!(matched[0].row, "whole", "{matched:?}");
        assert_eq!(matched[0].coloured, 76);
    }

    /// Ubiquitin twice: chains A and B, the same protein.
    fn dimer() -> Fold {
        let atoms: Vec<&str> = UBIQUITIN
            .lines()
            .filter(|l| l.starts_with("ATOM"))
            .collect();
        let mut text = String::new();
        for line in &atoms {
            text.push_str(line);
            text.push('\n');
        }
        text.push_str("TER\n");
        for line in &atoms {
            let mut bytes = line.as_bytes().to_vec();
            bytes[21] = b'B';
            text.push_str(&String::from_utf8(bytes).expect("text"));
            text.push('\n');
        }
        Fold::new(parse(&text).expect("a structure"))
    }

    #[test]
    fn the_chains_of_a_homodimer_are_lined_up_with_the_same_row_and_coloured_alike() {
        let mut f = dimer();
        assert_eq!(f.chain_sequences().len(), 2);
        let rows = vec![("ubiquitin".to_string(), ubiquitin_row(&[]))];
        let matched = f.apply_conservation(&rows, &scores_for(76, &[]));
        assert_eq!(
            matched.iter().map(|m| m.chain).collect::<Vec<_>>(),
            vec!['A', 'B']
        );
        assert!(matched
            .iter()
            .all(|m| m.row == "ubiquitin" && m.coloured == 76));
        let (a, b): (Vec<_>, Vec<_>) = f.residues.iter().partition(|r| r.chain == 'A');
        assert_eq!(a.len(), 76);
        assert!(a
            .iter()
            .zip(&b)
            .all(|(x, y)| x.conservation == y.conservation));
    }

    #[test]
    fn the_conservation_can_be_taken_off() {
        let mut f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        f.apply_conservation(
            &[("x".to_string(), ubiquitin_row(&[]))],
            &scores_for(76, &[]),
        );
        assert!(f.residues.iter().any(|r| r.conservation.is_some()));
        assert!(f.residues.iter().all(|r| r.column.is_some()));
        f.clear_conservation();
        assert!(f.residues.iter().all(|r| r.conservation.is_none()));
        // the columns go with the scores: the two were made together, and one without the other would say too much
        assert!(f.residues.iter().all(|r| r.column.is_none()));
    }

    #[test]
    fn a_short_chain_with_no_run_of_four_in_common_is_still_matched_to_a_row_of_its_length() {
        // too short for the quick look at shared runs to find anything: the rows of about its length are all tried
        let mut f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        let first: Vec<u8> = f.residues.iter().take(30).map(|r| r.letter).collect();
        // change every fourth residue: no run of four is left whole, but 73 % of the residues are the same
        f.residues.truncate(30);
        let mut row = first.clone();
        for k in (0..30).step_by(4) {
            row[k] = b'W';
        }
        let matched = f.apply_conservation(&[("changed".to_string(), row)], &[Some(0.5); 30]);
        assert_eq!(matched.len(), 1, "{matched:?}");
    }

    #[test]
    fn an_alignment_whose_columns_all_score_alike_is_coloured_in_the_middle() {
        let mut f = Fold::new(parse(UBIQUITIN).expect("a structure"));
        f.apply_conservation(
            &[("ubiquitin".to_string(), ubiquitin_row(&[]))],
            &[Some(0.8); 76],
        );
        assert!(f.residues.iter().all(|r| r.conservation == Some(0.5)));
    }
}
