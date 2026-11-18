//! The interface between two groups of chains of a structure: how much surface each residue buries when the groups meet,
//! which residues touch, and which pairs of atoms could be hydrogen bonds. Everything is measured on the coordinates that
//! are given; nothing is predicted, so for a complex that was itself predicted the numbers are only as good as the model.
//!
//! The surface is the solvent-accessible one of Lee and Richards, found with the points method of Shrake and Rupley (1973): a
//! probe of 1.4 angstroms rolled over atoms of the radii of Bondi (1964), hydrogens left out. An atom's area is the share of
//! points on its sphere (the atom's radius and the probe's) that no other atom's sphere holds. What a group buries is its area
//! alone less its area in the complex.

use super::par::par_map;
use super::pdb::Structure;
use std::collections::HashMap;
use std::f64::consts::PI;

/// The radius of the water probe, in angstroms.
pub const PROBE: f64 = 1.4;
/// Points on each atom's sphere: enough for the area of an atom to be right to a fraction of a percent.
pub const POINTS: usize = 960;
/// Residues are in contact when two of their atoms (hydrogens left out) are this close, in angstroms.
pub const CONTACT: f64 = 4.5;
/// A donor and an acceptor are a candidate hydrogen bond when their heavy atoms are this close, in angstroms.
pub const HBOND: f64 = 3.5;
/// A residue is at the interface when it buries at least this much, in square angstroms.
pub const MIN_BURIED: f64 = 1.0;

/// The van der Waals radius of an element (Bondi 1964 and Mantina et al. 2009 for the halogens and selenium), in angstroms;
/// `None` for hydrogen, which is left out. **An element not listed has the radius of carbon**, which is wrong for a metal ion
/// (sodium, potassium and calcium are larger, zinc smaller): the surface of a structure with ions is only approximate there.
pub fn radius(element: &str) -> Option<f64> {
    Some(match element.trim().to_ascii_uppercase().as_str() {
        "H" | "D" => return None,
        "N" => 1.55,
        "O" => 1.52,
        "S" => 1.80,
        "P" => 1.80,
        "SE" => 1.90,
        "F" => 1.47,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        _ => 1.70,
    })
}

/// `n` points spread evenly over the unit sphere, along a spiral that turns by the golden angle at each step.
pub fn sphere_points(n: usize) -> Vec<[f64; 3]> {
    let golden = PI * (3.0 - 5f64.sqrt());
    (0..n)
        .map(|k| {
            let z = 1.0 - (2.0 * k as f64 + 1.0) / n as f64;
            let r = (1.0 - z * z).max(0.0).sqrt();
            let phi = golden * k as f64;
            [r * phi.cos(), r * phi.sin(), z]
        })
        .collect()
}

fn dist2(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
}

/// The cell of a grid with cells `size` wide that holds `p`.
fn cell(p: &[f64; 3], size: f64) -> (i32, i32, i32) {
    let c = |x: f64| (x / size).floor() as i32;
    (c(p[0]), c(p[1]), c(p[2]))
}

fn grid(points: &[[f64; 3]], size: f64) -> HashMap<(i32, i32, i32), Vec<usize>> {
    let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
    for (i, p) in points.iter().enumerate() {
        grid.entry(cell(p, size)).or_default().push(i);
    }
    grid
}

/// The indices of the points of `grid` that lie in the cells next to the one of `p` (the 27 cells about it).
fn near<'a>(
    grid: &'a HashMap<(i32, i32, i32), Vec<usize>>,
    p: &[f64; 3],
    size: f64,
) -> impl Iterator<Item = usize> + 'a {
    let (cx, cy, cz) = cell(p, size);
    (-1..=1).flat_map(move |dx| {
        (-1..=1).flat_map(move |dy| {
            (-1..=1).flat_map(move |dz| {
                grid.get(&(cx + dx, cy + dy, cz + dz))
                    .into_iter()
                    .flatten()
                    .copied()
            })
        })
    })
}

/// The solvent-accessible area of each atom, in square angstroms: `coords` and `radii` are the atoms' (the probe is added to
/// each radius here), `points` is how many points are put on each sphere.
pub fn atom_areas(coords: &[[f64; 3]], radii: &[f64], probe: f64, points: usize) -> Vec<f64> {
    let spots = sphere_points(points.max(1));
    let reach: Vec<f64> = radii.iter().map(|r| r + probe).collect();
    let widest = reach.iter().copied().fold(0.0, f64::max);
    if coords.is_empty() {
        return Vec::new();
    }
    // atoms whose spheres can meet are less than two of the widest apart, so the cells are that wide
    let size = (2.0 * widest).max(1e-6);
    let cells = grid(coords, size);
    let all: Vec<usize> = (0..coords.len()).collect();
    par_map(&all, |&i| {
        let neighbours: Vec<usize> = near(&cells, &coords[i], size)
            .filter(|&j| j != i && dist2(&coords[i], &coords[j]) < (reach[i] + reach[j]).powi(2))
            .collect();
        let mut free = 0usize;
        // the neighbour that held the last point is tried first for the next: points beside each other are held by the same one
        let mut last = 0usize;
        'points: for u in &spots {
            let p = [
                coords[i][0] + reach[i] * u[0],
                coords[i][1] + reach[i] * u[1],
                coords[i][2] + reach[i] * u[2],
            ];
            if !neighbours.is_empty() {
                let held = |j: usize| dist2(&p, &coords[j]) < reach[j] * reach[j];
                if held(neighbours[last]) {
                    continue 'points;
                }
                for (k, &j) in neighbours.iter().enumerate() {
                    if held(j) {
                        last = k;
                        continue 'points;
                    }
                }
            }
            free += 1;
        }
        4.0 * PI * reach[i] * reach[i] * free as f64 / spots.len() as f64
    })
}

// --- hydrogen bonds -----------------------------------------------------------------------------------------------------------

/// Whether an atom of a residue can give a hydrogen to a bond (a donor) and can take one (an acceptor), by names. Only what the
/// atoms' names say: no hydrogens are looked at, and the angle is not measured. In a ligand, a sugar or a modified group any
/// nitrogen or oxygen (by its element, so a sodium ion named NA is not one) may be either.
fn bond_roles(residue: &str, atom: &str, element: &str, hetero: bool) -> (bool, bool) {
    if hetero {
        let polar = matches!(element.trim().to_ascii_uppercase().as_str(), "N" | "O");
        return (polar, polar);
    }
    match (residue, atom) {
        (_, "N") => (residue != "PRO", false),
        (_, "O" | "OXT") => (false, true),
        ("ARG", "NE" | "NH1" | "NH2") => (true, false),
        ("ASN", "OD1") | ("ASP", "OD1" | "OD2") | ("GLN", "OE1") | ("GLU", "OE1" | "OE2") => {
            (false, true)
        }
        ("ASN", "ND2") | ("GLN", "NE2") | ("LYS", "NZ") | ("TRP", "NE1") => (true, false),
        ("HIS", "ND1" | "NE2") | ("SER", "OG") | ("THR", "OG1") | ("TYR", "OH") => (true, true),
        _ => (false, false),
    }
}

// --- the interface --------------------------------------------------------------------------------------------------------------

/// The surface of one residue, alone in its group and in the complex, in square angstroms.
#[derive(Clone, Debug, PartialEq)]
pub struct Buried {
    /// Index into `Structure::residues`.
    pub residue: usize,
    pub alone: f64,
    pub bound: f64,
}

impl Buried {
    /// What the residue buries when its group meets the other.
    pub fn buried(&self) -> f64 {
        (self.alone - self.bound).max(0.0)
    }
}

/// Two residues, one of each group, that have atoms within `CONTACT` of each other.
#[derive(Clone, Debug, PartialEq)]
pub struct Contact {
    pub a: usize,
    pub b: usize,
    /// How many pairs of atoms are that close.
    pub atom_pairs: usize,
    pub closest: f64,
}

/// A donor and an acceptor, one in each group, within `HBOND`: a candidate, not a bond that was seen.
#[derive(Clone, Debug, PartialEq)]
pub struct HBond {
    /// Indices into `Structure::atoms`.
    pub a_atom: usize,
    pub b_atom: usize,
    pub distance: f64,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Interface {
    pub a_chains: Vec<char>,
    pub b_chains: Vec<char>,
    /// Every residue of the first group, with its surface alone and in the complex.
    pub a: Vec<Buried>,
    pub b: Vec<Buried>,
    pub contacts: Vec<Contact>,
    /// Closest first.
    pub hbonds: Vec<HBond>,
    /// The solvent-accessible area of each group alone, and of the two together.
    pub area_a: f64,
    pub area_b: f64,
    pub area_bound: f64,
}

impl Interface {
    /// What the two groups bury together: the area they have alone less the area they have as a complex. It is the sum of what
    /// every residue buries.
    pub fn buried_total(&self) -> f64 {
        self.area_a + self.area_b - self.area_bound
    }

    /// The area of the interface as it is usually quoted: half of what is buried, one half of it on each side.
    pub fn interface_area(&self) -> f64 {
        self.buried_total() / 2.0
    }

    /// The residues of a group that bury at least `MIN_BURIED`, the most buried first.
    pub fn interface_residues<'a>(&self, side: &'a [Buried]) -> Vec<&'a Buried> {
        let mut at: Vec<&Buried> = side.iter().filter(|r| r.buried() >= MIN_BURIED).collect();
        at.sort_by(|x, y| y.buried().total_cmp(&x.buried()));
        at
    }
}

/// One residue of the interface, for a table.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    /// Which group it is in: 1 for the first, 2 for the second.
    pub side: u8,
    /// Index into `Structure::residues`.
    pub residue: usize,
    pub chain: char,
    /// The three-letter name.
    pub name: String,
    pub number: i32,
    /// The insertion code that goes with the number (52A in antibody numbering), a space when there is none.
    pub insertion: char,
    /// What it buries, in square angstroms.
    pub buried: f64,
    /// With how many residues of the other group it is in contact.
    pub contacts: usize,
    /// How many candidate hydrogen bonds its atoms are in.
    pub hbonds: usize,
}

impl Row {
    /// The number as it is written, with its insertion code: `52A`.
    pub fn number_text(&self) -> String {
        if self.insertion.is_whitespace() {
            self.number.to_string()
        } else {
            format!("{}{}", self.number, self.insertion)
        }
    }
}

impl Interface {
    /// The residues at the interface, first group then second, the most buried first in each.
    pub fn rows(&self, s: &Structure) -> Vec<Row> {
        let in_contact = |residue: usize| {
            self.contacts
                .iter()
                .filter(|c| c.a == residue || c.b == residue)
                .count()
        };
        let in_hbonds = |residue: usize| {
            self.hbonds
                .iter()
                .filter(|h| {
                    s.atoms[h.a_atom].residue == residue || s.atoms[h.b_atom].residue == residue
                })
                .count()
        };
        let mut rows = Vec::new();
        for (side, group) in [(1u8, &self.a), (2u8, &self.b)] {
            for b in self.interface_residues(group) {
                let r = &s.residues[b.residue];
                rows.push(Row {
                    side,
                    residue: b.residue,
                    chain: r.chain,
                    name: r.name.clone(),
                    number: r.number,
                    insertion: r.insertion,
                    buried: b.buried(),
                    contacts: in_contact(b.residue),
                    hbonds: in_hbonds(b.residue),
                });
            }
        }
        rows
    }
}

/// The rows as text for a spreadsheet: tab separated, with a head line.
pub fn tsv(rows: &[Row]) -> String {
    let mut out =
        String::from("side\tchain\tresidue\tnumber\tburied_A2\tcontacts\thbond_candidates\n");
    for r in rows {
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\n",
            r.side,
            r.chain,
            r.name,
            r.number_text(),
            r.buried,
            r.contacts,
            r.hbonds
        ));
    }
    out
}

/// The interface of the chains `a` against the chains `b` of `s`. Water and hydrogens are left out; so are the chains that are in
/// neither group, which is what it is to look at two bodies.
pub fn analyse(s: &Structure, a: &[char], b: &[char]) -> Result<Interface, String> {
    if a.is_empty() || b.is_empty() {
        return Err("Choose at least one chain for each side of the interface.".into());
    }
    if let Some(c) = a.iter().find(|c| b.contains(c)) {
        return Err(format!(
            "Chain {c} is on both sides: a chain is on one side or the other."
        ));
    }
    let known: Vec<char> = s.residues.iter().map(|r| r.chain).collect();
    if let Some(c) = a.iter().chain(b).find(|c| !known.contains(c)) {
        return Err(format!("The structure has no chain {c}."));
    }
    // an atom with no coordinates (a file may say nan or inf) would be counted as free of every other, or as in contact with
    // everything: it is refused, naming it
    if let Some(r) = s.residues.iter().find(|r| {
        (a.contains(&r.chain) || b.contains(&r.chain))
            && r.atoms
                .iter()
                .any(|&ai| s.atoms[ai].pos.iter().any(|x| !x.is_finite()))
    }) {
        return Err(format!(
            "An atom of {} {}{} (chain {}) has no usable coordinates, so the surface cannot be measured.",
            r.name, r.number, r.insertion, r.chain
        ));
    }
    // the atoms of the two groups: where they are, how big, and whose
    struct Taken {
        atom: usize,
        residue: usize,
        side_a: bool,
        pos: [f64; 3],
        radius: f64,
    }
    let mut taken: Vec<Taken> = Vec::new();
    for (ri, r) in s.residues.iter().enumerate() {
        let side_a = a.contains(&r.chain);
        if !side_a && !b.contains(&r.chain) {
            continue;
        }
        for &ai in &r.atoms {
            let atom = &s.atoms[ai];
            if let Some(radius) = radius(&atom.element) {
                taken.push(Taken {
                    atom: ai,
                    residue: ri,
                    side_a,
                    pos: atom.pos.map(f64::from),
                    radius,
                });
            }
        }
    }
    let areas = |keep: &dyn Fn(&Taken) -> bool| -> Vec<(usize, f64)> {
        let (index, coords, radii): (Vec<usize>, Vec<[f64; 3]>, Vec<f64>) = taken
            .iter()
            .enumerate()
            .filter(|(_, t)| keep(t))
            .map(|(k, t)| (k, t.pos, t.radius))
            .collect();
        index
            .into_iter()
            .zip(atom_areas(&coords, &radii, PROBE, POINTS))
            .collect()
    };
    let bound = areas(&|_| true);
    let alone_a = areas(&|t| t.side_a);
    let alone_b = areas(&|t| !t.side_a);
    let mut per_residue: HashMap<usize, (f64, f64)> = HashMap::new();
    for (k, area) in &bound {
        per_residue.entry(taken[*k].residue).or_default().1 += area;
    }
    for (k, area) in alone_a.iter().chain(&alone_b) {
        per_residue.entry(taken[*k].residue).or_default().0 += area;
    }
    let side = |first: bool| -> Vec<Buried> {
        let mut rows: Vec<Buried> = per_residue
            .iter()
            .filter(|(r, _)| a.contains(&s.residues[**r].chain) == first)
            .map(|(&residue, &(alone, bound))| Buried {
                residue,
                alone,
                bound,
            })
            .collect();
        rows.sort_by_key(|r| r.residue);
        rows
    };
    let total = |v: &[(usize, f64)]| v.iter().map(|x| x.1).sum::<f64>();

    // residues in contact, and the candidate hydrogen bonds: each atom of the first group against the atoms of the second near it
    let side_b: Vec<usize> = (0..taken.len()).filter(|&k| !taken[k].side_a).collect();
    let b_positions: Vec<[f64; 3]> = side_b.iter().map(|&k| taken[k].pos).collect();
    let cells = grid(&b_positions, CONTACT);
    let mut pairs: HashMap<(usize, usize), (usize, f64)> = HashMap::new();
    let mut hbonds: Vec<HBond> = Vec::new();
    for ta in taken.iter().filter(|t| t.side_a) {
        let ra = &s.residues[ta.residue];
        let a_name = s.atoms[ta.atom].name.as_str();
        let (a_donor, a_acceptor) =
            bond_roles(&ra.name, a_name, &s.atoms[ta.atom].element, ra.hetero);
        for j in near(&cells, &ta.pos, CONTACT) {
            let tb = &taken[side_b[j]];
            let d = dist2(&ta.pos, &tb.pos).sqrt();
            // no coordinate is NaN here: `analyse` refused those above
            if d > CONTACT {
                continue;
            }
            let entry = pairs
                .entry((ta.residue, tb.residue))
                .or_insert((0, f64::MAX));
            entry.0 += 1;
            entry.1 = entry.1.min(d);
            if d <= HBOND {
                let rb = &s.residues[tb.residue];
                let (b_donor, b_acceptor) = bond_roles(
                    &rb.name,
                    &s.atoms[tb.atom].name,
                    &s.atoms[tb.atom].element,
                    rb.hetero,
                );
                if (a_donor && b_acceptor) || (a_acceptor && b_donor) {
                    hbonds.push(HBond {
                        a_atom: ta.atom,
                        b_atom: tb.atom,
                        distance: d,
                    });
                }
            }
        }
    }
    let mut contacts: Vec<Contact> = pairs
        .into_iter()
        .map(|((a, b), (atom_pairs, closest))| Contact {
            a,
            b,
            atom_pairs,
            closest,
        })
        .collect();
    contacts.sort_by_key(|c| (c.a, c.b));
    hbonds.sort_by(|x, y| {
        x.distance
            .total_cmp(&y.distance)
            .then(x.a_atom.cmp(&y.a_atom))
    });
    Ok(Interface {
        a_chains: a.to_vec(),
        b_chains: b.to_vec(),
        a: side(true),
        b: side(false),
        contacts,
        hbonds,
        area_a: total(&alone_a),
        area_b: total(&alone_b),
        area_bound: total(&bound),
    })
}

#[cfg(test)]
mod tests;
