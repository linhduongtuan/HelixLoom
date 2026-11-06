//! Reading protein structures in PDB format: the atoms with their coordinates, grouped into residues and chains,
//! with the helices and strands the file lists, if it lists any. Only the first model is read, and the first
//! conformer of an atom that has several. Water is left out.
//!
//! The fixed columns follow the wwPDB format description (version 3.3).

#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    pub name: String,
    pub element: String,
    pub pos: [f32; 3],
    /// The temperature factor column; predictors such as ESMFold and AlphaFold put their confidence (pLDDT) here.
    pub bfactor: f32,
    /// Index of the residue the atom belongs to.
    pub residue: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Residue {
    /// The three-letter name, upper case.
    pub name: String,
    pub chain: char,
    pub number: i32,
    pub insertion: char,
    /// Indices into `Structure::atoms`.
    pub atoms: Vec<usize>,
    /// Part of a ligand or a modified group rather than of the protein chain.
    pub hetero: bool,
}

/// A run of residues the file calls a helix or a strand: the chain, and the first and last residue numbers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Span {
    pub chain: char,
    pub from: i32,
    pub to: i32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Structure {
    pub title: String,
    pub atoms: Vec<Atom>,
    pub residues: Vec<Residue>,
    pub helices: Vec<Span>,
    pub strands: Vec<Span>,
    /// The file says how the structure was determined (`EXPDTA`): it is an experiment, not a prediction.
    pub experimental: bool,
}

/// The one-letter code of a residue name; modified residues that stand in for a standard one read as it.
pub fn letter(name: &str) -> u8 {
    match name {
        "ALA" => b'A',
        "ARG" => b'R',
        "ASN" => b'N',
        "ASP" => b'D',
        "CYS" => b'C',
        "GLN" => b'Q',
        "GLU" => b'E',
        "GLY" => b'G',
        "HIS" => b'H',
        "ILE" => b'I',
        "LEU" => b'L',
        "LYS" => b'K',
        "MET" | "MSE" => b'M',
        "PHE" => b'F',
        "PRO" => b'P',
        "SER" => b'S',
        "THR" => b'T',
        "TRP" => b'W',
        "TYR" => b'Y',
        "VAL" => b'V',
        _ => b'X',
    }
}

fn column(line: &str, from: usize, to: usize) -> &str {
    // columns are 1-based and inclusive; a short line gives what it has
    line.get(from - 1..to.min(line.len())).unwrap_or("")
}

fn number<T: std::str::FromStr>(line: &str, from: usize, to: usize) -> Option<T> {
    column(line, from, to).trim().parse().ok()
}

fn chain_of(line: &str, at: usize) -> char {
    line.as_bytes().get(at - 1).map_or(' ', |&c| char::from(c))
}

pub fn parse(text: &str) -> Result<Structure, String> {
    let mut s = Structure::default();
    let mut current: Option<(char, i32, char, String)> = None;
    for line in text.lines() {
        let record = column(line, 1, 6).trim_end();
        match record {
            "ENDMDL" => break,
            "TITLE" => {
                let part = column(line, 11, 80).trim();
                if !s.title.is_empty() {
                    s.title.push(' ');
                }
                s.title.push_str(part);
            }
            "EXPDTA" => s.experimental = true,
            "HELIX" => {
                if let (Some(from), Some(to)) = (number(line, 22, 25), number(line, 34, 37)) {
                    s.helices.push(Span {
                        chain: chain_of(line, 20),
                        from,
                        to,
                    });
                }
            }
            "SHEET" => {
                if let (Some(from), Some(to)) = (number(line, 23, 26), number(line, 34, 37)) {
                    s.strands.push(Span {
                        chain: chain_of(line, 22),
                        from,
                        to,
                    });
                }
            }
            "ATOM" | "HETATM" => {
                let hetero = record == "HETATM";
                let alt = chain_of(line, 17);
                let residue_name = column(line, 18, 20).trim().to_ascii_uppercase();
                if !matches!(alt, ' ' | 'A' | '1') || residue_name == "HOH" || residue_name == "WAT"
                {
                    continue;
                }
                let (Some(x), Some(y), Some(z)) = (
                    number::<f32>(line, 31, 38),
                    number::<f32>(line, 39, 46),
                    number::<f32>(line, 47, 54),
                ) else {
                    continue;
                };
                let key = (
                    chain_of(line, 22),
                    number(line, 23, 26).unwrap_or(0),
                    chain_of(line, 27),
                    residue_name.clone(),
                );
                if current.as_ref() != Some(&key) {
                    s.residues.push(Residue {
                        name: residue_name,
                        chain: key.0,
                        number: key.1,
                        insertion: key.2,
                        atoms: Vec::new(),
                        // a selenomethionine is part of the chain even though the file calls it HETATM
                        hetero: hetero && letter(&key.3) == b'X',
                    });
                    current = Some(key);
                }
                let name = column(line, 13, 16).trim().to_string();
                let element = {
                    let given = column(line, 77, 78).trim();
                    if given.is_empty() {
                        name.chars()
                            .find(char::is_ascii_alphabetic)
                            .map(String::from)
                            .unwrap_or_default()
                    } else {
                        given.to_string()
                    }
                };
                let residue = s.residues.len() - 1;
                s.residues[residue].atoms.push(s.atoms.len());
                s.atoms.push(Atom {
                    name,
                    element,
                    pos: [x, y, z],
                    bfactor: number(line, 61, 66).unwrap_or(0.0),
                    residue,
                });
            }
            _ => {}
        }
    }
    if s.atoms.is_empty() {
        return Err(
            "No atoms found: this does not look like a PDB file (it needs ATOM lines).".into(),
        );
    }
    Ok(s)
}

/// The backbone atoms of a residue, as the secondary-structure rules need them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Backbone {
    pub n: [f32; 3],
    pub ca: [f32; 3],
    pub c: [f32; 3],
    pub o: [f32; 3],
    pub proline: bool,
}

impl Structure {
    fn atom_of(&self, residue: &Residue, name: &str) -> Option<&Atom> {
        residue
            .atoms
            .iter()
            .map(|&a| &self.atoms[a])
            .find(|a| a.name == name)
    }

    /// The residues of the protein chains (no ligands) that have a C-alpha, with their index.
    pub fn amino_acids(&self) -> Vec<usize> {
        self.residues
            .iter()
            .enumerate()
            .filter(|(_, r)| !r.hetero && self.atom_of(r, "CA").is_some())
            .map(|(i, _)| i)
            .collect()
    }

    /// The chain letters in the order they first appear.
    pub fn chains(&self) -> Vec<char> {
        let mut seen = Vec::new();
        for &i in &self.amino_acids() {
            let c = self.residues[i].chain;
            if !seen.contains(&c) {
                seen.push(c);
            }
        }
        seen
    }

    /// The sequence of a chain, one letter per residue.
    #[cfg(test)]
    pub fn sequence(&self, chain: char) -> Vec<u8> {
        self.amino_acids()
            .into_iter()
            .map(|i| &self.residues[i])
            .filter(|r| r.chain == chain)
            .map(|r| letter(&r.name))
            .collect()
    }

    /// N, CA, C and O of each residue of the protein chains that has all four, in file order, and whether each
    /// one follows the one before it in the chain (their C and N are bonded).
    pub fn backbone(&self) -> Vec<(usize, Backbone, bool)> {
        let mut out: Vec<(usize, Backbone, bool)> = Vec::new();
        for i in self.amino_acids() {
            let r = &self.residues[i];
            let get = |n: &str| self.atom_of(r, n).map(|a| a.pos);
            let (Some(n), Some(ca), Some(c), Some(o)) = (get("N"), get("CA"), get("C"), get("O"))
            else {
                continue;
            };
            let linked = out
                .last()
                .is_some_and(|(_, prev, _)| distance(prev.c, n) < 2.0);
            out.push((
                i,
                Backbone {
                    n,
                    ca,
                    c,
                    o,
                    proline: r.name == "PRO",
                },
                linked,
            ));
        }
        out
    }

    /// The middle and the radius of a sphere that holds every atom (a little larger than needed).
    pub fn centre_and_radius(&self) -> ([f32; 3], f32) {
        let n = self.atoms.len().max(1) as f32;
        let mut centre = [0.0f32; 3];
        for a in &self.atoms {
            for (c, p) in centre.iter_mut().zip(a.pos) {
                *c += p / n;
            }
        }
        let radius = self
            .atoms
            .iter()
            .map(|a| distance(a.pos, centre))
            .fold(1.0, f32::max);
        (centre, radius)
    }

    /// How sure a predictor was about a residue, 0 to 100, read from the temperature factor of its C-alpha. Files
    /// that give the confidence as a fraction (0 to 1) are scaled up.
    pub fn confidence(&self, residue: usize) -> Option<f32> {
        let b = self.atom_of(self.residues.get(residue)?, "CA")?.bfactor;
        let fraction = self
            .residues
            .iter()
            .filter_map(|r| self.atom_of(r, "CA"))
            .all(|a| (0.0..=1.0).contains(&a.bfactor));
        Some(if fraction { b * 100.0 } else { b })
    }

    /// Whether the temperature factors look like predictor confidences rather than an experiment's.
    pub fn looks_predicted(&self) -> bool {
        !self.experimental
            && self
                .residues
                .iter()
                .filter_map(|r| self.atom_of(r, "CA"))
                .all(|a| (0.0..=100.0).contains(&a.bfactor))
    }
}

pub fn distance(a: [f32; 3], b: [f32; 3]) -> f32 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

#[cfg(test)]
mod tests;
