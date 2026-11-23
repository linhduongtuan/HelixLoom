//! A typing scheme as a folder holds it, in the layout that mlst's README gives for a scheme that is added by hand: a tab-separated
//! profile table named `NAME.txt` (a header, then a row for each sequence type: the type, one column of allele numbers for each locus, and
//! optionally columns such as `clonal_complex`), and one nucleotide FASTA file for each locus, `LOCUS.tfa`, whose records are `>LOCUS_ID` or
//! `>LOCUS-ID`. The scheme is named for its folder. The folder is one that the person downloaded (PubMLST's data may not be bundled).

use crate::pipeline::steps::read_fasta;
use std::path::Path;

/// Columns of a profile table that are about the type and not about a locus (mlst's README lists them as ignored).
const NOT_LOCI: [&str; 6] = [
    "mlst_clade",
    "clonal_complex",
    "species",
    "cc",
    "lineage",
    "st",
];

/// Extensions that an allele file may have.
const EXTENSIONS: [&str; 4] = ["tfa", "fa", "fasta", "fna"];

#[derive(Clone, Debug, PartialEq)]
pub struct Allele {
    /// What follows the locus name and its `_` or `-`: a number, or in some schemes letters.
    pub id: String,
    /// Upper case.
    pub seq: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Locus {
    pub name: String,
    pub alleles: Vec<Allele>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Scheme {
    pub name: String,
    pub loci: Vec<Locus>,
    /// A sequence type and the id of its allele at each locus, in the order of `loci`.
    pub profiles: Vec<(String, Vec<String>)>,
}

/// The records of one locus's FASTA text: `LOCUS_ID` and `LOCUS-ID` both give the id, and a header that does not begin with the locus's
/// name is refused with its place.
fn alleles_of(locus: &str, text: &str) -> Result<Vec<Allele>, String> {
    let mut out = Vec::new();
    for (n, (header, bases)) in read_fasta(text).into_iter().enumerate() {
        let name = header.split_whitespace().next().unwrap_or("");
        let id = name
            .strip_prefix(locus)
            .and_then(|rest| rest.strip_prefix(['_', '-']))
            .filter(|id| !id.is_empty())
            .ok_or_else(|| {
                format!("Record {} of {locus}'s alleles is named {name:?}: it should be {locus}_ID or {locus}-ID.", n + 1)
            })?;
        let seq: Vec<u8> = bases
            .iter()
            .filter(|b| !b.is_ascii_whitespace())
            .map(u8::to_ascii_uppercase)
            .collect();
        if seq.is_empty() {
            return Err(format!("The allele {name} has no sequence."));
        }
        out.push(Allele {
            id: id.to_string(),
            seq,
        });
    }
    if out.is_empty() {
        return Err(format!("There is no allele of {locus}."));
    }
    Ok(out)
}

impl Scheme {
    /// A scheme from its texts: the profile table, and for each locus the text of its alleles' FASTA file (looked up by the locus's name).
    pub fn from_texts(
        name: &str,
        profile_table: &str,
        fastas: &dyn Fn(&str) -> Option<String>,
    ) -> Result<Scheme, String> {
        let mut lines = profile_table
            .lines()
            .map(|l| l.trim_end_matches('\r'))
            .filter(|l| !l.trim().is_empty());
        let header: Vec<&str> = lines
            .next()
            .ok_or("The profile table is empty.")?
            .split('\t')
            .collect();
        if !header.first().is_some_and(|h| h.eq_ignore_ascii_case("ST")) {
            return Err(format!(
                "The first column of the profile table should be ST, not {:?}.",
                header.first().unwrap_or(&"")
            ));
        }
        // the columns that are loci, each with its alleles
        let mut columns = Vec::new();
        let mut loci = Vec::new();
        for (c, &h) in header.iter().enumerate().skip(1) {
            if NOT_LOCI.contains(&h.to_ascii_lowercase().as_str()) {
                continue;
            }
            let text = fastas(h).ok_or_else(|| format!("The profile table has a column {h}, but there is no file of alleles for it ({h}.tfa)."))?;
            loci.push(Locus {
                name: h.to_string(),
                alleles: alleles_of(h, &text)?,
            });
            columns.push(c);
        }
        if loci.is_empty() {
            return Err("The profile table has no locus column.".into());
        }
        let mut profiles = Vec::new();
        for (n, line) in lines.enumerate() {
            let cells: Vec<&str> = line.split('\t').collect();
            let id = cells[0].trim();
            let alleles: Option<Vec<String>> = columns
                .iter()
                .map(|&c| cells.get(c).map(|a| a.trim().to_string()))
                .collect();
            let Some(a) = alleles else {
                return Err(format!(
                    "Row {} of the profile table has too few columns.",
                    n + 2
                ));
            };
            if id.is_empty() {
                return Err(format!(
                    "Row {} of the profile table has no sequence type in its first column.",
                    n + 2
                ));
            }
            profiles.push((id.to_string(), a));
        }
        Ok(Scheme {
            name: name.to_string(),
            loci,
            profiles,
        })
    }

    /// Read the scheme in `folder`: its name is the folder's name, its profile table `NAME.txt` (or the one `.txt` file there is).
    pub fn load(folder: &Path) -> Result<Scheme, String> {
        let name = folder
            .file_name()
            .and_then(|n| n.to_str())
            .ok_or_else(|| format!("{} is not a folder with a name.", folder.display()))?
            .to_string();
        let own = folder.join(format!("{name}.txt"));
        let table = if own.is_file() {
            own
        } else {
            let mut txts: Vec<_> = std::fs::read_dir(folder)
                .map_err(|e| format!("Could not read {}: {e}", folder.display()))?
                .flatten()
                .map(|e| e.path())
                .filter(|p| p.extension().is_some_and(|x| x == "txt"))
                .collect();
            txts.sort();
            match txts.as_slice() {
                [one] => one.clone(),
                [] => {
                    return Err(format!(
                        "{} has no profile table: expected {name}.txt.",
                        folder.display()
                    ))
                }
                _ => {
                    return Err(format!(
                        "{} has several .txt files and none is {name}.txt.",
                        folder.display()
                    ))
                }
            }
        };
        let text = std::fs::read_to_string(&table)
            .map_err(|e| format!("Could not read {}: {e}", table.display()))?;
        let fastas = |locus: &str| {
            EXTENSIONS
                .iter()
                .map(|x| folder.join(format!("{locus}.{x}")))
                .find(|p| p.is_file())
                .and_then(|p| std::fs::read_to_string(p).ok())
        };
        Scheme::from_texts(&name, &text, &fastas)
    }

    /// The sequence type that has exactly these alleles (one id for each locus, in the order of `loci`), if there is one.
    pub fn sequence_type(&self, alleles: &[&str]) -> Option<&str> {
        self.profiles
            .iter()
            .find(|(_, a)| a.iter().map(String::as_str).eq(alleles.iter().copied()))
            .map(|(st, _)| st.as_str())
    }
}

#[cfg(test)]
mod tests;
