//! Searching protein sequences for domains by running HMMER 3's `hmmscan` against a database of profile HMMs (Pfam), and reading the
//! table it writes (`seq::hmmer`). The program and the database are other people's: they are found on the computer that runs the
//! search, never named by a request (`docs/DOMAINS.md`).

use super::tools::{self, Installed};
use crate::seq::hmmer::{parse_domtblout, DomainHit, MAX_LENGTH};
use std::path::Path;

/// How a hit is judged significant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Cutoff {
    /// The curated gathering thresholds of each family (`--cut_ga`): what Pfam itself uses to say a sequence belongs to a family.
    Curated,
    /// One E-value for the sequence and the domain alike.
    Evalue(f64),
}

impl Cutoff {
    /// The most permissive E-value that is taken.
    pub const MAX_EVALUE: f64 = 100.0;

    pub fn check(self) -> Result<(), String> {
        match self {
            Cutoff::Evalue(x) if !(x > 0.0 && x <= Self::MAX_EVALUE) => Err(format!(
                "The E-value cutoff must be above 0 and at most {}, not {x}.",
                Self::MAX_EVALUE
            )),
            _ => Ok(()),
        }
    }

    fn args(self) -> Vec<String> {
        match self {
            Cutoff::Curated => vec!["--cut_ga".into()],
            Cutoff::Evalue(x) => ["-E", &x.to_string(), "--domE", &x.to_string()]
                .map(String::from)
                .to_vec(),
        }
    }
}

/// The files `hmmpress` makes beside a database, which `hmmscan` cannot do without.
pub const PRESSED: [&str; 4] = [".h3f", ".h3i", ".h3m", ".h3p"];

/// Whether `db` has been pressed, or what to do about it.
pub fn check_pressed(db: &Path) -> Result<(), String> {
    let missing: Vec<String> = PRESSED
        .iter()
        .map(|suffix| format!("{}{suffix}", db.display()))
        .filter(|p| !Path::new(p).is_file())
        .collect();
    if missing.is_empty() {
        return Ok(());
    }
    Err(format!(
        "The database {} has not been pressed (missing {}): run `hmmpress {}` once.",
        db.display(),
        missing.join(", "),
        db.display()
    ))
}

/// The arguments of the `hmmscan` that writes the domain table `table` for the sequences in `fasta` against `db`. The main report,
/// which is long and not read, goes to `report`.
pub fn scan_args(
    cutoff: Cutoff,
    db: &Path,
    fasta: &Path,
    table: &Path,
    report: &Path,
) -> Result<Vec<String>, String> {
    let text = |p: &Path| {
        p.to_str()
            .map(String::from)
            .ok_or_else(|| format!("The path {} is not text.", p.display()))
    };
    let mut args: Vec<String> = ["--noali", "-o"].map(String::from).to_vec();
    args.push(text(report)?);
    args.extend(["--domtblout".to_string(), text(table)?]);
    args.extend(cutoff.args());
    args.extend([text(db)?, text(fasta)?]);
    Ok(args)
}

/// Search `queries` (name and residues, in capitals) with `hmmscan` against `db`. The sequences go to the program under names of
/// their own (`q1`, `q2`, …), so that nothing about a name can matter to it, and the answer is matched back by those.
pub fn scan(
    hmmscan: &Installed,
    db: &Path,
    queries: &[(String, Vec<u8>)],
    cutoff: Cutoff,
) -> Result<Vec<DomainHit>, String> {
    cutoff.check()?;
    if queries.is_empty() {
        return Err("There is no sequence to search.".into());
    }
    if let Some((name, s)) = queries.iter().find(|(_, s)| {
        s.is_empty() || s.len() > MAX_LENGTH || !s.iter().all(u8::is_ascii_alphabetic)
    }) {
        return Err(format!("The sequence {name} cannot be searched: it needs letters only, from 1 to {MAX_LENGTH}, not {}.", s.len()));
    }
    check_pressed(db)?;
    let dir = tools::WorkFolder::new("hmmscan")?;
    scan_in(&dir, hmmscan, db, queries, cutoff)
}

fn scan_in(
    dir: &Path,
    hmmscan: &Installed,
    db: &Path,
    queries: &[(String, Vec<u8>)],
    cutoff: Cutoff,
) -> Result<Vec<DomainHit>, String> {
    let (fasta, table, report) = (
        dir.join("queries.fa"),
        dir.join("domains.tbl"),
        dir.join("report.txt"),
    );
    let text: String = queries
        .iter()
        .enumerate()
        .map(|(k, (_, s))| format!(">q{}\n{}\n", k + 1, String::from_utf8_lossy(s)))
        .collect();
    std::fs::write(&fasta, text).map_err(|e| format!("Could not write the sequences: {e}"))?;
    let args = scan_args(cutoff, db, &fasta, &table, &report)?;
    let args: Vec<&str> = args.iter().map(String::as_str).collect();
    tools::run(&hmmscan.path, &args, dir)?;
    let written = std::fs::read_to_string(&table)
        .map_err(|_| "hmmscan finished but wrote no table of domains.".to_string())?;
    let names: Vec<String> = (1..=queries.len()).map(|k| format!("q{k}")).collect();
    let known: Vec<(&str, usize)> = names
        .iter()
        .map(String::as_str)
        .zip(queries.iter().map(|(_, s)| s.len()))
        .collect();
    parse_domtblout(&written, &known)
}

#[cfg(test)]
pub(crate) mod testing;
#[cfg(test)]
mod tests;
