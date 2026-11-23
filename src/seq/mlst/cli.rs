//! `helixloom mlst`: type genomes against schemes that the person has downloaded.

use super::call::{type_contigs, Params, Typing};
use super::output;
use super::scheme::Scheme;
use super::search::Contigs;
use crate::pipeline::steps::read_fasta;
use std::path::{Path, PathBuf};

pub const USAGE: &str = "\
Usage: helixloom mlst (--scheme FOLDER | --schemes FOLDER) [options] CONTIGS...

  Multi-locus sequence typing in the manner of mlst (Torsten Seemann), written from its README (docs/MLST.md). A scheme is a folder that you
  downloaded (PubMLST's data is not included): NAME.txt with the profile table, and LOCUS.tfa with the alleles of each locus. CONTIGS are
  FASTA files of the genome (not compressed), one genome to a file. One line is written for each:
  the file, the scheme, the sequence type, and the call for each locus: n exact, ~n novel (closest n), n? partial, - missing, n,m several.

  --scheme FOLDER    type against this scheme
  --schemes FOLDER   a folder of scheme folders: the one that scores best for each genome is used (the first in alphabetical order when
                     two score alike)
  --minid N          how alike (percent) a stretch must be to an allele to be novel or partial (default 95)
  --mincov N         how much (percent) of an allele a partial match must cover (default 10)
  --minscore N       the score (1 to 100) a scheme needs when schemes are tried (default 50)
  --full             the table with a header: FILE, SCHEME, ST, STATUS, SCORE, ALLELES
  --csv              commas in place of tabs
  --json FILE        also write JSON
  --novel FILE       also write the novel alleles as FASTA, named scheme.locus-MD5
  --nopath           write the name of the file without its folder
  --label NAME       the name to write for the file (one genome only)
  -o FILE            where to write the lines (the standard output when left out)
  -h, --help         this text
";

struct Args {
    schemes: Vec<PathBuf>,
    scan: Option<PathBuf>,
    params: Params,
    full: bool,
    csv: bool,
    json: Option<String>,
    novel: Option<String>,
    nopath: bool,
    label: Option<String>,
    output: Option<String>,
    files: Vec<String>,
}

fn number(flag: &str, value: &str) -> Result<f64, String> {
    value
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("{flag} needs a number, not {value}."))
}

fn parse(args: &[String]) -> Result<Args, String> {
    let mut a = Args {
        schemes: Vec::new(),
        scan: None,
        params: Params::default(),
        full: false,
        csv: false,
        json: None,
        novel: None,
        nopath: false,
        label: None,
        output: None,
        files: Vec::new(),
    };
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value."))
        };
        match flag.as_str() {
            "--scheme" => a.schemes.push(PathBuf::from(value("--scheme")?)),
            "--schemes" => a.scan = Some(PathBuf::from(value("--schemes")?)),
            "--minid" => a.params.min_identity = number("--minid", &value("--minid")?)?,
            "--mincov" => a.params.min_coverage = number("--mincov", &value("--mincov")?)?,
            "--minscore" => a.params.min_score = number("--minscore", &value("--minscore")?)?,
            "--full" => a.full = true,
            "--csv" => a.csv = true,
            "--json" => a.json = Some(value("--json")?),
            "--novel" => a.novel = Some(value("--novel")?),
            "--nopath" => a.nopath = true,
            "--label" => a.label = Some(value("--label")?),
            "-o" => a.output = Some(value("-o")?),
            other if other.starts_with('-') && other.len() > 1 => {
                return Err(format!(
                    "{other} is not an option (`helixloom mlst --help`)."
                ))
            }
            file => a.files.push(file.to_string()),
        }
    }
    if a.schemes.is_empty() && a.scan.is_none() {
        return Err(
            "Name a scheme: --scheme FOLDER, or --schemes FOLDER for a folder of them.".into(),
        );
    }
    if a.files.is_empty() {
        return Err("Name the FASTA file of a genome to type.".into());
    }
    if a.label.is_some() && a.files.len() > 1 {
        return Err("--label names one genome: give one file.".into());
    }
    for (flag, v) in [
        ("--minid", a.params.min_identity),
        ("--mincov", a.params.min_coverage),
        ("--minscore", a.params.min_score),
    ] {
        if !(0.0..=100.0).contains(&v) {
            return Err(format!("{flag} is a number from 0 to 100, not {v}."));
        }
    }
    Ok(a)
}

/// The scheme folders in a folder of them, in alphabetical order.
fn scan(folder: &Path) -> Result<Vec<PathBuf>, String> {
    let mut found: Vec<PathBuf> = std::fs::read_dir(folder)
        .map_err(|e| format!("Could not read {}: {e}", folder.display()))?
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.is_dir())
        .collect();
    found.sort();
    if found.is_empty() {
        return Err(format!("{} holds no scheme folder.", folder.display()));
    }
    Ok(found)
}

/// The best typing of a genome among the schemes: the highest score, the first when two are alike, none under `min_score`.
pub fn best(schemes: &[Scheme], contigs: &Contigs, p: &Params) -> Option<Typing> {
    let mut best: Option<Typing> = None;
    for s in schemes {
        let t = type_contigs(s, contigs, p);
        if t.score >= p.min_score && best.as_ref().is_none_or(|b| t.score > b.score) {
            best = Some(t);
        }
    }
    best
}

pub fn run(args: &[String]) -> i32 {
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{USAGE}");
        return 0;
    }
    match run_inner(args) {
        Ok(()) => 0,
        Err(why) => {
            eprintln!("{why}");
            1
        }
    }
}

fn run_inner(args: &[String]) -> Result<(), String> {
    let a = parse(args)?;
    let mut folders = a.schemes.clone();
    if let Some(dir) = &a.scan {
        folders.extend(scan(dir)?);
    }
    let schemes: Vec<Scheme> = folders
        .iter()
        .map(|f| Scheme::load(f))
        .collect::<Result<_, _>>()?;
    let separator = if a.csv { ',' } else { '\t' };
    let mut text = String::new();
    if a.full {
        text += &output::FULL_HEADER.join(&separator.to_string());
        text.push('\n');
    }
    let (mut rows, mut novel) = (Vec::new(), String::new());
    for path in &a.files {
        let contigs =
            std::fs::read_to_string(path).map_err(|e| format!("Could not read {path}: {e}"))?;
        let sequences = read_fasta(&contigs);
        if sequences.is_empty() {
            return Err(format!(
                "{path} holds no sequences: it should be a FASTA file."
            ));
        }
        let typing = best(&schemes, &Contigs::new(sequences)?, &a.params);
        let shown = match (&a.label, a.nopath) {
            (Some(label), _) => label.clone(),
            (None, true) => Path::new(path)
                .file_name()
                .map_or(path.clone(), |n| n.to_string_lossy().into_owned()),
            (None, false) => path.clone(),
        };
        text += &if a.full {
            output::full_line(&shown, typing.as_ref(), separator)
        } else {
            output::line(&shown, typing.as_ref(), separator)
        };
        text.push('\n');
        if let Some(t) = &typing {
            novel += &output::novel_fasta(&shown, t);
        }
        rows.push((shown.clone(), path.clone(), typing));
    }
    if let Some(path) = &a.json {
        std::fs::write(path, output::json(&rows))
            .map_err(|e| format!("Could not write {path}: {e}"))?;
    }
    if let Some(path) = &a.novel {
        std::fs::write(path, novel).map_err(|e| format!("Could not write {path}: {e}"))?;
    }
    match &a.output {
        Some(path) => {
            std::fs::write(path, text).map_err(|e| format!("Could not write {path}: {e}"))
        }
        None => {
            print!("{text}");
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests;
