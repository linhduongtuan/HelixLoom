//! `flye` and `medaka`: the two programs that stay outside (`docs/ROUTES.md`). They are called by name and found like MAFFT and IQ-TREE are; a
//! step that asks for one that is not installed fails and says so. **Neither was run in building this**: neither is installed on the computer
//! it was written on, so the command lines follow what each program's documentation gives, and the glue is tried on stand-ins that write the
//! files the programs write.

use super::{tag, AsmStep};
use crate::pipeline::flow::Value;
use crate::pipeline::steps::{fasta_text, read_fasta, Named};
use crate::pipeline::tools::{self, Program};
use crate::seq::fastq::{to_fastq, FastqRead};
use std::path::Path;

fn program(step: &AsmStep) -> Option<&'static Program> {
    match step {
        AsmStep::Flye { .. } => Some(&tools::FLYE),
        AsmStep::Medaka { .. } => Some(&tools::MEDAKA),
        _ => None,
    }
}

/// The installed program's version line (for the key of an answer), or nothing for a step that calls none, or that finds none.
pub fn version(step: &AsmStep) -> String {
    program(step)
        .and_then(tools::detect)
        .map_or(String::new(), |i| i.version)
}

pub fn ready(step: &AsmStep) -> Result<(), String> {
    match program(step) {
        Some(p) if tools::detect(p).is_none() => Err(format!(
            "{} is not installed here; the step asks for it by name, so nothing else is used instead (`helixloom tools` says how to install it)",
            p.label
        )),
        _ => Ok(()),
    }
}

fn threads() -> String {
    std::thread::available_parallelism()
        .map_or(1, |n| n.get())
        .to_string()
}

fn write(dir: &Path, file: &str, text: &str) -> Result<(), String> {
    std::fs::write(dir.join(file), text).map_err(|e| format!("Could not write {file}: {e}"))
}

/// The contigs of a Flye run in `out`: `assembly.fasta`, each said to be `circular=true` when `assembly_info.txt` says so (its fourth
/// column, `Y` or `N`).
fn read_flye(out: &Path) -> Result<Vec<Named>, String> {
    let text = std::fs::read_to_string(out.join("assembly.fasta"))
        .map_err(|e| format!("Flye wrote no assembly.fasta: {e}"))?;
    let info = std::fs::read_to_string(out.join("assembly_info.txt")).unwrap_or_default();
    let circular: Vec<&str> = info
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter_map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            (f.get(3) == Some(&"Y")).then(|| f[0])
        })
        .collect();
    let rows = read_fasta(&text);
    if rows.is_empty() {
        return Err("Flye made no contig".into());
    }
    Ok(rows
        .into_iter()
        .map(|(name, seq)| {
            let id = name.split_whitespace().next().unwrap_or("").to_string();
            let name = if circular.contains(&id.as_str()) {
                tag(&name, "circular=true")
            } else {
                name
            };
            (name, seq)
        })
        .collect())
}

/// Flye on `reads`, with the program at `exe`: `flye --nano-raw reads.fastq --out-dir flye-out --threads N [--genome-size S]`.
pub fn flye_with(
    exe: &Path,
    mode: &str,
    genome_size: &str,
    reads: &[FastqRead],
    dir: &Path,
) -> Result<Vec<Named>, String> {
    write(dir, "reads.fastq", &to_fastq(reads))?;
    let flag = format!("--{mode}");
    let n = threads();
    let mut args = vec![
        flag.as_str(),
        "reads.fastq",
        "--out-dir",
        "flye-out",
        "--threads",
        n.as_str(),
    ];
    if genome_size != "auto" {
        args.extend(["--genome-size", genome_size]);
    }
    tools::run(exe, &args, dir)?;
    read_flye(&dir.join("flye-out"))
}

pub fn flye(step: &AsmStep, reads: &[FastqRead]) -> Result<(Value, String), String> {
    let AsmStep::Flye { mode, genome_size } = step else {
        return Err("not a flye step".into());
    };
    ready(step)?;
    let installed = tools::detect(&tools::FLYE).ok_or("Flye is not installed here")?;
    let dir = tools::WorkFolder::new("flow")?;
    let rows = flye_with(&installed.path, mode, genome_size, reads, &dir)?;
    let circles = rows
        .iter()
        .filter(|r| super::has_tag(&r.0, "circular=true"))
        .count();
    let note = format!(
        "{} contigs, {circles} circular ({})",
        rows.len(),
        installed.describe("Flye")
    );
    Ok((Value::Sequences(rows), note))
}

/// Medaka with the program at `exe`: the draft written as `draft.fasta` with plain names, and
/// `medaka_consensus -i READS -d draft.fasta -o medaka-out -t N`; the polished contigs (`consensus.fasta`) take their old names back.
pub fn medaka_with(
    exe: &Path,
    draft: &[Named],
    reads: &Path,
    dir: &Path,
) -> Result<Vec<Named>, String> {
    let plain: Vec<Named> = draft
        .iter()
        .enumerate()
        .map(|(i, r)| (format!("c{i}"), r.1.clone()))
        .collect();
    write(dir, "draft.fasta", &fasta_text(&plain))?;
    let reads = std::fs::canonicalize(reads)
        .map_err(|e| format!("Could not find the reads {}: {e}", reads.display()))?;
    let reads = reads.to_string_lossy();
    let n = threads();
    tools::run(
        exe,
        &[
            "-i",
            &reads,
            "-d",
            "draft.fasta",
            "-o",
            "medaka-out",
            "-t",
            &n,
        ],
        dir,
    )?;
    let text = std::fs::read_to_string(dir.join("medaka-out").join("consensus.fasta"))
        .map_err(|e| format!("Medaka wrote no consensus.fasta: {e}"))?;
    let polished = read_fasta(&text);
    if polished.len() != draft.len() {
        return Err(format!(
            "Medaka gave {} contigs for the {} it was given",
            polished.len(),
            draft.len()
        ));
    }
    // by name where it kept ours, else in order
    Ok(draft
        .iter()
        .enumerate()
        .map(|(i, (name, _))| {
            let at = polished
                .iter()
                .position(|p| p.0.split_whitespace().next() == Some(format!("c{i}").as_str()))
                .unwrap_or(i);
            (name.clone(), polished[at].1.clone())
        })
        .collect())
}

pub fn medaka(
    step: &AsmStep,
    draft: &[Named],
    reads: Option<&Path>,
) -> Result<(Value, String), String> {
    ready(step)?;
    let installed = tools::detect(&tools::MEDAKA).ok_or("Medaka is not installed here")?;
    let reads = reads.ok_or("medaka has no reads to polish with")?;
    let dir = tools::WorkFolder::new("flow")?;
    let rows = medaka_with(&installed.path, draft, reads, &dir)?;
    let note = format!(
        "{} contigs polished ({})",
        rows.len(),
        installed.describe("Medaka")
    );
    Ok((Value::Sequences(rows), note))
}

#[cfg(test)]
mod tests;
