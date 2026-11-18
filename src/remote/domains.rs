//! The domain-search job: protein sequences go in, the domains HMMER finds in them come out (`seq/hmmer`, `pipeline/domains`). The
//! program (`hmmscan`) and the database (Pfam) are found on the computer that runs the job, never named by the request. The
//! answer carries the names and lengths of the sequences, so it can be shown without them.

use super::tasks::{bit, flag, lines, number};
use super::wire::{meta, Container};
use super::{clean_line, one_line};
use crate::models::{self, PFAM};
use crate::pipeline::domains::{scan, Cutoff};
use crate::pipeline::tools::{self, Installed, HMMSCAN};
use crate::seq::hmmer::{self, DomainHit, MAX_LENGTH, MAX_RESIDUES, MAX_SEQUENCES};
use std::collections::HashMap;
use std::path::Path;
use web_time::Instant;

/// What to run: the protein sequences (name and residues) and how a hit is judged significant.
#[derive(Clone, Debug, PartialEq)]
pub struct DomainsJob {
    pub queries: Vec<(String, Vec<u8>)>,
    pub cutoff: Cutoff,
}

impl DomainsJob {
    pub fn encode(&self) -> Container {
        let (curated, evalue) = match self.cutoff {
            Cutoff::Curated => (true, 0.0),
            Cutoff::Evalue(x) => (false, x),
        };
        Container::new()
            .with(
                "meta",
                lines(&[("curated", bit(curated)), ("evalue", evalue.to_string())]),
            )
            .with(
                "names",
                self.queries
                    .iter()
                    .map(|(n, _)| format!("{}\n", one_line(n)))
                    .collect::<String>(),
            )
            .with(
                "sequences",
                self.queries
                    .iter()
                    .map(|(_, s)| format!("{}\n", String::from_utf8_lossy(s)))
                    .collect::<String>(),
            )
    }

    /// A request comes from the network: the sequences are counted and measured, and the cutoff checked, before anything is run.
    pub fn decode(c: &Container, m: &HashMap<&str, &str>) -> Result<DomainsJob, String> {
        let cutoff = if flag(m, "curated") {
            Cutoff::Curated
        } else {
            Cutoff::Evalue(number(m, "evalue")?)
        };
        cutoff.check()?;
        let names: Vec<&str> = c.require("names")?.lines().collect();
        let residues: Vec<&str> = c.require("sequences")?.lines().collect();
        if names.len() != residues.len() {
            return Err("The message has a different number of names and sequences.".into());
        }
        if residues.is_empty() || residues.len() > MAX_SEQUENCES {
            return Err(format!(
                "A job takes from 1 to {MAX_SEQUENCES} sequences; this one has {}.",
                residues.len()
            ));
        }
        let total: usize = residues.iter().map(|s| s.len()).sum();
        if total > MAX_RESIDUES || residues.iter().any(|s| s.len() > MAX_LENGTH) {
            return Err(format!(
                "A job takes at most {MAX_RESIDUES} residues in all and {MAX_LENGTH} in a sequence."
            ));
        }
        if residues
            .iter()
            .any(|s| s.is_empty() || !s.bytes().all(|b| b.is_ascii_uppercase()))
        {
            return Err(
                "A sequence in the message is empty or has something that is not a capital letter."
                    .into(),
            );
        }
        Ok(DomainsJob {
            queries: names
                .iter()
                .zip(&residues)
                .map(|(n, s)| (one_line(n), s.as_bytes().to_vec()))
                .collect(),
            cutoff,
        })
    }
}

/// What a search found: for each sequence its name and length, and the domains.
#[derive(Clone, Debug, PartialEq)]
pub struct DomainsResult {
    pub names: Vec<String>,
    pub lengths: Vec<usize>,
    pub hits: Vec<DomainHit>,
    /// What made the answer: the program with its version ("HMMER 3.4") and the database file.
    pub program: String,
    pub database: String,
    pub seconds: f64,
}

/// Run the job here: with the `hmmscan` found on this computer and Pfam from the models folder, or why not.
pub fn run(job: &DomainsJob) -> Result<DomainsResult, String> {
    let program = tools::detect(&HMMSCAN).ok_or(
        "HMMER's hmmscan was not found on this computer: install HMMER 3 (for example `brew install hmmer`, \
         `apt install hmmer` or `conda install -c bioconda hmmer`) so that `hmmscan` is on the PATH.",
    )?;
    let db = models::installed(&PFAM)
        .then(|| models::path_of(&PFAM, "Pfam-A.hmm"))
        .flatten()
        .ok_or_else(|| {
            "The Pfam database is not installed here: download Pfam-A.hmm.gz from \
             https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/, unpack it, run `hmmpress Pfam-A.hmm`, and put the \
             five files (Pfam-A.hmm and the four it makes) in the pfam folder of the models folder (docs/MODELS.md)."
                .to_string()
        })?;
    run_with(job, &program, &db)
}

/// The same with the program and the database named.
pub fn run_with(job: &DomainsJob, program: &Installed, db: &Path) -> Result<DomainsResult, String> {
    let started = Instant::now();
    let hits = scan(program, db, &job.queries, job.cutoff)?;
    Ok(DomainsResult {
        names: job.queries.iter().map(|(n, _)| one_line(n)).collect(),
        lengths: job.queries.iter().map(|(_, s)| s.len()).collect(),
        hits,
        program: program.describe("HMMER"),
        database: db
            .file_name()
            .map_or_else(String::new, |n| n.to_string_lossy().into_owned()),
        seconds: started.elapsed().as_secs_f64(),
    })
}

impl DomainsResult {
    pub fn encode(&self) -> Container {
        let names: String = self
            .names
            .iter()
            .zip(&self.lengths)
            .map(|(n, l)| format!("{}\t{l}\n", clean_line(n)))
            .collect();
        let domains: String = self
            .hits
            .iter()
            .map(|h| {
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    h.sequence,
                    h.family,
                    h.accession,
                    h.family_length,
                    h.score,
                    h.i_evalue,
                    h.c_evalue,
                    h.hmm_from,
                    h.hmm_to,
                    h.ali_from,
                    h.ali_to,
                    h.env_from,
                    h.env_to,
                    h.accuracy,
                    clean_line(&h.description)
                )
            })
            .collect();
        Container::new()
            .with(
                "meta",
                lines(&[
                    ("seconds", self.seconds.to_string()),
                    ("program", clean_line(&self.program)),
                    ("database", clean_line(&self.database)),
                ]),
            )
            .with("names", names)
            .with("domains", domains)
    }

    /// A result comes from a server that is not trusted to keep to what it was asked: the sequences, their lengths and every
    /// coordinate are checked, so that nothing drawn from it can run off the end of a sequence or a profile.
    pub fn decode(c: &Container) -> Result<DomainsResult, String> {
        let m = meta(c.require("meta")?);
        let bad = |why: String| format!("The result's domains are damaged: {why}");
        let (mut names, mut lengths) = (Vec::new(), Vec::new());
        for row in c.require("names")?.lines() {
            let (name, length) = row
                .rsplit_once('\t')
                .ok_or_else(|| bad("a sequence has no length.".into()))?;
            let length: usize = length
                .parse()
                .ok()
                .filter(|l| (1..=MAX_LENGTH).contains(l))
                .ok_or_else(|| bad("a sequence has a length that is out of range.".into()))?;
            names.push(name.to_string());
            lengths.push(length);
            if names.len() > MAX_SEQUENCES {
                return Err(bad("there are too many sequences.".into()));
            }
        }
        let mut hits = Vec::new();
        for (k, row) in c.require("domains")?.lines().enumerate() {
            let f: Vec<&str> = row.split('\t').collect();
            let [sequence, family, accession, flen, score, ie, ce, hf, ht, af, at, ef, et, acc, description] =
                f[..]
            else {
                return Err(bad(format!(
                    "line {} has the wrong number of fields.",
                    k + 1
                )));
            };
            let number = |what: &str, x: &str| -> Result<f64, String> {
                x.parse::<f64>()
                    .map_err(|_| bad(format!("line {}: {what} is no number.", k + 1)))
            };
            let count = |what: &str, x: &str| -> Result<usize, String> {
                x.parse::<usize>()
                    .map_err(|_| bad(format!("line {}: {what} is no count.", k + 1)))
            };
            let hit = DomainHit {
                sequence: count("the sequence", sequence)?,
                family: family.to_string(),
                accession: accession.to_string(),
                family_length: count("the profile length", flen)?,
                score: number("the score", score)?,
                i_evalue: number("the E-value", ie)?,
                c_evalue: number("the E-value", ce)?,
                hmm_from: count("a coordinate", hf)?,
                hmm_to: count("a coordinate", ht)?,
                ali_from: count("a coordinate", af)?,
                ali_to: count("a coordinate", at)?,
                env_from: count("a coordinate", ef)?,
                env_to: count("a coordinate", et)?,
                accuracy: number("the probability", acc)?,
                description: description.to_string(),
            };
            let length = *lengths.get(hit.sequence).ok_or_else(|| {
                bad(format!(
                    "line {} is about a sequence that is not there.",
                    k + 1
                ))
            })?;
            hmmer::check(&hit, length, k + 1)?;
            hits.push(hit);
        }
        Ok(DomainsResult {
            names,
            lengths,
            hits,
            program: m.get("program").copied().unwrap_or_default().to_string(),
            database: m.get("database").copied().unwrap_or_default().to_string(),
            seconds: number(&m, "seconds")?,
        })
    }
}

#[cfg(test)]
mod real;
#[cfg(test)]
mod tests;
