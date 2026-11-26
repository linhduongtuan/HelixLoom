//! The steps of the assembly routes, as steps of a pipeline file (`docs/ROUTES.md`): cleaning short reads, selecting long reads, assembling short
//! reads into a graph, joining the graph with long reads, the checks that a long-read assembly of a plasmid wants (a concatemer, a circle,
//! where it starts), and the two programs that stay outside (Flye, which makes the long-read assembly, and Medaka, which polishes it).
//!
//! | step | takes | gives | what does it |
//! |------|-------|-------|--------------|
//! | `trim` | reads | reads | `seq::fastq::clean` (fastp's rules) |
//! | `filter` | reads | reads | `seq::readfilter` (the multiplex-MinION paper's selection) |
//! | `assemble` | reads | gfa | `helixloom dbg` (`seq::debruijn`) |
//! | `hybrid` | gfa | sequences | `helixloom hybrid` (`seq::hybrid`), with the long reads of `long=` |
//! | `flye` | reads | sequences | the installed Flye |
//! | `concatemer` | sequences | sequences | `seq::circular::find_concatemer` |
//! | `circularize` | sequences | sequences | `seq::circular::end_overlap` |
//! | `rotate` | sequences | sequences | `seq::hybrid::rotate`, with the genes of `genes=` |
//! | `medaka` | sequences | sequences | the installed Medaka, with the reads of `reads=` |
//!
//! A step that needs more than what flows to it (the long reads of `hybrid`, the genes of `rotate`, the reads of `medaka`) names a file or a
//! folder in a setting: a folder holds, for each input file, a file of the same name before the dot. The names of contigs carry what was
//! found as words after the name: `circular=true`, `concatemer=3`, `rotated=true`.

mod contigs;
mod external;
mod graph;
mod reads;
pub mod routes;

use super::flow::{digest, Kind, Value};
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// Where the quality cut of `trim` is made.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cut {
    /// From the 3' end, until a window reaches the quality.
    Tail,
    /// From the first window under the quality to the end.
    Right,
    Off,
}

impl Cut {
    fn word(self) -> &'static str {
        match self {
            Cut::Tail => "tail",
            Cut::Right => "right",
            Cut::Off => "off",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AsmStep {
    Trim {
        /// The adapter to cut from the 3' end; `none` for none.
        adapter: String,
        quality: u8,
        min_length: usize,
        cut: Cut,
    },
    Filter {
        min_length: usize,
        target_bases: u64,
        window: usize,
    },
    Assemble {
        /// `auto` or a list such as `21,33,55`.
        ks: String,
        /// `auto`, `off` or a number.
        cutoff: String,
        min_qual: u8,
    },
    Hybrid {
        long: String,
        mode: String,
        /// `auto` or a number from 0 to 100.
        min_quality: String,
        /// `auto` or a number.
        single_depth: String,
        overlap: usize,
    },
    Flye {
        mode: String,
        /// `auto` or a size such as `5m`.
        genome_size: String,
    },
    Concatemer {
        min_identity: u8,
        min_unit: usize,
    },
    Circularize {
        min_overlap: usize,
        min_identity: u8,
    },
    Rotate {
        genes: String,
    },
    Medaka {
        reads: String,
    },
}

fn number<T: FromStr>(step: &str, key: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{step}: {key} must be a number, not {value}"))
}

fn percent(step: &str, key: &str, value: &str) -> Result<u8, String> {
    let n: u8 = number(step, key, value)?;
    if n > 100 {
        return Err(format!("{step}: {key} is a share from 0 to 100, not {n}"));
    }
    Ok(n)
}

/// A setting that is `auto` or a number that `ok` accepts.
fn auto_or(
    step: &str,
    key: &str,
    value: &str,
    ok: impl Fn(f64) -> bool,
    said: &str,
) -> Result<String, String> {
    if value == "auto" || value.parse::<f64>().is_ok_and(|n| n.is_finite() && ok(n)) {
        Ok(value.to_string())
    } else {
        Err(format!("{step}: {key} is auto or {said}, not {value}"))
    }
}

fn one_of(step: &str, key: &str, value: &str, words: &[&str]) -> Result<String, String> {
    if words.contains(&value) {
        Ok(value.to_string())
    } else {
        Err(format!(
            "{step}: {key} is {}, not {value}",
            words.join(", ")
        ))
    }
}

impl AsmStep {
    pub fn name(&self) -> &'static str {
        match self {
            AsmStep::Trim { .. } => "trim",
            AsmStep::Filter { .. } => "filter",
            AsmStep::Assemble { .. } => "assemble",
            AsmStep::Hybrid { .. } => "hybrid",
            AsmStep::Flye { .. } => "flye",
            AsmStep::Concatemer { .. } => "concatemer",
            AsmStep::Circularize { .. } => "circularize",
            AsmStep::Rotate { .. } => "rotate",
            AsmStep::Medaka { .. } => "medaka",
        }
    }

    pub fn takes(&self) -> Kind {
        match self {
            AsmStep::Trim { .. }
            | AsmStep::Filter { .. }
            | AsmStep::Assemble { .. }
            | AsmStep::Flye { .. } => Kind::Reads,
            AsmStep::Hybrid { .. } => Kind::Gfa,
            _ => Kind::Sequences,
        }
    }

    pub fn gives(&self) -> Kind {
        match self {
            AsmStep::Trim { .. } | AsmStep::Filter { .. } => Kind::Reads,
            AsmStep::Assemble { .. } => Kind::Gfa,
            _ => Kind::Sequences,
        }
    }

    /// All the settings, in the order they are written.
    pub fn params(&self) -> Vec<(&'static str, String)> {
        let s = |v: &str| v.to_string();
        match self {
            AsmStep::Trim {
                adapter,
                quality,
                min_length,
                cut,
            } => vec![
                ("adapter", s(adapter)),
                ("quality", quality.to_string()),
                ("min-length", min_length.to_string()),
                ("cut", s(cut.word())),
            ],
            AsmStep::Filter {
                min_length,
                target_bases,
                window,
            } => vec![
                ("min-length", min_length.to_string()),
                ("target-bases", target_bases.to_string()),
                ("window", window.to_string()),
            ],
            AsmStep::Assemble {
                ks,
                cutoff,
                min_qual,
            } => vec![
                ("ks", s(ks)),
                ("cov-cutoff", s(cutoff)),
                ("min-qual", min_qual.to_string()),
            ],
            AsmStep::Hybrid {
                long,
                mode,
                min_quality,
                single_depth,
                overlap,
            } => vec![
                ("long", s(long)),
                ("mode", s(mode)),
                ("min-quality", s(min_quality)),
                ("single-depth", s(single_depth)),
                ("overlap", overlap.to_string()),
            ],
            AsmStep::Flye { mode, genome_size } => {
                vec![("mode", s(mode)), ("genome-size", s(genome_size))]
            }
            AsmStep::Concatemer {
                min_identity,
                min_unit,
            } => vec![
                ("min-identity", min_identity.to_string()),
                ("min-unit", min_unit.to_string()),
            ],
            AsmStep::Circularize {
                min_overlap,
                min_identity,
            } => vec![
                ("min-overlap", min_overlap.to_string()),
                ("min-identity", min_identity.to_string()),
            ],
            AsmStep::Rotate { genes } => vec![("genes", s(genes))],
            AsmStep::Medaka { reads } => vec![("reads", s(reads))],
        }
    }

    /// The step of this name with its settings; `None` when no step here has the name.
    pub fn parse(name: &str, params: &[(String, String)]) -> Result<Option<AsmStep>, String> {
        let auto = || "auto".to_string();
        let mut step = match name {
            "trim" => AsmStep::Trim {
                adapter: "AGATCGGAAGAGC".into(),
                quality: 20,
                min_length: 20,
                cut: Cut::Tail,
            },
            "filter" => AsmStep::Filter {
                min_length: 2000,
                target_bases: 500_000_000,
                window: 100,
            },
            "assemble" => AsmStep::Assemble {
                ks: auto(),
                cutoff: auto(),
                min_qual: 10,
            },
            "hybrid" => AsmStep::Hybrid {
                long: String::new(),
                mode: "normal".into(),
                min_quality: auto(),
                single_depth: auto(),
                overlap: 0,
            },
            "flye" => AsmStep::Flye {
                mode: "nano-raw".into(),
                genome_size: auto(),
            },
            "concatemer" => AsmStep::Concatemer {
                min_identity: 95,
                min_unit: 200,
            },
            "circularize" => AsmStep::Circularize {
                min_overlap: 30,
                min_identity: 95,
            },
            "rotate" => AsmStep::Rotate {
                genes: String::new(),
            },
            "medaka" => AsmStep::Medaka {
                reads: String::new(),
            },
            _ => return Ok(None),
        };
        for (key, value) in params {
            step.set(name, key, value)?;
        }
        match &step {
            AsmStep::Hybrid { long, .. } if long.is_empty() => {
                return Err("hybrid needs long=FILE or long=FOLDER: the long reads".into())
            }
            AsmStep::Rotate { genes } if genes.is_empty() => {
                return Err("rotate needs genes=FILE: a FASTA of the start genes".into())
            }
            AsmStep::Medaka { reads } if reads.is_empty() => {
                return Err(
                    "medaka needs reads=FILE or reads=FOLDER: the reads to polish with".into(),
                )
            }
            _ => {}
        }
        Ok(Some(step))
    }

    fn set(&mut self, name: &str, key: &str, v: &str) -> Result<(), String> {
        match (self, key) {
            (AsmStep::Trim { adapter, .. }, "adapter") => {
                if v != "none"
                    && (v.len() < 4
                        || !v
                            .bytes()
                            .all(|b| b"ACGTN".contains(&b.to_ascii_uppercase())))
                {
                    return Err(format!(
                        "trim: adapter is none or at least four bases of A, C, G, T, not {v}"
                    ));
                }
                *adapter = v.to_string();
            }
            (AsmStep::Trim { quality, .. }, "quality") => *quality = number(name, key, v)?,
            (AsmStep::Trim { min_length, .. }, "min-length")
            | (AsmStep::Filter { min_length, .. }, "min-length") => {
                *min_length = number(name, key, v)?
            }
            (AsmStep::Trim { cut, .. }, "cut") => {
                *cut = match one_of(name, key, v, &["tail", "right", "off"])?.as_str() {
                    "tail" => Cut::Tail,
                    "right" => Cut::Right,
                    _ => Cut::Off,
                }
            }
            (AsmStep::Filter { target_bases, .. }, "target-bases") => {
                *target_bases = number(name, key, v)?
            }
            (AsmStep::Filter { window, .. }, "window") => *window = number(name, key, v)?,
            (AsmStep::Assemble { ks, .. }, "ks") => {
                if v != "auto" {
                    graph::parse_ks(v).map_err(|why| format!("assemble: {why}"))?;
                }
                *ks = v.to_string();
            }
            (AsmStep::Assemble { cutoff, .. }, "cov-cutoff") => {
                *cutoff = if v == "off" {
                    v.to_string()
                } else {
                    auto_or(name, key, v, |n| n >= 0.0, "off or a number of 0 or more")?
                }
            }
            (AsmStep::Assemble { min_qual, .. }, "min-qual") => *min_qual = number(name, key, v)?,
            (AsmStep::Hybrid { long, .. }, "long") => *long = v.to_string(),
            (AsmStep::Hybrid { mode, .. }, "mode") => {
                *mode = one_of(name, key, v, &["conservative", "normal", "bold"])?
            }
            (AsmStep::Hybrid { min_quality, .. }, "min-quality") => {
                *min_quality = auto_or(
                    name,
                    key,
                    v,
                    |n| (0.0..=100.0).contains(&n),
                    "a number from 0 to 100",
                )?
            }
            (AsmStep::Hybrid { single_depth, .. }, "single-depth") => {
                *single_depth = auto_or(name, key, v, |n| n > 0.0, "a number above 0")?
            }
            (AsmStep::Hybrid { overlap, .. }, "overlap") => *overlap = number(name, key, v)?,
            (AsmStep::Flye { mode, .. }, "mode") => {
                *mode = one_of(name, key, v, &["nano-raw", "nano-corr", "nano-hq"])?
            }
            (AsmStep::Flye { genome_size, .. }, "genome-size") => {
                let digits = v.trim_end_matches(['k', 'K', 'm', 'M', 'g', 'G']);
                if v != "auto"
                    && !(digits.parse::<f64>().is_ok_and(|n| n > 0.0)
                        && v.len() - digits.len() <= 1)
                {
                    return Err(format!(
                        "flye: genome-size is auto or a size such as 5m or 60k, not {v}"
                    ));
                }
                *genome_size = v.to_string();
            }
            (AsmStep::Concatemer { min_identity, .. }, "min-identity")
            | (AsmStep::Circularize { min_identity, .. }, "min-identity") => {
                *min_identity = percent(name, key, v)?
            }
            (AsmStep::Concatemer { min_unit, .. }, "min-unit") => *min_unit = number(name, key, v)?,
            (AsmStep::Circularize { min_overlap, .. }, "min-overlap") => {
                *min_overlap = number(name, key, v)?
            }
            (AsmStep::Rotate { genes }, "genes") => *genes = v.to_string(),
            (AsmStep::Medaka { reads }, "reads") => *reads = v.to_string(),
            (_, other) => return Err(format!("{name} has no setting called {other}")),
        }
        Ok(())
    }

    /// The file or folder that a setting names for the step to read besides what flows to it.
    fn side(&self) -> Option<(&str, &[&'static str])> {
        match self {
            AsmStep::Hybrid { long, .. } => Some((long, &["fastq", "fq", "fasta", "fa"])),
            AsmStep::Medaka { reads, .. } => Some((reads, &["fastq", "fq"])),
            AsmStep::Rotate { genes } => Some((genes, &["fasta", "fa", "fna"])),
            _ => None,
        }
    }

    /// The file that the step reads for the input called `name`: the file the setting names, or in the folder it names the file with the
    /// name of the input before its ending.
    pub fn side_file(&self, name: &str) -> Result<Option<PathBuf>, String> {
        let Some((setting, endings)) = self.side() else {
            return Ok(None);
        };
        let path = Path::new(setting);
        if path.is_file() {
            return Ok(Some(path.to_path_buf()));
        }
        if path.is_dir() {
            return endings
                .iter()
                .map(|e| path.join(format!("{name}.{e}")))
                .find(|p| p.is_file())
                .map(Some)
                .ok_or_else(|| {
                    format!(
                        "{} holds no file called {name} with an ending of {} for {}",
                        path.display(),
                        endings.join(", "),
                        self.name()
                    )
                });
        }
        Err(format!(
            "{setting}, which {} is to read, is not a file or a folder",
            self.name()
        ))
    }

    /// What is added to the key of an answer besides the step's settings and its input: a checksum of the side file (so that a changed
    /// file is not read from the cache), and the version of the installed program.
    pub fn context(&self, name: &str) -> String {
        let side = match self.side_file(name) {
            Ok(Some(p)) => std::fs::read(&p).map_or("unreadable".to_string(), |b| digest(&b)),
            Ok(None) => String::new(),
            Err(_) => "missing".to_string(),
        };
        format!("{}|{side}", external::version(self))
    }

    /// Do the step on `input`, which is the input called `name`.
    pub fn work(&self, input: &Value, name: &str) -> Result<(Value, String), String> {
        let side = self.side_file(name)?;
        match (self, input) {
            (
                AsmStep::Trim {
                    adapter,
                    quality,
                    min_length,
                    cut,
                },
                Value::Reads(r),
            ) => reads::trim(r, adapter, *quality, *min_length, *cut),
            (
                AsmStep::Filter {
                    min_length,
                    target_bases,
                    window,
                },
                Value::Reads(r),
            ) => reads::filter(r, *min_length, *target_bases, *window),
            (
                AsmStep::Assemble {
                    ks,
                    cutoff,
                    min_qual,
                },
                Value::Reads(r),
            ) => graph::assemble(r, ks, cutoff, *min_qual),
            (
                AsmStep::Hybrid {
                    mode,
                    min_quality,
                    single_depth,
                    overlap,
                    ..
                },
                Value::Gfa(text),
            ) => graph::hybrid(
                text,
                side.as_deref(),
                mode,
                min_quality,
                single_depth,
                *overlap,
            ),
            (AsmStep::Flye { .. }, Value::Reads(r)) => external::flye(self, r),
            (
                AsmStep::Concatemer {
                    min_identity,
                    min_unit,
                },
                Value::Sequences(rows),
            ) => contigs::concatemer(rows, *min_identity, *min_unit),
            (
                AsmStep::Circularize {
                    min_overlap,
                    min_identity,
                },
                Value::Sequences(rows),
            ) => contigs::circularize(rows, *min_overlap, *min_identity),
            (AsmStep::Rotate { .. }, Value::Sequences(rows)) => {
                contigs::rotate(rows, side.as_deref())
            }
            (AsmStep::Medaka { .. }, Value::Sequences(rows)) => {
                external::medaka(self, rows, side.as_deref())
            }
            (step, input) => Err(format!(
                "{} cannot take {}",
                step.name(),
                input.kind().word()
            )),
        }
    }

    /// Whether the program that the step calls is there: an error that says so when it is not.
    pub fn ready(&self) -> Result<(), String> {
        external::ready(self)
    }
}

/// `name` with a word of what was found put after it, once.
pub fn tag(name: &str, word: &str) -> String {
    if has_tag(name, word) {
        name.to_string()
    } else {
        format!("{name} {word}")
    }
}

pub fn has_tag(name: &str, word: &str) -> bool {
    name.split_whitespace().any(|w| w == word)
}

#[cfg(test)]
mod tests;
