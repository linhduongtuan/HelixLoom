//! The motif job: DNA sequences go in, the motifs the search found and the places they were found at come out (`seq/motif`). Like the
//! other kinds it runs the same function on the workstation and on this computer, and travels as a container of sections: the
//! settings as `key=value` lines, the sequences one to a line, the motifs as text.

use super::tasks::{bit, flag, lines, number};
use super::wire::{meta, Container};
use crate::seq::motif::{check_size, discover, Found, Motif, Settings, Site, EFFORT_RANGE};
use std::collections::HashMap;
use web_time::Instant;

/// The most sequences a request may hold, and the most motifs, seeds and rounds it may ask for.
pub const MAX_SEQUENCES: usize = 100_000;
pub const MAX_MOTIFS: usize = 50;
const MAX_SEEDS: usize = 1_000;
const MAX_ITERATIONS: usize = 10_000;

/// What to run: the sequences and the settings of the search (see `motif::Settings`).
#[derive(Clone, Debug, PartialEq)]
pub struct MotifJob {
    pub seqs: Vec<Vec<u8>>,
    pub settings: Settings,
}

impl MotifJob {
    pub fn encode(&self) -> Container {
        let s = &self.settings;
        Container::new()
            .with(
                "meta",
                lines(&[
                    ("width", s.width.to_string()),
                    ("motifs", s.motifs.to_string()),
                    ("both_strands", bit(s.both_strands)),
                    ("seeds", s.seeds.to_string()),
                    ("iterations", s.iterations.to_string()),
                    ("effort", s.effort.to_string()),
                ]),
            )
            .with(
                "sequences",
                self.seqs
                    .iter()
                    .map(|q| format!("{}\n", String::from_utf8_lossy(q)))
                    .collect::<String>(),
            )
    }

    /// A request comes from the network: every setting is checked, and the sequences are counted and measured before any window
    /// of them is made.
    pub fn decode(c: &Container, m: &HashMap<&str, &str>) -> Result<MotifJob, String> {
        let settings = Settings {
            width: number(m, "width")?,
            motifs: number(m, "motifs")?,
            both_strands: flag(m, "both_strands"),
            seeds: number(m, "seeds")?,
            iterations: number(m, "iterations")?,
            effort: number(m, "effort")?,
        };
        if !(3..=40).contains(&settings.width)
            || !(1..=MAX_MOTIFS).contains(&settings.motifs)
            || !(1..=MAX_SEEDS).contains(&settings.seeds)
            || !(1..=MAX_ITERATIONS).contains(&settings.iterations)
            || !(EFFORT_RANGE.0..=EFFORT_RANGE.1).contains(&settings.effort)
        {
            return Err("The message's search settings are out of range.".into());
        }
        let text = c.require("sequences")?;
        let count = text.lines().count();
        if count > MAX_SEQUENCES {
            return Err(format!(
                "A job may hold at most {MAX_SEQUENCES} sequences; this one has {count}."
            ));
        }
        let seqs: Vec<Vec<u8>> = text.lines().map(|l| l.as_bytes().to_vec()).collect();
        if seqs.iter().flatten().any(|b| !b.is_ascii_alphabetic()) {
            return Err("A sequence in the message has a character that is no letter.".into());
        }
        check_size(
            seqs.iter().map(Vec::len).sum(),
            settings.width,
            settings.both_strands,
        )?;
        Ok(MotifJob { seqs, settings })
    }
}

/// What a search found.
#[derive(Clone, Debug, PartialEq)]
pub struct MotifResult {
    pub found: Vec<Found>,
    pub seconds: f64,
}

/// Run the job here.
pub fn run(job: &MotifJob) -> Result<MotifResult, String> {
    let started = Instant::now();
    let found = discover(&job.seqs, &job.settings)?;
    Ok(MotifResult {
        found,
        seconds: started.elapsed().as_secs_f64(),
    })
}

impl MotifResult {
    pub fn encode(&self) -> Container {
        let mut text = String::new();
        for f in &self.found {
            text.push_str(&format!(
                "motif\t{}\t{}\t{}\t{}\n",
                f.motif.name.replace(['\t', '\n', '\r'], " "),
                f.motif.width(),
                f.gain,
                f.sites.len()
            ));
            for c in &f.motif.probs {
                text.push_str(&format!("{}\t{}\t{}\t{}\n", c[0], c[1], c[2], c[3]));
            }
            for s in &f.sites {
                text.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    s.sequence,
                    s.position,
                    u8::from(s.reverse),
                    s.score,
                    s.probability
                ));
            }
        }
        Container::new()
            .with(
                "meta",
                lines(&[
                    ("seconds", self.seconds.to_string()),
                    ("motifs", self.found.len().to_string()),
                ]),
            )
            .with("motifs", text)
    }

    /// A result comes from a server that is not trusted to keep to what it was asked: everything is bounded and checked before
    /// anything is made of it.
    pub fn decode(c: &Container) -> Result<MotifResult, String> {
        let m = meta(c.require("meta")?);
        let bad = |why: &str| format!("The result's motifs are damaged: {why}.");
        let seconds: f64 = number(&m, "seconds")?;
        let asked: usize = number(&m, "motifs")?;
        if asked > MAX_MOTIFS {
            return Err(bad("there are too many"));
        }
        let mut rows = c.require("motifs")?.lines();
        let mut found = Vec::new();
        while let Some(head) = rows.next() {
            if found.len() == asked {
                return Err(bad("there are more than it says"));
            }
            let h: Vec<&str> = head.split('\t').collect();
            let ["motif", name, width, gain, sites] = h[..] else {
                return Err(bad("a motif has no heading"));
            };
            let (width, sites): (usize, usize) = (
                width.parse().map_err(|_| bad("a width"))?,
                sites.parse().map_err(|_| bad("a number of sites"))?,
            );
            let gain: f64 = gain.parse().map_err(|_| bad("a gain"))?;
            if !(3..=40).contains(&width) || sites > MAX_SEQUENCES || !gain.is_finite() {
                return Err(bad("a motif is out of range"));
            }
            if name.chars().count() > 100 || name.chars().any(char::is_control) {
                return Err(bad("a name"));
            }
            let mut probs = Vec::with_capacity(width);
            for _ in 0..width {
                let row = rows.next().ok_or_else(|| bad("a column is missing"))?;
                let p: Vec<f64> = row
                    .split('\t')
                    .map(|x| x.parse::<f64>().map_err(|_| bad("a probability")))
                    .collect::<Result<_, _>>()?;
                let [a, b, c, d] = p[..] else {
                    return Err(bad("a column is not four numbers"));
                };
                let column = [a, b, c, d];
                if column.iter().any(|x| !(0.0..=1.0).contains(x))
                    || (column.iter().sum::<f64>() - 1.0).abs() > 1e-6
                {
                    return Err(bad("a column is no distribution"));
                }
                probs.push(column);
            }
            let mut places = Vec::with_capacity(sites);
            for _ in 0..sites {
                let row = rows.next().ok_or_else(|| bad("a site is missing"))?;
                let s: Vec<&str> = row.split('\t').collect();
                let [sequence, position, reverse, score, probability] = s[..] else {
                    return Err(bad("a site has the wrong fields"));
                };
                let site = Site {
                    sequence: sequence.parse().map_err(|_| bad("a sequence number"))?,
                    position: position.parse().map_err(|_| bad("a position"))?,
                    reverse: match reverse {
                        "0" => false,
                        "1" => true,
                        _ => return Err(bad("a strand")),
                    },
                    score: score.parse().map_err(|_| bad("a score"))?,
                    probability: probability.parse().map_err(|_| bad("a probability"))?,
                };
                if site.sequence >= MAX_SEQUENCES
                    || !site.score.is_finite()
                    || !(0.0..=1.0).contains(&site.probability)
                {
                    return Err(bad("a site is out of range"));
                }
                places.push(site);
            }
            found.push(Found {
                motif: Motif {
                    name: name.to_string(),
                    probs,
                },
                sites: places,
                gain,
            });
        }
        if found.len() != asked {
            return Err(bad("there are fewer than it says"));
        }
        Ok(MotifResult { found, seconds })
    }

    /// Whether the result is about these sequences: every site is a window inside a sequence that exists. A server's answer is
    /// not taken for one until it is.
    pub fn fits(&self, seqs: &[Vec<u8>]) -> Result<(), String> {
        for f in &self.found {
            for s in &f.sites {
                let long = seqs.get(s.sequence).map_or(0, Vec::len);
                if s.position + f.motif.width() > long {
                    return Err(
                        "The workstation's motifs do not fit the sequences that were sent.".into(),
                    );
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
