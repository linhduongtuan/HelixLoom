//! `helixloom correct`: correct long reads with the reads that overlap them.

use super::{correct, Params};
use crate::seq::variants::reads::{self, Read};

pub const USAGE: &str = "\
Usage: helixloom correct [options] READS...

  Correct long noisy reads (nanopore) with the reads that overlap them, in two rounds, after NECAT (Chen et al., Nature Communications 2021),
  written from its paper (docs/CORRECT.md). READS are FASTA or FASTQ files (not compressed). The answer is FASTA. For plasmids and small
  genomes: every read is held in memory. NECAT was not run, and its accuracy is not claimed.

  -o FILE             where to write the corrected reads (the standard output when left out)
  --rounds N          rounds of correction, 1 to 5 (default 2); one round keeps every read whole, two or more write only the corrected
                      stretches of the last round
  --coverage N        pass over a candidate whose whole span has more reads than this (default 12)
  --candidates N      the best candidates looked at for each read (default 200)
  --min-support N     the fewest supporting reads to correct a position (default 4)
  --min-length N      the shortest corrected stretch written by the last round (default 500)
  -k N, -w N          the minimizers of the overlaps and alignments (default 15 and 5)
  -h, --help          this text
";

struct Args {
    params: Params,
    output: Option<String>,
    inputs: Vec<String>,
}

fn number<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} needs a whole number, not {value}."))
}

fn parse(args: &[String]) -> Result<Args, String> {
    let mut a = Args {
        params: Params::default(),
        output: None,
        inputs: Vec::new(),
    };
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value."))
        };
        match flag.as_str() {
            "-o" | "--output" => a.output = Some(value("-o")?),
            "--rounds" => a.params.rounds = number("--rounds", &value("--rounds")?)?,
            "--coverage" => a.params.coverage = number("--coverage", &value("--coverage")?)?,
            "--candidates" => {
                a.params.max_candidates = number("--candidates", &value("--candidates")?)?
            }
            "--min-support" => {
                a.params.min_support = number("--min-support", &value("--min-support")?)?
            }
            "--min-length" => a.params.min_piece = number("--min-length", &value("--min-length")?)?,
            "-k" => a.params.k = number("-k", &value("-k")?)?,
            "-w" => a.params.w = number("-w", &value("-w")?)?,
            other if other.starts_with('-') && other.len() > 1 => {
                return Err(format!(
                    "{other} is not an option (`helixloom correct --help`)."
                ))
            }
            path => a.inputs.push(path.to_string()),
        }
    }
    if a.inputs.is_empty() {
        return Err("Name the reads to correct.".into());
    }
    if !(1..=5).contains(&a.params.rounds) {
        return Err("--rounds must be from 1 to 5.".into());
    }
    if a.params.max_candidates == 0 || a.params.min_support == 0 {
        return Err("--candidates and --min-support must be at least 1.".into());
    }
    if !(1..=32).contains(&a.params.k) {
        return Err(format!("k must be from 1 to 32, not {}.", a.params.k));
    }
    if a.params.w == 0 {
        return Err("w must be 1 or more.".into());
    }
    Ok(a)
}

/// FASTA of the reads, in lines of 80.
fn fasta(reads: &[Read]) -> String {
    let mut out = String::new();
    for r in reads {
        out.push('>');
        out.push_str(&r.name);
        out.push('\n');
        for line in r.seq.chunks(80) {
            out.push_str(&String::from_utf8_lossy(line));
            out.push('\n');
        }
    }
    out
}

fn summary(before: &[Read], after: &[Read]) -> String {
    let bases = |r: &[Read]| r.iter().map(|x| x.seq.len()).sum::<usize>();
    format!(
        "{} reads of {} bases in, {} of {} bases out.\n",
        before.len(),
        bases(before),
        after.len(),
        bases(after)
    )
}

fn run_inner(args: &[String]) -> Result<(), String> {
    let a = parse(args)?;
    let input = reads::load(&a.inputs.join(","))?;
    let corrected = correct(&input, &a.params)?;
    let text = fasta(&corrected);
    match &a.output {
        Some(path) => {
            std::fs::write(path, text).map_err(|e| format!("Could not write {path}: {e}"))?
        }
        None => print!("{text}"),
    }
    eprint!("{}", summary(&input, &corrected));
    Ok(())
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

#[cfg(test)]
mod tests;
