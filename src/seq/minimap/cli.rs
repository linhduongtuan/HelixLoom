//! `helixloom map`: map reads to a reference, or find the overlaps between reads, and write PAF.

use super::{map_text, overlap_text, Options, MAX_K};
use crate::seq::kmer::cli::read_input;

/// The sequences of an input, each named by the first word of its header (a name in PAF has no spaces, as minimap2 writes it).
fn sequences(spec: &str) -> Result<Vec<(String, Vec<u8>)>, String> {
    Ok(read_input(spec)?
        .into_iter()
        .map(|(name, bases)| {
            (
                name.split_whitespace()
                    .next()
                    .unwrap_or("unnamed")
                    .to_string(),
                bases,
            )
        })
        .collect())
}

pub const USAGE: &str = "\
Usage: helixloom map [options] TARGET [READS]

  Map reads to a reference, in the manner of minimap2 (Li, Bioinformatics 2018), written from its paper and manual (docs/MAP.md). TARGET and
  READS are FASTA or FASTQ files (not compressed); several files that are to be one input are joined with commas. The answer is PAF.

  -x PRESET     map-ont (default: noisy long reads to a reference) or ava-ont (overlaps between reads: give one file, TARGET, which is
                mapped to itself, each pair once)
  -k N          the length of a k-mer, 1 to 32 (default 15)
  -w N          the number of k-mers in a window (default 10; 5 for ava-ont)
  -f FLOAT      ignore this share of the most frequent minimizers (default 0.0002)
  -g N          the most that two anchors of a chain may be apart (default 10000)
  -r N          the most that the gaps of two anchors may differ (default 500; 2000 for ava-ont)
  -n N          the fewest anchors in a chain (default 3)
  -m N          the lowest score of a chain (default 40; 100 for ava-ont)
  -p FLOAT      a secondary mapping must score this share of its primary (default 0.8)
  -N N          keep at most this many secondary mappings (default 5)
  -M FLOAT      a mapping is secondary to one it overlaps on the read by this share of the shorter (default 0.5)
  -c            align the bases along each chain: the CIGAR, NM, and exact matches (slower; the alignment runs from the first anchor to the last)
  -o FILE       where to write the PAF (the standard output when left out)
  -h, --help    this text
";

fn number<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} needs a number, not {value}."))
}

struct Args {
    options: Options,
    output: Option<String>,
    inputs: Vec<String>,
}

fn parse(args: &[String]) -> Result<Args, String> {
    // the preset first, so that the flags that follow change what it sets
    let mut options = Options::default();
    if let Some(at) = args.iter().position(|a| a == "-x") {
        let name = args
            .get(at + 1)
            .ok_or("-x needs a preset: map-ont or ava-ont.")?;
        options = Options::preset(name).ok_or_else(|| {
            format!("{name} is not a preset this program has: map-ont or ava-ont.")
        })?;
    }
    let (mut output, mut inputs) = (None, Vec::new());
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value."))
        };
        match flag.as_str() {
            "-x" => {
                value("-x")?;
            }
            "-k" => options.k = number("-k", &value("-k")?)?,
            "-w" => options.w = number("-w", &value("-w")?)?,
            "-f" => options.occ_fraction = number("-f", &value("-f")?)?,
            "-g" => options.max_gap = number("-g", &value("-g")?)?,
            "-r" => options.bandwidth = number("-r", &value("-r")?)?,
            "-n" => options.min_anchors = number("-n", &value("-n")?)?,
            "-m" => options.min_score = number("-m", &value("-m")?)?,
            "-p" => options.secondary_ratio = number("-p", &value("-p")?)?,
            "-N" => options.max_secondary = number("-N", &value("-N")?)?,
            "-M" => options.mask_level = number("-M", &value("-M")?)?,
            "-c" => options.cigar = true,
            "-o" => output = Some(value("-o")?),
            other if other.starts_with('-') && other.len() > 1 => {
                return Err(format!(
                    "{other} is not an option (`helixloom map --help`)."
                ))
            }
            path => inputs.push(path.to_string()),
        }
    }
    if !(1..=MAX_K).contains(&options.k) {
        return Err(format!("k must be from 1 to {MAX_K}, not {}.", options.k));
    }
    if options.w == 0 {
        return Err("w must be 1 or more.".into());
    }
    Ok(Args {
        options,
        output,
        inputs,
    })
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
    let paf = match (a.options.overlaps, a.inputs.as_slice()) {
        (true, [reads]) => overlap_text(&sequences(reads)?, &a.options)?,
        (true, _) => return Err("Overlaps (-x ava-ont) take one file of reads, which is mapped to itself.".into()),
        (false, [target, reads]) => map_text(&sequences(target)?, &sequences(reads)?, &a.options)?,
        (false, [_]) => return Err("Name the reads after the target: helixloom map TARGET READS (or -x ava-ont with one file).".into()),
        (false, _) => return Err("Name a target and a file of reads: helixloom map TARGET READS.".into()),
    };
    match &a.output {
        Some(path) => std::fs::write(path, paf).map_err(|e| format!("Could not write {path}: {e}")),
        None => {
            print!("{paf}");
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests;
