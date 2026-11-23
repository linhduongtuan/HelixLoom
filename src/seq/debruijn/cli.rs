//! `helixloom dbg`: assemble short reads with a de Bruijn graph.

use super::assemble::{assemble, Cutoff, Options, DEFAULT_KS};
use super::output::{fasta, gfa, number, summary};
use super::reads::{load, pieces};
use std::path::Path;

pub const USAGE: &str = "\
Usage: helixloom dbg [options] READS...

  Assemble short reads (FASTA or FASTQ, not compressed) with a de Bruijn graph, in the manner of SPAdes (Bankevich et al. 2012; docs/DEBRUIJN.md),
  for plasmids and small genomes. Reads that come in pairs are used as reads that do not: -1 and -2 name the two files of a pair, -s one of
  unpaired reads, and files without a flag are reads too. What is written is the contigs, as FASTA, named NODE_<number>_length_<bases>_cov_<coverage>
  (to the standard output unless --fasta or -o says where), and with --gfa or -o the graph, as GFA. A summary of the work goes to the error output.

  -1 FILE, -2 FILE, -s FILE   reads (any number of each; several files may be joined with commas)
  -k LIST        the k-mer lengths, odd and rising, from 3 to 63, as 21,33,55: each is made from the reads and the contigs of the one before
                 (default: those of 21, 33 and 55 that fit in the longest read)
  --cov-cutoff X  take away what is seen fewer times than X: a number, off, or auto (default auto: a tenth of the usual coverage, at least 1.5)
  --min-count N  drop the k-mers that are seen fewer than N times before the graph is made (default 1: none)
  --tip-length N the longest tip that is clipped, in k-mers (default: the length of the longest read)
  --tip-ratio X  a tip is kept if its coverage is more than X times that alongside it (default 2.5)
  --bulge-length N  the longest path that is removed as a bulge, in k-mers (default 150)
  --bulge-diff X the share of its bases by which the two paths of a bulge may differ (default 0.1)
  --isolated-length N  at the last k, a segment that is joined to nothing and is shorter than N bases is removed (default 200; 0: none is)
  --min-qual N   a base of lower quality is not used (FASTQ; default 10)
  --min-length N the shortest contig that is written to the FASTA (default 200)
  -o DIR         write DIR/contigs.fasta and DIR/assembly_graph.gfa (DIR is made if it is not there)
  --fasta FILE   where to write the contigs
  --gfa FILE     where to write the graph
  -q             no summary
  -h, --help     this text
";

struct Args {
    options: Options,
    files: Vec<String>,
    min_qual: u8,
    min_length: usize,
    fasta: Option<String>,
    gfa: Option<String>,
    quiet: bool,
}

fn number_of<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} needs a number, not {value}."))
}

fn ratio(flag: &str, value: &str) -> Result<f64, String> {
    let v: f64 = number_of(flag, value)?;
    if v.is_finite() && v >= 0.0 {
        Ok(v)
    } else {
        Err(format!("{flag} needs a number of 0 or more, not {value}."))
    }
}

fn parse(args: &[String]) -> Result<Args, String> {
    let mut a = Args {
        options: Options::default(),
        files: Vec::new(),
        min_qual: 10,
        min_length: 200,
        fasta: None,
        gfa: None,
        quiet: false,
    };
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value."))
        };
        match flag.as_str() {
            "-1" | "-2" | "-s" => a.files.push(value(flag)?),
            "-k" => {
                a.options.ks = value("-k")?
                    .split(',')
                    .map(|k| number_of("-k", k.trim()))
                    .collect::<Result<_, _>>()?;
            }
            "--cov-cutoff" => {
                let v = value("--cov-cutoff")?;
                a.options.cutoff = match v.as_str() {
                    "off" => Cutoff::Off,
                    "auto" => Cutoff::Auto,
                    other => Cutoff::Value(ratio("--cov-cutoff", other)?),
                };
            }
            "--min-count" => {
                a.options.min_count = number_of("--min-count", &value("--min-count")?)?
            }
            "--tip-length" => {
                a.options.tip_length = Some(number_of("--tip-length", &value("--tip-length")?)?)
            }
            "--tip-ratio" => a.options.tip_ratio = ratio("--tip-ratio", &value("--tip-ratio")?)?,
            "--bulge-length" => {
                a.options.bulge_length =
                    Some(number_of("--bulge-length", &value("--bulge-length")?)?)
            }
            "--bulge-diff" => {
                a.options.bulge_diff = ratio("--bulge-diff", &value("--bulge-diff")?)?
            }
            "--isolated-length" => {
                a.options.isolated_length =
                    number_of("--isolated-length", &value("--isolated-length")?)?
            }
            "--min-qual" => a.min_qual = number_of("--min-qual", &value("--min-qual")?)?,
            "--min-length" => a.min_length = number_of("--min-length", &value("--min-length")?)?,
            "-o" => {
                let dir = value("-o")?;
                a.fasta = Some(format!("{dir}/contigs.fasta"));
                a.gfa = Some(format!("{dir}/assembly_graph.gfa"));
            }
            "--fasta" => a.fasta = Some(value("--fasta")?),
            "--gfa" => a.gfa = Some(value("--gfa")?),
            "-q" => a.quiet = true,
            other if other.starts_with('-') && other.len() > 1 => {
                return Err(format!(
                    "{other} is not an option (`helixloom dbg --help`)."
                ))
            }
            file => a.files.push(file.to_string()),
        }
    }
    if a.files.is_empty() {
        return Err("Name the FASTA or FASTQ file of the reads to assemble.".into());
    }
    if a.options.min_count == 0 {
        return Err("--min-count is 1 or more.".into());
    }
    Ok(a)
}

fn write(path: &str, text: &str) -> Result<(), String> {
    if let Some(dir) = Path::new(path)
        .parent()
        .filter(|d| !d.as_os_str().is_empty())
    {
        std::fs::create_dir_all(dir)
            .map_err(|e| format!("Could not make {}: {e}", dir.display()))?;
    }
    std::fs::write(path, text).map_err(|e| format!("Could not write {path}: {e}"))
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
    let mut reads = Vec::new();
    for spec in &a.files {
        reads.extend(load(spec)?);
    }
    let shortest = a.options.ks.first().copied().unwrap_or(DEFAULT_KS[0]);
    let used = pieces(&reads, a.min_qual, shortest);
    if used.is_empty() {
        return Err(format!(
            "No run of {shortest} bases or more is left in the reads once the unknown bases and those under quality {} are cut out.",
            a.min_qual
        ));
    }
    let assembly = assemble(&used, &a.options)?;
    let contigs = number(&assembly.graph);
    let text = fasta(&contigs, a.min_length);
    match &a.fasta {
        Some(path) => write(path, &text)?,
        None => print!("{text}"),
    }
    if let Some(path) = &a.gfa {
        write(path, &gfa(&assembly.graph, &contigs))?;
    }
    if !a.quiet {
        eprint!("{}", summary(&assembly, &contigs, a.min_length));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
