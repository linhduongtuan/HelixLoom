//! `helixloom hybrid`: join an assembly graph with long reads.

use super::apply::Mode;
use super::assemble::{run as assemble, Options};
use super::input;
use super::output::{fasta, replicons};
use super::rotate::Options as RotateOptions;
use crate::seq::debruijn::output::{gfa, number};
use crate::seq::kmer::cli::read_input;
use std::path::Path;

pub const USAGE: &str = "\
Usage: helixloom hybrid [options] GRAPH.gfa LONGREADS...

  Join the segments of a short-read assembly graph with long reads, in the manner of Unicycler (Wick et al. 2017; docs/HYBRID.md). GRAPH is a
  GFA (from `helixloom dbg`, or SPAdes) whose segments have their bases and depth, and whose links overlap by k-1 bases; LONGREADS are FASTA
  or FASTQ files (not compressed). Each segment is given a number of copies from its depth and its links; the long reads that join two
  single-copy segments, directly or through repeats, make bridges, which are scored and put in, the best first; what is then unbranched
  is made one segment; circles are written once round. The replicons are written as FASTA, to the standard output unless --fasta or -o says
  where, and with --gfa or -o the graph, as GFA. A summary goes to the error output.

  --mode M         conservative, normal (default) or bold: how good a bridge must be (25, 10, 1) and what is joined up afterwards
  --min-quality N  the lowest quality of a bridge, 0 to 100 (replaces the mode's)
  --min-reads N    the fewest long reads that make a bridge (default 1)
  --single-depth X the depth of a single copy (default: the median depth of the bases of the graph)
  --depth-tolerance X  how far from the single-copy depth a segment may be to start as one copy, as a share (default 0.25)
  --min-align N    the shortest alignment of a read to a segment that counts, in bases (default 300)
  --end-slack N    how far from the end of a segment an alignment may stop to join it to the next (default 200)
  --start-genes F  FASTA of nucleotide genes (dnaA, repA of the species): a circle is turned to start at the best match
  --overlap N      the overlap of the links, k-1, for a graph that has no link to say it
  -o DIR           write DIR/assembly.fasta and DIR/assembly.gfa (DIR is made if it is not there)
  --fasta FILE     where to write the replicons
  --gfa FILE       where to write the graph
  -q               no summary
  -h, --help       this text
";

struct Args {
    options: Options,
    graph: String,
    reads: Vec<String>,
    genes: Option<String>,
    overlap: usize,
    fasta: Option<String>,
    gfa: Option<String>,
    quiet: bool,
}

fn number_of<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} needs a number, not {value}."))
}

fn positive(flag: &str, value: &str) -> Result<f64, String> {
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
        graph: String::new(),
        reads: Vec::new(),
        genes: None,
        overlap: 0,
        fasta: None,
        gfa: None,
        quiet: false,
    };
    let mut files: Vec<String> = Vec::new();
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value."))
        };
        match flag.as_str() {
            "--mode" => {
                a.options.mode = match value("--mode")?.as_str() {
                    "conservative" => Mode::Conservative,
                    "normal" => Mode::Normal,
                    "bold" => Mode::Bold,
                    other => {
                        return Err(format!(
                            "--mode is conservative, normal or bold, not {other}."
                        ))
                    }
                }
            }
            "--min-quality" => {
                let q = positive("--min-quality", &value("--min-quality")?)?;
                if q > 100.0 {
                    return Err(format!("--min-quality is from 0 to 100, not {q}."));
                }
                a.options.min_quality = Some(q);
            }
            "--min-reads" => {
                a.options.min_reads = number_of("--min-reads", &value("--min-reads")?)?
            }
            "--single-depth" => {
                let d = positive("--single-depth", &value("--single-depth")?)?;
                if d == 0.0 {
                    return Err("--single-depth must be more than 0.".into());
                }
                a.options.copies.single_depth = Some(d);
            }
            "--depth-tolerance" => {
                a.options.copies.tolerance =
                    positive("--depth-tolerance", &value("--depth-tolerance")?)?
            }
            "--min-align" => {
                a.options.bridges.min_align = number_of("--min-align", &value("--min-align")?)?
            }
            "--end-slack" => {
                a.options.bridges.end_slack = number_of("--end-slack", &value("--end-slack")?)?
            }
            "--start-genes" => a.genes = Some(value("--start-genes")?),
            "--overlap" => a.overlap = number_of("--overlap", &value("--overlap")?)?,
            "-o" => {
                let dir = value("-o")?;
                a.fasta = Some(format!("{dir}/assembly.fasta"));
                a.gfa = Some(format!("{dir}/assembly.gfa"));
            }
            "--fasta" => a.fasta = Some(value("--fasta")?),
            "--gfa" => a.gfa = Some(value("--gfa")?),
            "-q" => a.quiet = true,
            other if other.starts_with('-') && other.len() > 1 => {
                return Err(format!(
                    "{other} is not an option (`helixloom hybrid --help`)."
                ))
            }
            file => files.push(file.to_string()),
        }
    }
    let mut files = files.into_iter();
    a.graph = files
        .next()
        .ok_or("Name the GFA file of the assembly graph.")?;
    a.reads = files.collect();
    if a.reads.is_empty() {
        return Err("Name the FASTA or FASTQ file of the long reads.".into());
    }
    if a.options.min_reads == 0 {
        return Err("--min-reads is 1 or more.".into());
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
    let text = std::fs::read_to_string(&a.graph)
        .map_err(|e| format!("Could not read {}: {e}", a.graph))?;
    let input = input::read(&text, a.overlap)?;
    let mut reads = Vec::new();
    for spec in &a.reads {
        reads.extend(read_input(spec)?);
    }
    let genes = match &a.genes {
        Some(path) => read_input(path)?,
        None => Vec::new(),
    };
    let outcome = assemble(input.graph, &reads, &a.options)?;
    let found = replicons(&outcome.graph, &genes, &RotateOptions::default());
    let text = fasta(&found);
    match &a.fasta {
        Some(path) => write(path, &text)?,
        None => print!("{text}"),
    }
    if let Some(path) = &a.gfa {
        let contigs = number(&outcome.graph);
        write(path, &gfa(&outcome.graph, &contigs))?;
    }
    if !a.quiet {
        let r = &outcome.report;
        if !input.has_depth {
            eprintln!("The graph has no depth (DP or KC) on some segments: they were taken as one copy deep.");
        }
        eprintln!(
            "{} segments: {} single copy, {} repeats ({} by depth alone); a single copy is {:.2} deep",
            r.segments, r.single_copy, r.repeats, r.by_depth, r.single_depth
        );
        eprintln!(
            "{} joins by {} long reads; {} bridges, {} good enough, {} made; {} segments left, {} circular, {} turned to a start gene",
            r.joins_seen,
            reads.len(),
            r.bridges_found,
            r.bridges_good_enough,
            r.bridges_made,
            r.segments_after,
            r.circular,
            found.iter().filter(|f| f.rotated).count()
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests;
