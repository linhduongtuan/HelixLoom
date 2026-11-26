//! `assemble` and `hybrid`: the de Bruijn assembler and the joining of its graph with long reads (`seq::debruijn`, `seq::hybrid`).

use crate::pipeline::flow::Value;
use crate::seq::debruijn::assemble::{assemble as run_assembly, Cutoff, Options, DEFAULT_KS};
use crate::seq::debruijn::output::{gfa, number};
use crate::seq::debruijn::pack::Pack;
use crate::seq::debruijn::reads::{pieces, Read};
use crate::seq::fastq::FastqRead;
use crate::seq::hybrid::assemble::{run as join, Options as JoinOptions};
use crate::seq::hybrid::input;
use crate::seq::hybrid::output::replicons;
use crate::seq::hybrid::rotate::Options as RotateOptions;
use crate::seq::hybrid::Mode;
use crate::seq::kmer::cli::read_input;
use std::path::Path;

/// `21,33,55` as numbers: each odd, from 3 to 63, and rising.
pub fn parse_ks(text: &str) -> Result<Vec<usize>, String> {
    let ks: Vec<usize> = text
        .split(',')
        .map(|k| {
            k.trim()
                .parse::<usize>()
                .map_err(|_| format!("ks is a list of numbers such as 21,33,55, not {text}"))
        })
        .collect::<Result<_, _>>()?;
    ks.iter().try_for_each(|&k| Pack::new(k).map(|_| ()))?;
    if ks.windows(2).any(|w| w[0] >= w[1]) {
        return Err(format!("the k must rise, none twice: {text}"));
    }
    Ok(ks)
}

/// Short reads into the graph of the de Bruijn assembler, as GFA.
pub fn assemble(
    reads: &[FastqRead],
    ks: &str,
    cutoff: &str,
    min_qual: u8,
) -> Result<(Value, String), String> {
    let list = if ks == "auto" {
        Vec::new()
    } else {
        parse_ks(ks)?
    };
    let shortest = list.first().copied().unwrap_or(DEFAULT_KS[0]);
    let as_reads: Vec<Read> = reads
        .iter()
        .map(|r| Read {
            seq: r.seq.clone(),
            qual: r.qual.clone(),
        })
        .collect();
    let used = pieces(&as_reads, min_qual, shortest);
    if used.is_empty() {
        return Err(format!(
            "no run of {shortest} bases or more is left in the reads once the unknown bases and those under quality {min_qual} are cut out"
        ));
    }
    let o = Options {
        ks: list,
        cutoff: match cutoff {
            "auto" => Cutoff::Auto,
            "off" => Cutoff::Off,
            n => Cutoff::Value(
                n.parse()
                    .map_err(|_| format!("cov-cutoff is not a number: {n}"))?,
            ),
        },
        ..Options::default()
    };
    let a = run_assembly(&used, &o)?;
    let contigs = number(&a.graph);
    let last = a.steps.last().map_or(0, |s| s.k);
    let note = format!(
        "{} segments, {} bases, k = {last}",
        a.graph.count(),
        a.graph.total_length()
    );
    Ok((Value::Gfa(gfa(&a.graph, &contigs)), note))
}

/// A graph joined with the long reads of `long` into replicons, as sequences named `NUMBER length=… depth=…x [circular=true]`.
pub fn hybrid(
    text: &str,
    long: Option<&Path>,
    mode: &str,
    min_quality: &str,
    single_depth: &str,
    overlap: usize,
) -> Result<(Value, String), String> {
    let long = long.ok_or("hybrid has no long reads to read")?;
    let graph = input::read(text, overlap)?.graph;
    let reads = read_input(&long.to_string_lossy())?;
    let mut o = JoinOptions {
        mode: match mode {
            "conservative" => Mode::Conservative,
            "bold" => Mode::Bold,
            _ => Mode::Normal,
        },
        ..JoinOptions::default()
    };
    o.min_quality = min_quality.parse().ok();
    o.copies.single_depth = single_depth.parse().ok();
    let done = join(graph, &reads, &o)?;
    let found = replicons(&done.graph, &[], &RotateOptions::default());
    let rows = found
        .iter()
        .map(|r| {
            (
                format!(
                    "{} length={} depth={:.2}x{}",
                    r.number,
                    r.seq.len(),
                    r.depth,
                    if r.circular { " circular=true" } else { "" }
                ),
                r.seq.clone(),
            )
        })
        .collect();
    let r = &done.report;
    let note = format!(
        "{} segments, {} bridges made from {} joins by {} long reads; {} replicons, {} circular",
        r.segments,
        r.bridges_made,
        r.joins_seen,
        reads.len(),
        found.len(),
        r.circular
    );
    Ok((Value::Sequences(rows), note))
}

#[cfg(test)]
mod tests;
