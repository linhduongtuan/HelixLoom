//! `concatemer`, `circularize` and `rotate`: sequences in, sequences out, the finding written as a word after the name.

use super::{has_tag, tag};
use crate::pipeline::flow::Value;
use crate::pipeline::steps::Named;
use crate::seq::circular::{end_overlap, find_concatemer, Concatemer, Overlap};
use crate::seq::hybrid::rotate::{rotate as turn, Options as RotateOptions};
use crate::seq::kmer::cli::read_input;
use std::path::Path;

/// A contig that is several copies of one unit is cut to the first copy, and says `concatemer=N`.
pub fn concatemer(
    rows: &[Named],
    min_identity: u8,
    min_unit: usize,
) -> Result<(Value, String), String> {
    let o = Concatemer {
        min_unit,
        min_identity: f64::from(min_identity) / 100.0,
        ..Concatemer::default()
    };
    let mut collapsed = 0;
    let out: Vec<Named> = rows
        .iter()
        .map(|(name, seq)| match find_concatemer(seq, &o) {
            Some(f) => {
                collapsed += 1;
                (
                    tag(name, &format!("concatemer={}", f.copies)),
                    seq[..f.unit].to_vec(),
                )
            }
            None => (name.clone(), seq.clone()),
        })
        .collect();
    Ok((
        Value::Sequences(out),
        format!("{collapsed} of {} contigs were concatemers", rows.len()),
    ))
}

/// A contig whose end overlaps its start has the overlap cut off, and says `circular=true`; one that says it already is left alone.
pub fn circularize(
    rows: &[Named],
    min_overlap: usize,
    min_identity: u8,
) -> Result<(Value, String), String> {
    let o = Overlap {
        min_length: min_overlap,
        min_identity: f64::from(min_identity) / 100.0,
        ..Overlap::default()
    };
    let mut made = 0;
    let out: Vec<Named> = rows
        .iter()
        .map(|(name, seq)| {
            if has_tag(name, "circular=true") {
                return (name.clone(), seq.clone());
            }
            match end_overlap(seq, &o) {
                Some(l) => {
                    made += 1;
                    (tag(name, "circular=true"), seq[..seq.len() - l].to_vec())
                }
                None => (name.clone(), seq.clone()),
            }
        })
        .collect();
    Ok((
        Value::Sequences(out),
        format!("{made} of {} contigs are circles", rows.len()),
    ))
}

/// A circle (one that says `circular=true`) is turned to start at the best match of the genes in `genes`, and says `rotated=true`.
pub fn rotate(rows: &[Named], genes: Option<&Path>) -> Result<(Value, String), String> {
    let genes = read_input(
        &genes
            .ok_or("rotate has no genes to read")?
            .to_string_lossy(),
    )?;
    let mut turned = 0;
    let out: Vec<Named> = rows
        .iter()
        .map(|(name, seq)| {
            if !has_tag(name, "circular=true") {
                return (name.clone(), seq.clone());
            }
            match turn(seq, &genes, &RotateOptions::default()) {
                Some(s) => {
                    turned += 1;
                    (tag(name, "rotated=true"), s)
                }
                None => (name.clone(), seq.clone()),
            }
        })
        .collect();
    let circles = rows
        .iter()
        .filter(|r| has_tag(&r.0, "circular=true"))
        .count();
    Ok((
        Value::Sequences(out),
        format!("{turned} of {circles} circles turned to a start gene"),
    ))
}

#[cfg(test)]
mod tests;
