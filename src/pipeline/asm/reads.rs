//! `trim` and `filter`: reads in, reads out.

use super::Cut;
use crate::pipeline::flow::Value;
use crate::seq::fastq::{clean, Cleaning, FastqRead};
use crate::seq::readfilter::{filter_long, Filter};

/// Short reads cleaned by fastp's rules (`seq::fastq`): the adapter cut from the 3' end, the end of low quality cut off, and the reads that are
/// then too short, or have too many unknown or poor bases, dropped.
pub fn trim(
    reads: &[FastqRead],
    adapter: &str,
    quality: u8,
    min_length: usize,
    cut: Cut,
) -> Result<(Value, String), String> {
    let c = Cleaning {
        adapter: if adapter == "none" {
            Vec::new()
        } else {
            adapter.as_bytes().to_vec()
        },
        cut_tail: cut == Cut::Tail,
        cut_right: cut == Cut::Right,
        window_quality: usize::from(quality),
        min_length,
        ..Cleaning::default()
    };
    let (kept, report) = clean(reads, &c);
    if kept.is_empty() {
        return Err(format!(
            "no read is left of {} after trimming ({} too short, {} with too many poor bases, {} with too many N)",
            reads.len(),
            report.too_short,
            report.too_many_low_quality,
            report.too_many_n
        ));
    }
    let note = format!(
        "{} of {} reads kept; {} had the adapter cut, {} too short, {} too poor, {} too many N",
        kept.len(),
        reads.len(),
        report.adapter_trimmed,
        report.too_short,
        report.too_many_low_quality,
        report.too_many_n
    );
    Ok((Value::Reads(kept), note))
}

/// Long reads selected as the multiplex-MinION paper does (`seq::readfilter`).
pub fn filter(
    reads: &[FastqRead],
    min_length: usize,
    target_bases: u64,
    window: usize,
) -> Result<(Value, String), String> {
    let f = filter_long(
        reads,
        &Filter {
            min_length,
            target_bases,
            window,
        },
    );
    if f.reads.is_empty() {
        return Err(format!(
            "no read of {} is {min_length} bases or longer",
            reads.len()
        ));
    }
    let note = format!(
        "{} of {} reads kept, {} to {} bases; {} too short, {} dropped for quality",
        f.reads.len(),
        reads.len(),
        f.bases_in,
        f.bases_out,
        f.too_short,
        f.low_quality
    );
    Ok((Value::Reads(f.reads), note))
}

#[cfg(test)]
mod tests;
