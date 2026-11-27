//! The count table and the summary, laid out as featureCounts writes them.

use super::assign::{Index, Status, Tally};

/// One counted file.
pub struct Column<'a> {
    pub name: &'a str,
    pub tally: &'a Tally,
}

fn count_text(value: f64, fraction: bool) -> String {
    if fraction {
        format!("{value:.2}")
    } else {
        format!("{value:.0}")
    }
}

/// The table: a comment line, then `Geneid Chr Start End Strand Length` and a column of counts for each file. A meta-feature lists the
/// sequence, start (from 1), end and strand of each of its features, joined by `;`, and `Length` is the bases of its features counted once;
/// at feature level each row is one feature.
pub fn table(
    index: &Index<'_>,
    columns: &[Column<'_>],
    feature_level: bool,
    fraction: bool,
) -> String {
    let annotation = index.annotation;
    let mut out = format!(
        "# Program: helixloom {} count (after featureCounts)\nGeneid\tChr\tStart\tEnd\tStrand\tLength",
        env!("CARGO_PKG_VERSION")
    );
    for c in columns {
        out.push('\t');
        out.push_str(c.name);
    }
    out.push('\n');
    let rows: Vec<(usize, Vec<usize>)> = if feature_level {
        (0..annotation.features.len())
            .map(|i| (annotation.features[i].meta, vec![i]))
            .collect()
    } else {
        let mut parts = vec![Vec::new(); annotation.names.len()];
        for (i, f) in annotation.features.iter().enumerate() {
            parts[f.meta].push(i);
        }
        parts.into_iter().enumerate().collect()
    };
    for (target, (meta, features)) in rows.into_iter().enumerate() {
        let joined = |of: &dyn Fn(usize) -> String| {
            features
                .iter()
                .map(|&i| of(i))
                .collect::<Vec<_>>()
                .join(";")
        };
        out.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            annotation.names[meta],
            joined(&|i| annotation.features[i].seqname.clone()),
            joined(&|i| (annotation.features[i].start + 1).to_string()),
            joined(&|i| annotation.features[i].end.to_string()),
            joined(&|i| annotation.features[i].strand.symbol().to_string()),
            index.lengths[target],
        ));
        for c in columns {
            out.push('\t');
            out.push_str(&count_text(c.tally.counts[target], fraction));
        }
        out.push('\n');
    }
    out
}

/// `Status` and a column for each file: how many alignments had each status.
pub fn summary(columns: &[Column<'_>]) -> String {
    let mut out = String::from("Status");
    for c in columns {
        out.push('\t');
        out.push_str(c.name);
    }
    out.push('\n');
    for status in Status::ALL {
        out.push_str(status.name());
        for c in columns {
            out.push_str(&format!("\t{}", c.tally.of(status)));
        }
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests;
