//! What a reader of an alignment wants beside the letters: for each column the consensus symbol and how
//! conserved it is, and a summary of the whole alignment.

const GAP: u8 = b'-';

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Column {
    /// The commonest symbol (upper case), or a gap when gaps are at least half of the rows.
    pub consensus: u8,
    /// The share of the rows that carry the consensus symbol, 0 to 1 (gaps never count as carrying it).
    pub conservation: f32,
    pub gaps: usize,
}

fn is_gap(c: u8) -> bool {
    matches!(c, b'-' | b'.')
}

/// The consensus and conservation of every column of `rows` (rows of one length, gaps written `-`).
pub fn columns(rows: &[Vec<u8>]) -> Vec<Column> {
    let width = rows.iter().map(Vec::len).max().unwrap_or(0);
    (0..width)
        .map(|k| {
            let mut counts = [0usize; 256];
            let mut gaps = 0;
            for row in rows {
                match row.get(k) {
                    Some(&c) if !is_gap(c) => counts[usize::from(c.to_ascii_uppercase())] += 1,
                    _ => gaps += 1,
                }
            }
            // the commonest symbol; on a tie the one that sorts first, so the answer does not depend on row order
            let (symbol, top) = counts.iter().enumerate().fold((GAP, 0), |best, (c, &n)| {
                if n > best.1 {
                    (c as u8, n)
                } else {
                    best
                }
            });
            let consensus = if 2 * gaps >= rows.len() { GAP } else { symbol };
            Column {
                consensus,
                conservation: if rows.is_empty() {
                    0.0
                } else {
                    top as f32 / rows.len() as f32
                },
                gaps,
            }
        })
        .collect()
}

/// Whether the rows hold amino acids rather than nucleotides: more than a tenth of the letters are not
/// nucleotide codes.
pub fn is_protein(rows: &[Vec<u8>]) -> bool {
    let (mut letters, mut other) = (0usize, 0usize);
    for &c in rows.iter().flatten() {
        if c.is_ascii_alphabetic() {
            letters += 1;
            if !b"ACGTUNRYKMSWBDHV".contains(&c.to_ascii_uppercase()) {
                other += 1;
            }
        }
    }
    letters > 0 && other * 10 > letters
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Summary {
    pub rows: usize,
    pub columns: usize,
    /// The mean conservation of the columns, in percent.
    pub identity: f64,
    /// Cells that are gaps, in percent.
    pub gaps: f64,
    /// Columns with a symbol shared by every row and no gap.
    pub identical: usize,
    /// Columns where at least two different symbols occur.
    pub variable: usize,
    /// Columns where at least two symbols each occur in at least two rows (the sites a parsimony tree learns from).
    pub informative: usize,
}

pub fn summary(rows: &[Vec<u8>]) -> Summary {
    let cols = columns(rows);
    let width = cols.len();
    let (mut identical, mut variable, mut informative) = (0, 0, 0);
    for (k, column) in cols.iter().enumerate() {
        let mut counts = std::collections::HashMap::new();
        for row in rows {
            if let Some(&c) = row.get(k).filter(|&&c| !is_gap(c)) {
                *counts.entry(c.to_ascii_uppercase()).or_insert(0usize) += 1;
            }
        }
        identical += usize::from(counts.len() == 1 && column.gaps == 0);
        variable += usize::from(counts.len() > 1);
        informative += usize::from(counts.values().filter(|&&n| n >= 2).count() >= 2);
    }
    let cells = (rows.len() * width).max(1) as f64;
    Summary {
        rows: rows.len(),
        columns: width,
        identity: 100.0 * cols.iter().map(|c| f64::from(c.conservation)).sum::<f64>()
            / width.max(1) as f64,
        gaps: 100.0 * cols.iter().map(|c| c.gaps).sum::<usize>() as f64 / cells,
        identical,
        variable,
        informative,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(text: &[&str]) -> Vec<Vec<u8>> {
        text.iter().map(|r| r.as_bytes().to_vec()).collect()
    }

    #[test]
    fn the_consensus_is_the_commonest_symbol_and_a_gap_where_gaps_are_half() {
        let a = rows(&["ACGT-", "ACGA-", "ACTA-", "ACTA-"]);
        let cols = columns(&a);
        let consensus: Vec<u8> = cols.iter().map(|c| c.consensus).collect();
        // the third column is G, G, T, T: a tie goes to the symbol that sorts first
        assert_eq!(consensus, b"ACGA-");
        assert_eq!(cols[0].conservation, 1.0);
        assert_eq!(cols[2].conservation, 0.5);
        // a column that is half gaps is a gap in the consensus, and a gap never counts as conserved
        let b = columns(&rows(&["A", "A", "-", "-"]));
        assert_eq!(
            (b[0].consensus, b[0].gaps, b[0].conservation),
            (b'-', 2, 0.5)
        );
    }

    #[test]
    fn a_tie_does_not_depend_on_the_order_of_the_rows() {
        let one = columns(&rows(&["A", "C"]));
        let two = columns(&rows(&["C", "A"]));
        assert_eq!(one[0].consensus, two[0].consensus);
        assert_eq!(one[0].consensus, b'A');
    }

    #[test]
    fn case_and_dots_are_read_like_upper_case_and_dashes() {
        let cols = columns(&rows(&["aC.", "AC-", "Ac-"]));
        assert_eq!(cols[0].conservation, 1.0);
        assert_eq!(cols[1].consensus, b'C');
        assert_eq!(cols[2].gaps, 3);
    }

    #[test]
    fn the_summary_counts_identical_variable_and_informative_columns() {
        // columns: identical, identical, one change (variable only), two symbols twice each (informative), a gap
        let s = summary(&rows(&["AAAAG", "AAAAG", "AACC-", "AACC-"]));
        assert_eq!((s.rows, s.columns), (4, 5));
        assert_eq!(s.identical, 2);
        assert_eq!(s.variable, 2);
        assert_eq!(s.informative, 2);
        assert!((s.gaps - 10.0).abs() < 1e-9, "{}", s.gaps);
        let single_change = summary(&rows(&["AAA", "AAA", "AAC"]));
        assert_eq!((single_change.variable, single_change.informative), (1, 0));
    }

    #[test]
    fn identical_rows_are_a_hundred_percent_and_nothing_is_not_a_number() {
        assert_eq!(summary(&rows(&["ACGT", "ACGT"])).identity, 100.0);
        let empty = summary(&[]);
        assert_eq!((empty.rows, empty.columns), (0, 0));
        assert!(empty.identity.is_finite() && empty.gaps.is_finite());
        assert!(columns(&[]).is_empty());
    }

    #[test]
    fn protein_is_told_from_dna() {
        assert!(!is_protein(&rows(&["ACGTNACGT", "ACG-TRYAC"])));
        assert!(is_protein(&rows(&["MKTAYIAKQR", "MKT-YIAKQR"])));
        assert!(!is_protein(&rows(&["----"])));
    }
}
