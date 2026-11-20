//! The Clustal alignment format (`.aln`), the one ClustalW and ClustalX write and that most alignment programs can read and write:
//! a first line that says CLUSTAL, then blocks of 60 columns, each row a name and its residues (and the count of residues so far),
//! and under each block a line of conservation marks.
//!
//! The marks are ClustalW's own: `*` where every row has the same residue, `:` where all the residues belong to one of the nine
//! *strong* groups (STA, NEQK, NHQK, NDEQ, QHRK, MILV, MILF, HY, FYW), `.` where they belong to one of the eleven *weak* groups (CSA,
//! ATV, SAG, STNK, STPA, SGND, SNDEQK, NDEQHK, NEQHRK, FVLIM, HFY), and a blank for the rest and for any column with a gap. The
//! groups are those of the PAM 250 matrix as ClustalW documents them. Nucleotides get `*` only.

use crate::seq::alnstats::is_protein;

/// Columns in a block.
pub const BLOCK: usize = 60;
/// A name is padded to at least this wide, and a space is left after it.
const NAME_WIDTH: usize = 16;
const STRONG: [&str; 9] = [
    "STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW",
];
const WEAK: [&str; 11] = [
    "CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM", "HFY",
];

fn is_gap(c: u8) -> bool {
    matches!(c, b'-' | b'.')
}

/// The first words of the headers of the programs that write this format (ClustalW and Omega, MAFFT, MUSCLE version 3, ProbCons,
/// MSAProbs, Kalign), in capitals.
const HEADERS: [&str; 5] = ["CLUSTAL", "MUSCLE", "PROBCONS", "MSAPROBS", "KALIGN"];

/// Whether a text is a Clustal alignment: its first line that is not blank begins with one of those words and says more after it, as
/// every program's header does; or the keyword stands alone and the line under it is a row (a name and its residues). A protein that
/// is only the letters of one of those words is not taken for one.
pub fn is_clustal(text: &str) -> bool {
    let text = text.trim_start_matches('\u{feff}');
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let Some(first) = lines.next() else {
        return false;
    };
    let mut words = first.split_whitespace();
    let keyword = words.next().map(str::to_ascii_uppercase);
    if !keyword.as_deref().is_some_and(|k| HEADERS.contains(&k)) {
        return false;
    }
    words.next().is_some()
        || lines
            .next()
            .is_some_and(|l| l.split_whitespace().count() >= 2)
}

/// The mark of one column (upper case letters, at least one row).
fn mark(column: &[u8], protein: bool) -> char {
    if column.iter().any(|&c| is_gap(c)) {
        return ' ';
    }
    if column.iter().all(|&c| c == column[0]) {
        return '*';
    }
    if !protein {
        return ' ';
    }
    let within = |groups: &[&str]| {
        groups
            .iter()
            .any(|g| column.iter().all(|c| g.as_bytes().contains(c)))
    };
    if within(&STRONG) {
        ':'
    } else if within(&WEAK) {
        '.'
    } else {
        ' '
    }
}

/// One mark for each column of the rows (a row that is shorter than the others counts as gaps).
pub fn conservation(rows: &[(String, Vec<u8>)]) -> String {
    let columns = rows.iter().map(|r| r.1.len()).max().unwrap_or(0);
    let letters: Vec<Vec<u8>> = rows.iter().map(|r| r.1.clone()).collect();
    let protein = is_protein(&letters);
    (0..columns)
        .map(|k| {
            let column: Vec<u8> = rows
                .iter()
                .map(|r| r.1.get(k).map_or(b'-', u8::to_ascii_uppercase))
                .collect();
            mark(&column, protein)
        })
        .collect()
}

/// The common length of the rows, or why they have none.
fn common_length(rows: &[(String, Vec<u8>)]) -> Result<usize, String> {
    let Some(first) = rows.first() else {
        return Err("There are no sequences to write.".into());
    };
    if let Some(odd) = rows.iter().find(|r| r.1.len() != first.1.len()) {
        return Err(format!(
            "The rows are not of one length ({} has {} columns, {} has {}): they are not an alignment.",
            first.0,
            first.1.len(),
            odd.0,
            odd.1.len()
        ));
    }
    match first.1.len() {
        0 => Err("The rows are empty.".into()),
        n => Ok(n),
    }
}

/// The alignment as a Clustal text. A name cannot hold a space in this format, so each run of white space in one becomes an
/// underscore.
pub fn write(rows: &[(String, Vec<u8>)]) -> Result<String, String> {
    let columns = common_length(rows)?;
    let names: Vec<String> = rows
        .iter()
        .map(|r| r.0.split_whitespace().collect::<Vec<_>>().join("_"))
        .map(|n| {
            if n.is_empty() {
                "sequence".to_string()
            } else {
                n
            }
        })
        .collect();
    let width = names
        .iter()
        .map(|n| n.chars().count())
        .max()
        .unwrap_or(0)
        .max(NAME_WIDTH)
        + 1;
    let marks = conservation(rows);
    let mut out = String::from("CLUSTAL format alignment by HelixLoom\n\n");
    let mut residues = vec![0usize; rows.len()];
    for start in (0..columns).step_by(BLOCK) {
        let end = (start + BLOCK).min(columns);
        out.push('\n');
        for (i, (_, row)) in rows.iter().enumerate() {
            let piece = &row[start..end];
            residues[i] += piece.iter().filter(|&&c| !is_gap(c)).count();
            out += &format!(
                "{:<width$}{} {}\n",
                names[i],
                String::from_utf8_lossy(piece),
                residues[i]
            );
        }
        out += &format!("{}{}\n", " ".repeat(width), &marks[start..end]);
    }
    Ok(out)
}

/// The rows of a Clustal text: the rows of its first block name the sequences, and each later block adds to them in the same order.
/// The counts at the ends of the lines and the lines of marks are left out. The rows must come out of one length.
pub fn read(text: &str) -> Result<Vec<(String, Vec<u8>)>, String> {
    if !is_clustal(text) {
        return Err("This is not a Clustal alignment: its first line should say CLUSTAL.".into());
    }
    let text = text.trim_start_matches('\u{feff}');
    let mut rows: Vec<(String, Vec<u8>)> = Vec::new();
    let (mut place, mut first_block_over, mut header_seen) = (0usize, false, false);
    for (n, line) in text.lines().enumerate().map(|(n, l)| (n + 1, l)) {
        if !header_seen {
            header_seen = !line.trim().is_empty();
            continue;
        }
        // a line of nothing but marks and blanks ends a block, indented as the programs write it or not (an editor may have trimmed it)
        if line.chars().all(|c| " \t*:.".contains(c)) {
            first_block_over |= !rows.is_empty();
            place = 0;
            continue;
        }
        if line.starts_with([' ', '\t']) {
            return Err(format!(
                "Line {n}: a line that begins with a space should be the marks of conservation."
            ));
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        let (name, residues, count) = match parts.as_slice() {
            [name, residues] => (*name, *residues, None),
            [name, residues, count] if count.bytes().all(|c| c.is_ascii_digit()) => {
                (*name, *residues, Some(*count))
            }
            _ => return Err(format!("Line {n}: expected a name and its residues.")),
        };
        if !residues
            .bytes()
            .all(|c| c.is_ascii_alphabetic() || matches!(c, b'-' | b'.' | b'*'))
        {
            return Err(format!(
                "Line {n}: \"{residues}\" has a character that is neither a residue nor a gap."
            ));
        }
        if let Some(row) = rows.get_mut(place) {
            if row.0 != name {
                return Err(format!(
                    "Line {n}: {name} is not the sequence that the first block has at this place ({}).",
                    row.0
                ));
            }
            row.1.extend(residues.bytes());
        } else if first_block_over {
            return Err(format!(
                "Line {n}: {name} is not one of the sequences of the first block."
            ));
        } else {
            rows.push((name.to_string(), residues.bytes().collect()));
        }
        if let (Some(count), Some(row)) = (count, rows.get(place)) {
            // the count at the end of the line is the residues so far (or the columns, as some programs count): a block that is
            // missing leaves it ahead of them
            let residues_so_far = row.1.iter().filter(|&&c| !is_gap(c)).count();
            if count.parse() != Ok(residues_so_far) && count.parse() != Ok(row.1.len()) {
                return Err(format!(
                    "Line {n}: the count {count} at the end of the line is not the {residues_so_far} residues of {name} so far: a block of the alignment may be missing."
                ));
            }
        }
        place += 1;
    }
    if rows.is_empty() {
        return Err("The Clustal alignment holds no sequences.".into());
    }
    common_length(&rows).map_err(|why| why.replace("to write", "in it"))?;
    Ok(rows)
}

#[cfg(test)]
mod tests;
