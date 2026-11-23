//! The ways to write a typing: tab- or comma-separated lines, the full table with a header, JSON, and the novel alleles as FASTA.

use super::call::{Call, Typing};
use super::md5;

/// One line for a genome: the file, the scheme, the sequence type, and each locus as `locus(call)`; `-` for the scheme and the type when
/// no scheme was chosen.
pub fn line(file: &str, t: Option<&Typing>, separator: char) -> String {
    let mut cells = vec![file.to_string()];
    match t {
        Some(t) => {
            cells.push(t.scheme.clone());
            cells.push(t.sequence_type.clone().unwrap_or_else(|| "-".into()));
            cells.extend(t.calls.iter().map(|(l, c)| format!("{l}({})", c.text())));
        }
        None => cells.extend(["-".to_string(), "-".to_string()]),
    }
    cells.join(&separator.to_string())
}

pub const FULL_HEADER: [&str; 6] = ["FILE", "SCHEME", "ST", "STATUS", "SCORE", "ALLELES"];

/// The line of the full table: FILE, SCHEME, ST, STATUS, SCORE, and the calls joined with `;`.
pub fn full_line(file: &str, t: Option<&Typing>, separator: char) -> String {
    let cells: Vec<String> = match t {
        Some(t) => vec![
            file.to_string(),
            t.scheme.clone(),
            t.sequence_type.clone().unwrap_or_else(|| "-".into()),
            t.status.word().to_string(),
            format!("{:.0}", t.score),
            t.calls
                .iter()
                .map(|(l, c)| format!("{l}({})", c.text()))
                .collect::<Vec<_>>()
                .join(";"),
        ],
        None => vec![
            file.to_string(),
            "-".into(),
            "-".into(),
            "NONE".into(),
            "0".into(),
            String::new(),
        ],
    };
    cells.join(&separator.to_string())
}

fn json_text(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// A JSON array with an object for each genome: `id`, `filename`, `scheme`, `sequence_type` and `alleles` (locus to call).
pub fn json(rows: &[(String, String, Option<Typing>)]) -> String {
    let objects: Vec<String> = rows
        .iter()
        .map(|(id, file, t)| {
            let (scheme, st, alleles) = match t {
                Some(t) => (
                    t.scheme.clone(),
                    t.sequence_type.clone().unwrap_or_else(|| "-".into()),
                    t.calls.iter().map(|(l, c)| format!("{}: {}", json_text(l), json_text(&c.text()))).collect::<Vec<_>>().join(", "),
                ),
                None => ("-".to_string(), "-".to_string(), String::new()),
            };
            format!(
                "  {{\"id\": {}, \"filename\": {}, \"scheme\": {}, \"sequence_type\": {}, \"alleles\": {{{alleles}}}}}",
                json_text(id),
                json_text(file),
                json_text(&scheme),
                json_text(&st)
            )
        })
        .collect();
    format!("[\n{}\n]\n", objects.join(",\n"))
}

/// The novel alleles of a typing as FASTA: `>scheme.locus-MD5 file`, the MD5 of the sequence in lower-case hexadecimal.
pub fn novel_fasta(file: &str, t: &Typing) -> String {
    let mut out = String::new();
    for (locus, call) in &t.calls {
        if let Call::Novel(found) = call {
            for (_, seq) in found {
                out += &format!(
                    ">{}.{locus}-{} {file}\n{}\n",
                    t.scheme,
                    md5::hex(seq),
                    String::from_utf8_lossy(seq)
                );
            }
        }
    }
    out
}

#[cfg(test)]
mod tests;
