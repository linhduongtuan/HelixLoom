//! The original's Translate window (`translation_window`, AppMain.tcl line 36460, file offset 0x1a8803):
//! the options it offers, their defaults, and the exact text it builds: a header with the range and the
//! molecular weight, then the protein lines (and optionally the DNA lines) with their numbering.
//! The original's own code was run to produce the cases this is tested against (`translation/golden.rs`).

use super::dna::{revcom, three_letter, translate};

/// Mass of one amino-acid residue, as the molecular-weight line adds them up (line 36616, 0x1aaa55).
/// `X` and the stop `*` count as 0.
const RESIDUE_MASSES: [(u8, f64); 20] = [
    (b'G', 57.05),
    (b'A', 71.09),
    (b'S', 87.08),
    (b'T', 101.11),
    (b'C', 103.15),
    (b'V', 99.14),
    (b'L', 113.16),
    (b'I', 113.16),
    (b'M', 131.19),
    (b'P', 97.12),
    (b'F', 147.18),
    (b'Y', 163.18),
    (b'W', 186.21),
    (b'D', 115.09),
    (b'E', 129.12),
    (b'N', 114.11),
    (b'Q', 128.14),
    (b'H', 137.14),
    (b'K', 128.17),
    (b'R', 156.19),
];
/// The water added to the chain's mass.
const WATER_MASS: f64 = 18.02;
/// With the line-width entry left empty the original uses this many residues per line.
const UNLIMITED_WIDTH: usize = 200000;
const CODON: usize = 3;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Letters {
    One,
    Three,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Numbers {
    None,
    Left,
    Right,
    Both,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum DnaLine {
    None,
    Above,
    Below,
}

/// The window's settings. `Default` holds the original's defaults (`translation_window_*` in
/// `Accessory Files/ApE_Defaults.txt`, lines 46-54): three-letter codes, no codon spacing, no reverse
/// complement, 20 residues per line, no line numbers, no DNA line.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Options {
    pub letters: Letters,
    pub codon_spacing: bool,
    pub reverse_complement: bool,
    /// Residues per line; 0 for no wrapping.
    pub width: usize,
    pub numbers: Numbers,
    pub dna_line: DnaLine,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            letters: Letters::Three,
            codon_spacing: false,
            reverse_complement: false,
            width: 20,
            numbers: Numbers::None,
            dna_line: DnaLine::None,
        }
    }
}

/// What is being translated.
#[derive(Clone, Copy, Debug)]
pub enum Source<'a> {
    /// A stretch of the sequence (all of it, or the selection) whose first base is number `first`.
    Range { first: usize },
    /// A CDS feature, read along its own strand and numbered from 1; its name heads the output.
    Feature { name: &'a str },
}

pub struct Output {
    /// From the range (or the feature name) down to the blank line after the molecular weight.
    pub header: Vec<String>,
    pub lines: Vec<String>,
}

/// A number as Tcl prints a float: the shortest text that reads back the same, with `.0` on a whole number.
fn tcl_number(x: f64) -> String {
    let text = x.to_string();
    if text.contains('.') || text.contains('e') {
        text
    } else {
        format!("{text}.0")
    }
}

/// Mass of the residues plus water, rounded to 0.01 (`round($mw*100)/100.0`).
pub fn molecular_weight(residues: &[u8]) -> f64 {
    let chain: f64 = residues
        .iter()
        .map(|r| {
            RESIDUE_MASSES
                .iter()
                .find(|(a, _)| *a == r.to_ascii_uppercase())
                .map_or(0.0, |&(_, mass)| mass)
        })
        .sum();
    ((chain + WATER_MASS) * 100.0).round() / 100.0
}

/// `text[lo..hi]` with both ends held inside the text (all of it is ASCII).
fn part(text: &str, lo: usize, hi: usize) -> &str {
    let len = text.len();
    &text[lo.min(len)..hi.min(len)]
}

/// Build what the original's Translate window shows for `bases` (in document order; a feature's bases
/// already along its own strand).
pub fn build(bases: &[u8], table: &str, source: Source, o: &Options) -> Output {
    let (first, last) = match source {
        Source::Range { first } => (first, (first + bases.len()).saturating_sub(1)),
        Source::Feature { .. } => (1, bases.len()),
    };
    let reverse = o.reverse_complement && matches!(source, Source::Range { .. });
    let extra = bases.len() % CODON;
    // the first and last base number of the translated stretch, and the bases read
    let (start, end, bases) = match (source, reverse) {
        (Source::Range { .. }, true) => (first + extra, last, revcom(bases)),
        (Source::Range { .. }, false) => (first, last.saturating_sub(extra), bases.to_vec()),
        (Source::Feature { .. }, _) => (first, last, bases.to_vec()),
    };
    let whole = &bases[..bases.len() / CODON * CODON];

    // the molecular weight covers the protein up to the first stop; what follows is counted
    let protein = translate(whole, table);
    let stop = protein.iter().position(|&a| a == b'*');
    let chain = &protein[..stop.unwrap_or(protein.len())];
    let after_stop = stop.map_or(0, |s| protein.len() - s - 1);
    let mut header = vec![match source {
        Source::Range { .. } if reverse => format!("From {end} to {start}."),
        Source::Range { .. } => format!("From {start} to {end}."),
        Source::Feature { name } => format!("{name} "),
    }];
    header.push(format!(
        "Translation {} a.a. MW={}",
        chain.len(),
        tcl_number(molecular_weight(chain))
    ));
    header.push(if after_stop > 0 {
        format!("???({after_stop} extra codons after stop)")
    } else {
        String::new()
    });

    // a one-letter protein under a DNA line is spread so each letter sits under its codon's middle base
    let one = o.letters == Letters::One;
    let dna_shown = o.dna_line != DnaLine::None;
    let spacing = match (o.codon_spacing, one && dna_shown) {
        (true, true) => 3,
        (true, false) => 1,
        (false, true) => 2,
        (false, false) => 0,
    };
    let residues: Vec<String> = protein
        .iter()
        .map(|&a| {
            if one {
                char::from(a).to_string()
            } else {
                three_letter(a).to_string()
            }
        })
        .collect();
    let mut protein_text = residues.join(&" ".repeat(spacing));
    let mut dna_text = String::from_utf8_lossy(whole).into_owned();
    if o.codon_spacing {
        let codons: Vec<String> = whole
            .chunks(CODON)
            .map(|c| String::from_utf8_lossy(c).into_owned())
            .collect();
        dna_text = format!("{} ", codons.join(" "));
        protein_text.push(' ');
    }
    if one && dna_shown {
        protein_text = format!(" {protein_text} ");
    }
    // characters taken by one residue in both lines
    let unit = if one { 1 } else { CODON } + spacing;

    let width = if o.width == 0 {
        UNLIMITED_WIDTH
    } else {
        o.width
    };
    let (left, right) = (
        matches!(o.numbers, Numbers::Left | Numbers::Both),
        matches!(o.numbers, Numbers::Right | Numbers::Both),
    );
    let numbered = |before: usize, text: &str, after: usize| {
        let mut row = String::new();
        if left {
            row.push_str(&format!("{before:>7} "));
        }
        row.push_str(text);
        if right {
            row.push_str(&format!(" {after:<7}"));
        }
        row
    };
    let mut lines = Vec::new();
    let mut x = 0;
    while x * unit < protein_text.len() {
        let (lo, hi) = (x * unit, (width + x) * unit);
        let (pline, dline) = (part(&protein_text, lo, hi), part(&dna_text, lo, hi));
        let line_start = if reverse {
            end - CODON * x
        } else {
            start + CODON * x
        };
        let line_end = match (dline.len() == width * unit, reverse) {
            (true, true) => (end + 1).saturating_sub(CODON * (x + width)),
            (true, false) => start + CODON * (x + width) - 1,
            (false, true) => start,
            (false, false) => end,
        };
        let dna_row = numbered(line_start, dline, line_end);
        if o.dna_line == DnaLine::Above {
            lines.push(dna_row.clone());
        }
        lines.push(numbered(x + 1, pline, x + pline.len() / unit));
        if o.dna_line == DnaLine::Below {
            lines.push(dna_row);
        }
        if dna_shown {
            lines.push(String::new());
        }
        x += width;
    }
    Output { header, lines }
}

#[cfg(test)]
mod golden;
#[cfg(test)]
mod tests;
