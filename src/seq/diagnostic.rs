//! Find Diagnostic Site (Enzymes menu): where an enzyme's site is one, two or three bases away from the DNA, so a primer
//! with those changes makes a site that tells a mutated copy from the original.
//!
//! The rules are the original's `do_diagnostic_site` (AppMain.tcl L36062 @0x1a3e5b) and its dialog
//! `add_daignostic_site_dialog` (L35993 @0x1a32b3): every enzyme but a list that needs methylated DNA, both strands of
//! its site (a palindrome once), every place the DNA differs from the site in at least one and at most the chosen number of
//! bases. `tests` checks the reports against the original's own, run on random sequences.

use super::dna::revcom;
use super::enzymes::Enzyme;
use super::sitesearch::{flat_pattern, header, pattern_covers, upper};

#[cfg(test)]
mod tests;

/// Enzymes the search leaves out because they only cut methylated DNA (the list at L36099).
pub const SKIPPED: [&str; 16] = [
    "AbaSI", "AoxI", "BisI", "BlsI", "DpnI", "FspEI", "GlaI", "GluI", "KroI", "LpnPI", "MalI",
    "MspJI", "MteI", "PcsI", "PkrI", "SgeI",
];
/// The most mismatches the dialog offers.
pub const MAX_MISMATCHES: usize = 3;
/// Most places listed before the search stops (HelixLoom's own limit; with three mismatches a four-base site fits almost
/// everywhere).
pub const MAX_FOUND: usize = 100_000;

/// One place where an enzyme's site, or its reverse complement, nearly is.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Found {
    /// Where it starts in the DNA searched, counting from 0.
    pub index: usize,
    pub enzyme: String,
    /// The site as read along this strand.
    pub pattern: Vec<u8>,
    pub mismatches: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Search {
    pub found: Vec<Found>,
    /// The search stopped at `MAX_FOUND`.
    pub truncated: bool,
}

/// The mismatches between `pattern` and the DNA at `at`, or more than `most` when it stops counting.
fn mismatches(pattern: &[u8], text: &[u8], at: usize, most: usize) -> usize {
    let mut count = 0;
    for (&p, &c) in pattern.iter().zip(&text[at..]) {
        if !pattern_covers(p, c) {
            count += 1;
            if count > most {
                break;
            }
        }
    }
    count
}

fn strands(flat: &[u8]) -> Vec<Vec<u8>> {
    let back = revcom(flat);
    if back == flat {
        vec![flat.to_vec()]
    } else {
        vec![flat.to_vec(), back]
    }
}

/// The places in `dna` where an enzyme's site differs in 1 to `most` bases, in order of position (equal positions keep the
/// order of the enzymes, forward strand first).
pub fn find(dna: &[u8], enzymes: &[&Enzyme], most: usize) -> Search {
    let text = dna.to_ascii_uppercase();
    let mut found = Vec::new();
    let mut truncated = false;
    'enzymes: for e in enzymes
        .iter()
        .filter(|e| !SKIPPED.contains(&e.name.as_str()))
    {
        for pattern in strands(&flat_pattern(e)) {
            let windows = (text.len() + 1).saturating_sub(pattern.len());
            for at in 0..windows {
                let wrong = mismatches(&pattern, &text, at, most);
                if wrong == 0 || wrong > most {
                    continue;
                }
                if found.len() >= MAX_FOUND {
                    truncated = true;
                    break 'enzymes;
                }
                found.push(Found {
                    index: at,
                    enzyme: e.name.clone(),
                    pattern: pattern.clone(),
                    mismatches: wrong,
                });
            }
        }
    }
    found.sort_by_key(|f| f.index);
    Search { found, truncated }
}

/// Each letter written under the DNA: the DNA's own letter where the site allows it, otherwise the site's letter, flagged.
pub fn marks(f: &Found, dna: &[u8]) -> Vec<(char, bool)> {
    let window = dna.get(f.index..).unwrap_or_default();
    f.pattern
        .iter()
        .enumerate()
        .map(|(i, &p)| match window.get(i) {
            Some(&c) if pattern_covers(p, c.to_ascii_uppercase()) => (char::from(c), false),
            _ => (upper(p), true),
        })
        .collect()
}

/// What the original's window shows for one place: the position from 1, the DNA, the DNA with the changes marked
/// between `<` and `>`, the enzyme and the site. `from` is where the searched stretch starts in the sequence (0-based).
pub fn render_found(f: &Found, dna: &[u8], from: usize) -> String {
    let end = (f.index + f.pattern.len()).min(dna.len());
    let shown: String = dna[f.index..end].iter().map(|&c| char::from(c)).collect();
    let mut text = format!("\t{}\t{shown}\n\t\t", f.index + from + 1);
    for (c, changed) in marks(f, dna) {
        if changed {
            text.push_str(&format!("<{c}>"));
        } else {
            text.push(c);
        }
    }
    text.push_str(&format!(
        "\n\t{}\t{}\n\n",
        f.enzyme,
        String::from_utf8_lossy(&f.pattern)
    ));
    text
}

/// The whole text of the original's Diagnostic Sites window. `shown` is the stretch as the heading writes it (the original
/// counts it from 0 here, unlike the positions below it, which count from 1); `start` is where `dna` starts in the
/// sequence, counting from 0.
pub fn report(
    (date, name): (&str, &str),
    shown: (usize, usize),
    start: usize,
    dna: &[u8],
    found: &[Found],
) -> String {
    let mut text = header(date, name, shown.0, shown.1, "Diagnostic Sites");
    for f in found {
        text.push_str(&render_found(f, dna, start));
    }
    if found.is_empty() {
        text.push_str("No sites Found");
    }
    text
}
