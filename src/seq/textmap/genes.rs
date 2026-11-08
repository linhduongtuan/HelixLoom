//! The translations a Text Map writes under the DNA, and the "Uppercase as Genes" line: `translate` (AppMain.tcl L14518
//! @0xa51b0, with `gen_translation_map` L14477) and `uppercase_cdna` (L20180).

use super::tcl::{self, chars, range, spaces};

pub use crate::seq::dna::STANDARD_CODE;

/// 0 to 3 for t (or u), c, a, g in either case; none for any other letter.
fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_lowercase() {
        b't' | b'u' => Some(0),
        b'c' => Some(1),
        b'a' => Some(2),
        b'g' => Some(3),
        _ => None,
    }
}

/// The three-letter name of an amino acid (the original's `string map`, L14497).
fn three_letters(aa: u8) -> String {
    let name = match aa {
        b'X' => "Xxx",
        b'F' => "Phe",
        b'L' => "Leu",
        b'S' => "Ser",
        b'Y' => "Tyr",
        b'C' => "Cys",
        b'*' => "End",
        b'W' => "Trp",
        b'P' => "Pro",
        b'H' => "His",
        b'Q' => "Gln",
        b'R' => "Arg",
        b'I' => "Ile",
        b'M' => "Met",
        b'T' => "Thr",
        b'N' => "Asn",
        b'K' => "Lys",
        b'V' => "Val",
        b'A' => "Ala",
        b'D' => "Asp",
        b'E' => "Glu",
        b'G' => "Gly",
        other => return char::from(other).to_string(),
    };
    name.to_string()
}

/// `translate dna chars direction spacing`: every whole codon of `dna` read as it lies (`forward`) or as the codon of the
/// other strand that pairs with it, one letter or three, with `spacing` spaces between the amino acids. A codon with a letter
/// that is not a, c, g, t or u is `X`. The reverse reading does not turn the order round: each codon stays where it is.
pub(super) fn translate(
    dna: &[u8],
    code: &[u8],
    three: bool,
    forward: bool,
    spacing: usize,
) -> String {
    let code = if code.len() == 64 {
        code
    } else {
        STANDARD_CODE.as_bytes()
    };
    let words: Vec<String> = dna[..dna.len() / 3 * 3]
        .chunks(3)
        .map(|codon| {
            let bases: Option<Vec<usize>> = codon.iter().map(|&b| base_index(b)).collect();
            let aa = match bases.as_deref() {
                None => b'X',
                Some(&[a, b, c]) if forward => code[a * 16 + b * 4 + c],
                // the codon of the other strand reads the pairs of these bases backwards: t and a, c and g trade places
                Some(&[a, b, c]) => code[((c + 2) % 4) * 16 + ((b + 2) % 4) * 4 + (a + 2) % 4],
                Some(_) => b'X',
            };
            if three {
                three_letters(aa)
            } else {
                char::from(aa).to_string()
            }
        })
        .collect();
    words.join(&" ".repeat(spacing))
}

fn small(b: u8) -> bool {
    (b'a'..=b'y').contains(&b)
}

fn capital(b: u8) -> bool {
    (b'A'..=b'Y').contains(&b)
}

/// The stretches of small letters that `uppercase_cdna` takes out, as (first, last) places: the first run of small letters
/// that a capital follows (`([a-y]+)[A-Y]`, searched once), the introns (`[A-Y](g[tc][a-y]*?ag)[A-Y]`, all of them: the capital
/// after one is used up, so two introns need two capitals between them), and the small letters after a stop codon
/// (`(?:TGA|TAA|TAG)([a-y]+)`).
fn stretches(dna: &[u8]) -> Vec<(i64, i64)> {
    let n = dna.len();
    let mut found = Vec::new();
    let mut i = 0;
    while i < n {
        if !small(dna[i]) {
            i += 1;
            continue;
        }
        let end = (i..n).find(|&k| !small(dna[k])).unwrap_or(n);
        if end < n && capital(dna[end]) {
            found.push((i as i64, end as i64 - 1));
            break;
        }
        i = end;
    }
    let mut i = 0;
    while i < n {
        match intron_at(dna, i) {
            Some((first, last, next)) => {
                found.push((first as i64, last as i64));
                i = next;
            }
            None => i += 1,
        }
    }
    let mut i = 0;
    while i + 3 < n {
        let stop = matches!(&dna[i..i + 3], b"TGA" | b"TAA" | b"TAG");
        let end = (i + 3..n).find(|&k| !small(dna[k])).unwrap_or(n);
        if stop && end > i + 3 {
            found.push(((i + 3) as i64, end as i64 - 1));
            i = end;
        } else {
            i += 1;
        }
    }
    found
}

/// An intron that follows the capital at `i`: its first and last place, and where the search goes on (after the capital that
/// ends it). The body is as short as it can be.
fn intron_at(dna: &[u8], i: usize) -> Option<(usize, usize, usize)> {
    let n = dna.len();
    if !(capital(dna[i]) && i + 2 < n && dna[i + 1] == b'g' && matches!(dna[i + 2], b't' | b'c')) {
        return None;
    }
    let mut p = i + 3;
    while p + 2 < n {
        if dna[p] == b'a' && dna[p + 1] == b'g' && capital(dna[p + 2]) {
            return Some((i + 1, p + 1, p + 3));
        }
        if !small(dna[p]) {
            return None;
        }
        p += 1;
    }
    None
}

/// The DNA of a gene's exons, and where the introns were: for each, the place in the shortened DNA and its size.
pub(super) struct Cdna {
    pub dna: Vec<u8>,
    pub introns: Vec<(i64, i64)>,
}

/// `uppercase_cdna`: take the introns out (in the order of their last places, each taken from the DNA as it is by then).
pub(super) fn uppercase_cdna(dna: &[u8]) -> Cdna {
    let mut list = stretches(dna);
    list.sort_by_key(|&(_, last)| last);
    // `lsort -unique` keeps the last of the stretches that end in the same place
    let mut kept: Vec<(i64, i64)> = Vec::new();
    for pair in list {
        match kept.last_mut() {
            Some(last) if last.1 == pair.1 => *last = pair,
            _ => kept.push(pair),
        }
    }
    let mut text = chars(&String::from_utf8_lossy(dna));
    let mut introns = Vec::new();
    let mut taken = 0i64;
    for (first, last) in kept {
        text = tcl::replace(&text, first - taken, last - taken, &[]);
        introns.push((first - taken, last - first + 1));
        taken += last - first + 1;
    }
    let dna = text.iter().map(|&c| c as u8).collect();
    Cdna { dna, introns }
}

/// The line of amino acids under a gene: its exons translated, with room left where its introns are.
pub(super) fn genes_line(text: &[u8], code: &[u8]) -> Vec<char> {
    let cdna = uppercase_cdna(text);
    let mut line = chars(&translate(&cdna.dna, code, false, true, 2));
    let mut introns = cdna.introns;
    // `lsort -decreasing -integer -index 0` is stable
    introns.sort_by_key(|&(at, _)| std::cmp::Reverse(at));
    for (at, size) in introns {
        let mut next = range(&line, 0, at - 1);
        next.extend(spaces(size));
        next.extend(range(&line, at, line.len() as i64));
        line = next;
    }
    line
}
