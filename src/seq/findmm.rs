//! The original's Find (`do_find_mm`, AppMain.tcl line 13981), the search it uses for every nucleotide pattern:
//! IUPAC letters, several patterns at once, the reverse complement, mismatches and "minimum 3' bases", and
//! exactly where Find Next, Find Prev and Highlight All look, quirks included. The original's own code was
//! run on 300 random cases to make `findmm/golden.txt`, which the tests here compare with.

use super::dna::revcom;

/// The letters a Find pattern keeps (`regsub -all -nocase {[^abcdghkmnrstuvwy\*\|\;\,]}`): the base and
/// ambiguity letters (not E F I J L O P Q X Z), `*`, and the separators `|` `;` `,` between patterns.
const PATTERN_LETTERS: &str = "abcdghkmnrstuvwy*|;,";

/// What the Find dialog holds. `Default` is the original's settings (`info(find_*)`, AppMain.tcl L66966-66973
/// and L66920): not case sensitive, the reverse complement too, no mismatches, all bases of the pattern.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Query {
    pub pattern: String,
    pub case_sensitive: bool,
    /// "also find rev-com of string"
    pub rev_com: bool,
    /// "as literal": ambiguity letters in the pattern stand only for themselves
    pub literal: bool,
    /// "Allow mismatches" (the dialog offers 0 to 25)
    pub mismatches: usize,
    /// "Minimum 3' bases": only this many bases from the end of each pattern are used; `None` is "All"
    pub min_3_bases: Option<usize>,
    /// The original always wraps (`info(find_wrap)` is 1 and its checkbox is commented out)
    pub wrap: bool,
}

impl Default for Query {
    fn default() -> Self {
        Query {
            pattern: String::new(),
            case_sensitive: false,
            rev_com: true,
            literal: false,
            mismatches: 0,
            min_3_bases: None,
            wrap: true,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Action {
    Next,
    Prev,
    Highlight,
}

/// One highlighted stretch (`foundf` forward, `foundr` reverse): bases `start` up to `end`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Tag {
    pub reverse: bool,
    pub start: usize,
    pub end: usize,
}

#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Outcome {
    /// Select bases `start` up to `end` (which may run past the end of the sequence).
    Selected {
        start: usize,
        end: usize,
        message: String,
    },
    NotFound {
        message: String,
    },
    Highlighted {
        tags: Vec<Tag>,
        message: String,
    },
}

/// The bases a pattern letter stands for (`string map` in `do_find_mm`); other letters stand for themselves.
fn expand(c: u8) -> Vec<u8> {
    let set: &[u8] = match c {
        b'B' => b"BYSKCGT",
        b'D' => b"DRWKAGT",
        b'H' => b"HMWYACT",
        b'K' => b"KGT",
        b'M' => b"MAC",
        b'N' => b"NBDHVKMRYWSACGT",
        b'R' => b"RAG",
        b'S' => b"SGC",
        b'V' => b"VMSRACG",
        b'W' => b"WAT",
        b'Y' => b"YCT",
        b'b' => b"byskcgt",
        b'd' => b"drwkagt",
        b'h' => b"hmwyact",
        b'k' => b"kgt",
        b'm' => b"mac",
        b'n' => b"nbdhvkmrywsacgt",
        b'r' => b"rag",
        b's' => b"sgc",
        b'v' => b"vmsracg",
        b'w' => b"wat",
        b'y' => b"yct",
        other => return vec![other],
    };
    set.to_vec()
}

/// Every pattern to look for, each as the bases allowed at each place, and the length of the longest.
/// With the reverse complement each pattern is followed by its reverse complement.
fn prepare(q: &Query) -> (Vec<Vec<Vec<u8>>>, usize) {
    let kept: String = q
        .pattern
        .chars()
        .filter(|c| c.is_ascii() && PATTERN_LETTERS.contains(c.to_ascii_lowercase()))
        .collect();
    let pattern = if q.case_sensitive {
        kept
    } else {
        kept.to_ascii_uppercase()
    };
    // Tcl's `split` of an empty string gives no pieces at all
    let pieces: Vec<&str> = if pattern.is_empty() {
        Vec::new()
    } else {
        pattern.split(['|', ',', ';']).collect()
    };
    let (mut all, mut longest) = (Vec::new(), 0);
    for piece in pieces {
        let bytes = piece.as_bytes();
        let used = match q.min_3_bases {
            Some(n) => &bytes[bytes.len().saturating_sub(n)..],
            None => bytes,
        };
        all.push(used.to_vec());
        if q.rev_com {
            all.push(revcom(used));
        }
        longest = longest.max(used.len());
    }
    let allowed = all
        .into_iter()
        .map(|p| {
            p.into_iter()
                .map(|c| if q.literal { vec![c] } else { expand(c) })
                .collect()
        })
        .collect();
    (allowed, longest)
}

/// Tcl's `string range`: first to last, both included, held inside the text; nothing when `last < first`.
fn range(text: &[u8], first: i64, last: i64) -> Vec<u8> {
    let first = first.max(0);
    let last = last.min(text.len() as i64 - 1);
    if first > last {
        Vec::new()
    } else {
        text[first as usize..=last as usize].to_vec()
    }
}

/// Whether `pattern` matches the bases at the start of `window` with at most `allowed` mismatches. A base
/// that is missing (the window is short) counts as a mismatch.
fn matches_here(pattern: &[Vec<u8>], window: &[u8], allowed: usize) -> bool {
    let mut mismatches = 0;
    for (k, set) in pattern.iter().enumerate() {
        if !window.get(k).is_some_and(|c| set.contains(c)) {
            mismatches += 1;
            if mismatches > allowed {
                return false;
            }
        }
    }
    true
}

/// The first of the patterns that matches at `at`, by its place in the list.
fn first_match(patterns: &[Vec<Vec<u8>>], list: &[u8], at: usize, allowed: usize) -> Option<usize> {
    patterns.iter().position(|p| {
        let end = (at + p.len()).min(list.len());
        matches_here(p, &list[at.min(list.len())..end], allowed)
    })
}

fn direction(q: &Query, pattern_index: usize) -> &'static str {
    if q.rev_com && pattern_index % 2 == 1 {
        "reverse"
    } else {
        "forward"
    }
}

/// Run one of Find Next, Find Prev and Highlight All on `text`. The search starts from the selection's first
/// base when something is selected, otherwise from the caret.
pub fn find(
    text: &[u8],
    circular: bool,
    q: &Query,
    action: Action,
    caret: usize,
    selection_start: Option<usize>,
) -> Outcome {
    let (patterns, longest) = prepare(q);
    let upper_text = if q.case_sensitive {
        text.to_vec()
    } else {
        text.to_ascii_uppercase()
    };
    let text = upper_text.as_slice();
    let n = text.len() as i64;
    let none = || Outcome::NotFound {
        message: "Found no match.".into(),
    };
    if n == 0 {
        return none();
    }
    let longest = longest as i64;
    let allowed = q.mismatches;
    match action {
        Action::Next => {
            let ins = selection_start.map_or(caret as i64, |s| s as i64 + 1);
            let mut list = range(text, ins, n - 1);
            if circular {
                list.extend(range(text, 0, ins - 1));
                let twice = [text, &range(text, 0, longest)[..]].concat();
                list.extend(range(&twice, ins, ins + longest - 1));
            } else if q.wrap {
                list.push(b'\t');
                list.extend(range(text, 0, ins - 1));
                list.extend(range(text, ins, ins + longest - 1));
            }
            let found = (0..=(list.len() as i64 - longest))
                .find_map(|i| first_match(&patterns, &list, i as usize, allowed).map(|p| (i, p)));
            let Some((mut at, p)) = found else {
                return none();
            };
            if !circular && q.wrap && at + ins > n {
                at -= 1;
            }
            let start = (at + ins).rem_euclid(n);
            Outcome::Selected {
                start: start as usize,
                end: start as usize + patterns[p].len(),
                message: format!("Found {} match at {}.", direction(q, p), start + 1),
            }
        }
        Action::Prev => {
            let ins = selection_start.map_or(caret as i64, |s| s as i64);
            let list = if circular {
                let thrice = [text, text, &range(text, 0, longest)[..]].concat();
                range(&thrice, ins, n + ins + longest - 1)
            } else if q.wrap {
                let mut both = range(text, ins, n - 1);
                both.push(b'\t');
                both.extend(range(text, 0, longest + ins - 2));
                both
            } else {
                range(text, 0, longest + ins - 1)
            };
            let found = (0..=(list.len() as i64 - 1 - longest))
                .rev()
                .find_map(|i| first_match(&patterns, &list, i as usize, allowed).map(|p| (i, p)));
            let Some((at, p)) = found else {
                return none();
            };
            let len = patterns[p].len();
            let (start, shown) = if circular {
                let start = (at + ins).rem_euclid(n);
                (start, start + 1)
            } else if q.wrap {
                let gap = list.iter().position(|&c| c == b'\t').unwrap_or(0) as i64;
                if at > gap {
                    (at - gap - 1, at - gap - 1)
                } else {
                    (at + ins, at + ins + 1)
                }
            } else {
                (at, at + 1)
            };
            Outcome::Selected {
                start: start as usize,
                end: start as usize + len,
                message: format!("Found {} match at {shown}.", direction(q, p)),
            }
        }
        Action::Highlight => {
            let mut list = text.to_vec();
            if circular {
                list.extend(range(text, 0, longest - 1));
            }
            let (mut tags, mut count) = (Vec::new(), 0usize);
            for at in 0..list.len() {
                for (p, st) in patterns.iter().enumerate() {
                    let end = (at + st.len()).min(list.len());
                    if matches_here(st, &list[at..end], allowed) {
                        let reverse = q.rev_com && p % 2 == 1;
                        let to = (at + st.len()).min(n as usize);
                        if at < to {
                            tags.push(Tag {
                                reverse,
                                start: at,
                                end: to,
                            });
                        }
                        let wrapped = (at + st.len()) % n as usize;
                        if at + st.len() > n as usize && wrapped > 0 {
                            tags.push(Tag {
                                reverse,
                                start: 0,
                                end: wrapped,
                            });
                        }
                        count += 1;
                    }
                }
            }
            if tags.is_empty() {
                none()
            } else {
                Outcome::Highlighted {
                    tags,
                    message: format!("Found {count} matches."),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests;
