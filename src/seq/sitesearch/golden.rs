//! Reads the reference files made by running the original's Silent Sites and Find Diagnostic Site code, unchanged, on
//! random sequences (the harness is described in MECHANICS.md): one block per case, from `CASE kind n` to `END`.

use crate::seq::enzymes::{default_enzymes, find_sites};

pub struct Case {
    pub circular: bool,
    pub seq: Vec<u8>,
    pub selection: Option<(usize, usize)>,
    /// The enzymes chosen in the Enzyme Selector; none at all is empty.
    pub chosen: Vec<String>,
    pub mode: String,
    pub mismatches: usize,
    /// The 1-based starts of the sites the original found, for each enzyme that has any.
    pub sites: Vec<(String, Vec<usize>)>,
    /// What the original wrote into its result window (the date as `DATE`, the title as `TITLE`).
    pub text: String,
    /// The lines of the case that this reader has no field for, as (word, the rest of the line): the choices of List Enzymes.
    pub other: Vec<(String, String)>,
}

impl Case {
    fn new() -> Case {
        Case {
            circular: false,
            seq: Vec::new(),
            selection: None,
            chosen: Vec::new(),
            mode: String::new(),
            mismatches: 0,
            sites: Vec::new(),
            text: String::new(),
            other: Vec::new(),
        }
    }

    /// The stretch searched, `start..end`.
    pub fn region(&self) -> (usize, usize) {
        self.selection.unwrap_or((0, self.seq.len()))
    }

    /// What a line with this word at its start says after it.
    pub fn say(&self, word: &str) -> &str {
        let found = self.other.iter().find(|(w, _)| w == word);
        found.map_or("", |(_, rest)| rest.as_str())
    }
}

fn numbers(words: &[&str]) -> Vec<usize> {
    words
        .iter()
        .map(|w| w.parse().expect("a number in the reference file"))
        .collect()
}

fn read_line(case: &mut Case, line: &str) {
    let words: Vec<&str> = line.split(' ').collect();
    let rest = &words[1..];
    match words[0] {
        "circular" => case.circular = rest[0] == "1",
        "seq" => case.seq = rest[0].as_bytes().to_vec(),
        "sel" if rest[0] != "-" => {
            let n = numbers(rest);
            case.selection = Some((n[0], n[1]));
        }
        "enz" if rest[0] != "-" => case.chosen = rest.iter().map(|s| s.to_string()).collect(),
        "mode" => case.mode = rest[0].to_string(),
        "mismatches" => case.mismatches = numbers(rest)[0],
        "sites" => case.sites.push((rest[0].to_string(), numbers(&rest[1..]))),
        other => case.other.push((other.to_string(), rest.join(" "))),
    }
}

pub fn parse(file: &str) -> Vec<Case> {
    let mut cases = Vec::new();
    let mut open: Option<(Case, Vec<&str>, bool)> = None;
    for line in file.lines() {
        match (&mut open, line) {
            (_, l) if l.starts_with("CASE ") => open = Some((Case::new(), Vec::new(), false)),
            (Some((case, lines, _)), "END") => {
                case.text = lines.join("\n");
                cases.extend(open.take().map(|(c, _, _)| c));
            }
            (Some((_, _, reading)), "out") => *reading = true,
            (Some((_, lines, true)), l) => lines.push(&l[1..]),
            (Some((case, _, false)), l) => read_line(case, l),
            (None, _) => {}
        }
    }
    cases
}

/// The sites of every default enzyme as the build finds them, against the lists the original's `findenzymes` made.
pub fn check_sites(cases: &[Case]) {
    let all = default_enzymes();
    for (n, case) in cases.iter().enumerate() {
        for e in &all {
            let mut mine: Vec<usize> = find_sites(e, &case.seq, case.circular)
                .iter()
                .map(|s| s.start + 1)
                .collect();
            mine.sort_unstable();
            let mut theirs = case
                .sites
                .iter()
                .find(|(name, _)| *name == e.name)
                .map(|(_, at)| at.clone())
                .unwrap_or_default();
            theirs.sort_unstable();
            assert_eq!(mine, theirs, "case {n}: sites of {}", e.name);
        }
    }
}
