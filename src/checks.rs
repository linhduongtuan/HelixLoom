//! How each analysis was checked (`docs/ROADMAP.md` 2.4): what it was compared with, how close it came, what was not checked, and
//! the tests that make the comparison on every change. The page of the Help menu and `helixloom checks` show this table; a test
//! (`checks/tests.rs`) fails when the table names a test that is not there or does not run, so the page cannot claim a check that
//! nobody makes.

mod table;

pub use table::CHECKS;

/// What an analysis was compared with, strongest first.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Standard {
    Reference,
    Independent,
    Known,
    Untested,
}

impl Standard {
    pub const ALL: [Standard; 4] = [
        Standard::Reference,
        Standard::Independent,
        Standard::Known,
        Standard::Untested,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Standard::Reference => "reference program",
            Standard::Independent => "second implementation",
            Standard::Known => "known answers",
            Standard::Untested => "not checked",
        }
    }

    /// What the grade means, in a sentence.
    pub fn meaning(self) -> &'static str {
        match self {
            Standard::Reference => {
                "Compared with the program it follows, or with a reference program, on the same inputs."
            }
            Standard::Independent => {
                "Compared with a separate, plain implementation of the same method, or with its definition worked out the slow way: this shows the arithmetic is right, not that the method suits your data."
            }
            Standard::Known => {
                "Tried on cases whose answer is known: worked by hand, a closed form, or something planted and found again."
            }
            Standard::Untested => {
                "Not compared with anything here: either nobody has yet, or the answer is another program's own."
            }
        }
    }
}

/// How one analysis was checked.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Check {
    /// The heading it is listed under, and its name as the menus and pages say it.
    pub group: &'static str,
    pub analysis: &'static str,
    pub standard: Standard,
    /// What it was compared with, how close it came, and what was not checked.
    pub against: &'static str,
    pub result: &'static str,
    pub limits: &'static str,
    /// The tests that make the comparison: the file, from the top of the repository, and the name of the test.
    pub tests: &'static [(&'static str, &'static str)],
}

impl Check {
    /// Whether `word` (in small letters) is in any of its texts.
    fn holds(&self, word: &str) -> bool {
        let has = |text: &str| text.to_lowercase().contains(word);
        has(self.analysis)
            || has(self.group)
            || has(self.standard.label())
            || has(self.against)
            || has(self.result)
            || has(self.limits)
            || self.tests.iter().any(|(file, test)| has(file) || has(test))
    }
}

/// The checks that hold every word of `query`, whatever the case, in the order of the table; all of them for no word.
pub fn matching(query: &str) -> Vec<&'static Check> {
    let words: Vec<String> = query.split_whitespace().map(str::to_lowercase).collect();
    CHECKS
        .iter()
        .filter(|c| words.iter().all(|w| c.holds(w)))
        .collect()
}

/// The check of the method that made a result, from the record of the result: its kind and its settings, as the pages and
/// the windows write them (each has a `described`, and names its analysis where it starts a run). None for a result that the
/// table does not know, and for a run of the Pipelines page's fixed route, which is several methods with a check each.
pub fn of_result(analysis: &str, parameters: &[(String, String)]) -> Option<&'static Check> {
    let setting = |key: &str| {
        parameters
            .iter()
            .find(|(k, _)| k == key)
            .map_or("", |(_, v)| v.as_str())
    };
    let name = match analysis {
        "Alignment" if setting("engine") == "MAFFT" => {
            "MAFFT, IQ-TREE and the other programs HelixLoom runs"
        }
        "Alignment" if setting("engine") == "built in" => "Alignment (built in)",
        "Alignment" if setting("engine") == "progressive (ClustalW-style)" => {
            "Alignment (progressive, after ClustalW)"
        }
        "Tree estimation" => match setting("method") {
            "Neighbour-joining" | "UPGMA" => "Distance trees (neighbour-joining, UPGMA)",
            "Maximum likelihood" => "Maximum-likelihood trees",
            "Bayesian" => "Bayesian trees",
            _ => return None,
        },
        "Species tree" if setting("gene families") == "yes" => {
            "Gene families (duplications and losses)"
        }
        "Species tree" if !matches!(setting("weights"), "" | "none") => {
            "Weighted species trees (by support and branch length)"
        }
        "Species tree" => "Species tree from gene trees",
        "Motif search" => "Motif Finder (scan, discovery, comparison)",
        "Protein domain search" => "Protein domains (HMMER, Pfam)",
        "Structure prediction" => "Structure prediction (a folding program)",
        "PhyloCNN estimates" => "PhyloCNN (epidemic parameters from a tree)",
        "Pipeline file" => "Pipeline files (run, cache, export)",
        _ => return None,
    };
    CHECKS.iter().find(|c| c.analysis == name)
}

/// `checks` under their headings, the headings in the order they first come.
pub fn grouped<'a>(checks: &[&'a Check]) -> Vec<(&'static str, Vec<&'a Check>)> {
    let mut groups: Vec<(&'static str, Vec<&'a Check>)> = Vec::new();
    for &check in checks {
        match groups.iter_mut().find(|(group, _)| *group == check.group) {
            Some((_, list)) => list.push(check),
            None => groups.push((check.group, vec![check])),
        }
    }
    groups
}

/// How many of `checks` have each grade.
pub fn tally(checks: &[&Check]) -> [(Standard, usize); 4] {
    Standard::ALL.map(|s| (s, checks.iter().filter(|c| c.standard == s).count()))
}

/// `checks` as text, for the clipboard and the command line.
pub fn to_text(checks: &[&Check]) -> String {
    let counts: Vec<String> = tally(checks)
        .iter()
        .map(|(s, n)| format!("{n} {}", s.label()))
        .collect();
    let mut text = format!(
        "How each analysis was checked: {} {} ({})\n",
        checks.len(),
        if checks.len() == 1 {
            "analysis"
        } else {
            "analyses"
        },
        counts.join(", ")
    );
    for (group, list) in grouped(checks) {
        text.push_str(&format!("\n{group}\n"));
        for c in list {
            let tests: Vec<String> = c
                .tests
                .iter()
                .map(|(file, test)| format!("{test} ({file})"))
                .collect();
            text.push_str(&format!(
                "\n{} [{}]\n  compared with: {}\n  result: {}\n  not checked: {}\n  tests: {}\n",
                c.analysis,
                c.standard.label(),
                c.against,
                c.result,
                c.limits,
                if tests.is_empty() {
                    "none".to_string()
                } else {
                    tests.join("; ")
                }
            ));
        }
    }
    text
}

#[cfg(test)]
mod tests;
