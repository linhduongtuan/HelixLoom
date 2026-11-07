//! The engine against ASTRAL 5.7.8, run as a black box on ten generated data sets (`golden/*.txt`): the gene trees, the species
//! trees ASTRAL was asked to score (`-q`, with the full annotation `-t 2`), what it printed, and the best scores its
//! searches found (default, `-p 0`, and `-x` where the data are small). The data are simulated under the multispecies
//! coalescent, with polytomies and missing taxa; `golden/README` says how they were made.

use super::{infer, score_tree, Options};
use std::collections::{BTreeSet, HashMap};

const FILES: [(&str, &str); 10] = [
    ("a1", include_str!("golden/a1.txt")),
    ("a2", include_str!("golden/a2.txt")),
    ("a3", include_str!("golden/a3.txt")),
    ("a4", include_str!("golden/a4.txt")),
    ("a5", include_str!("golden/a5.txt")),
    ("a6", include_str!("golden/a6.txt")),
    ("b1", include_str!("golden/b1.txt")),
    ("b2", include_str!("golden/b2.txt")),
    ("b3", include_str!("golden/b3.txt")),
    ("b4", include_str!("golden/b4.txt")),
];

#[derive(Default)]
struct Case {
    genes: String,
    /// label, the species tree, what ASTRAL printed for it, its score and normalized score
    scored: Vec<(String, String, String, u64, f64)>,
    /// kind (default, p0, x), score, normalized score
    best: HashMap<String, (u64, f64)>,
}

fn parse(text: &str) -> Case {
    let mut case = Case::default();
    let (mut section, mut label) = ("", String::new());
    let (mut species, mut annotated) = (String::new(), String::new());
    for line in text.lines() {
        let words: Vec<&str> = line.split(' ').collect();
        match words[0] {
            "GENES" => section = "genes",
            "SPECIES" => {
                section = "species";
                label = words[1].to_string();
            }
            "ASTRAL" => section = "astral",
            "TREE" => section = "tree",
            "SCORE" => case.scored.push((
                label.clone(),
                species.clone(),
                annotated.clone(),
                words[2].parse().expect("a score"),
                words[3].parse().expect("a number"),
            )),
            "BEST" => {
                let score = words[2].parse().expect("a score");
                case.best.insert(
                    words[1].to_string(),
                    (score, words[3].parse().expect("a number")),
                );
            }
            _ => match section {
                "genes" => {
                    case.genes.push_str(line);
                    case.genes.push('\n');
                }
                "species" => species = line.to_string(),
                "astral" => annotated = line.to_string(),
                _ => {}
            },
        }
    }
    case
}

/// One annotated inner branch of ASTRAL's output: the taxa below it, its numbers, and its length.
struct Printed {
    clade: BTreeSet<String>,
    numbers: HashMap<String, f64>,
    length: Option<f64>,
}

struct Reader<'a> {
    text: &'a [u8],
    at: usize,
    found: Vec<Printed>,
}

impl Reader<'_> {
    fn node(&mut self) -> BTreeSet<String> {
        if self.text[self.at] == b'(' {
            self.at += 1;
            let mut clade = BTreeSet::new();
            loop {
                clade.extend(self.node());
                if self.text[self.at] == b',' {
                    self.at += 1;
                } else {
                    break;
                }
            }
            self.at += 1; // the closing bracket
            let numbers = self.annotation();
            let length = self.length();
            if !numbers.is_empty() {
                self.found.push(Printed {
                    clade: clade.clone(),
                    numbers,
                    length,
                });
            }
            clade
        } else {
            let start = self.at;
            while !b",():;".contains(&self.text[self.at]) {
                self.at += 1;
            }
            let name = String::from_utf8_lossy(&self.text[start..self.at]).into_owned();
            self.length();
            BTreeSet::from([name])
        }
    }

    /// The quoted `[q1=..;f1=..]` after a bracket, as numbers by name.
    fn annotation(&mut self) -> HashMap<String, f64> {
        let mut numbers = HashMap::new();
        if self.text[self.at] == b'\'' {
            let end = self.text[self.at + 1..]
                .iter()
                .position(|&c| c == b'\'')
                .expect("a closing quote");
            let inside = String::from_utf8_lossy(&self.text[self.at + 2..self.at + 1 + end - 1])
                .into_owned();
            for part in inside.split(';') {
                if let Some((k, v)) = part.split_once('=') {
                    numbers.insert(k.to_string(), v.parse().expect("a number"));
                }
            }
            self.at += end + 2;
        }
        numbers
    }

    fn length(&mut self) -> Option<f64> {
        if self.text.get(self.at) != Some(&b':') {
            return None;
        }
        self.at += 1;
        let start = self.at;
        while !b",();".contains(&self.text[self.at]) {
            self.at += 1;
        }
        String::from_utf8_lossy(&self.text[start..self.at])
            .parse()
            .ok()
    }
}

fn printed(tree: &str) -> Vec<Printed> {
    let mut reader = Reader {
        text: tree.as_bytes(),
        at: 0,
        found: Vec::new(),
    };
    reader.node();
    reader.found
}

/// The side of a split without the first taxon, so the same split reads the same however a tree is hung.
fn canonical(clade: &BTreeSet<String>, taxa: &[String]) -> BTreeSet<String> {
    if clade.contains(&taxa[0]) {
        taxa.iter()
            .filter(|t| !clade.contains(*t))
            .cloned()
            .collect()
    } else {
        clade.clone()
    }
}

fn close(a: f64, b: f64, relative: f64) -> bool {
    (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-300) || (a - b).abs() < 1e-12
}

fn sorted_pair(a: f64, b: f64) -> [f64; 2] {
    if a <= b {
        [a, b]
    } else {
        [b, a]
    }
}

#[test]
fn scores_and_annotations_equal_astrals_on_every_species_tree_it_scored() {
    let mut compared = 0;
    for (name, text) in FILES {
        let case = parse(text);
        for (label, species, ast, score, normalized) in &case.scored {
            let o = Options {
                annotate: 2,
                ..Options::default()
            };
            let mine = score_tree(&case.genes, species, &o).expect("scored");
            assert_eq!(mine.score, *score, "{name} {label}: quartet score");
            assert!(
                close(mine.normalized, *normalized, 1e-12),
                "{name} {label}: normalized {} against {normalized}",
                mine.normalized
            );
            // a rooted tree has its root branch annotated from both sides: the same split twice
            let theirs = printed(ast);
            let distinct: BTreeSet<BTreeSet<String>> = theirs
                .iter()
                .map(|p| canonical(&p.clade, &mine.taxa))
                .collect();
            assert_eq!(
                distinct.len(),
                mine.taxa.len() - 3,
                "{name} {label}: inner branches"
            );
            let mut seen = BTreeSet::new();
            let by_split: HashMap<BTreeSet<String>, &super::BranchOut> = mine
                .branches
                .iter()
                .map(|b| (canonical(&b.taxa.iter().cloned().collect(), &mine.taxa), b))
                .collect();
            for p in &theirs {
                if !seen.insert(canonical(&p.clade, &mine.taxa)) {
                    continue;
                }
                let b = by_split
                    .get(&canonical(&p.clade, &mine.taxa))
                    .unwrap_or_else(|| panic!("{name} {label}: no branch for {:?}", p.clade));
                let get = |k: &str| p.numbers[k];
                let where_ = format!("{name} {label} {:?}", p.clade);
                assert!(
                    close(b.f[0], get("f1"), 1e-9),
                    "{where_}: f1 {} against {}",
                    b.f[0],
                    get("f1")
                );
                assert_eq!(
                    sorted_pair(b.f[1], b.f[2]).map(|x| (x * 1e6).round()),
                    sorted_pair(get("f2"), get("f3")).map(|x| (x * 1e6).round()),
                    "{where_}: f2 f3"
                );
                assert!(close(b.en, get("EN"), 1e-9), "{where_}: EN");
                assert_eq!(b.quartets as f64, get("QC"), "{where_}: QC");
                assert!(
                    close(b.pp[0], get("pp1"), 1e-8),
                    "{where_}: pp1 {} against {}",
                    b.pp[0],
                    get("pp1")
                );
                assert!(close(b.q[0], get("q1"), 1e-9), "{where_}: q1");
                let length = p.length.expect("a length");
                assert!(
                    close(b.length, length, 1e-9),
                    "{where_}: length {} against {length}",
                    b.length
                );
                compared += 1;
            }
        }
    }
    assert!(compared > 100, "{compared} branches compared");
}

#[test]
fn the_exact_search_finds_the_optimum_that_astrals_exact_search_found() {
    let mut checked = 0;
    for (name, text) in FILES {
        let case = parse(text);
        let Some(&(score, normalized)) = case.best.get("x") else {
            continue;
        };
        let o = Options {
            exact: true,
            polish: false,
            ..Options::default()
        };
        let mine = infer(&case.genes, &o).expect("a tree");
        assert_eq!(mine.score, score, "{name}: the best score over every tree");
        assert!(
            close(mine.normalized, normalized, 1e-12),
            "{name}: normalized"
        );
        checked += 1;
    }
    assert_eq!(
        checked, 6,
        "the data sets small enough for the exact search"
    );
}

#[test]
fn the_search_reaches_astrals_scores_with_and_without_the_extra_clusters() {
    for (name, text) in FILES {
        let case = parse(text);
        let run = |o: Options| infer(&case.genes, &o).expect("a tree");
        let default = run(Options::default());
        let (best, normalized) = case.best["default"];
        assert!(
            default.score >= best,
            "{name}: {} against ASTRAL's {best}",
            default.score
        );
        if default.score == best {
            assert!(
                close(default.normalized, normalized, 1e-12),
                "{name}: normalized"
            );
        }
        // with only the gene trees' own clusters, as with ASTRAL's -p 0
        let plain = run(Options {
            extra: 0,
            polish: false,
            ..Options::default()
        });
        assert!(plain.score >= case.best["p0"].0, "{name}: no extras");
        assert!(plain.score <= default.score);
    }
}
