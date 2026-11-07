//! The species-tree job: gene trees go in, the species tree and what its branches show come out (`seq/astral`). Like the other
//! kinds it runs the same function on the workstation and on this computer, and travels as a container of sections: the
//! settings as `key=value` lines, the trees as Newick, the branches as tab-separated lines.

use super::tasks::{bit, flag, lines, number, MAX_TREES};
use super::wire::{meta, Container};
use super::{clean_line, one_line};
use crate::seq::astral::{self, BranchOut, Options};
use std::collections::HashMap;
use web_time::Instant;

/// The most taxa a result may name, and a species tree to score may have: far more than the search could ever finish with.
const MAX_TAXA: usize = 100_000;
/// The most leaves a request may hold in all its gene trees (about the species times the genes): the trees are read into a
/// node for every few characters, so this keeps what one request can ask of the workstation to a few gigabytes.
const MAX_LEAVES: usize = 20_000_000;

/// What to run: the gene trees and the settings of the search (see `astral::Options`).
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesJob {
    pub genes: String,
    /// A binary species tree to score and annotate; empty to search.
    pub score_tree: String,
    pub exact: bool,
    pub extra: u8,
    pub polish: bool,
    pub lambda: f64,
    /// The taxon to hang the tree from; empty for the first.
    pub root: String,
    pub min_leaves: usize,
}

impl SpeciesJob {
    /// A search with the defaults of ASTRAL.
    pub fn search(genes: &str) -> SpeciesJob {
        let o = Options::default();
        SpeciesJob {
            genes: genes.to_string(),
            score_tree: String::new(),
            exact: o.exact,
            extra: o.extra,
            polish: o.polish,
            lambda: o.lambda,
            root: String::new(),
            min_leaves: o.min_leaves,
        }
    }

    pub fn options(&self) -> Options {
        Options {
            exact: self.exact,
            extra: self.extra,
            polish: self.polish,
            lambda: self.lambda,
            root: (!self.root.is_empty()).then(|| self.root.clone()),
            min_leaves: self.min_leaves,
            ..Options::default()
        }
    }

    pub fn encode(&self) -> Container {
        Container::new()
            .with(
                "meta",
                lines(&[
                    ("exact", bit(self.exact)),
                    ("extra", self.extra.to_string()),
                    ("polish", bit(self.polish)),
                    ("lambda", self.lambda.to_string()),
                    ("root", clean_line(&self.root)),
                    ("min_leaves", self.min_leaves.to_string()),
                ]),
            )
            .with("genes", self.genes.clone())
            .with("species", self.score_tree.clone())
    }

    /// A request comes from the network: the gene trees are counted by their semicolons and their leaves by their commas (no
    /// tree is read before the counts are within limits), and every setting is checked.
    pub fn decode(c: &Container, m: &HashMap<&str, &str>) -> Result<SpeciesJob, String> {
        let genes = c.require("genes")?;
        let trees = genes.matches(';').count();
        if trees > MAX_TREES {
            return Err(format!(
                "A job may hold at most {MAX_TREES} gene trees; this one has {trees}."
            ));
        }
        let leaves = genes.matches([',', ';']).count();
        if leaves > MAX_LEAVES {
            return Err(format!(
                "A job may hold at most {MAX_LEAVES} leaves in all its gene trees; this one has about {leaves}."
            ));
        }
        let score_tree = c.text("species").unwrap_or_default();
        if score_tree.matches(',').count() > MAX_TAXA {
            return Err(format!(
                "The tree to score may have at most {MAX_TAXA} species."
            ));
        }
        let (extra, lambda): (u8, f64) = (number(m, "extra")?, number(m, "lambda")?);
        if extra > 1 || !(lambda >= 0.0 && lambda.is_finite()) {
            return Err("The message's search settings are out of range.".into());
        }
        Ok(SpeciesJob {
            genes: genes.to_string(),
            score_tree: score_tree.to_string(),
            exact: flag(m, "exact"),
            extra,
            polish: flag(m, "polish"),
            lambda,
            root: m.get("root").copied().unwrap_or_default().to_string(),
            min_leaves: number(m, "min_leaves")?,
        })
    }
}

/// What a run found.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesResult {
    pub newick: String,
    pub score: u64,
    pub normalized: f64,
    pub taxa: Vec<String>,
    pub genes: usize,
    pub skipped: usize,
    /// How many groups the search chose from; none when a tree was only scored.
    pub clusters: Option<usize>,
    pub branches: Vec<BranchOut>,
    pub seconds: f64,
}

/// Run the job here.
pub fn run(job: &SpeciesJob) -> Result<SpeciesResult, String> {
    let started = Instant::now();
    let o = job.options();
    let out = if job.score_tree.trim().is_empty() {
        astral::infer(&job.genes, &o)?
    } else {
        astral::score_tree(&job.genes, &job.score_tree, &o)?
    };
    Ok(SpeciesResult {
        newick: out.newick,
        score: out.score,
        normalized: out.normalized,
        taxa: out.taxa,
        genes: out.genes,
        skipped: out.skipped,
        clusters: out.clusters,
        branches: out.branches,
        seconds: started.elapsed().as_secs_f64(),
    })
}

fn numbers(values: &[f64]) -> String {
    values
        .iter()
        .map(f64::to_string)
        .collect::<Vec<_>>()
        .join("\t")
}

impl SpeciesResult {
    pub fn encode(&self) -> Container {
        let index: HashMap<&str, usize> = self
            .taxa
            .iter()
            .enumerate()
            .map(|(i, t)| (t.as_str(), i))
            .collect();
        let rows: String = self
            .branches
            .iter()
            .map(|b| {
                let side: Vec<String> = b
                    .taxa
                    .iter()
                    .filter_map(|t| index.get(t.as_str()).map(usize::to_string))
                    .collect();
                let values = [b.f, b.q, b.pp].concat();
                let all = [values, vec![b.en, b.length]].concat();
                format!("{}\t{}\t{}\n", side.join(","), b.quartets, numbers(&all))
            })
            .collect();
        let taxa: String = self
            .taxa
            .iter()
            .map(|t| format!("{}\n", one_line(t)))
            .collect();
        Container::new()
            .with(
                "meta",
                lines(&[
                    ("score", self.score.to_string()),
                    ("normalized", self.normalized.to_string()),
                    ("genes", self.genes.to_string()),
                    ("skipped", self.skipped.to_string()),
                    (
                        "clusters",
                        self.clusters.map(|c| c.to_string()).unwrap_or_default(),
                    ),
                    ("seconds", self.seconds.to_string()),
                ]),
            )
            .with("newick", self.newick.clone())
            .with("taxa", taxa)
            .with("branches", rows)
    }

    pub fn decode(c: &Container) -> Result<SpeciesResult, String> {
        let m = meta(c.require("meta")?);
        let taxa: Vec<String> = c.require("taxa")?.lines().map(str::to_string).collect();
        if taxa.len() > MAX_TAXA {
            return Err("The result names too many taxa.".into());
        }
        let branches = c
            .require("branches")?
            .lines()
            .map(|row| branch(row, &taxa))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(SpeciesResult {
            newick: c.require("newick")?.to_string(),
            score: number(&m, "score")?,
            normalized: number(&m, "normalized")?,
            genes: number(&m, "genes")?,
            skipped: number(&m, "skipped")?,
            clusters: m.get("clusters").and_then(|v| v.parse().ok()),
            seconds: number(&m, "seconds")?,
            taxa,
            branches,
        })
    }
}

/// One line of the branch table: the taxa on one side, the number of quartets, and eleven numbers.
fn branch(row: &str, taxa: &[String]) -> Result<BranchOut, String> {
    let bad = || "A branch in the result is damaged.".to_string();
    let fields: Vec<&str> = row.split('\t').collect();
    if fields.len() != 13 {
        return Err(bad());
    }
    let side = fields[0]
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<usize>().ok().and_then(|i| taxa.get(i)).cloned())
        .collect::<Option<Vec<String>>>()
        .ok_or_else(bad)?;
    let v: Vec<f64> = fields[2..]
        .iter()
        .map(|x| x.parse().map_err(|_| bad()))
        .collect::<Result<_, _>>()?;
    Ok(BranchOut {
        taxa: side,
        quartets: fields[1].parse().map_err(|_| bad())?,
        f: [v[0], v[1], v[2]],
        q: [v[3], v[4], v[5]],
        pp: [v[6], v[7], v[8]],
        en: v[9],
        length: v[10],
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const GENES: &str =
        "((A,B),(C,(D,E)));\n((A,B),(C,(D,E)));\n((A,C),(B,(D,E)));\n(((A,B),C),(D,E));\n";

    #[test]
    fn the_job_goes_through_bytes_and_comes_back_as_it_was() {
        let mut job = SpeciesJob::search(GENES);
        job.exact = true;
        job.lambda = 0.25;
        job.root = "C".into();
        job.min_leaves = 4;
        let bytes = job.encode().encode();
        let c = Container::decode(&bytes).expect("a container");
        let m = meta(c.require("meta").expect("meta"));
        assert_eq!(SpeciesJob::decode(&c, &m).expect("a job"), job);
    }

    #[test]
    fn a_job_with_settings_out_of_range_or_too_many_trees_is_refused() {
        let decode = |meta_text: &str, genes: &str| {
            let c = Container::new()
                .with("meta", meta_text.to_string())
                .with("genes", genes.to_string());
            SpeciesJob::decode(&c, &meta(c.require("meta").expect("meta")))
        };
        let good = "extra=1\nlambda=0.5\nmin_leaves=0\nexact=0\npolish=1\n";
        assert!(decode(good, GENES).is_ok());
        assert!(decode("extra=2\nlambda=0.5\nmin_leaves=0\n", GENES).is_err());
        assert!(decode("extra=1\nlambda=-1\nmin_leaves=0\n", GENES).is_err());
        assert!(decode("extra=1\nlambda=NaN\nmin_leaves=0\n", GENES).is_err());
        assert!(
            decode("extra=1\nlambda=0.5\n", GENES).is_err(),
            "min_leaves is missing"
        );
        assert!(decode(good, &";".repeat(MAX_TREES + 1)).is_err());
        // a few trees can hold a great many leaves: they are counted by their commas before any tree is read
        assert!(decode(good, &",".repeat(MAX_LEAVES)).is_ok());
        let why = decode(good, &",".repeat(MAX_LEAVES + 1)).expect_err("too many leaves");
        assert!(why.contains("leaves"), "{why}");
        // and so is the tree to score
        let c = Container::new()
            .with("meta", good.to_string())
            .with("genes", GENES.to_string())
            .with("species", ",".repeat(MAX_TAXA + 1));
        let m = meta(c.require("meta").expect("meta"));
        assert!(SpeciesJob::decode(&c, &m).is_err());
    }

    #[test]
    fn a_result_goes_through_bytes_and_the_run_gives_what_the_engine_gives() {
        let found = run(&SpeciesJob::search(GENES)).expect("a species tree");
        assert_eq!(found.taxa, ["A", "B", "C", "D", "E"]);
        assert_eq!(found.branches.len(), 2);
        assert!(found.score > 0 && found.clusters.is_some());
        let c = Container::decode(&found.encode().encode()).expect("a container");
        let back = SpeciesResult::decode(&c).expect("a result");
        assert_eq!(back, found);
        // scoring a given tree has no search to count
        let mut job = SpeciesJob::search(GENES);
        job.score_tree = "((A,B),C,(D,E));".into();
        let scored = run(&job).expect("a scored tree");
        assert!(scored.clusters.is_none());
        assert_eq!(scored.score, found.score);
        assert!(run(&SpeciesJob::search("(A,B,C);")).is_err());
    }

    /// The container of `found` with other rows in its branch table (`with` adds a section, it does not replace one).
    fn with_rows(found: &SpeciesResult, rows: &str) -> Container {
        let good = found.encode();
        let section = |name| good.text(name).expect("a section").to_string();
        Container::new()
            .with("meta", section("meta"))
            .with("newick", section("newick"))
            .with("taxa", section("taxa"))
            .with("branches", rows)
    }

    #[test]
    fn a_damaged_result_is_refused() {
        let found = run(&SpeciesJob::search(GENES)).expect("a species tree");
        let rows = found.encode().text("branches").expect("rows").to_string();
        let ok = |rows: &str| SpeciesResult::decode(&with_rows(&found, rows));
        assert_eq!(ok(&rows).expect("the genuine rows"), found);
        // a row cut short, a taxon number past the end, a word where a number should be, a taxon named by a word
        let row = "0\t1\t1\t0\t0\t1\t0\t0\t1\t0\t0\t2\t0.5\n";
        assert!(ok(row).is_ok(), "the model row is whole");
        assert!(ok("9\t1\t0\t0\t0\n").is_err());
        assert!(
            ok(&row.replacen('0', "99", 1)).is_err(),
            "a taxon number past the end"
        );
        assert!(ok(&row.replace("0.5", "slow")).is_err());
        assert!(ok(&row.replacen('0', "first", 1)).is_err());
        // a section missing, and a table of more taxa than a result may name
        let meta_only = Container::new().with(
            "meta",
            found.encode().text("meta").expect("meta").to_string(),
        );
        assert!(SpeciesResult::decode(&meta_only).is_err());
        let many = "x\n".repeat(MAX_TAXA + 1);
        let c = with_rows(&found, &rows);
        let c = Container::new()
            .with("meta", c.text("meta").expect("meta").to_string())
            .with("newick", c.text("newick").expect("newick").to_string())
            .with("taxa", many)
            .with("branches", "");
        assert!(SpeciesResult::decode(&c).is_err());
    }
}
