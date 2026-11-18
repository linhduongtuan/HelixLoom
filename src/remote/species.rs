//! The species-tree job: gene trees go in, the species tree and what its branches show come out (`seq/astral`). Like the other
//! kinds it runs the same function on the workstation and on this computer, and travels as a container of sections: the
//! settings as `key=value` lines, the trees as Newick, the branches as tab-separated lines.

use super::tasks::{bit, flag, lines, number, MAX_TREES};
use super::wire::{meta, Container};
use super::{clean_line, one_line, MAX_TEXT};
use crate::seq::astral::bootstrap::{DEFAULT_SEED, MAX_REPLICATES};
use crate::seq::astral::polytomy::MIN_GENES;
use crate::seq::astral::{self, BootInfo, BranchOut, Options, Polytomy, Resample, Weights};
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
    /// How the quartets of the gene trees are weighted. A server from before the weights ignores it, which is why the result
    /// says how it was weighted.
    pub weights: Weights,
    /// The gene trees are gene family trees (ASTRAL-Pro); a server from before them ignores this too, and says in the result.
    pub families: bool,
    /// Bootstrap replicates (none when 0), the seed that draws them, and what they resample; a server from before the bootstrap
    /// ignores them too, and says in the result.
    pub replicates: usize,
    pub seed: u64,
    pub resample: Resample,
    /// For resampling the sites: the bootstrap trees of each gene tree of the input, one text for each, in its order.
    pub site_trees: Vec<String>,
    /// Test every branch for being a polytomy; a server from before the test ignores this, and says in the result.
    pub polytomy: bool,
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
            weights: o.weights,
            families: o.families,
            replicates: o.bootstrap.replicates,
            seed: o.bootstrap.seed,
            resample: o.bootstrap.resample,
            site_trees: Vec::new(),
            polytomy: o.polytomy,
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
            bootstrap: astral::Boot {
                replicates: self.replicates,
                seed: self.seed,
                resample: self.resample,
            },
            site_trees: self.site_trees.clone(),
            families: self.families,
            weights: self.weights,
            polytomy: self.polytomy,
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
                    ("weights", self.weights.name().to_string()),
                    ("families", bit(self.families)),
                    ("replicates", self.replicates.to_string()),
                    ("seed", self.seed.to_string()),
                    ("resample", self.resample.name().to_string()),
                    ("polytomy", bit(self.polytomy)),
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
            .with(
                "site_sizes",
                self.site_trees
                    .iter()
                    .map(|t| format!("{}\n", t.len()))
                    .collect::<String>(),
            )
            .with("site_trees", self.site_trees.concat())
    }

    /// What this job asks of the server beyond the search itself, to check that the answer did it.
    pub fn asked(&self) -> Asked {
        Asked {
            weights: self.weights,
            families: self.families,
            bootstrap: (self.replicates > 0).then_some((self.resample, self.replicates)),
            polytomy: self.polytomy,
        }
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
        let site_trees = site_trees(c)?;
        let leaves = leaves
            + site_trees
                .iter()
                .map(|t| t.matches([',', ';']).count())
                .sum::<usize>();
        if leaves > MAX_LEAVES {
            return Err(format!(
                "A job may hold at most {MAX_LEAVES} leaves in all its trees, the bootstrap trees too; this one has about {leaves}."
            ));
        }
        let score_tree = c.text("species").unwrap_or_default();
        if score_tree.matches(',').count() > MAX_TAXA {
            return Err(format!(
                "The tree to score may have at most {MAX_TAXA} species."
            ));
        }
        // a message that says nothing of the weights (from a program from before them) asks for none
        let weights = match m.get("weights") {
            None => Weights::None,
            Some(w) => Weights::parse(w)
                .ok_or("The message's weights are not none, support, length or hybrid.")?,
        };
        let (extra, lambda): (u8, f64) = (number(m, "extra")?, number(m, "lambda")?);
        if extra > 1 || !(lambda >= 0.0 && lambda.is_finite()) {
            return Err("The message's search settings are out of range.".into());
        }
        // a message that says nothing of the bootstrap (from a program from before it) asks for none
        let replicates: usize = match m.get("replicates") {
            None => 0,
            Some(_) => number(m, "replicates")?,
        };
        if replicates > MAX_REPLICATES {
            return Err(format!(
                "A job may ask for at most {MAX_REPLICATES} bootstrap replicates; this one asks for {replicates}."
            ));
        }
        let seed: u64 = match m.get("seed") {
            None => DEFAULT_SEED,
            Some(_) => number(m, "seed")?,
        };
        let resample = match m.get("resample") {
            None => Resample::Genes,
            Some(r) => {
                Resample::parse(r).ok_or("The message's resampling is not genes, sites or both.")?
            }
        };
        Ok(SpeciesJob {
            genes: genes.to_string(),
            score_tree: score_tree.to_string(),
            weights,
            families: flag(m, "families"),
            replicates,
            seed,
            resample,
            site_trees,
            polytomy: flag(m, "polytomy"),
            exact: flag(m, "exact"),
            extra,
            polish: flag(m, "polish"),
            lambda,
            root: m.get("root").copied().unwrap_or_default().to_string(),
            min_leaves: number(m, "min_leaves")?,
        })
    }
}

/// The most bytes of bootstrap trees a request may carry.
const MAX_SITE_BYTES: usize = 256 << 20;

/// The bootstrap trees of a request: the texts of the genes, cut out by their sizes.
fn site_trees(c: &Container) -> Result<Vec<String>, String> {
    let (sizes, body) = (
        c.text("site_sizes").unwrap_or_default(),
        c.text("site_trees").unwrap_or_default(),
    );
    let bad = |why: &str| format!("The bootstrap trees in the message do not add up: {why}.");
    if body.len() > MAX_SITE_BYTES {
        return Err(format!(
            "A job may carry at most {} MB of bootstrap trees.",
            MAX_SITE_BYTES >> 20
        ));
    }
    let sizes: Vec<&str> = sizes.lines().collect();
    if sizes.len() > MAX_TREES {
        return Err(format!(
            "A job may hold the bootstrap trees of at most {MAX_TREES} genes; this one has {}.",
            sizes.len()
        ));
    }
    let mut start = 0usize;
    let mut texts = Vec::with_capacity(sizes.len());
    for size in sizes {
        let size: usize = size.parse().map_err(|_| bad("a size is not a number"))?;
        let end = start
            .checked_add(size)
            .ok_or_else(|| bad("a size is too large"))?;
        let text = body
            .get(start..end)
            .ok_or_else(|| bad("a size runs past the end or cuts a character"))?;
        texts.push(text.to_string());
        start = end;
    }
    if start != body.len() {
        return Err(bad("bytes are left over"));
    }
    Ok(texts)
}

/// What a job asked of the server beyond the search itself.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Asked {
    pub weights: Weights,
    pub families: bool,
    pub bootstrap: Option<(Resample, usize)>,
    pub polytomy: bool,
}

/// The result of a job that asked for these things, or why it cannot be used: a server from before the weights, the gene
/// families, the bootstrap and the polytomy test ignores what it does not know and says nothing of it, so its answer was not
/// weighted, not read as gene families, not bootstrapped and not tested, whatever was asked, and must not be shown as if it
/// were.
pub fn as_asked(asked: Asked, result: SpeciesResult) -> Result<SpeciesResult, String> {
    let Asked {
        weights,
        families,
        bootstrap,
        polytomy,
    } = asked;
    if result.polytomy != polytomy {
        return Err("The workstation did not test the branches for polytomies as asked: it is probably running an older version of the program.".into());
    }
    if result.bootstrap.map(|b| (b.resample, b.asked)) != bootstrap {
        return Err("The workstation did not do the bootstrap as asked: it is probably running an older version of the program.".into());
    }
    if result.weights != weights {
        return Err(format!(
            "The workstation answered without the weights that were asked for (\"{}\", not \"{}\"): it is probably running an older version of the program.",
            result.weights.name(),
            weights.name()
        ));
    }
    if result.families != families {
        return Err(format!(
            "The workstation did not read the trees as {} as asked: it is probably running an older version of the program.",
            if families { "gene families" } else { "single-copy trees" }
        ));
    }
    Ok(result)
}

/// What a run found.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesResult {
    pub newick: String,
    /// The quartets of the tree that the gene trees share, each counted one.
    pub score: u64,
    /// How the quartets were weighted, and the weights of those that the tree shares added up (when they were).
    pub weights: Weights,
    pub weighted_score: Option<f64>,
    /// The trees were read as gene families (ASTRAL-Pro): the score counts classes of quartets and no branch has a support.
    pub families: bool,
    /// What the bootstrap did, when one was asked for.
    pub bootstrap: Option<BootInfo>,
    /// The branches were tested for being polytomies (those with enough gene trees have the result in their row).
    pub polytomy: bool,
    /// The score (the weighted one in a weighted run) as a share of the quartets the gene trees resolve.
    pub normalized: f64,
    pub taxa: Vec<String>,
    pub genes: usize,
    pub skipped: usize,
    /// How many groups the search chose from; none when a tree was only scored.
    pub clusters: Option<usize>,
    /// The inner branches with what the gene trees say of them. A table too large to send (a tree as deep as it is wide
    /// lists every species below every branch) comes with its first rows only; `branch_count` says how many there are.
    pub branches: Vec<BranchOut>,
    pub branch_count: usize,
    pub seconds: f64,
}

/// What a branch table may take, which a reply is not trusted to keep to: the entries (a species named on one side of a
/// branch) and the bytes of names they copy. A tree that is a ladder of n species has n²/2 entries.
const MAX_ENTRIES: usize = 1_000_000;
const MAX_COPIED_BYTES: usize = 64 << 20;
/// The longest name a reply may give a species: `MAX_TEXT` characters of up to four bytes.
const MAX_NAME_BYTES: usize = 4 * MAX_TEXT;

/// How much of the table's allowance is left.
struct Budget {
    entries: usize,
    bytes: usize,
}

impl Budget {
    fn new() -> Budget {
        Budget {
            entries: MAX_ENTRIES,
            bytes: MAX_COPIED_BYTES,
        }
    }

    /// Take what a row naming species of these name lengths costs; false, and nothing taken, when it is more than is left.
    fn spend(&mut self, name_lengths: &[usize]) -> bool {
        let bytes: usize = name_lengths.iter().sum();
        if name_lengths.len() > self.entries || bytes > self.bytes {
            return false;
        }
        self.entries -= name_lengths.len();
        self.bytes -= bytes;
        true
    }
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
        weights: out.weights,
        weighted_score: out.weighted_score,
        families: out.families,
        bootstrap: out.bootstrap,
        polytomy: out.polytomy,
        normalized: out.normalized,
        taxa: out.taxa,
        genes: out.genes,
        skipped: out.skipped,
        clusters: out.clusters,
        branch_count: out.branches.len(),
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
    /// The branch table as rows of tab-separated fields, cut after the rows that fit the allowance of a table (`names` are the
    /// species names as the other end will read them, which is what the allowance is counted in); `branch_count` still says
    /// how many branches there are.
    fn table(&self, names: &[String]) -> String {
        let index: HashMap<&str, usize> = self
            .taxa
            .iter()
            .enumerate()
            .map(|(i, t)| (t.as_str(), i))
            .collect();
        let mut allowance = Budget::new();
        self.branches
            .iter()
            .map_while(|b| {
                let side: Vec<usize> = b
                    .taxa
                    .iter()
                    .filter_map(|t| index.get(t.as_str()).copied())
                    .collect();
                let lengths: Vec<usize> = side.iter().map(|&i| names[i].len()).collect();
                allowance.spend(&lengths).then(|| {
                    let values = [b.f, b.q, b.pp].concat();
                    let all = [values, vec![b.en, b.length]].concat();
                    let ids: Vec<String> = side.iter().map(usize::to_string).collect();
                    let boot = b.bootstrap.map(|x| x.to_string()).unwrap_or_default();
                    let test = b
                        .polytomy
                        .map(|t| format!("{},{},{}", t.genes, t.chi2, t.p))
                        .unwrap_or_default();
                    format!(
                        "{}\t{}\t{}\t{boot}\t{test}\n",
                        ids.join(","),
                        b.quartets,
                        numbers(&all)
                    )
                })
            })
            .collect()
    }

    pub fn encode(&self) -> Container {
        let names: Vec<String> = self.taxa.iter().map(|t| one_line(t)).collect();
        let taxa: String = names.iter().map(|t| format!("{t}\n")).collect();
        Container::new()
            .with(
                "meta",
                lines(&[
                    ("score", self.score.to_string()),
                    ("weights", self.weights.name().to_string()),
                    ("families", bit(self.families)),
                    ("polytomy", bit(self.polytomy)),
                    (
                        "weighted_score",
                        self.weighted_score
                            .map(|w| w.to_string())
                            .unwrap_or_default(),
                    ),
                    ("normalized", self.normalized.to_string()),
                    (
                        "boot_resample",
                        self.bootstrap
                            .map(|b| b.resample.name().to_string())
                            .unwrap_or_default(),
                    ),
                    (
                        "boot_asked",
                        self.bootstrap
                            .map(|b| b.asked.to_string())
                            .unwrap_or_default(),
                    ),
                    (
                        "boot_done",
                        self.bootstrap
                            .map(|b| b.done.to_string())
                            .unwrap_or_default(),
                    ),
                    (
                        "boot_failed",
                        self.bootstrap
                            .map(|b| b.failed.to_string())
                            .unwrap_or_default(),
                    ),
                    ("genes", self.genes.to_string()),
                    ("skipped", self.skipped.to_string()),
                    (
                        "clusters",
                        self.clusters.map(|c| c.to_string()).unwrap_or_default(),
                    ),
                    ("branch_count", self.branch_count.to_string()),
                    ("seconds", self.seconds.to_string()),
                ]),
            )
            .with("newick", self.newick.clone())
            .with("taxa", taxa)
            .with("branches", self.table(&names))
    }

    /// A result comes from a server that is not trusted to keep to what it was asked: the names, the rows of the branch table
    /// and what they copy are all bounded before anything is made of them.
    pub fn decode(c: &Container) -> Result<SpeciesResult, String> {
        let m = meta(c.require("meta")?);
        let taxa: Vec<String> = c.require("taxa")?.lines().map(str::to_string).collect();
        if taxa.len() > MAX_TAXA || taxa.iter().any(|t| t.len() > MAX_NAME_BYTES) {
            return Err("The result names too many taxa, or one with too long a name.".into());
        }
        let rows: Vec<&str> = c.require("branches")?.lines().collect();
        // an unrooted binary tree of n species has n - 3 inner branches
        if rows.len() > taxa.len() {
            return Err("The result has more branches than species.".into());
        }
        let mut allowance = Budget::new();
        let branches = rows
            .iter()
            .map(|row| branch(row, &taxa, &mut allowance))
            .collect::<Result<Vec<_>, _>>()?;
        let branch_count: usize = number(&m, "branch_count")?;
        if branch_count < branches.len() || branch_count > taxa.len() {
            return Err("The result counts its branches wrongly.".into());
        }
        // a result that says nothing of the test (from a program from before it) did none, and a test cannot have counted
        // more gene trees than the result has
        let genes: usize = number(&m, "genes")?;
        let polytomy = flag(&m, "polytomy");
        if branches
            .iter()
            .any(|b| b.polytomy.is_some_and(|t| !polytomy || t.genes > genes))
        {
            return Err(
                "The result's polytomy tests do not agree with what it says it did.".into(),
            );
        }
        // a result that says nothing of the weights (from a program from before them) was not weighted
        let weights = match m.get("weights") {
            None => Weights::None,
            Some(w) => Weights::parse(w)
                .ok_or("The result's weights are not none, support, length or hybrid.")?,
        };
        let weighted_score = match (weights, m.get("weighted_score").filter(|v| !v.is_empty())) {
            (Weights::None, None) => None,
            (w, Some(v)) if w != Weights::None => Some(
                v.parse::<f64>()
                    .ok()
                    .filter(|x| x.is_finite() && *x >= 0.0)
                    .ok_or("The result's weighted score is not a number.")?,
            ),
            _ => return Err("The result's weights and its weighted score do not agree.".into()),
        };
        // a result that says nothing of a bootstrap (from a program from before it) had none
        let bootstrap = match m.get("boot_asked").filter(|v| !v.is_empty()) {
            None => None,
            Some(_) => {
                let (asked, done, failed): (usize, usize, usize) = (
                    number(&m, "boot_asked")?,
                    number(&m, "boot_done")?,
                    number(&m, "boot_failed")?,
                );
                let resample = Resample::parse(m.get("boot_resample").copied().unwrap_or_default())
                    .ok_or("The result's resampling is not genes, sites or both.")?;
                if asked == 0
                    || asked > MAX_REPLICATES
                    || done == 0
                    || done.checked_add(failed) != Some(asked)
                {
                    return Err("The result counts its bootstrap replicates wrongly.".into());
                }
                Some(BootInfo {
                    resample,
                    asked,
                    done,
                    failed,
                })
            }
        };
        Ok(SpeciesResult {
            newick: c.require("newick")?.to_string(),
            score: number(&m, "score")?,
            weights,
            weighted_score,
            families: flag(&m, "families"),
            bootstrap,
            polytomy,
            normalized: number(&m, "normalized")?,
            genes,
            skipped: number(&m, "skipped")?,
            clusters: m.get("clusters").and_then(|v| v.parse().ok()),
            seconds: number(&m, "seconds")?,
            taxa,
            branches,
            branch_count,
        })
    }
}

/// One line of the branch table: the taxa on one side (numbers, each once and in order), the number of quartets, and eleven
/// numbers; then the bootstrap support and the test for a polytomy (`genes,chi2,p`) when the branch has them.
fn branch(row: &str, taxa: &[String], allowance: &mut Budget) -> Result<BranchOut, String> {
    let bad = || "A branch in the result is damaged.".to_string();
    let fields: Vec<&str> = row.split('\t').collect();
    // thirteen fields, a fourteenth that holds the bootstrap support (from a program that has it) and a fifteenth with the test
    if !(13..=15).contains(&fields.len()) {
        return Err(bad());
    }
    let ids: Vec<usize> = fields[0]
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse::<usize>().map_err(|_| bad()))
        .collect::<Result<_, _>>()?;
    if ids.iter().any(|&i| i >= taxa.len()) || ids.windows(2).any(|w| w[0] >= w[1]) {
        return Err(bad());
    }
    let lengths: Vec<usize> = ids.iter().map(|&i| taxa[i].len()).collect();
    if !allowance.spend(&lengths) {
        return Err("The result's table of branches is too large to take.".into());
    }
    let side: Vec<String> = ids.iter().map(|&i| taxa[i].clone()).collect();
    let v: Vec<f64> = fields[2..13]
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
        bootstrap: match fields.get(13).filter(|x| !x.is_empty()) {
            None => None,
            Some(x) => Some(
                x.parse::<f64>()
                    .ok()
                    .filter(|s| (0.0..=100.0).contains(s))
                    .ok_or_else(bad)?,
            ),
        },
        polytomy: match fields.get(14).filter(|x| !x.is_empty()) {
            None => None,
            Some(x) => Some(polytomy_test(x).ok_or_else(bad)?),
        },
    })
}

/// The test for a polytomy in a row: `genes,chi2,p`, from at least `MIN_GENES` gene trees, a statistic that is not negative
/// and a p-value that is a probability.
fn polytomy_test(text: &str) -> Option<Polytomy> {
    let mut part = text.split(',');
    let (genes, chi2, p) = (part.next()?, part.next()?, part.next()?);
    if part.next().is_some() {
        return None;
    }
    let (genes, chi2, p): (usize, f64, f64) =
        (genes.parse().ok()?, chi2.parse().ok()?, p.parse().ok()?);
    (genes >= MIN_GENES && chi2.is_finite() && chi2 >= 0.0 && (0.0..=1.0).contains(&p))
        .then_some(Polytomy { genes, chi2, p })
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

    /// The same sections with the lines of `meta` that `edit` changes.
    fn meta_edited(c: &Container, edit: impl Fn(&str) -> String) -> Container {
        let meta = edit(c.text("meta").expect("meta"));
        ["genes", "species", "newick", "taxa", "branches"]
            .into_iter()
            .filter_map(|name| c.bytes(name).map(|b| (name, b.to_vec())))
            .fold(Container::new().with("meta", meta), |c, (n, b)| {
                c.with(n, b)
            })
    }

    fn without_line(text: &str, key: &str) -> String {
        text.lines()
            .filter(|l| !l.starts_with(&format!("{key}=")))
            .map(|l| format!("{l}\n"))
            .collect()
    }

    #[test]
    fn the_weights_of_a_job_travel_and_a_message_that_says_nothing_of_them_asks_for_none() {
        for weights in Weights::ALL {
            let mut job = SpeciesJob::search(GENES);
            job.weights = weights;
            let c = Container::decode(&job.encode().encode()).expect("a container");
            let m = meta(c.require("meta").expect("meta"));
            assert_eq!(SpeciesJob::decode(&c, &m).expect("a job").weights, weights);
        }
        let mut job = SpeciesJob::search(GENES);
        job.weights = Weights::Hybrid;
        // a program from before the weights sends no such line: the job is unweighted
        let old = meta_edited(&job.encode(), |t| without_line(t, "weights"));
        let m = meta(old.require("meta").expect("meta"));
        assert_eq!(
            SpeciesJob::decode(&old, &m).expect("a job").weights,
            Weights::None
        );
        let odd = meta_edited(&job.encode(), |t| {
            t.replace("weights=hybrid", "weights=heavy")
        });
        let m = meta(odd.require("meta").expect("meta"));
        assert!(SpeciesJob::decode(&odd, &m)
            .unwrap_err()
            .contains("not none, support"));
    }

    #[test]
    fn a_weighted_result_goes_through_bytes_and_its_weights_and_score_must_agree() {
        let weighted = {
            let mut job =
                SpeciesJob::search("(A:1,B:1,(C:1,D:1)0.9:1);\n(A:1,B:1,(C:1,D:1)0.8:1);");
            job.weights = Weights::Support;
            run(&job).expect("a result")
        };
        assert_eq!(
            (weighted.weights, weighted.weighted_score.is_some()),
            (Weights::Support, true)
        );
        let c = weighted.encode();
        assert_eq!(SpeciesResult::decode(&c).expect("a result"), weighted);
        // an unweighted result has no weighted score
        let plain = run(&SpeciesJob::search(GENES)).expect("a result");
        assert_eq!(
            SpeciesResult::decode(&plain.encode()).expect("a result"),
            plain
        );
        assert_eq!((plain.weights, plain.weighted_score), (Weights::None, None));
        // a result from before the weights says nothing of them: it was not weighted
        let old = meta_edited(&plain.encode(), |t| {
            without_line(&without_line(t, "weights"), "weighted_score")
        });
        assert_eq!(
            SpeciesResult::decode(&old).expect("a result").weights,
            Weights::None
        );
        // the two must agree, and the score must be a number
        let say = |weights: &str, score: &str| {
            let c = meta_edited(&weighted.encode(), |t| {
                without_line(&without_line(t, "weights"), "weighted_score")
                    + &format!("weights={weights}\nweighted_score={score}\n")
            });
            SpeciesResult::decode(&c)
        };
        assert!(say("support", "1.5").is_ok());
        for (weights, score) in [
            ("support", ""),
            ("none", "1.5"),
            ("support", "NaN"),
            ("support", "-1"),
            ("support", "inf"),
            ("support", "x"),
            ("heavy", "1"),
        ] {
            assert!(say(weights, score).is_err(), "{weights} {score}");
        }
    }

    const FAMILIES: &str =
        "(((A,B),(C,(D,E))),((A,B),(C,(D,E))));\n(((A,B),(C,(D,E))),((A,B),(D,E)));\n";

    #[test]
    fn the_families_flag_travels_and_a_result_says_whether_it_was_read_as_families() {
        let mut job = SpeciesJob::search(FAMILIES);
        assert!(!job.families);
        job.families = true;
        let c = Container::decode(&job.encode().encode()).expect("a container");
        let m = meta(c.require("meta").expect("meta"));
        assert!(SpeciesJob::decode(&c, &m).expect("a job").families);
        // a program from before them sends no such line: single-copy trees
        let old = meta_edited(&job.encode(), |t| without_line(t, "families"));
        let m = meta(old.require("meta").expect("meta"));
        assert!(!SpeciesJob::decode(&old, &m).expect("a job").families);
        // the result says it, and goes through bytes
        let found = run(&job).expect("a result");
        assert!(found.families && found.branches.is_empty());
        assert_eq!(
            SpeciesResult::decode(&found.encode()).expect("a result"),
            found
        );
        let plain = run(&SpeciesJob::search(GENES)).expect("a result");
        assert!(!plain.families);
        // a result is taken only as it was asked for
        let asked = |weights, families| Asked {
            weights,
            families,
            bootstrap: None,
            polytomy: false,
        };
        assert!(as_asked(asked(Weights::None, true), found.clone()).is_ok());
        let why = as_asked(asked(Weights::None, true), plain).unwrap_err();
        assert!(
            why.contains("gene families") && why.contains("older version"),
            "{why}"
        );
        let why = as_asked(asked(Weights::None, false), found).unwrap_err();
        assert!(why.contains("single-copy"), "{why}");
    }

    #[test]
    fn the_bootstrap_of_a_job_travels_with_the_trees_of_the_sites_and_old_messages_ask_for_none() {
        let mut job = SpeciesJob::search(GENES);
        job.replicates = 12;
        job.seed = 5;
        job.resample = Resample::Both;
        job.site_trees = vec![
            "((A,B),(C,D));\n".into(),
            "((A,C),(B,D));\n((é,B),(C,D));\n".into(),
            String::new(),
        ];
        let c = Container::decode(&job.encode().encode()).expect("a container");
        let m = meta(c.require("meta").expect("meta"));
        assert_eq!(SpeciesJob::decode(&c, &m).expect("a job"), job);
        // a program from before the bootstrap sends none of it: none is asked for
        let old = meta_edited(&SpeciesJob::search(GENES).encode(), |t| {
            without_line(
                &without_line(&without_line(t, "replicates"), "seed"),
                "resample",
            )
        });
        let m = meta(old.require("meta").expect("meta"));
        let asked = SpeciesJob::decode(&old, &m).expect("a job");
        assert_eq!(
            (asked.replicates, asked.seed, asked.resample),
            (0, 692, Resample::Genes)
        );
        assert!(asked.site_trees.is_empty());
        // what is out of range or does not add up is refused
        let with = |edit: &dyn Fn(&str) -> String| {
            let c = meta_edited(&job.encode(), edit);
            SpeciesJob::decode(&c, &meta(c.require("meta").expect("meta")))
        };
        assert!(with(&|t| t.replace("replicates=12", "replicates=1001"))
            .unwrap_err()
            .contains("at most 1000"));
        assert!(with(&|t| t.replace("resample=both", "resample=some"))
            .unwrap_err()
            .contains("not genes, sites or both"));
        assert!(with(&|t| t.replace("seed=5", "seed=x")).is_err());
        let sized = |sizes: &str| {
            let c = Container::new()
                .with("meta", job.encode().text("meta").expect("meta").to_string())
                .with("genes", GENES)
                .with("site_sizes", sizes.to_string())
                .with("site_trees", "abcdef");
            SpeciesJob::decode(&c, &meta(c.require("meta").expect("meta")))
        };
        assert!(sized("3\n3\n").is_ok());
        for bad in ["3\n2\n", "3\n9\n", "x\n", "3\n3\n1\n"] {
            assert!(sized(bad).unwrap_err().contains("do not add up"), "{bad}");
        }
    }

    #[test]
    fn a_result_with_a_bootstrap_goes_through_bytes_and_one_that_counts_wrongly_is_refused() {
        let mut job = SpeciesJob::search("((A,B),(C,D));\n".repeat(6).as_str());
        job.replicates = 8;
        let found = run(&job).expect("a result");
        let info = found.bootstrap.expect("what the bootstrap did");
        assert_eq!((info.asked, info.done, info.failed), (8, 8, 0));
        assert_eq!(found.branches[0].bootstrap, Some(100.0));
        assert_eq!(
            SpeciesResult::decode(&found.encode()).expect("a result"),
            found
        );
        let say = |done: &str, failed: &str, resample: &str| {
            let c = meta_edited(&found.encode(), |t| {
                ["boot_asked", "boot_done", "boot_failed", "boot_resample"]
                    .iter()
                    .fold(t.to_string(), |t, k| without_line(&t, k))
                    + &format!("boot_asked=8\nboot_done={done}\nboot_failed={failed}\nboot_resample={resample}\n")
            });
            SpeciesResult::decode(&c)
        };
        assert!(say("6", "2", "genes").is_ok());
        for (done, failed, resample) in [
            ("5", "2", "genes"),
            ("0", "8", "genes"),
            ("9", "0", "genes"),
            ("8", "0", "most"),
        ] {
            assert!(
                say(done, failed, resample).is_err(),
                "{done} {failed} {resample}"
            );
        }
        // a branch row from a program from before it has thirteen fields, and a support must be a share
        let row = |extra: &str| {
            let t = found.table(&found.taxa.iter().map(|t| one_line(t)).collect::<Vec<_>>());
            let line = t.lines().next().expect("a row");
            let thirteen = line.split('\t').take(13).collect::<Vec<_>>().join("\t");
            branch(
                &format!("{thirteen}{extra}"),
                &found.taxa,
                &mut Budget::new(),
            )
        };
        assert_eq!(row("").expect("thirteen fields").bootstrap, None);
        assert_eq!(row("\t87.5").expect("fourteen").bootstrap, Some(87.5));
        assert!(row("\t101").is_err() && row("\t-1").is_err() && row("\tx").is_err());
    }

    /// Ten gene trees of one topology, ten of another and ten of the third, repeated so that every branch has enough to be tested.
    const MANY: &str = "((A,B),(C,(D,E)));\n((A,C),(B,(D,E)));\n(((A,B),C),(D,E));\n";

    fn many() -> String {
        MANY.repeat(10)
    }

    #[test]
    fn the_polytomy_test_of_a_job_travels_and_an_old_message_asks_for_none() {
        let mut job = SpeciesJob::search(GENES);
        assert!(!job.polytomy && !job.options().polytomy);
        job.polytomy = true;
        assert!(job.options().polytomy && job.asked().polytomy);
        let c = Container::decode(&job.encode().encode()).expect("a container");
        let m = meta(c.require("meta").expect("meta"));
        assert_eq!(SpeciesJob::decode(&c, &m).expect("a job"), job);
        // a message from a program from before the test says nothing of it
        let old = meta_edited(&job.encode(), |t| without_line(t, "polytomy"));
        let m = meta(old.require("meta").expect("meta"));
        assert!(!SpeciesJob::decode(&old, &m).expect("a job").polytomy);
    }

    #[test]
    fn a_result_with_the_test_goes_through_bytes_and_one_that_disagrees_with_itself_is_refused() {
        let mut job = SpeciesJob::search(&many());
        job.polytomy = true;
        let found = run(&job).expect("a result");
        assert!(found.polytomy && !found.branches.is_empty());
        assert!(
            found.branches.iter().all(|b| b.polytomy.is_some()),
            "thirty gene trees are enough for every branch"
        );
        let back = SpeciesResult::decode(&found.encode()).expect("a result that decodes");
        assert!(back.polytomy);
        assert_eq!(back.branches, found.branches);
        // the job that did not ask for the test gets rows with nothing in the last field
        let plain = run(&SpeciesJob::search(&many())).expect("a result");
        assert!(!plain.polytomy && plain.branches.iter().all(|b| b.polytomy.is_none()));
        let plain_back = SpeciesResult::decode(&plain.encode()).expect("a result");
        assert!(!plain_back.polytomy);
        // a test in a row of a result that says it did none, or from more gene trees than the result has, cannot be right
        let said = |edit: &dyn Fn(&str) -> String| {
            let bytes = found.encode();
            let c = meta_edited(&bytes, |t| edit(t));
            let c = c.clone().with(
                "branches",
                bytes.text("branches").expect("rows").to_string(),
            );
            SpeciesResult::decode(&c)
        };
        assert!(said(&|t| t.to_string()).is_ok());
        assert!(said(&|t| t.replace("polytomy=1", "polytomy=0")).is_err());
        assert!(said(&|t| t.replace(&format!("genes={}", found.genes), "genes=5")).is_err());
    }

    #[test]
    fn a_row_of_the_table_takes_a_polytomy_test_only_if_it_is_one() {
        let found = run(&SpeciesJob::search(GENES)).expect("a result");
        let row = |last: &str| {
            let t = found.table(&found.taxa.iter().map(|t| one_line(t)).collect::<Vec<_>>());
            let line = t.lines().next().expect("a row");
            let thirteen = line.split('\t').take(13).collect::<Vec<_>>().join("\t");
            branch(
                &format!("{thirteen}\t\t{last}"),
                &found.taxa,
                &mut Budget::new(),
            )
        };
        let ok = row("30,15,0.0005530843701478336").expect("a test");
        let t = ok.polytomy.expect("the test");
        assert_eq!((t.genes, t.chi2, t.p), (30, 15.0, 0.0005530843701478336));
        assert_eq!(row("").expect("none").polytomy, None);
        for bad in [
            "9,15,0.0005", // fewer gene trees than the test is made for
            "30,-1,0.5",   // a statistic is not negative
            "30,15,1.5",   // a p-value is a probability
            "30,15,-0.1",
            "30,NaN,0.5",
            "30,inf,0.5",
            "30,15",       // a field short
            "30,15,0.5,1", // a field over
            "x,15,0.5",
            "30.5,15,0.5",
        ] {
            assert!(row(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn an_answer_that_did_not_test_the_branches_as_asked_is_refused() {
        let plain = run(&SpeciesJob::search(&many())).expect("a result");
        let mut job = SpeciesJob::search(&many());
        job.polytomy = true;
        let why = as_asked(job.asked(), plain.clone()).unwrap_err();
        assert!(
            why.contains("polytomies") && why.contains("older version"),
            "{why}"
        );
        let tested = run(&job).expect("a result");
        assert!(as_asked(job.asked(), tested.clone()).is_ok());
        // and the other way: a test that was not asked for is not shown as if it were
        assert!(as_asked(SpeciesJob::search(&many()).asked(), tested).is_err());
        assert!(as_asked(SpeciesJob::search(&many()).asked(), plain).is_ok());
    }

    #[test]
    fn an_answer_with_no_bootstrap_to_a_job_that_asked_for_one_is_refused() {
        let plain = run(&SpeciesJob::search(GENES)).expect("a result");
        let mut job = SpeciesJob::search(GENES);
        job.replicates = 5;
        let why = as_asked(job.asked(), plain.clone()).unwrap_err();
        assert!(
            why.contains("bootstrap") && why.contains("older version"),
            "{why}"
        );
        assert!(as_asked(SpeciesJob::search(GENES).asked(), plain).is_ok());
        let found = run(&{
            let mut j = SpeciesJob::search("((A,B),(C,D));\n".repeat(5).as_str());
            j.replicates = 5;
            j
        })
        .expect("a result");
        assert!(as_asked(job.asked(), found.clone()).is_ok());
        // asked for another way of resampling: not what was asked
        job.resample = Resample::Both;
        assert!(as_asked(job.asked(), found).is_err());
    }

    #[test]
    fn a_result_not_weighted_as_asked_is_not_taken_for_one_that_was() {
        let plain = run(&SpeciesJob::search(GENES)).expect("a result");
        let none = Asked {
            weights: Weights::None,
            families: false,
            bootstrap: None,
            polytomy: false,
        };
        assert!(as_asked(none, plain.clone()).is_ok());
        let why = as_asked(
            Asked {
                weights: Weights::Support,
                ..none
            },
            plain,
        )
        .unwrap_err();
        assert!(
            why.contains("\"none\", not \"support\"") && why.contains("older version"),
            "{why}"
        );
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

    /// A reply with these names, rows and branch count in the form a server sends them.
    fn reply(taxa: &str, rows: &str, branch_count: usize) -> Container {
        let found = run(&SpeciesJob::search(GENES)).expect("a species tree");
        let good = found.encode();
        let meta = good
            .text("meta")
            .expect("meta")
            .lines()
            .map(|l| match l.split_once('=') {
                Some(("branch_count", _)) => format!("branch_count={branch_count}\n"),
                _ => format!("{l}\n"),
            })
            .collect::<String>();
        Container::new()
            .with("meta", meta)
            .with("newick", good.text("newick").expect("newick").to_string())
            .with("taxa", taxa.to_string())
            .with("branches", rows.to_string())
    }

    #[test]
    fn a_reply_cannot_make_the_app_hold_more_than_a_table_may_take() {
        let row = |ids: &str| format!("{ids}\t1\t1\t0\t0\t1\t0\t0\t1\t0\t0\t2\t0.5\n");
        // a name too long, more rows than species, a species twice or out of order in a row, one past the end
        let long = format!("{}\nb\nc\nd\n", "x".repeat(MAX_NAME_BYTES + 1));
        assert!(SpeciesResult::decode(&reply(&long, &row("0"), 1)).is_err());
        let four = "a\nb\nc\nd\n";
        assert!(SpeciesResult::decode(&reply(four, &row("0,1"), 1)).is_ok());
        let five_rows = row("0").repeat(5);
        assert!(SpeciesResult::decode(&reply(four, &five_rows, 5)).is_err());
        for bad in ["0,0", "1,0", "0,4", "0,1,1"] {
            assert!(
                SpeciesResult::decode(&reply(four, &row(bad), 1)).is_err(),
                "{bad}"
            );
        }
        // a branch count that is less than the rows, or more than the species there are
        assert!(SpeciesResult::decode(&reply(four, &row("0"), 0)).is_err());
        assert!(SpeciesResult::decode(&reply(four, &row("0"), 5)).is_err());

        // the whole table is bounded, whatever each row costs: 1,500 species, 1,500 rows of 1,000 each is 1.5 million
        // entries against the 1 million that are allowed (a hostile reply that asked for far more is refused the same way)
        let names: String = (0..1_500).map(|i| format!("s{i}\n")).collect();
        let ids = (0..1_000)
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let table = row(&ids).repeat(1_500);
        let why = SpeciesResult::decode(&reply(&names, &table, 1_500)).expect_err("too large");
        assert!(why.contains("too large"), "{why}");
        let fits = row(&ids).repeat(1_000);
        assert!(SpeciesResult::decode(&reply(&names, &fits, 1_500)).is_ok());
    }

    #[test]
    fn a_table_too_large_to_send_is_cut_to_the_rows_that_fit_and_says_how_many_there_were() {
        // a tree as deep as it is wide: 1,500 species, and every branch lists 1,000 of them
        let taxa: Vec<String> = (0..1_500).map(|i| format!("s{i}")).collect();
        let branch = |k: usize| BranchOut {
            taxa: taxa[..1_000].to_vec(),
            f: [1.0, 0.0, 0.0],
            q: [1.0, 0.0, 0.0],
            pp: [1.0, 0.0, 0.0],
            en: 1.0,
            quartets: k as u64,
            length: 1.0,
            bootstrap: None,
            polytomy: None,
        };
        let big = SpeciesResult {
            newick: "(s0,s1,s2);".into(),
            score: 1,
            weights: Weights::None,
            weighted_score: None,
            families: false,
            bootstrap: None,
            polytomy: false,
            normalized: 1.0,
            taxa: taxa.clone(),
            genes: 1,
            skipped: 0,
            clusters: None,
            branches: (0..1_400).map(branch).collect(),
            branch_count: 1_400,
            seconds: 0.0,
        };
        let back = SpeciesResult::decode(&big.encode()).expect("a result that decodes");
        assert_eq!(back.branch_count, 1_400, "the real number of branches");
        assert_eq!(
            back.branches.len(),
            1_000,
            "the rows that fit in a million entries"
        );
        assert_eq!(back.branches[..], big.branches[..1_000]);
        // a result within the limits goes through whole
        let found = run(&SpeciesJob::search(GENES)).expect("a species tree");
        assert_eq!(found.branch_count, found.branches.len());
    }
}
