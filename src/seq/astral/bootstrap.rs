//! Bootstrap support for a species tree, as ASTRAL gives it (`docs/SPECIES.md`): the search is run again on many sets of
//! gene trees made from the genes by resampling, and the support of a branch of the main species tree is the share of those
//! species trees (in percent) that have the branch.
//!
//! Three ways to make a replicate:
//!
//! * **genes**: the genes are drawn with replacement, as many as there are, from the gene trees;
//! * **sites**: every gene has bootstrap trees of its own (from bootstrapping its alignment); replicate `r` takes the `r`-th
//!   of each gene's, so each gene needs at least as many as there are replicates;
//! * **both** (the multilocus bootstrap of Seo 2008): the genes are drawn with replacement, and for every gene drawn one of
//!   its bootstrap trees is drawn.
//!
//! The source of random numbers is seeded, so the same input and seed give the same support every time.
//!
//! The replicates keep the weights of a weighted run: a gene drawn takes its supports and lengths with it. Bootstrap trees of the
//! sites have none, so the sites cannot be resampled with weights.

use super::bits::Bits;
use super::genes::{check_size, leaf_names, GTree, Gene, Genes};
use super::report;
use super::search;
use super::species::Sp;
use crate::seq::newick::NwTree;
use std::collections::HashMap;

/// The most replicates a run may ask for: every one is a whole search.
pub const MAX_REPLICATES: usize = 1_000;
/// The seed ASTRAL uses when none is given.
pub const DEFAULT_SEED: u64 = 692;

/// What the replicates are made of.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Resample {
    /// The genes, drawn with replacement.
    #[default]
    Genes,
    /// The bootstrap trees of the genes, one of each for every replicate.
    Sites,
    /// Genes drawn with replacement and, for each, one of its bootstrap trees.
    Both,
}

impl Resample {
    pub const ALL: [Resample; 3] = [Resample::Genes, Resample::Sites, Resample::Both];

    pub fn name(self) -> &'static str {
        match self {
            Resample::Genes => "genes",
            Resample::Sites => "sites",
            Resample::Both => "both",
        }
    }

    pub fn parse(word: &str) -> Option<Resample> {
        let word = word.trim().to_ascii_lowercase();
        Resample::ALL.into_iter().find(|r| r.name() == word)
    }

    /// Whether the gene trees' own bootstrap trees are used.
    pub fn uses_site_trees(self) -> bool {
        self != Resample::Genes
    }
}

/// What to ask of the bootstrap: no replicates (the default) is no bootstrap.
#[derive(Clone, Debug, PartialEq)]
pub struct Boot {
    pub replicates: usize,
    pub seed: u64,
    pub resample: Resample,
}

impl Default for Boot {
    fn default() -> Self {
        Boot {
            replicates: 0,
            seed: DEFAULT_SEED,
            resample: Resample::Genes,
        }
    }
}

/// What a bootstrap did, to say with the result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BootInfo {
    pub resample: Resample,
    pub asked: usize,
    pub done: usize,
    /// Replicates that gave no species tree (a species that the genes drawn do not have at all can leave no tree to build).
    pub failed: usize,
}

/// The bootstrap trees of every gene, read on the same species as the gene trees.
pub struct SiteTrees {
    per_gene: Vec<Vec<(GTree, Gene)>>,
}

impl SiteTrees {
    /// `texts` are the bootstrap trees of each gene tree of the input, in its order (one text for each tree of the input,
    /// the ones that were left out as too small included), each tree on species that the gene trees have.
    pub fn read(genes: &Genes, texts: &[String]) -> Result<SiteTrees, String> {
        if texts.len() != genes.input_trees {
            return Err(format!(
                "There are bootstrap trees for {} genes, and {} gene trees: each gene needs a file of bootstrap trees, in the order of the gene trees.",
                texts.len(),
                genes.input_trees
            ));
        }
        let n = genes.taxa.len();
        let index: HashMap<&str, usize> = genes
            .taxa
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        let mut per_gene = Vec::with_capacity(genes.origin.len());
        let (mut trees, mut nodes) = (0usize, 0u128);
        for &k in &genes.origin {
            let parsed = NwTree::parse_all(&texts[k - 1])
                .map_err(|why| format!("The bootstrap trees of gene {k}: {why}"))?;
            let mut list = Vec::with_capacity(parsed.len());
            for t in &parsed {
                for name in leaf_names(t, k)? {
                    if !index.contains_key(name) {
                        return Err(format!(
                            "A bootstrap tree of gene {k} has the species {name}, which the gene trees do not have."
                        ));
                    }
                }
                trees += 1;
                nodes += t.nodes.len() as u128;
                let tree = GTree::from_newick(t, &index);
                let gene = tree.gene(n);
                list.push((tree, gene));
            }
            per_gene.push(list);
        }
        check_size(trees, nodes, n, 8 << 30)?;
        Ok(SiteTrees { per_gene })
    }
}

/// The most files a list of bootstrap files may name, and the most bytes of trees in all: a file for each gene, as many genes
/// as a species-tree job may hold.
const MAX_LISTED_FILES: usize = 20_000;
const MAX_LISTED_BYTES: u64 = 512 << 20;

/// The bootstrap trees of the genes from a list of files, as ASTRAL's `-b` takes them: one path on each line (blank lines and
/// lines that start with `#` are left out), a relative path being taken from the folder of the list, each file holding the
/// bootstrap trees of one gene, in the order of the gene trees. The text of each file, in the order of the list.
pub fn read_list(list: &std::path::Path) -> Result<Vec<String>, String> {
    let text = std::fs::read_to_string(list)
        .map_err(|e| format!("Could not read the list {}: {e}", list.display()))?;
    let folder = list.parent().unwrap_or(std::path::Path::new("."));
    let paths: Vec<std::path::PathBuf> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(|l| folder.join(l))
        .collect();
    if paths.is_empty() {
        return Err(format!("The list {} names no file.", list.display()));
    }
    if paths.len() > MAX_LISTED_FILES {
        return Err(format!(
            "The list {} names {} files; at most {MAX_LISTED_FILES} are taken.",
            list.display(),
            paths.len()
        ));
    }
    let mut total = 0u64;
    paths
        .iter()
        .map(|p| {
            let size = std::fs::metadata(p)
                .map_err(|e| format!("Could not read {}: {e}", p.display()))?
                .len();
            total += size;
            if total > MAX_LISTED_BYTES {
                return Err(format!(
                    "The bootstrap files hold more than {} MB of trees: use fewer genes or fewer trees for each.",
                    MAX_LISTED_BYTES >> 20
                ));
            }
            std::fs::read_to_string(p).map_err(|e| format!("Could not read {}: {e}", p.display()))
        })
        .collect()
}

/// A source of random numbers that is the same on every run.
struct Random(u64);

impl Random {
    fn below(&mut self, n: usize) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 33) % n as u64) as usize
    }
}

/// The support of the branches of the main tree, by the taxa on the side that does not have the root taxon.
pub struct Support {
    pub of: HashMap<Bits, f64>,
    pub info: BootInfo,
    /// When asked for: the species tree of every replicate that gave one, as Newick, and its inner splits.
    pub trees: Vec<String>,
    pub splits: Vec<Vec<Bits>>,
}

/// The inner branches of a species tree: the taxa under every node that has children, but the root.
pub fn branches(sp: &Sp) -> impl Iterator<Item = &Bits> {
    sp.nodes
        .iter()
        .enumerate()
        .filter(|&(i, n)| n.kids.is_some() && i != sp.root)
        .map(|(_, n)| &n.set)
}

/// Run `b.replicates` searches on resamplings of `genes` (with the options `o` of the main search) and count how often each
/// branch of `main` is in the tree they find. With `keep`, the tree of every replicate is kept too.
pub fn run(
    genes: &Genes,
    sites: Option<&SiteTrees>,
    o: &search::Options,
    main: &Sp,
    b: &Boot,
    keep: bool,
) -> Result<Support, String> {
    if b.replicates == 0 || b.replicates > MAX_REPLICATES {
        return Err(format!(
            "The replicates must be from 1 to {MAX_REPLICATES}, not {}.",
            b.replicates
        ));
    }
    let n = genes.trees.len();
    if !genes.wtrees.is_empty() && b.resample.uses_site_trees() {
        return Err("Weights and a bootstrap of the sites cannot be used together: the bootstrap trees have no supports or lengths to weigh by.".into());
    }
    let site_trees = match (b.resample.uses_site_trees(), sites) {
        (false, _) => None,
        (true, Some(s)) => Some(s),
        (true, None) => {
            return Err("Resampling the sites needs the bootstrap trees of the genes.".into())
        }
    };
    if let (Resample::Sites, Some(s)) = (b.resample, site_trees) {
        if let Some((g, list)) = s
            .per_gene
            .iter()
            .enumerate()
            .find(|(_, l)| l.len() < b.replicates)
        {
            return Err(format!(
                "Gene {} has {} bootstrap trees, and {} replicates need that many each.",
                genes.origin[g],
                list.len(),
                b.replicates
            ));
        }
    }
    if let Some(s) = site_trees {
        if let Some((g, _)) = s.per_gene.iter().enumerate().find(|(_, l)| l.is_empty()) {
            return Err(format!("Gene {} has no bootstrap trees.", genes.origin[g]));
        }
    }
    let mut counts: HashMap<Bits, usize> = branches(main).map(|s| (s.clone(), 0)).collect();
    let mut random = Random(b.seed);
    let (mut done, mut failed) = (0usize, 0usize);
    let (mut kept_trees, mut kept_splits) = (Vec::new(), Vec::new());
    for r in 0..b.replicates {
        // the genes drawn, when they are drawn from the gene trees themselves
        let mut picked: Vec<usize> = Vec::new();
        let (trees, gs): (Vec<GTree>, Vec<Gene>) = match (b.resample, site_trees) {
            (Resample::Sites, Some(s)) => (0..n).map(|g| s.per_gene[g][r].clone()).unzip(),
            (Resample::Both, Some(s)) => (0..n)
                .map(|_| {
                    let list = &s.per_gene[random.below(n)];
                    list[random.below(list.len())].clone()
                })
                .unzip(),
            _ => {
                picked = (0..n).map(|_| random.below(n)).collect();
                picked
                    .iter()
                    .map(|&g| (genes.trees[g].clone(), genes.genes[g].clone()))
                    .unzip()
            }
        };
        let wtrees = if genes.wtrees.is_empty() {
            Vec::new()
        } else {
            picked.iter().map(|&g| genes.wtrees[g].clone()).collect()
        };
        let replicate = Genes {
            taxa: genes.taxa.clone(),
            trees,
            genes: gs,
            wtrees,
            notes: Vec::new(),
            skipped: 0,
            input_trees: n,
            origin: Vec::new(),
        };
        match search::search(&replicate, o) {
            Ok(found) => {
                done += 1;
                for set in branches(&found.sp) {
                    if let Some(c) = counts.get_mut(set) {
                        *c += 1;
                    }
                }
                if keep {
                    kept_trees.push(report::newick(&found.sp, &genes.taxa, &[], 0));
                    kept_splits.push(branches(&found.sp).cloned().collect());
                }
            }
            Err(_) => failed += 1,
        }
    }
    if done == 0 {
        return Err("None of the bootstrap replicates gave a species tree: the genes drawn left some species with nothing to place them by.".into());
    }
    Ok(Support {
        of: counts
            .into_iter()
            .map(|(set, c)| (set, 100.0 * c as f64 / done as f64))
            .collect(),
        info: BootInfo {
            resample: b.resample,
            asked: b.replicates,
            done,
            failed,
        },
        trees: kept_trees,
        splits: kept_splits,
    })
}
