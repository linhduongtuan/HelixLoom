//! Species trees from gene trees, in the way of ASTRAL (Mirarab et al. 2014; Zhang et al. 2018) and the ASTER family
//! (Zhang, Nielsen and Mirarab 2025): the species tree is the one whose quartets agree with the most quartets of the gene
//! trees, which under the multispecies coalescent is a consistent estimate of the true species tree despite incomplete
//! lineage sorting. Each branch also gets a length in coalescent units and a local posterior probability.
//!
//! This is written from the published descriptions (the papers and the tutorials), not from the programs' code, and it is
//! checked against ASTRAL 5.7.8 as a black box: see `tests`. What it does and does not cover is listed in MECHANICS.md.

pub mod bits;
pub mod bootstrap;
pub mod cli;
pub mod complete;
pub mod consensus;
pub mod family;
pub mod genes;
pub mod polytomy;
pub mod pscore;
pub mod report;
pub mod search;
pub mod special;
pub mod species;
pub mod weight;
pub mod wgenes;
pub mod wscore;

pub use bootstrap::{Boot, BootInfo, Resample};
pub use polytomy::Polytomy;
use species::{Branch, Sp};
pub use wgenes::Weights;

#[cfg(test)]
mod btests;
#[cfg(test)]
mod ctests;
#[cfg(test)]
mod ftests;
#[cfg(test)]
mod golden_tests;
#[cfg(test)]
mod tests;
#[cfg(test)]
mod wtests;

/// What to ask of a run. The defaults are ASTRAL's.
#[derive(Clone, Debug)]
pub struct Options {
    /// Bootstrap replicates for the support of the branches (none by default; `bootstrap.rs`).
    pub bootstrap: Boot,
    /// Keep the species tree of every bootstrap replicate and write their greedy consensus (needs replicates).
    pub keep_replicates: bool,
    /// The bootstrap trees of each gene tree of the input (one text of trees for each, in the order of the gene trees), for the
    /// replicates that resample the sites.
    pub site_trees: Vec<String>,
    /// The gene trees are gene family trees (ASTRAL-Pro, Zhang et al. 2020): a species may occur in a tree more than once, the
    /// leaf labels being the species names.
    pub families: bool,
    /// How the gene trees' quartets are weighted (the weighted ASTRAL of Zhang and Mirarab 2022); the default counts each one.
    pub weights: Weights,
    /// Look at every group of taxa instead of the enriched set from the gene trees (at most 12 taxa).
    pub exact: bool,
    /// How much to add to the set of groups: 0 nothing, 1 (the default) the ways to resolve polytomies.
    pub extra: u8,
    /// Most children of a polytomy whose pairs are joined when `extra` is 1.
    pub polylimit: usize,
    /// Swap neighbouring groups after the search while the score rises.
    pub polish: bool,
    /// The prior rate of the Yule model behind the branch lengths; 0 gives the maximum likelihood length.
    pub lambda: f64,
    /// What each branch carries in the Newick text: ASTRAL's `-t` values 0, 1, 2, 3, 4 and 8.
    pub annotate: u8,
    /// Test every branch for being a polytomy (Sayyari and Mirarab 2018; ASTRAL's `-t 10`): the p-value is the label of the
    /// branch unless there is a bootstrap support or labels are off.
    pub polytomy: bool,
    /// The taxon to hang the tree from when it is written.
    pub root: Option<String>,
    /// Gene trees with fewer leaves than this are left out.
    pub min_leaves: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            bootstrap: Boot::default(),
            keep_replicates: false,
            site_trees: Vec::new(),
            families: false,
            weights: Weights::None,
            exact: false,
            extra: 1,
            polylimit: complete::POLYLIMIT,
            polish: true,
            lambda: 0.5,
            annotate: 3,
            polytomy: false,
            root: None,
            min_leaves: 0,
        }
    }
}

/// What the gene trees say about one inner branch, with the taxa on one side of it.
#[derive(Clone, Debug, PartialEq)]
pub struct BranchOut {
    pub taxa: Vec<String>,
    pub f: [f64; 3],
    pub q: [f64; 3],
    pub pp: [f64; 3],
    pub en: f64,
    pub quartets: u64,
    pub length: f64,
    /// The share of the bootstrap replicates (in percent) that have the branch, when there were any.
    pub bootstrap: Option<f64>,
    /// The test for a polytomy, when it was asked for and at least `polytomy::MIN_GENES` gene trees count.
    pub polytomy: Option<Polytomy>,
}

#[derive(Clone, Debug)]
pub struct Output {
    /// The species tree in Newick.
    pub newick: String,
    /// The quartets of the species tree that the gene trees share, summed over the gene trees (each counted one, whatever the
    /// weights of the run).
    pub score: u64,
    /// With weights: the same quartets with their weights added up.
    pub weighted_score: Option<f64>,
    /// How the quartets were weighted.
    pub weights: Weights,
    /// The gene trees were read as gene families (ASTRAL-Pro), whose score counts classes of quartets, and which has no branch
    /// annotations.
    pub families: bool,
    /// What the bootstrap did, when one was asked for.
    pub bootstrap: Option<BootInfo>,
    /// When asked for: the species tree of each bootstrap replicate that gave one, and the greedy consensus of them.
    pub replicate_trees: Vec<String>,
    pub consensus: Option<String>,
    /// The branches were tested for being polytomies.
    pub polytomy: bool,
    /// The score (the weighted one in a weighted run) as a share of the quartets the gene trees resolve.
    pub normalized: f64,
    pub taxa: Vec<String>,
    pub genes: usize,
    /// Gene trees left out for being too small.
    pub skipped: usize,
    /// How many groups the search could choose from (none when a tree was only scored).
    pub clusters: Option<usize>,
    pub branches: Vec<BranchOut>,
    pub log: Vec<String>,
}

/// The taxon the tree is hung from when it is written: the one named, which the gene trees must have, else the first.
fn root_taxon(taxa: &[String], o: &Options) -> Result<usize, String> {
    match &o.root {
        Some(name) => taxa
            .iter()
            .position(|t| t == name)
            .ok_or_else(|| format!("The root taxon {name} is not in the gene trees.")),
        None => Ok(0),
    }
}

/// Find the species tree for gene trees given as Newick text, one tree after the other.
pub fn infer(gene_trees: &str, o: &Options) -> Result<Output, String> {
    let boot = o.bootstrap.replicates > 0;
    if o.families {
        if boot {
            return Err("Bootstrap replicates and gene families cannot be used together: a bootstrap of gene families is not built.".into());
        }
        return infer_families(gene_trees, o);
    }
    if boot && o.weights != Weights::None && o.bootstrap.resample.uses_site_trees() {
        return Err("Bootstrap replicates that resample the sites and weights cannot be used together: the bootstrap trees have no supports or lengths to weigh by (resample the genes instead).".into());
    }
    if o.keep_replicates && !boot {
        return Err(
            "The replicate trees are those of a bootstrap: ask for replicates as well.".into(),
        );
    }
    let genes = genes::read_weighted(gene_trees, o.min_leaves, o.weights)?;
    let root = root_taxon(&genes.taxa, o)?;
    let so = search::Options {
        weights: o.weights,
        exact: o.exact,
        extra: o.extra,
        polylimit: o.polylimit,
        polish: o.polish,
        root,
    };
    // the bootstrap trees are read, and the settings checked, before the search that may take long
    let sites = (boot && o.bootstrap.resample.uses_site_trees())
        .then(|| bootstrap::SiteTrees::read(&genes, &o.site_trees))
        .transpose()?;
    let found = search::search(&genes, &so)?;
    let support = boot
        .then(|| {
            bootstrap::run(
                &genes,
                sites.as_ref(),
                &so,
                &found.sp,
                &o.bootstrap,
                o.keep_replicates,
            )
        })
        .transpose()?;
    let mut out = finish(
        &genes,
        &found.sp,
        (found.score, found.weighted),
        Some(found.clusters),
        support.as_ref(),
        o,
    );
    if !found.completed {
        out.log.push(format!(
            "More than {} species, so no guide tree: the gene trees that lack species were not completed, and the search chose from the groups the gene trees give.",
            complete::GUIDE_LIMIT
        ));
    }
    Ok(out)
}

/// Score a given binary species tree against gene trees and annotate its branches (ASTRAL's `-q`).
pub fn score_tree(gene_trees: &str, species_tree: &str, o: &Options) -> Result<Output, String> {
    if o.families {
        return score_tree_families(gene_trees, species_tree, o);
    }
    let genes = genes::read_weighted(gene_trees, o.min_leaves, o.weights)?;
    let sp = Sp::read_hung(species_tree, &genes.taxa, root_taxon(&genes.taxa, o)?)?;
    if o.bootstrap.replicates > 0 {
        return Err(
            "Bootstrap replicates are for a search: a tree that is given is only scored.".into(),
        );
    }
    let score = sp.score(&genes);
    let weighted = (o.weights != Weights::None)
        .then(|| sp.score_with(|a, b| wscore::tripartition(&genes.wtrees, a, b, o.weights)));
    Ok(finish(&genes, &sp, (score, weighted), None, None, o))
}

/// What a run on gene families may not be asked for, and the families read.
fn families_for(gene_trees: &str, o: &Options) -> Result<family::Families, String> {
    if o.polytomy {
        return Err("The polytomy test and gene families cannot be used together: the species tree of gene families has no table of quartet counts to test.".into());
    }
    if o.weights != Weights::None {
        return Err(
            "Weights and gene families cannot be used together: weighted ASTRAL-Pro is not built."
                .into(),
        );
    }
    family::read(gene_trees, o.min_leaves)
}

/// Find the species tree for gene family trees (ASTRAL-Pro): the leaf labels are species names, which may repeat.
fn infer_families(gene_trees: &str, o: &Options) -> Result<Output, String> {
    let f = families_for(gene_trees, o)?;
    let root = root_taxon(&f.taxa, o)?;
    let found = search::search_families(
        &f,
        &search::Options {
            weights: Weights::None,
            exact: o.exact,
            extra: o.extra,
            polylimit: o.polylimit,
            polish: o.polish,
            root,
        },
    )?;
    Ok(finish_families(
        &f,
        &found.sp,
        found.score,
        Some(found.clusters),
    ))
}

/// Score a given binary species tree against gene family trees.
fn score_tree_families(
    gene_trees: &str,
    species_tree: &str,
    o: &Options,
) -> Result<Output, String> {
    let f = families_for(gene_trees, o)?;
    let sp = Sp::read_hung(species_tree, &f.taxa, root_taxon(&f.taxa, o)?)?;
    let score = sp.sum_with(|a, b| pscore::tripartition(&f.nodes, a, b));
    Ok(finish_families(&f, &sp, score, None))
}

fn finish_families(f: &family::Families, sp: &Sp, score: u64, clusters: Option<usize>) -> Output {
    let normalized = if f.classes == 0 {
        0.0
    } else {
        score as f64 / f.classes as f64
    };
    let mut first = format!(
        "{} gene family trees on {} species, {} duplications",
        f.families,
        f.taxa.len(),
        f.duplications
    );
    if f.rooted_by_cost > 0 {
        first.push_str(&format!(
            "; {} rooted where duplications and losses cost least",
            f.rooted_by_cost
        ));
    }
    if f.skipped > 0 {
        first.push_str(&format!("; {} with too few species, left out", f.skipped));
    }
    Output {
        newick: report::newick(sp, &f.taxa, &[], 0),
        score,
        weighted_score: None,
        weights: Weights::None,
        families: true,
        bootstrap: None,
        replicate_trees: Vec::new(),
        consensus: None,
        polytomy: false,
        normalized,
        taxa: f.taxa.clone(),
        genes: f.families,
        skipped: f.skipped,
        clusters,
        branches: Vec::new(),
        log: vec![
            first,
            format!("Final quartet score is: {score}"),
            format!("Final normalized quartet score is: {normalized}"),
            "The score is ASTRAL-Pro's: the quartet classes of the speciation nodes of the family trees (Zhang et al. 2020). The branch support and the lengths of the species tree are not given for gene families: ASTRAL-Pro's are not built.".into(),
        ],
    }
}

/// `scores` are the quartets of the tree that the gene trees share, and their weights in a weighted run.
fn finish(
    genes: &genes::Genes,
    sp: &Sp,
    scores: (u64, Option<f64>),
    clusters: Option<usize>,
    support: Option<&bootstrap::Support>,
    o: &Options,
) -> Output {
    let (score, weighted) = scores;
    let resolved: u128 = genes.genes.iter().map(genes::Gene::resolved_quartets).sum();
    let normalized = if resolved == 0 {
        0.0
    } else {
        weighted.unwrap_or(score as f64) / resolved as f64
    };
    let mut branches: Vec<Branch> = sp.annotate(genes, o.lambda);
    if let Some(s) = support {
        for b in &mut branches {
            b.bootstrap = s.of.get(&sp.nodes[b.node].set).copied();
        }
    }
    if o.polytomy {
        let tests = sp.polytomy_tests(genes);
        for b in &mut branches {
            b.polytomy = tests.iter().find(|t| t.0 == b.node).and_then(|t| t.1);
        }
    }
    let level = if o.polytomy && o.annotate != 0 {
        10
    } else {
        o.annotate
    };
    let newick = report::newick(sp, &genes.taxa, &branches, level);
    let named = |b: &Branch| BranchOut {
        taxa: sp.nodes[b.node]
            .set
            .iter()
            .map(|t| genes.taxa[t].clone())
            .collect(),
        f: b.f,
        q: b.q,
        pp: b.pp,
        en: b.en,
        quartets: b.quartets,
        length: b.length,
        bootstrap: b.bootstrap,
        polytomy: b.polytomy,
    };
    let mut log = vec![
        format!(
            "{} gene trees on {} taxa{}",
            genes.genes.len(),
            genes.taxa.len(),
            if genes.skipped > 0 {
                format!(" ({} too small, left out)", genes.skipped)
            } else {
                String::new()
            }
        ),
        format!("Final quartet score is: {score}"),
    ];
    if let Some(w) = weighted {
        log.push(format!("Weights: {}", o.weights.describe()));
        log.extend(genes.notes.iter().cloned());
        log.push(format!("Final weighted quartet score is: {w:.6}"));
        log.push(format!(
            "The branch support below is ASTRAL's, counted on the quartets without their weights: the weighted support of wastral is not built.{}",
            if o.polytomy { " The polytomy test is counted without the weights too." } else { "" }
        ));
    }
    log.push(format!("Final normalized quartet score is: {normalized}"));
    if let Some(s) = support {
        let i = s.info;
        log.push(format!(
            "Bootstrap: {} replicates resampling the {}, seed {}; {} gave a species tree{}. The numbers on the branches are the share of those trees (in percent) that have the branch.",
            i.asked,
            match i.resample {
                Resample::Genes => "genes",
                Resample::Sites => "sites (one bootstrap tree of every gene for each replicate)",
                Resample::Both => "genes and their sites (genes drawn with replacement, and one bootstrap tree of each)",
            },
            o.bootstrap.seed,
            i.done,
            if i.failed > 0 {
                format!(", {} did not", i.failed)
            } else {
                String::new()
            }
        ));
    }
    if o.polytomy {
        log.push(polytomy_line(&branches, support.is_some(), o.annotate));
    }
    Output {
        newick,
        score,
        weighted_score: weighted,
        weights: o.weights,
        families: false,
        bootstrap: support.map(|s| s.info),
        replicate_trees: support
            .filter(|_| o.keep_replicates)
            .map_or_else(Vec::new, |s| s.trees.clone()),
        consensus: support.filter(|_| o.keep_replicates).map(|s| {
            consensus::newick(
                &consensus::greedy(&s.splits),
                s.info.done,
                &genes.taxa,
                sp.root_taxon,
            )
        }),
        polytomy: o.polytomy,
        normalized,
        taxa: genes.taxa.clone(),
        genes: genes.genes.len(),
        skipped: genes.skipped,
        clusters,
        branches: branches.iter().map(named).collect(),
        log,
    }
}

/// What the test for a polytomy found, in a line for the log: how many branches were tested, how many it did not reject, and
/// what the numbers on the branches are.
fn polytomy_line(branches: &[Branch], bootstrapped: bool, annotate: u8) -> String {
    let tested: Vec<&Polytomy> = branches
        .iter()
        .filter_map(|b| b.polytomy.as_ref())
        .collect();
    let kept = tested.iter().filter(|t| !t.rejected()).count();
    let label = match (annotate, bootstrapped) {
        (0, _) => "No branch is labelled (-t 0).",
        (_, true) => "The numbers on the branches are the bootstrap support; the p-values are in the table of branches.",
        _ => "The numbers on the branches are the p-values; a branch with only its length was not tested.",
    };
    format!(
        "Polytomy test (Sayyari and Mirarab 2018): {} of {} branches tested, each with at least {} gene trees that resolve a quartet around it; {} of them cannot be told from a polytomy at p of {} or more. {label}",
        tested.len(),
        branches.len(),
        polytomy::MIN_GENES,
        kept,
        polytomy::ALPHA,
    )
}
