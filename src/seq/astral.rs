//! Species trees from gene trees, in the way of ASTRAL (Mirarab et al. 2014; Zhang et al. 2018) and the ASTER family
//! (Zhang, Nielsen and Mirarab 2025): the species tree is the one whose quartets agree with the most quartets of the gene
//! trees, which under the multispecies coalescent is a consistent estimate of the true species tree despite incomplete
//! lineage sorting. Each branch also gets a length in coalescent units and a local posterior probability.
//!
//! This is written from the published descriptions (the papers and the tutorials), not from the programs' code, and it is
//! checked against ASTRAL 5.7.8 as a black box: see `tests`. What it does and does not cover is listed in MECHANICS.md.

pub mod bits;
pub mod cli;
pub mod complete;
pub mod genes;
pub mod report;
pub mod search;
pub mod special;
pub mod species;
pub mod weight;

use species::{Branch, Sp};

#[cfg(test)]
mod golden_tests;
#[cfg(test)]
mod tests;

/// What to ask of a run. The defaults are ASTRAL's.
#[derive(Clone, Debug)]
pub struct Options {
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
    /// The taxon to hang the tree from when it is written.
    pub root: Option<String>,
    /// Gene trees with fewer leaves than this are left out.
    pub min_leaves: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            exact: false,
            extra: 1,
            polylimit: complete::POLYLIMIT,
            polish: true,
            lambda: 0.5,
            annotate: 3,
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
}

#[derive(Clone, Debug)]
pub struct Output {
    /// The species tree in Newick.
    pub newick: String,
    /// The quartets of the species tree that the gene trees share, summed over the gene trees.
    pub score: u64,
    /// The score as a share of the quartets the gene trees resolve.
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
fn root_taxon(genes: &genes::Genes, o: &Options) -> Result<usize, String> {
    match &o.root {
        Some(name) => genes
            .taxa
            .iter()
            .position(|t| t == name)
            .ok_or_else(|| format!("The root taxon {name} is not in the gene trees.")),
        None => Ok(0),
    }
}

/// Find the species tree for gene trees given as Newick text, one tree after the other.
pub fn infer(gene_trees: &str, o: &Options) -> Result<Output, String> {
    let genes = genes::read(gene_trees, o.min_leaves)?;
    let root = root_taxon(&genes, o)?;
    let found = search::search(
        &genes,
        &search::Options {
            exact: o.exact,
            extra: o.extra,
            polylimit: o.polylimit,
            polish: o.polish,
            root,
        },
    )?;
    let mut out = finish(&genes, &found.sp, found.score, Some(found.clusters), o);
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
    let genes = genes::read(gene_trees, o.min_leaves)?;
    let sp = Sp::read_hung(species_tree, &genes.taxa, root_taxon(&genes, o)?)?;
    let score = sp.score(&genes);
    Ok(finish(&genes, &sp, score, None, o))
}

fn finish(
    genes: &genes::Genes,
    sp: &Sp,
    score: u64,
    clusters: Option<usize>,
    o: &Options,
) -> Output {
    let resolved: u128 = genes.genes.iter().map(genes::Gene::resolved_quartets).sum();
    let normalized = if resolved == 0 {
        0.0
    } else {
        score as f64 / resolved as f64
    };
    let branches: Vec<Branch> = sp.annotate(genes, o.lambda);
    let newick = report::newick(sp, &genes.taxa, &branches, o.annotate);
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
    };
    let log = vec![
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
        format!("Final normalized quartet score is: {normalized}"),
    ];
    Output {
        newick,
        score,
        normalized,
        taxa: genes.taxa.clone(),
        genes: genes.genes.len(),
        skipped: genes.skipped,
        clusters,
        branches: branches.iter().map(named).collect(),
        log,
    }
}
