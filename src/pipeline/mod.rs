//! Workflows: chains of the program's tools (alignment, tree, primers), run with the installed MAFFT, IQ-TREE and
//! primer3 when they are there and with the built-in engines when they are not, and exported as Snakemake and
//! Nextflow files for running elsewhere.

pub mod export;
pub mod predict;
pub mod steps;
pub mod tools;

use crate::seq::primer3::Settings;
use steps::{Aligned, Named, Primers, Tree};
use tools::Installed;

/// The installed programs the workflow may call (`None`: use the built-in engine).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Available {
    pub mafft: Option<Installed>,
    pub iqtree: Option<Installed>,
    pub primer3: Option<Installed>,
}

/// Which steps to run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Choice {
    pub align: bool,
    pub tree: bool,
    pub primers: bool,
    /// Use the installed programs where there are some.
    pub use_installed: bool,
    pub bootstrap: usize,
}

#[derive(Default)]
pub struct Outcome {
    pub log: Vec<String>,
    pub aligned: Option<Aligned>,
    pub tree: Option<Tree>,
    pub primers: Option<Primers>,
}

/// Run the chosen steps on `seqs` in order: a tree needs the alignment, so choosing a tree aligns first.
/// Primers are picked on the first sequence around `target`.
pub fn run(
    seqs: &[Named],
    target: Option<(usize, usize)>,
    settings: &Settings,
    choice: &Choice,
    available: &Available,
) -> Result<Outcome, String> {
    let dir = tools::work_folder("workflow")?;
    let pick = |found: &Option<Installed>| found.clone().filter(|_| choice.use_installed);
    let mut out = Outcome::default();
    if choice.align || choice.tree {
        let aligned = steps::align(seqs, pick(&available.mafft).as_ref(), &dir)?;
        out.log.push(format!(
            "Align: {} sequences into {} columns with {} in {:.2} s.",
            aligned.rows.len(),
            aligned.columns(),
            aligned.engine,
            aligned.seconds
        ));
        out.aligned = Some(aligned);
    }
    if choice.tree {
        let aligned = out
            .aligned
            .as_ref()
            .ok_or("There is no alignment to build a tree from.")?;
        let tree = steps::tree(
            aligned,
            pick(&available.iqtree).as_ref(),
            choice.bootstrap,
            &dir,
        )?;
        out.log.push(format!(
            "Tree: {}{}{} in {:.2} s.",
            tree.engine,
            if tree.model.is_empty() {
                String::new()
            } else {
                format!(", model {}", tree.model)
            },
            tree.ln_likelihood
                .map_or(String::new(), |l| format!(", ln L {l:.2}")),
            tree.seconds
        ));
        out.tree = Some(tree);
    }
    if choice.primers {
        let first = seqs
            .first()
            .ok_or("There is no sequence to pick primers on.")?;
        let primers = steps::primers(
            &first.1,
            target,
            settings,
            pick(&available.primer3).as_ref(),
            &dir,
        )?;
        out.log.push(format!(
            "Primers: {} pairs on {} with {} in {:.2} s.",
            primers.pairs.len(),
            first.0,
            primers.engine,
            primers.seconds
        ));
        out.primers = Some(primers);
    }
    Ok(out)
}

#[cfg(test)]
mod tests;
