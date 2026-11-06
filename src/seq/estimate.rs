//! Building a tree from an alignment by the method the person chose: neighbour joining or UPGMA from distances (with
//! bootstrap support), maximum likelihood, or Bayesian sampling. Used by the Phylogenetic Tree window, the Phylogeny
//! page and the workstation server alike.

use super::likelihood::bayes::{bayes_tree, Settings as BayesSettings, CONVERGED_BELOW};
use super::likelihood::ml_tree;
use super::msa::Msa;
use super::tree::{add_bootstrap, msa_distances, DistanceModel, Method, Node};

/// Bootstrap resamples are repeatable: the same seed gives the same support values.
const BOOTSTRAP_SEED: u64 = 2024;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Analysis {
    NeighborJoining,
    Upgma,
    Likelihood,
    Bayesian,
}

impl Analysis {
    /// The name the workstation messages use.
    pub fn name(self) -> &'static str {
        match self {
            Analysis::NeighborJoining => "nj",
            Analysis::Upgma => "upgma",
            Analysis::Likelihood => "ml",
            Analysis::Bayesian => "bayes",
        }
    }

    pub fn from_name(name: &str) -> Option<Analysis> {
        [
            Analysis::NeighborJoining,
            Analysis::Upgma,
            Analysis::Likelihood,
            Analysis::Bayesian,
        ]
        .into_iter()
        .find(|a| a.name() == name)
    }
}

/// What a tree estimation needs: the method and its settings.
#[derive(Clone, Copy)]
pub struct EstimateOptions {
    pub analysis: Analysis,
    pub model: DistanceModel,
    /// Bootstrap replicates for the distance methods; 0 for none.
    pub bootstrap: usize,
    /// Generations of each Bayesian run.
    pub generations: usize,
}

/// A tree and the lines that say how it was made.
pub struct Estimate {
    pub tree: Node,
    pub details: Vec<String>,
}

/// Build the chosen kind of tree from an alignment.
pub fn estimate(msa: &Msa, o: &EstimateOptions) -> Result<Estimate, String> {
    let mut details = Vec::new();
    let method = match o.analysis {
        Analysis::Upgma => Method::Upgma,
        _ => Method::NeighborJoining,
    };
    let tree = if o.analysis == Analysis::Likelihood {
        ml_tree(msa).map(|r| {
            details.extend(r.fits.iter().map(|f| f.describe()));
            details.push(format!("Best model by BIC: {}.", r.best.describe()));
            r.tree
        })
    } else if o.analysis == Analysis::Bayesian {
        let settings = BayesSettings {
            generations: o.generations,
            ..BayesSettings::default()
        };
        bayes_tree(msa, &settings).map(|r| {
            details.push(format!(
                "{} runs of {} generations (first {:.0} % discarded), {} trees sampled; mean ln L {:.2}, mean kappa {:.2}.",
                settings.runs,
                settings.generations,
                settings.burn_in * 100.0,
                r.samples,
                r.mean_ln_l,
                r.mean_kappa
            ));
            details.push(match r.asdsf {
                Some(a) if a < CONVERGED_BELOW => format!(
                    "Average standard deviation of split frequencies {a:.4}: the runs agree."
                ),
                Some(a) => format!(
                    "Average standard deviation of split frequencies {a:.4} is above {CONVERGED_BELOW}: run more generations."
                ),
                None => "A single run: convergence can't be checked.".to_string(),
            });
            details.push("Branch labels: posterior probability, %. Majority-rule consensus.".into());
            r.tree
        })
    } else {
        msa_distances(msa, o.model).and_then(|d| {
            let mut tree = method
                .tree(&msa.names, &d)
                .ok_or("The distances can't make a tree.")?;
            if o.bootstrap > 0 {
                add_bootstrap(&mut tree, msa, o.model, method, o.bootstrap, BOOTSTRAP_SEED);
                details.push(format!(
                    "Branch labels: bootstrap support, % of {} resamples of the columns.",
                    o.bootstrap
                ));
            }
            Ok(tree)
        })
    }?;
    Ok(Estimate { tree, details })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn analyses_have_names_that_read_back() {
        for a in [
            Analysis::NeighborJoining,
            Analysis::Upgma,
            Analysis::Likelihood,
            Analysis::Bayesian,
        ] {
            assert_eq!(Analysis::from_name(a.name()), Some(a));
        }
        assert_eq!(Analysis::from_name("parsimony"), None);
    }
}
