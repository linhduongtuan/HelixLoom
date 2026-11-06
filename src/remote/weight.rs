//! How big a job is, for the choice of where it runs: with the choice set to Automatic a big job goes to the workstation
//! and a small one stays on this computer, where there is nothing to upload.

use crate::pipeline::Choice;
use crate::seq::estimate::Analysis;

/// A job is big once it passes any of these. They are estimates of where a job stops being quick, not measurements.
pub const HEAVY_SEQUENCES: usize = 100;
pub const HEAVY_LETTERS: usize = 1_000_000;
pub const HEAVY_REPLICATES: usize = 200;
pub const HEAVY_GENERATIONS: usize = 100_000;
pub const HEAVY_TREES: usize = 100;

/// How big a job is: the things that make it slow, counted from what a page has set up (before there is a request to
/// count). A new kind of analysis adds the counts that make it slow.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Weight {
    pub sequences: usize,
    /// Letters of sequence in all.
    pub letters: usize,
    /// Bootstrap resamples.
    pub replicates: usize,
    /// Bayesian generations.
    pub generations: usize,
    pub trees: usize,
    /// The job would run on a GPU of the workstation, whatever its size.
    pub gpu: bool,
}

impl Weight {
    /// `count` sequences of `letters` letters in all, to align.
    pub fn of_sequences(count: usize, letters: usize) -> Weight {
        Weight {
            sequences: count,
            letters,
            ..Weight::default()
        }
    }

    /// A tree from an alignment of `rows` rows and `columns` columns. Only the distance methods bootstrap and only the
    /// Bayesian one runs generations; the other settings are not used by the method chosen.
    pub fn of_tree(
        rows: usize,
        columns: usize,
        analysis: Analysis,
        bootstrap: usize,
        generations: usize,
    ) -> Weight {
        let distance = matches!(analysis, Analysis::NeighborJoining | Analysis::Upgma);
        Weight {
            replicates: if distance { bootstrap } else { 0 },
            generations: if analysis == Analysis::Bayesian {
                generations
            } else {
                0
            },
            ..Weight::of_sequences(rows, rows.saturating_mul(columns))
        }
    }

    /// A pipeline over `count` sequences of `letters` letters in all: it bootstraps only when it builds a tree.
    pub fn of_pipeline(count: usize, letters: usize, choice: &Choice) -> Weight {
        Weight {
            replicates: if choice.tree { choice.bootstrap } else { 0 },
            ..Weight::of_sequences(count, letters)
        }
    }

    /// The trees in `newick`, and whether they would be put through the networks on a GPU of the workstation.
    pub fn of_trees(newick: &str, gpu: bool) -> Weight {
        Weight {
            trees: newick.matches(';').count(),
            gpu,
            ..Weight::default()
        }
    }

    /// A job that would run on a GPU of the workstation, whatever its size (a structure the workstation folds itself).
    pub fn on_gpu() -> Weight {
        Weight {
            gpu: true,
            ..Weight::default()
        }
    }

    /// Whether the job is big enough for the workstation when the choice is Automatic.
    pub fn is_heavy(&self) -> bool {
        self.gpu
            || self.sequences > HEAVY_SEQUENCES
            || self.letters > HEAVY_LETTERS
            || self.replicates > HEAVY_REPLICATES
            || self.generations > HEAVY_GENERATIONS
            || self.trees > HEAVY_TREES
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_job_is_big_only_past_each_limit_and_only_for_what_its_method_uses() {
        assert!(!Weight::default().is_heavy());
        // each limit is where "more than" starts
        assert!(!Weight::of_sequences(HEAVY_SEQUENCES, 0).is_heavy());
        assert!(Weight::of_sequences(HEAVY_SEQUENCES + 1, 0).is_heavy());
        assert!(!Weight::of_sequences(1, HEAVY_LETTERS).is_heavy());
        assert!(Weight::of_sequences(1, HEAVY_LETTERS + 1).is_heavy());
        // a tree counts the letters of its alignment: rows times columns
        let tree = Weight::of_tree(50, 20_001, Analysis::Likelihood, 0, 0);
        assert_eq!((tree.sequences, tree.letters), (50, 1_000_050));
        assert!(tree.is_heavy());
        // settings that the method does not use are not counted
        let nj = |bootstrap, generations| {
            Weight::of_tree(10, 100, Analysis::NeighborJoining, bootstrap, generations)
        };
        assert!(!nj(HEAVY_REPLICATES, 1_000_000_000).is_heavy());
        assert!(nj(HEAVY_REPLICATES + 1, 0).is_heavy());
        let bayes = |bootstrap, generations| {
            Weight::of_tree(10, 100, Analysis::Bayesian, bootstrap, generations)
        };
        assert!(!bayes(1_000_000, HEAVY_GENERATIONS).is_heavy());
        assert!(bayes(0, HEAVY_GENERATIONS + 1).is_heavy());
        // a pipeline bootstraps only when it builds a tree
        let mut choice = Choice {
            align: true,
            tree: false,
            primers: false,
            use_installed: false,
            bootstrap: HEAVY_REPLICATES + 1,
        };
        assert!(!Weight::of_pipeline(6, 6000, &choice).is_heavy());
        choice.tree = true;
        assert!(Weight::of_pipeline(6, 6000, &choice).is_heavy());
        // trees are counted by their semicolons, and the GPU makes any job big
        assert_eq!(Weight::of_trees("(A,B);(A,B);", false).trees, 2);
        assert!(!Weight::of_trees(&";".repeat(HEAVY_TREES), false).is_heavy());
        assert!(Weight::of_trees(&";".repeat(HEAVY_TREES + 1), false).is_heavy());
        assert!(Weight::of_trees("(A,B);", true).is_heavy() && Weight::on_gpu().is_heavy());
    }
}
