//! The record of a run of the fixed route (`crate::provenance`): what it was made from and with, and what made it.

use super::{target_text, Source};
use crate::pipeline::steps::Named;
use crate::pipeline::{Choice, Outcome};
use crate::provenance;
use crate::remote::tasks::settings_text;

/// What a run was ticked with: where the sequences come from, and the boxes for aligning, the tree, the primers, the
/// installed programs and the bootstrap.
pub type Kept = (Source, [bool; 5]);

/// The inputs and the settings of a run, as its record says them: every sequence by its name, the steps, and for the
/// primers where they go and the primer settings (Tools → Primer3), which are an input like the sequences.
pub fn described(
    seqs: &[Named],
    target: Option<(usize, usize)>,
    settings: &crate::seq::primer3::Settings,
    choice: &Choice,
) -> (Vec<provenance::Input>, Vec<(String, String)>) {
    let mut inputs: Vec<provenance::Input> = seqs
        .iter()
        .map(|(name, bases)| provenance::Input::of(name, bases))
        .collect();
    let steps: Vec<&str> = [
        (choice.align, "align"),
        (choice.tree, "tree"),
        (choice.primers, "primers"),
    ]
    .into_iter()
    .filter_map(|(on, name)| on.then_some(name))
    .collect();
    let mut parameters = vec![
        ("steps".to_string(), steps.join(", ")),
        (
            "installed programs".to_string(),
            if choice.use_installed {
                "used where found"
            } else {
                "not used"
            }
            .to_string(),
        ),
    ];
    if choice.tree {
        parameters.push((
            "bootstrap replicates".to_string(),
            choice.bootstrap.to_string(),
        ));
    }
    if choice.primers {
        inputs.push(provenance::Input::of(
            "primer settings",
            settings_text(settings).as_bytes(),
        ));
        parameters.push(("primers go".to_string(), target_text(target)));
    }
    (inputs, parameters)
}

/// The programs that made the parts of an outcome (an aligner, a tree builder and a primer picker: never two the same), and
/// the time they took together.
pub(super) fn made_by(outcome: &Outcome) -> (Vec<String>, f64) {
    let parts = [
        outcome.aligned.as_ref().map(|a| (&a.engine, a.seconds)),
        outcome.tree.as_ref().map(|t| (&t.engine, t.seconds)),
        outcome.primers.as_ref().map(|p| (&p.engine, p.seconds)),
    ];
    let made: Vec<(&String, f64)> = parts.into_iter().flatten().collect();
    (
        made.iter().map(|(engine, _)| (*engine).clone()).collect(),
        made.iter().map(|(_, seconds)| seconds).sum(),
    )
}
