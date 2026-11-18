//! What the pages show, kept with the session and put back with it: the Alignment page (its list of sequences and its alignment),
//! the tree that the Trees and Phylogeny pages show, the structure, and the runs of the Pipelines page, each with the record of
//! what made it. A result is kept when there is room for it; the sequences never wait for one (`session::Store::keep_session`).
//!
//! Not kept: what the tool windows show, and how a page is scrolled, turned or zoomed.

use crate::kept::{fields, line, value, word};
use crate::provenance::Provenance;
use crate::seq::estimate::{Analysis, EstimateOptions};
use crate::seq::newick::NwTree;
use crate::seq::tree::DistanceModel;
use crate::ui::app::App;
use crate::ui::platform::phylogeny::{GENERATION_RANGE, MAX_REPLICATES};

pub const ALIGNMENT: &str = "alignment";
pub const TREE: &str = "tree";
pub const STRUCTURE: &str = "structure";
pub const PIPELINES: &str = "pipelines";

/// Said on the Trees page when the tree of the last session is put back.
pub const TREE_IS_BACK: &str = "The tree of your last session.";

const MODELS: [(DistanceModel, &str); 3] = [
    (DistanceModel::P, "p"),
    (DistanceModel::JukesCantor, "jc"),
    (DistanceModel::Kimura2, "k2p"),
];

/// What a kind of result is called in a sentence; none for a kind that this version does not know.
pub fn called(kind: &str) -> Option<&'static str> {
    match kind {
        ALIGNMENT => Some("the Alignment page"),
        TREE => Some("the tree"),
        STRUCTURE => Some("the structure"),
        PIPELINES => Some("the Pipelines page"),
        _ => None,
    }
}

/// The kinds of `kinds` that this version knows, as a sentence names them: `the Alignment page, the tree and the structure`.
pub fn named<'a>(kinds: impl IntoIterator<Item = &'a str>) -> String {
    let names: Vec<String> = kinds
        .into_iter()
        .filter_map(called)
        .map(str::to_string)
        .collect();
    crate::report::listed(&names)
}

impl App {
    /// The results that there are to keep, each with a number that is another one whenever the result has changed. Asked at every
    /// frame: it reads no sequence, tree or structure.
    pub(super) fn result_stamps(&self) -> Vec<(&'static str, u64)> {
        let (trees, made) = (&self.tree_win, &self.platform.phylogeny.made);
        let tree = trees
            .tree
            .as_ref()
            .map(|_| trees.changed * 2 + u64::from(made.is_some()));
        [
            (ALIGNMENT, self.platform.alignment.stamp()),
            (TREE, tree),
            (STRUCTURE, self.platform.structure.stamp()),
            (PIPELINES, self.platform.pipelines.stamp()),
        ]
        .into_iter()
        .filter_map(|(kind, stamp)| Some((kind, stamp?)))
        .collect()
    }

    /// One of them as the text that `put_result_back` reads; none when its page shows nothing.
    pub(super) fn result_text(&self, kind: &str) -> Option<String> {
        let alignment = &self.platform.alignment;
        match kind {
            ALIGNMENT => alignment.stamp().map(|_| alignment.kept()),
            TREE => self.tree_kept(),
            STRUCTURE => self.platform.structure.kept(),
            PIPELINES => self.platform.pipelines.kept_text(),
            _ => None,
        }
    }

    /// Show a kept result on its page again. One that cannot be read changes nothing: why is said.
    pub(super) fn put_result_back(&mut self, kind: &str, text: &str) -> Result<(), String> {
        match kind {
            ALIGNMENT => self.platform.alignment.put_back(text),
            TREE => self.put_tree_back(text),
            STRUCTURE => self.platform.structure.put_back(text),
            PIPELINES => self.platform.pipelines.put_back(text),
            _ => Err(format!(
                "This version does not know what a kept {kind:?} is."
            )),
        }
    }

    /// The tree with the lines that say how it was made and, for one that the Phylogeny page estimated, its record and settings.
    fn tree_kept(&self) -> Option<String> {
        let tree = self.tree_win.tree.as_ref()?;
        let mut text = line(&["newick", &tree.newick_exact()]);
        for detail in &self.tree_win.details {
            text += &line(&["detail", detail]);
        }
        if let Some((made, o)) = &self.platform.phylogeny.made {
            let (bootstrap, generations) = (o.bootstrap.to_string(), o.generations.to_string());
            text += &line(&[
                "made",
                o.analysis.name(),
                word(&MODELS, o.model),
                &bootstrap,
                &generations,
                &made.kept(),
            ]);
        }
        Some(text)
    }

    fn put_tree_back(&mut self, text: &str) -> Result<(), String> {
        let (mut tree, mut details, mut made) = (None, Vec::new(), None);
        for kept in text.lines() {
            let kept = fields(kept);
            let kept: Vec<&str> = kept.iter().map(String::as_str).collect();
            match kept.as_slice() {
                ["newick", held] => tree = NwTree::parse_all(held)?.first().map(NwTree::to_node),
                ["detail", held] => details.push(held.to_string()),
                ["made", analysis, model, bootstrap, generations, record] => {
                    made = settings_of(analysis, model, bootstrap, generations)
                        .and_then(|options| Some((Provenance::from_kept(record)?, options)));
                }
                _ => {}
            }
        }
        let tree = tree.ok_or("Nothing of the tree that was kept can be read.")?;
        self.show_tree(tree, TREE_IS_BACK);
        self.tree_win.details = details;
        self.platform.phylogeny.made = made;
        Ok(())
    }
}

/// The settings of an estimation as they were kept; none when a word is not one, or a number is one that the page does not allow
/// (the settings are what "Run again" runs with).
fn settings_of(
    analysis: &str,
    model: &str,
    bootstrap: &str,
    generations: &str,
) -> Option<EstimateOptions> {
    let bootstrap: usize = bootstrap.parse().ok().filter(|n| *n <= MAX_REPLICATES)?;
    let generations: usize = generations.parse().ok()?;
    Some(EstimateOptions {
        analysis: Analysis::from_name(analysis)?,
        model: value(&MODELS, model)?,
        bootstrap,
        generations: GENERATION_RANGE
            .contains(&generations)
            .then_some(generations)?,
    })
}

#[cfg(test)]
mod tests;
