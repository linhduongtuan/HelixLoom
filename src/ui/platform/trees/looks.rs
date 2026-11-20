//! What the Trees page keeps beside the tree: the tree with its annotations, how it is dressed (colours, labels, an axis), and the clades
//! that were collapsed or highlighted. All of it is about the tree that is on the page, and goes when another is put there.

use crate::seq::newick::NwTree;
use crate::seq::treelayout::Laid;
use crate::seq::treestyle::{Axis, Rgb, StyleOptions};

/// The tree on the page as the tree tools see it, and which tree (`TreeWin::changed`) of the app it is.
pub struct Annotated {
    pub nw: NwTree,
    pub for_change: u64,
}

/// How a values from other trees are taken.
#[derive(Clone, Copy, PartialEq, Eq, Default, Debug)]
pub enum SupportFrom {
    /// The share of the trees that have the branch (bootstrap replicates, posterior samples).
    #[default]
    Share,
    /// The number written on the same branch of one other tree.
    Written,
}

#[derive(Default)]
pub struct Looks {
    pub style: StyleOptions,
    pub tip_colours: bool,
    pub axis_on: bool,
    pub axis: Axis,
    /// Clades drawn as triangles, and clades with a colour behind them, by their place in the tree (`Placed::origin`).
    pub collapsed: Vec<usize>,
    pub highlights: Vec<(usize, Rgb)>,
    pub find: String,
    /// Collapsing weak branches: the annotation (empty: the label) and the value below which a branch goes.
    pub weak_key: String,
    pub weak_below: f64,
    pub names_paste: String,
    pub support_label: String,
    pub support_from: SupportFrom,
    pub support_key: String,
    /// The next tree file that arrives is the source of support values for the tree on the page, not a tree to show.
    pub awaiting_support: bool,
}

impl Looks {
    /// The clades marked on a tree that is gone.
    pub fn forget_marks(&mut self) {
        self.collapsed.clear();
        self.highlights.clear();
    }
}

/// The place in the layout of the node that was number `origin` in the tree.
pub fn laid_index(laid: &Laid, origin: usize) -> Option<usize> {
    laid.nodes.iter().position(|n| n.origin == origin)
}

/// The node of the annotated tree that node `i` of the layout is.
pub fn node_id(nw: &NwTree, laid: &Laid, i: usize) -> Option<usize> {
    nw.preorder().get(laid.nodes.get(i)?.origin).copied()
}

/// The names of the tips below node `id` of the annotated tree, in order.
pub fn tip_names_under(nw: &NwTree, id: usize) -> Vec<String> {
    let mut out = Vec::new();
    let mut stack = vec![id];
    while let Some(k) = stack.pop() {
        if nw.is_leaf(k) {
            out.push(nw.nodes[k].name.clone());
        }
        stack.extend(nw.nodes[k].children.iter().rev());
    }
    out
}

/// The tips of the layout (not collapsed clades) whose names hold `text`, in any case; none for an empty text.
pub fn found_tips(laid: &Laid, text: &str) -> Vec<usize> {
    let text = text.trim().to_lowercase();
    if text.is_empty() {
        return Vec::new();
    }
    laid.nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.is_tip() && n.cartoon.is_none() && n.name.to_lowercase().contains(&text))
        .map(|(i, _)| i)
        .collect()
}

/// For each node of the layout, what its annotations say, for the text that shows when the pointer rests on it.
pub fn notes(nw: &NwTree, laid: &Laid) -> Vec<String> {
    let order = nw.preorder();
    laid.nodes
        .iter()
        .map(|n| {
            let Some(&id) = order.get(n.origin) else {
                return String::new();
            };
            nw.nodes[id]
                .attrs
                .iter()
                .take(12)
                .map(|(k, v)| format!("{k} = {v}"))
                .collect::<Vec<_>>()
                .join("\n")
        })
        .collect()
}
