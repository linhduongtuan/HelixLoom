//! What the buttons of the Trees page do to the tree: the edits (root here, rotate, order, collapse weak branches, rename tips, cut out a
//! clade), the support values taken from other trees, and the marks on clades (collapsed, highlighted).

use super::looks::{laid_index, node_id, Annotated, Looks, SupportFrom};
use super::{Clicks, TreesPage};
use crate::seq::newick::NwTree;
use crate::seq::treeedit;
use crate::seq::treelayout::{lay_out_with, Laid};
use crate::seq::treestyle::PALETTE;
use crate::seq::treetools;
use crate::ui::app::App;

/// The share of the trees as a number from 0 to 100 with no more digits than a person can use.
fn percent(share: f64) -> f64 {
    (share * 10_000.0).round() / 100.0
}

impl App {
    /// Make sure that the page's annotated tree is the app's tree: when another tree has been put there (an estimate, a file) the page
    /// reads its text again, without annotations, and forgets what it knew of the other one.
    pub(in crate::ui) fn sync_annotated(&mut self) {
        let changed = self.tree_win.changed;
        let page = &mut self.platform.trees;
        if page.doc.as_ref().is_some_and(|d| d.for_change == changed) {
            return;
        }
        page.doc = self
            .tree_win
            .tree
            .as_ref()
            .and_then(|t| NwTree::parse_all(&t.newick_exact()).ok())
            .and_then(|mut v| (!v.is_empty()).then(|| v.remove(0)))
            .map(|nw| Annotated {
                nw,
                for_change: changed,
            });
        page.looks = Looks::default();
        page.selected = None;
    }

    /// Put `nw` on the page as the tree, as an edit of the one that was there (what was marked on that one is gone with it).
    pub(super) fn commit_tree(&mut self, nw: NwTree, message: &str) {
        // it is no longer the tree that the estimation made: running that again would not give it
        self.platform.phylogeny.made = None;
        self.tree_win.set_tree(Some(nw.to_node()));
        let page = &mut self.platform.trees;
        page.looks.forget_marks();
        page.selected = None;
        page.message = message.to_string();
        page.doc = Some(Annotated {
            nw,
            for_change: self.tree_win.changed,
        });
    }

    /// Show a tree read from a file or pasted, with its annotations; the first of `trees` is drawn.
    pub(in crate::ui) fn show_annotated(&mut self, mut trees: Vec<NwTree>, message: &str) {
        if trees.is_empty() {
            return;
        }
        let nw = trees.remove(0);
        self.show_tree(nw.to_node(), message);
        let numbered = nw
            .nodes
            .iter()
            .any(|n| !n.children.is_empty() && n.name.parse::<f64>().is_ok());
        let page = &mut self.platform.trees;
        page.looks = Looks::default();
        // a tree from BEAST says its support as an annotation: show it unless the nodes are labelled already
        if nw.attr_keys().iter().any(|k| k == "posterior") && !numbered {
            page.looks.style.labels = vec!["posterior".to_string()];
        }
        page.doc = Some(Annotated {
            nw,
            for_change: self.tree_win.changed,
        });
    }

    /// The value of every branch of the tree on the page that other trees give it, as annotations (taken from `trees`).
    pub(in crate::ui) fn add_support(&mut self, trees: &[NwTree], shown: &str) {
        let page = &mut self.platform.trees;
        page.looks.awaiting_support = false;
        let Some(nw) = page.doc.as_ref().map(|d| d.nw.clone()) else {
            page.message = "There is no tree on the page to add support to.".into();
            return;
        };
        let label = {
            let t = page.looks.support_label.trim();
            if t.is_empty() { "support" } else { t }.to_string()
        };
        let made = match (page.looks.support_from, trees.first()) {
            (SupportFrom::Written, Some(one)) => {
                treeedit::values_from(&nw, one, page.looks.support_key.trim()).map(|v| {
                    (
                        treeedit::with_values(&nw, &label, &v),
                        v.iter().flatten().count(),
                    )
                })
            }
            _ => treeedit::frequencies(&nw, trees).map(|s| {
                let shares: Vec<Option<f64>> =
                    s.iter().map(|x| x.map(|x| percent(x.frequency))).collect();
                let rivals: Vec<Option<f64>> =
                    s.iter().map(|x| x.map(|x| percent(x.conflict))).collect();
                let with = treeedit::with_values(&nw, &label, &shares);
                let count = shares.iter().flatten().count();
                (
                    treeedit::with_values(&with, &format!("{label}_conflict"), &rivals),
                    count,
                )
            }),
        };
        match made {
            Ok((with, count)) => {
                if let Some(d) = &mut page.doc {
                    d.nw = with;
                }
                if !page.looks.style.labels.contains(&label) {
                    page.looks.style.labels.push(label.clone());
                }
                page.message = format!(
                    "Added \"{label}\" from {shown} ({} tree{}): {count} branches have a value.",
                    trees.len(),
                    if trees.len() == 1 { "" } else { "s" }
                );
            }
            Err(why) => page.message = why,
        }
    }

    /// The buttons that change the tree or mark its clades.
    pub(super) fn apply_tree_edits(&mut self, c: &Clicks) {
        self.sync_annotated();
        let Some(tree) = self.tree_win.tree.clone() else {
            return;
        };
        let page = &self.platform.trees;
        let Some(nw) = page.doc.as_ref().map(|d| d.nw.clone()) else {
            return;
        };
        let laid = lay_out_with(&tree, page.options.layout, &page.looks.collapsed);
        let id = page.selected.and_then(|i| node_id(&nw, &laid, i));
        let wants_node = c.reroot || c.rotate || c.extract;
        if wants_node && id.is_none() {
            self.platform.trees.message = "Pick a node first: click it on the tree.".into();
            return;
        }
        self.edit_with(&nw, &laid, id, c);
        self.mark_clades(&laid, c);
    }

    fn edit_with(&mut self, nw: &NwTree, laid: &Laid, id: Option<usize>, c: &Clicks) {
        let looks = &self.platform.trees.looks;
        let done: Option<(Result<NwTree, String>, &str)> = if let (true, Some(id)) = (c.reroot, id)
        {
            Some((
                treeedit::reroot_above(nw, id, 0.5),
                "Rooted on the branch of the picked node.",
            ))
        } else if let (true, Some(id)) = (c.rotate, id) {
            Some((
                Ok(treeedit::rotate(nw, id)),
                "Turned the branches of the picked node round.",
            ))
        } else if let (true, Some(id)) = (c.extract, id) {
            Some((
                treeedit::extract_clade(nw, id),
                "Showing the picked clade as a tree of its own.",
            ))
        } else if c.order_up || c.order_down || c.ladderize {
            let up = c.order_up || c.ladderize;
            let text = if up {
                "Smaller clades first."
            } else {
                "Larger clades first."
            };
            Some((Ok(treeedit::order(nw, up)), text))
        } else if c.midpoint {
            Some((
                treetools::midpoint_root(nw).map(|(t, _)| t),
                "Rooted at the midpoint (annotations other than node labels are not kept).",
            ))
        } else if c.weak {
            let (tree, gone) =
                treeedit::collapse_below(nw, looks.weak_key.trim(), looks.weak_below);
            self.platform.trees.message = format!(
                "Took out {gone} branch{} below {}.",
                if gone == 1 { "" } else { "es" },
                looks.weak_below
            );
            Some((Ok(tree), ""))
        } else if c.rename {
            let table = treeedit::parse_name_table(&looks.names_paste);
            let (tree, done) = treeedit::rename_tips(nw, &table);
            self.platform.trees.message = format!(
                "Renamed {done} tip{} ({} in the table).",
                if done == 1 { "" } else { "s" },
                table.len()
            );
            Some((Ok(tree), ""))
        } else {
            None
        };
        let _ = laid;
        if let Some((result, text)) = done {
            match result {
                Ok(tree) => {
                    let said = self.platform.trees.message.clone();
                    self.commit_tree(tree, if text.is_empty() { &said } else { text });
                }
                Err(why) => self.platform.trees.message = why,
            }
        }
    }

    /// Collapse, expand and highlight the picked clade.
    fn mark_clades(&mut self, laid: &Laid, c: &Clicks) {
        let page = &mut self.platform.trees;
        let Some(origin) = page
            .selected
            .and_then(|i| laid.nodes.get(i))
            .map(|n| n.origin)
        else {
            return;
        };
        if c.fold {
            if let Some(at) = page.looks.collapsed.iter().position(|&o| o == origin) {
                page.looks.collapsed.remove(at);
            } else {
                page.looks.collapsed.push(origin);
            }
        }
        if let Some(k) = c.colour {
            page.looks.highlights.retain(|(o, _)| *o != origin);
            page.looks
                .highlights
                .push((origin, PALETTE[k % PALETTE.len()]));
        }
        if c.unmark {
            page.looks.highlights.retain(|(o, _)| *o != origin);
        }
        // a clade that was collapsed and is now inside another collapsed clade is not drawn: keep the picked one in range
        let _ = laid_index(laid, origin);
    }
}

impl TreesPage {
    /// The page's tree with its annotations, when it is the app's tree.
    pub(in crate::ui) fn annotated(&self) -> Option<&NwTree> {
        self.doc.as_ref().map(|d| &d.nw)
    }
}
