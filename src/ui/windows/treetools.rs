//! Tools → Tree Tools: read, root, sort and compare Newick trees, the jobs of ETE and DendroPy.

use crate::seq::newick::{number, NwTree};
use crate::seq::treetools::{distance_table, ladderize, midpoint_root, robinson_foulds, summary};
use crate::ui::app::App;
use crate::ui::dialog::{self, Wanted};
use crate::ui::platform::route::Source;
use eframe::egui::{self, RichText, ScrollArea};

const TEXT_BOX_MAX_HEIGHT: f32 = 140.0;
const LINES_MAX_HEIGHT: f32 = 150.0;
/// Most trees compared pair by pair; with more, only the first against each of the others.
const MAX_ALL_PAIRS: usize = 6;
/// Two small random trees on the same ten tips, so the first run shows something.
const DEMO_TREES: &str = "((((B:1.277,A:1.568):1.702,J:0.051):0.831,(I:0.459,H:1.825):1.701):0.187,((G:0.966,F:1.962):0.804,(E:0.825,(D:0.576,C:0.22):0.192):1.918):1.544);\n\
(((F:1.13,(E:0.263,D:1.049):1.39):1.24,(C:0.807,(B:0.548,A:1.493):0.991):1.465):0.656,(J:0.359,(I:0.203,(H:0.814,G:0.871):0.142):1.277):0.495);";

#[derive(Default)]
pub struct TreeToolsWin {
    pub open: bool,
    /// One or more Newick trees.
    pub newick: String,
    pub message: String,
    /// Details of the last operation, one per line.
    pub lines: Vec<String>,
}

/// Ask for a Newick file and read it: `None` when the dialog is cancelled, else its name and text,
/// or what went wrong. (In a browser the file arrives later, for `wanted`: see `dialog`.)
pub fn read_newick_file(wanted: Wanted) -> Option<Result<(String, String), String>> {
    let ask = dialog::FileDialog::new()
        .set_title("Open a Newick tree")
        .add_filter("Newick", &["nwk", "newick", "tree", "trees", "tre", "txt"]);
    let path = dialog::pick_one(wanted, ask)?;
    Some(newick_text(Source::Path(&path)))
}

/// A Newick file's name and text, or what went wrong.
pub fn newick_text(file: Source) -> Result<(String, String), String> {
    let name = file
        .path()
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("tree")
        .to_string();
    file.text()
        .map(|text| (name, text))
        .map_err(|e| format!("Couldn't read {}: {e}", file.shown()))
}

impl TreeToolsWin {
    fn say(&mut self, message: impl Into<String>, lines: Vec<String>) {
        self.message = message.into();
        self.lines = lines;
    }

    fn trees(&self) -> Result<Vec<NwTree>, String> {
        NwTree::parse_all(&self.newick)
    }

    /// Replace every tree in the box by what `change` makes of it.
    fn rewrite(&mut self, change: impl Fn(&NwTree) -> Result<NwTree, String>, done: String) {
        let outcome = self.trees().and_then(|trees| {
            trees
                .iter()
                .map(|t| change(t).map(|t| t.newick()))
                .collect::<Result<Vec<_>, _>>()
        });
        match outcome {
            Ok(texts) => {
                self.newick = texts.join("\n");
                self.say(done, Vec::new());
            }
            Err(e) => self.say(e, Vec::new()),
        }
    }

    pub fn midpoint(&mut self) {
        let first = self
            .trees()
            .and_then(|t| midpoint_root(&t[0]).map(|(_, found)| found));
        match first {
            Ok(found) => {
                let done = format!(
                    "Rooted halfway along the longest path, {} to {} ({}).",
                    found.tip_a,
                    found.tip_b,
                    number(found.diameter)
                );
                self.rewrite(|t| midpoint_root(t).map(|(rooted, _)| rooted), done);
            }
            Err(e) => self.say(e, Vec::new()),
        }
    }

    pub fn ladderize(&mut self) {
        self.rewrite(
            |t| Ok(ladderize(t)),
            "Smaller clades now come first at every node.".into(),
        );
    }

    pub fn describe(&mut self) {
        match self.trees() {
            Ok(trees) => {
                let lines = trees
                    .iter()
                    .enumerate()
                    .map(|(i, t)| {
                        let s = summary(t);
                        format!(
                            "Tree {}: {} tips, {} inner nodes, total length {}, height {}, {}",
                            i + 1,
                            s.tips,
                            s.inner,
                            number(s.length),
                            number(s.height),
                            if s.binary {
                                "binary"
                            } else {
                                "has multifurcations"
                            }
                        )
                    })
                    .collect();
                self.say(format!("{} tree(s) in the box.", trees.len()), lines);
            }
            Err(e) => self.say(e, Vec::new()),
        }
    }

    /// The first tree's distance table, on the clipboard.
    pub fn copy_distances(&mut self, ctx: &egui::Context) {
        match self.trees() {
            Ok(trees) => {
                let tips = summary(&trees[0]).tips;
                ctx.copy_text(distance_table(&trees[0]));
                self.say(
                    format!(
                        "Copied the {tips} × {tips} table of distances along the branches of the first tree (tab-separated)."
                    ),
                    Vec::new(),
                );
            }
            Err(e) => self.say(e, Vec::new()),
        }
    }

    pub fn compare(&mut self) {
        let trees = match self.trees() {
            Ok(t) if t.len() >= 2 => t,
            Ok(_) => {
                return self.say(
                    "Put two or more trees in the box, one after the other.",
                    Vec::new(),
                )
            }
            Err(e) => return self.say(e, Vec::new()),
        };
        let pairs: Vec<(usize, usize)> = if trees.len() <= MAX_ALL_PAIRS {
            (0..trees.len())
                .flat_map(|i| (i + 1..trees.len()).map(move |j| (i, j)))
                .collect()
        } else {
            (1..trees.len()).map(|j| (0, j)).collect()
        };
        let mut lines = Vec::new();
        for (i, j) in pairs {
            match robinson_foulds(&trees[i], &trees[j]) {
                Ok(rf) => lines.push(format!(
                    "Tree {} and tree {}: Robinson–Foulds distance {} of at most {}{}",
                    i + 1,
                    j + 1,
                    rf.distance,
                    rf.max,
                    if rf.max == 0 {
                        String::new()
                    } else {
                        format!(" ({:.2})", rf.distance as f64 / rf.max as f64)
                    }
                )),
                Err(e) => {
                    return self.say(
                        format!("Tree {} and tree {}: {e}", i + 1, j + 1),
                        Vec::new(),
                    )
                }
            }
        }
        self.say(
            "Splits found in only one tree, counted without regard to where each tree is rooted.",
            lines,
        );
    }
}

impl App {
    /// The two demo trees, rooted and compared.
    pub fn load_treetools_demo(&mut self) {
        let w = &mut self.treetools_win;
        w.newick = DEMO_TREES.to_string();
        w.midpoint();
        w.compare();
    }

    /// Put a Newick file that was read into the window's box, or say why it could not be read.
    fn use_newick_file(&mut self, read: Result<(String, String), String>) {
        let w = &mut self.treetools_win;
        match read {
            Ok((name, text)) => {
                w.newick = text;
                w.say(format!("Loaded {name}."), Vec::new());
            }
            Err(e) => w.say(e, Vec::new()),
        }
    }

    /// A Newick file for the window, from a source the browser brought.
    pub(in crate::ui) fn load_treetools_tree(&mut self, file: Source) {
        self.use_newick_file(newick_text(file));
    }

    pub fn treetools_window(&mut self, ctx: &egui::Context) {
        if !self.treetools_win.open {
            return;
        }
        let from_tree_window = self.tree_win.tree.as_ref().map(|t| t.newick());
        let (mut open, mut load, mut use_built) = (true, false, false);
        crate::ui::design::window("Tree Tools")
            .open(&mut open)
            .default_width(600.0)
            .show(ctx, |ui| {
                let w = &mut self.treetools_win;
                ui.label(
                    "Newick trees, one after the other. Root them halfway along the longest path, sort \
                     them, measure them, or compare two of them.",
                );
                ui.horizontal(|ui| {
                    load = ui.button("Open tree file…").clicked();
                    use_built = ui
                        .add_enabled(
                            from_tree_window.is_some(),
                            egui::Button::new("Tree from the Phylogenetic Tree window"),
                        )
                        .clicked();
                    if ui.button("Clear").clicked() {
                        w.newick.clear();
                        w.say(String::new(), Vec::new());
                    }
                });
                ScrollArea::vertical()
                    .max_height(TEXT_BOX_MAX_HEIGHT)
                    .show(ui, |ui| {
                        ui.add(
                            egui::TextEdit::multiline(&mut w.newick)
                                .font(egui::TextStyle::Monospace)
                                .hint_text("Paste Newick trees such as ((A:1,B:1):2,C:3);")
                                .desired_rows(5)
                                .desired_width(f32::INFINITY),
                        );
                    });
                ui.horizontal_wrapped(|ui| {
                    if ui.button("Midpoint root").clicked() {
                        w.midpoint();
                    }
                    if ui.button("Ladderize").clicked() {
                        w.ladderize();
                    }
                    if ui.button("Describe").clicked() {
                        w.describe();
                    }
                    if ui.button("Compare trees").clicked() {
                        w.compare();
                    }
                    if ui.button("Copy distance table").clicked() {
                        w.copy_distances(ui.ctx());
                    }
                    if ui.button("Copy Newick").clicked() {
                        ui.ctx().copy_text(w.newick.clone());
                    }
                });
                if !w.message.is_empty() {
                    ui.label(&w.message);
                }
                ScrollArea::vertical()
                    .id_salt("treetools_lines")
                    .max_height(LINES_MAX_HEIGHT)
                    .show(ui, |ui| {
                        for line in &w.lines {
                            ui.label(RichText::new(line).monospace());
                        }
                    });
            });
        if load {
            if let Some(read) = read_newick_file(Wanted::TreeToolsTree) {
                self.use_newick_file(read);
            }
        }
        if let (true, Some(text)) = (use_built, from_tree_window) {
            let w = &mut self.treetools_win;
            w.newick = text;
            w.say(
                "Took the tree from the Phylogenetic Tree window.",
                Vec::new(),
            );
        }
        self.treetools_win.open = open;
    }
}
