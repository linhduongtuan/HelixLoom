//! Tools → Phylogenetic Tree: align the open sequences, then build a neighbour-joining, UPGMA or
//! maximum-likelihood tree, with bootstrap support for the distance methods.

use crate::seq::msa::center_star;
use crate::seq::tree::{DistanceModel, Node};
use crate::ui::app::App;
use crate::ui::state::DocView;
use eframe::egui::{self, pos2, vec2, Align2, Color32, FontId, Painter, Sense, Stroke, Ui};
use std::collections::HashMap;

// the estimation itself is shared with the Phylogeny page and the workstation server
use crate::seq::estimate::{estimate, Analysis, EstimateOptions};

const ROW_HEIGHT: f32 = 22.0;
const LABEL_WIDTH: f32 = 150.0;
const TREE_MARGIN: f32 = 8.0;
const MIN_TREE_WIDTH: f32 = 320.0;
const TREE_MAX_HEIGHT: f32 = 320.0;
const MAX_REPLICATES: usize = 1000;

/// Generations a Bayesian run starts with, and the range the window allows.
const DEFAULT_GENERATIONS: usize = 10_000;
const GENERATION_RANGE: std::ops::RangeInclusive<usize> = 1_000..=1_000_000;

pub struct TreeWin {
    pub open: bool,
    pub analysis: Analysis,
    pub model: DistanceModel,
    /// Bootstrap replicates for the distance methods; 0 for none.
    pub bootstrap: usize,
    /// Generations of each Bayesian run.
    pub generations: usize,
    pub tree: Option<Node>,
    pub details: Vec<String>,
    pub message: String,
}

impl Default for TreeWin {
    fn default() -> Self {
        TreeWin {
            open: false,
            analysis: Analysis::NeighborJoining,
            model: DistanceModel::Kimura2,
            bootstrap: 100,
            generations: DEFAULT_GENERATIONS,
            tree: None,
            details: Vec::new(),
            message: String::new(),
        }
    }
}

/// Draws a tree left to right, branch lengths to scale, leaves in order from the top.
struct Plot {
    painter: Painter,
    scale: f32,
    top: f32,
    stroke: Stroke,
    text: Color32,
}

impl Plot {
    /// Draw `node` and everything below it; `x0` is where its branch starts. Returns the node's y.
    fn node(&self, node: &Node, x0: f32, next_row: &mut usize) -> f32 {
        let x = x0 + node.length as f32 * self.scale;
        let y = if node.is_leaf() {
            let y = self.top + ROW_HEIGHT * (*next_row as f32 + 0.5);
            *next_row += 1;
            let font = FontId::proportional(13.0);
            let at = pos2(x + 6.0, y);
            self.painter
                .text(at, Align2::LEFT_CENTER, &node.name, font, self.text);
            y
        } else {
            let ys: Vec<f32> = node
                .children
                .iter()
                .map(|c| self.node(c, x, next_row))
                .collect();
            let (lo, hi) = ys
                .iter()
                .fold((f32::MAX, f32::MIN), |(lo, hi), &y| (lo.min(y), hi.max(y)));
            self.painter
                .line_segment([pos2(x, lo), pos2(x, hi)], self.stroke);
            let y = (lo + hi) / 2.0;
            if let Some(s) = node.support {
                let font = FontId::proportional(10.0);
                // above the branch that arrives at the node, so it never runs into a leaf name
                let at = pos2(x - 3.0, y - 2.0);
                let label = format!("{s:.0}");
                self.painter.text(
                    at,
                    Align2::RIGHT_BOTTOM,
                    label,
                    font,
                    self.text.gamma_multiply(0.7),
                );
            }
            y
        };
        self.painter
            .line_segment([pos2(x0, y), pos2(x, y)], self.stroke);
        y
    }
}

fn draw_tree(ui: &mut Ui, root: &Node) {
    let width = ui.available_width().max(MIN_TREE_WIDTH);
    let height = ROW_HEIGHT * root.leaf_count() as f32 + 2.0 * TREE_MARGIN;
    let (rect, _) = ui.allocate_exact_size(vec2(width, height), Sense::hover());
    let depth = root.depth().max(f64::EPSILON) as f32;
    let color = ui.visuals().text_color();
    let plot = Plot {
        painter: ui.painter_at(rect),
        scale: (width - LABEL_WIDTH - 2.0 * TREE_MARGIN) / depth,
        top: rect.top() + TREE_MARGIN,
        stroke: Stroke::new(1.5, color),
        text: color,
    };
    plot.node(root, rect.left() + TREE_MARGIN, &mut 0);
}

/// Name and bases of every open sequence; a repeated name gets " #2", " #3" so leaves stay distinct.
pub fn open_sequences(docs: &[DocView]) -> Vec<(String, Vec<u8>)> {
    let mut seen: HashMap<&str, usize> = HashMap::new();
    docs.iter()
        .map(|d| {
            let count = seen.entry(d.doc.name.as_str()).or_insert(0);
            *count += 1;
            let name = match *count {
                1 => d.doc.name.clone(),
                n => format!("{} #{n}", d.doc.name),
            };
            (name, d.doc.seq.clone())
        })
        .collect()
}

impl App {
    /// Align the open sequences and build the chosen kind of tree.
    // ponytail: runs on the UI thread (alignments and the likelihood search use all cores, but the
    // window waits); the Phylogeny page runs the same estimation on a worker thread.
    pub fn build_tree(&mut self) {
        let seqs = open_sequences(&self.docs);
        let w = &mut self.tree_win;
        (w.tree, w.details, w.message) = (None, Vec::new(), String::new());
        if seqs.len() < 2 {
            w.message = "Open at least two sequences to build a tree.".into();
            return;
        }
        let msa = match center_star(&seqs) {
            Ok(m) => m,
            Err(e) => {
                w.message = e;
                return;
            }
        };
        w.details.push(format!(
            "{} sequences aligned into {} columns (centre-star alignment).",
            msa.rows.len(),
            msa.columns()
        ));
        let options = EstimateOptions {
            analysis: w.analysis,
            model: w.model,
            bootstrap: w.bootstrap,
            generations: w.generations,
        };
        match estimate(&msa, &options) {
            Ok(done) => {
                w.details.extend(done.details);
                w.tree = Some(done.tree);
            }
            Err(e) => w.message = e,
        }
    }

    pub fn tree_window(&mut self, ctx: &egui::Context) {
        if !self.tree_win.open {
            return;
        }
        let mut open = true;
        let mut build = false;
        crate::ui::design::window("Phylogenetic Tree")
            .open(&mut open)
            .default_size([580.0, 480.0])
            .show(ctx, |ui| {
                ui.label(format!(
                    "Uses the {} open sequences: they are aligned together first.",
                    self.docs.len()
                ));
                let w = &mut self.tree_win;
                ui.horizontal(|ui| {
                    ui.radio_value(&mut w.analysis, Analysis::NeighborJoining, "Neighbor-Joining");
                    ui.radio_value(&mut w.analysis, Analysis::Upgma, "UPGMA");
                    ui.radio_value(&mut w.analysis, Analysis::Likelihood, "Maximum likelihood");
                    ui.radio_value(&mut w.analysis, Analysis::Bayesian, "Bayesian");
                });
                if w.analysis == Analysis::Likelihood {
                    ui.label("Models JC69, K80, F81 and HKY85 are compared by BIC; the topology is improved by swaps.");
                } else if w.analysis == Analysis::Bayesian {
                    ui.horizontal(|ui| {
                        ui.label("HKY85, two independent runs of");
                        ui.add(
                            egui::DragValue::new(&mut w.generations)
                                .range(GENERATION_RANGE)
                                .speed(100.0),
                        );
                        ui.label("generations");
                    });
                } else {
                    ui.horizontal(|ui| {
                        for m in DistanceModel::ALL {
                            ui.radio_value(&mut w.model, m, m.label());
                        }
                        ui.label("Bootstrap");
                        ui.add(egui::DragValue::new(&mut w.bootstrap).range(0..=MAX_REPLICATES));
                    });
                }
                build = ui.button("Build tree").clicked();
                if !w.message.is_empty() {
                    ui.label(&w.message);
                }
                w.details.iter().for_each(|line| {
                    ui.label(line);
                });
                if let Some(tree) = &w.tree {
                    egui::ScrollArea::vertical()
                        .max_height(TREE_MAX_HEIGHT)
                        .show(ui, |ui| draw_tree(ui, tree));
                    ui.separator();
                    let mut text = tree.newick();
                    ui.add(
                        egui::TextEdit::multiline(&mut text)
                            .font(egui::TextStyle::Monospace)
                            .desired_rows(2)
                            .desired_width(f32::INFINITY),
                    );
                    if ui.button("Copy Newick").clicked() {
                        ui.ctx().copy_text(tree.newick());
                    }
                }
            });
        if build {
            self.build_tree();
        }
        self.tree_win.open = open;
    }
}
