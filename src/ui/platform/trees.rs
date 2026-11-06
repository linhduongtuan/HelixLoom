//! Trees: draw a tree as a phylogram, a cladogram or a circle, with support values and a scale bar, pick a clade
//! to see its tips, and save the picture as SVG. The tree is the one the Phylogeny page estimated, or one opened
//! from a Newick file or pasted.

use super::route::Source;
use super::treeview::{self, Camera, Options};
use super::Workspace;
use crate::seq::newick::NwTree;
use crate::seq::tree::Node;
use crate::seq::treelayout::{lay_out, to_svg, Laid, Layout, SvgOptions};
use crate::seq::treetools;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::dialog::{self, Wanted};
use eframe::egui::{self, RichText, ScrollArea, Ui};

const OPTIONS_WIDTH: f32 = 270.0;
const INFO_WIDTH: f32 = 260.0;
const NAMES_HEIGHT: f32 = 260.0;

#[derive(Default)]
pub struct TreesPage {
    pub options: Options,
    camera: Camera,
    selected: Option<usize>,
    paste: String,
    pub message: String,
}

impl TreesPage {
    /// Back to the whole tree, nothing picked.
    pub fn reset_view(&mut self) {
        self.camera = Camera::default();
        self.selected = None;
    }
}

#[derive(Default)]
struct Clicks {
    open: bool,
    example: bool,
    draw_pasted: bool,
    midpoint: bool,
    ladderize: bool,
    save_svg: bool,
    save_newick: bool,
    copy: bool,
    fit: bool,
    to_phylogeny: bool,
}

/// Rebuild a tree with one of the tree tools, which work on Newick trees.
fn through_tools(
    tree: &Node,
    tool: impl Fn(&NwTree) -> Result<NwTree, String>,
) -> Result<Node, String> {
    let parsed = NwTree::parse_all(&tree.newick())?.remove(0);
    tool(&parsed).map(|t| t.to_node())
}

impl App {
    pub(super) fn trees_page(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let look = design::look(self.prefs.dark);
        let tree = self.tree_win.tree.clone();
        let mut c = Clicks::default();
        let page = &mut self.platform.trees;
        let laid = tree.as_ref().map(|t| lay_out(t, page.options.layout));
        design::page(ui, |ui| {
            design::page_header(
                ui,
                &look,
                "Trees",
                "Draw and explore phylogenetic trees",
                |ui| {
                    ui.add_enabled_ui(tree.is_some(), |ui| {
                        c.save_svg = design::primary_button(ui, &look, "Save SVG…").clicked();
                        c.save_newick = ui.button("Save Newick…").clicked();
                        c.copy = ui.button("Copy Newick").clicked();
                        c.fit = ui.button("Fit to window").clicked();
                    });
                },
            );
            egui::Panel::left("trees_options")
                .exact_size(OPTIONS_WIDTH)
                .resizable(false)
                .frame(egui::Frame::new())
                .show(ui, |ui| {
                    ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            design::card(ui, &look, |ui| {
                                options_panel(ui, &look, page, tree.is_some(), &mut c)
                            });
                        });
                });
            if let Some(laid) = &laid {
                egui::Panel::right("trees_info")
                    .exact_size(INFO_WIDTH)
                    .resizable(false)
                    .frame(egui::Frame::new().inner_margin(egui::Margin {
                        left: design::PAGE_MARGIN,
                        ..egui::Margin::ZERO
                    }))
                    .show(ui, |ui| {
                        design::card(ui, &look, |ui| info_panel(ui, &look, page, laid));
                    });
            }
            egui::CentralPanel::default()
                .frame(egui::Frame::new().inner_margin(egui::Margin {
                    left: design::PAGE_MARGIN,
                    ..egui::Margin::ZERO
                }))
                .show(ui, |ui| canvas(ui, &look, page, laid.as_ref()));
        });
        self.apply_tree_clicks(c, ui.ctx(), out);
    }

    fn apply_tree_clicks(&mut self, c: Clicks, ctx: &egui::Context, out: &mut Vec<Action>) {
        if c.open {
            self.open_tree_file();
        }
        if c.example {
            self.show_tree(
                crate::sample::demo_tree(),
                "The example tree: 24 made-up strains in three groups.",
            );
        }
        if c.draw_pasted {
            let text = std::mem::take(&mut self.platform.trees.paste);
            match NwTree::parse_all(&text) {
                Ok(mut trees) => {
                    let n = trees.len();
                    self.show_tree(
                        trees.remove(0).to_node(),
                        &format!("Drew the first of {n} pasted tree(s)."),
                    );
                }
                Err(e) => self.platform.trees.message = e,
            }
        }
        let tree = self.tree_win.tree.clone();
        let Some(tree) = tree else { return };
        if c.midpoint || c.ladderize {
            let done = if c.midpoint {
                through_tools(&tree, |t| {
                    treetools::midpoint_root(t).map(|(rooted, _)| rooted)
                })
            } else {
                through_tools(&tree, |t| Ok(treetools::ladderize(t)))
            };
            match done {
                Ok(node) => {
                    self.tree_win.tree = Some(node);
                    self.platform.trees.selected = None;
                    self.platform.trees.message = if c.midpoint {
                        "Rooted at the midpoint."
                    } else {
                        "Ladderized."
                    }
                    .into();
                }
                Err(e) => self.platform.trees.message = e,
            }
        }
        let page = &mut self.platform.trees;
        if c.fit {
            page.reset_view();
        }
        if c.copy {
            ctx.copy_text(tree.newick());
            page.message = "Copied the tree as Newick.".into();
        }
        if c.save_newick {
            page.message = save_text("Save the tree", "tree.nwk", &tree.newick());
        }
        if c.save_svg {
            let laid = lay_out(&tree, page.options.layout);
            let o = SvgOptions {
                labels: page.options.labels,
                support: page.options.support,
                scale_bar: page.options.scale_bar,
            };
            page.message = save_text("Save the picture", "tree.svg", &to_svg(&laid, o));
        }
        if c.to_phylogeny {
            out.push(Action::Go(Workspace::Phylogeny));
        }
    }

    /// Show a tree on the page: it becomes the platform's current tree.
    pub(in crate::ui) fn show_tree(&mut self, tree: Node, message: &str) {
        self.tree_win.tree = Some(tree);
        self.platform.trees.reset_view();
        self.platform.trees.message = message.to_string();
    }

    fn open_tree_file(&mut self) {
        let ask = dialog::FileDialog::new().set_title("Open a tree");
        if let Some(path) = dialog::pick_one(Wanted::Tree, ask) {
            self.load_tree(Source::Path(&path));
        }
    }

    /// Read a Newick file and draw its first tree.
    pub(in crate::ui) fn load_tree(&mut self, file: Source) {
        let shown = file.shown();
        let text = match file.text() {
            Ok(t) => t,
            Err(e) => {
                self.platform.trees.message = format!("Could not read {shown}: {e}");
                return;
            }
        };
        match NwTree::parse_all(&text) {
            Ok(mut trees) => {
                let n = trees.len();
                let note = format!(
                    "Opened {shown} ({n} tree{}; the first is drawn).",
                    if n == 1 { "" } else { "s" }
                );
                self.show_tree(trees.remove(0).to_node(), &note);
            }
            Err(e) => self.platform.trees.message = format!("{shown}: {e}"),
        }
    }
}

/// Ask where to save `text` and write it; what to tell the person (nothing when they cancelled).
fn save_text(title: &str, name: &str, text: &str) -> String {
    let ask = dialog::FileDialog::new().set_title(title);
    dialog::save_bytes(ask, name, text.as_bytes()).map_or_else(String::new, |saved| saved.message)
}

fn options_panel(ui: &mut Ui, look: &Look, page: &mut TreesPage, have_tree: bool, c: &mut Clicks) {
    design::section(ui, look, "Tree");
    ui.horizontal_wrapped(|ui| {
        c.open = ui.button("Open Newick…").clicked();
        c.example = ui.button("Example").clicked();
    });
    ui.add(
        egui::TextEdit::multiline(&mut page.paste)
            .hint_text("Paste a Newick tree")
            .font(egui::TextStyle::Monospace)
            .desired_rows(2)
            .desired_width(f32::INFINITY),
    );
    c.draw_pasted = ui
        .add_enabled(
            !page.paste.trim().is_empty(),
            egui::Button::new("Draw pasted tree"),
        )
        .clicked();
    ui.add_enabled_ui(have_tree, |ui| {
        ui.horizontal_wrapped(|ui| {
            c.midpoint = ui.button("Root at midpoint").clicked();
            c.ladderize = ui.button("Ladderize").clicked();
        });
    });
    c.to_phylogeny = ui.small_button("Estimate a tree instead").clicked();

    design::section(ui, look, "Layout");
    let o = &mut page.options;
    ui.radio_value(
        &mut o.layout,
        Layout::Phylogram,
        "Phylogram: branch lengths to scale",
    );
    ui.radio_value(
        &mut o.layout,
        Layout::Cladogram,
        "Cladogram: branching only",
    );
    ui.radio_value(&mut o.layout, Layout::Circular, "Circular");

    design::section(ui, look, "Show");
    ui.checkbox(&mut o.labels, "Tip names");
    ui.checkbox(&mut o.support, "Support values");
    ui.add_enabled(
        o.layout != Layout::Cladogram && o.layout != Layout::Circular,
        egui::Checkbox::new(&mut o.scale_bar, "Scale bar"),
    );
    ui.checkbox(&mut o.groups, "Colour tips by the start of their names")
        .on_hover_text(
            "East_01 and East_02 share a colour; used when there are 2 to 8 such groups",
        );
    ui.add(egui::Slider::new(&mut o.font, 8.0..=22.0).text("Text"));
    ui.add(egui::Slider::new(&mut o.line, 0.8..=4.0).text("Lines"));
    ui.add_space(6.0);
    ui.label(
        RichText::new("Scroll to zoom, drag to move, click a node to pick its clade.")
            .size(12.0)
            .color(look.dim),
    );
    if !page.message.is_empty() {
        ui.add_space(4.0);
        ui.label(RichText::new(&page.message).size(12.0).color(look.dim));
    }
}

fn info_panel(ui: &mut Ui, look: &Look, page: &mut TreesPage, laid: &Laid) {
    let root = &laid.nodes[0];
    design::section(ui, look, "This tree");
    ui.horizontal_wrapped(|ui| {
        design::chip(ui, look, &format!("{} tips", root.tips), Tone::Accent);
        design::chip(
            ui,
            look,
            &format!("{} nodes", laid.nodes.len()),
            Tone::Neutral,
        );
        let rooted = if root.children.len() == 2 {
            "rooted"
        } else {
            "unrooted"
        };
        design::chip(ui, look, rooted, Tone::Neutral);
    });
    let depth = laid.nodes.iter().map(|n| n.radius).fold(0.0, f32::max);
    ui.label(
        RichText::new(format!("Longest root-to-tip distance {depth:.4}"))
            .size(12.0)
            .color(look.dim),
    );

    design::section(ui, look, "Picked clade");
    let Some(i) = page.selected.filter(|&i| i < laid.nodes.len()) else {
        ui.label(
            RichText::new("Click a node, or a tip name.")
                .size(12.5)
                .color(look.dim),
        );
        return;
    };
    let node = &laid.nodes[i];
    let tips = laid.tips_under(i);
    ui.horizontal_wrapped(|ui| {
        design::chip(ui, look, &format!("{} tips", tips.len()), Tone::Accent);
        if let Some(s) = node.support {
            design::chip(ui, look, &format!("support {s:.0}"), Tone::Good);
        }
    });
    ScrollArea::vertical()
        .max_height(NAMES_HEIGHT)
        .auto_shrink([false, true])
        .show(ui, |ui| {
            for &t in &tips {
                ui.label(
                    RichText::new(&laid.nodes[t].name)
                        .size(12.5)
                        .color(look.text),
                );
            }
        });
    if ui.button("Copy these names").clicked() {
        let names: Vec<&str> = tips.iter().map(|&t| laid.nodes[t].name.as_str()).collect();
        ui.ctx().copy_text(names.join("\n"));
        page.message = format!("Copied {} names.", names.len());
    }
}

fn canvas(ui: &mut Ui, look: &Look, page: &mut TreesPage, laid: Option<&Laid>) {
    let Some(laid) = laid else {
        design::card(ui, look, |ui| {
            ui.label(
                RichText::new("No tree yet")
                    .size(17.0)
                    .strong()
                    .color(look.text),
            );
            ui.label(
                RichText::new(
                    "Open a Newick file, paste a tree, or try the example on the left. A tree estimated on the Phylogeny \
                     page shows up here too.",
                )
                .color(look.dim),
            );
        });
        return;
    };
    egui::Frame::new()
        .fill(look.surface)
        .stroke(egui::Stroke::new(1.0, look.border))
        .corner_radius(egui::CornerRadius::same(design::CARD_RADIUS))
        .show(ui, |ui| {
            let (options, camera, selected) = (page.options, &mut page.camera, &mut page.selected);
            treeview::show(ui, look, laid, &options, camera, selected, true);
        });
}
