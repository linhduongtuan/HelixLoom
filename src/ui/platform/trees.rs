//! Trees: draw a tree as a phylogram, a cladogram or a circle, with support values and a scale bar, pick a clade
//! to see its tips, and save the picture as SVG. The tree is the one the Phylogeny page estimated, or one opened
//! from a Newick file or pasted.

use super::route::Source;
use super::treeview::{self, Camera, Dress, Options};
use super::Workspace;
use crate::seq::newick::NwTree;
use crate::seq::nexus;
use crate::seq::tree::Node;
use crate::seq::treelayout::{lay_out_with, to_svg, Laid, Layout, SvgOptions};
use crate::seq::treestyle::{decorate, Decor};
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::dialog::{self, Wanted};
use eframe::egui::{self, RichText, ScrollArea, Ui};

mod edits;
mod looks;
mod panels;
#[cfg(test)]
mod tests;

pub(in crate::ui) use looks::Looks;
use looks::{found_tips, laid_index, node_id, notes, tip_names_under, Annotated};

const INFO_WIDTH: f32 = 260.0;
/// The least width a tree is drawn in beside the card that says what it is.
const TREE_MIN_WIDTH: f32 = 360.0;
const NAMES_HEIGHT: f32 = 260.0;
/// The card that says what the tree is, when it is under the tree: its share of the height, and the least and the most.
const UNDER_SHARE: f32 = 0.4;
const UNDER_LEAST: f32 = 110.0;
const UNDER_MOST: f32 = 200.0;

#[derive(Default)]
pub struct TreesPage {
    pub options: Options,
    camera: Camera,
    pub(in crate::ui) selected: Option<usize>,
    paste: String,
    pub message: String,
    /// The tree with its annotations, and how it is dressed and marked.
    doc: Option<Annotated>,
    pub(in crate::ui) looks: Looks,
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
    save_nexus: bool,
    copy: bool,
    fit: bool,
    to_phylogeny: bool,
    reroot: bool,
    rotate: bool,
    extract: bool,
    order_up: bool,
    order_down: bool,
    weak: bool,
    rename: bool,
    support: bool,
    /// Collapse or expand the picked clade.
    fold: bool,
    /// The palette colour (by its place) to highlight the picked clade with, and to take the highlight off.
    colour: Option<usize>,
    unmark: bool,
}

impl App {
    pub(super) fn trees_page(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let look = self.prefs.look();
        self.sync_annotated();
        let tree = self.tree_win.tree.clone();
        let mut c = Clicks::default();
        let page = &mut self.platform.trees;
        let laid = tree
            .as_ref()
            .map(|t| lay_out_with(t, page.options.layout, &page.looks.collapsed));
        let decor = match (&laid, page.annotated()) {
            (Some(laid), Some(nw)) => Some(decorate(nw, laid, &page.looks.style)),
            _ => None,
        };
        let drawn = laid.as_ref().zip(page.annotated()).map(|(laid, nw)| Drawn {
            notes: notes(nw, laid),
            found: found_tips(laid, &page.looks.find),
            highlights: page
                .looks
                .highlights
                .iter()
                .filter_map(|&(o, c)| Some((laid_index(laid, o)?, c)))
                .collect(),
        });
        let keys = page.annotated().map(NwTree::attr_keys).unwrap_or_default();
        let has_labels = page.annotated().is_some_and(|nw| {
            nw.nodes
                .iter()
                .any(|n| !n.children.is_empty() && !n.name.is_empty())
        });
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
                        c.save_nexus = ui
                            .button("Save NEXUS…")
                            .on_hover_text("With the annotations of the nodes, as BEAST and FigTree write them")
                            .clicked();
                        c.copy = ui.button("Copy Newick").clicked();
                        c.fit = ui.button("Fit to window").clicked();
                        crate::ui::figure::button(ui);
                    });
                },
            );
            // a window too narrow for the tree between two cards: the card that says what the tree is, and what was picked
            // in it, goes under the tree, where a click on the tree still shows its answer without scrolling for it
            let for_the_tree = ui.available_width()
                - design::CONTROLS_WIDTH
                - INFO_WIDTH
                - 2.0 * f32::from(design::PAGE_MARGIN);
            let under_the_tree = for_the_tree < TREE_MIN_WIDTH;
            let found = drawn.as_ref().map_or(0, |d| d.found.len());
            design::controls_panel(ui, &look, "trees_options", |ui| {
                options_panel(
                    ui,
                    &look,
                    page,
                    (&keys, has_labels, found),
                    tree.is_some(),
                    &mut c,
                )
            });
            let apart = |left: i8, top: i8| {
                egui::Frame::new().inner_margin(egui::Margin {
                    left,
                    top,
                    ..egui::Margin::ZERO
                })
            };
            match &laid {
                Some(laid) if under_the_tree => {
                    // of one height whatever it holds: picking a clade must not make the tree smaller under the pointer
                    // that picked it, and a long list of names scrolls in the card, it does not take the tree's room
                    let high = (ui.available_height() * UNDER_SHARE).clamp(UNDER_LEAST, UNDER_MOST);
                    egui::Panel::bottom("trees_info_under")
                        .exact_size(high)
                        .resizable(false)
                        .frame(apart(design::PAGE_MARGIN, design::CARD_GAP as i8))
                        .show(ui, |ui| {
                            design::card(ui, &look, |ui| {
                                ui.set_min_width(ui.available_width());
                                ScrollArea::vertical()
                                    .id_salt("trees_info_under_scroll")
                                    .auto_shrink([false, false])
                                    .show(ui, |ui| info_panel(ui, &look, page, laid, &mut c));
                            });
                        });
                }
                Some(laid) => {
                    egui::Panel::right("trees_info")
                        .exact_size(INFO_WIDTH)
                        .resizable(false)
                        .frame(apart(design::PAGE_MARGIN, 0))
                        .show(ui, |ui| {
                            ScrollArea::vertical()
                                .id_salt("trees_info_scroll")
                                .auto_shrink([false, false])
                                .show(ui, |ui| {
                                    design::card(ui, &look, |ui| {
                                        info_panel(ui, &look, page, laid, &mut c)
                                    });
                                });
                        });
                }
                None => {}
            }
            design::results_area(ui, |ui| {
                canvas(
                    ui,
                    &look,
                    page,
                    laid.as_ref(),
                    (decor.as_ref(), drawn.as_ref()),
                    &mut c,
                )
            });
        });
        self.apply_tree_clicks(c, ui.ctx(), out);
    }

    fn apply_tree_clicks(&mut self, c: Clicks, ctx: &egui::Context, out: &mut Vec<Action>) {
        if c.open {
            self.open_tree_file();
        }
        if c.support {
            self.platform.trees.looks.awaiting_support = true;
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
            match nexus::parse_any(&text) {
                Ok(trees) => {
                    let n = trees.len();
                    if self.platform.trees.looks.awaiting_support {
                        self.add_support(&trees, "the pasted text");
                    } else {
                        self.show_annotated(
                            trees,
                            &format!("Drew the first of {n} pasted tree(s)."),
                        );
                    }
                }
                Err(e) => self.platform.trees.message = e,
            }
        }
        self.apply_tree_edits(&c);
        let tree = self.tree_win.tree.clone();
        let Some(tree) = tree else { return };
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
        if c.save_nexus {
            if let Some(nw) = page.annotated() {
                let text = nexus::write(nw, "tree", true);
                page.message = save_text("Save the tree", "tree.nex", &text);
            }
        }
        if c.save_svg {
            let laid = lay_out_with(&tree, page.options.layout, &page.looks.collapsed);
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
        // a tree from elsewhere is not the one the Phylogeny page estimated: its record goes with it, and so do the lines that say
        // how that one was made (a caller that has such lines for this tree sets them after this)
        self.platform.phylogeny.made = None;
        self.tree_win.details.clear();
        self.tree_win.set_tree(Some(tree));
        self.platform.trees.reset_view();
        self.platform.trees.message = message.to_string();
    }

    fn open_tree_file(&mut self) {
        let ask = dialog::FileDialog::new().set_title("Open a tree");
        if let Some(path) = dialog::pick_one(Wanted::Tree, ask) {
            self.load_tree(Source::Path(&path));
        }
    }

    /// Read a Newick or NEXUS file and draw its first tree (or, when the page waits for the trees that give support, take the support
    /// from them).
    pub(in crate::ui) fn load_tree(&mut self, file: Source) {
        let shown = file.shown();
        let text = match file.text() {
            Ok(t) => t,
            Err(e) => {
                self.platform.trees.message = format!("Could not read {shown}: {e}");
                self.platform.trees.looks.awaiting_support = false;
                return;
            }
        };
        match nexus::parse_any(&text) {
            Ok(trees) if self.platform.trees.looks.awaiting_support => {
                self.add_support(&trees, &shown)
            }
            Ok(trees) => {
                let n = trees.len();
                let note = format!(
                    "Opened {shown} ({n} tree{}; the first is drawn).",
                    if n == 1 { "" } else { "s" }
                );
                self.show_annotated(trees, &note);
            }
            Err(e) => {
                self.platform.trees.looks.awaiting_support = false;
                self.platform.trees.message = format!("{shown}: {e}");
            }
        }
    }
}

/// What the viewer needs besides the layout and the decoration, worked out for the frame.
struct Drawn {
    notes: Vec<String>,
    found: Vec<usize>,
    highlights: Vec<(usize, crate::seq::treestyle::Rgb)>,
}

/// Ask where to save `text` and write it; what to tell the person (nothing when they cancelled).
pub(in crate::ui) fn save_text(title: &str, name: &str, text: &str) -> String {
    let ask = dialog::FileDialog::new().set_title(title);
    dialog::save_bytes(ask, name, text.as_bytes()).map_or_else(String::new, |saved| saved.message)
}

fn options_panel(
    ui: &mut Ui,
    look: &Look,
    page: &mut TreesPage,
    facts: (&[String], bool, usize),
    have_tree: bool,
    c: &mut Clicks,
) {
    let (keys, has_labels, found) = facts;
    design::section(ui, look, "Tree");
    ui.horizontal_wrapped(|ui| {
        c.open = ui.button("Open Newick or NEXUS…").clicked();
        c.example = ui.button("Example").clicked();
    });
    ui.add(
        egui::TextEdit::multiline(&mut page.paste)
            .hint_text("Paste a Newick or NEXUS tree")
            .font(egui::TextStyle::Monospace)
            .desired_rows(2)
            .desired_width(f32::INFINITY),
    );
    c.draw_pasted = ui
        .add_enabled(
            !page.paste.trim().is_empty(),
            egui::Button::new(if page.looks.awaiting_support {
                "Use the pasted trees"
            } else {
                "Draw pasted tree"
            }),
        )
        .clicked();
    ui.add_enabled_ui(have_tree, |ui| {
        ui.horizontal_wrapped(|ui| {
            c.midpoint = ui.button("Root at midpoint").clicked();
            c.ladderize = ui.button("Ladderize").clicked();
        });
    });
    c.to_phylogeny = ui.button("Estimate a tree instead").clicked();

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
    ui.radio_value(
        &mut o.layout,
        Layout::Radial,
        "Radial: unrooted, branches in a fan",
    );

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
    if have_tree {
        panels::annotations(ui, look, page, keys, has_labels);
        panels::axis(ui, look, page);
        panels::find(ui, look, page, found);
        panels::edit(ui, look, page, keys, has_labels, c);
        panels::support(ui, look, page, c);
    }
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

fn info_panel(ui: &mut Ui, look: &Look, page: &mut TreesPage, laid: &Laid, c: &mut Clicks) {
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
    // the names come from the tree itself: a collapsed clade is one node of the layout, and its tips are all there
    let names: Vec<String> = page
        .annotated()
        .and_then(|nw| node_id(nw, laid, i).map(|id| tip_names_under(nw, id)))
        .unwrap_or_else(|| {
            laid.tips_under(i)
                .iter()
                .map(|&t| laid.nodes[t].name.clone())
                .collect()
        });
    ui.horizontal_wrapped(|ui| {
        design::chip(ui, look, &format!("{} tips", names.len()), Tone::Accent);
        if let Some(s) = node.support {
            design::chip(ui, look, &format!("support {s:.0}"), Tone::Good);
        }
    });
    panels::picked(ui, look, page, c, node.cartoon.is_some(), i != 0);
    ScrollArea::vertical()
        .max_height(NAMES_HEIGHT)
        .auto_shrink([false, true])
        .show(ui, |ui| {
            for name in &names {
                ui.label(RichText::new(name).size(12.5).color(look.text));
            }
        });
    if ui.button("Copy these names").clicked() {
        ui.ctx().copy_text(names.join("\n"));
        page.message = format!("Copied {} names.", names.len());
    }
}

fn canvas(
    ui: &mut Ui,
    look: &Look,
    page: &mut TreesPage,
    laid: Option<&Laid>,
    dressed: (Option<&Decor>, Option<&Drawn>),
    c: &mut Clicks,
) {
    let Some(laid) = laid else {
        let pick = design::empty_state(
            ui,
            look,
            "No tree yet",
            "Open a Newick or NEXUS file, paste a tree on the left, or try the example. A tree estimated on the Phylogeny \
             page shows up here too.",
            &["Open Newick or NEXUS…", "Try the example"],
        );
        match pick {
            Some(0) => c.open = true,
            Some(1) => c.example = true,
            _ => {}
        }
        return;
    };
    let (decor, drawn) = dressed;
    let empty = Drawn {
        notes: Vec::new(),
        found: Vec::new(),
        highlights: Vec::new(),
    };
    let drawn = drawn.unwrap_or(&empty);
    let dress = Dress {
        decor,
        highlights: &drawn.highlights,
        found: &drawn.found,
        axis: page.looks.axis_on.then_some(page.looks.axis),
        notes: &drawn.notes,
        tip_colours: page.looks.tip_colours,
    };
    if let Some(format) = crate::ui::figure::wanted() {
        // what the window shows: the same room, zoomed and moved as it is, with what is picked
        let (options, mut camera, mut selected) = (page.options, page.camera, page.selected);
        page.message = crate::ui::figure::save(format, "tree", ui.available_size(), |ui| {
            treeview::show_dressed(
                ui,
                look,
                laid,
                &options,
                &dress,
                &mut camera,
                &mut selected,
                false,
            );
        });
    }
    design::view_frame(ui, look, |ui| {
        let (options, camera, selected) = (page.options, &mut page.camera, &mut page.selected);
        treeview::show_dressed(ui, look, laid, &options, &dress, camera, selected, true);
    });
}
