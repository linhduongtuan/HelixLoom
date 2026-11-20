//! The controls that the Trees page has for annotated trees: what to colour, widen and label by, an axis, a search, and the edits.

use super::looks::SupportFrom;
use super::{Clicks, TreesPage};
use crate::seq::treestyle::{Axis, AxisMode};
use crate::ui::design::{self, Look};
use eframe::egui::{self, RichText, Ui};

/// A section that can be folded away, so that the long list of controls stays short; it is open when `open`.
fn fold(ui: &mut Ui, look: &Look, title: &str, open: bool, body: impl FnOnce(&mut Ui)) {
    ui.add_space(4.0);
    egui::CollapsingHeader::new(
        RichText::new(title.to_uppercase())
            .size(11.0)
            .strong()
            .color(look.dim),
    )
    .id_salt(("trees fold", title))
    .default_open(open)
    .show(ui, body);
}

/// A choice of one annotation (or none) from `keys`.
fn pick(ui: &mut Ui, id: &str, label: &str, current: &mut Option<String>, keys: &[String]) {
    ui.horizontal(|ui| {
        ui.label(label);
        egui::ComboBox::from_id_salt(id)
            .selected_text(current.as_deref().unwrap_or("nothing"))
            .show_ui(ui, |ui| {
                ui.selectable_value(current, None, "nothing");
                for k in keys {
                    ui.selectable_value(current, Some(k.clone()), k);
                }
            });
    });
}

/// How the tree's annotations show.
pub(super) fn annotations(
    ui: &mut Ui,
    look: &Look,
    page: &mut TreesPage,
    keys: &[String],
    has_labels: bool,
) {
    design::section(ui, look, "Annotations");
    if keys.is_empty() && !has_labels {
        ui.label(
            RichText::new("This tree carries none. A NEXUS or Newick file from BEAST, MrBayes or FigTree has them in square brackets, [&posterior=0.98].")
                .size(12.0)
                .color(look.dim),
        );
        return;
    }
    let s = &mut page.looks.style;
    pick(
        ui,
        "colour_by",
        "Colour branches by",
        &mut s.colour_by,
        keys,
    );
    pick(ui, "width_by", "Widen branches by", &mut s.width_by, keys);
    pick(ui, "shape_by", "Node shapes by", &mut s.shape_by, keys);
    pick(ui, "bars_by", "Bars across nodes", &mut s.bars_by, keys);
    ui.label(RichText::new("A bar is a range of heights, {least,greatest}, as a 95 % HPD interval is written; it shows on a phylogram.").size(11.5).color(look.dim));
    ui.label("Beside the inner nodes (several are joined with a slash):");
    ui.horizontal_wrapped(|ui| {
        let mut all: Vec<String> = Vec::new();
        if has_labels {
            all.push(String::new());
        }
        all.extend(keys.iter().cloned());
        for k in all {
            let mut on = s.labels.contains(&k);
            let name = if k.is_empty() {
                "node label"
            } else {
                k.as_str()
            };
            if ui.checkbox(&mut on, name).changed() {
                if on {
                    s.labels.push(k);
                } else {
                    s.labels.retain(|x| *x != k);
                }
            }
        }
    });
    ui.checkbox(
        &mut page.looks.tip_colours,
        "Colour tip names like their branches",
    );
}

/// The axis under a phylogram.
pub(super) fn axis(ui: &mut Ui, look: &Look, page: &mut TreesPage) {
    let open = page.looks.axis_on;
    fold(ui, look, "Axis", open, |ui| axis_body(ui, look, page));
}

fn axis_body(ui: &mut Ui, look: &Look, page: &mut TreesPage) {
    let on_phylogram = page.options.layout == crate::seq::treelayout::Layout::Phylogram;
    ui.add_enabled_ui(on_phylogram, |ui| {
        ui.checkbox(&mut page.looks.axis_on, "Show an axis under the tree")
            .on_hover_text(
                "A scale of distance or time along the branches, with ticks at round numbers",
            );
        if page.looks.axis_on {
            let a: &mut Axis = &mut page.looks.axis;
            ui.radio_value(&mut a.mode, AxisMode::FromRoot, "Distance from the root");
            ui.radio_value(
                &mut a.mode,
                AxisMode::BeforeTips,
                "Time before the tips (0 at the tips)",
            );
            ui.radio_value(&mut a.mode, AxisMode::Dates, "Dates (the tips at a date)");
            ui.horizontal(|ui| {
                ui.label("Units per length");
                ui.add(egui::DragValue::new(&mut a.scale).speed(0.05));
            });
            if a.mode == AxisMode::Dates {
                ui.horizontal(|ui| {
                    ui.label("Date of the tips");
                    ui.add(egui::DragValue::new(&mut a.tip_date).speed(0.01));
                });
            }
            ui.checkbox(&mut a.grid, "Lines up from the ticks");
        }
    });
    if !on_phylogram {
        ui.label(
            RichText::new("The axis is for the phylogram layout.")
                .size(12.0)
                .color(look.dim),
        );
    }
}

/// Search the tip names.
pub(super) fn find(ui: &mut Ui, look: &Look, page: &mut TreesPage, found: usize) {
    design::section(ui, look, "Find tips");
    ui.add(egui::TextEdit::singleline(&mut page.looks.find).hint_text("Part of a name"));
    if !page.looks.find.trim().is_empty() {
        ui.label(
            RichText::new(format!(
                "{found} tip{} match",
                if found == 1 { "" } else { "s" }
            ))
            .size(12.0)
            .color(look.dim),
        );
    }
}

/// The edits that work on the whole tree.
pub(super) fn edit(
    ui: &mut Ui,
    look: &Look,
    page: &mut TreesPage,
    keys: &[String],
    has_labels: bool,
    c: &mut Clicks,
) {
    design::section(ui, look, "Edit the tree");
    ui.horizontal_wrapped(|ui| {
        c.order_up = ui.button("Order: smaller first").clicked();
        c.order_down = ui.button("Order: larger first").clicked();
    });
    ui.label(
        RichText::new("Take out branches whose value is below:")
            .size(12.0)
            .color(look.dim),
    );
    ui.horizontal(|ui| {
        egui::ComboBox::from_id_salt("weak_key")
            .selected_text(if page.looks.weak_key.is_empty() {
                "node label"
            } else {
                page.looks.weak_key.as_str()
            })
            .show_ui(ui, |ui| {
                if has_labels {
                    ui.selectable_value(&mut page.looks.weak_key, String::new(), "node label");
                }
                for k in keys {
                    ui.selectable_value(&mut page.looks.weak_key, k.clone(), k);
                }
            });
        ui.add(egui::DragValue::new(&mut page.looks.weak_below).speed(0.5));
        c.weak = ui
            .button("Collapse")
            .on_hover_text("Their branches go to the parent, with the length added")
            .clicked();
    });
    ui.label(
        RichText::new("Rename tips (old name, a tab or comma, new name; a line each):")
            .size(12.0)
            .color(look.dim),
    );
    ui.add(
        egui::TextEdit::multiline(&mut page.looks.names_paste)
            .font(egui::TextStyle::Monospace)
            .desired_rows(2)
            .desired_width(f32::INFINITY),
    );
    c.rename = ui
        .add_enabled(
            !page.looks.names_paste.trim().is_empty(),
            egui::Button::new("Rename"),
        )
        .clicked();
}

/// Support values taken from other trees.
pub(super) fn support(ui: &mut Ui, look: &Look, page: &mut TreesPage, c: &mut Clicks) {
    let open = page.looks.awaiting_support;
    fold(ui, look, "Support from other trees", open, |ui| {
        support_body(ui, look, page, c)
    });
}

fn support_body(ui: &mut Ui, look: &Look, page: &mut TreesPage, c: &mut Clicks) {
    ui.horizontal(|ui| {
        ui.label("Name");
        ui.add(
            egui::TextEdit::singleline(&mut page.looks.support_label)
                .hint_text("support")
                .desired_width(110.0),
        );
    });
    ui.radio_value(
        &mut page.looks.support_from,
        SupportFrom::Share,
        "Share of the trees that have the branch",
    );
    ui.radio_value(
        &mut page.looks.support_from,
        SupportFrom::Written,
        "The value on the same branch of one tree",
    );
    if page.looks.support_from == SupportFrom::Written {
        ui.horizontal(|ui| {
            ui.label("Annotation");
            ui.add(
                egui::TextEdit::singleline(&mut page.looks.support_key)
                    .hint_text("empty: the node label")
                    .desired_width(150.0),
            );
        });
    }
    c.support = ui
        .button("Open the other trees…")
        .on_hover_text("Bootstrap replicates, posterior samples, or the tree of another analysis (Newick or NEXUS)")
        .clicked();
    if page.looks.awaiting_support {
        ui.label(
            RichText::new("Waiting for a file of trees.")
                .size(12.0)
                .color(look.dim),
        );
    }
}

/// What can be done to the picked clade.
pub(super) fn picked(
    ui: &mut Ui,
    look: &Look,
    page: &TreesPage,
    c: &mut Clicks,
    is_collapsed: bool,
    can_root: bool,
) {
    ui.horizontal_wrapped(|ui| {
        c.reroot = ui
            .add_enabled(can_root, egui::Button::new("Root here"))
            .clicked();
        c.rotate = ui.button("Rotate").clicked();
        c.fold = ui
            .button(if is_collapsed {
                "Expand"
            } else {
                "Collapse to a triangle"
            })
            .clicked();
        c.extract = ui.button("Show only this clade").clicked();
    });
    ui.horizontal_wrapped(|ui| {
        ui.label(RichText::new("Highlight").size(12.0).color(look.dim));
        for (k, rgb) in crate::seq::treestyle::PALETTE.iter().enumerate() {
            let colour = egui::Color32::from_rgb(rgb[0], rgb[1], rgb[2]);
            let button = egui::Button::new("  ")
                .fill(colour)
                .min_size(egui::vec2(18.0, 18.0));
            if ui.add(button).clicked() {
                c.colour = Some(k);
            }
        }
        c.unmark = ui.small_button("none").clicked();
    });
    let _ = page;
}
