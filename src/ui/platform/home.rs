//! The Home page: where to start. A search box that opens the command palette, the things people do most, a card for each
//! analysis, the routes a project usually takes, and which programs were found on this computer.

use super::{icons, Workspace};
use crate::pipeline::tools::{short_version, Installed};
use crate::remote::settings::Compute;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone, PAGE_MARGIN};
use crate::ui::{brand, topbar};
use eframe::egui::{
    self, pos2, vec2, Frame, Margin, Rect, RichText, ScrollArea, Sense, Shape, Stroke, Ui,
};

const CONTENT_WIDTH: f32 = 1180.0;
const CARD_GAP: f32 = 14.0;
const CARD_HEIGHT: f32 = 196.0;
const MIN_CARD_WIDTH: f32 = 270.0;
const ICON_TILE: f32 = 40.0;
/// Room kept for the buttons at the bottom of a card.
const BUTTONS_HEIGHT: f32 = 34.0;
/// The mark on the hero, and the widest the search box of the hero is.
const HERO_MARK: f32 = 132.0;
const HERO_SEARCH_WIDTH: f32 = 560.0;

/// The analyses the page offers, in the order of a typical project: the workspace, a title, what it does, and
/// whether it can start from an example.
const CARDS: [(Workspace, &str, &str, bool); 8] = [
    (
        Workspace::Signals,
        "Raw signals",
        "Open Sanger chromatograms (.ab1), check the quality of every base call, trim the weak ends and send the \
         clean sequence on. Sequencing reads (FASTQ) are summarised and cleaned here too.",
        true,
    ),
    (
        Workspace::Alignment,
        "Multiple alignment",
        "Align many sequences with the built-in aligner or MAFFT and read them in a colour-coded viewer with \
         conservation, consensus and identity.",
        true,
    ),
    (
        Workspace::Phylogeny,
        "Phylogenetic estimation",
        "Neighbour-joining, UPGMA, maximum likelihood with model choice, and Bayesian sampling, with bootstrap \
         support, run on a worker so the window stays free.",
        true,
    ),
    (
        Workspace::Trees,
        "Tree visualisation",
        "Draw trees as phylograms, cladograms or circles, with support values and a scale bar; explore clades and \
         save the picture as SVG.",
        true,
    ),
    (
        Workspace::Structure,
        "Protein 3D structure",
        "Predict a structure from an amino-acid sequence with ESMFold, or open a PDB file, and turn it around in \
         3D coloured by confidence, secondary structure or chain.",
        true,
    ),
    (
        Workspace::Sequences,
        "Sequences and plasmids",
        "The plasmid editor: maps, restriction enzymes, ORFs, features, primers and PCR, and the contig assembler.",
        false,
    ),
    (
        Workspace::Pipelines,
        "Pipelines",
        "Align, build a tree and pick primers in one run with the programs you have installed, or export the same \
         steps as Snakemake and Nextflow files.",
        true,
    ),
    (
        Workspace::Workstation,
        "Workstation",
        "Send the heavy work to a computer with GPUs: alignments, trees, pipelines and structure predictions are \
         computed there, and the results come back here.",
        false,
    ),
];

/// The routes a project usually takes, as steps to click through.
const ROUTES: [(&str, &[Workspace]); 2] = [
    (
        "From signal to tree",
        &[
            Workspace::Signals,
            Workspace::Alignment,
            Workspace::Phylogeny,
            Workspace::Trees,
        ],
    ),
    (
        "From gene to fold",
        &[Workspace::Sequences, Workspace::Structure],
    ),
];

/// What the hero offers to do at once: the words on the button and what it does.
const QUICK: [(&str, QuickAction); 4] = [
    ("Open a file…", QuickAction::Open),
    (
        "Look at a Sanger trace",
        QuickAction::Example(Workspace::Signals),
    ),
    (
        "Align some sequences",
        QuickAction::Example(Workspace::Alignment),
    ),
    ("Fold a protein", QuickAction::Example(Workspace::Structure)),
];

#[derive(Clone, Copy)]
enum QuickAction {
    Open,
    Example(Workspace),
}

/// A chip for a program: its name and version when it was found, else `missing`.
pub(super) fn program_chip(
    ui: &mut Ui,
    look: &Look,
    label: &str,
    installed: Option<&Installed>,
    missing: &str,
) {
    let Some(p) = installed else {
        return design::chip(ui, look, missing, Tone::Neutral);
    };
    let version = short_version(&p.version);
    let text = if version.is_empty() {
        format!("{label} found")
    } else {
        format!("{label} {version}")
    };
    design::chip(ui, look, &text, Tone::Good);
}

impl App {
    pub(super) fn home_page(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let look = design::look(self.prefs.dark);
        self.detect_programs(ui.ctx());
        ScrollArea::vertical()
            .auto_shrink([false, false])
            .show(ui, |ui| {
                let side = ((ui.available_width() - CONTENT_WIDTH) / 2.0)
                    .clamp(f32::from(PAGE_MARGIN), 120.0) as i8;
                Frame::new()
                    .inner_margin(Margin {
                        left: side,
                        right: side,
                        top: 22,
                        bottom: 22,
                    })
                    .show(ui, |ui| {
                        hero(ui, &look, out);
                        ui.add_space(22.0);
                        heading(ui, &look, "Start with an analysis");
                        ui.add_space(8.0);
                        cards(ui, &look, out);
                        ui.add_space(8.0);
                        routes(ui, &look, out);
                        ui.add_space(CARD_GAP);
                        self.this_computer(ui, &look);
                    });
            });
    }

    /// Which programs were found here, and where the heavy jobs go.
    fn this_computer(&self, ui: &mut Ui, look: &Look) {
        design::card(ui, look, |ui| {
            ui.label(
                RichText::new("On this computer")
                    .size(15.0)
                    .strong()
                    .color(look.text),
            );
            ui.add_space(4.0);
            ui.horizontal_wrapped(|ui| {
                match &self.platform.programs {
                    None => {
                        design::chip(ui, look, "Looking for installed programs…", Tone::Neutral)
                    }
                    Some(found) => {
                        for (label, installed) in found {
                            let missing = format!("{label} not found");
                            program_chip(ui, look, label, installed.as_ref(), &missing);
                        }
                    }
                }
                match self.platform.workstation.compute {
                    Compute::Local => {}
                    Compute::Workstation => {
                        design::chip(ui, look, "Heavy jobs go to the workstation", Tone::Accent);
                    }
                    Compute::Automatic => design::chip(
                        ui,
                        look,
                        "Big jobs go to the workstation, small ones stay here",
                        Tone::Accent,
                    ),
                }
            });
            ui.label(
                RichText::new("Optional: the built-in methods are used when a program is missing.")
                    .size(12.0)
                    .color(look.dim),
            );
        });
    }
}

/// A heading above a group of cards.
fn heading(ui: &mut Ui, look: &Look, text: &str) {
    ui.label(RichText::new(text).size(19.0).strong().color(look.text));
}

/// The top of the page: what the program is for, the search box that opens the command palette, and what people do most.
fn hero(ui: &mut Ui, look: &Look, out: &mut Vec<Action>) {
    design::card(ui, look, |ui| {
        ui.set_width(ui.available_width());
        let area = ui.max_rect();
        let mark = Rect::from_min_size(
            pos2(area.right() - HERO_MARK - 10.0, area.top() + 6.0),
            vec2(HERO_MARK, HERO_MARK),
        );
        brand::paint(ui.painter(), mark, look);
        ui.add_space(6.0);
        ui.label(
            RichText::new("A WORKBENCH FOR SEQUENCES, SIGNALS AND STRUCTURES")
                .size(11.5)
                .strong()
                .color(look.accent),
        );
        ui.label(
            RichText::new("What will you work on today?")
                .size(32.0)
                .strong()
                .color(look.text),
        );
        ui.label(
            RichText::new(
                "Edit plasmids, check Sanger traces, align, build trees and predict structures, on this computer or on \
                 your workstation.",
            )
            .size(15.0)
            .color(look.dim),
        );
        ui.add_space(14.0);
        let width = (ui.available_width() - HERO_MARK - 24.0).min(HERO_SEARCH_WIDTH);
        if topbar::search_field(ui, look, vec2(width, 44.0), true).clicked() {
            out.push(Action::Palette);
        }
        ui.add_space(12.0);
        ui.horizontal_wrapped(|ui| {
            for (text, what) in QUICK {
                if design::ghost_button(ui, look, text).clicked() {
                    out.push(match what {
                        QuickAction::Open => Action::Open,
                        QuickAction::Example(ws) => Action::Example(ws),
                    });
                }
            }
        });
        ui.add_space(6.0);
    });
}

fn cards(ui: &mut Ui, look: &Look, out: &mut Vec<Action>) {
    let available = ui.available_width();
    let fits = |n: f32| available >= n * MIN_CARD_WIDTH + (n - 1.0) * CARD_GAP;
    let columns = [4, 3, 2].into_iter().find(|&n| fits(n as f32)).unwrap_or(1);
    let width = (available - (columns as f32 - 1.0) * CARD_GAP) / columns as f32;
    for row in CARDS.chunks(columns) {
        ui.horizontal_top(|ui| {
            ui.spacing_mut().item_spacing.x = CARD_GAP;
            for &(ws, title, text, example) in row {
                let down = egui::Layout::top_down(egui::Align::Min);
                ui.allocate_ui_with_layout(vec2(width, CARD_HEIGHT), down, |ui| {
                    design::card(ui, look, |ui| {
                        ui.set_min_size(vec2(width - 32.0, CARD_HEIGHT - 32.0));
                        ui.set_max_width(width - 32.0);
                        card(ui, look, out, (ws, title, text, example));
                    });
                });
            }
        });
        ui.add_space(CARD_GAP);
    }
}

fn card(
    ui: &mut Ui,
    look: &Look,
    out: &mut Vec<Action>,
    (ws, title, text, example): (Workspace, &str, &str, bool),
) {
    ui.horizontal(|ui| {
        design::icon_tile(ui, look, ICON_TILE, |painter, rect| {
            icons::paint(painter, ws, rect, look.text, look.accent, look.dark);
        });
        ui.label(RichText::new(title).size(16.5).strong().color(look.text));
    });
    ui.add_space(2.0);
    ui.label(RichText::new(text).size(13.0).color(look.dim));
    let rest = ui.available_height() - BUTTONS_HEIGHT;
    if rest > 0.0 {
        ui.add_space(rest);
    }
    ui.horizontal(|ui| {
        if design::primary_button(ui, look, "Open").clicked() {
            out.push(Action::Go(ws));
        }
        if example && design::ghost_button(ui, look, "Try an example").clicked() {
            out.push(Action::Example(ws));
        }
    });
}

fn routes(ui: &mut Ui, look: &Look, out: &mut Vec<Action>) {
    design::card(ui, look, |ui| {
        ui.label(
            RichText::new("Follow a project")
                .size(16.5)
                .strong()
                .color(look.text),
        );
        for (name, steps) in ROUTES {
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                ui.label(RichText::new(name).size(13.0).color(look.dim));
                ui.add_space(8.0);
                for (k, &ws) in steps.iter().enumerate() {
                    if k > 0 {
                        chevron(ui, look);
                    }
                    let step =
                        egui::Button::new(RichText::new(ws.label()).strong().color(look.accent))
                            .fill(look.accent_soft)
                            .stroke(Stroke::NONE)
                            .corner_radius(design::BUTTON_RADIUS);
                    if ui.add(step).on_hover_text(ws.blurb()).clicked() {
                        out.push(Action::Go(ws));
                    }
                }
            });
        }
        ui.add_space(6.0);
        ui.label(
            RichText::new("Drop a sequence, tree, trace (.ab1) or structure (.pdb) file on the window and it opens in the right place.")
                .size(12.0)
                .color(look.dim),
        );
    });
}

/// A small arrow head between two steps.
fn chevron(ui: &mut Ui, look: &Look) {
    let (r, _) = ui.allocate_exact_size(vec2(14.0, 18.0), Sense::hover());
    let c = r.center();
    ui.painter().add(Shape::line(
        vec![
            c + vec2(-3.0, -5.0),
            c + vec2(3.0, 0.0),
            c + vec2(-3.0, 5.0),
        ],
        Stroke::new(1.8, look.dim),
    ));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_workspace_but_home_is_offered_once() {
        let offered: Vec<Workspace> = CARDS.iter().map(|c| c.0).collect();
        for ws in Workspace::ALL {
            let expected = usize::from(ws != Workspace::Home);
            assert_eq!(
                offered.iter().filter(|&&o| o == ws).count(),
                expected,
                "{ws:?}"
            );
        }
    }

    #[test]
    fn the_quick_actions_lead_somewhere_that_exists() {
        for (_, what) in QUICK {
            if let QuickAction::Example(ws) = what {
                assert!(
                    CARDS.iter().any(|c| c.0 == ws && c.3),
                    "{ws:?} has no example"
                );
            }
        }
    }
}
