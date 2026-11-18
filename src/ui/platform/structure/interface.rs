//! The Structure page's Interface section: choose which chains meet, measure what they bury, and see it on the structure.
//! The measuring is `seq::interface`; this is its window.

use super::viewer3d::Scheme;
use crate::seq::interface::{self, Interface, Row};
use crate::ui::app::App;
use crate::ui::design::{self, Look, Tone};
use crate::ui::fonts::FontWeight;
use crate::ui::job::{Job, Poll};
use eframe::egui::{self, RichText, ScrollArea, Ui};
use std::time::Duration;

const WORKING_REDRAW: Duration = Duration::from_millis(200);
/// How tall the table of residues is before it scrolls.
const TABLE_HEIGHT: f32 = 230.0;

/// Which side of the interface a chain is on.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Side {
    /// Left out: not looked at.
    #[default]
    Neither,
    First,
    Second,
}

/// What was measured, with the rows of its table.
pub struct Measured {
    pub interface: Interface,
    pub rows: Vec<Row>,
}

#[derive(Default)]
pub struct InterfacePage {
    /// The side of each chain of the structure shown, in the order the chains appear.
    pub sides: Vec<(char, Side)>,
    pub measured: Option<Measured>,
    pub message: String,
    job: Option<Job<Result<Interface, String>>>,
}

impl InterfacePage {
    pub fn is_running(&self) -> bool {
        self.job.is_some()
    }

    /// For a structure just shown: its first chain on the first side and the others on the second, which is where most
    /// complexes start (one protein and what binds it); the person changes it.
    pub fn reset_for(&mut self, chains: &[char]) {
        self.sides = chains
            .iter()
            .enumerate()
            .map(|(k, &c)| (c, if k == 0 { Side::First } else { Side::Second }))
            .collect();
        self.measured = None;
        self.job = None;
        self.message.clear();
    }

    /// Forget what was measured, and any measuring that is going on; the structure goes back from the colours of the interface
    /// (which would be all one grey) to those of its chains.
    pub fn take_back(&mut self, fold: &mut crate::seq::fold::Fold, style: &mut super::Style) {
        self.measured = None;
        self.job = None;
        fold.clear_interface();
        if style.scheme == Scheme::Interface {
            style.scheme = Scheme::Chain;
        }
    }

    /// The chains in `first` on the first side and every other chain on the second.
    pub fn set_first_side(&mut self, first: &[char]) {
        for (chain, side) in &mut self.sides {
            *side = if first.contains(chain) {
                Side::First
            } else {
                Side::Second
            };
        }
    }

    pub fn chains_on(&self, side: Side) -> Vec<char> {
        self.sides
            .iter()
            .filter(|(_, s)| *s == side)
            .map(|(c, _)| *c)
            .collect()
    }
}

#[derive(Default)]
pub struct InterfaceClicks {
    analyse: bool,
    colour: bool,
    copy: bool,
    /// A chain's side was changed: what was measured is not the answer to the new sides.
    changed: bool,
    /// A residue of the table was clicked: its index in the structure's residues.
    pick: Option<usize>,
}

pub(super) fn panel(ui: &mut Ui, look: &Look, page: &mut InterfacePage, c: &mut InterfaceClicks) {
    design::section(ui, look, "Interface");
    if page.sides.len() < 2 {
        ui.label(
            RichText::new("An interface needs at least two chains; this structure has one.")
                .size(12.0)
                .color(look.dim),
        );
        return;
    }
    ui.label(
        RichText::new("Which chains meet: put each on side 1, side 2, or leave it out.")
            .size(12.0)
            .color(look.dim),
    );
    for (chain, side) in &mut page.sides {
        ui.horizontal(|ui| {
            ui.label(format!("Chain {chain}"));
            for (value, label) in [
                (Side::First, "1"),
                (Side::Second, "2"),
                (Side::Neither, "left out"),
            ] {
                c.changed |= ui.selectable_value(side, value, label).changed();
            }
        });
    }
    let ready = !page.chains_on(Side::First).is_empty()
        && !page.chains_on(Side::Second).is_empty()
        && !page.is_running();
    c.analyse = ui
        .add_enabled(ready, egui::Button::new("Analyse the interface"))
        .on_disabled_hover_text("Put at least one chain on each side")
        .clicked();
    if !page.message.is_empty() {
        ui.label(RichText::new(&page.message).size(12.0).color(look.dim));
    }
    let Some(done) = &page.measured else { return };
    let i = &done.interface;
    ui.add_space(4.0);
    ui.horizontal_wrapped(|ui| {
        design::chip(
            ui,
            look,
            &format!("interface {:.0} Å²", i.interface_area()),
            Tone::Accent,
        );
        design::chip(
            ui,
            look,
            &format!("{} contacts", i.contacts.len()),
            Tone::Neutral,
        );
        design::chip(
            ui,
            look,
            &format!("{} H-bond candidates", i.hbonds.len()),
            Tone::Neutral,
        );
    });
    ui.label(
        RichText::new(format!(
            "Buried: {:.0} Å² in all, {:.0} on side 1 and {:.0} on side 2 (the interface is half of the total).",
            i.buried_total(),
            i.a.iter().map(|b| b.buried()).sum::<f64>(),
            i.b.iter().map(|b| b.buried()).sum::<f64>(),
        ))
        .size(12.0)
        .color(look.dim),
    );
    ScrollArea::vertical()
        .max_height(TABLE_HEIGHT)
        .show(ui, |ui| {
            egui::Grid::new("interface_rows")
                .striped(true)
                .spacing([10.0, 2.0])
                .show(ui, |ui| {
                    for h in ["Residue", "Buried Å²", "Contacts", "H-bonds"] {
                        ui.label(RichText::new(h).size(12.0).semibold());
                    }
                    ui.end_row();
                    for r in &done.rows {
                        if ui
                            .selectable_label(
                                false,
                                format!("{} {} {}", r.chain, r.name, r.number_text()),
                            )
                            .on_hover_text("Show it in the structure")
                            .clicked()
                        {
                            c.pick = Some(r.residue);
                        }
                        ui.label(RichText::new(format!("{:.1}", r.buried)).monospace());
                        ui.label(r.contacts.to_string());
                        ui.label(r.hbonds.to_string());
                        ui.end_row();
                    }
                });
        });
    ui.horizontal_wrapped(|ui| {
        c.colour = ui
            .button("Colour by buried area")
            .on_hover_text("Deep red where a residue buries much surface, pale where little")
            .clicked();
        c.copy = ui
            .button("Copy the table")
            .on_hover_text("As text for a spreadsheet")
            .clicked();
    });
    ui.label(RichText::new(EXPLANATION).size(11.5).color(look.dim));
}

/// What the numbers are, under the table: the counts above are of what this says.
const EXPLANATION: &str = "Measured on these coordinates, not predicted: solvent-accessible surface by the Shrake–Rupley method, \
    1.4 Å probe, hydrogens left out. A contact is a pair of residues, one from each side, with two atoms within 4.5 Å. An H-bond \
    candidate is a nitrogen or oxygen that can give and one that can take, within 3.5 Å: no angle is measured, so some are not \
    real bonds.";

impl InterfacePage {
    /// For a picture: measure the chains as they are set at once, and colour the structure by it.
    pub(super) fn measure_now(&mut self, loaded: &mut super::Loaded, style: &mut super::Style) {
        let (a, b) = (self.chains_on(Side::First), self.chains_on(Side::Second));
        match interface::analyse(&loaded.fold.structure, &a, &b) {
            Ok(interface) => {
                let rows = interface.rows(&loaded.fold.structure);
                self.message = format!("{} residues at the interface.", rows.len());
                loaded.fold.apply_interface(&interface);
                style.scheme = Scheme::Interface;
                self.measured = Some(Measured { interface, rows });
            }
            Err(why) => self.message = why,
        }
    }
}

impl App {
    /// Start a measurement, colour the structure by it, copy its table, or pick a residue of it.
    pub(super) fn apply_interface_clicks(&mut self, c: InterfaceClicks, ctx: &egui::Context) {
        let page = &mut self.platform.structure;
        let Some(loaded) = page.loaded.as_mut() else {
            return;
        };
        if c.changed {
            // what was measured, or is being measured, is not the answer to the new sides
            page.interface.take_back(&mut loaded.fold, &mut page.style);
            page.interface.message.clear();
        }
        if c.analyse {
            let (first, second) = (
                page.interface.chains_on(Side::First),
                page.interface.chains_on(Side::Second),
            );
            let structure = loaded.fold.structure.clone();
            page.interface.take_back(&mut loaded.fold, &mut page.style);
            page.interface.message = "Measuring…".into();
            page.interface.job = Some(Job::spawn(ctx, move || {
                interface::analyse(&structure, &first, &second)
            }));
        }
        if c.colour {
            if let Some(done) = &page.interface.measured {
                loaded.fold.apply_interface(&done.interface);
                page.style.scheme = Scheme::Interface;
            }
        }
        if c.copy {
            if let Some(done) = &page.interface.measured {
                ctx.copy_text(interface::tsv(&done.rows));
                page.interface.message = format!("Copied {} rows.", done.rows.len());
            }
        }
        if let Some(residue) = c.pick {
            // a ligand or a sugar is listed but not drawn as a chain: it cannot be picked, and the pick that was is kept
            match loaded.fold.residues.iter().position(|r| r.index == residue) {
                Some(k) => page.style.picked = Some(k),
                None => {
                    page.interface.message =
                        "That group is not part of the drawn protein chain, so it cannot be shown."
                            .into();
                }
            }
        }
    }

    /// Collect a finished measurement; once a frame.
    pub(super) fn interface_tick(&mut self, ctx: &egui::Context) {
        let page = &mut self.platform.structure;
        let Some(polled) = page.interface.job.as_ref().map(Job::poll) else {
            return;
        };
        match polled {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(done) => {
                page.interface.job = None;
                match (done, page.loaded.as_ref()) {
                    (Ok(interface), Some(loaded)) => {
                        let rows = interface.rows(&loaded.fold.structure);
                        page.interface.message =
                            format!("{} residues at the interface.", rows.len());
                        page.interface.measured = Some(Measured { interface, rows });
                    }
                    (Err(why), _) => page.interface.message = why,
                    (Ok(_), None) => {}
                }
            }
            Poll::Died => {
                page.interface.job = None;
                page.interface.message = "The measuring stopped without an answer.".into();
            }
        }
    }
}

#[cfg(test)]
mod tests;
