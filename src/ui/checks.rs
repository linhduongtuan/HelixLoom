//! The page that says how each analysis was checked (`crate::checks`, `docs/ROADMAP.md` 2.4): Help → How Each Analysis Was
//! Checked. A count for each grade, a filter, and under each heading the analyses with their grade; a click on one shows what it was
//! compared with, how close it came, what was not checked and the tests that make the comparison.

use super::app::App;
use super::design::{self, Look, Tone};
use crate::checks::{self, Check, Standard};
use eframe::egui::{self, Context, RichText, Ui};
use std::collections::BTreeSet;

pub const TITLE: &str = "How each analysis was checked";
pub const COPY: &str = "Copy as text";
pub const NONE: &str = "No analysis holds these words.";
const INTRO: &str = "What each analysis was compared with, how close it came, and what nobody has compared yet. The tests named run on every change to the program.";
const GRADES: &str = "What the grades mean";
const HINT: &str = "Filter: a method, a program, a word";
const WIDTH: f32 = 640.0;
/// Tall enough that an analysis can be opened without the list scrolling away (egui's own default is 420).
const HEIGHT: f32 = 600.0;
const FILTER_WIDTH: f32 = 300.0;
/// The room a grade takes before the name of an analysis, so that the names stand in one column.
const GRADE_WIDTH: f32 = 150.0;

#[derive(Default)]
pub struct ChecksWin {
    pub open: bool,
    /// The words of the filter.
    pub query: String,
    /// The analyses whose details are shown.
    pub shown: BTreeSet<&'static str>,
}

/// The colour of a grade: green for a reference program, amber for what was not checked.
fn tone(standard: Standard) -> Tone {
    match standard {
        Standard::Reference => Tone::Good,
        Standard::Independent => Tone::Accent,
        Standard::Known => Tone::Neutral,
        Standard::Untested => Tone::Warn,
    }
}

impl App {
    pub(super) fn checks_window(&mut self, ctx: &Context) {
        if !self.checks_win.open {
            return;
        }
        let look = self.prefs.look();
        let win = &mut self.checks_win;
        let mut open = true;
        let (mut copy, mut copied) = (false, None);
        design::window(TITLE)
            .open(&mut open)
            .default_size([WIDTH, HEIGHT])
            .show(ctx, |ui| {
                ui.label(RichText::new(INTRO).color(look.dim));
                ui.add_space(4.0);
                ui.horizontal_wrapped(|ui| {
                    ui.add(
                        egui::TextEdit::singleline(&mut win.query)
                            .hint_text(HINT)
                            .desired_width(FILTER_WIDTH),
                    );
                    copy = ui.button(COPY).clicked();
                });
                let found = checks::matching(&win.query);
                if copy {
                    copied = Some(checks::to_text(&found));
                }
                ui.horizontal_wrapped(|ui| {
                    for (standard, n) in checks::tally(&found) {
                        design::chip(
                            ui,
                            &look,
                            &format!("{n} {}", standard.label()),
                            tone(standard),
                        );
                    }
                });
                egui::CollapsingHeader::new(GRADES).show(ui, |ui| {
                    for standard in Standard::ALL {
                        ui.horizontal_wrapped(|ui| {
                            design::chip(ui, &look, standard.label(), tone(standard));
                            ui.label(standard.meaning());
                        });
                    }
                });
                ui.separator();
                // the list takes what is left of the window, so that a taller window shows more of it
                egui::ScrollArea::vertical()
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        if found.is_empty() {
                            ui.label(RichText::new(NONE).color(look.dim));
                        }
                        for (group, list) in checks::grouped(&found) {
                            design::section(ui, &look, group);
                            for check in list {
                                row(ui, &look, check, &mut win.shown);
                            }
                        }
                    });
            });
        if let Some(text) = copied {
            ctx.copy_text(text);
            self.say("Copied how the analyses were checked.");
        }
        self.checks_win.open = open;
    }
}

/// One analysis: its grade and its name, which a click opens and closes, and under it what the table says.
fn row(ui: &mut Ui, look: &Look, check: &'static Check, shown: &mut BTreeSet<&'static str>) {
    let is_shown = shown.contains(check.analysis);
    ui.horizontal(|ui| {
        let left = ui.cursor().min.x;
        design::chip(ui, look, check.standard.label(), tone(check.standard));
        ui.add_space((left + GRADE_WIDTH - ui.cursor().min.x).max(0.0));
        if ui.selectable_label(is_shown, check.analysis).clicked() && !shown.remove(check.analysis)
        {
            shown.insert(check.analysis);
        }
    });
    if !is_shown {
        return;
    }
    ui.indent(("check", check.analysis), |ui| {
        for (heading, text) in [
            ("Compared with", check.against),
            ("Result", check.result),
            ("Not checked", check.limits),
        ] {
            ui.label(RichText::new(heading).small().color(look.dim));
            ui.label(text);
        }
        ui.label(RichText::new("Tests").small().color(look.dim));
        if check.tests.is_empty() {
            ui.label("none");
        }
        for (file, test) in check.tests {
            ui.label(RichText::new(*test).small());
            ui.label(RichText::new(*file).small().color(look.dim));
        }
        ui.add_space(4.0);
    });
}

#[cfg(test)]
mod tests;
