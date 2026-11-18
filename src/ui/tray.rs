//! The job tray: what is running, for how long, and what ended lately, with a way to go to each. The top bar shows how many analyses
//! are running; a click on it opens this panel. When an analysis ends there is a toast that says which, and for one that ran long a notice
//! from the operating system if the window is not the one the person is looking at.
//!
//! The analyses are those of the pages and windows (`App::running_sources`). They run on threads or wait for the workstation, and none can
//! be stopped half way, so the tray offers no stop button: it offers the way to the page (the Workstation page stops a job there).

use super::app::{Action, App};
use super::design::{self, Look, Tone};
use super::platform::Workspace;
use super::topbar;
use crate::ui::fonts::FontWeight;
use eframe::egui::{self, Align2, Frame, Id, Margin, Order, RichText, Stroke};
use std::time::Duration;
// the browser has no `std::time::Instant` (it panics there); `web_time` is the same on a computer and works in a page
use web_time::Instant;

/// An analysis that ran at least this long is told of by a notice outside the window when it ends in the background.
pub const LONG: Duration = Duration::from_secs(20);
/// How many ended analyses the tray keeps, and for how long.
const REMEMBERED: usize = 8;
const REMEMBERED_FOR: Duration = Duration::from_secs(15 * 60);
const WIDTH: f32 = 340.0;

/// An analysis that is running now.
#[derive(Clone, Debug)]
pub struct Source {
    pub key: &'static str,
    pub label: &'static str,
    /// What goes to its page.
    pub open: Action,
}

#[derive(Clone, Debug)]
pub struct Active {
    pub source: Source,
    pub started: Instant,
}

#[derive(Clone, Debug)]
pub struct Finished {
    pub source: Source,
    pub took: Duration,
    pub ended: Instant,
}

#[derive(Default)]
pub struct Tray {
    pub open: bool,
    active: Vec<Active>,
    history: Vec<Finished>,
}

impl Tray {
    /// Take note of what is running at `now`: a source that was not before has started, one that is not now has ended. Returns those that ended.
    pub fn observe(&mut self, running: &[Source], now: Instant) -> Vec<Finished> {
        for source in running {
            if !self.active.iter().any(|a| a.source.key == source.key) {
                self.active.push(Active {
                    source: source.clone(),
                    started: now,
                });
            }
        }
        let (still, ended): (Vec<Active>, Vec<Active>) = self
            .active
            .drain(..)
            .partition(|a| running.iter().any(|s| s.key == a.source.key));
        self.active = still;
        let finished: Vec<Finished> = ended
            .into_iter()
            .map(|a| Finished {
                took: now.saturating_duration_since(a.started),
                source: a.source,
                ended: now,
            })
            .collect();
        for f in &finished {
            self.history.insert(0, f.clone());
        }
        self.history.truncate(REMEMBERED);
        // what ended long ago is of no news any more (and the top bar's "Jobs" goes when nothing is left to show)
        self.history
            .retain(|f| now.saturating_duration_since(f.ended) < REMEMBERED_FOR);
        finished
    }

    pub fn running(&self) -> &[Active] {
        &self.active
    }

    pub fn history(&self) -> &[Finished] {
        &self.history
    }
}

/// A time span as people say it: "7 s", "2 min 5 s", "1 h 3 min".
pub fn took(span: Duration) -> String {
    let seconds = span.as_secs();
    match seconds {
        0 => "under a second".into(),
        1..=59 => format!("{seconds} s"),
        60..=3599 => match seconds % 60 {
            0 => format!("{} min", seconds / 60),
            s => format!("{} min {s} s", seconds / 60),
        },
        _ => match (seconds % 3600) / 60 {
            0 => format!("{} h", seconds / 3600),
            m => format!("{} h {m} min", seconds / 3600),
        },
    }
}

impl App {
    /// The analyses that are running now, with where each is shown.
    pub(super) fn running_sources(&self) -> Vec<Source> {
        let p = &self.platform;
        [
            (
                p.alignment.is_running(),
                "alignment",
                "Alignment",
                Action::Go(Workspace::Alignment),
            ),
            (
                p.phylogeny.is_running(),
                "phylogeny",
                "Phylogeny",
                Action::Go(Workspace::Phylogeny),
            ),
            (
                p.structure.is_running(),
                "structure",
                "Structure prediction",
                Action::Go(Workspace::Structure),
            ),
            (
                p.pipelines.is_running(),
                "pipelines",
                "Pipeline",
                Action::Go(Workspace::Pipelines),
            ),
            (
                self.phylocnn_win.is_running(),
                "phylocnn",
                "PhyloCNN estimate",
                Action::PhyloCnn,
            ),
            (
                self.species_win.is_running(),
                "species",
                "Species tree",
                Action::SpeciesTree,
            ),
            (
                self.motif_win.running(),
                "motifs",
                "Motif search",
                Action::Motifs,
            ),
            (
                self.domains_win.running(),
                "domains",
                "Domain search",
                Action::Domains,
            ),
        ]
        .into_iter()
        .filter(|(running, ..)| *running)
        .map(|(_, key, label, open)| Source { key, label, open })
        .collect()
    }

    /// Notice analyses that start and end: a toast for each that ends, and for one that ran long a notice outside the window if the
    /// window is not the one in front.
    pub(super) fn tray_step(&mut self, ctx: &egui::Context) {
        let sources = self.running_sources();
        let ended = self.tray.observe(&sources, Instant::now());
        let in_front = ctx.input(|i| i.viewport().focused.unwrap_or(true));
        for f in ended {
            // "ended", not "finished": the tray cannot tell an answer from a failure, and the page says which it was
            let said = format!("{} ended after {}.", f.source.label, took(f.took));
            self.toasts.push(said.clone(), Tone::Neutral);
            if f.took >= LONG && !in_front && self.keeper.has_settings() {
                super::notify::desktop("HelixLoom", &said);
            }
        }
        if !sources.is_empty() {
            // the clock in the panel goes on while something runs
            ctx.request_repaint_after(Duration::from_secs(1));
        }
    }

    /// The panel under the top bar.
    pub(super) fn tray_panel(&mut self, ctx: &egui::Context, out: &mut Vec<Action>) {
        if !self.tray.open {
            return;
        }
        let look = self.prefs.look();
        let now = Instant::now();
        // Escape closes the panel, unless the command palette is open and is the one that Escape is for
        let mut close = !self.palette.is_open()
            && ctx.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape));
        // the "ago" of the rows goes on while it is looked at
        ctx.request_repaint_after(Duration::from_secs(1));
        let area = egui::Area::new(Id::new("job tray"))
            .order(Order::Foreground)
            .anchor(Align2::RIGHT_TOP, egui::vec2(-12.0, topbar::HEIGHT + 4.0))
            .show(ctx, |ui| {
                Frame::new()
                    .fill(look.surface)
                    .stroke(Stroke::new(1.0, look.border))
                    .corner_radius(design::CARD_RADIUS)
                    .inner_margin(Margin::same(12))
                    .show(ui, |ui| {
                        ui.set_width(WIDTH);
                        self.tray_rows(ui, &look, now, out, &mut close);
                    });
            });
        // a click anywhere else closes it (the click that opened it was on the bar, which is not this area)
        let outside = ctx.input(|i| i.pointer.any_pressed())
            && !area.response.contains_pointer()
            && !ctx
                .input(|i| i.pointer.hover_pos())
                .is_some_and(|p| p.y < topbar::HEIGHT);
        if close || outside {
            self.tray.open = false;
        }
    }

    fn tray_rows(
        &self,
        ui: &mut egui::Ui,
        look: &Look,
        now: Instant,
        out: &mut Vec<Action>,
        close: &mut bool,
    ) {
        let (running, history) = (self.tray.running(), self.tray.history());
        if running.is_empty() && history.is_empty() {
            ui.label(RichText::new("Nothing is running. Analyses that you start show here, and for a while after they end.").color(look.dim));
            return;
        }
        let mut row =
            |ui: &mut egui::Ui, label: &str, detail: String, open: &Action, done: bool| {
                ui.horizontal(|ui| {
                    if done {
                        // the bundled font has no tick: a dot (not coloured "good": it may have failed)
                        let (rect, _) =
                            ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
                        ui.painter().circle_filled(rect.center(), 4.0, look.accent2);
                    } else {
                        ui.add(egui::Spinner::new().size(12.0).color(look.accent));
                    }
                    ui.label(RichText::new(label).semibold());
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui
                            .small_button(if done { "Open" } else { "Show" })
                            .clicked()
                        {
                            out.push(open.clone());
                            *close = true;
                        }
                        ui.label(RichText::new(detail).color(look.dim));
                    });
                });
            };
        if !running.is_empty() {
            ui.label(RichText::new("Running").color(look.dim));
            for a in running {
                row(
                    ui,
                    a.source.label,
                    took(now.saturating_duration_since(a.started)),
                    &a.source.open,
                    false,
                );
            }
        }
        if !history.is_empty() {
            if !running.is_empty() {
                ui.add_space(6.0);
            }
            ui.label(RichText::new("Ended").color(look.dim));
            for f in history {
                row(
                    ui,
                    f.source.label,
                    format!(
                        "{} · {} ago",
                        took(f.took),
                        took(now.saturating_duration_since(f.ended))
                    ),
                    &f.source.open,
                    true,
                );
            }
        }
        if self.platform.workstation.answering() {
            ui.add_space(6.0);
            if ui.link("Jobs on the workstation…").clicked() {
                out.push(Action::Go(Workspace::Workstation));
                *close = true;
            }
        }
    }
}

#[cfg(test)]
mod tests;
