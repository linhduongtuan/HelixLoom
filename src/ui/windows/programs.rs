//! Tools → Programs HelixLoom Calls: the tool registry in a window (`docs/ROADMAP.md` 2.1). Each program the platform can call,
//! what it is for, whether it is found here and in which version, and its licence; one that is not found can be installed, by
//! the package manager that is on the computer, after the command has been shown and the person has said yes
//! (`pipeline/install.rs`). Nothing is installed without that yes, and nothing here is about the workstation, whose owner
//! installs what it has.

use crate::pipeline::install::{self, Found, Plan};
use crate::pipeline::tools::{self, Installed, Program};
use crate::ui::app::App;
use crate::ui::design::{self, Look, Tone};
use crate::ui::dialog;
use crate::ui::job::{Job, Poll};
use crate::ui::platform::workstation::Place;
use eframe::egui::{self, RichText, ScrollArea};
use std::sync::{Arc, Mutex};

pub const TITLE: &str = "Programs HelixLoom Calls";
pub const INSTALL: &str = "Install…";
pub const RUN: &str = "Run it";
pub const COPY: &str = "Copy the command";
pub const CANCEL: &str = "Cancel";
pub const AGAIN: &str = "Look again";
pub const LOOKING: &str = "Looking for installed programs…";
const INTRO: &str = "The programs of other projects that HelixLoom starts for some of its work. It looks for them on this computer and never installs one unless you ask it to here.";
const IN_BROWSER: &str = "A page in a browser starts no program: the programs are those of the workstation it is connected to (Workstation page).";
const WIDTH: f32 = 660.0;
const HEIGHT: f32 = 560.0;
const LOG_HEIGHT: f32 = 130.0;
/// The most lines of a package manager's output that are kept for the window.
const MOST_LINES: usize = 2_000;

/// Runs a plan, each line the package manager prints going to `say`.
type Installer = fn(&Plan, &mut dyn FnMut(&str)) -> Result<(), String>;

fn run_plan(plan: &Plan, say: &mut dyn FnMut(&str)) -> Result<(), String> {
    install::run(plan, say)
}

pub struct ProgramsWin {
    pub open: bool,
    /// The plan that waits for the person's yes.
    asked: Option<Plan>,
    /// The plan that is running, and its work.
    running: Option<(Plan, Job<Result<(), String>>)>,
    /// What the package manager printed.
    log: Arc<Mutex<Vec<String>>>,
    pub message: String,
    /// Finds the package managers, and runs a plan: the tests put stand-ins here.
    managers: fn() -> Vec<Found>,
    installer: Installer,
}

impl Default for ProgramsWin {
    fn default() -> Self {
        ProgramsWin {
            open: false,
            asked: None,
            running: None,
            log: Arc::default(),
            message: String::new(),
            managers: install::managers,
            installer: run_plan,
        }
    }
}

impl ProgramsWin {
    pub fn is_installing(&self) -> bool {
        self.running.is_some()
    }

    /// Make the plan to install `program` and show it, or say why there is none.
    fn ask(&mut self, program: &Program) {
        let active = std::env::var("CONDA_PREFIX").ok();
        match install::plan(program, &(self.managers)(), None, active.as_deref()) {
            Ok(plan) => (self.asked, self.message) = (Some(plan), String::new()),
            Err(why) => (self.asked, self.message) = (None, why),
        }
    }

    /// Start the plan that was shown.
    fn start(&mut self, ctx: &egui::Context) {
        let Some(plan) = self.asked.take() else {
            return;
        };
        if let Ok(mut log) = self.log.lock() {
            log.clear();
        }
        self.message = format!("Running {}…", plan.line());
        let (installer, log, work) = (self.installer, Arc::clone(&self.log), plan.clone());
        let job = Job::spawn(ctx, move || {
            installer(&work, &mut |line: &str| {
                if let Ok(mut log) = log.lock() {
                    if log.len() < MOST_LINES {
                        log.push(line.to_string());
                    }
                }
            })
        });
        self.running = Some((plan, job));
    }

    /// Whether an install has just ended well, so that the programs must be looked for again.
    fn poll(&mut self, ctx: &egui::Context) -> bool {
        let Some((plan, job)) = &self.running else {
            return false;
        };
        let (line, program) = (plan.line(), plan.program);
        match job.poll() {
            Poll::Running => {
                ctx.request_repaint_after(std::time::Duration::from_millis(150));
                false
            }
            Poll::Done(Ok(())) => {
                self.running = None;
                self.message = format!("{line} ended well: looking for {program} again.");
                true
            }
            Poll::Done(Err(why)) => {
                self.running = None;
                self.message = why;
                false
            }
            Poll::Died => {
                self.running = None;
                self.message = format!("{line} stopped without an answer.");
                false
            }
        }
    }
}

/// What the window says of one program that was found.
fn found_line(found: &Installed) -> String {
    let version = tools::short_version(&found.version);
    if version.is_empty() {
        found.path.display().to_string()
    } else {
        format!("{}  (version {version})", found.path.display())
    }
}

/// One program: its name, what it is for, whether it is here, its licence. True when its Install button was pressed.
fn row(
    ui: &mut egui::Ui,
    look: &Look,
    program: &Program,
    found: Option<&Installed>,
    can_install: bool,
) -> bool {
    let mut install = false;
    design::card(ui, look, |ui| {
        // every program's card is as wide as the list, whatever it says
        ui.set_min_width(ui.available_width());
        ui.horizontal(|ui| {
            ui.label(RichText::new(program.label).strong());
            match found {
                Some(_) => design::chip(ui, look, "found", Tone::Good),
                None => design::chip(ui, look, "not found", Tone::Warn),
            }
            if found.is_none() {
                install = ui
                    .add_enabled(can_install, egui::Button::new(INSTALL))
                    .clicked();
            }
        });
        ui.label(RichText::new(program.role).color(look.dim));
        if let Some(found) = found {
            ui.label(RichText::new(found_line(found)).monospace().small());
        } else {
            let ways = tools::install_ways(program);
            if !ways.is_empty() {
                ui.label(
                    RichText::new(format!("By hand: {}", ways.join("   or   ")))
                        .small()
                        .color(look.dim),
                );
            }
        }
        ui.label(
            RichText::new(format!("Licence: {}", program.licence))
                .small()
                .color(look.dim),
        );
    });
    install
}

impl App {
    pub fn programs_window(&mut self, ctx: &egui::Context) {
        if self.programs_win.poll(ctx) {
            self.platform.look_again();
        }
        if !self.programs_win.open {
            return;
        }
        let look = self.prefs.look();
        let found = self.platform.programs_at(&Place::Here);
        let mut open = true;
        let (mut again, mut copied) = (false, None);
        design::window(TITLE)
            .open(&mut open)
            .default_size([WIDTH, HEIGHT])
            .show(ctx, |ui| {
                let w = &mut self.programs_win;
                ui.label(RichText::new(INTRO).color(look.dim));
                if dialog::IN_BROWSER {
                    ui.label(IN_BROWSER);
                    return;
                }
                ui.horizontal(|ui| {
                    again = ui
                        .add_enabled(
                            found.is_some() && !w.is_installing(),
                            egui::Button::new(AGAIN),
                        )
                        .clicked();
                    if w.is_installing() {
                        ui.spinner();
                    }
                });
                if !w.message.is_empty() {
                    ui.label(&w.message);
                }
                if let Some(plan) = w.asked.clone() {
                    // nothing runs before this is read and answered
                    let program = tools::ALL.iter().find(|p| p.label == plan.program);
                    design::card(ui, &look, |ui| {
                        if let Some(program) = program {
                            ui.label(plan.notice(program));
                        }
                        ui.horizontal(|ui| {
                            if design::primary_button(ui, &look, RUN).clicked() {
                                w.start(ui.ctx());
                            }
                            if ui.button(COPY).clicked() {
                                copied = Some(plan.line());
                            }
                            if ui.button(CANCEL).clicked() {
                                w.asked = None;
                            }
                        });
                    });
                }
                let lines: Vec<String> = w.log.lock().map(|log| log.clone()).unwrap_or_default();
                if !lines.is_empty() {
                    ScrollArea::vertical()
                        .id_salt("programs_log")
                        .max_height(LOG_HEIGHT)
                        .stick_to_bottom(true)
                        .show(ui, |ui| {
                            for line in &lines {
                                ui.label(RichText::new(line).monospace().small());
                            }
                        });
                }
                ui.separator();
                let Some(found) = &found else {
                    design::chip(ui, &look, LOOKING, Tone::Neutral);
                    return;
                };
                let can_install = !w.is_installing() && w.asked.is_none();
                ScrollArea::vertical()
                    .id_salt("programs_list")
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        for program in tools::ALL {
                            let here = found
                                .iter()
                                .find(|(label, _)| *label == program.label)
                                .and_then(|(_, found)| found.as_ref());
                            if row(ui, &look, program, here, can_install) {
                                w.ask(program);
                            }
                            ui.add_space(4.0);
                        }
                    });
            });
        if again {
            self.programs_win.message.clear();
            self.platform.look_again();
        }
        if let Some(line) = copied {
            ctx.copy_text(line);
            self.say("Copied the command.");
        }
        self.programs_win.open = open;
    }
}

#[cfg(test)]
mod tests;
