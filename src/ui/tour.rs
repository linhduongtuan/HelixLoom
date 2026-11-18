//! The first-run tour: three cards, each pointing at a part of the window that a newcomer would not find by looking: the
//! search box, the pages on the left, and where the heavy jobs go. It is shown once; a file in the settings folder says that
//! it was (in a browser, a mark in the storage the browser gives the page). Help > Show the Tour brings it back.

use super::app::{Action, App, Tab};
use super::design::{self, Look};
use super::platform::Workspace;
use super::session::Store;
use crate::ui::fonts::FontWeight;
use eframe::egui::{
    self, pos2, vec2, CornerRadius, Id, Key, LayerId, Modifiers, Order, Pos2, Rect, Stroke,
    StrokeKind, Vec2,
};
use std::path::Path;

/// The file in the settings folder that says the tour was seen.
const SEEN: &str = "tour-seen";
/// How wide a card is, and how tall it is taken to be when it is put near the edge of the window.
const CARD: Vec2 = vec2(310.0, 170.0);
const MARGIN: f32 = 8.0;
const GAP: f32 = 12.0;

/// What a step points at.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Target {
    Search,
    Rail,
    Jobs,
}

pub struct Step {
    pub target: Target,
    pub title: &'static str,
    pub text: &'static str,
}

pub const STEPS: [Step; 3] = [
    Step {
        target: Target::Search,
        title: "Search for anything",
        text: "Click here, or press the shortcut shown in the box, and type what you want to do: align, tree, primer, \
               save. Every tool and page is in the list.",
    },
    Step {
        target: Target::Rail,
        title: "The pages",
        text: "Sequences edits and maps DNA. Signals reads traces. Alignment, Phylogeny and Trees work on many sequences; \
               Structure folds a protein; Pipelines chains the steps; Workstation is your big computer.",
    },
    Step {
        target: Target::Jobs,
        title: "Where the jobs run",
        text: "Long analyses run in the background, on this computer or on your workstation, and the window stays \
               free. While one runs, a Jobs button appears here with the ones that ended lately.",
    },
];

/// Where the tour is, and where the things it points at are on the screen.
#[derive(Default)]
pub struct Tour {
    step: Option<usize>,
    /// Where it is kept that the tour was seen; none: it is not kept (and the tour is not shown by itself).
    store: Option<Store>,
    targets: [Option<Rect>; STEPS.len()],
}

impl Tour {
    /// The tour of a first run: on when there is a settings folder to remember it in, it holds no mark that the tour was seen
    /// and the program was not asked for something particular (a file, a picture).
    pub fn first_run(folder: Option<&Path>, asked_for_something: bool) -> Tour {
        // with no folder this is the browser's storage in a page, and nothing anywhere else
        Tour::with(Store::open(folder), asked_for_something)
    }

    /// The tour of a first run for what `store` keeps (none: nothing is kept, and the tour is not shown by itself).
    pub(super) fn with(store: Option<Store>, asked_for_something: bool) -> Tour {
        let unseen = store.as_ref().is_some_and(|s| !s.has_mark(SEEN));
        Tour {
            step: (unseen && !asked_for_something).then_some(0),
            store,
            targets: Default::default(),
        }
    }

    pub fn replay(&mut self) {
        self.step = Some(0);
    }

    pub fn current(&self) -> Option<&'static Step> {
        self.step.map(|k| &STEPS[k])
    }

    /// 1-based number of the step shown, and how many there are.
    pub fn number(&self) -> Option<(usize, usize)> {
        self.step.map(|k| (k + 1, STEPS.len()))
    }

    pub fn is_last(&self) -> bool {
        self.step == Some(STEPS.len() - 1)
    }

    /// Where the window drew `target` this frame.
    pub fn aim(&mut self, target: Target, rect: Rect) {
        self.targets[target as usize] = Some(rect);
    }

    pub fn rect(&self, target: Target) -> Option<Rect> {
        self.targets[target as usize]
    }

    /// The next card; after the last one the tour ends.
    pub fn advance(&mut self) -> Result<(), String> {
        match self.step {
            Some(k) if k + 1 < STEPS.len() => {
                self.step = Some(k + 1);
                Ok(())
            }
            Some(_) => self.finish(),
            None => Ok(()),
        }
    }

    /// End the tour, and keep that it was seen. It ends for this run even if the mark cannot be written.
    pub fn finish(&mut self) -> Result<(), String> {
        self.step = None;
        match &self.store {
            Some(store) => store.keep_mark(SEEN),
            None => Ok(()),
        }
    }
}

/// Where a card goes: below a target at the top of the window, aligned to its right edge; beside the pages; and always inside
/// the window.
pub fn card_position(target: Option<Rect>, which: Target, screen: Rect) -> Pos2 {
    let wanted = match (target, which) {
        (Some(t), Target::Rail) => pos2(t.right() + GAP, t.top() + 30.0),
        (Some(t), _) => pos2(t.right() - CARD.x, t.bottom() + GAP),
        (None, _) => screen.center() - CARD / 2.0,
    };
    // `min` then `max`, not `clamp`: in a window smaller than the card the left and top edges win, and nothing panics
    pos2(
        wanted
            .x
            .min(screen.right() - CARD.x - MARGIN)
            .max(screen.left() + MARGIN),
        wanted
            .y
            .min(screen.bottom() - CARD.y - MARGIN)
            .max(screen.top() + MARGIN),
    )
}

/// Help > Open the Sample Project: the demo plasmid in Sequences, the demo trace in Signals and the demo tree in Trees, so that
/// there is something on every page to look at.
impl App {
    pub(super) fn open_sample_project(&mut self) {
        match self
            .docs
            .iter()
            .position(|d| d.doc.name == crate::sample::DEMO_NAME)
        {
            Some(i) => self.tab = Tab::Doc(i),
            None => self.add_doc(crate::sample::demo_plasmid()),
        }
        self.platform
            .signals
            .load("example trace", crate::sample::demo_trace());
        self.show_tree(
            crate::sample::demo_tree(),
            "The example tree: 24 made-up strains in three groups.",
        );
        self.workspace = Workspace::Sequences;
        self.say("Sample project opened: the plasmid is here, the trace is in Signals and the tree in Trees.");
    }
}

enum Choice {
    Next,
    Skip,
    Sample,
}

impl App {
    /// Draw the card of the step the tour is at, with a ring round what it points at.
    pub(super) fn tour_window(&mut self, ctx: &egui::Context, out: &mut Vec<Action>) {
        // a saved session is offered first
        if self.keeper.offer.is_some() {
            return;
        }
        let (Some(step), Some((number, total))) = (self.tour.current(), self.tour.number()) else {
            return;
        };
        let look = self.prefs.look();
        let target = self.tour.rect(step.target);
        if let Some(rect) = target {
            ctx.layer_painter(LayerId::new(Order::Foreground, Id::new("tour ring")))
                .rect_stroke(
                    rect.expand(3.0),
                    CornerRadius::same(8),
                    Stroke::new(2.5, look.accent),
                    StrokeKind::Outside,
                );
        }
        let pos = card_position(target, step.target, ctx.content_rect());
        let last = self.tour.is_last();
        let choice = egui::Area::new(Id::new("tour card"))
            .order(Order::Foreground)
            .fixed_pos(pos)
            .show(ctx, |ui| {
                design::card(ui, &look, |ui| {
                    ui.set_width(CARD.x - 2.0 * design::CARD_MARGIN as f32);
                    tour_card(ui, &look, step, (number, total), last)
                })
            })
            .inner;
        // Esc ends the tour unless it is meant for something else: the palette, or a text field that has the keyboard
        let escape = !self.palette.is_open()
            && !ctx.egui_wants_keyboard_input()
            && ctx.input_mut(|i| i.consume_key(Modifiers::NONE, Key::Escape));
        let choice = if escape { Some(Choice::Skip) } else { choice };
        let result = match choice {
            Some(Choice::Next) => self.tour.advance(),
            Some(Choice::Skip) => self.tour.finish(),
            Some(Choice::Sample) => {
                out.push(Action::SampleProject);
                self.tour.finish()
            }
            None => Ok(()),
        };
        if let Err(why) = result {
            self.say(format!("Could not keep that the tour was seen: {why}"));
        }
    }
}

fn tour_card(
    ui: &mut egui::Ui,
    look: &Look,
    step: &Step,
    (number, total): (usize, usize),
    last: bool,
) -> Option<Choice> {
    design::section(ui, look, &format!("Step {number} of {total}"));
    ui.label(
        egui::RichText::new(step.title)
            .size(16.5)
            .semibold()
            .color(look.text),
    );
    ui.label(egui::RichText::new(step.text).color(look.dim));
    ui.add_space(8.0);
    let mut choice = None;
    ui.horizontal_wrapped(|ui| {
        if design::primary_button(ui, look, if last { "Done" } else { "Next" }).clicked() {
            choice = Some(if last { Choice::Skip } else { Choice::Next });
        }
        if last && design::ghost_button(ui, look, "Open the sample project").clicked() {
            choice = Some(Choice::Sample);
        }
        if !last && design::ghost_button(ui, look, "Skip the tour").clicked() {
            choice = Some(Choice::Skip);
        }
    });
    choice
}

#[cfg(test)]
mod tests;
