//! Toasts: short notes that appear at the bottom right when something has happened (a file saved, text copied, a job
//! done) and fade away by themselves. `App::say` makes one of everything the program says in its status line, so every
//! message that was only a line of small text at the bottom is now seen, wherever the person is looking.

use super::app::App;
use super::design::{self, Look, Tone, CARD_RADIUS};
use eframe::egui::{
    self, vec2, Align2, Context, CornerRadius, Frame, Id, Margin, Order, RichText, Stroke,
};

/// How long a toast stays, how long its fading takes, how many are shown at once, and how wide one is.
const LIFETIME: f64 = 4.5;
const FADE: f64 = 0.6;
const MOST: usize = 4;
const WIDTH: f32 = 340.0;
/// How far from the corner of the window the stack of toasts sits (the status bar is under it).
const MARGIN: egui::Vec2 = vec2(-18.0, -46.0);

struct Toast {
    text: String,
    tone: Tone,
    /// When it first showed, in the seconds of the window's clock; none until it does.
    born: Option<f64>,
}

#[derive(Default)]
pub struct Toasts {
    items: Vec<Toast>,
    /// How many analyses were running at the last look.
    running: usize,
}

/// What a message says about how things went, from how it starts.
pub fn tone_of(text: &str) -> Tone {
    let lower = text.to_lowercase();
    let failed = ["could not", "couldn't", "cannot", "can't", "no "]
        .iter()
        .any(|w| lower.starts_with(w));
    if failed || lower.contains("failed") || lower.contains("error") {
        Tone::Bad
    } else if [
        "saved",
        "copied",
        "exported",
        "done",
        "found",
        "loaded",
        "installed",
    ]
    .iter()
    .any(|w| lower.starts_with(w))
    {
        Tone::Good
    } else {
        Tone::Neutral
    }
}

impl Toasts {
    pub fn push(&mut self, text: impl Into<String>, tone: Tone) {
        let text = text.into();
        if text.trim().is_empty() || self.items.last().is_some_and(|t| t.text == text) {
            return;
        }
        self.items.push(Toast {
            text,
            tone,
            born: None,
        });
        if self.items.len() > MOST {
            self.items.remove(0);
        }
    }

    /// Say so when fewer analyses are running than a moment ago: one has finished, wherever the person was looking.
    pub fn follow_jobs(&mut self, running: usize) {
        if running < self.running {
            self.push("An analysis has finished.", Tone::Good);
        }
        self.running = running;
    }

    /// Start the clock of the new toasts at `now` and drop the ones that have run out.
    fn tick(&mut self, now: f64) {
        for toast in &mut self.items {
            toast.born.get_or_insert(now);
        }
        self.items
            .retain(|t| t.born.is_none_or(|b| now - b < LIFETIME));
    }

    /// How opaque a toast is at `now`: whole while it lives, fading in the last moments.
    fn opacity(born: f64, now: f64) -> f32 {
        (((born + LIFETIME - now) / FADE).clamp(0.0, 1.0)) as f32
    }

    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }
}

impl App {
    /// Show the toasts in the corner of the window.
    pub(super) fn toasts_ui(&mut self, ctx: &Context) {
        let running = self.running_analyses();
        self.toasts.follow_jobs(running);
        let now = ctx.input(|i| i.time);
        self.toasts.tick(now);
        if self.toasts.is_empty() {
            return;
        }
        let look = design::look(self.prefs.dark);
        egui::Area::new(Id::new("toasts"))
            .order(Order::Tooltip)
            .anchor(Align2::RIGHT_BOTTOM, MARGIN)
            .interactable(false)
            .show(ctx, |ui| {
                ui.spacing_mut().item_spacing.y = 8.0;
                for toast in &self.toasts.items {
                    let alpha = Toasts::opacity(toast.born.unwrap_or(now), now);
                    ui.scope(|ui| {
                        ui.set_opacity(alpha);
                        show_toast(ui, &look, toast);
                    });
                }
            });
        // the clock must go on while a toast is up, so that it fades and goes
        ctx.request_repaint_after(std::time::Duration::from_millis(100));
    }
}

fn show_toast(ui: &mut egui::Ui, look: &Look, toast: &Toast) {
    let dot = match toast.tone {
        Tone::Good => look.good,
        Tone::Bad => look.bad,
        Tone::Accent => look.accent,
        Tone::Neutral => look.accent2,
    };
    Frame::new()
        .fill(look.surface)
        .stroke(Stroke::new(1.0, look.border))
        .corner_radius(CornerRadius::same(CARD_RADIUS - 2))
        .inner_margin(Margin::symmetric(12, 9))
        .shadow(ui.style().visuals.popup_shadow)
        .show(ui, |ui| {
            ui.set_max_width(WIDTH);
            ui.horizontal(|ui| {
                let (r, _) = ui.allocate_exact_size(vec2(8.0, 8.0), egui::Sense::hover());
                ui.painter().circle_filled(r.center(), 4.0, dot);
                ui.label(RichText::new(&toast.text).size(13.5).color(look.text));
            });
        });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_message_says_how_it_went_from_its_first_words() {
        assert_eq!(tone_of("Saved demo.gb."), Tone::Good);
        assert_eq!(tone_of("Copied 120 characters."), Tone::Good);
        assert_eq!(tone_of("Could not open x.gb: no such file"), Tone::Bad);
        assert_eq!(tone_of("Couldn't save /tmp/x.gb: denied"), Tone::Bad);
        assert_eq!(tone_of("The run failed"), Tone::Bad);
        assert_eq!(tone_of("Click a feature first."), Tone::Neutral);
    }

    #[test]
    fn a_toast_comes_once_lives_for_its_time_and_goes_and_the_same_message_comes_again_later() {
        let mut toasts = Toasts::default();
        toasts.push("Saved a.gb.", Tone::Good);
        toasts.push("Saved a.gb.", Tone::Good);
        assert_eq!(
            toasts.items.len(),
            1,
            "while it is up the same message is one toast"
        );
        toasts.tick(10.0);
        assert_eq!(toasts.items[0].born, Some(10.0));
        toasts.tick(10.0 + LIFETIME - 0.1);
        assert_eq!(toasts.items.len(), 1);
        toasts.tick(10.0 + LIFETIME + 0.1);
        assert!(toasts.is_empty());
        // an empty message is nothing to say
        toasts.push("  ", Tone::Neutral);
        assert!(toasts.is_empty());
        // copying twice is two toasts, a while apart
        toasts.push("Saved a.gb.", Tone::Good);
        assert_eq!(toasts.items.len(), 1);
    }

    #[test]
    fn a_toast_says_when_an_analysis_has_finished() {
        let mut toasts = Toasts::default();
        toasts.follow_jobs(0);
        toasts.follow_jobs(2);
        assert!(toasts.is_empty(), "starting is not news");
        toasts.follow_jobs(1);
        assert_eq!(toasts.items.len(), 1);
        assert_eq!(toasts.items[0].tone, Tone::Good);
    }

    #[test]
    fn only_the_newest_few_are_kept_and_they_fade_at_the_end() {
        let mut toasts = Toasts::default();
        for k in 0..MOST + 3 {
            toasts.push(format!("note {k}"), Tone::Neutral);
        }
        assert_eq!(toasts.items.len(), MOST);
        assert_eq!(toasts.items[0].text, "note 3");
        assert_eq!(Toasts::opacity(0.0, 1.0), 1.0);
        assert!(Toasts::opacity(0.0, LIFETIME - FADE / 2.0) < 1.0);
        assert_eq!(Toasts::opacity(0.0, LIFETIME + 1.0), 0.0);
    }
}
