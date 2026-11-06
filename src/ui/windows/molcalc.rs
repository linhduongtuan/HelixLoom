//! Tools → Molecular Reaction Calculator: the volumes of a molar ratio reaction (Gibson, HiFi, Golden Gate, ligation), a
//! master mix multiplier, and conversions between mass, moles, concentration and volume with a stock dilution. The three
//! tabs are those of the original's `mol_calc` (AppMain.tcl:34255, 0x18dde0); the arithmetic is in `seq/molcalc`.

mod mix;
mod molar;
mod molarity;

use crate::seq::molcalc::format::number;
use crate::ui::app::App;
use eframe::egui::{self, Ui};

#[derive(Clone, Copy, Default, PartialEq, Eq)]
enum Tab {
    #[default]
    Molar,
    Mix,
    Molarity,
}

impl Tab {
    const ALL: [Tab; 3] = [Tab::Molar, Tab::Mix, Tab::Molarity];

    fn label(self) -> &'static str {
        match self {
            Tab::Molar => "Molar Ratio Reaction",
            Tab::Mix => "Master Mix Multiplier",
            Tab::Molarity => "Molarity",
        }
    }
}

/// What the window keeps between frames (and while it is closed: the original starts again from the last values).
#[derive(Default)]
pub struct MolCalcWin {
    pub open: bool,
    tab: Tab,
    molar: molar::MolarTab,
    mix: mix::MixTab,
    molarity: molarity::MolarityTab,
}

impl MolCalcWin {
    /// Fill the first tab with two fragments, for `--view molcalc`.
    pub fn load_example(&mut self) {
        self.molar.example();
    }

    /// Show a tab by the name `--view molcalc-NAME` gives it (the first one for any other).
    pub fn show_tab(&mut self, name: &str) {
        self.tab = match name {
            "mix" => Tab::Mix,
            "molarity" => Tab::Molarity,
            _ => Tab::Molar,
        };
    }

    fn show(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            for tab in Tab::ALL {
                ui.selectable_value(&mut self.tab, tab, tab.label());
            }
        });
        ui.separator();
        match self.tab {
            Tab::Molar => self.molar.show(ui),
            Tab::Mix => self.mix.show(ui),
            Tab::Molarity => self.molarity.show(ui),
        }
    }
}

impl App {
    pub fn molcalc_window(&mut self, ctx: &egui::Context) {
        if !self.molcalc_win.open {
            return;
        }
        let mut open = true;
        crate::ui::design::window("Molecular Reaction Calculator")
            .open(&mut open)
            .default_width(640.0)
            .show(ctx, |ui| self.molcalc_win.show(ui));
        self.molcalc_win.open = open;
    }
}

/// What the original's spin boxes let into a field (`spinbox_validate`): nothing, a point, or a number that is not negative.
fn acceptable(text: &str) -> bool {
    text.is_empty() || text == "." || number(text).is_some_and(|v| v >= 0.0)
}

/// A field for a number. What is typed that the original would not take is put back; true when the text changed.
fn number_field(ui: &mut Ui, text: &mut String, width: f32) -> bool {
    let before = text.clone();
    ui.add(egui::TextEdit::singleline(text).desired_width(width));
    if *text == before {
        return false;
    }
    if acceptable(text) {
        true
    } else {
        *text = before;
        false
    }
}

/// A drop-down list of `options`; true when another one was chosen.
fn combo<T: Copy + PartialEq>(
    ui: &mut Ui,
    id: &str,
    value: &mut T,
    options: &[T],
    label: impl Fn(T) -> &'static str,
) -> bool {
    let before = *value;
    egui::ComboBox::from_id_salt(id)
        .selected_text(label(before))
        .show_ui(ui, |ui| {
            for &option in options {
                ui.selectable_value(value, option, label(option));
            }
        });
    *value != before
}

/// A bold result, or a dash while there is none.
fn result(ui: &mut Ui, text: &str) {
    ui.label(egui::RichText::new(if text.is_empty() { "-" } else { text }).strong());
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ui::app::Startup;
    use crate::ui::headless::frame_input;

    fn frame(app: &mut App, ctx: &egui::Context, n: u32) {
        let mut output = ctx.run_ui(frame_input(n), |ui| app.draw(ui));
        output.textures_delta.clear();
    }

    #[test]
    fn a_field_takes_what_the_spin_boxes_of_the_original_take() {
        for ok in ["", ".", "0", "1.5", "5.", ".5", "1e3", " 2 "] {
            assert!(acceptable(ok), "{ok:?}");
        }
        for not in ["-1", "abc", "1,5", "1e", "inf", "--", "1.2.3"] {
            assert!(!acceptable(not), "{not:?}");
        }
    }

    #[test]
    fn the_window_draws_every_tab_with_the_first_values_and_with_an_example_in_both_themes() {
        let ctx = egui::Context::default();
        for dark in [false, true] {
            for example in [false, true] {
                let mut app = App::new(Startup {
                    dark,
                    ..Startup::default()
                });
                app.molcalc_win.open = true;
                if example {
                    app.molcalc_win.load_example();
                }
                for (n, tab) in ["", "mix", "molarity"].into_iter().enumerate() {
                    app.molcalc_win.show_tab(tab);
                    frame(&mut app, &ctx, n as u32);
                }
            }
        }
    }
}
