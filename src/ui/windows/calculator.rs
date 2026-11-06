//! Tools → Calculator: a line of arithmetic and its result, as in the original's `calc` window (AppMain.tcl:34184). A
//! whole number is made real by a point after it, so `1/3` is 0.3333333333333333 and not 0; the rules are in
//! `seq/molcalc/expr.rs`.

use crate::seq::molcalc::expr::calculate;
use crate::ui::app::App;
use eframe::egui::{self, RichText};

#[derive(Default)]
pub struct CalculatorWin {
    pub open: bool,
    pub text: String,
    /// The window was just opened: the line gets the keyboard.
    pub fresh: bool,
    /// Something was typed: from then on an empty line is a mistake, as in the original.
    typed: bool,
}

/// What the result says: nothing before anything is typed (the original works nothing out until then), and after that what
/// the line comes to, which for an empty line is the mistake line.
fn shown(text: &str, typed: bool) -> String {
    if text.is_empty() && !typed {
        String::new()
    } else {
        calculate(text)
    }
}

impl App {
    pub fn calculator_window(&mut self, ctx: &egui::Context) {
        if !self.calculator_win.open {
            return;
        }
        let mut open = true;
        crate::ui::design::window("Calculator: Enter a calculation")
            .open(&mut open)
            .default_width(420.0)
            // the original can be made wider, not taller
            .resizable([true, false])
            .show(ctx, |ui| {
                let win = &mut self.calculator_win;
                let line = ui.add(
                    egui::TextEdit::singleline(&mut win.text)
                        .hint_text("for example 2 * (3 + 4) / 7, sqrt(2), $pi * 2 ** 3")
                        .desired_width(f32::INFINITY),
                );
                if std::mem::take(&mut win.fresh) {
                    line.request_focus();
                }
                win.typed |= line.changed();
                let result = shown(&win.text, win.typed);
                ui.horizontal(|ui| {
                    ui.label(RichText::new(&result).size(18.0).strong());
                    if !result.is_empty() && ui.small_button("Copy").clicked() {
                        ui.ctx().copy_text(result.clone());
                    }
                });
            });
        self.calculator_win.open = open;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ui::app::Startup;
    use crate::ui::headless::frame_input;

    #[test]
    fn an_empty_line_is_blank_until_something_is_typed_and_a_mistake_after() {
        assert_eq!(shown("", false), "");
        assert_eq!(shown("", true), "------");
        assert_eq!(shown("1/3", false), "0.3333333333333333");
    }

    #[test]
    fn the_window_shows_the_result_of_what_is_typed_and_nothing_before() {
        let ctx = egui::Context::default();
        let mut app = App::new(Startup::default());
        app.calculator_win.open = true;
        for (n, text) in ["", "1/3", "2 +", "0x10 + $pi"].into_iter().enumerate() {
            app.calculator_win.text = text.into();
            let mut output = ctx.run_ui(frame_input(n as u32), |ui| app.draw(ui));
            output.textures_delta.clear();
        }
    }
}
