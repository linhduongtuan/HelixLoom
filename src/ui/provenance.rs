//! The section under a result that says what made it (`crate::provenance`): the day, where it ran, the programs, the settings and the
//! inputs, then what has changed since, with buttons to copy the record and to run it again with the settings it was made with.

use super::design::Look;
use crate::provenance::{Input, Provenance};
use eframe::egui::{self, RichText, Ui};

/// The most inputs named in the hover text of the inputs line.
const MOST_NAMED: usize = 20;
pub const TITLE: &str = "How this was made";
pub const COPY: &str = "Copy the record";
pub const AGAIN: &str = "Run again with these settings";
pub const CHECKED: &str = "How this kind of analysis was checked";
pub const UNCHANGED: &str = "Nothing has changed since: running it again gives this result.";

/// What was pressed in the section.
#[derive(Default, Clone, Copy, PartialEq, Eq, Debug)]
pub struct Asks {
    pub copy: bool,
    pub again: bool,
    /// Show the page that says how each analysis was checked.
    pub checks: bool,
}

/// The section, closed until it is clicked. `now` gives the inputs and settings as they are, to say what has changed since; it is
/// asked only while the section is open. `can_run` is false while a run is going.
// ponytail: `now` makes the checksums of every input on each frame the section is open; cache them by a revision of the page if a
// big set makes the page slow
pub fn section(
    ui: &mut Ui,
    look: &Look,
    made: &Provenance,
    now: impl FnOnce() -> (Vec<Input>, Vec<(String, String)>),
    can_run: bool,
) -> Asks {
    let mut asks = Asks::default();
    egui::CollapsingHeader::new(TITLE)
        .id_salt(("provenance", made.analysis.as_str()))
        .show(ui, |ui| {
            // the buttons first: a long list of changes must not push them out of sight
            ui.horizontal_wrapped(|ui| {
                asks.copy = ui.button(COPY).clicked();
                asks.again = ui.add_enabled(can_run, egui::Button::new(AGAIN)).clicked();
            });
            ui.label(format!(
                "{} on {}, in {:.2} s",
                made.date, made.ran_on, made.seconds
            ));
            ui.label(format!("Programs: {}", made.programs.join("; ")));
            let settings: Vec<String> = made
                .parameters
                .iter()
                .map(|(k, v)| format!("{k}={v}"))
                .collect();
            ui.label(format!(
                "Settings: {}",
                if settings.is_empty() {
                    "none".to_string()
                } else {
                    settings.join("; ")
                }
            ));
            let named: Vec<String> = made
                .inputs
                .iter()
                .take(MOST_NAMED)
                .map(|i| format!("{}  ({} bytes, checksum {}…)", i.name, i.size, i.short()))
                .collect();
            ui.label(format!("Inputs: {}", made.inputs.len()))
                .on_hover_text(named.join("\n"));
            let (inputs, parameters) = now();
            let changes = made.changed_since(&inputs, &parameters);
            if changes.is_empty() {
                ui.label(RichText::new(UNCHANGED).color(look.good));
            } else {
                ui.label(RichText::new("Changed since:").color(look.warn));
                for change in &changes {
                    ui.label(format!("• {change}"));
                }
            }
            asks.checks = ui.link(CHECKED).clicked();
        });
    asks
}

#[cfg(test)]
mod tests;
