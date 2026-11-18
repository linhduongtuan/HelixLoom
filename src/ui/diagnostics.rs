//! Help → Copy Diagnostics (`crate::diagnostics`, `docs/ROADMAP.md` 1.8): what the program knows of itself and of the computer,
//! gathered for a report of a problem and put on the clipboard. Nothing of the person's work is in it, and nothing is sent.

use super::app::App;
use crate::diagnostics::{self, Facts};
use crate::pipeline::{install, tools};
use crate::ui::platform::workstation::Place;
use eframe::egui;

pub const COPIED: &str = "Copied the diagnostics: paste them into a report of a problem. They hold nothing of your work.";

/// Whether the settings folder is one the person named, not the usual one.
fn named_settings_folder() -> bool {
    std::env::var_os("HELIXLOOM_CONFIG_DIR").is_some_and(|dir| !dir.is_empty())
}

impl App {
    /// What a report of a problem needs, as text.
    pub(super) fn diagnostics_text(&self) -> String {
        let programs = self.platform.programs_at(&Place::Here).map(|found| {
            tools::ALL
                .iter()
                .map(|p| {
                    let here = found
                        .iter()
                        .find(|(label, _)| *label == p.label)
                        .and_then(|(_, found)| found.as_ref());
                    diagnostics::looked(
                        p.label,
                        here.map(|f| (f.path.as_path(), tools::short_version(&f.version))),
                    )
                })
                .collect()
        });
        let features = [
            ("gpu", cfg!(feature = "gpu")),
            ("cuda", cfg!(feature = "cuda")),
        ]
        .into_iter()
        .filter(|(_, built)| *built)
        .map(|(name, _)| name.to_string())
        .collect();
        let facts = Facts {
            version: env!("CARGO_PKG_VERSION").to_string(),
            os: std::env::consts::OS.to_string(),
            arch: std::env::consts::ARCH.to_string(),
            features,
            programs,
            managers: install::managers()
                .iter()
                .map(|found| found.manager.label().to_string())
                .collect(),
            models: [&crate::models::PHYLOCNN, &crate::models::PFAM]
                .iter()
                .map(|group| (group.folder.to_string(), crate::models::installed(group)))
                .collect(),
            // in words: the path of a folder that the person named can hold a name
            settings: diagnostics::settings_in_words(
                self.settings_dir.is_some(),
                named_settings_folder(),
            )
            .to_string(),
            workstation: self.platform.workstation.for_diagnostics(),
            appearance: format!(
                "{}, {} %",
                if self.prefs.dark { "dark" } else { "light" },
                self.prefs.appearance.scale
            ),
            last_crash: self
                .settings_dir
                .as_deref()
                .and_then(diagnostics::last_crash),
            // the home folder by the name each system gives it
            home: ["HOME", "USERPROFILE"]
                .iter()
                .find_map(|name| std::env::var(name).ok().filter(|home| home.len() > 1)),
        };
        diagnostics::text(&facts)
    }

    pub(super) fn copy_diagnostics(&mut self, ctx: &egui::Context) {
        ctx.copy_text(self.diagnostics_text());
        self.say(COPIED);
    }
}

#[cfg(test)]
mod tests;
