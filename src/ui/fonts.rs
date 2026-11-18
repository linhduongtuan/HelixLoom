//! The program's fonts: Inter for text and JetBrains Mono for sequences and code, both under the SIL Open Font Licence 1.1
//! (`assets/fonts`, `CREDITS.md`). Each is one variable file, so every weight comes from it: `FontWeight::semibold` asks for 600
//! where a heading or a label needs it. egui's own fonts stay behind them for the signs that these two do not have.

use eframe::egui::{self, FontData, FontDefinitions, FontFamily, RichText};
use std::sync::Arc;

const INTER: &[u8] = include_bytes!("../../assets/fonts/InterVariable.ttf");
const MONO: &[u8] = include_bytes!("../../assets/fonts/JetBrainsMono-wght.ttf");

/// The weight of text that stands out: a heading, a label, a name.
pub const SEMIBOLD: f32 = 600.0;

/// egui's defaults with Inter in front of the proportional fonts and JetBrains Mono in front of the monospace ones.
pub fn definitions() -> FontDefinitions {
    let mut fonts = FontDefinitions::default();
    for (name, bytes, family) in [
        ("Inter", INTER, FontFamily::Proportional),
        ("JetBrains Mono", MONO, FontFamily::Monospace),
    ] {
        fonts
            .font_data
            .insert(name.to_string(), Arc::new(FontData::from_static(bytes)));
        fonts
            .families
            .entry(family)
            .or_default()
            .insert(0, name.to_string());
    }
    fonts
}

/// Give a window its fonts; at the start of the program, before the first frame.
pub fn install(ctx: &egui::Context) {
    ctx.set_fonts(definitions());
}

/// Text with a real weight: the same font, drawn at another point of its weight axis.
pub trait FontWeight {
    /// Stands out: the stronger colour of `strong`, and the semibold weight.
    fn semibold(self) -> Self;
}

impl FontWeight for RichText {
    fn semibold(self) -> RichText {
        self.strong().variation(b"wght", SEMIBOLD)
    }
}

#[cfg(test)]
mod tests;
