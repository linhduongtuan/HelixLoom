//! Colours and fonts. Values come from the original's defaults (ApE_Defaults.txt and
//! the colour settings in AppMain.tcl); light mode is the default, dark mode follows its dark set.

use crate::seq::features::Rgb;
use eframe::egui::{self, Color32, FontFamily, FontId, TextStyle};

#[derive(Clone, Copy)]
pub struct Palette {
    pub dark: bool,
    pub window_bg: Color32,
    pub text_bg: Color32,
    pub fg: Color32,
    pub dim: Color32,
    pub select_bg: Color32,
    pub select_fg: Color32,
    pub find_fwd: Color32,
    pub find_rev: Color32,
    pub enzyme: Color32,
    pub line: Color32,
    pub accent: Color32,
}

pub fn palette(dark: bool) -> Palette {
    if dark {
        Palette {
            dark,
            // the platform redesign replaced the original's #334 behind the editor
            window_bg: super::design::look(true).bg,
            text_bg: Color32::from_rgb(0x1a, 0x1d, 0x21),
            fg: Color32::WHITE,
            dim: Color32::from_gray(150),
            select_bg: Color32::from_rgb(0x86, 0x86, 0xe8),
            select_fg: Color32::BLACK,
            find_fwd: Color32::from_gray(110),
            find_rev: Color32::from_gray(80),
            enzyme: Color32::from_rgb(255, 90, 90),
            line: Color32::from_gray(90),
            accent: Color32::from_rgb(90, 150, 255),
        }
    } else {
        Palette {
            dark,
            window_bg: super::design::look(false).bg,
            text_bg: Color32::WHITE,
            fg: Color32::BLACK,
            dim: Color32::from_gray(110),
            // text_select2_bg_color / text_select2_fg_color
            select_bg: Color32::from_rgb(0x60, 0x60, 0xd0),
            select_fg: Color32::WHITE,
            // find_forward_highlight gray75, find_reverse_highlight gray50
            find_fwd: Color32::from_gray(191),
            find_rev: Color32::from_gray(127),
            // enzyme_highlight_color
            enzyme: Color32::from_rgb(255, 0, 0),
            line: Color32::from_gray(190),
            accent: Color32::from_rgb(20, 100, 240),
        }
    }
}

pub fn rgb(c: Rgb) -> Color32 {
    Color32::from_rgb(c[0], c[1], c[2])
}

/// Black or white, whichever reads better on `bg`.
pub fn text_on(bg: Color32) -> Color32 {
    let lum = 0.299 * bg.r() as f32 + 0.587 * bg.g() as f32 + 0.114 * bg.b() as f32;
    if lum > 140.0 {
        Color32::BLACK
    } else {
        Color32::WHITE
    }
}

/// Blend `a` toward `b` by `pct` percent of `a` (the original's mix_color).
pub fn mix(a: Color32, b: Color32, pct: f32) -> Color32 {
    let t = pct / 100.0;
    let f = |x: u8, y: u8| (x as f32 * t + y as f32 * (1.0 - t)).round() as u8;
    Color32::from_rgb(f(a.r(), b.r()), f(a.g(), b.g()), f(a.b(), b.b()))
}

/// The visuals of one of the two themes: the platform's colours on egui's own light or dark base.
fn visuals_for(dark: bool) -> egui::Visuals {
    let pal = palette(dark);
    let mut visuals = if dark {
        egui::Visuals::dark()
    } else {
        egui::Visuals::light()
    };
    visuals.panel_fill = pal.window_bg;
    visuals.window_fill = pal.window_bg;
    visuals.extreme_bg_color = pal.text_bg;
    visuals.selection.bg_fill = pal.accent.gamma_multiply(0.6);
    visuals.window_corner_radius = 6.0.into();
    super::design::decorate(&mut visuals, &super::design::look(dark));
    visuals
}

/// Put the colours on the window. Both themes are set, and when `follow_system` the window is left to follow the system's own
/// setting: pinning it to light or dark would also pin the window's title bar, and on a Mac the system would then stop
/// telling the program when it changes. Where the system says nothing, the light colours are used.
pub fn apply(ctx: &egui::Context, dark: bool, follow_system: bool) {
    ctx.set_visuals_of(egui::Theme::Dark, visuals_for(true));
    ctx.set_visuals_of(egui::Theme::Light, visuals_for(false));
    ctx.options_mut(|options| options.fallback_theme = egui::Theme::Light);
    ctx.set_theme(match (follow_system, dark) {
        (true, _) => egui::ThemePreference::System,
        (false, true) => egui::ThemePreference::Dark,
        (false, false) => egui::ThemePreference::Light,
    });
    ctx.all_styles_mut(|s| {
        s.text_styles = [
            (
                TextStyle::Heading,
                FontId::new(17.0, FontFamily::Proportional),
            ),
            (TextStyle::Body, FontId::new(13.5, FontFamily::Proportional)),
            (
                TextStyle::Button,
                FontId::new(13.5, FontFamily::Proportional),
            ),
            (
                TextStyle::Small,
                FontId::new(11.0, FontFamily::Proportional),
            ),
            (
                TextStyle::Monospace,
                FontId::new(13.0, FontFamily::Monospace),
            ),
        ]
        .into();
    });
    super::design::space(ctx);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn text_colour_follows_background_brightness() {
        assert_eq!(text_on(Color32::WHITE), Color32::BLACK);
        assert_eq!(text_on(Color32::from_rgb(0x34, 0x6e, 0xe0)), Color32::WHITE);
        assert_eq!(text_on(Color32::from_rgb(0xe9, 0xd0, 0x24)), Color32::BLACK);
    }

    #[test]
    fn mixing_matches_the_original_formula() {
        let m = mix(Color32::BLACK, Color32::WHITE, 25.0);
        assert_eq!((m.r(), m.g(), m.b()), (191, 191, 191));
    }
}
