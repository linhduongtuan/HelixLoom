//! What a person can set to suit their eyes, their screen and their habits: the size of everything, how roomy it is, the colour of the
//! program, a high-contrast theme, and colours for bases and residues that stay apart for people who see colour differently.
//!
//! The settings are kept with the theme (light, dark or the system's) in `appearance.txt` of the settings folder, a line to a setting
//! (`scale 125`), so that a line that is not understood spoils nothing else. `theme::apply` puts them on the window.

pub mod safe;
mod window;

use super::design::{self, contrast, Look};
use super::state::ThemeMode;
use super::theme::{mix, Palette};
use eframe::egui::{vec2, Color32, Vec2};

#[cfg(test)]
mod screen_tests;
#[cfg(test)]
mod tests;

/// The sizes the whole window can be drawn at, in per cent. 100 is as it was designed.
pub const SCALES: [u16; 8] = [80, 90, 100, 110, 125, 150, 175, 200];
pub const NORMAL_SCALE: u16 = 100;

/// The shortest a button or a row that can be clicked may be, in points: the least that the Web Content Accessibility Guidelines
/// (2.2, success criterion 2.5.8) allow.
pub const MIN_TARGET: f32 = 24.0;

/// How thick outlines are with the high-contrast setting (1 otherwise).
const HIGH_OUTLINE: f32 = 1.5;
/// What the text and the lines of the high-contrast theme have to reach against what is behind them (WCAG level AAA for text, and
/// the 3 that edges of boxes need).
const HIGH_TEXT: f32 = 7.0;
const HIGH_EDGE: f32 = 3.0;

/// A setting that is one of a few named choices: the first is the one that nothing has been changed from.
macro_rules! choices {
    ($(#[$meta:meta])* $name:ident { $first:ident => $first_key:literal, $first_label:literal; $($variant:ident => $key:literal, $label:literal;)* }) => {
        $(#[$meta])*
        #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
        pub enum $name {
            #[default]
            $first,
            $($variant),*
        }

        impl $name {
            pub const ALL: &'static [$name] = &[$name::$first, $($name::$variant),*];

            /// The word it is written as in `appearance.txt`.
            pub fn key(self) -> &'static str {
                match self {
                    $name::$first => $first_key,
                    $($name::$variant => $key),*
                }
            }

            pub fn label(self) -> &'static str {
                match self {
                    $name::$first => $first_label,
                    $($name::$variant => $label),*
                }
            }

            fn from_key(key: &str) -> Option<$name> {
                Self::ALL.iter().copied().find(|choice| choice.key() == key)
            }
        }
    };
}

choices! {
    /// How much room there is between things.
    Density {
        Comfortable => "comfortable", "Comfortable";
        Compact => "compact", "Compact";
    }
}

choices! {
    /// The colour of the program: of the main buttons, of what is chosen and of links.
    Accent {
        Indigo => "indigo", "Indigo";
        Blue => "blue", "Blue";
        Teal => "teal", "Teal";
        Green => "green", "Green";
        Purple => "purple", "Purple";
        Graphite => "graphite", "Graphite";
    }
}

choices! {
    /// Whether text, lines and colours are pushed apart from what is behind them.
    Contrast {
        Standard => "standard", "Standard";
        High => "high", "High contrast";
    }
}

choices! {
    /// The colours of bases, residues and trace letters.
    Colours {
        Standard => "standard", "Standard";
        Safe => "safe", "Colour-blind safe";
    }
}

/// The spaces that a density sets (see `design::space`).
pub struct Spacing {
    pub item: Vec2,
    pub button_padding: Vec2,
    pub interact_height: f32,
    pub window_margin: i8,
    pub menu_margin: i8,
}

impl Density {
    pub fn spacing(self) -> Spacing {
        match self {
            Density::Comfortable => Spacing {
                item: vec2(8.0, 6.0),
                button_padding: vec2(10.0, 4.0),
                interact_height: 26.0,
                window_margin: 12,
                menu_margin: 6,
            },
            Density::Compact => Spacing {
                item: vec2(6.0, 3.0),
                button_padding: vec2(8.0, 2.0),
                interact_height: MIN_TARGET,
                window_margin: 9,
                menu_margin: 4,
            },
        }
    }
}

impl Accent {
    /// The colour on a page of the light or dark theme, none for the indigo that the program was designed in.
    fn colour(self, dark: bool) -> Option<Color32> {
        let (light, deep) = match self {
            Accent::Indigo => return None,
            Accent::Blue => ((0x1a, 0x5f, 0xd1), (0x6b, 0xa4, 0xff)),
            Accent::Teal => ((0x0b, 0x7a, 0x73), (0x3d, 0xd1, 0xc5)),
            Accent::Green => ((0x1a, 0x7a, 0x3a), (0x4c, 0xc3, 0x8a)),
            Accent::Purple => ((0x80, 0x35, 0xc4), (0xc5, 0x97, 0xff)),
            Accent::Graphite => ((0x45, 0x4e, 0x66), (0xae, 0xb8, 0xcc)),
        };
        let (r, g, b) = if dark { deep } else { light };
        Some(Color32::from_rgb(r, g, b))
    }

    /// The colour as the chrome draws it on the page of the light or dark theme.
    pub fn swatch(self, dark: bool) -> Color32 {
        self.colour(dark)
            .unwrap_or_else(|| design::look(dark).accent)
    }

    fn tint(self, l: &mut Look) {
        let Some(accent) = self.colour(l.dark) else {
            return;
        };
        if self == Accent::Teal {
            // the second colour is a teal too: it takes the indigo's place
            l.accent2 = design::look(l.dark).accent;
        }
        l.accent = accent;
        l.accent_soft = mix(accent, l.surface, if l.dark { 22.0 } else { 10.0 });
        // the accent shown through at a third (the designed indigo's way) is too pale a ground for the text of another accent
        l.selection_fill = l.accent_soft;
        l.on_accent = ink_for(accent);
    }
}

/// Black or white, whichever reads better on the accent (the program's own near-black on a dark page).
fn ink_for(accent: Color32) -> Color32 {
    let dark_ink = Color32::from_rgb(0x0b, 0x0e, 0x15);
    if contrast(accent, Color32::WHITE) >= contrast(accent, dark_ink) {
        Color32::WHITE
    } else {
        dark_ink
    }
}

/// `colour`, moved toward `toward` (black on a light page, white on a dark one) just far enough to reach `least` against every colour
/// of `behind`.
fn lifted(colour: Color32, toward: Color32, behind: &[Color32], least: f32) -> Color32 {
    (0..=100)
        .step_by(2)
        .map(|pct| mix(toward, colour, pct as f32))
        .find(|c| behind.iter().all(|bg| contrast(*c, *bg) >= least))
        .unwrap_or(toward)
}

/// Push the colours of the chrome apart from what is behind them: the page, the cards, and the grounds that the controls of a card are
/// drawn on (a box, a button under the pointer, a chosen button).
fn sharpen(l: &mut Look) {
    let toward = if l.dark {
        Color32::WHITE
    } else {
        Color32::BLACK
    };
    // an opaque tint of the accent is what a chosen button is on, whatever the accent
    l.selection_fill = l.accent_soft;
    let behind = [
        l.surface,
        l.bg,
        l.surface_alt,
        l.well,
        l.hover,
        l.accent_soft,
    ];
    let edge_behind = [l.surface, l.bg, l.surface_alt, l.well];
    let text_like = |c: Color32| lifted(c, toward, &behind, HIGH_TEXT);
    l.text = toward;
    l.dim = text_like(l.dim);
    l.accent = text_like(l.accent);
    l.accent2 = text_like(l.accent2);
    l.good = text_like(l.good);
    l.warn = text_like(l.warn);
    l.bad = text_like(l.bad);
    l.border = lifted(l.border, toward, &edge_behind, HIGH_EDGE);
    l.on_accent = ink_for(l.accent);
    l.outline = HIGH_OUTLINE;
    l.high_contrast = true;
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Appearance {
    /// How large everything is drawn, in per cent (one of `SCALES`).
    pub scale: u16,
    pub density: Density,
    pub accent: Accent,
    pub contrast: Contrast,
    pub colours: Colours,
}

impl Default for Appearance {
    fn default() -> Self {
        Appearance {
            scale: NORMAL_SCALE,
            density: Density::default(),
            accent: Accent::default(),
            contrast: Contrast::default(),
            colours: Colours::default(),
        }
    }
}

impl Appearance {
    /// The factor that the whole window is zoomed by.
    pub fn zoom(self) -> f32 {
        f32::from(self.scale) / 100.0
    }

    /// The colours of the chrome on the light or dark theme with these settings.
    pub fn look(self, dark: bool) -> Look {
        let mut look = design::look(dark);
        self.accent.tint(&mut look);
        if self.contrast == Contrast::High {
            sharpen(&mut look);
        }
        look.safe_colours = self.colours == Colours::Safe;
        look
    }

    /// The colours of the sequence editor with these settings.
    pub fn palette(self, mut palette: Palette) -> Palette {
        palette.safe_colours = self.colours == Colours::Safe;
        if self.contrast == Contrast::High {
            let toward = if palette.dark {
                Color32::WHITE
            } else {
                Color32::BLACK
            };
            palette.fg = toward;
            palette.dim = lifted(palette.dim, toward, &[palette.text_bg], HIGH_TEXT);
            palette.line = lifted(palette.line, toward, &[palette.text_bg], HIGH_EDGE);
        }
        palette
    }

    /// The next larger size (`by` of 1) or smaller one (-1) in `SCALES`, and the normal one for 0. The ends stay where they are.
    pub fn stepped(self, by: i8) -> Appearance {
        let scale = if by == 0 {
            NORMAL_SCALE
        } else {
            let at = SCALES
                .iter()
                .position(|&s| s == self.scale)
                .or_else(|| SCALES.iter().position(|&s| s == NORMAL_SCALE))
                .unwrap_or(0);
            let to = (at as i32 + i32::from(by)).clamp(0, SCALES.len() as i32 - 1);
            SCALES[to as usize]
        };
        Appearance { scale, ..self }
    }

    /// The text that `appearance.txt` holds: the theme and each setting on a line of its own.
    pub fn encode(&self, theme: ThemeMode) -> String {
        format!(
            "theme {}\nscale {}\ndensity {}\naccent {}\ncontrast {}\ncolours {}\n",
            theme.key(),
            self.scale,
            self.density.key(),
            self.accent.key(),
            self.contrast.key(),
            self.colours.key()
        )
    }

    /// The theme and settings that `text` says. A line that is not understood (an unknown word, or a value that is not one of the
    /// choices: a size that is not in `SCALES`) leaves that setting as it was; the last of two lines for a setting counts.
    pub fn decode(text: &str) -> (ThemeMode, Appearance) {
        Appearance::default().decode_over(ThemeMode::default(), text)
    }

    /// `decode`, with what `text` does not name left as `self` and `theme` have it.
    fn decode_over(self, theme: ThemeMode, text: &str) -> (ThemeMode, Appearance) {
        let (mut theme, mut a) = (theme, self);
        // an editor on some systems starts a file with a byte-order mark, which is not white space to Rust but is not a part of the word
        let text = text.strip_prefix('\u{feff}').unwrap_or(text);
        for line in text.lines() {
            let mut words = line.split_whitespace();
            let (Some(key), Some(value), None) = (words.next(), words.next(), words.next()) else {
                continue;
            };
            match key {
                "theme" => theme = ThemeMode::from_key(value).unwrap_or(theme),
                "scale" => {
                    a.scale = value
                        .parse()
                        .ok()
                        .filter(|s| SCALES.contains(s))
                        .unwrap_or(a.scale)
                }
                "density" => a.density = Density::from_key(value).unwrap_or(a.density),
                "accent" => a.accent = Accent::from_key(value).unwrap_or(a.accent),
                "contrast" => a.contrast = Contrast::from_key(value).unwrap_or(a.contrast),
                "colours" => a.colours = Colours::from_key(value).unwrap_or(a.colours),
                _ => {}
            }
        }
        (theme, a)
    }

    /// `self` with the settings written on a command line put on it: `scale=150,contrast=high` (what `--appearance` takes). What the
    /// line does not name stays as it is.
    pub fn with_spec(self, spec: &str) -> Appearance {
        self.decode_over(
            ThemeMode::default(),
            &spec.replace(',', "\n").replace('=', " "),
        )
        .1
    }

    /// The settings that a command line names, the defaults for what it does not.
    #[cfg(test)]
    pub fn from_spec(spec: &str) -> Appearance {
        Appearance::default().with_spec(spec)
    }
}
