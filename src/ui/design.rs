//! The platform's look: one palette, one set of spacings, and the few widgets every page shares (cards,
//! page headers, chips, buttons), so the whole program reads as a single interface. The base and
//! amino-acid colours used by the alignment viewer, the chromatogram and the structure viewer live here too.

use eframe::egui::{
    self, epaint::Shadow, Align, Color32, CornerRadius, FontId, Frame, Layout, Margin, Painter,
    Rect, Response, RichText, Sense, Stroke, Ui, Visuals,
};

/// Rounded corners of cards and windows.
pub const CARD_RADIUS: u8 = 12;
/// Small enough that a check box stays a square and does not read as a radio button.
pub const WIDGET_RADIUS: u8 = 4;
/// The corners of the buttons the helpers here draw: rounder than a plain widget, which has to keep its check boxes square.
pub const BUTTON_RADIUS: u8 = 8;
/// Space around the content of a page.
pub const PAGE_MARGIN: i8 = 18;
const CARD_MARGIN: i8 = 14;
/// A chip's text size, the space between the text and the edge, and the rounding of its corners.
const CHIP_TEXT: f32 = 11.5;
const CHIP_PAD_X: f32 = 8.0;
const CHIP_PAD_Y: f32 = 2.0;
const CHIP_RADIUS: u8 = 10;

/// Every colour the chrome uses, for one of the two themes.
#[derive(Clone, Copy)]
pub struct Look {
    pub dark: bool,
    /// The page behind the cards.
    pub bg: Color32,
    pub surface: Color32,
    pub surface_alt: Color32,
    /// The empty part of a slider and the inside of an unchecked box: darker than a card, so it shows on one.
    pub well: Color32,
    pub border: Color32,
    pub text: Color32,
    pub dim: Color32,
    /// The colour of the program: a deep indigo (a lighter one on the dark theme).
    pub accent: Color32,
    pub accent_soft: Color32,
    pub on_accent: Color32,
    /// The second colour of the mark and of what is new or alive: a teal.
    pub accent2: Color32,
    pub good: Color32,
    pub warn: Color32,
    pub bad: Color32,
    /// What a row or a plain button is filled with under the pointer.
    pub hover: Color32,
    /// The soft shadow under cards, windows and menus.
    pub shadow: Color32,
    /// What dims the window behind a dialog.
    pub scrim: Color32,
}

pub fn look(dark: bool) -> Look {
    let rgb = Color32::from_rgb;
    let rgba = Color32::from_rgba_unmultiplied;
    if dark {
        Look {
            dark,
            bg: rgb(0x0b, 0x0e, 0x15),
            surface: rgb(0x13, 0x18, 0x25),
            surface_alt: rgb(0x18, 0x1e, 0x2d),
            well: rgb(0x23, 0x2c, 0x40),
            border: rgb(0x25, 0x2e, 0x42),
            text: rgb(0xe9, 0xed, 0xf6),
            dim: rgb(0x9a, 0xa4, 0xb9),
            accent: rgb(0x8f, 0x8b, 0xff),
            accent_soft: rgb(0x1f, 0x21, 0x46),
            on_accent: rgb(0x0b, 0x0b, 0x24),
            accent2: rgb(0x3d, 0xd1, 0xc5),
            good: rgb(0x4c, 0xc3, 0x8a),
            warn: rgb(0xe0, 0xa0, 0x30),
            bad: rgb(0xf0, 0x71, 0x78),
            hover: rgb(0x1d, 0x24, 0x37),
            shadow: rgba(0, 0, 0, 90),
            scrim: rgba(0, 0, 0, 140),
        }
    } else {
        Look {
            dark,
            bg: rgb(0xf5, 0xf6, 0xfa),
            surface: rgb(0xff, 0xff, 0xff),
            surface_alt: rgb(0xf0, 0xf2, 0xf8),
            well: rgb(0xe4, 0xe8, 0xf1),
            border: rgb(0xe1, 0xe5, 0xef),
            text: rgb(0x13, 0x18, 0x26),
            dim: rgb(0x5a, 0x63, 0x76),
            accent: rgb(0x50, 0x46, 0xe5),
            accent_soft: rgb(0xec, 0xeb, 0xfc),
            on_accent: Color32::WHITE,
            accent2: rgb(0x0f, 0x9d, 0x94),
            good: rgb(0x1b, 0x9a, 0x59),
            warn: rgb(0xb2, 0x6a, 0x00),
            bad: rgb(0xd1, 0x3d, 0x3d),
            hover: rgb(0xee, 0xf0, 0xf7),
            shadow: rgba(19, 24, 38, 22),
            scrim: rgba(15, 18, 30, 96),
        }
    }
}

/// The shadow of a floating window or a menu, and the lighter one of a card.
fn shadow(l: &Look, offset: i8, blur: u8) -> Shadow {
    Shadow {
        offset: [0, offset],
        blur,
        spread: 0,
        color: l.shadow,
    }
}

/// Put the platform's colours, rounding and spacing on a set of visuals.
pub fn decorate(v: &mut Visuals, l: &Look) {
    let widget = CornerRadius::same(WIDGET_RADIUS);
    v.panel_fill = l.bg;
    v.window_fill = l.surface;
    v.extreme_bg_color = l.surface;
    v.faint_bg_color = l.surface_alt;
    v.window_stroke = Stroke::new(1.0, l.border);
    v.window_corner_radius = CornerRadius::same(CARD_RADIUS);
    v.window_shadow = shadow(l, 10, 28);
    v.popup_shadow = shadow(l, 6, 18);
    v.menu_corner_radius = CornerRadius::same(BUTTON_RADIUS);
    v.hyperlink_color = l.accent;
    v.warn_fg_color = l.warn;
    v.error_fg_color = l.bad;
    v.selection.bg_fill = l.accent.gamma_multiply(if l.dark { 0.45 } else { 0.30 });
    v.selection.stroke = Stroke::new(1.0, l.accent);
    let w = &mut v.widgets;
    w.noninteractive.bg_fill = l.surface;
    w.noninteractive.weak_bg_fill = l.surface;
    w.noninteractive.bg_stroke = Stroke::new(1.0, l.border);
    w.noninteractive.fg_stroke = Stroke::new(1.0, l.text);
    for state in [
        &mut w.noninteractive,
        &mut w.inactive,
        &mut w.hovered,
        &mut w.active,
        &mut w.open,
    ] {
        state.corner_radius = widget;
    }
    w.inactive.bg_fill = l.well;
    w.inactive.weak_bg_fill = l.surface_alt;
    w.inactive.bg_stroke = Stroke::new(1.0, l.border);
    w.inactive.fg_stroke = Stroke::new(1.0, l.text);
    w.hovered.bg_fill = l.hover;
    w.hovered.weak_bg_fill = l.hover;
    w.hovered.bg_stroke = Stroke::new(1.0, l.accent.gamma_multiply(0.6));
    w.hovered.fg_stroke = Stroke::new(1.5, l.text);
    w.active.bg_fill = l.accent_soft;
    w.active.weak_bg_fill = l.accent_soft;
    w.active.bg_stroke = Stroke::new(1.0, l.accent);
    w.active.fg_stroke = Stroke::new(1.5, l.text);
    w.open.bg_fill = l.accent_soft;
    w.open.weak_bg_fill = l.accent_soft;
}

/// Spacing for every widget, for both themes.
pub fn space(ctx: &egui::Context) {
    ctx.all_styles_mut(|s| {
        s.spacing.item_spacing = egui::vec2(8.0, 6.0);
        s.spacing.button_padding = egui::vec2(10.0, 4.0);
        s.spacing.interact_size.y = 26.0;
        s.spacing.window_margin = Margin::same(12);
        s.spacing.menu_margin = Margin::same(6);
    });
}

/// The margin every workspace page keeps around its content.
pub fn page<R>(ui: &mut Ui, add: impl FnOnce(&mut Ui) -> R) -> R {
    Frame::new()
        .inner_margin(Margin::same(PAGE_MARGIN))
        .show(ui, add)
        .inner
}

/// Where the first tool window opens: beside the rail (92 wide) and under the top bar (46 high), with a little air. Each
/// window is shifted a few steps by its title, so that several open at once do not hide each other.
const WINDOW_ORIGIN: (f32, f32) = (118.0, 74.0);
const WINDOW_STEP: f32 = 28.0;
const WINDOW_STEPS: u32 = 6;

/// How far down and to the right of the first place a window with this title opens: the same every time.
fn window_shift(title: &str) -> f32 {
    let step = title
        .bytes()
        .fold(0u32, |h, b| h.wrapping_mul(31).wrapping_add(u32::from(b)))
        % WINDOW_STEPS;
    step as f32 * WINDOW_STEP
}

/// A floating tool window, placed inside the working area and not over the bar and the rail.
pub fn window<'a>(title: impl Into<egui::WidgetText>) -> egui::Window<'a> {
    let title: egui::WidgetText = title.into();
    let shift = window_shift(title.text());
    egui::Window::new(title)
        .default_pos(egui::pos2(WINDOW_ORIGIN.0 + shift, WINDOW_ORIGIN.1 + shift))
}

/// A framed surface with rounded corners and a soft shadow.
pub fn card<R>(ui: &mut Ui, l: &Look, add: impl FnOnce(&mut Ui) -> R) -> R {
    Frame::new()
        .fill(l.surface)
        .stroke(Stroke::new(1.0, l.border))
        .corner_radius(CornerRadius::same(CARD_RADIUS))
        .inner_margin(Margin::same(CARD_MARGIN))
        .shadow(shadow(l, 1, 6))
        .show(ui, add)
        .inner
}

/// The strip every page starts with: a title, one line saying what the page is for, and the page's
/// actions on the right.
pub fn page_header(
    ui: &mut Ui,
    l: &Look,
    title: &str,
    subtitle: &str,
    actions: impl FnOnce(&mut Ui),
) {
    ui.horizontal(|ui| {
        ui.vertical(|ui| {
            ui.label(RichText::new(title).size(22.0).strong().color(l.text));
            ui.label(RichText::new(subtitle).size(13.0).color(l.dim));
        });
        ui.with_layout(Layout::right_to_left(Align::Center), actions);
    });
    ui.add_space(6.0);
}

/// A small heading above a group of controls.
pub fn section(ui: &mut Ui, l: &Look, text: &str) {
    ui.add_space(4.0);
    ui.label(
        RichText::new(text.to_uppercase())
            .size(11.0)
            .strong()
            .color(l.dim),
    );
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tone {
    Neutral,
    Accent,
    Good,
    Bad,
}

/// A small rounded label: a status, a count, a program that was found. It is drawn by hand, as one block that
/// moves to the next row of a wrapping row when it does not fit, never as text that wraps inside it.
pub fn chip(ui: &mut Ui, l: &Look, text: &str, tone: Tone) {
    let (fg, bg) = match tone {
        Tone::Neutral => (l.dim, l.surface_alt),
        Tone::Accent => (l.accent, l.accent_soft),
        Tone::Good => (l.good, l.good.gamma_multiply(0.15)),
        Tone::Bad => (l.bad, l.bad.gamma_multiply(0.15)),
    };
    let galley = ui
        .painter()
        .layout_no_wrap(text.to_string(), FontId::proportional(CHIP_TEXT), fg);
    let pad = egui::vec2(CHIP_PAD_X, CHIP_PAD_Y);
    let (rect, _) = ui.allocate_exact_size(galley.size() + 2.0 * pad, Sense::hover());
    ui.painter()
        .rect_filled(rect, CornerRadius::same(CHIP_RADIUS), bg);
    ui.painter().galley(rect.min + pad, galley, fg);
}

/// The one button on a page that does what the page is for.
pub fn primary_button(ui: &mut Ui, l: &Look, text: &str) -> Response {
    ui.add(
        egui::Button::new(RichText::new(text).strong().color(l.on_accent))
            .fill(l.accent)
            .stroke(Stroke::NONE)
            .corner_radius(CornerRadius::same(BUTTON_RADIUS))
            .min_size(egui::vec2(0.0, 30.0)),
    )
}

/// A button next to the main one: an outline, and the text in the ordinary colour.
pub fn ghost_button(ui: &mut Ui, l: &Look, text: &str) -> Response {
    ui.add(
        egui::Button::new(RichText::new(text).color(l.text))
            .fill(Color32::TRANSPARENT)
            .stroke(Stroke::new(1.0, l.border))
            .corner_radius(CornerRadius::same(BUTTON_RADIUS))
            .min_size(egui::vec2(0.0, 30.0)),
    )
}

/// A square with rounded corners in a tint of the accent, with an icon drawn inside it by `paint`: the head of a card.
pub fn icon_tile(ui: &mut Ui, l: &Look, size: f32, paint: impl FnOnce(&Painter, Rect)) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(size, size), Sense::hover());
    ui.painter()
        .rect_filled(rect, CornerRadius::same(BUTTON_RADIUS + 2), l.accent_soft);
    paint(ui.painter(), rect.shrink(size * 0.16));
}

/// The name of the key that the shortcuts of this program are made with.
pub fn command_key() -> &'static str {
    if cfg!(target_os = "macos") {
        "Cmd"
    } else {
        "Ctrl"
    }
}

/// A shortcut drawn as a little key: `Cmd K`.
pub fn keycap(ui: &mut Ui, l: &Look, text: &str) {
    let galley =
        ui.painter()
            .layout_no_wrap(text.to_string(), FontId::proportional(CHIP_TEXT), l.dim);
    let pad = egui::vec2(6.0, 2.0);
    let (rect, _) = ui.allocate_exact_size(galley.size() + 2.0 * pad, Sense::hover());
    ui.painter()
        .rect_filled(rect, CornerRadius::same(5), l.surface_alt);
    ui.painter().rect_stroke(
        rect,
        CornerRadius::same(5),
        Stroke::new(1.0, l.border),
        egui::StrokeKind::Inside,
    );
    ui.painter().galley(rect.min + pad, galley, l.dim);
}

// ---------------------------------------------------------------------------
// Colours of bases and residues
// ---------------------------------------------------------------------------

/// The colour of a trace, or of a base called along it, for base `b`, assigned as in the original
/// (AppMain.tcl:67217, `abi_info(colors,*)`): A green, C blue, G the foreground colour (black, or light on a dark
/// page), T red, and magenta for N, every ambiguity code and the gap. The shades are a little deeper than Tk's pure
/// green, blue, red and magenta so that the letters can be read on a white card.
pub fn trace_color(b: u8, dark: bool) -> Color32 {
    match b.to_ascii_uppercase() {
        b'A' => Color32::from_rgb(0x1b, 0x9e, 0x3e),
        b'C' => Color32::from_rgb(0x1e, 0x50, 0xd6),
        b'G' if dark => Color32::from_gray(0xe0),
        b'G' => Color32::from_gray(0x22),
        b'T' => Color32::from_rgb(0xd6, 0x28, 0x28),
        _ => Color32::from_rgb(0xc0, 0x26, 0xd3),
    }
}

/// The background of an alignment cell for a nucleotide, if it has a colour.
pub fn base_fill(b: u8, dark: bool) -> Option<Color32> {
    let c = match (b.to_ascii_uppercase(), dark) {
        (b'A', false) => [0x9a, 0xd9, 0x9a],
        (b'C', false) => [0x9d, 0xbd, 0xf2],
        (b'G', false) => [0xf2, 0xc9, 0x86],
        (b'T' | b'U', false) => [0xf0, 0xa4, 0xa4],
        (b'A', true) => [0x2e, 0x7d, 0x4f],
        (b'C', true) => [0x2f, 0x5f, 0xa8],
        (b'G', true) => [0xa8, 0x73, 0x21],
        (b'T' | b'U', true) => [0xa8, 0x3a, 0x3a],
        _ => return None,
    };
    Some(Color32::from_rgb(c[0], c[1], c[2]))
}

/// The background of an alignment cell for an amino acid (the Clustal groups: hydrophobic blue, positive red,
/// negative magenta, polar green, cysteine pink, glycine orange, proline yellow, aromatic cyan).
pub fn residue_fill(aa: u8, dark: bool) -> Option<Color32> {
    let c = match aa.to_ascii_uppercase() {
        b'A' | b'I' | b'L' | b'M' | b'F' | b'W' | b'V' => [0x80, 0xa0, 0xf0],
        b'K' | b'R' => [0xf0, 0x60, 0x55],
        b'D' | b'E' => [0xc0, 0x60, 0xc0],
        b'N' | b'Q' | b'S' | b'T' => [0x4c, 0xc8, 0x4c],
        b'C' => [0xf0, 0x90, 0x90],
        b'G' => [0xf0, 0x90, 0x48],
        b'P' => [0xd0, 0xd0, 0x30],
        b'H' | b'Y' => [0x30, 0xb4, 0xb4],
        _ => return None,
    };
    let c = Color32::from_rgb(c[0], c[1], c[2]);
    // on a light page the colours are softened so the letters stay easy to read, on a dark one darkened
    Some(if dark {
        c.gamma_multiply(0.55)
    } else {
        super::theme::mix(c, Color32::WHITE, 55.0)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn both_themes_have_readable_text_on_their_surfaces() {
        for dark in [false, true] {
            let l = look(dark);
            let lum = |c: Color32| {
                0.299 * f32::from(c.r()) + 0.587 * f32::from(c.g()) + 0.114 * f32::from(c.b())
            };
            assert!(
                (lum(l.text) - lum(l.surface)).abs() > 120.0,
                "text on the surface, dark={dark}"
            );
            assert!(
                (lum(l.dim) - lum(l.surface)).abs() > 60.0,
                "dim text on the surface, dark={dark}"
            );
            assert!(
                (lum(l.on_accent) - lum(l.accent)).abs() > 60.0,
                "text on the accent, dark={dark}"
            );
        }
    }

    #[test]
    fn windows_open_in_a_few_places_inside_the_working_area_and_the_same_one_every_time() {
        let titles = [
            "Find",
            "Digest",
            "Jump To",
            "Translate",
            "Tree Tools",
            "Enzyme Selector",
            "Keyboard shortcuts",
            "PCR and Primers",
            "Calculator: Enter a calculation",
        ];
        let shifts: Vec<f32> = titles.iter().map(|t| window_shift(t)).collect();
        assert!(shifts
            .iter()
            .all(|s| (0.0..=(WINDOW_STEPS - 1) as f32 * WINDOW_STEP).contains(s)));
        assert_eq!(window_shift("Find"), window_shift("Find"));
        let mut distinct = shifts.clone();
        distinct.sort_by(f32::total_cmp);
        distinct.dedup();
        assert!(distinct.len() >= 4, "{shifts:?}");
        // the first place is clear of the rail and the top bar
        assert!(WINDOW_ORIGIN.0 > 92.0 && WINDOW_ORIGIN.1 > 46.0);
    }

    #[test]
    fn the_accent_and_the_second_colour_read_on_the_surface_in_both_themes() {
        for dark in [false, true] {
            let l = look(dark);
            let lum = |c: Color32| {
                0.299 * f32::from(c.r()) + 0.587 * f32::from(c.g()) + 0.114 * f32::from(c.b())
            };
            assert!(
                (lum(l.accent) - lum(l.surface)).abs() > 60.0,
                "accent, dark={dark}"
            );
            assert!(
                (lum(l.accent2) - lum(l.surface)).abs() > 60.0,
                "teal, dark={dark}"
            );
            assert_ne!(l.accent, l.accent2);
        }
    }

    #[test]
    fn bases_and_residues_get_colours_and_gaps_do_not() {
        for dark in [false, true] {
            assert!(b"ACGTU".iter().all(|&b| base_fill(b, dark).is_some()));
            assert!(base_fill(b'-', dark).is_none() && base_fill(b'N', dark).is_none());
            assert!(b"ACDEFGHIKLMNPQRSTVWY"
                .iter()
                .all(|&b| residue_fill(b, dark).is_some()));
            assert!(residue_fill(b'-', dark).is_none() && residue_fill(b'X', dark).is_none());
        }
        // lower case reads like upper case
        assert_eq!(base_fill(b'a', false), base_fill(b'A', false));
        assert_eq!(residue_fill(b'k', true), residue_fill(b'K', true));
    }

    #[test]
    fn the_four_trace_colours_differ() {
        let colours: std::collections::HashSet<_> =
            b"ACGT".iter().map(|&b| trace_color(b, false)).collect();
        assert_eq!(colours.len(), 4);
        assert_ne!(trace_color(b'G', false), trace_color(b'G', true));
    }

    #[test]
    fn n_the_ambiguity_codes_and_the_gap_are_all_magenta_as_in_the_original() {
        let magenta = trace_color(b'N', false);
        for b in *b"BDHKMRSVWYU-n*" {
            assert_eq!(trace_color(b, false), magenta, "{}", char::from(b));
        }
        assert!(!b"ACGT".iter().any(|&b| trace_color(b, false) == magenta));
        // lower case reads like upper case
        assert_eq!(trace_color(b'a', false), trace_color(b'A', false));
    }
}
