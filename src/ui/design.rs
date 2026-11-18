//! The platform's look: one palette, one set of spacings, and the few widgets every page shares (cards,
//! page headers, chips, buttons), so the whole program reads as a single interface. The base and
//! amino-acid colours used by the alignment viewer, the chromatogram and the structure viewer live here too.

use crate::ui::appearance::{safe, Density};
use crate::ui::fonts::FontWeight;
use eframe::egui::{
    self, epaint::Shadow, Align, Color32, CornerRadius, FontId, Frame, Layout, Margin, Painter,
    Rect, Response, RichText, ScrollArea, Sense, Stroke, Ui, Visuals,
};

/// Rounded corners of cards and windows.
pub const CARD_RADIUS: u8 = 12;
/// Small enough that a check box stays a square and does not read as a radio button.
pub const WIDGET_RADIUS: u8 = 4;
/// The corners of the buttons the helpers here draw: rounder than a plain widget, which has to keep its check boxes square.
pub const BUTTON_RADIUS: u8 = 8;
/// Space around the content of a page.
pub const PAGE_MARGIN: i8 = 18;
pub const CARD_MARGIN: i8 = 14;
/// A chip's text size, the space between the text and the edge, and the rounding of its corners.
const CHIP_TEXT: f32 = 11.5;
const CHIP_PAD_X: f32 = 8.0;
const CHIP_PAD_Y: f32 = 2.0;
const CHIP_RADIUS: u8 = 10;

/// Every colour the chrome uses, for one of the two themes.
#[derive(Clone, Copy, Debug, PartialEq)]
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
    /// What a chosen button, a selected row and a text selection are filled with (the accent shown through at a third, as it always was;
    /// an accent of the person's choosing, and the high-contrast setting, put an opaque tint there, so that the text of the accent
    /// reads on it).
    pub selection_fill: Color32,
    /// How thick the outlines of cards, buttons and boxes are (thicker with the high-contrast setting).
    pub outline: f32,
    /// Whether the high-contrast setting is on.
    pub high_contrast: bool,
    /// Whether bases, residues and trace letters are drawn in colours that stay apart for people who see colour differently.
    pub safe_colours: bool,
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
            selection_fill: rgb(0x8f, 0x8b, 0xff).gamma_multiply(0.45),
            outline: 1.0,
            high_contrast: false,
            safe_colours: false,
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
            selection_fill: rgb(0x50, 0x46, 0xe5).gamma_multiply(0.30),
            outline: 1.0,
            high_contrast: false,
            safe_colours: false,
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
    v.window_stroke = Stroke::new(l.outline, l.border);
    v.window_corner_radius = CornerRadius::same(CARD_RADIUS);
    v.window_shadow = shadow(l, 10, 28);
    v.popup_shadow = shadow(l, 6, 18);
    v.menu_corner_radius = CornerRadius::same(BUTTON_RADIUS);
    v.hyperlink_color = l.accent;
    v.warn_fg_color = l.warn;
    v.error_fg_color = l.bad;
    v.selection.bg_fill = l.selection_fill;
    v.selection.stroke = Stroke::new(l.outline, l.accent);
    let w = &mut v.widgets;
    w.noninteractive.bg_fill = l.surface;
    w.noninteractive.weak_bg_fill = l.surface;
    w.noninteractive.bg_stroke = Stroke::new(l.outline, l.border);
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
    w.inactive.bg_stroke = Stroke::new(l.outline, l.border);
    w.inactive.fg_stroke = Stroke::new(1.0, l.text);
    w.hovered.bg_fill = l.hover;
    w.hovered.weak_bg_fill = l.hover;
    w.hovered.bg_stroke = Stroke::new(l.outline, l.accent.gamma_multiply(0.6));
    w.hovered.fg_stroke = Stroke::new(1.5, l.text);
    w.active.bg_fill = l.accent_soft;
    w.active.weak_bg_fill = l.accent_soft;
    w.active.bg_stroke = Stroke::new(l.outline, l.accent);
    w.active.fg_stroke = Stroke::new(1.5, l.text);
    w.open.bg_fill = l.accent_soft;
    w.open.weak_bg_fill = l.accent_soft;
}

/// Spacing for every widget, for both themes, as roomy as `density` says.
pub fn space(ctx: &egui::Context, density: Density) {
    let spacing = density.spacing();
    ctx.all_styles_mut(|s| {
        s.spacing.item_spacing = spacing.item;
        s.spacing.button_padding = spacing.button_padding;
        s.spacing.interact_size.y = spacing.interact_height;
        s.spacing.window_margin = Margin::same(spacing.window_margin);
        s.spacing.menu_margin = Margin::same(spacing.menu_margin);
    });
}

/// The margin every workspace page keeps around its content.
pub fn page<R>(ui: &mut Ui, add: impl FnOnce(&mut Ui) -> R) -> R {
    Frame::new()
        .inner_margin(Margin::same(PAGE_MARGIN))
        .show(ui, add)
        .inner
}

/// How wide the controls of a workspace page are: the same on every page, so that the results begin at the same place
/// whichever page is open.
pub const CONTROLS_WIDTH: f32 = 320.0;

/// The controls of a workspace page: a column of one width at the left, which scrolls by itself when the window is short.
/// They are in one card, as wide as the column whatever is in it. `id` tells the pages' columns apart.
pub fn controls_panel(ui: &mut Ui, l: &Look, id: &'static str, add: impl FnOnce(&mut Ui)) {
    egui::Panel::left(id)
        .exact_size(CONTROLS_WIDTH)
        .resizable(false)
        .frame(Frame::new())
        .show(ui, |ui| {
            ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    card(ui, l, |ui| {
                        ui.set_min_width(ui.available_width());
                        add(ui);
                    });
                });
        });
}

/// The results of a workspace page: the room that is left beside its controls, a page's margin away from them.
pub fn results_area<R>(ui: &mut Ui, add: impl FnOnce(&mut Ui) -> R) -> R {
    egui::CentralPanel::default()
        .frame(Frame::new().inner_margin(Margin {
            left: PAGE_MARGIN,
            ..Margin::ZERO
        }))
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

thread_local! {
    /// Where the tool windows were left, as the settings folder kept them, to open them there again (set once, at the start).
    static SAVED: std::cell::RefCell<std::collections::HashMap<String, egui::Pos2>> = Default::default();
    /// The titles of the windows drawn since they were last taken: the places of these are noted after the frame.
    static SEEN: std::cell::RefCell<Vec<String>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Where the window with this title opens when no place was kept for it.
pub fn default_place(title: &str) -> egui::Pos2 {
    let shift = window_shift(title);
    egui::pos2(WINDOW_ORIGIN.0 + shift, WINDOW_ORIGIN.1 + shift)
}

/// The id a tool window with this title is drawn under, so that where it is can be read from egui's memory by the title.
pub fn window_id(title: &str) -> egui::Id {
    egui::Id::new(("tool window", title))
}

/// The places the tool windows were left at, for `window` to open them at.
pub fn set_saved_places(places: std::collections::HashMap<String, egui::Pos2>) {
    SAVED.with(|s| *s.borrow_mut() = places);
}

/// The titles of the windows that were drawn since this was last called (each once, in the order they were first drawn).
pub fn take_seen_windows() -> Vec<String> {
    SEEN.with(|s| std::mem::take(&mut *s.borrow_mut()))
}

/// A floating tool window, placed inside the working area and not over the bar and the rail: where it was left, when a place
/// was kept for it (never above or left of `places::FLOOR`), else where its title says.
pub fn window<'a>(title: impl Into<egui::WidgetText>) -> egui::Window<'a> {
    let title: egui::WidgetText = title.into();
    let name = title.text().to_string();
    SEEN.with(|s| {
        let mut seen = s.borrow_mut();
        if !seen.contains(&name) {
            seen.push(name.clone());
        }
    });
    let place = SAVED.with(|s| s.borrow().get(&name).copied()).map_or_else(
        || default_place(&name),
        |p| {
            let (x, y) = super::places::floored(p.x, p.y);
            egui::pos2(x, y)
        },
    );
    egui::Window::new(title)
        .id(window_id(&name))
        .default_pos(place)
}

/// A framed surface with rounded corners and a soft shadow.
pub fn card<R>(ui: &mut Ui, l: &Look, add: impl FnOnce(&mut Ui) -> R) -> R {
    Frame::new()
        .fill(l.surface)
        .stroke(Stroke::new(l.outline, l.border))
        .corner_radius(CornerRadius::same(CARD_RADIUS))
        .inner_margin(Margin::same(CARD_MARGIN))
        .shadow(shadow(l, 1, 6))
        .show(ui, add)
        .inner
}

/// The room between two cards of a page's results.
pub const CARD_GAP: f32 = 10.0;

/// A card of a page's results with its title, as wide as the results whatever is in it.
pub fn titled_card(ui: &mut Ui, l: &Look, title: &str, add: impl FnOnce(&mut Ui)) {
    card(ui, l, |ui| {
        ui.set_min_width(ui.available_width());
        ui.label(RichText::new(title).size(15.0).semibold().color(l.text));
        ui.add_space(4.0);
        add(ui);
    });
}

/// The frame of a view (an alignment, a tree, a structure, a trace): the surface and the outline of a card around a picture
/// that fills it from edge to edge, so without a card's margin, and flat, since the picture is the page's result and not
/// one card among others.
pub fn view_frame<R>(ui: &mut Ui, l: &Look, add: impl FnOnce(&mut Ui) -> R) -> R {
    Frame::new()
        .fill(l.surface)
        .stroke(Stroke::new(l.outline, l.border))
        .corner_radius(CornerRadius::same(CARD_RADIUS))
        .show(ui, add)
        .inner
}

/// What a page says where its result will be, until there is one: what will appear, and the buttons that get it started
/// (the first is the main one). Returns the number of the button that was clicked.
pub fn empty_state(
    ui: &mut Ui,
    l: &Look,
    title: &str,
    text: &str,
    buttons: &[&str],
) -> Option<usize> {
    card(ui, l, |ui| {
        ui.label(RichText::new(title).size(17.0).semibold().color(l.text));
        ui.label(RichText::new(text).color(l.dim));
        ui.add_space(8.0);
        let mut hit = None;
        ui.horizontal_wrapped(|ui| {
            for (k, label) in buttons.iter().enumerate() {
                let button = if k == 0 {
                    primary_button(ui, l, label)
                } else {
                    ghost_button(ui, l, label)
                };
                if button.clicked() {
                    hit = Some(k);
                }
            }
        });
        hit
    })
}

/// The strip every page starts with: a title, one line saying what the page is for, and the page's
/// actions on the right.
///
/// The actions keep their place at the right whatever the width of the window, and the title and the line of purpose have
/// the room that is left: a line that is longer than that is cut short (the pointer shows it whole), never drawn under the
/// buttons, which is what a small window did until November 2026. The title is never cut: where it has no room beside
/// the actions, they go on a row of their own under it, and on to the next when that is full.
pub fn page_header(
    ui: &mut Ui,
    l: &Look,
    title: &str,
    subtitle: &str,
    actions: impl FnOnce(&mut Ui),
) {
    // how wide the actions were when they were last drawn beside the title: nothing is known of them before they are
    // drawn, so the first frame of a page gives the title the whole row and the next puts it right
    let id = egui::Id::new(("page_header", title));
    let taken: f32 = ui.data(|d| d.get_temp(id)).unwrap_or(0.0);
    let row = ui.available_width();
    let beside = row - taken - HEADER_GAP;
    let title = RichText::new(title).size(22.0).semibold().color(l.text);
    // the title decides whether the actions have room beside it: it is shown whole
    let title_wide = egui::WidgetText::from(title.clone())
        .into_galley(
            ui,
            Some(egui::TextWrapMode::Extend),
            f32::INFINITY,
            egui::TextStyle::Heading,
        )
        .size()
        .x;
    let words = |ui: &mut Ui, room: f32| {
        ui.vertical(|ui| {
            ui.set_max_width(room);
            ui.label(title);
            let purpose = RichText::new(subtitle).size(13.0).color(l.dim);
            ui.add(egui::Label::new(purpose).truncate());
        });
    };
    if beside >= title_wide {
        ui.horizontal(|ui| {
            words(ui, beside);
            let drawn = ui.with_layout(Layout::right_to_left(Align::Center), actions);
            let used = drawn.response.rect.width();
            if (used - taken).abs() > 0.5 {
                ui.data_mut(|d| d.insert_temp(id, used));
                ui.ctx().request_repaint();
            }
        });
    } else {
        // how wide the actions are in one row is what was measured when they were beside the title: it is kept, so that
        // a wider window puts them back there
        words(ui, row);
        // a row as high as a button, and as many more as the actions need: asked for with its size, since a layout that
        // is only given the room that is left takes all of it, and the page under the header would have none
        let wrapping = Layout::right_to_left(Align::Center).with_main_wrap(true);
        let a_row = egui::vec2(row, ui.spacing().interact_size.y);
        ui.allocate_ui_with_layout(a_row, wrapping, actions);
    }
    ui.add_space(6.0);
}

/// The room between the title of a page and its actions.
const HEADER_GAP: f32 = 12.0;

/// How large the letters of a heading above a group of controls are.
pub const SECTION_SIZE: f32 = 11.0;

/// A small heading above a group of controls.
pub fn section(ui: &mut Ui, l: &Look, text: &str) {
    ui.add_space(4.0);
    ui.label(
        RichText::new(text.to_uppercase())
            .size(SECTION_SIZE)
            .semibold()
            .color(l.dim),
    );
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tone {
    Neutral,
    Accent,
    Good,
    Warn,
    Bad,
}

/// A small rounded label: a status, a count, a program that was found. It is drawn by hand, as one block that
/// moves to the next row of a wrapping row when it does not fit, never as text that wraps inside it.
/// The line that says what is selected on a page (`docs/DESIGN.md`, *One selection on every page*): `words`, or `hint`
/// while nothing is selected, with `links` at the right (each its words and what the pointer is told over it). It is one
/// row whatever it says, so that what lies under the pointer stays where it is when something is clicked: the links keep
/// their room, and words too long for the rest are cut short there and read whole under the pointer. Returns the link
/// that was clicked, by its place in `links`.
pub fn selection_line(
    ui: &mut Ui,
    l: &Look,
    hint: &str,
    words: Option<&str>,
    links: &[(&str, &str)],
) -> Option<usize> {
    let mut clicked = None;
    ui.horizontal(|ui| {
        let Some(words) = words else {
            ui.label(RichText::new(hint).size(12.0).color(l.dim));
            return;
        };
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            // from the right: the last link first
            for (k, (label, told)) in links.iter().enumerate().rev() {
                if ui.link(*label).on_hover_text(*told).clicked() {
                    clicked = Some(k);
                }
            }
            ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
                ui.add(egui::Label::new(RichText::new(words).size(12.5).semibold()).truncate())
                    .on_hover_text(words);
            });
        });
    });
    clicked
}

pub fn chip(ui: &mut Ui, l: &Look, text: &str, tone: Tone) {
    let (fg, bg) = match tone {
        Tone::Neutral => (l.dim, l.surface_alt),
        Tone::Accent => (l.accent, l.accent_soft),
        Tone::Good => (l.good, l.good.gamma_multiply(0.15)),
        Tone::Warn => (l.warn, l.warn.gamma_multiply(0.15)),
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
        egui::Button::new(RichText::new(text).semibold().color(l.on_accent))
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
            .stroke(Stroke::new(l.outline, l.border))
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
        Stroke::new(l.outline, l.border),
        egui::StrokeKind::Inside,
    );
    ui.painter().galley(rect.min + pad, galley, l.dim);
}

/// The least contrast that letters need against what is behind them (WCAG level AA for text).
pub const READABLE: f32 = 4.5;

/// How well two colours can be told apart as the Web Content Accessibility Guidelines measure it: 1 for the same colour and 21 for
/// black on white. Text wants 4.5 or more, the edge of a box or button 3.
pub fn contrast(a: Color32, b: Color32) -> f32 {
    let light = |c: Color32| {
        let linear = |v: u8| {
            let v = f32::from(v) / 255.0;
            if v <= 0.04045 {
                v / 12.92
            } else {
                ((v + 0.055) / 1.055).powf(2.4)
            }
        };
        0.2126 * linear(c.r()) + 0.7152 * linear(c.g()) + 0.0722 * linear(c.b())
    };
    let (x, y) = (light(a), light(b));
    (x.max(y) + 0.05) / (x.min(y) + 0.05)
}

/// The colour of letters on `fill`: the theme's own text colour when that reads well there, otherwise black or white, whichever reads better.
pub fn ink_on(l: &Look, fill: Color32) -> Color32 {
    if contrast(l.text, fill) >= READABLE {
        return l.text;
    }
    let better = |a: &Color32, b: &Color32| contrast(*a, fill).total_cmp(&contrast(*b, fill));
    [Color32::BLACK, Color32::WHITE]
        .into_iter()
        .max_by(better)
        .unwrap_or(l.text)
}

/// The colour of the letter in a cell that is filled with `fill` (or has no fill): the theme's text colour, as it always was, except
/// where the person asked for the safe colours or high contrast, which choose black or white by the cell (`ink_on`). The one rule
/// for the alignment, the sample in the Appearance window and anything else that writes letters on these fills.
pub fn letter_ink(l: &Look, fill: Option<Color32>) -> Color32 {
    match fill {
        Some(fill) if l.safe_colours || l.high_contrast => ink_on(l, fill),
        _ => l.text,
    }
}

// ---------------------------------------------------------------------------
// Colours of bases and residues
// ---------------------------------------------------------------------------

/// The colour of a trace, or of a base called along it, for base `b`, assigned as in the original
/// (AppMain.tcl:67217, `abi_info(colors,*)`): A green, C blue, G the foreground colour (black, or light on a dark
/// page), T red, and magenta for N, every ambiguity code and the gap. The shades are a little deeper than Tk's pure
/// green, blue, red and magenta so that the letters can be read on a white card.
pub fn trace_color(b: u8, l: &Look) -> Color32 {
    let dark = l.dark;
    if l.safe_colours {
        return safe::trace(b, dark);
    }
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
pub fn base_fill(b: u8, l: &Look) -> Option<Color32> {
    let dark = l.dark;
    if l.safe_colours {
        return safe::base(b, dark);
    }
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
pub fn residue_fill(aa: u8, l: &Look) -> Option<Color32> {
    let dark = l.dark;
    if l.safe_colours {
        return safe::residue(aa, dark);
    }
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

    /// What a header drew: its texts with the room each takes, and the frames of its buttons.
    struct Header {
        texts: Vec<(String, Rect)>,
        buttons: Vec<Rect>,
        /// How much of the height of its room the header took: the page under it begins there.
        high: f32,
    }

    impl Header {
        fn room_of(&self, text: &str) -> Rect {
            let found = self.texts.iter().find(|(t, _)| t == text);
            found
                .unwrap_or_else(|| panic!("no {text:?} among {:?}", self.words()))
                .1
        }

        fn words(&self) -> Vec<&str> {
            self.texts.iter().map(|(t, _)| t.as_str()).collect()
        }
    }

    /// The header of a page drawn in a room `wide` points wide, with a button for each of `actions`, as it is once it has
    /// been drawn a few times (the first frame does not know how wide its buttons are).
    fn header_in(wide: f32, title: &str, purpose: &str, actions: &[&str]) -> Header {
        use egui::epaint::Shape;
        fn walk(shape: &Shape, out: &mut Header) {
            match shape {
                Shape::Text(t) => out.texts.push((
                    t.galley.text().to_string(),
                    t.galley.rect.translate(t.pos.to_vec2()),
                )),
                Shape::Rect(r) if r.fill.a() > 0 || r.stroke.width > 0.0 => {
                    out.buttons.push(r.rect)
                }
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, out)),
                _ => {}
            }
        }
        let ctx = egui::Context::default();
        let l = look(false);
        let mut drawn = Header {
            texts: Vec::new(),
            buttons: Vec::new(),
            high: 0.0,
        };
        let mut high = 0.0;
        for n in 0..4 {
            let mut input = crate::ui::headless::frame_input(n);
            input.events = vec![];
            let mut out = ctx.run_ui(input, |ui| {
                ui.allocate_ui(egui::vec2(wide, 300.0), |ui| {
                    let top = ui.cursor().top();
                    page_header(ui, &l, title, purpose, |ui| {
                        for action in actions {
                            let _ = ui.button(*action);
                        }
                    });
                    high = ui.cursor().top() - top;
                });
            });
            out.textures_delta.clear();
            drawn = Header {
                texts: Vec::new(),
                buttons: Vec::new(),
                high,
            };
            out.shapes.iter().for_each(|c| walk(&c.shape, &mut drawn));
        }
        drawn
    }

    const PURPOSE: &str =
        "A line that says at some length what the page is for, and so does not fit in a small window";
    const ACTIONS: [&str; 4] = [
        "Save figure",
        "Copy the result",
        "Save the result…",
        "Go on",
    ];

    #[test]
    fn a_header_with_room_shows_its_line_whole_and_its_actions_at_the_right() {
        let wide = 1100.0;
        let h = header_in(wide, "Trees", PURPOSE, &ACTIONS);
        let (title, purpose) = (h.room_of("Trees"), h.room_of(PURPOSE));
        assert_eq!(h.buttons.len(), ACTIONS.len());
        // whole: as wide as in a room with far more room than it needs
        let at_ease = header_in(3000.0, "Trees", PURPOSE, &ACTIONS).room_of(PURPOSE);
        assert_eq!(purpose.width(), at_ease.width());
        // two lines of words and the room under them: the header is no higher
        assert!(h.high < 70.0, "the header is {} high", h.high);
        let right = h.buttons.iter().map(|b| b.max.x).fold(0.0, f32::max);
        assert!(
            (right - (title.min.x + wide)).abs() < 2.0,
            "the last button ends where the room does"
        );
        // the actions are beside the title and the line, not under them
        assert!(h
            .buttons
            .iter()
            .all(|b| b.min.y < purpose.max.y && b.min.x > purpose.max.x));
    }

    #[test]
    fn a_line_of_purpose_that_has_no_room_beside_the_actions_is_cut_short_and_never_drawn_under_them(
    ) {
        let wide = 620.0;
        let h = header_in(wide, "Trees", PURPOSE, &ACTIONS);
        let title = h.room_of("Trees");
        // the words are all there for the pointer to show; what is drawn of them is narrower than in a room that has room
        let purpose = h.room_of(PURPOSE);
        let whole = header_in(1100.0, "Trees", PURPOSE, &ACTIONS).room_of(PURPOSE);
        assert!(
            purpose.width() < whole.width() - 20.0,
            "{} wide here, {} where it is whole",
            purpose.width(),
            whole.width()
        );
        // one row, and it ends before the first button begins
        assert_eq!(h.buttons.len(), ACTIONS.len());
        assert!(purpose.height() < 26.0, "{purpose:?}");
        assert!(h.high < 70.0, "the header is {} high", h.high);
        let first = h
            .buttons
            .iter()
            .map(|b| b.min.x)
            .fold(f32::INFINITY, f32::min);
        assert!(
            purpose.max.x < first,
            "{purpose:?} runs under a button that begins at {first}"
        );
        // the actions kept their place: beside the title, ending where the room ends
        assert!(h.buttons.iter().all(|b| b.min.y < purpose.max.y));
        let right = h.buttons.iter().map(|b| b.max.x).fold(0.0, f32::max);
        assert!((right - (title.min.x + wide)).abs() < 2.0);
    }

    #[test]
    fn actions_with_no_room_beside_the_title_go_under_it_and_on_to_another_row_when_one_is_full() {
        let wide = 260.0;
        let h = header_in(wide, "Phylogenetic estimation", PURPOSE, &ACTIONS);
        // the title is whole, and the line under it is cut to the room
        let title = h.room_of("Phylogenetic estimation");
        let (left, right) = (title.min.x, title.min.x + wide);
        let (_, purpose) = h
            .texts
            .iter()
            .find(|(t, _)| t.starts_with("A line that says"))
            .expect("the line of purpose");
        assert!(
            purpose.max.x <= right + 1.0,
            "{purpose:?} in a room that ends at {right}"
        );
        assert!(purpose.height() < 26.0, "one row: {purpose:?}");
        // every button is under the line, inside the room, and they are on more rows than one
        assert_eq!(h.buttons.len(), ACTIONS.len());
        for button in &h.buttons {
            assert!(
                button.min.y >= purpose.max.y - 1.0,
                "{button:?} is not under {purpose:?}"
            );
            assert!(
                button.min.x >= left - 1.0 && button.max.x <= right + 1.0,
                "{button:?} is out of the room from {left} to {right}"
            );
        }
        let mut rows: Vec<i32> = h.buttons.iter().map(|b| b.min.y.round() as i32).collect();
        rows.sort_unstable();
        rows.dedup();
        assert!(
            rows.len() > 1,
            "four buttons in one row of {wide} points: {:?}",
            h.buttons
        );
        // and the header is as high as it needs to be, no more: the page under it has the rest of the window
        let bottom = h.buttons.iter().map(|b| b.max.y).fold(0.0, f32::max);
        let drawn = bottom - title.min.y;
        assert!(
            h.high >= drawn && h.high < drawn + 20.0,
            "the header takes {} of its room, and what it draws ends {drawn} down",
            h.high
        );
        assert!(h.high < 160.0, "the header is {} high", h.high);
    }

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
    fn a_window_opens_where_it_was_left_never_above_the_floor_and_every_window_drawn_is_noted() {
        let ctx = egui::Context::default();
        let draw = |n: u32, titles: &[&str]| {
            let mut input = crate::ui::headless::frame_input(n);
            input.events = vec![];
            let mut out = ctx.run_ui(input, |ui| {
                for t in titles {
                    window(*t).show(ui.ctx(), |ui| {
                        ui.label("x");
                    });
                }
            });
            out.textures_delta.clear();
        };
        let corner = |t: &str| {
            ctx.memory(|m| m.area_rect(window_id(t)))
                .expect("the window was drawn")
                .min
        };
        set_saved_places(
            [
                ("Find".to_string(), egui::pos2(300.0, 220.0)),
                ("Low".to_string(), egui::pos2(10.0, 5.0)),
            ]
            .into(),
        );
        draw(0, &["Find", "Low", "Other"]);
        draw(1, &["Find", "Low", "Other"]);
        assert_eq!(corner("Find"), egui::pos2(300.0, 220.0));
        // a place above or left of the floor is brought to it
        assert_eq!(
            corner("Low"),
            egui::pos2(crate::ui::places::FLOOR.0, crate::ui::places::FLOOR.1)
        );
        // a window with no place kept opens where it always did
        assert_eq!(corner("Other"), default_place("Other"));
        // every window that was drawn is noted, once, and the notes are taken
        let mut seen = take_seen_windows();
        seen.sort();
        assert_eq!(seen, ["Find", "Low", "Other"].map(String::from));
        assert!(take_seen_windows().is_empty());
        set_saved_places(Default::default());
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
            let l = look(dark);
            assert!(b"ACGTU".iter().all(|&b| base_fill(b, &l).is_some()));
            assert!(base_fill(b'-', &l).is_none() && base_fill(b'N', &l).is_none());
            assert!(b"ACDEFGHIKLMNPQRSTVWY"
                .iter()
                .all(|&b| residue_fill(b, &l).is_some()));
            assert!(residue_fill(b'-', &l).is_none() && residue_fill(b'X', &l).is_none());
        }
        // lower case reads like upper case
        assert_eq!(base_fill(b'a', &look(false)), base_fill(b'A', &look(false)));
        assert_eq!(
            residue_fill(b'k', &look(true)),
            residue_fill(b'K', &look(true))
        );
    }

    #[test]
    fn the_four_trace_colours_differ() {
        let (light, dark) = (look(false), look(true));
        let colours: std::collections::HashSet<_> =
            b"ACGT".iter().map(|&b| trace_color(b, &light)).collect();
        assert_eq!(colours.len(), 4);
        assert_ne!(trace_color(b'G', &light), trace_color(b'G', &dark));
    }

    #[test]
    fn n_the_ambiguity_codes_and_the_gap_are_all_magenta_as_in_the_original() {
        let l = look(false);
        let magenta = trace_color(b'N', &l);
        for b in *b"BDHKMRSVWYU-n*" {
            assert_eq!(trace_color(b, &l), magenta, "{}", char::from(b));
        }
        assert!(!b"ACGT".iter().any(|&b| trace_color(b, &l) == magenta));
        // lower case reads like upper case
        assert_eq!(trace_color(b'a', &l), trace_color(b'A', &l));
    }
}
