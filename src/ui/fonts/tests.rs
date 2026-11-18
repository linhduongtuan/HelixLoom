use super::*;
use crate::ui::headless::frame_input;
use eframe::egui::{TextStyle, TextWrapMode, WidgetText};

fn axis(bytes: &'static [u8]) -> egui::epaint::text::FontVariationAxis {
    FontData::from_static(bytes)
        .variation_axes()
        .into_iter()
        .find(|a| &a.tag.to_be_bytes() == b"wght")
        .expect("the font has a weight axis")
}

#[test]
fn both_fonts_have_a_weight_axis_that_reaches_semibold() {
    for (name, bytes) in [("Inter", INTER), ("JetBrains Mono", MONO)] {
        let a = axis(bytes);
        assert!(
            a.range.min <= 400.0 && a.range.max >= SEMIBOLD,
            "{name}: {:?}",
            a.range
        );
    }
}

#[test]
fn the_two_fonts_lead_their_families_and_egui_s_own_stay_behind_for_other_signs() {
    let d = definitions();
    let first = |f: FontFamily| d.families[&f].first().cloned();
    assert_eq!(first(FontFamily::Proportional).as_deref(), Some("Inter"));
    assert_eq!(
        first(FontFamily::Monospace).as_deref(),
        Some("JetBrains Mono")
    );
    assert!(d.families[&FontFamily::Proportional].len() > 1);
    assert!(d.families[&FontFamily::Monospace].len() > 1);
}

/// The width of `text` as the window draws it, in the first frame after the fonts are in.
fn width(text: impl Fn() -> RichText) -> f32 {
    let ctx = egui::Context::default();
    install(&ctx);
    // the fonts are taken up as a frame begins
    ctx.run_ui(frame_input(0), |_| {}).textures_delta.clear();
    let mut width = 0.0;
    let mut out = ctx.run_ui(frame_input(1), |ui| {
        let galley = WidgetText::from(text()).into_galley(
            ui,
            Some(TextWrapMode::Extend),
            f32::INFINITY,
            TextStyle::Body,
        );
        width = galley.size().x;
    });
    out.textures_delta.clear();
    width
}

#[test]
fn semibold_text_is_drawn_wider_than_regular_text_of_the_same_words() {
    let words = "Hamburgefonstiv Multiple alignment";
    let regular = width(|| RichText::new(words));
    let semibold = width(|| RichText::new(words).semibold());
    assert!(regular > 50.0, "text has a width: {regular}");
    assert!(
        semibold > regular * 1.01,
        "semibold {semibold} against regular {regular}"
    );
}

#[test]
fn the_sequence_font_has_one_width_for_every_base() {
    let one = |c: char| width(|| RichText::new(c.to_string()).monospace());
    let a = one('A');
    for c in ['C', 'G', 'T', 'N', '-'] {
        assert!((one(c) - a).abs() < 0.01, "{c}");
    }
}
