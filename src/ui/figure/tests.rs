use super::*;
use crate::ui::fonts::FontWeight;
use eframe::egui::epaint::text::{LayoutJob, TextFormat};
use eframe::egui::epaint::TextShape;
use eframe::egui::{pos2, vec2, Align2, FontId, RichText, Stroke};

const SIZE: Vec2 = vec2(300.0, 120.0);
const RED: Color32 = Color32::from_rgb(255, 0, 0);
const BLUE: Color32 = Color32::from_rgb(0, 128, 255);

fn looks() -> Looks {
    Looks {
        dark: false,
        appearance: Appearance::default(),
    }
}

/// The lines of the SVG of what `add` draws, without its first line, its ground and its last line.
fn drawn(add: impl FnMut(&mut Ui)) -> Vec<String> {
    let svg = draw(looks(), SIZE, 1.0, add).svg();
    let lines: Vec<String> = svg.lines().map(str::to_string).collect();
    assert!(
        lines[0].contains("width=\"300\" height=\"120\" viewBox=\"0 0 300 120\""),
        "{}",
        lines[0]
    );
    assert!(
        lines[1].starts_with("<rect width=\"300\" height=\"120\" fill="),
        "{}",
        lines[1]
    );
    assert_eq!(lines.last().map(String::as_str), Some("</svg>"));
    lines[2..lines.len() - 1].to_vec()
}

/// The number that `attribute` of an element has.
fn number(element: &str, attribute: &str) -> f32 {
    let from = element
        .find(&format!(" {attribute}=\""))
        .unwrap_or_else(|| panic!("no {attribute} in {element}"))
        + attribute.len()
        + 3;
    let to = from + element[from..].find('"').expect("the end of the value");
    element[from..to].parse().expect("a number")
}

#[test]
fn a_figure_holds_what_the_view_draws_and_nothing_else_and_is_the_same_every_time() {
    let view = |ui: &mut Ui| {
        ui.painter().rect_filled(
            Rect::from_min_size(pos2(10.0, 10.0), vec2(20.0, 20.0)),
            0.0,
            RED,
        );
        ui.painter()
            .line_segment([pos2(0.0, 50.0), pos2(300.0, 50.0)], Stroke::new(1.0, BLUE));
    };
    let lines = drawn(view);
    assert_eq!(
        lines,
        [
            "<rect x=\"10\" y=\"10\" width=\"20\" height=\"20\" fill=\"#ff0000\"/>",
            "<line x1=\"0\" y1=\"50\" x2=\"300\" y2=\"50\" stroke=\"#0080ff\" stroke-width=\"1\"/>",
        ]
    );
    assert_eq!(
        draw(looks(), SIZE, 1.0, view).svg(),
        draw(looks(), SIZE, 1.0, view).svg()
    );
}

#[test]
fn text_stays_text_with_its_font_its_size_its_colour_and_the_width_it_has_on_screen() {
    let lines = drawn(|ui| {
        ui.painter().text(
            pos2(20.0, 30.0),
            Align2::LEFT_TOP,
            "ACGT <&>",
            FontId::monospace(14.0),
            Color32::from_rgb(10, 20, 30),
        );
    });
    assert_eq!(lines.len(), 1, "{lines:?}");
    let text = &lines[0];
    assert!(
        text.starts_with("<text x=\"20\" y=\"") && text.ends_with(">ACGT &lt;&amp;&gt;</text>"),
        "{text}"
    );
    assert!(
        text.contains(" font-family=\"'JetBrains Mono', Menlo, Consolas, monospace\" font-size=\"14\" fill=\"#0a141e\""),
        "{text}"
    );
    // the baseline is under the top of the text by less than its size, and eight letters of one width are that wide
    let baseline = number(text, "y");
    assert!(baseline > 38.0 && baseline < 46.0, "{text}");
    let width = number(text, "textLength");
    assert!(width > 8.0 * 7.0 && width < 8.0 * 10.0, "{text}");
}

#[test]
fn a_text_of_two_colours_is_two_stretches_side_by_side_and_a_heavy_text_says_its_weight() {
    let lines = drawn(|ui| {
        let mut job = LayoutJob::default();
        for (letters, colour) in [("AB", RED), ("CD", BLUE)] {
            job.append(
                letters,
                0.0,
                TextFormat {
                    font_id: FontId::proportional(20.0),
                    color: colour,
                    ..Default::default()
                },
            );
        }
        let galley = ui.painter().layout_job(job);
        ui.painter()
            .galley(pos2(10.0, 10.0), galley, Color32::BLACK);
        ui.put(
            Rect::from_min_size(pos2(10.0, 60.0), vec2(200.0, 30.0)),
            eframe::egui::Label::new(RichText::new("Heavy").semibold()),
        );
    });
    assert_eq!(lines.len(), 3, "{lines:?}");
    assert!(
        lines[0].contains("fill=\"#ff0000\"") && lines[0].ends_with(">AB</text>"),
        "{}",
        lines[0]
    );
    assert!(
        lines[1].contains("fill=\"#0080ff\"") && lines[1].ends_with(">CD</text>"),
        "{}",
        lines[1]
    );
    assert!(
        lines[0].contains(" font-family=\"Inter, 'Helvetica Neue', Helvetica, Arial, sans-serif\" font-size=\"20\""),
        "{}",
        lines[0]
    );
    // the second starts where the first ends (to the pixel that a letter is put on), on the same line
    let end = number(&lines[0], "x") + number(&lines[0], "textLength");
    assert!((number(&lines[1], "x") - end).abs() < 1.0, "{lines:?}");
    assert_eq!(number(&lines[0], "y"), number(&lines[1], "y"));
    assert!(
        lines[2].contains(" font-weight=\"600\"") && lines[2].ends_with(">Heavy</text>"),
        "{}",
        lines[2]
    );
    assert!(!lines[0].contains("font-weight"), "{}", lines[0]);
    // a label that takes its colour from the window has the window's colour for text
    let [r, g, b, _] = looks().appearance.look(false).text.to_srgba_unmultiplied();
    assert!(
        lines[2].contains(&format!(" fill=\"#{r:02x}{g:02x}{b:02x}\" ")),
        "{}",
        lines[2]
    );
}

#[test]
fn a_figure_has_the_sizes_of_the_window_at_its_plain_size_however_far_the_window_is_zoomed() {
    let view = |ui: &mut Ui| {
        ui.label("A label of the window");
        ui.painter().circle_filled(pos2(100.0, 80.0), 12.0, BLUE);
    };
    let zoomed = Looks {
        dark: false,
        appearance: Appearance {
            scale: 150,
            ..Appearance::default()
        },
    };
    assert_eq!(
        draw(zoomed, SIZE, 1.0, view).svg(),
        draw(looks(), SIZE, 1.0, view).svg()
    );
}

#[test]
fn a_turned_text_is_turned_about_the_place_it_was_put_at_and_spaces_alone_draw_nothing() {
    let lines = drawn(|ui| {
        let galley = ui
            .painter()
            .layout_no_wrap("up".to_string(), FontId::proportional(12.0), RED);
        ui.painter().add(
            TextShape::new(pos2(100.0, 80.0), galley, RED).with_angle(-std::f32::consts::FRAC_PI_2),
        );
        ui.painter().text(
            pos2(5.0, 5.0),
            Align2::LEFT_TOP,
            "   ",
            FontId::monospace(12.0),
            RED,
        );
    });
    assert_eq!(lines.len(), 3, "{lines:?}");
    assert_eq!(lines[0], "<g transform=\"rotate(-90 100 80)\">");
    assert!(
        lines[1].starts_with("<text x=\"100\"") && lines[1].ends_with(">up</text>"),
        "{}",
        lines[1]
    );
    assert_eq!(lines[2], "</g>");
}

#[test]
fn a_text_with_no_colour_of_its_own_takes_the_colour_it_is_drawn_with() {
    let lines = drawn(|ui| {
        // egui lays out a label's text with a colour that stands for "the one I am painted with"
        let galley = ui.painter().layout_no_wrap(
            "plain".to_string(),
            FontId::proportional(12.0),
            Color32::PLACEHOLDER,
        );
        ui.painter()
            .add(TextShape::new(pos2(10.0, 10.0), galley, RED));
    });
    assert_eq!(lines.len(), 1, "{lines:?}");
    assert!(lines[0].contains(" fill=\"#ff0000\" "), "{}", lines[0]);
    assert!(!lines[0].contains("fill-opacity"), "{}", lines[0]);
}

#[test]
fn a_view_that_measures_itself_in_one_frame_and_lays_itself_out_in_the_next_is_drawn_laid_out() {
    // the columns of a grid are as wide as their widest cell, which the grid knows from the frame before
    let lines = drawn(|ui| {
        eframe::egui::Grid::new("columns").show(ui, |ui| {
            ui.label("a");
            ui.label("first");
            ui.end_row();
            ui.label("a much longer cell");
            ui.label("second");
            ui.end_row();
        });
    });
    let x_of = |word: &str| {
        let line = lines
            .iter()
            .find(|l| l.ends_with(&format!(">{word}</text>")))
            .unwrap_or_else(|| panic!("{word} in {lines:?}"));
        number(line, "x")
    };
    assert_eq!(x_of("first"), x_of("second"), "{lines:?}");
    assert!(
        x_of("first") > x_of("a much longer cell") + 60.0,
        "{lines:?}"
    );
}

#[test]
fn a_png_has_the_points_times_the_scale_in_pixels_with_the_colours_of_the_view_on_the_ground() {
    let figure = draw(looks(), vec2(100.0, 60.0), 4.0, |ui| {
        ui.painter().rect_filled(
            Rect::from_min_size(pos2(10.0, 10.0), vec2(20.0, 20.0)),
            0.0,
            RED,
        );
    });
    let image = figure.image();
    assert_eq!(image.size, [400, 240]);
    let at = |x: usize, y: usize| image.pixels[y * 400 + x];
    assert_eq!(at(80, 80), RED, "the middle of the box");
    let ground = looks().appearance.look(false).bg;
    assert_eq!(at(300, 200), ground, "beside the box");
    assert_eq!(at(20, 20), ground, "before the box begins at 40 pixels");
    // and the bytes are a PNG of that size
    let bytes = headless::png_bytes(&image).expect("a PNG");
    assert_eq!(&bytes[..8], b"\x89PNG\r\n\x1a\n");
    assert_eq!(&bytes[16..24], [0, 0, 1, 144, 0, 0, 0, 240], "400 by 240");
}

#[test]
fn a_figure_is_never_smaller_or_larger_than_it_can_be_and_a_large_one_has_fewer_pixels_to_a_point()
{
    assert_eq!(room(vec2(0.0, -5.0)), vec2(16.0, 16.0));
    assert_eq!(room(vec2(f32::NAN, f32::INFINITY)), vec2(16.0, 16.0));
    assert_eq!(room(vec2(800.4, 99_999.0)), vec2(800.0, 6000.0));
    assert_eq!(png_scale(vec2(800.0, 600.0)), 4.0);
    // 3000 by 3000 points at four pixels each would be 144 million pixels
    let scale = png_scale(vec2(3000.0, 3000.0));
    assert!(scale > 1.6 && scale < 1.7, "{scale}");
    assert!(3000.0 * 3000.0 * scale * scale <= MOST_PIXELS * 1.001);
    // and never less than one pixel to a point
    assert_eq!(png_scale(vec2(6000.0, 6000.0)), 1.0);
}

#[test]
fn saving_draws_the_view_as_the_kind_asked_for_under_its_name() {
    let view = |ui: &mut Ui| {
        ui.painter().circle_filled(pos2(30.0, 30.0), 10.0, BLUE);
        ui.painter().text(
            pos2(20.0, 60.0),
            Align2::LEFT_TOP,
            "ACGTACGTAC",
            FontId::monospace(13.0),
            BLUE,
        );
    };
    let _ = take_saved();
    assert_eq!(save(Format::Svg, "map", SIZE, view), "Saved map.svg.");
    assert_eq!(save(Format::Png, "map", SIZE, view), "Saved map.png.");
    let saved = take_saved();
    assert_eq!(saved.len(), 2);
    assert_eq!(saved[0].0, "map.svg");
    let svg = String::from_utf8(saved[0].1.clone()).expect("text");
    assert!(
        svg.contains("<circle cx=\"30\" cy=\"30\" r=\"10\" fill=\"#0080ff\"/>"),
        "{svg}"
    );
    // letters are as wide as they are, not put on whole points: ten letters of 7.8 points are 78 wide, not 80
    let text = svg
        .lines()
        .find(|l| l.ends_with(">ACGTACGTAC</text>"))
        .expect("the letters");
    let width = number(text, "textLength");
    assert!((width - 78.0).abs() < 0.6, "{text}");
    assert_eq!(saved[1].0, "map.png");
    assert_eq!(&saved[1].1[..8], b"\x89PNG\r\n\x1a\n");
    // 300 by 120 points, four pixels each
    assert_eq!(&saved[1].1[16..24], [0, 0, 4, 176, 0, 0, 1, 224]);
    assert!(take_saved().is_empty());
}

#[test]
fn a_png_keeps_a_colour_and_how_see_through_it_is_apart() {
    // a picture keeps them multiplied together: red at half is (128, 0, 0, 128) there, and (255, 0, 0, 128) in a PNG
    let half_red = Color32::from_rgba_unmultiplied(255, 0, 0, 128);
    let image = ColorImage {
        size: [2, 1],
        source_size: vec2(2.0, 1.0),
        pixels: vec![half_red, Color32::TRANSPARENT],
    };
    let bytes = headless::png_bytes(&image).expect("a PNG");
    let mut reader = png::Decoder::new(std::io::Cursor::new(bytes))
        .read_info()
        .expect("a PNG that reads");
    let mut pixels = vec![0; reader.output_buffer_size().expect("its size")];
    reader.next_frame(&mut pixels).expect("its pixels");
    assert_eq!(pixels[..4], [255, 0, 0, 128]);
    assert_eq!(pixels[7], 0, "nothing where there is nothing");
}
