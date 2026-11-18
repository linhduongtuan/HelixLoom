use super::*;
use eframe::egui::epaint::{
    CornerRadius, CubicBezierShape, EllipseShape, QuadraticBezierShape, RectShape,
};
use eframe::egui::{pos2, vec2};

const SIZE: Vec2 = vec2(200.0, 100.0);
const RED: Color32 = Color32::from_rgb(255, 0, 0);
const BLUE: Color32 = Color32::from_rgb(0, 128, 255);

fn whole() -> Rect {
    Rect::from_min_size(Pos2::ZERO, SIZE)
}

/// The lines of the picture of `shapes`, each cut to `cut`, between its first line and its last.
fn drawn_in(shapes: Vec<(Rect, Shape)>) -> Vec<String> {
    let clipped: Vec<ClippedShape> = shapes
        .into_iter()
        .map(|(clip_rect, shape)| ClippedShape { clip_rect, shape })
        .collect();
    let svg = to_svg(&clipped, SIZE, Color32::TRANSPARENT);
    let lines: Vec<String> = svg.lines().map(str::to_string).collect();
    assert_eq!(
        lines[0],
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"200\" height=\"100\" viewBox=\"0 0 200 100\" xml:space=\"preserve\">"
    );
    assert_eq!(lines.last().map(String::as_str), Some("</svg>"));
    lines[1..lines.len() - 1].to_vec()
}

fn drawn(shapes: Vec<Shape>) -> Vec<String> {
    drawn_in(shapes.into_iter().map(|s| (whole(), s)).collect())
}

#[test]
fn a_number_has_two_decimals_at_most_and_no_zeros_after_them() {
    assert_eq!(n(10.0), "10");
    assert_eq!(n(0.5), "0.5");
    assert_eq!(n(1.256), "1.26");
    assert_eq!(n(-2.5), "-2.5");
    assert_eq!(n(-0.001), "0");
    assert_eq!(n(100.0), "100");
    // what is no number would make a file that nothing opens
    assert_eq!(n(f32::NAN), "0");
    assert_eq!(n(f32::INFINITY), "0");
    assert_eq!(n(f32::NEG_INFINITY), "0");
}

#[test]
fn a_shape_with_a_place_that_is_no_number_still_makes_a_file_of_numbers() {
    let got = drawn(vec![Shape::circle_filled(
        pos2(f32::NAN, 5.0),
        f32::INFINITY,
        RED,
    )]);
    assert_eq!(
        got,
        ["<circle cx=\"0\" cy=\"5\" r=\"0\" fill=\"#ff0000\"/>"]
    );
}

#[test]
fn characters_that_a_file_may_not_hold_are_left_out_of_its_text_and_the_signs_of_markup_are_written_out(
) {
    assert_eq!(escaped("a\u{1}b\tc\u{1b}d\u{0}\u{7f}\u{FFFF}"), "ab\tcd");
    assert_eq!(escaped("<a & b>"), "&lt;a &amp; b&gt;");
    assert_eq!(escaped("plain é 中"), "plain é 中");
}

#[test]
fn a_glyph_that_draws_no_character_of_the_text_keeps_its_neighbours_in_their_sections() {
    use eframe::egui::epaint::text::{LayoutJob, TextFormat};
    let mut job = LayoutJob::default();
    job.append("ab", 0.0, TextFormat::default());
    // another colour: two stretches of one format are one section
    job.append(
        "cd",
        0.0,
        TextFormat {
            color: RED,
            ..Default::default()
        },
    );
    assert_eq!(job.sections.len(), 2);
    let mut characters = job.text.char_indices();
    assert_eq!(section_of('a', &mut characters, &job, 9), 0);
    // the dots of a text cut short are no character of it: the section before stays, and so does the place in the text
    assert_eq!(section_of('…', &mut characters, &job, 0), 0);
    assert_eq!(section_of('b', &mut characters, &job, 0), 0);
    assert_eq!(section_of('c', &mut characters, &job, 0), 1);
    assert_eq!(section_of('d', &mut characters, &job, 1), 1);
    // nothing is left of the text
    assert_eq!(section_of('d', &mut characters, &job, 1), 1);

    // a letter that comes again is looked for after the one before it, not from the start
    let mut twice = LayoutJob::default();
    twice.append("ab", 0.0, TextFormat::default());
    twice.append(
        "ab",
        0.0,
        TextFormat {
            color: RED,
            ..Default::default()
        },
    );
    let mut characters = twice.text.char_indices();
    let sections: Vec<usize> = "abab"
        .chars()
        .map(|c| section_of(c, &mut characters, &twice, 0))
        .collect();
    assert_eq!(sections, [0, 0, 1, 1]);
}

#[test]
fn a_line_a_circle_and_an_ellipse_are_the_elements_that_draw_them() {
    let stroke = Stroke::new(1.5, RED);
    let got = drawn(vec![
        Shape::line_segment([pos2(1.0, 2.0), pos2(30.5, 40.25)], stroke),
        Shape::circle_filled(pos2(50.0, 50.0), 10.0, BLUE),
        Shape::circle_stroke(pos2(50.0, 50.0), 10.0, stroke),
        Shape::Ellipse(EllipseShape {
            center: pos2(20.0, 30.0),
            radius: vec2(8.0, 4.0),
            fill: BLUE,
            stroke: Stroke::NONE,
            angle: std::f32::consts::FRAC_PI_2,
        }),
        // a line of no width, and one of no colour, draw nothing
        Shape::line_segment([pos2(0.0, 0.0), pos2(9.0, 9.0)], Stroke::new(0.0, RED)),
        Shape::line_segment(
            [pos2(0.0, 0.0), pos2(9.0, 9.0)],
            Stroke::new(2.0, Color32::TRANSPARENT),
        ),
        Shape::Noop,
    ]);
    assert_eq!(
        got,
        [
            "<line x1=\"1\" y1=\"2\" x2=\"30.5\" y2=\"40.25\" stroke=\"#ff0000\" stroke-width=\"1.5\"/>",
            "<circle cx=\"50\" cy=\"50\" r=\"10\" fill=\"#0080ff\"/>",
            "<circle cx=\"50\" cy=\"50\" r=\"10\" fill=\"none\" stroke=\"#ff0000\" stroke-width=\"1.5\"/>",
            "<ellipse cx=\"20\" cy=\"30\" rx=\"8\" ry=\"4\" fill=\"#0080ff\" transform=\"rotate(90 20 30)\"/>",
        ]
    );
}

#[test]
fn a_colour_that_is_see_through_says_how_much_and_shapes_in_a_group_come_in_their_order() {
    let half = Color32::from_rgba_unmultiplied(255, 0, 0, 128);
    let got = drawn(vec![Shape::Vec(vec![
        Shape::circle_filled(pos2(5.0, 5.0), 2.0, half),
        Shape::circle_filled(pos2(6.0, 6.0), 3.0, BLUE),
    ])]);
    assert_eq!(
        got,
        [
            "<circle cx=\"5\" cy=\"5\" r=\"2\" fill=\"#ff0000\" fill-opacity=\"0.5\"/>",
            "<circle cx=\"6\" cy=\"6\" r=\"3\" fill=\"#0080ff\"/>",
        ]
    );
}

#[test]
fn a_rectangle_keeps_its_round_corners_and_an_outline_beside_its_edge_is_moved_onto_it() {
    let rect = Rect::from_min_size(pos2(10.0, 20.0), vec2(40.0, 30.0));
    let line = Stroke::new(2.0, RED);
    let got = drawn(vec![
        Shape::rect_filled(rect, 0.0, BLUE),
        Shape::rect_filled(rect, 4.0, BLUE),
        Shape::Rect(RectShape::new(rect, 0.0, BLUE, line, StrokeKind::Inside)),
        Shape::Rect(RectShape::new(rect, 0.0, BLUE, line, StrokeKind::Outside)),
        Shape::Rect(RectShape::new(rect, 0.0, BLUE, line, StrokeKind::Middle)),
        Shape::Rect(RectShape::new(
            rect,
            CornerRadius {
                nw: 5,
                ne: 0,
                sw: 0,
                se: 5,
            },
            BLUE,
            Stroke::NONE,
            StrokeKind::Middle,
        )),
        // a soft shadow is not drawn as a hard box
        Shape::Rect(RectShape::filled(rect, 0.0, BLUE).with_blur_width(6.0)),
        // an outline that cannot be seen moves nothing
        Shape::Rect(RectShape::new(
            rect,
            0.0,
            BLUE,
            Stroke::new(2.0, Color32::TRANSPARENT),
            StrokeKind::Inside,
        )),
    ]);
    let outlined = "fill=\"#0080ff\" stroke=\"#ff0000\" stroke-width=\"2\"/>";
    assert_eq!(
        got,
        [
            "<rect x=\"10\" y=\"20\" width=\"40\" height=\"30\" fill=\"#0080ff\"/>".to_string(),
            "<rect x=\"10\" y=\"20\" width=\"40\" height=\"30\" rx=\"4\" fill=\"#0080ff\"/>".to_string(),
            format!("<rect x=\"11\" y=\"21\" width=\"38\" height=\"28\" {outlined}"),
            format!("<rect x=\"9\" y=\"19\" width=\"42\" height=\"32\" {outlined}"),
            format!("<rect x=\"10\" y=\"20\" width=\"40\" height=\"30\" {outlined}"),
            "<path d=\"M15 20H50A0 0 0 0 1 50 20V45A5 5 0 0 1 45 50H10A0 0 0 0 1 10 50V25A5 5 0 0 1 15 20Z\" fill=\"#0080ff\"/>"
                .to_string(),
            "<rect x=\"10\" y=\"20\" width=\"40\" height=\"30\" fill=\"#0080ff\"/>".to_string(),
        ]
    );
}

#[test]
fn an_open_path_is_a_line_through_its_points_and_a_closed_one_is_filled() {
    let points = vec![pos2(1.0, 1.0), pos2(20.0, 1.5), pos2(10.0, 30.0)];
    let got = drawn(vec![
        Shape::line(points.clone(), Stroke::new(1.0, RED)),
        Shape::convex_polygon(points.clone(), BLUE, Stroke::NONE),
        Shape::closed_line(points, Stroke::new(1.0, RED)),
        // one point is no path
        Shape::line(vec![pos2(3.0, 3.0)], Stroke::new(1.0, RED)),
    ]);
    assert_eq!(
        got,
        [
            "<polyline points=\"1,1 20,1.5 10,30\" fill=\"none\" stroke=\"#ff0000\" stroke-width=\"1\"/>",
            "<polygon points=\"1,1 20,1.5 10,30\" fill=\"#0080ff\"/>",
            "<polygon points=\"1,1 20,1.5 10,30\" fill=\"none\" stroke=\"#ff0000\" stroke-width=\"1\"/>",
        ]
    );
}

#[test]
fn curves_are_paths_by_their_control_points() {
    let line = Stroke::new(1.0, RED);
    let got = drawn(vec![
        Shape::QuadraticBezier(QuadraticBezierShape::from_points_stroke(
            [pos2(0.0, 0.0), pos2(10.0, 20.0), pos2(30.0, 0.0)],
            false,
            Color32::TRANSPARENT,
            line,
        )),
        Shape::CubicBezier(CubicBezierShape::from_points_stroke(
            [
                pos2(0.0, 0.0),
                pos2(10.0, 20.0),
                pos2(20.0, 20.0),
                pos2(30.0, 0.0),
            ],
            true,
            BLUE,
            line,
        )),
    ]);
    assert_eq!(
        got,
        [
            "<path d=\"M0 0 Q10 20 30 0\" fill=\"none\" stroke=\"#ff0000\" stroke-width=\"1\"/>",
            "<path d=\"M0 0 C10 20 20 20 30 0Z\" fill=\"#0080ff\" stroke=\"#ff0000\" stroke-width=\"1\"/>",
        ]
    );
}

#[test]
fn the_triangles_of_a_mesh_that_follow_each_other_in_one_colour_are_one_area() {
    let mut mesh = Mesh::default();
    for (colour, x) in [(RED, 0.0), (RED, 10.0), (BLUE, 20.0)] {
        let first = mesh.vertices.len() as u32;
        for corner in [pos2(x, 0.0), pos2(x + 10.0, 0.0), pos2(x, 10.0)] {
            mesh.colored_vertex(corner, colour);
        }
        mesh.add_triangle(first, first + 1, first + 2);
    }
    let seam = |colour: &str| {
        format!("stroke=\"{colour}\" stroke-width=\"0.3\" stroke-linejoin=\"round\"/>")
    };
    assert_eq!(
        drawn(vec![Shape::mesh(mesh.clone())]),
        [
            format!(
                "<path d=\"M0 0L10 0L0 10ZM10 0L20 0L10 10Z\" fill=\"#ff0000\" {}",
                seam("#ff0000")
            ),
            format!(
                "<path d=\"M20 0L30 0L20 10Z\" fill=\"#0080ff\" {}",
                seam("#0080ff")
            ),
        ]
    );
    // a mesh with a picture on it is left out
    mesh.texture_id = TextureId::User(7);
    assert!(drawn(vec![Shape::mesh(mesh)]).is_empty());
}

#[test]
fn a_shape_its_view_cut_to_a_part_of_the_room_is_cut_to_that_part() {
    let part = Rect::from_min_max(pos2(10.0, 10.0), pos2(50.0, 60.0));
    // a cut that reaches past the room is the part of it that is in the room
    let past = Rect::from_min_max(pos2(150.0, -20.0), pos2(400.0, 60.0));
    let dot = |x: f32| Shape::circle_filled(pos2(x, 5.0), 1.0, RED);
    let got = drawn_in(vec![
        (part, dot(1.0)),
        (part, dot(2.0)),
        (whole(), dot(3.0)),
        (past, dot(4.0)),
        (part, dot(5.0)),
        // nothing drawn in a cut makes no group
        (past, Shape::Noop),
        // a cut that holds the whole room cuts nothing
        (whole().expand(30.0), dot(6.0)),
    ]);
    let circle = |x: u8| format!("<circle cx=\"{x}\" cy=\"5\" r=\"1\" fill=\"#ff0000\"/>");
    // the cuts are named after the figure: `cut`, ten letters of its checksum, and their number
    let mark = got[1]
        .strip_prefix("<clipPath id=\"cut")
        .and_then(|rest| rest.split_once("-0\">"))
        .map(|(mark, _)| mark.to_string())
        .unwrap_or_else(|| panic!("no mark in {}", got[1]));
    assert_eq!(mark.len(), 10, "{mark}");
    assert!(mark.chars().all(|c| c.is_ascii_hexdigit()), "{mark}");
    assert_eq!(
        got,
        [
            "<defs>".to_string(),
            format!(
                "<clipPath id=\"cut{mark}-0\"><rect x=\"10\" y=\"10\" width=\"40\" height=\"50\"/></clipPath>"
            ),
            format!(
                "<clipPath id=\"cut{mark}-1\"><rect x=\"150\" y=\"0\" width=\"50\" height=\"60\"/></clipPath>"
            ),
            "</defs>".to_string(),
            format!("<g clip-path=\"url(#cut{mark}-0)\">"),
            circle(1),
            circle(2),
            "</g>".to_string(),
            circle(3),
            format!("<g clip-path=\"url(#cut{mark}-1)\">"),
            circle(4),
            "</g>".to_string(),
            format!("<g clip-path=\"url(#cut{mark}-0)\">"),
            circle(5),
            "</g>".to_string(),
            circle(6),
        ]
    );
}

#[test]
fn the_cuts_of_two_figures_have_different_names_so_that_a_page_can_hold_both() {
    let part = Rect::from_min_max(pos2(10.0, 10.0), pos2(50.0, 60.0));
    let dot = |x: f32| Shape::circle_filled(pos2(x, 5.0), 1.0, RED);
    let names = |shapes: Vec<(Rect, Shape)>| -> Vec<String> {
        drawn_in(shapes)
            .iter()
            .filter_map(|line| line.strip_prefix("<clipPath id=\""))
            .filter_map(|rest| rest.split('"').next())
            .map(str::to_string)
            .collect()
    };
    // the same cut, and another thing drawn in it: another name
    let one = names(vec![(part, dot(1.0))]);
    let other = names(vec![(part, dot(2.0))]);
    assert_eq!((one.len(), other.len()), (1, 1));
    assert_ne!(one, other);
    // the same figure made again has the same name: a figure is the same file every time it is saved
    assert_eq!(one, names(vec![(part, dot(1.0))]));
    // a figure that cuts nothing has no names at all
    assert!(names(vec![(whole(), dot(1.0))]).is_empty());
}

#[test]
fn the_ground_is_a_box_as_large_as_the_picture_under_everything() {
    let dot = ClippedShape {
        clip_rect: whole(),
        shape: Shape::circle_filled(pos2(5.0, 5.0), 1.0, RED),
    };
    let svg = to_svg(&[dot], SIZE, Color32::WHITE);
    let lines: Vec<&str> = svg.lines().collect();
    assert_eq!(
        lines[1],
        "<rect width=\"200\" height=\"100\" fill=\"#ffffff\"/>"
    );
    assert!(lines[2].starts_with("<circle"));
    assert_eq!(lines.len(), 4);
}
