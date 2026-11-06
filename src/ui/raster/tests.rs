use super::*;
use eframe::egui::epaint::{ImageDelta, Mesh};
use eframe::egui::{pos2, TextureOptions};

/// What egui makes of a texture it has just created.
fn created(image: ColorImage, options: TextureOptions) -> TexturesDelta {
    let mut delta = TexturesDelta::default();
    delta.push(TextureId::default(), ImageDelta::full(image, options));
    delta
}

/// Texture 0, where egui keeps its glyphs and the white texel that plain fills sample.
fn white_texture() -> Textures {
    let mut textures = Textures::default();
    textures.apply(&mut created(
        ColorImage::filled([1, 1], Color32::WHITE),
        TextureOptions::NEAREST,
    ));
    textures
}

fn rect(left: f32, top: f32, right: f32, bottom: f32) -> Rect {
    Rect::from_min_max(pos2(left, top), pos2(right, bottom))
}

fn filled(rect: Rect, color: Color32) -> Mesh {
    let mut mesh = Mesh::default();
    mesh.add_colored_rect(rect, color);
    mesh
}

/// `mesh` drawn on black, `size` pixels across, clipped to `clip`.
fn paint(mesh: Mesh, clip: Rect, size: [usize; 2], scale: f32, textures: &Textures) -> ColorImage {
    let mut canvas = Canvas::new(size, Color32::BLACK);
    let job = ClippedPrimitive {
        clip_rect: clip,
        primitive: Primitive::Mesh(mesh),
    };
    canvas.paint(&[job], textures, scale);
    canvas.into_image()
}

fn pixel(image: &ColorImage, x: usize, y: usize) -> Color32 {
    image.pixels[y * image.size[0] + x]
}

fn everywhere() -> Rect {
    rect(0.0, 0.0, 100.0, 100.0)
}

#[test]
fn a_rectangle_fills_exactly_the_pixels_it_covers() {
    let mesh = filled(rect(2.0, 3.0, 7.0, 9.0), Color32::RED);
    let image = paint(mesh, everywhere(), [10, 12], 1.0, &white_texture());
    for y in 0..12 {
        for x in 0..10 {
            let covered = (2..7).contains(&x) && (3..9).contains(&y);
            let expected = if covered {
                Color32::RED
            } else {
                Color32::BLACK
            };
            assert_eq!(pixel(&image, x, y), expected, "pixel ({x}, {y})");
        }
    }
}

#[test]
fn a_see_through_rectangle_is_blended_once_everywhere_even_along_the_diagonal() {
    let half_white = Color32::from_rgba_premultiplied(128, 128, 128, 128);
    let mesh = filled(rect(0.0, 0.0, 8.0, 8.0), half_white);
    let image = paint(mesh, everywhere(), [8, 8], 1.0, &white_texture());
    // twice over would give 191 where the two triangles meet
    assert!(image
        .pixels
        .iter()
        .all(|&p| p == Color32::from_rgb(128, 128, 128)));
}

#[test]
fn a_triangle_that_goes_the_other_way_round_is_drawn_the_same() {
    let forward = filled(rect(1.0, 1.0, 5.0, 4.0), Color32::GREEN);
    let mut backward = forward.clone();
    backward.indices.reverse();
    let textures = white_texture();
    assert_eq!(
        paint(forward, everywhere(), [6, 5], 1.0, &textures).pixels,
        paint(backward, everywhere(), [6, 5], 1.0, &textures).pixels
    );
}

#[test]
fn drawing_stays_inside_the_clip_rectangle() {
    let mesh = filled(rect(0.0, 0.0, 10.0, 4.0), Color32::RED);
    let image = paint(
        mesh,
        rect(3.0, 1.0, 6.0, 3.0),
        [10, 4],
        1.0,
        &white_texture(),
    );
    for y in 0..4 {
        for x in 0..10 {
            let shown = (3..6).contains(&x) && (1..3).contains(&y);
            let expected = if shown { Color32::RED } else { Color32::BLACK };
            assert_eq!(pixel(&image, x, y), expected, "pixel ({x}, {y})");
        }
    }
}

#[test]
fn points_become_pixels_at_the_given_scale() {
    // one point is two pixels: (1, 1) to (3, 2) points is pixels 2..6 across and 2..4 down
    let mesh = filled(rect(1.0, 1.0, 3.0, 2.0), Color32::BLUE);
    let image = paint(mesh, everywhere(), [8, 6], 2.0, &white_texture());
    let blue = image.pixels.iter().filter(|&&p| p == Color32::BLUE).count();
    assert_eq!(blue, 4 * 2);
    assert_eq!(pixel(&image, 2, 2), Color32::BLUE);
    assert_eq!(pixel(&image, 5, 3), Color32::BLUE);
    assert_eq!(pixel(&image, 6, 3), Color32::BLACK);
}

#[test]
fn a_texture_is_looked_up_pixel_by_pixel_with_either_filter() {
    for options in [TextureOptions::NEAREST, TextureOptions::LINEAR] {
        let mut textures = Textures::default();
        textures.apply(&mut created(
            ColorImage::new([2, 1], vec![Color32::RED, Color32::BLUE]),
            options,
        ));
        let mut mesh = Mesh::default();
        mesh.add_rect_with_uv(
            rect(0.0, 0.0, 2.0, 1.0),
            rect(0.0, 0.0, 1.0, 1.0),
            Color32::WHITE,
        );
        let image = paint(mesh, everywhere(), [2, 1], 1.0, &textures);
        assert_eq!(pixel(&image, 0, 0), Color32::RED, "{options:?}");
        assert_eq!(pixel(&image, 1, 0), Color32::BLUE, "{options:?}");
    }
}

#[test]
fn a_smooth_texture_blends_between_texels_when_magnified() {
    let mut textures = Textures::default();
    let black_and_white = ColorImage::new([2, 1], vec![Color32::BLACK, Color32::WHITE]);
    textures.apply(&mut created(black_and_white, TextureOptions::LINEAR));
    let mut mesh = Mesh::default();
    // four pixels across for two texels: the ones in between are grey
    mesh.add_rect_with_uv(
        rect(0.0, 0.0, 4.0, 1.0),
        rect(0.0, 0.0, 1.0, 1.0),
        Color32::WHITE,
    );
    let image = paint(mesh, everywhere(), [4, 1], 1.0, &textures);
    let reds: Vec<u8> = (0..4).map(|x| pixel(&image, x, 0).r()).collect();
    assert_eq!(reds[0], 0);
    assert_eq!(reds[3], 255);
    assert!(
        reds[1] > 0 && reds[1] < reds[2] && reds[2] < 255,
        "{reds:?}"
    );
}

#[test]
fn a_patch_replaces_only_its_own_part_of_a_texture() {
    let mut textures = Textures::default();
    textures.apply(&mut created(
        ColorImage::filled([2, 1], Color32::WHITE),
        TextureOptions::NEAREST,
    ));
    let mut patch = TexturesDelta::default();
    patch.push(
        TextureId::default(),
        ImageDelta::partial(
            [1, 0],
            ColorImage::filled([1, 1], Color32::RED),
            TextureOptions::NEAREST,
        ),
    );
    textures.apply(&mut patch);
    let mut mesh = Mesh::default();
    mesh.add_rect_with_uv(
        rect(0.0, 0.0, 2.0, 1.0),
        rect(0.0, 0.0, 1.0, 1.0),
        Color32::WHITE,
    );
    let image = paint(mesh, everywhere(), [2, 1], 1.0, &textures);
    assert_eq!(pixel(&image, 0, 0), Color32::WHITE);
    assert_eq!(pixel(&image, 1, 0), Color32::RED);
}

#[test]
fn a_patch_that_hangs_over_the_edge_is_cut_off_instead_of_failing() {
    let mut textures = white_texture();
    let mut patch = TexturesDelta::default();
    for at in [[0, 0], [1, 0], [0, 1], [5, 5]] {
        patch.push(
            TextureId::default(),
            ImageDelta::partial(
                at,
                ColorImage::filled([3, 3], Color32::RED),
                TextureOptions::NEAREST,
            ),
        );
    }
    textures.apply(&mut patch);
    assert!(patch.is_empty());
}

#[test]
fn freed_textures_stay_until_they_are_forgotten() {
    let mut textures = white_texture();
    let mut delta = TexturesDelta::default();
    delta.free(TextureId::default());
    let freed = textures.apply(&mut delta);
    assert_eq!(freed, vec![TextureId::default()]);
    let mesh = filled(rect(0.0, 0.0, 1.0, 1.0), Color32::RED);
    let shown = paint(mesh.clone(), everywhere(), [1, 1], 1.0, &textures);
    assert_eq!(shown.pixels, vec![Color32::RED]);
    textures.forget(&freed);
    let gone = paint(mesh, everywhere(), [1, 1], 1.0, &textures);
    assert_eq!(gone.pixels, vec![Color32::BLACK]);
}

#[test]
fn broken_triangles_draw_nothing_and_do_not_fail() {
    let mut mesh = filled(rect(1.0, 1.0, 3.0, 3.0), Color32::RED);
    mesh.vertices[0].pos = pos2(f32::NAN, 1.0);
    mesh.vertices[1].pos = pos2(f32::INFINITY, 1.0);
    // a triangle with no area, and an index past the end of the vertices
    mesh.indices.extend([0, 0, 0, 2, 3, 99]);
    let image = paint(mesh, everywhere(), [4, 4], 1.0, &white_texture());
    assert!(image
        .pixels
        .iter()
        .all(|&p| p == Color32::BLACK || p == Color32::RED));
    let mut missing = filled(rect(0.0, 0.0, 2.0, 2.0), Color32::RED);
    missing.texture_id = TextureId::Managed(7);
    let image = paint(missing, everywhere(), [2, 2], 1.0, &white_texture());
    assert!(image.pixels.iter().all(|&p| p == Color32::BLACK));
}
