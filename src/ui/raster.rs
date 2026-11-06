//! A small software rasterizer for egui's triangles: what the graphics card does when it draws the
//! window, done on the processor so a picture can be made with no window at all (see `headless`).
//!
//! egui blends in gamma space with premultiplied alpha (`out = src + dst * (1 - src.a)`) and
//! multiplies the vertex colour by the texture colour for every pixel; this does the same.

use eframe::egui::epaint::{Primitive, Vertex};
use eframe::egui::{
    ClippedPrimitive, Color32, ColorImage, ImageData, Pos2, Rect, TextureFilter, TextureId,
    TexturesDelta,
};
use std::collections::HashMap;

/// One texture the triangles can sample: egui's glyph atlas, and any image the program loaded.
struct Texture {
    size: [usize; 2],
    pixels: Vec<Color32>,
    smooth: bool,
}

impl Texture {
    /// The texel at `(x, y)`, premultiplied, 0..=255; outside the image the nearest edge texel counts.
    fn texel(&self, x: i64, y: i64) -> [f32; 4] {
        let [w, h] = self.size;
        let x = x.clamp(0, w as i64 - 1) as usize;
        let y = y.clamp(0, h as i64 - 1) as usize;
        self.pixels
            .get(y * w + x)
            .map_or([0.0; 4], |c| c.to_array().map(f32::from))
    }

    /// The colour under texture coordinate `uv`, filtered the way the texture asks for.
    // ponytail: filters in gamma space and always clamps at the edges, where the graphics card
    // filters in linear space and can repeat; only visible for an image drawn scaled or tiled,
    // and the glyph atlas and the plain fills are read texel for texel.
    fn sample(&self, uv: Pos2) -> [f32; 4] {
        // texel centres sit half a texel in from the edges
        let x = uv.x * self.size[0] as f32 - 0.5;
        let y = uv.y * self.size[1] as f32 - 0.5;
        if !self.smooth {
            return self.texel((x + 0.5).floor() as i64, (y + 0.5).floor() as i64);
        }
        let (left, top) = (x.floor(), y.floor());
        let (fx, fy) = (x - left, y - top);
        let (left, top) = (left as i64, top as i64);
        let corners = [
            self.texel(left, top),
            self.texel(left + 1, top),
            self.texel(left, top + 1),
            self.texel(left + 1, top + 1),
        ];
        std::array::from_fn(|k| {
            let upper = corners[0][k] * (1.0 - fx) + corners[1][k] * fx;
            let lower = corners[2][k] * (1.0 - fx) + corners[3][k] * fx;
            upper * (1.0 - fy) + lower * fy
        })
    }

    /// Overwrite the part of the texture at `at` with `image`; what lies outside is left out.
    fn patch(&mut self, at: [usize; 2], image: &ColorImage) {
        let [w, h] = self.size;
        if at[0] >= w {
            return;
        }
        for (row, source) in image.pixels.chunks_exact(image.size[0].max(1)).enumerate() {
            let y = at[1] + row;
            if y >= h {
                break;
            }
            let width = source.len().min(w - at[0]);
            let start = y * w + at[0];
            self.pixels[start..start + width].copy_from_slice(&source[..width]);
        }
    }
}

/// Every texture egui has made so far, kept up to date from what it reports each frame.
#[derive(Default)]
pub struct Textures(HashMap<TextureId, Texture>);

impl Textures {
    /// Take in the textures egui made or changed this frame. The ids it freed are handed back: they
    /// stay usable until the frame that may still draw with them has been painted.
    pub fn apply(&mut self, delta: &mut TexturesDelta) -> Vec<TextureId> {
        for (id, patches) in delta.set.drain() {
            for patch in patches {
                let ImageData::Color(image) = patch.image;
                let smooth = patch.options.magnification == TextureFilter::Linear;
                match patch.pos {
                    None if image.size.contains(&0) => {}
                    None => {
                        let texture = Texture {
                            size: image.size,
                            pixels: image.pixels.clone(),
                            smooth,
                        };
                        self.0.insert(id, texture);
                    }
                    Some(at) => {
                        if let Some(texture) = self.0.get_mut(&id) {
                            texture.patch(at, &image);
                        }
                    }
                }
            }
        }
        delta.free.drain().collect()
    }

    pub fn forget(&mut self, ids: &[TextureId]) {
        for id in ids {
            self.0.remove(id);
        }
    }
}

/// Twice the area of triangle `a`, `b`, `p`, positive when `p` is on the left of the edge `a → b`.
fn edge(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])
}

/// Is a pixel whose edge value is `w` inside the edge `from → to`? A pixel exactly on an edge
/// belongs to one of the two triangles that share it and never to both, so a see-through fill
/// shows no seam along the diagonal of a rectangle.
fn inside(w: f64, from: [f64; 2], to: [f64; 2]) -> bool {
    let (dx, dy) = (to[0] - from[0], to[1] - from[1]);
    w > 0.0 || (w == 0.0 && (dy > 0.0 || (dy == 0.0 && dx < 0.0)))
}

/// The colour at barycentric weights `t` of the three vertices: their colours blended, times the
/// texture colour under their blended texture coordinate. Premultiplied, 0..=255.
fn shade(v: [&Vertex; 3], t: [f32; 3], texture: &Texture) -> [f32; 4] {
    let uv = Pos2::new(
        t[0] * v[0].uv.x + t[1] * v[1].uv.x + t[2] * v[2].uv.x,
        t[0] * v[0].uv.y + t[1] * v[1].uv.y + t[2] * v[2].uv.y,
    );
    let texel = texture.sample(uv);
    let colors = v.map(|vertex| vertex.color.to_array().map(f32::from));
    std::array::from_fn(|k| {
        (t[0] * colors[0][k] + t[1] * colors[1][k] + t[2] * colors[2][k]) * texel[k] / 255.0
    })
}

/// The picture being drawn.
pub struct Canvas {
    size: [usize; 2],
    pixels: Vec<Color32>,
}

impl Canvas {
    pub fn new(size: [usize; 2], fill: Color32) -> Canvas {
        Canvas {
            size,
            pixels: vec![fill; size[0] * size[1]],
        }
    }

    pub fn into_image(self) -> ColorImage {
        ColorImage::new(self.size, self.pixels)
    }

    /// Draw what egui's `tessellate` made, front to back as listed, at `pixels_per_point`.
    pub fn paint(&mut self, jobs: &[ClippedPrimitive], textures: &Textures, pixels_per_point: f32) {
        for job in jobs {
            // a paint callback is custom drawing for the graphics card; this program has none
            let Primitive::Mesh(mesh) = &job.primitive else {
                continue;
            };
            let Some(texture) = textures.0.get(&mesh.texture_id) else {
                continue;
            };
            let clip = self.clip(job.clip_rect, pixels_per_point);
            for &[i, j, k] in mesh.indices.as_chunks::<3>().0 {
                let vertex = |index: u32| mesh.vertices.get(index as usize);
                if let (Some(a), Some(b), Some(c)) = (vertex(i), vertex(j), vertex(k)) {
                    self.fill_triangle([a, b, c], texture, clip, pixels_per_point);
                }
            }
        }
    }

    /// The pixels a clip rectangle (in points) lets through, as `[left, top, right, bottom]`.
    fn clip(&self, rect: Rect, pixels_per_point: f32) -> [usize; 4] {
        let to_pixel =
            |v: f32, limit: usize| (v * pixels_per_point).round().clamp(0.0, limit as f32) as usize;
        [
            to_pixel(rect.min.x, self.size[0]),
            to_pixel(rect.min.y, self.size[1]),
            to_pixel(rect.max.x, self.size[0]),
            to_pixel(rect.max.y, self.size[1]),
        ]
    }

    fn fill_triangle(
        &mut self,
        corners: [&Vertex; 3],
        texture: &Texture,
        clip: [usize; 4],
        pixels_per_point: f32,
    ) {
        let at = |v: &Vertex| {
            [
                f64::from(v.pos.x * pixels_per_point),
                f64::from(v.pos.y * pixels_per_point),
            ]
        };
        let signed_area = edge(at(corners[0]), at(corners[1]), at(corners[2]));
        if signed_area.is_nan() || signed_area == 0.0 {
            return;
        }
        // the tests below want the corners going one way round: turn a triangle that goes the other
        let [a, b, c] = if signed_area < 0.0 {
            [corners[0], corners[2], corners[1]]
        } else {
            corners
        };
        let (pa, pb, pc) = (at(a), at(b), at(c));
        if ![pa, pb, pc].iter().flatten().all(|v| v.is_finite()) {
            return;
        }
        let area = signed_area.abs();
        // ponytail: tests every pixel of the triangle's bounding box, which wastes work on a long
        // thin diagonal; walk the edges instead if such a triangle ever makes a picture slow
        let reach = |pick: fn(f64, f64) -> f64, k: usize| pick(pick(pa[k], pb[k]), pc[k]);
        let first = |low: f64, limit: usize| low.floor().max(limit as f64) as usize;
        let past = |high: f64, limit: usize| high.ceil().min(limit as f64) as usize;
        let (x_from, x_to) = (
            first(reach(f64::min, 0), clip[0]),
            past(reach(f64::max, 0), clip[2]),
        );
        let (y_from, y_to) = (
            first(reach(f64::min, 1), clip[1]),
            past(reach(f64::max, 1), clip[3]),
        );
        for y in y_from..y_to {
            for x in x_from..x_to {
                let p = [x as f64 + 0.5, y as f64 + 0.5];
                let w = [edge(pb, pc, p), edge(pc, pa, p), edge(pa, pb, p)];
                if inside(w[0], pb, pc) && inside(w[1], pc, pa) && inside(w[2], pa, pb) {
                    let weights = w.map(|v| (v / area) as f32);
                    self.blend(x, y, shade([a, b, c], weights, texture));
                }
            }
        }
    }

    /// Put `src` (premultiplied, 0..=255) over the pixel at `(x, y)`.
    fn blend(&mut self, x: usize, y: usize, src: [f32; 4]) {
        let Some(dst) = self.pixels.get_mut(y * self.size[0] + x) else {
            return;
        };
        let keep = 1.0 - src[3] / 255.0;
        let over = |s: f32, d: u8| (s + f32::from(d) * keep).round().clamp(0.0, 255.0) as u8;
        *dst = Color32::from_rgba_premultiplied(
            over(src[0], dst.r()),
            over(src[1], dst.g()),
            over(src[2], dst.b()),
            over(src[3], dst.a()),
        );
    }
}

#[cfg(test)]
mod tests;
