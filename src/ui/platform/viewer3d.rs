//! The 3D viewer: a protein drawn as a ribbon, as its backbone or as all its bonds, turned with the mouse, with
//! the colour chosen by confidence, secondary structure, chain or position along the chain. Everything is painted
//! with lines and discs, far parts first and faded a little, so no graphics card code is needed.

use crate::seq::fold::Fold;
use crate::seq::pdb::Atom;
use crate::seq::secstruct::Ss;
use crate::ui::design::Look;
use crate::ui::theme::mix;
use eframe::egui::{
    pos2, vec2, Align2, Color32, FontId, Pos2, Rect, Response, Sense, Shape, Stroke, Ui, Vec2,
};

/// How far the camera is from the middle of the structure, in radii.
const CAMERA_DISTANCE: f32 = 3.0;
/// The structure fills this much of the smaller side of the view at zoom 1.
const FILL: f32 = 0.86;
const SUBDIVISIONS: usize = 5;
const MAX_PITCH: f32 = 1.5;
const TURN_PER_PIXEL: f32 = 0.01;
const MIN_ZOOM: f32 = 0.3;
const MAX_ZOOM: f32 = 12.0;
/// Tube widths in angstroms: a loop, a helix and a strand.
const COIL_WIDTH: f32 = 0.55;
const HELIX_WIDTH: f32 = 1.5;
const STRAND_WIDTH: f32 = 1.3;
/// Widths in angstroms of the backbone and of the bonds of the all-bonds drawing.
const BACKBONE_WIDTH: f32 = 0.35;
const STICK_WIDTH: f32 = 0.28;
/// Every line is drawn on a slightly wider darker one, this many points wider.
const OUTLINE: f32 = 1.6;
const MIN_LINE: f32 = 1.4;
const PICK_RADIUS: f32 = 14.0;
/// Stick drawings of structures with more atoms than this fall back to the backbone.
pub const MAX_STICK_ATOMS: usize = 40_000;
/// What the pointer does, in the corner of the picture.
const HINT: &str = "drag to turn · scroll to zoom · shift-drag to move · click a residue to pick it · double-click to start again";

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Representation {
    Ribbon,
    Backbone,
    Sticks,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Scheme {
    /// pLDDT bands as AlphaFold draws them.
    Confidence,
    Secondary,
    Chain,
    /// Blue at the start of the chain to red at its end.
    Rainbow,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Camera {
    pub yaw: f32,
    pub pitch: f32,
    pub zoom: f32,
    pub pan: Vec2,
}

impl Default for Camera {
    fn default() -> Self {
        Camera {
            yaw: 0.5,
            pitch: -0.35,
            zoom: 1.0,
            pan: Vec2::ZERO,
        }
    }
}

/// How the structure's coordinates become points on the screen.
pub struct Projection {
    rotation: [[f32; 3]; 3],
    centre: [f32; 3],
    origin: Pos2,
    /// Points on the screen for one angstrom at the middle depth.
    scale: f32,
    distance: f32,
}

impl Projection {
    /// The view of `fold` in `rect`: at zoom 1 what is drawn fills the view, the ribbon and the backbone fitted by
    /// the alpha carbons and the all-bonds drawing by every atom.
    pub fn new(
        camera: &Camera,
        fold: &Fold,
        rect: Rect,
        representation: Representation,
    ) -> Projection {
        let (sy, cy) = camera.yaw.sin_cos();
        let (sp, cp) = camera.pitch.sin_cos();
        // a turn about the vertical axis, then a tilt about the horizontal one
        let rotation = [
            [cy, 0.0, sy],
            [sp * sy, cp, -sp * cy],
            [-cp * sy, sp, cp * cy],
        ];
        let fit = match representation {
            Representation::Sticks => fold.radius,
            Representation::Ribbon | Representation::Backbone => fold.backbone_radius,
        }
        .max(1.0);
        Projection {
            rotation,
            centre: fold.centre,
            origin: rect.center() + camera.pan,
            scale: rect.width().min(rect.height()) / 2.0 * FILL / fit * camera.zoom,
            distance: CAMERA_DISTANCE * fit,
        }
    }

    /// The point on the screen, and how near the viewer it is (larger is nearer).
    pub fn project(&self, p: [f32; 3]) -> (Pos2, f32) {
        let v = [
            p[0] - self.centre[0],
            p[1] - self.centre[1],
            p[2] - self.centre[2],
        ];
        let r = self.rotation;
        let at = |row: usize| r[row][0] * v[0] + r[row][1] * v[1] + r[row][2] * v[2];
        let (x, y, z) = (at(0), at(1), at(2));
        let k = self.distance / (self.distance - z).max(self.distance * 0.1);
        (
            self.origin + vec2(x * k * self.scale, -y * k * self.scale),
            z,
        )
    }

    /// Points on the screen for one angstrom at the depth of `z`.
    fn size_at(&self, z: f32) -> f32 {
        self.scale * self.distance / (self.distance - z).max(self.distance * 0.1)
    }
}

/// pLDDT as AlphaFold colours it: very high blue, confident light blue, low yellow, very low orange.
pub fn confidence_colour(c: f32) -> Color32 {
    match c {
        c if c > 90.0 => Color32::from_rgb(0x00, 0x53, 0xd6),
        c if c > 70.0 => Color32::from_rgb(0x65, 0xcb, 0xf3),
        c if c > 50.0 => Color32::from_rgb(0xff, 0xdb, 0x13),
        _ => Color32::from_rgb(0xff, 0x7d, 0x45),
    }
}

/// A rainbow from blue (0) through green to red (1).
pub fn rainbow(t: f32) -> Color32 {
    let hue = (1.0 - t.clamp(0.0, 1.0)) * 240.0;
    let (c, x) = (1.0, 1.0 - ((hue / 60.0) % 2.0 - 1.0).abs());
    let (r, g, b) = match hue {
        h if h < 60.0 => (c, x, 0.0),
        h if h < 120.0 => (x, c, 0.0),
        h if h < 180.0 => (0.0, c, x),
        _ => (0.0, x, c),
    };
    let level = |v: f32| (255.0 * (0.15 + 0.7 * v)).round() as u8;
    Color32::from_rgb(level(r), level(g), level(b))
}

const CHAIN_COLOURS: [[u8; 3]; 8] = [
    [0x4a, 0x8f, 0xe7],
    [0xe0, 0x7a, 0x5f],
    [0x5f, 0xb8, 0x7e],
    [0xc2, 0x7b, 0xd6],
    [0xe3, 0xb5, 0x4a],
    [0x4a, 0xb5, 0xc4],
    [0xd6, 0x5f, 0x8f],
    [0x8a, 0x93, 0x4a],
];

pub fn secondary_colour(ss: Ss) -> Color32 {
    match ss {
        Ss::Helix => Color32::from_rgb(0xe0, 0x40, 0x5a),
        Ss::Strand => Color32::from_rgb(0xf0, 0xb0, 0x20),
        Ss::Coil => Color32::from_rgb(0x8f, 0x9d, 0xb2),
    }
}

/// The colour of residue `k` of `fold` (an index into `fold.residues`).
pub fn residue_colour(fold: &Fold, k: usize, scheme: Scheme) -> Color32 {
    let r = &fold.residues[k];
    match scheme {
        Scheme::Confidence => r
            .confidence
            .map_or(Color32::from_gray(0x9a), confidence_colour),
        Scheme::Secondary => secondary_colour(r.ss),
        Scheme::Chain => {
            let c = CHAIN_COLOURS[r.chain_index % CHAIN_COLOURS.len()];
            Color32::from_rgb(c[0], c[1], c[2])
        }
        Scheme::Rainbow => rainbow(k as f32 / (fold.residues.len().max(2) - 1) as f32),
    }
}

/// The colour of a ligand atom by its element.
fn element_colour(a: &Atom) -> Color32 {
    match a.element.as_str() {
        "N" => Color32::from_rgb(0x4a, 0x6f, 0xe0),
        "O" => Color32::from_rgb(0xe0, 0x4a, 0x4a),
        "S" => Color32::from_rgb(0xe0, 0xc0, 0x30),
        "P" => Color32::from_rgb(0xe0, 0x90, 0x30),
        "FE" | "ZN" | "MG" | "CA" | "MN" | "CU" => Color32::from_rgb(0xb0, 0x70, 0xd0),
        _ => Color32::from_gray(0x90),
    }
}

/// A thick line with round ends. What is drawn is first collected as these, then painted from the far side to
/// the near.
struct Drawn {
    depth: f32,
    a: Pos2,
    b: Pos2,
    width: f32,
    colour: Color32,
}

/// A point on the Catmull-Rom curve through `p1` and `p2`, `t` of the way from one to the other.
fn spline(p: [[f32; 3]; 4], t: f32) -> [f32; 3] {
    let (t2, t3) = (t * t, t * t * t);
    std::array::from_fn(|i| {
        0.5 * (2.0 * p[1][i]
            + (-p[0][i] + p[2][i]) * t
            + (2.0 * p[0][i] - 5.0 * p[1][i] + 4.0 * p[2][i] - p[3][i]) * t2
            + (-p[0][i] + 3.0 * p[1][i] - 3.0 * p[2][i] + p[3][i]) * t3)
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Style {
    pub representation: Representation,
    pub scheme: Scheme,
    /// The picked residue, an index into `fold.residues`.
    pub picked: Option<usize>,
}

impl Default for Style {
    fn default() -> Self {
        Style {
            representation: Representation::Ribbon,
            scheme: Scheme::Secondary,
            picked: None,
        }
    }
}

/// What the pointer did to the picture.
#[derive(Default)]
pub struct Touched {
    pub picked: Option<usize>,
    pub cleared: bool,
}

fn collect(fold: &Fold, proj: &Projection, style: &Style, bonds: &[(usize, usize)]) -> Vec<Drawn> {
    let mut out = Vec::new();
    let colour = |k: usize| residue_colour(fold, k, style.scheme);
    let sticks = style.representation == Representation::Sticks
        && fold.structure.atoms.len() <= MAX_STICK_ATOMS;
    match style.representation {
        Representation::Ribbon => {
            for run in &fold.runs {
                let at = |i: usize| fold.residues[run[i.min(run.len() - 1)]].ca;
                for s in 0..run.len() {
                    let last = s + 1 == run.len();
                    let steps = if last { 0 } else { SUBDIVISIONS };
                    let mut previous = proj.project(at(s));
                    for sub in 1..=steps {
                        let t = sub as f32 / SUBDIVISIONS as f32;
                        let p = spline([at(s.saturating_sub(1)), at(s), at(s + 1), at(s + 2)], t);
                        let here = proj.project(p);
                        // the first half of the way belongs to this residue, the second to the next
                        let k = run[if t < 0.5 { s } else { s + 1 }];
                        let width = match fold.residues[k].ss {
                            Ss::Helix => HELIX_WIDTH,
                            Ss::Strand => STRAND_WIDTH,
                            Ss::Coil => COIL_WIDTH,
                        };
                        let depth = (previous.1 + here.1) / 2.0;
                        out.push(Drawn {
                            depth,
                            a: previous.0,
                            b: here.0,
                            width: (width * proj.size_at(depth)).max(MIN_LINE),
                            colour: colour(k),
                        });
                        previous = here;
                    }
                }
            }
        }
        Representation::Backbone => {
            for run in &fold.runs {
                for (i, &k) in run.iter().enumerate() {
                    let b = &fold.backbone[k];
                    let mut link = |a: [f32; 3], c: [f32; 3]| {
                        let (pa, pc) = (proj.project(a), proj.project(c));
                        let depth = (pa.1 + pc.1) / 2.0;
                        out.push(Drawn {
                            depth,
                            a: pa.0,
                            b: pc.0,
                            width: (BACKBONE_WIDTH * proj.size_at(depth)).max(MIN_LINE),
                            colour: colour(k),
                        });
                    };
                    link(b.n, b.ca);
                    link(b.ca, b.c);
                    link(b.c, b.o);
                    if let Some(&next) = run.get(i + 1) {
                        link(b.c, fold.backbone[next].n);
                    }
                }
            }
        }
        Representation::Sticks if sticks => {
            let residue_colour: std::collections::HashMap<usize, usize> = fold
                .residues
                .iter()
                .enumerate()
                .map(|(k, r)| (r.index, k))
                .collect();
            let atoms = &fold.structure.atoms;
            let tone = |i: usize| match residue_colour.get(&atoms[i].residue) {
                Some(&k) => colour(k),
                None => element_colour(&atoms[i]),
            };
            for &(a, b) in bonds {
                let (pa, pb) = (proj.project(atoms[a].pos), proj.project(atoms[b].pos));
                let depth = (pa.1 + pb.1) / 2.0;
                out.push(Drawn {
                    depth,
                    a: pa.0,
                    b: pb.0,
                    width: (STICK_WIDTH * proj.size_at(depth)).max(MIN_LINE),
                    colour: tone(a),
                });
            }
        }
        Representation::Sticks => {}
    }
    out.sort_by(|a, b| a.depth.total_cmp(&b.depth));
    out
}

/// Draw the fold in the free space of `ui`; turn it by dragging, zoom with the wheel, move it with shift or the
/// right button, and double-click to start again.
pub fn show(
    ui: &mut Ui,
    look: &Look,
    fold: &Fold,
    camera: &mut Camera,
    style: &Style,
    bonds: &[(usize, usize)],
    spin: bool,
) -> (Response, Touched) {
    let size = ui.available_size_before_wrap().max(vec2(120.0, 120.0));
    let (response, painter) = ui.allocate_painter(size, Sense::click_and_drag());
    let rect = response.rect;
    if response.dragged() {
        let d = response.drag_delta();
        let move_it = ui.input(|i| i.modifiers.shift)
            || response.dragged_by(eframe::egui::PointerButton::Secondary);
        if move_it {
            camera.pan += d;
        } else {
            camera.yaw += d.x * TURN_PER_PIXEL;
            camera.pitch = (camera.pitch + d.y * TURN_PER_PIXEL).clamp(-MAX_PITCH, MAX_PITCH);
        }
    }
    if response.hovered() {
        let (scroll, pinch) = ui.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
        camera.zoom = (camera.zoom * pinch * (scroll * 0.004).exp()).clamp(MIN_ZOOM, MAX_ZOOM);
    }
    if response.double_clicked() {
        *camera = Camera::default();
    }
    if spin {
        camera.yaw += 0.008;
        ui.ctx().request_repaint();
    }

    let proj = Projection::new(camera, fold, rect, style.representation);
    let drawn = collect(fold, &proj, style, bonds);
    let (near, far) = drawn.iter().fold((f32::MIN, f32::MAX), |(hi, lo), d| {
        (hi.max(d.depth), lo.min(d.depth))
    });
    let span = (near - far).max(f32::EPSILON);
    let background = look.surface;
    for d in &drawn {
        // farther parts fade towards the background
        let nearness = (d.depth - far) / span;
        let fade = |c: Color32| mix(c, background, 58.0 + 42.0 * nearness);
        let outline = mix(d.colour, Color32::BLACK, 55.0);
        painter.line_segment([d.a, d.b], Stroke::new(d.width + OUTLINE, fade(outline)));
        painter.line_segment([d.a, d.b], Stroke::new(d.width, fade(d.colour)));
        painter.circle_filled(d.b, d.width / 2.0, fade(d.colour));
    }

    // the picked residue
    if let Some(k) = style.picked.filter(|&k| k < fold.residues.len()) {
        let r = &fold.residues[k];
        let (at, _) = proj.project(r.ca);
        painter.add(Shape::circle_stroke(at, 9.0, Stroke::new(2.2, look.accent)));
        let kind = match r.ss {
            Ss::Helix => "helix",
            Ss::Strand => "strand",
            Ss::Coil => "loop",
        };
        let sure = r
            .confidence
            .map_or(String::new(), |c| format!(" · pLDDT {c:.0}"));
        let label = format!("{} {} · {kind}{sure}", char::from(r.letter), r.number);
        painter.text(
            at + vec2(12.0, -12.0),
            Align2::LEFT_BOTTOM,
            label,
            FontId::proportional(13.0),
            look.text,
        );
    }
    legend(&painter, rect, look, style.scheme, fold);
    painter.text(
        rect.right_bottom() - vec2(14.0, 12.0),
        Align2::RIGHT_BOTTOM,
        HINT,
        FontId::proportional(11.5),
        look.dim,
    );

    // picking: the alpha carbon nearest to the click
    let mut touched = Touched::default();
    if response.clicked() {
        if let Some(pointer) = response.interact_pointer_pos() {
            let nearest = fold
                .residues
                .iter()
                .enumerate()
                .map(|(k, r)| (k, (proj.project(r.ca).0 - pointer).length()))
                .filter(|&(_, d)| d <= PICK_RADIUS)
                .min_by(|a, b| a.1.total_cmp(&b.1));
            match nearest {
                Some((k, _)) => touched.picked = Some(k),
                None => touched.cleared = true,
            }
        }
    }
    (response, touched)
}

/// What the colours mean, in the corner.
fn legend(painter: &eframe::egui::Painter, rect: Rect, look: &Look, scheme: Scheme, fold: &Fold) {
    let entries: Vec<(Color32, String)> = match scheme {
        Scheme::Confidence if !fold.predicted => return,
        Scheme::Confidence => vec![
            (
                confidence_colour(95.0),
                "very high (pLDDT > 90)".to_string(),
            ),
            (confidence_colour(80.0), "confident (70 to 90)".to_string()),
            (confidence_colour(60.0), "low (50 to 70)".to_string()),
            (confidence_colour(30.0), "very low (< 50)".to_string()),
        ],
        Scheme::Secondary => vec![
            (secondary_colour(Ss::Helix), "helix".to_string()),
            (secondary_colour(Ss::Strand), "strand".to_string()),
            (secondary_colour(Ss::Coil), "loop".to_string()),
        ],
        Scheme::Chain => fold
            .structure
            .chains()
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let tone = CHAIN_COLOURS[i % CHAIN_COLOURS.len()];
                (
                    Color32::from_rgb(tone[0], tone[1], tone[2]),
                    format!("chain {c}"),
                )
            })
            .collect(),
        Scheme::Rainbow => vec![
            (rainbow(0.0), "N terminus".to_string()),
            (rainbow(1.0), "C terminus".to_string()),
        ],
    };
    let mut y = rect.bottom() - 12.0 - 17.0 * entries.len() as f32;
    for (colour, label) in entries {
        painter.rect_filled(
            Rect::from_min_size(pos2(rect.left() + 14.0, y + 2.0), vec2(12.0, 12.0)),
            3.0,
            colour,
        );
        painter.text(
            pos2(rect.left() + 32.0, y + 8.0),
            Align2::LEFT_CENTER,
            label,
            FontId::proportional(12.0),
            look.dim,
        );
        y += 17.0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::pdb::parse;

    fn fold() -> Fold {
        Fold::new(parse(include_str!("../../../assets/structure/1ubq.pdb")).expect("a structure"))
    }

    #[test]
    fn the_middle_of_the_structure_lands_in_the_middle_of_the_view() {
        let f = fold();
        let rect = Rect::from_min_size(pos2(100.0, 50.0), vec2(800.0, 600.0));
        for camera in [
            Camera::default(),
            Camera {
                yaw: 2.0,
                pitch: 0.7,
                zoom: 3.0,
                pan: Vec2::ZERO,
            },
        ] {
            let p = Projection::new(&camera, &f, rect, Representation::Ribbon);
            let (at, depth) = p.project(f.centre);
            assert!((at - rect.center()).length() < 1e-3, "{at:?}");
            assert!(depth.abs() < 1e-3);
        }
        // panning moves everything the same way
        let panned = Camera {
            pan: vec2(30.0, -10.0),
            ..Camera::default()
        };
        let (at, _) = Projection::new(&panned, &f, rect, Representation::Ribbon).project(f.centre);
        assert!((at - (rect.center() + vec2(30.0, -10.0))).length() < 1e-3);
    }

    #[test]
    fn turning_keeps_distances_and_the_whole_structure_fits_the_view() {
        let f = fold();
        let rect = Rect::from_min_size(Pos2::ZERO, vec2(700.0, 500.0));
        let camera = Camera {
            yaw: 1.1,
            pitch: -0.4,
            zoom: 1.0,
            pan: Vec2::ZERO,
        };
        let p = Projection::new(&camera, &f, rect, Representation::Sticks);
        // the rotation matrix is orthonormal
        for i in 0..3 {
            for j in 0..3 {
                let dot: f32 = (0..3).map(|k| p.rotation[i][k] * p.rotation[j][k]).sum();
                assert!((dot - if i == j { 1.0 } else { 0.0 }).abs() < 1e-5);
            }
        }
        // every atom is inside the view at zoom 1 when all the bonds are drawn
        for a in &f.structure.atoms {
            let (at, _) = p.project(a.pos);
            assert!(rect.expand(1.0).contains(at), "{at:?}");
        }
        // the ribbon is fitted to the alpha carbons: they fill the view, and it is bigger than the all-bonds one
        let ribbon = Projection::new(&camera, &f, rect, Representation::Ribbon);
        for r in &f.residues {
            let (at, _) = ribbon.project(r.ca);
            assert!(rect.expand(1.0).contains(at), "{at:?}");
        }
        let (a, b) = (f.structure.atoms[0].pos, f.structure.atoms[1].pos);
        let length = |p: &Projection| (p.project(a).0 - p.project(b).0).length();
        assert!(length(&ribbon) > length(&p));
        // zooming in makes the same bond longer on the screen
        let zoomed = |zoom: f32| {
            Projection::new(&Camera { zoom, ..camera }, &f, rect, Representation::Sticks)
        };
        assert!(length(&zoomed(2.0)) > 1.9 * length(&zoomed(1.0)));
    }

    #[test]
    fn a_nearer_point_has_the_larger_depth_and_looks_larger() {
        let f = fold();
        let rect = Rect::from_min_size(Pos2::ZERO, vec2(700.0, 500.0));
        let straight = Camera {
            yaw: 0.0,
            pitch: 0.0,
            zoom: 1.0,
            pan: Vec2::ZERO,
        };
        let p = Projection::new(&straight, &f, rect, Representation::Ribbon);
        // with no turn the viewer looks along z: a bigger z is nearer
        let near = [f.centre[0], f.centre[1], f.centre[2] + 10.0];
        let far = [f.centre[0], f.centre[1], f.centre[2] - 10.0];
        assert!(p.project(near).1 > p.project(far).1);
        assert!(p.size_at(p.project(near).1) > p.size_at(p.project(far).1));
    }

    #[test]
    fn colours_follow_the_bands_the_ends_and_the_groups() {
        assert_eq!(confidence_colour(95.0), Color32::from_rgb(0x00, 0x53, 0xd6));
        assert_eq!(confidence_colour(90.0), Color32::from_rgb(0x65, 0xcb, 0xf3));
        assert_eq!(confidence_colour(55.0), Color32::from_rgb(0xff, 0xdb, 0x13));
        assert_eq!(confidence_colour(10.0), Color32::from_rgb(0xff, 0x7d, 0x45));
        assert_ne!(rainbow(0.0), rainbow(1.0));
        assert!(rainbow(0.0).b() > rainbow(0.0).r(), "blue at the start");
        assert!(rainbow(1.0).r() > rainbow(1.0).b(), "red at the end");
        assert_eq!(rainbow(-3.0), rainbow(0.0));
        let f = fold();
        let secondary: std::collections::HashSet<_> = (0..f.residues.len())
            .map(|k| residue_colour(&f, k, Scheme::Secondary))
            .collect();
        assert_eq!(secondary.len(), 3, "helix, strand and loop");
        // an experiment has no confidence to show
        assert_eq!(
            residue_colour(&f, 0, Scheme::Confidence),
            Color32::from_gray(0x9a)
        );
        assert_eq!(
            residue_colour(&f, 0, Scheme::Chain),
            residue_colour(&f, 75, Scheme::Chain)
        );
    }

    #[test]
    fn the_curve_passes_through_its_ends_and_a_straight_run_stays_straight() {
        let p = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        assert_eq!(spline(p, 0.0), [1.0, 0.0, 0.0]);
        assert_eq!(spline(p, 1.0), [2.0, 0.0, 0.0]);
        let mid = spline(p, 0.5);
        assert!((mid[0] - 1.5).abs() < 1e-6 && mid[1] == 0.0 && mid[2] == 0.0);
    }
}
