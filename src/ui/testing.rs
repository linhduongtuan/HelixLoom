//! For the tests that draw the whole window: a screen that runs frames, reads the text that they drew and clicks a button by its
//! words, as a mouse would.

use super::app::App;
use super::headless::frame_input_at;
use eframe::egui::{
    self,
    epaint::{ClippedShape, Shape},
    Event, Modifiers, PointerButton, Pos2,
};

/// A text of the last frame: where it begins, what it says, how large its letters are and the room it takes.
pub struct Words {
    pub at: Pos2,
    pub text: String,
    pub size: f32,
    pub room: egui::Rect,
}

/// How near the edge of its column a button may be and still be clicked where it is: nearer, and it is cut or about to be.
const OUT_OF_REACH: f32 = 20.0;
/// The frames a column is given to glide to where the wheel sent it.
const GLIDE: usize = 40;

pub struct Screen {
    ctx: egui::Context,
    frame: u32,
    shapes: Vec<ClippedShape>,
    /// The text the program last put on the clipboard.
    copied: Option<String>,
}

impl Screen {
    pub fn new() -> Screen {
        Screen {
            ctx: egui::Context::default(),
            frame: 0,
            shapes: Vec::new(),
            copied: None,
        }
    }

    /// Draw one frame of `app`, with `events` arriving in it.
    pub fn frame(&mut self, app: &mut App, events: Vec<Event>) {
        // a window that is zoomed holds fewer points, as it does when it is drawn for a picture
        let mut input = frame_input_at(self.frame, self.ctx.zoom_factor());
        input.events = events;
        self.frame += 1;
        let mut out = self.ctx.run_ui(input, |ui| app.draw(ui));
        out.textures_delta.clear();
        for command in &out.platform_output.commands {
            if let egui::OutputCommand::CopyText(text) = command {
                self.copied = Some(text.clone());
            }
        }
        self.shapes = out.shapes;
    }

    /// The text the program last put on the clipboard, in any frame so far.
    pub fn copied(&self) -> Option<&str> {
        self.copied.as_deref()
    }

    /// The context that the frames are drawn with: what the window did to it (the zoom, the spacing, the colours) is read here.
    pub fn ctx(&self) -> &egui::Context {
        &self.ctx
    }

    /// Draw a few frames, for a page that takes one to lay itself out.
    pub fn settle(&mut self, app: &mut App) {
        for _ in 0..3 {
            self.frame(app, vec![]);
        }
    }

    /// All the text of the last frame, one piece to a line.
    pub fn text(&self) -> String {
        fn walk(shape: &Shape, out: &mut String) {
            match shape {
                Shape::Text(t) => {
                    out.push_str(t.galley.text());
                    out.push('\n');
                }
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, out)),
                _ => {}
            }
        }
        let mut out = String::new();
        self.shapes.iter().for_each(|c| walk(&c.shape, &mut out));
        out
    }

    /// The middle of the text `label` in the last frame, which is on the button that holds it.
    pub fn find(&self, label: &str) -> Option<Pos2> {
        fn walk(shape: &Shape, label: &str) -> Option<Pos2> {
            match shape {
                Shape::Text(t) if t.galley.text() == label => {
                    Some(t.pos + t.galley.rect.center().to_vec2())
                }
                Shape::Vec(v) => v.iter().find_map(|s| walk(s, label)),
                _ => None,
            }
        }
        self.shapes.iter().find_map(|c| walk(&c.shape, label))
    }

    /// Every text of the last frame that begins inside `room`, from the top down (and from the left, on one line).
    pub fn words_in(&self, room: egui::Rect) -> Vec<Words> {
        fn walk(shape: &Shape, room: egui::Rect, out: &mut Vec<Words>) {
            match shape {
                Shape::Text(t) if room.contains(t.pos) => out.push(Words {
                    at: t.pos,
                    room: t.galley.rect.translate(t.pos.to_vec2()),
                    text: t.galley.text().to_string(),
                    size: t
                        .galley
                        .job
                        .sections
                        .first()
                        .map_or(0.0, |s| s.format.font_id.size),
                }),
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, room, out)),
                _ => {}
            }
        }
        let mut out = Vec::new();
        self.shapes
            .iter()
            .for_each(|c| walk(&c.shape, room, &mut out));
        out.sort_by(|a, b| (a.at.y, a.at.x).partial_cmp(&(b.at.y, b.at.x)).unwrap());
        out
    }

    /// Whether the first text `label` of the last frame is drawn pale, as the words of a button that cannot be pressed
    /// are: seen through, not in the full colour. None when there is no such text.
    pub fn is_pale(&self, label: &str) -> Option<bool> {
        fn walk(shape: &Shape, label: &str) -> Option<bool> {
            match shape {
                Shape::Text(t) if t.galley.text() == label => {
                    // letters that were given no colour of their own are drawn in the colour the text falls back on
                    let own = t
                        .galley
                        .rows
                        .iter()
                        .flat_map(|r| r.row.visuals.mesh.vertices.iter())
                        .map(|v| v.color)
                        .find(|c| *c != egui::Color32::PLACEHOLDER);
                    let colour = t.override_text_color.or(own).unwrap_or(t.fallback_color);
                    Some(colour.a() < 255)
                }
                Shape::Vec(v) => v.iter().find_map(|s| walk(s, label)),
                _ => None,
            }
        }
        self.shapes.iter().find_map(|c| walk(&c.shape, label))
    }

    /// How far left anything is seen that begins inside `room`: a text, or a rectangle that is filled or outlined (a
    /// shadow is not: it reaches beyond what casts it). What a view cuts away at its edge is not seen: a label of a tree
    /// that begins left of the tree's frame counts from the frame. None when nothing does.
    pub fn leftmost_in(&self, room: egui::Rect) -> Option<f32> {
        fn walk(shape: &Shape, room: egui::Rect, out: &mut Vec<f32>) {
            match shape {
                Shape::Text(t) if room.contains(t.pos) => out.push(t.pos.x),
                Shape::Rect(r)
                    if r.blur_width == 0.0
                        && (r.fill.a() > 0 || r.stroke.width > 0.0)
                        && room.contains(r.rect.min) =>
                {
                    out.push(r.rect.min.x);
                }
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, room, out)),
                _ => {}
            }
        }
        let mut seen = Vec::new();
        for clipped in &self.shapes {
            let mut out = Vec::new();
            walk(&clipped.shape, room, &mut out);
            seen.extend(out.into_iter().map(|x| x.max(clipped.clip_rect.min.x)));
        }
        seen.into_iter().reduce(f32::min)
    }

    /// The middle of every text `label` in the last frame, in the order they were drawn.
    pub fn find_all(&self, label: &str) -> Vec<Pos2> {
        fn walk(shape: &Shape, label: &str, out: &mut Vec<Pos2>) {
            match shape {
                Shape::Text(t) if t.galley.text() == label => {
                    out.push(t.pos + t.galley.rect.center().to_vec2());
                }
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, label, out)),
                _ => {}
            }
        }
        let mut out = Vec::new();
        self.shapes
            .iter()
            .for_each(|c| walk(&c.shape, label, &mut out));
        out
    }

    /// The middle of the first circle in the last frame that is filled with each of `colours` (those that were not drawn are left out).
    pub fn circles_of(&self, colours: &[egui::Color32]) -> Vec<Pos2> {
        fn walk(shape: &Shape, colour: egui::Color32) -> Option<Pos2> {
            match shape {
                Shape::Circle(c) if c.fill == colour => Some(c.center),
                Shape::Vec(v) => v.iter().find_map(|s| walk(s, colour)),
                _ => None,
            }
        }
        colours
            .iter()
            .filter_map(|&colour| self.shapes.iter().find_map(|c| walk(&c.shape, colour)))
            .collect()
    }

    /// The middle of every circle in the last frame that has `radius` (a ring that something is marked with).
    pub fn circles_of_radius(&self, radius: f32) -> Vec<Pos2> {
        fn walk(shape: &Shape, radius: f32, out: &mut Vec<Pos2>) {
            match shape {
                Shape::Circle(c) if (c.radius - radius).abs() < 0.01 => out.push(c.center),
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, radius, out)),
                _ => {}
            }
        }
        let mut out = Vec::new();
        self.shapes
            .iter()
            .for_each(|c| walk(&c.shape, radius, &mut out));
        out
    }

    /// Every rectangle in the last frame that is outlined in `colour` (a mark round something), in the order they were drawn.
    pub fn rects_outlined(&self, colour: egui::Color32) -> Vec<egui::Rect> {
        fn walk(shape: &Shape, colour: egui::Color32, out: &mut Vec<egui::Rect>) {
            match shape {
                Shape::Rect(r) if r.stroke.width > 0.0 && r.stroke.color == colour => {
                    out.push(r.rect);
                }
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, colour, out)),
                _ => {}
            }
        }
        let mut out = Vec::new();
        self.shapes
            .iter()
            .for_each(|c| walk(&c.shape, colour, &mut out));
        out
    }

    /// Every rectangle in the last frame that is filled with `colour`, in the order they were drawn.
    pub fn rects_filled(&self, colour: egui::Color32) -> Vec<egui::Rect> {
        fn walk(shape: &Shape, colour: egui::Color32, out: &mut Vec<egui::Rect>) {
            match shape {
                Shape::Rect(r) if r.fill == colour => out.push(r.rect),
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, colour, out)),
                _ => {}
            }
        }
        let mut out = Vec::new();
        self.shapes
            .iter()
            .for_each(|c| walk(&c.shape, colour, &mut out));
        out
    }

    /// Every filled rectangle in the last frame that `at` lies in, with its colour, in the order they were drawn: what
    /// stands behind a text.
    pub fn fills_at(&self, at: Pos2) -> Vec<(egui::Rect, egui::Color32)> {
        fn walk(shape: &Shape, at: Pos2, out: &mut Vec<(egui::Rect, egui::Color32)>) {
            match shape {
                Shape::Rect(r) if r.fill.a() > 0 && r.rect.contains(at) => {
                    out.push((r.rect, r.fill));
                }
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, at, out)),
                _ => {}
            }
        }
        let mut out = Vec::new();
        self.shapes
            .iter()
            .for_each(|c| walk(&c.shape, at, &mut out));
        out
    }

    /// The room of every outline in the last frame that is drawn as a path (a shape with corners, not a plain
    /// rectangle) with a line of `colour` that is `width` heavy, in the order they were drawn.
    pub fn paths_stroked(&self, colour: egui::Color32, width: f32) -> Vec<egui::Rect> {
        fn walk(shape: &Shape, colour: egui::Color32, width: f32, out: &mut Vec<egui::Rect>) {
            match shape {
                Shape::Path(p)
                    if (p.stroke.width - width).abs() < 0.01
                        && matches!(p.stroke.color, egui::epaint::ColorMode::Solid(c) if c == colour) =>
                {
                    out.push(egui::Rect::from_points(&p.points));
                }
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, colour, width, out)),
                _ => {}
            }
        }
        let mut out = Vec::new();
        self.shapes
            .iter()
            .for_each(|c| walk(&c.shape, colour, width, &mut out));
        out
    }

    /// The room of every path in the last frame that is filled with `colour`, in the order they were drawn.
    pub fn paths_filled(&self, colour: egui::Color32) -> Vec<egui::Rect> {
        fn walk(shape: &Shape, colour: egui::Color32, out: &mut Vec<egui::Rect>) {
            match shape {
                Shape::Path(p) if p.fill == colour => out.push(egui::Rect::from_points(&p.points)),
                Shape::Vec(v) => v.iter().for_each(|s| walk(s, colour, out)),
                _ => {}
            }
        }
        let mut out = Vec::new();
        self.shapes
            .iter()
            .for_each(|c| walk(&c.shape, colour, &mut out));
        out
    }

    /// The middle of the first text `label` of the last frame, and the room it may be drawn in: what lies outside that
    /// room is cut away, as what has been scrolled out of a column is.
    fn find_in_its_room(&self, label: &str) -> Option<(Pos2, egui::Rect)> {
        fn walk(shape: &Shape, label: &str) -> Option<Pos2> {
            match shape {
                Shape::Text(t) if t.galley.text() == label => {
                    Some(t.pos + t.galley.rect.center().to_vec2())
                }
                Shape::Vec(v) => v.iter().find_map(|s| walk(s, label)),
                _ => None,
            }
        }
        self.shapes
            .iter()
            .find_map(|c| walk(&c.shape, label).map(|at| (at, c.clip_rect)))
    }

    /// Click the button that says `label`: move onto it, press, let go. A button that is below or above what its column
    /// shows is first brought into sight by turning the wheel over the column, as a person would.
    pub fn click(&mut self, app: &mut App, label: &str) {
        let missing =
            |screen: &Screen| -> ! { panic!("no \"{label}\" on the screen:\n{}", screen.text()) };
        let Some((mut at, room)) = self.find_in_its_room(label) else {
            missing(self)
        };
        if at.y < room.min.y + OUT_OF_REACH || at.y > room.max.y - OUT_OF_REACH {
            let to = if at.y > room.center().y {
                room.max.y - 2.0 * OUT_OF_REACH
            } else {
                room.min.y + 2.0 * OUT_OF_REACH
            };
            let wheel = Event::MouseWheel {
                unit: egui::MouseWheelUnit::Point,
                delta: egui::vec2(0.0, to - at.y),
                phase: egui::TouchPhase::Move,
                modifiers: Modifiers::NONE,
            };
            self.frame(app, vec![Event::PointerMoved(room.center())]);
            self.frame(app, vec![wheel]);
            // the column glides to where the wheel sent it
            for _ in 0..GLIDE {
                self.frame(app, vec![]);
            }
            at = self.find(label).unwrap_or_else(|| missing(self));
        }
        self.click_at(app, at);
    }

    /// Click at a place: move onto it, press, let go.
    pub fn click_at(&mut self, app: &mut App, at: Pos2) {
        let button = |pressed| Event::PointerButton {
            pos: at,
            button: PointerButton::Primary,
            pressed,
            modifiers: Modifiers::NONE,
        };
        self.frame(app, vec![Event::PointerMoved(at)]);
        self.frame(app, vec![button(true)]);
        self.frame(app, vec![button(false)]);
        self.frame(app, vec![]);
    }
}
