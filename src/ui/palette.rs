//! The command palette: press Cmd+K (or click the search box of the top bar), type a few letters of what you want to do, and
//! press Enter. It finds every command of the program (`commands.rs`), so nobody has to know which menu a thing is in.

use super::app::{Action, App};
use super::commands::{self, Command};
use super::design::{self, Look, BUTTON_RADIUS, CARD_RADIUS};
use eframe::egui::{
    self, vec2, Align2, Context, CornerRadius, Frame, Id, Key, Margin, Modifiers, Order, RichText,
    Sense, Stroke, TextEdit, Ui,
};

/// How many answers the palette lists, and how wide and how far from the top it is.
const LIMIT: usize = 9;
const WIDTH: f32 = 580.0;
const TOP_OFFSET: f32 = 84.0;
const ROW_HEIGHT: f32 = 34.0;

#[derive(Default)]
pub struct CommandPalette {
    open: bool,
    query: String,
    /// Which of the listed answers is chosen.
    chosen: usize,
    /// The commands, as listed when the palette opened.
    commands: &'static [Command],
}

impl CommandPalette {
    pub fn is_open(&self) -> bool {
        self.open
    }

    /// Open it with nothing typed, or close it when it is open.
    pub fn toggle(&mut self) {
        if self.open {
            self.open = false;
            return;
        }
        *self = CommandPalette {
            open: true,
            commands: commands::all(),
            ..CommandPalette::default()
        };
    }

    /// The commands that answer what is typed.
    fn answers(&self) -> Vec<usize> {
        commands::search(self.commands, &self.query, LIMIT)
    }
}

/// What the keys of the keyboard asked of the list this frame.
struct Keys {
    /// How many times each arrow was pressed in the frame (a held key presses it again and again).
    up: usize,
    down: usize,
    enter: bool,
    escape: bool,
}

fn read_keys(ctx: &Context) -> Keys {
    ctx.input_mut(|i| Keys {
        up: i.count_and_consume_key(Modifiers::NONE, Key::ArrowUp),
        down: i.count_and_consume_key(Modifiers::NONE, Key::ArrowDown),
        enter: i.consume_key(Modifiers::NONE, Key::Enter),
        escape: i.consume_key(Modifiers::NONE, Key::Escape),
    })
}

impl App {
    /// Draw the palette over the window while it is open, and put the command that was chosen in `out`.
    pub(super) fn palette_ui(&mut self, ctx: &Context, out: &mut Vec<Action>) {
        if !self.palette.open {
            return;
        }
        let look = design::look(self.prefs.dark);
        let keys = read_keys(ctx);
        let mut close = keys.escape || self.scrim(ctx, &look);
        let mut run: Option<usize> = None;
        egui::Area::new(Id::new("command_palette"))
            .order(Order::Tooltip)
            .anchor(Align2::CENTER_TOP, vec2(0.0, TOP_OFFSET))
            .show(ctx, |ui| {
                Frame::new()
                    .fill(look.surface)
                    .stroke(Stroke::new(1.0, look.border))
                    .corner_radius(CornerRadius::same(CARD_RADIUS + 2))
                    .inner_margin(Margin::same(10))
                    .shadow(ctx.global_style().visuals.window_shadow)
                    .show(ui, |ui| {
                        // a click on text must not count as a click on nothing: that would take the keyboard from the box
                        ui.style_mut().interaction.selectable_labels = false;
                        ui.set_width(WIDTH);
                        run = self.palette_body(ui, &look, &keys);
                    });
            });
        if let Some(i) = run {
            out.push(self.palette.commands[i].action.clone());
            close = true;
        }
        if close {
            self.palette.open = false;
        }
    }

    /// The dimmed window behind the palette: a click on it closes the palette. True when it was clicked.
    fn scrim(&self, ctx: &Context, look: &Look) -> bool {
        let screen = ctx.content_rect();
        egui::Area::new(Id::new("command_palette_scrim"))
            .order(Order::Foreground)
            .fixed_pos(screen.min)
            .show(ctx, |ui| {
                ui.painter().rect_filled(screen, 0.0, look.scrim);
                ui.allocate_rect(screen, Sense::click()).clicked()
            })
            .inner
    }

    /// The box and the list inside the palette's frame; gives the place in the list of the command to run, if one was chosen.
    fn palette_body(&mut self, ui: &mut Ui, look: &Look, keys: &Keys) -> Option<usize> {
        let search = TextEdit::singleline(&mut self.palette.query)
            .id(Id::new("palette_query"))
            .lock_focus(true)
            .hint_text("Type a command, a tool or a place…")
            .font(egui::FontId::proportional(17.0))
            .frame(Frame::new())
            .desired_width(f32::INFINITY)
            .margin(Margin::symmetric(8, 8));
        let field = ui.add(search);
        // the box keeps the keyboard for as long as the palette is open, whatever was clicked, so that nothing typed reaches
        // the page underneath
        field.request_focus();
        if field.changed() {
            self.palette.chosen = 0;
        }
        ui.add_space(2.0);
        ui.separator();
        let answers = self.palette.answers();
        if answers.is_empty() {
            ui.add_space(10.0);
            ui.label(
                RichText::new("Nothing matches. Try a word such as align, tree, primer or save.")
                    .color(look.dim),
            );
            ui.add_space(10.0);
            return None;
        }
        let last = answers.len() - 1;
        self.palette.chosen = (self.palette.chosen + keys.down).min(last);
        self.palette.chosen = self.palette.chosen.saturating_sub(keys.up);
        self.palette.chosen = self.palette.chosen.min(last);
        let mut run = keys.enter.then(|| answers[self.palette.chosen]);
        ui.add_space(4.0);
        for (place, &index) in answers.iter().enumerate() {
            let command = &self.palette.commands[index];
            let (clicked, hovered) = row(ui, look, command, place, place == self.palette.chosen);
            if hovered {
                self.palette.chosen = place;
            }
            if clicked {
                run = Some(index);
            }
        }
        ui.add_space(6.0);
        ui.label(
            RichText::new("Up and Down to choose, Enter to run, Esc to close")
                .size(11.5)
                .color(look.dim),
        );
        run
    }
}

/// The id of the answer `place` places down the list.
fn row_id(place: usize) -> Id {
    Id::new(("palette_row", place))
}

/// One answer: its group, its name and its key. Gives (clicked, hovered).
fn row(ui: &mut Ui, look: &Look, command: &Command, place: usize, chosen: bool) -> (bool, bool) {
    let rect = ui.allocate_space(vec2(ui.available_width(), ROW_HEIGHT)).1;
    let response = ui.interact(rect, row_id(place), Sense::click());
    if chosen {
        ui.painter()
            .rect_filled(rect, CornerRadius::same(BUTTON_RADIUS), look.accent_soft);
    }
    let ink = if chosen { look.accent } else { look.text };
    let mut inner = ui.new_child(
        egui::UiBuilder::new()
            .max_rect(rect.shrink2(vec2(10.0, 0.0)))
            .layout(egui::Layout::left_to_right(egui::Align::Center)),
    );
    inner.label(RichText::new(command.title).size(14.5).color(ink));
    inner.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
        if !command.keys.is_empty() {
            let keys = command.keys.replace("Cmd", design::command_key());
            design::keycap(ui, look, &keys);
        }
        ui.label(RichText::new(command.group).size(12.0).color(look.dim));
    });
    (
        response.clicked(),
        response.hovered() && ui.input(|i| i.pointer.delta() != vec2(0.0, 0.0)),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ui::app::Startup;
    use crate::ui::headless::frame_input;
    use crate::ui::platform::Workspace;
    use eframe::egui::Event;

    fn key(key: Key) -> Event {
        Event::Key {
            key,
            physical_key: Some(key),
            pressed: true,
            repeat: false,
            modifiers: Modifiers::NONE,
        }
    }

    /// One frame of the whole window with `events` arriving in it.
    fn frame(app: &mut App, ctx: &Context, n: u32, events: Vec<Event>) {
        let mut input = frame_input(n);
        input.events = events;
        let mut output = ctx.run_ui(input, |ui| app.draw(ui));
        output.textures_delta.clear();
    }

    /// A click on the name of the answer `place` places down the list, in the middle of the text.
    fn click_on_the_name(app: &mut App, ctx: &Context, n: u32, place: usize) {
        use eframe::egui::PointerButton;
        let rect = ctx
            .read_response(row_id(place))
            .expect("the row was drawn")
            .rect;
        let at = eframe::egui::pos2(rect.left() + 30.0, rect.center().y);
        let button = |pressed| Event::PointerButton {
            pos: at,
            button: PointerButton::Primary,
            pressed,
            modifiers: Modifiers::NONE,
        };
        frame(app, ctx, n, vec![Event::PointerMoved(at)]);
        frame(app, ctx, n + 1, vec![button(true), button(false)]);
    }

    #[test]
    fn a_click_on_the_name_of_an_answer_runs_it_and_a_click_elsewhere_leaves_the_keyboard_with_the_box(
    ) {
        let (ctx, mut app) = (Context::default(), App::new(Startup::default()));
        app.palette.toggle();
        frame(&mut app, &ctx, 0, vec![]);
        frame(&mut app, &ctx, 1, vec![Event::Text("tr".into())]);
        click_on_the_name(&mut app, &ctx, 2, 0);
        assert!(!app.palette.is_open(), "the click ran the answer");
        assert_eq!(
            app.workspace,
            Workspace::Trees,
            "Trees is the first answer for tr"
        );
        // opened again, a click on the hint line (no answer, no widget) must not hand the keyboard to the editor below
        app.workspace = Workspace::Sequences;
        app.palette.toggle();
        frame(&mut app, &ctx, 10, vec![]);
        frame(&mut app, &ctx, 11, vec![]);
        let at = eframe::egui::pos2(640.0, 94.0 + 6.0);
        frame(&mut app, &ctx, 12, vec![Event::PointerMoved(at)]);
        let press = |pressed| Event::PointerButton {
            pos: at,
            button: eframe::egui::PointerButton::Primary,
            pressed,
            modifiers: Modifiers::NONE,
        };
        frame(&mut app, &ctx, 13, vec![press(true), press(false)]);
        let before = app.docs[0].doc.seq.clone();
        frame(&mut app, &ctx, 14, vec![Event::Text("a".into())]);
        frame(&mut app, &ctx, 15, vec![Event::Text("c".into())]);
        assert_eq!(
            app.docs[0].doc.seq, before,
            "typed letters did not reach the sequence"
        );
        assert_eq!(app.palette.query, "ac", "they went into the palette's box");
    }

    #[test]
    fn typing_narrows_the_list_and_enter_runs_the_answer_that_is_chosen() {
        let (ctx, mut app) = (Context::default(), App::new(Startup::default()));
        app.palette.toggle();
        frame(&mut app, &ctx, 0, vec![]);
        frame(&mut app, &ctx, 1, vec![Event::Text("sanger".into())]);
        assert_eq!(app.palette.query, "sanger");
        frame(&mut app, &ctx, 2, vec![key(Key::Enter)]);
        assert!(
            !app.palette.is_open(),
            "running a command closes the palette"
        );
        assert_eq!(app.workspace, Workspace::Signals);
    }

    #[test]
    fn the_arrow_keys_move_the_choice_and_escape_closes_without_running_anything() {
        let (ctx, mut app) = (Context::default(), App::new(Startup::default()));
        let before = app.workspace;
        app.palette.toggle();
        frame(&mut app, &ctx, 0, vec![]);
        frame(&mut app, &ctx, 1, vec![Event::Text("tree".into())]);
        assert_eq!(app.palette.chosen, 0);
        frame(
            &mut app,
            &ctx,
            2,
            vec![key(Key::ArrowDown), key(Key::ArrowDown)],
        );
        assert_eq!(app.palette.chosen, 2);
        frame(&mut app, &ctx, 3, vec![key(Key::ArrowUp)]);
        assert_eq!(app.palette.chosen, 1);
        frame(&mut app, &ctx, 4, vec![key(Key::Escape)]);
        assert!(!app.palette.is_open());
        assert_eq!(app.workspace, before);
    }

    #[test]
    fn it_opens_empty_with_all_the_commands_and_closes_when_toggled_again() {
        let mut palette = CommandPalette::default();
        assert!(!palette.is_open());
        palette.query = "left over".into();
        palette.toggle();
        assert!(palette.is_open() && palette.query.is_empty());
        assert_eq!(palette.answers().len(), LIMIT);
        palette.toggle();
        assert!(!palette.is_open());
    }

    #[test]
    fn what_is_typed_narrows_the_answers() {
        let mut palette = CommandPalette::default();
        palette.toggle();
        palette.query = "sanger".into();
        let answers = palette.answers();
        assert!(!answers.is_empty() && answers.len() < LIMIT);
        let titles: Vec<&str> = answers.iter().map(|&i| palette.commands[i].title).collect();
        // a word of the name beats a word that only describes the command
        assert_eq!(titles[0], "Try the Sanger trace example");
        assert!(titles.contains(&"Signals"), "{titles:?}");
    }
}
