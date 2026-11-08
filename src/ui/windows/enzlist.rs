//! List Enzymes and the Quick Lists (Enzymes menu). List Enzymes asks what the original's dialog asks (the format, the
//! columns, which enzymes, which of the DNA); a Quick List is the same code with its own choices (L1186-1188). Either makes
//! one result window holding the text of `seq::enzlist`, which Copy as text puts on the clipboard.

use super::sitetools::PROGRAM;
use crate::seq::enzlist::{self, Entry, Format, Span, Which, DEFAULT_COLUMNS, MAX_COLUMNS};
use crate::seq::enzymes::Enzyme;
use crate::ui::app::App;
use crate::ui::design::{self, Look, Tone};
use crate::ui::map::commify;
use crate::ui::state::DocView;
use eframe::egui::{self, Id, Key, Modifiers, RichText, ScrollArea, TextWrapMode, Ui};
use std::collections::BTreeSet;

#[cfg(test)]
mod tests;

/// The three Quick Lists of the Enzymes menu, each a format and a list fixed by the original: Site Usage is
/// `list_enzymes_window $w 1 all`, Unique Sites `... 0 unique` and Absent Sites `... 0 absent`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Quick {
    SiteUsage,
    UniqueSites,
    AbsentSites,
}

impl Quick {
    pub const ALL: [Quick; 3] = [Quick::SiteUsage, Quick::UniqueSites, Quick::AbsentSites];

    pub fn label(self) -> &'static str {
        match self {
            Quick::SiteUsage => "Site Usage",
            Quick::UniqueSites => "Unique Sites",
            Quick::AbsentSites => "Absent Sites",
        }
    }

    fn format_and_list(self) -> (Format, Which) {
        match self {
            Quick::SiteUsage => (Format::Counts, Which::All),
            Quick::UniqueSites => (Format::Names, Which::Unique),
            Quick::AbsentSites => (Format::Names, Which::Absent),
        }
    }
}

/// What a list was asked for: the dialog's four choices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Choice {
    format: Format,
    which: Which,
    columns: usize,
    /// "Selected DNA only": the bases of the selection, when there is one (the whole sequence otherwise).
    selection_only: bool,
}

/// What a list was made from. When any of it changes the list is out of date: the sequence (told apart by `DocView::id`, and
/// edited: `rev`), its topology and its Dam/Dcm switch, which change the sites without an edit, and for a list of the selected
/// enzymes the choice in the Enzyme Selector.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Basis {
    doc: u64,
    rev: u64,
    circular: bool,
    methylated: bool,
    chosen: Option<BTreeSet<String>>,
}

impl Basis {
    fn of(dv: &DocView, which: Which, selected: &BTreeSet<String>) -> Basis {
        Basis {
            doc: dv.id,
            rev: dv.rev,
            circular: dv.doc.circular,
            methylated: dv.doc.methylated,
            chosen: (which == Which::Selected).then(|| selected.clone()),
        }
    }
}

/// A list that was made: what it was made from, what it was asked for, how many enzymes it has, the text of the original's
/// window and its size in lines and characters. `layout` stays the same while the choices do, so that the window keeps the
/// place and the size the person gave it; another kind of list opens at the size of its own text.
struct ListRun {
    basis: Basis,
    choice: Choice,
    layout: u32,
    shown: usize,
    /// The lines of the text and the characters of the longest.
    shape: (usize, usize),
    title: String,
    report: String,
}

/// The frame, scroll bars and buttons around the text of a result window, and the least it opens at.
const WINDOW_EXTRA: [f32; 2] = [56.0, 140.0];
const MIN_WINDOW: [f32; 2] = [440.0, 260.0];
/// How much of the screen an opening window leaves free.
const SCREEN_ROOM: f32 = 80.0;

pub struct EnzListWin {
    pub dialog_open: bool,
    pub result_open: bool,
    pub format: Format,
    pub which: Which,
    pub columns: usize,
    pub selection_only: bool,
    run: Option<ListRun>,
}

impl Default for EnzListWin {
    fn default() -> Self {
        EnzListWin {
            dialog_open: false,
            result_open: false,
            format: Format::Names,
            which: Which::All,
            columns: DEFAULT_COLUMNS,
            selection_only: false,
            run: None,
        }
    }
}

/// "1 enzyme", "116 enzymes".
fn count_of(n: usize) -> String {
    let word = if n == 1 { "enzyme" } else { "enzymes" };
    format!("{} {word}", commify(n))
}

/// The lines of a text and the characters of its longest line.
fn shape_of(text: &str) -> (usize, usize) {
    let longest = text.lines().map(|l| l.chars().count()).max();
    (text.lines().count(), longest.unwrap_or(0))
}

/// The enzymes with their sites on the sequence as the Enzyme Selector counts them (Dam/Dcm switch included), numbered from 1.
fn entries_of(dv: &mut DocView, enzymes: &[Enzyme]) -> Vec<Entry> {
    let sites = dv.sites(enzymes);
    enzymes
        .iter()
        .zip(sites)
        .map(|(e, s)| Entry {
            name: e.name.clone(),
            pattern: e.pattern_text(),
            sites: s.iter().map(|site| site.start + 1).collect(),
        })
        .collect()
}

impl App {
    /// List Enzymes…: the dialog, with the format, columns and list of the last time. "DNA" starts from the Selection Only
    /// switch each time (L16339). With no sequence open there is nothing to list, and the dialog does not open.
    pub fn open_list_dialog(&mut self) {
        if self.current_doc().is_none() {
            return self.say("Open a sequence first.");
        }
        self.enz_list.selection_only = self.prefs.selection_only;
        self.enz_list.dialog_open = true;
    }

    /// A Quick List: its own format and enzymes, the columns of the dialog, and the DNA the Selection Only switch says.
    pub fn run_quick_list(&mut self, quick: Quick) {
        let (format, which) = quick.format_and_list();
        let columns = self.enz_list.columns;
        self.make_list(Choice {
            format,
            which,
            columns,
            selection_only: self.prefs.selection_only,
        });
    }

    /// What the dialog's OK does: list with its choices and close it (the original destroys it after listing, L16365-16368).
    pub(super) fn list_dialog_ok(&mut self) {
        let w = &self.enz_list;
        let choice = Choice {
            format: w.format,
            which: w.which,
            columns: w.columns,
            selection_only: w.selection_only,
        };
        self.enz_list.dialog_open = false;
        self.make_list(choice);
    }

    /// Make the list for the sequence in front and show it.
    fn make_list(&mut self, choice: Choice) {
        let Some(i) = self.current_doc() else {
            return self.say("Open a sequence first.");
        };
        let dv = &mut self.docs[i];
        // no selection is the same as all of the DNA (L16401)
        let span = dv
            .selection()
            .filter(|_| choice.selection_only)
            .map(|(start, end)| Span { start, end });
        let entries = entries_of(dv, &self.enzymes);
        let chosen = enzlist::choose(&entries, choice.which, &self.selected, span);
        let body = enzlist::text(&chosen, choice.format, choice.columns, span);
        let (len, name) = (dv.doc.len(), dv.doc.name.clone());
        let report = enzlist::report(PROGRAM, &name, choice.which, span, len, &body);
        let shown = chosen.len();
        let title = format!("{name} {}", choice.which.title());
        let message = format!("{}: {}.", choice.which.title(), count_of(shown));
        let layout = match &self.enz_list.run {
            Some(before) if before.choice == choice => before.layout,
            Some(before) => before.layout + 1,
            None => 0,
        };
        self.enz_list.run = Some(ListRun {
            basis: Basis::of(dv, choice.which, &self.selected),
            choice,
            layout,
            shown,
            shape: shape_of(&report),
            title,
            report,
        });
        self.enz_list.result_open = true;
        self.say(message);
    }

    /// Both windows; `typing` is whether a text field had the keyboard when the frame began.
    pub fn enz_list_windows(&mut self, ctx: &egui::Context, typing: bool) {
        self.list_dialog_window(ctx, typing);
        self.list_result_window(ctx);
    }

    fn list_dialog_window(&mut self, ctx: &egui::Context, typing: bool) {
        if !self.enz_list.dialog_open {
            return;
        }
        let look = design::look(self.prefs.dark);
        let (mut open, mut ok, mut cancel) = (true, false, false);
        design::window("List Enzymes")
            .open(&mut open)
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                (ok, cancel) = self.list_dialog_body(ui, &look, typing);
            });
        self.enz_list.dialog_open = open && !cancel;
        if ok {
            self.list_dialog_ok();
        }
    }

    /// The dialog's four choices and its buttons; (OK, Cancel) as pressed, or as Return and Escape ask for them.
    fn list_dialog_body(&mut self, ui: &mut Ui, look: &Look, typing: bool) -> (bool, bool) {
        let w = &mut self.enz_list;
        ui.horizontal(|ui| {
            ui.label("Format:");
            egui::ComboBox::from_id_salt("enz_list_format")
                .selected_text(w.format.label())
                .show_ui(ui, |ui| {
                    for f in Format::ALL {
                        ui.selectable_value(&mut w.format, f, f.label());
                    }
                });
            // the two formats that write a line to an enzyme have no columns
            ui.add_enabled_ui(w.format.uses_columns(), |ui| {
                ui.label("Columns:");
                egui::ComboBox::from_id_salt("enz_list_columns")
                    .selected_text(w.columns.to_string())
                    .show_ui(ui, |ui| {
                        for c in 1..=MAX_COLUMNS {
                            ui.selectable_value(&mut w.columns, c, c.to_string());
                        }
                    });
            });
        });
        ui.horizontal(|ui| {
            ui.label("Enzymes:");
            for which in Which::ALL {
                ui.radio_value(&mut w.which, which, which.label());
            }
        });
        ui.horizontal(|ui| {
            ui.label("DNA:");
            ui.radio_value(&mut w.selection_only, false, "All");
            ui.radio_value(&mut w.selection_only, true, "Selected DNA only");
        });
        let (mut ok, mut cancel) = (false, false);
        ui.horizontal(|ui| {
            ok = design::primary_button(ui, look, "OK").clicked();
            cancel = design::ghost_button(ui, look, "Cancel").clicked();
        });
        if self.dialog_has_the_keys(ui, typing) {
            // the original binds Return to OK and its <<Cancel>> (Escape) to Cancel (L16312-16313, L16347)
            ok |= ui.input_mut(|i| i.consume_key(Modifiers::NONE, Key::Enter));
            cancel |= ui.input_mut(|i| i.consume_key(Modifiers::NONE, Key::Escape));
        }
        (ok, cancel)
    }

    /// Whether Return and Escape are meant for the dialog: not while the command palette is open, nor while a text field had
    /// the keyboard (Find's, say), and only when the dialog is the window in front.
    pub(super) fn dialog_has_the_keys(&self, ui: &Ui, typing: bool) -> bool {
        !self.palette.is_open() && !typing && ui.ctx().top_layer_id() == Some(ui.layer_id())
    }

    fn list_result_window(&mut self, ctx: &egui::Context) {
        if !self.enz_list.result_open {
            return;
        }
        let Some(run) = &self.enz_list.run else {
            self.enz_list.result_open = false;
            return;
        };
        let look = design::look(self.prefs.dark);
        let live = self.current_doc().is_some_and(|i| {
            Basis::of(&self.docs[i], run.choice.which, &self.selected) == run.basis
        });
        let (mut open, mut again, mut copy) = (true, false, false);
        design::window(run.title.as_str())
            .id(Id::new(("enz_list_result", run.layout)))
            .open(&mut open)
            .default_size(natural_size(ctx, run.shape))
            .show(ctx, |ui| {
                (again, copy) = list_head(ui, &look, run, live);
                list_text(ui, &run.report);
            });
        if copy {
            ctx.copy_text(run.report.clone());
        }
        let choice = run.choice;
        self.enz_list.result_open = open;
        if again {
            self.make_list(choice);
        }
    }
}

/// The size the original gives its window: as wide as the longest line and as high as the text (`textbox_geometry`, L16817),
/// but not more than the screen has, so that every column of a list is there when the window opens. `shape` is the lines of
/// the text and the characters of the longest.
pub(super) fn natural_size(ctx: &egui::Context, shape: (usize, usize)) -> [f32; 2] {
    let style = ctx.global_style();
    let font = egui::TextStyle::Monospace.resolve(&style);
    let (advance, row) = ctx.fonts_mut(|f| (f.glyph_width(&font, 'M'), f.row_height(&font)));
    let screen = ctx.content_rect().size();
    let fit =
        |wanted: f32, least: f32, room: f32| wanted.clamp(least, (room - SCREEN_ROOM).max(least));
    let row = row + style.spacing.item_spacing.y;
    [
        fit(
            shape.1 as f32 * advance + WINDOW_EXTRA[0],
            MIN_WINDOW[0],
            screen.x,
        ),
        fit(
            shape.0 as f32 * row + WINDOW_EXTRA[1],
            MIN_WINDOW[1],
            screen.y,
        ),
    ]
}

/// How many enzymes, the buttons and the note that the sequence has changed; (list again, copy as text).
fn list_head(ui: &mut Ui, look: &Look, run: &ListRun, live: bool) -> (bool, bool) {
    let mut pressed = (false, false);
    ui.horizontal(|ui| {
        design::chip(ui, look, &count_of(run.shown), Tone::Accent);
        pressed.0 = design::primary_button(ui, look, "List again").clicked();
        pressed.1 = design::ghost_button(ui, look, "Copy as text").clicked();
    });
    if !live {
        let text = "The sequence has changed since this list was made; list again to refresh it.";
        ui.label(RichText::new(text).small().color(look.warn));
    }
    ui.separator();
    pressed
}

/// The text of the original's window, a line to a row, scrolling both ways: the lines are not wrapped, as the original's are not.
fn list_text(ui: &mut Ui, report: &str) {
    let lines: Vec<&str> = report.lines().collect();
    let row = ui.text_style_height(&egui::TextStyle::Monospace);
    ScrollArea::both()
        .auto_shrink([false, false])
        .show_rows(ui, row, lines.len(), |ui, rows| {
            for line in &lines[rows] {
                // an empty line still takes a row
                let shown = if line.is_empty() { " " } else { line };
                ui.add(
                    egui::Label::new(RichText::new(shown).monospace())
                        .wrap_mode(TextWrapMode::Extend),
                );
            }
        });
}
