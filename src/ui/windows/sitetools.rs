//! Silent Sites and Find Diagnostic Site (Enzymes menu): two result windows. Each takes the sequence in front, the stretch
//! of it that is selected (or all of it) and the enzymes chosen in the Enzyme Selector, runs the search of `seq::silent` /
//! `seq::diagnostic` once, and lists what it found; a click on a row selects those bases in the sequence.

use crate::seq::diagnostic::{self, Found, MAX_MISMATCHES};
use crate::seq::dna::GENETIC_CODES;
use crate::seq::enzymes::Enzyme;
use crate::seq::silent::{self, Hit, Query, MAX_PLACES};
use crate::ui::app::App;
use crate::ui::design::{self, Look, Tone};
use crate::ui::map::commify;
use crate::ui::platform::Workspace;
use eframe::egui::{
    self, text::LayoutJob, Color32, CursorIcon, FontId, Id, Label, Rect, RichText, ScrollArea,
    Sense, TextFormat, Ui,
};

#[cfg(test)]
mod tests;

/// Rows a list draws; "Copy as text" has all of them. (HelixLoom's own limit.)
const ROWS_DRAWN: usize = 1000;
const MONO: f32 = 13.0;
/// The width of the first column of a row, in characters.
const LEAD: usize = 10;
const CLEAR: Color32 = Color32::TRANSPARENT;
/// The heading of the copied text where the original writes the date.
const PROGRAM: &str = "HelixLoom";

/// The sequence a search was made on (told apart by `DocView::id`, not by its name), and when (`rev`: a later edit makes its
/// rows stale), and the stretch of it that was searched, `from..to` from 0.
struct Stretch {
    doc: u64,
    name: String,
    rev: u64,
    from: usize,
    to: usize,
}

struct SilentRun {
    at: Stretch,
    enzymes: usize,
    all: bool,
    result: silent::Search,
}

struct DiagnosticRun {
    at: Stretch,
    enzymes: usize,
    absent: bool,
    mismatches: usize,
    dna: Vec<u8>,
    result: diagnostic::Search,
}

pub struct SiteToolsWin {
    pub silent_open: bool,
    pub diagnostic_open: bool,
    /// Find Diagnostic Site: search the enzymes that are absent from the sequence, not the selected ones.
    pub absent: bool,
    pub mismatches: usize,
    /// The bases a click on a row selected last. A search that finds exactly these selected does not take them for a stretch
    /// the person chose: it goes on with the stretch of the search whose row was clicked.
    picked: Option<(usize, usize)>,
    silent: Option<SilentRun>,
    diagnostic: Option<DiagnosticRun>,
}

impl Default for SiteToolsWin {
    fn default() -> Self {
        SiteToolsWin {
            silent_open: false,
            diagnostic_open: false,
            absent: false,
            mismatches: 1,
            picked: None,
            silent: None,
            diagnostic: None,
        }
    }
}

fn line(pieces: Vec<(String, Color32, Color32)>) -> LayoutJob {
    let mut job = LayoutJob::default();
    for (text, color, background) in pieces {
        let format = TextFormat {
            font_id: FontId::monospace(MONO),
            color,
            background,
            ..TextFormat::default()
        };
        job.append(&text, 0.0, format);
    }
    job
}

fn lead(text: &str) -> String {
    format!("{text:<LEAD$}")
}

/// The bases of a variant or a site with the changed ones in colour, in threes when `threes`.
fn marked(look: &Look, first: String, marks: &[(char, bool)], threes: bool) -> LayoutJob {
    let mut pieces = vec![(first, look.text, CLEAR)];
    for (i, &(c, changed)) in marks.iter().enumerate() {
        if threes && i > 0 && i % 3 == 0 {
            pieces.push((" ".into(), look.text, CLEAR));
        }
        let (color, background) = if changed {
            (look.bad, look.bad.gamma_multiply(0.15))
        } else {
            (look.text, CLEAR)
        };
        pieces.push((c.to_string(), color, background));
    }
    line(pieces)
}

/// One line of a row; where it is, so that the row can be clicked as a whole.
fn text_line(ui: &mut Ui, job: LayoutJob) -> Rect {
    ui.add(Label::new(job).selectable(false)).rect
}

/// A row as one target, under an id of its own so that a test can find it; true when it was clicked, and `live`.
fn row_target(ui: &mut Ui, id: Id, area: Rect, live: bool) -> bool {
    if !live {
        return false;
    }
    let response = ui
        .interact(area, id, Sense::click())
        .on_hover_cursor(CursorIcon::PointingHand)
        .on_hover_text("Select these bases in the sequence");
    response.clicked()
}

/// The ids of a row of each list.
pub fn silent_row_id(index: usize) -> Id {
    Id::new(("silent_row", index))
}

pub fn diagnostic_row_id(index: usize) -> Id {
    Id::new(("diagnostic_row", index))
}

fn hit_row(ui: &mut Ui, look: &Look, (index, h): (usize, &Hit), live: bool) -> bool {
    let codons = line(vec![
        (lead(&h.position.to_string()), look.dim, CLEAR),
        (silent::in_threes(&h.codons), look.text, CLEAR),
    ]);
    let mut area = text_line(ui, codons);
    for v in &h.variants {
        let marks = silent::marks(v, &h.codons);
        area = area.union(text_line(ui, marked(look, lead(""), &marks, true)));
    }
    let site = vec![
        (lead(&h.enzyme), look.accent, CLEAR),
        (h.layout.clone(), look.dim, CLEAR),
    ];
    area = area.union(text_line(ui, line(site)));
    let residues = vec![
        (lead(""), look.dim, CLEAR),
        (format!(" {}", h.residues), look.dim, CLEAR),
    ];
    area = area.union(text_line(ui, line(residues)));
    row_target(ui, silent_row_id(index), area, live)
}

fn found_row(
    ui: &mut Ui,
    look: &Look,
    (index, f): (usize, &Found),
    run: &DiagnosticRun,
    live: bool,
) -> bool {
    let end = (f.index + f.pattern.len()).min(run.dna.len());
    let shown: String = run.dna[f.index..end]
        .iter()
        .map(|&c| char::from(c))
        .collect();
    let dna = line(vec![
        (
            lead(&(f.index + run.at.from + 1).to_string()),
            look.dim,
            CLEAR,
        ),
        (shown, look.text, CLEAR),
    ]);
    let mut area = text_line(ui, dna);
    let marks = diagnostic::marks(f, &run.dna);
    area = area.union(text_line(ui, marked(look, lead(""), &marks, false)));
    let site = vec![
        (lead(&f.enzyme), look.accent, CLEAR),
        (
            String::from_utf8_lossy(&f.pattern).into_owned(),
            look.dim,
            CLEAR,
        ),
    ];
    area = area.union(text_line(ui, line(site)));
    row_target(ui, diagnostic_row_id(index), area, live)
}

fn note(ui: &mut Ui, text: &str, color: Color32) {
    ui.label(RichText::new(text).small().color(color));
}

/// The rows of a list, one block each, with room between; the position of the bases of the row that was clicked.
fn rows<T>(
    ui: &mut Ui,
    items: &[T],
    mut draw: impl FnMut(&mut Ui, (usize, &T)) -> Option<(usize, usize)>,
) -> Option<(usize, usize)> {
    let mut picked = None;
    ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.spacing_mut().item_spacing.y = 1.0;
            for item in items.iter().enumerate().take(ROWS_DRAWN) {
                if let Some(range) = draw(ui, item) {
                    picked = Some(range);
                }
                ui.add_space(8.0);
            }
        });
    picked
}

fn drawn_note(ui: &mut Ui, look: &Look, total: usize) {
    if total > ROWS_DRAWN {
        let text = format!(
            "The first {} of {} are drawn; Copy as text has all of them.",
            commify(ROWS_DRAWN),
            commify(total)
        );
        note(ui, &text, look.dim);
    }
}

impl App {
    /// The sequence in front and the stretch the tools look at: the selection, or all of it. A selection that is only the
    /// bases a click on a row selected is not one the person chose: the stretch of the search it came from (`previous`)
    /// goes on, as long as that search is still current.
    fn site_target(&self, previous: Option<&Stretch>) -> Option<(usize, Stretch)> {
        let i = self.current_doc()?;
        let dv = &self.docs[i];
        let len = dv.doc.len();
        let selection = dv.selection();
        let (from, to) = match previous {
            Some(p)
                if selection.is_some()
                    && selection == self.site_win.picked
                    && self.is_current(previous) =>
            {
                (p.from, p.to)
            }
            _ => selection.map_or((0, len), |(a, b)| (a.min(len), b.min(len))),
        };
        let at = Stretch {
            doc: dv.id,
            name: dv.doc.name.clone(),
            rev: dv.rev,
            from,
            to,
        };
        Some((i, at))
    }

    /// The enzymes chosen in the Enzyme Selector, as places in the enzyme list.
    fn chosen_enzymes(&self) -> Vec<usize> {
        let named = |e: &Enzyme| self.selected.contains(&e.name);
        let all = self.enzymes.iter().enumerate();
        all.filter(|(_, e)| named(e)).map(|(k, _)| k).collect()
    }

    /// Whether `at` is the sequence in front, unchanged since: its rows can still be selected.
    fn is_current(&self, at: Option<&Stretch>) -> bool {
        let (Some(at), Some(i)) = (at, self.current_doc()) else {
            return false;
        };
        self.docs[i].id == at.doc && self.docs[i].rev == at.rev
    }

    /// Run Silent Sites on the sequence in front and show the result.
    pub fn run_silent_sites(&mut self) {
        let previous = self.site_win.silent.as_ref().map(|r| &r.at);
        let Some((i, at)) = self.site_target(previous) else {
            return self.say("Open a sequence first.");
        };
        let chosen = self.chosen_enzymes();
        let all = chosen.is_empty();
        let indexes: Vec<usize> = if all {
            (0..self.enzymes.len()).collect()
        } else {
            chosen
        };
        let sites = self.docs[i].sites(&self.enzymes);
        // in order, as the search looks them up by halving
        let existing: Vec<Vec<usize>> = indexes
            .iter()
            .map(|&k| {
                let mut starts: Vec<usize> = sites[k].iter().map(|s| s.start + 1).collect();
                starts.sort_unstable();
                starts.dedup();
                starts
            })
            .collect();
        let picked: Vec<&Enzyme> = indexes.iter().map(|&k| &self.enzymes[k]).collect();
        let result = silent::search(&Query {
            seq: &self.docs[i].doc.seq,
            region: (at.from, at.to),
            table: GENETIC_CODES[self.prefs.code].table,
            enzymes: &picked,
            existing: &existing,
        });
        let message = format!("Silent sites: {} places.", result.hits.len());
        self.site_win.silent = Some(SilentRun {
            at,
            enzymes: indexes.len(),
            all,
            result,
        });
        self.site_win.silent_open = true;
        self.say(message);
    }

    /// Run Find Diagnostic Site with the window's choices on the sequence in front and show the result.
    pub fn run_diagnostic_sites(&mut self) {
        let previous = self.site_win.diagnostic.as_ref().map(|r| &r.at);
        let Some((i, at)) = self.site_target(previous) else {
            return self.say("Open a sequence first.");
        };
        let (absent, mismatches) = (self.site_win.absent, self.site_win.mismatches);
        let sites = self.docs[i].sites(&self.enzymes);
        let indexes: Vec<usize> = if absent {
            (0..self.enzymes.len())
                .filter(|&k| sites[k].is_empty())
                .collect()
        } else {
            self.chosen_enzymes()
        };
        let picked: Vec<&Enzyme> = indexes.iter().map(|&k| &self.enzymes[k]).collect();
        let dna = self.docs[i].doc.slice(at.from, at.to).to_vec();
        let result = diagnostic::find(&dna, &picked, mismatches);
        let message = format!("Diagnostic sites: {} places.", result.found.len());
        self.site_win.diagnostic = Some(DiagnosticRun {
            at,
            enzymes: indexes.len(),
            absent,
            mismatches,
            dna,
            result,
        });
        self.say(message);
    }

    /// Select `from..to` in the sequence in front and show it.
    fn select_found(&mut self, from: usize, to: usize) {
        let Some(i) = self.current_doc() else { return };
        self.workspace = Workspace::Sequences;
        let dv = &mut self.docs[i];
        dv.select(from, to);
        dv.scroll_to = Some(from);
        let picked = dv.selection();
        self.site_win.picked = picked;
    }

    pub fn site_windows(&mut self, ctx: &egui::Context) {
        self.silent_window(ctx);
        self.diagnostic_window(ctx);
    }

    fn silent_window(&mut self, ctx: &egui::Context) {
        if !self.site_win.silent_open {
            return;
        }
        let look = design::look(self.prefs.dark);
        let live = self.is_current(self.site_win.silent.as_ref().map(|r| &r.at));
        let (mut open, mut again, mut copy, mut pick) = (true, false, false, None);
        design::window("Silent Sites")
            .open(&mut open)
            .default_size([600.0, 540.0])
            .show(ctx, |ui| {
                ui.label(
                    "Where a restriction site could be made without changing the protein. A click selects the codons.",
                );
                let Some(run) = &self.site_win.silent else {
                    return;
                };
                (again, copy) = silent_head(ui, &look, run, live);
                pick = rows(ui, &run.result.hits, |ui, (n, h)| {
                    hit_row(ui, &look, (n, h), live).then_some(h.codons_at)
                });
            });
        self.site_win.silent_open = open;
        if copy {
            ctx.copy_text(self.silent_text());
        }
        if again {
            self.run_silent_sites();
        }
        if let Some((a, b)) = pick {
            self.select_found(a, b);
        }
    }

    fn silent_text(&self) -> String {
        let Some(run) = &self.site_win.silent else {
            return String::new();
        };
        let (from, to) = (run.at.from + 1, run.at.to);
        silent::report(PROGRAM, &run.at.name, from, to, &run.result.hits)
    }

    fn diagnostic_window(&mut self, ctx: &egui::Context) {
        if !self.site_win.diagnostic_open {
            return;
        }
        let look = design::look(self.prefs.dark);
        let live = self.is_current(self.site_win.diagnostic.as_ref().map(|r| &r.at));
        let (mut open, mut find, mut copy, mut pick) = (true, false, false, None);
        design::window("Diagnostic Sites")
            .open(&mut open)
            .default_size([600.0, 540.0])
            .show(ctx, |ui| {
                ui.label(
                    "Where an enzyme's site is a few bases away from the DNA: a primer with those changes makes a site that tells the changed copy from the original.",
                );
                find = self.diagnostic_options(ui, &look);
                let Some(run) = &self.site_win.diagnostic else {
                    return;
                };
                copy = diagnostic_head(ui, &look, run, live);
                pick = rows(ui, &run.result.found, |ui, (n, f)| {
                    let range = (run.at.from + f.index, run.at.from + f.index + f.pattern.len());
                    found_row(ui, &look, (n, f), run, live).then_some(range)
                });
            });
        self.site_win.diagnostic_open = open;
        if copy {
            ctx.copy_text(self.diagnostic_text());
        }
        if find {
            self.run_diagnostic_sites();
        }
        if let Some((a, b)) = pick {
            self.select_found(a, b);
        }
    }

    /// The window's two choices and the button that searches; true when it was pressed.
    fn diagnostic_options(&mut self, ui: &mut Ui, look: &Look) -> bool {
        ui.horizontal(|ui| {
            ui.label("Enzymes:");
            ui.radio_value(&mut self.site_win.absent, false, "Selected enzymes");
            ui.radio_value(&mut self.site_win.absent, true, "Absent from this sequence");
            if ui.small_button("Enzyme Selector…").clicked() {
                self.enzyme_win.open = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("Mismatches:");
            for n in 1..=MAX_MISMATCHES {
                ui.radio_value(&mut self.site_win.mismatches, n, n.to_string());
            }
        });
        let find = design::primary_button(ui, look, "Find").clicked();
        ui.separator();
        find
    }

    fn diagnostic_text(&self) -> String {
        let Some(run) = &self.site_win.diagnostic else {
            return String::new();
        };
        let shown = (run.at.from + 1, run.at.to);
        let head = (PROGRAM, run.at.name.as_str());
        diagnostic::report(head, shown, run.at.from, &run.dna, &run.result.found)
    }
}

/// What a Silent Sites search covered, its buttons, and the notes; (search again, copy).
fn silent_head(ui: &mut Ui, look: &Look, run: &SilentRun, live: bool) -> (bool, bool) {
    let covered = if run.all {
        "No enzymes are selected, so all of them were searched.".to_string()
    } else {
        format!("Searched the {} selected enzymes.", run.enzymes)
    };
    let stretch = format!(
        "{}: bases {}–{}. {covered}",
        run.at.name,
        commify(run.at.from + 1),
        commify(run.at.to)
    );
    ui.label(RichText::new(stretch).small().color(look.dim));
    let (again, copy) = head_buttons(ui, look, "Search again", run.result.hits.len(), "places");
    if !live {
        let text = "The sequence has changed since this search; search again to refresh it.";
        note(ui, text, look.warn);
    }
    if run.result.truncated {
        let text = format!(
            "Stopped after {} places; select fewer enzymes or a shorter stretch.",
            commify(MAX_PLACES)
        );
        note(ui, &text, look.warn);
    }
    drawn_note(ui, look, run.result.hits.len());
    if run.result.hits.is_empty() {
        ui.label("No sites found.");
    }
    ui.separator();
    (again, copy)
}

/// The same for Find Diagnostic Site; true when Copy as text was pressed.
fn diagnostic_head(ui: &mut Ui, look: &Look, run: &DiagnosticRun, live: bool) -> bool {
    let which = if run.absent {
        "that are absent from the sequence"
    } else {
        "that were selected"
    };
    let stretch = format!(
        "{}: bases {}–{}. Searched the {} enzymes {which}, up to {} mismatch{}.",
        run.at.name,
        commify(run.at.from + 1),
        commify(run.at.to),
        run.enzymes,
        run.mismatches,
        if run.mismatches == 1 { "" } else { "es" }
    );
    ui.label(RichText::new(stretch).small().color(look.dim));
    let (_, copy) = head_buttons(ui, look, "", run.result.found.len(), "places");
    if run.enzymes == 0 {
        let text =
            "No enzymes to search: choose some in the Enzyme Selector, or search the absent ones.";
        note(ui, text, look.warn);
    }
    if !live {
        let text = "The sequence has changed since this search; press Find to refresh it.";
        note(ui, text, look.warn);
    }
    if run.result.truncated {
        let text = "Stopped at the limit; use fewer mismatches or fewer enzymes.";
        note(ui, text, look.warn);
    }
    drawn_note(ui, look, run.result.found.len());
    if run.result.found.is_empty() {
        ui.label("No sites found.");
    }
    ui.separator();
    copy
}

/// The count of what was found and the buttons; (the first button, Copy as text). A blank `again` draws no first button.
fn head_buttons(ui: &mut Ui, look: &Look, again: &str, count: usize, what: &str) -> (bool, bool) {
    let mut pressed = (false, false);
    ui.horizontal(|ui| {
        design::chip(
            ui,
            look,
            &format!("{} {what}", commify(count)),
            Tone::Accent,
        );
        if !again.is_empty() {
            pressed.0 = design::primary_button(ui, look, again).clicked();
        }
        pressed.1 = design::ghost_button(ui, look, "Copy as text").clicked();
    });
    pressed
}
