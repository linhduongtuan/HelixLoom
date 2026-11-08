//! Text Map (Enzymes menu, ⇧⌘T): the dialog that asks what the page holds and in which order, and the window that shows the page
//! `seq::textmap` makes, with the colours of the features on it. A click on a line of DNA puts the caret on that base of the
//! sequence, a Shift-click selects up to it (the original's `alignlink2`, L20148-20151).

mod dialog;
mod view;

#[cfg(test)]
mod tests;

use crate::seq::features::Feature;
use crate::seq::textmap::{self, EnzymeSites, MapFeature, Region, Request, Settings, TextMap};
use crate::ui::app::App;
use crate::ui::design;
use crate::ui::state::DocView;
use eframe::egui::{self, Id};
use std::collections::BTreeSet;
use view::Painted;

#[cfg(test)]
use dialog::move_frame;

/// What the dialog holds, and the last page made.
#[derive(Default)]
pub struct TextMapWin {
    pub dialog_open: bool,
    pub result_open: bool,
    pub settings: Settings,
    /// "Selected Enzymes Only": the enzymes picked in the Enzyme Selector, not all of them (`info(enz_text_enz)`).
    pub selected_only: bool,
    /// "Selected DNA only": the bases of the selection, when there is one.
    pub selection_only: bool,
    run: Option<MapRun>,
}

/// What a page was made from. When any of it changes the page is out of date: the sequence (told apart by `DocView::id`, and
/// edited: `rev`), its topology and Dam/Dcm switch, which change the sites without an edit, and the choice in the Enzyme
/// Selector when the page is about the selected enzymes.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Basis {
    doc: u64,
    rev: u64,
    circular: bool,
    methylated: bool,
    chosen: Option<BTreeSet<String>>,
}

impl Basis {
    fn of(dv: &DocView, selected_only: bool, selected: &BTreeSet<String>) -> Basis {
        Basis {
            doc: dv.id,
            rev: dv.rev,
            circular: dv.doc.circular,
            methylated: dv.doc.methylated,
            chosen: selected_only.then(|| selected.clone()),
        }
    }
}

/// A page that was made, with what it was made from and what is needed to show it.
struct MapRun {
    basis: Basis,
    /// Stays the same while the choices do, so that the window keeps the place and size the person gave it.
    layout: u32,
    title: String,
    page: TextMap,
    painted: Painted,
    selected_only: bool,
    /// How many bases the sequence has: a click is on a base of it.
    seq_len: usize,
}

/// The features as the page needs them. HelixLoom has no hidden features, so none is.
fn map_features(features: &[Feature]) -> Vec<MapFeature> {
    features
        .iter()
        .map(|f| MapFeature {
            name: f.name.clone(),
            kind: f.kind.clone(),
            spans: f.spans.clone(),
            reverse: f.reverse,
            hidden: false,
            fwd: f.fwd_color,
            rev: f.rev_color,
        })
        .collect()
}

/// What to do about an error the original stops with, where there is something to do.
fn hint(why: &str) -> &'static str {
    if why.contains("textline") {
        " (The second strand is made from the line of DNA before it: put the 2nd Strand box below the DNA box.)"
    } else {
        ""
    }
}

/// The date and the time the way the original writes them (`%a %b %d, %Y %k:%M %Z`); the program knows no time zone, so it is UTC.
fn clock_line() -> String {
    let secs = web_time::SystemTime::now()
        .duration_since(web_time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    let days = (secs / 86_400) as i64;
    let (hour, minute) = ((secs % 86_400) / 3_600, (secs % 3_600) / 60);
    // 1 January 1970 was a Thursday
    const DAYS: [&str; 7] = ["Thu", "Fri", "Sat", "Sun", "Mon", "Tue", "Wed"];
    const MONTHS: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    let weekday = DAYS[days.rem_euclid(7) as usize];
    format!(
        "{weekday} {} {day:02}, {year} {hour:>2}:{minute:02} UTC",
        MONTHS[(month - 1) as usize]
    )
}

impl App {
    /// Text Map…: the dialog, with what it held last time. "DNA" starts from the Selection Only switch each time, as the
    /// original's does (L19440). With no sequence open there is nothing to map.
    pub fn open_text_map_dialog(&mut self) {
        if self.current_doc().is_none() {
            return self.say("Open a sequence first.");
        }
        self.text_map.selection_only = self.prefs.selection_only;
        self.text_map.dialog_open = true;
    }

    /// Make the page for the sequence in front with what the dialog holds, and show it.
    pub fn make_text_map(&mut self) {
        let Some(i) = self.current_doc() else {
            return self.say("Open a sequence first.");
        };
        let w = &self.text_map;
        let (settings, selected_only, selection_only) =
            (w.settings.clone(), w.selected_only, w.selection_only);
        let dv = &mut self.docs[i];
        let enzymes = enzyme_sites(dv, &self.enzymes, &self.selected, selected_only);
        let sites: Vec<EnzymeSites> = enzymes
            .iter()
            .map(|(name, sites)| EnzymeSites { name, sites })
            .collect();
        let features = map_features(&dv.doc.features);
        let selection = dv.selection();
        let region = match selection.filter(|_| selection_only) {
            Some((a, b)) => Region::Part(a, b),
            None => Region::All,
        };
        let request = Request {
            seq: &dv.doc.seq,
            circular: dv.doc.circular,
            title: &dv.doc.name,
            date: &clock_line(),
            region,
            selection,
            row_width: dv.per_row,
            enzymes: &sites,
            features: &features,
            code: textmap::STANDARD_CODE,
            settings: &settings,
        };
        let page = match textmap::text_map(&request) {
            Ok(Some(page)) => page,
            Ok(None) => return self.say("There is nothing to map."),
            Err(why) => return self.say(format!("Text Map: {why}{}", hint(&why.0))),
        };
        let title = format!("{} Text Map", dv.doc.name);
        let basis = Basis::of(dv, selected_only, &self.selected);
        let layout = self.text_map.run.as_ref().map_or(0, |r| r.layout);
        let painted = Painted::of(&page);
        self.text_map.run = Some(MapRun {
            basis,
            layout,
            title,
            page,
            painted,
            selected_only,
            seq_len: self.docs[i].doc.len(),
        });
        self.text_map.result_open = true;
        self.say("Text Map made.");
    }

    /// Both windows; `typing` is whether a text field had the keyboard when the frame began.
    pub fn text_map_windows(&mut self, ctx: &egui::Context, typing: bool) {
        self.text_map_dialog_window(ctx, typing);
        self.text_map_result_window(ctx);
    }

    fn text_map_result_window(&mut self, ctx: &egui::Context) {
        if !self.text_map.result_open {
            return;
        }
        let Some(run) = &self.text_map.run else {
            self.text_map.result_open = false;
            return;
        };
        let look = design::look(self.prefs.dark);
        let live = self
            .docs
            .iter()
            .any(|d| Basis::of(d, run.selected_only, &self.selected) == run.basis);
        let (mut open, mut again, mut copy) = (true, false, false);
        let mut clicked = None;
        let shape = run.painted.shape();
        design::window(run.title.as_str())
            .id(Id::new(("text_map_result", run.layout)))
            .open(&mut open)
            .default_size(super::enzlist::natural_size(ctx, shape))
            .show(ctx, |ui| {
                (again, copy) = view::head(ui, &look, live);
                clicked = view::rows(ui, &look, run);
            });
        if copy {
            ctx.copy_text(run.page.text.clone());
        }
        let jump = clicked.and_then(|c| view::base_at(run, c.line, c.col).map(|b| (b, c.shift)));
        let doc = run.basis.doc;
        self.text_map.result_open = open;
        if let Some((base, shift)) = jump {
            self.put_caret_on(doc, base, shift);
        }
        if again {
            self.make_text_map();
        }
    }

    /// Put the caret on a base of the sequence a page was made from, or select up to it.
    fn put_caret_on(&mut self, doc: u64, base: usize, shift: bool) {
        let Some(dv) = self.docs.iter_mut().find(|d| d.id == doc) else {
            return;
        };
        if shift {
            let from = dv.anchor.unwrap_or(dv.caret);
            dv.select(from, base);
        } else {
            dv.place_caret(base);
        }
        dv.scroll_to = Some(base);
    }
}

/// The enzymes the page is about, each with its sites numbered from 1: all of them in the order of the enzyme file, or those
/// the Enzyme Selector has picked, by name (the order `enzlist` lists them in).
fn enzyme_sites(
    dv: &mut DocView,
    enzymes: &[crate::seq::enzymes::Enzyme],
    selected: &BTreeSet<String>,
    selected_only: bool,
) -> Vec<(String, Vec<usize>)> {
    let found: Vec<(String, Vec<usize>)> = enzymes
        .iter()
        .zip(dv.sites(enzymes))
        .map(|(e, s)| (e.name.clone(), s.iter().map(|x| x.start + 1).collect()))
        .collect();
    if !selected_only {
        return found;
    }
    let mut picked: Vec<(String, Vec<usize>)> = found
        .into_iter()
        .filter(|(n, _)| selected.contains(n))
        .collect();
    picked.sort_by(|a, b| a.0.cmp(&b.0));
    picked
}
