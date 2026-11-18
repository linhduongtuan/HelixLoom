//! File > Fetch from a Database…: a record of NCBI, UniProt, the PDB or the AlphaFold DB by its accession, kept in a folder of
//! the program's settings, with where it came from and its licence, and opened where it belongs.

use crate::connect::curl::Curl;
use crate::connect::{self, Fetched, Fetcher, Format, Provenance, Request, Source, Store};
use crate::pipeline::tools;
use crate::ui::app::App;
use crate::ui::design::{self, Tone};
use crate::ui::job::{Job, Poll};
use crate::ui::platform::Workspace;
use eframe::egui::{self, RichText};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Duration;

const WORKING_REDRAW: Duration = Duration::from_millis(200);

pub struct FetchWin {
    pub open: bool,
    pub source: Source,
    pub text: String,
    /// Ask the database even when a copy is kept.
    pub refresh: bool,
    pub message: String,
    /// Where the record last opened came from.
    pub last: Option<(Provenance, bool)>,
    /// The folder the records are kept in; none keeps nothing.
    pub cache_dir: Option<PathBuf>,
    /// What gets the record. None is `curl`; the tests put a table here.
    pub fetcher: Option<Arc<dyn Fetcher + Send + Sync>>,
    job: Option<Job<Result<Fetched, String>>>,
    /// The window was just opened: the accession box gets the keyboard.
    pub fresh: bool,
}

impl Default for FetchWin {
    fn default() -> Self {
        FetchWin {
            open: false,
            source: Source::NcbiNucleotide,
            text: String::new(),
            refresh: false,
            message: String::new(),
            last: None,
            cache_dir: None,
            fetcher: None,
            job: None,
            fresh: false,
        }
    }
}

impl FetchWin {
    pub fn is_running(&self) -> bool {
        self.job.is_some()
    }
}

impl App {
    pub fn fetch_window(&mut self, ctx: &egui::Context) {
        self.fetch_tick(ctx);
        if !self.fetch_win.open {
            return;
        }
        let look = self.prefs.look();
        let mut open = true;
        let mut go = false;
        design::window("Fetch from a database")
            .open(&mut open)
            .default_width(460.0)
            .resizable([true, false])
            .show(ctx, |ui| {
                let win = &mut self.fetch_win;
                ui.horizontal(|ui| {
                    ui.label("Database");
                    egui::ComboBox::from_id_salt("fetch_source")
                        .selected_text(win.source.label())
                        .show_ui(ui, |ui| {
                            for s in Source::ALL {
                                ui.selectable_value(&mut win.source, s, s.label());
                            }
                        });
                });
                let line = ui.add(
                    egui::TextEdit::singleline(&mut win.text)
                        .hint_text(format!("accession, for example {}", win.source.example()))
                        .desired_width(f32::INFINITY),
                );
                if std::mem::take(&mut win.fresh) {
                    line.request_focus();
                }
                let pressed_enter = line.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
                ui.checkbox(&mut win.refresh, "Get a fresh copy, even if one is kept here");
                ui.horizontal(|ui| {
                    let working = win.is_running();
                    let button = ui.add_enabled_ui(!working, |ui| {
                        design::primary_button(ui, &look, if working { "Asking…" } else { "Fetch" })
                    });
                    go = button.inner.clicked() || (pressed_enter && !working);
                });
                if !win.message.is_empty() {
                    ui.label(RichText::new(&win.message).color(look.dim));
                }
                if let Some((p, from_cache)) = &win.last {
                    ui.add_space(6.0);
                    design::section(ui, &look, "Where it came from");
                    ui.horizontal_wrapped(|ui| {
                        design::chip(ui, &look, &p.source, Tone::Accent);
                        design::chip(ui, &look, &p.id, Tone::Neutral);
                        let when = if *from_cache {
                            format!("kept since {}", p.fetched)
                        } else {
                            format!("fetched {}", p.fetched)
                        };
                        design::chip(ui, &look, &when, Tone::Good);
                    });
                    ui.label(RichText::new(&p.url).monospace().size(11.5).color(look.dim));
                    ui.label(RichText::new(format!("Licence: {}", p.licence)).size(12.0).color(look.dim));
                }
                ui.add_space(4.0);
                let kept = match &win.cache_dir {
                    Some(dir) => format!("Records are kept in {}.", dir.display()),
                    None => "Nothing is kept: there is no settings folder.".to_string(),
                };
                ui.label(RichText::new(kept).size(11.5).color(look.dim));
                ui.label(
                    RichText::new("Reached with curl, which most computers have. The data belong to their databases: see the licence.")
                        .size(11.5)
                        .color(look.dim),
                );
            });
        self.fetch_win.open = open;
        if go {
            self.start_fetch(ctx);
        }
    }

    fn start_fetch(&mut self, ctx: &egui::Context) {
        let win = &mut self.fetch_win;
        let request = match Request::check(win.source, &win.text) {
            Ok(r) => r,
            Err(why) => {
                win.message = why;
                return;
            }
        };
        let (dir, refresh, injected) = (win.cache_dir.clone(), win.refresh, win.fetcher.clone());
        win.message = format!("Asking {} for {}…", request.source.label(), request.id);
        win.job = Some(Job::spawn(ctx, move || {
            // a folder that cannot be made is a record not kept, not a record not fetched
            let (store, no_store) = match dir.as_deref().map(Store::open) {
                Some(Ok(s)) => (Some(s), None),
                Some(Err(why)) => (None, Some(why)),
                None => (None, None),
            };
            let mut got = match injected {
                Some(fetcher) => {
                    connect::fetch(&request, fetcher.as_ref(), store.as_ref(), refresh)
                }
                None => {
                    let curl = tools::detect(&tools::CURL).ok_or(
                        "curl is not installed here, and the databases are reached with it. Install curl and try again.",
                    )?;
                    connect::fetch(&request, &Curl::new(curl.path), store.as_ref(), refresh)
                }
            }?;
            got.warning = got.warning.or(no_store);
            Ok(got)
        }));
    }

    /// Collect an answer; once a frame.
    fn fetch_tick(&mut self, ctx: &egui::Context) {
        let Some(polled) = self.fetch_win.job.as_ref().map(Job::poll) else {
            return;
        };
        match polled {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(done) => {
                self.fetch_win.job = None;
                match done {
                    Ok(fetched) => self.open_fetched(fetched),
                    Err(why) => self.fetch_win.message = why,
                }
            }
            Poll::Died => {
                self.fetch_win.job = None;
                self.fetch_win.message = "The request stopped without an answer.".into();
            }
        }
    }

    /// Put a record where it belongs: a nucleotide record in Sequences, protein sequences in the list of the Alignment page, a
    /// structure in the Structure page.
    fn open_fetched(&mut self, f: Fetched) {
        let how = if f.from_cache {
            "from the copy kept here"
        } else {
            "downloaded"
        };
        let id = f.request.id.clone();
        let said = match f.request.source.format() {
            Format::GenBank => f.documents().map(|docs| {
                let bases: usize = docs.iter().map(|d| d.len()).sum();
                docs.into_iter().for_each(|d| self.add_doc(d));
                self.workspace = Workspace::Sequences;
                format!("Opened {id} ({bases} bases), {how}.")
            }),
            Format::Fasta => f.documents().map(|docs| {
                let residues: usize = docs.iter().map(|d| d.len()).sum();
                for d in &docs {
                    self.platform.alignment.add(&d.name, &d.seq);
                }
                format!(
                    "Added {id} ({residues} residues) to the list of the Alignment page, {how}."
                )
            }),
            Format::Pdb => {
                let name = format!("{} {id}", f.request.source.label());
                self.platform.structure.show_pdb(&name, &f.body).map(|()| {
                    self.workspace = Workspace::Structure;
                    format!("Opened the structure {id}, {how}.")
                })
            }
        };
        let mut message = said.unwrap_or_else(|why| why);
        if let Some(why) = f.warning {
            message.push_str(&format!(" It was not kept for next time: {why}"));
        }
        self.fetch_win.message = message;
        self.fetch_win.last = Some((f.provenance, f.from_cache));
    }
}

#[cfg(test)]
mod tests;
