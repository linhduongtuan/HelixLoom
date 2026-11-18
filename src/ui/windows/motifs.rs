//! Tools → Motif Finder: look for motifs shared by DNA sequences (the way of MEME, `seq::motif`), see them as logos, scan the sequences
//! for more places that match one, and copy a motif in MEME's text format for other programs.

#[cfg(test)]
mod tests;

use crate::provenance::{Asked, Input, Provenance};
use crate::remote::client::Remote;
use crate::remote::motif::{self as remote_motif, MotifJob, MotifResult};
use crate::remote::tasks::{Request, Response};
use crate::remote::weight::Weight;
use crate::seq::motif::{
    background, compare, discover, scan, Found, Motif, ScoreDistribution, Settings, BASES,
};
use crate::seq::motif::{check_size, EFFORT_RANGE};
use crate::ui::app::App;
use crate::ui::design::{self, Look};
use crate::ui::job::{Job, Poll};
use crate::ui::platform::workstation::{compute_picker, Place};
use crate::ui::provenance as made_ui;
use eframe::egui::{self, Align2, Color32, FontId, RichText, ScrollArea, Sense, Stroke, Vec2};
use std::time::Duration;
// the browser has no `std::time::Instant` (it panics there); `web_time` is the same on a computer and works in a page
use web_time::Instant;

const MOST_SEQUENCES: usize = 2000;
/// The most bases searched at once (the search is cut down to a budget of work, but the sequences are all held, both strands, as windows;
/// a page has less memory and no second thread, so it is given less).
#[cfg(not(target_arch = "wasm32"))]
const MOST_BASES: usize = 300_000;
#[cfg(target_arch = "wasm32")]
const MOST_BASES: usize = 60_000;
/// The most rows of sites or hits written out: the counts say how many there are.
const MOST_ROWS: usize = 200;
const WORKING_REDRAW: Duration = Duration::from_millis(200);
const LOGO_COLUMN: f32 = 22.0;
const LOGO_HEIGHT: f32 = 90.0;
/// In a report: the room of the line above a logo and the gap under it, and the least width (for that line).
const LOGO_ROOM: f32 = 40.0;
const LOGOS_MIN_WIDTH: f32 = 360.0;
/// The most bits a column can say, which is the height of the logo.
const MOST_BITS: f32 = 2.0;
/// A block is given its letter when it is at least this tall, and the letter is at most this big.
const LETTER_FROM: f32 = 12.0;
const LETTER_SIZE: f32 = 15.0;

/// What a search was made on and found.
pub struct Run {
    pub names: Vec<String>,
    pub seqs: Vec<Vec<u8>>,
    pub found: Vec<Found>,
    pub settings: Settings,
    pub seconds: f64,
}

/// One place a motif matches, for the list.
pub struct Row {
    pub sequence: String,
    pub position: usize,
    pub reverse: bool,
    pub score: f64,
    pub p_value: f64,
}

pub struct MotifWin {
    pub open: bool,
    /// FASTA, or one sequence on each line.
    pub text: String,
    pub settings: Settings,
    /// A scan keeps the windows whose p-value is at most this.
    pub p_value: f64,
    job: Option<Job<Result<Run, String>>>,
    /// Where the search that is going (or the last one) ran, for the message.
    ran_on: &'static str,
    pub message: String,
    pub run: Option<Run>,
    /// The motif shown, of those found.
    pub chosen: usize,
    pub hits: Vec<Row>,
    pub hit_total: usize,
    /// The record of what made the motifs that are shown; none for a result that came from the workstation's list, whose
    /// settings are not known here.
    pub made: Option<Provenance>,
}

/// The inputs and the settings of a search, as its record says them: every sequence by its name, and what the search was told.
pub fn described(
    names: &[String],
    seqs: &[Vec<u8>],
    s: &Settings,
) -> (Vec<Input>, Vec<(String, String)>) {
    let inputs = names
        .iter()
        .zip(seqs)
        .map(|(name, bases)| Input::of(name, bases))
        .collect();
    let settings = [
        ("width", s.width.to_string()),
        ("motifs", s.motifs.to_string()),
        (
            "both strands",
            if s.both_strands { "yes" } else { "no" }.to_string(),
        ),
        ("starting words", s.seeds.to_string()),
        ("rounds", s.iterations.to_string()),
        ("effort", s.effort.to_string()),
    ];
    (
        inputs,
        settings
            .into_iter()
            .map(|(k, v)| (k.to_string(), v))
            .collect(),
    )
}

impl Default for MotifWin {
    fn default() -> MotifWin {
        MotifWin {
            open: false,
            text: String::new(),
            settings: Settings::default(),
            p_value: 1e-4,
            job: None,
            ran_on: "this computer",
            message: String::new(),
            run: None,
            chosen: 0,
            hits: Vec::new(),
            made: None,
            hit_total: 0,
        }
    }
}

/// The names and bases of the sequences in `text`: FASTA, or when there is no `>` one sequence for each line.
pub fn read_sequences(text: &str) -> Result<(Vec<String>, Vec<Vec<u8>>), String> {
    let (names, seqs) = crate::seq::motif::read_sequences(text)?;
    if seqs.len() > MOST_SEQUENCES {
        return Err(format!(
            "{} sequences are more than the {MOST_SEQUENCES} that can be searched at once.",
            seqs.len()
        ));
    }
    let bases: usize = seqs.iter().map(Vec::len).sum();
    if bases > MOST_BASES {
        return Err(format!(
            "{bases} bases are more than the {MOST_BASES} that can be searched at once."
        ));
    }
    Ok((names, seqs))
}

/// For each column of `motif`, the four letters with the height each is drawn at (a share of the column's two bits, so that the
/// letters of a column that says little are short), the tallest last, so that it is drawn on top.
pub fn logo(motif: &Motif) -> Vec<[(char, f32); 4]> {
    motif
        .probs
        .iter()
        .map(|column| {
            let bits = 2.0
                + column
                    .iter()
                    .filter(|&&p| p > 0.0)
                    .map(|&p| p * p.log2())
                    .sum::<f64>();
            let mut letters: [(char, f32); 4] =
                std::array::from_fn(|b| (BASES[b] as char, (column[b] * bits) as f32 / MOST_BITS));
            letters.sort_by(|a, b| a.1.total_cmp(&b.1));
            letters
        })
        .collect()
}

/// The colour of the block of `letter` in the logo: the letters of the trace when the safe colours are chosen.
fn letter_colour(letter: char, look: &Look) -> Color32 {
    if look.safe_colours {
        return design::trace_color(letter as u8, look);
    }
    match letter {
        'A' => Color32::from_rgb(0, 160, 60),
        'C' => Color32::from_rgb(40, 90, 220),
        'G' => Color32::from_rgb(240, 150, 0),
        _ => Color32::from_rgb(210, 40, 40),
    }
}

/// The room the logos of `found` take, one under the other with a line above each: for a report.
pub(in crate::ui) fn logos_room(found: &[Found]) -> Vec2 {
    let widest = found.iter().map(|f| f.motif.probs.len()).max().unwrap_or(0);
    Vec2::new(
        (LOGO_COLUMN * widest as f32 + 4.0).max(LOGOS_MIN_WIDTH),
        (LOGO_HEIGHT + LOGO_ROOM) * found.len() as f32,
    )
}

/// The logos of `found`, one under the other, each under its number and its likeliest letters.
pub(in crate::ui) fn draw_logos(ui: &mut egui::Ui, found: &[Found], look: &Look) {
    for (k, one) in found.iter().enumerate() {
        ui.label(format!("Motif {}: {}", k + 1, one.motif.consensus()));
        draw_logo(ui, &one.motif, look);
    }
}

fn draw_logo(ui: &mut egui::Ui, motif: &Motif, look: &Look) {
    let columns = logo(motif);
    let (rect, _) = ui.allocate_exact_size(
        Vec2::new(LOGO_COLUMN * columns.len() as f32 + 4.0, LOGO_HEIGHT + 4.0),
        Sense::hover(),
    );
    let painter = ui.painter_at(rect);
    let floor = rect.bottom() - 2.0;
    // Each letter is a block of its colour as wide as the column and as tall as its share of the column's two bits (a glyph cannot be
    // stretched to fit, so the letter sits in the block), the likeliest on top.
    for (k, letters) in columns.iter().enumerate() {
        let left = rect.left() + 2.0 + LOGO_COLUMN * k as f32;
        let mut bottom = floor;
        for &(letter, share) in letters {
            let height = share * LOGO_HEIGHT;
            if height >= 1.0 {
                let block = egui::Rect::from_min_max(
                    egui::pos2(left + 1.0, bottom - height),
                    egui::pos2(left + LOGO_COLUMN - 1.0, bottom),
                );
                let fill = letter_colour(letter, look);
                painter.rect_filled(block, 2.0, fill);
                if height >= LETTER_FROM {
                    // white on the blocks of old; the safe ones include a light grey, which needs black
                    let ink = if look.safe_colours {
                        design::ink_on(look, fill)
                    } else {
                        Color32::WHITE
                    };
                    painter.text(
                        block.center(),
                        Align2::CENTER_CENTER,
                        letter,
                        FontId::monospace((height * 0.8).min(LETTER_SIZE)),
                        ink,
                    );
                }
            }
            bottom -= height;
        }
    }
    painter.line_segment(
        [
            egui::pos2(rect.left(), floor),
            egui::pos2(rect.right(), floor),
        ],
        Stroke::new(1.0, ui.visuals().widgets.noninteractive.bg_stroke.color),
    );
}

/// What the workstation found, as a run of this window: only if its motifs are about the sequences that were sent (an answer is not
/// trusted to be), so that no site names a sequence or a place that is not there.
fn from_workstation(
    result: MotifResult,
    names: Vec<String>,
    seqs: Vec<Vec<u8>>,
    settings: Settings,
    seconds: f64,
) -> Result<Run, String> {
    result.fits(&seqs)?;
    Ok(Run {
        names,
        seqs,
        found: result.found,
        settings,
        seconds,
    })
}

/// The search itself, which a window runs on a thread of its own.
fn search(names: Vec<String>, seqs: Vec<Vec<u8>>, settings: Settings) -> Result<Run, String> {
    let start = Instant::now();
    discover(&seqs, &settings).map(|found| Run {
        names,
        seqs,
        found,
        settings,
        seconds: start.elapsed().as_secs_f64(),
    })
}

impl MotifWin {
    /// The example, searched at once (for a picture of the window): the sequences, a search over them, and what it found.
    pub fn show_example(&mut self) {
        self.open = true;
        self.text = crate::sample::motif_demo();
        self.settings.width = 10;
        if let Ok((names, seqs)) = read_sequences(&self.text) {
            let outcome = search(names, seqs, self.settings);
            self.take(outcome, "this computer", true);
        }
    }

    /// Start the search that the text and the settings ask for, or say why not.
    pub fn start(&mut self, ctx: &egui::Context) {
        if self.job.is_some() {
            return;
        }
        let (names, seqs) = match read_sequences(&self.text) {
            Ok(read) => read,
            Err(why) => {
                self.message = why;
                return;
            }
        };
        let settings = self.settings;
        (self.ran_on, self.message) = ("this computer", "Searching…".into());
        self.job = Some(Job::spawn(ctx, move || search(names, seqs, settings)));
    }

    /// Start the search as a job of the workstation, or say why not. The sequences are sent as they are; the names stay here.
    pub fn start_there(
        &mut self,
        ctx: &egui::Context,
        remote: Remote,
        names: Vec<String>,
        seqs: Vec<Vec<u8>>,
    ) {
        if self.job.is_some() {
            return;
        }
        let settings = self.settings;
        if seqs.len() > remote_motif::MAX_SEQUENCES {
            self.message = format!(
                "{} sequences are more than the {} that a job takes.",
                seqs.len(),
                remote_motif::MAX_SEQUENCES
            );
            return;
        }
        if let Err(why) = check_size(
            seqs.iter().map(Vec::len).sum(),
            settings.width,
            settings.both_strands,
        ) {
            self.message = why;
            return;
        }
        (self.ran_on, self.message) = ("the workstation", "Searching on the workstation…".into());
        let label = format!("Motif search in {} sequences", seqs.len());
        let request = Request::Motif(MotifJob {
            seqs: seqs.clone(),
            settings,
        });
        self.job = Some(Job::spawn_async(ctx, async move {
            let started = web_time::Instant::now();
            let result = remote
                .run_async(&request, &label)
                .await
                .and_then(Response::into_motif)?;
            from_workstation(
                result,
                names,
                seqs,
                settings,
                started.elapsed().as_secs_f64(),
            )
        }));
    }

    pub fn running(&self) -> bool {
        self.job.is_some()
    }

    /// Take the answer of the search when it has come (once a frame, as every window with a job does).
    pub fn poll(&mut self, ctx: &egui::Context) {
        let Some(job) = &self.job else {
            return;
        };
        match job.poll() {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(outcome) => {
                self.job = None;
                // a search that this window started: its sequences and its settings are the ones of the run
                self.take(outcome, self.ran_on, true);
            }
            Poll::Died => {
                self.job = None;
                self.message = "The search stopped without an answer.".into();
            }
        }
    }

    /// The inputs and the settings as they are in the window now, to say what has changed since the motifs were found.
    pub fn as_it_is(&self) -> (Vec<Input>, Vec<(String, String)>) {
        match read_sequences(&self.text) {
            Ok((names, seqs)) => described(&names, &seqs, &self.settings),
            Err(_) => described(&[], &[], &self.settings),
        }
    }

    /// Take in what a search found, or why it could not. `ran_on` is where it ran; `known` says that the search was made with
    /// the sequences and the settings the run holds, so that a record of it can be kept.
    fn take(&mut self, outcome: Result<Run, String>, ran_on: &'static str, known: bool) {
        self.hits.clear();
        self.hit_total = 0;
        self.made = None;
        match outcome {
            Ok(run) if run.found.is_empty() => {
                self.message =
                    "No motif was found: the sequences hold no stretch of A, C, G and T that long."
                        .into();
                self.run = None;
            }
            Ok(run) => {
                self.message = format!(
                    "{} motif{} in {} sequences, in {:.1} s{}.",
                    run.found.len(),
                    if run.found.len() == 1 { "" } else { "s" },
                    run.seqs.len(),
                    run.seconds,
                    if ran_on == "the workstation" {
                        " on the workstation"
                    } else {
                        ""
                    }
                );
                self.chosen = 0;
                if known {
                    let (inputs, parameters) = described(&run.names, &run.seqs, &run.settings);
                    let asked = Asked {
                        analysis: "Motif search",
                        inputs,
                        parameters,
                        ran_on,
                    };
                    self.made = Some(asked.made(Vec::new(), run.seconds));
                }
                self.run = Some(run);
            }
            Err(why) => {
                self.message = why;
                self.run = None;
            }
        }
    }

    /// The windows of the sequences that match the chosen motif with a p-value of at most `p_value` (the most rows are kept).
    pub fn scan_chosen(&mut self) {
        self.hits.clear();
        self.hit_total = 0;
        let Some(run) = &self.run else {
            return;
        };
        let Some(found) = run.found.get(self.chosen) else {
            return;
        };
        let bg = background(&run.seqs);
        let distribution = ScoreDistribution::new(&found.motif, &bg);
        let least = distribution.score_for(self.p_value.clamp(1e-12, 1.0));
        for (name, seq) in run.names.iter().zip(&run.seqs) {
            for hit in scan(seq, &found.motif, &bg, least, run.settings.both_strands) {
                self.hit_total += 1;
                if self.hits.len() < MOST_ROWS {
                    self.hits.push(Row {
                        sequence: name.clone(),
                        position: hit.position,
                        reverse: hit.reverse,
                        score: hit.score,
                        p_value: distribution.p_value(hit.score),
                    });
                }
            }
        }
        self.message = format!(
            "{} window{} match {} with a p-value of at most {:e}.",
            self.hit_total,
            if self.hit_total == 1 { "" } else { "s" },
            found.motif.name,
            self.p_value
        );
    }

    /// The chosen motif as MEME's text.
    pub fn chosen_meme(&self) -> Option<String> {
        let run = self.run.as_ref()?;
        let found = run.found.get(self.chosen)?;
        Some(
            found
                .motif
                .to_meme(&background(&run.seqs), found.sites.len()),
        )
    }
}

/// The bases in the sequences of `text` (a FASTA header is not counted), roughly: for the choice of where the search runs.
fn bases_in(text: &str) -> usize {
    text.lines()
        .filter(|l| !l.starts_with('>'))
        .map(|l| l.bytes().filter(u8::is_ascii_alphabetic).count())
        .sum()
}

impl App {
    /// Show what the workstation found, and open the window: the sequences are not kept with a result, so it is shown only if the
    /// window still holds the ones it is about.
    pub fn show_motif_result(&mut self, result: MotifResult) {
        let w = &mut self.motif_win;
        w.open = true;
        let held = crate::seq::motif::read_sequences(&w.text)
            .and_then(|(names, seqs)| result.fits(&seqs).map(|()| (names, seqs)));
        match held {
            Ok((names, seqs)) if !result.found.is_empty() => {
                let run = Run {
                    names,
                    seqs,
                    found: result.found,
                    settings: w.settings,
                    seconds: result.seconds,
                };
                // fetched from the workstation's list: the settings it was made with are not known here. Where a search of
                // this window that is going runs is left as it is
                w.take(Ok(run), "the workstation", false);
            }
            _ => {
                w.message = format!(
                    "The workstation found {} motif{}, for sequences this window no longer holds: put those sequences in the box and open the result again.",
                    result.found.len(),
                    if result.found.len() == 1 { "" } else { "s" }
                );
            }
        }
    }

    /// Start the search where the compute choice says: on this computer, or as a job of the workstation.
    pub(in crate::ui) fn start_motifs(&mut self, ctx: &egui::Context) {
        let (names, seqs) = match crate::seq::motif::read_sequences(&self.motif_win.text) {
            Ok(read) => read,
            Err(why) => {
                self.motif_win.message = why;
                return;
            }
        };
        let bases = seqs.iter().map(Vec::len).sum();
        match self.platform.place(&Weight::of_motifs(bases)) {
            Place::There(remote) => self.motif_win.start_there(ctx, remote, names, seqs),
            Place::Here => self.motif_win.start(ctx),
        }
    }

    pub fn motif_window(&mut self, ctx: &egui::Context) {
        self.motif_win.poll(ctx);
        if !self.motif_win.open {
            return;
        }
        let look = self.prefs.look();
        let view = self
            .platform
            .workstation
            .view(&Weight::of_motifs(bases_in(&self.motif_win.text)));
        let mut start = false;
        let mut open = true;
        let mut asks = made_ui::Asks::default();
        crate::ui::design::window("Motif Finder")
            .open(&mut open)
            .default_width(640.0)
            .show(ctx, |ui| {
                let w = &mut self.motif_win;
                ui.label(
                    "Look for a motif that the DNA sequences share, as MEME does: at most one site in a sequence, on either \
                     strand. Then scan them for more places that match it.",
                );
                ui.horizontal(|ui| {
                    if ui.button("Use an example").clicked() {
                        w.text = crate::sample::motif_demo();
                        w.settings.width = 10;
                    }
                    if ui.button("Clear").clicked() {
                        *w = MotifWin::default();
                        w.open = true;
                    }
                });
                ScrollArea::vertical().max_height(100.0).show(ui, |ui| {
                    ui.add(
                        egui::TextEdit::multiline(&mut w.text)
                            .font(egui::TextStyle::Monospace)
                            .hint_text("Paste DNA sequences: FASTA, or one on each line")
                            .desired_rows(4)
                            .desired_width(f32::INFINITY),
                    );
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Width");
                    ui.add(egui::DragValue::new(&mut w.settings.width).range(3..=30));
                    ui.label("Motifs");
                    ui.add(egui::DragValue::new(&mut w.settings.motifs).range(1..=5));
                    ui.checkbox(&mut w.settings.both_strands, "Both strands");
                    ui.label("Effort");
                    ui.add(
                        egui::DragValue::new(&mut w.settings.effort)
                            .range(EFFORT_RANGE.0..=EFFORT_RANGE.1)
                            .speed(0.05)
                            .max_decimals(2),
                    )
                    .on_hover_text(
                        "How much work the search may spend over a big set, 1 being the usual budget of about five seconds. A small \
                         set is searched in full whatever this is.",
                    );
                });
                compute_picker(ui, &look, view);
                ui.horizontal(|ui| {
                    let idle = !w.running();
                    start = ui
                        .add_enabled(idle, egui::Button::new("Find motifs"))
                        .clicked();
                    if !idle {
                        ui.spinner();
                    }
                });
                if !w.message.is_empty() {
                    ui.label(RichText::new(w.message.as_str()).weak());
                }
                if let Some(made) = &w.made {
                    asks = made_ui::section(ui, &look, made, || w.as_it_is(), !w.running());
                }
                let Some(run) = &w.run else {
                    return;
                };
                ui.separator();
                for (k, found) in run.found.iter().enumerate() {
                    ui.horizontal(|ui| {
                        if ui.radio_value(&mut w.chosen, k, "").clicked() {
                            w.hits.clear();
                            w.hit_total = 0;
                        }
                        ui.label(RichText::new(format!(
                            "{}: {}  ({:.1} bits, {} sites, gain {:.0} bits)",
                            found.motif.name,
                            found.motif.consensus(),
                            found.motif.information(),
                            found.sites.len(),
                            found.gain
                        )).monospace());
                    });
                }
                let Some(found) = run.found.get(w.chosen) else {
                    return;
                };
                if let [first, second, ..] = run.found.as_slice() {
                    // the second motif that is only the first again is a sign that there was one motif to find
                    ui.label(
                        RichText::new(format!(
                            "{} and {} are alike by {:.2} (1: the same motif, 0: unrelated).",
                            first.motif.name,
                            second.motif.name,
                            compare(&first.motif, &second.motif, 5)
                        ))
                        .weak(),
                    );
                }
                draw_logo(ui, &found.motif, &look);
                ScrollArea::vertical().id_salt("motif sites").max_height(90.0).show(ui, |ui| {
                    for site in found.sites.iter().take(MOST_ROWS) {
                        ui.label(
                            RichText::new(format!(
                                "{:<16} {:>6}  {}  score {:>6.2}  p {:.2}",
                                run.names[site.sequence],
                                site.position + 1,
                                if site.reverse { "−" } else { "+" },
                                site.score,
                                site.probability
                            ))
                            .monospace(),
                        );
                    }
                });
                let mut copy = false;
                let mut scan_now = false;
                ui.horizontal_wrapped(|ui| {
                    copy = ui.button("Copy as MEME text").clicked();
                    ui.label("Scan for p-value at most");
                    ui.add(
                        egui::DragValue::new(&mut w.p_value)
                            .speed(1e-5)
                            .range(1e-9..=0.1)
                            .custom_formatter(|v, _| format!("{v:e}")),
                    );
                    scan_now = ui.button("Scan").clicked();
                });
                if copy {
                    if let Some(text) = w.chosen_meme() {
                        ctx.copy_text(text);
                    }
                }
                if scan_now {
                    w.scan_chosen();
                }
                if !w.hits.is_empty() {
                    ScrollArea::vertical().id_salt("motif hits").max_height(140.0).show(ui, |ui| {
                        for row in &w.hits {
                            ui.label(
                                RichText::new(format!(
                                    "{:<16} {:>6}  {}  score {:>6.2}  p {:.1e}",
                                    row.sequence,
                                    row.position + 1,
                                    if row.reverse { "−" } else { "+" },
                                    row.score,
                                    row.p_value
                                ))
                                .monospace(),
                            );
                        }
                    });
                }
            });
        self.motif_win.open = open;
        if let (true, Some(made)) = (asks.copy, &self.motif_win.made) {
            ctx.copy_text(made.to_text());
            self.motif_win.message = "Copied the record of how the motifs were found.".into();
        }
        self.checks_win.open |= asks.checks;
        if asks.again {
            // with the settings the motifs that are shown were found with
            if let Some(run) = &self.motif_win.run {
                self.motif_win.settings = run.settings;
            }
        }
        if start || asks.again {
            self.start_motifs(ctx);
        }
    }
}
