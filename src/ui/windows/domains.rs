//! Tools → Protein Domains: which domains a protein has, from HMMER's search of Pfam (`pipeline::domains`, `seq::hmmer`), drawn along
//! the sequence. HMMER and Pfam are other projects' programs and data: the search runs where the compute choice says (this computer
//! or the workstation), wherever they are installed.

#[cfg(test)]
mod tests;

use crate::pipeline::domains::Cutoff;
use crate::provenance::{Asked, Input, Provenance};
use crate::remote::domains::{self as remote_domains, DomainsJob, DomainsResult};
use crate::remote::tasks::{Request, Response};
use crate::remote::weight::Weight;
use crate::seq::hmmer::{along, read_proteins, to_table, DomainHit};
use crate::ui::app::App;
use crate::ui::design::{self, Look};
use crate::ui::fonts::FontWeight;
use crate::ui::job::{Job, Poll};
use crate::ui::platform::workstation::{compute_picker, Place};
use crate::ui::provenance as made_ui;
use eframe::egui::{self, Align2, Color32, FontId, RichText, ScrollArea, Sense, Stroke, Vec2};
use std::sync::Arc;
use std::time::Duration;

const WORKING_REDRAW: Duration = Duration::from_millis(200);
/// The most sequences drawn and the most rows of the table, whatever the search found: the copied table has them all.
const MOST_DRAWN: usize = 50;
const MOST_ROWS: usize = 500;
const LANE_HEIGHT: f32 = 20.0;
const BAR_HEIGHT: f32 = 6.0;
const LABEL_WIDTH: f32 = 150.0;
/// In a report: the width the sequences are drawn at, and the room between two of them.
const TRACKS_WIDTH: f32 = 860.0;
const TRACK_GAP: f32 = 8.0;
/// The E-value that is offered when the person leaves Pfam's own thresholds.
const DEFAULT_EVALUE: f64 = 1e-3;
/// The window opens tall enough for a result and no taller than the screen.
const WINDOW_ROOM: f32 = 120.0;
const MIN_WINDOW_HEIGHT: f32 = 320.0;
const MAX_WINDOW_HEIGHT: f32 = 900.0;
/// Ubiquitin (the sequence of the structure 1UBQ), a protein with one small, well-known domain, for a first try.
/// The chains of the two structures that come with the program (`assets/structure`): ubiquitin (PDB 1UBQ), and the two variable
/// domains of the antibody D1.3 with the lysozyme it binds (PDB 1VFB). With Pfam 39.0 they are in the ubiquitin family, the
/// immunoglobulin V-set and the C-type lysozymes (`remote/domains/real.rs` checks it where HMMER and Pfam are installed).
pub(crate) const EXAMPLE: &str = "\
>ubiquitin
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
>D1.3_light_chain variable domain of the antibody D1.3
DIVLTQSPASLSASVGETVTITCRASGNIHNYLAWYQQKQGKSPQLLVYYTTTLADGVPSRFSGSGSGTQYSLKINSLQPEDFGSYYCQHFWSTPRTFGGGTKLEIK
>D1.3_heavy_chain variable domain of the antibody D1.3
QVQLQESGPGLVAPSQSLSITCTVSGFSLTGYGVNWVRQPPGKGLEWLGMIWGDGNTDYNSALKSRLSISKDNSKSQVFLKMNSLHTDDTARYYCARERDYRLDYWGQGTTLTVSS
>lysozyme of hen egg white
KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL
";

/// A colour-blind-safe set (Okabe and Ito), for telling the families apart.
const PALETTE: [(u8, u8, u8); 8] = [
    (230, 159, 0),
    (86, 180, 233),
    (0, 158, 115),
    (240, 228, 66),
    (0, 114, 178),
    (213, 94, 0),
    (204, 121, 167),
    (153, 153, 153),
];

/// Which colour a family is drawn in: the same family always has the same one, in every sequence and every run.
pub fn colour_index(family: &str) -> usize {
    let hash = family.bytes().fold(0xcbf2_9ce4_8422_2325_u64, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    });
    (hash % PALETTE.len() as u64) as usize
}

/// The lane (a row under the bar) of each domain of `hits` as `along` orders them: the first lane in which the one before ends
/// before this one starts, so that domains that overlap are drawn one under the other.
pub fn lanes(hits: &[&DomainHit]) -> Vec<usize> {
    let mut ends: Vec<usize> = Vec::new();
    hits.iter()
        .map(|h| match ends.iter().position(|&end| end < h.ali_from) {
            Some(lane) => {
                ends[lane] = h.ali_to;
                lane
            }
            None => {
                ends.push(h.ali_to);
                ends.len() - 1
            }
        })
        .collect()
}

/// The horizontal place of residue `position` (counted from 1, its left edge) in a bar `width` wide for a sequence of `length`.
pub fn x_of(position: usize, length: usize, width: f32) -> f32 {
    width * (position.saturating_sub(1)) as f32 / length.max(1) as f32
}

pub struct DomainsWin {
    pub open: bool,
    /// FASTA, or one protein on each line.
    pub text: String,
    /// Pfam's own thresholds for each family, or one E-value for all.
    pub curated: bool,
    pub evalue: f64,
    job: Option<Job<Result<DomainsResult, String>>>,
    /// Where the search that is going (or the last one) ran, for the message.
    ran_on: &'static str,
    /// Runs the search on this computer; the tests replace it with a program that stands in for HMMER.
    runner: fn(&DomainsJob) -> Result<DomainsResult, String>,
    pub message: String,
    pub result: Option<Arc<DomainsResult>>,
    /// What the search that is going was asked, with the threshold to bring back for *Run again*.
    asked: Option<(Asked, Kept)>,
    /// The record of what made the domains that are shown; none for a result that came from the workstation's list, whose
    /// sequences and threshold are not known here.
    pub made: Option<(Provenance, Kept)>,
}

/// The threshold of a search as the window holds it: Pfam's own, or else the E-value.
pub type Kept = (bool, f64);

/// What a search is asked, for its record.
fn asking(inputs: Vec<Input>, parameters: Vec<(String, String)>, ran_on: &'static str) -> Asked {
    Asked {
        analysis: "Protein domain search",
        inputs,
        parameters,
        ran_on,
    }
}

/// The inputs and the settings of a search, as its record says them: every protein by its name, and the threshold.
pub fn described(job: &DomainsJob) -> (Vec<Input>, Vec<(String, String)>) {
    let inputs = job
        .queries
        .iter()
        .map(|(name, residues)| Input::of(name, residues))
        .collect();
    let threshold = match job.cutoff {
        Cutoff::Curated => "Pfam's own for each family".to_string(),
        Cutoff::Evalue(most) => format!("E-value at most {most:e}"),
    };
    (inputs, vec![("threshold".to_string(), threshold)])
}

impl Default for DomainsWin {
    fn default() -> Self {
        DomainsWin {
            open: false,
            text: String::new(),
            curated: true,
            evalue: DEFAULT_EVALUE,
            job: None,
            ran_on: "this computer",
            runner: remote_domains::run,
            message: String::new(),
            result: None,
            asked: None,
            made: None,
        }
    }
}

/// What a search found, in a sentence.
pub fn summary(r: &DomainsResult, ran_on: &str) -> String {
    let with = r
        .hits
        .iter()
        .map(|h| h.sequence)
        .collect::<std::collections::HashSet<_>>()
        .len();
    let what = if r.program.is_empty() {
        "HMMER"
    } else {
        r.program.as_str()
    };
    let against = if r.database.is_empty() {
        String::new()
    } else {
        format!(" against {}", r.database)
    };
    if r.hits.is_empty() {
        return format!(
            "No domain found in {} sequence{} by {what}{against}, in {:.1} s on {ran_on}.",
            r.names.len(),
            if r.names.len() == 1 { "" } else { "s" },
            r.seconds
        );
    }
    format!(
        "{} domain{} in {with} of {} sequence{} by {what}{against}, in {:.1} s on {ran_on}.",
        r.hits.len(),
        if r.hits.len() == 1 { "" } else { "s" },
        r.names.len(),
        if r.names.len() == 1 { "" } else { "s" },
        r.seconds
    )
}

impl DomainsWin {
    pub fn running(&self) -> bool {
        self.job.is_some()
    }

    fn cutoff(&self) -> Cutoff {
        if self.curated {
            Cutoff::Curated
        } else {
            Cutoff::Evalue(self.evalue)
        }
    }

    /// The residues of the sequences in the text, for the choice of where the search runs.
    pub(in crate::ui) fn weight(&self) -> Weight {
        let residues = self
            .text
            .lines()
            .filter(|l| !l.starts_with('>'))
            .map(|l| l.bytes().filter(u8::is_ascii_alphabetic).count())
            .sum();
        let count = self
            .text
            .lines()
            .filter(|l| l.starts_with('>'))
            .count()
            .max(1);
        Weight::of_domains(count, residues)
    }

    fn job_asked_for(&self) -> Result<DomainsJob, String> {
        let queries = read_proteins(&self.text)?;
        let job = DomainsJob {
            queries,
            cutoff: self.cutoff(),
        };
        job.cutoff.check()?;
        Ok(job)
    }

    /// Start the search where it is to run, or say why not.
    pub fn start(&mut self, ctx: &egui::Context, place: Place) {
        if self.job.is_some() {
            return;
        }
        let job = match self.job_asked_for() {
            Ok(job) => job,
            Err(why) => {
                self.message = why;
                return;
            }
        };
        (self.result, self.made) = (None, None);
        let (inputs, parameters) = described(&job);
        match place {
            Place::Here => {
                (self.ran_on, self.message) = ("this computer", "Searching…".into());
                let runner = self.runner;
                self.job = Some(Job::spawn(ctx, move || runner(&job)));
            }
            Place::There(remote) => {
                (self.ran_on, self.message) =
                    ("the workstation", "Searching on the workstation…".into());
                let label = format!("Domain search of {} sequences", job.queries.len());
                let request = Request::Domains(job);
                self.job = Some(Job::spawn_async(ctx, async move {
                    remote
                        .run_async(&request, &label)
                        .await
                        .and_then(Response::into_domains)
                }));
            }
        }
        let asked = asking(inputs, parameters, self.ran_on);
        self.asked = Some((asked, (self.curated, self.evalue)));
    }

    /// The inputs and the threshold as they are in the window now, to say what has changed since the domains were found.
    pub fn as_it_is(&self) -> (Vec<Input>, Vec<(String, String)>) {
        described(&DomainsJob {
            queries: read_proteins(&self.text).unwrap_or_default(),
            cutoff: self.cutoff(),
        })
    }

    pub fn poll(&mut self, ctx: &egui::Context) {
        let Some(job) = &self.job else {
            return;
        };
        match job.poll() {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(outcome) => {
                self.job = None;
                let asked = self.asked.take();
                self.take(outcome, self.ran_on, asked);
            }
            Poll::Died => {
                self.job = None;
                self.message = "The search stopped without an answer.".into();
            }
        }
    }

    /// The window open with its example searched, for a picture of it: this waits for HMMER (a second or two).
    pub fn show_example(&mut self) {
        self.open = true;
        self.text = EXAMPLE.to_string();
        match self.job_asked_for() {
            Ok(job) => {
                let (inputs, parameters) = described(&job);
                // the example is searched on this computer, wherever the last search ran
                let asked = asking(inputs, parameters, "this computer");
                let outcome = (self.runner)(&job);
                self.take(
                    outcome,
                    "this computer",
                    Some((asked, (self.curated, self.evalue))),
                );
            }
            Err(why) => self.message = why,
        }
    }

    /// Take in what a search found, or why it could not. `ran_on` is where it ran. `asked` is what the search was asked, when
    /// that is known here: the record of the domains is made from it, with the program and the database the answer names
    /// (an answer of an older workstation names none).
    fn take(
        &mut self,
        outcome: Result<DomainsResult, String>,
        ran_on: &'static str,
        asked: Option<(Asked, Kept)>,
    ) {
        match outcome {
            Ok(result) => {
                self.message = summary(&result, ran_on);
                self.made = asked.map(|(asked, kept)| {
                    let programs = [&result.program, &result.database]
                        .into_iter()
                        .filter(|named| !named.is_empty())
                        .cloned()
                        .collect();
                    (asked.made(programs, result.seconds), kept)
                });
                self.result = Some(Arc::new(result));
            }
            Err(why) => {
                self.message = why;
                (self.result, self.made) = (None, None);
            }
        }
    }
}

fn family_colour(family: &str) -> Color32 {
    let (r, g, b) = PALETTE[colour_index(family)];
    Color32::from_rgb(r, g, b)
}

/// The sequences that are drawn of a result: all of them, up to `MOST_DRAWN`.
pub(in crate::ui) fn drawn(r: &DomainsResult) -> usize {
    r.names.len().min(MOST_DRAWN)
}

/// The room the drawing of a result takes at a width of `TRACKS_WIDTH`: for a report.
pub(in crate::ui) fn tracks_room(r: &DomainsResult) -> Vec2 {
    let height: f32 = (0..drawn(r))
        .map(|k| {
            let lanes_used = lanes(&along(&r.hits, k)).iter().max().map_or(1, |m| m + 1);
            BAR_HEIGHT + 4.0 + LANE_HEIGHT * lanes_used as f32 + TRACK_GAP
        })
        .sum();
    Vec2::new(TRACKS_WIDTH, height + TRACK_GAP)
}

/// Every sequence of a result (the first `MOST_DRAWN`) with its domains along it.
pub(in crate::ui) fn draw_tracks(ui: &mut egui::Ui, look: &Look, r: &DomainsResult) {
    for (k, (name, length)) in r.names.iter().zip(&r.lengths).enumerate().take(drawn(r)) {
        draw_sequence(ui, look, name, *length, &along(&r.hits, k));
    }
}

/// One sequence: its name and length, and its domains as boxes under a bar the length of the sequence.
fn draw_sequence(ui: &mut egui::Ui, look: &Look, name: &str, length: usize, hits: &[&DomainHit]) {
    let lane_of = lanes(hits);
    let lanes_used = lane_of.iter().max().map_or(1, |m| m + 1);
    let width = (ui.available_width() - LABEL_WIDTH - 12.0).max(120.0);
    let height = BAR_HEIGHT + 4.0 + LANE_HEIGHT * lanes_used as f32;
    ui.horizontal(|ui| {
        ui.allocate_ui(Vec2::new(LABEL_WIDTH, height), |ui| {
            // the names take the same room whatever their length, so that the sequences begin in one column; a name too long
            // for it is cut, and read whole under the pointer
            ui.set_min_width(LABEL_WIDTH);
            ui.set_max_width(LABEL_WIDTH);
            ui.vertical(|ui| {
                ui.add(egui::Label::new(RichText::new(name).monospace()).truncate())
                    .on_hover_text(name);
                ui.label(RichText::new(format!("{length} residues")).small().weak());
            });
        });
        let (rect, _) = ui.allocate_exact_size(Vec2::new(width, height), Sense::hover());
        let painter = ui.painter_at(rect);
        let bar = egui::Rect::from_min_size(rect.min, Vec2::new(width, BAR_HEIGHT));
        painter.rect_filled(bar, 2.0, ui.visuals().widgets.noninteractive.bg_stroke.color);
        for (h, lane) in hits.iter().zip(&lane_of) {
            let left = rect.left() + x_of(h.ali_from, length, width);
            let right = (rect.left() + x_of(h.ali_to + 1, length, width)).max(left + 3.0);
            let top = rect.top() + BAR_HEIGHT + 4.0 + LANE_HEIGHT * *lane as f32;
            let block = egui::Rect::from_min_max(egui::pos2(left, top), egui::pos2(right, top + LANE_HEIGHT - 3.0));
            let fill = family_colour(&h.family);
            painter.rect_filled(block, 3.0, fill);
            painter.rect_stroke(block, 3.0, Stroke::new(1.0, Color32::from_black_alpha(90)), egui::StrokeKind::Inside);
            let galley_width = h.family.chars().count() as f32 * 7.0 + 6.0;
            if block.width() >= galley_width {
                painter.text(
                    block.center(),
                    Align2::CENTER_CENTER,
                    &h.family,
                    FontId::proportional(12.0),
                    design::ink_on(look, fill),
                );
            }
            let id = ui.id().with(("domain", h.ali_from, h.ali_to, &h.family));
            ui.interact(block, id, Sense::hover()).on_hover_text(format!(
                "{} {}\n{}\nresidues {}-{} (envelope {}-{})\nscore {:.1}, E-value {:e}\nprofile {}-{} of {}",
                h.family,
                h.accession,
                h.description,
                h.ali_from,
                h.ali_to,
                h.env_from,
                h.env_to,
                h.score,
                h.i_evalue,
                h.hmm_from,
                h.hmm_to,
                h.family_length
            ));
        }
    });
}

fn table(ui: &mut egui::Ui, r: &DomainsResult) {
    egui::Grid::new("domains_table")
        .striped(true)
        .show(ui, |ui| {
            for h in [
                "Sequence",
                "Domain",
                "Accession",
                "From–To",
                "Score",
                "E-value",
                "Profile",
            ] {
                ui.label(RichText::new(h).semibold());
            }
            ui.end_row();
            for hit in r.hits.iter().take(MOST_ROWS) {
                ui.label(&r.names[hit.sequence]);
                ui.label(RichText::new(&hit.family).semibold())
                    .on_hover_text(&hit.description);
                ui.label(RichText::new(&hit.accession).monospace());
                ui.label(RichText::new(format!("{}–{}", hit.ali_from, hit.ali_to)).monospace());
                ui.label(RichText::new(format!("{:.1}", hit.score)).monospace());
                ui.label(RichText::new(format!("{:.1e}", hit.i_evalue)).monospace());
                ui.label(
                    RichText::new(format!(
                        "{}–{} of {}",
                        hit.hmm_from, hit.hmm_to, hit.family_length
                    ))
                    .monospace(),
                );
                ui.end_row();
            }
        });
}

impl App {
    /// Show what the workstation found, and open the window.
    pub fn show_domains_result(&mut self, result: DomainsResult) {
        let w = &mut self.domains_win;
        w.open = true;
        // fetched from the workstation's list: the sequences and the threshold it was made with are not known here. Where
        // a search of this window that is going runs is left as it is
        w.take(Ok(result), "the workstation", None);
    }

    pub fn domains_window(&mut self, ctx: &egui::Context) {
        self.domains_win.poll(ctx);
        if !self.domains_win.open {
            return;
        }
        let look = self.prefs.look();
        let view = self.platform.workstation.view(&self.domains_win.weight());
        let place = self.platform.place(&self.domains_win.weight());
        let mut start = false;
        let mut open = true;
        let mut asks = made_ui::Asks::default();
        let height =
            (ctx.content_rect().height() - WINDOW_ROOM).clamp(MIN_WINDOW_HEIGHT, MAX_WINDOW_HEIGHT);
        design::window("Protein Domains")
            .open(&mut open)
            .default_width(680.0)
            .default_height(height)
            .vscroll(true)
            .show(ctx, |ui| {
                let w = &mut self.domains_win;
                ui.label(
                    "Which domains a protein has: HMMER searches the sequence against Pfam, a database of protein families, and the \
                     domains it finds are drawn along the sequence. HMMER and Pfam are other projects' programs and data; they must \
                     be installed on the computer that runs the search (docs/DOMAINS.md).",
                );
                ui.horizontal(|ui| {
                    if ui.button("Use an example").clicked() {
                        w.text = EXAMPLE.to_string();
                    }
                    if ui.button("Clear").clicked() {
                        let runner = w.runner;
                        *w = DomainsWin { open: true, runner, ..DomainsWin::default() };
                    }
                });
                ScrollArea::vertical().max_height(110.0).show(ui, |ui| {
                    ui.add(
                        egui::TextEdit::multiline(&mut w.text)
                            .font(egui::TextStyle::Monospace)
                            .hint_text("Paste protein sequences: FASTA, or one on each line")
                            .desired_rows(4)
                            .desired_width(f32::INFINITY),
                    );
                });
                ui.horizontal_wrapped(|ui| {
                    ui.radio_value(&mut w.curated, true, "Pfam's own thresholds")
                        .on_hover_text("The curated gathering threshold of each family, which Pfam uses to say that a sequence belongs to it.");
                    ui.radio_value(&mut w.curated, false, "E-value at most");
                    ui.add_enabled(
                        !w.curated,
                        egui::DragValue::new(&mut w.evalue)
                            .range(1e-30..=Cutoff::MAX_EVALUE)
                            .speed(1e-4)
                            .custom_formatter(|v, _| format!("{v:e}")),
                    );
                });
                compute_picker(ui, &look, view);
                ui.horizontal(|ui| {
                    start = ui
                        .add_enabled(!w.running(), egui::Button::new(RichText::new("Find domains").semibold()))
                        .clicked();
                    if w.running() {
                        ui.spinner();
                    }
                });
                if !w.message.is_empty() {
                    ui.label(&w.message);
                }
                if let Some((made, _)) = &w.made {
                    asks = made_ui::section(ui, &look, made, || w.as_it_is(), !w.running());
                }
                let Some(result) = w.result.clone() else {
                    return;
                };
                ui.separator();
                if ui.button("Copy as a table").clicked() {
                    ui.ctx().copy_text(to_table(&result.hits, &result.names));
                }
                draw_tracks(ui, &look, &result);
                if result.names.len() > MOST_DRAWN {
                    ui.label(format!("The first {MOST_DRAWN} of {} sequences are drawn; the table to copy has them all.", result.names.len()));
                }
                if result.hits.iter().enumerate().any(|(i, a)| result.hits[i + 1..].iter().any(|b| a.overlaps(b))) {
                    ui.label(
                        RichText::new(
                            "Some domains overlap: HMMER reports every family that matches, and Pfam's rule that families of one \
                             clan may not overlap is not applied. They are drawn one under the other.",
                        )
                        .small(),
                    );
                }
                if !result.hits.is_empty() {
                    ScrollArea::vertical().id_salt("domains_rows").max_height(260.0).show(ui, |ui| table(ui, &result));
                    if result.hits.len() > MOST_ROWS {
                        ui.label(format!("The first {MOST_ROWS} of {} domains; the table to copy has them all.", result.hits.len()));
                    }
                }
            });
        self.domains_win.open = open;
        let w = &mut self.domains_win;
        if let (true, Some((made, _))) = (asks.copy, &w.made) {
            ctx.copy_text(made.to_text());
            w.message = "Copied the record of how the domains were found.".into();
        }
        self.checks_win.open |= asks.checks;
        if let (true, Some((_, kept))) = (asks.again, &w.made) {
            // with the threshold the domains that are shown were found with
            (w.curated, w.evalue) = *kept;
        }
        if start || asks.again {
            w.start(ctx, place);
        }
    }
}
