//! Tools → Read Quality and Trimming: FastQC-style statistics and fastp-style cleaning of FASTQ
//! reads, with the clean reads handed on to the contig assembler.

use crate::seq::fastq::{clean, parse_fastq, stats, to_fastq, Cleaning, FastqRead, Report, Stats};
use crate::ui::app::{App, Tab};
use crate::ui::dialog::{self, Wanted};
use crate::ui::platform::route::Source;
use eframe::egui::{self, pos2, vec2, Color32, RichText, Sense, Shape, Stroke, Ui};

const PLOT_HEIGHT: f32 = 80.0;
const PLOT_MAX_WIDTH: f32 = 420.0;
/// The top of the quality plot, in Phred units.
const PLOT_TOP_QUALITY: f64 = 40.0;

#[derive(Default)]
pub struct ReadsWin {
    pub open: bool,
    pub pasted: String,
    pub reads: Vec<FastqRead>,
    pub source: String,
    pub cleaning: Cleaning,
    pub before: Stats,
    pub cleaned: Vec<FastqRead>,
    pub after: Option<(Stats, Report)>,
    pub message: String,
}

fn summary(s: &Stats) -> String {
    format!(
        "{} reads, {} bases, lengths {}–{}, GC {:.1} %, mean quality {:.1}, Q20 {:.1} %, Q30 {:.1} %",
        s.reads, s.bases, s.min_length, s.max_length, s.gc_percent, s.mean_quality, s.q20_percent, s.q30_percent
    )
}

/// Mean quality along the reads: one line per series, guides at Q20 and Q30.
fn quality_plot(ui: &mut Ui, series: &[(&[f64], Color32)]) {
    let width = ui.available_width().min(PLOT_MAX_WIDTH);
    let (rect, _) = ui.allocate_exact_size(vec2(width, PLOT_HEIGHT), Sense::hover());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 2.0, ui.visuals().extreme_bg_color);
    let y_of =
        |q: f64| rect.bottom() - (q / PLOT_TOP_QUALITY).clamp(0.0, 1.0) as f32 * rect.height();
    for guide in [20.0, 30.0] {
        let stroke = Stroke::new(1.0, ui.visuals().weak_text_color().gamma_multiply(0.5));
        painter.line_segment(
            [
                pos2(rect.left(), y_of(guide)),
                pos2(rect.right(), y_of(guide)),
            ],
            stroke,
        );
    }
    for (values, color) in series {
        let last = values.len().saturating_sub(1).max(1) as f32;
        let points: Vec<_> = values
            .iter()
            .enumerate()
            .map(|(i, &q)| pos2(rect.left() + i as f32 / last * rect.width(), y_of(q)))
            .collect();
        painter.add(Shape::line(points, Stroke::new(1.5, *color)));
    }
}

impl App {
    fn load_fastq_text(&mut self, text: &str, source: &str) {
        let w = &mut self.reads_win;
        match parse_fastq(text) {
            Ok(reads) => {
                w.before = stats(&reads);
                w.message = format!("Loaded {} reads from {source}.", reads.len());
                (w.reads, w.source, w.cleaned, w.after) =
                    (reads, source.to_string(), Vec::new(), None);
            }
            Err(e) => w.message = e,
        }
    }

    fn open_fastq(&mut self) {
        let ask = dialog::FileDialog::new()
            .set_title("Open FASTQ reads")
            .add_filter("FASTQ", &["fastq", "fq", "txt"]);
        if let Some(path) = dialog::pick_one(Wanted::Fastq, ask) {
            self.load_fastq(Source::Path(&path));
        }
    }

    /// Read a FASTQ file into the window.
    pub(in crate::ui) fn load_fastq(&mut self, file: Source) {
        match file.text() {
            Ok(text) => {
                let name = file.path().file_name().and_then(|n| n.to_str());
                self.load_fastq_text(&text, name.unwrap_or("reads"));
            }
            Err(e) => {
                self.reads_win.message = format!("Couldn't read {}: {e}", file.shown());
            }
        }
    }

    pub fn clean_reads(&mut self) {
        let w = &mut self.reads_win;
        let (cleaned, report) = clean(&w.reads, &w.cleaning);
        w.after = Some((stats(&cleaned), report));
        w.message = format!("{} of {} reads kept.", cleaned.len(), w.reads.len());
        w.cleaned = cleaned;
    }

    fn reads_to_assembler(&mut self) {
        let items: Vec<(String, Vec<u8>)> = self
            .reads_win
            .cleaned
            .iter()
            .map(|r| {
                (
                    r.name
                        .split_whitespace()
                        .next()
                        .unwrap_or("read")
                        .to_string(),
                    r.seq.clone(),
                )
            })
            .collect();
        let added = self.assembly.add_reads(items);
        self.tab = Tab::Assembly;
        self.say(format!(
            "Added {added} cleaned reads to the Contig Assembler."
        ));
    }

    fn save_cleaned(&mut self) {
        let ask = dialog::FileDialog::new().set_title("Save cleaned reads");
        let text = to_fastq(&self.reads_win.cleaned);
        if let Some(saved) = dialog::save_bytes(ask, "cleaned.fastq", text.as_bytes()) {
            self.reads_win.message = saved.message;
        }
    }

    /// Open the window with FASTQ text already loaded (the demo, or a pasted record set).
    pub fn show_reads(&mut self, text: &str, source: &str) {
        self.reads_win.open = true;
        self.load_fastq_text(text, source);
    }

    pub fn reads_window(&mut self, ctx: &egui::Context) {
        if !self.reads_win.open {
            return;
        }
        let mut open = true;
        let mut clicks = ReadsClicks::default();
        crate::ui::design::window("Read Quality and Trimming")
            .open(&mut open)
            .default_size([560.0, 560.0])
            .show(ctx, |ui| clicks = self.reads_body(ui));
        self.apply_reads_clicks(clicks);
        self.reads_win.open = open;
    }

    /// The controls and results of the read tools, for the window and for the Signals page.
    pub(in crate::ui) fn reads_body(&mut self, ui: &mut Ui) -> ReadsClicks {
        let mut clicks = ReadsClicks::default();
        let w = &mut self.reads_win;
        ui.horizontal(|ui| {
            clicks.load = ui.button("Open FASTQ…").clicked();
            ui.label("FASTQ with Phred+33 qualities, four lines to a read");
        });
        ui.collapsing("Paste FASTQ text", |ui| {
            ui.add(
                egui::TextEdit::multiline(&mut w.pasted)
                    .font(egui::TextStyle::Monospace)
                    .desired_rows(4)
                    .desired_width(f32::INFINITY),
            );
            clicks.use_pasted = ui.button("Use pasted text").clicked();
        });
        if !w.message.is_empty() {
            ui.label(&w.message);
        }
        if w.reads.is_empty() {
            return clicks;
        }
        ui.separator();
        ui.label(RichText::new(format!("{}: before", w.source)).strong());
        ui.label(summary(&w.before));
        let blue = Color32::from_rgb(60, 110, 200);
        let green = Color32::from_rgb(40, 160, 80);
        let after = w.after.as_ref().map(|(s, _)| s.per_position.clone());
        let mut lines: Vec<(&[f64], Color32)> = vec![(&w.before.per_position, blue)];
        if let Some(a) = &after {
            lines.push((a, green));
        }
        quality_plot(ui, &lines);
        ui.label(
            RichText::new(
                "Mean quality along the reads (blue: before, green: after); lines at Q20 and Q30",
            )
            .small(),
        );
        ui.separator();
        ui.label(RichText::new("Cleaning (fastp's defaults)").strong());
        let c = &mut w.cleaning;
        ui.horizontal(|ui| {
            ui.label("Adapter");
            let mut text = String::from_utf8_lossy(&c.adapter).into_owned();
            if ui
                .add(egui::TextEdit::singleline(&mut text).desired_width(260.0))
                .changed()
            {
                c.adapter = text
                    .bytes()
                    .filter(u8::is_ascii_alphabetic)
                    .map(|b| b.to_ascii_uppercase())
                    .collect();
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.checkbox(&mut c.cut_front, "cut front");
            ui.checkbox(&mut c.cut_tail, "cut tail");
            ui.checkbox(&mut c.cut_right, "cut right");
            ui.label("window");
            ui.add(egui::DragValue::new(&mut c.window).range(1..=1000));
            ui.label("mean quality");
            ui.add(egui::DragValue::new(&mut c.window_quality).range(1..=30));
        });
        ui.horizontal_wrapped(|ui| {
            ui.label("low quality below Q");
            ui.add(egui::DragValue::new(&mut c.qualified_quality).range(0..=93));
            ui.label("drop above");
            ui.add(
                egui::DragValue::new(&mut c.unqualified_percent)
                    .range(0..=100)
                    .suffix(" % low"),
            );
            ui.label("N limit");
            ui.add(egui::DragValue::new(&mut c.n_base_limit).range(0..=50));
            ui.label("mean Q at least");
            ui.add(egui::DragValue::new(&mut c.average_quality).range(0..=93));
            ui.label("length at least");
            ui.add(egui::DragValue::new(&mut c.min_length).range(1..=10_000));
        });
        ui.horizontal(|ui| {
            clicks.run = ui.button(RichText::new("Clean reads").strong()).clicked();
            if ui.button("fastp defaults").clicked() {
                *c = Cleaning::default();
            }
        });
        if let Some((s, r)) = &w.after {
            ui.separator();
            ui.label(RichText::new("After cleaning").strong());
            ui.label(summary(s));
            ui.label(format!(
                "Dropped: {} low quality, {} too many N, {} too short. Adapter found in {} reads.",
                r.too_many_low_quality, r.too_many_n, r.too_short, r.adapter_trimmed
            ));
            ui.horizontal(|ui| {
                clicks.assemble = ui.button("Add to Contig Assembler").clicked();
                clicks.save = ui.button("Save cleaned FASTQ…").clicked();
            });
        }
        clicks
    }

    pub(in crate::ui) fn apply_reads_clicks(&mut self, c: ReadsClicks) {
        if c.load {
            self.open_fastq();
        }
        if c.use_pasted {
            let text = std::mem::take(&mut self.reads_win.pasted);
            self.load_fastq_text(&text, "pasted text");
        }
        if c.run {
            self.clean_reads();
        }
        if c.assemble {
            self.reads_to_assembler();
        }
        if c.save {
            self.save_cleaned();
        }
    }
}

/// What was pressed in the read tools this frame.
#[derive(Default)]
pub struct ReadsClicks {
    load: bool,
    use_pasted: bool,
    run: bool,
    assemble: bool,
    save: bool,
}
