//! Signals: raw sequencing signals. Sanger chromatograms (.ab1 and .scf) are read and drawn with their base calls and
//! quality by the rules of the original's trace window, trimmed at the weak ends and sent on to Sequences or
//! Alignment or exported as SCF; sequencing reads (FASTQ) are summarised and cleaned with the read tools.

use super::chromatogram::{self, shown_base, View};
use super::route::{self, Source};
use super::Workspace;
use crate::pipeline::steps::fasta_text;
use crate::seq::abif::{self, Trace};
use crate::seq::dna::revcom;
use crate::seq::doc::Doc;
use crate::seq::scf;
use crate::ui::app::{Action, App};
use crate::ui::design::{self, Look, Tone};
use crate::ui::dialog::{self, Wanted};
use eframe::egui::{self, RichText, ScrollArea, Ui};

const OPTIONS_WIDTH: f32 = 300.0;
/// The quality that trimming starts from: a base this good has a 1 in 20 chance of being wrong.
const DEFAULT_TRIM_QUALITY: f32 = 13.0;
/// The original's scales (`abi_x_scale`, `abi_scale`): across 0.5 to 5 times and up 1 to 5 times, both 1 at first.
const X_SCALE_RANGE: std::ops::RangeInclusive<f32> = 0.5..=5.0;
const Y_SCALE_RANGE: std::ops::RangeInclusive<f32> = 1.0..=5.0;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Tab {
    Sanger,
    Reads,
}

pub struct Loaded {
    pub name: String,
    pub trace: Trace,
}

pub struct SignalsPage {
    tab: Tab,
    pub loaded: Option<Loaded>,
    x_scale: f32,
    y_scale: f32,
    /// Which of the A, C, G, T traces are drawn (by the letter the file gives the channel).
    shown: [bool; 4],
    /// The graph of how sure each call was, behind the traces (the original's "Toggle Phred", on when a trace opens).
    quality_graph: bool,
    trim: bool,
    trim_quality: f32,
    /// Mirror the picture and complement the letters, and hand on the reverse complement of the bases.
    pub reverse: bool,
    selection: Option<(usize, usize)>,
    pub message: String,
}

impl Default for SignalsPage {
    fn default() -> Self {
        SignalsPage {
            tab: Tab::Sanger,
            loaded: None,
            x_scale: 1.0,
            y_scale: 1.0,
            shown: [true; 4],
            quality_graph: true,
            trim: true,
            trim_quality: DEFAULT_TRIM_QUALITY,
            reverse: false,
            selection: None,
            message: String::new(),
        }
    }
}

impl SignalsPage {
    /// Show a trace and forget what was picked on the one before.
    pub fn load(&mut self, name: &str, trace: Trace) {
        self.message = format!(
            "{name}: {} bases, {} scans.",
            trace.bases.len(),
            trace.scans()
        );
        self.loaded = Some(Loaded {
            name: name.to_string(),
            trace,
        });
        (self.selection, self.tab) = (None, Tab::Sanger);
    }

    /// Save the whole trace as an SCF file, as the original's "Export as SCF" does (`abi_save_scf`): not only the
    /// bases that are kept or picked.
    fn export_scf(&mut self) {
        let Some(loaded) = &self.loaded else {
            return;
        };
        let stem = std::path::Path::new(&loaded.name)
            .file_stem()
            .map_or("trace".to_string(), |s| s.to_string_lossy().into_owned());
        let ask = dialog::FileDialog::new()
            .set_title("Save SCF File")
            .add_filter("Sequence Coloring Format", &["scf"]);
        let name = format!("{stem}.scf");
        if let Some(saved) = dialog::save_bytes(ask, &name, &scf::to_scf(&loaded.trace)) {
            self.message = saved.message;
        }
    }

    /// The error chance that matches the trimming quality.
    fn cutoff(&self) -> f64 {
        10f64.powf(-f64::from(self.trim_quality) / 10.0)
    }

    /// The bases trimming keeps, when it is on and the file has qualities.
    fn kept(&self) -> Option<(usize, usize)> {
        let loaded = self.loaded.as_ref()?;
        (self.trim && loaded.trace.has_quality())
            .then(|| loaded.trace.good_stretch(self.cutoff()))
            .flatten()
    }

    /// The bases to use: what was picked, else what trimming keeps, else all of them.
    pub fn range(&self) -> Option<(usize, usize)> {
        let n = self.loaded.as_ref()?.trace.bases.len();
        Some(self.selection.or_else(|| self.kept()).unwrap_or((0, n)))
    }

    /// The name and bases to hand on.
    pub fn sequence(&self) -> Option<(String, Vec<u8>)> {
        let loaded = self.loaded.as_ref()?;
        let (from, to) = self.range()?;
        let n = loaded.trace.bases.len();
        let mut bases = loaded.trace.bases.get(from..to.min(n))?.to_vec();
        let mut name = if loaded.trace.sample.is_empty() {
            loaded.name.clone()
        } else {
            loaded.trace.sample.clone()
        };
        if (from, to) != (0, n) {
            name.push_str(&format!(" bases {}-{to}", from + 1));
        }
        if self.reverse {
            bases = revcom(&bases);
            name.push_str(" (reverse complement)");
        }
        Some((name, bases))
    }
}

#[derive(Default)]
struct Clicks {
    open: bool,
    example: bool,
    use_in_sequences: bool,
    to_alignment: bool,
    copy: bool,
    save: bool,
    export_scf: bool,
}

impl App {
    pub(super) fn signals_page(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let look = design::look(self.prefs.dark);
        let mut c = Clicks::default();
        let mut reads_clicks = None;
        design::page(ui, |ui| {
            let mut tab = self.platform.signals.tab;
            design::page_header(
                ui,
                &look,
                "Raw signals",
                "Sanger chromatograms and sequencing reads: check, trim and send on",
                |ui| {
                    ui.selectable_value(&mut tab, Tab::Reads, "Reads (FASTQ)");
                    ui.selectable_value(&mut tab, Tab::Sanger, "Sanger traces");
                },
            );
            self.platform.signals.tab = tab;
            match tab {
                Tab::Sanger => sanger_tab(ui, &look, &mut self.platform.signals, &mut c),
                Tab::Reads => {
                    ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            design::card(ui, &look, |ui| reads_clicks = Some(self.reads_body(ui)));
                        });
                }
            }
        });
        if let Some(clicks) = reads_clicks {
            self.apply_reads_clicks(clicks);
        }
        self.apply_signal_clicks(c, ui.ctx(), out);
    }

    fn apply_signal_clicks(&mut self, c: Clicks, ctx: &egui::Context, out: &mut Vec<Action>) {
        if c.open {
            self.open_trace_file();
        }
        if c.example {
            self.platform
                .signals
                .load("example trace", crate::sample::demo_trace());
        }
        if c.export_scf {
            self.platform.signals.export_scf();
        }
        let Some((name, bases)) = self.platform.signals.sequence() else {
            return;
        };
        if c.use_in_sequences {
            let mut doc = Doc::new(name.clone(), bases.clone());
            doc.definition = format!("From a Sanger trace: {name}.");
            doc.dirty = true;
            self.add_doc(doc);
            out.push(Action::Go(Workspace::Sequences));
            self.platform.signals.message = format!("Opened {} bases in Sequences.", bases.len());
        }
        if c.to_alignment {
            self.platform.alignment.add(&name, &bases);
            let count = self.platform.alignment.inputs.len();
            self.platform.signals.message =
                format!("Sent to Alignment, which now has {count} sequences.");
        }
        if c.copy {
            ctx.copy_text(fasta_text(&[(name.clone(), bases.clone())]));
            self.platform.signals.message = "Copied as FASTA.".into();
        }
        if c.save {
            let ask = dialog::FileDialog::new().set_title("Save the sequence");
            let fasta = fasta_text(&[(name, bases)]);
            if let Some(saved) = dialog::save_bytes(ask, "trace.fasta", fasta.as_bytes()) {
                self.platform.signals.message = saved.message;
            }
        }
    }

    fn open_trace_file(&mut self) {
        let ask = dialog::FileDialog::new()
            .set_title("Open a sequencing trace")
            .add_filter("Sanger traces", &route::TRACE_EXTENSIONS);
        if let Some(path) = dialog::pick_one(Wanted::Trace, ask) {
            self.load_trace(Source::Path(&path));
        }
    }

    /// Read a trace file (ABIF or SCF) and show it.
    pub(in crate::ui) fn load_trace(&mut self, file: Source) {
        let name = file
            .path()
            .file_name()
            .map_or("trace".to_string(), |n| n.to_string_lossy().into_owned());
        let page = &mut self.platform.signals;
        match file
            .bytes()
            .map_err(|e| format!("Could not read {}: {e}", file.shown()))
            .and_then(|b| abif::read(&b))
        {
            Ok(trace) => page.load(&name, trace),
            Err(why) => page.message = format!("{name}: {why}"),
        }
    }
}

fn sanger_tab(ui: &mut Ui, look: &Look, page: &mut SignalsPage, c: &mut Clicks) {
    egui::Panel::left("signals_options")
        .exact_size(OPTIONS_WIDTH)
        .resizable(false)
        .frame(egui::Frame::new())
        .show(ui, |ui| {
            ScrollArea::vertical()
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    design::card(ui, look, |ui| options_panel(ui, look, page, c));
                });
        });
    egui::CentralPanel::default()
        .frame(egui::Frame::new().inner_margin(egui::Margin {
            left: design::PAGE_MARGIN,
            ..egui::Margin::ZERO
        }))
        .show(ui, |ui| viewer_panel(ui, look, page));
}

/// What the original's File Info window says about a trace (`abi_file_info`, AppMain.tcl:31840, 0x176398): the sample,
/// the signal strength of each channel and when the run began and ended; the instrument is added. Only what the file
/// gives is listed.
fn file_info(t: &Trace) -> Vec<(&'static str, String)> {
    let mut lines = Vec::new();
    if !t.sample.is_empty() {
        lines.push(("Sample", t.sample.clone()));
    }
    if !t.instrument.is_empty() {
        lines.push(("Instrument", t.instrument.clone()));
    }
    if t.strength.iter().any(|&s| s != 0) {
        let each: Vec<String> = t
            .order
            .iter()
            .zip(t.strength)
            .map(|(&b, s)| format!("{}:{s}", char::from(b)))
            .collect();
        lines.push(("Scaling", each.join(", ")));
    }
    if let Some(s) = t.started {
        lines.push(("Run began", s.long_text()));
    }
    if let Some(s) = t.ended {
        lines.push(("Run ended", s.long_text()));
    }
    lines
}

fn options_panel(ui: &mut Ui, look: &Look, page: &mut SignalsPage, c: &mut Clicks) {
    design::section(ui, look, "Trace");
    ui.horizontal_wrapped(|ui| {
        c.open = ui.button("Open trace…").clicked();
        c.example = ui.button("Example").clicked();
    });
    let Some(loaded) = &page.loaded else {
        ui.add_space(4.0);
        ui.label(
            RichText::new(
                "Open a trace from an Applied Biosystems sequencer (.ab1) or an SCF file, or try the example.",
            )
            .size(12.5)
            .color(look.dim),
        );
        if !page.message.is_empty() {
            ui.label(RichText::new(&page.message).size(12.0).color(look.bad));
        }
        return;
    };
    let t = &loaded.trace;
    ui.horizontal_wrapped(|ui| {
        design::chip(ui, look, &format!("{} bases", t.bases.len()), Tone::Accent);
        design::chip(ui, look, &format!("{} scans", t.scans()), Tone::Neutral);
    });
    design::section(ui, look, "File info");
    for (label, text) in file_info(t) {
        ui.label(
            RichText::new(format!("{label}: {text}"))
                .size(12.5)
                .color(look.dim),
        );
    }
    design::section(ui, look, "Quality");
    if t.has_quality() {
        let mean =
            t.quality.iter().map(|&q| f64::from(q)).sum::<f64>() / t.quality.len().max(1) as f64;
        ui.horizontal_wrapped(|ui| {
            design::chip(ui, look, &format!("mean Q {mean:.1}"), Tone::Neutral);
            design::chip(
                ui,
                look,
                &format!("Q20 {:.0} %", t.percent_at_least(20)),
                Tone::Good,
            );
            design::chip(
                ui,
                look,
                &format!("Q30 {:.0} %", t.percent_at_least(30)),
                Tone::Good,
            );
        });
    } else {
        ui.label(RichText::new("This instrument gives no quality values, so the ends cannot be trimmed by quality.").size(12.5).color(look.dim));
    }
    ui.add_enabled_ui(t.has_quality(), |ui| {
        ui.checkbox(&mut page.trim, "Trim the weak ends");
        ui.add_enabled(
            page.trim,
            egui::Slider::new(&mut page.trim_quality, 5.0..=30.0)
                .text("Q")
                .fixed_decimals(0),
        )
        .on_hover_text(
            "A base counts as good when its chance of being wrong is below this quality's",
        );
    });
    let n = t.bases.len();
    if let Some((from, to)) = page.kept() {
        ui.label(
            RichText::new(format!("Keeping bases {}–{to} of {n}", from + 1))
                .size(12.5)
                .color(look.dim),
        );
    }

    design::section(ui, look, "View");
    ui.add(egui::Slider::new(&mut page.x_scale, X_SCALE_RANGE).text("x scale"))
        .on_hover_text("Across: 0.5 to 5 times");
    ui.add(egui::Slider::new(&mut page.y_scale, Y_SCALE_RANGE).text("y scale"))
        .on_hover_text(
            "Up: at 1 the tallest peak of the four traces reaches the top of the trace area",
        );
    let reverse = page.reverse;
    ui.horizontal(|ui| {
        for (k, base) in b"ACGT".iter().enumerate() {
            // in a reverse-complemented view a channel is drawn as the base it pairs with
            let shown = shown_base(*base, reverse);
            let colour = design::trace_color(shown, look.dark);
            ui.checkbox(
                &mut page.shown[k],
                RichText::new(char::from(shown).to_string())
                    .strong()
                    .color(colour),
            );
        }
    });
    ui.checkbox(&mut page.quality_graph, "Quality graph")
        .on_hover_text(
        "How sure the machine was of each base, behind the traces (Toggle Phred in the original)",
    );
    ui.label(
        RichText::new("Drag across bases to pick them.")
            .size(12.0)
            .color(look.dim),
    );

    design::section(ui, look, "Use the sequence");
    ui.checkbox(&mut page.reverse, "Reverse-complement")
        .on_hover_text(
        "Mirrors the picture and complements the letters, for a read made with the reverse primer",
    );
    if let Some((from, to)) = page.range() {
        let what = if page.selection.is_some() {
            "the picked bases"
        } else if page.kept().is_some() {
            "the kept bases"
        } else {
            "all the bases"
        };
        ui.label(
            RichText::new(format!("{what}: {}–{to} ({} bases)", from + 1, to - from))
                .size(12.5)
                .color(look.dim),
        );
    }
    if page.selection.is_some() && ui.small_button("Clear the pick").clicked() {
        page.selection = None;
    }
    c.use_in_sequences = ui.button("Open in Sequences").clicked();
    c.to_alignment = ui.button("Send to Alignment").clicked();
    ui.horizontal_wrapped(|ui| {
        c.copy = ui.button("Copy FASTA").clicked();
        c.save = ui.button("Save FASTA…").clicked();
    });
    c.export_scf = ui
        .button("Export as SCF…")
        .on_hover_text("The whole trace, as the original's Export as SCF writes it")
        .clicked();
    if !page.message.is_empty() {
        ui.add_space(4.0);
        ui.label(RichText::new(&page.message).size(12.0).color(look.dim));
    }
}

fn viewer_panel(ui: &mut Ui, look: &Look, page: &mut SignalsPage) {
    let kept = page.kept();
    let Some(loaded) = &page.loaded else {
        design::card(ui, look, |ui| {
            ui.label(
                RichText::new("No trace yet")
                    .size(17.0)
                    .strong()
                    .color(look.text),
            );
            ui.label(
                RichText::new(
                    "Open a .ab1 or .scf file, or drop one on the window. The four colour traces appear here with the \
                     bases the machine called, a graph of how sure it was of each, and the weak ends shaded.",
                )
                .color(look.dim),
            );
        });
        return;
    };
    let view = View {
        trace: &loaded.trace,
        look,
        x_scale: page.x_scale,
        y_scale: page.y_scale,
        shown: page.shown,
        kept,
        quality_graph: page.quality_graph,
        reverse: page.reverse,
    };
    let touched = egui::Frame::new()
        .fill(look.surface)
        .stroke(egui::Stroke::new(1.0, look.border))
        .corner_radius(egui::CornerRadius::same(design::CARD_RADIUS))
        .show(ui, |ui| chromatogram::show(ui, &view, page.selection))
        .inner;
    if touched.picked.is_some() {
        page.selection = touched.picked;
    } else if touched.cleared {
        page.selection = None;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn page() -> SignalsPage {
        let mut p = SignalsPage::default();
        p.load("demo.ab1", crate::sample::demo_trace());
        p
    }

    #[test]
    fn trimming_keeps_the_good_middle_and_a_pick_overrides_it() {
        let mut p = page();
        let n = p.loaded.as_ref().expect("loaded").trace.bases.len();
        let (from, to) = p.range().expect("a range");
        assert!(
            from > 3 && to < n,
            "the weak ends are trimmed: {from}..{to} of {n}"
        );
        p.trim = false;
        assert_eq!(p.range(), Some((0, n)));
        p.selection = Some((10, 20));
        assert_eq!(p.range(), Some((10, 20)));
    }

    #[test]
    fn the_sequence_to_hand_on_is_named_for_what_was_taken() {
        let mut p = page();
        p.trim = false;
        let (name, bases) = p.sequence().expect("a sequence");
        assert_eq!(name, "pDemo_T7_forward (example)");
        assert_eq!(bases.len(), 480);
        p.selection = Some((4, 12));
        p.reverse = true;
        let (name, bases) = p.sequence().expect("a sequence");
        assert_eq!(
            name,
            "pDemo_T7_forward (example) bases 5-12 (reverse complement)"
        );
        let whole = p.loaded.as_ref().expect("loaded").trace.bases[4..12].to_vec();
        assert_eq!(bases, revcom(&whole));
        assert!(SignalsPage::default().sequence().is_none());
    }

    #[test]
    fn the_trimming_cutoff_is_the_chance_of_a_wrong_base() {
        let p = SignalsPage::default();
        assert!((p.cutoff() - 0.05).abs() < 0.001, "{}", p.cutoff());
    }
}
