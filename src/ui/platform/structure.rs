//! Structure: predict the 3D structure of a protein from its amino-acid sequence with ESMFold (through the ESM
//! Metagenomic Atlas web service) or open a PDB file, and turn it round in a 3D viewer coloured by confidence,
//! secondary structure, chain or position along the chain.

use super::route::{self, Source};
use super::viewer3d::{self, Camera, Representation, Scheme, Style};
use super::workstation::{run_button, ComputeView, Place};
use super::Workspace;
use crate::pipeline::predict::{self, ESM_ATLAS_URL, MAX_LENGTH};
use crate::pipeline::steps::fasta_text;
use crate::pipeline::tools::Installed;
use crate::provenance::{Asked, Input, Provenance};
use crate::remote::tasks::Request;
use crate::remote::weight::Weight;
use crate::seq::dna::{find_orfs, revcom, translate, StartRule, Strand, GENETIC_CODES};
use crate::seq::doc::Doc;
use crate::seq::fold::Fold;
use crate::seq::pdb::parse;
use crate::ui::app::App;
use crate::ui::design::{self, Look, Tone};
use crate::ui::dialog::{self, Wanted};
use crate::ui::job::{Job, Poll};
use crate::ui::provenance as made_ui;
use eframe::egui::{self, Align2, Color32, FontId, RichText, ScrollArea, Sense, Ui};
use std::path::Path;
use std::time::Duration;
use web_time::Instant;

const STRIP_HEIGHT: f32 = 62.0;
const LETTER_WIDTH: f32 = 15.0;
const WORKING_REDRAW: Duration = Duration::from_millis(250);
/// The ORFs the open sequence is searched for must be at least this many bases long.
const MIN_ORF_BASES: usize = 90;
/// The viewer never gets less room than this, however small the window.
const MIN_VIEWER_HEIGHT: f32 = 160.0;
/// The example structure: ubiquitin, PDB entry 1UBQ.
const EXAMPLE_PDB: &str = include_str!("../../../assets/structure/1ubq.pdb");
/// The antibody D1.3 Fv (chains A and B) with hen egg-white lysozyme (chain C): PDB 1VFB, an example of an interface.
const ANTIBODY_PDB: &str = include_str!("../../../assets/structure/1vfb.pdb");

pub struct Loaded {
    pub name: String,
    /// The file as it was read, so saving gives back exactly what came in.
    text: String,
    pub fold: Fold,
    /// The bonds of every atom, found when the all-bonds drawing is first asked for.
    bonds: Option<Vec<(usize, usize)>>,
}

/// What a prediction gives: the structure as PDB text and what made it, or why there is none.
type Predicted = Result<(String, String), String>;

#[derive(Default)]
pub struct StructurePage {
    pub loaded: Option<Loaded>,
    /// Counted up whenever a structure is shown: the session keeps the structure, and tells by this that it is another one.
    shown: u64,
    /// The interface between two groups of chains: which chains, and what was measured.
    interface: interface::InterfacePage,
    sequence: String,
    style: Style,
    camera: Camera,
    spin: bool,
    /// A prediction that is going, and when it was asked for.
    job: Option<(Job<Predicted>, Instant)>,
    /// What the prediction that is going was asked.
    asked: Option<Asked>,
    /// The record of the prediction that is shown; none for a structure that was opened, which nothing here made.
    pub made: Option<Provenance>,
    pub message: String,
    /// The version of the alignment's scores (`AlignmentPage::revision`) that the colours by conservation were made from, if
    /// there are such colours: they are taken off when the scores change.
    conservation_from: Option<u64>,
}

impl StructurePage {
    /// For a picture: show the example of two chains with its interface measured and drawn.
    pub fn show_example_interface(&mut self) {
        self.load_antibody_example();
        if let Some(loaded) = self.loaded.as_mut() {
            self.interface.measure_now(loaded, &mut self.style);
        }
    }

    pub fn is_running(&self) -> bool {
        self.job.is_some()
    }

    /// Show a structure given as PDB text; a predicted one starts coloured by confidence.
    pub fn show_pdb(&mut self, name: &str, text: &str) -> Result<(), String> {
        let fold = Fold::new(parse(text)?);
        if fold.residues.is_empty() {
            return Err(
                "The file has atoms but no protein chain (no residue with N, CA, C and O).".into(),
            );
        }
        let s = fold.summary();
        self.style = Style {
            representation: match self.style.representation {
                Representation::Sticks if s.atoms > viewer3d::MAX_STICK_ATOMS => {
                    Representation::Ribbon
                }
                kept => kept,
            },
            scheme: if fold.predicted {
                Scheme::Confidence
            } else {
                Scheme::Secondary
            },
            picked: None,
        };
        self.camera = Camera::default();
        self.conservation_from = None;
        self.made = None;
        self.interface.reset_for(&fold.structure.chains());
        self.message = format!(
            "{name}: {} residues, {} chain(s), {} atoms.",
            s.residues, s.chains, s.atoms
        );
        self.loaded = Some(Loaded {
            name: name.to_string(),
            text: text.to_string(),
            fold,
            bonds: None,
        });
        self.shown += 1;
        Ok(())
    }

    /// Read the PDB file at `path` and show it.
    pub fn open_path(&mut self, path: &Path) {
        self.open_source(Source::Path(path));
    }

    /// Read a PDB file and show it.
    pub fn open_source(&mut self, file: Source) {
        let name = file
            .path()
            .file_name()
            .map_or("structure".to_string(), |n| {
                n.to_string_lossy().into_owned()
            });
        match file.text() {
            Ok(text) => {
                if let Err(e) = self.show_pdb(&name, &text) {
                    self.message = format!("{name}: {e}");
                }
            }
            Err(e) => self.message = format!("Could not read {}: {e}", file.shown()),
        }
    }

    /// Show the example structure, ubiquitin.
    pub fn load_example(&mut self) {
        self.sequence.clear();
        if let Err(e) = self.show_pdb("ubiquitin (PDB 1UBQ, example)", EXAMPLE_PDB) {
            self.message = e;
        }
    }

    /// Show the antibody example: D1.3 against lysozyme, with the antibody's two chains on one side and the antigen on the other.
    pub fn load_antibody_example(&mut self) {
        self.sequence.clear();
        match self.show_pdb(
            "antibody D1.3 with lysozyme (PDB 1VFB, example)",
            ANTIBODY_PDB,
        ) {
            Ok(()) => {
                self.interface.set_first_side(&['A', 'B']);
                self.message
                    .push_str(" The antibody is chains A and B, lysozyme chain C: see Interface.");
            }
            Err(e) => self.message = e,
        }
    }

    /// Draw the structure that is shown, turned and coloured as it is in the window and standing still: for a report.
    pub(in crate::ui) fn draw_still(&self, ui: &mut Ui, look: &Look) {
        let Some(loaded) = &self.loaded else {
            return;
        };
        let mut camera = self.camera;
        let bonds = loaded.bonds.as_deref().unwrap_or(&[]);
        viewer3d::show(
            ui,
            look,
            &loaded.fold,
            &mut camera,
            &self.style,
            bonds,
            false,
        );
    }

    /// What is written in the sequence box, checked.
    fn checked(&self) -> Result<Vec<u8>, String> {
        predict::clean_sequence(&self.sequence)
    }

    /// The sequence in the box and who would fold it now (`via`, as for `predictor`), to say what has changed since the
    /// prediction that is shown was made.
    pub fn as_it_is(&self, via: Option<bool>) -> (Vec<Input>, Vec<(String, String)>) {
        described(&self.checked().unwrap_or_default(), via)
    }
}

/// Who folds a sequence: the ESM Atlas web service, asked from this computer (`None`) or by the workstation (`Some(false)`),
/// or the workstation's own folding program (`Some(true)`).
pub fn predictor(via: Option<bool>) -> &'static str {
    match via {
        Some(true) => "the workstation's folding program",
        Some(false) | None => "ESMFold, through the ESM Atlas web service",
    }
}

/// Whether the prediction of this record was folded by who would fold one now.
fn made_by(made: &Provenance, via: Option<bool>) -> bool {
    made.parameters
        .iter()
        .any(|(key, who)| key == "predictor" && who == predictor(via))
}

/// The input and the setting of a prediction, as its record says them: the residues that are sent, and who folds them.
pub fn described(sequence: &[u8], via: Option<bool>) -> (Vec<Input>, Vec<(String, String)>) {
    (
        vec![Input::of("protein sequence", sequence)],
        vec![("predictor".to_string(), predictor(via).to_string())],
    )
}

/// The protein an open document codes for: the picked bases read from their first base, or else its longest ORF.
pub fn protein_of(doc: &Doc, selection: Option<(usize, usize)>, table: &str) -> Option<Vec<u8>> {
    let bases: Vec<u8> = match selection {
        Some((a, b)) if b > a => doc.seq.get(a..b)?.to_vec(),
        _ => {
            let longest = find_orfs(
                &doc.seq,
                doc.circular,
                Strand::Both,
                MIN_ORF_BASES,
                StartRule::Met,
                false,
            )
            .into_iter()
            .max_by_key(|o| o.length)?;
            // an ORF of a circular sequence may run past the end, so read it from two turns round
            let doubled = [doc.seq.as_slice(), doc.seq.as_slice()].concat();
            let top = doubled.get(longest.start..longest.start + longest.length)?;
            if longest.reverse {
                revcom(top)
            } else {
                top.to_vec()
            }
        }
    };
    let protein = translate(&bases, table);
    let upto: Vec<u8> = protein.into_iter().take_while(|&c| c != b'*').collect();
    (!upto.is_empty()).then_some(upto)
}

#[derive(Default)]
struct Clicks {
    predict: bool,
    from_document: bool,
    open: bool,
    example: bool,
    antibody_example: bool,
    save: bool,
    reset: bool,
    /// Colour the structure by the conservation of the alignment on the Alignment page.
    conservation: bool,
    interface: interface::InterfaceClicks,
    /// What was pressed under the record of the prediction.
    made: made_ui::Asks,
    /// Go to the Alignment page, where the column of the picked residue is selected.
    to_alignment: bool,
}

impl App {
    /// Collect a finished prediction.
    pub(super) fn structure_tick(&mut self, ctx: &egui::Context) {
        self.interface_tick(ctx);
        // colours made from scores that have since changed (another alignment, another score) say what is no longer shown
        let revision = self.platform.alignment.revision();
        let page = &mut self.platform.structure;
        if page.conservation_from.is_some_and(|from| from != revision) {
            page.conservation_from = None;
            if let Some(loaded) = page.loaded.as_mut() {
                loaded.fold.clear_conservation();
                if page.style.scheme == Scheme::Conservation {
                    page.style.scheme = if loaded.fold.predicted {
                        Scheme::Confidence
                    } else {
                        Scheme::Secondary
                    };
                }
            }
            page.message = "The alignment or its score changed, so the colours by conservation were taken off: press \"Colour by conservation\" again."
                .into();
        }
        let Some(polled) = page.job.as_ref().map(|(job, _)| job.poll()) else {
            return;
        };
        match polled {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(Ok((text, engine))) => {
                let seconds = page
                    .job
                    .take()
                    .map_or(0.0, |(_, started)| started.elapsed().as_secs_f64());
                let asked = page.asked.take();
                match page.show_pdb("ESMFold prediction", &text) {
                    Ok(()) => {
                        page.message = format!("Predicted in {seconds:.1} s. {}", page.message);
                        // the time is from the asking to the answer, the wait for the service or the queue with it
                        page.made = asked.map(|asked| asked.made(vec![engine], seconds));
                    }
                    Err(e) => {
                        page.message = format!(
                            "The service answered, but the answer is not a protein structure: {e}"
                        );
                    }
                }
            }
            Poll::Done(Err(why)) => (page.message, page.job) = (why, None),
            Poll::Died => {
                (page.message, page.job) =
                    ("The prediction stopped without an answer.".into(), None);
            }
        }
    }

    /// How big a prediction is for the choice of where it runs: the workstation is worth it when it folds on a GPU of its
    /// own; when it would only ask the web service, going there gains nothing (the choice of Workstation sends it anyway).
    fn fold_weight(&self) -> Weight {
        if self.platform.workstation.folds_on_gpu() {
            Weight::on_gpu()
        } else {
            Weight::default()
        }
    }

    pub(super) fn structure_page(&mut self, ui: &mut Ui) {
        let look = self.prefs.look();
        let weight = self.fold_weight();
        let place = self.platform.place(&weight);
        let curl = self.platform.installed_at("curl", &place);
        let searched = self.platform.programs_known_at(&place);
        // where the sequence would go: the workstation (which may fold it on its own GPU) or the web service from here
        let via = place
            .is_there()
            .then(|| self.platform.workstation.folds_on_gpu());
        let view = self.platform.workstation.view(&weight);
        let mut c = Clicks::default();
        let aligned = self.platform.alignment.conservation().is_some();
        let page = &mut self.platform.structure;
        let picked_before = page.picked();
        design::page(ui, |ui| {
            design::page_header(
                ui,
                &look,
                "Protein 3D structure",
                "Predict a structure from a sequence, or open one, and explore it in 3D",
                |ui| {
                    ui.add_enabled_ui(page.loaded.is_some(), |ui| {
                        c.save = ui.button("Save PDB…").clicked();
                        c.reset = ui.button("Reset view").clicked();
                        crate::ui::figure::button(ui);
                    });
                },
            );
            design::controls_panel(ui, &look, "structure_options", |ui| {
                let curl = Curl {
                    found: curl.is_some() || via == Some(true),
                    searched,
                    via,
                    view,
                };
                options_panel(ui, &look, page, curl, aligned, &mut c);
            });
            design::results_area(ui, |ui| viewer_panel(ui, &look, page, &mut c));
        });
        let coloured = c.conservation;
        self.apply_structure_clicks(c, curl, ui.ctx());
        // another residue, or none, by a click in the picture or on the strip, or by Reset view: the alignment follows.
        // Colouring by the alignment has settled the selection itself: the column that is selected there stays, whatever
        // it does to the residue that was picked here
        if !coloured && self.platform.structure.picked() != picked_before {
            self.platform.residue_picked();
        }
    }

    /// Colour the structure by the conservation of the columns of the alignment, through the rows that its chains match.
    pub(super) fn colour_structure_by_alignment(&mut self) {
        let alignment = &self.platform.alignment;
        let page = &mut self.platform.structure;
        let revision = alignment.revision();
        let Some((rows, scores, score)) = alignment.conservation() else {
            return;
        };
        let Some(loaded) = page.loaded.as_mut() else {
            return;
        };
        let matched = loaded.fold.apply_conservation(rows, scores);
        if matched.is_empty() {
            page.message = "No sequence of the alignment is this protein or a close relative of it (at least 70 % the same), so nothing was coloured."
                .into();
            return;
        }
        page.conservation_from = Some(revision);
        page.style.scheme = Scheme::Conservation;
        let parts: Vec<String> = matched
            .iter()
            .map(|m| {
                format!(
                    "chain {} from \"{}\", {:.0} % the same, {} of {} residues",
                    m.chain,
                    m.row,
                    100.0 * m.identity,
                    m.coloured,
                    m.residues
                )
            })
            .collect();
        page.message = format!("Conservation ({}): {}.", score.label(), parts.join("; "));
        // from now on the one selection is shown on both pages: the column that is selected there already, or else the
        // residue that is picked here
        match self.platform.alignment.selected_column() {
            Some(column) => self.platform.structure.pick_column(Some(column)),
            None => self.platform.residue_picked(),
        }
    }

    fn apply_structure_clicks(&mut self, c: Clicks, curl: Option<Installed>, ctx: &egui::Context) {
        self.apply_interface_clicks(c.interface, ctx);
        if c.from_document {
            self.protein_from_open_document();
        }
        if c.conservation {
            self.colour_structure_by_alignment();
        }
        let place = self.platform.place(&self.fold_weight());
        let gpu_fold = self.platform.workstation.folds_on_gpu();
        let page = &mut self.platform.structure;
        if let (true, Some(made)) = (c.made.copy, &page.made) {
            ctx.copy_text(made.to_text());
            page.message = "Copied the record of how the structure was predicted.".into();
        }
        self.checks_win.open |= c.made.checks;
        if c.predict || c.made.again {
            let asking = |sequence: &[u8], via: Option<bool>| {
                let (inputs, parameters) = described(sequence, via);
                Asked {
                    analysis: "Structure prediction",
                    inputs,
                    parameters,
                    ran_on: if via.is_some() {
                        "the workstation"
                    } else {
                        "this computer"
                    },
                }
            };
            match (page.checked(), place, curl) {
                (Err(why), _, _) => page.message = why,
                (Ok(sequence), Place::There(remote), _) => {
                    page.message = "Asking the workstation to fold it…".into();
                    page.asked = Some(asking(&sequence, Some(gpu_fold)));
                    let label = format!("Fold of {} residues", sequence.len());
                    let request = Request::Fold {
                        sequence,
                        local: gpu_fold,
                    };
                    page.job = Some((
                        Job::spawn_async(ctx, async move {
                            let folded = remote.run_async(&request, &label).await?.into_fold()?;
                            Ok((folded.pdb, folded.engine))
                        }),
                        Instant::now(),
                    ));
                }
                (Ok(sequence), Place::Here, Some(curl)) => {
                    page.message = "Asking the ESM Atlas service to fold it…".into();
                    page.asked = Some(asking(&sequence, None));
                    let path = curl.path;
                    page.job = Some((
                        Job::spawn(ctx, move || {
                            predict::predict(&path, ESM_ATLAS_URL, &sequence)
                                .map(|pdb| (pdb, "the ESM Atlas web service".to_string()))
                        }),
                        Instant::now(),
                    ));
                }
                (Ok(_), Place::Here, None) => {
                    page.message = "curl was not found; it is needed to ask the service.".into();
                }
            }
        }
        if c.example {
            page.load_example();
        }
        if c.antibody_example {
            page.load_antibody_example();
        }
        if c.open {
            let ask = dialog::FileDialog::new()
                .set_title("Open a protein structure")
                .add_filter("PDB files", &route::STRUCTURE_EXTENSIONS);
            if let Some(path) = dialog::pick_one(Wanted::Structure, ask) {
                page.open_path(&path);
            }
        }
        if c.reset {
            page.camera = Camera::default();
            page.style.picked = None;
        }
        if c.to_alignment {
            self.workspace = Workspace::Alignment;
        }
        let page = &mut self.platform.structure;
        if c.save {
            if let Some(loaded) = &page.loaded {
                let ask = dialog::FileDialog::new().set_title("Save the structure");
                if let Some(saved) =
                    dialog::save_bytes(ask, "structure.pdb", loaded.text.as_bytes())
                {
                    page.message = saved.message;
                }
            }
        }
    }

    /// Put the protein of the open sequence in the sequence box.
    fn protein_from_open_document(&mut self) {
        let Some(i) = self.current_doc() else {
            self.platform.structure.message = "No sequence is open.".into();
            return;
        };
        let dv = &self.docs[i];
        let table = GENETIC_CODES[self.prefs.code].table;
        let page = &mut self.platform.structure;
        match protein_of(&dv.doc, dv.selection(), table) {
            Some(protein) => {
                page.message = format!("Took {} residues from {}.", protein.len(), dv.doc.name);
                page.sequence = String::from_utf8_lossy(&protein).into_owned();
            }
            None => {
                page.message = "No reading frame starts with ATG and runs 90 bases or more; \
                                select the coding bases and try again."
                    .into();
            }
        }
    }
}

/// What the page knows about the program that asks the service.
#[derive(Clone, Copy)]
struct Curl {
    found: bool,
    /// The search for programs has finished, so a missing one is really missing.
    searched: bool,
    /// The sequence goes to the workstation, and whether that folds it on a GPU of its own.
    via: Option<bool>,
    /// Where it goes, and the choice that decides it.
    view: ComputeView,
}

/// The options, with what is being worked on first: a structure that is shown has its own sections above the ones for starting
/// another.
fn options_panel(
    ui: &mut Ui,
    look: &Look,
    page: &mut StructurePage,
    curl: Curl,
    aligned: bool,
    c: &mut Clicks,
) {
    if page.loaded.is_some() {
        open_section(ui, look, c);
        structure_section(ui, look, page, curl, aligned, c);
        predict_section(ui, look, page, curl, c);
    } else {
        predict_section(ui, look, page, curl, c);
        open_section(ui, look, c);
    }
}

fn predict_section(ui: &mut Ui, look: &Look, page: &mut StructurePage, curl: Curl, c: &mut Clicks) {
    design::section(ui, look, "Predict from a sequence");
    ui.add(
        egui::TextEdit::multiline(&mut page.sequence)
            .hint_text("Amino acids, or a FASTA protein")
            .font(egui::TextStyle::Monospace)
            .desired_rows(4)
            .desired_width(f32::INFINITY),
    );
    let checked = page.checked();
    match &checked {
        Ok(s) => design::chip(ui, look, &format!("{} residues", s.len()), Tone::Good),
        Err(_) if page.sequence.trim().is_empty() => {}
        Err(why) => {
            ui.label(RichText::new(why).size(12.0).color(look.bad));
        }
    }
    c.from_document = ui
        .button("Use the open sequence")
        .on_hover_text("The selected bases, or else the longest ORF, translated")
        .clicked();
    let state = (checked.is_ok() && curl.found, page.is_running());
    c.predict = run_button(ui, look, curl.view, "Predict structure", state);
    // right under the button, where a failure (the workstation cannot be reached, say) is seen
    if !page.message.is_empty() {
        ui.label(RichText::new(&page.message).size(12.0).color(look.dim));
    }
    if curl.searched && !curl.found {
        let place = if curl.via.is_some() {
            "the workstation"
        } else {
            "this computer"
        };
        ui.label(
            RichText::new(format!(
                "curl was not found on {place}, so the service cannot be asked."
            ))
            .size(12.0)
            .color(look.bad),
        );
    }
    let how = match curl.via {
        None => format!(
            "ESMFold runs on the ESM Metagenomic Atlas web service: the sequence is sent over the internet to \
             api.esmatlas.com (up to {MAX_LENGTH} residues; needs curl). Nothing else leaves this computer."
        ),
        Some(false) => format!(
            "The workstation asks the ESM Metagenomic Atlas web service: the sequence goes to the workstation, and \
             from there over the internet to api.esmatlas.com (up to {MAX_LENGTH} residues; needs curl there)."
        ),
        Some(true) => "The sequence goes to the workstation, which folds it on a GPU of its own with the program \
                       its owner set up. It is not sent anywhere else."
            .to_string(),
    };
    ui.label(RichText::new(how).size(11.5).color(look.warn));
}

fn open_section(ui: &mut Ui, look: &Look, c: &mut Clicks) {
    design::section(ui, look, "Open a structure");
    ui.horizontal_wrapped(|ui| {
        c.open = ui.button("Open PDB…").clicked();
        c.example = ui.button("Example").clicked();
        c.antibody_example = ui
            .button("Antibody")
            .on_hover_text("Example: an antibody (chains A and B) bound to lysozyme (chain C), PDB 1VFB, to try the Interface on")
            .clicked();
    });
}

fn structure_section(
    ui: &mut Ui,
    look: &Look,
    page: &mut StructurePage,
    curl: Curl,
    aligned: bool,
    c: &mut Clicks,
) {
    let Some(loaded) = &page.loaded else {
        return;
    };
    let s = loaded.fold.summary();
    design::section(ui, look, "This structure");
    ui.label(RichText::new(&loaded.name).size(12.5).color(look.text));
    if !loaded.fold.structure.title.is_empty() {
        ui.label(
            RichText::new(&loaded.fold.structure.title)
                .size(11.5)
                .color(look.dim),
        );
    }
    ui.horizontal_wrapped(|ui| {
        design::chip(ui, look, &format!("{} residues", s.residues), Tone::Accent);
        design::chip(ui, look, &format!("{} chain(s)", s.chains), Tone::Neutral);
        design::chip(ui, look, &format!("{} atoms", s.atoms), Tone::Neutral);
        design::chip(ui, look, &format!("helix {:.0} %", s.helix), Tone::Neutral);
        design::chip(
            ui,
            look,
            &format!("strand {:.0} %", s.strand),
            Tone::Neutral,
        );
        if let Some(confidence) = s.confidence {
            design::chip(ui, look, &format!("mean pLDDT {confidence:.0}"), Tone::Good);
        }
    });
    if let Some(made) = &page.made {
        // Run again folds what is in the box with who would fold it now: it is offered while that is who made this one
        let same = made_by(made, curl.via);
        let ready = page.checked().is_ok() && curl.found && !page.is_running();
        c.made = made_ui::section(ui, look, made, || page.as_it_is(curl.via), ready && same);
    }
    if ui
        .button("Copy sequence")
        .on_hover_text("As FASTA, one record for each chain")
        .clicked()
    {
        let stem = loaded.name.split_whitespace().next().unwrap_or("structure");
        let records: Vec<(String, Vec<u8>)> = loaded
            .fold
            .chain_sequences()
            .into_iter()
            .map(|(chain, letters)| (format!("{stem}_{chain}"), letters))
            .collect();
        ui.ctx().copy_text(fasta_text(&records));
        page.message = "Copied the sequence as FASTA.".into();
    }
    ui.add_space(6.0);
    c.conservation = ui
            .add_enabled(aligned, egui::Button::new("Colour by conservation"))
            .on_hover_text(
                "Line this structure up with the sequences of the alignment on the Alignment page, and colour each residue \
                 by how conserved its column is there (by the score chosen there)",
            )
            .on_disabled_hover_text(
                "Align sequences first, on the Alignment page: the family of this protein",
            )
            .clicked();
    interface::panel(ui, look, &mut page.interface, &mut c.interface);
}

/// A choice that is always for choosing.
const ALWAYS: (bool, &str) = (true, "");

/// One of a row of choices, which may not be for choosing now (`why_not` is told to the pointer then). The row wraps in a
/// narrow window, and a choice has to go to the next row whole, with the others: so it is one widget, never a room of its
/// own (`add_enabled_ui`), which stays where it began, pushes what comes after it out of the window and takes the height
/// of the page with it. That is what hid the structure in a small window until November 2026.
fn choice<V: PartialEq>(
    ui: &mut Ui,
    current: &mut V,
    value: V,
    text: &str,
    (enabled, why_not): (bool, &str),
) {
    let button = egui::Button::selectable(*current == value, text);
    let asked = ui.add_enabled(enabled, button);
    if asked.on_disabled_hover_text(why_not).clicked() {
        *current = value;
    }
}

/// The choices of how the structure is drawn and coloured, in a row above it; whether the link to the column of the picked
/// residue was clicked.
fn look_bar(ui: &mut Ui, look: &Look, style: &mut Style, spin: &mut bool, loaded: &Loaded) -> bool {
    let mut to_alignment = false;
    let (atoms, predicted) = (loaded.fold.structure.atoms.len(), loaded.fold.predicted);
    let heading =
        |ui: &mut Ui, text: &str| ui.label(RichText::new(text).size(12.0).color(look.dim));
    ui.horizontal_wrapped(|ui| {
        heading(ui, "Show");
        let how = &mut style.representation;
        choice(ui, how, Representation::Ribbon, "Ribbon", ALWAYS);
        choice(ui, how, Representation::Backbone, "Backbone", ALWAYS);
        choice(
            ui,
            &mut style.representation,
            Representation::Sticks,
            "All bonds",
            (
                atoms <= viewer3d::MAX_STICK_ATOMS,
                "Too many atoms to draw every bond",
            ),
        );
        // the picked residue's column, once the structure is lined up with the alignment: at the end of this row, which
        // has room for it, and moved to the next whole (never cut in two) where it has not
        let picked = style.picked.and_then(|k| loaded.fold.residues.get(k));
        if let Some((r, column)) = picked.and_then(|r| r.column.map(|k| (r, k))) {
            ui.add_space(10.0);
            ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
            to_alignment = ui
                .link(format!(
                    "{} {} is column {} of the alignment",
                    char::from(r.letter),
                    r.number,
                    column + 1
                ))
                .on_hover_text("Go to the Alignment page, where the column is selected")
                .clicked();
            ui.style_mut().wrap_mode = None;
        }
        // the colours are a row of their own: with the interface and the conservation there are too many for one
        ui.end_row();
        heading(ui, "Colour by");
        let only_predicted = (predicted, "Only a prediction says how sure it is");
        choice(
            ui,
            &mut style.scheme,
            Scheme::Confidence,
            "Confidence",
            only_predicted,
        );
        choice(ui, &mut style.scheme, Scheme::Secondary, "Structure", ALWAYS);
        choice(ui, &mut style.scheme, Scheme::Chain, "Chain", ALWAYS);
        choice(ui, &mut style.scheme, Scheme::Rainbow, "N to C", ALWAYS);
        let has_interface = loaded.fold.residues.iter().any(|r| r.buried.is_some());
        let after_analysing = (
            has_interface,
            "Use \"Colour by buried area\" under Interface in the options, after analysing two groups of chains",
        );
        choice(
            ui,
            &mut style.scheme,
            Scheme::Interface,
            "Interface",
            after_analysing,
        );
        let has_conservation = loaded.fold.residues.iter().any(|r| r.conservation.is_some());
        let after_colouring = (
            has_conservation,
            "Use \"Colour by conservation\" in the options: it needs an alignment of this protein's family",
        );
        choice(
            ui,
            &mut style.scheme,
            Scheme::Conservation,
            "Conservation",
            after_colouring,
        );
        ui.add_space(10.0);
        ui.checkbox(spin, "Spin");
    });
    to_alignment
}

fn viewer_panel(ui: &mut Ui, look: &Look, page: &mut StructurePage, c: &mut Clicks) {
    let Some(loaded) = &mut page.loaded else {
        let pick = design::empty_state(
            ui,
            look,
            "No structure yet",
            "Paste an amino-acid sequence and press Predict structure, or open a PDB file, or try the example. \
             The structure appears here, ready to turn round.",
            &["Open PDB…", "Try the example"],
        );
        match pick {
            Some(0) => c.open = true,
            Some(1) => c.example = true,
            _ => {}
        }
        return;
    };
    c.to_alignment = look_bar(ui, look, &mut page.style, &mut page.spin, loaded);
    ui.add_space(6.0);
    if page.style.representation == Representation::Sticks && loaded.bonds.is_none() {
        loaded.bonds = Some(loaded.fold.bonds());
    }
    let bonds = loaded.bonds.as_deref().unwrap_or(&[]);
    let mut touched = viewer3d::Touched::default();
    design::view_frame(ui, look, |ui| {
        let free = ui.available_height() - STRIP_HEIGHT - ui.spacing().item_spacing.y * 2.0;
        let room = egui::vec2(ui.available_width(), free.max(MIN_VIEWER_HEIGHT));
        if let Some(format) = crate::ui::figure::wanted() {
            // what the window shows: the same room, turned as it is, standing still
            let mut camera = page.camera;
            page.message = crate::ui::figure::save(format, "structure", room, |ui| {
                viewer3d::show(
                    ui,
                    look,
                    &loaded.fold,
                    &mut camera,
                    &page.style,
                    bonds,
                    false,
                );
            });
        }
        ui.allocate_ui(room, |ui| {
            (_, touched) = viewer3d::show(
                ui,
                look,
                &loaded.fold,
                &mut page.camera,
                &page.style,
                bonds,
                page.spin,
            );
        });
    });
    if touched.picked.is_some() {
        page.style.picked = touched.picked;
    } else if touched.cleared {
        page.style.picked = None;
    }
    sequence_strip(
        ui,
        look,
        &loaded.fold,
        page.style.scheme,
        &mut page.style.picked,
    );
}

/// The sequence in a row, each letter on the colour of its residue; click one to pick it in the 3D view.
fn sequence_strip(
    ui: &mut Ui,
    look: &Look,
    fold: &Fold,
    scheme: Scheme,
    picked: &mut Option<usize>,
) {
    design::view_frame(ui, look, |ui| {
        ScrollArea::horizontal()
            .auto_shrink([false, false])
            .max_height(STRIP_HEIGHT)
            .show(ui, |ui| {
                let n = fold.residues.len();
                let size = egui::vec2(n as f32 * LETTER_WIDTH + 16.0, STRIP_HEIGHT - 8.0);
                let (rect, response) = ui.allocate_exact_size(size, Sense::click());
                let p = ui.painter_at(rect);
                for (k, r) in fold.residues.iter().enumerate() {
                    let x = rect.left() + 8.0 + k as f32 * LETTER_WIDTH;
                    let cell = egui::Rect::from_min_size(
                        egui::pos2(x, rect.top() + 14.0),
                        egui::vec2(LETTER_WIDTH - 1.0, 22.0),
                    );
                    let colour = viewer3d::residue_colour(fold, k, scheme);
                    let tint = if look.dark {
                        colour.gamma_multiply(0.7)
                    } else {
                        crate::ui::theme::mix(colour, Color32::WHITE, 55.0)
                    };
                    p.rect_filled(cell, 2.0, tint);
                    if *picked == Some(k) {
                        p.rect_stroke(
                            cell.expand(1.0),
                            2.0,
                            egui::Stroke::new(2.0, look.accent),
                            egui::StrokeKind::Middle,
                        );
                    }
                    p.text(
                        cell.center(),
                        Align2::CENTER_CENTER,
                        char::from(r.letter).to_string(),
                        FontId::monospace(12.5),
                        look.text,
                    );
                    if (k + 1) % 10 == 0 {
                        p.text(
                            egui::pos2(cell.center().x, rect.top() + 6.0),
                            Align2::CENTER_CENTER,
                            r.number.to_string(),
                            FontId::proportional(9.5),
                            look.dim,
                        );
                    }
                }
                if response.clicked() {
                    if let Some(pointer) = response.interact_pointer_pos() {
                        let k = ((pointer.x - rect.left() - 8.0) / LETTER_WIDTH).floor();
                        *picked = (k >= 0.0 && (k as usize) < n).then_some(k as usize);
                    }
                }
            });
    });
}

mod interface;
mod kept;
mod link;

#[cfg(test)]
mod tests;
