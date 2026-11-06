//! Structure: predict the 3D structure of a protein from its amino-acid sequence with ESMFold (through the ESM
//! Metagenomic Atlas web service) or open a PDB file, and turn it round in a 3D viewer coloured by confidence,
//! secondary structure, chain or position along the chain.

use super::route::{self, Source};
use super::viewer3d::{self, Camera, Representation, Scheme, Style};
use super::workstation::{compute_picker, ComputeView, Place};
use crate::pipeline::predict::{self, ESM_ATLAS_URL, MAX_LENGTH};
use crate::pipeline::steps::fasta_text;
use crate::pipeline::tools::Installed;
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
use eframe::egui::{self, Align2, Color32, FontId, RichText, ScrollArea, Sense, Ui};
use std::path::Path;
use std::time::Duration;
use web_time::Instant;

const OPTIONS_WIDTH: f32 = 310.0;
const STRIP_HEIGHT: f32 = 62.0;
const LETTER_WIDTH: f32 = 15.0;
const WORKING_REDRAW: Duration = Duration::from_millis(250);
/// The ORFs the open sequence is searched for must be at least this many bases long.
const MIN_ORF_BASES: usize = 90;
/// The viewer never gets less room than this, however small the window.
const MIN_VIEWER_HEIGHT: f32 = 160.0;
/// The example structure: ubiquitin, PDB entry 1UBQ.
const EXAMPLE_PDB: &str = include_str!("../../../assets/structure/1ubq.pdb");

pub struct Loaded {
    pub name: String,
    /// The file as it was read, so saving gives back exactly what came in.
    text: String,
    pub fold: Fold,
    /// The bonds of every atom, found when the all-bonds drawing is first asked for.
    bonds: Option<Vec<(usize, usize)>>,
}

#[derive(Default)]
pub struct StructurePage {
    pub loaded: Option<Loaded>,
    sequence: String,
    style: Style,
    camera: Camera,
    spin: bool,
    job: Option<(Job<Result<String, String>>, Instant)>,
    pub message: String,
}

impl StructurePage {
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

    /// What is written in the sequence box, checked.
    fn checked(&self) -> Result<Vec<u8>, String> {
        predict::clean_sequence(&self.sequence)
    }
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
    save: bool,
    reset: bool,
}

impl App {
    /// Collect a finished prediction.
    pub(super) fn structure_tick(&mut self, ctx: &egui::Context) {
        let page = &mut self.platform.structure;
        let Some(polled) = page.job.as_ref().map(|(job, _)| job.poll()) else {
            return;
        };
        match polled {
            Poll::Running => ctx.request_repaint_after(WORKING_REDRAW),
            Poll::Done(Ok(text)) => {
                let seconds = page
                    .job
                    .take()
                    .map_or(0.0, |(_, started)| started.elapsed().as_secs_f64());
                match page.show_pdb("ESMFold prediction", &text) {
                    Ok(()) => {
                        page.message = format!("Predicted in {seconds:.1} s. {}", page.message);
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
        let look = design::look(self.prefs.dark);
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
        let page = &mut self.platform.structure;
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
                    });
                },
            );
            egui::Panel::left("structure_options")
                .exact_size(OPTIONS_WIDTH)
                .resizable(false)
                .frame(egui::Frame::new())
                .show(ui, |ui| {
                    ScrollArea::vertical()
                        .auto_shrink([false, false])
                        .show(ui, |ui| {
                            design::card(ui, &look, |ui| {
                                let curl = Curl {
                                    found: curl.is_some() || via == Some(true),
                                    searched,
                                    via,
                                    view,
                                };
                                options_panel(ui, &look, page, curl, &mut c);
                            });
                        });
                });
            egui::CentralPanel::default()
                .frame(egui::Frame::new().inner_margin(egui::Margin {
                    left: design::PAGE_MARGIN,
                    ..egui::Margin::ZERO
                }))
                .show(ui, |ui| viewer_panel(ui, &look, page));
        });
        self.apply_structure_clicks(c, curl, ui.ctx());
    }

    fn apply_structure_clicks(&mut self, c: Clicks, curl: Option<Installed>, ctx: &egui::Context) {
        if c.from_document {
            self.protein_from_open_document();
        }
        let place = self.platform.place(&self.fold_weight());
        let gpu_fold = self.platform.workstation.folds_on_gpu();
        let page = &mut self.platform.structure;
        if c.predict {
            match (page.checked(), place, curl) {
                (Err(why), _, _) => page.message = why,
                (Ok(sequence), Place::There(remote), _) => {
                    page.message = "Asking the workstation to fold it…".into();
                    let label = format!("Fold of {} residues", sequence.len());
                    let request = Request::Fold {
                        sequence,
                        local: gpu_fold,
                    };
                    page.job = Some((
                        Job::spawn(ctx, move || {
                            let folded = remote.run(&request, &label)?.into_fold()?;
                            Ok(folded.pdb)
                        }),
                        Instant::now(),
                    ));
                }
                (Ok(sequence), Place::Here, Some(curl)) => {
                    page.message = "Asking the ESM Atlas service to fold it…".into();
                    let path = curl.path;
                    page.job = Some((
                        Job::spawn(ctx, move || {
                            predict::predict(&path, ESM_ATLAS_URL, &sequence)
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

fn options_panel(ui: &mut Ui, look: &Look, page: &mut StructurePage, curl: Curl, c: &mut Clicks) {
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
    let ready = checked.is_ok() && curl.found && !page.is_running();
    compute_picker(ui, look, curl.view);
    ui.horizontal(|ui| {
        c.predict = ui
            .add_enabled(
                ready,
                egui::Button::new(RichText::new("Predict structure").strong()),
            )
            .clicked();
        if page.is_running() {
            ui.spinner();
        }
    });
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

    design::section(ui, look, "Or open a structure");
    ui.horizontal_wrapped(|ui| {
        c.open = ui.button("Open PDB…").clicked();
        c.example = ui.button("Example").clicked();
    });

    if let Some(loaded) = &page.loaded {
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
    }
}

/// The choices of how the structure is drawn and coloured, in a row above it.
fn look_bar(ui: &mut Ui, look: &Look, style: &mut Style, spin: &mut bool, loaded: &Loaded) {
    let (atoms, predicted) = (loaded.fold.structure.atoms.len(), loaded.fold.predicted);
    let heading =
        |ui: &mut Ui, text: &str| ui.label(RichText::new(text).size(12.0).color(look.dim));
    ui.horizontal_wrapped(|ui| {
        heading(ui, "Show");
        ui.selectable_value(&mut style.representation, Representation::Ribbon, "Ribbon");
        ui.selectable_value(
            &mut style.representation,
            Representation::Backbone,
            "Backbone",
        );
        ui.add_enabled_ui(atoms <= viewer3d::MAX_STICK_ATOMS, |ui| {
            ui.selectable_value(
                &mut style.representation,
                Representation::Sticks,
                "All bonds",
            )
            .on_disabled_hover_text("Too many atoms to draw every bond")
        });
        ui.add_space(10.0);
        heading(ui, "Colour by");
        ui.add_enabled_ui(predicted, |ui| {
            ui.selectable_value(&mut style.scheme, Scheme::Confidence, "Confidence")
                .on_disabled_hover_text("Only a prediction says how sure it is")
        });
        ui.selectable_value(&mut style.scheme, Scheme::Secondary, "Structure");
        ui.selectable_value(&mut style.scheme, Scheme::Chain, "Chain");
        ui.selectable_value(&mut style.scheme, Scheme::Rainbow, "N to C");
        ui.add_space(10.0);
        ui.checkbox(spin, "Spin");
    });
}

fn viewer_panel(ui: &mut Ui, look: &Look, page: &mut StructurePage) {
    let Some(loaded) = &mut page.loaded else {
        design::card(ui, look, |ui| {
            ui.label(
                RichText::new("No structure yet")
                    .size(17.0)
                    .strong()
                    .color(look.text),
            );
            ui.label(
                RichText::new(
                    "Paste an amino-acid sequence and press Predict structure, or open a PDB file, or try the example. \
                     The structure appears here, ready to turn round.",
                )
                .color(look.dim),
            );
        });
        return;
    };
    look_bar(ui, look, &mut page.style, &mut page.spin, loaded);
    ui.add_space(6.0);
    if page.style.representation == Representation::Sticks && loaded.bonds.is_none() {
        loaded.bonds = Some(loaded.fold.bonds());
    }
    let bonds = loaded.bonds.as_deref().unwrap_or(&[]);
    let mut touched = viewer3d::Touched::default();
    egui::Frame::new()
        .fill(look.surface)
        .stroke(egui::Stroke::new(1.0, look.border))
        .corner_radius(egui::CornerRadius::same(design::CARD_RADIUS))
        .show(ui, |ui| {
            let free = ui.available_height() - STRIP_HEIGHT - ui.spacing().item_spacing.y * 2.0;
            ui.allocate_ui(
                egui::vec2(ui.available_width(), free.max(MIN_VIEWER_HEIGHT)),
                |ui| {
                    (_, touched) = viewer3d::show(
                        ui,
                        look,
                        &loaded.fold,
                        &mut page.camera,
                        &page.style,
                        bonds,
                        page.spin,
                    );
                },
            );
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
    egui::Frame::new()
        .fill(look.surface)
        .stroke(egui::Stroke::new(1.0, look.border))
        .corner_radius(egui::CornerRadius::same(design::CARD_RADIUS))
        .show(ui, |ui| {
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

#[cfg(test)]
mod tests;
