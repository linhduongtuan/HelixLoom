//! Opening, saving and closing sequences.

use super::app::{App, Tab};
use super::dialog::{self, FileDialog, Wanted};
use super::platform::route::{self, Source};
use super::state::DocView;
use crate::seq::doc::Doc;
use crate::seq::io;
use eframe::egui::{self};
use std::path::Path;

impl App {
    // -----------------------------------------------------------------------
    // Files
    // -----------------------------------------------------------------------

    /// Open the sequences in a file, each in a tab of its own.
    pub(super) fn load_sequences(&mut self, file: Source) {
        let shown = file.shown();
        let bytes = match file.bytes() {
            Ok(b) => b,
            Err(e) => return self.say(format!("Couldn't open {shown}: {e}")),
        };
        let stem = file
            .path()
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("Untitled")
            .to_string();
        match io::parse_any(&String::from_utf8_lossy(&bytes), &stem) {
            Ok(docs) => {
                let count = docs.len();
                for mut d in docs {
                    d.path = file.on_disk().map(Path::to_path_buf);
                    d.dirty = false;
                    self.docs.push(DocView::new(d));
                }
                self.tab = Tab::Doc(self.docs.len() - 1);
                self.say(format!(
                    "Opened {} sequence{} from {shown}",
                    count,
                    if count == 1 { "" } else { "s" },
                ));
            }
            Err(e) => self.say(format!("{shown}: {e}")),
        }
    }

    /// Open pasted GenBank or FASTA text as new, unsaved documents.
    pub fn open_text(&mut self, text: &str, name: &str) {
        match io::parse_any(text, name) {
            Ok(docs) => {
                let count = docs.len();
                for mut d in docs {
                    d.dirty = true;
                    self.docs.push(DocView::new(d));
                }
                self.tab = Tab::Doc(self.docs.len() - 1);
                self.say(format!(
                    "Opened {count} sequence{} from the clipboard",
                    if count == 1 { "" } else { "s" }
                ));
            }
            Err(e) => self.say(format!("The clipboard text isn't a sequence record: {e}")),
        }
    }

    /// The Open dialog: sequences, and the other files the workspaces read, which go where they belong.
    pub(super) fn pick_open(&mut self) {
        let ask = FileDialog::new()
            .set_title("Open")
            .add_filter(
                "Sequences",
                &[
                    "gb", "gbk", "ape", "gen", "fasta", "fa", "fas", "seq", "txt", "dna",
                ],
            )
            .add_filter("Sanger traces", &route::TRACE_EXTENSIONS)
            .add_filter("Protein structures", &route::STRUCTURE_EXTENSIONS)
            .add_filter("Trees", &route::TREE_EXTENSIONS);
        for f in dialog::pick_many(Wanted::Anything, ask) {
            self.open_file(&f);
        }
    }

    /// In a browser a sequence is saved by downloading it: as FASTA, or else as GenBank. There is no file to save over.
    #[cfg(target_arch = "wasm32")]
    pub(super) fn save_doc(&mut self, i: usize, _ask: bool, fasta: bool) {
        let doc = &self.docs[i].doc;
        let (ext, text) = if fasta {
            ("fasta", io::to_fasta(doc, &doc.name))
        } else {
            ("gb", io::to_genbank(doc))
        };
        let name = format!("{}.{ext}", doc.name);
        let Some(saved) = dialog::save_bytes(FileDialog::new(), &name, text.as_bytes()) else {
            return;
        };
        if saved.ok && !fasta {
            self.docs[i].doc.dirty = false;
        }
        self.say(saved.message);
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub(super) fn save_doc(&mut self, i: usize, ask: bool, fasta: bool) {
        let doc = &self.docs[i].doc;
        let ext = if fasta { "fasta" } else { "gb" };
        let existing = doc.path.clone().filter(|_| !ask && !fasta);
        let path = existing.or_else(|| {
            FileDialog::new()
                .set_title(if fasta {
                    "Export FASTA"
                } else {
                    "Save sequence"
                })
                .set_file_name(format!("{}.{ext}", doc.name))
                .save_file()
        });
        let Some(path) = path else { return };
        let is_fasta = matches!(
            path.extension().and_then(|e| e.to_str()),
            Some("fasta" | "fa" | "fas")
        );
        let text = if is_fasta {
            let stem = path.file_stem().and_then(|s| s.to_str());
            io::to_fasta(doc, stem.unwrap_or(&doc.name))
        } else {
            io::to_genbank(doc)
        };
        match std::fs::write(&path, text) {
            Ok(()) => {
                let dv = &mut self.docs[i];
                if !fasta {
                    dv.doc.path = Some(path.clone());
                    dv.doc.dirty = false;
                    if let Some(stem) = path.file_stem().and_then(|s| s.to_str()) {
                        dv.doc.name = stem.to_string();
                    }
                }
                self.say(format!("Saved {}", path.display()));
            }
            Err(e) => self.say(format!("Couldn't save {}: {e}", path.display())),
        }
    }

    pub fn add_doc(&mut self, doc: Doc) {
        self.docs.push(DocView::new(doc));
        self.tab = Tab::Doc(self.docs.len() - 1);
    }

    pub(super) fn close_tab(&mut self, i: usize) {
        if i >= self.docs.len() {
            return;
        }
        self.docs.remove(i);
        if self.docs.is_empty() {
            self.docs
                .push(DocView::new(Doc::new("Untitled", Vec::new())));
        }
        self.tab = match self.tab {
            Tab::Doc(t) if t > i => Tab::Doc(t - 1),
            Tab::Doc(t) => Tab::Doc(t.min(self.docs.len() - 1)),
            Tab::Assembly => Tab::Assembly,
        };
    }

    pub(super) fn request_close(&mut self, i: usize) {
        if self.docs.get(i).is_some_and(|d| d.doc.dirty) {
            self.confirm_close = Some(i);
        } else {
            self.close_tab(i);
        }
    }

    pub(super) fn close_prompt(&mut self, ctx: &egui::Context) {
        let Some(i) = self.confirm_close else { return };
        let name = self
            .docs
            .get(i)
            .map(|d| d.doc.name.clone())
            .unwrap_or_default();
        let mut decision: Option<bool> = None;
        crate::ui::design::window("Unsaved changes")
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.label(format!("“{name}” has changes that haven't been saved."));
                ui.horizontal(|ui| {
                    if ui.button("Save").clicked() {
                        self.save_doc(i, false, false);
                        decision = Some(!self.docs[i].doc.dirty);
                    }
                    if ui.button("Discard").clicked() {
                        decision = Some(true);
                    }
                    if ui.button("Cancel").clicked() {
                        decision = Some(false);
                    }
                });
            });
        if let Some(close) = decision {
            self.confirm_close = None;
            if close {
                self.close_tab(i);
            }
        }
    }
}
