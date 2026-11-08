//! Menu bar, toolbar, document tabs and the status line.

use super::app::{Action, App, Tab};
use super::clipboard::Special;
use super::icons::{self, Icon};
use super::platform::Workspace;
use super::state::{ThemeMode, TransMode};
use super::windows::enzlist::Quick;
use super::windows::orfs::{rule_label, ORF_MIN_CHOICES};
use crate::seq::copyspecial::Kind;
use crate::seq::dna::{StartRule, Strand};
use eframe::egui::{self, RichText, Ui};

/// Enzymes per submenu in "Insert Recognition Site" (the original's insert_enzyme_fill_cascade).
const INSERT_MENU_GROUP: usize = 15;
/// The sizes the original's Font Size menu offers (AppMain.tcl line 1382, 0x13a0e).
const FONT_SIZES: [u8; 15] = [8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 24, 28, 32, 36, 40];

impl App {
    pub(super) fn menu_bar(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let (has_sel, can_undo, can_redo, circular) = match self.current_doc() {
            Some(i) => {
                let dv = &self.docs[i];
                (
                    dv.selection().is_some(),
                    dv.can_undo(),
                    dv.can_redo(),
                    dv.doc.circular,
                )
            }
            None => (false, false, false, false),
        };
        egui::MenuBar::new().ui(ui, |ui| {
            let mut item = |ui: &mut Ui, label: &str, key: &str, enabled: bool, action: Action| {
                // the bundled font has no ⌘ or ⇧, so shortcuts are spelled out as the original does ("Cmd+S")
                let key = key
                    .replace("⇧⌘", "Cmd+Shift+")
                    .replace("⌥⌘", "Cmd+Alt+")
                    .replace('⌘', "Cmd+");
                let button = egui::Button::new(label).shortcut_text(key);
                if ui.add_enabled(enabled, button).clicked() {
                    out.push(action);
                    ui.close();
                }
            };
            ui.menu_button("File", |ui| {
                item(ui, "New", "⌘N", true, Action::New);
                item(ui, "Open…", "⌘O", true, Action::Open);
                item(
                    ui,
                    "Duplicate Selection",
                    "",
                    true,
                    Action::DuplicateSelection,
                );
                item(
                    ui,
                    "Duplicate Full Sequence",
                    "",
                    true,
                    Action::DuplicateFull,
                );
                ui.separator();
                item(ui, "Close", "⌘W", true, Action::CloseCurrent);
                item(ui, "Save", "⌘S", true, Action::Save);
                item(ui, "Save As…", "⇧⌘S", true, Action::SaveAs);
                item(ui, "Export as FASTA…", "", true, Action::ExportFasta);
            });
            ui.menu_button("Edit", |ui| {
                item(ui, "Undo", "⌘Z", can_undo, Action::Undo);
                item(ui, "Redo", "⇧⌘Z", can_redo, Action::Redo);
                ui.separator();
                item(ui, "Cut", "⌘X", has_sel, Action::Cut);
                item(ui, "Copy", "⌘C", has_sel, Action::Copy);
                item(ui, "Paste", "⌘V", true, Action::Paste);
                ui.separator();
                item(ui, "Cut Rev-Com", "", has_sel, Action::CutRevCom);
                item(ui, "Copy Rev-Com", "", has_sel, Action::CopyRevCom);
                item(ui, "Paste Rev-Com", "", true, Action::PasteRevCom);
                ui.separator();
                // the original's Copy Special menu, in its order (AppMain.tcl lines 1117-1139)
                ui.menu_button("Copy Special", |ui| {
                    let convert = |kind, reverse_complement| {
                        Action::CopySpecial(Special::Convert {
                            kind,
                            reverse_complement,
                        })
                    };
                    item(
                        ui,
                        "Copy All as GenBank",
                        "⌘B",
                        true,
                        Action::CopySpecial(Special::Genbank),
                    );
                    ui.separator();
                    item(
                        ui,
                        "Copy Uppercase",
                        "",
                        true,
                        convert(Kind::Uppercase, false),
                    );
                    item(
                        ui,
                        "Copy Uppercase Rev-Com",
                        "",
                        true,
                        convert(Kind::Uppercase, true),
                    );
                    ui.separator();
                    item(
                        ui,
                        "Copy Translated",
                        "",
                        true,
                        convert(Kind::Translated, false),
                    );
                    item(
                        ui,
                        "Copy Uppercase Translated",
                        "",
                        true,
                        convert(Kind::TranslatedUppercase, false),
                    );
                    item(
                        ui,
                        "Copy Translated Rev-Com",
                        "",
                        true,
                        convert(Kind::Translated, true),
                    );
                    item(
                        ui,
                        "Copy Uppercase Translated Rev-Com",
                        "",
                        true,
                        convert(Kind::TranslatedUppercase, true),
                    );
                    ui.separator();
                    item(
                        ui,
                        "Copy as FASTA",
                        "",
                        true,
                        Action::CopySpecial(Special::Fasta),
                    );
                    ui.separator();
                    item(ui, "Copy as RNA", "", true, convert(Kind::Rna, false));
                    ui.separator();
                    item(
                        ui,
                        "Copy Reverse (NOT A GOOD IDEA)",
                        "",
                        true,
                        convert(Kind::Reverse, false),
                    );
                    item(
                        ui,
                        "Copy Complement (NOT A GOOD IDEA)",
                        "",
                        true,
                        convert(Kind::Complement, false),
                    );
                    ui.separator();
                    item(
                        ui,
                        "Copy Features as NCBI Bankit table",
                        "",
                        true,
                        Action::CopySpecial(Special::Bankit),
                    );
                });
                ui.separator();
                item(ui, "Select All", "⌘A", true, Action::SelectAll);
                item(ui, "Select From-To…", "⇧⌘A", true, Action::SelectRange);
                item(ui, "Jump To…", "⌘J", true, Action::Jump);
                ui.separator();
                item(ui, "Find…", "⌘F", true, Action::Find);
                item(ui, "Find Again", "⌘G", true, Action::FindAgain);
                item(
                    ui,
                    "Clear Find Highlighting",
                    "⇧⌘F",
                    true,
                    Action::ClearFind,
                );
                ui.separator();
                item(ui, "Convert to UPPERCASE", "", true, Action::Upper);
                item(ui, "Convert to lowercase", "", true, Action::Lower);
                item(ui, "UPPER<->lower", "", true, Action::SwapCase);
                ui.separator();
                item(ui, "Reverse-Complement", "⌘/", true, Action::RevCom);
                ui.separator();
                item(ui, "Circular<->Linear", "", true, Action::ToggleCircular);
                item(ui, "Set Origin", "⌘I", circular, Action::SetOrigin);
                item(
                    ui,
                    "Linearize @ insert site",
                    "",
                    circular,
                    Action::Linearize,
                );
            });
            ui.menu_button("Enzymes", |ui| {
                item(ui, "Enzyme Selector…", "⌘E", true, Action::EnzymeSelector);
                ui.separator();
                // the original's check item (L1170); the state is in its words, as an item has no check mark here
                let state = if self.prefs.selection_only {
                    "on"
                } else {
                    "off"
                };
                let label = format!("Selection Only: {state}");
                item(ui, &label, "⇧⌘E", true, Action::ToggleSelectionOnly);
                ui.separator();
                item(ui, "Digestion", "⌘R", true, Action::Digest);
                item(ui, "Highlight Sites", "", true, Action::ToggleHighlight);
                ui.separator();
                item(ui, "List Enzymes…", "", true, Action::ListEnzymes);
                item(ui, "Text Map…", "⇧⌘T", true, Action::TextMap);
                ui.menu_button("Quick Lists", |ui| {
                    for quick in Quick::ALL {
                        item(ui, quick.label(), "", true, Action::QuickList(quick));
                    }
                });
                ui.separator();
                ui.menu_button("Insert Recognition Site", |ui| {
                    for (g, group) in self.enzymes.chunks(INSERT_MENU_GROUP).enumerate() {
                        let (Some(first), Some(last)) = (group.first(), group.last()) else {
                            continue;
                        };
                        ui.menu_button(format!("{}-{}", first.name, last.name), |ui| {
                            for (k, e) in group.iter().enumerate() {
                                let enzyme = g * INSERT_MENU_GROUP + k;
                                let site = |reverse| Action::InsertSite { enzyme, reverse };
                                if e.is_palindromic() {
                                    item(ui, &e.name, "", true, site(false));
                                } else {
                                    item(
                                        ui,
                                        &format!("{} (Forward)", e.name),
                                        "",
                                        true,
                                        site(false),
                                    );
                                    item(
                                        ui,
                                        &format!("{} (Reverse)", e.name),
                                        "",
                                        true,
                                        site(true),
                                    );
                                }
                            }
                        });
                    }
                });
                ui.separator();
                item(ui, "Silent Sites…", "", true, Action::SilentSites);
                item(
                    ui,
                    "Find Diagnostic Site…",
                    "",
                    true,
                    Action::DiagnosticSites,
                );
            });
            ui.menu_button("ORFs", |ui| {
                item(ui, "Find Next", "⇧⌘.", true, Action::OrfNext);
                item(ui, "Find Previous", "⇧⌘,", true, Action::OrfPrev);
                ui.separator();
                ui.menu_button("ORF Starts With", |ui| {
                    for r in [StartRule::Met, StartRule::AfterStop] {
                        ui.radio_value(&mut self.prefs.orf_rule, r, rule_label(r));
                    }
                });
                ui.menu_button("Search Strand", |ui| {
                    for (s, name) in [
                        (Strand::Top, "Top"),
                        (Strand::Bottom, "Bottom"),
                        (Strand::Both, "Both"),
                    ] {
                        ui.radio_value(&mut self.prefs.orf_strand, s, name);
                    }
                });
                ui.menu_button("Minimum bp", |ui| {
                    for m in ORF_MIN_CHOICES {
                        ui.radio_value(&mut self.prefs.orf_min, m, m.to_string());
                    }
                });
                ui.separator();
                item(ui, "Translate…", "⌘T", true, Action::Translate);
                item(ui, "ORF Map and List…", "", true, Action::OrfList);
                ui.separator();
                item(
                    ui,
                    "Selection Translation: 3 frames",
                    "",
                    true,
                    Action::Translation(TransMode::Forward),
                );
                item(
                    ui,
                    "Selection Translation: 6 frames",
                    "",
                    true,
                    Action::Translation(TransMode::Six),
                );
            });
            ui.menu_button("Features", |ui| {
                item(ui, "New Feature…", "⌘.", true, Action::NewFeature);
                item(
                    ui,
                    "New Features From Found Highlighting…",
                    "",
                    true,
                    Action::FeaturesFromFound,
                );
                item(ui, "Edit Selected Feature…", "", true, Action::EditFeature);
                item(
                    ui,
                    "Delete Selected Feature",
                    "",
                    true,
                    Action::DeleteFeature,
                );
                item(ui, "Clear Features", "", true, Action::ClearFeatures);
                ui.separator();
                item(
                    ui,
                    "Annotate Features Using Library",
                    "⇧⌘K",
                    true,
                    Action::Annotate,
                );
                item(
                    ui,
                    "Open Feature Library and Annotate…",
                    "",
                    true,
                    Action::OpenLibrary,
                );
                ui.separator();
                item(ui, "Show Feature Table", "", true, Action::ToggleFeatures);
            });
            ui.menu_button("Tools", |ui| {
                item(ui, "Align Sequences…", "⌘L", true, Action::Align);
                ui.separator();
                item(ui, "Contig Assembler…", "", true, Action::OpenAssembler);
                item(
                    ui,
                    "Gibson-style Assembly…",
                    "",
                    true,
                    Action::GibsonAssembler,
                );
                ui.separator();
                item(ui, "PCR Products and Primers…", "", true, Action::Pcr);
                item(ui, "Primer3 Pair Design…", "", true, Action::Primer3);
                item(ui, "Read Quality and Trimming…", "", true, Action::Reads);
                item(ui, "Calculator…", "⌘1", true, Action::Calculator);
                item(
                    ui,
                    "Molecular Reaction Calculator…",
                    "⌘2",
                    true,
                    Action::MolecularCalculator,
                );
                item(ui, "Tm and GC Calculator…", "", true, Action::TmCalculator);
                ui.separator();
                item(ui, "Phylogenetic Tree…", "", true, Action::Tree);
                item(
                    ui,
                    "PhyloCNN Epidemic Estimates…",
                    "",
                    true,
                    Action::PhyloCnn,
                );
                item(
                    ui,
                    "Species Tree from Gene Trees…",
                    "",
                    true,
                    Action::SpeciesTree,
                );
                item(ui, "Tree Tools (Newick)…", "", true, Action::TreeTools);
                item(
                    ui,
                    "Pipelines (MAFFT, IQ-TREE, primer3)",
                    "⌥⌘8",
                    true,
                    Action::Workflow,
                );
            });
            ui.menu_button("Assembly", |ui| {
                item(ui, "Contig Assembler…", "", true, Action::OpenAssembler);
                item(
                    ui,
                    "Add Open Sequences to Assembly",
                    "",
                    true,
                    Action::AddOpenToAssembly,
                );
                item(ui, "Assemble", "", true, Action::AssembleNow);
                ui.separator();
                item(
                    ui,
                    "Gibson-style Assembly…",
                    "",
                    true,
                    Action::GibsonAssembler,
                );
                item(ui, "Load Demo Reads", "", true, Action::DemoReads);
            });
            ui.menu_button("Go", |ui| {
                for (k, ws) in Workspace::ALL.into_iter().enumerate() {
                    let key = format!("⌥⌘{}", k + 1);
                    item(ui, ws.label(), &key, true, Action::Go(ws));
                }
            });
            ui.menu_button("View", |ui| {
                item(ui, "Show Map", "", true, Action::ToggleMap);
                item(ui, "Show Feature Table", "", true, Action::ToggleFeatures);
                item(
                    ui,
                    "Show Complement Strand",
                    "",
                    true,
                    Action::ToggleSecondStrand,
                );
                ui.separator();
                item(
                    ui,
                    "Translation: 3 frames",
                    "",
                    true,
                    Action::Translation(TransMode::Forward),
                );
                item(
                    ui,
                    "Translation: 6 frames",
                    "",
                    true,
                    Action::Translation(TransMode::Six),
                );
                ui.separator();
                ui.menu_button("Font Size", |ui| {
                    for size in FONT_SIZES {
                        let value = f32::from(size);
                        ui.radio_value(&mut self.prefs.font_size, value, size.to_string());
                    }
                });
                ui.menu_button("Theme", |ui| {
                    for mode in ThemeMode::ALL {
                        ui.radio_value(&mut self.prefs.theme, mode, mode.label());
                    }
                });
            });
            ui.menu_button("Help", |ui| {
                item(ui, "Command Palette…", "⌘K", true, Action::Palette);
                item(ui, "Keyboard Shortcuts", "", true, Action::Shortcuts);
                ui.separator();
                item(ui, "About HelixLoom", "", true, Action::About);
            });
        });
    }

    pub(super) fn toolbar(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        let pal = self.palette();
        let selection_only = self.prefs.selection_only;
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 2.0;
            let shift = ui.input(|i| i.modifiers.shift);
            let mut b = |ui: &mut Ui, icon: Icon, tip: &str, action: Action| {
                if icons::button(ui, icon, tip, &pal).clicked() {
                    out.push(action);
                }
            };
            b(ui, Icon::New, "New", Action::New);
            b(ui, Icon::Open, "Open…", Action::Open);
            b(ui, Icon::Save, "Save", Action::Save);
            ui.separator();
            b(ui, Icon::Undo, "Undo", Action::Undo);
            b(ui, Icon::Redo, "Redo", Action::Redo);
            ui.separator();
            b(ui, Icon::Cut, "Cut", Action::Cut);
            b(ui, Icon::Copy, "Copy", Action::Copy);
            b(ui, Icon::Paste, "Paste", Action::Paste);
            ui.separator();
            b(ui, Icon::Find, "Find…", Action::Find);
            ui.separator();
            // the original's enzyme button opens the selector, and with Shift turns Selection Only on or off (L4370); its
            // picture shows which
            let (icon, action) = enzyme_button(selection_only, shift);
            b(
                ui,
                icon,
                "Enzyme Selector…   Shift-click: Selection Only on or off",
                action,
            );
            b(
                ui,
                Icon::Map,
                "Show / hide the graphic map",
                Action::ToggleMap,
            );
            b(ui, Icon::Digest, "Digestion", Action::Digest);
            ui.separator();
            b(
                ui,
                Icon::Translate,
                "Show translation (3 frames)",
                Action::Translation(TransMode::Forward),
            );
            b(ui, Icon::Orf, "ORF map and list", Action::OrfList);
            b(
                ui,
                Icon::Features,
                "Show / hide the feature table",
                Action::ToggleFeatures,
            );
            ui.separator();
            b(ui, Icon::Primers, "PCR products and primers", Action::Pcr);
            b(ui, Icon::Align, "Align sequences", Action::Align);
            ui.separator();
            b(
                ui,
                Icon::Assemble,
                "Contig Assembler: join reads or fragments into contigs",
                Action::OpenAssembler,
            );
        });
    }

    pub(super) fn tab_strip(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        ui.horizontal_wrapped(|ui| {
            let current = self.current_doc();
            let on_assembler = self.tab == Tab::Assembly;
            for (i, dv) in self.docs.iter().enumerate() {
                let selected = !on_assembler && current == Some(i);
                if ui.selectable_label(selected, dv.doc.title()).clicked() {
                    out.push(Action::SwitchTab(Tab::Doc(i)));
                }
                if ui.small_button("×").on_hover_text("Close").clicked() {
                    out.push(Action::CloseTab(i));
                }
                ui.add_space(6.0);
            }
            ui.separator();
            if ui
                .selectable_label(on_assembler, "Contig Assembler")
                .clicked()
            {
                out.push(Action::SwitchTab(Tab::Assembly));
            }
        });
    }

    pub(super) fn status_bar(&self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.label(RichText::new(&self.status).small());
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(
                    RichText::new(format!("{} v{}", "HelixLoom", env!("CARGO_PKG_VERSION")))
                        .small()
                        .weak(),
                );
            });
        });
    }
}

/// What the toolbar's enzyme button shows and does: the selector, or with Shift the Selection Only switch; its picture says
/// whether the switch is on.
fn enzyme_button(selection_only: bool, shift: bool) -> (Icon, Action) {
    let icon = if selection_only {
        Icon::EnzymesSelectionOnly
    } else {
        Icon::Enzymes
    };
    let action = if shift {
        Action::ToggleSelectionOnly
    } else {
        Action::EnzymeSelector
    };
    (icon, action)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_enzyme_button_opens_the_selector_and_with_shift_turns_the_switch_over() {
        assert!(matches!(
            enzyme_button(false, false),
            (Icon::Enzymes, Action::EnzymeSelector)
        ));
        assert!(matches!(
            enzyme_button(false, true),
            (Icon::Enzymes, Action::ToggleSelectionOnly)
        ));
    }

    #[test]
    fn the_picture_of_the_enzyme_button_shows_whether_selection_only_is_on() {
        assert!(matches!(
            enzyme_button(true, false),
            (Icon::EnzymesSelectionOnly, Action::EnzymeSelector)
        ));
        assert!(matches!(
            enzyme_button(true, true),
            (Icon::EnzymesSelectionOnly, Action::ToggleSelectionOnly)
        ));
    }
}
