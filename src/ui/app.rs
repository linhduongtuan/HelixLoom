//! The application shell: the top bar and its menus, the rail, toolbar, tabs, the sequence view with its map and feature
//! table, the command palette, and the screenshot hook used by the build checks.

use super::clipboard::{self, Special};
use super::design;
use super::headless::{claim_picture, save_png, SCREENSHOT_FRAMES};
use super::palette::CommandPalette;
use super::platform::{Platform, Workspace, RAIL_WIDTH};
use super::state::{DocView, Prefs, ThemeMode, TransMode};
use super::theme::{self, Palette};
use super::toasts::Toasts;
use super::windows::assembly::{AssemblyWin, Preset};
use super::windows::calculator::CalculatorWin;
use super::windows::digest::DigestWin;
use super::windows::enzlist::{EnzListWin, Quick};
use super::windows::enzymes::EnzymeWin;
use super::windows::features::FeatureWin;
use super::windows::find::{Dir, FindWin};
use super::windows::molcalc::MolCalcWin;
use super::windows::orfs::{OrfStep, OrfWin};
use super::windows::phylocnn::PhyloCnnWin;
use super::windows::primer3::Primer3Win;
use super::windows::reads::ReadsWin;
use super::windows::sitetools::SiteToolsWin;
use super::windows::species::SpeciesWin;
use super::windows::tools::{AlignWin, CalcWin, PcrWin};
use super::windows::tree::TreeWin;
use super::windows::treetools::TreeToolsWin;
use crate::seq::dna::{revcom, swap_case, Strand, GENETIC_CODES};
use crate::seq::doc::Doc;
use crate::seq::enzymes::{default_enzymes, Enzyme};
use crate::seq::estimate::Analysis;
use crate::seq::features::{builtin_library, LibraryEntry};
use eframe::egui::{self, Key, KeyboardShortcut, Modifiers, Ui};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Tab {
    Doc(usize),
    Assembly,
}

/// What the person asked the program to do at startup (command line).
#[derive(Default, Clone)]
pub struct Startup {
    pub screenshot: Option<String>,
    pub view: String,
    pub files: Vec<String>,
    /// Draw the screenshot without opening a window at all.
    pub no_window: bool,
    /// Start in the dark theme.
    pub dark: bool,
    /// The folder the program keeps its settings in (the workstation's address and token); none reads and keeps
    /// nothing.
    pub settings_dir: Option<PathBuf>,
}

#[derive(Clone, Debug)]
pub enum Action {
    /// Go to a workspace, and go to a workspace with its example loaded.
    Go(Workspace),
    Example(Workspace),
    New,
    Open,
    DuplicateSelection,
    DuplicateFull,
    Save,
    SaveAs,
    ExportFasta,
    CloseTab(usize),
    CloseCurrent,
    Undo,
    Redo,
    Cut,
    Copy,
    Paste,
    CutRevCom,
    CopyRevCom,
    PasteRevCom,
    CopySpecial(Special),
    SelectAll,
    Find,
    FindAgain,
    Jump,
    SelectRange,
    ClearFind,
    Upper,
    Lower,
    SwapCase,
    RevCom,
    ToggleCircular,
    SetOrigin,
    Linearize,
    EnzymeSelector,
    Digest,
    /// List Enzymes… (its dialog) and the three Quick Lists of the Enzymes menu.
    ListEnzymes,
    QuickList(Quick),
    SilentSites,
    DiagnosticSites,
    ToggleHighlight,
    /// Put an enzyme's recognition sequence at the caret (index into the enzyme list).
    InsertSite {
        enzyme: usize,
        reverse: bool,
    },
    Translation(TransMode),
    ToggleMap,
    ToggleFeatures,
    ToggleSecondStrand,
    SwitchTab(Tab),
    OpenAssembler,
    GibsonAssembler,
    AssembleNow,
    AddOpenToAssembly,
    DemoReads,
    OrfNext,
    OrfPrev,
    OrfList,
    Translate,
    NewFeature,
    FeaturesFromFound,
    EditFeature,
    DeleteFeature,
    ClearFeatures,
    Annotate,
    OpenLibrary,
    Align,
    Pcr,
    /// The arithmetic Calculator of the original.
    Calculator,
    TmCalculator,
    MolecularCalculator,
    Tree,
    PhyloCnn,
    SpeciesTree,
    TreeTools,
    Workflow,
    Primer3,
    Reads,
    About,
    /// Open or close the command palette.
    Palette,
    /// Show the sheet of keyboard shortcuts.
    Shortcuts,
    /// Choose how the colours are decided: by the system, light or dark.
    Theme(ThemeMode),
}

pub struct App {
    pub docs: Vec<DocView>,
    pub tab: Tab,
    pub(super) last_doc: usize,
    pub enzymes: Vec<Enzyme>,
    pub selected: BTreeSet<String>,
    pub plus_unique: bool,
    pub library: Vec<LibraryEntry>,
    pub prefs: Prefs,
    pub status: String,
    pub map_linear: bool,
    pub enzyme_win: EnzymeWin,
    pub digest_win: DigestWin,
    pub enz_list: EnzListWin,
    pub find_win: FindWin,
    pub assembly: AssemblyWin,
    pub orf_win: OrfWin,
    pub site_win: SiteToolsWin,
    pub feature_win: FeatureWin,
    pub align_win: AlignWin,
    pub pcr_win: PcrWin,
    pub calc_win: CalcWin,
    pub calculator_win: CalculatorWin,
    pub molcalc_win: MolCalcWin,
    pub tree_win: TreeWin,
    pub phylocnn_win: PhyloCnnWin,
    pub species_win: SpeciesWin,
    pub treetools_win: TreeToolsWin,
    /// The workspace showing beside the rail, and what the platform keeps for them.
    pub workspace: Workspace,
    pub platform: Platform,
    pub p3_win: Primer3Win,
    pub reads_win: ReadsWin,
    pub about_open: bool,
    pub confirm_close: Option<usize>,
    pub palette: CommandPalette,
    pub toasts: Toasts,
    pub shortcuts_open: bool,
    paste_rc_pending: bool,
    shot: Option<String>,
    shot_asked: bool,
    frames: u32,
    /// The colours last put on the window: whether dark, and whether they follow the system.
    theme_applied: Option<(bool, bool)>,
}

impl App {
    pub fn new(start: Startup) -> App {
        let mut app = App {
            docs: Vec::new(),
            tab: Tab::Doc(0),
            last_doc: 0,
            enzymes: default_enzymes(),
            selected: BTreeSet::new(),
            plus_unique: true,
            library: builtin_library(),
            prefs: Prefs {
                dark: start.dark,
                theme: if start.dark {
                    ThemeMode::Dark
                } else {
                    ThemeMode::System
                },
                ..Prefs::default()
            },
            status: String::new(),
            map_linear: false,
            enzyme_win: EnzymeWin::default(),
            digest_win: DigestWin::default(),
            enz_list: EnzListWin::default(),
            find_win: FindWin::default(),
            assembly: AssemblyWin::default(),
            orf_win: OrfWin::default(),
            site_win: SiteToolsWin::default(),
            feature_win: FeatureWin::default(),
            align_win: AlignWin::default(),
            pcr_win: PcrWin::default(),
            calc_win: CalcWin::default(),
            calculator_win: CalculatorWin::default(),
            molcalc_win: MolCalcWin::default(),
            tree_win: TreeWin::default(),
            phylocnn_win: PhyloCnnWin::default(),
            species_win: SpeciesWin::default(),
            treetools_win: TreeToolsWin::default(),
            workspace: if start.view.is_empty() && start.files.is_empty() {
                Workspace::Home
            } else {
                Workspace::Sequences
            },
            platform: Platform::new(start.settings_dir.clone()),
            p3_win: Primer3Win::default(),
            reads_win: ReadsWin::default(),
            about_open: false,
            confirm_close: None,
            palette: CommandPalette::default(),
            toasts: Toasts::default(),
            shortcuts_open: false,
            paste_rc_pending: false,
            shot: start.screenshot.clone(),
            shot_asked: false,
            frames: 0,
            theme_applied: None,
        };
        for f in &start.files {
            app.open_file(Path::new(f));
        }
        if app.docs.is_empty() {
            app.docs.push(DocView::new(crate::sample::demo_plasmid()));
        }
        app.apply_startup_view(&start.view);
        app
    }

    /// The sequence tab being worked on; while the assembler tab is showing, the last one viewed.
    pub fn current_doc(&self) -> Option<usize> {
        if self.docs.is_empty() {
            return None;
        }
        match self.tab {
            Tab::Doc(i) if i < self.docs.len() => Some(i),
            _ => Some(self.last_doc.min(self.docs.len() - 1)),
        }
    }

    pub(super) fn say(&mut self, text: impl Into<String>) {
        let text = text.into();
        // every message is a toast too, so the same message twice is seen twice (while the first is still up it is one)
        self.toasts.push(&text, super::toasts::tone_of(&text));
        self.status = text;
    }

    // -----------------------------------------------------------------------
    // Actions
    // -----------------------------------------------------------------------

    fn run(&mut self, action: Action, ctx: &egui::Context) {
        // the palette and the menus are on every page: what works on the open sequence takes the person to it first
        if super::commands::acts_on_the_sequence(&action) && self.workspace != Workspace::Sequences
        {
            self.workspace = Workspace::Sequences;
        }
        let cur = self.current_doc();
        match action {
            Action::Go(ws) => self.workspace = ws,
            Action::Example(ws) => self.load_example(ws),
            Action::New => self.add_doc(Doc::new("Untitled", Vec::new())),
            Action::Open => self.pick_open(),
            Action::DuplicateSelection => self.duplicate_current(false),
            Action::DuplicateFull => self.duplicate_current(true),
            Action::Save => cur.into_iter().for_each(|i| self.save_doc(i, false, false)),
            Action::SaveAs => cur.into_iter().for_each(|i| self.save_doc(i, true, false)),
            Action::ExportFasta => cur.into_iter().for_each(|i| self.save_doc(i, true, true)),
            Action::CloseCurrent => cur.into_iter().for_each(|i| self.request_close(i)),
            Action::CloseTab(i) => self.request_close(i),
            Action::SwitchTab(t) => self.tab = t,
            Action::Find => self.find_win.open = true,
            Action::FindAgain => self.find_step(Dir::Next),
            Action::Jump => self.open_jump(),
            Action::SelectRange => self.open_range(),
            Action::EnzymeSelector => self.enzyme_win.open = true,
            Action::Digest => self.digest_win.open = true,
            Action::ListEnzymes => self.open_list_dialog(),
            Action::QuickList(quick) => self.run_quick_list(quick),
            Action::SilentSites => self.run_silent_sites(),
            Action::DiagnosticSites => self.site_win.diagnostic_open = true,
            Action::ToggleMap => self.prefs.show_map = !self.prefs.show_map,
            Action::ToggleFeatures => self.prefs.show_features = !self.prefs.show_features,
            Action::Theme(mode) => self.prefs.theme = mode,
            Action::Palette => self.palette.toggle(),
            Action::Shortcuts => self.shortcuts_open = true,
            Action::ToggleSecondStrand => self.prefs.second_strand = !self.prefs.second_strand,
            Action::OpenAssembler => self.tab = Tab::Assembly,
            Action::GibsonAssembler => {
                self.tab = Tab::Assembly;
                self.assembly.set_preset(Preset::Gibson);
                if self.assembly.reads.is_empty() {
                    self.add_open_sequences();
                }
            }
            Action::AssembleNow => {
                self.tab = Tab::Assembly;
                self.assembly.run();
            }
            Action::AddOpenToAssembly => {
                self.tab = Tab::Assembly;
                self.add_open_sequences();
            }
            Action::DemoReads => {
                self.tab = Tab::Assembly;
                self.assembly.load_demo();
            }
            Action::OrfNext => self.orf_step(OrfStep::Next),
            Action::OrfPrev => self.orf_step(OrfStep::Prev),
            Action::OrfList => self.orf_win.list_open = true,
            Action::CopySpecial(what) => {
                if let Some(i) = cur {
                    let table = GENETIC_CODES[self.prefs.code].table;
                    let text = clipboard::copy_special(&self.docs[i], what, table);
                    self.say(format!("Copied {} characters.", text.chars().count()));
                    ctx.copy_text(text);
                }
            }
            Action::Translate => self.open_translate(),
            Action::NewFeature => self.open_feature_window(None),
            Action::FeaturesFromFound => self.open_found_feature_window(),
            Action::EditFeature => match cur.and_then(|i| self.docs[i].selected_feature) {
                Some(k) => self.open_feature_window(Some(k)),
                None => self
                    .say("Click a feature (in the text or the feature table) to choose it first."),
            },
            Action::DeleteFeature => self.delete_selected_feature(),
            Action::ClearFeatures => self.clear_features(),
            Action::Annotate => self.annotate_with_library(),
            Action::OpenLibrary => self.open_feature_library(),
            Action::Align => self.align_win.open = true,
            Action::Pcr => self.pcr_win.open = true,
            Action::Calculator => {
                self.calculator_win.fresh = !self.calculator_win.open;
                self.calculator_win.open = true;
            }
            Action::TmCalculator => self.calc_win.open = true,
            Action::MolecularCalculator => self.molcalc_win.open = true,
            Action::Tree => self.tree_win.open = true,
            Action::PhyloCnn => self.phylocnn_win.open = true,
            Action::SpeciesTree => self.species_win.open = true,
            Action::TreeTools => self.treetools_win.open = true,
            Action::Workflow => self.workspace = Workspace::Pipelines,
            Action::Primer3 => self.p3_win.open = true,
            Action::Reads => self.reads_win.open = true,
            Action::About => self.about_open = true,
            other => {
                if let Some(i) = cur {
                    self.edit_action(i, other, ctx);
                }
            }
        }
    }

    fn edit_action(&mut self, i: usize, action: Action, ctx: &egui::Context) {
        let dv = &mut self.docs[i];
        match action {
            Action::Undo => dv.undo(),
            Action::Redo => dv.redo(),
            Action::Copy => clipboard::copy(ctx, dv, false),
            Action::CopyRevCom => clipboard::copy(ctx, dv, true),
            Action::Cut | Action::CutRevCom => {
                clipboard::copy(ctx, dv, matches!(action, Action::CutRevCom));
                if dv.selection().is_some() {
                    dv.replace_selection(b"");
                }
            }
            Action::Paste | Action::PasteRevCom => {
                let rc = matches!(action, Action::PasteRevCom);
                ctx.send_viewport_cmd(egui::ViewportCommand::RequestPaste);
                self.paste_rc_pending = rc;
            }
            Action::SelectAll => {
                let n = dv.doc.len();
                dv.select(0, n);
            }
            Action::ClearFind => dv.find_hits.clear(),
            Action::Upper => map_selection(dv, |c| c.to_ascii_uppercase()),
            Action::Lower => map_selection(dv, |c| c.to_ascii_lowercase()),
            Action::SwapCase => map_selection(dv, swap_case),
            Action::RevCom => {
                let n = dv.doc.len();
                match dv.selection() {
                    // part of the sequence: only that stretch is turned around
                    Some((a, b)) if (a, b) != (0, n) => {
                        dv.edit(|d| d.reverse_complement_range(a, b))
                    }
                    // everything selected: the whole sequence turns, the selection stays
                    Some(_) => dv.edit(|d| d.reverse_complement()),
                    // nothing selected: the caret keeps its distance from the other end
                    None => {
                        let caret = dv.caret;
                        dv.edit(|d| d.reverse_complement());
                        dv.place_caret(n - caret);
                    }
                }
            }
            Action::InsertSite { enzyme, reverse } => {
                if let Some(e) = self.enzymes.get(enzyme) {
                    let (a, b) = dv.selection().unwrap_or((dv.caret, dv.caret));
                    let bases = e.site_to_insert(reverse);
                    dv.edit(|d| {
                        d.delete(a, b);
                        d.insert_site(a, &e.name, &bases, reverse);
                    });
                    dv.select(a, a + bases.len());
                    dv.scroll_to = Some(a);
                }
            }
            Action::ToggleCircular => {
                dv.doc.circular = !dv.doc.circular;
                dv.doc.dirty = true;
            }
            Action::SetOrigin if dv.doc.circular => {
                let at = dv.caret;
                dv.edit(|d| d.set_origin(at));
                dv.place_caret(0);
            }
            Action::Linearize if dv.doc.circular => {
                let at = dv.caret;
                dv.edit(|d| {
                    d.set_origin(at);
                    d.circular = false;
                });
                dv.place_caret(0);
            }
            Action::ToggleHighlight => dv.highlight_sites = !dv.highlight_sites,
            Action::Translation(mode) => {
                dv.trans = if dv.trans == mode {
                    TransMode::Off
                } else {
                    mode
                }
            }
            _ => {}
        }
    }

    // -----------------------------------------------------------------------
    // Frame
    // -----------------------------------------------------------------------

    fn shortcuts(&self, ctx: &egui::Context, out: &mut Vec<Action>) {
        // Cmd+Alt+1 to 9 go to the workspaces. They come before the table, where Cmd+1 and Cmd+2 would also take them
        // (an Alt that is pressed too is ignored by a shortcut without one).
        let numbers = [
            Key::Num1,
            Key::Num2,
            Key::Num3,
            Key::Num4,
            Key::Num5,
            Key::Num6,
            Key::Num7,
            Key::Num8,
            Key::Num9,
        ];
        for (key, ws) in numbers.into_iter().zip(Workspace::ALL) {
            let shortcut = KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::ALT, key);
            if ctx.input_mut(|i| i.consume_shortcut(&shortcut)) {
                out.push(Action::Go(ws));
            }
        }
        // Shortcuts with more modifiers come first: a plain Cmd+S would otherwise also swallow Cmd+Shift+S.
        let table = [
            (
                Modifiers::COMMAND | Modifiers::SHIFT,
                Key::S,
                Action::SaveAs,
            ),
            (
                Modifiers::COMMAND | Modifiers::SHIFT,
                Key::K,
                Action::Annotate,
            ),
            (
                Modifiers::COMMAND | Modifiers::SHIFT,
                Key::A,
                Action::SelectRange,
            ),
            (
                Modifiers::COMMAND | Modifiers::SHIFT,
                Key::F,
                Action::ClearFind,
            ),
            (
                Modifiers::COMMAND | Modifiers::SHIFT,
                Key::Period,
                Action::OrfNext,
            ),
            (
                Modifiers::COMMAND | Modifiers::SHIFT,
                Key::Comma,
                Action::OrfPrev,
            ),
            (
                Modifiers::COMMAND,
                Key::B,
                Action::CopySpecial(Special::Genbank),
            ),
            (Modifiers::COMMAND, Key::T, Action::Translate),
            (Modifiers::COMMAND, Key::L, Action::Align),
            (Modifiers::COMMAND, Key::K, Action::Palette),
            (Modifiers::COMMAND, Key::Num1, Action::Calculator),
            (Modifiers::COMMAND, Key::Num2, Action::MolecularCalculator),
            (Modifiers::COMMAND, Key::Period, Action::NewFeature),
            (Modifiers::COMMAND, Key::N, Action::New),
            (Modifiers::COMMAND, Key::O, Action::Open),
            (Modifiers::COMMAND, Key::S, Action::Save),
            (Modifiers::COMMAND, Key::W, Action::CloseCurrent),
            (Modifiers::COMMAND, Key::F, Action::Find),
            (Modifiers::COMMAND, Key::G, Action::FindAgain),
            (Modifiers::COMMAND, Key::J, Action::Jump),
            (Modifiers::COMMAND, Key::E, Action::EnzymeSelector),
            (Modifiers::COMMAND, Key::R, Action::Digest),
            (Modifiers::COMMAND, Key::I, Action::SetOrigin),
            (Modifiers::COMMAND, Key::Slash, Action::RevCom),
        ];
        for (mods, key, action) in table {
            if ctx.input_mut(|i| i.consume_shortcut(&KeyboardShortcut::new(mods, key))) {
                out.push(action);
            }
        }
    }

    pub fn palette(&self) -> Palette {
        theme::palette(self.prefs.dark)
    }

    // -----------------------------------------------------------------------
    // Screenshot hook
    // -----------------------------------------------------------------------

    fn screenshot_step(&mut self, ctx: &egui::Context) {
        let Some(path) = self.shot.clone() else {
            return;
        };
        self.frames += 1;
        ctx.request_repaint();
        // a window that isn't in front can be throttled by the system; ask for focus while drawing
        if self.frames == 1 || self.frames == SCREENSHOT_FRAMES / 2 {
            ctx.send_viewport_cmd(egui::ViewportCommand::Focus);
        }
        if self.frames == SCREENSHOT_FRAMES && !self.shot_asked {
            self.shot_asked = true;
            ctx.send_viewport_cmd(egui::ViewportCommand::Screenshot(egui::UserData::default()));
        }
        let image = ctx.input(|i| {
            i.events.iter().find_map(|e| match e {
                egui::Event::Screenshot { image, .. } => Some(image.clone()),
                _ => None,
            })
        });
        if let Some(img) = image {
            // if the window was too slow, the picture is being made without it, and that ends the
            // program: only the one that takes the picture closes the window
            if claim_picture() {
                if let Err(e) = save_png(&path, &img) {
                    eprintln!("could not save the screenshot: {e}");
                }
                ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            }
        }
    }
}

fn map_selection(dv: &mut DocView, f: impl Fn(u8) -> u8) {
    let (a, b) = dv.selection().unwrap_or((0, dv.doc.len()));
    let mapped: Vec<u8> = dv.doc.slice(a, b).iter().map(|&c| f(c)).collect();
    dv.edit(|d| d.overwrite(a, &mapped));
}

impl eframe::App for App {
    /// While taking a screenshot, ignore the mouse so the picture doesn't depend on where the
    /// pointer happens to rest (hover read-outs, stray scroll events).
    fn raw_input_hook(&mut self, _ctx: &egui::Context, raw_input: &mut egui::RawInput) {
        if self.shot.is_some() {
            raw_input.events.retain(|e| {
                !matches!(
                    e,
                    egui::Event::PointerMoved(_)
                        | egui::Event::PointerButton { .. }
                        | egui::Event::PointerGone
                        | egui::Event::MouseWheel { .. }
                        | egui::Event::Touch { .. }
                )
            });
        }
    }

    fn ui(&mut self, ui: &mut Ui, _frame: &mut eframe::Frame) {
        self.draw(ui);
        let ctx = ui.ctx().clone();
        self.screenshot_step(&ctx);
    }
}

impl App {
    /// Everything on screen for one frame. Kept apart from the window so that the picture can be
    /// drawn with no window at all (see `headless`).
    pub(super) fn draw(&mut self, ui: &mut Ui) {
        let ctx = ui.ctx().clone();
        // whether a text field has the keyboard as the frame begins: a key it takes (Return in Find's) is not for a dialog
        let typing = ctx.text_edit_focused();
        // the colours follow the system until the person chooses light or dark
        self.prefs.dark = self.prefs.theme.is_dark(ctx.system_theme());
        let colours = (self.prefs.dark, self.prefs.theme == ThemeMode::System);
        if self.theme_applied != Some(colours) {
            theme::apply(&ctx, colours.0, colours.1);
            self.theme_applied = Some(colours);
        }
        self.take_reverse_paste(&ctx);
        let mut actions: Vec<Action> = Vec::new();
        self.shortcuts(&ctx, &mut actions);
        // while the palette is open it has the keyboard: only the key that closes it goes through
        if self.palette.is_open() {
            actions.retain(|a| matches!(a, Action::Palette));
        }
        let look = design::look(self.prefs.dark);
        egui::Panel::top("topbar")
            .exact_size(super::topbar::HEIGHT)
            .frame(
                egui::Frame::new()
                    .fill(look.surface)
                    .inner_margin(egui::Margin::symmetric(12, 0)),
            )
            .show(ui, |ui| self.top_bar(ui, &mut actions));
        egui::Panel::left("rail")
            .exact_size(RAIL_WIDTH)
            .resizable(false)
            .frame(egui::Frame::new().fill(look.surface))
            .show(ui, |ui| self.rail(ui, &mut actions));
        if self.workspace == Workspace::Sequences {
            egui::Panel::top("toolbar").show(ui, |ui| self.toolbar(ui, &mut actions));
            egui::Panel::top("tabs").show(ui, |ui| self.tab_strip(ui, &mut actions));
        }
        egui::Panel::bottom("status").show(ui, |ui| self.status_bar(ui));
        egui::CentralPanel::default()
            .frame(egui::Frame::new().fill(look.bg))
            .show(ui, |ui| self.workspace_content(ui, &mut actions));
        self.enzyme_window(&ctx);
        self.digest_window(&ctx);
        self.find_window(&ctx);
        self.assembly_windows(&ctx);
        self.orf_window(&ctx);
        self.site_windows(&ctx);
        self.enz_list_windows(&ctx, typing);
        self.feature_window(&ctx);
        self.align_window(&ctx);
        self.pcr_window(&ctx);
        self.calc_window(&ctx);
        self.calculator_window(&ctx);
        self.molcalc_window(&ctx);
        self.tree_window(&ctx);
        self.phylocnn_window(&ctx);
        self.species_window(&ctx);
        self.treetools_window(&ctx);
        self.primer3_window(&ctx);
        self.reads_window(&ctx);
        self.about_window(&ctx);
        self.shortcuts_window(&ctx);
        self.close_prompt(&ctx);
        self.paste_follow_up(&ctx);
        self.take_dropped_files(&ctx);
        self.palette_ui(&ctx, &mut actions);
        self.toasts_ui(&ctx);
        for a in actions {
            self.run(a, &ctx);
        }
    }

    /// "Paste Rev-Com" asks for the clipboard; when it arrives, insert it reverse-complemented
    /// and hide it from the editor so it isn't pasted twice.
    fn take_reverse_paste(&mut self, ctx: &egui::Context) {
        if !self.paste_rc_pending {
            return;
        }
        let text = ctx.input_mut(|i| {
            let mut found = None;
            i.events.retain(|e| match e {
                egui::Event::Paste(s) => {
                    found = Some(s.clone());
                    false
                }
                _ => true,
            });
            found
        });
        if let Some(t) = text {
            self.paste_rc_pending = false;
            if let Some(i) = self.current_doc() {
                clipboard::paste(ctx, &mut self.docs[i], &t, true);
            }
        }
    }
}

mod startup;

#[cfg(test)]
mod tests;
