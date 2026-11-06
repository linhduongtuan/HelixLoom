//! State shared by the windows: open documents (with selection, undo and cached enzyme sites)
//! and preferences.

use super::clipboard::PasteJob;
use crate::seq::dna::{Hit, StartRule, Strand};
use crate::seq::doc::Doc;
use crate::seq::enzymes::{all_sites, Enzyme, Site};
use crate::seq::primers::PrimerParams;
use std::sync::atomic::{AtomicU64, Ordering};

const UNDO_LIMIT: usize = 100;

/// The next number a sequence opened in this run is told apart by.
static NEXT_VIEW_ID: AtomicU64 = AtomicU64::new(1);

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum TransMode {
    Off,
    Forward,
    Six,
}

impl TransMode {
    pub fn rows(self) -> usize {
        match self {
            TransMode::Off => 0,
            TransMode::Forward => 3,
            TransMode::Six => 6,
        }
    }
}

struct SiteCache {
    rev: u64,
    methylated: bool,
    circular: bool,
    enzyme_count: usize,
    sites: Vec<Vec<Site>>,
}

pub struct DocView {
    pub doc: Doc,
    /// Told apart from every other sequence opened in this run (not its name, which two can share): a result that was
    /// worked out from this sequence remembers it.
    pub id: u64,
    undo: Vec<Doc>,
    redo: Vec<Doc>,
    /// Bumped on every edit; caches key off it.
    pub rev: u64,
    pub caret: usize,
    pub anchor: Option<usize>,
    pub scroll_to: Option<usize>,
    pub find_hits: Vec<Hit>,
    /// A paste that still needs the application (open a record, or answer a question).
    pub paste_job: Option<PasteJob>,
    pub highlight_sites: bool,
    pub trans: TransMode,
    pub hover: Option<usize>,
    pub selected_feature: Option<usize>,
    cache: Option<SiteCache>,
}

impl DocView {
    pub fn new(doc: Doc) -> DocView {
        DocView {
            doc,
            id: NEXT_VIEW_ID.fetch_add(1, Ordering::Relaxed),
            undo: Vec::new(),
            redo: Vec::new(),
            rev: 0,
            caret: 0,
            anchor: None,
            scroll_to: None,
            find_hits: Vec::new(),
            paste_job: None,
            highlight_sites: false,
            trans: TransMode::Off,
            hover: None,
            selected_feature: None,
            cache: None,
        }
    }

    /// Selected range `(start, end)` with `start < end`.
    pub fn selection(&self) -> Option<(usize, usize)> {
        self.anchor
            .map(|a| (a.min(self.caret), a.max(self.caret)))
            .filter(|(a, b)| a != b)
    }

    pub fn select(&mut self, a: usize, b: usize) {
        let n = self.doc.len();
        self.anchor = Some(a.min(n));
        self.caret = b.min(n);
    }

    pub fn place_caret(&mut self, pos: usize) {
        self.anchor = None;
        self.caret = pos.min(self.doc.len());
    }

    /// Selected bases, or the whole sequence when nothing is selected.
    pub fn selected_or_all(&self) -> &[u8] {
        match self.selection() {
            Some((a, b)) => self.doc.slice(a, b),
            None => &self.doc.seq,
        }
    }

    /// Run an edit with undo support.
    pub fn edit<R>(&mut self, f: impl FnOnce(&mut Doc) -> R) -> R {
        self.undo.push(self.doc.clone());
        if self.undo.len() > UNDO_LIMIT {
            self.undo.remove(0);
        }
        self.redo.clear();
        let out = f(&mut self.doc);
        self.after_change();
        out
    }

    fn after_change(&mut self) {
        self.rev += 1;
        let n = self.doc.len();
        self.caret = self.caret.min(n);
        self.anchor = self.anchor.map(|a| a.min(n));
        self.find_hits.clear();
        self.selected_feature = self
            .selected_feature
            .filter(|&i| i < self.doc.features.len());
    }

    pub fn can_undo(&self) -> bool {
        !self.undo.is_empty()
    }

    pub fn can_redo(&self) -> bool {
        !self.redo.is_empty()
    }

    pub fn undo(&mut self) {
        if let Some(prev) = self.undo.pop() {
            self.redo.push(std::mem::replace(&mut self.doc, prev));
            self.after_change();
        }
    }

    pub fn redo(&mut self) {
        if let Some(next) = self.redo.pop() {
            self.undo.push(std::mem::replace(&mut self.doc, next));
            self.after_change();
        }
    }

    /// Replace the selection (or insert at the caret) and put the caret after the new bases.
    pub fn replace_selection(&mut self, bases: &[u8]) {
        let (a, b) = self.selection().unwrap_or((self.caret, self.caret));
        self.edit(|d| d.replace(a, b, bases));
        self.place_caret(a + bases.len());
        self.scroll_to = Some(self.caret);
    }

    /// Sites of every enzyme, recomputed only after edits or a methylation / topology change.
    pub fn sites(&mut self, enzymes: &[Enzyme]) -> &[Vec<Site>] {
        let (m, c) = (self.doc.methylated, self.doc.circular);
        let stale = self.cache.as_ref().is_none_or(|k| {
            k.rev != self.rev
                || k.methylated != m
                || k.circular != c
                || k.enzyme_count != enzymes.len()
        });
        if stale {
            let sites = all_sites(enzymes, &self.doc.seq, c, m);
            self.cache = Some(SiteCache {
                rev: self.rev,
                methylated: m,
                circular: c,
                enzyme_count: enzymes.len(),
                sites,
            });
        }
        &self.cache.as_ref().expect("cache was just filled").sites
    }
}

#[derive(Clone)]
pub struct Prefs {
    /// Whether the dark colours are in use right now: what `theme` comes to.
    pub dark: bool,
    pub theme: ThemeMode,
    pub font_size: f32,
    pub second_strand: bool,
    pub orf_min: usize,
    pub orf_rule: StartRule,
    pub orf_strand: Strand,
    pub code: usize,
    pub three_letter: bool,
    pub dna_nm: f64,
    pub salt_mm: f64,
    pub primer: PrimerParams,
    pub show_map: bool,
    pub show_features: bool,
}

/// How the colours are chosen: by the system (light by day, dark by night, as the person set it), or always one of them.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ThemeMode {
    #[default]
    System,
    Light,
    Dark,
}

impl ThemeMode {
    pub const ALL: [ThemeMode; 3] = [ThemeMode::System, ThemeMode::Light, ThemeMode::Dark];

    pub fn label(self) -> &'static str {
        match self {
            ThemeMode::System => "Same as the system",
            ThemeMode::Light => "Light",
            ThemeMode::Dark => "Dark",
        }
    }

    /// Whether the dark colours are wanted, given what the system says (nothing, when it does not say).
    pub fn is_dark(self, system: Option<eframe::egui::Theme>) -> bool {
        match self {
            ThemeMode::System => system == Some(eframe::egui::Theme::Dark),
            ThemeMode::Light => false,
            ThemeMode::Dark => true,
        }
    }

    /// The next one in the order System, Light, Dark.
    pub fn next(self) -> ThemeMode {
        match self {
            ThemeMode::System => ThemeMode::Light,
            ThemeMode::Light => ThemeMode::Dark,
            ThemeMode::Dark => ThemeMode::System,
        }
    }
}

impl Default for Prefs {
    fn default() -> Self {
        Prefs {
            dark: false,
            theme: ThemeMode::System,
            font_size: 12.0, // dnafontsize in ApE_Defaults.txt
            second_strand: true,
            orf_min: 90,
            orf_rule: StartRule::Met,
            orf_strand: Strand::Top,
            code: 0,
            three_letter: false,
            dna_nm: 250.0,
            salt_mm: 50.0,
            primer: PrimerParams::default(),
            show_map: true,
            show_features: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::enzymes::default_enzymes;

    fn view() -> DocView {
        DocView::new(Doc::new("t", b"AAAAGAATTCAAAA".to_vec()))
    }

    #[test]
    fn the_theme_follows_the_system_only_when_asked_to() {
        use eframe::egui::Theme;
        assert!(ThemeMode::System.is_dark(Some(Theme::Dark)));
        assert!(!ThemeMode::System.is_dark(Some(Theme::Light)));
        assert!(
            !ThemeMode::System.is_dark(None),
            "no word from the system is the light colours"
        );
        assert!(ThemeMode::Dark.is_dark(Some(Theme::Light)));
        assert!(!ThemeMode::Light.is_dark(Some(Theme::Dark)));
        let mut mode = ThemeMode::System;
        for expected in [ThemeMode::Light, ThemeMode::Dark, ThemeMode::System] {
            mode = mode.next();
            assert_eq!(mode, expected);
        }
    }

    #[test]
    fn selections_are_normalised_and_empty_ones_vanish() {
        let mut v = view();
        v.select(8, 3);
        assert_eq!(v.selection(), Some((3, 8)));
        v.select(5, 5);
        assert_eq!(v.selection(), None);
        v.select(0, 99);
        assert_eq!(v.selection(), Some((0, 14)));
    }

    #[test]
    fn edits_can_be_undone_and_redone() {
        let mut v = view();
        v.select(4, 10);
        v.replace_selection(b"GG");
        assert_eq!(v.doc.seq, b"AAAAGGAAAA");
        assert_eq!(v.caret, 6);
        v.undo();
        assert_eq!(v.doc.seq, b"AAAAGAATTCAAAA");
        assert!(v.can_redo());
        v.redo();
        assert_eq!(v.doc.seq, b"AAAAGGAAAA");
        assert!(!v.can_redo());
    }

    #[test]
    fn site_cache_refreshes_after_an_edit() {
        let enzymes = default_enzymes();
        let eco = enzymes.iter().position(|e| e.name == "EcoRI").unwrap();
        let mut v = view();
        assert_eq!(v.sites(&enzymes)[eco].len(), 1);
        v.select(4, 10);
        v.replace_selection(b"CCCCCC");
        assert_eq!(v.sites(&enzymes)[eco].len(), 0);
        v.undo();
        assert_eq!(v.sites(&enzymes)[eco].len(), 1);
    }

    #[test]
    fn undo_history_is_bounded() {
        let mut v = view();
        for _ in 0..UNDO_LIMIT + 20 {
            v.edit(|d| d.insert(0, b"A"));
        }
        assert_eq!(v.undo.len(), UNDO_LIMIT);
    }
}
