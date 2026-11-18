//! What the command palette can answer from the work in hand, besides its fixed list of commands: a position typed as a command
//! ("go to 1200"), a search ("find GAATTC"), and a few letters of the name of an open sequence, of a feature of the sequence in view
//! or of a recent file. A pure function of the query and a description of the work, so that it can be tested without a window;
//! `palette.rs` shows the answers and runs them.

use super::app::Action;

/// The most answers of each kind, so that a common word does not fill the list, and the longest pattern a search takes.
pub const MOST_SEQUENCES: usize = 3;
pub const MOST_FEATURES: usize = 5;
pub const MOST_RECENT: usize = 3;
pub const MOST_PATTERN: usize = 200;
/// Fewer letters than this say too little to look for a name.
const LEAST_LETTERS: usize = 2;

/// A feature of the sequence in view, as the palette says it: its name and kind, and its first and last base counted from 1.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FeatureInfo {
    pub name: String,
    pub kind: String,
    pub from: usize,
    pub to: usize,
}

/// The work in hand: the open sequences by name, which one is in view and how long it is, its features, and the recent files by name.
pub struct Work<'a> {
    pub sequences: &'a [String],
    pub current: Option<usize>,
    pub bases: usize,
    pub features: &'a [FeatureInfo],
    pub recent: &'a [String],
}

/// One answer: what it is called, the group it is listed under, and what running it does.
#[derive(Clone, Debug)]
pub struct Entry {
    pub group: &'static str,
    pub title: String,
    pub action: Action,
}

/// The answers: those that take a command's place at the top (a position, a search), and those that follow the commands (names).
#[derive(Clone, Debug, Default)]
pub struct Answers {
    pub lead: Vec<Entry>,
    pub trail: Vec<Entry>,
}

/// The base a query asks to go to: "go to 1200", "goto 1200", "jump (to) 1200", "base 1200", "pos(ition) 1200" or the number alone.
fn position(query: &str) -> Option<usize> {
    let q = query.trim().to_lowercase();
    let words = [
        "go to", "goto", "jump to", "jump", "base", "position", "pos",
    ];
    let rest = words.iter().find_map(|w| q.strip_prefix(w)).unwrap_or(&q);
    let digits = rest.trim();
    if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_digit()) {
        return None;
    }
    let n = digits.parse::<usize>().unwrap_or(usize::MAX);
    (n >= 1).then_some(n)
}

/// The pattern a query asks to find: "find GAATTC" (or search, seek), as typed.
fn pattern(query: &str) -> Option<&str> {
    let q = query.trim();
    let lower = q.to_lowercase();
    let word = ["find ", "search ", "seek "]
        .into_iter()
        .find(|w| lower.starts_with(w))?;
    let pattern = q[word.len()..].trim();
    (!pattern.is_empty() && pattern.chars().count() <= MOST_PATTERN).then_some(pattern)
}

pub fn answers(query: &str, work: &Work) -> Answers {
    let mut out = Answers::default();
    let in_view = work.current.is_some();
    if in_view && work.bases > 0 {
        if let Some(n) = position(query) {
            out.lead.push(if n > work.bases {
                Entry {
                    group: "Sequence",
                    title: format!("Jump to the end of the sequence ({} bases)", work.bases),
                    action: Action::JumpTo(work.bases),
                }
            } else {
                Entry {
                    group: "Sequence",
                    title: format!("Jump to base {n}"),
                    action: Action::JumpTo(n),
                }
            });
        }
    }
    if in_view {
        if let Some(p) = pattern(query) {
            out.lead.push(Entry {
                group: "Sequence",
                title: format!("Find {p} in the sequence"),
                action: Action::FindText(p.to_string()),
            });
        }
    }
    // a position or a pattern is not a name, and a query with a keyword before it is for the sequence, not for a search of names
    let q = query.trim().to_lowercase();
    if !out.lead.is_empty()
        || q.chars().count() < LEAST_LETTERS
        || position(query).is_some()
        || pattern(query).is_some()
    {
        return out;
    }
    let has = |text: &str| text.to_lowercase().contains(&q);
    out.trail.extend(
        work.sequences
            .iter()
            .enumerate()
            .filter(|(i, name)| Some(*i) != work.current && has(name))
            .take(MOST_SEQUENCES)
            .map(|(i, name)| Entry {
                group: "Open sequences",
                title: format!("Switch to {name}"),
                action: Action::SwitchDoc(i),
            }),
    );
    if in_view {
        out.trail.extend(
            work.features
                .iter()
                .enumerate()
                .filter(|(_, f)| has(&f.name) || has(&f.kind))
                .take(MOST_FEATURES)
                .map(|(i, f)| Entry {
                    group: "Features",
                    title: format!("{} — {} {}–{}", f.name, f.kind, f.from, f.to),
                    action: Action::SelectFeature(i),
                }),
        );
    }
    out.trail.extend(
        work.recent
            .iter()
            .enumerate()
            .filter(|(_, name)| has(name))
            .take(MOST_RECENT)
            .map(|(i, name)| Entry {
                group: "Recent files",
                title: format!("Open {name}"),
                action: Action::OpenRecent(i),
            }),
    );
    out
}

use super::app::{App, Tab};
use super::platform::Workspace;
use super::windows::find::Dir;
use crate::seq::goto::jump_to;

impl App {
    /// Go to base `n` (counted from 1) of the sequence in view and show it.
    pub(super) fn jump_to_base(&mut self, n: usize) {
        let Some(i) = self.current_doc() else { return };
        let dv = &mut self.docs[i];
        dv.place_caret(jump_to(&n.to_string()));
        dv.scroll_to = Some(dv.caret);
        self.workspace = Workspace::Sequences;
        self.say(format!("Jumped to base {n}."));
    }

    /// Find the next place of `pattern` in the sequence in view, with the settings of the Find window, and say what was found.
    pub(super) fn find_text(&mut self, pattern: String) {
        if self.current_doc().is_none() {
            return;
        }
        self.find_win.text = pattern;
        self.workspace = Workspace::Sequences;
        self.find_step(Dir::Next);
        let said = self.find_win.message.clone();
        self.say(said);
    }

    /// Show the open sequence `i`.
    pub(super) fn switch_doc(&mut self, i: usize) {
        if i < self.docs.len() {
            self.tab = Tab::Doc(i);
            self.workspace = Workspace::Sequences;
        }
    }

    /// Select feature `fi` of the sequence in view, as a click on it does, and show it.
    pub(super) fn select_feature(&mut self, fi: usize) {
        let Some(i) = self.current_doc() else { return };
        let dv = &mut self.docs[i];
        let Some(f) = dv.doc.features.get(fi) else {
            return;
        };
        let (a, b) = (f.start(), f.end());
        dv.selected_feature = Some(fi);
        dv.select(a, b);
        dv.scroll_to = Some(a);
        self.workspace = Workspace::Sequences;
    }
}

#[cfg(test)]
mod tests;
