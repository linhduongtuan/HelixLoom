//! The program's side of `session`: the open sequences are saved a moment after they change (and when the window is closed), the
//! last session is offered when the program starts, the files that were opened are remembered for File ▸ Open Recent, and the window
//! comes back where it was. What the Alignment, Trees, Structure and Pipelines pages show is kept with the sequences (`results`): each result when
//! there is room for it, and never in the way of the sequences.
//!
//! What this is careful about, because it keeps the person's unsaved work:
//! * a change is written when it has been quiet for two seconds, or has waited fifteen, whichever comes first; a sequence that has not
//!   changed since the last save is not written again (its file is kept), so a big one costs nothing while another is edited;
//! * the last session waiting on offer is **not written over**: if the person goes on working before answering, it is put aside
//!   (`session.previous`) and the new work is kept as usual, so closing the window at any moment loses nothing;
//! * two copies of the program do not share a session: the second does not keep one (a lock file);
//! * a save that fails is said once, in red, and tried again in a few seconds.

use super::app::{Action, App, Tab};
use super::appearance::Appearance;
use super::places::Places;
use super::platform::Workspace;
use super::results;
use super::session::{
    push_recent, Lock, SavedDoc, SavedResult, Session, Store, WindowPlace, MOST_RECENT,
};
use super::state::{DocView, ThemeMode};
use crate::seq::doc::Doc;
use crate::seq::io;
use eframe::egui::{self, Id, RichText};
use std::path::{Path, PathBuf};
use std::time::Duration;
// the browser has no `std::time::Instant` (it panics there); `web_time` is the same on a computer and works in a page
use web_time::Instant;

/// A change is written when nothing has changed for this long...
const QUIET: Duration = Duration::from_secs(2);
/// ...or when the first change has waited this long, however much is still going on.
const MOST_WAIT: Duration = Duration::from_secs(15);
/// After a save that failed, the next try is this long after.
const RETRY_AFTER: Duration = Duration::from_secs(5);
/// The same for a move or a resize of the window.
const WINDOW_SAVE_AFTER: Duration = Duration::from_secs(2);
/// The most sequences the offer names.
const MOST_NAMED: usize = 6;
/// What a page says when a page that was opened after it keeps the session now.
pub const TAKEN_OVER: &str =
    "Another tab of HelixLoom keeps the session now, so what you change in this tab is not kept: \
                              save what you need (File menu).";
/// What the File menu calls the switch that says whether the work is remembered from one run to the next: the open
/// sequences and what the pages show (`results`), which are offered again, and the list of recent files.
pub const KEEP_SWITCH: &str = "Remember My Work for Next Time";
/// What the person is told when they turn that off, and on again.
pub const NOT_KEPT_NOW: &str =
    "Your work is no longer remembered: the open sequences, what the pages showed and the recent files that were kept are removed.";
pub const KEPT_AGAIN: &str = "Your work is remembered for next time again.";
/// What a page says when the remembering was turned off in another page of the site.
pub const TURNED_OFF_ELSEWHERE: &str =
    "Remembering was turned off in another tab: the work in this tab is no longer kept.";
/// What a copy of the program says when another copy has the session.
pub const OTHER_COPY: &str =
    "Another copy of HelixLoom is open: this one does not keep a session, so that the two do not write over each other's.";
/// What the command says where there is nowhere to keep anything.
pub const NOWHERE_TO_KEEP: &str =
    "Nothing is remembered here: this copy has no place to keep anything.";
/// What the command says while the last session is on offer: it would remove that session without the question about
/// work that was not saved.
pub const ANSWER_FIRST: &str = "Answer the question about your last session first.";
/// A window of which less than this (in points, each way) is on the screen is moved back onto it.
const LEAST_VISIBLE: f32 = 100.0;
/// Where a window that was off every screen is put.
const BACK_ON_SCREEN: egui::Pos2 = egui::pos2(40.0, 40.0);

/// What one open sequence is, for telling whether it changed since it was saved.
#[derive(Clone, Debug, PartialEq)]
struct DocSig {
    id: u64,
    rev: u64,
    len: usize,
    name: String,
    path: Option<PathBuf>,
    dirty: bool,
    circular: bool,
    methylated: bool,
}

/// What the session is made of at a moment: its sequences, what the pages show, and which workspace and tab were showing.
#[derive(Clone, Debug, PartialEq)]
struct Signature {
    docs: Vec<DocSig>,
    /// The results of the pages (`results`): the kind of each, and a number that is another one when it has changed.
    results: Vec<(&'static str, u64)>,
    workspace: &'static str,
    active: usize,
}

/// Whether the session is to be written now.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Due {
    /// Nothing has changed, or this copy keeps no session.
    Never,
    /// Something has changed, and is waiting for things to go quiet.
    Wait,
    Now,
}

/// What `restore_session` did.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Restored {
    pub opened: usize,
    /// The sequences that could not be read back.
    pub lost: usize,
    /// The sequences whose file was open already, so that they were not opened again.
    pub already_open: usize,
    /// The kinds of the results that are on their pages again, and of those that could not be put back (a text that cannot be
    /// read, or a page that would not take it).
    pub results: Vec<String>,
    pub results_lost: Vec<String>,
}

pub struct Keeper {
    store: Option<Store>,
    /// The lock of the session folder, held as long as the program runs.
    _lock: Option<std::fs::File>,
    /// Whether this copy keeps a session (not when another copy of the program does, or the person said not to).
    sessions: bool,
    /// Whether the person wants the open sequences kept (File menu): they do until they say not to.
    wanted: bool,
    /// Another copy of the program keeps the session (said once; the tests ask).
    #[cfg_attr(not(test), allow(dead_code))]
    pub other_copy: bool,
    saved: Option<Signature>,
    /// What each sequence of the last save was, and the file that holds its text.
    written: Vec<(DocSig, Option<String>)>,
    /// The same for each result that the last save kept.
    written_results: Vec<((&'static str, u64), String)>,
    /// The kinds of the results that the last save left out (said when it is not what the save before left out).
    left_out: Vec<&'static str>,
    seen: Option<Signature>,
    first_change: Option<Instant>,
    last_change: Option<Instant>,
    /// After a failed save, nothing is tried before this.
    blocked_until: Option<Instant>,
    /// The session of the last run, waiting for the person to say whether to reopen it.
    pub offer: Option<Session>,
    /// The offered session was put aside (`session.previous`), because the person went on working before answering.
    archived: bool,
    /// "Start fresh" was clicked for a session with unsaved work, and is asked about once more.
    pub confirm_fresh: bool,
    /// The recent files, the latest first.
    pub recent: Vec<PathBuf>,
    window: Option<WindowPlace>,
    window_changed_at: Option<Instant>,
    window_saved: Option<WindowPlace>,
    /// Where the tool windows were left, as they are now and as they were last written.
    places: Places,
    places_saved: Places,
    places_changed_at: Option<Instant>,
    /// The theme and appearance settings that the settings folder holds.
    appearance_kept: Option<(ThemeMode, Appearance)>,
    /// After a failed write of them, nothing is tried before this.
    appearance_retry_at: Option<Instant>,
    /// The last write of them failed and has not been made good yet.
    pub appearance_failed: bool,
    /// Whether the window was looked at to see that it is on a screen.
    on_screen_checked: bool,
    /// In a browser, the mark this page took the session over with (`session::Lock::Taken`).
    owner: Option<String>,
    /// Something to tell the person, once.
    problem: Option<String>,
    said_problem: Option<String>,
    /// What the larger table's refusals were told as in this visit: each reason is told once, however many saves it spoils.
    said_late: Option<String>,
}

/// Take the session of `store` for this copy of the program: the lock that says so, whether this copy keeps the session,
/// whether another copy does, and in a browser the mark this page took it over with. With no store nothing is kept.
fn taken(store: Option<&Store>) -> (Option<std::fs::File>, bool, bool, Option<String>) {
    match store.map(Store::lock_session) {
        Some(Lock::Held(file)) => (Some(file), true, false, None),
        Some(Lock::Busy) => (None, false, true, None),
        Some(Lock::Unavailable) => (None, true, false, None),
        Some(Lock::Taken(mark)) => (None, true, false, Some(mark)),
        None => (None, false, false, None),
    }
}

impl Keeper {
    /// A keeper for the settings folder `dir`, with the last session on offer if it left one. With no folder it keeps things
    /// in the browser's storage when the program runs in a page, and nothing anywhere else.
    pub fn new(dir: Option<&Path>) -> Keeper {
        Keeper::with(Store::open(dir))
    }

    /// A keeper for `store` (none: nothing is read or kept).
    pub(super) fn with(store: Option<Store>) -> Keeper {
        // a copy that is to keep nothing does not take the session from one that may
        let wanted = store.as_ref().is_none_or(Store::sessions_wanted);
        let (lock, sessions, other_copy, owner) = taken(store.as_ref().filter(|_| wanted));
        // nothing is to be remembered: a session that is there all the same (a page that had not heard yet wrote it) and
        // a list of recent files are removed, not left where the next person could find them
        let left_over = store
            .as_ref()
            .filter(|_| !wanted)
            .and_then(|store| store.forget_sessions().and(store.clear_recent()).err())
            .map(|why| format!("What was not to be remembered could not be removed: {why}"));
        let offer = if sessions {
            store.as_ref().and_then(Store::load_session)
        } else {
            None
        };
        let recent = store
            .as_ref()
            .filter(|_| wanted)
            .map(Store::recent)
            .unwrap_or_default();
        let window = store.as_ref().and_then(Store::window);
        let places = store.as_ref().map(Store::places).unwrap_or_default();
        Keeper {
            store,
            _lock: lock,
            sessions,
            wanted,
            other_copy,
            saved: None,
            written: Vec::new(),
            written_results: Vec::new(),
            left_out: Vec::new(),
            seen: None,
            first_change: None,
            last_change: None,
            blocked_until: None,
            offer,
            archived: false,
            confirm_fresh: false,
            recent,
            window,
            window_changed_at: None,
            window_saved: window,
            places_saved: places.clone(),
            places,
            places_changed_at: None,
            appearance_kept: None,
            appearance_retry_at: None,
            appearance_failed: false,
            on_screen_checked: false,
            problem: other_copy.then(|| OTHER_COPY.to_string()).or(left_over),
            said_problem: None,
            said_late: None,
            owner,
        }
    }

    /// Whether the session is still this copy's to write. In a browser a page that was opened after this one has taken it
    /// over: this page then stops keeping it, and says so once.
    fn still_keeps(&mut self) -> bool {
        if !self.sessions {
            return false;
        }
        // another page of the site was told to remember nothing: that goes for this one too
        if self.store.as_ref().is_some_and(|s| !s.sessions_wanted()) {
            self.stop_remembering();
            self.problem = Some(TURNED_OFF_ELSEWHERE.to_string());
            return false;
        }
        let ours = self
            .store
            .as_ref()
            .is_some_and(|store| store.owns(self.owner.as_deref()));
        if !ours {
            self.sessions = false;
            self.problem = Some(TAKEN_OVER.to_string());
        }
        ours
    }

    /// Where the tool windows were left, for `design::window` to open them at.
    pub fn places_for_windows(&self) -> std::collections::HashMap<String, egui::Pos2> {
        self.places
            .iter()
            .map(|(title, &(x, y))| (title.clone(), egui::pos2(x, y)))
            .collect()
    }

    /// The theme and appearance settings that the settings folder holds, none when it holds none (or there is no folder).
    pub fn saved_appearance(&self) -> Option<(ThemeMode, Appearance)> {
        self.store.as_ref().and_then(Store::appearance)
    }

    /// Say what the theme and appearance settings are that the window starts with (what the settings folder holds, or the defaults when
    /// it holds none), so that nothing is written until the person changes one.
    pub fn appearance_starts_as(&mut self, theme: ThemeMode, appearance: Appearance) {
        self.appearance_kept = Some((theme, appearance));
    }

    /// Write the theme and appearance settings when they are not what the settings folder holds. A write that fails is told once, and
    /// tried again every `RETRY_AFTER` (`appearance_failed` is set for as long as it has not gone through). Returns how long to wait
    /// before the next try, when there is one to make.
    pub fn keep_appearance(
        &mut self,
        theme: ThemeMode,
        appearance: &Appearance,
        now: Instant,
    ) -> Option<Duration> {
        let store = self.store.as_ref()?;
        let wanted = (theme, *appearance);
        if self.appearance_kept == Some(wanted) {
            // what the folder holds is what is wanted again (a choice that was changed and changed back)
            self.appearance_failed = false;
            self.appearance_retry_at = None;
            return None;
        }
        if let Some(at) = self.appearance_retry_at.filter(|at| now < *at) {
            return Some(at - now);
        }
        match store.save_appearance(theme, appearance) {
            Ok(()) => {
                self.appearance_kept = Some(wanted);
                self.appearance_retry_at = None;
                self.appearance_failed = false;
                None
            }
            Err(why) => {
                self.problem = Some(why);
                self.appearance_retry_at = Some(now + RETRY_AFTER);
                self.appearance_failed = true;
                Some(RETRY_AFTER)
            }
        }
    }

    /// Whether there is a settings folder at all (a picture or a test has none, and makes no notices outside the window either).
    pub fn has_settings(&self) -> bool {
        self.store.is_some()
    }

    /// Whether the work is remembered for the next run (the open sequences and the recent files): the person's choice,
    /// as it stands now, so that a page shows at once what another page of the site was told. None where there is nothing
    /// to choose: nothing can be kept, or another copy of the program keeps the session.
    pub fn keeps_work(&self) -> Option<bool> {
        let store = self.store.as_ref().filter(|_| !self.other_copy)?;
        Some(self.wanted && store.sessions_wanted())
    }

    /// What the command says where the switch has nothing to choose.
    pub fn why_nothing_to_choose(&self) -> &'static str {
        if self.other_copy {
            OTHER_COPY
        } else {
            NOWHERE_TO_KEEP
        }
    }

    /// This copy remembers nothing from here on, and holds nothing of a session: no lock, no mark, no offer, and nothing
    /// that it believes to be written.
    fn stop_remembering(&mut self) {
        // the lock is a file of the folder that is removed: it is let go before that
        (self._lock, self.sessions, self.wanted, self.owner) = (None, false, false, None);
        (self.offer, self.archived, self.confirm_fresh) = (None, false, false);
        self.forget_what_was_written();
        self.recent.clear();
    }

    /// What this copy believes to be in the session is no longer so: the next save writes everything.
    fn forget_what_was_written(&mut self) {
        (self.saved, self.seen, self.first_change, self.last_change) = (None, None, None, None);
        (self.blocked_until, self.said_problem) = (None, None);
        self.written.clear();
        self.written_results.clear();
        self.left_out.clear();
    }

    /// Turn the remembering on or off, for this run and the next. Off removes what was kept: the session, one that is on
    /// offer, one that was put aside, one of another version, and the list of recent files. What to tell the person;
    /// nothing where there is nothing to choose, or nothing changes.
    pub fn set_keeping(&mut self, wanted: bool) -> Option<String> {
        if self.keeps_work()? == wanted {
            return None;
        }
        if wanted {
            return Some(self.remember_again());
        }
        self.stop_remembering();
        let store = self.store.as_ref()?;
        // what was kept goes first, whatever comes of noting the choice: a storage that is full has room for the note
        // only after that, and the person is leaving the computer
        let removed = store.forget_sessions().and(store.clear_recent());
        let noted = store.want_sessions(false);
        Some(match (removed, noted) {
            (Ok(()), Ok(())) => NOT_KEPT_NOW.to_string(),
            (Err(why), Ok(())) => format!(
                "Your work is no longer remembered, but what was kept could not be removed: {why}"
            ),
            (Ok(()), Err(why)) => format!(
                "What was kept is removed and nothing more is kept in this run, but that could not be noted for the next: {why}"
            ),
            (Err(removing), Err(noting)) => format!(
                "Nothing more is kept in this run, but what was kept could not be removed ({removing}) and that could not be noted for the next ({noting})"
            ),
        })
    }

    /// Take the session for this copy again and note that the work is to be remembered; what to tell the person.
    fn remember_again(&mut self) -> String {
        let Some(store) = &self.store else {
            return NOWHERE_TO_KEEP.to_string();
        };
        // whether the site was set to remember nothing, or only this copy was (its choice could not be noted)
        let was_off = !store.sessions_wanted();
        if let Err(why) = store.want_sessions(true) {
            return why;
        }
        self.wanted = true;
        (self._lock, self.sessions, self.other_copy, self.owner) = taken(Some(store));
        if self.other_copy {
            return OTHER_COPY.to_string();
        }
        // a session that would not go when the remembering was turned off is not to come back as an offer. Where only
        // this copy had stopped, the session is the work of a copy that went on, and stays
        let left_over = if was_off {
            store.forget_sessions()
        } else {
            Ok(())
        };
        self.recent = store.recent();
        match left_over {
            Ok(()) => KEPT_AGAIN.to_string(),
            Err(why) => format!(
                "Your work is remembered again, but a session that was to be removed is still there: {why}"
            ),
        }
    }

    /// Whether the settings are kept by a browser for the page's site, and not in a folder of the computer.
    pub fn kept_in_a_browser(&self) -> bool {
        self.store.as_ref().is_some_and(Store::in_a_browser)
    }

    /// Where the window was left, to start it there.
    pub fn window_place(dir: Option<&Path>) -> Option<WindowPlace> {
        Store::open(dir).and_then(|s| s.window())
    }

    /// Remember that `path` was opened or saved: it goes to the top of the recent files.
    pub fn remember(&mut self, path: &Path) {
        // not while the person asked for nothing to be remembered: the name of a file says what they work on
        let Some(store) = self
            .store
            .as_ref()
            .filter(|s| self.wanted && s.sessions_wanted())
        else {
            return;
        };
        self.recent = push_recent(&self.recent, path, MOST_RECENT);
        if let Err(why) = store.add_recent(path) {
            self.problem = Some(why);
        }
    }

    pub fn forget_recent(&mut self) {
        self.recent.clear();
        if let Some(Err(why)) = self.store.as_ref().map(Store::clear_recent) {
            self.problem = Some(why);
        }
    }

    /// Whether `signature` is something that is not in the session yet.
    fn pending(&self, signature: &Signature) -> bool {
        if !self.sessions || self.store.is_none() || self.saved.as_ref() == Some(signature) {
            return false;
        }
        // a session of nothing, when nothing was saved before it, has nothing to write
        !(self.saved.is_none() && signature.docs.is_empty() && signature.results.is_empty())
    }

    /// Whether the session is to be written now, given what it is made of at `now`: a change waits until things have been quiet for
    /// a moment, or until it has waited a long time.
    fn due(&mut self, signature: &Signature, now: Instant) -> Due {
        if !self.pending(signature) {
            self.seen = None;
            self.first_change = None;
            self.last_change = None;
            return Due::Never;
        }
        if self.seen.as_ref() != Some(signature) {
            self.seen = Some(signature.clone());
            self.last_change = Some(now);
            self.first_change.get_or_insert(now);
        }
        if self.blocked_until.is_some_and(|t| now < t) {
            return Due::Wait;
        }
        let waited = |t: Option<Instant>, most: Duration| {
            t.is_some_and(|t| now.saturating_duration_since(t) >= most)
        };
        if waited(self.last_change, QUIET) || waited(self.first_change, MOST_WAIT) {
            Due::Now
        } else {
            Due::Wait
        }
    }

    /// The page's larger table (`session/large.rs`) answers a write later than it is asked, and may answer that it did not happen:
    /// that is told like any failure of a save (once), and the results are written again at the next save, which a change of
    /// the work brings: not at once, or a database that refuses would be asked for ever.
    fn note_late_failure(&mut self) {
        let Some(why) = self.store.as_ref().and_then(Store::shelf_problem) else {
            return;
        };
        // nothing is being kept (the person said not to, or another copy keeps the session): a refused removal is not a refused
        // keeping, and is cleaned up at the next visit with the remembering off
        if !self.wanted || !self.sessions {
            return;
        }
        // what a save believed was written is not so
        self.written_results.clear();
        let said = format!(
            "Could not keep what the pages show for next time: the browser's database refused it ({why})."
        );
        if self.said_late.replace(said.clone()).as_deref() != Some(said.as_str()) {
            self.problem = Some(said);
        }
    }

    /// A failure to tell the person: once, whatever number of times it happens.
    fn fail(&mut self, text: String, now: Instant) {
        self.problem = Some(text);
        self.blocked_until = Some(now + RETRY_AFTER);
        // what was reused from the last save may be what is wrong: the next try writes everything
        self.written.clear();
        self.written_results.clear();
    }

    /// What there is to say, if it is new.
    fn take_problem(&mut self) -> Option<String> {
        let text = self.problem.take()?;
        if self.said_problem.as_deref() == Some(text.as_str()) {
            return None;
        }
        self.said_problem = Some(text.clone());
        Some(text)
    }

    /// Write `session`, which was worked out from `signature`.
    fn write(&mut self, session: Session, signature: Signature, now: Instant) {
        if !self.still_keeps() {
            return;
        }
        let Some(store) = &self.store else {
            return;
        };
        // the last session is waiting on offer and the person has gone on: it is put aside, and what is written is the new work
        if self.offer.is_some() && !self.archived {
            if let Err(why) = store.archive_session() {
                self.fail(format!("Could not save the session: {why}"), now);
                return;
            }
            self.archived = true;
        }
        match store.keep_session(&session) {
            Ok(files) => {
                let sequences = !signature.docs.is_empty();
                let results = signature.results.iter().zip(files.results);
                let (kept, left_out): (Vec<_>, Vec<_>) =
                    results.partition(|(_, file)| file.is_some());
                self.written_results = kept
                    .into_iter()
                    .filter_map(|(stamp, file)| Some((*stamp, file?)))
                    .collect();
                let left_out: Vec<&'static str> =
                    left_out.iter().map(|(stamp, _)| stamp.0).collect();
                // a result that there is no room for is tried at every save, and said when it is news
                let news = !left_out.is_empty() && left_out != self.left_out;
                self.left_out = left_out;
                self.written = signature.docs.iter().cloned().zip(files.docs).collect();
                self.saved = Some(signature);
                (
                    self.seen,
                    self.first_change,
                    self.last_change,
                    self.blocked_until,
                ) = (None, None, None, None);
                self.said_problem = None;
                // what was put aside is of no use once the person has answered the offer and what they reopened is in this session
                // (when it cannot be removed it stays marked as put aside, and the next save tries again)
                if self.archived && self.offer.is_none() && store.clear_previous().is_ok() {
                    self.archived = false;
                }
                if news {
                    self.problem = Some(not_kept(&self.left_out, sequences));
                }
            }
            Err(why) => self.fail(format!("Could not save the session: {why}"), now),
        }
    }
}

/// What is said when results were left out of the session; `sequences`: there are open sequences, which were kept.
fn not_kept(kinds: &[&'static str], sequences: bool) -> String {
    let (it, was) = if kinds.len() == 1 {
        ("it", "it")
    } else {
        ("them", "they")
    };
    // it begins as a failure does, so that the notice is the colour of one (`toasts::tone_of`)
    format!(
        "Could not keep {} for next time: there is no room for {it} here, or {was} could not be written.{}",
        results::named(kinds.iter().copied()),
        if sequences {
            " The open sequences are kept."
        } else {
            ""
        }
    )
}

/// What is said after the last session was reopened: what came back, and what did not.
fn reopened(session: &Session, done: &Restored) -> String {
    let one = |n: usize| n == 1;
    let back = results::named(done.results.iter().map(String::as_str));
    // a kind that this version does not know is one that another version kept: it is not passed over in silence
    let unread: Vec<String> = done
        .results_lost
        .iter()
        .map(|kind| results::called(kind).unwrap_or("a result of another version"))
        .map(str::to_string)
        .collect();
    // "could not be put back" is so whatever the reason: a text that cannot be read, or a page that would not take it now
    let unread = crate::report::listed(&unread);
    let no_sequences = session.docs.is_empty() && session.lost == 0;
    if no_sequences && back.is_empty() {
        return format!(
            "Nothing was put back from the last session: {unread} could not be put back"
        );
    }
    let mut said = if no_sequences {
        format!("Put back {back} from the last session")
    } else {
        format!(
            "Reopened {} sequence{} from the last session{}{back}",
            done.opened,
            if one(done.opened) { "" } else { "s" },
            if back.is_empty() { "" } else { ", with " }
        )
    };
    if done.already_open > 0 {
        let were = if one(done.already_open) {
            "was"
        } else {
            "were"
        };
        said.push_str(&format!("; {} {were} open already", done.already_open));
    }
    if done.lost > 0 {
        let files = if one(done.lost) {
            "its file has"
        } else {
            "their files have"
        };
        said.push_str(&format!(
            "; {} could not be reopened ({files} gone or cannot be read)",
            done.lost
        ));
    }
    if !unread.is_empty() {
        said.push_str(&format!("; {unread} could not be put back"));
    }
    said
}

impl App {
    /// The sequences the session keeps: all but those that nobody has touched and that are not worth coming back to: the example
    /// plasmid that the program opens with (it would be opened again beside the new one at every start), and an empty new one.
    pub(super) fn kept_docs(&self) -> Vec<(usize, &DocView)> {
        self.docs
            .iter()
            .enumerate()
            .filter(|(_, d)| !is_pristine(d))
            .collect()
    }

    fn session_signature(&self) -> Signature {
        let kept = self.kept_docs();
        let active = match self.tab {
            Tab::Doc(i) => i,
            Tab::Assembly => self.last_doc,
        };
        Signature {
            docs: kept
                .iter()
                .map(|(_, d)| DocSig {
                    id: d.id,
                    rev: d.rev,
                    len: d.doc.len(),
                    name: d.doc.name.clone(),
                    path: d.doc.path.clone(),
                    dirty: d.doc.dirty,
                    circular: d.doc.circular,
                    methylated: d.doc.methylated,
                })
                .collect(),
            results: self.result_stamps(),
            workspace: self.workspace.view_name(),
            // the number of the tab among those kept (the first, if the one in front is not kept)
            active: kept.iter().position(|(i, _)| *i == active).unwrap_or(0),
        }
    }

    /// The session as it is now. A sequence that has not changed since the last save is not written again: its file is kept.
    fn session_now(&self, signature: &Signature) -> Session {
        let docs = self
            .kept_docs()
            .iter()
            .zip(&signature.docs)
            .map(|((_, view), sig)| {
                let before = self.keeper.written.iter().find(|(s, _)| s == sig);
                match before {
                    Some((_, Some(file))) if SavedDoc::needs_text(&view.doc) => {
                        let mut saved = SavedDoc::without_text(&view.doc);
                        saved.keep = Some(file.clone());
                        saved
                    }
                    _ => SavedDoc::of(&view.doc),
                }
            })
            .collect();
        // a result that has not changed since the last save is not written again either; one that the last save left out is tried
        let written = &self.keeper.written_results;
        let results = signature
            .results
            .iter()
            .map(|stamp| {
                let file = written.iter().find(|(before, _)| before == stamp);
                SavedResult {
                    kind: stamp.0.to_string(),
                    text: file.is_none().then(|| self.result_text(stamp.0)).flatten(),
                    keep: file.map(|(_, file)| file.clone()),
                }
            })
            .collect();
        Session {
            workspace: signature.workspace.to_string(),
            active: signature.active,
            docs,
            lost: 0,
            results,
        }
    }

    /// Write the session now if there is something in it that is not written (used when the window closes).
    pub fn flush_session(&mut self) {
        // a refusal that came since the last frame is known before what is written is worked out
        self.keeper.note_late_failure();
        let signature = self.session_signature();
        if self.keeper.pending(&signature) {
            let session = self.session_now(&signature);
            self.keeper.write(session, signature, Instant::now());
        }
    }

    /// Keep the session, and the place of the window, up to date; once a frame, and cheap when nothing has changed.
    pub(super) fn keep_work(&mut self, ctx: &egui::Context) {
        // a refusal that came since the last frame is known before what is written is worked out
        self.keeper.note_late_failure();
        let signature = self.session_signature();
        let now = Instant::now();
        let closing = ctx.input(|i| i.viewport().close_requested());
        match self.keeper.due(&signature, now) {
            Due::Now => {
                let session = self.session_now(&signature);
                self.keeper.write(session, signature, now);
                // a save that failed is tried again without waiting for the mouse to move
                if self.keeper.blocked_until.is_some() {
                    ctx.request_repaint_after(RETRY_AFTER);
                }
            }
            Due::Wait => ctx.request_repaint_after(QUIET),
            Due::Never => {}
        }
        // the window is closing: whatever is waiting is written now
        if closing {
            self.flush_session();
        }
        self.keep_window(ctx, now);
        self.keep_places(ctx, now);
        if let Some(why) = self.keeper.take_problem() {
            self.say(why);
        }
    }

    /// Notice where the tool windows drawn in this frame are, and write the places a moment after the last of them stops moving. A
    /// window that is where it opens by itself has no place kept. Nothing is kept without a settings folder.
    pub(super) fn keep_places(&mut self, ctx: &egui::Context, now: Instant) {
        let seen = crate::ui::design::take_seen_windows();
        let keeper = &mut self.keeper;
        if keeper.store.is_none() {
            return;
        }
        for title in seen {
            let Some(rect) = ctx.memory(|m| m.area_rect(crate::ui::design::window_id(&title)))
            else {
                continue;
            };
            let opens = crate::ui::design::default_place(&title);
            let at_default =
                (rect.min.x - opens.x).abs() <= 1.0 && (rect.min.y - opens.y).abs() <= 1.0;
            let changed = if at_default {
                keeper.places.remove(&title)
            } else {
                keeper.places.set(&title, rect.min.x, rect.min.y)
            };
            if changed {
                keeper.places_changed_at = Some(now);
            }
        }
        let (Some(at), Some(store)) = (keeper.places_changed_at, keeper.store.as_ref()) else {
            return;
        };
        if keeper.places == keeper.places_saved {
            keeper.places_changed_at = None;
        } else if now.duration_since(at) >= WINDOW_SAVE_AFTER {
            match store.save_places(&keeper.places) {
                Ok(()) => keeper.places_saved = keeper.places.clone(),
                Err(why) => keeper.problem = Some(why),
            }
            keeper.places_changed_at = None;
        } else {
            ctx.request_repaint_after(WINDOW_SAVE_AFTER);
        }
    }

    /// Notice where the window is and write it a moment after it stops moving; and see, once, that it is on a screen.
    fn keep_window(&mut self, ctx: &egui::Context, now: Instant) {
        self.bring_window_back(ctx);
        let zoom = ctx.zoom_factor();
        let place = ctx.input(|i| {
            let v = i.viewport();
            // not while it is maximised, minimised or full screen: those places are not the one to come back to
            if v.maximized == Some(true) || v.minimized == Some(true) || v.fullscreen == Some(true)
            {
                return None;
            }
            let (outer, inner) = (v.outer_rect?, v.inner_rect?);
            Some(
                WindowPlace {
                    x: outer.min.x,
                    y: outer.min.y,
                    width: inner.width(),
                    height: inner.height(),
                }
                .scaled(zoom),
            )
        });
        let Some(place) = place.filter(WindowPlace::sane) else {
            return;
        };
        let keeper = &mut self.keeper;
        let moved = |a: &WindowPlace, b: &WindowPlace| {
            (a.x - b.x).abs() > 1.0
                || (a.y - b.y).abs() > 1.0
                || (a.width - b.width).abs() > 1.0
                || (a.height - b.height).abs() > 1.0
        };
        if keeper.window.as_ref().is_none_or(|w| moved(w, &place)) {
            keeper.window = Some(place);
            keeper.window_changed_at = Some(now);
        }
        let waiting = keeper
            .window_saved
            .as_ref()
            .is_none_or(|w| keeper.window.as_ref().is_some_and(|c| moved(w, c)));
        if let (true, Some(at), Some(current), Some(store)) = (
            waiting,
            keeper.window_changed_at,
            keeper.window,
            keeper.store.as_ref(),
        ) {
            if now.duration_since(at) >= WINDOW_SAVE_AFTER {
                match store.save_window(&current) {
                    Ok(()) => keeper.window_saved = Some(current),
                    Err(why) => keeper.problem = Some(why),
                }
                keeper.window_changed_at = None;
            } else {
                ctx.request_repaint_after(WINDOW_SAVE_AFTER);
            }
        }
    }

    /// A window that was left on a screen that is not there any more starts where no screen shows it: once, at the start, when the
    /// size of the screen is known, it is moved back to a corner of it.
    fn bring_window_back(&mut self, ctx: &egui::Context) {
        if self.keeper.on_screen_checked || self.keeper.store.is_none() {
            return;
        }
        let looked = ctx.input(|i| {
            let v = i.viewport();
            Some((v.monitor_size?, v.outer_rect?))
        });
        let Some((screen, outer)) = looked else {
            return;
        };
        self.keeper.on_screen_checked = true;
        let visible = egui::Rect::from_min_max(egui::Pos2::ZERO, screen.to_pos2()).intersect(outer);
        if visible.width() < LEAST_VISIBLE || visible.height() < LEAST_VISIBLE {
            ctx.send_viewport_cmd(egui::ViewportCommand::OuterPosition(BACK_ON_SCREEN));
        }
    }

    /// Put the sequences of `session` back, each as it was (with the changes that were not saved), and the tab and workspace that were
    /// showing. A sequence whose file is open already is not opened again (one with unsaved changes comes back as a copy of its own,
    /// not tied to the file, so that saving it cannot write over what is open).
    pub fn restore_session(&mut self, session: &Session) -> Restored {
        let mut done = Restored {
            lost: session.lost,
            ..Restored::default()
        };
        let mut open: Vec<PathBuf> = self
            .docs
            .iter()
            .filter_map(|d| d.doc.path.clone())
            .collect();
        let mut back: Vec<Doc> = Vec::new();
        // where each saved sequence went among those that came back
        let mut went: Vec<Option<usize>> = Vec::new();
        for saved in &session.docs {
            let Some(mut doc) = read_saved(saved) else {
                done.lost += 1;
                went.push(None);
                continue;
            };
            doc.name = saved.name.clone();
            let already = saved.path.as_ref().is_some_and(|p| open.contains(p));
            if already && !saved.dirty {
                done.already_open += 1;
                went.push(None);
                continue;
            }
            if already {
                doc.name = format!("{} (recovered)", saved.name);
                doc.path = None;
            } else {
                doc.path = saved.path.clone();
                open.extend(saved.path.clone());
            }
            doc.dirty = saved.dirty;
            went.push(Some(back.len()));
            back.push(doc);
        }
        done.opened = back.len();
        for kept in &session.results {
            let put = kept
                .text
                .as_deref()
                .map(|text| self.put_result_back(&kept.kind, text));
            match put {
                Some(Ok(())) => done.results.push(kept.kind.clone()),
                _ => done.results_lost.push(kept.kind.clone()),
            }
        }
        if back.is_empty() && done.results.is_empty() {
            return done;
        }
        self.workspace = Workspace::from_view(&session.workspace).unwrap_or(Workspace::Sequences);
        if back.is_empty() {
            return done;
        }
        // the example plasmid that the program opened with is not wanted beside what the person was working on; and what points at a
        // tab by its number points at nothing now
        self.docs.retain(|d| !is_pristine(d));
        self.confirm_close = None;
        self.last_doc = 0;
        let first = self.docs.len();
        // the tab that was in front, or the next that came back
        let front = went
            .iter()
            .skip(session.active)
            .flatten()
            .next()
            .or_else(|| went.iter().flatten().last())
            .copied()
            .unwrap_or(0);
        self.docs.extend(back.into_iter().map(DocView::new));
        self.tab = Tab::Doc(first + front);
        done
    }

    /// What was said to the offer of the last session: reopen it, or start with nothing.
    fn answer_offer(&mut self, reopen: bool) {
        let Some(session) = self.keeper.offer.take() else {
            return;
        };
        self.keeper.confirm_fresh = false;
        if reopen {
            let done = self.restore_session(&session);
            self.say(reopened(&session, &done));
        } else if !self.keeper.still_keeps() {
            // the session is another page's now: what was offered here is not this page's to remove
        } else if let Some(store) = &self.keeper.store {
            // what was offered is in the session folder, or was put aside from it
            let result = if self.keeper.archived {
                store.clear_previous().map_err(|e| e.to_string())
            } else {
                store.clear_session()
            };
            self.keeper.archived = false;
            if let Err(why) = result {
                self.keeper.problem = Some(why);
            }
        }
    }

    /// The question that the last session is asked with, until it is answered: nothing else in the window can be used meanwhile.
    pub(super) fn session_offer_window(&mut self, ctx: &egui::Context) {
        let Some(session) = &self.keeper.offer else {
            return;
        };
        let count = session.docs.len();
        let unsaved = session.docs.iter().filter(|d| d.dirty).count();
        let named: Vec<String> = session
            .docs
            .iter()
            .take(MOST_NAMED)
            .map(|d| {
                format!(
                    "{}{}",
                    d.name,
                    if d.dirty { "  (changes not saved)" } else { "" }
                )
            })
            .collect();
        let more = count.saturating_sub(MOST_NAMED);
        let lost = session.lost;
        let results = results::named(session.results.iter().map(|r| r.kind.as_str()));
        let confirming = self.keeper.confirm_fresh;
        let (mut reopen, mut fresh, mut sure, mut back) = (false, false, false, false);
        let look = self.prefs.look();
        egui::Modal::new(Id::new("session offer")).show(ctx, |ui| {
            ui.set_max_width(420.0);
            if confirming {
                ui.heading("Forget the unsaved work?");
                ui.label(format!(
                    "{unsaved} sequence{} had changes that were not saved. They are only in the last session; starting fresh removes them.",
                    if unsaved == 1 { "" } else { "s" }
                ));
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    sure = ui.button("Forget them").clicked();
                    back = ui.button("Go back").clicked();
                });
                return;
            }
            ui.heading("Reopen your last session?");
            if count == 0 {
                let none = if lost > 0 {
                    "The sequences that were open last time cannot be read any more: their files have gone."
                } else {
                    "No sequence was open last time."
                };
                if results.is_empty() {
                    ui.label(none);
                    ui.add_space(6.0);
                    fresh = ui.button("OK").clicked();
                    return;
                }
                // what the pages showed is still something to come back to
                ui.label(format!("{none} Kept from it: {results}."));
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    reopen = ui.button("Reopen").clicked();
                    fresh = ui.button("Start fresh").clicked();
                });
                return;
            }
            ui.label(format!(
                "Last time {count} sequence{} {} open{}.",
                if count == 1 { "" } else { "s" },
                if count == 1 { "was" } else { "were" },
                if unsaved > 0 {
                    format!(", {unsaved} with changes that were not saved")
                } else {
                    String::new()
                }
            ));
            for name in &named {
                ui.label(RichText::new(format!("•  {name}")).monospace());
            }
            if more > 0 {
                ui.label(RichText::new(format!("and {more} more")).color(look.dim));
            }
            if lost > 0 {
                ui.label(
                    RichText::new(format!(
                        "{lost} more cannot be reopened: {} cannot be read.",
                        if lost == 1 { "it" } else { "they" }
                    ))
                        .color(look.warn),
                );
            }
            if !results.is_empty() {
                ui.label(format!("Also kept: {results}."));
            }
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                reopen = ui.button("Reopen").clicked();
                fresh = ui.button("Start fresh").clicked();
            });
        });
        if fresh && unsaved > 0 {
            // forgetting work that was never saved is asked about once more
            self.keeper.confirm_fresh = true;
        } else if reopen {
            self.answer_offer(true);
        } else if fresh || sure {
            self.answer_offer(false);
        } else if back {
            self.keeper.confirm_fresh = false;
        }
    }

    /// File ▸ Open Recent: open the `i`th recent file.
    pub(super) fn open_recent(&mut self, i: usize) {
        let Some(path) = self.keeper.recent.get(i).cloned() else {
            return;
        };
        if path.is_file() {
            self.open_file(&path);
        } else {
            self.keeper.recent.retain(|p| *p != path);
            self.say(format!("{} is not there any more", path.display()));
        }
    }

    /// The submenu of File ▸ Open Recent: what it lists, as (label, action).
    pub(super) fn recent_items(&self) -> Vec<(String, Action)> {
        self.keeper
            .recent
            .iter()
            .enumerate()
            .map(|(i, p)| (recent_label(p), Action::OpenRecent(i)))
            .collect()
    }
}

/// A recent file as the menu shows it: its name, and its folder when that helps tell two of the same name apart.
pub fn recent_label(path: &Path) -> String {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("file");
    match path
        .parent()
        .and_then(|p| p.file_name())
        .and_then(|n| n.to_str())
    {
        Some(folder) => format!("{name}   —   {folder}"),
        None => name.to_string(),
    }
}

/// Whether nobody has touched this sequence and it is not worth keeping: it is the example plasmid that the program opens with, or
/// it is a new one that is still empty.
fn is_pristine(d: &DocView) -> bool {
    d.rev == 0
        && !d.doc.dirty
        && d.doc.path.is_none()
        && (d.doc.name == crate::sample::DEMO_NAME || d.doc.is_empty())
}

/// A sequence of the session as a document: from its own text, or from its file, taking the record of its name (the only record of a
/// file that has one). `None` when it cannot be read, or its record is not in the file.
fn read_saved(saved: &SavedDoc) -> Option<Doc> {
    let docs = match (&saved.text, &saved.path) {
        (Some(text), _) => io::parse_any(text, &saved.name).ok()?,
        (None, Some(path)) => {
            let text = std::fs::read(path).ok()?;
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("Untitled");
            io::parse_any(&String::from_utf8_lossy(&text), stem).ok()?
        }
        (None, None) => return None,
    };
    let at = if docs.len() == 1 {
        0
    } else {
        docs.iter().position(|d| d.name == saved.name)?
    };
    docs.into_iter().nth(at)
}

#[cfg(test)]
mod tests;
