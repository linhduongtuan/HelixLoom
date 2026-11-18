//! Keeping the person's work and place between runs, in the program's settings folder: the **session** (the sequences that are open,
//! saved as they change, so that a crash or a power cut loses at most a few seconds), the **recent files**, and **where the window was**.
//!
//! The session is a folder: `manifest.txt`, which says what is open, and one GenBank file for each sequence that cannot be read back
//! from its own file (one that has no file, or has unsaved changes). A save writes new files under new names, makes them durable
//! (`fsync`), replaces the manifest in one step, and then removes the old files, so that whatever moment the program stops at, even a
//! power cut, the folder holds a whole session: the one before or the one after. The folder and its files are private to the person
//! (mode 0700 and 0600 on Unix), and only one copy of the program at a time keeps a session in it (a lock file).
//!
//! In a browser there is no settings folder, and the same things are kept in the storage the browser gives the page
//! (`session/browser.rs`), through the same code: what is written here is written to a `Medium`, which is the folder or that
//! storage (`session/medium.rs`). There the texts of the sequences are kept first and the list of what is open last, as here,
//! so that a page that is closed in the middle of a save comes back with the session before it.

use super::appearance::Appearance;
use super::places::Places;
use super::state::ThemeMode;
use crate::kept::{escape, unescape};
use crate::seq::doc::Doc;
use crate::seq::io;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

mod browser;
#[cfg(any(target_arch = "wasm32", test))]
pub(in crate::ui) mod large;
mod medium;

/// The table a browser gives one tab, which is gone when the tab is closed (`browser.rs`).
#[cfg(target_arch = "wasm32")]
pub use browser::{open_database, tab_shelf, wake_on_refusal};
#[cfg(test)]
use medium::write_atomic;
use medium::Medium;
#[cfg(test)]
pub use medium::Memory;
pub use medium::Shelf;

const SESSION: &str = "session";
const MANIFEST: &str = "manifest.txt";
const RECENT: &str = "recent.txt";
const WINDOW: &str = "window.txt";
const PLACES: &str = "places.txt";
/// How much of the file of the places is read: far more than the most places a list holds.
const MAX_PLACES_BYTES: u64 = 32 * 1024;
const APPEARANCE: &str = "appearance.txt";
/// The most that is read of `appearance.txt`: it holds six short lines.
const MAX_APPEARANCE_BYTES: u64 = 4096;
/// Where a session waiting on offer is put when the person goes on working without answering, so that the new work can be kept without
/// writing over it; and where one that this version cannot read is put before it is replaced.
const PREVIOUS: &str = "session.previous";
/// How many times removing the put-aside session is tried, and how long to wait between tries.
const CLEAR_TRIES: usize = 3;
const CLEAR_PAUSE: std::time::Duration = std::time::Duration::from_millis(25);
const UNREADABLE: &str = "session.unreadable";
const LOCK: &str = "session.lock";
/// In a browser, the mark of the page that keeps the session (see `Lock::Taken`).
const OWNER: &str = "session.owner";
/// The mark that says the person asked for the open sequences not to be kept.
const NOT_WANTED: &str = "session.off";
/// The first line of a manifest: what it is, and its version.
const MAGIC: &str = "helixloom-session\t1";
/// The ending of a file that holds what a page showed.
const RESULT: &str = ".result";
/// How many recent files are kept.
pub const MOST_RECENT: usize = 12;

/// What taking the session lock came to.
#[derive(Debug)]
pub enum Lock {
    /// This copy has it (the file holds it until the program ends).
    Held(fs::File),
    /// Another copy of the program has it.
    Busy,
    /// The system could not lock a file here: go on without.
    Unavailable,
    /// In a browser, where there is no file to lock and no telling whether another page that is open on the site is still
    /// there: this page has taken the session over, under this mark. A page that was open before it finds another mark than
    /// its own when it next wants to write, and stops keeping the session (`Store::owns`), so that two pages never write
    /// into one session.
    Taken(String),
}

/// Where the settings are kept, and what is read and written there.
#[derive(Clone, Debug)]
pub struct Store {
    medium: Medium,
}

/// The name of a file of the session folder, as a medium takes it.
fn in_session(file: &str) -> String {
    format!("{SESSION}/{file}")
}

/// One open sequence as the session keeps it.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct SavedDoc {
    pub name: String,
    /// The file it was opened from or saved to.
    pub path: Option<PathBuf>,
    /// It has changes that are not in its file.
    pub dirty: bool,
    /// The sequence itself (GenBank text), kept when it cannot be read back from `path`: it has none, or has unsaved changes.
    pub text: Option<String>,
    /// The file of the last save that holds this sequence's text, when it has not changed since: the save leaves it where it is
    /// instead of writing the text again.
    pub keep: Option<String>,
}

impl SavedDoc {
    /// A sequence as the session keeps it: with its text unless its own file holds it.
    pub fn of(doc: &Doc) -> SavedDoc {
        let mut saved = SavedDoc::without_text(doc);
        if !(doc.path.is_some() && !doc.dirty) {
            saved.text = Some(io::to_genbank(doc));
        }
        saved
    }

    /// The same without its text (for one whose text is already in a file of the session).
    pub fn without_text(doc: &Doc) -> SavedDoc {
        SavedDoc {
            name: doc.name.clone(),
            path: doc.path.clone(),
            dirty: doc.dirty,
            text: None,
            keep: None,
        }
    }

    /// Whether the session must hold this sequence's text (its own file does not).
    pub fn needs_text(doc: &Doc) -> bool {
        !(doc.path.is_some() && !doc.dirty)
    }
}

/// What a page shows, as the session keeps it: an alignment, a tree or a structure, in the text its page gives (`ui::results`).
#[derive(Clone, Debug, PartialEq, Default)]
pub struct SavedResult {
    /// The page it is of: a word of small letters.
    pub kind: String,
    /// The result itself; none when the file of the last save holds it as it is (`keep`).
    pub text: Option<String>,
    pub keep: Option<String>,
}

/// The files that a save put the texts of a session in.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Kept {
    /// For each sequence, the file that holds its text (none for one that is read back from its own file).
    pub docs: Vec<Option<String>>,
    /// For each result, the file that holds it; none for one that was left out, because there was no room for it or it could not be
    /// written. A result is kept when it can be: one that cannot never stands in the way of the sequences.
    pub results: Vec<Option<String>>,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct Session {
    /// The workspace showing (its name, as `Workspace::name` gives it).
    pub workspace: String,
    /// The tab that was in front: the number of a sequence.
    pub active: usize,
    pub docs: Vec<SavedDoc>,
    /// How many sequences the folder named and could not give back (a file that is gone).
    pub lost: usize,
    /// What the pages showed.
    pub results: Vec<SavedResult>,
}

/// Where the window was and how big.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindowPlace {
    pub x: f32,
    pub y: f32,
    pub width: f32,
    pub height: f32,
}

impl WindowPlace {
    /// The same place with every number multiplied by `factor`. The rectangles that egui reports for the window are in points, which
    /// are the window's own units divided by the zoom (see `appearance`); the window is made at the size that is kept as it stands, with
    /// no zoom yet, so what is kept has to be in the window's own units: the places go in at the zoom they were seen at.
    pub fn scaled(self, factor: f32) -> WindowPlace {
        WindowPlace {
            x: self.x * factor,
            y: self.y * factor,
            width: self.width * factor,
            height: self.height * factor,
        }
    }

    /// Whether this is a place a window could be: a size that fits a screen, a corner that is not far away (a monitor that has gone leaves a
    /// place that no window should be put at, and these limits cannot tell that; they only keep out numbers that are nonsense).
    pub fn sane(&self) -> bool {
        let ok = |v: f32, least: f32, most: f32| v.is_finite() && (least..=most).contains(&v);
        ok(self.width, 400.0, 12_000.0)
            && ok(self.height, 300.0, 12_000.0)
            && ok(self.x, -2_000.0, 20_000.0)
            && ok(self.y, -200.0, 20_000.0)
    }
}

/// What the list says of a result, to tell that the text read back is the one that was written: its length and its checksum.
fn check_of(text: &str) -> String {
    format!("{}\t{:016x}", text.len(), crate::kept::hash(text))
}

/// Whether `name` is the name of a file in a folder, and not a way to somewhere else.
fn plain_name(name: &str) -> bool {
    !name.is_empty() && !name.starts_with('.') && !name.contains(['/', '\\', '\0', ':'])
}

/// Whether `kind` can name a result: a word of small letters, which is all a page calls its own.
fn is_kind(kind: &str) -> bool {
    !kind.is_empty() && kind.bytes().all(|c| c.is_ascii_lowercase())
}

/// The kind of the result that the file `name` of the session folder holds (`generation-kind.result`), if it is one.
fn result_kind(name: &str) -> Option<&str> {
    let (generation, kind) = name.strip_suffix(RESULT)?.split_once('-')?;
    let numbered = !generation.is_empty() && generation.bytes().all(|c| c.is_ascii_digit());
    (numbered && is_kind(kind)).then_some(kind)
}

/// Whether `name` is a file this program writes in the session folder: a sequence's `generation-number.gb`, a result's
/// `generation-kind.result`, or a temporary file.
fn is_ours(name: &str) -> bool {
    let sequence = name.strip_suffix(".gb").is_some_and(|stem| {
        stem.split_once('-').is_some_and(|(a, b)| {
            !a.is_empty()
                && !b.is_empty()
                && a.bytes().all(|c| c.is_ascii_digit())
                && b.bytes().all(|c| c.is_ascii_digit())
        })
    });
    sequence || result_kind(name).is_some() || (name.starts_with('.') && name.ends_with(".tmp"))
}

/// `path` put first in `list`, once, with at most `most` kept.
pub fn push_recent(list: &[PathBuf], path: &Path, most: usize) -> Vec<PathBuf> {
    std::iter::once(path.to_path_buf())
        .chain(list.iter().filter(|p| p.as_path() != path).cloned())
        .take(most)
        .collect()
}

impl Store {
    /// The store in `dir`. With no settings folder it is the browser's storage when the program runs in a page, and none
    /// anywhere else (a picture is being made, say): nothing is then read or kept.
    pub fn open(dir: Option<&Path>) -> Option<Store> {
        match dir {
            Some(dir) => Some(Store {
                medium: Medium::Folder(dir.to_path_buf()),
            }),
            None => Store::of_the_page(),
        }
    }

    #[cfg(target_arch = "wasm32")]
    fn of_the_page() -> Option<Store> {
        browser::shelf().map(Store::on)
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn of_the_page() -> Option<Store> {
        None
    }

    /// A store that keeps things in a table: the browser's, or a test's.
    #[cfg(any(target_arch = "wasm32", test))]
    pub(in crate::ui) fn on(shelf: std::rc::Rc<dyn Shelf>) -> Store {
        Store {
            medium: Medium::Shelf(shelf),
        }
    }

    /// Whether what is kept is kept by a browser, for this page's site, and not in a folder of the computer: the window says
    /// where its settings are.
    pub fn in_a_browser(&self) -> bool {
        matches!(self.medium, Medium::Shelf(_))
    }

    // ---- marks ---------------------------------------------------------------------------------------------------------------------

    /// Whether the mark `name` was kept: something that is done once, such as the tour of the first run.
    pub fn has_mark(&self, name: &str) -> bool {
        self.medium.holds(name)
    }

    pub fn keep_mark(&self, name: &str) -> Result<(), String> {
        self.medium.write(name, "seen\n", false)
    }

    /// Why a write that was not waited for did not happen (the page's larger table), once; none otherwise.
    pub fn shelf_problem(&self) -> Option<String> {
        match &self.medium {
            Medium::Shelf(shelf) => shelf.take_problem(),
            Medium::Folder(_) => None,
        }
    }

    /// Whether the open sequences are to be kept from one run to the next: they are, until the person says not to (on a
    /// computer that others use, say).
    pub fn sessions_wanted(&self) -> bool {
        !self.has_mark(NOT_WANTED)
    }

    /// Remove everything that is kept of the open sequences: the session, one that was put aside, and one of another
    /// version of the program that could not be read and was moved out of the way. Each is tried, whatever the others say.
    pub fn forget_sessions(&self) -> Result<(), String> {
        let removed = [SESSION, PREVIOUS, UNREADABLE].map(|folder| {
            self.remove_kept(folder)
                .map_err(|e| format!("Could not remove {folder}: {e}"))
        });
        // whatever another page of the site put in the larger table is taken away too, by what the page does not know the names of
        self.medium.clear_apart();
        removed.into_iter().collect()
    }

    /// Remove a folder of the settings, trying again for a moment when it will not go (a program that looks through the
    /// files of the computer may hold one of them open); one that is not there is removed.
    fn remove_kept(&self, folder: &str) -> std::io::Result<()> {
        let mut last = Ok(());
        for attempt in 0..CLEAR_TRIES {
            match self.medium.remove_folder(folder) {
                Ok(()) => return Ok(()),
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(()),
                Err(e) => last = Err(e),
            }
            if attempt + 1 < CLEAR_TRIES {
                std::thread::sleep(CLEAR_PAUSE);
            }
        }
        last
    }

    /// Keep the person's answer to that, for this run and the next.
    pub fn want_sessions(&self, wanted: bool) -> Result<(), String> {
        if !wanted {
            return self.keep_mark(NOT_WANTED);
        }
        self.medium.remove(NOT_WANTED);
        if self.has_mark(NOT_WANTED) {
            return Err("Could not turn the keeping on again: the mark that turns it off cannot be removed.".into());
        }
        Ok(())
    }

    // ---- recent files --------------------------------------------------------------------------------------------------------------

    /// The recent files, the latest first; those that are no longer there are left out.
    pub fn recent(&self) -> Vec<PathBuf> {
        let text = self.medium.read(RECENT, None).unwrap_or_default();
        text.lines()
            .map(|l| PathBuf::from(unescape(l)))
            .filter(|p| p.is_file())
            .take(MOST_RECENT)
            .collect()
    }

    /// Put `path` at the top of the recent files.
    pub fn add_recent(&self, path: &Path) -> Result<(), String> {
        let list = push_recent(&self.recent(), path, MOST_RECENT);
        self.write_recent(&list)
    }

    pub fn clear_recent(&self) -> Result<(), String> {
        self.write_recent(&[])
    }

    fn write_recent(&self, list: &[PathBuf]) -> Result<(), String> {
        let text: String = list
            .iter()
            .map(|p| format!("{}\n", escape(&p.to_string_lossy())))
            .collect();
        self.medium.write(RECENT, &text, false)
    }

    // ---- the window ----------------------------------------------------------------------------------------------------------------

    /// Where the window was left, if that was kept and makes sense.
    pub fn window(&self) -> Option<WindowPlace> {
        let text = self.medium.read(WINDOW, None)?;
        let numbers: Vec<f32> = text
            .split_whitespace()
            .filter_map(|n| n.parse().ok())
            .collect();
        let [x, y, width, height] = numbers.as_slice() else {
            return None;
        };
        let place = WindowPlace {
            x: *x,
            y: *y,
            width: *width,
            height: *height,
        };
        place.sane().then_some(place)
    }

    pub fn save_window(&self, place: &WindowPlace) -> Result<(), String> {
        let text = format!("{} {} {} {}\n", place.x, place.y, place.width, place.height);
        self.medium.write(WINDOW, &text, false)
    }

    // ---- the places of the tool windows ---------------------------------------------------------------------------------------------

    /// Where the tool windows were left (none kept, or a file that is no file, gives none).
    pub fn places(&self) -> Places {
        match self.medium.read(PLACES, Some(MAX_PLACES_BYTES)) {
            Some(text) => Places::parse(&text),
            None => Places::default(),
        }
    }

    pub fn save_places(&self, places: &Places) -> Result<(), String> {
        self.medium.write(PLACES, &places.to_text(), false)
    }

    // ---- the appearance ------------------------------------------------------------------------------------------------------------

    /// The theme and the appearance settings that were kept, none when nothing was.
    pub fn appearance(&self) -> Option<(ThemeMode, Appearance)> {
        // what this program writes is a few lines of text in a file: not a pipe or a device that something else put there (reading
        // one could wait for ever), and no more than a few thousand bytes of it are read, whatever its size
        let text = self.medium.read(APPEARANCE, Some(MAX_APPEARANCE_BYTES))?;
        Some(Appearance::decode(&text))
    }

    pub fn save_appearance(&self, theme: ThemeMode, appearance: &Appearance) -> Result<(), String> {
        self.medium
            .write(APPEARANCE, &appearance.encode(theme), false)
    }

    // ---- the session ---------------------------------------------------------------------------------------------------------------

    /// Take the lock that says this copy of the program keeps the session: `Some` (hold it as long as the program runs) when nobody else
    /// has it, none when another copy does. Where the system cannot lock a file the lock is taken as given: a session that is
    /// shared is better than none.
    pub fn lock_session(&self) -> Lock {
        let Medium::Folder(dir) = &self.medium else {
            return match self.take_over() {
                Some(mark) => Lock::Taken(mark),
                // the mark could not be written (the storage is full, or switched off): nothing else will be either
                None => Lock::Unavailable,
            };
        };
        if fs::create_dir_all(dir).is_err() {
            return Lock::Unavailable;
        }
        let opened = fs::OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(dir.join(LOCK));
        let Ok(file) = opened else {
            return Lock::Unavailable;
        };
        match file.try_lock() {
            Ok(()) => Lock::Held(file),
            Err(fs::TryLockError::WouldBlock) => Lock::Busy,
            Err(fs::TryLockError::Error(_)) => Lock::Unavailable,
        }
    }

    /// Write this page's mark over the one that is there: a number greater than any before it (the time, or one more than
    /// the last mark when two pages are opened in the same thousandth of a second).
    fn take_over(&self) -> Option<String> {
        let last: u64 = self
            .medium
            .read(OWNER, None)
            .and_then(|kept| kept.trim().parse().ok())
            .unwrap_or(0);
        let now = web_time::SystemTime::now()
            .duration_since(web_time::UNIX_EPOCH)
            .map_or(0, |since| since.as_millis() as u64);
        let mark = now.max(last.saturating_add(1)).to_string();
        self.medium.write(OWNER, &mark, false).ok()?;
        Some(mark)
    }

    /// Whether the session is still this copy's to write: always in a folder, where the lock file says so for as long as the
    /// program runs; in a browser, as long as the mark that is kept is the one this page took the session over with (`mark`;
    /// a page that could write no mark has none to compare).
    pub fn owns(&self, mark: Option<&str>) -> bool {
        match (&self.medium, mark) {
            (Medium::Folder(_), _) | (Medium::Shelf(_), None) => true,
            (Medium::Shelf(_), Some(mark)) => self
                .medium
                .read(OWNER, None)
                .is_some_and(|kept| kept.trim() == mark),
        }
    }

    /// Keep `session`: the sequences, and what the pages show when there is room for it. With neither there is nothing to come back
    /// to, and the folder is removed. The sequences come first: when a text of theirs or the list cannot be written while results
    /// are kept, there may be no room, so the results that are written where the texts are are removed and the save is tried once
    /// more without them (a result that is kept apart, in a page's larger table, stays and is written). A save that fails for
    /// another reason leaves the results where they are: making room would not mend it.
    pub fn keep_session(&self, session: &Session) -> Result<Kept, String> {
        if session.docs.is_empty() && session.results.is_empty() {
            self.clear_session()?;
            return Ok(Kept::default());
        }
        // the results that make room when a text cannot be written: those that are written where the texts are
        let results = |store: &Store| -> Vec<String> {
            let names = store.medium.names_in(SESSION).into_iter();
            let shared = |name: &String| store.medium.shares_room(&in_session(name));
            names
                .filter(|name| result_kind(name).is_some() && shared(name))
                .collect()
        };
        let done = match self.write_session(session, true) {
            Err((_, true)) if !results(self).is_empty() => {
                for name in results(self) {
                    self.medium.remove(&in_session(&name));
                }
                self.write_session(session, false)
            }
            done => done,
        };
        done.map_err(|(why, _)| why)
    }

    /// The same, giving only the files of the sequences.
    #[cfg(test)]
    pub fn save_session(&self, session: &Session) -> Result<Vec<Option<String>>, String> {
        self.keep_session(session).map(|kept| kept.docs)
    }

    /// Write the session; when it cannot be written, why, and whether it was a text that could not be written (room may help).
    /// `shared_too`: the results that are written where the texts are are written too (not in the retry that makes room).
    fn write_session(&self, session: &Session, shared_too: bool) -> Result<Kept, (String, bool)> {
        self.set_aside_unreadable().map_err(|why| (why, false))?;
        self.medium
            .make_folder(SESSION)
            .map_err(|why| (why, false))?;
        let generation = self.generation() + 1;
        let mut manifest = format!(
            "{MAGIC}\ngeneration\t{generation}\nworkspace\t{}\nactive\t{}\n",
            escape(&session.workspace),
            session.active
        );
        let (mut kept, mut files) = (Vec::new(), Kept::default());
        for (n, doc) in session.docs.iter().enumerate() {
            let file = match (&doc.text, &doc.keep) {
                (Some(text), _) => {
                    let file = format!("{generation}-{n}.gb");
                    self.medium
                        .write(&in_session(&file), text, true)
                        .map_err(|why| (why, true))?;
                    Some(file)
                }
                // unchanged since the last save: its file is kept as it is, if it is still there
                (None, Some(file)) if plain_name(file) && self.medium.holds(&in_session(file)) => {
                    Some(file.clone())
                }
                (None, Some(file)) => {
                    let why = format!("{file}, kept from the last save, is not there");
                    return Err((why, false));
                }
                _ => None,
            };
            if let Some(file) = &file {
                kept.push(file.clone());
            }
            let path = doc
                .path
                .as_ref()
                .map_or(String::new(), |p| p.to_string_lossy().into_owned());
            manifest.push_str(&format!(
                "doc\t{}\t{}\t{}\t{}\n",
                u8::from(doc.dirty),
                escape(&doc.name),
                escape(&path),
                file.as_deref().unwrap_or("")
            ));
            files.docs.push(file);
        }
        let before = self.checks_of_the_list();
        for result in &session.results {
            let wanted = |file: &str| shared_too || !self.medium.shares_room(&in_session(file));
            let new_file = format!("{generation}-{}{RESULT}", result.kind);
            let (file, check) = match (&result.text, &result.keep) {
                _ if !is_kind(&result.kind) => (None, None),
                (Some(text), _) if wanted(&new_file) => (
                    self.write_result(generation, &result.kind, text),
                    Some(check_of(text)),
                ),
                // unchanged since the last save: what that list said of it is said again
                (None, Some(file)) if wanted(file) => {
                    let there = result_kind(file) == Some(result.kind.as_str())
                        && self.medium.holds(&in_session(file));
                    (there.then(|| file.clone()), before.get(file).cloned())
                }
                _ => (None, None),
            };
            if let Some(file) = &file {
                kept.push(file.clone());
                let check = check.map_or(String::new(), |check| format!("\t{check}"));
                manifest.push_str(&format!("result\t{}\t{file}{check}\n", result.kind));
            }
            files.results.push(file);
        }
        // the manifest is what makes the new files the session: only now
        self.medium
            .write(&in_session(MANIFEST), &manifest, true)
            .map_err(|why| (why, true))?;
        // and what the old manifest named is of no use now (only files this program wrote are removed)
        for name in self.medium.names_in(SESSION) {
            if is_ours(&name) && !kept.contains(&name) {
                self.medium.remove(&in_session(&name));
            }
        }
        Ok(files)
    }

    /// Write what a page shows, and give the file it is in; none when it cannot be written. The result of the last save is still there
    /// then, and its page no longer shows it: where the two do not fit together, the old one makes room and the new one is tried
    /// once more. (A sequence is not treated so: until the new manifest is written, the old texts are the session.)
    fn write_result(&self, generation: u64, kind: &str, text: &str) -> Option<String> {
        let file = format!("{generation}-{kind}{RESULT}");
        let write = || self.medium.write(&in_session(&file), text, true).is_ok();
        if write() {
            return Some(file);
        }
        for name in self.medium.names_in(SESSION) {
            if name != file && result_kind(&name) == Some(kind) {
                self.medium.remove(&in_session(&name));
            }
        }
        write().then_some(file)
    }

    /// What the list that is there says of each of its results: the length and the checksum of the file, by the file's name.
    /// A result that is kept as it is gets the same said of it in the next list, without being read again.
    fn checks_of_the_list(&self) -> HashMap<String, String> {
        let Some(text) = self.medium.read(&in_session(MANIFEST), None) else {
            return HashMap::new();
        };
        text.lines()
            .filter_map(
                |line| match line.split('\t').collect::<Vec<_>>().as_slice() {
                    ["result", _, file, len, sum] => {
                        Some((file.to_string(), format!("{len}\t{sum}")))
                    }
                    _ => None,
                },
            )
            .collect()
    }

    /// The generation of the manifest that is there (0 when there is none).
    fn generation(&self) -> u64 {
        self.medium
            .read(&in_session(MANIFEST), None)
            .and_then(|text| {
                text.lines().find_map(|l| {
                    l.strip_prefix("generation\t")
                        .and_then(|n| n.trim().parse().ok())
                })
            })
            .unwrap_or(0)
    }

    /// A session of another version of the program, or not a session at all, is not written over: its folder is put aside, once.
    /// When it cannot be put aside nothing is saved into it: its files may be named as this program names its own.
    fn set_aside_unreadable(&self) -> Result<(), String> {
        let Some(text) = self.medium.read(&in_session(MANIFEST), None) else {
            return Ok(());
        };
        if text.lines().next() == Some(MAGIC) {
            return Ok(());
        }
        self.medium
            .move_folder(SESSION, UNREADABLE)
            .map_err(|e| format!("a session of another version could not be put aside: {e}"))
    }

    /// Put the session that is there aside, in `session.previous`, so that a new session can be kept without writing over it (the
    /// person has gone on working before answering the offer of the last one).
    pub fn archive_session(&self) -> Result<(), String> {
        if !self.medium.has_folder(SESSION) {
            return Ok(());
        }
        self.medium
            .move_folder(SESSION, PREVIOUS)
            .map_err(|e| format!("Could not put the last session aside: {e}"))
    }

    /// Remove the session that was put aside. A folder may fail to go for a moment (something else, such as the system's indexer, has
    /// a file of it open, and "not empty" is said), so it is tried a few times; a folder that is not there is removed.
    pub fn clear_previous(&self) -> std::io::Result<()> {
        self.remove_kept(PREVIOUS)
    }

    /// The session that was kept, if there is one with something in it. A sequence whose file is gone is counted in `lost`.
    pub fn load_session(&self) -> Option<Session> {
        let text = self.medium.read(&in_session(MANIFEST), None)?;
        let mut lines = text.lines();
        if lines.next()? != MAGIC {
            return None;
        }
        let mut session = Session::default();
        for line in lines {
            let fields: Vec<&str> = line.split('\t').collect();
            match fields.as_slice() {
                ["workspace", name] => session.workspace = unescape(name),
                ["active", n] => session.active = n.parse().unwrap_or(0),
                ["doc", dirty, name, path, file] => {
                    let path = (!path.is_empty()).then(|| PathBuf::from(unescape(path)));
                    let text = if file.is_empty() {
                        None
                    } else if !plain_name(file) {
                        // a name that leads out of the folder is not one this program wrote
                        session.lost += 1;
                        continue;
                    } else {
                        match self.medium.read(&in_session(file), None) {
                            Some(text) => Some(text),
                            None => {
                                session.lost += 1;
                                continue;
                            }
                        }
                    };
                    session.docs.push(SavedDoc {
                        name: unescape(name),
                        path,
                        dirty: *dirty == "1",
                        text,
                        keep: None,
                    });
                }
                // a result that cannot be read back is left out: nothing depends on it. Nor is one that is not the text that was
                // listed (a list of an earlier version has no length and checksum to tell by)
                ["result", kind, file, check @ ..]
                    if result_kind(file) == Some(*kind) && matches!(check.len(), 0 | 2) =>
                {
                    let text = self.medium.read(&in_session(file), None);
                    let listed =
                        |text: &str| check.is_empty() || check.join("\t") == check_of(text);
                    if let Some(text) = text.filter(|text| listed(text)) {
                        session.results.push(SavedResult {
                            kind: kind.to_string(),
                            text: Some(text),
                            keep: None,
                        });
                    }
                }
                _ => {}
            }
        }
        let something = !session.docs.is_empty() || session.lost > 0 || !session.results.is_empty();
        something.then_some(session)
    }

    pub fn clear_session(&self) -> Result<(), String> {
        match self.medium.remove_folder(SESSION) {
            Ok(()) => Ok(()),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
            Err(e) => Err(format!("Could not remove the saved session: {e}")),
        }
    }
}

#[cfg(test)]
mod browser_tests;
#[cfg(test)]
mod tests;
