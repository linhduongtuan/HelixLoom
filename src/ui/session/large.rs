//! A larger table for what the pages show, beside the small one (`browser.rs`): the browser's `localStorage` gives a site a few
//! megabytes and takes every write at once, and a database (IndexedDB) gives hundreds but answers only later, and a store reads and
//! writes synchronously (`medium.rs`). So the database is read whole before the page starts into a table in memory, which the store
//! reads from and writes to at once; every change is also queued, and the queue is handed to the database one change at a time, in
//! the order it was made. Nothing here knows IndexedDB: it is given something that applies a change (`Persist`), which in a page is
//! the database (`browser.rs`) and in a test is a stand-in.
//!
//! **What goes here and what does not.** Only the files that hold a result (`Roomy`, `is_large`). The open sequences stay in the
//! small table, because that one is written when the write is asked for, and a page that is closed at that moment keeps them; the
//! database may not have finished by then. A result is something that can be made again, so it may be lost in that moment; the
//! sequences may not. A result whose write has not finished when the page is closed is simply not there at the next visit, which
//! the session already says (`Store::load_session` leaves out a result that cannot be read back).
//!
//! **When the database says no.** A change that it refuses is not kept: what it was to put is taken out of the table, so that
//! nothing says it is there, and why is kept to be told (`take_failure`). What the database may still hold under that name (an older
//! text that the refused one was to replace, an entry that a refused removal was to take away) is taken away with the next change
//! that is asked: not at once, or a database that refuses everything would be asked for ever. And the session does not take a
//! result on trust: its list gives the length and a checksum of each result, and a text that is not the one listed is left out
//! (`session.rs`), so that a text that is left behind by a refusal can never come back as the result of another session.

use super::medium::Shelf;
use std::cell::RefCell;
use std::collections::{BTreeMap, BTreeSet};
use std::rc::Rc;

/// What is said when a text would not fit in the table in memory.
pub const FULL_LARGE: &str =
    "There is no more room for results in this page's memory: they are kept only while they fit.";

/// The most that the table in memory holds, in bytes of names and texts. The database could hold more; this is what a page can
/// hold twice, here and on the page itself, without a phone running out of memory.
pub const ROOM: usize = 160 << 20;

/// One change to the database.
#[derive(Clone, Debug, PartialEq)]
pub enum Op {
    Put(String, Rc<str>),
    Delete(String),
    /// Everything the database holds is taken away (the person said not to remember anything).
    Clear,
}

impl Op {
    /// The name the change is about; none for a change of everything.
    fn name(&self) -> Option<&str> {
        match self {
            Op::Put(name, _) | Op::Delete(name) => Some(name),
            Op::Clear => None,
        }
    }
}

/// What keeps the changes: a database in a page, a stand-in in a test.
pub trait Persist: std::fmt::Debug {
    /// Apply `op`, and call `done` when it is over (at once or later) with how it went. The next change is not handed over before
    /// `done` is called. Each change is a step of its own: one that is refused takes no other with it.
    fn apply(&self, op: Op, done: Box<dyn FnOnce(Result<(), String>)>);
}

/// What is called when the database has refused a change, to have the page look at it (the person is told on the next frame, and
/// no frame comes by itself).
struct Waker(Rc<dyn Fn()>);

impl std::fmt::Debug for Waker {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("a waker")
    }
}

#[derive(Debug)]
struct Inner {
    held: BTreeMap<String, Rc<str>>,
    bytes: usize,
    room: usize,
    /// Changes that are waiting for the one that is out; at most one for a name.
    queue: Vec<Op>,
    busy: bool,
    /// Why a change was refused, until it is taken (`Shelf::take_problem`).
    failed: Option<String>,
    /// Names that the database may hold although the table does not: taken away with the next change that is asked.
    stale: BTreeSet<String>,
    /// A change of everything was refused: it is asked again with the next change.
    clear_again: bool,
    wake: Option<Waker>,
}

impl Inner {
    /// Queue `op`. A change to a name takes the place of those that still wait for that name (what was to be put is no longer
    /// wanted, and what is to be put is put over what was there); a change of everything takes the place of all. What an earlier
    /// refusal left in the database goes first.
    fn enqueue(&mut self, op: Op) {
        if matches!(op, Op::Clear) {
            self.queue.clear();
            self.stale.clear();
            self.clear_again = false;
            self.queue.push(op);
            return;
        }
        let name = op.name().map(str::to_string).unwrap_or_default();
        self.queue.retain(|queued| queued.name() != Some(&name));
        let mut first = Vec::new();
        if std::mem::take(&mut self.clear_again) {
            self.stale.clear();
            first.push(Op::Clear);
        }
        for stale in std::mem::take(&mut self.stale) {
            let asked = stale == name || self.queue.iter().any(|q| q.name() == Some(&stale));
            if !asked {
                first.push(Op::Delete(stale));
            }
        }
        self.queue.splice(0..0, first);
        self.queue.push(op);
    }
}

/// The table in memory and the queue behind it.
#[derive(Clone, Debug)]
pub struct Large {
    inner: Rc<RefCell<Inner>>,
    persist: Rc<dyn Persist>,
}

fn size(name: &str, text: &str) -> usize {
    name.len() + text.len()
}

impl Large {
    /// A table that holds what the database held (`loaded`), whatever its size, and takes no more than `room` bytes in all.
    pub fn new(persist: Rc<dyn Persist>, loaded: Vec<(String, String)>, room: usize) -> Large {
        let held: BTreeMap<String, Rc<str>> = loaded
            .into_iter()
            .map(|(name, text)| (name, Rc::from(text)))
            .collect();
        let bytes = held.iter().map(|(name, text)| size(name, text)).sum();
        Large {
            inner: Rc::new(RefCell::new(Inner {
                held,
                bytes,
                room,
                queue: Vec::new(),
                busy: false,
                failed: None,
                stale: BTreeSet::new(),
                clear_again: false,
                wake: None,
            })),
            persist,
        }
    }

    /// Have `wake` called when the database refuses a change.
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    pub fn set_waker(&self, wake: Rc<dyn Fn()>) {
        self.inner.borrow_mut().wake = Some(Waker(wake));
    }

    pub fn get(&self, name: &str) -> Option<String> {
        self.inner
            .borrow()
            .held
            .get(name)
            .map(|text| text.to_string())
    }

    pub fn holds(&self, name: &str) -> bool {
        self.inner.borrow().held.contains_key(name)
    }

    pub fn names(&self) -> Vec<String> {
        self.inner.borrow().held.keys().cloned().collect()
    }

    /// Keep `text` under `name`, at once; and have the database keep it when it is its turn. Refused when it would not fit.
    pub fn set(&self, name: &str, text: &str) -> Result<(), String> {
        let text: Rc<str> = Rc::from(text);
        {
            let mut inner = self.inner.borrow_mut();
            let before = inner.held.get(name).map_or(0, |held| size(name, held));
            let after = inner.bytes - before + size(name, &text);
            if after > inner.room {
                return Err(FULL_LARGE.to_string());
            }
            inner.bytes = after;
            inner.held.insert(name.to_string(), text.clone());
            inner.enqueue(Op::Put(name.to_string(), text));
        }
        self.pump();
        Ok(())
    }

    /// Take `name` away at once, and from the database when it is its turn. A name that is not held asks nothing of it.
    pub fn remove(&self, name: &str) {
        {
            let mut inner = self.inner.borrow_mut();
            let Some(held) = inner.held.remove(name) else {
                return;
            };
            inner.bytes -= size(name, &held);
            inner.enqueue(Op::Delete(name.to_string()));
        }
        self.pump();
    }

    /// The text under `from` is put under `to` and `from` is taken away, with no second copy of the text in memory and the room
    /// of one text, not two. Nothing is done when there is no text under `from`; refused when `to` would not fit.
    pub fn rename(&self, from: &str, to: &str) -> Result<(), String> {
        {
            let mut inner = self.inner.borrow_mut();
            let Some(text) = inner.held.get(from).cloned() else {
                return Ok(());
            };
            let replaced = inner.held.get(to).map_or(0, |old| size(to, old));
            let after = inner.bytes - size(from, &text) - replaced + size(to, &text);
            if after > inner.room {
                return Err(FULL_LARGE.to_string());
            }
            inner.held.remove(from);
            inner.held.insert(to.to_string(), text.clone());
            inner.bytes = after;
            inner.enqueue(Op::Put(to.to_string(), text));
            inner.enqueue(Op::Delete(from.to_string()));
        }
        self.pump();
        Ok(())
    }

    /// Everything is taken away, from memory at once and from the database when it is its turn: not only what this page put there,
    /// but what another page of the site did.
    pub fn clear(&self) {
        {
            let mut inner = self.inner.borrow_mut();
            inner.held.clear();
            inner.bytes = 0;
            inner.enqueue(Op::Clear);
        }
        self.pump();
    }

    /// Why the database refused a change, once: none when all went well, or when it was told already.
    pub fn take_failure(&self) -> Option<String> {
        self.inner.borrow_mut().failed.take()
    }

    /// Hand the database the change that is next, if it is not busy with one.
    fn pump(&self) {
        let op = {
            let mut inner = self.inner.borrow_mut();
            if inner.busy || inner.queue.is_empty() {
                return;
            }
            inner.busy = true;
            inner.queue.remove(0)
        };
        let table = self.clone();
        let asked = op.clone();
        self.persist
            .apply(op, Box::new(move |outcome| table.finished(&asked, outcome)));
    }

    /// A change is over. When it was refused, what it was to keep is not kept (unless a newer text has been put under its name
    /// since): it is taken out of the table, what the database may hold under its name is left to be taken away, and why is kept.
    fn finished(&self, op: &Op, outcome: Result<(), String>) {
        let wake = {
            let mut inner = self.inner.borrow_mut();
            inner.busy = false;
            match outcome {
                Ok(()) => None,
                Err(why) => {
                    match op {
                        Op::Put(name, text) => {
                            let now = inner.held.get(name).cloned();
                            if now.as_ref().is_some_and(|now| Rc::ptr_eq(now, text)) {
                                inner.held.remove(name);
                                inner.bytes -= size(name, text);
                            }
                            // a newer text that waits to be put will take its place in the database
                            if now.is_none_or(|now| Rc::ptr_eq(&now, text)) {
                                inner.stale.insert(name.clone());
                            }
                        }
                        Op::Delete(name) => {
                            if !inner.held.contains_key(name) {
                                inner.stale.insert(name.clone());
                            }
                        }
                        Op::Clear => inner.clear_again = true,
                    }
                    inner.failed.get_or_insert(why);
                    inner.wake.as_ref().map(|waker| waker.0.clone())
                }
            }
        };
        if let Some(wake) = wake {
            wake();
        }
        self.pump();
    }
}

/// Whether the file `name` holds a result, which is kept in the larger table (`session.rs` names them `generation-kind.result`).
pub fn is_large(name: &str) -> bool {
    name.ends_with(super::RESULT)
}

/// The two tables as one: the results in the larger one, everything else in the small one. Where there is no larger one (the
/// page has no database, or could not open it) the results are in the small one, as they were before it.
#[derive(Debug)]
pub struct Roomy {
    small: Rc<dyn Shelf>,
    large: Option<Large>,
}

impl Roomy {
    pub fn new(small: Rc<dyn Shelf>, large: Option<Large>) -> Roomy {
        Roomy { small, large }
    }

    /// The larger table, for a result's name, when there is one.
    fn larger(&self, name: &str) -> Option<&Large> {
        self.large.as_ref().filter(|_| is_large(name))
    }
}

impl Shelf for Roomy {
    /// A result that an earlier version kept in the small table is read from there when the larger one has none.
    fn get(&self, name: &str) -> Option<String> {
        match self.larger(name) {
            Some(large) => large.get(name).or_else(|| self.small.get(name)),
            None => self.small.get(name),
        }
    }

    fn holds(&self, name: &str) -> bool {
        self.larger(name).is_some_and(|large| large.holds(name)) || self.small.holds(name)
    }

    fn set(&self, name: &str, text: &str) -> Result<(), String> {
        match self.larger(name) {
            Some(large) => large.set(name, text),
            None => self.small.set(name, text),
        }
    }

    /// From both: a result that an earlier version kept in the small table goes too, when it is no longer wanted.
    fn remove(&self, name: &str) {
        if let Some(large) = self.larger(name) {
            large.remove(name);
        }
        self.small.remove(name);
    }

    /// A result in the larger table is moved inside it, with no second copy; any other (a result of an earlier version that is
    /// in the small table, a text of the session) is read, put under the new name, and taken away.
    fn rename(&self, from: &str, to: &str) -> Result<(), String> {
        if let Some(large) = self.larger(from).filter(|large| large.holds(from)) {
            if self.larger(to).is_some() {
                return large.rename(from, to);
            }
        }
        let Some(text) = self.get(from) else {
            return Ok(());
        };
        self.set(to, &text)?;
        self.remove(from);
        Ok(())
    }

    fn names(&self) -> Vec<String> {
        let mut names: BTreeSet<String> = self.small.names().into_iter().collect();
        if let Some(large) = &self.large {
            names.extend(large.names());
        }
        names.into_iter().collect()
    }

    /// A result takes the room of the sequences when it is in the small table (a page without a database, or a result of an
    /// earlier version that is still there); one in the larger table does not.
    fn shares_room(&self, name: &str) -> bool {
        self.larger(name).is_none() || self.small.holds(name)
    }

    fn clear_apart(&self) {
        if let Some(large) = &self.large {
            large.clear();
        }
    }

    fn take_problem(&self) -> Option<String> {
        self.large.as_ref()?.take_failure()
    }
}

/// A database that is in memory, for the tests: it applies a change at once, or holds it until the test says, or refuses it.
#[cfg(test)]
pub mod stand_in {
    use super::{Large, Op, Persist};
    use std::cell::{Cell, RefCell};
    use std::collections::BTreeMap;
    use std::rc::Rc;

    #[derive(Clone, Copy, Debug, PartialEq)]
    pub enum Mode {
        /// A change is applied and answered as it is handed over.
        Now,
        /// A change waits until the test lets it go (`finish`, `fail`).
        Later,
        /// A change is refused as it is handed over (the database is full, or is not allowed).
        Refuse,
    }

    type Done = Box<dyn FnOnce(Result<(), String>)>;

    pub struct Idb {
        pub stored: RefCell<BTreeMap<String, String>>,
        /// Every change that was handed over, in the order it was.
        pub ops: RefCell<Vec<Op>>,
        waiting: RefCell<Vec<(Op, Done)>>,
        pub mode: Cell<Mode>,
        /// A text longer than this is refused (the database is full for it), in the mode `Now`.
        pub quota: Cell<Option<usize>>,
    }

    impl std::fmt::Debug for Idb {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "a database of {} texts", self.stored.borrow().len())
        }
    }

    impl Idb {
        pub fn new(mode: Mode) -> Rc<Idb> {
            Rc::new(Idb {
                stored: RefCell::new(BTreeMap::new()),
                ops: RefCell::new(Vec::new()),
                waiting: RefCell::new(Vec::new()),
                mode: Cell::new(mode),
                quota: Cell::new(None),
            })
        }

        fn apply_to_stored(&self, op: &Op) {
            let mut stored = self.stored.borrow_mut();
            match op {
                Op::Put(name, text) => stored.insert(name.clone(), text.to_string()),
                Op::Delete(name) => stored.remove(name),
                Op::Clear => {
                    stored.clear();
                    None
                }
            };
        }

        /// How many changes wait for the test.
        pub fn waiting(&self) -> usize {
            self.waiting.borrow().len()
        }

        /// The oldest change that waits is applied, and its `done` is called with success.
        pub fn finish(&self) {
            let (op, done) = self.waiting.borrow_mut().remove(0);
            self.apply_to_stored(&op);
            done(Ok(()));
        }

        /// The oldest change that waits is not applied, and its `done` is called with `why`.
        pub fn fail(&self, why: &str) {
            let (_, done) = self.waiting.borrow_mut().remove(0);
            done(Err(why.to_string()));
        }

        /// What the database holds, as a page that is opened finds it.
        pub fn entries(&self) -> Vec<(String, String)> {
            self.stored
                .borrow()
                .iter()
                .map(|(name, text)| (name.clone(), text.clone()))
                .collect()
        }

        /// The names of what the database holds.
        pub fn names(&self) -> Vec<String> {
            self.stored.borrow().keys().cloned().collect()
        }

        /// The table of a page that is opened now, with `room` bytes.
        pub fn open(self: &Rc<Idb>, room: usize) -> Large {
            Large::new(self.clone(), self.entries(), room)
        }
    }

    impl Persist for Idb {
        fn apply(&self, op: Op, done: Done) {
            self.ops.borrow_mut().push(op.clone());
            let too_long =
                matches!(&op, Op::Put(_, text) if self.quota.get().is_some_and(|q| text.len() > q));
            match self.mode.get() {
                Mode::Now if too_long => done(Err("QuotaExceededError".to_string())),
                Mode::Now => {
                    self.apply_to_stored(&op);
                    done(Ok(()));
                }
                Mode::Later => self.waiting.borrow_mut().push((op, done)),
                Mode::Refuse => done(Err("QuotaExceededError".to_string())),
            }
        }
    }
}

#[cfg(test)]
mod tests;
