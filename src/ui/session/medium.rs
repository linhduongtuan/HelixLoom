//! Where a store keeps what it keeps: the settings folder of a computer, or, in a browser, the storage the browser gives the
//! page (`docs/WEB.md`). The store asks for the same few things of both, by a name that is a file of the folder (`places.txt`)
//! or a file of a folder in it (`session/3-0.gb`), so that what is kept, in which order and what is cleaned away is written
//! once (`session.rs`) and is the same in both.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::rc::Rc;

/// A table of names and texts that outlives the run: what a browser gives a page to keep things in (`browser.rs`), and what a
/// test puts in its place.
pub trait Shelf: std::fmt::Debug {
    fn get(&self, name: &str) -> Option<String>;
    /// Keep `text` under `name`; why not, when the table is full or may not be written.
    fn set(&self, name: &str, text: &str) -> Result<(), String>;
    fn remove(&self, name: &str);
    /// Every name the table holds.
    fn names(&self) -> Vec<String>;
    /// Whether a text is kept under `name`, without reading it (a result may be megabytes).
    fn holds(&self, name: &str) -> bool {
        self.get(name).is_some()
    }
    /// Put the text under `from` under `to`, and take `from` away: nothing is done when there is none; refused (and `from` as it
    /// was) when `to` cannot be written. A table that can move a text without copying it does.
    fn rename(&self, from: &str, to: &str) -> Result<(), String> {
        let Some(text) = self.get(from) else {
            return Ok(());
        };
        self.set(to, &text)?;
        self.remove(from);
        Ok(())
    }
    /// Whether writing `name` takes the room that the open sequences are written in. A result that is kept apart (a page's larger
    /// table) does not: when a sequence does not fit, taking it away would not make room.
    fn shares_room(&self, _name: &str) -> bool {
        true
    }
    /// Take away everything that the table keeps apart (a page's larger table), whoever put it there. The names of the table are
    /// removed by name; this is for what the page does not know the names of (another page of the site wrote them).
    fn clear_apart(&self) {}
    /// Why a write that was not waited for did not happen, once; none for a table that waits for every write.
    fn take_problem(&self) -> Option<String> {
        None
    }
}

#[derive(Clone, Debug)]
pub enum Medium {
    Folder(PathBuf),
    // only a page, or a test, has a table to keep things in
    #[cfg_attr(not(any(target_arch = "wasm32", test)), allow(dead_code))]
    Shelf(Rc<dyn Shelf>),
}

/// Make a folder that only the person can read (on Unix; elsewhere the system's default).
fn make_private_dir(dir: &Path) -> Result<(), String> {
    let mut builder = fs::DirBuilder::new();
    builder.recursive(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::DirBuilderExt;
        builder.mode(0o700);
    }
    builder
        .create(dir)
        .map_err(|e| format!("Could not make {}: {e}", dir.display()))
}

/// Ask the system to put the folder's entries on disk, so that a file renamed into it is still there after a power cut. Best effort:
/// some systems cannot open a folder to do it.
fn sync_dir(dir: &Path) {
    if let Ok(handle) = fs::File::open(dir) {
        let _ = handle.sync_all();
    }
}

/// Write `bytes` to `path` so that it is whole or not there: to a file beside it, put on disk, then renamed over it. `private` makes the
/// file readable by the person only (on Unix), and puts the folder's entry on disk too.
pub fn write_atomic(path: &Path, bytes: &[u8], private: bool) -> Result<(), String> {
    use std::io::Write;
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("file");
    let temporary = path.with_file_name(format!(".{name}.tmp"));
    let written = (|| -> std::io::Result<()> {
        let mut options = fs::OpenOptions::new();
        options.write(true).create(true).truncate(true);
        #[cfg(unix)]
        if private {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        let mut file = options.open(&temporary)?;
        file.write_all(bytes)?;
        file.sync_all()
    })();
    if let Err(e) = written {
        let _ = fs::remove_file(&temporary);
        return Err(format!("Could not write {}: {e}", temporary.display()));
    }
    fs::rename(&temporary, path).map_err(|e| {
        let _ = fs::remove_file(&temporary);
        format!("Could not replace {}: {e}", path.display())
    })?;
    if private {
        if let Some(parent) = path.parent() {
            sync_dir(parent);
        }
    }
    Ok(())
}

/// The names of a shelf that are in `folder`, each with the folder taken off.
fn in_folder(shelf: &dyn Shelf, folder: &str) -> Vec<String> {
    let prefix = format!("{folder}/");
    shelf
        .names()
        .into_iter()
        .filter_map(|name| name.strip_prefix(&prefix).map(str::to_string))
        .collect()
}

impl Medium {
    /// The text kept under `name`. With `most`, what is read is a short plain file and at most so many bytes of it (not a pipe
    /// or a device that something else put there, which reading could wait on for ever; of a table's text, so many bytes cut
    /// back to where a letter ends); none when there is none.
    pub fn read(&self, name: &str, most: Option<u64>) -> Option<String> {
        match (self, most) {
            (Medium::Folder(dir), None) => fs::read_to_string(dir.join(name)).ok(),
            (Medium::Folder(dir), Some(most)) => {
                let path = dir.join(name);
                if !fs::metadata(&path).ok()?.is_file() {
                    return None;
                }
                let mut bytes = Vec::new();
                fs::File::open(&path)
                    .ok()?
                    .take(most)
                    .read_to_end(&mut bytes)
                    .ok()?;
                Some(String::from_utf8_lossy(&bytes).into_owned())
            }
            (Medium::Shelf(shelf), most) => {
                let text = shelf.get(name)?;
                // a text that is read whole is not copied again (a result may be megabytes)
                let Some(most) = most else {
                    return Some(text);
                };
                let most = usize::try_from(most).unwrap_or(usize::MAX);
                let end = (0..=most.min(text.len()))
                    .rev()
                    .find(|&at| text.is_char_boundary(at))
                    .unwrap_or(0);
                Some(text[..end].to_string())
            }
        }
    }

    /// Keep `text` under `name`, whole or not at all. `private` is for the files of the session folder, which only the person
    /// may read and whose folder is made beforehand (`make_folder`); any other file is in the settings folder itself, which is
    /// made if it is not there.
    pub fn write(&self, name: &str, text: &str, private: bool) -> Result<(), String> {
        match self {
            Medium::Folder(dir) => {
                if !private {
                    fs::create_dir_all(dir)
                        .map_err(|e| format!("Could not make {}: {e}", dir.display()))?;
                }
                write_atomic(&dir.join(name), text.as_bytes(), private)
            }
            Medium::Shelf(shelf) => shelf.set(name, text),
        }
    }

    /// Whether a text is kept under `name` (in a folder: a plain file).
    pub fn holds(&self, name: &str) -> bool {
        match self {
            Medium::Folder(dir) => dir.join(name).is_file(),
            Medium::Shelf(shelf) => shelf.holds(name),
        }
    }

    /// Whether writing `name` takes the room that the other files of the session are written in (a folder: it does; a table may
    /// keep some of its names apart, `Shelf::shares_room`).
    pub fn shares_room(&self, name: &str) -> bool {
        match self {
            Medium::Folder(_) => true,
            Medium::Shelf(shelf) => shelf.shares_room(name),
        }
    }

    /// Take away what a table keeps apart, whoever put it there (a folder has nothing apart).
    pub fn clear_apart(&self) {
        if let Medium::Shelf(shelf) = self {
            shelf.clear_apart();
        }
    }

    pub fn remove(&self, name: &str) {
        match self {
            Medium::Folder(dir) => {
                let _ = fs::remove_file(dir.join(name));
            }
            Medium::Shelf(shelf) => shelf.remove(name),
        }
    }

    /// Make `folder`, private to the person, if it is not there. A table has no folders to make: a name says which it is in.
    pub fn make_folder(&self, folder: &str) -> Result<(), String> {
        match self {
            Medium::Folder(dir) => make_private_dir(&dir.join(folder)),
            Medium::Shelf(_) => Ok(()),
        }
    }

    pub fn has_folder(&self, folder: &str) -> bool {
        match self {
            Medium::Folder(dir) => dir.join(folder).exists(),
            Medium::Shelf(shelf) => !in_folder(shelf.as_ref(), folder).is_empty(),
        }
    }

    /// The names of what `folder` holds, without the folder.
    pub fn names_in(&self, folder: &str) -> Vec<String> {
        match self {
            Medium::Folder(dir) => fs::read_dir(dir.join(folder))
                .map(|entries| {
                    entries
                        .flatten()
                        .map(|e| e.file_name().to_string_lossy().into_owned())
                        .collect()
                })
                .unwrap_or_default(),
            Medium::Shelf(shelf) => in_folder(shelf.as_ref(), folder),
        }
    }

    /// Remove `folder` and all it holds. One that is not there is an error of the kind `NotFound` in a folder of the computer,
    /// as the system says it, and nothing at all in a table.
    pub fn remove_folder(&self, folder: &str) -> std::io::Result<()> {
        match self {
            Medium::Folder(dir) => fs::remove_dir_all(dir.join(folder)),
            Medium::Shelf(shelf) => {
                for name in in_folder(shelf.as_ref(), folder) {
                    shelf.remove(&format!("{folder}/{name}"));
                }
                Ok(())
            }
        }
    }

    /// Put what `from` holds in `to` instead, in the place of whatever `to` held.
    pub fn move_folder(&self, from: &str, to: &str) -> Result<(), String> {
        let _ = self.remove_folder(to);
        match self {
            Medium::Folder(dir) => {
                fs::rename(dir.join(from), dir.join(to)).map_err(|e| e.to_string())
            }
            // a table moves one text at a time (`Shelf::rename`: it needs room for one text twice, not for all; a table that keeps
            // some apart moves them without a copy). When one cannot be moved (the table is full) what was moved is moved back,
            // and `from` is as it was: a move is whole or it has not happened, as a rename is. A page runs this within a frame,
            // so there is no moment in which it is half done for anyone else.
            Medium::Shelf(shelf) => {
                let names = in_folder(shelf.as_ref(), from);
                let mut moved: Vec<&String> = Vec::new();
                for name in &names {
                    if let Err(why) =
                        shelf.rename(&format!("{from}/{name}"), &format!("{to}/{name}"))
                    {
                        for done in moved.iter().rev() {
                            let _ =
                                shelf.rename(&format!("{to}/{done}"), &format!("{from}/{done}"));
                        }
                        return Err(why);
                    }
                    moved.push(name);
                }
                Ok(())
            }
        }
    }
}

/// A table in memory, for the tests: it stands where a browser's storage does. It can be full at some size, as a browser's is,
/// and it can refuse one name.
#[cfg(test)]
#[derive(Debug, Default)]
pub struct Memory {
    pub held: std::cell::RefCell<std::collections::BTreeMap<String, String>>,
    /// The most letters it holds in all its texts and names; none: any number.
    pub room: std::cell::Cell<Option<usize>>,
    /// A name that cannot be written.
    pub refused: std::cell::RefCell<Option<String>>,
    /// The names are given last first: a browser gives the names of its table in no order that can be counted on.
    pub backwards: std::cell::Cell<bool>,
}

#[cfg(test)]
impl Shelf for Memory {
    fn get(&self, name: &str) -> Option<String> {
        self.held.borrow().get(name).cloned()
    }

    fn set(&self, name: &str, text: &str) -> Result<(), String> {
        if self.refused.borrow().as_deref() == Some(name) {
            return Err(format!("{name} may not be written"));
        }
        let mut held = self.held.borrow_mut();
        let others: usize = held
            .iter()
            .filter(|(key, _)| key.as_str() != name)
            .map(|(key, value)| key.len() + value.len())
            .sum();
        if self
            .room
            .get()
            .is_some_and(|room| others + name.len() + text.len() > room)
        {
            return Err(super::browser::FULL.to_string());
        }
        held.insert(name.to_string(), text.to_string());
        Ok(())
    }

    fn remove(&self, name: &str) {
        self.held.borrow_mut().remove(name);
    }

    fn names(&self) -> Vec<String> {
        let mut names: Vec<String> = self.held.borrow().keys().cloned().collect();
        if self.backwards.get() {
            names.reverse();
        }
        names
    }
}
