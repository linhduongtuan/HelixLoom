//! Which workspace a file belongs in, told by the end of its name, so a file opened from the menu, dropped on the
//! window or named on the command line lands where it can be used.

use super::Workspace;
use crate::ui::app::App;
use std::borrow::Cow;
use std::io;
use std::path::Path;

/// A file to open. On a computer it is a path, read when it is needed; in a browser a page cannot see where a file is, so
/// the file comes as its name and its bytes. Every loader takes a `Source`, so a file opens the same way whichever it is.
#[derive(Clone, Copy, Debug)]
pub enum Source<'a> {
    Path(&'a Path),
    Bytes { name: &'a str, bytes: &'a [u8] },
}

impl<'a> Source<'a> {
    /// The path, or the bare name of a file from a browser: enough to tell its kind and to take its stem.
    pub fn path(&self) -> &'a Path {
        match *self {
            Source::Path(path) => path,
            Source::Bytes { name, .. } => Path::new(name),
        }
    }

    /// Where the file is on this computer, if it is on one.
    pub fn on_disk(&self) -> Option<&'a Path> {
        match *self {
            Source::Path(path) => Some(path),
            Source::Bytes { .. } => None,
        }
    }

    /// How messages name the file.
    pub fn shown(&self) -> String {
        self.path().display().to_string()
    }

    pub fn bytes(&self) -> io::Result<Cow<'a, [u8]>> {
        match *self {
            Source::Path(path) => std::fs::read(path).map(Cow::Owned),
            Source::Bytes { bytes, .. } => Ok(Cow::Borrowed(bytes)),
        }
    }

    /// The file as text; bytes that are not UTF-8 are an error, as they are for `std::fs::read_to_string`.
    pub fn text(&self) -> io::Result<String> {
        match *self {
            Source::Path(path) => std::fs::read_to_string(path),
            Source::Bytes { bytes, .. } => String::from_utf8(bytes.to_vec())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    Trace,
    Structure,
    Tree,
    Sequence,
}

/// File name endings of each kind, as the open dialog offers them.
/// The original's list of file endings it opens as traces (`possible_extensions`, AppMain.tcl:67221): `.abi`, `.ab1`
/// and `.scf`.
pub const TRACE_EXTENSIONS: [&str; 3] = ["ab1", "abi", "scf"];
pub const STRUCTURE_EXTENSIONS: [&str; 2] = ["pdb", "ent"];
pub const TREE_EXTENSIONS: [&str; 5] = ["nwk", "newick", "tre", "tree", "treefile"];

pub fn kind_of(path: &Path) -> Kind {
    let ending = path
        .extension()
        .and_then(|e| e.to_str())
        .map(str::to_ascii_lowercase);
    match ending.as_deref() {
        Some(e) if TRACE_EXTENSIONS.contains(&e) => Kind::Trace,
        Some(e) if STRUCTURE_EXTENSIONS.contains(&e) => Kind::Structure,
        Some(e) if TREE_EXTENSIONS.contains(&e) => Kind::Tree,
        _ => Kind::Sequence,
    }
}

impl App {
    /// Open the file at `path` where it belongs (see `open_source`).
    pub fn open_file(&mut self, path: &Path) {
        self.open_source(Source::Path(path));
    }

    /// Open a file where it belongs and go there: a trace in Signals, a structure in Structure, a tree in Trees,
    /// sequences in Sequences, or added to the list when Alignment or Phylogeny is showing.
    pub fn open_source(&mut self, file: Source) {
        match kind_of(file.path()) {
            Kind::Trace => {
                self.load_trace(file);
                self.workspace = Workspace::Signals;
            }
            Kind::Structure => {
                self.platform.structure.open_source(file);
                self.workspace = Workspace::Structure;
            }
            Kind::Tree => {
                self.load_tree(file);
                self.workspace = Workspace::Trees;
            }
            Kind::Sequence
                if matches!(self.workspace, Workspace::Alignment | Workspace::Phylogeny) =>
            {
                self.platform.alignment.add_source(file);
            }
            Kind::Sequence => {
                self.load_sequences(file);
                self.workspace = Workspace::Sequences;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_end_of_the_name_tells_the_kind_whatever_its_case() {
        let kind = |name: &str| kind_of(Path::new(name));
        assert_eq!(kind("run/sample_A01.ab1"), Kind::Trace);
        assert_eq!(kind("SAMPLE.ABI"), Kind::Trace);
        assert_eq!(kind("run/read.scf"), Kind::Trace);
        assert_eq!(kind("1ubq.pdb"), Kind::Structure);
        assert_eq!(kind("pdb1ubq.ent"), Kind::Structure);
        assert_eq!(kind("out.treefile"), Kind::Tree);
        assert_eq!(kind("strains.NWK"), Kind::Tree);
        assert_eq!(kind("pUC19.gb"), Kind::Sequence);
        assert_eq!(kind("reads.fasta"), Kind::Sequence);
        assert_eq!(kind("notes"), Kind::Sequence);
        assert_eq!(kind("archive.ab1.bak"), Kind::Sequence);
    }

    #[test]
    fn every_listed_ending_is_in_lower_case() {
        for e in TRACE_EXTENSIONS
            .iter()
            .chain(&STRUCTURE_EXTENSIONS)
            .chain(&TREE_EXTENSIONS)
        {
            assert_eq!(*e, e.to_ascii_lowercase());
        }
    }
}
