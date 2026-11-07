//! The file dialogs of a computer and of a browser, behind the same few calls.
//!
//! A computer has the system's dialogs, which answer at once with paths. A browser has no such thing: the page's file
//! chooser answers later, with each file's name and bytes (a page cannot see where a file is), and a file is saved by
//! handing it to the browser as a download. So an open dialog is asked for together with what its files are for
//! (`Wanted`): on a computer the paths come back from the call, in a browser the files come back through `take_arrived`,
//! which the window asks once a frame.

#[cfg(target_arch = "wasm32")]
mod browser;
#[cfg(not(target_arch = "wasm32"))]
mod desktop;

#[cfg(target_arch = "wasm32")]
pub use browser::{dropped, pick_many, pick_one, save_bytes, take_arrived, FileDialog};
#[cfg(not(target_arch = "wasm32"))]
pub use desktop::{dropped, pick_many, pick_one, save_bytes, take_arrived, FileDialog};

/// Whether this is the browser, where no folder can be chosen and nothing is written to a disk.
pub const IN_BROWSER: bool = cfg!(target_arch = "wasm32");

/// What an open dialog is for, which says where its files go when a browser brings them.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Wanted {
    /// Whatever they are: each goes where its name says (File → Open, files dropped on the window).
    Anything,
    /// Sequences for the Contig Assembler's list.
    AssemblyReads,
    /// Sequences for the Alignment page.
    Alignment,
    /// A sequencing trace for the Signals page.
    Trace,
    /// A Newick tree for the Trees page.
    Tree,
    /// A PDB file for the Structure page.
    Structure,
    /// Reads for the read-quality window.
    Fastq,
    /// A feature library for annotating.
    FeatureLibrary,
    /// A Newick tree for the Tree Tools window.
    TreeToolsTree,
    /// A Newick tree for the PhyloCNN window.
    PhyloCnnTree,
    /// A file of gene trees for the species-tree window.
    SpeciesTrees,
}

/// Files a browser brought: each one's name and bytes, or why it could not be read.
// only a browser makes one of these
#[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
pub struct Arrived {
    pub wanted: Wanted,
    pub files: Vec<(String, Result<Vec<u8>, String>)>,
}

/// What came of saving a file.
pub struct Saved {
    /// What to tell the person.
    pub message: String,
    /// Whether the file was written (in a browser: handed to the browser).
    // only the browser's way of saving a sequence looks at it
    #[cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]
    pub ok: bool,
}
