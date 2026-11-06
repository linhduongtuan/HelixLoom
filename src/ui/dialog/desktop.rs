//! The computer's file dialogs: the system's own, answering at once.

use super::{Arrived, Saved, Wanted};
use eframe::egui;
use std::path::PathBuf;

pub use rfd::FileDialog;

/// Ask for one file: its path, or none when the person cancelled. `wanted` matters only to a browser.
pub fn pick_one(_wanted: Wanted, dialog: FileDialog) -> Option<PathBuf> {
    dialog.pick_file()
}

/// Ask for any number of files: their paths.
pub fn pick_many(_wanted: Wanted, dialog: FileDialog) -> Vec<PathBuf> {
    dialog.pick_files().unwrap_or_default()
}

/// Ask where to save `bytes` as `name`, and write them there: none when the person cancelled.
pub fn save_bytes(dialog: FileDialog, name: &str, bytes: &[u8]) -> Option<Saved> {
    let path = dialog.set_file_name(name).save_file()?;
    Some(match std::fs::write(&path, bytes) {
        Ok(()) => Saved {
            message: format!("Saved {}.", path.display()),
            ok: true,
        },
        Err(e) => Saved {
            message: format!("Could not save {}: {e}", path.display()),
            ok: false,
        },
    })
}

/// The paths of the files dropped on the window.
pub fn dropped(files: &[egui::DroppedFileHandle], _wanted: Wanted) -> Vec<PathBuf> {
    files.iter().map(|f| f.path().to_path_buf()).collect()
}

/// A computer's dialogs answer at once, so nothing comes later.
pub fn take_arrived(_ctx: &egui::Context) -> Vec<Arrived> {
    Vec::new()
}
