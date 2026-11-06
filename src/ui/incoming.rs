//! Files that reach the window without being asked for by name: dropped on it, or, in a browser, chosen in the page's file
//! chooser, which answers later with each file's name and bytes (see `dialog`).

use super::app::{App, Tab};
use super::dialog::{self, Arrived, Wanted};
use super::platform::route::Source;
use super::platform::Workspace;
use eframe::egui;
use std::path::PathBuf;

impl App {
    /// Files dropped on the window open in the workspace they belong to, or join the assembly when that tab is showing.
    /// What a browser brought since the last frame is taken in too.
    pub(super) fn take_dropped_files(&mut self, ctx: &egui::Context) {
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        if !dropped.is_empty() {
            let wanted = if self.workspace == Workspace::Sequences && self.tab == Tab::Assembly {
                Wanted::AssemblyReads
            } else {
                Wanted::Anything
            };
            let paths = dialog::dropped(&dropped, wanted);
            self.open_paths(wanted, &paths);
        }
        for arrived in dialog::take_arrived(ctx) {
            self.take_arrived(arrived);
        }
    }

    /// Files on this computer, for what they are wanted for.
    fn open_paths(&mut self, wanted: Wanted, paths: &[PathBuf]) {
        if wanted == Wanted::AssemblyReads {
            self.add_read_files(paths);
        } else {
            paths.iter().for_each(|p| self.open_file(p));
        }
    }

    /// Files a browser brought: each is read from its bytes by what it is wanted for, like a file from a path.
    fn take_arrived(&mut self, arrived: Arrived) {
        let mut files = Vec::new();
        for (name, bytes) in &arrived.files {
            match bytes {
                Ok(bytes) => files.push(Source::Bytes { name, bytes }),
                Err(why) => self.say(format!("Couldn't read {name}: {why}")),
            }
        }
        if arrived.wanted == Wanted::AssemblyReads {
            self.add_read_sources(&files);
        } else {
            for file in files {
                self.take_file(arrived.wanted, file);
            }
        }
    }

    fn take_file(&mut self, wanted: Wanted, file: Source) {
        match wanted {
            Wanted::Anything => self.open_source(file),
            Wanted::Alignment => self.platform.alignment.add_source(file),
            Wanted::Trace => self.load_trace(file),
            Wanted::Tree => self.load_tree(file),
            Wanted::Structure => self.platform.structure.open_source(file),
            Wanted::Fastq => self.load_fastq(file),
            Wanted::FeatureLibrary => self.load_feature_library(file),
            Wanted::TreeToolsTree => self.load_treetools_tree(file),
            Wanted::PhyloCnnTree => self.load_phylocnn_tree(file),
            // a batch, taken by `take_arrived`
            Wanted::AssemblyReads => self.add_read_sources(&[file]),
        }
    }
}

#[cfg(test)]
mod tests;
