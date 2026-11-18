//! The platform shell: a rail that directs to the analysis workspaces, a Home page, and the state each
//! workspace keeps. Every workspace is one page in its own file.

pub mod alignment;
mod alnview;
mod chromatogram;
#[cfg(test)]
mod empty_states;
mod home;
mod icons;
mod linked;
pub mod phylogeny;
pub mod pipelines;
mod rail;
mod report;
/// The texts of a report that other parts of the window check.
#[cfg(test)]
pub mod report_texts {
    pub use super::report::NOT_DRAWN;
}
pub mod route;
#[cfg(test)]
mod same_look;
pub mod signals;
pub mod structure;
pub mod trees;
mod treeview;
mod viewer3d;
pub mod workstation;

use crate::pipeline::tools::{self, Installed};
use crate::remote::weight::Weight;
use crate::ui::app::{Action, App};
use crate::ui::job::{Job, Poll};
use eframe::egui::{self, Ui};
use std::path::PathBuf;
use workstation::Place;

/// Width of the rail, in points.
pub const RAIL_WIDTH: f32 = 92.0;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Workspace {
    Home,
    Sequences,
    Signals,
    Alignment,
    Phylogeny,
    Trees,
    Structure,
    Pipelines,
    Workstation,
}

impl Workspace {
    pub const ALL: [Workspace; 9] = [
        Workspace::Home,
        Workspace::Sequences,
        Workspace::Signals,
        Workspace::Alignment,
        Workspace::Phylogeny,
        Workspace::Trees,
        Workspace::Structure,
        Workspace::Pipelines,
        Workspace::Workstation,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Workspace::Home => "Home",
            Workspace::Sequences => "Sequences",
            Workspace::Signals => "Signals",
            Workspace::Alignment => "Alignment",
            Workspace::Phylogeny => "Phylogeny",
            Workspace::Trees => "Trees",
            Workspace::Structure => "Structure",
            Workspace::Pipelines => "Pipelines",
            Workspace::Workstation => "Workstation",
        }
    }

    /// What the workspace is for, in one line: the rail's tooltip.
    pub fn blurb(self) -> &'static str {
        match self {
            Workspace::Home => "Where to start: every analysis in one place",
            Workspace::Sequences => "Read, edit and annotate DNA: maps, enzymes, ORFs, primers",
            Workspace::Signals => "Raw sequencing signals: Sanger chromatograms and reads",
            Workspace::Alignment => "Multiple sequence alignment",
            Workspace::Phylogeny => "Estimate phylogenetic trees",
            Workspace::Trees => "Draw and explore phylogenetic trees",
            Workspace::Structure => "Predict and view protein 3D structures",
            Workspace::Pipelines => "Run several steps in one go, or export a workflow",
            Workspace::Workstation => "Compute the heavy work on a workstation with GPUs",
        }
    }

    /// The name `--view` takes.
    pub fn view_name(self) -> &'static str {
        match self {
            Workspace::Home => "home",
            Workspace::Sequences => "sequences",
            Workspace::Signals => "signals",
            Workspace::Alignment => "alignment",
            Workspace::Phylogeny => "phylogeny",
            Workspace::Trees => "trees",
            Workspace::Structure => "structure",
            Workspace::Pipelines => "pipelines",
            Workspace::Workstation => "workstation",
        }
    }

    pub fn from_view(name: &str) -> Option<Workspace> {
        Workspace::ALL.into_iter().find(|w| w.view_name() == name)
    }
}

type Detected = Vec<(&'static str, Option<Installed>)>;

/// What the platform keeps between frames, beyond what each workspace keeps for itself.
#[derive(Default)]
pub struct Platform {
    programs: Option<Detected>,
    detecting: Option<Job<Detected>>,
    /// The programs were asked to be looked for again while a search was going: what that search finds is from before.
    stale: bool,
    pub alignment: alignment::AlignmentPage,
    pub phylogeny: phylogeny::PhylogenyPage,
    pub trees: trees::TreesPage,
    pub signals: signals::SignalsPage,
    pub structure: structure::StructurePage,
    pub pipelines: pipelines::PipelinesPage,
    pub workstation: workstation::WorkstationPage,
}

impl Platform {
    /// A platform whose workstation connection is kept in `settings_dir` (none: nothing is read or kept).
    pub fn new(settings_dir: Option<PathBuf>) -> Platform {
        Platform {
            workstation: workstation::WorkstationPage::open(settings_dir),
            ..Platform::default()
        }
    }

    /// How many analyses the pages are running now, here or on the workstation (a page that waits for a workstation's answer
    /// is running too). The windows keep their own.
    pub fn busy_count(&self) -> usize {
        [
            self.alignment.is_running(),
            self.phylogeny.is_running(),
            self.structure.is_running(),
            self.pipelines.is_running(),
        ]
        .into_iter()
        .filter(|running| *running)
        .count()
    }

    /// Where a job of this weight runs, by the choice made on the Workstation page. Every analysis asks here before it
    /// starts, so a new one gets the choice (and the workstation) by describing how big its jobs are.
    pub fn place(&self, weight: &Weight) -> Place {
        self.workstation.place(weight)
    }

    /// The programs a job at `place` would use, as (name, found): the workstation's for a job that goes there (known once
    /// it has said what it has), else this computer's (known once the search ends).
    pub fn programs_at(&self, place: &Place) -> Option<Detected> {
        if !place.is_there() {
            return self.programs.clone();
        }
        let theirs = self.workstation.tools()?;
        Some(
            tools::ALL
                .iter()
                .map(|p| {
                    let found = theirs.iter().find(|(name, _)| name == p.label);
                    // only the version is known of a program on another computer
                    let found = found.map(|(_, version)| Installed {
                        path: PathBuf::new(),
                        version: version.clone(),
                    });
                    (p.label, found)
                })
                .collect(),
        )
    }

    /// The program called `label` at `place`, if it was found.
    /// For a test: the search for programs on this computer is over, and it found none.
    #[cfg(test)]
    pub(in crate::ui) fn no_programs_here(&mut self) {
        (self.programs, self.detecting) = (Some(Vec::new()), None);
    }

    pub fn installed_at(&self, label: &str, place: &Place) -> Option<Installed> {
        self.programs_at(place)?
            .into_iter()
            .find(|(l, _)| *l == label)
            .and_then(|(_, found)| found)
    }

    /// Forget what was found on this computer, so that the programs are looked for again (one was installed).
    pub fn look_again(&mut self) {
        // a search that is going began before the install: its answer is thrown away, and another search made
        self.stale = self.detecting.is_some();
        self.programs = None;
    }

    /// Take the answer of a search for the programs: kept, unless the programs were asked for again while it went.
    fn found_programs(&mut self, found: Detected) {
        self.detecting = None;
        self.programs = if std::mem::take(&mut self.stale) {
            None
        } else {
            Some(found)
        };
    }

    /// The program called `label` on this computer, if it was found.
    pub fn installed(&self, label: &str) -> Option<Installed> {
        self.installed_at(label, &Place::Here)
    }

    /// Whether the search for the programs at `place` has finished, so a program that `installed_at` does not give is
    /// missing.
    pub fn programs_known_at(&self, place: &Place) -> bool {
        self.programs_at(place).is_some()
    }
}

impl App {
    /// How many analyses are running now: those of the pages and the PhyloCNN and species-tree windows.
    pub fn running_analyses(&self) -> usize {
        self.platform.busy_count()
            + usize::from(self.phylocnn_win.is_running())
            + usize::from(self.species_win.is_running())
            + usize::from(self.motif_win.running())
            + usize::from(self.domains_win.running())
    }

    /// Look for the programs the platform can call, once, without making the window wait.
    fn detect_programs(&mut self, ctx: &egui::Context) {
        let p = &mut self.platform;
        if p.programs.is_none() && p.detecting.is_none() {
            p.detecting = Some(Job::spawn(ctx, || {
                tools::ALL
                    .iter()
                    .map(|p| (p.label, tools::detect(p)))
                    .collect()
            }));
        }
        if let Some(job) = &p.detecting {
            match job.poll() {
                Poll::Done(found) => p.found_programs(found),
                Poll::Died => (p.programs, p.detecting) = (Some(Vec::new()), None),
                Poll::Running => ctx.request_repaint_after(std::time::Duration::from_millis(200)),
            }
        }
    }

    /// What fills the window beside the rail.
    pub(in crate::ui) fn workspace_content(&mut self, ui: &mut Ui, out: &mut Vec<Action>) {
        // work that was started goes on, and finishes, whichever page is showing
        self.detect_programs(ui.ctx());
        self.alignment_tick(ui.ctx());
        self.phylogeny_tick(ui.ctx());
        self.structure_tick(ui.ctx());
        self.pipelines_tick(ui.ctx());
        self.workstation_tick(ui.ctx());
        match self.workspace {
            Workspace::Home => self.home_page(ui, out),
            Workspace::Sequences => self.content(ui),
            Workspace::Alignment => self.alignment_page(ui, out),
            Workspace::Phylogeny => self.phylogeny_page(ui, out),
            Workspace::Trees => self.trees_page(ui, out),
            Workspace::Signals => self.signals_page(ui, out),
            Workspace::Structure => self.structure_page(ui),
            Workspace::Pipelines => self.pipelines_page(ui, out),
            Workspace::Workstation => self.workstation_page(ui, out),
        }
    }

    /// Go to a workspace with its example loaded, so the page shows what it does.
    pub(in crate::ui) fn load_example(&mut self, ws: Workspace) {
        self.workspace = ws;
        match ws {
            Workspace::Alignment => self.load_alignment_example(),
            Workspace::Phylogeny => {
                self.load_alignment_example();
                self.platform.phylogeny.run_when_ready = true;
            }
            Workspace::Trees => self.show_tree(
                crate::sample::demo_tree(),
                "The example tree: 24 made-up strains in three groups.",
            ),
            Workspace::Signals => self
                .platform
                .signals
                .load("example trace", crate::sample::demo_trace()),
            Workspace::Structure => self.platform.structure.load_example(),
            Workspace::Pipelines => {
                let page = &mut self.platform.pipelines;
                page.source = pipelines::Source::Example;
                page.run_requested = true;
            }
            Workspace::Home | Workspace::Sequences | Workspace::Workstation => {}
        }
    }

    /// The six demo isolates, ticked and set to be aligned as soon as the page shows.
    fn load_alignment_example(&mut self) {
        let page = &mut self.platform.alignment;
        page.clear_inputs();
        for doc in crate::sample::tree_demo() {
            page.add(&doc.name, &doc.seq);
        }
        page.run_when_shown = true;
    }
}

#[cfg(test)]
mod look_again_tests {
    use super::*;

    #[test]
    fn a_search_that_was_going_when_a_program_was_installed_is_not_believed() {
        let mut p = Platform::default();
        let before: Detected = vec![("MAFFT", None)];
        // no search is going: what was found is forgotten, and the next frame searches
        p.programs = Some(before.clone());
        p.look_again();
        assert!(p.programs.is_none() && !p.stale);
        // a search is going (it began before the install): its answer is from before, so it is thrown away
        let ctx = egui::Context::default();
        p.detecting = Some(Job::spawn(&ctx, Vec::new));
        p.look_again();
        assert!(p.stale);
        p.found_programs(before.clone());
        assert!(p.programs.is_none() && p.detecting.is_none() && !p.stale);
        // and the search after it is believed
        p.found_programs(before.clone());
        assert_eq!(p.programs, Some(before));
    }
}
