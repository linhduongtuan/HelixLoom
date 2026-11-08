//! Where a job runs: the choice of this computer, the workstation, or each job by its size, and the small control that every
//! page that starts a job shows beside its button.

use super::{dim, Clicks, WorkstationPage};
use crate::remote::client::Remote;
use crate::remote::settings::{Compute, GpusPerJob};
use crate::remote::weight::{
    HEAVY_GENERATIONS, HEAVY_LETTERS, HEAVY_REPLICATES, HEAVY_SEQUENCES, HEAVY_TREES,
};
use crate::ui::design::{self, Look};
use crate::ui::map::commify;
use eframe::egui::{self, RichText, Ui};

/// Where one job runs: on this computer, or on the workstation (with the connection to it).
pub enum Place {
    Here,
    There(Remote),
}

impl Place {
    pub fn is_there(&self) -> bool {
        matches!(self, Place::There(_))
    }
}

/// What a page shows about where its job would run, and the choice that decides it.
#[derive(Clone, Copy)]
pub(in crate::ui) struct ComputeView {
    pub mode: Compute,
    /// A workstation is connected: the choices of Workstation and Automatic need one.
    pub connected: bool,
    /// The job as the page has it set up would run on the workstation.
    pub there: bool,
}

/// What each choice of where jobs run means.
pub(super) fn explain(choice: Compute) -> String {
    match choice {
        Compute::Local => "Every job runs on this computer. Nothing is sent anywhere.".into(),
        Compute::Workstation => "Alignments, tree estimates, pipelines, PhyloCNN estimates and structure predictions are \
                                 computed there. The results come back here and open as usual, and they wait on the \
                                 workstation until you fetch them."
            .into(),
        Compute::Automatic => format!(
            "Big jobs go to the workstation and small ones stay here, with nothing to upload. A job is big when it has \
             more than {} sequences or {} letters, asks for more than {} bootstrap resamples or {} Bayesian generations, \
             has more than {} trees, or can use a GPU there.",
            commify(HEAVY_SEQUENCES),
            commify(HEAVY_LETTERS),
            commify(HEAVY_REPLICATES),
            commify(HEAVY_GENERATIONS),
            commify(HEAVY_TREES),
        ),
    }
}

/// What each choice of how many GPUs one job is given means.
pub(super) fn explain_gpus(choice: GpusPerJob) -> String {
    match choice {
        GpusPerJob::Count(1) => "Each job that needs GPUs is given one, so several jobs run side by side, a GPU each.".into(),
        GpusPerJob::Count(n) => format!(
            "A job that can use several GPUs (PhyloCNN with many trees) is given {n} at once and shares its trees out over \
             them; a job that can use fewer is given what it can use. The jobs behind it wait for GPUs in the order they \
             were sent."
        ),
        GpusPerJob::All => "A job that can use several GPUs is given every GPU the workstation has: one job at a time, as \
                            fast as the GPUs together allow."
            .into(),
    }
}

/// The choices of how many GPUs one job is given on a workstation with `have` GPUs: one up to all of them, and every GPU.
pub(super) fn gpu_choices(have: usize) -> Vec<GpusPerJob> {
    let counts = (1..=have.min(crate::remote::settings::MOST_GPUS_PER_JOB)).map(GpusPerJob::Count);
    counts.chain([GpusPerJob::All]).collect()
}

/// How many GPUs one job is given, when the workstation has several. It is the choice that a job that can use several GPUs is
/// sent with; with one GPU there is nothing to choose.
fn gpus_section(ui: &mut Ui, look: &Look, page: &WorkstationPage, c: &mut Clicks) {
    if page.active.is_none() || page.gpu_count() < 2 {
        return;
    }
    design::section(ui, look, "GPUs for one job");
    let mut chosen = page.gpus;
    egui::ComboBox::from_id_salt("workstation_gpus_per_job")
        .selected_text(chosen.label())
        .show_ui(ui, |ui| {
            for choice in gpu_choices(page.gpu_count()) {
                ui.selectable_value(&mut chosen, choice, choice.label());
            }
        });
    if chosen != page.gpus {
        c.gpus = Some(chosen);
    }
    ui.label(dim(look, explain_gpus(page.gpus)));
}

/// Where the heavy jobs are computed: this computer, the workstation, or each job by its size, and how many GPUs of the
/// workstation one job is given. The same choice of where is on every page that starts a job.
pub(super) fn heavy_jobs_section(ui: &mut Ui, look: &Look, page: &WorkstationPage, c: &mut Clicks) {
    design::section(ui, look, "Where jobs run");
    let connected = page.active.is_some();
    let mut mode = page.compute;
    for choice in Compute::ALL {
        // this computer can always be chosen; the others need a workstation to send to
        let usable = connected || choice == Compute::Local;
        ui.add_enabled_ui(usable, |ui| {
            ui.radio_value(&mut mode, choice, choice.label());
        });
        ui.label(dim(look, explain(choice)));
    }
    if mode != page.compute {
        c.choose = Some(mode);
    }
    if !connected {
        ui.label(dim(look, "Connect a workstation to send jobs there."));
    }
    gpus_section(ui, look, page, c);
}

/// Where a click on the choice waits for the next frame: the pages that show the choice are drawn while the page that
/// keeps it is borrowed.
pub(super) fn choice_id() -> egui::Id {
    egui::Id::new("workstation_compute_choice")
}

/// The choice of where jobs run, drawn where a page starts its work, with where this job would run. A click is applied at
/// the start of the next frame. With no workstation connected there is nothing to choose, and the line says how to get one.
pub(in crate::ui) fn compute_picker(ui: &mut Ui, look: &Look, view: ComputeView) {
    if !view.connected {
        ui.label(dim(
            look,
            "Runs on this computer. To compute on a workstation, connect one on the Workstation page.",
        ));
        return;
    }
    // small and close together, so that the three fit in the narrowest options panel (and wrap, not spill, if one does not)
    ui.horizontal_wrapped(|ui| {
        ui.spacing_mut().item_spacing.x = 4.0;
        let mut mode = view.mode;
        for choice in Compute::ALL {
            ui.selectable_value(&mut mode, choice, RichText::new(choice.label()).size(12.0))
                .on_hover_text(explain(choice));
        }
        if mode != view.mode {
            ui.ctx().data_mut(|d| d.insert_temp(choice_id(), mode));
            ui.ctx().request_repaint();
        }
    });
    let said = match (view.mode, view.there) {
        (Compute::Local, _) => "Runs on this computer.",
        (Compute::Workstation, _) => "Runs on the workstation.",
        (Compute::Automatic, true) => "Big: it runs on the workstation.",
        (Compute::Automatic, false) => "Small: it runs on this computer.",
    };
    let colour = if view.there { look.accent } else { look.dim };
    ui.label(RichText::new(said).size(12.0).color(colour));
}
