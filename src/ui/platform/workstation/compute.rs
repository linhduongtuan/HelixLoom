//! Where a job runs: the choice of this computer, the workstation, or each job by its size, and the small control that every
//! page that starts a job shows beside its button.

use super::dim;
use crate::remote::client::Remote;
use crate::remote::settings::Compute;
use crate::remote::weight::{
    HEAVY_GENERATIONS, HEAVY_LETTERS, HEAVY_REPLICATES, HEAVY_SEQUENCES, HEAVY_TREES,
};
use crate::ui::design::Look;
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
