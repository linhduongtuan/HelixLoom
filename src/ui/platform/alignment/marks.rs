//! Stretches of columns marked on the alignment, for the primers and the probe of a design: a bar at the foot of the header and an
//! outline over the rows.

use super::AlignmentPage;
use crate::ui::design::Look;
use eframe::egui::Color32;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    /// A primer that reads along the alignment (a left primer).
    Left,
    /// A primer on the other strand (a right primer).
    Right,
    Probe,
}

impl Role {
    pub fn colour(self, look: &Look) -> Color32 {
        match self {
            Role::Left => look.good,
            Role::Right => look.warn,
            Role::Probe => look.accent2,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Role::Left => "left primer",
            Role::Right => "right primer",
            Role::Probe => "probe",
        }
    }
}

/// Columns `first` to `last`, both included, counted from 0.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Mark {
    pub first: usize,
    pub last: usize,
    pub role: Role,
}

impl AlignmentPage {
    /// Mark these stretches, and select the first column of the first so that the viewer scrolls to it.
    pub fn set_marks(&mut self, marks: Vec<Mark>) {
        let first = marks.iter().map(|m| m.first).min();
        self.marks = marks;
        if first.is_some() {
            self.select_column(first, true);
        }
    }

    #[cfg(test)]
    pub fn marks(&self) -> &[Mark] {
        &self.marks
    }
}
