//! The sheet of keyboard shortcuts: every command that has a key, in its group, as it is written in the command palette.

use super::app::App;
use super::commands::{self, Command};
use super::design::{self, Look};
use eframe::egui::{self, Context, RichText, Ui};

/// How wide the sheet is.
const WIDTH: f32 = 560.0;

/// The commands that have a key, grouped, in the order the groups first appear.
fn with_keys(commands: &[Command]) -> Vec<(&'static str, Vec<&Command>)> {
    let mut groups: Vec<(&'static str, Vec<&Command>)> = Vec::new();
    for command in commands.iter().filter(|c| !c.keys.is_empty()) {
        match groups.iter_mut().find(|(group, _)| *group == command.group) {
            Some((_, list)) => list.push(command),
            None => groups.push((command.group, vec![command])),
        }
    }
    groups
}

impl App {
    pub(super) fn shortcuts_window(&mut self, ctx: &Context) {
        if !self.shortcuts_open {
            return;
        }
        let look = design::look(self.prefs.dark);
        let mut open = true;
        crate::ui::design::window("Keyboard shortcuts")
            .open(&mut open)
            .default_width(WIDTH)
            .collapsible(false)
            .show(ctx, |ui| {
                ui.label(
                    RichText::new(format!(
                        "Press {} K to search every command, so that none of these has to be remembered.",
                        design::command_key()
                    ))
                    .color(look.dim),
                );
                ui.add_space(6.0);
                egui::ScrollArea::vertical().max_height(460.0).show(ui, |ui| {
                    for (group, list) in with_keys(commands::all()) {
                        sheet_group(ui, &look, group, &list);
                    }
                });
            });
        self.shortcuts_open = open;
    }
}

fn sheet_group(ui: &mut Ui, look: &Look, group: &str, list: &[&Command]) {
    design::section(ui, look, group);
    egui::Grid::new(("shortcut_sheet", group))
        .num_columns(2)
        .spacing([24.0, 6.0])
        .show(ui, |ui| {
            for command in list {
                ui.label(command.title);
                ui.horizontal(|ui| {
                    design::keycap(
                        ui,
                        look,
                        &command.keys.replace("Cmd", design::command_key()),
                    );
                });
                ui.end_row();
            }
        });
    ui.add_space(4.0);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_sheet_has_a_group_for_each_kind_of_command_that_has_a_key_and_none_without() {
        let groups = with_keys(commands::all());
        let names: Vec<&str> = groups.iter().map(|(g, _)| *g).collect();
        assert_eq!(names[0], "Go to");
        assert!(names.contains(&"File") && names.contains(&"Tools"));
        assert!(groups
            .iter()
            .all(|(_, list)| list.iter().all(|c| !c.keys.is_empty())));
        // a command without a key is not on the sheet
        assert!(!groups
            .iter()
            .flat_map(|(_, list)| list)
            .any(|c| c.title == "Export as FASTA…"));
    }
}
