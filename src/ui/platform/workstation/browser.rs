//! The Workstation page in a browser. A page cannot open the plain TCP connection the workstation server speaks, so the page
//! says that plainly; the roadmap reaches the workstation through the browser's fetch (docs/WEB.md).

use eframe::egui::Ui;

pub(super) fn not_available(ui: &mut Ui) {
    ui.heading("Workstation");
    ui.add_space(8.0);
    ui.label("Connecting to a workstation is not available in the browser yet.");
    ui.add_space(4.0);
    ui.label(
        "Every analysis in the browser runs on this computer, inside the page. To send heavy jobs to a computer with \
         GPUs, use the desktop app for now; the browser version will reach a workstation through the browser's fetch.",
    );
}
