//! Help ▸ Check for Updates…: asks GitHub for the latest release when the person presses the button (never before), says whether a later
//! version is out, and downloads the file for this computer, checked against the release's checksums (`crate::update`). The program does
//! not replace itself: the file is saved and shown, and what to do with it is said. In a browser the page is the program, and it updates
//! when it is loaded again.

use crate::ui::app::App;
use crate::ui::design;
use crate::update::{Downloaded, Outcome, Release};
use eframe::egui;
#[cfg(not(target_arch = "wasm32"))]
use eframe::egui::RichText;

#[derive(Clone, Debug, PartialEq, Default)]
pub enum State {
    #[default]
    Idle,
    Checking,
    UpToDate(String),
    Available(Release),
    Failed(String),
    Downloading(Release),
    Downloaded(Release, Downloaded),
}

/// What a piece of work that ran on a thread of its own came to.
pub enum Msg {
    Checked(Result<Outcome, String>),
    Downloaded(Result<Downloaded, String>),
}

impl State {
    /// The state after `msg` came back; a message that no work of the present state is waiting for changes nothing.
    pub fn after(self, msg: Msg) -> State {
        match (self, msg) {
            (State::Checking, Msg::Checked(Ok(Outcome::UpToDate { latest }))) => {
                State::UpToDate(latest)
            }
            (State::Checking, Msg::Checked(Ok(Outcome::Available(release)))) => {
                State::Available(release)
            }
            (State::Checking, Msg::Checked(Err(why))) => State::Failed(why),
            (State::Downloading(release), Msg::Downloaded(Ok(file))) => {
                State::Downloaded(release, file)
            }
            (State::Downloading(_), Msg::Downloaded(Err(why))) => State::Failed(why),
            (state, _) => state,
        }
    }
}

#[derive(Default)]
pub struct UpdateWin {
    pub open: bool,
    pub state: State,
    #[cfg(not(target_arch = "wasm32"))]
    job: Option<crate::ui::job::Job<Msg>>,
}

/// What a press in the window asks for.
enum Press {
    Check,
    Download(Release),
    Page(String),
    Show(std::path::PathBuf),
    #[cfg(target_arch = "wasm32")]
    Reload,
}

fn megabytes(bytes: u64) -> String {
    format!("{:.1} MB", bytes as f64 / 1_000_000.0)
}

impl App {
    /// Collect what the thread of a check or a download has to say.
    #[cfg(not(target_arch = "wasm32"))]
    fn poll_update(&mut self) {
        use crate::ui::job::Poll;
        let Some(job) = &self.update_win.job else {
            return;
        };
        match job.poll() {
            Poll::Running => {}
            Poll::Done(msg) => {
                self.update_win.job = None;
                self.update_win.state = std::mem::take(&mut self.update_win.state).after(msg);
            }
            Poll::Died => {
                self.update_win.job = None;
                self.update_win.state = State::Failed("The work stopped without an answer.".into());
            }
        }
    }

    #[cfg(not(target_arch = "wasm32"))]
    fn start_update_work(&mut self, ctx: &egui::Context, press: Press) {
        use crate::connect::curl::Curl;
        use crate::ui::job::Job;
        use crate::update;
        let Some(curl) = crate::pipeline::tools::locate(&["curl"]) else {
            self.update_win.state = State::Failed(
                "curl was not found on this computer: it is what asks GitHub.".into(),
            );
            return;
        };
        match press {
            Press::Check => {
                self.update_win.state = State::Checking;
                self.update_win.job = Some(Job::spawn(ctx, move || {
                    Msg::Checked(update::check(&Curl::new(curl), env!("CARGO_PKG_VERSION")))
                }));
            }
            Press::Download(release) => {
                let Some(asset) =
                    update::asset_for(&release, std::env::consts::OS, std::env::consts::ARCH)
                        .cloned()
                else {
                    return;
                };
                self.update_win.state = State::Downloading(release.clone());
                self.update_win.job = Some(Job::spawn(ctx, move || {
                    let got = update::download(
                        &release,
                        &asset,
                        &update::downloads_dir(),
                        &|url, file| update::curl_file(&curl, url, file),
                        &Curl::new(curl.clone()),
                    );
                    Msg::Downloaded(got)
                }));
            }
            _ => {}
        }
    }

    pub fn update_window(&mut self, ctx: &egui::Context) {
        if !self.update_win.open {
            return;
        }
        #[cfg(not(target_arch = "wasm32"))]
        self.poll_update();
        let (mut open, mut press) = (true, None::<Press>);
        let state = self.update_win.state.clone();
        design::window("Check for Updates")
            .open(&mut open)
            .default_width(480.0)
            .show(ctx, |ui| {
                ui.label(format!("This is HelixLoom {}.", env!("CARGO_PKG_VERSION")));
                #[cfg(target_arch = "wasm32")]
                {
                    ui.label("This is the web version: it updates when the page is loaded again.");
                    if ui.button("Reload the page").clicked() {
                        press = Some(Press::Reload);
                    }
                    let _ = &state;
                }
                #[cfg(not(target_arch = "wasm32"))]
                {
                    ui.label(
                        RichText::new("Nothing is asked of any server until you press the button; then one request goes to GitHub.")
                            .small(),
                    );
                    let busy = matches!(state, State::Checking | State::Downloading(_));
                    if ui.add_enabled(!busy, egui::Button::new("Check for updates")).clicked() {
                        press = Some(Press::Check);
                    }
                    ui.separator();
                    match &state {
                        State::Idle => {}
                        State::Checking => {
                            ui.label("Asking GitHub…");
                        }
                        State::UpToDate(latest) => {
                            ui.label(format!("This is the latest version (the latest release is {latest})."));
                        }
                        State::Failed(why) => {
                            ui.label(RichText::new(format!("Could not do that: {why}")).color(ui.visuals().error_fg_color));
                        }
                        State::Available(release) | State::Downloading(release) | State::Downloaded(release, _) => {
                            ui.label(RichText::new(format!("A later version is out: {}.", release.tag)).strong());
                            egui::ScrollArea::vertical().max_height(140.0).show(ui, |ui| {
                                ui.label(RichText::new(&release.notes).small());
                            });
                            let asset = crate::update::asset_for(release, std::env::consts::OS, std::env::consts::ARCH);
                            ui.horizontal(|ui| {
                                match asset {
                                    Some(a) => {
                                        let label = format!("Download {} ({})", a.name, megabytes(a.size));
                                        let can = matches!(state, State::Available(_));
                                        if ui.add_enabled(can, egui::Button::new(label)).clicked() {
                                            press = Some(Press::Download(release.clone()));
                                        }
                                    }
                                    None => {
                                        ui.label("This release has no file for this computer's system: build it from the source (docs/UBUNTU.md).");
                                    }
                                }
                                if ui.button("Open the release page").clicked() {
                                    press = Some(Press::Page(release.page.clone()));
                                }
                            });
                            match &state {
                                State::Downloading(_) => {
                                    ui.label("Downloading…");
                                }
                                State::Downloaded(_, got) => {
                                    ui.label(format!("Saved {}", got.path.display()));
                                    ui.label(if got.verified {
                                        "Its checksum matches the one the release publishes.".to_string()
                                    } else {
                                        "The release has no checksum for it, so it could not be checked.".to_string()
                                    });
                                    ui.label(crate::update::next_step(&got.path));
                                    if ui.button("Show it in its folder").clicked() {
                                        press = Some(Press::Show(got.path.clone()));
                                    }
                                }
                                _ => {}
                            }
                        }
                    }
                }
            });
        match press {
            #[cfg(not(target_arch = "wasm32"))]
            Some(p @ (Press::Check | Press::Download(_))) => self.start_update_work(ctx, p),
            Some(Press::Page(url)) => ctx.open_url(egui::OpenUrl::new_tab(url)),
            Some(Press::Show(path)) => crate::update::reveal(&path),
            #[cfg(target_arch = "wasm32")]
            Some(Press::Reload) => {
                if let Some(w) = web_sys::window() {
                    let _ = w.location().reload();
                }
            }
            _ => {}
        }
        self.update_win.open = open;
    }
}
