//! HelixLoom: a bioinformatics platform with a plasmid editor (enzyme maps, ORFs, features, a contig assembler) and
//! workspaces for sequencing signals, alignment, phylogeny, trees, protein structure and pipelines, and for sending the
//! heavy work to a workstation.
//!
//! Command line: `helixloom [files…] [--screenshot out.png [--no-window]] [--view <name>] [--dark]`. Files go to the
//! workspace they belong in.
//! With `--screenshot` the window is drawn for about 120 frames, saved as a PNG, and closed. When the
//! window cannot draw (the screen is locked, or there is no display) the same frames are drawn
//! without a window instead; `--no-window` does that from the start.
//!
//! On the workstation: `helixloom serve` takes the jobs the app sends and computes them (`helixloom help` lists it all).
//!
//! In a browser (the `wasm32` build, see `docs/WEB.md`) the same app starts on the page's canvas; the address takes
//! `?view=<name>` and `?dark` where the command line takes `--view` and `--dark`.

// The web build has no workstation server and no command line, so what only they use is not called there.
#![cfg_attr(target_arch = "wasm32", allow(dead_code))]

mod models;
mod pipeline;
mod remote;
mod sample;
mod seq;
mod ui;
#[cfg(target_arch = "wasm32")]
mod web;

#[cfg(not(target_arch = "wasm32"))]
use eframe::egui;
use ui::app::Startup;
#[cfg(not(target_arch = "wasm32"))]
use ui::{app::App, headless};

#[cfg(not(target_arch = "wasm32"))]
fn parse_args(args: &[String]) -> Startup {
    let mut start = Startup::default();
    let mut it = args.iter().skip(1);
    while let Some(a) = it.next() {
        match a.as_str() {
            "--screenshot" => start.screenshot = it.next().cloned(),
            "--view" => start.view = it.next().cloned().unwrap_or_default(),
            "--no-window" => start.no_window = true,
            "--dark" => start.dark = true,
            other if !other.starts_with("--") => start.files.push(other.to_string()),
            _ => {}
        }
    }
    start
}

/// What the page's address asks for, as the command line's options do (see `parse_args`): `?view=trees` opens a workspace
/// with its example, or any other screen `--view` knows, and `?dark` starts in the dark theme.
#[cfg(any(target_arch = "wasm32", test))]
fn parse_query(search: &str) -> Startup {
    let mut start = Startup::default();
    for pair in search.trim_start_matches('?').split('&') {
        match pair.split_once('=').unwrap_or((pair, "")) {
            ("view", name) => start.view = name.to_string(),
            ("dark", value) => start.dark = !matches!(value, "0" | "false"),
            _ => {}
        }
    }
    start
}

/// The browser runs `main` when the page has loaded the module; the app starts on the page's canvas.
#[cfg(target_arch = "wasm32")]
fn main() {
    web::run();
}

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result {
    let args: Vec<String> = std::env::args().collect();
    // the workstation's commands (serve, token, status, run-job) run without a window
    if let Some(code) = remote::cli::dispatch(&args) {
        std::process::exit(code);
    }
    let mut start = parse_args(&args);
    // a picture is made without the person's saved workstation (no calls over the network, no pixels that depend on it)
    // unless a folder of settings is named for it with HELIXLOOM_CONFIG_DIR
    let explicit = std::env::var_os("HELIXLOOM_CONFIG_DIR").is_some_and(|d| !d.is_empty());
    if start.screenshot.is_none() || explicit {
        start.settings_dir = remote::settings::dir();
    }
    if let Some(path) = &start.screenshot {
        if start.no_window {
            headless::finish(start.clone(), path);
        }
        headless::watch(start.clone(), path.clone());
    }
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size(headless::WINDOW_POINTS)
            .with_min_inner_size([760.0, 480.0])
            .with_title("HelixLoom"),
        ..Default::default()
    };
    let again = start.clone();
    let shown = eframe::run_native(
        "HelixLoom",
        options,
        Box::new(move |_cc| Ok(Box::new(App::new(start)))),
    );
    if let (Err(why), Some(path)) = (&shown, again.screenshot.clone()) {
        picture_without_window(again, &path, why);
    }
    shown
}

/// The window could not run: make the picture without it, or wait for the watchdog that is already
/// doing so (it ends the program once the picture is saved).
#[cfg(not(target_arch = "wasm32"))]
fn picture_without_window(start: Startup, path: &str, why: &eframe::Error) -> ! {
    if headless::claim_picture() {
        eprintln!("The window could not be opened ({why}); drawing the picture without it.");
        headless::finish(start, path);
    }
    loop {
        std::thread::park();
    }
}

#[cfg(all(test, not(target_arch = "wasm32")))]
mod tests {
    use super::*;

    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn the_page_address_asks_for_a_view_and_a_theme() {
        let s = parse_query("?view=trees&dark");
        assert_eq!((s.view.as_str(), s.dark), ("trees", true));
        let s = parse_query("?dark=0&other=1&view=digest");
        assert_eq!((s.view.as_str(), s.dark), ("digest", false));
        let s = parse_query("");
        assert_eq!((s.view.as_str(), s.dark), ("", false));
    }

    #[test]
    fn command_line_options_are_parsed() {
        let s = parse_args(&args(&[
            "helixloom",
            "a.gb",
            "--screenshot",
            "out.png",
            "--view",
            "digest",
            "b.fa",
        ]));
        assert_eq!(s.screenshot.as_deref(), Some("out.png"));
        assert_eq!(s.view, "digest");
        assert_eq!(s.files, vec!["a.gb", "b.fa"]);
        assert!(!s.no_window);
        assert!(parse_args(&args(&["helixloom"])).screenshot.is_none());
    }

    #[test]
    fn a_picture_without_a_window_is_asked_for_with_no_window() {
        let s = parse_args(&args(&[
            "helixloom",
            "--screenshot",
            "out.png",
            "--no-window",
        ]));
        assert!(s.no_window);
        assert_eq!(s.screenshot.as_deref(), Some("out.png"));
        assert!(s.files.is_empty());
    }
}
