use super::*;
use crate::ui::app::Startup;
use crate::ui::headless::frame_input;

fn source(key: &'static str) -> Source {
    Source {
        key,
        label: key,
        open: Action::Go(Workspace::Alignment),
    }
}

#[test]
fn a_span_is_said_the_way_people_say_it() {
    for (seconds, said) in [
        (0, "under a second"),
        (1, "1 s"),
        (59, "59 s"),
        (60, "1 min"),
        (125, "2 min 5 s"),
        (3599, "59 min 59 s"),
        (3600, "1 h"),
        (3780, "1 h 3 min"),
    ] {
        assert_eq!(took(Duration::from_secs(seconds)), said, "{seconds}");
    }
    assert_eq!(took(Duration::from_millis(900)), "under a second");
}

#[test]
fn a_job_is_seen_to_start_to_go_on_and_to_end_with_the_time_it_took_and_the_newest_end_is_first() {
    let mut tray = Tray::default();
    let t0 = Instant::now();
    let at = |s: u64| t0 + Duration::from_secs(s);
    assert!(tray.observe(&[], at(0)).is_empty());
    assert!(tray.observe(&[source("a")], at(1)).is_empty());
    assert!(tray.observe(&[source("a"), source("b")], at(5)).is_empty());
    assert_eq!(tray.running().len(), 2);
    // going on is not news, and the clock of a job runs from when it was first seen
    assert!(tray.observe(&[source("a"), source("b")], at(9)).is_empty());
    assert_eq!(tray.running()[0].started, at(1));
    // a ends
    let ended = tray.observe(&[source("b")], at(11));
    assert_eq!(ended.len(), 1);
    assert_eq!(
        (ended[0].source.key, ended[0].took),
        ("a", Duration::from_secs(10))
    );
    assert_eq!(tray.running().len(), 1);
    // a starts again: it is a new job, not the old one going on
    tray.observe(&[source("a"), source("b")], at(12));
    assert_eq!(
        tray.running()
            .iter()
            .find(|r| r.source.key == "a")
            .map(|r| r.started),
        Some(at(12))
    );
    // both end at once: both are told, and the history has the newest first
    let ended = tray.observe(&[], at(20));
    assert_eq!(ended.len(), 2);
    assert_eq!(tray.history().len(), 3);
    assert!(tray.running().is_empty());
    assert_eq!(tray.history()[2].source.key, "a");
    assert_eq!(tray.history()[2].took, Duration::from_secs(10));
}

#[test]
fn the_tray_remembers_only_the_last_few_that_ended() {
    let mut tray = Tray::default();
    let t0 = Instant::now();
    for k in 0..REMEMBERED + 5 {
        let now = t0 + Duration::from_secs(k as u64 * 2);
        tray.observe(&[source("a")], now);
        tray.observe(&[], now + Duration::from_secs(1));
    }
    assert_eq!(tray.history().len(), REMEMBERED);
    assert!(tray.history().windows(2).all(|w| w[0].ended >= w[1].ended));
}

fn frame(app: &mut App, ctx: &egui::Context, n: u32, events: Vec<egui::Event>) {
    let mut input = frame_input(n);
    input.events = events;
    let mut output = ctx.run_ui(input, |ui| app.draw(ui));
    output.textures_delta.clear();
}

#[test]
fn a_search_that_runs_shows_in_the_tray_and_its_end_is_a_toast_that_names_it() {
    let ctx = egui::Context::default();
    let mut app = App::new(Startup::default());
    app.motif_win.text = crate::sample::motif_demo();
    app.motif_win.settings.width = 10;
    app.motif_win.open = true;
    app.motif_win.start(&ctx);
    assert!(app.running_sources().iter().any(|s| s.key == "motifs"));
    let mut n = 0;
    let until = Instant::now() + Duration::from_secs(120);
    while app.motif_win.running() && Instant::now() < until {
        frame(&mut app, &ctx, n, Vec::new());
        n += 1;
        std::thread::sleep(Duration::from_millis(10));
    }
    // one more frame, to see that it is over
    frame(&mut app, &ctx, n, Vec::new());
    assert!(!app.motif_win.running());
    let done = app.tray.history().first().expect("it ended");
    assert_eq!(done.source.key, "motifs");
    assert!(
        app.toasts
            .texts()
            .iter()
            .any(|t| t.starts_with("Motif search ended after ")),
        "{:?}",
        app.toasts.texts()
    );
    assert!(app.tray.running().is_empty());
}

#[test]
fn the_panel_opens_from_the_bar_shows_what_runs_and_what_ran_and_closes_with_escape() {
    let ctx = egui::Context::default();
    let mut app = App::new(Startup::default());
    // nothing has ever run: no pill, and an open panel says so
    frame(&mut app, &ctx, 0, Vec::new());
    app.run(Action::JobTray, &ctx);
    assert!(app.tray.open);
    frame(&mut app, &ctx, 1, Vec::new());
    // a search runs: the panel is drawn with it in the list
    app.motif_win.text = crate::sample::motif_demo();
    app.motif_win.settings.width = 10;
    app.motif_win.settings.seeds = 200;
    app.motif_win.start(&ctx);
    frame(&mut app, &ctx, 2, Vec::new());
    assert_eq!(app.tray.running().len(), 1);
    // and when it has ended, in the list of what ended
    let until = Instant::now() + Duration::from_secs(120);
    let mut n = 3;
    while app.motif_win.running() && Instant::now() < until {
        frame(&mut app, &ctx, n, Vec::new());
        n += 1;
        std::thread::sleep(Duration::from_millis(10));
    }
    frame(&mut app, &ctx, n, Vec::new());
    assert_eq!((app.tray.running().len(), app.tray.history().len()), (0, 1));
    assert!(app.tray.open, "the panel stays open while it is looked at");
    // Escape closes it
    let escape = egui::Event::Key {
        key: egui::Key::Escape,
        physical_key: Some(egui::Key::Escape),
        pressed: true,
        repeat: false,
        modifiers: egui::Modifiers::NONE,
    };
    frame(&mut app, &ctx, n + 1, vec![escape]);
    assert!(!app.tray.open);
    // and the bar's action toggles it
    app.run(Action::JobTray, &ctx);
    assert!(app.tray.open);
    app.run(Action::JobTray, &ctx);
    assert!(!app.tray.open);
}

fn press(at: egui::Pos2) -> Vec<egui::Event> {
    vec![
        egui::Event::PointerMoved(at),
        egui::Event::PointerButton {
            pos: at,
            button: egui::PointerButton::Primary,
            pressed: true,
            modifiers: egui::Modifiers::NONE,
        },
    ]
}

#[test]
fn a_click_inside_the_panel_leaves_it_open_and_a_click_anywhere_else_closes_it() {
    let ctx = egui::Context::default();
    let mut app = App::new(Startup::default());
    let width = frame_input(0).screen_rect.expect("a screen").width();
    app.run(Action::JobTray, &ctx);
    // the panel is there for a frame, with the pointer over it, before the click
    let inside = egui::pos2(width - 100.0, topbar::HEIGHT + 20.0);
    frame(&mut app, &ctx, 0, vec![egui::Event::PointerMoved(inside)]);
    frame(&mut app, &ctx, 1, vec![egui::Event::PointerMoved(inside)]);
    frame(&mut app, &ctx, 2, press(inside));
    assert!(app.tray.open, "a click inside the panel");
    let outside = egui::pos2(60.0, 500.0);
    frame(&mut app, &ctx, 3, vec![egui::Event::PointerMoved(outside)]);
    frame(&mut app, &ctx, 4, press(outside));
    assert!(!app.tray.open, "a click elsewhere");
}

#[test]
fn what_ended_long_ago_is_forgotten_and_the_bar_stops_offering_it() {
    let mut tray = Tray::default();
    let t0 = Instant::now();
    tray.observe(&[source("a")], t0);
    tray.observe(&[], t0 + Duration::from_secs(5));
    assert_eq!(tray.history().len(), 1);
    // a while later, with nothing running, a look is enough to drop it
    tray.observe(&[], t0 + REMEMBERED_FOR + Duration::from_secs(10));
    assert!(tray.history().is_empty());
}

#[test]
fn escape_closes_the_panel_but_when_the_palette_is_open_it_is_the_palettes() {
    let ctx = egui::Context::default();
    let mut app = App::new(Startup::default());
    let escape = || {
        vec![egui::Event::Key {
            key: egui::Key::Escape,
            physical_key: Some(egui::Key::Escape),
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        }]
    };
    app.run(Action::JobTray, &ctx);
    frame(&mut app, &ctx, 0, Vec::new());
    app.palette.toggle();
    assert!(app.palette.is_open());
    frame(&mut app, &ctx, 1, escape());
    assert!(app.tray.open, "the palette took the key");
    assert!(!app.palette.is_open());
    frame(&mut app, &ctx, 2, escape());
    assert!(!app.tray.open);
}
