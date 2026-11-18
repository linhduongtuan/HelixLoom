//! The web entry point: the same app, started by `eframe::WebRunner` on the canvas of `index.html`.

use crate::ui::app::App;
use wasm_bindgen::JsCast as _;

/// The canvas the app draws on, in `index.html`.
const CANVAS_ID: &str = "helixloom_canvas";
/// The text that says "Loading…" until the app runs; it takes the place of the page's message if the app cannot start.
const STATUS_ID: &str = "loading";

/// Start the app on the page. The start is asynchronous (the graphics are set up first), so this returns at once.
pub fn run() {
    // what eframe says through `log` goes to the browser's console
    let _ = eframe::WebLogger::init(log::LevelFilter::Info);
    wasm_bindgen_futures::spawn_local(async {
        let outcome = start().await;
        let status = web_sys::window()
            .and_then(|window| window.document())
            .and_then(|document| document.get_element_by_id(STATUS_ID));
        match (outcome, status) {
            (Ok(()), Some(status)) => status.remove(),
            (Ok(()), None) => {}
            (Err(why), status) => {
                log::error!("HelixLoom could not start: {why}");
                if let Some(status) = status {
                    status.set_text_content(Some(&format!("HelixLoom could not start: {why}")));
                }
            }
        }
    });
}

async fn start() -> Result<(), String> {
    let window = web_sys::window().ok_or("there is no browser window")?;
    let document = window.document().ok_or("there is no page")?;
    let canvas = document
        .get_element_by_id(CANVAS_ID)
        .ok_or_else(|| format!("the page has no canvas called {CANVAS_ID}"))?
        .dyn_into::<web_sys::HtmlCanvasElement>()
        .map_err(|_| format!("{CANVAS_ID} is not a canvas"))?;
    // the address may ask for dark (`?dark`), for this visit; otherwise the theme is the one that is kept, which follows the
    // browser's own light or dark setting until the person chooses (a visit that is given a theme keeps none, as a run that is
    // given one on the command line: so the browser's setting is not put here as if it had been asked for)
    let startup = crate::parse_query(&window.location().search().unwrap_or_default());
    // the database that keeps the results of the pages is read whole before the app starts (it is asked at once for them then)
    crate::ui::session::open_database().await;
    let runner = eframe::WebRunner::new();
    runner
        .start(
            canvas,
            eframe::WebOptions::default(),
            Box::new(move |cc| {
                crate::ui::fonts::install(&cc.egui_ctx);
                crate::ui::session::wake_on_refusal(&cc.egui_ctx);
                Ok(Box::new(App::new(startup)))
            }),
        )
        .await
        .map_err(|why| format!("{why:?}"))?;
    // The open sequences are kept a moment after the last change, from a frame (`ui/autosave.rs`), and a page that is closed
    // or reloaded draws no more frames: what is waiting is written when the page is left, and when it is hidden (a tab that
    // is hidden may be closed without another word).
    let leaving: [(&web_sys::EventTarget, &'static str); 2] = [
        (window.as_ref(), "pagehide"),
        (document.as_ref(), "visibilitychange"),
    ];
    for (target, event) in leaving {
        let listening = runner.add_event_listener(target, event, |_: web_sys::Event, runner| {
            runner.app_mut::<App>().flush_session();
        });
        // the app runs without it: the last seconds of work are then kept only by the frames
        if let Err(why) = listening {
            log::warn!("HelixLoom cannot tell when the page is left ({event}): {why:?}");
        }
    }
    Ok(())
}
