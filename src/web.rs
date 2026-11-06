//! The web entry point: the same app, started by `eframe::WebRunner` on the canvas of `index.html`.

use crate::ui::app::App;
use eframe::egui;
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
    let mut startup = crate::parse_query(&window.location().search().unwrap_or_default());
    // the browser's own light or dark setting, unless the address asked for dark
    startup.dark |= eframe::web::system_theme() == Some(egui::Theme::Dark);
    eframe::WebRunner::new()
        .start(
            canvas,
            eframe::WebOptions::default(),
            Box::new(move |_cc| Ok(Box::new(App::new(startup)))),
        )
        .await
        .map_err(|why| format!("{why:?}"))
}
