//! The page's way to the workstation: one request through the browser's `fetch`, and a wait that does not stop the page.
//! The program reaches the server over a plain TCP connection (`http::call`), which a page cannot open; the server speaks
//! HTTP, so the same requests, with the same bearer token, go through `fetch` instead.

use super::http::{Allowed, Reply, MIN_UPLOAD_RATE};
use js_sys::{Promise, Uint8Array};
use std::time::Duration;
use wasm_bindgen::{closure::Closure, JsCast, JsValue};
use wasm_bindgen_futures::JsFuture;
use web_sys::{AbortController, Headers, Request, RequestInit, RequestMode, Response, Window};

/// The words for what the browser threw.
fn said(error: &JsValue) -> String {
    error
        .dyn_ref::<js_sys::Error>()
        .map(|e| String::from(e.message()))
        .or_else(|| error.as_string())
        .unwrap_or_else(|| "the browser gave no reason".into())
}

/// Whether the browser threw because the request was stopped (by the deadline).
fn was_stopped(error: &JsValue) -> bool {
    error
        .dyn_ref::<js_sys::Error>()
        .is_some_and(|e| e.name() == "AbortError")
}

/// Stop `abort` after `time`, since a request has no time limit of its own. The timer to cancel when the request is over,
/// with the closure it calls: the closure is freed when it is dropped, so it must be kept until the timer is cancelled (one
/// that is handed over to the browser is freed only if it is called, and this one is nearly always cancelled first).
fn set_deadline(
    window: &Window,
    abort: &AbortController,
    time: Duration,
) -> Option<(i32, Closure<dyn FnMut()>)> {
    let abort = abort.clone();
    let stop = Closure::once(move || abort.abort());
    let ms = i32::try_from(time.as_millis()).unwrap_or(i32::MAX);
    let timer = window
        .set_timeout_with_callback_and_timeout_and_arguments_0(stop.as_ref().unchecked_ref(), ms)
        .ok()?;
    Some((timer, stop))
}

/// Send one request to `base` (`http://host:7878`, no slash at the end) and read the whole reply, as long as it is no
/// larger than `allowed` says and arrives in the time it says.
pub async fn call(
    base: &str,
    method: &str,
    target: &str,
    token: &str,
    body: &[u8],
    allowed: Allowed,
) -> Result<Reply, String> {
    let window = web_sys::window().ok_or("There is no browser window.")?;
    let abort = AbortController::new().map_err(|e| said(&e))?;
    // the time `allowed` gives is for the reply, as in a program (where it starts once the body is written): an upload has the
    // time its size needs at the slowest rate the server waits for, besides
    let time = allowed.time + Duration::from_secs((body.len() / MIN_UPLOAD_RATE) as u64);
    let deadline = set_deadline(&window, &abort, time);
    let outcome = request(&window, &abort, base, method, target, token, body, allowed).await;
    if let Some((timer, _stop)) = deadline {
        window.clear_timeout_with_handle(timer);
    }
    outcome
}

#[allow(clippy::too_many_arguments)]
async fn request(
    window: &Window,
    abort: &AbortController,
    base: &str,
    method: &str,
    target: &str,
    token: &str,
    body: &[u8],
    allowed: Allowed,
) -> Result<Reply, String> {
    let thrown = |e: JsValue| said(&e);
    let headers = Headers::new().map_err(thrown)?;
    headers
        .set("Authorization", &format!("Bearer {token}"))
        .map_err(thrown)?;
    headers
        .set("Content-Type", "application/x-helixloom")
        .map_err(thrown)?;
    let init = RequestInit::new();
    init.set_method(method);
    init.set_headers(&headers);
    init.set_mode(RequestMode::Cors);
    init.set_signal(Some(&abort.signal()));
    if !body.is_empty() {
        init.set_body(&Uint8Array::from(body));
    }
    let request =
        Request::new_with_str_and_init(&format!("{base}{target}"), &init).map_err(thrown)?;
    let response = JsFuture::from(window.fetch_with_request(&request))
        .await
        .map_err(|e| {
            if was_stopped(&e) {
                format!("{base} did not answer in time.")
            } else {
                format!(
                    "Cannot reach {base} ({}). Check that the server is running and the address is right; a page opened \
                     over https cannot call an http address, and a page that the server does not serve itself must be \
                     named to it with --allow-origin.",
                    said(&e)
                )
            }
        })?;
    let response: Response = response
        .dyn_into()
        .map_err(|_| format!("{base} did not answer with a reply."))?;
    let too_large = || {
        format!(
            "The reply of {base} is larger than this program takes ({} MB).",
            allowed.bytes >> 20
        )
    };
    // the size the server announces, before any of the reply is read
    // ponytail: a reply that announces no length (a proxy that compresses or chunks it) is read whole before the size is
    // checked, where a program stops at the limit; only a broken or hostile server sends none. Upgrade: read the body as a
    // stream and count (`response.body().get_reader()`).
    let announced = response
        .headers()
        .get("content-length")
        .ok()
        .flatten()
        .and_then(|v| v.parse::<usize>().ok());
    if announced.is_some_and(|n| n > allowed.bytes) {
        // stop the download too: the browser would go on fetching the rest
        abort.abort();
        return Err(too_large());
    }
    let buffer = JsFuture::from(response.array_buffer().map_err(thrown)?)
        .await
        .map_err(|e| {
            if was_stopped(&e) {
                format!("The reply of {base} took too long.")
            } else {
                format!("The reply of {base} was cut short ({}).", said(&e))
            }
        })?;
    let bytes = Uint8Array::new(&buffer).to_vec();
    if bytes.len() > allowed.bytes {
        return Err(too_large());
    }
    Ok(Reply {
        status: response.status(),
        body: bytes,
    })
}

/// Wait `time` without stopping the page: the timer of the browser, not a loop.
pub async fn sleep(time: Duration) {
    let ms = i32::try_from(time.as_millis()).unwrap_or(i32::MAX);
    let waited = Promise::new(&mut |resolve, _| match web_sys::window() {
        Some(window) => {
            let _ = window.set_timeout_with_callback_and_timeout_and_arguments_0(&resolve, ms);
        }
        // no window to keep time: do not wait for ever
        None => {
            let _ = resolve.call0(&JsValue::NULL);
        }
    });
    let _ = JsFuture::from(waited).await;
}
