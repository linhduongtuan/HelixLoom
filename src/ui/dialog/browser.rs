//! The browser's file dialogs: the page's own file chooser (through `rfd`'s async dialog), the files dropped on the page, and
//! downloads for saving. A browser answers later, and with names and bytes, so the answers wait in an inbox that the window
//! empties once a frame (`take_arrived`).

use super::{Arrived, Saved, Wanted};
use eframe::egui;
use std::cell::RefCell;
use std::path::PathBuf;
use wasm_bindgen::{closure::Closure, JsCast, JsValue};

/// How long a download's address is kept: the browser has been given the file by then.
const KEEP_DOWNLOAD_MS: i32 = 30_000;

thread_local! {
    /// What the browser brought and the window has not taken yet.
    static INBOX: RefCell<Vec<Arrived>> = const { RefCell::new(Vec::new()) };
    /// The window, to ask for a new frame when something is brought.
    static WINDOW: RefCell<Option<egui::Context>> = const { RefCell::new(None) };
}

/// What the desktop's file dialog is asked, kept so that the browser's chooser can be asked the same: a browser has no title
/// to show and no folders to choose, but it can offer only the file endings that were listed.
#[derive(Default)]
pub struct FileDialog {
    extensions: Vec<String>,
}

impl FileDialog {
    pub fn new() -> FileDialog {
        FileDialog::default()
    }

    pub fn set_title(self, _title: impl Into<String>) -> FileDialog {
        self
    }

    pub fn add_filter(
        mut self,
        _name: impl Into<String>,
        extensions: &[impl ToString],
    ) -> FileDialog {
        self.extensions
            .extend(extensions.iter().map(ToString::to_string));
        self
    }

    /// A browser cannot choose a folder.
    pub fn pick_folder(self) -> Option<PathBuf> {
        None
    }
}

/// Open the page's file chooser for one file. Nothing comes back now: the file arrives later, for `wanted`.
pub fn pick_one(wanted: Wanted, dialog: FileDialog) -> Option<PathBuf> {
    choose(wanted, dialog, false);
    None
}

/// The same for any number of files.
pub fn pick_many(wanted: Wanted, dialog: FileDialog) -> Vec<PathBuf> {
    choose(wanted, dialog, true);
    Vec::new()
}

fn choose(wanted: Wanted, dialog: FileDialog, many: bool) {
    let mut chooser = rfd::AsyncFileDialog::new();
    if !dialog.extensions.is_empty() {
        chooser = chooser.add_filter("Files", &dialog.extensions);
    }
    wasm_bindgen_futures::spawn_local(async move {
        let handles = if many {
            chooser.pick_files().await.unwrap_or_default()
        } else {
            chooser.pick_file().await.into_iter().collect()
        };
        let mut files = Vec::with_capacity(handles.len());
        for handle in &handles {
            files.push((handle.file_name(), Ok(handle.read().await)));
        }
        if !files.is_empty() {
            deliver(Arrived { wanted, files });
        }
    });
}

/// Read the files dropped on the page, one after the other, and bring them for `wanted`. They have no path to give.
pub fn dropped(files: &[egui::DroppedFileHandle], wanted: Wanted) -> Vec<PathBuf> {
    let files = files.to_vec();
    wasm_bindgen_futures::spawn_local(async move {
        let mut read = Vec::with_capacity(files.len());
        for file in &files {
            let name = file.path().to_string_lossy().into_owned();
            read.push((name, file.bytes_async().await));
        }
        deliver(Arrived {
            wanted,
            files: read,
        });
    });
    Vec::new()
}

/// What the browser brought since the last frame. The window is remembered, to be woken when more comes.
pub fn take_arrived(ctx: &egui::Context) -> Vec<Arrived> {
    WINDOW.with_borrow_mut(|window| {
        window.get_or_insert_with(|| ctx.clone());
    });
    INBOX.with_borrow_mut(std::mem::take)
}

fn deliver(arrived: Arrived) {
    INBOX.with_borrow_mut(|inbox| inbox.push(arrived));
    WINDOW.with_borrow(|window| {
        if let Some(ctx) = window {
            ctx.request_repaint();
        }
    });
}

/// Hand `bytes` to the browser as a download named `name`; there is no dialog, the browser decides where it goes.
pub fn save_bytes(_dialog: FileDialog, name: &str, bytes: &[u8]) -> Option<Saved> {
    Some(match download(name, bytes) {
        Ok(()) => Saved {
            message: format!("Downloaded {name}."),
            ok: true,
        },
        Err(why) => Saved {
            message: format!("Could not download {name}: {why}"),
            ok: false,
        },
    })
}

fn describe(error: JsValue) -> String {
    error.as_string().unwrap_or_else(|| format!("{error:?}"))
}

/// A hidden link to a blob of the bytes, clicked once.
fn download(name: &str, bytes: &[u8]) -> Result<(), String> {
    let window = web_sys::window().ok_or("there is no window")?;
    let document = window.document().ok_or("there is no page")?;
    let body = document.body().ok_or("the page has no body")?;
    let parts = js_sys::Array::of1(&js_sys::Uint8Array::from(bytes));
    let kind = web_sys::BlobPropertyBag::new();
    kind.set_type("application/octet-stream");
    let blob =
        web_sys::Blob::new_with_u8_array_sequence_and_options(&parts, &kind).map_err(describe)?;
    let url = web_sys::Url::create_object_url_with_blob(&blob).map_err(describe)?;
    let link = document
        .create_element("a")
        .map_err(describe)?
        .dyn_into::<web_sys::HtmlAnchorElement>()
        .map_err(|_| "the page made no link")?;
    link.set_href(&url);
    link.set_download(name);
    link.set_hidden(true);
    body.append_child(&link).map_err(describe)?;
    link.click();
    link.remove();
    release_later(&window, url);
    Ok(())
}

/// Let go of the download's address after a while.
fn release_later(window: &web_sys::Window, url: String) {
    let release = Closure::once_into_js(move || {
        if let Err(why) = web_sys::Url::revoke_object_url(&url) {
            log::warn!("could not release a download: {}", describe(why));
        }
    });
    let asked = window.set_timeout_with_callback_and_timeout_and_arguments_0(
        release.unchecked_ref(),
        KEEP_DOWNLOAD_MS,
    );
    if let Err(why) = asked {
        log::warn!("could not schedule releasing a download: {}", describe(why));
    }
}
