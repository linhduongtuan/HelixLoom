//! Work that takes seconds (the GPU's first run, MAFFT, IQ-TREE) runs on its own thread, so the window keeps
//! drawing while it goes. The window asks `poll` once a frame.
//!
//! A browser page has no threads to start (yet), so there the work runs at once, on the window's own thread, and is done
//! by the first `poll`: the window waits for it, and what takes long makes the page stall.

use eframe::egui;
use std::sync::mpsc::{self, Receiver, TryRecvError};

pub struct Job<T> {
    results: Receiver<T>,
}

pub enum Poll<T> {
    Running,
    Done(T),
    /// The work stopped without an answer (it panicked).
    Died,
}

impl<T: Send + 'static> Job<T> {
    /// Start `work` on a new thread; the window is asked to redraw when it is done.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn spawn(ctx: &egui::Context, work: impl FnOnce() -> T + Send + 'static) -> Job<T> {
        let (send, results) = mpsc::channel();
        let ctx = ctx.clone();
        std::thread::spawn(move || {
            // a panic drops `send` unused, which `poll` reports as `Died`
            let _ = send.send(work());
            ctx.request_repaint();
        });
        Job { results }
    }

    /// Run `work` now: there are no threads in a browser page. The answer is there for the first `poll`.
    #[cfg(target_arch = "wasm32")]
    pub fn spawn(ctx: &egui::Context, work: impl FnOnce() -> T + Send + 'static) -> Job<T> {
        let (send, results) = mpsc::channel();
        // the receiver is held by the job, so sending cannot fail
        let _ = send.send(work());
        ctx.request_repaint();
        Job { results }
    }

    pub fn poll(&self) -> Poll<T> {
        match self.results.try_recv() {
            Ok(done) => Poll::Done(done),
            Err(TryRecvError::Empty) => Poll::Running,
            Err(TryRecvError::Disconnected) => Poll::Died,
        }
    }
}
