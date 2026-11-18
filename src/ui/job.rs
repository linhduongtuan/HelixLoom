//! Work that takes seconds (the GPU's first run, MAFFT, IQ-TREE) runs on its own thread, so the window keeps
//! drawing while it goes. The window asks `poll` once a frame.
//!
//! A browser page has no threads to start (yet), so there the work runs at once, on the window's own thread, and is done
//! by the first `poll`: the window waits for it, and what takes long makes the page stall. Work that only waits for
//! another machine (a job on the workstation) is written as an `async` block and started with `spawn_async`: in a page it
//! goes on in the background, between frames, and the page does not stall.

use eframe::egui;
use std::future::Future;
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

    /// Start `work`, which waits for another machine, on a thread of its own; the window is asked to redraw when it is done.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn spawn_async(
        ctx: &egui::Context,
        work: impl Future<Output = T> + Send + 'static,
    ) -> Job<T> {
        Job::spawn(ctx, move || crate::remote::block::block_on(work))
    }

    /// Start `work`, which waits for another machine, in the background of the page: it goes on between frames, and the
    /// window is asked to redraw when it is done. (A page's futures need not be `Send`: there is one thread.)
    #[cfg(target_arch = "wasm32")]
    pub fn spawn_async(ctx: &egui::Context, work: impl Future<Output = T> + 'static) -> Job<T> {
        let (send, results) = mpsc::channel();
        let ctx = ctx.clone();
        wasm_bindgen_futures::spawn_local(async move {
            // the receiver is held by the job; if the job was dropped, nobody wants the answer
            let _ = send.send(work.await);
            ctx.request_repaint();
        });
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

#[cfg(all(test, not(target_arch = "wasm32")))]
mod tests {
    use super::*;
    use std::time::{Duration, Instant};

    /// Poll until the job is over, as a window does once a frame.
    fn finish<T: Send + 'static>(job: &Job<T>) -> Poll<T> {
        let until = Instant::now() + Duration::from_secs(10);
        loop {
            match job.poll() {
                Poll::Running if Instant::now() < until => {
                    std::thread::sleep(Duration::from_millis(5))
                }
                other => return other,
            }
        }
    }

    #[test]
    fn work_that_waits_for_another_machine_is_done_on_its_own_thread_and_polled() {
        let ctx = egui::Context::default();
        let main = std::thread::current().id();
        let job = Job::spawn_async(&ctx, async move {
            // a wait that is not ready at once, as a request is not
            std::thread::sleep(Duration::from_millis(30));
            (6 * 7, std::thread::current().id() != main)
        });
        // the window is not held up: the answer is not there yet
        assert!(matches!(job.poll(), Poll::Running));
        match finish(&job) {
            Poll::Done(answer) => assert_eq!(answer, (42, true)),
            _ => panic!("the job did not finish"),
        }
    }

    #[test]
    fn async_work_that_panics_is_reported_as_died_not_waited_for_ever() {
        let ctx = egui::Context::default();
        let job: Job<u8> = Job::spawn_async(&ctx, async { panic!("the work failed") });
        assert!(matches!(finish(&job), Poll::Died));
    }
}
