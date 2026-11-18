//! Waiting for an `async fn` on a thread of its own. The client is written once as async functions, because a page can only
//! ask the workstation through `fetch`, which never blocks; a program has threads, so it just waits here for the end. No
//! executor is needed: a future that is not ready is waited for by parking the thread until it wakes it.

use std::future::Future;
use std::pin::pin;
use std::sync::Arc;
use std::task::{Context, Poll, Wake, Waker};
use std::thread::{self, Thread};

struct Unpark(Thread);

impl Wake for Unpark {
    fn wake(self: Arc<Self>) {
        self.0.unpark();
    }
}

/// Run `future` to its end on this thread.
pub fn block_on<F: Future>(future: F) -> F::Output {
    let mut future = pin!(future);
    let waker = Waker::from(Arc::new(Unpark(thread::current())));
    let mut context = Context::from_waker(&waker);
    loop {
        match future.as_mut().poll(&mut context) {
            Poll::Ready(output) => return output,
            Poll::Pending => thread::park(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::pin::Pin;

    /// Not ready the first `pending` times it is asked, and wakes whoever asks from another thread, as a timer does.
    struct After {
        pending: usize,
    }

    impl Future for After {
        type Output = &'static str;

        fn poll(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Self::Output> {
            if self.pending == 0 {
                return Poll::Ready("ready");
            }
            self.pending -= 1;
            let waker = cx.waker().clone();
            thread::spawn(move || {
                thread::sleep(std::time::Duration::from_millis(5));
                waker.wake();
            });
            Poll::Pending
        }
    }

    #[test]
    fn a_future_that_is_ready_gives_its_answer_at_once() {
        assert_eq!(block_on(async { 6 * 7 }), 42);
    }

    #[test]
    fn a_future_that_waits_is_waited_for_and_woken() {
        assert_eq!(block_on(After { pending: 3 }), "ready");
        let chained = async { After { pending: 1 }.await.len() + After { pending: 1 }.await.len() };
        assert_eq!(block_on(chained), 10);
    }
}
