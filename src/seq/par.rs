//! Run independent jobs on every CPU core. The heavy loops of the platform (pairwise alignments,
//! likelihood searches) are plain functions over byte slices, so this is the one place a GPU
//! backend would have to replace later.

thread_local! {
    static WORKER: std::cell::Cell<Option<usize>> = const { std::cell::Cell::new(None) };
}

/// Which thread of `par_map` this is (the number of the part of the items it was given, from 0), or `None` on any other thread.
/// Work that goes to GPUs uses it to give each of the threads that run side by side a card of its own.
#[cfg_attr(
    not(all(feature = "gpu", not(target_arch = "wasm32"))),
    allow(dead_code)
)]
pub fn worker() -> Option<usize> {
    WORKER.with(std::cell::Cell::get)
}

/// Apply `f` to every item on all cores; results keep the order of `items`.
#[cfg(not(target_arch = "wasm32"))]
pub fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let workers = std::thread::available_parallelism().map_or(1, |w| w.get());
    let per_worker = items.len().div_ceil(workers).max(1);
    let f = &f;
    std::thread::scope(|s| {
        let jobs: Vec<_> = items
            .chunks(per_worker)
            .enumerate()
            .map(|(k, part)| {
                s.spawn(move || {
                    WORKER.with(|w| w.set(Some(k)));
                    part.iter().map(f).collect::<Vec<R>>()
                })
            })
            .collect();
        jobs.into_iter()
            .flat_map(|job| job.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    })
}

/// Apply `f` to every item, one after the other: a browser page has no threads to share the work among (yet).
#[cfg(target_arch = "wasm32")]
pub fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    items.iter().map(f).collect()
}

/// Apply `f` to every item on all cores, each thread taking the next item that nobody has taken, so that the cores that are
/// slower (the efficiency cores of a laptop, or a core that something else is using) do fewer items instead of making the others
/// wait for them: for items of about the same size that are several for each core, where `par_map` would cut them in equal parts.
/// Results keep the order of `items`; a thread's `worker()` is its own number, from 0.
#[cfg(not(target_arch = "wasm32"))]
pub fn par_queue<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let workers = std::thread::available_parallelism()
        .map_or(1, |w| w.get())
        .min(items.len());
    let next = std::sync::atomic::AtomicUsize::new(0);
    let (f, next) = (&f, &next);
    let mut done: Vec<(usize, R)> = std::thread::scope(|s| {
        let jobs: Vec<_> = (0..workers)
            .map(|k| {
                s.spawn(move || {
                    WORKER.with(|w| w.set(Some(k)));
                    let mut mine = Vec::new();
                    loop {
                        let at = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        let Some(item) = items.get(at) else { break };
                        mine.push((at, f(item)));
                    }
                    mine
                })
            })
            .collect();
        jobs.into_iter()
            .flat_map(|job| job.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    });
    done.sort_by_key(|&(at, _)| at);
    done.into_iter().map(|(_, r)| r).collect()
}

/// Apply `f` to every item, one after the other: a browser page has no threads to share the work among (yet).
#[cfg(target_arch = "wasm32")]
pub fn par_queue<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    items.iter().map(f).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_thread_knows_its_number_and_other_threads_do_not() {
        assert_eq!(worker(), None);
        let items: Vec<u32> = (0..4).collect();
        let numbers = par_map(&items, |_| worker());
        // every item was done by a worker, and the numbers are those of the parts (one each, or shared when there are fewer cores)
        assert!(numbers.iter().all(Option::is_some));
        assert_eq!(worker(), None);
    }

    #[test]
    fn the_queue_does_every_item_once_in_order_on_threads_that_know_their_number() {
        let items: Vec<u64> = (0..1000).collect();
        let done = par_queue(&items, |&x| (x * x, worker()));
        assert_eq!(done.len(), 1000);
        assert!(done
            .iter()
            .enumerate()
            .all(|(i, &(s, w))| s == (i * i) as u64 && w.is_some()));
        assert_eq!(worker(), None);
        assert!(par_queue(&[] as &[u32], |&x| x).is_empty());
        // a slow item does not hold up the others: with one item slow, the rest are still all done
        let items: Vec<u32> = (0..64).collect();
        let slow = par_queue(&items, |&x| {
            if x == 0 {
                std::thread::sleep(std::time::Duration::from_millis(30));
            }
            x + 1
        });
        assert_eq!(slow, (1..=64).collect::<Vec<u32>>());
    }

    #[test]
    fn results_keep_their_order_whatever_the_core_count() {
        let items: Vec<u32> = (0..1000).collect();
        let squares = par_map(&items, |&x| x * x);
        assert_eq!(squares.len(), 1000);
        assert!(squares
            .iter()
            .enumerate()
            .all(|(i, &s)| s == (i * i) as u32));
        assert!(par_map(&[] as &[u32], |&x| x).is_empty());
    }
}
