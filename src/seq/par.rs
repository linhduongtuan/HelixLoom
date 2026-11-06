//! Run independent jobs on every CPU core. The heavy loops of the platform (pairwise alignments,
//! likelihood searches) are plain functions over byte slices, so this is the one place a GPU
//! backend would have to replace later.

/// Apply `f` to every item on all cores; results keep the order of `items`.
#[cfg(not(target_arch = "wasm32"))]
pub fn par_map<T: Sync, R: Send>(items: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let workers = std::thread::available_parallelism().map_or(1, |w| w.get());
    let per_worker = items.len().div_ceil(workers).max(1);
    let f = &f;
    std::thread::scope(|s| {
        let jobs: Vec<_> = items
            .chunks(per_worker)
            .map(|part| s.spawn(move || part.iter().map(f).collect::<Vec<R>>()))
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

#[cfg(test)]
mod tests {
    use super::*;

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
