//! These need a Python with PyTorch: `HELIXLOOM_TORCH_PYTHON=/path/to/python cargo test`. Without one they say so and pass, as the
//! build the continuous tests run (no PyTorch there) must; the workers run on the CPU device, which runs the very code a card does.

use super::*;
use crate::accel::distance::tests::{alignment, msa};
use crate::seq::tree::{cpu_distances, DistanceModel, Rng};

/// A pool of PyTorch workers on the CPU, one for each weight, or `None` (and a note) where no Python is named. Pools are kept for
/// the tests that follow: starting one means importing PyTorch, seconds for each worker.
pub(crate) fn cpu_pool(weights: &[f64]) -> Option<Arc<TorchPool>> {
    // the weights as bits, with the pool they made
    type Kept = Vec<(Vec<u64>, Arc<TorchPool>)>;
    static POOLS: OnceLock<Mutex<Kept>> = OnceLock::new();
    let Some(python) = configured_python() else {
        eprintln!("skipped: {PYTHON_VARIABLE} names no Python with PyTorch");
        return None;
    };
    let key: Vec<u64> = weights.iter().map(|w| w.to_bits()).collect();
    let mut kept = POOLS
        .get_or_init(Mutex::default)
        .lock()
        .unwrap_or_else(PoisonError::into_inner);
    if let Some((_, pool)) = kept.iter().find(|(k, _)| *k == key) {
        return Some(Arc::clone(pool));
    }
    let devices = vec!["cpu".to_string(); weights.len()];
    let pool = TorchPool::start_weighted(&python, Some(devices), Some(weights.to_vec()))
        .expect("workers that start")
        .expect("a pool");
    let pool = Arc::new(pool);
    kept.push((key, Arc::clone(&pool)));
    Some(pool)
}

/// The pools to try: one worker, two with unequal shares, and three (so that the work does not divide).
pub(crate) fn pools() -> Vec<Arc<TorchPool>> {
    [vec![1.0], vec![2.0, 1.0], vec![5.0, 3.0, 1.0]]
        .iter()
        .filter_map(|w| cpu_pool(w))
        .collect()
}

/// Kill worker `k` of `pool`, as a card that dies in the middle of a job would be (use a pool of the test's own, not a shared one).
pub(crate) fn kill_worker(pool: &TorchPool, k: usize) {
    pool.workers()[k]
        .child
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .kill()
        .expect("killed");
}

#[test]
fn the_workers_say_what_they_are_and_the_shares_cover_the_work() {
    let Some(pool) = cpu_pool(&[2.0, 1.0]) else {
        return;
    };
    assert_eq!(pool.len(), 2);
    assert!(pool.versions.starts_with("PyTorch "), "{}", pool.versions);
    assert!(pool.workers()[0].name.starts_with("cpu "));
    assert_eq!(pool.distance_shares(300), vec![0..200, 200..300]);
    assert_eq!(
        pool.likelihood_shares(3)
            .iter()
            .map(|r| r.len())
            .sum::<usize>(),
        3
    );
    assert!(pool.likelihood_description().contains("67 %"));
}

#[test]
fn the_pytorch_distances_are_those_of_the_processor_to_the_last_digit() {
    let rows = alignment(13, 777, 0.15, 7);
    for pool in pools() {
        let engine = distance::Engine::new(&pool, &rows).expect("an engine");
        let counts = engine.counts(None).expect("counts");
        for model in DistanceModel::ALL {
            let want = cpu_distances(&msa(&rows), model).expect("distances");
            assert_eq!(
                counts.distances(model).expect("distances"),
                want,
                "{} on {} workers",
                model.label(),
                pool.len()
            );
        }
    }
}

#[test]
fn a_bootstrap_resample_as_weights_gives_the_distances_of_the_resampled_alignment_in_pytorch() {
    let rows = alignment(8, 400, 0.12, 21);
    let width = rows[0].len();
    let mut rng = Rng::new(99);
    for pool in pools() {
        let engine = distance::Engine::new(&pool, &rows).expect("an engine");
        for _ in 0..3 {
            let picks: Vec<usize> = (0..width).map(|_| rng.below(width)).collect();
            let mut weights = vec![0u32; width];
            picks.iter().for_each(|&c| weights[c] += 1);
            let resampled: Vec<Vec<u8>> = rows
                .iter()
                .map(|r| picks.iter().map(|&c| r[c]).collect())
                .collect();
            let counts = engine.counts(Some(&weights)).expect("counts");
            for model in DistanceModel::ALL {
                assert_eq!(
                    counts.distances(model),
                    Ok(cpu_distances(&msa(&resampled), model).expect("distances")),
                    "{} on {} workers",
                    model.label(),
                    pool.len()
                );
            }
        }
    }
}

#[test]
fn what_cannot_be_counted_is_refused_before_it_is_sent_and_a_pair_that_cannot_be_compared_is_named()
{
    let Some(pool) = cpu_pool(&[1.0, 1.0]) else {
        return;
    };
    assert!(distance::Engine::new(&pool, &[b"ACGT".to_vec()]).is_err());
    assert!(distance::Engine::new(&pool, &[b"ACGT".to_vec(), b"ACG".to_vec()]).is_err());
    let rows: Vec<Vec<u8>> = ["ACGT----", "----ACGT", "ACGTACGT"]
        .iter()
        .map(|r| r.as_bytes().to_vec())
        .collect();
    let engine = distance::Engine::new(&pool, &rows).expect("an engine");
    assert!(engine.counts(Some(&[1, 2, 3])).is_err());
    assert!(engine.counts(Some(&[2000; 8])).is_err());
    assert_eq!(
        engine
            .counts(None)
            .expect("counts")
            .distances(DistanceModel::P),
        Err((0, 1))
    );
}

#[test]
fn a_worker_that_dies_gives_an_error_and_is_not_asked_again_and_the_others_go_on() {
    let Some(python) = configured_python() else {
        eprintln!("skipped: {PYTHON_VARIABLE} names no Python with PyTorch");
        return;
    };
    let pool = Arc::new(
        TorchPool::start_weighted(
            &python,
            Some(vec!["cpu".into(), "cpu".into()]),
            Some(vec![1.0, 1.0]),
        )
        .expect("workers")
        .expect("a pool"),
    );
    let rows = alignment(6, 200, 0.1, 3);
    let engine = distance::Engine::new(&pool, &rows).expect("an engine");
    assert!(engine.counts(None).is_ok());
    pool.workers()[1]
        .child
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .kill()
        .expect("killed");
    let first = engine.counts(None).expect_err("a worker is dead");
    assert!(first.contains("stopped"), "{first}");
    let second = engine.counts(None).expect_err("still dead");
    assert!(second.contains("has stopped"), "{second}");
    // the other worker still answers for what it holds
    assert!(pool.workers()[0].call(&Message::new("speed")).is_ok());
}

#[test]
fn an_error_inside_a_worker_comes_back_as_text_and_the_worker_stays_up() {
    let Some(pool) = cpu_pool(&[1.0]) else {
        return;
    };
    let worker = &pool.workers()[0];
    let error = worker
        .call(&Message::new("counts").with_number("id", 999_999))
        .expect_err("no such set");
    assert!(error.contains("KeyError"), "{error}");
    assert!(worker.call(&Message::new("speed")).is_ok());
}

/// `worker.py` run inside `before`, which is Python code that runs first (it may change what PyTorch does): the worker's script is
/// handed over in a variable, because `from __future__` must stay the first thing in it.
fn wrapped(before: &str) -> String {
    format!("import os\n{before}\nexec(compile(os.environ['HX_SCRIPT'], 'worker.py', 'exec'))")
}

fn spawn_wrapped(
    before: &str,
    start_timeout: Duration,
    timeout: Duration,
) -> Option<Result<(Worker, Hello), String>> {
    let python = configured_python()?;
    Some(Worker::spawn_script(
        &python,
        &wrapped(before),
        &[("HX_SCRIPT", SCRIPT)],
        "cpu",
        start_timeout,
        timeout,
    ))
}

#[test]
fn a_library_that_writes_to_standard_output_does_not_break_the_messages() {
    // anything that prints, in Python or in C, would land in the middle of a reply if the protocol still used descriptor 1
    let noisy = "import torch\n_ones = torch.ones\ndef ones(*a, **k):\n    print('noise from print')\n    os.write(1, b'noise from a C library\\n')\n    return _ones(*a, **k)\ntorch.ones = ones";
    let Some(started) = spawn_wrapped(noisy, Duration::from_secs(120), Duration::from_secs(120))
    else {
        eprintln!("skipped: {PYTHON_VARIABLE} names no Python with PyTorch");
        return;
    };
    let (worker, _) = started.expect("a worker");
    for _ in 0..3 {
        let reply = worker.call(&Message::new("speed")).expect("an answer");
        assert!(reply.f64s("memory").is_some_and(|m| m[0] > 0.0));
    }
}

#[test]
fn a_worker_that_hangs_is_killed_after_the_time_it_was_given_and_is_not_asked_again() {
    let hang = "import torch, time\n_ones = torch.ones\ndef ones(*a, **k):\n    time.sleep(3600)\ntorch.ones = ones";
    let Some(started) = spawn_wrapped(hang, Duration::from_secs(120), Duration::from_secs(1))
    else {
        eprintln!("skipped: {PYTHON_VARIABLE} names no Python with PyTorch");
        return;
    };
    let (worker, _) = started.expect("a worker");
    let begun = std::time::Instant::now();
    let error = worker.call(&Message::new("speed")).expect_err("it hangs");
    assert!(error.contains("did not answer within 1 seconds"), "{error}");
    assert!(begun.elapsed() < Duration::from_secs(30));
    assert!(worker
        .call(&Message::new("speed"))
        .expect_err("stopped")
        .contains("has stopped"));
    // and the process is gone, not left sleeping
    let status = worker
        .child
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .wait()
        .expect("ended");
    assert!(!status.success());
}

#[test]
fn a_worker_that_does_not_start_in_time_or_cannot_import_pytorch_is_an_error_with_the_reason() {
    let Some(python) = configured_python() else {
        eprintln!("skipped: {PYTHON_VARIABLE} names no Python with PyTorch");
        return;
    };
    let slow = Worker::spawn_script(
        &python,
        &wrapped("import time\ntime.sleep(3600)"),
        &[("HX_SCRIPT", SCRIPT)],
        "cpu",
        Duration::from_secs(1),
        Duration::from_secs(1),
    )
    .err()
    .expect("it never says hello");
    assert!(slow.contains("did not start within 1 seconds"), "{slow}");
    let missing = Worker::spawn("/nonexistent/python", "cpu")
        .err()
        .expect("no such program");
    assert!(missing.contains("Cannot start"), "{missing}");
    // a Python that does not have PyTorch: a Python of the system, when it has none
    let without = std::process::Command::new("python3")
        .args(["-I", "-c", "import torch"])
        .output();
    if without.is_ok_and(|o| !o.status.success()) {
        let why = Worker::spawn("python3", "cpu").err().expect("no PyTorch");
        assert!(why.contains("is PyTorch (and NumPy) installed"), "{why}");
    }
}

#[test]
fn an_error_in_a_request_is_answered_once_and_the_next_request_still_works() {
    let Some(pool) = cpu_pool(&[1.0]) else {
        return;
    };
    let worker = &pool.workers()[0];
    for _ in 0..2 {
        let error = worker
            .call(&Message::new("no_such_operation"))
            .expect_err("unknown");
        assert!(error.contains("KeyError"), "{error}");
        assert!(worker.call(&Message::new("speed")).is_ok());
    }
}

#[test]
fn measured_shares_are_positive_add_up_to_one_and_the_first_worker_is_named_for_its_device() {
    let Some(python) = configured_python() else {
        eprintln!("skipped: {PYTHON_VARIABLE} names no Python with PyTorch");
        return;
    };
    // no weights given: each worker's speed is measured
    let pool = TorchPool::start(&python, Some(vec!["cpu".into(), "cpu".into()]))
        .expect("workers")
        .expect("a pool");
    let cards = pool.cards();
    assert_eq!(cards.len(), 2);
    for (name, likelihood, distance) in &cards {
        assert!(name.starts_with("cpu "), "{name}");
        assert!(*likelihood > 0.0 && *distance > 0.0);
    }
    assert!((cards.iter().map(|c| c.1).sum::<f64>() - 1.0).abs() < 1e-9);
    assert!((cards.iter().map(|c| c.2).sum::<f64>() - 1.0).abs() < 1e-9);
}

#[test]
fn a_likelihood_is_cut_among_workers_only_as_far_as_each_has_enough_patterns() {
    // little: the whole to the fastest, the others nothing
    assert_eq!(split(&[2.0, 1.0], 100_000, 200_000), vec![0..100_000, 0..0]);
    assert_eq!(split(&[1.0, 2.0], 100_000, 200_000), vec![0..0, 0..100_000]);
    assert_eq!(split(&[1.0], 5, 200_000), vec![0..5]);
    // enough for two of three: the two fastest, cut by their shares, in the order of the workers
    assert_eq!(
        split(&[5.0, 3.0, 1.0], 450_000, 200_000),
        vec![0..281_250, 281_250..450_000, 0..0]
    );
    assert_eq!(
        split(&[1.0, 3.0, 5.0], 450_000, 200_000),
        vec![0..0, 0..168_750, 168_750..450_000]
    );
    // plenty: all of them, as `shares` cuts
    assert_eq!(
        split(&[2.0, 1.0], 900_000, 200_000),
        shares(&[2.0 / 3.0, 1.0 / 3.0], 900_000)
    );
    // nothing is lost or counted twice, whatever the sizes
    for total in [
        0, 1, 199_999, 200_000, 399_999, 400_000, 600_001, 10_000_000,
    ] {
        let cut = split(&[5.0, 3.0, 1.0], total, 200_000);
        assert_eq!(cut.iter().map(|r| r.len()).sum::<usize>(), total, "{total}");
    }
}

/// Distances between `n` random points in the plane: a matrix with no ties, which the joins of the two sides must agree on.
fn plane_distances(n: usize, seed: u64) -> Vec<Vec<f64>> {
    let mut rng = Rng::new(seed);
    let points: Vec<(f64, f64)> = (0..n)
        .map(|_| (rng.unit() * 10.0, rng.unit() * 10.0))
        .collect();
    points
        .iter()
        .map(|a| {
            points
                .iter()
                .map(|b| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt())
                .collect()
        })
        .collect()
}

/// The lengths of all the branches of a tree added up.
fn total_length(node: &crate::seq::tree::Node) -> f64 {
    node.length + node.children.iter().map(total_length).sum::<f64>()
}

#[test]
fn the_pytorch_joins_are_those_of_the_processor_for_matrices_without_ties() {
    use crate::seq::tree::{cpu_joins, join_tree, splits};
    let Some(pool) = cpu_pool(&[1.0]) else {
        return;
    };
    for (n, seed) in [(4usize, 1u64), (5, 2), (9, 3), (40, 4), (150, 5)] {
        let dist = plane_distances(n, seed);
        let (want, want_last) = cpu_joins(&dist);
        let (got, got_last) = nj::joins(&pool, &dist).expect("joins");
        assert_eq!(got.len(), want.len(), "{n} sequences");
        // the last join leaves a choice for four clusters, whatever the distances: the two pairs that split them the same way score the
        // same, and which is taken is a matter of rounding. Every join before it is the processor's.
        let strict = want.len().saturating_sub(1);
        for (k, (g, w)) in got.iter().zip(&want).enumerate().take(strict) {
            assert_eq!((g.i, g.j), (w.i, w.j), "{n} sequences, join {k}");
            assert!(
                (g.li - w.li).abs() <= 1e-9 * w.li.abs().max(1.0),
                "{n}: join {k}"
            );
            assert!(
                (g.lj - w.lj).abs() <= 1e-9 * w.lj.abs().max(1.0),
                "{n}: join {k}"
            );
        }
        // and the trees are the same tree: the same splits and the same total length
        let names: Vec<String> = (0..n).map(|k| format!("s{k}")).collect();
        let on_cards = join_tree(&names, &got, &got_last);
        let on_cpu = join_tree(&names, &want, &want_last);
        assert_eq!(on_cards.leaf_count(), n);
        assert_eq!(
            splits(&on_cards, &names),
            splits(&on_cpu, &names),
            "{n} sequences"
        );
        let (a, b) = (total_length(&on_cards), total_length(&on_cpu));
        assert!(
            (a - b).abs() <= 1e-9 * b.abs().max(1.0),
            "{n} sequences: {a} against {b}"
        );
    }
}

#[test]
fn identical_sequences_tie_every_join_and_the_pytorch_joins_are_the_processors_and_make_a_complete_tree(
) {
    use crate::seq::tree::join_tree;
    let Some(pool) = cpu_pool(&[1.0]) else {
        return;
    };
    // five sequences at distance zero from each other, three apart: which of the tied pairs is joined first is not fixed
    let n = 8;
    let dist: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                // zero within a group (and from a sequence to itself), three between the groups
                .map(|j| if (i < 5) == (j < 5) { 0.0 } else { 3.0 })
                .collect()
        })
        .collect();
    let (joins, last) = nj::joins(&pool, &dist).expect("joins");
    // the sums are exact for such distances, so the tied pairs are taken in the order the processor takes them: the earliest in the list
    let (want, want_last) = crate::seq::tree::cpu_joins(&dist);
    assert_eq!(joins, want);
    assert_eq!(last, want_last);
    let names: Vec<String> = (0..n).map(|k| format!("s{k}")).collect();
    let tree = join_tree(&names, &joins, &last);
    assert_eq!(tree.leaf_count(), n);
    // every sequence once: the names are all in the Newick text
    let text = tree.newick();
    assert!(
        names
            .iter()
            .all(|name| text.matches(name.as_str()).count() == 1),
        "{text}"
    );
    assert!(nj::joins(
        &pool,
        &dist[..3]
            .iter()
            .map(|r| r[..3].to_vec())
            .collect::<Vec<_>>()
    )
    .is_err());
    // a distance that is not a finite number, and a matrix that is not square, are refused before anything is sent
    let mut broken = dist.clone();
    broken[2][5] = f64::NAN;
    assert!(nj::joins(&pool, &broken)
        .expect_err("not finite")
        .contains("finite"));
    broken[2][5] = f64::INFINITY;
    assert!(nj::joins(&pool, &broken).is_err());
    let mut ragged = dist.clone();
    ragged[1].pop();
    assert!(nj::joins(&pool, &ragged).is_err());
}

#[test]
fn a_job_shown_no_gpu_does_not_start_pytorch_unless_a_device_is_named() {
    // the server shows a job that was given no GPU `CUDA_VISIBLE_DEVICES=` (nothing): there is nothing for PyTorch to find
    assert!(!can_have_devices(Some(""), false));
    assert!(!can_have_devices(Some("  "), false));
    assert!(can_have_devices(Some("0,1"), false));
    // not set at all (a Mac, or a person's own run): ask PyTorch
    assert!(can_have_devices(None, false));
    // a device named on purpose (`cpu`, to try the code) is used whatever the job was shown
    assert!(can_have_devices(Some(""), true));
}

#[test]
fn neighbour_joining_with_ties_in_a_matrix_big_enough_to_cut_down_agrees_with_the_processor_but_for_the_last_join(
) {
    use crate::seq::tree::cpu_joins;
    let Some(pool) = cpu_pool(&[1.0]) else {
        return;
    };
    // sixty sequences in groups of five (distance 0 inside a group) at whole-number distances between groups: ties at every join, and
    // the matrix is cut down to the places still alive more than once on the way
    let n = 60;
    let dist: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let (a, b) = (i / 5, j / 5);
                    if a == b {
                        0.0
                    } else {
                        1.0 + ((a * 7 + b * 7 + a * b) % 5) as f64
                    }
                })
                .collect()
        })
        .collect();
    let dist: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| dist[i.min(j)][i.max(j)]).collect())
        .collect();
    let (want, _) = cpu_joins(&dist);
    let (got, _) = nj::joins(&pool, &dist).expect("joins");
    assert_eq!(got.len(), want.len());
    let strict = want.len() - 1;
    assert_eq!(got[..strict], want[..strict]);
}

#[test]
fn a_distance_too_large_to_multiply_is_refused_before_it_is_sent() {
    let Some(pool) = cpu_pool(&[1.0]) else {
        return;
    };
    let mut dist = plane_distances(8, 1);
    dist[1][4] = 1e150;
    dist[4][1] = 1e150;
    assert!(nj::joins(&pool, &dist).is_err());
}
