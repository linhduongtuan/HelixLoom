//! These need a Python with PyTorch (`HELIXLOOM_TORCH_PYTHON`); without one they say so and pass. The workers run on the CPU device.

use super::*;
use crate::accel::torch::tests::pools;
use crate::seq::likelihood::bayes::{bayes_tree_using, Settings};
use crate::seq::likelihood::check::{parameters, setup as make, simulate};
use crate::seq::likelihood::engine::Engine as Any;
use crate::seq::likelihood::{ml_tree_using, Kind};
use crate::seq::msa::Msa;
use crate::seq::par::par_map;
use crate::seq::tree::splits;

fn setup(msa: &Msa) -> (Sites, Tree) {
    make(msa).expect("sites and a tree")
}

fn near(got: f64, want: f64, relative: f64) -> bool {
    (got - want).abs() <= relative * want.abs().max(1.0)
}

#[test]
fn a_tree_has_the_same_likelihood_in_pytorch_as_on_the_processor_for_every_model() {
    let msa = simulate(14, 500, 0.12, 3);
    let (sites, tree) = setup(&msa);
    for pool in pools() {
        let engine = Engine::new(&pool, &sites).expect("an engine");
        for kind in Kind::ALL {
            let (q, freqs) = parameters(kind, &sites);
            let want = tree.ln_likelihood(&sites, &q, &freqs);
            let got = engine
                .ln_likelihood(&tree, &q, &freqs)
                .expect("a likelihood");
            assert!(
                near(got, want, 2e-6),
                "{} on {} workers: {got} against {want}",
                kind.label(),
                pool.len()
            );
        }
        assert!(!engine.broken());
    }
}

#[test]
fn a_tree_of_hundreds_of_tips_does_not_underflow_in_pytorch_either() {
    let msa = simulate(300, 120, 0.2, 9);
    let (sites, tree) = setup(&msa);
    let (q, freqs) = parameters(Kind::Hky85, &sites);
    let want = tree.ln_likelihood(&sites, &q, &freqs);
    assert!(want < -10_000.0, "{want}");
    for pool in pools() {
        let engine = Engine::new(&pool, &sites).expect("an engine");
        let got = engine
            .ln_likelihood(&tree, &q, &freqs)
            .expect("a likelihood");
        assert!(
            near(got, want, 5e-6),
            "{} workers: {got} against {want}",
            pool.len()
        );
    }
}

#[test]
fn the_threads_of_par_map_each_use_a_worker_and_the_main_thread_shares_the_columns_in_pytorch() {
    let msa = simulate(10, 400, 0.15, 17);
    let (sites, tree) = setup(&msa);
    let (q, freqs) = parameters(Kind::K80, &sites);
    let want = tree.ln_likelihood(&sites, &q, &freqs);
    for pool in pools() {
        let engine = Engine::new(&pool, &sites).expect("an engine");
        let shared = engine
            .ln_likelihood(&tree, &q, &freqs)
            .expect("a likelihood");
        assert!(near(shared, want, 2e-6), "{} workers", pool.len());
        let items: Vec<usize> = (0..12).collect();
        let each = par_map(&items, |_| engine.ln_likelihood(&tree, &q, &freqs));
        for value in each {
            let value = value.expect("a likelihood");
            assert!(near(value, want, 2e-6), "{} workers", pool.len());
        }
    }
}

#[test]
fn a_dead_worker_hands_the_likelihood_to_the_processor_and_the_report_says_why() {
    let Some(pool) = crate::accel::torch::tests::cpu_pool(&[1.0]) else {
        return;
    };
    let msa = simulate(8, 200, 0.1, 5);
    let (mut sites, tree) = setup(&msa);
    let (q, freqs) = parameters(Kind::F81, &sites);
    let want = tree.ln_likelihood(&sites, &q, &freqs);
    let engine = Engine::new(&pool, &sites).expect("an engine");
    // a tree the engine cannot take (two tips) is a failure as far as the caller can tell
    let two_tips = Tree {
        adj: vec![vec![(1, 0.1)], vec![(0, 0.1)]],
        tips: 2,
    };
    assert_eq!(engine.ln_likelihood(&two_tips, &q, &freqs), None);
    assert!(engine.broken());
    assert_eq!(engine.ln_likelihood(&tree, &q, &freqs), None);
    sites.accel = Some(Any::Torch(engine));
    assert_eq!(tree.ln_likelihood(&sites, &q, &freqs), want);
    let report = sites
        .accelerator(Some("on the cards".into()))
        .expect("a report");
    assert!(report.contains("failed"), "{report}");
}

#[test]
fn what_cannot_be_put_on_the_workers_is_refused_with_a_reason() {
    let Some(pool) = crate::accel::torch::tests::cpu_pool(&[1.0]) else {
        return;
    };
    let msa = simulate(8, 50, 0.1, 1);
    let (mut sites, _) = setup(&msa);
    sites.tips = 2;
    assert!(Engine::new(&pool, &sites).is_err());
    sites.tips = 8;
    sites.patterns.clear();
    assert!(Engine::new(&pool, &sites).is_err());
}

#[test]
fn the_search_in_pytorch_finds_the_tree_and_the_model_that_the_processor_finds() {
    let msa = simulate(9, 700, 0.1, 29);
    let on_cpu = ml_tree_using(&msa, |_| None).expect("a tree");
    for pool in pools() {
        let on_workers = ml_tree_using(&msa, |sites| {
            sites.accel = Engine::new(&pool, sites).ok().map(Any::Torch);
            sites.accel.as_ref().map(|e| e.describe())
        })
        .expect("a tree");
        assert_eq!(
            splits(&on_workers.tree, &msa.names),
            splits(&on_cpu.tree, &msa.names),
            "{} workers",
            pool.len()
        );
        assert_eq!(on_workers.best.kind, on_cpu.best.kind);
        assert!((on_workers.best.ln_l - on_cpu.best.ln_l).abs() < 0.5);
        assert!(on_workers
            .accelerator
            .as_deref()
            .unwrap_or("")
            .contains("PyTorch"));
    }
}

#[test]
fn the_sampler_in_pytorch_runs_its_chains_side_by_side_and_agrees_with_the_processor() {
    let msa = simulate(8, 400, 0.1, 41);
    let settings = Settings {
        generations: 400,
        sample_every: 5,
        ..Settings::default()
    };
    let on_cpu = bayes_tree_using(&msa, &settings, |_| None).expect("a sample");
    for pool in pools() {
        let on_workers = bayes_tree_using(&msa, &settings, |sites| {
            sites.accel = Engine::new(&pool, sites).ok().map(Any::Torch);
            sites.accel.as_ref().map(|e| e.describe())
        })
        .expect("a sample");
        assert_eq!(on_workers.samples, on_cpu.samples, "{} workers", pool.len());
        assert!(
            (on_workers.mean_ln_l - on_cpu.mean_ln_l).abs() < 3.0,
            "{} workers: {} against {}",
            pool.len(),
            on_workers.mean_ln_l,
            on_cpu.mean_ln_l
        );
        assert!(on_workers.accelerator.is_some());
    }
}

#[test]
fn the_likelihood_of_one_branch_for_every_length_on_pytorch_workers_is_that_of_the_whole_tree() {
    let msa = simulate(14, 500, 0.12, 3);
    let (sites, tree) = setup(&msa);
    for pool in pools() {
        let engine = Engine::new(&pool, &sites).expect("an engine");
        for kind in [Kind::Jc69, Kind::Hky85] {
            let (q, freqs) = parameters(kind, &sites);
            // every branch: to a tip and between inner nodes
            for (a, b) in tree.edges() {
                let scan = engine.edge_scan(&tree, a, b, &q, &freqs).expect("a scan");
                for length in [1e-6, 0.01, 0.1, 0.7, 5.0] {
                    let mut set = tree.clone();
                    set.set_length(a, b, length);
                    let want = set.ln_likelihood(&sites, &q, &freqs);
                    let got = scan.ln_likelihood(&q, length).expect("a likelihood");
                    assert!(
                        near(got, want, 5e-6),
                        "{} on {} workers, branch {a}-{b} at {length}: {got} against {want}",
                        kind.label(),
                        pool.len()
                    );
                }
            }
        }
        assert!(!engine.broken());
    }
}

#[test]
fn scans_made_by_threads_that_run_side_by_side_each_go_to_a_worker_and_agree_with_the_processor() {
    let msa = simulate(12, 400, 0.1, 17);
    let (sites, tree) = setup(&msa);
    let (q, freqs) = parameters(Kind::Hky85, &sites);
    let edges = tree.edges();
    for pool in pools() {
        let engine = Engine::new(&pool, &sites).expect("an engine");
        let got = par_map(&edges, |&(a, b)| {
            let scan = engine.edge_scan(&tree, a, b, &q, &freqs).expect("a scan");
            scan.ln_likelihood(&q, 0.2).expect("a likelihood")
        });
        for (&(a, b), got) in edges.iter().zip(got) {
            let mut set = tree.clone();
            set.set_length(a, b, 0.2);
            let want = set.ln_likelihood(&sites, &q, &freqs);
            assert!(
                near(got, want, 5e-6),
                "{} workers, branch {a}-{b}: {got} against {want}",
                pool.len()
            );
        }
        assert!(!engine.broken());
    }
}

#[test]
fn setting_the_branches_on_the_workers_reaches_the_likelihood_that_the_processor_reaches() {
    let msa = simulate(10, 500, 0.12, 5);
    let (sites, tree) = setup(&msa);
    let (q, freqs) = parameters(Kind::Hky85, &sites);
    let edges = tree.edges();
    let mut on_cpu = tree.clone();
    crate::seq::likelihood::optimize_edges(&mut on_cpu, &edges, 2, &sites, &q, &freqs);
    let want = on_cpu.ln_likelihood(&sites, &q, &freqs);
    for pool in pools() {
        let (mut carded, _) = setup(&msa);
        carded.accel = Engine::new(&pool, &carded).ok().map(Any::Torch);
        let mut on_workers = tree.clone();
        crate::seq::likelihood::optimize_edges(&mut on_workers, &edges, 2, &carded, &q, &freqs);
        let got = on_workers.ln_likelihood(&sites, &q, &freqs);
        // the workers sum 32-bit logarithms, so the likelihood is the processor's to about seven digits, and the branches a little apart
        assert!(
            near(got, want, 5e-5),
            "{} workers: {got} against {want}",
            pool.len()
        );
        assert!(!carded.accel.as_ref().is_some_and(|e| e.broken()));
    }
}

#[test]
#[ignore = "timing: HELIXLOOM_TORCH_PYTHON=... cargo test timing_of_setting_branches_on_pytorch -- --ignored --nocapture"]
fn timing_of_setting_branches_on_pytorch_workers() {
    use crate::seq::likelihood::{best_length_by_whole_tree, best_length_on_cards};
    for (tips, columns) in [(10usize, 500usize), (20, 1000), (40, 1500)] {
        let msa = simulate(tips, columns, 0.12, 9);
        let (mut sites, tree) = setup(&msa);
        let (q, freqs) = parameters(Kind::Hky85, &sites);
        let Some(pool) = pools().into_iter().next() else {
            return;
        };
        sites.accel = Engine::new(&pool, &sites).ok().map(Any::Torch);
        let edges: Vec<(usize, usize)> = tree.edges().into_iter().take(12).collect();
        let time = |on_cards: bool| {
            let mut set = tree.clone();
            let start = std::time::Instant::now();
            let best: Vec<f64> = edges
                .iter()
                .map(|&edge| {
                    if on_cards {
                        best_length_on_cards(&set, edge, &sites, &q, &freqs).expect("a length")
                    } else {
                        best_length_by_whole_tree(&mut set, edge, &sites, &q, &freqs)
                    }
                })
                .collect();
            (start.elapsed().as_secs_f64(), best)
        };
        let (whole, a) = time(false);
        let (scan, b) = time(true);
        let worst = a
            .iter()
            .zip(&b)
            .map(|(x, y)| (x - y).abs())
            .fold(0.0, f64::max);
        println!(
            "{tips} tips x {columns} columns, 12 branches: {whole:.2} s with the whole tree for every length, {scan:.2} s with the scan ({:.1} times); lengths agree to {worst:.1e} in the logarithm",
            whole / scan
        );
    }
}

#[test]
fn a_card_that_dies_in_the_middle_costs_a_share_and_the_other_card_goes_on_with_the_likelihood() {
    use crate::accel::torch::{configured_python, tests::kill_worker, TorchPool};
    let Some(python) = configured_python() else {
        return;
    };
    // a pool of its own: one worker of it is going to be killed
    let pool = std::sync::Arc::new(
        TorchPool::start_weighted(
            &python,
            Some(vec!["cpu".into(), "cpu".into()]),
            Some(vec![1.0, 1.0]),
        )
        .expect("workers")
        .expect("a pool"),
    );
    let msa = simulate(10, 400, 0.1, 7);
    let (mut sites, tree) = setup(&msa);
    let (q, freqs) = parameters(Kind::Hky85, &sites);
    let want = tree.ln_likelihood(&sites, &q, &freqs);
    let engine = Engine::new(&pool, &sites).expect("an engine");
    let before = engine
        .ln_likelihood(&tree, &q, &freqs)
        .expect("a likelihood");
    assert!(near(before, want, 2e-6));
    kill_worker(&pool, 1);
    // each thread of `par_map` has a worker of its own: the threads that were given the dead one have it done by the other, and the cards
    // are not given up (a likelihood of a few patterns is not cut among the workers: it goes whole to the fastest)
    let threads: Vec<usize> = (0..4).collect();
    for _ in 0..2 {
        let after = par_map(&threads, |_| engine.ln_likelihood(&tree, &q, &freqs));
        for value in after {
            let value = value.expect("the other worker did it");
            assert!(near(value, want, 2e-6), "{value} against {want}");
        }
    }
    assert!(!engine.broken());
    assert!(engine.retried().is_some());
    // a branch's scan may not be kept by a dead worker (the whole tree's likelihood, which retries, does that length instead), but it
    // never gives a wrong number
    let (a, b) = tree.edges()[3];
    let scans = par_map(&threads, |_| {
        engine
            .edge_scan(&tree, a, b, &q, &freqs)
            .and_then(|scan| scan.ln_likelihood(&q, 0.1))
    });
    let mut whole = tree.clone();
    whole.set_length(a, b, 0.1);
    let want_scan = whole.ln_likelihood(&sites, &q, &freqs);
    assert!(scans.iter().flatten().all(|&v| near(v, want_scan, 5e-6)));
    assert!(
        scans.iter().flatten().count() >= 1,
        "the live worker's threads keep theirs"
    );
    sites.accel = Some(Any::Torch(engine));
    let report = sites
        .accelerator(Some("on the cards".into()))
        .expect("a report");
    assert!(report.contains("another GPU did its share"), "{report}");
}
