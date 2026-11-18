use super::*;
use crate::accel::pool::Member;
use crate::seq::likelihood::bayes::{bayes_tree_using, Settings};
use crate::seq::likelihood::check::{parameters, simulate};
use crate::seq::likelihood::{ml_tree_using, Kind};
use crate::seq::msa::Msa;
use crate::seq::par::par_map;
use crate::seq::tree::splits;
use candle_core::Device;

fn cpu(weight: f64) -> Member {
    Member::new(Device::Cpu, "CPU".into(), weight)
}

/// The pools to try: one device, two unequal ones, three, and the real GPU of this machine (alone and beside a CPU) if it has one.
fn pools() -> Vec<Pool> {
    let mut pools = vec![
        Pool::new(vec![cpu(1.0)]),
        Pool::new(vec![cpu(2.0), cpu(1.0)]),
        Pool::new(vec![cpu(5.0), cpu(3.0), cpu(1.0)]),
    ];
    if let Some(real) = Pool::global() {
        pools.push(real.clone());
        let mut members: Vec<Member> = real.members().to_vec();
        members.push(cpu(1.0));
        pools.push(Pool::new(members));
    }
    pools
}

fn setup(msa: &Msa) -> (Sites, Tree) {
    crate::seq::likelihood::check::setup(msa).expect("sites and a tree")
}

fn near(got: f64, want: f64, relative: f64) -> bool {
    (got - want).abs() <= relative * want.abs().max(1.0)
}

#[test]
fn a_tree_has_the_same_likelihood_on_the_cards_as_on_the_processor_for_every_model() {
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
                "{} on {} cards: {got} against {want}",
                kind.label(),
                pool.len()
            );
        }
        assert!(!engine.broken());
    }
}

#[test]
fn a_tree_of_hundreds_of_tips_does_not_underflow_in_thirty_two_bits() {
    // the likelihood of a column is a product over three hundred tips: far below what 32 bits hold without the scaling
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
            "{} cards: {got} against {want}",
            pool.len()
        );
    }
}

#[test]
fn the_threads_of_par_map_each_use_a_card_and_the_main_thread_shares_the_columns_and_all_agree() {
    let msa = simulate(10, 400, 0.15, 17);
    let (sites, tree) = setup(&msa);
    let (q, freqs) = parameters(Kind::K80, &sites);
    let want = tree.ln_likelihood(&sites, &q, &freqs);
    for pool in pools() {
        let engine = Engine::new(&pool, &sites).expect("an engine");
        // the main thread: the columns are cut among the cards
        let shared = engine
            .ln_likelihood(&tree, &q, &freqs)
            .expect("a likelihood");
        assert!(
            near(shared, want, 2e-6),
            "{} cards: {shared} against {want}",
            pool.len()
        );
        // many threads at once, each on the card of its number: all get the same, and none waits for ever or gets a wrong answer
        let items: Vec<usize> = (0..12).collect();
        let each = par_map(&items, |_| engine.ln_likelihood(&tree, &q, &freqs));
        for value in each {
            let value = value.expect("a likelihood");
            assert!(
                near(value, want, 2e-6),
                "{} cards: {value} against {want}",
                pool.len()
            );
        }
    }
}

#[test]
fn a_card_that_fails_hands_the_work_to_the_processor_and_the_answer_is_still_right() {
    let msa = simulate(8, 200, 0.1, 5);
    let (mut sites, tree) = setup(&msa);
    let (q, freqs) = parameters(Kind::F81, &sites);
    let want = tree.ln_likelihood(&sites, &q, &freqs);
    let engine = Engine::new(&pools()[0], &sites).expect("an engine");
    // a tree the engine cannot take (two tips) is a failure of the card as far as the caller can tell
    let two_tips = Tree {
        adj: vec![vec![(1, 0.1)], vec![(0, 0.1)]],
        tips: 2,
    };
    assert_eq!(engine.ln_likelihood(&two_tips, &q, &freqs), None);
    assert!(engine.broken());
    // from then on it says nothing, and `Tree::ln_likelihood` goes on with the processor
    assert_eq!(engine.ln_likelihood(&tree, &q, &freqs), None);
    sites.accel = Some(crate::seq::likelihood::engine::Engine::Candle(engine));
    assert_eq!(tree.ln_likelihood(&sites, &q, &freqs), want);
    let report = sites
        .accelerator(Some("on the cards".into()))
        .expect("a report");
    assert!(report.contains("failed"), "{report}");
    // and the report says why, so that a person is not left guessing which card or what
    let why = sites
        .accel
        .as_ref()
        .and_then(|e| e.failure())
        .expect("a reason");
    assert!(report.contains(&why), "{report} / {why}");
}

#[test]
fn on_cards_a_change_of_topology_must_rise_by_more_than_the_noise_of_thirty_two_bits() {
    let msa = simulate(8, 200, 0.1, 5);
    let (mut sites, _) = setup(&msa);
    // on the processor, anything above rounding counts, whatever the size of the likelihood
    assert_eq!(sites.smallest_gain(-1.0e6), 1e-3);
    sites.accel = Some(crate::seq::likelihood::engine::Engine::Candle(
        Engine::new(&pools()[0], &sites).expect("an engine"),
    ));
    // on cards, a share of the likelihood, but never less than the processor's
    assert_eq!(sites.smallest_gain(-10.0), 1e-3);
    assert!((sites.smallest_gain(-1.0e6) - 1.0).abs() < 1e-9);
    assert!((sites.smallest_gain(1.0e6) - 1.0).abs() < 1e-9);
    // a card that has failed is not used any more, so the processor's rule applies again
    let two_tips = Tree {
        adj: vec![vec![(1, 0.1)], vec![(0, 0.1)]],
        tips: 2,
    };
    let (q, freqs) = parameters(Kind::F81, &sites);
    assert_eq!(
        sites
            .accel
            .as_ref()
            .and_then(|e| e.ln_likelihood(&two_tips, &q, &freqs)),
        None
    );
    assert_eq!(sites.smallest_gain(-1.0e6), 1e-3);
}

#[test]
fn what_cannot_be_put_on_the_cards_is_refused_with_a_reason() {
    let pool = &pools()[0];
    let msa = simulate(8, 50, 0.1, 1);
    let (mut sites, _) = setup(&msa);
    sites.tips = 2;
    assert!(Engine::new(pool, &sites).is_err());
    sites.tips = 8;
    sites.patterns.clear();
    assert!(Engine::new(pool, &sites).is_err());
}

#[test]
fn the_search_on_the_cards_finds_the_tree_and_the_model_that_the_processor_finds() {
    let msa = simulate(9, 700, 0.1, 29);
    let on_cpu = ml_tree_using(&msa, |_| None).expect("a tree");
    for pool in pools() {
        let on_cards = ml_tree_using(&msa, |sites| {
            sites.accel = Engine::new(&pool, sites)
                .ok()
                .map(crate::seq::likelihood::engine::Engine::Candle);
            sites.accel.as_ref().map(|e| e.describe())
        })
        .expect("a tree");
        assert_eq!(
            splits(&on_cards.tree, &msa.names),
            splits(&on_cpu.tree, &msa.names),
            "{} cards",
            pool.len()
        );
        assert_eq!(on_cards.best.kind, on_cpu.best.kind, "{} cards", pool.len());
        assert!(
            (on_cards.best.ln_l - on_cpu.best.ln_l).abs() < 0.5,
            "{} cards: {} against {}",
            pool.len(),
            on_cards.best.ln_l,
            on_cpu.best.ln_l
        );
        assert!(on_cards.accelerator.is_some());
    }
    assert!(on_cpu.accelerator.is_none());
}

#[test]
fn the_sampler_on_the_cards_runs_its_chains_side_by_side_and_agrees_with_the_processor() {
    let msa = simulate(8, 400, 0.1, 41);
    let settings = Settings {
        generations: 400,
        sample_every: 5,
        ..Settings::default()
    };
    let on_cpu = bayes_tree_using(&msa, &settings, |_| None).expect("a sample");
    for pool in pools() {
        let on_cards = bayes_tree_using(&msa, &settings, |sites| {
            sites.accel = Engine::new(&pool, sites)
                .ok()
                .map(crate::seq::likelihood::engine::Engine::Candle);
            sites.accel.as_ref().map(|e| e.describe())
        })
        .expect("a sample");
        assert_eq!(on_cards.samples, on_cpu.samples, "{} cards", pool.len());
        // the chains may part where two likelihoods that differ in the seventh digit decide a move differently: the posterior may not
        assert!(
            (on_cards.mean_ln_l - on_cpu.mean_ln_l).abs() < 3.0,
            "{} cards: {} against {}",
            pool.len(),
            on_cards.mean_ln_l,
            on_cpu.mean_ln_l
        );
        assert!(on_cards.accelerator.is_some());
    }
}
