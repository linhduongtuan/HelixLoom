use super::*;

fn zurich() -> NwTree {
    NwTree::parse_all(&zurich_tree().unwrap())
        .unwrap()
        .remove(0)
}

#[test]
fn the_encoding_reproduces_the_rescale_factor_of_the_authors_notebook() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let enc = encode(&zurich()).unwrap();
    // CI_HIV.ipynb: test_rescale = 5.174108653956917
    assert!(
        (enc.rescale - 5.174108653956917).abs() < 1e-9,
        "{}",
        enc.rescale
    );
    assert_eq!(enc.rows.len(), 399);
    let root = enc.rows[0];
    assert_eq!((root[1], root[2], root[3]), (0.0, 0.0, 200.0));
    assert_eq!((root[4], root[9]), (-1.0, -1.0));
    // every tip row has -1 for its children
    assert_eq!(
        enc.rows
            .iter()
            .filter(|r| r[3] == 1.0 && r[5] == -1.0 && r[10] == -1.0)
            .count(),
        200
    );
}

#[test]
fn the_network_input_has_two_padded_channels_and_the_sampling_probability() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let enc = encode(&zurich()).unwrap();
    let x = model_input(&enc, 0.25).unwrap();
    assert_eq!(x.len(), 500 * 19 * 2);
    let at = |row: usize, feature: usize, channel: usize| x[(row * 19 + feature) * 2 + channel];
    // 200 tips fill the first channel, 199 inner nodes and one copy fill the second
    assert_eq!((at(199, 3, 0), at(200, 3, 0)), (1.0, 0.0));
    assert_eq!(
        (at(198, 3, 1) > 1.0, at(199, 3, 1) > 1.0, at(200, 3, 1)),
        (true, true, 0.0)
    );
    assert_eq!(
        (at(0, 18, 0), at(0, 18, 1), at(300, 18, 0)),
        (0.25, 0.25, 0.0)
    );
}

// Rows of the authors' own encoder (Encoding/PhyloCNN_Encoding_PhyloDyn.py, ete3) for this tree, as it
// writes them to its float32 file: the root, two inner nodes, and three tips.
#[test]
fn the_encoding_matches_the_authors_encoder() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let enc = encode(&zurich()).unwrap();
    let theirs: [(usize, [f64; FEATURES]); 6] = [
        (
            0,
            [
                7.860106, 0.0, 0.0, 200.0, -1.0, 1.6097466, 0.18655863, 8.0, 192.0, -1.0,
                6.0723586, 2.4993448, 1.1124566, 0.8548161, 1.0, 7.0, 89.0, 103.0,
            ],
        ),
        (
            2,
            [
                6.2503595, 1.6097466, 1.6097466, 8.0, -1.0, 4.4626117, 0.88959837, 1.0, 7.0, -1.0,
                -1.0, -1.0, 2.8352785, 0.9085393, -1.0, -1.0, 1.0, 6.0,
            ],
        ),
        (
            3,
            [
                7.00529, 0.8548161, 0.6682574, 103.0, 0.8548161, 0.7770147, 0.24375501, 6.0, 97.0,
                -1.0, 3.549622, 1.7110063, 0.43974623, 0.56075287, 1.0, 5.0, 14.0, 83.0,
            ],
        ),
        (
            150,
            [
                3.3145807, 4.5455256, 1.9163785, 1.0, 2.010746, -1.0, -1.0, -1.0, -1.0, 2.8783963,
                -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
            ],
        ),
        (
            300,
            [
                1.6485932, 6.2115126, 1.3006096, 1.0, 1.6192442, -1.0, -1.0, -1.0, -1.0, 2.708724,
                -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
            ],
        ),
        (
            398,
            [
                1.6099391,
                6.250167,
                0.084330015,
                1.0,
                0.6062148,
                -1.0,
                -1.0,
                -1.0,
                -1.0,
                0.792053,
                -1.0,
                -1.0,
                -1.0,
                -1.0,
                -1.0,
                -1.0,
                -1.0,
                -1.0,
            ],
        ),
    ];
    for (at, row) in theirs {
        for (feature, want) in row.into_iter().enumerate() {
            let got = enc.rows[at][feature];
            assert!(
                (got - want).abs() < 1e-6 * want.abs().max(1.0),
                "row {at} feature {feature}: {got} vs {want}"
            );
        }
    }
}

// The output of the authors' Keras layers (weights from their .h5 files) on exactly the network input
// this module builds for the Zurich tree at sampling probability 0.25, one line per trained network.
#[test]
fn every_network_matches_keras_on_the_same_input() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let enc = encode(&zurich()).unwrap();
    let x = model_input(&enc, 0.25).unwrap();
    let keras: [(Task, Priors, &[f32]); 6] = [
        (Task::Bd, Priors::Simulated, &[1.2909807, 1.1191297]),
        (
            Task::Bdei,
            Priors::Simulated,
            &[1.2124276, 0.4007007, 0.8226939],
        ),
        (
            Task::Bdss,
            Priors::Simulated,
            &[1.1770825, 1.5029598, 9.531764, 0.08032211],
        ),
        (
            Task::ChooseModel,
            Priors::Simulated,
            &[4.3746075e-4, 3.9996526e-10, 0.99956256],
        ),
        (
            Task::Bdss,
            Priors::Hiv,
            &[1.4004343, 1.7627712, 19.707214, 0.12167474],
        ),
        (
            Task::ChooseModel,
            Priors::Hiv,
            &[2.8337555e-8, 8.479723e-12, 1.0],
        ),
    ];
    for (task, priors, want) in keras {
        let net = Network::parse(&network_bytes(task, priors).unwrap()).unwrap();
        let got = net.run(&x).unwrap();
        assert_eq!(got.len(), want.len(), "{task:?} {priors:?}");
        for (g, w) in got.iter().zip(want) {
            assert!(
                (g - w).abs() < 1e-6 + 1e-4 * w.abs(),
                "{task:?} {priors:?}: {got:?} vs {want:?}"
            );
        }
    }
}

// The authors' pipeline end to end (their encoder, then their Keras layers, as in their CI_HIV notebook
// but without the rescale row that training dropped) gives 1.4019606, 1.7617304 x 5.1741085, 19.73183 and
// 0.12139099. Equal-age tips can be ordered in more than one way, and the two channels are paired row by
// row, so a different order moves the answer by about a tenth of a percent.
#[test]
fn the_estimate_agrees_with_the_authors_pipeline_within_a_quarter_percent() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let got = estimate(&zurich(), 0.25, Task::Bdss, Priors::Hiv).unwrap();
    let want = [1.4019606, 1.7617304 * 5.1741085, 19.73183, 0.12139099];
    let values: Vec<f64> = got.iter().map(|e| e.value).collect();
    for (g, w) in values.iter().zip(want) {
        assert!((g - w).abs() < 3e-3 * w, "{values:?} vs {want:?}");
    }
}

#[test]
fn model_choice_gives_probabilities_and_other_priors_are_refused() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let p = estimate(&zurich(), 0.25, Task::ChooseModel, Priors::Hiv).unwrap();
    assert_eq!(p.len(), 3);
    assert!((p.iter().map(|e| e.value).sum::<f64>() - 1.0).abs() < 1e-5);
    assert!(estimate(&zurich(), 0.25, Task::Bd, Priors::Hiv).is_err());
    assert!(estimate(&zurich(), 0.0, Task::Bd, Priors::Simulated).is_err());
    for task in Task::ALL {
        assert!(
            estimate(&zurich(), 0.25, task, Priors::Simulated).is_ok(),
            "{task:?}"
        );
    }
    let polytomy = NwTree::parse_all("((A:1,B:1,C:1):1,D:2);")
        .unwrap()
        .remove(0);
    assert!(encode(&polytomy).is_err());
}

/// The Zurich tree with its branch lengths shaken a little, at several sampling probabilities: varied
/// inputs for comparing two ways of running a network.
#[cfg(feature = "gpu")]
fn varied_inputs() -> Vec<Vec<f32>> {
    let mut inputs = Vec::new();
    for (k, sampling) in [0.05, 0.1, 0.25, 0.5, 0.9, 1.0].into_iter().enumerate() {
        let mut tree = zurich();
        for (i, node) in tree.nodes.iter_mut().enumerate() {
            node.dist *= 1.0 + 0.3 * ((i + 7 * k) as f64).sin();
        }
        inputs.push(model_input(&encode(&tree).unwrap(), sampling).unwrap());
    }
    inputs
}

#[cfg(feature = "gpu")]
const ALL_NETWORKS: [(Task, Priors); 6] = [
    (Task::Bd, Priors::Simulated),
    (Task::Bdei, Priors::Simulated),
    (Task::Bdss, Priors::Simulated),
    (Task::ChooseModel, Priors::Simulated),
    (Task::Bdss, Priors::Hiv),
    (Task::ChooseModel, Priors::Hiv),
];

#[cfg(feature = "gpu")]
fn same_outputs(device: &candle_core::Device, tolerance: f32) {
    let inputs = varied_inputs();
    for (task, priors) in ALL_NETWORKS {
        let net = Network::parse(&network_bytes(task, priors).unwrap()).unwrap();
        let got = gpu::run_batch_on(&net, &inputs, device).unwrap();
        for (x, g) in inputs.iter().zip(&got) {
            let want = net.run(x).unwrap();
            assert_eq!(g.len(), want.len());
            for (a, b) in g.iter().zip(&want) {
                assert!(
                    (a - b).abs() <= tolerance * b.abs().max(1.0),
                    "{task:?} {priors:?}: {g:?} vs {want:?}"
                );
            }
        }
    }
}

#[cfg(feature = "gpu")]
#[test]
fn candle_on_the_cpu_gives_the_plain_networks_answers() {
    if crate::models::phylocnn_missing() {
        return;
    }
    same_outputs(&candle_core::Device::Cpu, 1e-4);
}

#[cfg(feature = "gpu")]
#[test]
#[ignore = "builds the GPU's shaders, about fifteen seconds; run with --ignored"]
fn the_gpu_gives_the_plain_networks_answers() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let Some((device, name)) = gpu::device() else {
        eprintln!("no GPU here");
        return;
    };
    eprintln!("checking {name}");
    same_outputs(&device, 1e-3);
}

#[test]
fn many_trees_give_each_trees_own_estimate_and_unusable_trees_are_reported() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let mut trees = vec![zurich(), zurich()];
    trees[1].nodes[3].dist *= 2.0;
    trees.push(
        NwTree::parse_all("((A:1,B:1,C:1):1,D:2);")
            .unwrap()
            .remove(0),
    );
    let batch = estimate_many(&trees, 0.25, Task::Bdss, Priors::Hiv, Engine::Cpu).unwrap();
    assert_eq!(batch.per_tree.len(), 2);
    assert_eq!(batch.skipped.len(), 1);
    assert_eq!(batch.skipped[0].0, 2);
    for (tree, got) in trees.iter().zip(&batch.per_tree) {
        let alone = estimate(tree, 0.25, Task::Bdss, Priors::Hiv).unwrap();
        assert_eq!(got, &alone);
    }
    assert!(estimate_many(&trees[2..], 0.25, Task::Bd, Priors::Simulated, Engine::Cpu).is_err());
}

#[test]
fn the_summary_of_many_estimates_is_the_mean_the_median_and_the_middle_95_percent() {
    let one = |v: f64| {
        vec![Estimate {
            name: "x",
            value: v,
        }]
    };
    let trees: Vec<_> = (1..=5).map(|v| one(f64::from(v))).collect();
    let s = &summarize(&trees)[0];
    assert_eq!((s.mean, s.median), (3.0, 3.0));
    assert!((s.low - 1.1).abs() < 1e-12 && (s.high - 4.9).abs() < 1e-12);
    assert!(summarize(&[]).is_empty());
}

#[cfg(feature = "gpu")]
#[test]
#[ignore = "a benchmark; run with --ignored --nocapture"]
fn how_long_a_thousand_trees_take_on_the_cpu_and_on_the_gpu() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let trees: Vec<NwTree> = (0..1000)
        .map(|k| {
            let mut t = zurich();
            for (i, n) in t.nodes.iter_mut().enumerate() {
                n.dist *= 1.0 + 0.3 * ((i + 7 * k) as f64).sin();
            }
            t
        })
        .collect();
    let run = |engine| estimate_many(&trees, 0.25, Task::Bdss, Priors::Hiv, engine).unwrap();
    let cpu = run(Engine::Cpu);
    eprintln!("1000 trees, CPU on every core: {:.2} s", cpu.seconds);
    let first = estimate_many(&trees[..1], 0.25, Task::Bdss, Priors::Hiv, Engine::Gpu).unwrap();
    eprintln!(
        "GPU, first use (builds its shaders): {:.2} s",
        first.seconds
    );
    let gpu = run(Engine::Gpu);
    eprintln!("1000 trees, GPU: {:.2} s", gpu.seconds);
    let worst = cpu
        .per_tree
        .iter()
        .flatten()
        .zip(gpu.per_tree.iter().flatten())
        .map(|(a, b)| (a.value - b.value).abs() / a.value.abs().max(1.0))
        .fold(0.0, f64::max);
    eprintln!("largest relative difference between them: {worst:.2e}");
}
