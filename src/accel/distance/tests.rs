use super::*;
use crate::accel::pool::Member;
use crate::seq::msa::Msa;
use crate::seq::tree::{cpu_distances, Rng};

fn cpu(weight: f64) -> Member {
    Member::new(Device::Cpu, "CPU".into(), weight)
}

/// The pools to try: one device, two unequal ones, three (so that the work does not divide), and the real GPU of this machine
/// beside a CPU, when there is one.
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

/// `n` sequences of `columns` columns: a common ancestor, each copy changed at about `divergence` of its columns (and with gaps,
/// ambiguity codes and lower case here and there, which the distances must treat as the processor does).
pub(crate) fn alignment(n: usize, columns: usize, divergence: f64, seed: u64) -> Vec<Vec<u8>> {
    let mut rng = Rng::new(seed);
    let letters = b"ACGT";
    let ancestor: Vec<u8> = (0..columns).map(|_| letters[rng.below(4)]).collect();
    (0..n)
        .map(|_| {
            ancestor
                .iter()
                .map(|&c| {
                    let roll = rng.below(1000);
                    let mut x = c;
                    if (roll as f64) < divergence * 1000.0 {
                        x = letters[rng.below(4)];
                    }
                    match rng.below(40) {
                        0 => b'-',
                        1 => b'N',
                        2 => b'R',
                        3 => x.to_ascii_lowercase(),
                        _ => x,
                    }
                })
                .collect()
        })
        .collect()
}

pub(crate) fn msa(rows: &[Vec<u8>]) -> Msa {
    Msa {
        names: (0..rows.len()).map(|k| format!("s{k}")).collect(),
        rows: rows.to_vec(),
    }
}

#[test]
fn the_gpu_distances_are_those_of_the_processor_to_the_last_digit() {
    let rows = alignment(13, 777, 0.15, 7);
    for pool in pools() {
        let engine = Engine::new(&pool, &rows).expect("an engine");
        let counts = engine.counts(None).expect("counts");
        for model in DistanceModel::ALL {
            let want = cpu_distances(&msa(&rows), model).expect("distances");
            let got = counts.distances(model).expect("distances");
            // the same whole numbers through the same formula: the same floats, not "close"
            assert_eq!(got, want, "{} on {} cards", model.label(), pool.len());
        }
    }
}

#[test]
fn columns_taken_in_chunks_of_any_size_give_the_same_counts() {
    let rows = alignment(9, 203, 0.2, 11);
    let pool = Pool::new(vec![cpu(2.0), cpu(1.0)]);
    let whole = Engine::new(&pool, &rows)
        .expect("an engine")
        .counts(None)
        .expect("counts");
    // one column at a time, chunks that do not divide the columns, one chunk longer than they are
    for chunk in [1, 2, 37, 202, 203, 204, 1000] {
        let counts = Engine::with_chunk(&pool, &rows, Some(chunk))
            .expect("an engine")
            .counts(None)
            .expect("counts");
        assert_eq!(counts.sites, whole.sites, "chunk {chunk}");
        assert_eq!(counts.transitions, whole.transitions, "chunk {chunk}");
        assert_eq!(counts.transversions, whole.transversions, "chunk {chunk}");
    }
}

#[test]
fn a_bootstrap_resample_is_a_weight_for_each_column_and_gives_the_distances_of_the_resampled_alignment(
) {
    let rows = alignment(8, 400, 0.12, 21);
    let width = rows[0].len();
    let mut rng = Rng::new(99);
    for pool in pools() {
        let engine = Engine::with_chunk(&pool, &rows, Some(64)).expect("an engine");
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
                    "{} on {} cards",
                    model.label(),
                    pool.len()
                );
            }
        }
    }
}

#[test]
fn a_pair_that_cannot_be_compared_is_named_as_the_processor_names_it() {
    let pool = Pool::new(vec![cpu(1.0), cpu(1.0)]);
    // no column has a base in both the first and the second
    let rows: Vec<Vec<u8>> = ["ACGT----", "----ACGT", "ACGTACGT"]
        .iter()
        .map(|r| r.as_bytes().to_vec())
        .collect();
    let counts = Engine::new(&pool, &rows)
        .expect("an engine")
        .counts(None)
        .expect("counts");
    assert_eq!(counts.distances(DistanceModel::P), Err((0, 1)));
    assert!(cpu_distances(&msa(&rows), DistanceModel::P)
        .unwrap_err()
        .contains("s0 and s1"));
    // the sequences differ too much for the Jukes-Cantor logarithm: p = 1 is past 3/4
    let far: Vec<Vec<u8>> = ["AAAA", "CCCC"]
        .iter()
        .map(|r| r.as_bytes().to_vec())
        .collect();
    let counts = Engine::new(&pool, &far)
        .expect("an engine")
        .counts(None)
        .expect("counts");
    assert_eq!(counts.distances(DistanceModel::JukesCantor), Err((0, 1)));
    assert_eq!(
        counts.distances(DistanceModel::P),
        Ok(vec![vec![0.0, 1.0], vec![1.0, 0.0]])
    );
}

#[test]
fn what_cannot_be_counted_is_refused_with_a_reason() {
    let pool = Pool::new(vec![cpu(1.0)]);
    let rows = alignment(4, 50, 0.1, 3);
    assert!(Engine::new(&pool, &rows[..1]).is_err());
    assert!(Engine::new(&pool, &[]).is_err());
    assert!(Engine::new(&pool, &[Vec::new(), Vec::new()]).is_err());
    let mut ragged = rows.clone();
    ragged[2].pop();
    assert!(Engine::new(&pool, &ragged)
        .unwrap_err()
        .contains("same length"));
    let engine = Engine::new(&pool, &rows).expect("an engine");
    assert!(engine
        .counts(Some(&[1; 49]))
        .unwrap_err()
        .contains("one weight"));
    let mut heavy = vec![1u32; 50];
    heavy[3] = MOST_WEIGHT + 1;
    assert!(engine
        .counts(Some(&heavy))
        .unwrap_err()
        .contains("more than"));
}

#[test]
fn equal_sequences_are_at_distance_nothing_and_u_is_not_a_base_as_for_the_processor() {
    let pool = Pool::new(vec![cpu(1.0), cpu(1.0)]);
    let rows: Vec<Vec<u8>> = ["ACGTACGT", "ACGTACGT", "ACGUACGU"]
        .iter()
        .map(|r| r.as_bytes().to_vec())
        .collect();
    let counts = Engine::new(&pool, &rows)
        .expect("an engine")
        .counts(None)
        .expect("counts");
    let got = counts.distances(DistanceModel::P).expect("distances");
    assert_eq!(got[0][1], 0.0);
    assert_eq!(Ok(got), cpu_distances(&msa(&rows), DistanceModel::P));
}

#[test]
fn many_threads_using_one_card_at_once_all_get_the_right_counts() {
    // candle's Metal backend gave wrong answers when several threads used one device together; the card of a pool is held by one
    // thread at a time, and this is the test that fails if that stops. (The real GPU if there is one, else the processor.)
    let pool = Pool::global()
        .cloned()
        .unwrap_or_else(|| Pool::new(vec![cpu(1.0)]));
    let rows = alignment(11, 641, 0.15, 5);
    let want = cpu_distances(&msa(&rows), DistanceModel::P).expect("distances");
    std::thread::scope(|s| {
        let jobs: Vec<_> = (0..8)
            .map(|_| {
                s.spawn(|| {
                    (0..4).all(|_| {
                        let engine =
                            Engine::with_chunk(&pool, &rows, Some(100)).expect("an engine");
                        engine
                            .counts(None)
                            .expect("counts")
                            .distances(DistanceModel::P)
                            == Ok(want.clone())
                    })
                })
            })
            .collect();
        for job in jobs {
            assert!(
                job.join().expect("a thread"),
                "a thread got other distances than the processor"
            );
        }
    });
}
