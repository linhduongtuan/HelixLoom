use super::*;
use crate::seq::msa::center_star;
use std::collections::HashMap;

struct Rng(u32);

impl Rng {
    fn next(&mut self) -> u32 {
        self.0 = self.0.wrapping_mul(1_103_515_245).wrapping_add(12_345);
        self.0 >> 16
    }

    fn below(&mut self, n: usize) -> usize {
        self.next() as usize % n
    }

    fn chance(&mut self, p: f64) -> bool {
        f64::from(self.next() % 10_000) < p * 10_000.0
    }
}

fn letters(kind: Kind) -> &'static [u8] {
    match kind {
        Kind::Dna => NUCLEO,
        Kind::Protein => AMINO,
    }
}

/// A sequence of residues, each with the number of the column of the ancestor it came from (new numbers for what was inserted).
type Evolved = Vec<(u32, u8)>;

/// What `residue` changes to: another nucleotide, any of the three; another amino acid, with the odds that BLOSUM62 stands for (a
/// change to a similar residue is likelier than one to a different residue), as real proteins change.
fn substitute(residue: u8, kind: Kind, rng: &mut Rng) -> u8 {
    let alphabet = letters(kind);
    let from = alphabet.iter().position(|&l| l == residue).unwrap_or(0);
    let weight = |to: usize| -> f64 {
        match kind {
            _ if to == from => 0.0,
            Kind::Dna => 1.0,
            Kind::Protein => 2f64.powf(f64::from(BLOSUM62[from][to]) / 2.0),
        }
    };
    let total: f64 = (0..alphabet.len()).map(weight).sum();
    let mut left = f64::from(rng.next() % 100_000) / 100_000.0 * total;
    for (to, &letter) in alphabet.iter().enumerate() {
        left -= weight(to);
        if left < 0.0 {
            return letter;
        }
    }
    alphabet[(from + 1) % alphabet.len()]
}

/// One branch of the tree: substitutions, and short insertions and deletions, at the rate `rate` for each residue.
fn mutate(parent: &Evolved, rate: f64, kind: Kind, rng: &mut Rng, next_id: &mut u32) -> Evolved {
    let alphabet = letters(kind);
    let mut out = Vec::with_capacity(parent.len());
    let mut skip = 0;
    for &(id, residue) in parent {
        if skip > 0 {
            skip -= 1;
            continue;
        }
        if rng.chance(rate / 4.0) {
            // a deletion of one to three residues
            skip = rng.below(3);
            continue;
        }
        let kept = if rng.chance(rate) {
            substitute(residue, kind, rng)
        } else {
            residue
        };
        out.push((id, kept));
        if rng.chance(rate / 4.0) {
            for _ in 0..=rng.below(3) {
                *next_id += 1;
                out.push((*next_id, alphabet[rng.below(alphabet.len())]));
            }
        }
    }
    out
}

/// The sequences that were made, and for each row the number of the column each residue belongs to (the true alignment).
type Simulated = (Vec<(String, Vec<u8>)>, Vec<Vec<u32>>);

/// Sequences that evolved along a balanced tree of `leaves` leaves from one ancestor, with the numbers that say which residues are
/// the same column (the true alignment).
fn simulate(kind: Kind, leaves: usize, length: usize, rate: f64, seed: u32) -> Simulated {
    let mut rng = Rng(seed);
    let alphabet = letters(kind);
    let ancestor: Evolved = (0..length)
        .map(|k| (k as u32, alphabet[rng.below(alphabet.len())]))
        .collect();
    let mut next_id = length as u32;
    let mut level = vec![ancestor];
    while level.len() < leaves {
        level = level
            .iter()
            .flat_map(|p| {
                [
                    mutate(p, rate, kind, &mut rng, &mut next_id),
                    mutate(p, rate, kind, &mut rng, &mut next_id),
                ]
            })
            .collect();
    }
    let seqs = level
        .iter()
        .enumerate()
        .map(|(k, e)| (format!("seq{k}"), e.iter().map(|&(_, r)| r).collect()))
        .collect();
    (
        seqs,
        level
            .iter()
            .map(|e| e.iter().map(|&(id, _)| id).collect())
            .collect(),
    )
}

/// The share of the pairs of residues that are one column in the true alignment which are one column in `rows` (how much of the
/// truth is found), and the share of the pairs that are one column in `rows` which are one in the truth (how much of what is
/// found is true).
fn accuracy(rows: &[Vec<u8>], truth: &[Vec<u32>]) -> (f64, f64) {
    let columns: Vec<Vec<usize>> = rows
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .filter(|(_, &b)| b != GAP)
                .map(|(c, _)| c)
                .collect()
        })
        .collect();
    let (mut true_pairs, mut found_true, mut found) = (0usize, 0usize, 0usize);
    for i in 0..rows.len() {
        let by_id: HashMap<u32, usize> = truth[i]
            .iter()
            .enumerate()
            .map(|(k, &id)| (id, k))
            .collect();
        for j in i + 1..rows.len() {
            let mut theirs: HashMap<usize, usize> = HashMap::new();
            for (k, &c) in columns[j].iter().enumerate() {
                theirs.insert(c, k);
            }
            for (k, &id) in truth[j].iter().enumerate() {
                if by_id.contains_key(&id) {
                    true_pairs += 1;
                    if columns[i][by_id[&id]] == columns[j][k] {
                        found_true += 1;
                    }
                }
            }
            for (a, &c) in columns[i].iter().enumerate() {
                if let Some(&b) = theirs.get(&c) {
                    found += 1;
                    let _ = (a, b);
                }
            }
        }
    }
    (
        found_true as f64 / true_pairs.max(1) as f64,
        found_true as f64 / found.max(1) as f64,
    )
}

fn strip(row: &[u8]) -> Vec<u8> {
    row.iter().copied().filter(|&b| b != GAP).collect()
}

fn check_alignment(seqs: &[(String, Vec<u8>)], msa: &Msa) {
    assert_eq!(
        msa.names,
        seqs.iter().map(|s| s.0.clone()).collect::<Vec<_>>()
    );
    assert_eq!(msa.rows.len(), seqs.len());
    let len = msa.rows[0].len();
    for (row, (_, original)) in msa.rows.iter().zip(seqs) {
        assert_eq!(row.len(), len, "rows of one length");
        assert_eq!(&strip(row), original, "a row is its sequence with gaps");
    }
    for c in 0..len {
        assert!(
            msa.rows.iter().any(|r| r[c] != GAP),
            "column {c} is all gaps"
        );
    }
}

#[test]
fn blosum62_is_the_published_matrix() {
    for a in 0..20 {
        for b in 0..20 {
            assert_eq!(
                BLOSUM62[a][b], BLOSUM62[b][a],
                "{} and {}",
                AMINO[a] as char, AMINO[b] as char
            );
        }
    }
    let diagonal: Vec<i8> = (0..20).map(|k| BLOSUM62[k][k]).collect();
    assert_eq!(
        diagonal,
        [4, 5, 6, 6, 9, 5, 5, 6, 8, 4, 4, 5, 5, 6, 7, 4, 5, 11, 7, 4]
    );
    let at = |a: u8, b: u8| {
        let place = |x: u8| AMINO.iter().position(|&l| l == x).expect("an amino acid");
        BLOSUM62[place(a)][place(b)]
    };
    // entries that the published matrix is known by
    for (a, b, want) in [
        (b'W', b'F', 1),
        (b'W', b'Y', 2),
        (b'I', b'V', 3),
        (b'I', b'L', 2),
        (b'L', b'V', 1),
        (b'K', b'R', 2),
        (b'D', b'E', 2),
        (b'F', b'Y', 3),
        (b'H', b'Y', 2),
        (b'C', b'W', -2),
        (b'P', b'W', -4),
        (b'D', b'W', -4),
        (b'A', b'R', -1),
        (b'S', b'T', 1),
        (b'N', b'D', 1),
        (b'E', b'Q', 2),
    ] {
        assert_eq!(at(a, b), want, "{} and {}", a as char, b as char);
    }
}

#[test]
fn every_alignment_is_its_sequences_with_gaps_in_rows_of_one_length_with_no_empty_column() {
    for (kind, seed) in [(Kind::Dna, 7), (Kind::Protein, 11), (Kind::Dna, 23)] {
        let (seqs, _) = simulate(kind, 8, 90, 0.05, seed);
        let report = align(&seqs, &Options::default()).expect("aligned");
        check_alignment(&seqs, &report.msa);
        assert!(!report.by_words);
        // the same input gives the same answer
        assert_eq!(
            align(&seqs, &Options::default()).expect("aligned").msa,
            report.msa
        );
    }
}

#[test]
fn a_gap_sign_in_a_sequence_is_ignored_and_what_cannot_be_aligned_is_said() {
    let seqs = vec![
        ("a".to_string(), b"AC-GTACGTAC".to_vec()),
        ("b".to_string(), b"ACGTACGTAC".to_vec()),
    ];
    let report = align(&seqs, &Options::default()).expect("aligned");
    assert_eq!(report.msa.rows[0], report.msa.rows[1]);
    assert!(align(&[], &Options::default()).is_err());
    assert!(align(
        &[("a".into(), b"ACGT".to_vec()), ("b".into(), Vec::new())],
        &Options::default()
    )
    .is_err());
    let one = align(
        &[("only".to_string(), b"ACGT".to_vec())],
        &Options::default(),
    )
    .expect("one");
    assert_eq!(one.msa.rows, [b"ACGT".to_vec()]);
    let many: Vec<(String, Vec<u8>)> = (0..=MAX_SEQUENCES)
        .map(|k| (format!("s{k}"), b"ACGTACGT".to_vec()))
        .collect();
    let why = align(&many, &Options::default()).expect_err("too many");
    assert!(why.contains("2000"), "{why}");
}

#[test]
fn identical_sequences_align_without_a_gap_and_one_deleted_base_makes_one_gap() {
    let seq = b"ACGTTGCAAGCTTGACCAGTA".to_vec();
    let same = vec![
        ("a".to_string(), seq.clone()),
        ("b".to_string(), seq.clone()),
        ("c".to_string(), seq.clone()),
    ];
    for guide in [Guide::NeighborJoining, Guide::Upgma] {
        let r = align(
            &same,
            &Options {
                guide,
                ..Options::default()
            },
        )
        .expect("aligned");
        assert!(r.msa.rows.iter().all(|row| row == &seq), "{guide:?}");
    }
    let mut shorter = seq.clone();
    shorter.remove(10);
    let three = vec![
        ("a".to_string(), seq.clone()),
        ("b".to_string(), shorter.clone()),
        ("c".to_string(), seq),
    ];
    let r = align(&three, &Options::default()).expect("aligned");
    check_alignment(&three, &r.msa);
    assert_eq!(r.msa.rows[1].iter().filter(|&&b| b == GAP).count(), 1);
    assert_eq!(r.msa.columns(), 21);
}

#[test]
fn a_fragment_of_a_sequence_sits_inside_it_with_no_gap_in_it_when_the_ends_are_free() {
    let whole = b"GATTACAGGCTTAACCGGTTAAGCTTGGATCCAAGCTT".to_vec();
    let part = whole[8..28].to_vec();
    let seqs = vec![
        ("whole".to_string(), whole.clone()),
        ("part".to_string(), part.clone()),
    ];
    let r = align(&seqs, &Options::default()).expect("aligned");
    let row = &r.msa.rows[1];
    let first = row.iter().position(|&b| b != GAP).expect("residues");
    let last = row.iter().rposition(|&b| b != GAP).expect("residues");
    assert!(
        row[first..=last].iter().all(|&b| b != GAP),
        "{}",
        String::from_utf8_lossy(row)
    );
    assert_eq!(first, 8);
    assert_eq!(&r.msa.rows[0][first..=last], &part[..]);
}

#[test]
fn the_alignment_is_close_to_the_true_one_for_sequences_that_evolved_with_substitutions_and_indels()
{
    // (kind, per-branch rate, the least share of the true pairs that must be found, and of the pairs found that must be true):
    // about 6 and 12 residues in 100 differ along a path of three branches; measured with other seeds than these, the plain method
    // finds 0.93 (nucleotides) and 0.88 (proteins) and the refined one 0.95 and 0.92
    let mut report = Vec::new();
    for (kind, rate, want) in [(Kind::Dna, 0.02, 0.88), (Kind::Protein, 0.04, 0.82)] {
        let (mut found, mut true_share, mut refined, mut centre) = (0.0, 0.0, 0.0, 0.0);
        let seeds = [3u32, 17, 29, 41];
        let n = seeds.len() as f64;
        for &seed in &seeds {
            let (seqs, truth) = simulate(kind, 8, 110, rate, seed);
            let r = align(&seqs, &Options::default()).expect("aligned");
            check_alignment(&seqs, &r.msa);
            let (f, t) = accuracy(&r.msa.rows, &truth);
            (found, true_share) = (found + f / n, true_share + t / n);
            let better = align(
                &seqs,
                &Options {
                    refine: 3,
                    ..Options::default()
                },
            )
            .expect("aligned");
            refined += accuracy(&better.msa.rows, &truth).0 / n;
            centre += accuracy(&center_star(&seqs).expect("aligned").rows, &truth).0 / n;
        }
        report.push(format!(
            "{kind:?}: found {found:.3}, true {true_share:.3}, refined {refined:.3}, centre-star {centre:.3}"
        ));
        assert!(found >= want, "{report:?}");
        assert!(true_share >= want, "{report:?}");
        // refining does not lose what the first alignment found, to within what the seeds can move
        assert!(refined + 0.02 >= found, "{report:?}");
        if kind == Kind::Dna {
            // and with refinement the method does about as well as aligning everything to one sequence
            assert!(refined + 0.03 >= centre, "{report:?}");
        }
    }
}

#[test]
fn refinement_never_lowers_the_score_and_never_changes_what_a_row_is() {
    for (kind, seed) in [(Kind::Protein, 5u32), (Kind::Dna, 13), (Kind::Protein, 31)] {
        let (seqs, _) = simulate(kind, 8, 80, 0.07, seed);
        let plain = align(&seqs, &Options::default()).expect("aligned");
        let refined = align(
            &seqs,
            &Options {
                refine: 3,
                ..Options::default()
            },
        )
        .expect("aligned");
        check_alignment(&seqs, &refined.msa);
        let (before, after) = (
            plain.score.expect("a score"),
            refined.score.expect("a score"),
        );
        assert_eq!(refined.unrefined_score, plain.score);
        assert!(after >= before - 1e-9, "{before} then {after}");
        assert!(!refined.refine_skipped);
    }
    // two sequences have nothing to refine, and too many say so
    let two = vec![
        ("a".to_string(), b"ACGTACGTAC".to_vec()),
        ("b".to_string(), b"ACGTACGTAC".to_vec()),
    ];
    assert!(
        align(
            &two,
            &Options {
                refine: 2,
                ..Options::default()
            }
        )
        .expect("aligned")
        .refine_skipped
    );
}

#[test]
fn near_copies_weigh_less_than_a_sequence_that_stands_alone_and_the_weights_add_up_to_their_number()
{
    let base = b"ACGTTGCAAGCTTGACCAGTAACGGATCCATTGACCTAGGCAT".to_vec();
    let mut copy1 = base.clone();
    copy1[5] = b'A';
    let mut copy2 = base.clone();
    copy2[20] = b'C';
    let alone: Vec<u8> = b"TTTTGGGGAAAACCCCTTGGAACCTTGGAACCGGTTAACCGGAA".to_vec();
    let seqs: Vec<Vec<u8>> = vec![base, copy1, copy2, alone];
    let sc = Scores::new(Kind::Dna);
    let (d, by_words) = distances(&seqs, &sc, false).expect("distances");
    assert!(!by_words);
    for guide in [Guide::NeighborJoining, Guide::Upgma] {
        let tree = guide_tree(guide, &d).expect("a tree");
        let w = weights_from(&tree, 4);
        assert!((w.iter().sum::<f32>() - 4.0).abs() < 1e-3, "{w:?}");
        assert!(w[3] > w[0] && w[3] > w[1] && w[3] > w[2], "{guide:?} {w:?}");
    }
    // a tree with no length gives every sequence the same weight
    let flat = vec![vec![0.0; 3]; 3];
    let tree = guide_tree(Guide::NeighborJoining, &flat).expect("a tree");
    let w = weights_from(&tree, 3);
    assert!(w.iter().all(|&x| (x - 1.0).abs() < 1e-3), "{w:?}");
}

#[test]
fn the_tree_is_rooted_halfway_along_its_longest_path() {
    // leaves A and B hang 1 from one join, C hangs 1 and D hangs 5 from the other, the joins are 1 apart: the longest path is
    // A to D (7), whose middle is 3.5 from D, on D's own branch
    let d = vec![
        vec![0.0, 2.0, 3.0, 7.0],
        vec![2.0, 0.0, 3.0, 7.0],
        vec![3.0, 3.0, 0.0, 6.0],
        vec![7.0, 7.0, 6.0, 0.0],
    ];
    let tree = guide_tree(Guide::NeighborJoining, &d).expect("a tree");
    assert_eq!(tree.kids[tree.root].len(), 2);
    let leaf_d = (0..tree.kids.len())
        .find(|&v| tree.leaf[v] == Some(3))
        .expect("leaf D");
    assert!(
        tree.kids[tree.root].contains(&leaf_d),
        "D hangs from the root"
    );
    assert!(
        (tree.len[leaf_d] - 3.5).abs() < 1e-6,
        "{}",
        tree.len[leaf_d]
    );
    let other = tree.kids[tree.root]
        .iter()
        .copied()
        .find(|&v| v != leaf_d)
        .expect("the rest");
    assert!((tree.len[other] - 1.5).abs() < 1e-6, "{}", tree.len[other]);
    assert_eq!(tree.leaves_below()[other], 3);
    // every join has two sequences or groups below it
    assert!(tree.kids.iter().all(|k| k.is_empty() || k.len() == 2));
}

#[test]
fn a_gap_costs_less_where_sequences_already_have_gaps_and_more_near_one() {
    let sc = Scores::new(Kind::Dna);
    // columns: 0-2 no gap, 3 a gap in one of two, 4.. no gap again; the second gap column is far from the first
    let rows: Vec<&[u8]> = vec![b"ACGT-ACGTACGTACGTAC", b"ACGTTACGTACGTACGTAC"];
    let group = Group {
        ids: vec![0, 1],
        rows: rows.iter().map(|r| r.to_vec()).collect(),
    };
    let p = profile(&group, &[1.0, 1.0], &sc);
    // the column with a gap: the share of the weight that has none
    assert!((p.opening[4] - 0.5).abs() < 1e-6, "{}", p.opening[4]);
    // near it: 2 + (8 - distance) * 2 / 8, the distance counted in columns from the gap
    for d in 1..GAP_DISTANCE {
        let want = 2.0 + (GAP_DISTANCE - d) as f32 * 2.0 / GAP_DISTANCE as f32;
        assert!(
            (p.opening[4 + d] - want).abs() < 1e-6,
            "{d}: {}",
            p.opening[4 + d]
        );
        assert!((p.opening[4 - d.min(4)] >= 2.0) || d > 4);
    }
    // and far from it nothing changes
    assert!((p.opening[4 + GAP_DISTANCE] - 1.0).abs() < 1e-6);
    // weights count: when the sequence with the gap weighs little, the gap is nearly as dear as a new one
    let light = profile(&group, &[0.1, 1.9], &sc);
    assert!(
        (light.opening[4] - 0.95).abs() < 1e-6,
        "{}",
        light.opening[4]
    );
}

#[test]
fn a_second_sequence_with_the_same_gap_puts_it_in_the_same_columns() {
    let a = b"ACGTACGTAAAACCCGGGTTTACGTACGTCATGCA".to_vec();
    let b: Vec<u8> = [&a[..8], &a[18..]].concat();
    let seqs = vec![
        ("a".to_string(), a.clone()),
        ("b".to_string(), b.clone()),
        ("c".to_string(), b),
        ("d".to_string(), a),
    ];
    let r = align(&seqs, &Options::default()).expect("aligned");
    check_alignment(&seqs, &r.msa);
    assert_eq!(r.msa.rows[1], r.msa.rows[2]);
    let gaps: Vec<usize> = r.msa.rows[1]
        .iter()
        .enumerate()
        .filter(|(_, &c)| c == GAP)
        .map(|(k, _)| k)
        .collect();
    assert_eq!(gaps.len(), 10);
    assert_eq!(
        gaps.last().expect("a gap") - gaps[0],
        9,
        "one block: {gaps:?}"
    );
}

#[test]
fn the_fast_distances_come_from_shared_words_and_still_give_an_alignment() {
    let (seqs, truth) = simulate(Kind::Dna, 16, 400, 0.02, 99);
    let sc = Scores::new(Kind::Dna);
    let clean: Vec<Vec<u8>> = seqs.iter().map(|s| s.1.clone()).collect();
    // the words of two close sequences are mostly shared
    let (set0, set1) = (words(&clean[0], 6, &sc), words(&clean[1], 6, &sc));
    let by_words = word_distance(&set0, &set1, 6, (clean[0].len(), clean[1].len()));
    assert!((0.0..=1.0).contains(&by_words));
    // identical sequences share every word, and sequences with nothing in common none
    let w = words(b"ACGTACGTACGT", 6, &sc);
    assert_eq!(word_distance(&w, &w, 6, (12, 12)), 0.0);
    let (x, y) = (
        words(b"AAAAAAAAAAAA", 6, &sc),
        words(b"CCCCCCCCCCCC", 6, &sc),
    );
    assert_eq!(word_distance(&x, &y, 6, (12, 12)), 1.0);
    // a sequence shorter than a word has none to share
    assert_eq!(word_distance(&[], &w, 6, (3, 12)), 1.0);
    // the whole alignment with the fast distances is a good one, and says that it was made so; the slow one says it was not
    let fast = align(
        &seqs,
        &Options {
            fast: true,
            ..Options::default()
        },
    )
    .expect("aligned");
    check_alignment(&seqs, &fast.msa);
    assert!(fast.by_words);
    let (found, _) = accuracy(&fast.msa.rows, &truth);
    let slow = align(&seqs, &Options::default()).expect("aligned");
    assert!(!slow.by_words);
    let (found_slow, _) = accuracy(&slow.msa.rows, &truth);
    // sixteen sequences of 400 at four branches from each other are harder than eight of 110: the two ways are compared with each
    // other, and the fast one may do a little worse than the slow one and no more
    assert!(
        found + 0.05 >= found_slow,
        "fast {found}, slow {found_slow}"
    );
    assert!(found > 0.7, "fast {found}, slow {found_slow}");
}

fn named(rows: &[(&str, &[u8])]) -> Vec<(String, Vec<u8>)> {
    rows.iter()
        .map(|(n, s)| (n.to_string(), s.to_vec()))
        .collect()
}

#[test]
fn unrelated_sequences_that_share_a_letter_at_their_ends_are_far_apart_not_identical() {
    // with free end gaps one T against one T would be the whole of the alignment, a distance of nothing
    let dna = Scores::new(Kind::Dna);
    let d = pair_distance(b"AAAAAAAT", b"TCCCCCCC", &dna).expect("a distance");
    assert!(d > 0.5, "{d}");
    let protein = Scores::new(Kind::Protein);
    let d = pair_distance(b"GGGGGGGW", b"WLLLLLLL", &protein).expect("a distance");
    assert!(d > 0.5, "{d}");
    // and the same sequences with one letter changed are close
    let d = pair_distance(b"ACGTTGCAAGCTTGACC", b"ACGTTGCAAGCTAGACC", &dna).expect("a distance");
    assert!(d > 0.0 && d < 0.3, "{d}");
}

#[test]
fn the_distance_of_proteins_is_kimuras_correction_of_the_share_that_differs() {
    let protein = Scores::new(Kind::Protein);
    let same = pair_distance(b"ACDEFGHIKL", b"ACDEFGHIKL", &protein).expect("a distance");
    assert_eq!(same, 0.0);
    // two of ten differ: -ln(1 - D - D^2/5) for D = 0.2
    let d = pair_distance(b"ACDEFGHIKL", b"ACDEFGHIPP", &protein).expect("a distance");
    let want = -(1.0f64 - 0.2 - 0.2 * 0.2 * 0.2).ln();
    assert!((d - want).abs() < 1e-9, "{d} against {want}");
}

#[test]
fn rna_and_lower_case_are_measured_as_dna_is() {
    let dna = Scores::new(Kind::Dna);
    let t = pair_distance(b"ACGTACGTACGTAC", b"ACGTACGAACGTAC", &dna).expect("a distance");
    let u = pair_distance(b"ACGUACGUACGUAC", b"ACGUACGAACGUAC", &dna).expect("a distance");
    let lower = pair_distance(b"acguacguacguac", b"acguacgaacguac", &dna).expect("a distance");
    assert!(t > 0.0 && t < 0.5, "{t}");
    assert_eq!((t, t), (u, lower));
    // identical RNA is not at the largest distance
    assert_eq!(pair_distance(b"UUUUACGU", b"UUUUACGU", &dna), Ok(0.0));
}

#[test]
fn the_score_of_an_alignment_is_what_the_costs_say_worked_out_by_hand() {
    let dna = Scores::new(Kind::Dna);
    let sop = |a: &str, b: &str, free: bool| {
        let rows = vec![a.as_bytes().to_vec(), b.as_bytes().to_vec()];
        sum_of_pairs(&rows, &[0, 1], &[1.0, 1.0], &dna, free)
    };
    assert_eq!(sop("ACGT", "ACGT", false), 20.0); // four matches of +5
    assert_eq!(sop("ACGT", "ACGA", false), 11.0); // three matches, one mismatch of -4
    assert_eq!(sop("AC-GT", "ACAGT", false), 10.0); // four matches, one gap opened at 10
    assert_eq!(sop("AC--T", "ACAGT", false), 3.0); // three matches, a gap of two: 10 and 2
    assert_eq!(sop("-ACGT", "TACGT", false), 10.0); // an end gap costs
    assert_eq!(sop("-ACGT", "TACGT", true), 20.0); // and is free when the ends are
    assert_eq!(sop("A-CGT", "TACGT", true), 1.0); // not at the end: -4 for A against T, -10 for the gap, +15
                                                  // weights multiply: a pair of weight 2 and 0.5 counts as 1, one of weight 2 and 2 as 4
    let rows = vec![b"ACGT".to_vec(), b"ACGT".to_vec()];
    assert_eq!(sum_of_pairs(&rows, &[0, 1], &[2.0, 0.5], &dna, true), 20.0);
    assert_eq!(sum_of_pairs(&rows, &[0, 1], &[2.0, 2.0], &dna, true), 80.0);
}

#[test]
fn a_fragment_is_placed_where_it_belongs_when_the_ends_cost_too() {
    let mut rng = Rng(5);
    let whole: Vec<u8> = (0..60).map(|_| b"ACGT"[rng.below(4)]).collect();
    let part = whole[10..40].to_vec();
    let seqs = named(&[("whole", &whole), ("part", &part)]);
    for free_end_gaps in [true, false] {
        let o = Options {
            free_end_gaps,
            ..Options::default()
        };
        let r = align(&seqs, &o).expect("aligned");
        check_alignment(&seqs, &r.msa);
        assert_eq!(r.msa.columns(), 60, "{free_end_gaps}");
        let row = &r.msa.rows[1];
        assert_eq!(&row[10..40], &part[..], "{free_end_gaps}");
        assert!(row[..10].iter().chain(&row[40..]).all(|&b| b == GAP));
    }
}

#[test]
fn refinement_improves_the_score_of_some_alignments_and_not_only_leaves_it_alone() {
    let (mut improved, mut tried) = (0, 0);
    for kind in [Kind::Protein, Kind::Dna] {
        for seed in 1..=12 {
            let (seqs, _) = simulate(
                kind,
                8,
                110,
                if kind == Kind::Dna { 0.04 } else { 0.06 },
                seed,
            );
            let r = align(
                &seqs,
                &Options {
                    refine: 3,
                    ..Options::default()
                },
            )
            .expect("aligned");
            tried += 1;
            let (after, before) = (
                r.score.expect("a score"),
                r.unrefined_score.expect("a score"),
            );
            assert!(
                after >= before - 1e-9,
                "{kind:?} seed {seed}: {before} then {after}"
            );
            if after > before + 1e-6 {
                improved += 1;
            }
        }
    }
    assert!(improved > 0, "none of {tried} was improved by refining");
}

#[test]
fn a_halving_of_the_tree_that_is_too_large_to_align_is_left_alone_and_the_alignment_is_kept() {
    // three rows of 7,000 columns: any one against the other two would be 7,001 x 7,001 cells, more than one piece may hold
    let group = Group {
        ids: vec![0, 1, 2],
        rows: vec![vec![b'A'; 7000]; 3],
    };
    let tree = Tree {
        kids: vec![vec![1, 2], vec![], vec![3, 4], vec![], vec![]],
        leaf: vec![None, Some(0), None, Some(1), Some(2)],
        len: vec![0.0; 5],
        root: 0,
    };
    let sc = Scores::new(Kind::Dna);
    let out = refine(group.clone(), &tree, &[1.0; 3], &sc, true, 2).expect("the alignment is kept");
    assert_eq!(out.rows, group.rows);
}

#[test]
fn what_is_said_of_refinement_is_what_happened_to_it() {
    let two = named(&[("a", b"ACGTACGTAC"), ("b", b"ACGTACGTAC")]);
    let one = named(&[("a", b"ACGTACGTAC")]);
    let o = Options {
        refine: 2,
        ..Options::default()
    };
    for seqs in [&two, &one] {
        let r = align(seqs, &o).expect("aligned");
        assert!(r.refine_skipped);
        let said = describe(&o, &r);
        assert!(said.contains("refinement left out"), "{said}");
        assert!(!said.contains("refined, up to"), "{said}");
    }
    // no more passes are made, or said, than the program allows
    let (seqs, _) = simulate(Kind::Dna, 4, 60, 0.02, 3);
    let many = Options {
        refine: 1000,
        ..Options::default()
    };
    let r = align(&seqs, &many).expect("aligned");
    let said = describe(&many, &r);
    assert!(
        said.contains(&format!("up to {MOST_REFINE} passes")),
        "{said}"
    );
    // asked for none, it says nothing of it
    let none = Options::default();
    assert!(!describe(&none, &align(&seqs, &none).expect("aligned")).contains("refine"));
}

#[test]
fn many_long_sequences_have_their_distances_from_words_and_a_few_short_ones_from_alignments() {
    let sc = Scores::new(Kind::Dna);
    let long: Vec<Vec<u8>> = (0..8u8)
        .map(|k| {
            (0..5000)
                .map(|i| b"ACGT"[(i * 7 + usize::from(k) * 3 + i / 5) % 4])
                .collect()
        })
        .collect();
    let (_, by_words) = distances(&long, &sc, false).expect("distances");
    assert!(by_words, "28 pairs of 5,000 are over the budget");
    let short: Vec<Vec<u8>> = long.iter().map(|s| s[..100].to_vec()).collect();
    let (_, by_words) = distances(&short, &sc, false).expect("distances");
    assert!(!by_words);
}

#[test]
fn lower_case_rna_ambiguity_codes_stops_and_one_letter_sequences_are_aligned_without_harm() {
    let cases: Vec<Vec<(String, Vec<u8>)>> = vec![
        named(&[
            ("a", b"acgtNNacgtacgt"),
            ("b", b"ACGTACACGTACGT"),
            ("c", b"ACGTNNACGTACGA"),
        ]),
        named(&[
            ("a", b"ACGUACGUACGU"),
            ("b", b"ACGUACGAACGU"),
            ("c", b"ACGUACGUACGA"),
        ]),
        named(&[
            ("p", b"MKVLAAGIX*GLL"),
            ("q", b"MKVLAAGIVGGLL"),
            ("r", b"MKVLAAGIVGGML"),
        ]),
        named(&[("a", b"A"), ("b", b"C")]),
        named(&[("a", b"A"), ("b", b"C"), ("c", b"G")]),
    ];
    for seqs in &cases {
        for guide in [Guide::NeighborJoining, Guide::Upgma] {
            let o = Options {
                guide,
                refine: 2,
                ..Options::default()
            };
            let r = align(seqs, &o).expect("aligned");
            check_alignment(seqs, &r.msa);
        }
    }
    // RNA gives the alignment of the same sequences as DNA
    let dna: Vec<(String, Vec<u8>)> = cases[1]
        .iter()
        .map(|(n, s)| {
            (
                n.clone(),
                s.iter()
                    .map(|&c| if c == b'U' { b'T' } else { c })
                    .collect(),
            )
        })
        .collect();
    let (r, d) = (
        align(&cases[1], &Options::default()).expect("aligned"),
        align(&dna, &Options::default()).expect("aligned"),
    );
    let as_dna: Vec<Vec<u8>> = r
        .msa
        .rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|&c| if c == b'U' { b'T' } else { c })
                .collect()
        })
        .collect();
    assert_eq!(as_dna, d.msa.rows);
}

#[test]
fn two_groups_that_are_too_long_to_align_in_one_piece_are_refused_with_the_reason() {
    let sc = Scores::new(Kind::Dna);
    let long = |n: usize| Group::single(0, vec![b'A'; n]);
    let why = merge(&long(7_000), &long(7_000), &[1.0], &sc, true)
        .err()
        .expect("too large");
    assert!(why.contains("too large"), "{why}");
}

/// How good the alignments of simulated sequences are with other gap costs, to choose the costs by: `cargo test tuning -- --ignored
/// --nocapture`. The seeds here are not the ones the tests above use.
#[test]
#[ignore]
fn tuning_the_gap_costs_on_simulated_sequences() {
    let seeds: Vec<u32> = (0..8).map(|k| 1000 + 37 * k).collect();
    let rows = |kind: Kind, rate: f64, open: f32, ext: f32, mismatch: f32| -> (f64, f64) {
        let (mut found, mut true_share) = (0.0, 0.0);
        for &seed in &seeds {
            let (seqs, truth) = simulate(kind, 8, 110, rate, seed);
            let r = align_scored(&seqs, &Options::default(), |k| {
                Scores::with(k, open, ext, mismatch)
            })
            .expect("aligned");
            let (f, t) = accuracy(&r.msa.rows, &truth);
            found += f / seeds.len() as f64;
            true_share += t / seeds.len() as f64;
        }
        (found, true_share)
    };
    let mut centre = Vec::new();
    for rate in [0.02, 0.04, 0.07] {
        let mut sum = 0.0;
        for &seed in &seeds {
            let (seqs, truth) = simulate(Kind::Dna, 8, 110, rate, seed);
            sum +=
                accuracy(&center_star(&seqs).expect("aligned").rows, &truth).0 / seeds.len() as f64;
        }
        centre.push(sum);
    }
    println!("DNA centre-star at 0.02 0.04 0.07: {centre:.3?}");
    for mismatch in [-4.0f32, -1.0, 0.0] {
        for open in [10.0f32, 15.0, 20.0, 30.0] {
            for ext in [0.5f32, 2.0, 4.0, 6.66] {
                let r: Vec<String> = [0.02, 0.04, 0.07]
                    .iter()
                    .map(|&rate| {
                        let (f, t) = rows(Kind::Dna, rate, open, ext, mismatch);
                        format!("{f:.3}/{t:.3}")
                    })
                    .collect();
                println!(
                    "DNA mismatch {mismatch:>4} open {open:>4} extend {ext:>4}: {}",
                    r.join("  ")
                );
            }
        }
    }
    for open in [6.0f32, 10.0, 14.0, 20.0] {
        for ext in [0.2f32, 0.5, 1.0, 2.0] {
            let r: Vec<String> = [0.04, 0.06, 0.10]
                .iter()
                .map(|&rate| {
                    let (f, t) = rows(Kind::Protein, rate, open, ext, 0.0);
                    format!("{f:.3}/{t:.3}")
                })
                .collect();
            println!("protein open {open:>4} extend {ext:>4}: {}", r.join("  "));
        }
    }
}

/// The three built-in methods on simulated sequences, with the costs that were chosen: `cargo test comparison -- --ignored --nocapture`.
#[test]
#[ignore]
fn comparison_of_the_methods_on_simulated_sequences() {
    let seeds: Vec<u32> = (0..8).map(|k| 1000 + 37 * k).collect();
    for (kind, rates) in [
        (Kind::Dna, [0.02, 0.04, 0.07]),
        (Kind::Protein, [0.04, 0.06, 0.10]),
    ] {
        for rate in rates {
            let (mut plain, mut refined, mut centre, mut upgma) = (0.0, 0.0, 0.0, 0.0);
            for &seed in &seeds {
                let (seqs, truth) = simulate(kind, 8, 110, rate, seed);
                let n = seeds.len() as f64;
                let q = |r: &Report| accuracy(&r.msa.rows, &truth).0 / n;
                plain += q(&align(&seqs, &Options::default()).expect("aligned"));
                refined += q(&align(
                    &seqs,
                    &Options {
                        refine: 3,
                        ..Options::default()
                    },
                )
                .expect("aligned"));
                upgma += q(&align(
                    &seqs,
                    &Options {
                        guide: Guide::Upgma,
                        ..Options::default()
                    },
                )
                .expect("aligned"));
                centre += accuracy(&center_star(&seqs).expect("aligned").rows, &truth).0 / n;
            }
            println!("{kind:?} rate {rate}: progressive {plain:.3}, refined {refined:.3}, UPGMA {upgma:.3}, centre-star {centre:.3}");
        }
    }
}

#[test]
fn a_gap_costs_nothing_at_the_ends_when_they_are_free_and_what_the_columns_say_elsewhere() {
    let dna = Scores::new(Kind::Dna);
    let opening = [1.0f32, 2.0, 0.5];
    let (open, extend) = (dna.open, dna.extend);
    // before the first column and after the last, free or not
    assert_eq!(gap_cost(&dna, &opening, 0, true), (0.0, 0.0));
    assert_eq!(gap_cost(&dna, &opening, 3, true), (0.0, 0.0));
    assert_eq!(gap_cost(&dna, &opening, 0, false), (open, extend));
    assert_eq!(gap_cost(&dna, &opening, 3, false), (open * 0.5, extend));
    // between columns, free ends or not: the opening is scaled by the column the gap lies before
    for free in [true, false] {
        assert_eq!(gap_cost(&dna, &opening, 1, free), (open * 2.0, extend));
        assert_eq!(gap_cost(&dna, &opening, 2, free), (open * 0.5, extend));
    }
}

/// The best score of two sequences by an affine-gap dynamic program written plainly (three tables), with the gaps at the ends of
/// a sequence free or not: the same method as `merge`, but by a separate piece of code that works on the two sequences alone.
fn best_score(a: &[u8], b: &[u8], sc: &Scores, free_ends: bool) -> f64 {
    let (m, n) = (a.len(), b.len());
    let neg = f64::NEG_INFINITY;
    let cost = |at: usize, len: usize, step: bool| -> f64 {
        if free_ends && (at == 0 || at == len) {
            0.0
        } else if step {
            f64::from(sc.open)
        } else {
            f64::from(sc.extend)
        }
    };
    let mut mm = vec![vec![neg; n + 1]; m + 1];
    let (mut xx, mut yy) = (mm.clone(), mm.clone());
    mm[0][0] = 0.0;
    for i in 0..=m {
        for j in 0..=n {
            if i > 0 && j > 0 {
                let before = mm[i - 1][j - 1].max(xx[i - 1][j - 1]).max(yy[i - 1][j - 1]);
                mm[i][j] = before + f64::from(sc.score(sc.place(a[i - 1]), sc.place(b[j - 1])));
            }
            if i > 0 {
                // a residue of a against a gap in b, with j residues of b before it
                xx[i][j] = (mm[i - 1][j] - cost(j, n, true))
                    .max(xx[i - 1][j] - cost(j, n, false))
                    .max(yy[i - 1][j] - cost(j, n, true));
            }
            if j > 0 {
                yy[i][j] = (mm[i][j - 1] - cost(i, m, true))
                    .max(yy[i][j - 1] - cost(i, m, false))
                    .max(xx[i][j - 1] - cost(i, m, true));
            }
        }
    }
    mm[m][n].max(xx[m][n]).max(yy[m][n])
}

#[test]
fn the_alignment_of_two_sequences_is_the_best_that_a_plain_dynamic_program_finds_with_the_ends_free_or_not(
) {
    let mut rng = Rng(11);
    for kind in [Kind::Dna, Kind::Protein] {
        let sc = Scores::new(kind);
        let alphabet = letters(kind);
        for round in 0..60 {
            let a: Vec<u8> = (0..14 + round % 7)
                .map(|_| alphabet[rng.below(alphabet.len())])
                .collect();
            // the other is the first with substitutions, deletions and insertions, and sometimes cut at an end
            let mut b = Vec::new();
            for &c in &a {
                match rng.below(10) {
                    0 => b.push(alphabet[rng.below(alphabet.len())]),
                    1 => {}
                    2 => b.extend([c, alphabet[rng.below(alphabet.len())]]),
                    _ => b.push(c),
                }
            }
            if round % 3 == 0 {
                b.drain(..3);
            }
            if round % 5 == 0 {
                b.truncate(b.len() - 2);
            }
            for free_ends in [true, false] {
                let both = merge(
                    &Group::single(0, a.clone()),
                    &Group::single(1, b.clone()),
                    &[1.0, 1.0],
                    &sc,
                    free_ends,
                )
                .expect("merged");
                assert_eq!(strip(&both.rows[0]), a);
                assert_eq!(strip(&both.rows[1]), b);
                let got = sum_of_pairs(&both.rows, &[0, 1], &[1.0, 1.0], &sc, free_ends);
                let want = best_score(&a, &b, &sc, free_ends);
                assert!(
                    (got - want).abs() < 1e-3,
                    "{kind:?} round {round} free {free_ends}: {got} against {want}\n{}\n{}",
                    String::from_utf8_lossy(&both.rows[0]),
                    String::from_utf8_lossy(&both.rows[1])
                );
            }
        }
    }
}

#[test]
fn a_group_of_copies_is_aligned_as_the_one_it_is_a_copy_of() {
    // two identical rows against a row, or one: the scores are averages, so the alignment is the same
    let mut rng = Rng(23);
    let sc = Scores::new(Kind::Dna);
    let mut differing = 0;
    for _ in 0..120 {
        let a: Vec<u8> = (0..12).map(|_| b"ACGT"[rng.below(4)]).collect();
        let mut b = a.clone();
        b.remove(rng.below(b.len()));
        b[rng.below(11)] = b"ACGT"[rng.below(4)];
        let one = Group::single(1, b.clone());
        let two = Group {
            ids: vec![1, 2],
            rows: vec![b.clone(), b.clone()],
        };
        let p = Group::single(0, a);
        let w = [1.0f32; 3];
        let alone = merge(&p, &one, &w, &sc, true).expect("merged");
        let copies = merge(&p, &two, &w, &sc, true).expect("merged");
        differing +=
            usize::from(alone.rows[0] != copies.rows[0] || alone.rows[1] != copies.rows[1]);
    }
    assert_eq!(
        differing, 0,
        "a group of copies was aligned another way {differing} times of 120"
    );
}

#[test]
fn upgma_makes_every_sequence_equally_far_from_the_root_and_neighbour_joining_does_not() {
    // the example distances of the neighbour-joining article: a tree that is not the same distance from its root to every leaf
    let d = vec![
        vec![0.0, 5.0, 9.0, 9.0, 8.0],
        vec![5.0, 0.0, 10.0, 10.0, 9.0],
        vec![9.0, 10.0, 0.0, 8.0, 7.0],
        vec![9.0, 10.0, 8.0, 0.0, 3.0],
        vec![8.0, 9.0, 7.0, 3.0, 0.0],
    ];
    let depths = |tree: &Tree| -> Vec<f64> {
        let mut depth = vec![0.0f64; tree.kids.len()];
        for v in tree.preorder() {
            for &k in &tree.kids[v] {
                depth[k] = depth[v] + tree.len[k];
            }
        }
        (0..tree.kids.len())
            .filter(|&v| tree.leaf[v].is_some())
            .map(|v| depth[v])
            .collect()
    };
    let spread = |d: &[f64]| {
        d.iter().cloned().fold(f64::MIN, f64::max) - d.iter().cloned().fold(f64::MAX, f64::min)
    };
    let upgma = depths(&guide_tree(Guide::Upgma, &d).expect("a tree"));
    let nj = depths(&guide_tree(Guide::NeighborJoining, &d).expect("a tree"));
    assert_eq!(upgma.len(), 5);
    assert!(spread(&upgma) < 1e-9, "{upgma:?}");
    assert!(spread(&nj) > 0.5, "{nj:?}");
}

#[test]
fn a_dot_in_a_sequence_is_a_gap_sign_too_and_too_many_sequences_for_refining_are_not_refined_and_say_so(
) {
    let seqs = named(&[("a", b"AC.GTAC-G"), ("b", b"ACGTACG")]);
    let r = align(&seqs, &Options::default()).expect("aligned");
    assert_eq!(strip(&r.msa.rows[0]), b"ACGTACG");
    assert!(r.msa.rows.iter().all(|row| !row.contains(&b'.')));
    // sixty-one sequences are one more than refinement goes to
    let mut rng = Rng(7);
    let base: Vec<u8> = (0..24).map(|_| b"ACGT"[rng.below(4)]).collect();
    let many: Vec<(String, Vec<u8>)> = (0..=REFINE_MOST)
        .map(|k| {
            let mut s = base.clone();
            s[k % 24] = b"ACGT"[rng.below(4)];
            (format!("s{k}"), s)
        })
        .collect();
    let o = Options {
        refine: 1,
        ..Options::default()
    };
    let r = align(&many, &o).expect("aligned");
    assert!(r.refine_skipped);
    assert!(describe(&o, &r).contains("refinement left out"));
    // and sixty are refined
    let r = align(&many[..REFINE_MOST], &o).expect("aligned");
    assert!(!r.refine_skipped);
}

#[test]
fn residues_of_the_shorter_sequence_that_have_no_partner_count_as_differences_in_the_distance() {
    let mut rng = Rng(31);
    let dna = Scores::new(Kind::Dna);
    let a: Vec<u8> = (0..40).map(|_| b"ACGT"[rng.below(4)]).collect();
    // the same stretch of it, with six bases put in that it has not: they have no place in the alignment, and the rest is as alike
    // as can be
    let plain = a[..28].to_vec();
    let mut with_insertion = a[..15].to_vec();
    with_insertion.extend(b"GGGGGG");
    with_insertion.extend(&a[15..28]);
    let alike = pair_distance(&a, &plain, &dna).expect("a distance");
    let apart = pair_distance(&a, &with_insertion, &dna).expect("a distance");
    assert_eq!(alike, 0.0);
    assert!(apart > 0.1, "{apart}");
}

#[test]
fn refinement_never_returns_a_worse_alignment_whatever_the_sequences() {
    for seed in 1..=150u32 {
        let kind = if seed % 2 == 0 {
            Kind::Dna
        } else {
            Kind::Protein
        };
        let rate = if kind == Kind::Dna { 0.10 } else { 0.12 };
        let (seqs, _) = simulate(kind, 5, 40, rate, seed);
        let r = align(
            &seqs,
            &Options {
                refine: 3,
                ..Options::default()
            },
        )
        .expect("aligned");
        let (after, before) = (
            r.score.expect("a score"),
            r.unrefined_score.expect("a score"),
        );
        assert!(
            after >= before - 1e-9,
            "{kind:?} seed {seed}: {before} then {after}"
        );
    }
}

#[test]
fn more_passes_of_refinement_find_better_alignments_that_one_pass_leaves() {
    let (mut better, mut tried) = (0, 0);
    for seed in 1..=60u32 {
        let kind = if seed % 2 == 0 {
            Kind::Dna
        } else {
            Kind::Protein
        };
        let rate = if kind == Kind::Dna { 0.07 } else { 0.10 };
        let (seqs, _) = simulate(kind, 8, 60, rate, seed);
        let score = |passes| {
            align(
                &seqs,
                &Options {
                    refine: passes,
                    ..Options::default()
                },
            )
            .expect("aligned")
            .score
            .expect("a score")
        };
        let (one, many) = (score(1), score(MOST_REFINE));
        assert!(
            many >= one - 1e-9,
            "{kind:?} seed {seed}: {one} then {many}"
        );
        tried += 1;
        better += usize::from(many > one + 1e-6);
    }
    assert!(better > 0, "none of {tried} gained from more than one pass");
}
