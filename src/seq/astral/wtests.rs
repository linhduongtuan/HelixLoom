//! Weighted quartet counting against the definition: the weight of every quartet found by walking the tree, one quartet at a
//! time, for small random trees with polytomies, missing species, supports and lengths.

use super::bits::Bits;
use super::genes::{read, read_weighted, Genes};
use super::weight::tripartition;
use super::wgenes::{WTree, Weights};
use super::wscore;

// --- random trees ---------------------------------------------------------------------------------------------------------------

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn unit(&mut self) -> f64 {
        (self.next() % 10_000) as f64 / 10_000.0
    }
}

/// A random unrooted binary tree on the leaves `names`, as neighbour lists: the leaves are 0.., then the inner nodes.
fn random_shape(names: usize, rng: &mut Rng) -> Vec<Vec<usize>> {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); names];
    // an inner node joins the first three leaves
    let first = adj.len();
    adj.push(vec![0, 1, 2]);
    for list in adj.iter_mut().take(3) {
        list.push(first);
    }
    let mut edges = vec![(0, first), (1, first), (2, first)];
    for leaf in 3..names {
        // the new leaf hangs from a new node on a random branch
        let k = rng.below(edges.len());
        let (u, v) = edges[k];
        let mid = adj.len();
        adj.push(vec![u, v, leaf]);
        for (a, b) in [(u, v), (v, u)] {
            let slot = adj[a].iter_mut().find(|x| **x == b).expect("neighbour");
            *slot = mid;
        }
        adj[leaf].push(mid);
        edges[k] = (u, mid);
        edges.push((mid, v));
        edges.push((mid, leaf));
    }
    adj
}

/// The Newick text of a random tree on the species `present` (indices into `letters`), hung from an inner node: with a
/// support (two decimals) and a length on each inner branch, a length on each leaf branch, and some inner branches taken out
/// (polytomies). `rooted` hangs it from a binary root instead.
struct Style {
    supports: bool,
    lengths: bool,
    collapse: bool,
    /// Every inner branch has the support 1 (instead of a random one).
    ones: bool,
}

fn newick(letters: &[&str], present: &[usize], style: &Style, rng: &mut Rng) -> String {
    let adj = random_shape(present.len(), rng);
    let first = present.len();
    fn print(
        v: usize,
        from: usize,
        adj: &[Vec<usize>],
        letters: &[&str],
        present: &[usize],
        style: &Style,
        rng: &mut Rng,
    ) -> String {
        let length = |rng: &mut Rng| {
            if style.lengths {
                format!(":{:.3}", 0.01 + rng.unit() * 0.4)
            } else {
                String::new()
            }
        };
        if v < present.len() {
            return format!("{}{}", letters[present[v]], length(rng));
        }
        let mut kids = Vec::new();
        for &c in adj[v].iter().filter(|&&c| c != from) {
            let inner = c >= present.len();
            if inner && style.collapse && rng.below(4) == 0 {
                // the branch is taken out: what hangs from it hangs from this node
                let text = print(c, v, adj, letters, present, style, rng);
                kids.push(text[1..text.rfind(')').expect("a close")].to_string());
            } else {
                kids.push(print(c, v, adj, letters, present, style, rng));
            }
        }
        let label = if style.ones {
            "1".to_string()
        } else if style.supports {
            format!("{:.2}", rng.unit())
        } else {
            String::new()
        };
        format!("({}){label}{}", kids.join(","), length(rng))
    }
    let body = print(first, usize::MAX, &adj, letters, present, style, rng);
    // the root has no branch above it: its own label and length are cut off
    let body = &body[..=body.rfind(')').expect("a close")];
    format!("{body};")
}

const LETTERS: [&str; 9] = ["A", "B", "C", "D", "E", "F", "G", "H", "I"];

/// Some gene trees on up to nine species, some of them missing some of the species.
fn random_genes(count: usize, style: &Style, seed: u64) -> String {
    let mut rng = Rng(seed);
    (0..count)
        .map(|_| {
            let size = 4 + rng.below(LETTERS.len() - 3);
            let mut present: Vec<usize> = (0..LETTERS.len()).collect();
            for k in (1..present.len()).rev() {
                present.swap(k, rng.below(k + 1));
            }
            present.truncate(size);
            newick(&LETTERS, &present, style, &mut rng)
        })
        .collect::<Vec<_>>()
        .join("\n")
}

// --- the definition ---------------------------------------------------------------------------------------------------------

/// The branches on the path between two nodes of a tree, each named by its lower node.
fn path(t: &WTree, a: usize, b: usize) -> Vec<usize> {
    let depth = |mut v: usize| {
        let mut d = 0;
        while let Some(p) = t.nodes[v].parent {
            v = p;
            d += 1;
        }
        d
    };
    let (mut x, mut y) = (a, b);
    let (mut dx, mut dy) = (depth(x), depth(y));
    let mut edges = Vec::new();
    while dx > dy {
        edges.push(x);
        x = t.nodes[x].parent.expect("a parent");
        dx -= 1;
    }
    while dy > dx {
        edges.push(y);
        y = t.nodes[y].parent.expect("a parent");
        dy -= 1;
    }
    while x != y {
        edges.push(x);
        edges.push(y);
        x = t.nodes[x].parent.expect("a parent");
        y = t.nodes[y].parent.expect("a parent");
    }
    edges
}

/// The weighted quartets of one tree for the clusters `a`, `b` and the rest, by looking at every four species.
fn reference_tree(t: &WTree, a: &Bits, b: &Bits, w: Weights) -> f64 {
    let leaves: Vec<(usize, usize)> = t
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(i, n)| n.taxon.map(|x| (i, x)))
        .collect();
    let class =
        |x: usize| usize::from(!a.contains(x)) + usize::from(!a.contains(x) && !b.contains(x));
    let mut total = 0.0;
    let n = leaves.len();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                for l in k + 1..n {
                    let four = [leaves[i], leaves[j], leaves[k], leaves[l]];
                    let mut counts = [0; 3];
                    for (_, x) in four {
                        counts[class(x)] += 1;
                    }
                    // two species of a cluster, one of each of the others
                    if counts.iter().max() != Some(&2) || counts.iter().min() != Some(&1) {
                        continue;
                    }
                    total += quartet(t, four.map(|(v, _)| v), four.map(|(_, x)| class(x)), w);
                }
            }
        }
    }
    total
}

/// The weight of the quartet of the nodes `v` (of the clusters `c`) if the tree resolves it with the two of one cluster
/// together; 0 otherwise.
fn quartet(t: &WTree, v: [usize; 4], c: [usize; 3 + 1], w: Weights) -> f64 {
    let paired = |p: usize, q: usize, r: usize, s: usize| {
        path(t, v[p], v[q]).len() + path(t, v[r], v[s]).len()
    };
    // the three pairings, by the edges their paths hold: the one the tree has is shorter than the other two
    let cost = [paired(0, 1, 2, 3), paired(0, 2, 1, 3), paired(0, 3, 1, 2)];
    let pairings = [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)];
    let Some(best) = (0..3).find(|&p| (0..3).all(|o| o == p || cost[p] < cost[o])) else {
        return 0.0;
    };
    let (p, q, r, s) = pairings[best];
    // the pair that is together must be the two of one cluster
    if c[p] != c[q] && c[r] != c[s] {
        return 0.0;
    }
    let (ab, cd) = (path(t, v[p], v[q]), path(t, v[r], v[s]));
    let between: Vec<usize> = path(t, v[p], v[r])
        .into_iter()
        .filter(|e| !ab.contains(e) && !cd.contains(e))
        .collect();
    let mut weight = 1.0;
    if w.uses_support() {
        weight *= 1.0
            - between
                .iter()
                .map(|&e| t.nodes[e].edge.miss)
                .product::<f64>();
    }
    if w.uses_length() {
        weight *= ab
            .iter()
            .chain(&cd)
            .map(|&e| t.nodes[e].edge.decay)
            .product::<f64>();
    }
    weight
}

fn reference(genes: &Genes, a: &Bits, b: &Bits, w: Weights) -> f64 {
    genes
        .wtrees
        .iter()
        .map(|t| reference_tree(t, a, b, w))
        .sum()
}

/// Some pairs of disjoint groups of the nine species.
fn groups(n: usize, rng: &mut Rng, count: usize) -> Vec<(Bits, Bits)> {
    (0..count)
        .map(|_| {
            let (mut a, mut b) = (Bits::new(n), Bits::new(n));
            for x in 0..n {
                match rng.below(3) {
                    0 => a.insert(x),
                    1 => b.insert(x),
                    _ => {}
                }
            }
            (a, b)
        })
        .collect()
}

fn close(x: f64, y: f64) -> bool {
    (x - y).abs() <= 1e-9 * (1.0 + x.abs().max(y.abs()))
}

// --- tests ----------------------------------------------------------------------------------------------------------------------

#[test]
fn with_every_quartet_counting_one_the_weighted_counts_are_the_plain_ones() {
    let mut counted = 0;
    // supports of 1 and no lengths: every resolved quartet has the weight 1
    for seed in 1..=12 {
        let text = random_genes(
            5,
            &Style {
                supports: false,
                lengths: false,
                collapse: true,
                ones: true,
            },
            seed,
        );
        let plain = read(&text, 0).expect("trees");
        let weighted = read_weighted(&text, 0, Weights::Support).expect("trees");
        let mut rng = Rng(seed * 7);
        for (a, b) in groups(plain.taxa.len(), &mut rng, 8) {
            let expected = tripartition(&plain.genes, &a, &b) as f64;
            let found = wscore::tripartition(&weighted.wtrees, &a, &b, Weights::Support);
            assert!(
                close(found, expected),
                "seed {seed}: {found} against {expected}"
            );
            counted += usize::from(expected > 0.0);
        }
    }
    // the sums are not all zero
    assert!(counted > 20, "{counted}");
}

#[test]
fn the_fast_weighted_counts_are_the_weights_of_the_quartets_one_by_one() {
    let mut counted = 0;
    for (w, supports, lengths) in [
        (Weights::Support, true, false),
        (Weights::Length, false, true),
        (Weights::Hybrid, true, true),
    ] {
        for seed in 1..=10 {
            let text = random_genes(
                4,
                &Style {
                    supports,
                    lengths,
                    collapse: seed % 2 == 0,
                    ones: false,
                },
                seed * 31,
            );
            let genes = read_weighted(&text, 0, w).expect("trees");
            let mut rng = Rng(seed);
            for (a, b) in groups(genes.taxa.len(), &mut rng, 6) {
                let (found, expected) = (
                    wscore::tripartition(&genes.wtrees, &a, &b, w),
                    reference(&genes, &a, &b, w),
                );
                assert!(
                    close(found, expected),
                    "{w:?}, seed {seed}: {found} against {expected}\n{text}"
                );
                counted += usize::from(expected > 0.0);
            }
        }
    }
    // the sums are not all zero: the comparison is of something
    assert!(counted > 60, "{counted}");
}

#[test]
fn a_quartet_has_the_weight_of_the_formula() {
    // A and B meet at the root, C and D at the other node, the branch between them has the support 0.5; the paths of the
    // pairs are 0.1 + 0.2 and 0.3 + 0.4 long
    let text = "(A:0.1,B:0.2,(C:0.3,D:0.4)0.5:0.6);";
    let at = |w| {
        let genes = read_weighted(text, 0, w).expect("a tree");
        let (a, b) = (Bits::of(4, [0, 1]), Bits::of(4, [2]));
        wscore::tripartition(&genes.wtrees, &a, &b, w)
    };
    assert!(close(at(Weights::Support), 0.5));
    assert!(close(at(Weights::Length), (-1.0f64).exp()));
    assert!(close(at(Weights::Hybrid), 0.5 * (-1.0f64).exp()));
    // a quartet the tree does not have counts nothing
    let genes = read_weighted(text, 0, Weights::Hybrid).expect("a tree");
    let (a, b) = (Bits::of(4, [0, 2]), Bits::of(4, [1]));
    assert_eq!(
        wscore::tripartition(&genes.wtrees, &a, &b, Weights::Hybrid),
        0.0
    );
}

#[test]
fn the_support_of_several_branches_is_one_minus_the_product_of_what_each_misses() {
    // A B | D E: the path between the pairs holds two branches, with the supports 0.5 and 0.8; A B | D C: one, with 0.5
    let text = "(A:1,B:1,(C:1,(D:1,E:1)0.8:1)0.5:1);";
    let genes = read_weighted(text, 0, Weights::Support).expect("a tree");
    let (a, b) = (Bits::of(5, [0, 1]), Bits::of(5, [3]));
    let found = wscore::tripartition(&genes.wtrees, &a, &b, Weights::Support);
    let expected = (1.0 - (1.0 - 0.5) * (1.0 - 0.8)) + 0.5;
    assert!(close(found, expected), "{found} against {expected}");
}

#[test]
fn a_rooted_tree_and_the_same_tree_without_a_root_weigh_the_same() {
    let unrooted = "(A:0.1,B:0.2,(C:0.3,D:0.4,E:0.2)0.8:0.5);";
    // the branch between the two inner nodes is cut in two by the root: its length is the sum and its support is the one
    // that is there
    let rooted = "((A:0.1,B:0.2):0.2,(C:0.3,D:0.4,E:0.2)0.8:0.3);";
    // with a support on both sides of the root, the larger is the support of the branch
    let both = "((A:0.1,B:0.2)0.6:0.2,(C:0.3,D:0.4,E:0.2)0.8:0.3);";
    for w in [Weights::Support, Weights::Length, Weights::Hybrid] {
        let (u, r) = (
            read_weighted(unrooted, 0, w).expect("a tree"),
            read_weighted(rooted, 0, w).expect("a tree"),
        );
        let both = read_weighted(both, 0, w).expect("a tree");
        let mut rng = Rng(5);
        for (a, b) in groups(5, &mut rng, 12) {
            let (x, y) = (
                wscore::tripartition(&u.wtrees, &a, &b, w),
                wscore::tripartition(&r.wtrees, &a, &b, w),
            );
            assert!(close(x, y), "{w:?}: {x} against {y}");
            let z = wscore::tripartition(&both.wtrees, &a, &b, w);
            assert!(close(x, z), "{w:?}, both supported: {x} against {z}");
        }
    }
}

#[test]
fn supports_in_percent_are_read_as_probabilities() {
    let percent = read_weighted(
        "(A:1,B:1,(C:1,D:1)50:1);\n(A:1,B:1,(C:1,D:1)98:1);",
        0,
        Weights::Support,
    )
    .expect("trees");
    let unit = read_weighted(
        "(A:1,B:1,(C:1,D:1)0.5:1);\n(A:1,B:1,(C:1,D:1)0.98:1);",
        0,
        Weights::Support,
    )
    .expect("trees");
    let (a, b) = (Bits::of(4, [0, 1]), Bits::of(4, [2]));
    let (x, y) = (
        wscore::tripartition(&percent.wtrees, &a, &b, Weights::Support),
        wscore::tripartition(&unit.wtrees, &a, &b, Weights::Support),
    );
    assert!(close(x, 0.5 + 0.98) && close(x, y), "{x} {y}");
}

#[test]
fn what_the_weights_cannot_read_is_said_with_the_gene_it_is_in() {
    let err = |text: &str, w| read_weighted(text, 0, w).expect_err("refused");
    // no supports in a tree
    let why = err(
        "(A:1,B:1,(C:1,D:1)0.9:1);\n(A:1,B:1,(C:1,D:1):1);",
        Weights::Support,
    );
    assert!(why.contains("Gene tree 2 has no branch support"), "{why}");
    // no lengths
    let why = err("(A,B,(C,D)0.9);", Weights::Length);
    assert!(why.contains("Gene tree 1 has no branch lengths"), "{why}");
    let why = err("(A:1,B:1,(C:1,D:1)95/100:1);", Weights::Hybrid);
    assert!(why.contains("more than one number"), "{why}");
    let why = err("(A:1,B:1,(C:1,D:1)high:1);", Weights::Support);
    assert!(why.contains("not a support value"), "{why}");
    let why = err("(A:1,B:1,(C:1,D:1)150:1);", Weights::Support);
    assert!(why.contains("probabilities"), "{why}");
    // the same trees are fine for weights that do not read what is missing
    assert!(read_weighted("(A,B,(C,D));", 0, Weights::None).is_ok());
    assert!(read_weighted("(A:1,B:1,(C:1,D:1):1);", 0, Weights::Length).is_ok());
    assert!(read_weighted("(A,B,(C,D)0.5);", 0, Weights::Support).is_ok());
}

#[test]
fn a_branch_without_support_next_to_ones_with_is_taken_as_fully_supported_and_said() {
    let genes =
        read_weighted("(A:1,B:1,(C:1,(D:1,E:1):1)0.5:1);", 0, Weights::Support).expect("a tree");
    assert_eq!(genes.notes.len(), 1, "{:?}", genes.notes);
    assert!(
        genes.notes[0].contains("1 inner branch"),
        "{:?}",
        genes.notes
    );
    // AB | DE passes the branch with 0.5 and the one with no support (taken as 1): 1 - 0.5 * 0 = 1; AB | DC passes only the
    // first: 0.5
    let (a, b) = (Bits::of(5, [0, 1]), Bits::of(5, [3]));
    let found = wscore::tripartition(&genes.wtrees, &a, &b, Weights::Support);
    assert!(close(found, 1.5), "{found}");
}

#[test]
fn weights_that_read_nothing_ask_nothing_of_the_trees() {
    let genes = read_weighted("((A,B),(C,D));", 0, Weights::None).expect("a tree");
    assert!(genes.wtrees.is_empty() && genes.notes.is_empty());
}

// --- the search --------------------------------------------------------------------------------------------------------------

use super::{infer, score_tree, Options};

fn with(weights: Weights) -> Options {
    Options {
        weights,
        ..Options::default()
    }
}

/// The two taxa that share the inner branch of a tree on four species.
fn cherry(text: &str, o: &Options) -> Vec<String> {
    let out = infer(text, o).expect("a species tree");
    assert_eq!(out.branches.len(), 1, "{}", out.newick);
    out.branches[0].taxa.clone()
}

#[test]
fn weights_of_one_give_the_tree_and_the_score_of_astral() {
    for seed in 1..=6 {
        let text = random_genes(
            12,
            &Style {
                supports: false,
                lengths: false,
                collapse: seed % 2 == 0,
                ones: true,
            },
            seed * 101,
        );
        let plain = infer(&text, &Options::default()).expect("a tree");
        let weighted = infer(&text, &with(Weights::Support)).expect("a tree");
        assert_eq!(weighted.newick, plain.newick, "seed {seed}");
        assert_eq!(weighted.score, plain.score);
        assert_eq!(weighted.weighted_score, Some(plain.score as f64));
        assert_eq!(weighted.weights, Weights::Support);
        assert_eq!(plain.weighted_score, None);
        assert!((weighted.normalized - plain.normalized).abs() < 1e-12);
    }
}

#[test]
fn well_supported_gene_trees_outweigh_more_of_the_poorly_supported() {
    // three gene trees that join A with B, each sure of it, against four that join A with C, each hardly so
    let ab = "(A:0.1,B:0.1,(C:0.1,D:0.1)0.95:0.1);\n".repeat(3);
    let ac = "(A:0.1,C:0.1,(B:0.1,D:0.1)0.05:0.1);\n".repeat(4);
    let text = format!("{ab}{ac}");
    // counted one each, the four win
    assert_eq!(cherry(&text, &Options::default()), ["B", "D"]);
    // weighted by their support, the three
    for w in [Weights::Support, Weights::Hybrid] {
        assert_eq!(cherry(&text, &with(w)), ["C", "D"], "{w:?}");
    }
    let out = infer(&text, &with(Weights::Support)).expect("a tree");
    assert_eq!(
        out.score, 3,
        "the three quartets it shares, each counted one"
    );
    assert!(close(out.weighted_score.expect("weighted"), 3.0 * 0.95));
    // the branch support is ASTRAL's, and the log says so
    assert!(
        out.log.iter().any(|l| l.contains("ASTRAL's")),
        "{:?}",
        out.log
    );
}

#[test]
fn trees_with_short_paths_outweigh_more_with_long_ones() {
    let short = "(A:0.01,B:0.01,(C:0.01,D:0.01)0.9:0.1);\n".repeat(3);
    let long = "(A:2,C:2,(B:2,D:2)0.9:0.1);\n".repeat(4);
    let text = format!("{short}{long}");
    assert_eq!(cherry(&text, &Options::default()), ["B", "D"]);
    for w in [Weights::Length, Weights::Hybrid] {
        assert_eq!(cherry(&text, &with(w)), ["C", "D"], "{w:?}");
    }
}

#[test]
fn a_given_tree_is_scored_with_the_weights_too() {
    let ab = "(A:0.1,B:0.1,(C:0.1,D:0.1)0.8:0.1);\n".repeat(2);
    let o = with(Weights::Support);
    let scored = score_tree(&ab, "((A,B),(C,D));", &o).expect("a score");
    assert_eq!(scored.score, 2);
    assert!(close(scored.weighted_score.expect("weighted"), 1.6));
    let other = score_tree(&ab, "((A,C),(B,D));", &o).expect("a score");
    assert_eq!((other.score, other.weighted_score), (0, Some(0.0)));
    assert!(
        scored.normalized > 0.7 && scored.normalized < 0.9,
        "{}",
        scored.normalized
    );
}

#[test]
fn what_the_trees_lack_stops_a_weighted_run_before_it_starts() {
    let why = infer("((A,B),(C,D));", &with(Weights::Support)).expect_err("no supports");
    assert!(why.contains("no branch support"), "{why}");
    let why = infer("((A,B)0.9,(C,D)0.9);", &with(Weights::Hybrid)).expect_err("no lengths");
    assert!(why.contains("no branch lengths"), "{why}");
    assert!(infer("((A,B),(C,D));", &Options::default()).is_ok());
}

#[test]
fn the_names_of_the_weights_are_read_back() {
    for w in Weights::ALL {
        assert_eq!(Weights::parse(w.name()), Some(w));
        assert_eq!(Weights::parse(&w.name().to_uppercase()), Some(w));
    }
    assert_eq!(Weights::parse("hybrid "), Some(Weights::Hybrid));
    assert_eq!(Weights::parse("heavy"), None);
    assert_eq!(Weights::default(), Weights::None);
}

/// How much dearer a weighted count is than a plain one, to set `WEIGHTED_WORDS`:
/// `cargo test --bin helixloom weighted_count_cost -- --ignored --nocapture`.
#[test]
#[ignore = "a measurement, not a check"]
fn weighted_count_cost() {
    for (species, genes) in [(25usize, 200usize), (40, 400), (60, 600)] {
        let mut rng = Rng(9);
        let names: Vec<String> = (0..species).map(|i| format!("t{i}")).collect();
        let letters: Vec<&str> = names.iter().map(String::as_str).collect();
        let present: Vec<usize> = (0..species).collect();
        let style = Style {
            supports: true,
            lengths: true,
            collapse: false,
            ones: false,
        };
        let text: String = (0..genes)
            .map(|_| newick(&letters, &present, &style, &mut rng))
            .collect::<Vec<_>>()
            .join("\n");
        let plain = read(&text, 0).expect("trees");
        let weighted = read_weighted(&text, 0, Weights::Hybrid).expect("trees");
        let pairs = groups(species, &mut rng, 300);
        let time = |f: &dyn Fn(&Bits, &Bits) -> f64| {
            let start = std::time::Instant::now();
            let sum: f64 = pairs.iter().map(|(a, b)| f(a, b)).sum();
            (start.elapsed().as_secs_f64() / pairs.len() as f64, sum)
        };
        let (p, _) = time(&|a, b| tripartition(&plain.genes, a, b) as f64);
        let (w, _) = time(&|a, b| wscore::tripartition(&weighted.wtrees, a, b, Weights::Hybrid));
        println!(
            "{species} species, {genes} genes: plain {:.1} us, weighted {:.1} us, {:.0} times",
            p * 1e6,
            w * 1e6,
            w / p
        );
    }
}
