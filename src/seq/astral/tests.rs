//! The counting is checked against counting every quartet one by one, on random trees with polytomies and missing taxa.

use super::bits::Bits;
use super::genes::{read, Genes};
use super::species::Sp;
use super::{infer, Options};
use crate::seq::tree::Rng;

#[derive(Clone)]
enum Shape {
    Leaf(String),
    Join(Vec<Shape>),
}

impl Shape {
    fn text(&self) -> String {
        match self {
            Shape::Leaf(s) => s.clone(),
            Shape::Join(kids) => format!(
                "({})",
                kids.iter().map(Shape::text).collect::<Vec<_>>().join(",")
            ),
        }
    }
}

/// A random binary tree over `leaves`, with each inner branch collapsed with probability `collapse`.
fn random_shape(rng: &mut Rng, leaves: &[String], collapse: f64) -> Shape {
    let mut items: Vec<Shape> = leaves.iter().map(|s| Shape::Leaf(s.clone())).collect();
    while items.len() > 1 {
        let a = items.swap_remove(rng.below(items.len()));
        let b = items.swap_remove(rng.below(items.len()));
        items.push(Shape::Join(vec![a, b]));
    }
    let root = items.pop().expect("at least one leaf");
    flatten(root, rng, collapse)
}

fn flatten(shape: Shape, rng: &mut Rng, collapse: f64) -> Shape {
    let Shape::Join(kids) = shape else {
        return shape;
    };
    let mut out = Vec::new();
    for kid in kids {
        match flatten(kid, rng, collapse) {
            Shape::Join(inner) if rng.unit() < collapse => out.extend(inner),
            other => out.push(other),
        }
    }
    Shape::Join(out)
}

fn names(n: usize) -> Vec<String> {
    (0..n).map(|i| format!("t{i:02}")).collect()
}

/// Gene trees over `n` taxa: each keeps every taxon with probability `1 - missing` (at least four stay).
fn random_genes(rng: &mut Rng, n: usize, count: usize, collapse: f64, missing: f64) -> String {
    let all = names(n);
    let mut text = String::new();
    for _ in 0..count {
        let mut kept: Vec<String> = all
            .iter()
            .filter(|_| rng.unit() >= missing)
            .cloned()
            .collect();
        while kept.len() < 4 {
            let extra = all[rng.below(n)].clone();
            if !kept.contains(&extra) {
                kept.push(extra);
            }
        }
        text.push_str(&random_shape(rng, &kept, collapse).text());
        text.push_str(";\n");
    }
    text
}

fn random_species(rng: &mut Rng, n: usize) -> String {
    format!("{};", random_shape(rng, &names(n), 0.0).text())
}

/// The way a tree resolves four taxa: 0 for `q0 q1 | q2 q3`, 1 for `q0 q2 | q1 q3`, 2 for `q0 q3 | q1 q2`; none if it does
/// not. A tree resolves them when some branch has exactly two of them on one side.
fn topology(clusters: &[Bits], q: [usize; 4]) -> Option<u8> {
    clusters.iter().find_map(|x| {
        let inside: Vec<usize> = (0..4).filter(|&i| x.contains(q[i])).collect();
        match inside[..] {
            [0, 1] | [2, 3] => Some(0),
            [0, 2] | [1, 3] => Some(1),
            [0, 3] | [1, 2] => Some(2),
            _ => None,
        }
    })
}

fn quartets(n: usize) -> Vec<[usize; 4]> {
    let mut out = Vec::new();
    for a in 0..n {
        for b in a + 1..n {
            for c in b + 1..n {
                for d in c + 1..n {
                    out.push([a, b, c, d]);
                }
            }
        }
    }
    out
}

fn species_clusters(sp: &Sp) -> Vec<Bits> {
    sp.nodes.iter().map(|n| n.set.clone()).collect()
}

/// Every quartet counted one by one: the quartets that gene tree `g` and the species tree resolve the same way.
fn brute_score(genes: &Genes, sp: &Sp) -> u64 {
    let species = species_clusters(sp);
    let mut score = 0;
    for tree in &genes.trees {
        let clusters = tree.clusters(genes.taxa.len());
        let present = clusters[tree.root].clone();
        for q in quartets(genes.taxa.len()) {
            if q.iter().all(|&t| present.contains(t)) {
                let (g, s) = (topology(&clusters, q), topology(&species, q));
                score += u64::from(g.is_some() && g == s);
            }
        }
    }
    score
}

/// For a branch with the groups `l r | s o`, the share of each gene's quartets (one taxon from each group) that each tree
/// resolves as `lr|so`, `rs|lo` and `ro|ls`, added over the genes.
fn brute_around(genes: &Genes, groups: [&Bits; 4]) -> ([f64; 3], f64) {
    let mut f = [0.0; 3];
    for tree in &genes.trees {
        let clusters = tree.clusters(genes.taxa.len());
        let present = &clusters[tree.root];
        let members: Vec<Vec<usize>> = groups
            .iter()
            .map(|g| g.iter().filter(|&t| present.contains(t)).collect())
            .collect();
        let total: usize = members.iter().map(Vec::len).product();
        if total == 0 {
            continue;
        }
        let mut counts = [0usize; 3];
        for &l in &members[0] {
            for &r in &members[1] {
                for &s in &members[2] {
                    for &o in &members[3] {
                        // [l, r, s, o] as q0 q1 q2 q3: pairing 0 is lr|so, 1 is ls|ro, 2 is lo|rs
                        match topology(&clusters, [l, r, s, o]) {
                            Some(0) => counts[0] += 1,
                            Some(1) => counts[2] += 1, // ls|ro is the second alternative
                            Some(2) => counts[1] += 1, // lo|rs is the first alternative
                            _ => {}
                        }
                    }
                }
            }
        }
        for (x, c) in f.iter_mut().zip(counts) {
            *x += c as f64 / total as f64;
        }
    }
    let en = f.iter().sum();
    (f, en)
}

#[test]
fn the_score_equals_counting_every_quartet_on_random_trees() {
    let mut rng = Rng::new(7);
    for n in [4usize, 5, 6, 8, 9] {
        for (collapse, missing) in [(0.0, 0.0), (0.3, 0.0), (0.0, 0.3), (0.3, 0.3)] {
            let genes = read(&random_genes(&mut rng, n, 12, collapse, missing), 0).expect("genes");
            let n = genes.taxa.len();
            for _ in 0..4 {
                let sp = Sp::from_newick(&random_species(&mut rng, n), &genes.taxa)
                    .expect("species tree");
                assert_eq!(
                    sp.score(&genes),
                    brute_score(&genes, &sp),
                    "n={n} collapse={collapse} missing={missing}"
                );
            }
        }
    }
}

#[test]
fn resolved_quartets_equal_counting_the_resolved_ones() {
    let mut rng = Rng::new(11);
    for (collapse, missing) in [(0.0, 0.0), (0.4, 0.0), (0.4, 0.3)] {
        let genes = read(&random_genes(&mut rng, 9, 15, collapse, missing), 0).expect("genes");
        for (tree, gene) in genes.trees.iter().zip(&genes.genes) {
            let clusters = tree.clusters(genes.taxa.len());
            let present = &clusters[tree.root];
            let by_hand = quartets(genes.taxa.len())
                .into_iter()
                .filter(|q| q.iter().all(|&t| present.contains(t)))
                .filter(|&q| topology(&clusters, q).is_some())
                .count() as u128;
            assert_eq!(gene.resolved_quartets(), by_hand);
        }
    }
}

#[test]
fn the_frequencies_around_every_branch_equal_counting_every_quartet() {
    let mut rng = Rng::new(23);
    for (collapse, missing) in [(0.0, 0.0), (0.3, 0.0), (0.0, 0.3), (0.3, 0.3)] {
        let genes = read(&random_genes(&mut rng, 8, 14, collapse, missing), 0).expect("genes");
        let n = genes.taxa.len();
        let sp = Sp::from_newick(&random_species(&mut rng, n), &genes.taxa).expect("species tree");
        let all = Bits::of(n, 0..n);
        let branches = sp.annotate(&genes, 0.5);
        assert_eq!(
            branches.len(),
            n - 3,
            "an unrooted binary tree has n - 3 inner branches"
        );
        for b in &branches {
            let node = &sp.nodes[b.node];
            let (l, r) = node.kids.expect("an inner branch");
            let parent = &sp.nodes[node.parent.expect("a parent")];
            let sib = match parent.kids.expect("two children") {
                (a, _) if a != b.node => a,
                (_, c) => c,
            };
            let o = all.minus(&node.set.union(&sp.nodes[sib].set));
            let groups = [&sp.nodes[l].set, &sp.nodes[r].set, &sp.nodes[sib].set, &o];
            let (f, en) = brute_around(&genes, groups);
            for (i, (mine, counted)) in b.f.iter().zip(f).enumerate() {
                assert!(
                    (mine - counted).abs() < 1e-9,
                    "f{i}: {mine} against {counted}"
                );
            }
            assert!((b.en - en).abs() < 1e-9);
            assert_eq!(
                b.quartets,
                groups.iter().map(|g| u64::from(g.count())).product::<u64>()
            );
            assert!((b.pp.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        }
    }
}

#[test]
fn a_species_tree_is_read_whatever_its_root_and_refused_when_it_is_not_binary() {
    let genes = read("((A,B),(C,(D,E)));", 0).expect("genes");
    let rooted = Sp::from_newick("(((A,B),C),(D,E));", &genes.taxa).expect("rooted");
    let unrooted = Sp::from_newick("((A,B),C,(D,E));", &genes.taxa).expect("unrooted");
    assert_eq!(rooted.score(&genes), unrooted.score(&genes));
    assert_eq!(rooted.score(&genes), 5);
    let star = Sp::from_newick("(A,B,C,D,E);", &genes.taxa);
    assert!(star.is_err());
    let polytomy = Sp::from_newick("((A,B,C),(D,E));", &genes.taxa);
    assert!(polytomy.is_err());
    assert!(
        Sp::from_newick("((A,B),(C,D));", &genes.taxa).is_err(),
        "a taxon is missing"
    );
    assert!(
        Sp::from_newick("((A,B),(C,(D,F)));", &genes.taxa).is_err(),
        "a taxon is unknown"
    );
    assert!(
        Sp::from_newick("((A,B),(C,(D,D)));", &genes.taxa).is_err(),
        "a name twice"
    );
}

/// Every rooted binary tree over the leaves, by adding the leaves one at a time to every branch of every tree so far.
fn all_shapes(leaves: &[String]) -> Vec<Shape> {
    let first = [
        Shape::Leaf(leaves[0].clone()),
        Shape::Leaf(leaves[1].clone()),
    ];
    let mut trees = vec![Shape::Join(first.to_vec())];
    for leaf in &leaves[2..] {
        trees = trees.iter().flat_map(|t| attach(t, leaf)).collect();
    }
    trees
}

fn attach(tree: &Shape, leaf: &str) -> Vec<Shape> {
    let mut out = vec![Shape::Join(vec![tree.clone(), Shape::Leaf(leaf.into())])];
    if let Shape::Join(kids) = tree {
        for (i, kid) in kids.iter().enumerate() {
            for new in attach(kid, leaf) {
                let mut changed = kids.clone();
                changed[i] = new;
                out.push(Shape::Join(changed));
            }
        }
    }
    out
}

#[test]
fn the_exact_search_finds_the_best_of_all_trees_counted_quartet_by_quartet() {
    let mut rng = Rng::new(31);
    for (n, collapse, missing) in [
        (5usize, 0.0, 0.0),
        (6, 0.3, 0.2),
        (7, 0.2, 0.0),
        (7, 0.3, 0.2),
    ] {
        let text = random_genes(&mut rng, n, 14, collapse, missing);
        let genes = read(&text, 0).expect("genes");
        let taxa = names(genes.taxa.len());
        let best = all_shapes(&taxa[1..])
            .iter()
            .map(|shape| {
                let tree = format!("({},{});", taxa[0], shape.text());
                brute_score(
                    &genes,
                    &Sp::from_newick(&tree, &genes.taxa).expect("a tree"),
                )
            })
            .max()
            .expect("trees");
        let o = Options {
            exact: true,
            polish: false,
            ..Options::default()
        };
        let found = infer(&text, &o).expect("a tree");
        assert_eq!(
            found.score, best,
            "n={n} collapse={collapse} missing={missing}"
        );
        // the tree it returns scores what it says
        let again = Sp::from_newick(&found.newick, &genes.taxa).expect("its own tree");
        assert_eq!(brute_score(&genes, &again), found.score);
    }
}

#[test]
fn polishing_never_lowers_the_score_and_the_root_taxon_only_changes_the_writing() {
    let mut rng = Rng::new(41);
    for (collapse, missing) in [(0.1, 0.0), (0.3, 0.0), (0.2, 0.3)] {
        let text = random_genes(&mut rng, 11, 25, collapse, missing);
        let plain = Options {
            polish: false,
            ..Options::default()
        };
        let without = infer(&text, &plain).expect("a tree");
        let with = infer(&text, &Options::default()).expect("a tree");
        assert!(
            with.score >= without.score,
            "{} against {}",
            with.score,
            without.score
        );
        let hung = Options {
            root: Some("t05".into()),
            polish: false,
            ..Options::default()
        };
        let other = infer(&text, &hung).expect("a tree");
        assert!(other.newick.starts_with("(t05,"), "{}", other.newick);
        if missing == 0.0 {
            // with every taxon in every tree the clusters do not depend on where the tree hangs
            assert_eq!(other.score, without.score, "collapse {collapse}");
        }
        assert!(infer(
            &text,
            &Options {
                root: Some("nobody".into()),
                ..Options::default()
            }
        )
        .is_err());
    }
}

#[test]
fn the_exact_search_refuses_more_taxa_than_it_can_take_and_a_tree_needs_four() {
    let mut rng = Rng::new(5);
    let text = random_genes(&mut rng, 14, 4, 0.0, 0.0);
    let o = Options {
        exact: true,
        ..Options::default()
    };
    assert!(infer(&text, &o).expect_err("too many").contains("at most"));
    assert!(infer("(A,B,C);", &Options::default()).is_err());
}

#[test]
fn the_same_gene_trees_give_the_same_tree_every_time_down_to_the_order_of_the_branches() {
    // the guide tree's splits come out of a hash set, whose order differs from one run to the next: the order the
    // clusters are added in must not decide which of several equal trees is built, nor how its branches are listed
    let text = crate::sample::SPECIES_GENE_TREES;
    for options in [
        Options::default(),
        Options {
            extra: 0,
            ..Options::default()
        },
        Options {
            polish: false,
            ..Options::default()
        },
    ] {
        let first = infer(text, &options).expect("a tree");
        for run in 1..24 {
            let again = infer(text, &options).expect("a tree");
            assert_eq!(again.newick, first.newick, "run {run}");
            assert_eq!(again.branches, first.branches, "run {run}");
            assert_eq!(again.clusters, first.clusters, "run {run}");
        }
    }
}
