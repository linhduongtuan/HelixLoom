//! Gene families against the definition: the score of a species tree found by counting the quartet classes of the trees one
//! at a time, for small random trees with several copies of some species.

use super::bits::Bits;
use super::family::{read, Families, Rooted};
use super::pscore;
use super::species::Sp;
use std::collections::BTreeSet;

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
}

/// A random rooted binary tree over `labels` (in this order, cut at random places).
fn random_rooted(labels: &[String], rng: &mut Rng) -> String {
    if labels.len() == 1 {
        return labels[0].clone();
    }
    let cut = 1 + rng.below(labels.len() - 1);
    format!(
        "({},{})",
        random_rooted(&labels[..cut], rng),
        random_rooted(&labels[cut..], rng)
    )
}

const SPECIES: [&str; 7] = ["A", "B", "C", "D", "E", "F", "G"];

/// A family of `leaves` copies of `species` species (each at least once when it can be), in a random shape.
fn random_family(species: usize, leaves: usize, rng: &mut Rng) -> String {
    let mut labels: Vec<String> = (0..species).map(|s| SPECIES[s].to_string()).collect();
    while labels.len() < leaves {
        labels.push(SPECIES[rng.below(species)].to_string());
    }
    for k in (1..labels.len()).rev() {
        labels.swap(k, rng.below(k + 1));
    }
    format!("{};", random_rooted(&labels, rng))
}

fn species_tree(names: &[String], rng: &mut Rng) -> String {
    let mut shuffled = names.to_vec();
    for k in (1..shuffled.len()).rev() {
        shuffled.swap(k, rng.below(k + 1));
    }
    format!("{};", random_rooted(&shuffled, rng))
}

// --- the definition ---------------------------------------------------------------------------------------------------------

struct Tree {
    parent: Vec<Option<usize>>,
    depth: Vec<usize>,
}

impl Tree {
    fn of(r: &Rooted) -> Tree {
        let mut parent = vec![None; r.kids.len()];
        let mut depth = vec![0; r.kids.len()];
        let mut stack = vec![r.root];
        while let Some(v) = stack.pop() {
            if let Some((a, b)) = r.kids[v] {
                for c in [a, b] {
                    parent[c] = Some(v);
                    depth[c] = depth[v] + 1;
                    stack.push(c);
                }
            }
        }
        Tree { parent, depth }
    }

    fn lca(&self, mut a: usize, mut b: usize) -> usize {
        while self.depth[a] > self.depth[b] {
            a = self.parent[a].expect("a parent");
        }
        while self.depth[b] > self.depth[a] {
            b = self.parent[b].expect("a parent");
        }
        while a != b {
            a = self.parent[a].expect("a parent");
            b = self.parent[b].expect("a parent");
        }
        a
    }
}

/// How a species tree resolves four taxa: the pair that is together with the smaller number first, as a mask.
fn pairing(sp: &Sp, quartet: [usize; 4]) -> u8 {
    for node in &sp.nodes {
        let inside: Vec<usize> = (0..4).filter(|&i| node.set.contains(quartet[i])).collect();
        if inside.len() == 2 {
            let mask: u8 = inside.iter().map(|&i| 1u8 << i).sum();
            // the pair that has the first of the four
            return if mask & 1 == 1 { mask } else { 15 - mask };
        }
    }
    unreachable!("some branch of a tree splits four taxa two and two");
}

/// The speciation-driven quartet classes of one family tree: (the four species, the node where the pairs meet, the pairing of
/// the species), by the definition.
fn classes(r: &Rooted, n: usize) -> BTreeSet<([usize; 4], usize, u8)> {
    let t = Tree::of(r);
    let sets = r.species(n);
    let (speciation, _) = r.speciations(&sets);
    let leaves: Vec<usize> = (0..r.kids.len()).filter(|&v| r.kids[v].is_none()).collect();
    let species = |v: usize| r.taxon[v].expect("a species");
    let mut found = BTreeSet::new();
    for (i, &a) in leaves.iter().enumerate() {
        for (j, &b) in leaves.iter().enumerate().skip(i + 1) {
            for (k, &c) in leaves.iter().enumerate().skip(j + 1) {
                for &d in leaves.iter().skip(k + 1) {
                    let q = [a, b, c, d];
                    let spp = q.map(species);
                    let distinct: BTreeSet<usize> = spp.iter().copied().collect();
                    if distinct.len() < 4 {
                        continue;
                    }
                    // every three of the leaves meet at a speciation
                    let triples = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
                    if !triples
                        .iter()
                        .all(|[x, y, z]| speciation[t.lca(t.lca(q[*x], q[*y]), q[*z])])
                    {
                        continue;
                    }
                    // the pair that is together: the pair that meets deepest
                    let pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
                    let deepest = pairs
                        .iter()
                        .max_by_key(|(x, y)| t.depth[t.lca(q[*x], q[*y])])
                        .copied()
                        .expect("a pair");
                    let together: u8 = (1 << deepest.0) | (1 << deepest.1);
                    // the anchors: where pairs meet, but not where all four do; the quartet's node is where the anchors meet
                    let all = t.lca(t.lca(q[0], q[1]), t.lca(q[2], q[3]));
                    let anchors: BTreeSet<usize> = pairs
                        .iter()
                        .map(|(x, y)| t.lca(q[*x], q[*y]))
                        .filter(|&v| v != all)
                        .collect();
                    let anchors: Vec<usize> = anchors.into_iter().collect();
                    let psi = t.lca(anchors[0], anchors[1]);
                    // sorted by species: the pairing in the order of the sorted species
                    let mut order = [0usize, 1, 2, 3];
                    order.sort_by_key(|&x| spp[x]);
                    let sorted = order.map(|x| spp[x]);
                    let mapped: u8 = (0..4)
                        .filter(|&p| together >> order[p] & 1 == 1)
                        .map(|p| 1u8 << p)
                        .sum();
                    let key = if mapped & 1 == 1 { mapped } else { 15 - mapped };
                    found.insert((sorted, psi, key));
                }
            }
        }
    }
    found
}

fn score_by_definition(f: &Families, sp: &Sp) -> (u64, u64) {
    let n = f.taxa.len();
    let (mut all, mut agreeing) = (0, 0);
    for r in &f.rooted {
        for (sorted, _, key) in classes(r, n) {
            all += 1;
            agreeing += u64::from(pairing(sp, sorted) == key);
        }
    }
    (agreeing, all)
}

fn score_by_formula(f: &Families, sp: &Sp) -> u64 {
    sp.sum_with(|a: &Bits, b: &Bits| pscore::tripartition(&f.nodes, a, b))
}

// --- tests ------------------------------------------------------------------------------------------------------------------

#[test]
fn the_score_of_a_species_tree_is_the_number_of_classes_it_agrees_with() {
    let mut nonzero = 0;
    let mut with_duplications = 0;
    for seed in 1..=40u64 {
        let mut rng = Rng(seed * 977);
        let species = 5 + rng.below(3);
        let families: Vec<String> = (0..3)
            .map(|_| {
                let leaves = species + rng.below(5);
                random_family(species, leaves, &mut rng)
            })
            .collect();
        let text = families.join("\n");
        let f = match read(&text, 0) {
            Ok(f) => f,
            // no speciation with species of four kinds at all
            Err(_) => continue,
        };
        with_duplications += usize::from(f.duplications > 0);
        for _ in 0..3 {
            let sp_text = species_tree(&f.taxa, &mut rng);
            let sp = Sp::from_newick(&sp_text, &f.taxa).expect("a species tree");
            let (agreeing, all) = score_by_definition(&f, &sp);
            let found = score_by_formula(&f, &sp);
            assert_eq!(found, agreeing, "seed {seed}\n{text}\n{sp_text}");
            assert_eq!(f.classes, u128::from(all), "the classes, seed {seed}");
            nonzero += usize::from(agreeing > 0);
        }
    }
    // the comparison is of something, and of families that had duplications
    assert!(nonzero > 40, "{nonzero}");
    assert!(with_duplications > 15, "{with_duplications}");
}

#[test]
fn trees_with_one_copy_of_each_species_score_as_in_astral() {
    for seed in 1..=10u64 {
        let mut rng = Rng(seed * 31);
        let names: Vec<String> = SPECIES[..6].iter().map(|s| s.to_string()).collect();
        let text: String = (0..5)
            .map(|_| {
                format!(
                    "{};\n",
                    random_rooted(&shuffled(&names, &mut rng), &mut rng)
                )
            })
            .collect();
        let f = read(&text, 0).expect("families");
        assert_eq!(f.duplications, 0);
        let genes = super::genes::read(&text, 0).expect("genes");
        for _ in 0..4 {
            let sp_text = species_tree(&f.taxa, &mut rng);
            let sp = Sp::from_newick(&sp_text, &f.taxa).expect("a species tree");
            assert_eq!(score_by_formula(&f, &sp), sp.score(&genes), "seed {seed}");
        }
    }
}

fn shuffled(names: &[String], rng: &mut Rng) -> Vec<String> {
    let mut v = names.to_vec();
    for k in (1..v.len()).rev() {
        v.swap(k, rng.below(k + 1));
    }
    v
}

#[test]
fn a_duplication_is_a_node_whose_branches_share_species() {
    // rooted as written: the duplication is not at the root
    let f = read("(((A,B),(A,B)),(C,D));", 0).expect("a family");
    assert_eq!((f.duplications, f.rooted_by_cost, f.families), (1, 0, 1));
    // the speciations are the root and the three nodes (A,B), (A,B) and (C,D)
    assert_eq!(f.nodes.len(), 4);
}

#[test]
fn an_unrooted_tree_is_rooted_where_the_duplications_and_losses_cost_least() {
    // the tree of two copies of the species tree ((A,B),(C,D)), a duplication at its root, drawn with no root
    let f = read("((A,B),(C,D),((A,B),(C,D)));", 0).expect("a family");
    assert_eq!((f.rooted_by_cost, f.duplications), (1, 1));
    let r = &f.rooted[0];
    let sets = r.species(f.taxa.len());
    let (a, b) = r.kids[r.root].expect("a root with two branches");
    // both branches of the root hold all four species: the duplication is at the root
    assert_eq!((sets[a].count(), sets[b].count()), (4, 4));
}

#[test]
fn what_cannot_be_read_as_gene_families_is_refused_or_left_out_and_said() {
    let why = |text: &str| read(text, 0).expect_err("refused");
    assert!(
        why("((A,B,C),(D,E));").contains("not binary"),
        "{}",
        why("((A,B,C),(D,E));")
    );
    assert!(why("((A,B),(C,));").contains("has no name"));
    assert!(why("((A,B),(A,B));").contains("at least four"));
    // copies of a species that sit together are read, and make a quartet of the four species where the branches share none
    let f = read("(((A,A),(B,B)),((C,C),(D,D)));", 0).expect("a family");
    assert_eq!((f.duplications, f.classes), (4, 1));
    // a family with fewer than four species is left out and counted
    let f = read("((A,B),(C,D));\n((A,B),(A,C));", 0).expect("families");
    assert_eq!((f.families, f.skipped), (1, 1));
}

#[test]
fn the_single_copy_trees_that_are_sampled_have_each_species_once() {
    let f = read(
        "(((A,B),(A,B)),((C,D),(C,E)));\n(((A,B),(C,D)),((A,B),(C,D)));",
        0,
    )
    .expect("families");
    assert!(!f.sample.trees.is_empty());
    for t in &f.sample.trees {
        let leaves: Vec<usize> = t.nodes.iter().filter_map(|n| n.taxon).collect();
        let distinct: BTreeSet<usize> = leaves.iter().copied().collect();
        assert_eq!(leaves.len(), distinct.len(), "{leaves:?}");
        assert!(leaves.len() >= 4);
    }
    // the same trees every time
    let again = read(
        "(((A,B),(A,B)),((C,D),(C,E)));\n(((A,B),(C,D)),((A,B),(C,D)));",
        0,
    )
    .expect("families");
    assert_eq!(again.sample.trees.len(), f.sample.trees.len());
}

// --- the search --------------------------------------------------------------------------------------------------------------

use super::{infer, score_tree, Options, Weights};
use crate::seq::newick::NwTree;
use crate::seq::treetools::robinson_foulds;

fn families() -> Options {
    Options {
        families: true,
        ..Options::default()
    }
}

fn tree(text: &str) -> NwTree {
    NwTree::parse_all(text).expect("a tree").remove(0)
}

/// How many splits two trees on the same species do not share.
fn apart(a: &str, b: &str) -> usize {
    robinson_foulds(&tree(a), &tree(b))
        .expect("the same species")
        .distance
}

#[test]
fn the_species_tree_is_found_from_families_that_have_duplications_and_losses() {
    // two copies of the species tree ((A,B),(C,(D,E))) after a duplication at the root, and some of the same with losses
    let truth = "((A,B),(C,(D,E)));";
    let text = "(((A,B),(C,(D,E))),((A,B),(C,(D,E))));\n\
                (((A,B),(C,(D,E))),((A,B),(D,E)));\n\
                (((B),(C,(D,E))),((A,B),(C,(D,E))));\n\
                ((A,(B,(C,(D,E)))),((A,B),(D,E)));\n"
        .replace("(B),", "B,");
    let out = infer(&text, &families()).expect("a species tree");
    assert!(out.families && out.branches.is_empty());
    assert_eq!(out.weights, Weights::None);
    assert_eq!(apart(&out.newick, truth), 0, "{}", out.newick);
    assert!(
        out.score > 0 && out.normalized > 0.5,
        "{} {}",
        out.score,
        out.normalized
    );
    assert!(
        out.log[0].contains("4 gene family trees on 5 species"),
        "{:?}",
        out.log
    );
    // the tree that it found scores best of the trees that were tried
    let wrong = score_tree(&text, "((A,C),(B,(D,E)));", &families()).expect("a score");
    assert!(wrong.score < out.score, "{} {}", wrong.score, out.score);
    let same = score_tree(&text, truth, &families()).expect("a score");
    assert_eq!(same.score, out.score);
}

#[test]
fn without_duplications_the_families_give_astrals_score_and_tree() {
    for seed in 1..=5u64 {
        let mut rng = Rng(seed * 7);
        let names: Vec<String> = SPECIES[..6].iter().map(|s| s.to_string()).collect();
        let text: String = (0..8)
            .map(|_| {
                format!(
                    "{};\n",
                    random_rooted(&shuffled(&names, &mut rng), &mut rng)
                )
            })
            .collect();
        let plain = infer(&text, &Options::default()).expect("a tree");
        let pro = infer(&text, &families()).expect("a tree");
        assert_eq!(pro.score, plain.score, "seed {seed}");
        assert!(
            (pro.normalized - plain.normalized).abs() < 1e-12,
            "seed {seed}"
        );
        assert_eq!(
            apart(&pro.newick, &plain.newick),
            0,
            "seed {seed}: {} against {}",
            pro.newick,
            plain.newick
        );
    }
}

#[test]
fn families_and_weights_are_not_used_together_and_say_so() {
    let o = Options {
        families: true,
        weights: Weights::Support,
        ..Options::default()
    };
    let why = infer("((A,B),(C,D));", &o).expect_err("refused");
    assert!(why.contains("cannot be used together"), "{why}");
    assert!(score_tree("((A,B),(C,D));", "((A,B),(C,D));", &o).is_err());
    // single-copy trees without families are still refused for a repeated leaf
    let why = infer("(((A,B),(A,B)),(C,D));", &Options::default()).expect_err("refused");
    assert!(why.contains("more than once"), "{why}");
}
