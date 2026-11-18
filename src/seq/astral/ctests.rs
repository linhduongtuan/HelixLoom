//! The greedy consensus of species trees: checked against a plain version written with sorted lists of numbers, on random trees.

use super::bits::Bits;
use super::consensus::{greedy, newick};
use crate::seq::newick::NwTree;
use crate::seq::tree::Rng;
use std::collections::HashMap;

/// The inner splits of a random binary tree on the taxa 0..n, hung from taxon 0: each is the set below a branch (taxon 0 is
/// never in one). Made by joining two random clades at a time.
fn random_tree(rng: &mut Rng, n: usize) -> Vec<Bits> {
    let mut clades: Vec<Vec<usize>> = (1..n).map(|t| vec![t]).collect();
    let mut splits = Vec::new();
    while clades.len() > 1 {
        let a = clades.swap_remove(rng.below(clades.len()));
        let b = clades.swap_remove(rng.below(clades.len()));
        let mut joined = [a, b].concat();
        joined.sort_unstable();
        if clades.is_empty() {
            break; // the last join is everything but taxon 0: not an inner branch
        }
        splits.push(Bits::of(n, joined.iter().copied()));
        clades.push(joined);
    }
    splits
}

fn listed(b: &Bits) -> Vec<usize> {
    b.iter().collect()
}

/// Two splits fit in one tree when one is inside the other or they share nothing.
fn compatible(a: &[usize], b: &[usize]) -> bool {
    let inside = |x: &[usize], y: &[usize]| x.iter().all(|t| y.contains(t));
    let shared = a.iter().any(|t| b.contains(t));
    !shared || inside(a, b) || inside(b, a)
}

/// The same choice as `greedy`, done on sorted lists: by how many trees have the split (the most first), the lists in order for
/// ties, a split kept when it fits all that are kept.
fn plain(trees: &[Vec<Bits>]) -> Vec<(Vec<usize>, usize)> {
    let mut count: HashMap<Vec<usize>, usize> = HashMap::new();
    for t in trees {
        for s in t {
            *count.entry(listed(s)).or_insert(0) += 1;
        }
    }
    let mut all: Vec<(Vec<usize>, usize)> = count.into_iter().collect();
    all.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let mut kept: Vec<(Vec<usize>, usize)> = Vec::new();
    for (s, c) in all {
        if kept.iter().all(|(k, _)| compatible(&s, k)) {
            kept.push((s, c));
        }
    }
    kept
}

#[test]
fn trees_that_are_all_the_same_give_that_tree_with_every_split_in_all_of_them() {
    let mut rng = Rng::new(5);
    let one = random_tree(&mut rng, 9);
    let got = greedy(&vec![one.clone(); 7]);
    assert_eq!(
        got.len(),
        one.len(),
        "n - 3 inner branches of a tree of 9 taxa hung from one of them"
    );
    assert!(got.iter().all(|(_, c)| *c == 7));
    let mut a: Vec<Vec<usize>> = got.iter().map(|(s, _)| listed(s)).collect();
    let mut b: Vec<Vec<usize>> = one.iter().map(listed).collect();
    a.sort();
    b.sort();
    assert_eq!(a, b);
}

#[test]
fn a_split_that_most_trees_have_is_kept_and_one_that_cannot_go_with_it_is_not() {
    // taxa 0 to 5, hung from 0. {1,2} in 3 of 4 trees; {2,3} (not fitting {1,2}) in 1; {1,2,3} in 2; {4,5} in 3
    let s = |items: &[usize]| Bits::of(6, items.iter().copied());
    let trees = vec![
        vec![s(&[1, 2]), s(&[1, 2, 3]), s(&[4, 5])],
        vec![s(&[1, 2]), s(&[1, 2, 3]), s(&[4, 5])],
        vec![s(&[1, 2]), s(&[4, 5])],
        vec![s(&[2, 3]), s(&[1, 2, 3])],
    ];
    let got = greedy(&trees);
    let as_lists: Vec<(Vec<usize>, usize)> = got.iter().map(|(b, c)| (listed(b), *c)).collect();
    // {1,2}, {1,2,3} and {4,5} are each in three trees and fit together; {2,3} is in one and does not fit {1,2}
    assert_eq!(as_lists.len(), 3, "{as_lists:?}");
    assert!(as_lists.contains(&(vec![1, 2], 3)));
    assert!(as_lists.contains(&(vec![4, 5], 3)));
    assert!(as_lists.contains(&(vec![1, 2, 3], 3)));
    assert!(
        !as_lists.iter().any(|(l, _)| l == &vec![2, 3]),
        "{as_lists:?}"
    );
}

#[test]
fn the_choice_is_the_same_as_the_plain_one_and_follows_the_rules_of_a_greedy_consensus_on_random_trees(
) {
    let mut rng = Rng::new(17);
    for (n, count) in [(6usize, 5usize), (8, 9), (10, 25), (12, 40), (20, 30)] {
        let trees: Vec<Vec<Bits>> = (0..count).map(|_| random_tree(&mut rng, n)).collect();
        let got = greedy(&trees);
        let mine: Vec<(Vec<usize>, usize)> = got.iter().map(|(b, c)| (listed(b), *c)).collect();
        assert_eq!(mine, plain(&trees), "n={n} trees={count}");
        // every pair fits in one tree
        for (i, (a, _)) in mine.iter().enumerate() {
            for (b, _) in &mine[..i] {
                assert!(compatible(a, b), "{a:?} {b:?}");
            }
        }
        // every split in more than half of the trees is kept (the majority rule's tree is inside the greedy one)
        let mut seen: HashMap<Vec<usize>, usize> = HashMap::new();
        for t in &trees {
            for s in t {
                *seen.entry(listed(s)).or_insert(0) += 1;
            }
        }
        for (s, c) in &seen {
            if 2 * c > count {
                assert!(
                    mine.iter().any(|(k, _)| k == s),
                    "the majority split {s:?} ({c} of {count}) is missing"
                );
            }
            // and one that is left out cannot be added: it does not fit one that is kept
            if !mine.iter().any(|(k, _)| k == s) {
                assert!(
                    mine.iter().any(|(k, _)| !compatible(s, k)),
                    "{s:?} could have been added"
                );
            }
        }
        // a tree of n taxa has at most n - 3 inner branches
        assert!(mine.len() <= n - 3);
    }
}

#[test]
fn the_result_does_not_depend_on_the_order_of_the_trees() {
    let mut rng = Rng::new(3);
    let mut trees: Vec<Vec<Bits>> = (0..12).map(|_| random_tree(&mut rng, 9)).collect();
    let before = greedy(&trees);
    trees.reverse();
    trees.rotate_left(5);
    assert_eq!(greedy(&trees), before);
    assert!(greedy(&[]).is_empty());
    assert!(greedy(&[Vec::new(), Vec::new()]).is_empty());
}

fn names(n: usize) -> Vec<String> {
    (0..n)
        .map(|i| ((b'A' + i as u8) as char).to_string())
        .collect()
}

#[test]
fn the_tree_is_written_hung_from_the_root_taxon_with_the_share_of_trees_on_each_branch() {
    let s = |items: &[usize]| Bits::of(5, items.iter().copied());
    let chosen = vec![(s(&[1, 2]), 8), (s(&[1, 2, 3]), 6)];
    assert_eq!(newick(&chosen, 10, &names(5), 0), "(A,(((B,C)80,D)60,E));");
    // hung from another taxon, the clades are of the others; the children of a node come in the order of their lowest taxon
    let from_c = vec![(Bits::of(5, [3, 4]), 5)];
    assert_eq!(newick(&from_c, 8, &names(5), 2), "(C,(A,B,(D,E)62.5));");
    // nothing chosen: a star
    assert_eq!(newick(&[], 10, &names(5), 0), "(A,(B,C,D,E));");
    // a share is written with one decimal at most
    assert_eq!(
        newick(&[(s(&[1, 2]), 1)], 3, &names(5), 0),
        "(A,((B,C)33.3,D,E));"
    );
}

#[test]
fn what_is_written_reads_back_as_the_splits_that_were_chosen_with_their_shares() {
    let mut rng = Rng::new(29);
    let n = 14;
    let trees: Vec<Vec<Bits>> = (0..30).map(|_| random_tree(&mut rng, n)).collect();
    let chosen = greedy(&trees);
    let text = newick(&chosen, trees.len(), &names(n), 0);
    let tree = &NwTree::parse_all(&text).expect("a tree")[0];
    // the splits of the text, by the leaves below each inner branch, on the side without taxon A, with the labels
    let mut found: Vec<(Vec<String>, String)> = Vec::new();
    for (i, node) in tree.nodes.iter().enumerate() {
        if node.children.is_empty() || node.parent.is_none() {
            continue;
        }
        let mut stack = vec![i];
        let mut below = Vec::new();
        while let Some(j) = stack.pop() {
            if tree.nodes[j].children.is_empty() {
                below.push(tree.nodes[j].name.clone());
            }
            stack.extend(tree.nodes[j].children.iter().copied());
        }
        below.sort();
        if !below.contains(&"A".to_string()) && !node.name.is_empty() {
            found.push((below, node.name.clone()));
        }
    }
    found.sort();
    let mut want: Vec<(Vec<String>, String)> = chosen
        .iter()
        .map(|(b, c)| {
            let mut leaves: Vec<String> = b.iter().map(|t| names(n)[t].clone()).collect();
            leaves.sort();
            let percent = 100.0 * *c as f64 / 30.0;
            let label = format!("{percent:.1}")
                .trim_end_matches('0')
                .trim_end_matches('.')
                .to_string();
            (leaves, label)
        })
        .collect();
    want.sort();
    assert_eq!(found, want);
}
