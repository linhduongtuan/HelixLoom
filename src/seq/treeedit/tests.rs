use super::*;
use std::collections::BTreeMap;

fn tree(text: &str) -> NwTree {
    NwTree::parse_all(text).unwrap().remove(0)
}

fn find(t: &NwTree, name: &str) -> usize {
    t.nodes
        .iter()
        .position(|n| n.name == name)
        .unwrap_or_else(|| panic!("no node {name}"))
}

fn tips_under(t: &NwTree, i: usize) -> Vec<String> {
    let mut out = Vec::new();
    let mut stack = vec![i];
    while let Some(k) = stack.pop() {
        if t.is_leaf(k) {
            out.push(t.nodes[k].name.clone());
        }
        stack.extend(t.nodes[k].children.iter());
    }
    out.sort();
    out
}

/// Every branch of the tree without regard to where its root is: the tips on the side that has not the least tip, with the length of the
/// branch and its `posterior`. The two branches of a root with two children are one branch.
fn edges(t: &NwTree) -> BTreeMap<Vec<String>, (f64, Option<f64>)> {
    let all = tips_under(t, t.root);
    let mut out: BTreeMap<Vec<String>, (f64, Option<f64>)> = BTreeMap::new();
    for i in 0..t.nodes.len() {
        if i == t.root {
            continue;
        }
        let mut side = tips_under(t, i);
        if side.contains(&all[0]) {
            side = all.iter().filter(|n| !side.contains(n)).cloned().collect();
        }
        let e = out.entry(side).or_insert((0.0, None));
        e.0 += t.nodes[i].dist;
        e.1 = e.1.or(t.attr_number(i, "posterior"));
    }
    out
}

fn same_edges(a: &NwTree, b: &NwTree) {
    let (ea, eb) = (edges(a), edges(b));
    assert_eq!(ea.keys().collect::<Vec<_>>(), eb.keys().collect::<Vec<_>>());
    for (k, (la, pa)) in &ea {
        let (lb, pb) = eb[k];
        assert!((la - lb).abs() < 1e-9, "{k:?}: {la} against {lb}");
        // a branch that has a value keeps it
        assert_eq!(pa.is_some(), pb.is_some(), "{k:?}");
        if let (Some(x), Some(y)) = (pa, pb) {
            assert!((x - y).abs() < 1e-12, "{k:?}");
        }
    }
}

const THREE: &str =
    "(((A:1,B:1)x[&posterior=0.9]:2,(C:1,D:1)y[&posterior=0.8]:3)z[&posterior=0.7]:4,E:5,F:6);";
const TWO: &str = "(((A:1,B:1)x[&posterior=0.9]:2,(C:1,D:1)y[&posterior=0.8]:3)z[&posterior=0.7]:4,(E:5,F:6)w[&posterior=0.7]:1);";

#[test]
fn rooting_on_a_branch_keeps_every_branch_with_its_length_and_value() {
    for text in [THREE, TWO] {
        let t = tree(text);
        for name in ["x", "y", "z", "A", "D", "E", "F"] {
            for f in [0.5, 0.25] {
                let Some(i) = t.nodes.iter().position(|n| n.name == name) else {
                    continue;
                };
                let r = reroot_above(&t, i, f).unwrap();
                same_edges(&t, &r);
                assert_eq!(r.leaf_count(), 6, "{name}");
                assert_eq!(
                    r.nodes.len(),
                    11,
                    "{name}: two nodes more than tips... one root, no unary node"
                );
                // the new root has two children and the clade is one of them
                let kids = &r.nodes[r.root].children;
                assert_eq!(kids.len(), 2, "{name}");
                let moved = find(&r, name);
                assert!(kids.contains(&moved), "{name}");
                assert!(
                    r.nodes
                        .iter()
                        .all(|n| n.parent.is_some() == (n != &r.nodes[r.root])),
                    "{name}"
                );
            }
        }
    }
}

#[test]
fn the_branch_above_is_shared_between_the_two_sides_of_the_new_root() {
    let t = tree(THREE);
    let x = find(&t, "x");
    let r = reroot_above(&t, x, 0.25).unwrap();
    let (below, above) = (find(&r, "x"), r.nodes[r.root].children[1]);
    assert!((r.nodes[below].dist - 0.5).abs() < 1e-12);
    assert!((r.nodes[above].dist - 1.5).abs() < 1e-12);
    // both halves carry the value written on the branch
    assert_eq!(r.attr_number(below, "posterior"), Some(0.9));
    assert_eq!(r.attr_number(above, "posterior"), Some(0.9));
    // a tip's name stays with the tip and is not the label of the other half
    let a = find(&r, "A");
    let r2 = reroot_above(&r, a, 0.5).unwrap();
    assert_eq!(r2.nodes.iter().filter(|n| n.name == "A").count(), 1);
    // root on the branch of the root's child of a tree with a two-child root: the two root branches are one, shared again
    let two = tree(TWO);
    let w = find(&two, "w");
    let back = reroot_above(&two, w, 0.5).unwrap();
    let (a, b) = (
        back.nodes[back.root].children[0],
        back.nodes[back.root].children[1],
    );
    assert!((back.nodes[a].dist - back.nodes[b].dist).abs() < 1e-12);
    assert!(
        (back.nodes[a].dist + back.nodes[b].dist - 5.0).abs() < 1e-12,
        "z's 4 and w's 1"
    );
}

#[test]
fn the_root_cannot_be_rooted_on_and_a_wild_fraction_is_kept_in_range() {
    let t = tree(THREE);
    assert!(reroot_above(&t, t.root, 0.5)
        .unwrap_err()
        .contains("no branch above"));
    assert!(reroot_above(&t, 99, 0.5)
        .unwrap_err()
        .contains("no such node"));
    let x = find(&t, "x");
    for f in [-3.0, 7.0, f64::NAN] {
        let r = reroot_above(&t, x, f).unwrap();
        same_edges(&t, &r);
        // the two sides of the new root share the branch's 2: neither has less than none or more than all
        let below = find(&r, "x");
        let above = r.nodes[r.root]
            .children
            .iter()
            .copied()
            .find(|&c| c != below)
            .unwrap();
        let (lo, hi) = (r.nodes[below].dist, r.nodes[above].dist);
        assert!(
            lo >= 0.0 && hi >= 0.0 && (lo + hi - 2.0).abs() < 1e-12,
            "{f}: {lo} and {hi}"
        );
        if f.is_nan() {
            assert_eq!((lo, hi), (1.0, 1.0));
        }
    }
}

#[test]
fn rotating_reverses_one_node_and_ordering_sorts_by_the_tips_beyond() {
    let t = tree("((A,B,C)x,D,(E,(F,G))y);");
    let r = rotate(&t, find(&t, "x"));
    let kids: Vec<&str> = r.nodes[find(&r, "x")]
        .children
        .iter()
        .map(|&c| r.nodes[c].name.as_str())
        .collect();
    assert_eq!(kids, ["C", "B", "A"]);
    assert_eq!(rotate(&t, 99), t);
    let names = |t: &NwTree, i: usize| -> Vec<String> {
        t.nodes[i]
            .children
            .iter()
            .map(|&c| tips_under(t, c).len().to_string())
            .collect()
    };
    assert_eq!(names(&order(&t, true), t.root), ["1", "3", "3"]);
    assert_eq!(names(&order(&t, false), t.root), ["3", "3", "1"]);
    // equal sizes keep their order
    let up = order(&t, true);
    assert_eq!(
        tips_under(&up, up.nodes[up.root].children[1]),
        ["A", "B", "C"]
    );
    // ordering is only of the order: the branches are the same
    same_edges(&t, &up);
}

#[test]
fn branches_below_a_value_are_collapsed_into_their_parent() {
    let t = tree(THREE);
    let (c, gone) = collapse_below(&t, "posterior", 0.85);
    assert_eq!(gone, 2, "y 0.8 and z 0.7 are below, x 0.9 is not");
    assert!(c.nodes.iter().all(|n| n.name != "y" && n.name != "z"));
    // C and D hang from the root now, with y's and z's lengths added: 1 + 3 + 4
    let d = find(&c, "D");
    assert!((c.nodes[d].dist - 8.0).abs() < 1e-12);
    assert_eq!(c.nodes[find(&c, "x")].parent, Some(c.root));
    assert!(
        (c.nodes[find(&c, "x")].dist - 6.0).abs() < 1e-12,
        "x's 2 and z's 4"
    );
    assert_eq!(c.leaf_count(), 6);
    // nothing is below a low value, and the tips and the root are never taken out
    assert_eq!(collapse_below(&t, "posterior", 0.0).1, 0);
    assert_eq!(collapse_below(&t, "posterior", 5.0).1, 3);
    // a branch with no value stays, and the label is a value when no key is given
    assert_eq!(collapse_below(&t, "rate", 5.0).1, 0);
    let labelled = tree("(((A,B)90,(C,D)40)70,E,F);");
    assert_eq!(collapse_below(&labelled, "", 75.0).1, 2);
}

#[test]
fn a_clade_can_be_cut_out_and_tips_renamed_from_a_table() {
    let t = tree(THREE);
    let x = extract_clade(&t, find(&t, "z")).unwrap();
    assert_eq!(tips_under(&x, x.root), ["A", "B", "C", "D"]);
    assert_eq!(x.nodes[x.root].dist, 0.0);
    assert_eq!(x.attr_number(x.root, "posterior"), Some(0.7));
    assert!(extract_clade(&t, 99).is_err());
    let table = parse_name_table("A\tHomo sapiens\n\nB, Pan\n\"C\";Gorilla\nonly one\nD\t\n");
    assert_eq!(table.len(), 3);
    let (r, n) = rename_tips(&t, &table);
    assert_eq!(n, 3);
    let mut names: Vec<&str> = r
        .nodes
        .iter()
        .filter(|n| n.children.is_empty())
        .map(|n| n.name.as_str())
        .collect();
    names.sort_unstable();
    assert_eq!(names, ["D", "E", "F", "Gorilla", "Homo sapiens", "Pan"]);
}

const REFERENCE: &str = "((A,B),(C,D),(E,F));";

fn replicates() -> Vec<NwTree> {
    [
        "((A,B),(C,D),(E,F));",
        "((A,B),(C,E),(D,F));",
        "(((A,B),(C,D)),E,F);",
        "((A,C),(B,D),(E,F));",
    ]
    .iter()
    .map(|t| tree(t))
    .collect()
}

#[test]
fn how_often_a_branch_is_among_the_replicates_and_how_often_a_rival_is() {
    let t = tree(REFERENCE);
    let s = frequencies(&t, &replicates()).unwrap();
    let of = |tips: &[&str]| {
        let i = (0..t.nodes.len())
            .find(|&i| tips_under(&t, i) == tips)
            .unwrap();
        s[i].unwrap()
    };
    let ab = of(&["A", "B"]);
    assert_eq!((ab.frequency, ab.conflict), (0.75, 0.25));
    let cd = of(&["C", "D"]);
    assert_eq!((cd.frequency, cd.conflict), (0.5, 0.25));
    let ef = of(&["E", "F"]);
    assert_eq!((ef.frequency, ef.conflict), (0.75, 0.25));
    // the root and the tips have none
    assert!(s[t.root].is_none());
    assert!(t
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| n.children.is_empty())
        .all(|(i, _)| s[i].is_none()));
}

#[test]
fn the_rooting_of_a_replicate_does_not_matter_and_a_branch_counts_once_per_tree() {
    let t = tree(REFERENCE);
    // the same tree rooted elsewhere, and written with its root branches as two
    let other = [
        tree("(A,(B,((C,D),(E,F))));"),
        tree("((E,F),((C,D),(A,B)));"),
    ];
    let s = frequencies(&t, &other).unwrap();
    assert!(
        s.iter()
            .flatten()
            .all(|x| x.frequency == 1.0 && x.conflict == 0.0),
        "{s:?}"
    );
}

#[test]
fn tips_that_a_tree_lacks_are_left_out_and_too_few_shared_tips_is_an_error() {
    let t = tree(REFERENCE);
    // F is missing and G is extra: A to E remain, and {E,F} is no branch of five
    let s = frequencies(&t, &[tree("((A,B),(C,D),(E,G));")]).unwrap();
    let ab = (0..t.nodes.len())
        .find(|&i| tips_under(&t, i) == ["A", "B"])
        .unwrap();
    assert_eq!(s[ab].unwrap().frequency, 1.0);
    assert!(frequencies(&t, &[tree("((A,B),(C,X));")])
        .unwrap_err()
        .contains("share 3 tips"));
    assert!(frequencies(&t, &[]).unwrap_err().contains("no other trees"));
}

#[test]
fn the_value_on_the_same_branch_of_another_tree_is_found_wherever_it_is_rooted() {
    let t = tree(REFERENCE);
    let other = tree("(((A,B)0.99,(C,D)0.5)0.77,E,F);");
    let v = values_from(&t, &other, "").unwrap();
    let of = |tips: &[&str]| {
        v[(0..t.nodes.len())
            .find(|&i| tips_under(&t, i) == tips)
            .unwrap()]
    };
    assert_eq!(of(&["A", "B"]), Some(0.99));
    assert_eq!(of(&["C", "D"]), Some(0.5));
    // {E,F} is the other side of the node that has A to D below it
    assert_eq!(of(&["E", "F"]), Some(0.77));
    assert_eq!(v[t.root], None);
    // an annotation as the value, and a branch the other tree lacks has none
    let noted = tree("((A,B)[&posterior=0.4],(C,E),(D,F));");
    let w = values_from(&t, &noted, "posterior").unwrap();
    let ab = (0..t.nodes.len())
        .find(|&i| tips_under(&t, i) == ["A", "B"])
        .unwrap();
    let cd = (0..t.nodes.len())
        .find(|&i| tips_under(&t, i) == ["C", "D"])
        .unwrap();
    assert_eq!((w[ab], w[cd]), (Some(0.4), None));
}

#[test]
fn values_are_written_as_annotations_and_replace_the_ones_there() {
    let t = tree(REFERENCE);
    let mut values = vec![None; t.nodes.len()];
    values[1] = Some(0.5);
    let once = with_values(&t, "boot", &values);
    assert_eq!(once.attr_number(1, "boot"), Some(0.5));
    values[1] = Some(0.75);
    values[2] = Some(1.0);
    let twice = with_values(&once, "boot", &values);
    assert_eq!(
        twice.nodes[1]
            .attrs
            .iter()
            .filter(|(k, _)| k == "boot")
            .count(),
        1
    );
    assert_eq!(
        (twice.attr_number(1, "boot"), twice.attr_number(2, "boot")),
        (Some(0.75), Some(1.0))
    );
}
