//! The species tree as Newick text, with each inner branch annotated the way ASTRAL's `-t` option asks: nothing (0), the
//! quartet support of the branch in percent (1), all the quantities (2), the local posterior probability (3, the usual one),
//! the three posteriors (4) or the three quartet supports (8). Lengths are in coalescent units; leaves have none.

use super::species::{Branch, Sp};

/// A number written as ASTRAL's Java does for the `#.##` format: at most two decimals, no trailing zeros.
fn two(x: f64) -> String {
    let text = format!("{x:.2}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    if text.is_empty() || text == "-0" {
        "0".into()
    } else {
        text.into()
    }
}

/// A number with all its digits, as a float is written in Java (`2.0`, `1.5E-6`).
fn full(x: f64) -> String {
    if x.is_infinite() {
        return "Infinity".into();
    }
    format!("{x:?}").replace('e', "E")
}

fn label(b: &Branch, level: u8) -> String {
    let length = full(b.length);
    let three = |v: [f64; 3], name: &str| {
        let parts: Vec<String> = v
            .iter()
            .enumerate()
            .map(|(i, x)| format!("{name}{}={}", i + 1, two(*x)))
            .collect();
        format!("'[{}]'", parts.join(";"))
    };
    match level {
        0 => String::new(),
        1 => format!("{}:{length}", two(b.q[0] * 100.0)),
        2 => {
            let each = |name: &str, v: [f64; 3]| {
                (0..3)
                    .map(|i| format!("{name}{}={}", i + 1, full(v[i])))
                    .collect::<Vec<_>>()
                    .join(";")
            };
            format!(
                "'[{};{};{};QC={};EN={}]':{length}",
                each("q", b.q),
                each("f", b.f),
                each("pp", b.pp),
                b.quartets,
                full(b.en)
            )
        }
        4 => format!("{}:{length}", three(b.pp, "pp")),
        8 => format!("{}:{length}", three(b.q, "q")),
        _ => format!("{}:{length}", two(b.pp[0])),
    }
}

/// The tree as Newick, hung from its root taxon, with the annotations for `level` on each branch that has one.
pub fn newick(sp: &Sp, taxa: &[String], branches: &[Branch], level: u8) -> String {
    let mut by_node: Vec<Option<&Branch>> = vec![None; sp.nodes.len()];
    for b in branches {
        by_node[b.node] = Some(b);
    }
    let mut out = format!("({},", name(&taxa[sp.root_taxon]));
    write(sp, taxa, &by_node, level, &mut out);
    out.push_str(");");
    out
}

/// A species name as Newick writes it: as it is, or in single quotes (a quote inside is written twice) when it holds a space,
/// a bracket, a comma, a colon, a semicolon or a quote, which would otherwise end or break the name.
fn name(species: &str) -> String {
    let plain = !species.is_empty()
        && !species
            .chars()
            .any(|c| c.is_whitespace() || "(),:;[]'".contains(c));
    if plain {
        species.to_string()
    } else {
        format!("'{}'", species.replace('\'', "''"))
    }
}

/// What is still to be written, in the order it is taken from the stack.
enum Step {
    /// A node: its brackets and children, or its name.
    Node(usize),
    Comma,
    /// The bracket that closes an inner node, and the label of its branch.
    Close(usize),
}

/// Write the tree below the root without recursion: a tree of thousands of species can be as deep as it is wide, and the
/// thread that writes it has a small stack.
fn write(sp: &Sp, taxa: &[String], by_node: &[Option<&Branch>], level: u8, out: &mut String) {
    let mut stack = vec![Step::Node(sp.root)];
    while let Some(step) = stack.pop() {
        match step {
            Step::Comma => out.push(','),
            Step::Close(i) => {
                out.push(')');
                if let Some(branch) = by_node[i] {
                    out.push_str(&label(branch, level));
                }
            }
            Step::Node(i) => match (sp.nodes[i].kids, sp.nodes[i].taxon) {
                (Some((a, b)), _) => {
                    out.push('(');
                    // taken in the reverse order: the first child, a comma, the second child, the closing bracket
                    stack.extend([Step::Close(i), Step::Node(b), Step::Comma, Step::Node(a)]);
                }
                (None, Some(t)) => out.push_str(&name(&taxa[t])),
                (None, None) => {}
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_are_written_as_javas_decimal_format_and_double_do() {
        assert_eq!(two(0.9247), "0.92");
        assert_eq!(two(1.0), "1");
        assert_eq!(two(63.75), "63.75");
        assert_eq!(two(0.0), "0");
        assert_eq!(two(0.004), "0");
        assert_eq!(full(2.0), "2.0");
        assert_eq!(full(0.2876820724517809), "0.2876820724517809");
        assert_eq!(full(1.3691758792688155e-6), "1.3691758792688155E-6");
        assert_eq!(full(f64::INFINITY), "Infinity");
    }

    #[test]
    fn a_tree_as_deep_as_it_is_wide_is_written_without_running_out_of_stack() {
        // a caterpillar of 100,000 leaves, all of the one species: only the shape matters to the writer, and nodes 2k + 1
        // are the new leaves while nodes 2k + 2 join the tree so far (node 2k) with them
        let leaves = 100_000;
        let mut kids = vec![None];
        let mut taxon = vec![Some(0)];
        for k in 0..leaves - 1 {
            kids.extend([None, Some((2 * k, 2 * k + 1))]);
            taxon.extend([Some(0), None]);
        }
        let sp = Sp::from_shape(4, 1, &kids, &taxon, kids.len() - 1);
        let names: Vec<String> = ["a", "b", "c", "d"].map(String::from).to_vec();
        let text = newick(&sp, &names, &[], 3);
        // the root species, the deepest cherry first, and the last leaf added at the top
        assert!(
            text.starts_with("(b,((((") && text.contains("(a,a),a),a)") && text.ends_with(",a));"),
            "{}",
            &text[..20]
        );
        assert_eq!(text.matches('(').count(), leaves);
        assert_eq!(text.matches(')').count(), leaves);
    }

    #[test]
    fn names_that_would_break_the_text_are_quoted_and_read_back_as_they_were() {
        assert_eq!(name("Homo_sapiens"), "Homo_sapiens");
        assert_eq!(name("Homo sapiens"), "'Homo sapiens'");
        assert_eq!(name("Nomascus (gibbon)"), "'Nomascus (gibbon)'");
        assert_eq!(name("a,b:c;d[e]"), "'a,b:c;d[e]'");
        assert_eq!(name("it's"), "'it''s'");
        assert_eq!(name(""), "''");
        // a whole tree of such names is read by the same reader that reads the gene trees, to the same names
        let kids = [None, None, Some((0, 1)), None, Some((2, 3))];
        let taxon = [Some(0), Some(1), None, Some(2), None];
        let sp = Sp::from_shape(4, 3, &kids, &taxon, 4);
        let names: Vec<String> = ["Homo sapiens", "Nomascus (gibbon)", "it's, odd:", "Pan"]
            .map(String::from)
            .to_vec();
        let text = newick(&sp, &names, &[], 3);
        let back = &crate::seq::newick::NwTree::parse_all(&text).expect("the text reads back")[0];
        let mut read: Vec<&str> = back
            .nodes
            .iter()
            .filter(|n| n.children.is_empty())
            .map(|n| n.name.as_str())
            .collect();
        read.sort_unstable();
        let mut wrote: Vec<&str> = names.iter().map(String::as_str).collect();
        wrote.sort_unstable();
        assert_eq!(read, wrote, "{text}");
    }

    #[test]
    fn the_children_and_labels_come_in_their_order() {
        // ((a,b)x,c)y hung from d: the labels sit after the closing brackets of their own nodes
        let kids = [None, None, Some((0, 1)), None, Some((2, 3))];
        let taxon = [Some(0), Some(1), None, Some(2), None];
        let sp = Sp::from_shape(4, 3, &kids, &taxon, 4);
        let names: Vec<String> = ["a", "b", "c", "d"].map(String::from).to_vec();
        let branch = |node, pp| Branch {
            node,
            f: [1.0, 0.0, 0.0],
            q: [1.0, 0.0, 0.0],
            pp: [pp, 0.0, 0.0],
            en: 1.0,
            quartets: 1,
            length: 2.0,
        };
        let text = newick(&sp, &names, &[branch(2, 0.5), branch(4, 0.25)], 3);
        assert_eq!(text, "(d,((a,b)0.5:2.0,c)0.25:2.0);");
        assert_eq!(newick(&sp, &names, &[], 3), "(d,((a,b),c));");
    }
}
