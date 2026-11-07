//! Newick trees: reading text into a tree of nodes with names and branch lengths, and writing it
//! back. The base of the tree tools (the ETE / DendroPy side of the platform) and of PhyloCNN.

use std::collections::VecDeque;

#[derive(Clone, Debug, PartialEq)]
pub struct NwNode {
    /// Leaf name, or the label of an inner node (often a support value).
    pub name: String,
    /// Length of the branch to the parent; 0 for a root without one.
    pub dist: f64,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct NwTree {
    pub nodes: Vec<NwNode>,
    pub root: usize,
}

/// The deepest nesting of parentheses a tree may have. The programs that read, draw and write a tree go down it by
/// recursion, and a tree nested far deeper than any real one would use up their stack (and end the program), so text
/// that nests deeper is refused. A tree of two thousand levels is a ladder of two thousand tips.
// ponytail: a ceiling on the depth instead of reading, drawing and writing trees without recursion; a ladder of more tips
// than this (a tree of thousands of near-identical sequences can be one) cannot be opened. Walking the tree with a
// stack of its own would lift it.
pub const MAX_DEPTH: usize = 2_000;

struct Reader<'a> {
    text: &'a [u8],
    pos: usize,
    nodes: Vec<NwNode>,
    /// How deep the node being read lies.
    depth: usize,
}

impl Reader<'_> {
    fn skip(&mut self) {
        loop {
            match self.text.get(self.pos) {
                Some(c) if c.is_ascii_whitespace() => self.pos += 1,
                // comments and NHX annotations in square brackets
                Some(b'[') => {
                    while self.text.get(self.pos).is_some_and(|&c| c != b']') {
                        self.pos += 1;
                    }
                    // a comment that is never closed runs to the end of the text
                    self.pos = (self.pos + 1).min(self.text.len());
                }
                _ => break,
            }
        }
    }

    fn peek(&mut self) -> Option<u8> {
        self.skip();
        self.text.get(self.pos).copied()
    }

    fn label(&mut self) -> String {
        self.skip();
        if self.text.get(self.pos) == Some(&b'\'') {
            self.pos += 1;
            // a quote inside a quoted name is written twice
            let mut name = Vec::new();
            while let Some(&c) = self.text.get(self.pos) {
                if c != b'\'' {
                    name.push(c);
                    self.pos += 1;
                } else if self.text.get(self.pos + 1) == Some(&b'\'') {
                    name.push(b'\'');
                    self.pos += 2;
                } else {
                    break;
                }
            }
            // past the closing quote, or to the end of the text when it is never closed
            self.pos = (self.pos + 1).min(self.text.len());
            return String::from_utf8_lossy(&name).into_owned();
        }
        let start = self.pos;
        while self
            .text
            .get(self.pos)
            .is_some_and(|&c| !b"(),:;[".contains(&c) && !c.is_ascii_whitespace())
        {
            self.pos += 1;
        }
        String::from_utf8_lossy(&self.text[start..self.pos]).into_owned()
    }

    /// One node and everything below it; returns its index.
    fn node(&mut self, parent: Option<usize>) -> Result<usize, String> {
        if self.depth >= MAX_DEPTH {
            return Err(format!(
                "The tree is nested more than {MAX_DEPTH} levels deep (character {}).",
                self.pos + 1
            ));
        }
        self.depth += 1;
        let read = self.node_inside(parent);
        self.depth -= 1;
        read
    }

    fn node_inside(&mut self, parent: Option<usize>) -> Result<usize, String> {
        let me = self.nodes.len();
        self.nodes.push(NwNode {
            name: String::new(),
            dist: f64::NAN,
            parent,
            children: Vec::new(),
        });
        if self.peek() == Some(b'(') {
            self.pos += 1;
            loop {
                let child = self.node(Some(me))?;
                self.nodes[me].children.push(child);
                match self.peek() {
                    Some(b',') => self.pos += 1,
                    Some(b')') => {
                        self.pos += 1;
                        break;
                    }
                    _ => {
                        return Err(format!(
                            "Expected ',' or ')' at character {}.",
                            self.pos + 1
                        ))
                    }
                }
            }
        }
        let name = self.label();
        self.nodes[me].name = name;
        if self.peek() == Some(b':') {
            self.pos += 1;
            self.skip();
            let start = self.pos;
            while self
                .text
                .get(self.pos)
                .is_some_and(|c| c.is_ascii_digit() || b"+-.eE".contains(c))
            {
                self.pos += 1;
            }
            let number = String::from_utf8_lossy(&self.text[start..self.pos]).into_owned();
            self.nodes[me].dist = number.parse().map_err(|_| {
                format!(
                    "'{number}' is not a branch length (character {}).",
                    start + 1
                )
            })?;
        }
        Ok(me)
    }
}

impl NwTree {
    /// Every tree in `text` (each ends with `;`). A missing branch length is 0 for the root and
    /// 1 for any other node, as ETE reads it.
    pub fn parse_all(text: &str) -> Result<Vec<NwTree>, String> {
        let mut reader = Reader {
            text: text.as_bytes(),
            pos: 0,
            nodes: Vec::new(),
            depth: 0,
        };
        let mut trees = Vec::new();
        while reader.peek().is_some() {
            reader.nodes.clear();
            let root = reader.node(None)?;
            if reader.peek() == Some(b';') {
                reader.pos += 1;
            } else if reader.peek().is_some() {
                return Err(format!("Expected ';' at character {}.", reader.pos + 1));
            }
            let mut nodes = std::mem::take(&mut reader.nodes);
            for (i, n) in nodes.iter_mut().enumerate() {
                if n.dist.is_nan() {
                    n.dist = if i == root { 0.0 } else { 1.0 };
                }
            }
            trees.push(NwTree { nodes, root });
        }
        if trees.is_empty() {
            return Err("No tree found: a Newick tree looks like (A:1,B:2);".into());
        }
        Ok(trees)
    }

    /// The tree as the phylogeny code draws and fits it: branch lengths, and an inner node's label as its
    /// support when it is a number.
    pub fn to_node(&self) -> super::tree::Node {
        fn below(t: &NwTree, i: usize) -> super::tree::Node {
            let n = &t.nodes[i];
            let mut node =
                super::tree::Node::inner(n.children.iter().map(|&c| below(t, c)).collect());
            if n.children.is_empty() {
                node.name.clone_from(&n.name);
            } else {
                node.support = n.name.parse().ok();
            }
            node.length = n.dist;
            node
        }
        below(self, self.root)
    }

    pub fn is_leaf(&self, i: usize) -> bool {
        self.nodes[i].children.is_empty()
    }

    pub fn leaf_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.children.is_empty()).count()
    }

    /// Nodes root first, then level by level, children in order (ETE's "levelorder").
    pub fn levelorder(&self) -> Vec<usize> {
        let mut out = Vec::with_capacity(self.nodes.len());
        let mut queue = VecDeque::from([self.root]);
        while let Some(i) = queue.pop_front() {
            out.push(i);
            queue.extend(self.nodes[i].children.iter().copied());
        }
        out
    }

    /// Children before their parent (ETE's "postorder").
    pub fn postorder(&self) -> Vec<usize> {
        let mut out = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![(self.root, 0usize)];
        while let Some((i, next)) = stack.pop() {
            match self.nodes[i].children.get(next) {
                Some(&child) => {
                    stack.push((i, next + 1));
                    stack.push((child, 0));
                }
                None => out.push(i),
            }
        }
        out
    }

    /// Newick text with a closing `;` and no length on the root. Characters Newick uses are
    /// replaced by `_` in names.
    pub fn newick(&self) -> String {
        format!("{};", self.write(self.root, false))
    }

    fn write(&self, i: usize, with_length: bool) -> String {
        let n = &self.nodes[i];
        let inside = if n.children.is_empty() {
            String::new()
        } else {
            let kids: Vec<String> = n.children.iter().map(|&c| self.write(c, true)).collect();
            format!("({})", kids.join(","))
        };
        let name: String = n
            .name
            .chars()
            .map(|c| {
                if c.is_whitespace() || "()[]:;,'".contains(c) {
                    '_'
                } else {
                    c
                }
            })
            .collect();
        if with_length {
            format!("{inside}{name}:{}", number(n.dist))
        } else {
            format!("{inside}{name}")
        }
    }
}

/// A length for text output: twelve decimals at most, without the noise a sum of floats carries
/// (1.5330000000000006 is written 1.533), and an exponent for very small values.
pub fn number(x: f64) -> String {
    if x == 0.0 || !x.is_finite() {
        return "0".into();
    }
    if x.abs() < 1e-6 {
        return format!("{x:e}");
    }
    let text = format!("{x:.12}");
    text.trim_end_matches('0').trim_end_matches('.').to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_tree_nested_far_deeper_than_any_real_one_is_refused_instead_of_overflowing_the_stack() {
        // a ladder of exactly MAX_DEPTH levels is read, one level more is not
        let ladder =
            |levels: usize| format!("{}A{};", "(".repeat(levels - 1), ",B)".repeat(levels - 1));
        let tree = &NwTree::parse_all(&ladder(MAX_DEPTH)).expect("a ladder")[0];
        assert_eq!(tree.leaf_count(), MAX_DEPTH);
        let why = NwTree::parse_all(&ladder(MAX_DEPTH + 1)).expect_err("too deep");
        assert!(why.contains("nested more than 2000 levels"), "{why}");
        // text that would have overflowed the stack answers at once, without recursing that far
        let hostile = "(".repeat(1_000_000);
        assert!(NwTree::parse_all(&hostile).is_err());
        // the depth counter is given back, so a long run of shallow trees is fine
        let many = "((A,B),(C,D));".repeat(5_000);
        assert_eq!(NwTree::parse_all(&many).expect("many trees").len(), 5_000);
    }

    #[test]
    fn text_that_stops_inside_a_comment_or_a_quote_is_refused_or_read_but_never_panics() {
        // the end of the text came before the closing bracket or quote
        for text in [
            "((A,B),(C,D))[",
            "(A,[",
            "((A,B),(C,D));\n((A,B),(C,D))[unterminated",
            "'",
            "(A,'b",
            "(A,B)'",
            "(A,B)[x",
            "[",
        ] {
            // either answer is fine; what matters is that the reader returns
            let _ = NwTree::parse_all(text);
        }
        // a comment after the last tree that is cut off leaves the tree before it as it was
        let trees = NwTree::parse_all("((A,B),(C,D));[cut").expect("a tree");
        assert_eq!((trees.len(), trees[0].leaf_count()), (1, 4));
    }

    #[test]
    fn trees_read_with_names_lengths_comments_and_quotes() {
        let t = &NwTree::parse_all("((A:1,'B c':2.5e-1)90:3,C:4)[&&NHX:x=1];").unwrap()[0];
        assert_eq!((t.nodes.len(), t.leaf_count()), (5, 3));
        let names: Vec<&str> = t
            .levelorder()
            .iter()
            .map(|&i| t.nodes[i].name.as_str())
            .collect();
        assert_eq!(names, ["", "90", "C", "A", "B c"]);
        let b = t.nodes.iter().find(|n| n.name == "B c").unwrap();
        assert!((b.dist - 0.25).abs() < 1e-12);
        assert_eq!(t.nodes[t.root].dist, 0.0); // no root length given
        let post: Vec<&str> = t
            .postorder()
            .iter()
            .map(|&i| t.nodes[i].name.as_str())
            .collect();
        assert_eq!(post, ["A", "B c", "90", "C", ""]);
    }

    #[test]
    fn several_trees_and_errors() {
        assert_eq!(NwTree::parse_all("(A:1,B:1);(C:2,D:2);").unwrap().len(), 2);
        assert!(NwTree::parse_all("(A:1,B:1").is_err());
        assert!(NwTree::parse_all("(A:x,B:1);").is_err());
        assert!(NwTree::parse_all("   ").is_err());
        let again =
            NwTree::parse_all(&NwTree::parse_all("((A:1,B:2):3,C:4);").unwrap()[0].newick())
                .unwrap();
        assert_eq!(again[0].leaf_count(), 3);
    }
}
