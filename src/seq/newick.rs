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
    /// What the file says in `[&key=value,...]` after the node (the way BEAST's TreeAnnotator and FigTree write it): the keys as
    /// written and the values as written, a list in braces (`{1.2,3.4}`) still a list.
    pub attrs: Vec<(String, String)>,
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
    /// The node that a `[&...]` comment read now belongs to: the one being read, or the last one that was finished.
    current: usize,
}

/// The `key=value` pairs of a comment that begins with `&` (`&posterior=0.98,height_95%_HPD={1.2,3.4},!color=#ff0000`), cut at the commas that
/// are neither inside braces nor inside quotes. Anything else (an ordinary comment, `&R`, `&U`, NHX's `&&NHX:`) has none.
fn annotations(comment: &str) -> Vec<(String, String)> {
    let Some(body) = comment.strip_prefix('&').filter(|b| !b.starts_with('&')) else {
        return Vec::new();
    };
    let (mut parts, mut depth, mut quoted, mut from) = (Vec::new(), 0i32, false, 0);
    for (i, c) in body.char_indices() {
        match c {
            '"' => quoted = !quoted,
            '{' if !quoted => depth += 1,
            '}' if !quoted => depth -= 1,
            ',' if !quoted && depth <= 0 => {
                parts.push(&body[from..i]);
                from = i + 1;
            }
            _ => {}
        }
    }
    parts.push(&body[from..]);
    parts
        .into_iter()
        .filter_map(|p| p.split_once('='))
        .map(|(k, v)| (k.trim().to_string(), v.trim().trim_matches('"').to_string()))
        .filter(|(k, _)| !k.is_empty())
        .collect()
}

impl Reader<'_> {
    fn skip(&mut self) {
        loop {
            match self.text.get(self.pos) {
                Some(c) if c.is_ascii_whitespace() => self.pos += 1,
                // comments and NHX annotations in square brackets
                Some(b'[') => {
                    let start = self.pos + 1;
                    let mut quoted = false;
                    while let Some(&c) = self.text.get(self.pos) {
                        if c == b'"' {
                            quoted = !quoted;
                        }
                        if c == b']' && !quoted {
                            break;
                        }
                        self.pos += 1;
                    }
                    let inside = String::from_utf8_lossy(&self.text[start.min(self.pos)..self.pos])
                        .into_owned();
                    if let Some(node) = self.nodes.get_mut(self.current) {
                        node.attrs.extend(annotations(&inside));
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
            attrs: Vec::new(),
        });
        self.current = me;
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
            // a comment after the closing bracket is this node's, not its last child's
            self.current = me;
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
        NwTree::parse_trees(text, true)
    }

    /// `parse_all`, but a branch length that is not written stays NaN, so that a tree with no lengths can be told from one whose
    /// lengths are all 1.
    pub fn parse_all_raw(text: &str) -> Result<Vec<NwTree>, String> {
        NwTree::parse_trees(text, false)
    }

    fn parse_trees(text: &str, fill_missing: bool) -> Result<Vec<NwTree>, String> {
        let mut reader = Reader {
            text: text.as_bytes(),
            pos: 0,
            nodes: Vec::new(),
            depth: 0,
            current: usize::MAX,
        };
        let mut trees = Vec::new();
        while reader.peek().is_some() {
            reader.nodes.clear();
            reader.current = usize::MAX;
            let root = reader.node(None)?;
            if reader.peek() == Some(b';') {
                reader.pos += 1;
            } else if reader.peek().is_some() {
                return Err(format!("Expected ';' at character {}.", reader.pos + 1));
            }
            let mut nodes = std::mem::take(&mut reader.nodes);
            for (i, n) in nodes.iter_mut().enumerate() {
                if fill_missing && n.dist.is_nan() {
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
        format!("{};", self.write(self.root, false, false))
    }

    /// `newick`, with each node's annotations written after its name as `[&key=value,...]`, where BEAST and FigTree put them.
    pub fn newick_annotated(&self) -> String {
        format!("{};", self.write(self.root, false, true))
    }

    /// The node `i`'s annotation `key`, as written.
    pub fn attr(&self, i: usize, key: &str) -> Option<&str> {
        self.nodes[i]
            .attrs
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// The annotation as a number (a value that is not one, or a list, is none).
    pub fn attr_number(&self, i: usize, key: &str) -> Option<f64> {
        self.attr(i, key)?
            .parse()
            .ok()
            .filter(|v: &f64| v.is_finite())
    }

    /// The annotation as a range: a list of two numbers in braces (`{1.2,3.4}`, as `height_95%_HPD` is written), the least first.
    pub fn attr_range(&self, i: usize, key: &str) -> Option<(f64, f64)> {
        let inside = self.attr(i, key)?.strip_prefix('{')?.strip_suffix('}')?;
        let (a, b) = inside.split_once(',')?;
        let (a, b) = (a.trim().parse::<f64>().ok()?, b.trim().parse::<f64>().ok()?);
        (a.is_finite() && b.is_finite()).then_some((a.min(b), a.max(b)))
    }

    /// The keys that some node has, in the order of the nodes (a parent before its children) in which they were first met.
    pub fn attr_keys(&self) -> Vec<String> {
        let mut keys: Vec<String> = Vec::new();
        for n in &self.nodes {
            for (k, _) in &n.attrs {
                if !keys.contains(k) {
                    keys.push(k.clone());
                }
            }
        }
        keys
    }

    /// The nodes root first, then every node before its children, children in order: the order that `treelayout` numbers them in.
    pub fn preorder(&self) -> Vec<usize> {
        let mut out = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![self.root];
        while let Some(i) = stack.pop() {
            out.push(i);
            stack.extend(self.nodes[i].children.iter().rev());
        }
        out
    }

    fn write(&self, i: usize, with_length: bool, annotated: bool) -> String {
        let n = &self.nodes[i];
        let inside = if n.children.is_empty() {
            String::new()
        } else {
            let kids: Vec<String> = n
                .children
                .iter()
                .map(|&c| self.write(c, true, annotated))
                .collect();
            format!("({})", kids.join(","))
        };
        let notes = if annotated && !n.attrs.is_empty() {
            let pairs: Vec<String> = n
                .attrs
                .iter()
                .map(|(k, v)| format!("{k}={}", written(v)))
                .collect();
            format!("[&{}]", pairs.join(","))
        } else {
            String::new()
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
            format!("{inside}{name}{notes}:{}", number(n.dist))
        } else {
            format!("{inside}{name}{notes}")
        }
    }
}

/// An annotation's value as it is written: text with a comma, a space or a bracket in it is quoted, a list or a number is not.
fn written(value: &str) -> String {
    let plain = value.starts_with('{')
        || value.parse::<f64>().is_ok()
        || value
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || "_.-%#".contains(c));
    if plain {
        value.to_string()
    } else {
        format!("\"{}\"", value.replace('"', "'"))
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
    fn annotations_in_square_brackets_belong_to_the_node_before_them() {
        let text = "((A[&rate=1.5]:1,B:2)[&posterior=0.98,height_95%_HPD={0.5,1.5},!color=#ff0000,note=\"a, b\"]:3,C:4)[&R];";
        let t = &NwTree::parse_all(text).unwrap()[0];
        let by = |name: &str| t.nodes.iter().position(|n| n.name == name).unwrap();
        assert_eq!(t.attr(by("A"), "rate"), Some("1.5"));
        assert!(t.nodes[by("B")].attrs.is_empty() && t.nodes[by("C")].attrs.is_empty());
        let inner = t.nodes[by("A")].parent.unwrap();
        assert_eq!(t.attr_number(inner, "posterior"), Some(0.98));
        assert_eq!(t.attr_range(inner, "height_95%_HPD"), Some((0.5, 1.5)));
        assert_eq!(t.attr(inner, "!color"), Some("#ff0000"));
        // a comma inside quotes does not cut the value
        assert_eq!(t.attr(inner, "note"), Some("a, b"));
        // `[&R]` marks the tree as rooted and is no annotation
        assert!(t.nodes[t.root].attrs.is_empty());
        assert_eq!(t.attr_number(inner, "note"), None);
        assert_eq!(t.attr_range(inner, "posterior"), None);
        assert_eq!(
            t.attr_keys(),
            ["posterior", "height_95%_HPD", "!color", "note", "rate"]
        );
        // other comments, and NHX's, are not annotations
        let plain = &NwTree::parse_all("(A[just a comment]:1,B[&&NHX:S=x]:1);").unwrap()[0];
        assert!(plain.nodes.iter().all(|n| n.attrs.is_empty()));
    }

    #[test]
    fn a_range_written_the_wrong_way_round_is_put_right_and_what_is_not_a_range_is_none() {
        let t = &NwTree::parse_all("(A[&r={3,1},s={1},t={1,x},u=1,v={1,2,3},w={1e400,2}]:1,B:1);")
            .unwrap()[0];
        let a = t.nodes.iter().position(|n| n.name == "A").unwrap();
        assert_eq!(t.attr_range(a, "r"), Some((1.0, 3.0)));
        for key in ["s", "t", "u", "v", "w", "missing"] {
            assert_eq!(t.attr_range(a, key), None, "{key}");
        }
    }

    #[test]
    fn a_tree_written_with_its_annotations_reads_back_the_same() {
        let text =
            "((A[&rate=1.5]:1,B:2)[&posterior=0.98,height_95%_HPD={0.5,1.5},note=\"a, b\"]:3,C:4);";
        let t = &NwTree::parse_all(text).unwrap()[0];
        let again = &NwTree::parse_all(&t.newick_annotated()).unwrap()[0];
        assert_eq!(again, t);
        // the plain writer leaves them out
        assert!(!t.newick().contains('['));
    }

    #[test]
    fn the_order_of_the_nodes_is_the_one_the_layout_numbers_them_in() {
        let t = &NwTree::parse_all("((A,B)x,(C,D)y)z;").unwrap()[0];
        let names: Vec<&str> = t
            .preorder()
            .iter()
            .map(|&i| t.nodes[i].name.as_str())
            .collect();
        assert_eq!(names, ["z", "x", "A", "B", "y", "C", "D"]);
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
