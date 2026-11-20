//! Where the parts of a tree go on a page: rectangular phylograms and cladograms and circular trees, as plain
//! coordinates, lines and arcs that the viewer and the SVG export both draw.

use super::tree::Node;
use std::f32::consts::TAU;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Layout {
    /// Branch lengths to scale, root at the left, tips in a column to the right of their own distance.
    Phylogram,
    /// Every tip in one column, whatever the lengths: only the branching shows.
    Cladogram,
    /// Branch lengths to scale, root in the middle, tips round the edge.
    Circular,
    /// Branch lengths to scale, no root: the equal-angle layout of an unrooted tree, branches leaving every node in a fan.
    Radial,
}

#[derive(Clone, Debug)]
pub struct Placed {
    /// The node's point in layout units. Rectangular layouts: x grows away from the root and y runs down
    /// the tips, one unit apart. Circular layout: the cartesian point round the root at the origin.
    pub x: f32,
    pub y: f32,
    /// Distance from the root, and the direction in turns of a circle (the circular layout only).
    pub radius: f32,
    pub angle: f32,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    pub name: String,
    pub support: Option<f64>,
    /// Tips at or below this node.
    pub tips: usize,
    /// Where the node is in the tree it was laid out from, counted in pre-order from the root as 0 (a clade that is collapsed counts all its nodes).
    pub origin: usize,
    /// Set for a node that stands for a clade collapsed into a triangle.
    pub cartoon: Option<Cartoon>,
}

/// A collapsed clade: the triangle reaches `far` (a distance from the root, in the layout's units) and is `rows` rows tall.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cartoon {
    pub far: f32,
    pub rows: f32,
}

impl Placed {
    pub fn is_tip(&self) -> bool {
        self.children.is_empty()
    }
}

/// One stroke of the drawing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Piece {
    Line {
        from: (f32, f32),
        to: (f32, f32),
    },
    /// An arc of the circle of `radius` round the root, from one angle to another (radians).
    Arc {
        radius: f32,
        from: f32,
        to: f32,
    },
}

pub struct Laid {
    pub kind: Layout,
    /// Root first, then every node before its children (pre-order); children in the tree's order.
    pub nodes: Vec<Placed>,
    /// The rows that the tips take together (a collapsed clade takes several): the angle of a circular layout is shared over them.
    pub rows_total: f32,
}

fn finite(length: f64) -> f32 {
    if length.is_finite() && length > 0.0 {
        length as f32
    } else {
        0.0
    }
}

/// How many rows the triangle of a collapsed clade of `tips` tips stands for: more tips, a taller triangle, up to a limit.
fn rows_for(tips: usize) -> f32 {
    (tips as f32).sqrt().ceil().clamp(1.0, MOST_CARTOON_ROWS)
}

/// The most rows a collapsed clade takes.
const MOST_CARTOON_ROWS: f32 = 8.0;

fn node_count(root: &Node) -> usize {
    let mut n = 0;
    let mut stack = vec![root];
    while let Some(node) = stack.pop() {
        n += 1;
        stack.extend(node.children.iter());
    }
    n
}

pub fn lay_out(root: &Node, kind: Layout) -> Laid {
    lay_out_with(root, kind, &[])
}

/// `lay_out`, with the nodes whose `origin` (their place in the tree's pre-order) is in `collapsed` drawn as one triangle each instead of
/// the clade below them (FigTree's cartoon). Such a node is a tip of the layout that stands for all the tips below it.
pub fn lay_out_with(root: &Node, kind: Layout, collapsed: &[usize]) -> Laid {
    let mut nodes: Vec<Placed> = Vec::new();
    let mut lengths: Vec<f32> = Vec::new();
    // for a collapsed node: the longest distance from it to a tip below
    let mut deepest: Vec<f32> = Vec::new();
    let mut folded: Vec<bool> = Vec::new();
    let mut counter = 0usize;
    let mut stack: Vec<(&Node, Option<usize>)> = vec![(root, None)];
    while let Some((node, parent)) = stack.pop() {
        let i = nodes.len();
        let origin = counter;
        counter += 1;
        let fold = !node.is_leaf() && collapsed.contains(&origin);
        nodes.push(Placed {
            x: 0.0,
            y: 0.0,
            radius: 0.0,
            angle: 0.0,
            parent,
            children: Vec::new(),
            name: node.name.clone(),
            support: node.support,
            tips: if fold { node.leaf_count() } else { 0 },
            origin,
            cartoon: None,
        });
        lengths.push(finite(node.length));
        deepest.push(if fold {
            finite(node.depth() - node.length)
        } else {
            0.0
        });
        folded.push(fold);
        if let Some(p) = parent {
            nodes[p].children.push(i);
        }
        if fold {
            // the nodes below are not drawn, but they are counted, so that the numbers of the nodes after this one stay the tree's own
            counter += node_count(node) - 1;
        } else {
            stack.extend(node.children.iter().rev().map(|c| (c, Some(i))));
        }
    }
    let n = nodes.len();
    // children come after their parents, so going backwards sees every child first
    let mut height = vec![0usize; n];
    for i in (0..n).rev() {
        if folded[i] {
            continue;
        }
        if nodes[i].children.is_empty() {
            nodes[i].tips = 1;
        } else {
            nodes[i].tips = nodes[i].children.iter().map(|&c| nodes[c].tips).sum();
            height[i] = 1 + nodes[i]
                .children
                .iter()
                .map(|&c| height[c])
                .max()
                .unwrap_or(0);
        }
    }
    let mut distance = vec![0.0f32; n];
    for i in 1..n {
        let parent = nodes[i].parent.unwrap_or(0);
        distance[i] = match kind {
            Layout::Cladogram => (height[0] - height[i]) as f32,
            _ => distance[parent] + lengths[i],
        };
    }
    // the rows a tip takes: one, and for a collapsed clade as many as its triangle is tall
    let weight: Vec<f32> = (0..n)
        .map(|i| {
            if folded[i] {
                rows_for(nodes[i].tips)
            } else {
                1.0
            }
        })
        .collect();
    let mut row = vec![0.0f32; n];
    let mut next = 0.0;
    for i in 0..n {
        if nodes[i].children.is_empty() {
            row[i] = next + (weight[i] - 1.0) / 2.0;
            next += weight[i];
        }
    }
    for i in (0..n).rev() {
        if let (Some(&first), Some(&last)) = (nodes[i].children.first(), nodes[i].children.last()) {
            let rows = nodes[i].children.iter().map(|&c| row[c]);
            let (lo, hi) = rows.fold((row[first], row[last]), |(lo, hi), r| {
                (lo.min(r), hi.max(r))
            });
            row[i] = (lo + hi) / 2.0;
        }
    }
    let rows_total = next.max(1.0);
    for i in 0..n {
        let turn = (row[i] + 0.5) / rows_total;
        let p = &mut nodes[i];
        p.radius = distance[i];
        p.angle = turn;
        (p.x, p.y) = match kind {
            Layout::Circular => (
                distance[i] * (TAU * turn).cos(),
                distance[i] * (TAU * turn).sin(),
            ),
            _ => (distance[i], row[i]),
        };
        if folded[i] {
            let far = if kind == Layout::Cladogram {
                height[0] as f32
            } else {
                distance[i] + deepest[i]
            };
            p.cartoon = Some(Cartoon {
                far,
                rows: weight[i],
            });
        }
    }
    if kind == Layout::Radial {
        equal_angle(&mut nodes, &lengths, &weight);
    }
    Laid {
        kind,
        nodes,
        rows_total,
    }
}

/// The unrooted layout of Felsenstein's equal-angle method: every node's branch leaves its parent in the middle of a share of the full
/// circle that is as large as the share of the tips that lie below the node; the root is at the centre.
fn equal_angle(nodes: &mut [Placed], lengths: &[f32], weight: &[f32]) {
    let n = nodes.len();
    let mut below = vec![0.0f32; n];
    for i in (0..n).rev() {
        below[i] = if nodes[i].children.is_empty() {
            weight[i]
        } else {
            nodes[i].children.iter().map(|&c| below[c]).sum()
        };
    }
    // the wedge of each node as turns of a circle: where it starts and where it ends
    let mut wedge = vec![(0.0f32, 1.0f32); n];
    for i in 0..n {
        let (start, end) = wedge[i];
        let mut at = start;
        for &c in &nodes[i].children.clone() {
            let share = (end - start) * below[c] / below[i].max(f32::EPSILON);
            wedge[c] = (at, at + share);
            at += share;
        }
        let middle = (start + end) / 2.0;
        nodes[i].angle = middle;
        let (px, py) = nodes[i]
            .parent
            .map_or((0.0, 0.0), |p| (nodes[p].x, nodes[p].y));
        (nodes[i].x, nodes[i].y) = if nodes[i].parent.is_none() {
            (0.0, 0.0)
        } else {
            (
                px + lengths[i] * (TAU * middle).cos(),
                py + lengths[i] * (TAU * middle).sin(),
            )
        };
    }
}

impl Laid {
    /// The tips at or below node `i`, in order.
    pub fn tips_under(&self, i: usize) -> Vec<usize> {
        let mut found = Vec::new();
        let mut stack = vec![i];
        while let Some(k) = stack.pop() {
            match self.nodes.get(k) {
                Some(n) if n.is_tip() => found.push(k),
                Some(n) => stack.extend(n.children.iter().rev()),
                None => {}
            }
        }
        found
    }

    /// Whether node `i` is `ancestor` or lies below it.
    pub fn is_within(&self, mut i: usize, ancestor: usize) -> bool {
        loop {
            if i == ancestor {
                return true;
            }
            match self.nodes.get(i).and_then(|n| n.parent) {
                Some(p) => i = p,
                None => return false,
            }
        }
    }

    /// Every stroke, each with the node it belongs to: the branch that arrives at a node, and for an inner node
    /// the line or arc that joins its children.
    pub fn pieces(&self) -> Vec<(usize, Piece)> {
        let mut out = Vec::with_capacity(2 * self.nodes.len());
        for (i, node) in self.nodes.iter().enumerate() {
            if let Some(p) = node.parent {
                let from = &self.nodes[p];
                out.push((
                    i,
                    match self.kind {
                        Layout::Radial => Piece::Line {
                            from: (from.x, from.y),
                            to: (node.x, node.y),
                        },
                        Layout::Circular => Piece::Line {
                            from: (
                                from.radius * (TAU * node.angle).cos(),
                                from.radius * (TAU * node.angle).sin(),
                            ),
                            to: (node.x, node.y),
                        },
                        _ => Piece::Line {
                            from: (from.x, node.y),
                            to: (node.x, node.y),
                        },
                    },
                ));
            }
            if self.kind != Layout::Radial
                && (node.children.len() > 1
                    || (node.children.len() == 1 && self.kind != Layout::Circular))
            {
                let first = &self.nodes[node.children[0]];
                let last = &self.nodes[*node.children.last().unwrap_or(&node.children[0])];
                out.push((
                    i,
                    match self.kind {
                        Layout::Circular => Piece::Arc {
                            radius: node.radius,
                            from: TAU * first.angle,
                            to: TAU * last.angle,
                        },
                        _ => Piece::Line {
                            from: (node.x, first.y),
                            to: (node.x, last.y),
                        },
                    },
                ));
            }
        }
        out
    }

    /// The smallest and largest x and y of the drawing, labels not counted.
    pub fn extent(&self) -> ((f32, f32), (f32, f32)) {
        let reach = |n: &Placed| n.cartoon.map_or(n.radius, |c| c.far.max(n.radius));
        let longest = self.nodes.iter().map(reach).fold(0.0, f32::max);
        match self.kind {
            Layout::Circular => ((-longest, -longest), (longest, longest)),
            Layout::Radial => {
                let (mut lo, mut hi) = ((0.0f32, 0.0f32), (0.0f32, 0.0f32));
                for n in &self.nodes {
                    let far = n.cartoon.map_or(0.0, |c| c.far - n.radius);
                    let (dx, dy) = ((TAU * n.angle).cos() * far, (TAU * n.angle).sin() * far);
                    lo = (lo.0.min(n.x + dx.min(0.0)), lo.1.min(n.y + dy.min(0.0)));
                    hi = (hi.0.max(n.x + dx.max(0.0)), hi.1.max(n.y + dy.max(0.0)));
                }
                (lo, hi)
            }
            _ => ((0.0, 0.0), (longest, (self.rows_total - 1.0).max(0.0))),
        }
    }

    /// A length for a scale bar: 1, 2 or 5 times a power of ten, at most `limit`.
    pub fn nice_length(limit: f32) -> f32 {
        if !(limit.is_finite() && limit > 0.0) {
            return 0.0;
        }
        let power = 10f32.powf(limit.log10().floor());
        [5.0, 2.0, 1.0]
            .into_iter()
            .map(|m| m * power)
            .find(|&l| l <= limit)
            .unwrap_or(power)
    }
}

/// The names of the tips, in order.
#[cfg(test)]
pub fn tip_names(laid: &Laid) -> Vec<&str> {
    laid.nodes
        .iter()
        .filter(|n| n.is_tip())
        .map(|n| n.name.as_str())
        .collect()
}

fn escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// What an SVG picture of the tree shows.
#[derive(Clone, Copy)]
pub struct SvgOptions {
    pub labels: bool,
    pub support: bool,
    pub scale_bar: bool,
}

const ROW: f32 = 18.0;
const MARGIN: f32 = 20.0;
const LABEL_SPACE: f32 = 190.0;
const WIDTH: f32 = 900.0;
const CIRCLE_SIZE: f32 = 920.0;
const FONT: f32 = 12.0;

/// The tree as an SVG picture: strokes in black, tip names beside the tips, support values small, and for the
/// layouts that keep branch lengths a scale bar.
pub fn to_svg(laid: &Laid, o: SvgOptions) -> String {
    let (min, max) = laid.extent();
    // the circular and the radial layout are drawn round a middle point, with one scale in both directions
    let circular = matches!(laid.kind, Layout::Circular | Layout::Radial);
    let centre = ((min.0 + max.0) / 2.0, (min.1 + max.1) / 2.0);
    let (width, height) = if circular {
        (CIRCLE_SIZE, CIRCLE_SIZE)
    } else {
        (
            WIDTH,
            2.0 * MARGIN + ROW * (max.1 + 1.0) + if o.scale_bar { 30.0 } else { 0.0 },
        )
    };
    let span = if circular {
        (max.0 - min.0).max(max.1 - min.1).max(f32::EPSILON)
    } else {
        max.0.max(f32::EPSILON)
    };
    let scale = if circular {
        (CIRCLE_SIZE - 2.0 * MARGIN - 2.0 * LABEL_SPACE * 0.6) / span
    } else {
        (WIDTH - 2.0 * MARGIN - LABEL_SPACE) / span
    };
    let at = |x: f32, y: f32| -> (f32, f32) {
        if circular {
            (
                CIRCLE_SIZE / 2.0 + (x - centre.0) * scale,
                CIRCLE_SIZE / 2.0 + (y - centre.1) * scale,
            )
        } else {
            (MARGIN + x * scale, MARGIN + ROW / 2.0 + y * ROW)
        }
    };
    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width:.0}\" height=\"{height:.0}\" viewBox=\"0 0 {width:.0} {height:.0}\">\n\
         <rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n\
         <g stroke=\"black\" stroke-width=\"1.4\" fill=\"none\" stroke-linecap=\"round\">\n"
    );
    for (_, piece) in laid.pieces() {
        match piece {
            Piece::Line { from, to } => {
                let ((x1, y1), (x2, y2)) = (at(from.0, from.1), at(to.0, to.1));
                svg.push_str(&format!(
                    "<line x1=\"{x1:.1}\" y1=\"{y1:.1}\" x2=\"{x2:.1}\" y2=\"{y2:.1}\"/>\n"
                ));
            }
            Piece::Arc { radius, from, to } => {
                let r = radius * scale;
                let ((x1, y1), (x2, y2)) = (
                    at(radius * from.cos(), radius * from.sin()),
                    at(radius * to.cos(), radius * to.sin()),
                );
                let large = u8::from((to - from).abs() > std::f32::consts::PI);
                let sweep = u8::from(to >= from);
                svg.push_str(&format!(
                    "<path d=\"M {x1:.1} {y1:.1} A {r:.1} {r:.1} 0 {large} {sweep} {x2:.1} {y2:.1}\"/>\n"
                ));
            }
        }
    }
    svg.push_str("</g>\n<g font-family=\"Helvetica, Arial, sans-serif\">\n");
    for node in &laid.nodes {
        let (x, y) = at(node.x, node.y);
        if node.is_tip() && o.labels {
            let side = if circular && node.x < 0.0 {
                ("end", -6.0)
            } else {
                ("start", 6.0)
            };
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{y:.1}\" font-size=\"{FONT}\" text-anchor=\"{}\" dominant-baseline=\"middle\">{}</text>\n",
                x + side.1,
                side.0,
                escape(&node.name)
            ));
        } else if let (Some(s), true) = (node.support, o.support && !node.is_tip()) {
            svg.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" font-size=\"9\" text-anchor=\"end\" fill=\"#666\">{s:.0}</text>\n",
                x - 3.0,
                y - 3.0
            ));
        }
    }
    if o.scale_bar && laid.kind != Layout::Circular {
        let length = Laid::nice_length(max.0 * 0.25);
        if length > 0.0 {
            let (x1, y) = (MARGIN, height - 14.0);
            let x2 = x1 + length * scale;
            svg.push_str(&format!(
                "<line x1=\"{x1:.1}\" y1=\"{y:.1}\" x2=\"{x2:.1}\" y2=\"{y:.1}\" stroke=\"black\" stroke-width=\"1.4\"/>\n\
                 <text x=\"{:.1}\" y=\"{:.1}\" font-size=\"10\" text-anchor=\"middle\">{length}</text>\n",
                (x1 + x2) / 2.0,
                y - 4.0
            ));
        }
    }
    svg.push_str("</g>\n</svg>\n");
    svg
}

#[cfg(test)]
mod folded_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::newick::NwTree;

    /// ((A:1,B:2)X:1,C:3): X joins A and B.
    fn small() -> Node {
        NwTree::parse_all("((A:1,B:2)90:1,C:3);")
            .expect("a tree")
            .remove(0)
            .to_node()
    }

    #[test]
    fn a_phylogram_puts_nodes_at_their_distance_and_tips_one_row_apart() {
        let laid = lay_out(&small(), Layout::Phylogram);
        let by_name = |name: &str| laid.nodes.iter().find(|n| n.name == name).expect("a tip");
        let (a, b, c) = (by_name("A"), by_name("B"), by_name("C"));
        assert_eq!((a.x, b.x, c.x), (2.0, 3.0, 3.0));
        assert_eq!((a.y, b.y, c.y), (0.0, 1.0, 2.0));
        // the inner node X is 1 from the root and halfway between its tips; the root halfway between its children
        let x = &laid.nodes[1];
        assert_eq!((x.x, x.y, x.tips), (1.0, 0.5, 2));
        assert_eq!(
            (laid.nodes[0].x, laid.nodes[0].y, laid.nodes[0].tips),
            (0.0, 1.25, 3)
        );
    }

    #[test]
    fn a_cladogram_lines_every_tip_up_whatever_the_lengths() {
        let laid = lay_out(&small(), Layout::Cladogram);
        let tips: Vec<f32> = laid
            .nodes
            .iter()
            .filter(|n| n.is_tip())
            .map(|n| n.x)
            .collect();
        assert_eq!(tips, vec![2.0, 2.0, 2.0]);
        assert_eq!(laid.nodes[1].x, 1.0);
        assert_eq!(laid.nodes[0].x, 0.0);
    }

    #[test]
    fn a_circular_tree_spaces_the_tips_round_the_circle_at_their_distance() {
        let laid = lay_out(&small(), Layout::Circular);
        let tips: Vec<&Placed> = laid.nodes.iter().filter(|n| n.is_tip()).collect();
        let turns: Vec<f32> = tips.iter().map(|t| t.angle).collect();
        let step = turns[1] - turns[0];
        assert!((step - 1.0 / 3.0).abs() < 1e-6 && (turns[2] - turns[1] - step).abs() < 1e-6);
        for t in &tips {
            let r = (t.x * t.x + t.y * t.y).sqrt();
            assert!((r - t.radius).abs() < 1e-5);
        }
        assert_eq!(tips[0].radius, 2.0);
        let ((x0, y0), (x1, y1)) = laid.extent();
        assert_eq!((x0, y0, x1, y1), (-3.0, -3.0, 3.0, 3.0));
    }

    #[test]
    fn every_branch_and_every_join_is_drawn_once() {
        // four branches and two joins, in either style
        for kind in [Layout::Phylogram, Layout::Cladogram, Layout::Circular] {
            let laid = lay_out(&small(), kind);
            let pieces = laid.pieces();
            assert_eq!(pieces.len(), 6, "{kind:?}");
            let arcs = pieces
                .iter()
                .filter(|(_, p)| matches!(p, Piece::Arc { .. }))
                .count();
            assert_eq!(arcs, if kind == Layout::Circular { 2 } else { 0 });
        }
    }

    #[test]
    fn a_clade_is_its_tips_and_its_inner_nodes() {
        let laid = lay_out(&small(), Layout::Phylogram);
        let names: Vec<&str> = laid
            .tips_under(1)
            .iter()
            .map(|&i| laid.nodes[i].name.as_str())
            .collect();
        assert_eq!(names, ["A", "B"]);
        assert_eq!(laid.tips_under(0).len(), 3);
        assert!(laid.is_within(2, 1) && laid.is_within(1, 1) && !laid.is_within(4, 1));
        assert_eq!(tip_names(&laid), ["A", "B", "C"]);
    }

    #[test]
    fn odd_trees_do_not_fail() {
        let single = Node {
            name: "only".into(),
            length: 1.0,
            support: None,
            children: Vec::new(),
        };
        for kind in [Layout::Phylogram, Layout::Cladogram, Layout::Circular] {
            let laid = lay_out(&single, kind);
            assert_eq!(laid.nodes.len(), 1);
            assert!(laid.pieces().is_empty());
            assert!(laid.nodes[0].x.is_finite() && laid.nodes[0].y.is_finite());
        }
        // a negative or missing length counts as none
        let tip = |name: &str, length: f64| Node {
            name: name.into(),
            length,
            support: None,
            children: Vec::new(),
        };
        let odd = Node::inner(vec![
            tip("A", -1.0),
            tip("B", f64::NAN),
            tip("C", f64::INFINITY),
        ]);
        let laid = lay_out(&odd, Layout::Phylogram);
        assert!(laid.nodes.iter().all(|n| n.x.is_finite() && n.x >= 0.0));
    }

    #[test]
    fn scale_bars_use_round_lengths() {
        assert_eq!(Laid::nice_length(0.37), 0.2);
        assert_eq!(Laid::nice_length(7.0), 5.0);
        assert_eq!(Laid::nice_length(1.0), 1.0);
        assert_eq!(Laid::nice_length(0.0), 0.0);
        assert_eq!(Laid::nice_length(f32::NAN), 0.0);
    }

    #[test]
    fn the_svg_has_the_strokes_the_names_and_a_scale_bar() {
        let laid = lay_out(&small(), Layout::Phylogram);
        let all = SvgOptions {
            labels: true,
            support: true,
            scale_bar: true,
        };
        let svg = to_svg(&laid, all);
        assert!(svg.starts_with("<svg") && svg.trim_end().ends_with("</svg>"));
        assert_eq!(
            svg.matches("<line").count(),
            6 + 1,
            "six strokes and the scale bar"
        );
        for name in ["A", "B", "C"] {
            assert!(svg.contains(&format!(">{name}</text>")), "{name}");
        }
        assert!(svg.contains(">90</text>"));
        let bare = to_svg(
            &laid,
            SvgOptions {
                labels: false,
                support: false,
                scale_bar: false,
            },
        );
        assert!(!bare.contains("<text"));
        let circle = to_svg(&lay_out(&small(), Layout::Circular), all);
        assert_eq!(circle.matches("<path").count(), 2);
    }

    #[test]
    fn names_are_escaped_in_svg() {
        let odd = Node {
            name: String::new(),
            length: 0.0,
            support: None,
            children: vec![
                Node {
                    name: "a<b & \"c\"".into(),
                    length: 1.0,
                    support: None,
                    children: Vec::new(),
                },
                Node {
                    name: "d".into(),
                    length: 1.0,
                    support: None,
                    children: Vec::new(),
                },
            ],
        };
        let svg = to_svg(
            &lay_out(&odd, Layout::Phylogram),
            SvgOptions {
                labels: true,
                support: false,
                scale_bar: false,
            },
        );
        assert!(svg.contains("a&lt;b &amp; &quot;c&quot;"));
    }
}
