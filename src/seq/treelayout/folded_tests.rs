use super::*;
use crate::seq::newick::NwTree;

fn node(text: &str) -> Node {
    NwTree::parse_all(text).expect("a tree").remove(0).to_node()
}

/// root 0, x 1 (A to D 2 to 5), y 6 (E 7, F 8)
const TWO_CLADES: &str = "((A:1,B:1,C:1,D:1)x:1,(E:2,F:2)y:1);";

#[test]
fn a_collapsed_clade_is_one_tip_that_stands_for_all_of_its_tips() {
    let laid = lay_out_with(&node(TWO_CLADES), Layout::Phylogram, &[1]);
    let origins: Vec<usize> = laid.nodes.iter().map(|n| n.origin).collect();
    // the nodes below x are not there, but they are counted: y is still 6
    assert_eq!(origins, [0, 1, 6, 7, 8]);
    let x = &laid.nodes[1];
    assert!(x.is_tip() && x.tips == 4);
    let c = x.cartoon.expect("a triangle");
    // it reaches the farthest tip below it (1 + 1) and is as tall as 4 tips need: two rows
    assert_eq!((c.far, c.rows), (2.0, 2.0));
    assert!(laid.nodes[0].cartoon.is_none() && laid.nodes[2].cartoon.is_none());
    // two rows for the triangle, then E and F: its middle is half a row from the first
    assert_eq!(x.y, 0.5);
    assert_eq!((laid.nodes[3].y, laid.nodes[4].y), (2.0, 3.0));
    assert_eq!(laid.rows_total, 4.0);
    assert_eq!(laid.extent(), ((0.0, 0.0), (3.0, 3.0)));
    // the root's tips are all of them
    assert_eq!(laid.nodes[0].tips, 6);
}

#[test]
fn nothing_collapsed_gives_the_layout_it_always_gave() {
    let tree = node(TWO_CLADES);
    for kind in [Layout::Phylogram, Layout::Cladogram, Layout::Circular] {
        let plain = lay_out(&tree, kind);
        let none = lay_out_with(&tree, kind, &[]);
        assert_eq!(plain.nodes.len(), 9);
        for (a, b) in plain.nodes.iter().zip(&none.nodes) {
            assert_eq!((a.x, a.y, a.tips, a.origin), (b.x, b.y, b.tips, b.origin));
            assert!(a.cartoon.is_none());
        }
        assert_eq!(
            plain.nodes.iter().map(|n| n.origin).collect::<Vec<_>>(),
            (0..9).collect::<Vec<_>>()
        );
    }
}

#[test]
fn a_cladogram_triangle_reaches_the_column_of_the_tips_and_a_tip_cannot_be_collapsed() {
    let laid = lay_out_with(&node(TWO_CLADES), Layout::Cladogram, &[1, 2]);
    assert_eq!(laid.nodes[1].cartoon.map(|c| c.far), Some(2.0));
    // origin 2 is the tip A, which is inside x and not drawn; asking for a tip alone changes nothing
    let tip = lay_out_with(&node(TWO_CLADES), Layout::Phylogram, &[2]);
    assert_eq!(tip.nodes.len(), 9);
    assert!(tip.nodes.iter().all(|n| n.cartoon.is_none()));
}

#[test]
fn collapsing_the_root_leaves_one_triangle_and_the_circle_is_shared_over_the_rows() {
    let all = lay_out_with(&node(TWO_CLADES), Layout::Phylogram, &[0]);
    assert_eq!(all.nodes.len(), 1);
    assert_eq!(all.nodes[0].tips, 6);
    let circle = lay_out_with(&node(TWO_CLADES), Layout::Circular, &[1]);
    // x holds 2 of the 4 rows: its direction is a quarter of the circle from the start, E and F follow
    assert!(
        (circle.nodes[1].angle - 0.5 / 4.0 * 1.0 - 0.0625).abs() < 1e-6
            || (circle.nodes[1].angle - 0.25 / 1.0).abs() < 1.0
    );
    assert!((circle.nodes[1].angle - (0.5 + 0.5) / 4.0).abs() < 1e-6);
}

#[test]
fn the_radial_layout_puts_the_root_in_the_middle_and_every_branch_its_own_length_long() {
    let tree = node("((A:1,B:2):1,(C:3,D:1):2);");
    let laid = lay_out(&tree, Layout::Radial);
    assert_eq!((laid.nodes[0].x, laid.nodes[0].y), (0.0, 0.0));
    for (i, n) in laid.nodes.iter().enumerate().skip(1) {
        let p = &laid.nodes[n.parent.expect("a parent")];
        let d = ((n.x - p.x).powi(2) + (n.y - p.y).powi(2)).sqrt();
        let length = n.radius - p.radius;
        assert!((d - length).abs() < 1e-5, "node {i}: {d} against {length}");
    }
    // four equal tips share the circle in quarters of a half: the directions are 45, 135, 225 and 315 degrees
    let tips: Vec<f32> = laid
        .nodes
        .iter()
        .filter(|n| n.is_tip())
        .map(|n| n.angle)
        .collect();
    assert_eq!(tips, [0.125, 0.375, 0.625, 0.875]);
    // a branch is a line from its parent, and there are no arcs or joins
    let pieces = laid.pieces();
    assert_eq!(pieces.len(), 6);
    for (i, piece) in &pieces {
        let (n, p) = (
            &laid.nodes[*i],
            &laid.nodes[laid.nodes[*i].parent.expect("a parent")],
        );
        assert_eq!(
            *piece,
            Piece::Line {
                from: (p.x, p.y),
                to: (n.x, n.y)
            }
        );
    }
    let ((x0, y0), (x1, y1)) = laid.extent();
    assert!(laid
        .nodes
        .iter()
        .all(|n| n.x >= x0 && n.x <= x1 && n.y >= y0 && n.y <= y1));
}

#[test]
fn a_bigger_clade_gets_a_bigger_wedge_in_the_radial_layout() {
    // three tips on one side and one on the other: the wedge of the three is three times as wide
    let laid = lay_out(&node("((A:1,B:1,C:1):1,D:1);"), Layout::Radial);
    let angle = |name: &str| {
        laid.nodes
            .iter()
            .find(|n| n.name == name)
            .expect("a tip")
            .angle
    };
    // D's wedge is the last quarter: its middle is at 0.875; A, B, C share the first three quarters
    assert!((angle("D") - 0.875).abs() < 1e-6);
    assert!(
        (angle("A") - 0.125).abs() < 1e-6
            && (angle("B") - 0.375).abs() < 1e-6
            && (angle("C") - 0.625).abs() < 1e-6
    );
}

#[test]
fn a_radial_tree_is_drawn_as_svg_inside_the_picture_with_every_tip_name() {
    let laid = lay_out(&node("((A:1,B:2):1,(C:3,D:1):2);"), Layout::Radial);
    let svg = to_svg(
        &laid,
        SvgOptions {
            labels: true,
            support: false,
            scale_bar: true,
        },
    );
    for name in ["A", "B", "C", "D"] {
        assert!(svg.contains(&format!(">{name}</text>")), "{svg}");
    }
    // the branches are straight lines (no arcs), and every point is inside the picture
    assert!(!svg.contains(" A "), "{svg}");
    for chunk in svg.split("x1=\"").skip(1) {
        let x: f32 = chunk.split('"').next().unwrap().parse().unwrap();
        assert!((0.0..=CIRCLE_SIZE).contains(&x), "{x}");
    }
    // a scale bar is there: the lengths are to scale
    assert!(svg.contains("text-anchor=\"middle\""));
}
