use super::*;
use crate::seq::treelayout::{lay_out, lay_out_with, Layout};

const TEXT: &str = "((A[&loc=\"Asia\",rate=1]:1,B[&loc=\"Europe\",rate=3]:1)[&posterior=0.9,height_95%_HPD={0.5,1.5}]:1,(C[&loc=\"Asia\",rate=2]:1,D:1)[&posterior=0.5]:1)[&posterior=1];";

fn setup() -> (NwTree, Laid) {
    let t = NwTree::parse_all(TEXT).unwrap().remove(0);
    let laid = lay_out(&t.to_node(), Layout::Phylogram);
    (t, laid)
}

fn index_of(laid: &Laid, name: &str) -> usize {
    laid.nodes.iter().position(|n| n.name == name).unwrap()
}

#[test]
fn a_number_gets_a_colour_from_the_gradient_by_its_place_between_the_least_and_the_greatest() {
    let (t, laid) = setup();
    let d = decorate(
        &t,
        &laid,
        &StyleOptions {
            colour_by: Some("rate".into()),
            ..StyleOptions::default()
        },
    );
    let (a, b, c) = (
        index_of(&laid, "A"),
        index_of(&laid, "B"),
        index_of(&laid, "C"),
    );
    assert_eq!(d.colour[a], Some(gradient(0.0)));
    assert_eq!(d.colour[b], Some(gradient(1.0)));
    assert_eq!(d.colour[c], Some(gradient(0.5)));
    // D has no rate, and no node above has one: they keep the usual colour
    assert_eq!(d.colour[index_of(&laid, "D")], None);
    assert_eq!(d.colour[0], None);
    assert_eq!(
        d.legend,
        Some(Legend::Gradient {
            key: "rate".into(),
            low: 1.0,
            high: 3.0
        })
    );
    assert_eq!(gradient(0.0), GRADIENT[0]);
    assert_eq!(gradient(0.5), GRADIENT[1]);
    assert_eq!(gradient(1.0), GRADIENT[2]);
    assert_eq!(gradient(7.0), GRADIENT[2]);
}

#[test]
fn text_gets_a_colour_for_each_value_in_the_order_they_are_met() {
    let (t, laid) = setup();
    let d = decorate(
        &t,
        &laid,
        &StyleOptions {
            colour_by: Some("loc".into()),
            ..StyleOptions::default()
        },
    );
    let (a, b, c) = (
        index_of(&laid, "A"),
        index_of(&laid, "B"),
        index_of(&laid, "C"),
    );
    assert_eq!(d.colour[a], Some(PALETTE[0]));
    assert_eq!(d.colour[b], Some(PALETTE[1]));
    assert_eq!(d.colour[c], d.colour[a]);
    let Some(Legend::Categories { items, .. }) = &d.legend else {
        panic!("a legend of categories")
    };
    assert_eq!(
        items,
        &[
            ("Asia".to_string(), PALETTE[0]),
            ("Europe".to_string(), PALETTE[1])
        ]
    );
    // an annotation nobody has colours nothing and has no legend
    let none = decorate(
        &t,
        &laid,
        &StyleOptions {
            colour_by: Some("nothing".into()),
            ..StyleOptions::default()
        },
    );
    assert!(none.colour.iter().all(Option::is_none) && none.legend.is_none());
}

#[test]
fn width_follows_a_number_and_is_left_alone_for_text() {
    let (t, laid) = setup();
    let wide = decorate(
        &t,
        &laid,
        &StyleOptions {
            width_by: Some("rate".into()),
            ..StyleOptions::default()
        },
    );
    assert_eq!(wide.width[index_of(&laid, "A")], 1.0);
    assert_eq!(wide.width[index_of(&laid, "B")], MOST_WIDTH);
    assert_eq!(wide.width[index_of(&laid, "D")], 1.0);
    let text = decorate(
        &t,
        &laid,
        &StyleOptions {
            width_by: Some("loc".into()),
            ..StyleOptions::default()
        },
    );
    assert!(text.width.iter().all(|&w| w == 1.0));
}

#[test]
fn labels_join_the_annotations_asked_for_and_leave_the_tips_to_their_names() {
    let (t, laid) = setup();
    let d = decorate(
        &t,
        &laid,
        &StyleOptions {
            labels: vec!["posterior".into(), "rate".into(), "none".into()],
            ..StyleOptions::default()
        },
    );
    // the root has a posterior of 1 (whole: "1"), the inner nodes 0.9 and 0.5
    assert_eq!(d.label[0], "1");
    assert_eq!(d.label[1], "0.90");
    assert!(d.label[index_of(&laid, "A")].is_empty());
    // the node's own label stands for an empty key
    let named = NwTree::parse_all("((A,B)90,(C,D)75);").unwrap().remove(0);
    let laid2 = lay_out(&named.to_node(), Layout::Phylogram);
    let d2 = decorate(
        &named,
        &laid2,
        &StyleOptions {
            labels: vec![String::new()],
            ..StyleOptions::default()
        },
    );
    assert_eq!((d2.label[1].as_str(), d2.label[4].as_str()), ("90", "75"));
    assert_eq!(format_number(0.98765), "0.99");
    assert_eq!(format_number(12.0), "12");
    assert_eq!(format_number(0.0), "0");
    assert_eq!(format_number(0.0000123), "1.23e-5");
}

#[test]
fn a_shape_is_on_every_node_that_has_the_annotation_and_is_larger_for_a_greater_number() {
    let (t, laid) = setup();
    let d = decorate(
        &t,
        &laid,
        &StyleOptions {
            shape_by: Some("posterior".into()),
            ..StyleOptions::default()
        },
    );
    assert_eq!(d.shape[0], Some(SHAPE_RADII.1));
    assert_eq!(
        d.shape[2.min(laid.nodes.len() - 1)],
        d.shape[2.min(laid.nodes.len() - 1)]
    );
    let inner: Vec<Option<f32>> = laid
        .nodes
        .iter()
        .enumerate()
        .filter(|(_, n)| !n.is_tip())
        .map(|(i, _)| d.shape[i])
        .collect();
    assert_eq!(inner.len(), 3);
    assert!(inner.iter().all(Option::is_some));
    assert!(inner[1].unwrap() > inner[2].unwrap(), "0.9 against 0.5");
    assert!(d.shape[index_of(&laid, "A")].is_none());
    let text = decorate(
        &t,
        &laid,
        &StyleOptions {
            shape_by: Some("loc".into()),
            ..StyleOptions::default()
        },
    );
    assert_eq!(text.shape[index_of(&laid, "A")], Some(SHAPE_PLAIN));
}

#[test]
fn a_range_of_heights_is_a_bar_measured_from_the_tips_backwards() {
    let (t, laid) = setup();
    let d = decorate(
        &t,
        &laid,
        &StyleOptions {
            bars_by: Some("height_95%_HPD".into()),
            ..StyleOptions::default()
        },
    );
    // the tips are 2 from the root; a node whose height lies between 0.5 and 1.5 lies between 0.5 and 1.5 from the root
    let (lo, hi) = d.bar[1].expect("a bar");
    assert!((lo - 0.5).abs() < 1e-6 && (hi - 1.5).abs() < 1e-6);
    assert!(d
        .bar
        .iter()
        .enumerate()
        .all(|(i, b)| (i == 1) == b.is_some()));
}

#[test]
fn a_collapsed_clade_is_decorated_with_the_annotations_of_its_own_node() {
    let t = NwTree::parse_all(TEXT).unwrap().remove(0);
    // the first inner node (origin 1) is collapsed: the nodes after it keep their own numbers
    let laid = lay_out_with(&t.to_node(), Layout::Phylogram, &[1]);
    let d = decorate(
        &t,
        &laid,
        &StyleOptions {
            labels: vec!["posterior".into()],
            colour_by: Some("rate".into()),
            ..StyleOptions::default()
        },
    );
    assert_eq!(d.label[1], "0.90");
    let c = index_of(&laid, "C");
    assert_eq!(laid.nodes[c].origin, 5);
    assert_eq!(d.colour[c], Some(gradient(0.5)));
}

#[test]
fn a_tree_without_nodes_or_a_layout_of_one_node_is_decorated_without_panicking() {
    let single = NwTree::parse_all("A;").unwrap().remove(0);
    let laid = lay_out(&single.to_node(), Layout::Phylogram);
    let d = decorate(
        &single,
        &laid,
        &StyleOptions {
            colour_by: Some("x".into()),
            labels: vec!["y".into()],
            ..StyleOptions::default()
        },
    );
    assert_eq!(d.colour.len(), 1);
}

#[test]
fn the_axis_counts_from_the_root_from_the_tips_or_from_a_date() {
    let deepest = 8.0;
    let root = Axis::default();
    assert_eq!(
        (root.value(0.0, deepest), root.value(8.0, deepest)),
        (0.0, 8.0)
    );
    let before = Axis {
        mode: AxisMode::BeforeTips,
        ..Axis::default()
    };
    assert_eq!(
        (before.value(0.0, deepest), before.value(8.0, deepest)),
        (8.0, 0.0)
    );
    // 8 units of length are 80 years: the root is 80 years before the tips, which were sampled in 2015.5
    let dates = Axis {
        mode: AxisMode::Dates,
        scale: 10.0,
        tip_date: 2015.5,
        grid: false,
    };
    assert_eq!(
        (dates.value(8.0, deepest), dates.value(0.0, deepest)),
        (2015.5, 1935.5)
    );
    for axis in [root, before, dates] {
        for v in [0.0, 1.0, 5.0] {
            let there = axis.distance_of(v, deepest);
            assert!(
                (axis.value(there, deepest) - v).abs() < 1e-4,
                "{axis:?} {v}"
            );
        }
    }
    // a scale of zero does not divide by it
    assert!(Axis {
        scale: 0.0,
        ..Axis::default()
    }
    .distance_of(3.0, deepest)
    .is_finite());
}

#[test]
fn ticks_are_at_round_numbers_in_the_range_of_the_axis_and_written_with_the_decimals_they_need() {
    let ticks = Axis::default().ticks(10.0, 5);
    let at: Vec<f32> = ticks.iter().map(|t| t.0).collect();
    assert_eq!(at, [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]);
    assert_eq!(ticks[1].1, "2");
    // fine steps get decimals; a tick at zero is written 0
    let fine = Axis::default().ticks(0.7, 7);
    assert_eq!(fine.first().map(|t| t.1.as_str()), Some("0.0"));
    assert_eq!(fine[1].1, "0.1");
    // reversed (time before the tips) the numbers grow to the left: 0 at the far end
    let before = Axis {
        mode: AxisMode::BeforeTips,
        ..Axis::default()
    }
    .ticks(10.0, 5);
    assert_eq!((before[0].0, before[0].1.as_str()), (10.0, "0"));
    assert_eq!(
        before.last().map(|t| (t.0, t.1.as_str())),
        Some((0.0, "10"))
    );
    // dates in a range that does not start at zero
    let dates = Axis {
        mode: AxisMode::Dates,
        scale: 1.0,
        tip_date: 2020.0,
        grid: true,
    }
    .ticks(100.0, 5);
    assert_eq!(dates[0].1, "1920");
    assert_eq!(dates.last().map(|t| t.1.as_str()), Some("2020"));
    // nothing to draw for a tree with no length
    assert!(Axis::default().ticks(0.0, 5).is_empty());
    assert!(Axis::default().ticks(f32::NAN, 5).is_empty());
    assert_eq!(tick_step(0.0, 5), 1.0);
    assert_eq!(tick_step(10.0, 5), 2.0);
    assert_eq!(tick_step(0.9, 5), 0.2);
}

#[test]
fn a_tree_whose_values_are_all_one_number_is_drawn_in_the_middle_of_the_scale() {
    let t = NwTree::parse_all("((A[&rate=2]:1,B[&rate=2]:1):1,C[&rate=2]:1);")
        .unwrap()
        .remove(0);
    let laid = lay_out(&t.to_node(), Layout::Phylogram);
    let d = decorate(
        &t,
        &laid,
        &StyleOptions {
            colour_by: Some("rate".into()),
            width_by: Some("rate".into()),
            ..StyleOptions::default()
        },
    );
    let a = index_of(&laid, "A");
    assert_eq!(d.colour[a], Some(gradient(0.5)));
    assert_eq!(d.width[a], 2.5);
}
