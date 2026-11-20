use super::*;
use crate::seq::treelayout::lay_out_with;
use crate::ui::app::Startup;
use crate::ui::testing::Screen;

const BEAST: &str = "#NEXUS\nbegin trees;\n translate 1 Alpha, 2 Beta, 3 Gamma, 4 Delta, 5 Epsilon;\n tree t = [&R] (((1[&rate=0.5]:1,2[&rate=1.5]:1)[&posterior=0.98,height_95%_HPD={0.4,1.6}]:1,3[&rate=1.0]:2)[&posterior=0.61]:1,(4:2,5:2)[&posterior=0.61]:2);\nend;\n";

fn trees_app() -> (App, Screen) {
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Trees;
    let trees = nexus::parse_any(BEAST).expect("trees");
    app.show_annotated(trees, "opened");
    let mut screen = Screen::new();
    screen.settle(&mut app);
    (app, screen)
}

/// The place in the layout of the node that has these tips below it (the layout the page draws now).
fn laid_node(app: &App, tips: &[&str]) -> usize {
    let page = &app.platform.trees;
    let laid = lay_out_with(
        app.tree_win.tree.as_ref().expect("a tree"),
        page.options.layout,
        &page.looks.collapsed,
    );
    (0..laid.nodes.len())
        .find(|&i| {
            let mut names: Vec<&str> = laid
                .tips_under(i)
                .iter()
                .map(|&t| laid.nodes[t].name.as_str())
                .collect();
            names.sort_unstable();
            let mut want = tips.to_vec();
            want.sort_unstable();
            names == want
        })
        .expect("a node with these tips")
}

fn click(app: &mut App, screen: &mut Screen, c: Clicks) {
    app.apply_tree_edits(&c);
    screen.settle(app);
}

#[test]
fn a_tree_from_beast_shows_its_posterior_beside_the_nodes_and_keeps_every_annotation() {
    let (app, screen) = trees_app();
    let page = &app.platform.trees;
    assert_eq!(page.looks.style.labels, ["posterior"]);
    let nw = page.annotated().expect("an annotated tree");
    assert_eq!(nw.attr_keys(), ["posterior", "height_95%_HPD", "rate"]);
    let text = screen.text();
    for want in ["Alpha", "Epsilon", "0.98", "0.61"] {
        assert!(text.contains(want), "{want} is missing:\n{text}");
    }
    assert!(!text.contains("[&"), "{text}");
}

#[test]
fn colouring_by_a_number_draws_its_scale_and_by_a_word_its_colours() {
    let (mut app, mut screen) = trees_app();
    app.platform.trees.looks.style.colour_by = Some("rate".into());
    screen.settle(&mut app);
    let text = screen.text();
    assert!(
        text.contains("rate") && text.contains("0.50") && text.contains("1.50"),
        "{text}"
    );
    app.platform.trees.looks.style.colour_by = Some("none".into());
    screen.settle(&mut app);
    assert!(!screen.text().contains("1.50"));
}

#[test]
fn every_layout_draws_every_tip_name() {
    let (mut app, mut screen) = trees_app();
    for layout in [
        Layout::Phylogram,
        Layout::Cladogram,
        Layout::Circular,
        Layout::Radial,
    ] {
        app.platform.trees.options.layout = layout;
        screen.settle(&mut app);
        let text = screen.text();
        for name in ["Alpha", "Beta", "Gamma", "Delta", "Epsilon"] {
            assert!(text.contains(name), "{layout:?} is missing {name}:\n{text}");
        }
    }
}

#[test]
fn collapsing_the_picked_clade_draws_it_as_a_triangle_and_expanding_puts_it_back() {
    let (mut app, mut screen) = trees_app();
    let ab = laid_node(&app, &["Alpha", "Beta"]);
    app.platform.trees.selected = Some(ab);
    // the card of the picked clade is drawn (and says its "2 tips")
    screen.settle(&mut app);
    let before = screen.text().matches("2 tips").count();
    assert!(before >= 1);
    click(
        &mut app,
        &mut screen,
        Clicks {
            fold: true,
            ..Clicks::default()
        },
    );
    assert_eq!(app.platform.trees.looks.collapsed.len(), 1);
    // the picked clade says "2 tips" in the card, and its triangle says it on the tree
    assert_eq!(screen.text().matches("2 tips").count(), before + 1);
    let laid = lay_out_with(
        app.tree_win.tree.as_ref().unwrap(),
        Layout::Phylogram,
        &app.platform.trees.looks.collapsed,
    );
    assert_eq!(laid.nodes.len(), 7);
    click(
        &mut app,
        &mut screen,
        Clicks {
            fold: true,
            ..Clicks::default()
        },
    );
    assert!(app.platform.trees.looks.collapsed.is_empty());
    assert_eq!(screen.text().matches("2 tips").count(), before);
}

#[test]
fn the_names_of_a_collapsed_clade_are_still_those_of_its_tips() {
    let (mut app, mut screen) = trees_app();
    let abc = laid_node(&app, &["Alpha", "Beta", "Gamma"]);
    app.platform.trees.selected = Some(abc);
    click(
        &mut app,
        &mut screen,
        Clicks {
            fold: true,
            ..Clicks::default()
        },
    );
    // the card lists the three, though the layout holds one node for them
    let text = screen.text();
    assert!(
        text.contains("3 tips") && text.contains("Alpha") && text.contains("Gamma"),
        "{text}"
    );
    // and a tree edit works on the real tree under it
    let node = app.platform.trees.selected.expect("still picked");
    assert!(node < 5);
}

#[test]
fn rooting_on_the_picked_branch_keeps_the_tips_and_the_annotations() {
    let (mut app, mut screen) = trees_app();
    let changed = app.tree_win.changed;
    let gamma = laid_node(&app, &["Gamma"]);
    app.platform.trees.selected = Some(gamma);
    click(
        &mut app,
        &mut screen,
        Clicks {
            reroot: true,
            ..Clicks::default()
        },
    );
    assert_eq!(app.tree_win.changed, changed + 1);
    let nw = app.platform.trees.annotated().expect("a tree");
    assert_eq!(nw.leaf_count(), 5);
    assert_eq!(
        nw.attr_number(
            nw.nodes.iter().position(|n| n.name == "Alpha").unwrap(),
            "rate"
        ),
        Some(0.5)
    );
    // the two inner branches that are not a tip's keep their posterior (the root's two branches were one, and say it twice)
    let mut posteriors: Vec<f64> = nw
        .nodes
        .iter()
        .filter_map(|n| {
            n.attrs
                .iter()
                .find(|(k, _)| k == "posterior")?
                .1
                .parse()
                .ok()
        })
        .collect();
    posteriors.sort_by(f64::total_cmp);
    assert_eq!(posteriors, [0.61, 0.98]);
    assert!(app.platform.trees.message.contains("Rooted"));
    assert!(app.platform.trees.selected.is_none());
    // with nothing picked it says what to do, and changes nothing
    click(
        &mut app,
        &mut screen,
        Clicks {
            rotate: true,
            ..Clicks::default()
        },
    );
    assert!(app.platform.trees.message.contains("Pick a node"));
    assert_eq!(app.tree_win.changed, changed + 1);
}

#[test]
fn rotating_and_ordering_change_only_the_order() {
    let (mut app, mut screen) = trees_app();
    let tips = |app: &App| -> Vec<String> {
        let nw = app.platform.trees.annotated().unwrap();
        nw.preorder()
            .into_iter()
            .filter(|&i| nw.is_leaf(i))
            .map(|i| nw.nodes[i].name.clone())
            .collect()
    };
    let first = tips(&app);
    app.platform.trees.selected = Some(0);
    click(
        &mut app,
        &mut screen,
        Clicks {
            rotate: true,
            ..Clicks::default()
        },
    );
    let rotated = tips(&app);
    assert_ne!(first, rotated);
    let (mut a, mut b) = (first.clone(), rotated);
    a.sort();
    b.sort();
    assert_eq!(a, b);
    click(
        &mut app,
        &mut screen,
        Clicks {
            order_down: true,
            ..Clicks::default()
        },
    );
    assert_eq!(tips(&app).len(), 5);
    assert!(app.platform.trees.message.contains("Larger"));
}

#[test]
fn weak_branches_are_collapsed_and_tips_renamed_from_a_table() {
    let (mut app, mut screen) = trees_app();
    app.platform.trees.looks.weak_key = "posterior".into();
    app.platform.trees.looks.weak_below = 0.7;
    click(
        &mut app,
        &mut screen,
        Clicks {
            weak: true,
            ..Clicks::default()
        },
    );
    assert!(
        app.platform
            .trees
            .message
            .starts_with("Took out 2 branches"),
        "{}",
        app.platform.trees.message
    );
    let nw = app.platform.trees.annotated().unwrap();
    // the two weak branches are gone: five tips, the one clade that stays, and the root
    assert_eq!(nw.nodes.len(), 5 + 2);
    app.platform.trees.looks.names_paste = "Alpha\tHomo sapiens\nnobody,Nobody\n".into();
    click(
        &mut app,
        &mut screen,
        Clicks {
            rename: true,
            ..Clicks::default()
        },
    );
    assert!(
        app.platform.trees.message.starts_with("Renamed 1 tip"),
        "{}",
        app.platform.trees.message
    );
    assert!(screen.text().contains("Homo sapiens"));
}

#[test]
fn support_from_other_trees_is_added_as_an_annotation_and_shown() {
    let (mut app, mut screen) = trees_app();
    let same = "(((Alpha,Beta),Gamma),(Delta,Epsilon));";
    let other = "(((Alpha,Gamma),Beta),(Delta,Epsilon));";
    let trees: Vec<NwTree> = [same, same, same, other]
        .iter()
        .map(|t| NwTree::parse_all(t).unwrap().remove(0))
        .collect();
    app.platform.trees.looks.support_label = "boot".into();
    app.platform.trees.looks.awaiting_support = true;
    app.add_support(&trees, "four trees");
    assert!(!app.platform.trees.looks.awaiting_support);
    assert!(
        app.platform.trees.message.contains("\"boot\"")
            && app.platform.trees.message.contains("4 trees"),
        "{}",
        app.platform.trees.message
    );
    let nw = app.platform.trees.annotated().unwrap();
    let ab = nw
        .nodes
        .iter()
        .position(|n| {
            !n.children.is_empty()
                && n.children
                    .iter()
                    .all(|&c| nw.nodes[c].name == "Alpha" || nw.nodes[c].name == "Beta")
        })
        .unwrap();
    assert_eq!(nw.attr_number(ab, "boot"), Some(75.0));
    assert_eq!(nw.attr_number(ab, "boot_conflict"), Some(25.0));
    assert!(app
        .platform
        .trees
        .looks
        .style
        .labels
        .contains(&"boot".to_string()));
    screen.settle(&mut app);
    assert!(screen.text().contains("75"), "{}", screen.text());
    // no tree to add to, or trees that share too little, says so
    app.platform.trees.doc = None;
    app.add_support(&trees, "x");
    assert!(app.platform.trees.message.contains("no tree on the page"));
}

#[test]
fn an_axis_with_dates_is_drawn_under_the_phylogram_and_find_marks_the_tips_that_match() {
    let (mut app, mut screen) = trees_app();
    let looks = &mut app.platform.trees.looks;
    looks.axis_on = true;
    looks.axis.mode = crate::seq::treestyle::AxisMode::Dates;
    looks.axis.tip_date = 2020.0;
    screen.settle(&mut app);
    let text = screen.text();
    // the tips are at 2020 and the root is three length units before them
    for date in ["2017", "2018", "2019", "2020"] {
        assert!(text.contains(date), "{date}:\n{text}");
    }
}

#[test]
fn find_marks_the_tips_whose_names_hold_the_text_in_any_case() {
    let (app, _) = trees_app();
    let page = &app.platform.trees;
    let laid = lay_out_with(
        app.tree_win.tree.as_ref().unwrap(),
        Layout::Phylogram,
        &page.looks.collapsed,
    );
    let found = |t: &str| -> Vec<String> {
        found_tips(&laid, t)
            .iter()
            .map(|&i| laid.nodes[i].name.clone())
            .collect()
    };
    assert_eq!(found("EPS"), ["Epsilon"]);
    assert_eq!(found("a"), ["Alpha", "Beta", "Gamma", "Delta"]);
    assert!(found("").is_empty() && found("  ").is_empty() && found("zzz").is_empty());
}

#[test]
fn another_tree_put_on_the_page_takes_the_place_of_the_annotated_one_with_all_it_had() {
    let (mut app, mut screen) = trees_app();
    app.platform.trees.looks.style.colour_by = Some("rate".into());
    app.tree_win.set_tree(Some(crate::sample::demo_tree()));
    screen.settle(&mut app);
    let page = &app.platform.trees;
    assert_eq!(page.annotated().map(NwTree::leaf_count), Some(24));
    assert!(page.looks.style.colour_by.is_none() && page.looks.collapsed.is_empty());
}

#[test]
fn a_nexus_file_that_is_opened_is_read_and_a_bad_one_says_why() {
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Trees;
    let dir = std::env::temp_dir().join(format!("helixloom-nexus-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let good = dir.join("beast.nex");
    std::fs::write(&good, BEAST).unwrap();
    app.load_tree(Source::Path(&good));
    assert_eq!(
        app.platform.trees.annotated().map(NwTree::leaf_count),
        Some(5)
    );
    assert!(
        app.platform.trees.message.starts_with("Opened"),
        "{}",
        app.platform.trees.message
    );
    let bad = dir.join("bad.nex");
    std::fs::write(&bad, "#NEXUS\nbegin taxa;\nend;\n").unwrap();
    app.load_tree(Source::Path(&bad));
    assert!(
        app.platform.trees.message.contains("no trees block"),
        "{}",
        app.platform.trees.message
    );
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn a_highlighted_clade_has_a_block_of_colour_behind_it_on_a_phylogram_and_none_when_taken_off() {
    let (mut app, mut screen) = trees_app();
    let block = |c: [u8; 3]| egui::Color32::from_rgba_unmultiplied(c[0], c[1], c[2], 46);
    let tone = crate::seq::treestyle::PALETTE[2];
    assert!(screen.rects_filled(block(tone)).is_empty());
    let ab = laid_node(&app, &["Alpha", "Beta"]);
    app.platform.trees.selected = Some(ab);
    click(
        &mut app,
        &mut screen,
        Clicks {
            colour: Some(2),
            ..Clicks::default()
        },
    );
    assert_eq!(app.platform.trees.looks.highlights.len(), 1);
    let blocks = screen.rects_filled(block(tone));
    assert_eq!(blocks.len(), 1);
    // it spans the two tips of the clade, and not the other three
    let rows: Vec<f32> = ["Alpha", "Beta", "Gamma"]
        .iter()
        .map(|n| screen.find_all(n).last().expect("the name on the tree").y)
        .collect();
    assert!(
        blocks[0].top() < rows[0] && blocks[0].bottom() > rows[1] && blocks[0].bottom() < rows[2]
    );
    // a second colour for the same clade replaces the first
    click(
        &mut app,
        &mut screen,
        Clicks {
            colour: Some(3),
            ..Clicks::default()
        },
    );
    assert_eq!(app.platform.trees.looks.highlights.len(), 1);
    assert!(screen.rects_filled(block(tone)).is_empty());
    click(
        &mut app,
        &mut screen,
        Clicks {
            unmark: true,
            ..Clicks::default()
        },
    );
    assert!(app.platform.trees.looks.highlights.is_empty());
    assert!(screen
        .rects_filled(block(crate::seq::treestyle::PALETTE[3]))
        .is_empty());
}

#[test]
fn node_shapes_follow_the_annotation_and_the_legend_names_what_the_colours_are() {
    let (mut app, mut screen) = trees_app();
    // radio buttons are round too: count what the shapes add
    let (big, small) = (
        crate::seq::treestyle::SHAPE_RADII.1,
        crate::seq::treestyle::SHAPE_RADII.0,
    );
    let (had_big, had_small) = (
        screen.circles_of_radius(big).len(),
        screen.circles_of_radius(small).len(),
    );
    app.platform.trees.looks.style.shape_by = Some("posterior".into());
    screen.settle(&mut app);
    // the greatest posterior (0.98) is on one node, the least (0.61) on two (the root's two branches); each shape is a disc and its ring
    assert_eq!(screen.circles_of_radius(big).len(), had_big + 2);
    assert_eq!(screen.circles_of_radius(small).len(), had_small + 4);
    app.platform.trees.looks.style.shape_by = None;
    screen.settle(&mut app);
    assert_eq!(screen.circles_of_radius(big).len(), had_big);
}
