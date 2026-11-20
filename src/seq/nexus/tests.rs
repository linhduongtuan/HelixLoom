use super::*;

const BEAST: &str = "#NEXUS\n\nBegin taxa;\n\tDimensions ntax=3;\n\tTaxlabels\n\t\tA\n\t\tB\n\t\tC\n\t\t;\nEnd;\n\nBegin trees;\n\tTranslate\n\t\t1 'Homo sapiens',\n\t\t2 Pan,\n\t\t3 Gorilla\n\t\t;\ntree TREE1 = [&R] ((1[&rate=1.1]:0.5,2:0.5)[&posterior=0.99,height_95%_HPD={0.4,0.6}]:1,3:1.5);\ntree * second = [&U] (1:1,2:1,3:1);\nEnd;\n";

#[test]
fn trees_are_read_with_their_names_translated_tip_names_and_annotations() {
    let trees = parse_trees(BEAST).unwrap();
    assert_eq!(
        trees.iter().map(|t| t.0.as_str()).collect::<Vec<_>>(),
        ["TREE1", "second"]
    );
    let t = &trees[0].1;
    let tips: Vec<&str> = t
        .preorder()
        .into_iter()
        .filter(|&i| t.is_leaf(i))
        .map(|i| t.nodes[i].name.as_str())
        .collect();
    assert_eq!(tips, ["Homo sapiens", "Pan", "Gorilla"]);
    let human = t
        .nodes
        .iter()
        .position(|n| n.name == "Homo sapiens")
        .unwrap();
    assert_eq!(t.attr_number(human, "rate"), Some(1.1));
    let inner = t.nodes[human].parent.unwrap();
    assert_eq!(t.attr_number(inner, "posterior"), Some(0.99));
    assert_eq!(t.attr_range(inner, "height_95%_HPD"), Some((0.4, 0.6)));
    assert!((t.nodes[inner].dist - 1.0).abs() < 1e-12);
    // the second tree is translated too
    assert!(trees[1].1.nodes.iter().any(|n| n.name == "Gorilla"));
}

#[test]
fn a_file_without_trees_or_that_is_not_nexus_says_so() {
    assert!(parse_trees("(A,B);").unwrap_err().contains("not a NEXUS"));
    assert!(parse_trees("#NEXUS\nbegin taxa;\nend;\n")
        .unwrap_err()
        .contains("no trees block"));
    assert!(parse_trees("#NEXUS\nbegin trees;\nend;\n")
        .unwrap_err()
        .contains("no tree statement"));
    assert!(parse_trees("#NEXUS\nbegin trees;\ntree t (A,B);\nend;\n")
        .unwrap_err()
        .contains("equals"));
    assert!(parse_trees("#NEXUS\nbegin trees;\ntree t = (A,B;\nend;\n").is_err());
}

#[test]
fn case_comments_quotes_and_a_semicolon_in_a_name_do_not_confuse_the_reader() {
    let text = "#nexus\n[a comment; with a semicolon]\nBEGIN TREES;\n[another; one]\n\tTRANSLATE 1 'a;b', 2 c;\n\tTREE t1 = (1:1,2:1);\nENDBLOCK;\n";
    let trees = parse_trees(text).unwrap();
    let names: Vec<&str> = trees[0]
        .1
        .nodes
        .iter()
        .filter(|n| n.children.is_empty())
        .map(|n| n.name.as_str())
        .collect();
    assert_eq!(names, ["a;b", "c"]);
    assert!(is_nexus("\u{feff}  #NEXUS"));
    assert!(!is_nexus("#NEX"));
}

#[test]
fn nexus_and_newick_text_are_both_read_by_parse_any() {
    assert_eq!(parse_any(BEAST).unwrap().len(), 2);
    assert_eq!(parse_any("(A:1,B:1);(C:1,D:1);").unwrap().len(), 2);
    assert!(parse_any("not a tree").is_err());
}

#[test]
fn written_nexus_reads_back_with_its_annotations() {
    let t = &parse_trees(BEAST).unwrap()[0].1;
    for annotated in [true, false] {
        let text = write(t, "my tree", annotated);
        assert!(text.contains("tree 'my tree' = [&R]"), "{text}");
        let back = parse_trees(&text).unwrap();
        assert_eq!(back[0].0, "my tree");
        let inner = back[0]
            .1
            .nodes
            .iter()
            .position(|n| !n.children.is_empty() && n.parent.is_some())
            .unwrap();
        assert_eq!(
            back[0].1.attr_number(inner, "posterior").is_some(),
            annotated
        );
        assert_eq!(back[0].1.leaf_count(), 3);
    }
}

#[test]
fn an_unrooted_tree_statement_is_a_tree_too() {
    let trees = parse_trees("#NEXUS\nbegin trees;\n utree u1 = (A:1,B:1,C:1);\nend;\n").unwrap();
    assert_eq!((trees[0].0.as_str(), trees[0].1.leaf_count()), ("u1", 3));
}
