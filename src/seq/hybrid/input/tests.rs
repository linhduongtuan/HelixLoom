use super::*;
use crate::seq::debruijn::graph::rc;

// k = 5, so the links overlap by 4: A ends with CGGG and B begins with it
const GFA: &str = "H\tVN:Z:1.0\nS\tA\tAAACCCGGG\tLN:i:9\tKC:i:50\tDP:f:10.0\nS\tB\tCGGGTTTT\tDP:f:20.5\nL\tA\t+\tB\t+\t4M\n";

#[test]
fn segments_links_and_depths_are_read() {
    let i = read(GFA, 0).unwrap();
    let g = &i.graph;
    assert!(i.has_depth);
    assert_eq!((g.k(), g.count()), (5, 2));
    assert_eq!(g.segment(0).seq, b"AAACCCGGG");
    assert_eq!(g.segment(1).seq, b"CGGGTTTT");
    // the depth is the DP tag
    assert!((g.coverage(0) - 10.0).abs() < 1e-9);
    assert!((g.coverage(1) - 20.5).abs() < 1e-9);
    assert_eq!(g.succ(handle(0, false)), [handle(1, false)]);
    assert_eq!(g.succ(handle(1, true)), [handle(0, true)]);
}

#[test]
fn the_depth_comes_from_the_k_mer_count_when_there_is_no_dp_and_is_one_when_there_is_neither() {
    let text = "S\tA\tAAACCCGGG\tKC:i:50\nS\tB\tCGGGTTTT\nL\tA\t+\tB\t+\t4M\n";
    let i = read(text, 0).unwrap();
    // 50 k-mers counted over the 5 k-mers of A
    assert!((i.graph.coverage(0) - 10.0).abs() < 1e-9);
    assert!((i.graph.coverage(1) - 1.0).abs() < 1e-9);
    assert!(!i.has_depth);
}

#[test]
fn an_orientation_is_read_and_a_reverse_link_joins_the_other_ends() {
    // A leads into the reverse of C, whose reverse complement is CGGGTTTT's: CCCCAAAC... C is B reversed
    let b = "CGGGTTTT";
    let c: String = b
        .bytes()
        .rev()
        .map(|x| match x {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            _ => 'C',
        })
        .collect();
    let text = format!("S\tA\tAAACCCGGG\tDP:f:1\nS\tC\t{c}\tDP:f:1\nL\tA\t+\tC\t-\t4M\n");
    let i = read(&text, 0).unwrap();
    assert_eq!(i.graph.succ(handle(0, false)), [handle(1, true)]);
    assert_eq!(i.graph.succ(rc(handle(1, true))), [rc(handle(0, false))]);
}

#[test]
fn a_graph_without_links_takes_the_overlap_that_is_given() {
    let text = "S\tA\tAAACCCGGG\nS\tB\tCGGGTTTT\n";
    assert_eq!(read(text, 0).unwrap().graph.k(), 1);
    assert_eq!(read(text, 4).unwrap().graph.k(), 5);
    // the links say it when there are some, whatever was given
    assert_eq!(read(GFA, 9).unwrap().graph.k(), 5);
}

#[test]
fn what_is_wrong_with_a_graph_is_said_with_its_line() {
    let why = |t: &str| read(t, 0).err().unwrap();
    assert!(why("H\tVN:Z:1.0\n").contains("no segment"));
    assert!(why("S\tA\n").contains("Line 1") && why("S\tA\n").contains("name and its bases"));
    assert!(why("S\tA\t*\n").contains("no bases"));
    assert!(
        why("S\tA\tACGT\nS\tA\tACGT\n").contains("Line 2")
            && why("S\tA\tACGT\nS\tA\tACGT\n").contains("twice")
    );
    assert!(why("S\tA\tACGTACGT\nL\tA\t+\tB\t+\t4M\n").contains("not (yet) in the file"));
    assert!(why("S\tA\tACGTACGT\nL\tA\t+\n").contains("two segments"));
    assert!(why("S\tA\tACGTACGT\nS\tB\tACGTACGT\nL\tA\tx\tB\t+\t4M\n").contains("+ or -"));
    assert!(why("S\tA\tACGTACGT\nS\tB\tACGTACGT\nL\tA\t+\tB\t+\t*\n").contains("such as 54M"));
}

#[test]
fn links_that_do_not_all_overlap_by_the_same_number_or_do_not_match_are_refused() {
    let mixed =
        "S\tA\tAAACCCGGG\nS\tB\tCGGGTTTT\nS\tC\tGGTTTTAA\nL\tA\t+\tB\t+\t4M\nL\tB\t+\tC\t+\t3M\n";
    let e = read(mixed, 0).err().unwrap();
    assert!(
        e.contains("Line 5") && e.contains("overlaps by 3") && e.contains("first overlaps by 4"),
        "{e}"
    );
    let wrong = "S\tA\tAAACCCGGG\nS\tB\tTTTTTTTT\nL\tA\t+\tB\t+\t4M\n";
    let e = read(wrong, 0).err().unwrap();
    assert!(e.contains("Line 3") && e.contains("do not match"), "{e}");
}

#[test]
fn a_segment_shorter_than_k_is_refused() {
    let text = "S\tA\tACG\nS\tB\tCGGGTTTT\nL\tA\t+\tB\t+\t4M\n";
    let e = read(text, 0).err().unwrap();
    assert!(e.contains("Segment A has 3 bases"), "{e}");
}

#[test]
fn bases_are_taken_in_upper_case_and_other_lines_are_left_alone() {
    let text = "H\tVN:Z:1.0\nP\tpath\tA+\t*\n# a comment\nS\tA\tacgtacgt\tDP:f:3\n";
    let i = read(text, 0).unwrap();
    assert_eq!(i.graph.segment(0).seq, b"ACGTACGT");
}

#[test]
fn a_link_with_no_overlap_field_is_refused_and_a_segment_of_exactly_k_bases_is_allowed() {
    let e = read("S\tA\tAAACCCGGG\nS\tB\tCGGGT\nL\tA\t+\tB\t+\n", 0)
        .err()
        .unwrap();
    assert!(
        e.contains("Line 3") && e.contains("two segments with their orientations"),
        "{e}"
    );
    // k = 5: B has exactly 5 bases, and its first 4 are the end of A
    let i = read("S\tA\tAAACCCGGG\nS\tB\tCGGGT\nL\tA\t+\tB\t+\t4M\n", 0).unwrap();
    assert_eq!(i.graph.length(1), 5);
    assert_eq!(i.graph.kmers(1), 1);
}
