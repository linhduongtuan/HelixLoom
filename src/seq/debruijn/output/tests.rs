use super::*;
use crate::seq::debruijn::graph::handle;

const K: usize = 4;

fn f(id: usize) -> usize {
    handle(id, false)
}

fn r(id: usize) -> usize {
    handle(id, true)
}

/// Three segments: a short one, a long one, and one of the length of the first (to see how the numbering breaks the tie), and a circle.
fn graph() -> Graph {
    let mut g = Graph::new(K);
    g.add(b"GGGGTT".to_vec(), 30.0); // 3 k-mers, cov 10
    g.add(b"ACGTACGTACGT".to_vec(), 90.0); // 9 k-mers, cov 10
    g.add(b"AAAATT".to_vec(), 6.0); // 3 k-mers, cov 2
    g
}

#[test]
fn the_longest_segment_is_number_one_and_a_tie_goes_by_the_bases() {
    let g = graph();
    let c = number(&g);
    assert_eq!(c.iter().map(|c| c.number).collect::<Vec<_>>(), [1, 2, 3]);
    assert_eq!(c.iter().map(|c| c.id).collect::<Vec<_>>(), [1, 2, 0]);
    assert_eq!(c[0].seq, b"ACGTACGTACGT");
    assert_eq!(c[0].coverage, 10.0);
    assert!(!c[0].circular);
}

#[test]
fn a_removed_segment_has_no_number() {
    let mut g = graph();
    g.remove(1);
    assert_eq!(number(&g).len(), 2);
}

#[test]
fn a_contig_is_written_with_the_name_the_spades_manual_gives_and_its_bases_in_lines_of_sixty() {
    let mut g = Graph::new(K);
    let seq: Vec<u8> = b"ACGT".iter().copied().cycle().take(130).collect();
    g.add(seq, 127.0 * 3.5);
    let c = number(&g);
    let text = fasta(&c, 1);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], ">NODE_1_length_130_cov_3.500000");
    assert_eq!(lines.len(), 1 + 3);
    assert_eq!(
        (lines[1].len(), lines[2].len(), lines[3].len()),
        (60, 60, 10)
    );
}

#[test]
fn only_the_contigs_of_the_length_asked_for_or_more_are_written_but_the_numbers_stay() {
    let c = number(&graph());
    let text = fasta(&c, 6);
    assert_eq!(text.matches('>').count(), 3);
    let text = fasta(&c, 7);
    assert_eq!(text.matches('>').count(), 1);
    assert!(text.starts_with(">NODE_1_length_12_"));
    let text = fasta(&c, 1);
    assert!(text.contains(">NODE_2_length_6_") && text.contains(">NODE_3_length_6_"));
    assert_eq!(fasta(&c, 100), "");
}

#[test]
fn a_circle_is_written_once_round_and_says_so() {
    let mut g = Graph::new(K);
    // 8 k-mers: 8 + k - 1 = 11 bases; the last 3 are the first 3
    let id = g.add(b"ACGGTCAACGG".to_vec(), 16.0);
    g.add_link(f(id), f(id));
    let c = number(&g);
    assert!(c[0].circular);
    assert_eq!(c[0].seq, b"ACGGTCAA");
    let text = fasta(&c, 1);
    assert!(text.starts_with(">NODE_1_length_8_cov_2.000000 circular\n"));
    assert!(text.contains("\nACGGTCAA\n"));
    // the graph keeps the whole segment
    let gfa = gfa(&g, &c);
    assert!(gfa.contains("S\t1\tACGGTCAACGG\tLN:i:11"));
    assert!(gfa.contains("L\t1\t+\t1\t+\t3M"));
}

#[test]
fn the_graph_is_gfa_with_a_segment_for_each_contig_and_the_links_overlapping_by_k_minus_one() {
    let mut g = graph();
    g.add_link(f(1), f(0));
    g.add_link(f(1), r(2));
    let c = number(&g);
    let text = gfa(&g, &c);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "H\tVN:Z:1.0");
    assert_eq!(
        lines[1],
        "S\t1\tACGTACGTACGT\tLN:i:12\tKC:i:90\tDP:f:10.000000"
    );
    assert_eq!(lines[2], "S\t2\tAAAATT\tLN:i:6\tKC:i:6\tDP:f:2.000000");
    assert_eq!(lines[3], "S\t3\tGGGGTT\tLN:i:6\tKC:i:30\tDP:f:10.000000");
    // contig 1 (id 1) leads to contig 3 (id 0), and to contig 2 (id 2) on the other strand; the lines are sorted
    assert_eq!(&lines[4..], ["L\t1\t+\t2\t-\t3M", "L\t1\t+\t3\t+\t3M"]);
}

#[test]
fn every_link_is_written_once_though_it_is_held_twice() {
    let mut g = graph();
    g.add_link(f(0), f(1));
    let text = gfa(&g, &number(&g));
    assert_eq!(text.lines().filter(|l| l.starts_with('L')).count(), 1);
}

#[test]
fn n50_is_the_length_at_which_half_the_bases_are_in_contigs_that_long_or_longer() {
    assert_eq!(n50(&[100]), 100);
    assert_eq!(n50(&[10, 20, 30, 40]), 30);
    assert_eq!(n50(&[40, 30, 20, 10]), 30);
    assert_eq!(n50(&[50, 50]), 50);
    assert_eq!(n50(&[2, 2, 2, 2, 2, 2, 2, 2, 2, 100]), 100);
    assert_eq!(n50(&[]), 0);
}

#[test]
fn the_k_mer_count_is_rounded_to_the_nearest_whole_number() {
    let mut g = Graph::new(K);
    g.add(b"ACGTACGT".to_vec(), 29.6);
    g.add(b"TTTTGGGGCC".to_vec(), 40.4);
    let text = gfa(&g, &number(&g));
    assert!(text.contains("\tKC:i:30\t"), "{text}");
    assert!(text.contains("\tKC:i:40\t"), "{text}");
}

#[test]
fn n50_counts_the_length_that_reaches_half_the_bases_itself() {
    // 30 + 20 is exactly half of 100: the length that gets there is 20, not the next
    assert_eq!(n50(&[30, 20, 15, 15, 10, 10]), 20);
}
