use super::*;
use crate::seq::debruijn::graph::handle;
use crate::seq::debruijn::sim::random;

const K: usize = 5;

fn graph() -> Graph {
    let mut g = Graph::new(K);
    g.add(random(80, 1), 76.0 * 10.0);
    g.add(random(40, 2), 36.0 * 2.5);
    g
}

#[test]
fn the_replicons_are_numbered_by_length_with_their_depth() {
    let r = replicons(&graph(), &[], &RotateOptions::default());
    assert_eq!(
        r.iter()
            .map(|r| (r.number, r.seq.len()))
            .collect::<Vec<_>>(),
        [(1, 80), (2, 40)]
    );
    assert_eq!((r[0].depth, r[1].depth), (10.0, 2.5));
    assert!(r.iter().all(|r| !r.circular && !r.rotated));
}

#[test]
fn the_fasta_name_has_the_number_the_length_the_depth_and_for_a_circle_circular() {
    let mut g = graph();
    // a circle: 60 bases, 56 k-mers, its last 4 are its first 4
    let ring = random(56, 3);
    let seq = [&ring[..], &ring[..4]].concat();
    let id = g.add(seq, 56.0 * 4.0);
    g.add_link(handle(id, false), handle(id, false));
    let r = replicons(&g, &[], &RotateOptions::default());
    let text = fasta(&r);
    let names: Vec<&str> = text.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(
        names,
        [
            ">1 length=80 depth=10.00x",
            ">2 length=56 depth=4.00x circular=true",
            ">3 length=40 depth=2.50x"
        ]
    );
    // the circle is once round
    assert_eq!(r[1].seq, ring);
}

#[test]
fn the_bases_are_in_lines_of_seventy() {
    let g = graph();
    let text = fasta(&replicons(&g, &[], &RotateOptions::default()));
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!((lines[1].len(), lines[2].len()), (70, 10));
}

#[test]
fn a_circle_is_turned_to_a_start_gene_and_says_so_while_a_line_is_not() {
    let mut g = Graph::new(K);
    let ring = random(500, 4);
    let id = g.add([&ring[..], &ring[..4]].concat(), 500.0 * 4.0);
    g.add_link(handle(id, false), handle(id, false));
    g.add(random(300, 5), 296.0 * 4.0);
    let gene = ring[200..380].to_vec();
    let r = replicons(
        &g,
        &[("dnaA".to_string(), gene.clone())],
        &RotateOptions::default(),
    );
    assert!(r[0].circular && r[0].rotated);
    assert_eq!(&r[0].seq[..180], &gene[..]);
    assert!(!r[1].rotated);
    // no gene: not turned
    let none = replicons(&g, &[], &RotateOptions::default());
    assert!(!none[0].rotated && none[0].seq == ring);
}

#[test]
fn a_line_is_not_turned_to_a_start_gene_even_when_the_gene_is_in_it() {
    let mut g = Graph::new(K);
    let line = random(300, 6);
    g.add(line.clone(), 296.0 * 4.0);
    let gene = line[100..250].to_vec();
    let r = replicons(&g, &[("dnaA".to_string(), gene)], &RotateOptions::default());
    assert!(!r[0].circular && !r[0].rotated);
    assert_eq!(r[0].seq, line);
}
