use super::*;
use crate::seq::debruijn::sim::random;
use crate::seq::hybrid::sim::is_circle_of;

fn turned(seq: &[u8], at: usize) -> Vec<u8> {
    [&seq[at..], &seq[..at]].concat()
}

fn genes(g: Vec<u8>) -> Vec<(String, Vec<u8>)> {
    vec![("dnaA".to_string(), g)]
}

fn substituted(g: &[u8], every: usize) -> Vec<u8> {
    g.iter()
        .enumerate()
        .map(|(i, &b)| {
            if i % every == every - 1 {
                if b == b'A' {
                    b'C'
                } else {
                    b'A'
                }
            } else {
                b
            }
        })
        .collect()
}

#[test]
fn a_gene_on_the_forward_strand_becomes_the_start() {
    let circle = random(3000, 1);
    let gene = circle[1000..1400].to_vec();
    let r = rotate(&circle, &genes(gene.clone()), &Options::default()).unwrap();
    assert_eq!(&r[..400], &gene[..]);
    assert_eq!(r, turned(&circle, 1000));
}

#[test]
fn a_gene_on_the_other_strand_flips_the_circle_and_starts_there() {
    let circle = random(3000, 2);
    let gene = revcom(&circle[2000..2400]);
    let r = rotate(&circle, &genes(gene.clone()), &Options::default()).unwrap();
    assert_eq!(&r[..400], &gene[..]);
    assert!(is_circle_of(&r, &circle));
    // the flipped circle, turned to where the gene starts
    let flipped = revcom(&circle);
    assert_eq!(r, turned(&flipped, 3000 - 2400));
}

#[test]
fn a_gene_across_the_join_is_found_whole() {
    let circle = random(3000, 3);
    // the circle turned so that the gene lies over the end
    let at = 2800;
    let gene = [&circle[2800..], &circle[..200]].concat();
    let r = rotate(&circle, &genes(gene.clone()), &Options::default()).unwrap();
    assert_eq!(&r[..400], &gene[..]);
    assert_eq!(r, turned(&circle, at));
}

#[test]
fn a_gene_with_a_base_in_twenty_changed_is_found_and_the_start_is_where_it_begins() {
    let circle = random(4000, 4);
    let gene = substituted(&circle[1500..2200], 20);
    let r = rotate(&circle, &genes(gene), &Options::default()).unwrap();
    // the matching bases are all on one line, so the start is exact whatever bases the alignment began at
    assert_eq!(r, turned(&circle, 1500));
}

#[test]
fn of_several_genes_the_best_match_decides() {
    let circle = random(3000, 5);
    let weak = substituted(&circle[500..800], 20);
    let strong = circle[2000..2300].to_vec();
    let r = rotate(
        &circle,
        &[
            ("weak".to_string(), weak),
            ("strong".to_string(), strong.clone()),
        ],
        &Options::default(),
    )
    .unwrap();
    assert_eq!(&r[..300], &strong[..]);
}

#[test]
fn a_gene_that_is_not_there_or_an_empty_list_leaves_the_circle_alone() {
    let circle = random(3000, 6);
    assert!(rotate(&circle, &genes(random(400, 99)), &Options::default()).is_none());
    assert!(rotate(&circle, &[], &Options::default()).is_none());
    assert!(rotate(&[], &genes(random(400, 99)), &Options::default()).is_none());
}

#[test]
fn a_gene_that_the_alignment_covers_too_little_of_is_not_taken() {
    let circle = random(3000, 7);
    // 150 bases of the circle in a gene of 1000 (the rest not from it): 15 % covered
    let gene = [&circle[1000..1150], &random(850, 98)[..]].concat();
    assert!(rotate(&circle, &genes(gene.clone()), &Options { min_cover: 0.3 }).is_none());
    assert!(rotate(&circle, &genes(gene), &Options { min_cover: 0.1 }).is_some());
}

#[test]
fn the_default_share_of_a_gene_that_must_be_found_is_three_tenths() {
    assert_eq!(Options::default(), Options { min_cover: 0.3 });
}
