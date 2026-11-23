//! The assembler on reads simulated from sequences that are known, so that what it gives back can be compared with what it was given. The
//! reads have no pairs, and the limits that come of that are tested as limits: a repeat longer than a read stays a repeat.

use super::assemble::{assemble, Assembly, Cutoff, Options};
use super::output::{gfa, number, Contig};
use super::sim::{random, reads, tiled};
use crate::seq::dna::revcom;

fn contains(hay: &[u8], needle: &[u8]) -> bool {
    hay.windows(needle.len()).any(|w| w == needle)
}

/// Is `contig` a stretch of `genome` or of its reverse complement (going round the join, for a circle)?
fn belongs(contig: &[u8], genome: &[u8], circular: bool) -> bool {
    let line: Vec<u8> = if circular {
        [genome, genome].concat()
    } else {
        genome.to_vec()
    };
    contains(&line, contig) || contains(&revcom(&line), contig)
}

fn options(ks: &[usize], cutoff: Cutoff) -> Options {
    Options {
        ks: ks.to_vec(),
        cutoff,
        ..Options::default()
    }
}

fn run(pieces: &[Vec<u8>], ks: &[usize]) -> (Assembly, Vec<Contig>) {
    let a = assemble(pieces, &options(ks, Cutoff::Auto)).unwrap();
    let c = number(&a.graph);
    (a, c)
}

fn is_genome(contig: &[u8], genome: &[u8]) -> bool {
    contig == genome || contig == revcom(genome)
}

#[test]
fn perfect_reads_of_a_line_give_the_line_back_at_one_k_and_at_several() {
    let genome = random(4000, 1);
    let pieces = reads(&genome, false, 100, 30.0, 0.0, 2);
    for ks in [vec![21], vec![31], vec![21, 33, 55], vec![]] {
        let (_, c) = run(&pieces, &ks);
        assert_eq!(c.len(), 1, "ks {ks:?}");
        assert!(is_genome(&c[0].seq, &genome), "ks {ks:?}");
        assert!(!c[0].circular);
    }
}

#[test]
fn the_coverage_written_is_the_kmer_coverage_at_the_last_k() {
    let genome = random(4000, 3);
    // 30 times, reads of 100: a k-mer of 21 is in about 30 x 80 / 100 = 24 reads, one of 55 in 30 x 46 / 100 = 13.8
    let pieces = reads(&genome, false, 100, 30.0, 0.0, 4);
    let (_, c) = run(&pieces, &[21]);
    assert!(
        c[0].coverage > 20.0 && c[0].coverage < 28.0,
        "{}",
        c[0].coverage
    );
    let (_, c) = run(&pieces, &[21, 55]);
    assert!(
        c[0].coverage > 11.0 && c[0].coverage < 17.0,
        "{}",
        c[0].coverage
    );
}

#[test]
fn reads_of_a_circle_give_one_contig_that_is_a_circle_of_the_same_length() {
    let genome = random(3000, 5);
    let pieces = reads(&genome, true, 100, 40.0, 0.0, 6);
    let (a, c) = run(&pieces, &[]);
    assert_eq!(c.len(), 1);
    assert!(c[0].circular);
    assert_eq!(c[0].seq.len(), 3000);
    assert!(belongs(&c[0].seq, &genome, true));
    assert!(a.graph.circular(c[0].id));
}

#[test]
fn a_circle_that_is_read_with_errors_is_still_one_circle() {
    let genome = random(3000, 7);
    let pieces = reads(&genome, true, 100, 60.0, 0.005, 8);
    let (_, c) = run(&pieces, &[]);
    assert_eq!(
        c.len(),
        1,
        "{:?}",
        c.iter().map(|c| c.seq.len()).collect::<Vec<_>>()
    );
    assert!(c[0].circular && c[0].seq.len() == 3000);
    assert!(belongs(&c[0].seq, &genome, true));
}

#[test]
fn errors_in_the_reads_are_cleaned_away_and_what_is_left_is_the_genome() {
    let genome = random(5000, 9);
    let pieces = reads(&genome, false, 100, 60.0, 0.01, 10);
    let (a, c) = run(&pieces, &[]);
    let lengths: Vec<usize> = c.iter().map(|c| c.seq.len()).collect();
    assert!(
        c.iter().all(|c| belongs(&c.seq, &genome, false)),
        "{lengths:?}"
    );
    assert_eq!(c.len(), 1, "{lengths:?}");
    assert!(lengths[0] >= 4900, "{lengths:?}");
    // the cleaning had something to do
    assert!(a.steps.iter().map(|s| s.tips + s.bulges).sum::<usize>() > 0);
}

#[test]
fn two_copies_that_differ_by_one_base_are_one_contig_after_the_bulge_is_removed() {
    let a = random(2000, 11);
    let mut b = a.clone();
    b[1000] = if b[1000] == b'A' { b'C' } else { b'A' };
    let mut pieces = reads(&a, false, 100, 25.0, 0.0, 12);
    pieces.extend(reads(&b, false, 100, 25.0, 0.0, 13));
    let (asm, c) = run(&pieces, &[21]);
    assert_eq!(c.len(), 1);
    assert!(is_genome(&c[0].seq, &a) || is_genome(&c[0].seq, &b));
    assert_eq!(asm.steps[0].bulges, 1);
    // and the one that stayed was seen as often as both
    assert!(c[0].coverage > 35.0, "{}", c[0].coverage);
}

fn with_repeat(repeat: usize) -> (Vec<u8>, Vec<u8>) {
    let r = random(repeat, 20);
    let parts = [
        random(1000, 21),
        r.clone(),
        random(1000, 22),
        r,
        random(1000, 23),
    ];
    (parts.concat(), parts[1].clone())
}

#[test]
fn a_repeat_shorter_than_the_larger_k_is_resolved_by_it_and_not_by_the_smaller() {
    let (genome, _) = with_repeat(28);
    let pieces = reads(&genome, false, 100, 40.0, 0.0, 24);
    // at k = 21 the repeat is 8 k-mers that are the same in two places: the graph branches
    let (_, small) = run(&pieces, &[21]);
    assert!(small.len() >= 3, "{}", small.len());
    assert!(small.iter().all(|c| belongs(&c.seq, &genome, false)));
    // at k = 33 no k-mer is in both: one contig
    let (_, both) = run(&pieces, &[21, 33]);
    assert_eq!(both.len(), 1);
    assert!(is_genome(&both[0].seq, &genome));
    let (_, large) = run(&pieces, &[33]);
    assert_eq!(large.len(), 1);
}

#[test]
fn a_repeat_longer_than_the_reads_is_not_resolved_without_pairs_and_the_graph_shows_it() {
    let (genome, _) = with_repeat(300);
    let pieces = reads(&genome, false, 100, 40.0, 0.0, 25);
    let (a, c) = run(&pieces, &[]);
    assert!(c.len() >= 3, "{}", c.len());
    assert!(c.iter().all(|c| belongs(&c.seq, &genome, false)));
    // the contigs are linked where they meet: the repeat leads on to both what follows it
    let text = gfa(&a.graph, &c);
    assert!(text.lines().filter(|l| l.starts_with('L')).count() >= 4);
    // and no contig is the whole genome
    assert!(c.iter().all(|c| c.seq.len() < genome.len()));
}

#[test]
fn what_is_seen_few_times_beside_a_genome_seen_often_is_taken_away_by_the_automatic_cutoff_and_kept_without(
) {
    let genome = random(4000, 31);
    let other = random(1500, 32);
    let mut pieces = tiled(&genome, 100, 2);
    pieces.extend(tiled(&other, 100, 50));
    let (_, kept) = {
        let a = assemble(&pieces, &options(&[21], Cutoff::Off)).unwrap();
        let c = number(&a.graph);
        (a, c)
    };
    assert_eq!(kept.len(), 2);
    let (a, c) = run(&pieces, &[21]);
    assert_eq!(c.len(), 1);
    assert!(is_genome(&c[0].seq, &genome));
    assert!(a.steps[0].low_coverage >= 1 && a.steps[0].cutoff > 3.0);
}

#[test]
fn a_chromosome_and_a_plasmid_seen_at_different_depths_are_two_circles() {
    let chromosome = random(3000, 41);
    let plasmid = random(1500, 42);
    let mut pieces = reads(&chromosome, true, 100, 30.0, 0.002, 43);
    pieces.extend(reads(&plasmid, true, 100, 90.0, 0.002, 44));
    let (_, c) = run(&pieces, &[]);
    assert_eq!(
        c.len(),
        2,
        "{:?}",
        c.iter().map(|c| c.seq.len()).collect::<Vec<_>>()
    );
    assert!(c.iter().all(|c| c.circular));
    assert_eq!((c[0].seq.len(), c[1].seq.len()), (3000, 1500));
    assert!(belongs(&c[0].seq, &chromosome, true) && belongs(&c[1].seq, &plasmid, true));
    // the plasmid is seen about three times as often
    assert!(
        c[1].coverage > 2.0 * c[0].coverage,
        "{} {}",
        c[1].coverage,
        c[0].coverage
    );
}

#[test]
fn reads_from_the_other_strand_give_the_same_contigs() {
    let genome = random(3000, 51);
    let pieces = reads(&genome, false, 100, 40.0, 0.003, 52);
    let flipped: Vec<Vec<u8>> = pieces.iter().map(|p| revcom(p)).collect();
    let canon = |c: &[Contig]| {
        let mut v: Vec<Vec<u8>> = c
            .iter()
            .map(|c| c.seq.clone().min(revcom(&c.seq)))
            .collect();
        v.sort();
        v
    };
    let (_, a) = run(&pieces, &[]);
    let (_, b) = run(&flipped, &[]);
    assert_eq!(canon(&a), canon(&b));
}

#[test]
fn the_same_reads_give_the_same_contigs_every_time() {
    let genome = random(3000, 53);
    let pieces = reads(&genome, false, 100, 40.0, 0.005, 54);
    let (_, a) = run(&pieces, &[]);
    let (_, b) = run(&pieces, &[]);
    assert_eq!(a, b);
}

#[test]
fn every_link_of_the_graph_overlaps_by_k_minus_one_bases_in_the_direction_it_goes() {
    let (genome, _) = with_repeat(300);
    let pieces = reads(&genome, false, 100, 40.0, 0.0, 55);
    let (a, _) = run(&pieces, &[]);
    let k = a.graph.k();
    assert!(!a.graph.links().is_empty());
    for (x, y) in a.graph.links() {
        let (ox, oy) = (a.graph.oriented(x), a.graph.oriented(y));
        assert_eq!(ox[ox.len() - (k - 1)..], oy[..k - 1]);
    }
}

#[test]
fn a_branch_seen_too_few_times_is_cut_off_and_the_path_it_left_is_one_contig_again() {
    let genome = random(3000, 81);
    // a branch that leaves the genome at 1100 and goes on for 300 bases of something else, seen 1.6 times beside 40
    let branch: Vec<u8> = [&genome[1000..1100], &random(300, 82)[..]].concat();
    let mut pieces = tiled(&genome, 100, 2);
    pieces.extend(tiled(&branch, 100, 50));
    let a = assemble(&pieces, &options(&[21], Cutoff::Auto)).unwrap();
    assert_eq!(a.steps[0].low_coverage, 1);
    let c = number(&a.graph);
    assert_eq!(
        c.len(),
        1,
        "{:?}",
        c.iter().map(|c| c.seq.len()).collect::<Vec<_>>()
    );
    assert!(is_genome(&c[0].seq, &genome));
}

#[test]
fn a_circle_whose_two_copies_differ_by_one_base_is_one_circle_after_the_bulge_goes() {
    let a = random(1500, 83);
    let mut b = a.clone();
    b[700] = if b[700] == b'A' { b'C' } else { b'A' };
    let mut pieces = reads(&a, true, 100, 30.0, 0.0, 84);
    pieces.extend(reads(&b, true, 100, 30.0, 0.0, 85));
    let (asm, c) = run(&pieces, &[21]);
    assert_eq!(c.len(), 1);
    assert!(c[0].circular && c[0].seq.len() == 1500);
    assert!(belongs(&c[0].seq, &a, true) || belongs(&c[0].seq, &b, true));
    assert_eq!(asm.steps[0].bulges, 1);
}
