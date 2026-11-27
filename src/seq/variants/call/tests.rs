use super::*;
use crate::seq::variants::pileup::{pile, Filter};
use crate::seq::variants::sam::Record;
use crate::seq::variants::sim::{other, record};

const REF: &[u8] = b"ACGTAGCTAGGCTTACGATCGGATTCGACCATGGCTAAGC";

fn calls_with(reference: &[u8], records: &[Record], p: &Params) -> Calls {
    let piles = pile(
        &[("t".to_string(), reference.len())],
        records,
        &Filter::default(),
    );
    call_target(reference, &piles[0], p)
}

fn calls(records: &[Record]) -> Calls {
    calls_with(REF, records, &Params::default())
}

/// `n` identical reads from `pos` with the CIGAR and the bases, each base of quality `q`.
fn copies(n: usize, pos: usize, cigar: &str, seq: &[u8], q: u8) -> Vec<Record> {
    (0..n)
        .map(|i| {
            record(
                &format!("r{i}"),
                pos,
                cigar,
                std::str::from_utf8(seq).unwrap(),
                q,
            )
        })
        .collect()
}

fn full(n: usize, sample: &[u8], q: u8) -> Vec<Record> {
    copies(n, 0, &format!("{}M", sample.len()), sample, q)
}

fn changed(at: usize) -> Vec<u8> {
    let mut s = REF.to_vec();
    s[at] = other(s[at]);
    s
}

fn near(a: f64, b: f64, tolerance: f64) -> bool {
    (a - b).abs() <= tolerance
}

#[test]
fn the_defaults_are_the_ones_the_documentation_names() {
    let p = Params::default();
    assert_eq!(
        (
            p.min_depth,
            p.min_quality,
            p.snp_prior,
            p.indel_prior,
            p.indel_error
        ),
        (5, 20.0, 0.001, 0.0001, 0.1)
    );
    assert!(p.check().is_ok());
    assert_eq!(MAX_QUALITY, 999.0);
}

#[test]
fn a_parameter_out_of_its_range_is_refused() {
    let bad = |f: &dyn Fn(&mut Params)| {
        let mut p = Params::default();
        f(&mut p);
        p.check().unwrap_err()
    };
    assert!(bad(&|p| p.min_depth = 0).contains("depth"));
    assert!(bad(&|p| p.min_quality = -1.0).contains("quality"));
    assert!(bad(&|p| p.min_quality = f64::NAN).contains("quality"));
    assert!(bad(&|p| p.snp_prior = 0.0).contains("priors"));
    assert!(bad(&|p| p.snp_prior = 1.0).contains("priors"));
    assert!(bad(&|p| p.indel_prior = 0.0).contains("priors"));
    assert!(bad(&|p| p.indel_prior = 1.5).contains("priors"));
    assert!(bad(&|p| p.indel_error = 0.0).contains("wrong indel"));
    assert!(bad(&|p| p.indel_error = 0.5).contains("wrong indel"));
    // and the ends that are fine
    let ok = Params {
        min_quality: 0.0,
        indel_error: 0.49,
        ..Params::default()
    };
    assert!(ok.check().is_ok());
}

#[test]
fn a_base_that_every_read_shows_is_a_substitution_with_the_quality_that_the_model_gives() {
    let sample = changed(20);
    let c = calls(&full(10, &sample, 20));
    assert_eq!(c.variants.len(), 1);
    let v = &c.variants[0];
    assert_eq!(v.kind, Kind::Snp);
    assert_eq!(
        (v.pos, v.reference.as_slice(), v.alt.as_slice()),
        (20, &REF[20..21], &sample[20..21])
    );
    assert_eq!((v.depth, v.ref_reads, v.alt_reads), (10, 0, 10));
    // worked by hand: ln(0.001 / 3) + 10 ln 297 against ln(0.999), with the two other bases at ln(0.001 / 3)
    assert!(near(v.quality, 212.51, 0.05), "{}", v.quality);
    assert_eq!(c.consensus, sample);
    assert_eq!(c.unknown, 0);
}

#[test]
fn a_base_with_enough_reads_against_it_is_not_a_difference_and_one_with_enough_for_it_is() {
    let sample = changed(20);
    let mixed = |alt: usize, same: usize| {
        let mut reads = full(alt, &sample, 30);
        reads.extend(full(same, REF, 30));
        calls(&reads)
    };
    let c = mixed(27, 3);
    assert_eq!(c.variants.len(), 1);
    assert_eq!((c.variants[0].ref_reads, c.variants[0].alt_reads), (3, 27));
    // two in the same read set that do not agree: 3 of 30 is an error
    let c = mixed(3, 27);
    assert!(c.variants.is_empty());
    assert_eq!(c.consensus, REF);
}

#[test]
fn a_site_needs_the_least_depth_and_below_it_the_consensus_is_n() {
    let sample = changed(20);
    // exactly five reads: called
    let five = calls(&full(5, &sample, 30));
    assert_eq!(five.variants.len(), 1);
    assert_eq!(
        (five.unknown, five.consensus.as_slice()),
        (0, sample.as_slice())
    );
    // four: nothing is called and nothing is known
    let four = calls(&full(4, &sample, 30));
    assert!(four.variants.is_empty());
    assert_eq!(four.unknown, REF.len());
    assert!(four.consensus.iter().all(|&b| b == b'N'));
    // the least depth is a choice
    let p = Params {
        min_depth: 4,
        ..Params::default()
    };
    let c = calls_with(REF, &full(4, &sample, 30), &p);
    assert_eq!(c.variants.len(), 1);
}

#[test]
fn a_call_must_have_the_lowest_quality_and_not_less() {
    // two reads of quality 13 at one position: the base is the best, with a quality of 3.18 by the model
    let reads = copies(
        2,
        20,
        "1M",
        &REF[20..21].iter().map(|&b| other(b)).collect::<Vec<u8>>(),
        13,
    );
    let at = |q: f64| {
        calls_with(
            REF,
            &reads,
            &Params {
                min_depth: 1,
                min_quality: q,
                ..Params::default()
            },
        )
    };
    assert!(at(20.0).variants.is_empty());
    assert!(at(3.2).variants.is_empty());
    let v = at(3.1).variants;
    assert_eq!(v.len(), 1);
    assert!(near(v[0].quality, 3.179, 0.01), "{}", v[0].quality);
    // one read of quality 13 is not enough to move off the reference at all
    let one = copies(
        1,
        20,
        "1M",
        &REF[20..21].iter().map(|&b| other(b)).collect::<Vec<u8>>(),
        13,
    );
    let p = Params {
        min_depth: 1,
        min_quality: 0.0,
        ..Params::default()
    };
    assert!(calls_with(REF, &one, &p).variants.is_empty());
}

#[test]
fn the_prior_that_a_base_differs_moves_the_call() {
    let reads = copies(
        2,
        20,
        "1M",
        &REF[20..21].iter().map(|&b| other(b)).collect::<Vec<u8>>(),
        13,
    );
    let p = Params {
        min_depth: 1,
        min_quality: 3.0,
        snp_prior: 1e-6,
        ..Params::default()
    };
    // a much smaller prior: two reads are no longer enough
    assert!(calls_with(REF, &reads, &p).variants.is_empty());
    let p = Params {
        snp_prior: 0.1,
        ..p
    };
    assert_eq!(calls_with(REF, &reads, &p).variants.len(), 1);
}

#[test]
fn the_quality_is_held_at_the_highest_that_is_written() {
    let sample = changed(20);
    let c = calls(&full(250, &sample, 30));
    assert_eq!(c.variants[0].quality, MAX_QUALITY);
}

#[test]
fn a_reference_in_lower_case_is_the_same_and_one_with_n_has_the_base_that_the_reads_give() {
    let lower: Vec<u8> = REF.to_ascii_lowercase();
    let sample = changed(20);
    let c = calls_with(&lower, &full(10, &sample, 30), &Params::default());
    assert_eq!(c.variants.len(), 1);
    assert_eq!(c.variants[0].reference, &REF[20..21]);
    assert_eq!(c.consensus, sample);
    // an N in the reference is no difference, and the consensus has what the reads say
    let mut with_n = REF.to_vec();
    with_n[10] = b'N';
    let c = calls_with(&with_n, &full(10, REF, 30), &Params::default());
    assert!(c.variants.is_empty());
    assert_eq!(c.consensus, REF);
}

#[test]
fn a_run_of_deleted_positions_is_one_deletion_written_with_the_base_before_it() {
    let sample = [&REF[..10], &REF[13..]].concat();
    let reads = copies(30, 0, "10M3D27M", &sample, 30);
    let c = calls(&reads);
    assert_eq!(c.variants.len(), 1);
    let v = &c.variants[0];
    assert_eq!(v.kind, Kind::Deletion);
    assert_eq!(
        (v.pos, v.reference.as_slice(), v.alt.as_slice()),
        (9, &REF[9..13], &REF[9..10])
    );
    assert_eq!((v.depth, v.ref_reads, v.alt_reads), (30, 0, 30));
    assert_eq!(c.consensus, sample);
}

#[test]
fn a_deletion_at_the_very_start_is_written_with_the_base_after_it() {
    let reads = copies(30, 0, "2D38M", &REF[2..], 30);
    let c = calls(&reads);
    assert_eq!(c.variants.len(), 1);
    let v = &c.variants[0];
    assert_eq!(
        (v.pos, v.reference.as_slice(), v.alt.as_slice()),
        (0, &REF[..3], &REF[2..3])
    );
    assert_eq!(c.consensus, &REF[2..]);
}

#[test]
fn a_deletion_that_few_reads_show_is_not_called_and_where_the_line_is_depends_on_the_reads_that_do_not(
) {
    let sample = [&REF[..10], &REF[12..]].concat();
    let mixed = |with: usize, without: usize| {
        let mut reads = copies(with, 0, "10M2D28M", &sample, 30);
        reads.extend(full(without, REF, 30));
        calls(&reads)
    };
    assert!(mixed(3, 27).variants.is_empty());
    // of twenty reads, 13 are not enough and 14 are: the quality is 17.3 and 36.3 by the model
    assert!(mixed(13, 7).variants.is_empty());
    let c = mixed(14, 6);
    assert_eq!(c.variants.len(), 1);
    let v = &c.variants[0];
    assert!(near(v.quality, 36.34, 0.05), "{}", v.quality);
    assert_eq!((v.depth, v.ref_reads, v.alt_reads), (20, 6, 14));
    // the consensus has it
    assert_eq!(c.consensus, sample);
    // a lower chance that a read errs makes 13 enough, a higher one does not make 14 so
    let p = Params {
        indel_error: 0.3,
        ..Params::default()
    };
    let mut reads = copies(14, 0, "10M2D28M", &sample, 30);
    reads.extend(full(6, REF, 30));
    assert!(calls_with(REF, &reads, &p).variants.is_empty());
}

#[test]
fn an_indel_needs_the_least_depth_over_the_place() {
    // five reads that all show it: the quality is 8.4 by the model, under the lowest; with the lowest lowered it is called at five and not at four
    let sample = [&REF[..10], &REF[12..]].concat();
    let p = Params {
        min_quality: 5.0,
        ..Params::default()
    };
    let five = calls_with(REF, &copies(5, 0, "10M2D28M", &sample, 30), &p);
    assert_eq!(five.variants.len(), 1);
    assert!(near(five.variants[0].quality, 8.39, 0.05));
    let four = calls_with(REF, &copies(4, 0, "10M2D28M", &sample, 30), &p);
    assert!(four.variants.is_empty());
    // and under the default quality, five are not enough
    assert!(calls(&copies(5, 0, "10M2D28M", &sample, 30))
        .variants
        .is_empty());
}

#[test]
fn a_run_takes_the_numbers_of_its_weakest_position() {
    // 15 reads delete positions 10 to 12, and 5 delete only 10 and 11
    let long = [&REF[..10], &REF[13..]].concat();
    let short = [&REF[..10], &REF[12..]].concat();
    let mut reads = copies(15, 0, "10M3D27M", &long, 30);
    reads.extend(copies(5, 0, "10M2D28M", &short, 30));
    let c = calls(&reads);
    assert_eq!(c.variants.len(), 1);
    let v = &c.variants[0];
    assert_eq!(v.kind, Kind::Deletion);
    assert_eq!(v.reference, &REF[9..13]);
    assert!(near(v.quality, 55.42, 0.05), "{}", v.quality);
    assert_eq!((v.depth, v.ref_reads, v.alt_reads), (20, 5, 15));
}

#[test]
fn two_deletions_apart_are_two_variants() {
    let sample = [&REF[..5], &REF[7..17], &REF[19..]].concat();
    let c = calls(&copies(30, 0, "5M2D10M2D21M", &sample, 30));
    let places: Vec<(Kind, usize)> = c.variants.iter().map(|v| (v.kind, v.pos)).collect();
    assert_eq!(places, [(Kind::Deletion, 4), (Kind::Deletion, 16)]);
    assert_eq!(c.consensus, sample);
}

#[test]
fn an_insertion_that_the_reads_show_is_written_after_the_base_before_it() {
    let sample = [&REF[..10], b"GT", &REF[10..]].concat();
    let c = calls(&copies(25, 0, "10M2I30M", &sample, 30));
    assert_eq!(c.variants.len(), 1);
    let v = &c.variants[0];
    assert_eq!(v.kind, Kind::Insertion);
    assert_eq!(v.pos, 9);
    assert_eq!(v.reference, &REF[9..10]);
    assert_eq!(v.alt, [&REF[9..10], b"GT".as_slice()].concat());
    assert_eq!((v.depth, v.ref_reads, v.alt_reads), (25, 0, 25));
    assert_eq!(c.consensus, sample);
}

#[test]
fn of_the_insertions_after_a_position_the_commonest_is_taken_and_the_others_count_against_it() {
    let with = |a: usize, b: usize| {
        let gt = [&REF[..10], b"GT", &REF[10..]].concat();
        let ac = [&REF[..10], b"AC", &REF[10..]].concat();
        let mut reads = copies(a, 0, "10M2I30M", &gt, 30);
        reads.extend(copies(b, 0, "10M2I30M", &ac, 30));
        calls(&reads)
    };
    let c = with(15, 5);
    assert_eq!(c.variants.len(), 1);
    assert_eq!(c.variants[0].alt, [&REF[9..10], b"GT".as_slice()].concat());
    assert_eq!(c.variants[0].alt_reads, 15);
    assert_eq!(c.variants[0].ref_reads, 5);
    // half and half is not a call, and neither is a minority
    assert!(with(10, 10).variants.is_empty());
    assert!(with(3, 0).variants.is_empty() || with(3, 0).variants[0].alt_reads == 3);
}

#[test]
fn an_insertion_of_bases_that_are_not_a_c_g_or_t_is_not_called() {
    let sample = [&REF[..10], b"NN", &REF[10..]].concat();
    let c = calls(&copies(25, 0, "10M2I30M", &sample, 30));
    assert!(c.variants.is_empty());
    assert_eq!(c.consensus, REF);
}

#[test]
fn an_insertion_after_a_position_that_is_deleted_is_left_out() {
    let deleted = [&REF[..10], &REF[13..]].concat();
    let inserted = [&REF[..12], b"G", &REF[12..]].concat();
    let mut reads = copies(30, 0, "10M3D27M", &deleted, 30);
    reads.extend(copies(12, 0, "12M1I28M", &inserted, 30));
    let c = calls(&reads);
    let kinds: Vec<Kind> = c.variants.iter().map(|v| v.kind).collect();
    assert_eq!(kinds, [Kind::Deletion]);
    assert_eq!(c.consensus, deleted);
}

#[test]
fn the_variants_are_in_the_order_of_the_target_and_each_kind_is_found_in_one_run() {
    let mut sample = REF.to_vec();
    sample[5] = other(sample[5]);
    let sample = [&sample[..20], &sample[22..30], b"A", &sample[30..]].concat();
    let c = calls(&copies(30, 0, "20M2D8M1I10M", &sample, 30));
    let found: Vec<(Kind, usize)> = c.variants.iter().map(|v| (v.kind, v.pos)).collect();
    assert_eq!(
        found,
        [(Kind::Snp, 5), (Kind::Deletion, 19), (Kind::Insertion, 29)]
    );
}

#[test]
fn the_depth_of_each_position_and_its_mean_are_reported() {
    let c = calls(&copies(10, 0, "20M", &REF[..20], 30));
    assert_eq!(c.depth.len(), 40);
    assert!(c.depth[..20].iter().all(|&d| d == 10));
    assert!(c.depth[20..].iter().all(|&d| d == 0));
    assert!(near(c.mean_depth, 5.0, 1e-9));
    // the positions that no read reaches are N
    assert_eq!(c.unknown, 20);
    assert_eq!(&c.consensus[..20], &REF[..20]);
    assert!(c.consensus[20..].iter().all(|&b| b == b'N'));
    let none = call_target(b"", &Pileup::new(0), &Params::default());
    assert_eq!((none.consensus.len(), none.mean_depth), (0, 0.0));
}

#[test]
fn an_indel_with_overwhelming_support_has_the_highest_quality_that_is_written() {
    let sample = [&REF[..10], &REF[12..]].concat();
    let c = calls(&copies(250, 0, "10M2D28M", &sample, 30));
    assert_eq!(c.variants.len(), 1);
    assert_eq!(c.variants[0].quality, MAX_QUALITY);
}

#[test]
fn of_two_other_bases_that_the_reads_show_equally_the_earlier_in_the_order_acgt_is_taken() {
    // the reference has A at 0; five reads show C and five show G, with the same qualities
    let mut reads = copies(5, 0, "1M", b"C", 30);
    reads.extend(copies(5, 0, "1M", b"G", 30));
    let c = calls_with(
        REF,
        &reads,
        &Params {
            min_depth: 1,
            ..Params::default()
        },
    );
    let first = c
        .variants
        .iter()
        .find(|v| v.pos == 0)
        .expect("a call at the first position");
    assert_eq!(first.alt, b"C");
}

#[test]
fn an_indel_whose_quality_is_exactly_the_lowest_is_called_and_one_under_it_is_not() {
    let p = Params::default();
    let q = indel(14, 20, &p).unwrap();
    assert!(near(q, 36.34, 0.05));
    let exactly = Params {
        min_quality: q,
        ..Params::default()
    };
    assert_eq!(indel(14, 20, &exactly), Some(q));
    let above = Params {
        min_quality: q + 1e-9,
        ..Params::default()
    };
    assert_eq!(indel(14, 20, &above), None);
}

#[test]
fn a_substitution_whose_quality_is_exactly_the_lowest_is_called() {
    let sample = changed(20);
    let reads = full(10, &sample, 20);
    let q = calls(&reads).variants[0].quality;
    let exactly = Params {
        min_quality: q,
        ..Params::default()
    };
    assert_eq!(calls_with(REF, &reads, &exactly).variants.len(), 1);
    let above = Params {
        min_quality: q + 1e-6,
        ..Params::default()
    };
    assert!(calls_with(REF, &reads, &above).variants.is_empty());
}

#[test]
fn an_insertion_with_one_base_that_is_not_a_c_g_or_t_among_good_ones_is_not_called() {
    let sample = [&REF[..10], b"AN", &REF[10..]].concat();
    let c = calls(&copies(25, 0, "10M2I30M", &sample, 30));
    assert!(c.variants.is_empty());
    assert_eq!(c.consensus, REF);
}

#[test]
fn the_depth_of_an_insertion_is_the_reads_that_have_a_base_before_it_and_not_those_that_delete_it()
{
    // 25 reads have the insertion after position 9, and 3 delete position 9 (too few to be called)
    let with = [&REF[..10], b"GT", &REF[10..]].concat();
    let without = [&REF[..9], &REF[10..]].concat();
    let mut reads = copies(25, 0, "10M2I30M", &with, 30);
    reads.extend(copies(3, 0, "9M1D30M", &without, 30));
    let c = calls(&reads);
    let ins = c
        .variants
        .iter()
        .find(|v| v.kind == Kind::Insertion)
        .expect("the insertion");
    assert_eq!((ins.depth, ins.alt_reads, ins.ref_reads), (25, 25, 0));
}

#[test]
fn a_pileup_longer_or_shorter_than_the_reference_is_used_as_far_as_both_go() {
    let sample = changed(20);
    // the reads cover 30 bases of a 40-base reference: the pileup is as long as the reference here
    let piles = pile(
        &[("t".to_string(), 30)],
        &full_to(30, &sample),
        &Filter::default(),
    );
    let c = call_target(REF, &piles[0], &Params::default());
    assert_eq!(c.consensus.len(), 30);
    assert_eq!(c.depth.len(), 30);
    // a reference shorter than the pileup: the reference's length
    let piles = pile(
        &[("t".to_string(), 40)],
        &full(10, REF, 30),
        &Filter::default(),
    );
    let c = call_target(&REF[..25], &piles[0], &Params::default());
    assert_eq!(c.consensus.len(), 25);
}

fn full_to(n: usize, sample: &[u8]) -> Vec<Record> {
    copies(10, 0, &format!("{n}M"), &sample[..n], 30)
}

#[test]
fn a_deletion_that_runs_from_the_start_to_the_end_is_not_written_and_does_not_stop_the_calls() {
    // reads that delete the whole target and have only soft-clipped bases
    let c = calls(&copies(30, 0, "40D2S", b"AC", 30));
    assert!(c.variants.is_empty());
    assert_eq!(c.consensus.len(), REF.len());
    assert!(c.consensus.iter().all(|&b| b == b'N'));
}

#[test]
fn the_depth_of_a_position_counts_the_reads_over_it_whatever_their_bases_and_the_quality_of_them() {
    // ten reads over positions 0 to 9, one of them with quality 5 (not counted as a base), and five that delete position 9
    let mut reads = copies(9, 0, "10M", &REF[..10], 30);
    reads.extend(copies(1, 0, "10M", &REF[..10], 5));
    reads.extend(copies(5, 0, "9M1D", &REF[..9], 30));
    let c = calls(&reads);
    assert_eq!(&c.depth[..9], &[15; 9]);
    assert_eq!(c.depth[9], 15);
}

#[test]
fn where_the_reference_is_unknown_a_split_vote_is_written_as_n() {
    let mut with_n = REF.to_vec();
    with_n[10] = b'N';
    // three reads show A there and three C, with the same qualities: no base is ahead
    let mut a = REF.to_vec();
    a[10] = b'A';
    let mut c = REF.to_vec();
    c[10] = b'C';
    let mut reads = full(3, &a, 20);
    reads.extend(full(3, &c, 20));
    let calls = calls_with(&with_n, &reads, &Params::default());
    assert_eq!(calls.consensus[10], b'N');
    assert_eq!(calls.unknown, 1);
    // five to one is far ahead: the base is written
    let mut reads = full(5, &c, 20);
    reads.extend(full(1, &a, 20));
    let calls = calls_with(&with_n, &reads, &Params::default());
    assert_eq!((calls.consensus[10], calls.unknown), (b'C', 0));
    // three reads against two at quality 13 are ahead by one read's weight, a lead of 17.5: under the lowest quality of 20
    let mut reads = full(3, &c, 13);
    reads.extend(full(2, &a, 13));
    let calls = calls_with(&with_n, &reads, &Params::default());
    assert_eq!(calls.consensus[10], b'N');
    let p = Params {
        min_quality: 17.0,
        ..Params::default()
    };
    assert_eq!(calls_with(&with_n, &reads, &p).consensus[10], b'C');
}
