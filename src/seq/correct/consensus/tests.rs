use super::*;
use crate::seq::variants::sam::parse_cigar;

fn record(pos: usize, cigar: &str, seq: &[u8]) -> Record {
    Record {
        name: "r".into(),
        flag: 0,
        target: Some("template".into()),
        pos,
        mapq: 60,
        cigar: parse_cigar(cigar).unwrap(),
        seq: seq.to_vec(),
        qual: Vec::new(),
        tags: Vec::new(),
    }
}

fn support(r: Record) -> Aligned {
    let length: usize = r
        .cigar
        .iter()
        .filter(|c| uses_target(c.0))
        .map(|c| c.1)
        .sum();
    Aligned {
        from: r.pos,
        to: r.pos + length,
        identity: 1.0,
        long: true,
        record: r,
    }
}

#[test]
fn a_part_of_a_record_keeps_the_bases_and_the_operations_that_lie_in_it() {
    let r = record(
        5,
        "2S4M2I3M2D3M",
        b"xxACGTttGGAcca".to_ascii_uppercase().as_slice(),
    );
    // the target positions: M 5..9, I after 8, M 9..12, D 12..14, M 14..17
    let c = clip(&r, 7, 15).unwrap();
    assert_eq!(c.pos, 7);
    assert_eq!(
        crate::seq::variants::sam::cigar_text(&c.cigar),
        "2M2I3M2D1M"
    );
    // the read bases that go with them: of the first block the last two, the insertion, the second block, one of the last
    let want: Vec<u8> = [&r.seq[4..6], &r.seq[6..8], &r.seq[8..11], &r.seq[11..12]].concat();
    assert_eq!(c.seq, want);
    // an insertion on the edge of the part is left out
    let edge = clip(&r, 9, 20).unwrap();
    assert_eq!(crate::seq::variants::sam::cigar_text(&edge.cigar), "3M2D3M");
    let edge = clip(&r, 0, 9).unwrap();
    assert_eq!(crate::seq::variants::sam::cigar_text(&edge.cigar), "4M");
    // nothing there
    assert!(clip(&r, 17, 30).is_none());
    assert!(clip(&r, 0, 5).is_none());
}

#[test]
fn the_qualities_go_with_the_bases() {
    let mut r = record(0, "6M", b"ACGTAC");
    r.qual = vec![1, 2, 3, 4, 5, 6];
    assert_eq!(clip(&r, 2, 4).unwrap().qual, [3, 4]);
}

/// A template of `blocks` blocks of 10 bases, and the template's bases with the bases of the chosen blocks changed, `wrong` of each ten.
fn blocks(n: usize, bad: &[usize], wrong: usize) -> (Vec<u8>, Record) {
    let template = crate::seq::minimap::sim::genome(10 * n, 31);
    let mut seq = template.clone();
    for &b in bad {
        for i in 0..wrong {
            let at = 10 * b + i;
            seq[at] = crate::seq::variants::sim::other(seq[at]);
        }
    }
    let r = record(0, &format!("{}M", 10 * n), &seq);
    (template, r)
}

#[test]
fn the_longest_run_of_blocks_with_half_the_errors_or_fewer_is_the_good_run() {
    let (t, r) = blocks(4, &[1], 6);
    assert_eq!(good_run(&r, &t, 10, 0.5), Some((20, 40)));
    // exactly half is still good
    let (t, r) = blocks(4, &[1], 5);
    assert_eq!(good_run(&r, &t, 10, 0.5), Some((0, 40)));
    // two runs of a size: the first
    let (t, r) = blocks(3, &[1], 6);
    assert_eq!(good_run(&r, &t, 10, 0.5), Some((0, 10)));
    // nothing good
    let (t, r) = blocks(2, &[0, 1], 6);
    assert_eq!(good_run(&r, &t, 10, 0.5), None);
}

#[test]
fn a_block_that_no_part_of_the_alignment_reaches_ends_a_run() {
    let template = crate::seq::minimap::sim::genome(45, 32);
    // aligned only on 0..10 and 20..45
    let mut seq = template[0..10].to_vec();
    seq.extend_from_slice(&template[20..45]);
    let r = record(0, "10M10N25M", &seq);
    assert_eq!(good_run(&r, &template, 10, 0.5), Some((20, 45)));
}

#[test]
fn deletions_and_insertions_count_as_errors_in_their_block() {
    let template = crate::seq::minimap::sim::genome(20, 33);
    // six bases of the second block deleted
    let mut seq = template[..10].to_vec();
    seq.extend_from_slice(&template[16..]);
    let r = record(0, "10M6D4M", &seq);
    assert_eq!(good_run(&r, &template, 10, 0.5), Some((0, 10)));
    // eleven inserted into the first block make it bad
    let mut seq = template[..5].to_vec();
    seq.extend(std::iter::repeat_n(b'A', 11));
    seq.extend_from_slice(&template[5..]);
    let r = record(0, "5M11I15M", &seq);
    assert_eq!(good_run(&r, &template, 10, 0.5), Some((10, 20)));
}

const T: &[u8] = b"ACGTACGTACGT";

fn voted(supports: &[Aligned], min: usize) -> Voted {
    let p = Params {
        min_support: min,
        ..Params::default()
    };
    vote(T, supports, false, &p)
}

fn reads(n: usize, pos: usize, cigar: &str, seq: &[u8]) -> Vec<Aligned> {
    (0..n).map(|_| support(record(pos, cigar, seq))).collect()
}

#[test]
fn a_position_with_fewer_supports_than_the_least_keeps_the_template() {
    let changed = b"ACGTACGAACGT";
    let three = voted(&reads(3, 0, "12M", changed), 4);
    assert_eq!(three.whole(), T);
    assert!(three.positions.iter().all(|p| !p.1));
    let four = voted(&reads(4, 0, "12M", changed), 4);
    assert_eq!(four.whole(), changed);
    assert!(four.positions.iter().all(|p| p.1));
}

#[test]
fn a_tie_between_the_template_and_one_read_keeps_the_template() {
    let changed = b"ACGTACGAACGT";
    let v = voted(&reads(1, 0, "12M", changed), 1);
    assert_eq!(v.whole(), T);
    let v = voted(&reads(2, 0, "12M", changed), 1);
    assert_eq!(v.whole(), changed);
    // and it keeps the template whichever of the four bases the template's is: here an A, the first, against a read's C
    let changed = b"ACGTCCGTACGT";
    let v = voted(&reads(1, 0, "12M", changed), 1);
    assert_eq!(v.whole(), T);
}

#[test]
fn a_position_is_left_out_when_more_reads_delete_it_than_have_a_base_there() {
    let without: Vec<u8> = [&T[..5], &T[6..]].concat();
    let mut five = reads(5, 0, "5M1D6M", &without);
    assert_eq!(voted(&five, 4).whole(), without);
    // three delete it, and two with the template have it: not more
    five.truncate(3);
    five.extend(reads(2, 0, "12M", T));
    assert_eq!(voted(&five, 4).whole(), T);
    // four delete it against two with the template
    five.extend(reads(1, 0, "5M1D6M", &without));
    assert_eq!(voted(&five, 4).whole(), without);
}

#[test]
fn an_insertion_is_written_when_over_half_of_those_with_a_base_there_have_it() {
    let with: Vec<u8> = [&T[..5], b"GG", &T[5..]].concat();
    let mut set = reads(4, 0, "5M2I7M", &with);
    set.extend(reads(2, 0, "12M", T));
    // seven with a base there, the template among them: four of them insert
    assert_eq!(voted(&set, 4).whole(), with);
    // exactly half of six do not
    let mut set = reads(3, 0, "5M2I7M", &with);
    set.extend(reads(2, 0, "12M", T));
    assert_eq!(voted(&set, 4).whole(), T);
    // three of seven do not
    let mut set = reads(3, 0, "5M2I7M", &with);
    set.extend(reads(3, 0, "12M", T));
    assert_eq!(voted(&set, 4).whole(), T);
}

#[test]
fn pieces_are_the_corrected_stretches_that_are_long_enough() {
    // reads over the first six and the last four positions; positions 6 and 7 have none
    let mut set = reads(4, 0, "6M", &T[..6]);
    set.extend(reads(4, 8, "4M", &T[8..]));
    let v = voted(&set, 4);
    assert_eq!(v.pieces(1), vec![T[..6].to_vec(), T[8..].to_vec()]);
    assert_eq!(v.pieces(5), vec![T[..6].to_vec()]);
    assert_eq!(v.pieces(6), vec![T[..6].to_vec()]);
    assert!(v.pieces(7).is_empty());
    assert_eq!(v.whole(), T);
}

#[test]
fn in_a_round_that_counts_blocks_an_alignment_counts_only_over_its_good_blocks() {
    let template = crate::seq::minimap::sim::genome(30, 34);
    // a support whose second block is all wrong
    let mut seq = template.clone();
    for base in &mut seq[10..20] {
        *base = crate::seq::variants::sim::other(*base);
    }
    let set: Vec<Aligned> = (0..4).map(|_| support(record(0, "30M", &seq))).collect();
    let p = Params {
        block: 10,
        min_support: 4,
        ..Params::default()
    };
    let blocks = vote(&template, &set, true, &p);
    // the good run is the first block (the first of two of a size): the second is untouched, and the third not counted
    assert_eq!(blocks.whole(), template);
    assert!(blocks.positions[..10].iter().all(|x| x.1));
    assert!(blocks.positions[10..].iter().all(|x| !x.1));
    // the whole alignments change the wrong block
    let all = vote(&template, &set, false, &p);
    assert_eq!(all.whole(), seq);
}

#[test]
fn with_no_supports_nothing_is_corrected() {
    let v = voted(&[], 1);
    assert_eq!(v.whole(), T);
    assert!(v.pieces(1).is_empty());
}
