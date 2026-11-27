use super::*;
use crate::seq::minimap::sim::{genome, read, Rng};
use crate::seq::minimap::{map_read, Index, Options};
use crate::seq::variants::sim::record;

#[test]
fn a_cigar_is_read_in_runs_and_written_back_and_star_is_none() {
    let c = parse_cigar("10M2I3D5S1H4N6P7=8X").unwrap();
    assert_eq!(
        c,
        [
            ('M', 10),
            ('I', 2),
            ('D', 3),
            ('S', 5),
            ('H', 1),
            ('N', 4),
            ('P', 6),
            ('=', 7),
            ('X', 8)
        ]
    );
    assert_eq!(cigar_text(&c), "10M2I3D5S1H4N6P7=8X");
    assert_eq!(parse_cigar("*").unwrap(), Vec::new());
    assert_eq!(cigar_text(&Vec::new()), "*");
    // an operation of no length is no operation
    assert_eq!(parse_cigar("0M5I").unwrap(), [('I', 5)]);
}

#[test]
fn a_cigar_that_is_not_runs_of_a_number_and_an_operation_is_refused() {
    for bad in ["", "10", "M", "10Z", "5M3", "M5", "99999999999999999999M"] {
        assert!(parse_cigar(bad).is_err(), "{bad}");
    }
    assert!(parse_cigar("10Z").unwrap_err().contains("not an operation"));
    assert!(parse_cigar("M5").unwrap_err().contains("no number"));
}

#[test]
fn the_operations_know_which_of_the_read_and_the_target_they_use_up() {
    let read: String = "MISXN=DHP".chars().filter(|&c| uses_read(c)).collect();
    let target: String = "MISXN=DHP".chars().filter(|&c| uses_target(c)).collect();
    assert_eq!(read, "MISX=");
    assert_eq!(target, "MXN=D");
}

#[test]
fn a_record_is_written_with_its_position_from_one_its_qualities_in_phred_plus_33_and_its_tags() {
    let mut r = record("r1", 9, "4M", "ACGT", 30);
    r.mapq = 42;
    r.tags = vec!["NM:i:0".into(), "tp:A:P".into()];
    assert_eq!(
        r.line(),
        "r1\t0\tt\t10\t42\t4M\t*\t0\t0\tACGT\t????\tNM:i:0\ttp:A:P"
    );
    // a quality above 93 is held there (the last printable character)
    r.qual = vec![200; 4];
    assert!(r.line().contains("\t~~~~\t"));
}

#[test]
fn a_read_that_did_not_map_is_one_line_with_stars_and_flag_four() {
    let read = Read {
        name: "lost".into(),
        seq: b"ACGT".to_vec(),
        qual: Vec::new(),
    };
    let u = Record::unmapped(&read);
    assert_eq!(u.line(), "lost\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*");
    assert!(!u.mapped() && !u.usable());
}

#[test]
fn only_a_mapped_primary_that_passed_its_checks_is_usable() {
    let r = record("r", 0, "4M", "ACGT", 30);
    assert!(r.mapped() && r.usable());
    for flag in [UNMAPPED, SECONDARY, DUPLICATE, QC_FAIL] {
        let mut other = r.clone();
        other.flag = flag;
        assert!(!other.usable(), "flag {flag}");
    }
    // the reverse strand is usable, and so is the other part of a split read
    for flag in [REVERSE, SUPPLEMENTARY] {
        let mut other = r.clone();
        other.flag = flag;
        assert!(other.usable(), "flag {flag}");
    }
    // no bases, or no CIGAR, or no target
    let mut no_seq = r.clone();
    no_seq.seq.clear();
    assert!(!no_seq.usable() && no_seq.mapped());
    let mut no_cigar = r.clone();
    no_cigar.cigar.clear();
    assert!(!no_cigar.mapped());
    let mut no_target = r;
    no_target.target = None;
    assert!(!no_target.mapped());
}

#[test]
fn a_sam_file_is_read_with_its_targets_and_written_back_the_same() {
    let text = "@HD\tVN:1.6\n@SQ\tSN:t\tLN:50\n@SQ\tSN:u\tLN:20\n\
                a\t0\tt\t5\t60\t4M\t*\t0\t0\tACGT\t!!II\tNM:i:1\n\
                b\t16\tu\t1\t7\t2M1D2M\t*\t0\t0\tAcGT\t*\n\
                c\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\t*\n";
    let sam = parse(text).unwrap();
    assert_eq!(sam.targets, [("t".into(), 50), ("u".into(), 20)]);
    assert_eq!(sam.records.len(), 3);
    let a = &sam.records[0];
    assert_eq!((a.pos, a.mapq, a.seq.as_slice()), (4, 60, &b"ACGT"[..]));
    assert_eq!(a.qual, [0, 0, 40, 40]);
    assert_eq!(a.tags, ["NM:i:1"]);
    // lower case is upper case, and * for qualities is none
    let b = &sam.records[1];
    assert_eq!(
        (b.flag, b.seq.as_slice(), b.qual.len()),
        (16, &b"ACGT"[..], 0)
    );
    assert!(!sam.records[2].mapped());
    // written and read again, the records are the same
    let again = parse(&write(&sam)).unwrap();
    assert_eq!(again.records, sam.records);
    assert_eq!(again.targets, sam.targets);
}

#[test]
fn a_header_names_the_targets_and_this_program() {
    let h = header(&[("a".into(), 10), ("b".into(), 20)]);
    let lines: Vec<&str> = h.lines().collect();
    assert_eq!(lines[0], "@HD\tVN:1.6\tSO:unsorted");
    assert_eq!(lines[1], "@SQ\tSN:a\tLN:10");
    assert_eq!(lines[2], "@SQ\tSN:b\tLN:20");
    assert!(lines[3].starts_with("@PG\tID:helixloom\tPN:helixloom\tVN:"));
    assert_eq!(lines.len(), 4);
}

#[test]
fn a_sam_with_a_fault_is_refused_with_the_line() {
    let two = parse("a\t0\tt\t1\t5\t4M\t*\t0\t0\tACGT\n").unwrap_err();
    assert!(
        two.contains("Line 1") && two.contains("10 columns"),
        "{two}"
    );
    let flag = parse("a\tx\tt\t1\t5\t4M\t*\t0\t0\tACGT\t*\n").unwrap_err();
    assert!(flag.contains("flag"), "{flag}");
    let pos = parse("a\t0\tt\t-1\t5\t4M\t*\t0\t0\tACGT\t*\n").unwrap_err();
    assert!(pos.contains("position"), "{pos}");
    let mapq = parse("a\t0\tt\t1\t300\t4M\t*\t0\t0\tACGT\t*\n").unwrap_err();
    assert!(mapq.contains("mapping quality"), "{mapq}");
    let cigar =
        parse("ok\t0\tt\t1\t5\t4M\t*\t0\t0\tACGT\t*\na\t0\tt\t1\t5\t5M\t*\t0\t0\tACGT\t*\n")
            .unwrap_err();
    assert!(
        cigar.contains("Line 2") && cigar.contains("5 bases"),
        "{cigar}"
    );
    let qual = parse("a\t0\tt\t1\t5\t4M\t*\t0\t0\tACGT\t!!!\n").unwrap_err();
    assert!(qual.contains("4 bases and 3 qualities"), "{qual}");
    let low = parse("a\t0\tt\t1\t5\t4M\t*\t0\t0\tACGT\t!! !\n").unwrap_err();
    assert!(low.contains("quality"), "{low}");
    let bad_cigar = parse("a\t0\tt\t1\t5\t4Q\t*\t0\t0\tACGT\t*\n").unwrap_err();
    assert!(bad_cigar.contains("Line 1") && bad_cigar.contains("not an operation"));
}

#[test]
fn the_header_lines_of_targets_need_a_name_and_a_length_and_the_reads_must_lie_on_them() {
    assert!(parse("@SQ\tLN:5\n").unwrap_err().contains("no SN"));
    assert!(parse("@SQ\tSN:t\n").unwrap_err().contains("no LN"));
    assert!(parse("@SQ\tSN:t\tLN:x\n").unwrap_err().contains("no LN"));
    let off = parse("@SQ\tSN:t\tLN:9\na\t0\tz\t1\t5\t4M\t*\t0\t0\tACGT\t*\n").unwrap_err();
    assert!(off.contains("target z") && off.contains("no @SQ"), "{off}");
    // with no @SQ line at all, any target will do
    assert!(parse("a\t0\tz\t1\t5\t4M\t*\t0\t0\tACGT\t*\n").is_ok());
    // blank lines and CRLF are let through, other header lines too
    let sam = parse("@CO\tnote\r\n\r\na\t0\tt\t1\t5\t4M\t*\t0\t0\tACGT\t*\r\n").unwrap();
    assert_eq!(sam.records.len(), 1);
    assert_eq!(sam.records[0].tags.len(), 0);
}

fn exact(rev: bool) -> (Vec<u8>, Read, Record) {
    let g = genome(6000, 91);
    // a read of two thousand bases of the genome with 70 bases on each side that belong to nothing
    let mut r = Rng::new(3);
    let flank =
        |n: usize, r: &mut Rng| -> Vec<u8> { (0..n).map(|_| b"ACGT"[r.below(4)]).collect() };
    let mut seq = flank(70, &mut r);
    seq.extend_from_slice(&read(&g, 1000, 2000, false, 0.0, &mut r));
    seq.extend(flank(70, &mut r));
    let seq = if rev {
        crate::seq::dna::revcom(&seq)
    } else {
        seq
    };
    let quals: Vec<u8> = (0..seq.len()).map(|i| (i % 40) as u8).collect();
    let o = Options {
        cigar: true,
        ..Options::default()
    };
    let index = Index::build(vec![("ref".to_string(), g.clone())], &o).unwrap();
    let rd = Read {
        name: "r".into(),
        seq: seq.clone(),
        qual: quals,
    };
    let hits = map_read(&index, &seq, &o, &|_| false);
    let rec = Record::from_hit(&rd, &index, &hits[0]).unwrap();
    (g, rd, rec)
}

#[test]
fn a_read_that_maps_forward_has_its_ends_that_were_not_aligned_as_soft_clips() {
    let (g, rd, rec) = exact(false);
    assert_eq!((rec.flag, rec.target.as_deref()), (0, Some("ref")));
    // the clips and the aligned bases account for every base of the read
    let used: usize = rec.cigar.iter().map(|c| c.1).sum();
    assert_eq!(used, rd.seq.len());
    let (first, last) = (rec.cigar[0], rec.cigar[rec.cigar.len() - 1]);
    assert!(first.0 == 'S' && first.1 >= 70, "{first:?}");
    assert!(last.0 == 'S' && last.1 >= 70, "{last:?}");
    // the read is as it was, and the bases after the clip are the genome's from the position
    assert_eq!(rec.seq, rd.seq);
    assert_eq!(rec.qual, rd.qual);
    let n = rec.cigar[1].1;
    assert_eq!(rec.cigar[1].0, 'M');
    assert_eq!(&rec.seq[first.1..first.1 + n], &g[rec.pos..rec.pos + n]);
    assert!(rec.usable());
    // the read has no errors: nothing differs
    assert_eq!(rec.tags[0], "NM:i:0");
}

#[test]
fn nm_counts_the_bases_that_differ_from_the_target() {
    let g = genome(6000, 92);
    let mut r = Rng::new(5);
    let seq = read(&g, 1000, 2500, false, 0.04, &mut r);
    let o = Options {
        cigar: true,
        ..Options::default()
    };
    let index = Index::build(vec![("ref".to_string(), g)], &o).unwrap();
    let rd = Read {
        name: "n".into(),
        seq: seq.clone(),
        qual: Vec::new(),
    };
    let hits = map_read(&index, &seq, &o, &|_| false);
    let rec = Record::from_hit(&rd, &index, &hits[0]).unwrap();
    let base = hits[0].base.as_ref().unwrap();
    assert_eq!(rec.tags[0], format!("NM:i:{}", base.differences));
    // a read with 12 % of its bases wrong differs in a good many, and in less than the bases it holds
    assert!(
        base.differences > 100 && base.differences < base.matches,
        "{} {}",
        base.differences,
        base.matches
    );
}

#[test]
fn a_secondary_on_the_other_strand_has_both_flags() {
    let g = genome(2000, 8);
    let o = Options {
        cigar: true,
        ..Options::default()
    };
    let index = Index::build(vec![("ref".to_string(), g.clone())], &o).unwrap();
    let seq = crate::seq::dna::revcom(&g[100..900]);
    let mut hits = map_read(&index, &seq, &o, &|_| false);
    hits[0].primary = false;
    let rd = Read {
        name: "x".into(),
        seq,
        qual: Vec::new(),
    };
    let rec = Record::from_hit(&rd, &index, &hits[0]).unwrap();
    assert_eq!(rec.flag, REVERSE | SECONDARY);
}

#[test]
fn a_read_that_maps_on_the_other_strand_is_written_as_its_reverse_complement_with_the_qualities_turned(
) {
    let (g, rd, rec) = exact(true);
    assert_eq!(rec.flag, REVERSE);
    assert_eq!(rec.seq, crate::seq::dna::revcom(&rd.seq));
    let turned: Vec<u8> = rd.qual.iter().rev().copied().collect();
    assert_eq!(rec.qual, turned);
    let used: usize = rec.cigar.iter().map(|c| c.1).sum();
    assert_eq!(used, rd.seq.len());
    let first = rec.cigar[0];
    let n = rec.cigar[1].1;
    assert!(first.0 == 'S' && first.1 >= 70);
    // against the genome the written bases are the genome's
    assert_eq!(&rec.seq[first.1..first.1 + n], &g[rec.pos..rec.pos + n]);
}

#[test]
fn a_secondary_mapping_has_flag_256_and_no_bases() {
    let (_, rd, rec) = exact(false);
    let g = genome(2000, 5);
    let o = Options {
        cigar: true,
        ..Options::default()
    };
    let index = Index::build(vec![("ref".to_string(), g.clone())], &o).unwrap();
    let mut hits = map_read(&index, &g[100..900], &o, &|_| false);
    assert!(!hits.is_empty());
    hits[0].primary = false;
    let sec = Record::from_hit(&rd, &index, &hits[0]).unwrap();
    assert_eq!(sec.flag & SECONDARY, SECONDARY);
    assert!(sec.seq.is_empty() && sec.qual.is_empty());
    assert!(sec.line().contains("\t*\t*\t") || sec.line().contains("\t*\ttp:A:S"));
    assert!(!sec.usable());
    assert!(rec.usable());
}

#[test]
fn a_hit_with_no_alignment_gives_no_record() {
    let g = genome(2000, 6);
    let o = Options {
        cigar: false,
        ..Options::default()
    };
    let index = Index::build(vec![("ref".to_string(), g.clone())], &o).unwrap();
    let hits = map_read(&index, &g[100..900], &o, &|_| false);
    let rd = Read {
        name: "x".into(),
        seq: g[100..900].to_vec(),
        qual: Vec::new(),
    };
    assert!(Record::from_hit(&rd, &index, &hits[0]).is_none());
}

#[test]
fn a_line_with_position_zero_lies_on_no_sequence() {
    let sam = parse("a\t0\tt\t0\t60\t4M\t*\t0\t0\tACGT\t*\n").unwrap();
    assert_eq!(sam.records[0].target, None);
    assert!(!sam.records[0].mapped());
}
