use super::*;
use crate::seq::variants::call::{Calls, Kind, Params, Variant};

fn variant(
    kind: Kind,
    pos: usize,
    reference: &str,
    alt: &str,
    quality: f64,
    counts: (usize, usize, usize),
) -> Variant {
    Variant {
        kind,
        pos,
        reference: reference.bytes().collect(),
        alt: alt.bytes().collect(),
        quality,
        depth: counts.0,
        ref_reads: counts.1,
        alt_reads: counts.2,
    }
}

fn calls() -> Calls {
    Calls {
        variants: vec![
            variant(Kind::Snp, 4, "A", "G", 212.508, (10, 0, 10)),
            variant(Kind::Deletion, 9, "CTT", "C", 36.34, (20, 6, 14)),
            variant(Kind::Insertion, 30, "G", "GAC", 55.0, (20, 5, 15)),
        ],
        consensus: (0..130).map(|i| b"ACGT"[i % 4]).collect(),
        depth: vec![7, 8, 9, 10],
        unknown: 2,
        mean_depth: 8.5,
    }
}

fn targets<'a>(c: &'a Calls) -> Vec<Target<'a>> {
    vec![Target {
        name: "plasmid",
        length: 130,
        calls: c,
    }]
}

#[test]
fn the_vcf_has_its_header_the_columns_and_a_line_to_each_difference() {
    let c = calls();
    let text = vcf(&targets(&c), "mysample", &Params::default());
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], "##fileformat=VCFv4.3");
    assert!(lines[1].starts_with("##source=helixloom "));
    assert!(lines.contains(&"##contig=<ID=plasmid,length=130>"));
    assert!(lines
        .iter()
        .any(|l| l.starts_with("##FORMAT=<ID=AD,Number=R")));
    assert!(lines
        .iter()
        .any(|l| l.starts_with("##INFO=<ID=INDEL,Number=0,Type=Flag")));
    assert!(lines
        .iter()
        .any(|l| l.contains("lowest quality 20") && l.contains("lowest depth 5")));
    let header = lines.iter().position(|l| l.starts_with("#CHROM")).unwrap();
    assert_eq!(
        lines[header],
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tmysample"
    );
    assert_eq!(
        &lines[header + 1..],
        [
            "plasmid\t5\t.\tA\tG\t212.5\tPASS\tDP=10\tGT:DP:AD\t1:10:0,10",
            "plasmid\t10\t.\tCTT\tC\t36.3\tPASS\tDP=20;INDEL\tGT:DP:AD\t1:20:6,14",
            "plasmid\t31\t.\tG\tGAC\t55.0\tPASS\tDP=20;INDEL\tGT:DP:AD\t1:20:5,15",
        ]
    );
}

#[test]
fn each_target_is_a_contig_and_its_differences_follow_in_the_order_of_the_targets() {
    let (a, b) = (
        calls(),
        Calls {
            variants: vec![variant(Kind::Snp, 0, "T", "C", 50.0, (9, 1, 8))],
            ..calls()
        },
    );
    let both = [
        Target {
            name: "one",
            length: 130,
            calls: &a,
        },
        Target {
            name: "two",
            length: 99,
            calls: &b,
        },
    ];
    let text = vcf(&both, "s", &Params::default());
    assert!(text.contains("##contig=<ID=one,length=130>\n##contig=<ID=two,length=99>\n"));
    let rows: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(rows.len(), 4);
    assert!(rows[0].starts_with("one\t5\t") && rows[3].starts_with("two\t1\t"));
}

#[test]
fn the_consensus_is_named_after_the_target_and_cut_into_lines_of_sixty() {
    let c = calls();
    let text = consensus(&targets(&c));
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines[0], ">plasmid consensus");
    assert_eq!(lines.len(), 4);
    assert_eq!(
        (lines[1].len(), lines[2].len(), lines[3].len()),
        (60, 60, 10)
    );
    assert_eq!(lines[1..].concat().as_bytes(), c.consensus.as_slice());
    // a consensus of exactly sixty is one line
    let sixty = Calls {
        consensus: vec![b'A'; 60],
        ..calls()
    };
    assert_eq!(consensus(&targets(&sixty)).lines().count(), 2);
}

#[test]
fn the_depth_table_counts_positions_from_one() {
    let c = calls();
    assert_eq!(
        depth_table(&targets(&c)),
        "target\tposition\tdepth\nplasmid\t1\t7\nplasmid\t2\t8\nplasmid\t3\t9\nplasmid\t4\t10\n"
    );
}

#[test]
fn the_summary_says_how_many_of_each_kind_and_how_much_was_not_known() {
    let c = calls();
    assert_eq!(
        summary(&targets(&c)),
        "plasmid: 130 bases, mean depth 8.5, 2 written as N, 1 substitutions, 1 insertions, 1 deletions\n"
    );
}
