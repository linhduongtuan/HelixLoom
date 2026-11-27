use super::*;
use crate::seq::count::annotation;
use crate::seq::count::assign::Tally;

fn ann() -> annotation::Annotation {
    annotation::parse(
        "GeneID\tChr\tStart\tEnd\tStrand\ng1\tc\t101\t200\t+\ng2\tc\t300\t310\t-\ng1\td\t151\t160\t.\n",
        &["exon".to_string()],
        "gene_id",
    )
    .unwrap()
}

fn tally(counts: &[f64]) -> Tally {
    let mut statuses = [0u64; 11];
    statuses[0] = 3;
    statuses[10] = 2;
    Tally {
        counts: counts.to_vec(),
        statuses,
    }
}

#[test]
fn the_table_has_a_row_to_a_gene_with_its_features_joined_and_its_length() {
    let a = ann();
    let index = Index::new(&a, false);
    let (one, two) = (tally(&[3.0, 0.0]), tally(&[1.0, 7.0]));
    let columns = [
        Column {
            name: "a.sam",
            tally: &one,
        },
        Column {
            name: "b.sam",
            tally: &two,
        },
    ];
    let text = table(&index, &columns, false, false);
    let lines: Vec<&str> = text.lines().collect();
    assert!(
        lines[0].starts_with("# Program: helixloom ") && lines[0].contains("after featureCounts")
    );
    assert_eq!(
        lines[1],
        "Geneid\tChr\tStart\tEnd\tStrand\tLength\ta.sam\tb.sam"
    );
    assert_eq!(lines[2], "g1\tc;d\t101;151\t200;160\t+;.\t110\t3\t1");
    assert_eq!(lines[3], "g2\tc\t300\t310\t-\t11\t0\t7");
    assert_eq!(lines.len(), 4);
}

#[test]
fn at_feature_level_each_feature_is_a_row_named_after_its_gene() {
    let a = ann();
    let index = Index::new(&a, true);
    let t = tally(&[1.0, 2.0, 3.0]);
    let text = table(
        &index,
        &[Column {
            name: "x",
            tally: &t,
        }],
        true,
        false,
    );
    let rows: Vec<&str> = text.lines().skip(2).collect();
    assert_eq!(
        rows,
        [
            "g1\tc\t101\t200\t+\t100\t1",
            "g2\tc\t300\t310\t-\t11\t2",
            "g1\td\t151\t160\t.\t10\t3",
        ]
    );
}

#[test]
fn fractions_are_written_with_two_places() {
    let a = ann();
    let index = Index::new(&a, false);
    let t = tally(&[0.125, 2.5]);
    let text = table(
        &index,
        &[Column {
            name: "x",
            tally: &t,
        }],
        false,
        true,
    );
    assert!(
        text.contains("\t110\t0.12\n") || text.contains("\t110\t0.13\n"),
        "{text}"
    );
    assert!(text.contains("\t11\t2.50\n"));
}

#[test]
fn the_summary_has_each_status_in_order_for_each_file() {
    let (one, two) = (tally(&[0.0, 0.0]), tally(&[0.0, 0.0]));
    let text = summary(&[
        Column {
            name: "a",
            tally: &one,
        },
        Column {
            name: "b",
            tally: &two,
        },
    ]);
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 12);
    assert_eq!(lines[0], "Status\ta\tb");
    assert_eq!(lines[1], "Assigned\t3\t3");
    assert_eq!(lines[2], "Unassigned_Unmapped\t0\t0");
    assert_eq!(lines[11], "Unassigned_Ambiguity\t2\t2");
}
