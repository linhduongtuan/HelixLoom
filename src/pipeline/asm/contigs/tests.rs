use super::*;
use crate::seq::debruijn::sim::random;

fn rows(v: Value) -> Vec<Named> {
    match v {
        Value::Sequences(r) => r,
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_contig_of_two_copies_is_cut_to_one_and_says_how_many() {
    let unit = random(1500, 1);
    let two = [&unit[..], &unit].concat();
    let (out, note) =
        concatemer(&[("a".into(), two), ("b".into(), random(2000, 2))], 95, 200).unwrap();
    let r = rows(out);
    assert_eq!(r[0], ("a concatemer=2".to_string(), unit));
    assert_eq!(r[1].0, "b");
    assert_eq!(r[1].1.len(), 2000);
    assert!(note.starts_with("1 of 2 contigs"), "{note}");
}

#[test]
fn the_identity_and_the_shortest_unit_are_the_ones_asked_for() {
    let unit = random(300, 3);
    let two = [&unit[..], &unit].concat();
    // a unit of 300 is under a shortest of 400
    let r = rows(concatemer(&[("a".into(), two.clone())], 95, 400).unwrap().0);
    assert_eq!(r[0].0, "a");
    let r = rows(concatemer(&[("a".into(), two)], 95, 200).unwrap().0);
    assert_eq!(r[0].0, "a concatemer=2");
}

#[test]
fn an_end_that_repeats_the_start_is_cut_and_the_contig_is_a_circle() {
    let circle = random(2000, 4);
    let contig = [&circle[..], &circle[..80]].concat();
    let (out, note) = circularize(
        &[("p".into(), contig), ("q".into(), random(1000, 5))],
        30,
        95,
    )
    .unwrap();
    let r = rows(out);
    assert_eq!(r[0], ("p circular=true".to_string(), circle));
    assert_eq!(r[1].0, "q");
    assert!(note.starts_with("1 of 2 contigs are circles"), "{note}");
}

#[test]
fn a_contig_that_is_a_circle_already_is_not_cut_again() {
    let circle = random(2000, 6);
    let contig = [&circle[..], &circle[..80]].concat();
    let r = rows(
        circularize(&[("p circular=true".into(), contig.clone())], 30, 95)
            .unwrap()
            .0,
    );
    assert_eq!(r[0], ("p circular=true".to_string(), contig));
}

#[test]
fn the_overlap_and_identity_asked_for_are_the_ones_used() {
    let circle = random(1000, 7);
    let contig = [&circle[..], &circle[..40]].concat();
    assert_eq!(
        rows(
            circularize(&[("p".into(), contig.clone())], 50, 95)
                .unwrap()
                .0
        )[0]
        .0,
        "p"
    );
    assert_eq!(
        rows(
            circularize(&[("p".into(), contig.clone())], 40, 95)
                .unwrap()
                .0
        )[0]
        .0,
        "p circular=true"
    );
    let mut worn = contig;
    let n = worn.len();
    for i in 0..8 {
        worn[n - 1 - i * 5] = if worn[n - 1 - i * 5] == b'A' {
            b'C'
        } else {
            b'A'
        };
    }
    // 8 of 40 changed: 80 % alike
    assert_eq!(
        rows(
            circularize(&[("p".into(), worn.clone())], 40, 95)
                .unwrap()
                .0
        )[0]
        .0,
        "p"
    );
    assert_eq!(
        rows(circularize(&[("p".into(), worn)], 40, 75).unwrap().0)[0].0,
        "p circular=true"
    );
}

#[test]
fn a_circle_is_turned_to_the_gene_and_a_line_is_not() {
    let d = crate::pipeline::asm::tests::folder("rotate");
    let circle = random(2500, 8);
    let genes = d.join("genes.fasta");
    std::fs::write(
        &genes,
        format!(">g\n{}\n", String::from_utf8_lossy(&circle[900..1300])),
    )
    .unwrap();
    let input = [
        ("c circular=true".to_string(), circle.clone()),
        ("l".to_string(), circle.clone()),
    ];
    let (out, note) = rotate(&input, Some(&genes)).unwrap();
    let r = rows(out);
    assert_eq!(r[0].0, "c circular=true rotated=true");
    assert_eq!(&r[0].1[..400], &circle[900..1300]);
    assert_eq!(r[1], ("l".to_string(), circle));
    assert!(note.starts_with("1 of 1 circles turned"), "{note}");
    // a gene that is not there leaves the circle as it was
    let other = d.join("other.fasta");
    std::fs::write(
        &other,
        format!(">g\n{}\n", String::from_utf8_lossy(&random(400, 99))),
    )
    .unwrap();
    let r = rows(rotate(&input[..1], Some(&other)).unwrap().0);
    assert_eq!(r[0].0, "c circular=true");
    assert!(rotate(&input, None).unwrap_err().contains("no genes"));
    std::fs::remove_dir_all(&d).ok();
}
