use super::*;
use crate::pipeline::asm::tests::{fastq_reads, folder};
use crate::seq::debruijn::sim::{random, reads as short_reads};
use crate::seq::hybrid::sim::{is_circle_of, long_reads, ring, short_graph};

#[test]
fn the_list_of_k_is_odd_rising_and_in_range() {
    assert_eq!(parse_ks("21,33,55").unwrap(), [21, 33, 55]);
    assert_eq!(parse_ks("21").unwrap(), [21]);
    assert_eq!(parse_ks(" 21 , 33 ").unwrap(), [21, 33]);
    assert!(parse_ks("22").unwrap_err().contains("odd"));
    assert!(parse_ks("65").is_err());
    assert!(parse_ks("33,21").unwrap_err().contains("must rise"));
    assert!(parse_ks("21,21").unwrap_err().contains("none twice"));
    assert!(parse_ks("a").unwrap_err().contains("list of numbers"));
    assert!(parse_ks("").is_err());
}

#[test]
fn short_reads_become_a_graph_that_the_next_step_can_read() {
    let genome = random(3000, 1);
    let reads = fastq_reads(&short_reads(&genome, false, 100, 40.0, 0.0, 2));
    let (value, note) = assemble(&reads, "auto", "auto", 10).unwrap();
    let Value::Gfa(text) = value else { panic!() };
    let graph = input::read(&text, 0).unwrap().graph;
    assert_eq!(graph.count(), 1);
    assert_eq!(graph.length(graph.ids()[0]), 3000);
    assert!(note.contains("1 segments, 3000 bases, k = 55"), "{note}");
}

#[test]
fn the_k_that_are_asked_for_are_the_ones_used_and_the_cutoff_too() {
    let genome = random(2000, 3);
    let reads = fastq_reads(&short_reads(&genome, false, 100, 40.0, 0.0, 4));
    let (_, note) = assemble(&reads, "21,31", "off", 10).unwrap();
    assert!(note.ends_with("k = 31"), "{note}");
    assert!(assemble(&reads, "21", "2.5", 10).is_ok());
    assert!(assemble(&reads, "21", "lots", 10)
        .unwrap_err()
        .contains("not a number"));
}

#[test]
fn reads_with_nothing_good_enough_in_them_fail_with_the_reason() {
    let mut r = fastq_reads(&[random(100, 5)]);
    r[0].qual = vec![2; 100];
    let e = assemble(&r, "auto", "auto", 10).unwrap_err();
    assert!(
        e.contains("no run of 21 bases or more is left") && e.contains("under quality 10"),
        "{e}"
    );
    let e = assemble(&fastq_reads(&[random(100, 6)]), "31", "auto", 10);
    assert!(e.is_ok());
    let e = assemble(&fastq_reads(&[random(40, 7)]), "45", "auto", 10).unwrap_err();
    assert!(e.contains("no run of 45 bases"), "{e}");
}

#[test]
fn a_graph_and_the_long_reads_become_the_plasmid_as_one_circle() {
    let d = folder("hybrid");
    let r = ring(3000, 1200, 1);
    let g = short_graph(&r.genome, true, 2);
    let text = gfa(&g, &number(&g));
    let reads = long_reads(&r.genome, true, 5000, 40.0, 0.02, 3);
    let file = d.join("long.fasta");
    let fasta: String = reads
        .iter()
        .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
        .collect();
    std::fs::write(&file, fasta).unwrap();
    let (value, note) = hybrid(&text, Some(&file), "normal", "auto", "auto", 0).unwrap();
    let Value::Sequences(rows) = value else {
        panic!()
    };
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].0.starts_with("1 length=11400 depth=") && rows[0].0.ends_with("x circular=true"),
        "{}",
        rows[0].0
    );
    assert!(is_circle_of(&rows[0].1, &r.genome));
    assert!(
        note.contains("2 bridges made") && note.contains("1 replicons, 1 circular"),
        "{note}"
    );
    // the settings are read: a quality that nothing reaches makes no bridge
    let (none, note) = hybrid(&text, Some(&file), "normal", "100", "auto", 0).unwrap();
    let Value::Sequences(rows) = none else {
        panic!()
    };
    assert!(rows.len() > 1 && note.contains("0 bridges made"), "{note}");
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn the_mode_and_the_depth_of_a_single_copy_are_passed_on() {
    let d = folder("hybrid-mode");
    let r = ring(3000, 1200, 1);
    let g = short_graph(&r.genome, true, 2);
    let text = gfa(&g, &number(&g));
    let reads = long_reads(&r.genome, true, 5000, 40.0, 0.02, 3);
    let file = d.join("long.fasta");
    let fasta: String = reads
        .iter()
        .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
        .collect();
    std::fs::write(&file, fasta).unwrap();
    // told that a single copy is half of what it is, the segments look like repeats and nothing is bridged
    let (_, note) = hybrid(&text, Some(&file), "bold", "auto", "9.0", 0).unwrap();
    assert!(note.contains("0 bridges made"), "{note}");
    let (_, note) = hybrid(&text, Some(&file), "conservative", "auto", "auto", 0).unwrap();
    assert!(note.contains("2 bridges made"), "{note}");
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn a_hybrid_with_no_long_reads_or_a_graph_that_is_not_one_says_so() {
    assert!(
        hybrid("S\tA\tACGTACGT\n", None, "normal", "auto", "auto", 0)
            .unwrap_err()
            .contains("no long reads")
    );
    let d = folder("hybrid-bad");
    let f = d.join("l.fasta");
    std::fs::write(&f, ">r\nACGT\n").unwrap();
    assert!(hybrid("not a graph", Some(&f), "normal", "auto", "auto", 0).is_err());
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn the_cutoff_that_is_asked_for_is_the_one_the_assembler_uses() {
    // a line at 40 times and a stretch seen about twice beside it: with the cutoff off the stretch stays as a segment, with auto it goes
    let genome = random(3000, 31);
    let mut pieces = short_reads(&genome, false, 100, 40.0, 0.0, 32);
    let island = random(300, 33);
    pieces.extend(short_reads(&island, false, 100, 2.0, 0.0, 34));
    let reads = fastq_reads(&pieces);
    let segments = |cutoff: &str| {
        let Value::Gfa(text) = assemble(&reads, "21", cutoff, 10).unwrap().0 else {
            panic!()
        };
        text.lines().filter(|l| l.starts_with("S\t")).count()
    };
    assert!(segments("off") > segments("auto"));
}

#[test]
fn conservative_is_the_mode_of_that_name() {
    // four reads cross the gap between two segments of 700 bases: a bridge of the reads' own bases, rated 19.6 (4/5 of the reads, 0.49 for the
    // lengths, half for the bases), which normal (10) makes and conservative (25) does not
    let d = folder("hybrid-conservative");
    let (x, z, gap) = (random(700, 41), random(700, 42), random(300, 43));
    let mut text = String::new();
    for (n, s) in [("A", &x), ("B", &z)] {
        text += &format!("S\t{n}\t{}\tDP:f:20\n", String::from_utf8_lossy(s));
    }
    let read = [&x[100..], &gap[..], &z[..600]].concat();
    let reads: String = (0..4)
        .map(|i| format!(">r{i}\n{}\n", String::from_utf8_lossy(&read)))
        .collect();
    let file = d.join("long.fasta");
    std::fs::write(&file, reads).unwrap();
    let normal = hybrid(&text, Some(&file), "normal", "auto", "auto", 0)
        .unwrap()
        .1;
    let conservative = hybrid(&text, Some(&file), "conservative", "auto", "auto", 0)
        .unwrap()
        .1;
    assert!(normal.contains("1 bridges made"), "{normal}");
    assert!(conservative.contains("0 bridges made"), "{conservative}");
    std::fs::remove_dir_all(&d).ok();
}
