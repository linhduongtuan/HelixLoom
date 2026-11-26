use super::*;
use crate::pipeline::asm::tests::{fastq_reads, folder};
use crate::seq::debruijn::sim::random;

#[cfg(unix)]
fn script(dir: &Path, name: &str, body: &str) -> std::path::PathBuf {
    use std::os::unix::fs::PermissionsExt;
    let path = dir.join(name);
    std::fs::write(&path, format!("#!/bin/sh\n{body}\n")).unwrap();
    std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
    path
}

fn flye_stand_in(dir: &Path) -> std::path::PathBuf {
    script(
        dir,
        "flye",
        "echo \"$@\" > args.txt\nmkdir -p flye-out\nprintf '>contig_1 extra\\nACGTACGT\\n>contig_2\\nGGGGCCCC\\n' > flye-out/assembly.fasta\nprintf '#seq_name\\tlength\\tcov.\\tcirc.\\trepeat\\tmult.\\talt_group\\tgraph_path\\ncontig_1\\t8\\t30\\tY\\tN\\t1\\t*\\t1\\ncontig_2\\t8\\t30\\tN\\tN\\t1\\t*\\t2\\n' > flye-out/assembly_info.txt",
    )
}

#[cfg(unix)]
#[test]
fn flye_is_called_with_the_reads_the_mode_and_the_size_and_its_circles_are_read() {
    let d = folder("flye");
    let exe = flye_stand_in(&d);
    let work = d.join("work");
    std::fs::create_dir_all(&work).unwrap();
    let reads = fastq_reads(&[random(100, 1), random(100, 2)]);
    let rows = flye_with(&exe, "nano-hq", "5m", &reads, &work).unwrap();
    let args = std::fs::read_to_string(work.join("args.txt")).unwrap();
    let words: Vec<&str> = args.split_whitespace().collect();
    assert_eq!(
        &words[..5],
        [
            "--nano-hq",
            "reads.fastq",
            "--out-dir",
            "flye-out",
            "--threads"
        ]
    );
    assert_eq!(&words[6..], ["--genome-size", "5m"]);
    // the reads were written for it
    let written = std::fs::read_to_string(work.join("reads.fastq")).unwrap();
    assert_eq!(written, crate::seq::fastq::to_fastq(&reads));
    assert_eq!(rows[0].0, "contig_1 extra circular=true");
    assert_eq!(rows[0].1, b"ACGTACGT");
    assert_eq!(rows[1], ("contig_2".to_string(), b"GGGGCCCC".to_vec()));
    // with no size, none is given
    let work2 = d.join("work2");
    std::fs::create_dir_all(&work2).unwrap();
    flye_with(&exe, "nano-raw", "auto", &reads, &work2).unwrap();
    assert!(!std::fs::read_to_string(work2.join("args.txt"))
        .unwrap()
        .contains("genome-size"));
    std::fs::remove_dir_all(&d).ok();
}

#[cfg(unix)]
#[test]
fn a_flye_that_stops_or_writes_nothing_is_an_error_that_says_so() {
    let d = folder("flye-bad");
    let reads = fastq_reads(&[random(100, 3)]);
    let stops = script(&d, "stops", "echo 'out of memory' >&2\nexit 3");
    let e = flye_with(&stops, "nano-raw", "auto", &reads, &d).unwrap_err();
    assert!(e.contains("stopped") && e.contains("out of memory"), "{e}");
    let empty = script(
        &d,
        "empty",
        "mkdir -p flye-out\n: > flye-out/assembly.fasta",
    );
    let e = flye_with(&empty, "nano-raw", "auto", &reads, &d).unwrap_err();
    assert!(e.contains("Flye made no contig"), "{e}");
    // a folder that the stand-ins before have not written into
    let clean = d.join("clean");
    std::fs::create_dir_all(&clean).unwrap();
    let silent = script(&d, "silent", "true");
    assert!(flye_with(&silent, "nano-raw", "auto", &reads, &clean)
        .unwrap_err()
        .contains("wrote no assembly.fasta"));
    std::fs::remove_dir_all(&d).ok();
}

#[cfg(unix)]
#[test]
fn medaka_is_called_with_the_draft_and_the_reads_and_the_names_come_back() {
    let d = folder("medaka");
    let exe = script(
        &d,
        "medaka",
        "echo \"$@\" > args.txt\nmkdir -p medaka-out\ncp draft.fasta medaka-out/consensus.fasta",
    );
    let reads = d.join("reads.fastq");
    std::fs::write(&reads, "@r\nACGT\n+\nIIII\n").unwrap();
    let work = d.join("work");
    std::fs::create_dir_all(&work).unwrap();
    let draft = vec![
        ("contig_1 circular=true".to_string(), b"ACGTACGT".to_vec()),
        ("contig_2".to_string(), b"GGGGCCCC".to_vec()),
    ];
    let out = medaka_with(&exe, &draft, &reads, &work).unwrap();
    assert_eq!(out, draft);
    let args = std::fs::read_to_string(work.join("args.txt")).unwrap();
    let words: Vec<&str> = args.split_whitespace().collect();
    assert_eq!(words[0], "-i");
    assert!(
        words[1].ends_with("reads.fastq") && words[1].starts_with('/'),
        "{args}"
    );
    assert_eq!(
        &words[2..7],
        ["-d", "draft.fasta", "-o", "medaka-out", "-t"]
    );
    // it was given plain names
    let given = std::fs::read_to_string(work.join("draft.fasta")).unwrap();
    assert!(
        given.starts_with(">c0\n") && given.contains(">c1\n"),
        "{given}"
    );
    std::fs::remove_dir_all(&d).ok();
}

#[cfg(unix)]
#[test]
fn a_medaka_that_gives_back_other_names_in_the_same_order_or_other_contigs_is_dealt_with() {
    let d = folder("medaka-names");
    let reads = d.join("reads.fastq");
    std::fs::write(&reads, "@r\nACGT\n+\nIIII\n").unwrap();
    let draft = vec![
        ("a".to_string(), b"AAAA".to_vec()),
        ("b".to_string(), b"CCCC".to_vec()),
    ];
    // other names, same order, and the polished bases differ
    let renames = script(
        &d,
        "renames",
        "mkdir -p medaka-out\nprintf '>x\\nAAAT\\n>y\\nCCCG\\n' > medaka-out/consensus.fasta",
    );
    let out = medaka_with(&renames, &draft, &reads, &d).unwrap();
    assert_eq!(
        out,
        [
            ("a".to_string(), b"AAAT".to_vec()),
            ("b".to_string(), b"CCCG".to_vec())
        ]
    );
    // our names, in the other order
    let swaps = script(
        &d,
        "swaps",
        "mkdir -p medaka-out\nprintf '>c1\\nCCCG\\n>c0\\nAAAT\\n' > medaka-out/consensus.fasta",
    );
    let out = medaka_with(&swaps, &draft, &reads, &d).unwrap();
    assert_eq!(
        out,
        [
            ("a".to_string(), b"AAAT".to_vec()),
            ("b".to_string(), b"CCCG".to_vec())
        ]
    );
    // fewer contigs than it was given
    let fewer = script(
        &d,
        "fewer",
        "mkdir -p medaka-out\nprintf '>c0\\nAAAT\\n' > medaka-out/consensus.fasta",
    );
    let e = medaka_with(&fewer, &draft, &reads, &d).unwrap_err();
    assert!(e.contains("1 contigs for the 2"), "{e}");
    let clean = d.join("clean");
    std::fs::create_dir_all(&clean).unwrap();
    let none = script(&d, "none", "true");
    assert!(medaka_with(&none, &draft, &reads, &clean)
        .unwrap_err()
        .contains("wrote no consensus.fasta"));
    assert!(medaka_with(&none, &draft, &d.join("missing.fastq"), &d)
        .unwrap_err()
        .contains("Could not find the reads"));
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn a_step_that_calls_no_program_has_no_version_and_is_always_ready() {
    let t = AsmStep::Trim {
        adapter: "none".into(),
        quality: 20,
        min_length: 20,
        cut: crate::pipeline::asm::Cut::Off,
    };
    assert_eq!(version(&t), "");
    assert!(ready(&t).is_ok());
}

#[test]
fn flye_and_medaka_are_found_by_name_and_say_so_when_they_are_not_there() {
    for (step, label) in [
        (
            AsmStep::Flye {
                mode: "nano-raw".into(),
                genome_size: "auto".into(),
            },
            "Flye",
        ),
        (AsmStep::Medaka { reads: "r".into() }, "Medaka"),
    ] {
        let there = tools::detect(program(&step).unwrap()).is_some();
        match ready(&step) {
            Ok(()) => assert!(there),
            Err(e) => {
                assert!(!there);
                assert!(
                    e.contains(&format!("{label} is not installed here"))
                        && e.contains("helixloom tools"),
                    "{e}"
                );
            }
        }
    }
}

#[test]
fn the_two_programs_are_in_the_list_of_programs_the_platform_looks_for() {
    assert!(tools::ALL
        .iter()
        .any(|p| p.label == "Flye" && p.names == ["flye"]));
    assert!(tools::ALL
        .iter()
        .any(|p| p.label == "Medaka" && p.names == ["medaka"]));
}

#[cfg(unix)]
#[test]
fn a_medaka_that_gives_more_contigs_than_it_was_given_is_also_an_error() {
    let d = folder("medaka-more");
    let reads = d.join("reads.fastq");
    std::fs::write(&reads, "@r\nACGT\n+\nIIII\n").unwrap();
    let draft = vec![("a".to_string(), b"AAAA".to_vec())];
    let more = script(
        &d,
        "more",
        "mkdir -p medaka-out\nprintf '>c0\\nAAAT\\n>c1\\nCCCG\\n' > medaka-out/consensus.fasta",
    );
    let e = medaka_with(&more, &draft, &reads, &d).unwrap_err();
    assert!(e.contains("2 contigs for the 1"), "{e}");
    std::fs::remove_dir_all(&d).ok();
}
