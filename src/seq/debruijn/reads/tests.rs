use super::*;

fn dir(name: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("helixloom-dbg-reads-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

fn read(seq: &str, qual: &[u8]) -> Read {
    Read {
        seq: seq.as_bytes().to_vec(),
        qual: qual.to_vec(),
    }
}

#[test]
fn a_read_is_cut_at_a_base_that_is_not_known_and_the_runs_that_are_long_enough_are_kept() {
    let r = read("ACGTACGNACGTAC", &[]);
    assert_eq!(
        pieces(std::slice::from_ref(&r), 0, 5),
        [b"ACGTACG".to_vec(), b"ACGTAC".to_vec()]
    );
    assert_eq!(
        pieces(std::slice::from_ref(&r), 0, 7),
        [b"ACGTACG".to_vec()]
    );
    assert_eq!(pieces(&[r], 0, 8), Vec::<Vec<u8>>::new());
}

#[test]
fn a_read_is_cut_at_a_base_under_the_quality_asked_for_and_not_at_one_of_exactly_it() {
    // qualities 40 40 40 9 10 40 40 40: the fourth is under 10
    let q = [40, 40, 40, 9, 10, 40, 40, 40];
    let r = read("ACGTACGT", &q);
    assert_eq!(
        pieces(std::slice::from_ref(&r), 10, 3),
        [b"ACG".to_vec(), b"ACGT".to_vec()]
    );
    assert_eq!(
        pieces(std::slice::from_ref(&r), 11, 3),
        [b"ACG".to_vec(), b"CGT".to_vec()]
    );
    assert_eq!(pieces(&[r], 0, 3), [b"ACGTACGT".to_vec()]);
}

#[test]
fn a_run_at_the_end_is_kept_lower_case_is_made_upper_and_u_is_a_base() {
    let r = read("acgtNacguacgt", &[]);
    assert_eq!(pieces(&[r], 0, 4), [b"ACGT".to_vec(), b"ACGUACGT".to_vec()]);
}

#[test]
fn reads_with_no_quality_are_all_good() {
    let r = read("ACGTACGT", &[]);
    assert_eq!(pieces(&[r], 40, 8), [b"ACGTACGT".to_vec()]);
}

#[test]
fn fasta_and_fastq_files_are_read_and_several_are_joined_with_commas() {
    let d = dir("load");
    let fa = d.join("a.fa");
    let fq = d.join("b.fq");
    std::fs::write(&fa, ">r1 first\nACGT\nACGT\n>r2\nTTTT\n").unwrap();
    std::fs::write(&fq, "@q1\nACGTA\n+\nIIII#\n").unwrap();
    let a = load(fa.to_str().unwrap()).unwrap();
    assert_eq!(a.len(), 2);
    assert_eq!(
        (a[0].seq.as_slice(), a[0].qual.len()),
        (&b"ACGTACGT"[..], 0)
    );
    let b = load(fq.to_str().unwrap()).unwrap();
    assert_eq!(b[0].seq, b"ACGTA");
    assert_eq!(b[0].qual, [40, 40, 40, 40, 2]);
    let both = load(&format!("{},{}", fa.display(), fq.display())).unwrap();
    assert_eq!(both.len(), 3);
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn what_cannot_be_read_says_why() {
    let d = dir("errors");
    let not = d.join("x.txt");
    std::fs::write(&not, "hello\n").unwrap();
    assert!(load(not.to_str().unwrap())
        .unwrap_err()
        .contains("neither FASTA nor FASTQ"));
    let empty = d.join("e.fa");
    std::fs::write(&empty, "").unwrap();
    assert!(load(empty.to_str().unwrap())
        .unwrap_err()
        .contains("neither FASTA nor FASTQ"));
    assert!(load(d.join("nothing").to_str().unwrap())
        .unwrap_err()
        .contains("Could not read"));
    assert!(load("").unwrap_err().contains("holds no reads"));
    let cut = d.join("c.fq");
    std::fs::write(&cut, "@r\nACGT\n").unwrap();
    assert!(load(cut.to_str().unwrap())
        .unwrap_err()
        .contains("cut short"));
    std::fs::remove_dir_all(&d).ok();
}
