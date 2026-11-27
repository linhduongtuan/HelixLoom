use super::*;

fn folder(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("helixloom-reads-{name}-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn fasta_reads_are_named_by_the_first_word_in_capitals_and_have_no_qualities() {
    let reads = from_text(">r1 some words\nacgt\nAC\n>r2\nTTTT\n").unwrap();
    assert_eq!(reads.len(), 2);
    assert_eq!(reads[0].name, "r1");
    assert_eq!(reads[0].seq, b"ACGTAC");
    assert!(reads[0].qual.is_empty());
    assert_eq!(reads[1].name, "r2");
}

#[test]
fn fastq_reads_keep_their_qualities_and_are_named_by_the_first_word() {
    let reads = from_text("@a b c\nacgt\n+\n!5I?\n@d\nGG\n+\nII\n").unwrap();
    assert_eq!(reads[0].name, "a");
    assert_eq!(reads[0].seq, b"ACGT");
    assert_eq!(reads[0].qual, [0, 20, 40, 30]);
    assert_eq!(reads[1].qual, [40, 40]);
}

#[test]
fn a_text_that_is_neither_gives_an_error_and_a_broken_fastq_says_which_record() {
    assert!(from_text("hello\n")
        .unwrap_err()
        .contains("neither FASTA nor FASTQ"));
    assert!(from_text("").unwrap_err().contains("neither"));
    let e = from_text("@a\nACGT\n+\n!!\n").unwrap_err();
    assert!(e.contains("Record 1"), "{e}");
}

#[test]
fn files_are_joined_with_commas_and_a_missing_or_empty_one_is_an_error() {
    let d = folder("load");
    std::fs::write(d.join("a.fa"), ">x\nACGT\n").unwrap();
    std::fs::write(d.join("b.fq"), "@y\nTT\n+\nII\n").unwrap();
    std::fs::write(d.join("empty.fa"), "").unwrap();
    let both = format!("{},{}", d.join("a.fa").display(), d.join("b.fq").display());
    let reads = load(&both).unwrap();
    assert_eq!(
        reads.iter().map(|r| r.name.as_str()).collect::<Vec<_>>(),
        ["x", "y"]
    );
    let missing = load(&d.join("none.fa").display().to_string()).unwrap_err();
    assert!(missing.contains("Could not read"), "{missing}");
    let empty = load(&d.join("empty.fa").display().to_string()).unwrap_err();
    assert!(empty.contains("neither"), "{empty}");
    assert!(load("").unwrap_err().contains("holds no sequences"));
    // a header line with nothing is named unnamed
    assert_eq!(from_text(">\nAC\n").unwrap()[0].name, "unnamed");
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn u_is_read_as_t_from_fasta_and_fastq() {
    let fasta = from_text(">r\nacguACGU\n").unwrap();
    assert_eq!(fasta[0].seq, b"ACGTACGT");
    let fastq = from_text("@r\nUuAg\n+\nIIII\n").unwrap();
    assert_eq!(fastq[0].seq, b"TTAG");
}
