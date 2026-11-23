use super::*;

fn fastas(locus: &str) -> Option<String> {
    match locus {
        "abcZ" => Some(
            ">abcZ_1\nACGTACGTAC\n>abcZ_2 a comment\nACGTACGTAA\n>abcZ-3\nacgtacgtag\n".to_string(),
        ),
        "adk" => Some(">adk_1\nGGGGCCCC\n>adk_7\nGGGGCCCA\n".to_string()),
        _ => None,
    }
}

const TABLE: &str = "ST\tabcZ\tadk\tclonal_complex\n1\t1\t7\tCC1\n2\t3\t1\t\n";

#[test]
fn a_scheme_is_read_with_its_loci_alleles_and_profiles() {
    let s = Scheme::from_texts("neisseria", TABLE, &fastas).unwrap();
    assert_eq!(s.name, "neisseria");
    assert_eq!(
        s.loci.iter().map(|l| l.name.as_str()).collect::<Vec<_>>(),
        ["abcZ", "adk"]
    );
    let ids: Vec<&str> = s.loci[0].alleles.iter().map(|a| a.id.as_str()).collect();
    assert_eq!(ids, ["1", "2", "3"]);
    // a lower-case sequence is upper case, and a header may be followed by a comment
    assert_eq!(s.loci[0].alleles[2].seq, b"ACGTACGTAG");
    assert_eq!(
        s.profiles,
        [
            ("1".to_string(), vec!["1".to_string(), "7".to_string()]),
            ("2".to_string(), vec!["3".to_string(), "1".to_string()])
        ]
    );
}

#[test]
fn a_sequence_type_is_found_by_its_alleles_and_only_by_all_of_them() {
    let s = Scheme::from_texts("n", TABLE, &fastas).unwrap();
    assert_eq!(s.sequence_type(&["1", "7"]), Some("1"));
    assert_eq!(s.sequence_type(&["3", "1"]), Some("2"));
    assert_eq!(s.sequence_type(&["1", "1"]), None);
    assert_eq!(s.sequence_type(&["1"]), None);
}

#[test]
fn ids_that_are_letters_and_the_other_columns_are_taken_as_they_come() {
    let table = "ST\tadk\tspecies\tCC\tLineage\n10\tA1\tx\ty\tz\n";
    let f = |l: &str| (l == "adk").then(|| ">adk_A1\nACGT\n>adk-B2\nACGA\n".to_string());
    let s = Scheme::from_texts("s", table, &f).unwrap();
    assert_eq!(s.loci.len(), 1);
    assert_eq!(
        s.loci[0]
            .alleles
            .iter()
            .map(|a| a.id.as_str())
            .collect::<Vec<_>>(),
        ["A1", "B2"]
    );
    assert_eq!(s.sequence_type(&["A1"]), Some("10"));
    // the header's ST may be in any case, and blank lines and carriage returns do no harm
    let crlf = "st\tadk\r\n\r\n10\tA1\r\n";
    assert_eq!(Scheme::from_texts("s", crlf, &f).unwrap().profiles.len(), 1);
}

#[test]
fn a_line_of_blanks_is_skipped_and_a_row_without_a_type_is_refused() {
    let table = "ST\tabcZ\tadk\n1\t1\t7\n  \t \n\t\n2\t3\t1\n";
    assert_eq!(
        Scheme::from_texts("s", table, &fastas)
            .unwrap()
            .profiles
            .len(),
        2
    );
    let why = Scheme::from_texts("s", "ST\tabcZ\tadk\n1\t1\t7\n\t3\t1\n", &fastas).unwrap_err();
    assert!(
        why.contains("Row 3") && why.contains("no sequence type"),
        "{why}"
    );
}

#[test]
fn what_is_wrong_with_a_scheme_is_said_with_the_place() {
    let why = |table: &str, f: &dyn Fn(&str) -> Option<String>| {
        Scheme::from_texts("s", table, f).unwrap_err()
    };
    let all = fastas;
    assert!(why("", &all).contains("empty"));
    assert!(why("locus\tabcZ\n", &all).contains("should be ST"));
    assert!(why("ST\tnolocus\n1\t1\n", &all).contains("no file of alleles for it (nolocus.tfa)"));
    assert!(why("ST\tclonal_complex\n1\tx\n", &all).contains("no locus column"));
    assert!(why("ST\tabcZ\tadk\n1\t1\n", &all)
        .contains("Row 2 of the profile table has too few columns"));
    let bad_name = |l: &str| (l == "adk").then(|| ">other_1\nACGT\n".to_string());
    assert!(why("ST\tadk\n1\t1\n", &bad_name)
        .contains("Record 1 of adk's alleles is named \"other_1\""));
    let no_id = |l: &str| (l == "adk").then(|| ">adk_\nACGT\n".to_string());
    assert!(why("ST\tadk\n1\t1\n", &no_id).contains("should be adk_ID"));
    let none = |l: &str| (l == "adk").then(String::new);
    assert!(why("ST\tadk\n1\t1\n", &none).contains("no allele of adk"));
    let empty_seq = |l: &str| (l == "adk").then(|| ">adk_1\n>adk_2\nACGT\n".to_string());
    assert!(why("ST\tadk\n1\t1\n", &empty_seq).contains("adk_1 has no sequence"));
}

#[test]
fn a_folder_is_a_scheme_named_for_it_and_its_files_are_found_by_name() {
    let dir = std::env::temp_dir().join(format!("helixloom-scheme-{}", std::process::id()));
    let folder = dir.join("neis");
    std::fs::create_dir_all(&folder).unwrap();
    std::fs::write(folder.join("neis.txt"), TABLE).unwrap();
    std::fs::write(folder.join("abcZ.tfa"), ">abcZ_1\nACGT\n>abcZ_3\nACGA\n").unwrap();
    std::fs::write(folder.join("adk.fasta"), ">adk_1\nGGGG\n>adk_7\nGGGA\n").unwrap();
    let s = Scheme::load(&folder).unwrap();
    assert_eq!(
        (s.name.as_str(), s.loci.len(), s.profiles.len()),
        ("neis", 2, 2)
    );
    // a table that is not named for the folder is taken when it is the only one
    let other = dir.join("other");
    std::fs::create_dir_all(&other).unwrap();
    std::fs::write(other.join("profiles.txt"), TABLE).unwrap();
    std::fs::write(other.join("abcZ.fa"), ">abcZ_1\nACGT\n").unwrap();
    std::fs::write(other.join("adk.fna"), ">adk_1\nGGGG\n").unwrap();
    assert_eq!(Scheme::load(&other).unwrap().name, "other");
    // none, or several, say so
    let empty = dir.join("empty");
    std::fs::create_dir_all(&empty).unwrap();
    assert!(Scheme::load(&empty)
        .unwrap_err()
        .contains("no profile table: expected empty.txt"));
    std::fs::write(empty.join("a.txt"), "x").unwrap();
    std::fs::write(empty.join("b.txt"), "y").unwrap();
    assert!(Scheme::load(&empty)
        .unwrap_err()
        .contains("several .txt files"));
    assert!(Scheme::load(&dir.join("missing")).is_err());
    std::fs::remove_dir_all(&dir).ok();
}
