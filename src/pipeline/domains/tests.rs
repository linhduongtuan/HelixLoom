use super::testing::*;
use super::*;
use std::path::PathBuf;

fn queries() -> Vec<(String, Vec<u8>)> {
    vec![
        ("first protein".to_string(), vec![b'M'; 60]),
        ("second".to_string(), vec![b'G'; 30]),
    ]
}

fn lines(path: PathBuf) -> Vec<String> {
    std::fs::read_to_string(path)
        .expect("a file")
        .lines()
        .map(str::to_string)
        .collect()
}

#[cfg(unix)]
#[test]
fn the_program_is_asked_for_the_domain_table_with_the_curated_cutoffs_and_the_answer_is_matched_to_the_sequences(
) {
    let folder = Folder::new("domains-run");
    let db = pressed_db(&folder);
    let program = stand_in(&folder, TABLE, "ok");
    let hits = scan(&program, &db, &queries(), Cutoff::Curated).expect("domains");
    assert_eq!(hits.len(), 2);
    assert_eq!((hits[0].sequence, hits[0].family.as_str()), (0, "Pkinase"));
    assert_eq!((hits[1].sequence, hits[1].family.as_str()), (1, "Globin"));
    // what the program was given: no alignments, the main report and the table in the work folder, the curated cutoff, the database,
    // and last the sequences
    let args = lines(folder.0.join("args.txt"));
    assert_eq!(args[0], "--noali");
    assert_eq!(args[1], "-o");
    assert_eq!(args[3], "--domtblout");
    assert_eq!(
        &args[5..],
        ["--cut_ga", &db.display().to_string(), &args[args.len() - 1]]
    );
    assert!(
        args[2].ends_with("report.txt") && args[4].ends_with("domains.tbl"),
        "{args:?}"
    );
    // the sequences went under names of the program's own, whatever the names were
    assert_eq!(
        std::fs::read_to_string(folder.0.join("seen.fa")).expect("the sequences"),
        format!(">q1\n{}\n>q2\n{}\n", "M".repeat(60), "G".repeat(30))
    );
    // and the work folder is gone afterwards
    let work = Path::new(&args[2])
        .parent()
        .expect("a folder")
        .to_path_buf();
    assert!(!work.exists(), "{} was left", work.display());
}

#[cfg(unix)]
#[test]
fn an_e_value_is_given_for_the_sequence_and_the_domain_alike() {
    let folder = Folder::new("domains-evalue");
    let db = pressed_db(&folder);
    let program = stand_in(&folder, TABLE, "ok");
    scan(&program, &db, &queries(), Cutoff::Evalue(0.001)).expect("domains");
    let args = lines(folder.0.join("args.txt"));
    assert_eq!(&args[5..9], ["-E", "0.001", "--domE", "0.001"]);
    assert!(!args.iter().any(|a| a == "--cut_ga"));
}

#[cfg(unix)]
#[test]
fn a_program_that_fails_or_writes_nothing_or_writes_what_does_not_fit_is_an_error_with_the_reason()
{
    let folder = Folder::new("domains-fail");
    let db = pressed_db(&folder);
    let why = scan(
        &stand_in(&folder, TABLE, "fail"),
        &db,
        &queries(),
        Cutoff::Curated,
    )
    .unwrap_err();
    assert!(
        why.contains("hmmscan stopped") && why.contains("database is damaged"),
        "{why}"
    );
    let why = scan(
        &stand_in(&folder, TABLE, "silent"),
        &db,
        &queries(),
        Cutoff::Curated,
    )
    .unwrap_err();
    assert!(why.contains("wrote no table"), "{why}");
    // a table about a sequence that was not sent, or of another length
    let other = TABLE.replace(" q2 ", " q7 ");
    let why = scan(
        &stand_in(&folder, &other, "ok"),
        &db,
        &queries(),
        Cutoff::Curated,
    )
    .unwrap_err();
    assert!(why.contains("not searched"), "{why}");
    let longer = TABLE.replace(" - 60 ", " - 61 ");
    let why = scan(
        &stand_in(&folder, &longer, "ok"),
        &db,
        &queries(),
        Cutoff::Curated,
    )
    .unwrap_err();
    assert!(why.contains("not 60"), "{why}");
    // no domain found is not an error
    let none = scan(
        &stand_in(&folder, "# nothing found\n", "ok"),
        &db,
        &queries(),
        Cutoff::Curated,
    );
    assert_eq!(none, Ok(Vec::new()));
}

#[cfg(unix)]
#[test]
fn what_cannot_be_searched_is_refused_before_the_program_is_started() {
    let folder = Folder::new("domains-refuse");
    let db = pressed_db(&folder);
    let program = stand_in(&folder, TABLE, "ok");
    let refused = |q: Vec<(String, Vec<u8>)>, cutoff: Cutoff, why: &str| {
        let said = scan(&program, &db, &q, cutoff).unwrap_err();
        assert!(said.contains(why), "{why}: {said}");
    };
    refused(Vec::new(), Cutoff::Curated, "no sequence");
    refused(
        vec![("e".into(), Vec::new())],
        Cutoff::Curated,
        "cannot be searched",
    );
    refused(
        vec![("d".into(), b"MK1".to_vec())],
        Cutoff::Curated,
        "letters only",
    );
    refused(queries(), Cutoff::Evalue(0.0), "above 0");
    refused(queries(), Cutoff::Evalue(1000.0), "at most 100");
    refused(queries(), Cutoff::Evalue(f64::NAN), "above 0");
    assert!(
        !folder.0.join("args.txt").exists(),
        "the program was started"
    );
    // an unpressed database says what to do
    let bare = folder.0.join("Other.hmm");
    std::fs::write(&bare, "HMMER3/f\n").expect("a file");
    let said = scan(&program, &bare, &queries(), Cutoff::Curated).unwrap_err();
    assert!(said.contains("hmmpress") && said.contains(".h3m"), "{said}");
    assert!(!folder.0.join("args.txt").exists());
}

#[test]
fn a_cutoff_is_checked_and_a_pressed_database_has_its_four_files() {
    assert!(Cutoff::Curated.check().is_ok());
    assert!(Cutoff::Evalue(1e-5).check().is_ok());
    assert!(Cutoff::Evalue(100.0).check().is_ok());
    assert!(Cutoff::Evalue(-1.0).check().is_err());
    let folder = Folder::new("domains-pressed");
    let db = pressed_db(&folder);
    assert!(check_pressed(&db).is_ok());
    std::fs::remove_file(format!("{}.h3p", db.display())).expect("removed");
    let said = check_pressed(&db).unwrap_err();
    assert!(said.contains(".h3p") && !said.contains(".h3m"), "{said}");
}
