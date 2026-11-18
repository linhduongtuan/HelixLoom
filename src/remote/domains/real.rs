//! The domain search with a real HMMER and with the real Pfam. These are left out of the usual run, because neither HMMER nor Pfam
//! is part of the repository; they are run by hand on a computer that has them (`docs/DOMAINS.md`):
//!
//! * `HELIXLOOM_HMMER_TUTORIAL=/opt/homebrew/share/doc/hmmer/tutorial cargo test --bin helixloom real_hmmer -- --ignored` needs
//!   HMMER 3.4 on the PATH and the tutorial folder that comes with it (`fn3.hmm`, `Pkinase.hmm`, `7LESS_DROME`);
//! * `cargo test --bin helixloom real_pfam -- --ignored` needs HMMER and Pfam-A, pressed, in the models folder.
//!
//! The numbers of the first are what `hmmscan` 3.4 itself writes for the same files when it is run directly
//! (`hmmscan --domtblout`), copied from its table; the families of the second are what is known of the proteins.

use super::*;
use std::path::{Path, PathBuf};
use std::process::Command;

/// The residues of a protein in UniProt's flat-file format: the lines between `SQ` and `//`.
fn flat_file_sequence(text: &str) -> Vec<u8> {
    text.lines()
        .skip_while(|l| !l.starts_with("SQ "))
        .skip(1)
        .take_while(|l| !l.starts_with("//"))
        .flat_map(str::bytes)
        .filter(u8::is_ascii_alphabetic)
        .collect()
}

/// The two Pfam families of HMMER's tutorial as one database, pressed, in a folder of its own.
fn tutorial_db(tutorial: &Path, hmmscan: &Path) -> PathBuf {
    let dir = tools::work_folder("real-hmmer").expect("a folder");
    let db = dir.join("minifam");
    let models: String = ["fn3.hmm", "Pkinase.hmm"]
        .iter()
        .map(|f| std::fs::read_to_string(tutorial.join(f)).expect("a model of the tutorial"))
        .collect();
    std::fs::write(&db, models).expect("the database");
    let pressed = Command::new(hmmscan.with_file_name("hmmpress"))
        .arg(&db)
        .output()
        .expect("hmmpress beside hmmscan");
    assert!(pressed.status.success(), "{pressed:?}");
    db
}

/// What is compared of a domain: the family, its place on the profile, on the sequence (the alignment and the envelope), its
/// score and the mean posterior probability.
type Seen = (&'static str, [usize; 6], f64, f64);

fn seen(hits: &[DomainHit]) -> Vec<(String, [usize; 6], f64, f64)> {
    hits.iter()
        .map(|h| {
            (
                h.family.clone(),
                [
                    h.hmm_from, h.hmm_to, h.ali_from, h.ali_to, h.env_from, h.env_to,
                ],
                h.score,
                h.accuracy,
            )
        })
        .collect()
}

fn assert_seen(hits: &[DomainHit], want: &[Seen]) {
    let got = seen(hits);
    assert_eq!(got.len(), want.len(), "{got:?}");
    for (g, w) in got.iter().zip(want) {
        assert_eq!(
            (g.0.as_str(), g.1, g.2, g.3),
            (w.0, w.1, w.2, w.3),
            "{got:?}"
        );
    }
}

#[test]
#[ignore = "needs HMMER 3.4 and its tutorial folder; see the top of the file"]
fn real_hmmer_gives_the_program_the_domains_it_writes_when_it_is_run_directly() {
    let tutorial = PathBuf::from(
        std::env::var("HELIXLOOM_HMMER_TUTORIAL")
            .expect("HELIXLOOM_HMMER_TUTORIAL names the folder"),
    );
    let program = tools::detect(&HMMSCAN).expect("hmmscan on the PATH");
    assert!(program.version.contains("3.4"), "{}", program.version);
    let db = tutorial_db(&tutorial, &program.path);
    let protein = flat_file_sequence(
        &std::fs::read_to_string(tutorial.join("7LESS_DROME"))
            .expect("the protein of the tutorial"),
    );
    assert_eq!(protein.len(), 2554);
    let job = |cutoff| DomainsJob {
        queries: vec![("sevenless".to_string(), protein.clone())],
        cutoff,
    };

    // the program's own thresholds are an E-value of 10 for the sequence and for a domain
    let all = run_with(&job(Cutoff::Evalue(10.0)), &program, &db).expect("a search");
    assert_seen(
        &all.hits,
        &[
            ("fn3", [61, 74, 396, 409, 395, 411], -1.3, 0.85),
            ("fn3", [2, 84, 439, 520, 437, 521], 40.7, 0.95),
            ("fn3", [13, 85, 836, 913, 826, 914], 14.4, 0.73),
            ("fn3", [10, 36, 1209, 1235, 1203, 1259], 5.1, 0.82),
            ("fn3", [14, 80, 1313, 1380, 1304, 1386], 24.3, 0.82),
            ("fn3", [58, 72, 1754, 1768, 1739, 1769], 0.0, 0.89),
            ("fn3", [1, 85, 1799, 1890, 1799, 1891], 47.2, 0.91),
            ("fn3", [6, 74, 1904, 1966, 1901, 1976], 17.8, 0.90),
            ("fn3", [1, 86, 1993, 2107, 1993, 2107], 12.8, 0.89),
            ("Pkinase", [2, 256, 2210, 2479, 2209, 2482], 136.5, 0.85),
        ],
    );
    let kinase = all.hits.last().expect("the kinase domain");
    assert_eq!(
        (
            kinase.accession.as_str(),
            kinase.family_length,
            kinase.description.as_str()
        ),
        ("PF00069.17", 260, "Protein kinase domain")
    );
    assert_eq!(
        (all.names.as_slice(), all.lengths.as_slice()),
        (&["sevenless".to_string()][..], &[2554][..])
    );
    assert_eq!(all.database, "minifam");
    assert!(all.program.contains("3.4"), "{}", all.program);

    // Pfam's own thresholds keep the six fibronectin domains that score enough, and the kinase
    let curated = run_with(&job(Cutoff::Curated), &program, &db).expect("a search");
    assert_seen(
        &curated.hits,
        &[
            ("fn3", [2, 84, 439, 520, 437, 521], 40.7, 0.95),
            ("fn3", [13, 85, 836, 913, 826, 914], 14.4, 0.73),
            ("fn3", [14, 80, 1313, 1380, 1304, 1386], 24.3, 0.82),
            ("fn3", [1, 85, 1799, 1890, 1799, 1891], 47.2, 0.91),
            ("fn3", [6, 74, 1904, 1966, 1901, 1976], 17.8, 0.90),
            ("fn3", [1, 86, 1993, 2107, 1993, 2107], 12.8, 0.89),
            ("Pkinase", [2, 256, 2210, 2479, 2209, 2482], 136.5, 0.85),
        ],
    );

    // and the answer survives the way to a workstation and back
    let back =
        DomainsResult::decode(&Container::decode(&all.encode().encode()).expect("a container"))
            .expect("a result");
    assert_eq!(back.hits, all.hits);
    let _ = std::fs::remove_dir_all(db.parent().expect("the folder of the database"));
}

/// The chains of the two structures that come with the program, as proteins to search: ubiquitin (PDB 1UBQ), and the light and
/// heavy variable domains of an antibody with the lysozyme it binds (PDB 1VFB).
fn chains() -> Vec<(String, Vec<u8>)> {
    let read = |text: &str| crate::seq::pdb::parse(text).expect("a structure");
    let ubiquitin = read(include_str!("../../../assets/structure/1ubq.pdb"));
    let complex = read(include_str!("../../../assets/structure/1vfb.pdb"));
    vec![
        ("ubiquitin".to_string(), ubiquitin.sequence('A')),
        ("antibody light chain".to_string(), complex.sequence('A')),
        ("antibody heavy chain".to_string(), complex.sequence('B')),
        ("lysozyme".to_string(), complex.sequence('C')),
    ]
}

#[test]
#[ignore = "needs HMMER and Pfam-A, pressed, in the models folder; see the top of the file"]
fn real_pfam_finds_the_families_that_are_known_for_the_proteins_that_come_with_the_program() {
    let queries = chains();
    assert_eq!(
        queries.iter().map(|(_, s)| s.len()).collect::<Vec<_>>(),
        [76, 107, 116, 129]
    );
    let found = run(&DomainsJob {
        queries,
        cutoff: Cutoff::Curated,
    })
    .expect("a search of Pfam");
    assert_eq!(found.database, "Pfam-A.hmm");
    assert!(found.program.contains("3.4"), "{}", found.program);
    // the family that scores best in each protein is the one it is known to belong to, and covers most of it
    let best = |sequence: usize| {
        found
            .hits
            .iter()
            .filter(|h| h.sequence == sequence)
            .max_by(|a, b| a.score.total_cmp(&b.score))
            .unwrap_or_else(|| panic!("no domain in {}", found.names[sequence]))
    };
    for (sequence, family, accession) in [
        (0, "ubiquitin", "PF00240."),
        (1, "V-set", "PF07686."),
        (2, "V-set", "PF07686."),
        (3, "Lys", "PF00062."),
    ] {
        let hit = best(sequence);
        let name = &found.names[sequence];
        assert_eq!(hit.family, family, "{name}");
        assert!(
            hit.accession.starts_with(accession),
            "{name}: {}",
            hit.accession
        );
        let covered = (hit.ali_to + 1 - hit.ali_from) as f64 / found.lengths[sequence] as f64;
        assert!(covered > 0.85, "{name}: {}-{}", hit.ali_from, hit.ali_to);
        assert!(hit.i_evalue < 1e-10, "{name}: {:e}", hit.i_evalue);
    }
    // and the example of the window is these four proteins: the bundled structures and the example cannot drift apart
    let example = hmmer::read_proteins(crate::ui::windows::domains::EXAMPLE).expect("the example");
    let residues =
        |list: &[(String, Vec<u8>)]| list.iter().map(|(_, s)| s.clone()).collect::<Vec<_>>();
    assert_eq!(residues(&example), residues(&chains()));
}
