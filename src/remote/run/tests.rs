use super::*;

fn isolates() -> Vec<steps::Named> {
    crate::sample::tree_demo()
        .into_iter()
        .map(|d| (d.name, d.seq))
        .collect()
}

#[test]
fn an_alignment_is_made_as_the_app_makes_it_itself() {
    let seqs = isolates();
    let Response::Align(here) = execute(
        Request::Align {
            seqs: seqs.clone(),
            mafft: false,
        },
        &Env::default(),
    )
    .expect("aligned") else {
        panic!("an alignment");
    };
    let dir = tools::work_folder("test-align").expect("a folder");
    let app = steps::align(&seqs, None, &dir).expect("aligned");
    assert_eq!((here.rows, here.engine), (app.rows, app.engine));
}

#[test]
fn one_tree_runs_on_the_cpu_and_the_estimates_have_their_names() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let tree = crate::seq::phylocnn::zurich_tree().expect("the Zurich tree");
    let r = run_phylocnn(tree.trim(), 0.25, Task::Bdss, Priors::Hiv, true).expect("estimates");
    assert_eq!(r.engine, "CPU", "a single tree never goes to the GPU");
    assert_eq!((r.trees, r.tips), (1, 200));
    let names: Vec<&str> = r.per_tree[0].iter().map(|e| e.name).collect();
    assert_eq!(names[0], "R0");
    assert!(run_phylocnn("not a tree", 0.25, Task::Bd, Priors::Simulated, false).is_err());
    assert!(run_phylocnn("(A:1,B:1);", 2.0, Task::Bd, Priors::Simulated, false).is_err());
}

#[test]
fn a_fold_with_no_program_configured_says_so() {
    let why = execute(
        Request::Fold {
            sequence: b"MQIFVKTLTGKTITLEVEP".to_vec(),
            local: true,
        },
        &Env::default(),
    )
    .err()
    .expect("an error");
    assert!(why.contains("no folding program"), "{why}");
}

#[cfg(unix)]
mod fold_program {
    use super::*;
    use std::path::{Path, PathBuf};

    const STRUCTURE: &str = include_str!("../../../assets/structure/1ubq.pdb");

    /// A shell script that is run as `/bin/sh script input output`.
    fn program(dir: &Path, body: &str) -> Vec<String> {
        let script: PathBuf = dir.join("fold.sh");
        std::fs::write(&script, format!("{body}\n")).expect("a script");
        vec!["/bin/sh".into(), script.display().to_string()]
    }

    fn fold(command: Vec<String>, seconds: u64) -> Result<Response, String> {
        execute(
            Request::Fold {
                sequence: b"MQIFVKTLTGKTITLEVEP".to_vec(),
                local: true,
            },
            &Env {
                fold_command: Some(command),
                fold_seconds: seconds,
            },
        )
    }

    #[test]
    fn the_program_gets_the_files_and_its_structure_comes_back() {
        let dir = tools::work_folder("test-fold-ok").expect("a folder");
        let canned = dir.join("canned.pdb");
        std::fs::write(&canned, STRUCTURE).expect("a structure");
        // the sequence arrives as FASTA in $1, the answer goes to $2
        let command = program(
            &dir,
            &format!(
                "grep -q '>query' \"$1\" && grep -q MQIFVK \"$1\" && cp '{}' \"$2\"",
                canned.display()
            ),
        );
        let Response::Fold(done) = fold(command, 30).expect("a structure") else {
            panic!("a structure");
        };
        assert_eq!(done.pdb, STRUCTURE);
        assert!(done.engine.contains("folding program"));
    }

    #[test]
    fn a_program_that_fails_or_writes_nothing_useful_is_reported() {
        let dir = tools::work_folder("test-fold-bad").expect("a folder");
        let why = |body: &str| fold(program(&dir, body), 30).err().expect("an error");
        assert!(why("exit 3").contains("stopped with"), "exit code");
        assert!(why("true").contains("wrote no structure"), "no file");
        assert!(
            why("echo nothing > \"$2\"").contains("no atoms"),
            "no atoms"
        );
        let missing = fold(vec!["/definitely/not/a/program".into()], 30)
            .err()
            .expect("an error");
        assert!(missing.contains("Could not start"), "{missing}");
        assert!(fold(Vec::new(), 30).is_err());
    }

    #[test]
    fn a_program_that_runs_too_long_is_stopped() {
        let dir = tools::work_folder("test-fold-slow").expect("a folder");
        let started = Instant::now();
        let why = fold(program(&dir, "sleep 30"), 1).err().expect("an error");
        assert!(why.contains("did not finish"), "{why}");
        assert!(
            started.elapsed() < Duration::from_secs(20),
            "it was stopped, not waited for"
        );
    }
}
