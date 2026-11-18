//! What the tests of the domain search share: a folder, a pressed database and a program that stands in for `hmmscan`.

use super::*;
use std::path::PathBuf;

/// A folder for one test, removed when the test is over.
pub(crate) struct Folder(pub PathBuf);

impl Folder {
    pub(crate) fn new(label: &str) -> Folder {
        let dir = tools::work_folder(label).expect("a folder");
        Folder(dir)
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A database that has been "pressed": its text file and the four files beside it, empty (the stand-in program reads none).
pub(crate) fn pressed_db(folder: &Folder) -> PathBuf {
    let db = folder.0.join("Pfam-A.hmm");
    std::fs::write(&db, "HMMER3/f\n").expect("a file");
    for suffix in PRESSED {
        std::fs::write(format!("{}{suffix}", db.display()), "x").expect("a file");
    }
    db
}

/// A table in the format of the HMMER 3.4 User's Guide, made by hand, for the sequences `q1` (60 residues) and `q2` (30).
pub(crate) const TABLE: &str = "\
# target name accession tlen query name accession qlen E-value score bias # of c-Evalue i-Evalue score bias from to from to from to acc description of target
Pkinase PF00069.30 264 q1 - 60 1e-20 70.0 0.1 1 1 1e-17 2e-17 60.0 0.0 1 264 5 55 3 58 0.97 Protein kinase domain
Globin PF00042.25 109 q2 - 30 1e-5 20.0 0.0 1 1 1e-3 2e-3 15.0 0.0 5 25 2 22 1 24 0.90 Globin
";

/// A program that stands in for `hmmscan` (HMMER is another project's program and is not part of the tests): it writes the
/// arguments it was given, one to a line, beside itself, copies the sequence file it was given, then writes `table` where it
/// was asked to put the domain table, or does what `mode` says.
#[cfg(unix)]
pub(crate) fn stand_in(folder: &Folder, table: &str, mode: &str) -> Installed {
    use std::os::unix::fs::PermissionsExt;
    let script = folder.0.join("hmmscan");
    let fixture = folder.0.join("fixture.tbl");
    std::fs::write(&fixture, table).expect("a file");
    let args = folder.0.join("args.txt");
    let seen = folder.0.join("seen.fa");
    let body = format!(
        "#!/bin/sh\nprintf '%s\\n' \"$@\" > '{args}'\nout=''\nprev=''\nfor a in \"$@\"; do\n  [ \"$prev\" = '--domtblout' ] && out=\"$a\"\n  prev=\"$a\"\n  last=\"$a\"\ndone\ncp \"$last\" '{seen}'\ncase '{mode}' in\n  fail) echo 'Error: the profile database is damaged' >&2; echo 'second line' >&2; exit 1;;\n  silent) exit 0;;\nesac\ncp '{fixture}' \"$out\"\n",
        args = args.display(),
        seen = seen.display(),
        fixture = fixture.display(),
    );
    std::fs::write(&script, body).expect("a script");
    std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).expect("executable");
    Installed {
        path: script,
        version: "# HMMER 3.4 (Aug 2023); http://hmmer.org/".into(),
    }
}
