//! Finding and running other programs (MAFFT, IQ-TREE, HMMER, primer3, fastp, curl, Snakemake, Nextflow). A program started
//! from the Finder gets a short PATH, so the folders Homebrew and MacPorts use are searched as well.

use std::path::{Path, PathBuf};
use std::process::Command;

/// Searched after PATH.
const EXTRA_FOLDERS: [&str; 4] = [
    "/opt/homebrew/bin",
    "/usr/local/bin",
    "/usr/bin",
    "/opt/local/bin",
];

/// A program the workflows can call, by the names it goes by.
pub struct Program {
    pub label: &'static str,
    pub names: &'static [&'static str],
    /// Arguments that make it print its version.
    pub version_args: &'static [&'static str],
    /// What this program uses it for.
    pub role: &'static str,
    /// Its licence as Homebrew's formula lists it (November 2026): check the program's own before passing it on.
    pub licence: &'static str,
    /// The name of its Homebrew formula and of its Bioconda package; empty where there is none to name.
    pub homebrew: &'static str,
    pub conda: &'static str,
}

pub const MAFFT: Program = Program {
    label: "MAFFT",
    names: &["mafft"],
    version_args: &["--version"],
    role: "multiple sequence alignment (Alignment, Pipelines)",
    licence: "BSD-3-Clause",
    homebrew: "mafft",
    conda: "mafft",
};
pub const IQTREE: Program = Program {
    label: "IQ-TREE",
    names: &["iqtree3", "iqtree2", "iqtree"],
    version_args: &["--version"],
    role: "maximum-likelihood trees with model finding and ultrafast bootstrap (Pipelines)",
    licence: "GPL-2.0-or-later",
    homebrew: "iqtree3",
    conda: "iqtree",
};
/// HMMER 3's `hmmscan`, which searches protein sequences against profile HMMs (`pipeline::domains`). It has no `--version`: `-h`
/// prints a banner whose second line is `# HMMER 3.4 (Aug 2023); http://hmmer.org/`.
pub const HMMSCAN: Program = Program {
    label: "HMMER",
    names: &["hmmscan"],
    version_args: &["-h"],
    role: "profile HMM searches against Pfam (Protein Domains)",
    licence: "BSD-3-Clause",
    homebrew: "hmmer",
    conda: "hmmer",
};
pub const PRIMER3: Program = Program {
    label: "primer3",
    names: &["primer3_core"],
    version_args: &["-about"],
    role: "primer design (Pipelines, Primer3)",
    licence: "GPL-2.0-or-later and GPL-3.0-or-later",
    homebrew: "primer3",
    conda: "primer3",
};
pub const FASTP: Program = Program {
    label: "fastp",
    names: &["fastp"],
    version_args: &["--version"],
    role: "cleaning sequencing reads (Pipelines)",
    licence: "MIT",
    homebrew: "fastp",
    conda: "fastp",
};
/// Used for the online structure prediction (a web request, so no network library is linked in).
pub const CURL: Program = Program {
    label: "curl",
    names: &["curl"],
    version_args: &["--version"],
    role: "fetching records from databases and the online structure prediction",
    licence: "curl (a permissive licence of its own)",
    homebrew: "curl",
    conda: "",
};
pub const SNAKEMAKE: Program = Program {
    label: "Snakemake",
    names: &["snakemake"],
    version_args: &["--version"],
    role: "running exported Snakemake workflows",
    licence: "MIT",
    homebrew: "snakemake",
    conda: "snakemake",
};
pub const NEXTFLOW: Program = Program {
    label: "Nextflow",
    names: &["nextflow"],
    version_args: &["-version"],
    role: "running exported Nextflow workflows",
    licence: "Apache-2.0",
    homebrew: "nextflow",
    conda: "nextflow",
};

/// Every program the platform looks for.
pub const ALL: [&Program; 8] = [
    &MAFFT, &IQTREE, &HMMSCAN, &PRIMER3, &FASTP, &CURL, &SNAKEMAKE, &NEXTFLOW,
];

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Installed {
    pub path: PathBuf,
    /// The first line the program prints for its version.
    pub version: String,
}

/// The part of a version line worth showing: the first word that looks like a number ("v7.526", "3.1.4").
pub fn short_version(version: &str) -> &str {
    version
        .split_whitespace()
        .find(|t| t.contains('.') && t.chars().any(|c| c.is_ascii_digit()))
        .unwrap_or("")
}

impl Installed {
    /// The program's name and short version, for saying what made a result: "MAFFT v7.526".
    pub fn describe(&self, label: &str) -> String {
        match short_version(&self.version) {
            "" => label.to_string(),
            v => format!("{label} {v}"),
        }
    }
}

#[cfg(unix)]
fn can_run(path: &Path) -> bool {
    use std::os::unix::fs::PermissionsExt;
    std::fs::metadata(path).is_ok_and(|m| m.is_file() && m.permissions().mode() & 0o111 != 0)
}

#[cfg(not(unix))]
fn can_run(path: &Path) -> bool {
    path.is_file()
}

/// Where the first of `names` is, looking in PATH, then Homebrew's and MacPorts' folders, then `~/.local/bin`.
pub fn locate(names: &[&str]) -> Option<PathBuf> {
    // a browser starts no program (and std's `split_paths` panics there): none is found
    if cfg!(target_arch = "wasm32") {
        return None;
    }
    let path = std::env::var_os("PATH").unwrap_or_default();
    let mut folders: Vec<PathBuf> = std::env::split_paths(&path).collect();
    folders.extend(EXTRA_FOLDERS.iter().map(PathBuf::from));
    if let Some(home) = std::env::var_os("HOME") {
        folders.push(Path::new(&home).join(".local/bin"));
    }
    names
        .iter()
        .flat_map(|name| folders.iter().map(move |dir| dir.join(name)))
        .find(|candidate| can_run(candidate))
}

/// The commands that install `p`, as a person types them: one for Homebrew and one for conda, where it has a package there.
pub fn install_ways(p: &Program) -> Vec<String> {
    let mut ways = Vec::new();
    if !p.homebrew.is_empty() {
        ways.push(format!("brew install {}", p.homebrew));
    }
    if !p.conda.is_empty() {
        ways.push(format!(
            "conda install -c bioconda -c conda-forge {}",
            p.conda
        ));
    }
    ways
}

/// The program of the registry that `name` means: its label or the name of one of its programs, in any case.
pub fn by_name(name: &str) -> Option<&'static Program> {
    let name = name.to_ascii_lowercase();
    ALL.into_iter()
        .find(|p| p.label.to_ascii_lowercase() == name || p.names.contains(&name.as_str()))
}

/// What the `tools` command and the workstation page say about one program: what it is for, whether it is found and where, its
/// licence, and how to install it. Nothing here installs anything.
pub fn describe_program(p: &Program, found: Option<&Installed>) -> String {
    let mut text = format!("{} — {}\n", p.label, p.role);
    match found {
        Some(f) => {
            let version = short_version(&f.version);
            text.push_str(&format!(
                "  found:   {}{}\n",
                f.path.display(),
                if version.is_empty() {
                    String::new()
                } else {
                    format!(" (version {version})")
                }
            ));
        }
        None => text.push_str("  found:   no\n"),
    }
    text.push_str(&format!("  licence: {}\n", p.licence));
    let ways = install_ways(p);
    if !ways.is_empty() {
        text.push_str(&format!("  install: {}\n", ways.join("   or   ")));
    }
    text
}

/// The programs HelixLoom looks for, as `describe_program` says them, with a line at the end that counts the ones found.
pub fn report(found: &[(&Program, Option<Installed>)]) -> String {
    let mut text = String::from(
        "HelixLoom looks for these programs on the PATH (and in Homebrew's folders); it installs one only when asked to \
         (`helixloom tools install NAME`).\n\n",
    );
    for (p, f) in found {
        text.push_str(&describe_program(p, f.as_ref()));
        text.push('\n');
    }
    text.push_str(&format!(
        "{} of {} found. Licences are as Homebrew lists them (November 2026): check a program's own before passing it on.\n",
        found.iter().filter(|(_, f)| f.is_some()).count(),
        found.len()
    ));
    text
}

/// The program, if it is installed, with its version.
pub fn detect(program: &Program) -> Option<Installed> {
    let path = locate(program.names)?;
    let output = Command::new(&path)
        .args(program.version_args)
        .output()
        .ok()?;
    // some programs print their version on the error stream
    let text = [output.stdout, output.stderr].concat();
    Some(Installed {
        path,
        version: version_line(&String::from_utf8_lossy(&text)),
    })
}

/// The line of a program's help or version output that says its version: the line that says "version" if there is one (Nextflow
/// prints a banner first), or HMMER's `# HMMER 3.4 (Aug 2023)` (its banner has no such word), else the first line.
fn version_line(text: &str) -> String {
    let lines: Vec<&str> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .collect();
    lines
        .iter()
        .find(|l| l.to_ascii_lowercase().contains("version") || l.starts_with("# HMMER"))
        .or(lines.first())
        .map_or(String::new(), |l| (*l).to_string())
}

/// Run `program` with `args` in `dir`; what it printed on its output stream, or why it failed.
pub fn run(program: &Path, args: &[&str], dir: &Path) -> Result<String, String> {
    let output = Command::new(program)
        .args(args)
        .current_dir(dir)
        .output()
        .map_err(|e| format!("Could not start {}: {e}", program.display()))?;
    if output.status.success() {
        return Ok(String::from_utf8_lossy(&output.stdout).into_owned());
    }
    let errors = String::from_utf8_lossy(&output.stderr);
    let tail: Vec<&str> = errors
        .lines()
        .rev()
        .filter(|l| !l.trim().is_empty())
        .take(3)
        .collect();
    Err(format!(
        "{} stopped ({}): {}",
        program
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("the program"),
        output.status,
        tail.into_iter().rev().collect::<Vec<_>>().join(" | ")
    ))
}

/// The folder of one run's files, which is removed with everything the run wrote into it (the person's sequences are among
/// it) when the run is over, however it ends: with its answer, with a failure, or because the program is stopping. A run
/// keeps this for as long as it reads or writes there, and nothing more is needed to clean up. (Until November 2026 the
/// runs were given a bare path and three of them never removed it: every alignment with MAFFT, every pipeline and every
/// pipeline file left its sequences in the temporary folder, readable by everyone on the computer.)
pub struct WorkFolder(PathBuf);

impl WorkFolder {
    /// A new empty folder under the system's temporary folder, which only the person who runs the program may look into.
    pub fn new(label: &str) -> Result<WorkFolder, String> {
        make_work_folder(label).map(WorkFolder)
    }

    /// For a run that writes no files: no folder is made, and none is removed.
    pub fn none() -> WorkFolder {
        WorkFolder(PathBuf::new())
    }
}

impl std::ops::Deref for WorkFolder {
    type Target = Path;

    fn deref(&self) -> &Path {
        &self.0
    }
}

impl Drop for WorkFolder {
    fn drop(&mut self) {
        // a run without files, and every run in a browser, has an empty path: there is nothing to remove, and removing
        // "" must never be tried
        if self.0.as_os_str().is_empty() {
            return;
        }
        // a folder that will not go cannot be told to anyone from here; the run's answer is what its caller is given
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A folder of the temporary folder that a test makes and looks after itself. A run uses `WorkFolder`, which is removed.
#[cfg(test)]
pub fn work_folder(label: &str) -> Result<PathBuf, String> {
    make_work_folder(label)
}

/// A browser has no folders and starts no program, so nothing is written where a run would work: an empty path stands in
/// (std's `temp_dir` and `process::id` panic there).
#[cfg(target_arch = "wasm32")]
fn make_work_folder(_label: &str) -> Result<PathBuf, String> {
    Ok(PathBuf::new())
}

/// A new empty folder under the system's temporary folder for one run's files.
#[cfg(not(target_arch = "wasm32"))]
fn make_work_folder(label: &str) -> Result<PathBuf, String> {
    use std::sync::atomic::{AtomicU64, Ordering};
    // a counter beside the clock: threads that ask in the same tick still get folders of their own, and `create_dir` (not
    // `create_dir_all`) fails rather than share one that is there
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let dir = std::env::temp_dir().join(format!(
        "helixloom-{label}-{}-{stamp}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    // not `create_dir_all`: a folder that is there is somebody else's. The temporary folder is everyone's on many
    // computers, and a run writes sequences here: only the person who runs it may look in
    let mut builder = std::fs::DirBuilder::new();
    #[cfg(unix)]
    {
        use std::os::unix::fs::DirBuilderExt;
        builder.mode(0o700);
    }
    builder
        .create(&dir)
        .map_err(|e| format!("Could not make {}: {e}", dir.display()))?;
    Ok(dir)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_folder_of_a_run_goes_with_what_the_run_wrote_when_the_run_is_over_however_it_ends() {
        let run = |fail: bool| -> Result<PathBuf, (PathBuf, String)> {
            let dir = WorkFolder::new("guard").expect("a folder");
            let path = dir.to_path_buf();
            assert!(path.is_dir());
            std::fs::write(dir.join("sequences.fasta"), ">someone's\nACGT\n").expect("written");
            std::fs::create_dir(dir.join("deeper")).expect("a folder in it");
            std::fs::write(dir.join("deeper").join("more.txt"), "more").expect("written");
            if fail {
                return Err((path, "the program failed".into()));
            }
            Ok(path)
        };
        let ended_well = run(false).expect("a run");
        assert!(
            !ended_well.exists(),
            "{} is left behind",
            ended_well.display()
        );
        let (failed, _) = run(true).expect_err("a failed run");
        assert!(!failed.exists(), "{} is left behind", failed.display());
    }

    #[test]
    fn a_run_that_writes_no_files_has_no_folder_and_removes_none() {
        let none = WorkFolder::none();
        assert!(none.as_os_str().is_empty());
        // dropping it must not remove the folder the program happens to stand in
        let here = std::env::current_dir().expect("a folder to stand in");
        drop(none);
        assert!(here.is_dir());
    }

    #[cfg(unix)]
    #[test]
    fn only_the_person_who_runs_it_may_look_into_the_folder_of_a_run() {
        use std::os::unix::fs::PermissionsExt;
        let dir = WorkFolder::new("private").expect("a folder");
        let mode = std::fs::metadata(&*dir)
            .expect("it is there")
            .permissions()
            .mode();
        assert_eq!(mode & 0o777, 0o700, "{mode:o}");
    }

    #[test]
    fn every_program_says_what_it_is_for_its_licence_and_how_to_get_it() {
        let mut labels: Vec<&str> = ALL.iter().map(|p| p.label).collect();
        labels.sort_unstable();
        labels.dedup();
        assert_eq!(labels.len(), ALL.len(), "two programs have one label");
        for p in ALL {
            assert!(!p.role.is_empty() && !p.licence.is_empty(), "{}", p.label);
            // a program with nothing to install it by must at least be one that is usually there
            assert!(!p.homebrew.is_empty() || p.label == "curl", "{}", p.label);
            // the names go into a line for people to copy: plain
            assert!(
                p.homebrew
                    .chars()
                    .chain(p.conda.chars())
                    .all(|c| c.is_ascii_alphanumeric() || c == '-'),
                "{}",
                p.label
            );
        }
    }

    #[test]
    fn a_program_is_described_found_or_not_with_its_licence_and_the_ways_to_install_it() {
        let found = Installed {
            path: PathBuf::from("/opt/homebrew/bin/hmmscan"),
            version: "# HMMER 3.4 (Aug 2023); http://hmmer.org/".into(),
        };
        assert_eq!(
            describe_program(&HMMSCAN, Some(&found)),
            "HMMER — profile HMM searches against Pfam (Protein Domains)\n  found:   /opt/homebrew/bin/hmmscan (version 3.4)\n  licence: BSD-3-Clause\n  install: brew install hmmer   or   conda install -c bioconda -c conda-forge hmmer\n"
        );
        let missing = describe_program(&FASTP, None);
        assert!(
            missing.contains("found:   no") && missing.contains("licence: MIT"),
            "{missing}"
        );
        // curl has no Bioconda package to name, so only the one way is given
        let curl = describe_program(&CURL, None);
        assert!(
            curl.contains("install: brew install curl\n") && !curl.contains("conda"),
            "{curl}"
        );
        // the report counts what is found
        let all: Vec<(&Program, Option<Installed>)> = vec![(&HMMSCAN, Some(found)), (&FASTP, None)];
        let text = report(&all);
        assert!(
            text.starts_with("HelixLoom looks for these programs")
                && text.contains("1 of 2 found."),
            "{text}"
        );
    }

    #[test]
    fn work_folders_asked_for_at_the_same_moment_are_all_different() {
        let made: Vec<PathBuf> = std::thread::scope(|s| {
            let jobs: Vec<_> = (0..16)
                .map(|_| s.spawn(|| work_folder("race").expect("a folder")))
                .collect();
            jobs.into_iter().map(|j| j.join().unwrap()).collect()
        });
        let mut sorted = made.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(
            sorted.len(),
            made.len(),
            "two threads were given one folder"
        );
        for dir in &made {
            assert!(dir.is_dir());
            let _ = std::fs::remove_dir(dir);
        }
    }

    #[test]
    fn a_version_shows_the_number_and_nothing_else() {
        assert_eq!(short_version("v7.526 (2024/Apr/26)"), "v7.526");
        assert_eq!(
            short_version(
                "IQ-TREE single-core version 3.1.4 for MacOS ARM 64-bit built Sep 10 2026"
            ),
            "3.1.4"
        );
        assert_eq!(
            short_version("curl 8.7.1 (x86_64-apple-darwin23.0) libcurl/8.7.1"),
            "8.7.1"
        );
        assert_eq!(short_version(""), "");
        assert_eq!(short_version("no number here"), "");
    }

    #[test]
    fn hmmer_says_its_version_on_the_second_line_of_its_help() {
        // what `hmmscan -h` begins with (HMMER 3.4 User's Guide): the version is found, not the line before it
        let banner = "# hmmscan :: search sequence(s) against a profile database\n# HMMER 3.4 (Aug 2023); http://hmmer.org/\n# Copyright (C) 2023 Howard Hughes Medical Institute.\n";
        let found = version_line(banner);
        assert_eq!(found, "# HMMER 3.4 (Aug 2023); http://hmmer.org/");
        assert_eq!(short_version(&found), "3.4");
        // the other programs' output is read as before: the line that says "version", else the first line
        assert_eq!(
            version_line("N E X T F L O W\n  version 24.10.0 build 5936\n"),
            "version 24.10.0 build 5936"
        );
        assert_eq!(
            version_line("\nv7.526 (2024/Apr/26)\n"),
            "v7.526 (2024/Apr/26)"
        );
        assert_eq!(version_line(""), "");
        assert_eq!(HMMSCAN.names, ["hmmscan"]);
        assert!(ALL.iter().any(|p| p.label == "HMMER"));
    }

    #[test]
    fn a_result_is_credited_to_the_program_and_its_short_version() {
        let found = |version: &str| Installed {
            path: PathBuf::from("/usr/bin/mafft"),
            version: version.to_string(),
        };
        assert_eq!(
            found("v7.526 (2024/Apr/26)").describe("MAFFT"),
            "MAFFT v7.526"
        );
        assert_eq!(found("").describe("MAFFT"), "MAFFT");
    }
}
