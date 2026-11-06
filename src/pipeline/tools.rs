//! Finding and running other programs (MAFFT, IQ-TREE, primer3, fastp, curl, Snakemake, Nextflow). A program started
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
}

pub const MAFFT: Program = Program {
    label: "MAFFT",
    names: &["mafft"],
    version_args: &["--version"],
};
pub const IQTREE: Program = Program {
    label: "IQ-TREE",
    names: &["iqtree3", "iqtree2", "iqtree"],
    version_args: &["--version"],
};
pub const PRIMER3: Program = Program {
    label: "primer3",
    names: &["primer3_core"],
    version_args: &["-about"],
};
pub const FASTP: Program = Program {
    label: "fastp",
    names: &["fastp"],
    version_args: &["--version"],
};
/// Used for the online structure prediction (a web request, so no network library is linked in).
pub const CURL: Program = Program {
    label: "curl",
    names: &["curl"],
    version_args: &["--version"],
};
pub const SNAKEMAKE: Program = Program {
    label: "Snakemake",
    names: &["snakemake"],
    version_args: &["--version"],
};
pub const NEXTFLOW: Program = Program {
    label: "Nextflow",
    names: &["nextflow"],
    version_args: &["-version"],
};

/// Every program the platform looks for.
pub const ALL: [&Program; 7] = [
    &MAFFT, &IQTREE, &PRIMER3, &FASTP, &CURL, &SNAKEMAKE, &NEXTFLOW,
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

/// The program, if it is installed, with its version.
pub fn detect(program: &Program) -> Option<Installed> {
    let path = locate(program.names)?;
    let output = Command::new(&path)
        .args(program.version_args)
        .output()
        .ok()?;
    // some programs print their version on the error stream
    let text = [output.stdout, output.stderr].concat();
    // the line that says "version" if there is one (Nextflow prints a banner first), else the first line
    let text = String::from_utf8_lossy(&text);
    let lines: Vec<&str> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty())
        .collect();
    let version = lines
        .iter()
        .find(|l| l.to_ascii_lowercase().contains("version"))
        .or(lines.first())
        .map_or(String::new(), |l| (*l).to_string());
    Some(Installed { path, version })
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

/// A browser has no folders and starts no program, so nothing is written where a run would work: an empty path stands in
/// (std's `temp_dir` and `process::id` panic there).
#[cfg(target_arch = "wasm32")]
pub fn work_folder(_label: &str) -> Result<PathBuf, String> {
    Ok(PathBuf::new())
}

/// A new empty folder under the system's temporary folder for one run's files.
#[cfg(not(target_arch = "wasm32"))]
pub fn work_folder(label: &str) -> Result<PathBuf, String> {
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let dir =
        std::env::temp_dir().join(format!("helixloom-{label}-{}-{stamp}", std::process::id()));
    std::fs::create_dir_all(&dir).map_err(|e| format!("Could not make {}: {e}", dir.display()))?;
    Ok(dir)
}

#[cfg(test)]
mod tests {
    use super::*;

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
