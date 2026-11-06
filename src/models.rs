//! Files that belong to other projects and are not part of this program: the trained PhyloCNN networks and the tree that
//! comes with them. Their authors have not given a licence that lets this repository pass them on, so the program looks for
//! them on the computer instead of carrying them (`docs/MODELS.md` says where they go). The same folder is where trained
//! models for later analyses (motifs, domains, epitopes) will live.

use std::path::{Path, PathBuf};

/// A set of files that belong together: the folder they are kept in, inside the models folder, and their names.
pub struct Group {
    pub folder: &'static str,
    /// The files the analysis cannot run without.
    pub files: &'static [&'static str],
    /// Files that only make it nicer (an example): installed with the others, but not waited for.
    pub extras: &'static [&'static str],
}

/// The PhyloCNN networks (Perez & Gascuel 2025) converted to plain binary files, and the Zurich HIV tree of its repository.
pub const PHYLOCNN: Group = Group {
    folder: "phylocnn",
    files: &[
        "Trained_2Generation_PhyDyn.bin",
        "Trained_2Generation_BD.bin",
        "Trained_2Generation_BDEI.bin",
        "Trained_2Generation_BDSS.bin",
        "Trained_2Generation_PhyDyn_HIV.bin",
        "Trained_2Generation_BDSS_HIV.bin",
    ],
    extras: &["ZurichHIV-tree.nwk"],
};

/// The folder that holds the model files the person installed: `models` inside the program's settings folder.
pub fn install_dir() -> Option<PathBuf> {
    crate::remote::settings::dir().map(|d| d.join("models"))
}

/// Where a model file is looked for, first to last: the folder named by `HELIXLOOM_MODELS`, the install folder, and, for
/// a program built from its source folder, the `assets` folder there (where the people who work on it keep the files).
pub fn search_dirs() -> Vec<PathBuf> {
    let named = std::env::var_os("HELIXLOOM_MODELS").filter(|d| !d.is_empty());
    named
        .map(PathBuf::from)
        .into_iter()
        .chain(install_dir())
        .chain(source_tree_assets())
        .collect()
}

/// The `assets` folder of the source tree the program was built from.
#[cfg(not(target_arch = "wasm32"))]
fn source_tree_assets() -> Option<PathBuf> {
    Some(Path::new(env!("CARGO_MANIFEST_DIR")).join("assets"))
}

/// A browser has no folders, and the page must not carry the path of the machine that built it.
#[cfg(target_arch = "wasm32")]
fn source_tree_assets() -> Option<PathBuf> {
    None
}

/// Where `file` of `group` is, in the first of `dirs` that has it; a file with nothing in it is not there.
fn find(dirs: &[PathBuf], group: &Group, file: &str) -> Option<PathBuf> {
    dirs.iter()
        .map(|dir| dir.join(group.folder).join(file))
        .find(|path| std::fs::metadata(path).is_ok_and(|m| m.is_file() && m.len() > 0))
}

/// The bytes of `file` of `group`, from the first of `dirs` that has it.
pub fn read_from(dirs: &[PathBuf], group: &Group, file: &str) -> Result<Vec<u8>, String> {
    if let Some(bytes) = find(dirs, group, file).and_then(|path| std::fs::read(path).ok()) {
        return Ok(bytes);
    }
    Err(format!(
        "The model file {file} is not installed. It is not part of HelixLoom's code (its authors have not allowed that); \
         put it in {} (see docs/MODELS.md).",
        dirs.first()
            .map_or_else(|| group.folder.to_string(), |d| d.join(group.folder).display().to_string())
    ))
}

/// The bytes of `file` of `group`, from the folders the program looks in.
pub fn read(group: &Group, file: &str) -> Result<Vec<u8>, String> {
    // a browser page has no folders to look in
    if cfg!(target_arch = "wasm32") {
        return Err(format!(
            "The model file {file} is not available in the browser yet."
        ));
    }
    read_from(&search_dirs(), group, file)
}

/// The text of `file` of `group`.
pub fn read_text(group: &Group, file: &str) -> Result<String, String> {
    String::from_utf8(read(group, file)?).map_err(|_| format!("{file} is not text."))
}

/// Whether every file that `group` needs can be found (looked at, not read).
pub fn installed(group: &Group) -> bool {
    let dirs = search_dirs();
    group
        .files
        .iter()
        .all(|file| find(&dirs, group, file).is_some())
}

/// Put the files of `group` that `source` has (directly, or in a folder named like the group inside it) into `dest`. Gives
/// how many are there now. A file that already is that very file (the folder chosen was the models folder) is left alone:
/// copying a file onto itself would empty it.
pub fn install_from(source: &Path, dest: &Path, group: &Group) -> Result<usize, String> {
    let inside = source.join(group.folder);
    let from = if inside.is_dir() {
        inside
    } else {
        source.to_path_buf()
    };
    let target = dest.join(group.folder);
    let wanted = group.files.iter().chain(group.extras);
    let found: Vec<&&str> = wanted.filter(|file| from.join(file).is_file()).collect();
    if found.is_empty() {
        return Err(format!(
            "None of the files ({}, …) are in {}.",
            group.files[0],
            from.display()
        ));
    }
    std::fs::create_dir_all(&target)
        .map_err(|e| format!("Could not make {}: {e}", target.display()))?;
    for file in &found {
        let (here, there) = (from.join(file), target.join(file));
        let same = std::fs::canonicalize(&here).ok() == std::fs::canonicalize(&there).ok();
        if there.exists() && same {
            continue;
        }
        std::fs::copy(&here, &there).map_err(|e| format!("Could not copy {file}: {e}"))?;
    }
    Ok(found.len())
}

/// For the tests that need the authors' files: true (with a note) when they are not installed, so that the test can end there.
#[cfg(test)]
pub fn phylocnn_missing() -> bool {
    // the example tree too: the tests that use it need it
    let missing = !installed(&PHYLOCNN) || read_text(&PHYLOCNN, "ZurichHIV-tree.nwk").is_err();
    if missing {
        eprintln!("skipped: the PhyloCNN files are not installed (see docs/MODELS.md)");
    }
    missing
}

#[cfg(test)]
mod tests {
    use super::*;

    const TWO: Group = Group {
        folder: "demo",
        files: &["a.bin", "b.txt"],
        extras: &["c.txt"],
    };

    fn scratch(label: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "helixloom-models-{label}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map_or(0, |d| d.as_nanos())
        ));
        std::fs::create_dir_all(&dir).expect("a scratch folder");
        dir
    }

    #[test]
    fn a_file_is_taken_from_the_first_folder_that_has_it_and_a_missing_one_says_where_to_put_it() {
        let (first, second) = (scratch("first"), scratch("second"));
        std::fs::create_dir_all(second.join("demo")).expect("a folder");
        std::fs::write(second.join("demo/a.bin"), b"two").expect("a file");
        let dirs = vec![first.clone(), second];
        assert_eq!(read_from(&dirs, &TWO, "a.bin").expect("found"), b"two");
        std::fs::create_dir_all(first.join("demo")).expect("a folder");
        std::fs::write(first.join("demo/a.bin"), b"one").expect("a file");
        assert_eq!(read_from(&dirs, &TWO, "a.bin").expect("found"), b"one");
        let why = read_from(&dirs, &TWO, "b.txt").expect_err("not there");
        assert!(why.contains("b.txt") && why.contains("MODELS.md"), "{why}");
        assert!(
            why.contains(&first.join("demo").display().to_string()),
            "{why}"
        );
    }

    #[test]
    fn files_are_installed_from_a_folder_or_from_the_folder_named_like_the_group() {
        let (source, dest) = (scratch("source"), scratch("dest"));
        assert!(install_from(&source, &dest, &TWO).is_err());
        std::fs::write(source.join("b.txt"), b"text").expect("a file");
        assert_eq!(install_from(&source, &dest, &TWO), Ok(1));
        assert_eq!(
            std::fs::read(dest.join("demo/b.txt")).expect("copied"),
            b"text"
        );
        // the folder named like the group wins when there is one
        std::fs::create_dir_all(source.join("demo")).expect("a folder");
        std::fs::write(source.join("demo/a.bin"), b"bin").expect("a file");
        assert_eq!(install_from(&source, &dest, &TWO), Ok(1));
        assert_eq!(
            std::fs::read(dest.join("demo/a.bin")).expect("copied"),
            b"bin"
        );
    }

    #[test]
    fn choosing_the_models_folder_itself_leaves_its_files_as_they_are() {
        let dest = scratch("same");
        let folder = dest.join("demo");
        std::fs::create_dir_all(&folder).expect("a folder");
        std::fs::write(folder.join("a.bin"), b"networks").expect("a file");
        std::fs::write(folder.join("b.txt"), b"text").expect("a file");
        // the folder that holds the files, and the folder that holds that folder
        for source in [folder.clone(), dest.clone()] {
            assert_eq!(install_from(&source, &dest, &TWO), Ok(2));
            assert_eq!(
                std::fs::read(folder.join("a.bin")).expect("kept"),
                b"networks"
            );
            assert_eq!(std::fs::read(folder.join("b.txt")).expect("kept"), b"text");
        }
    }

    #[test]
    fn a_file_with_nothing_in_it_is_not_installed_and_an_extra_is_not_needed() {
        let dir = scratch("empty");
        let (dirs, folder) = (vec![dir.clone()], dir.join("demo"));
        std::fs::create_dir_all(&folder).expect("a folder");
        std::fs::write(folder.join("a.bin"), b"").expect("a file");
        std::fs::write(folder.join("b.txt"), b"text").expect("a file");
        assert!(find(&dirs, &TWO, "a.bin").is_none(), "empty is missing");
        assert!(read_from(&dirs, &TWO, "a.bin").is_err());
        std::fs::write(folder.join("a.bin"), b"bin").expect("a file");
        assert!(TWO.files.iter().all(|f| find(&dirs, &TWO, f).is_some()));
        assert!(
            find(&dirs, &TWO, "c.txt").is_none(),
            "the extra is not there, and nothing waits for it"
        );
    }

    #[test]
    fn the_source_tree_is_the_last_place_looked_in() {
        let dirs = search_dirs();
        assert!(
            dirs.last().is_some_and(|d| d.ends_with("assets")),
            "{dirs:?}"
        );
    }
}
