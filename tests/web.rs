//! The web version as `scripts/build-web.sh` packages it, and its service worker as a browser would run it
//! (`docs/WEB.md`, *Installed, and without a network*). The script is given two small files that stand for what
//! wasm-bindgen writes, so nothing is built; the worker is run by node in `tests/web/sw_harness.mjs`.
#![cfg(unix)]

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const ROOT: &str = env!("CARGO_MANIFEST_DIR");

/// A folder of the temporary folder for one test, removed when the test ends.
struct Scratch(PathBuf);

impl Scratch {
    fn new(label: &str) -> Scratch {
        let dir =
            std::env::temp_dir().join(format!("helixloom-web-test-{}-{label}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("a folder");
        Scratch(dir)
    }

    /// The two files that stand for what wasm-bindgen writes, the program being `wasm`; the folder they are in.
    fn prebuilt(&self, wasm: &[u8]) -> PathBuf {
        let prebuilt = self.0.join("prebuilt");
        std::fs::create_dir_all(&prebuilt).expect("a folder");
        // the script that loads the program names it, as the one wasm-bindgen writes does
        let loader = "export default function init(o) { return fetch(new URL('helixloom_bg.wasm', import.meta.url)); }\n";
        std::fs::write(prebuilt.join("helixloom.js"), loader).expect("a script");
        std::fs::write(prebuilt.join("helixloom_bg.wasm"), wasm).expect("a program");
        prebuilt
    }

    /// A copy of the script and of the files it reads, to try what the script must refuse without the repository's own
    /// files being at stake; the root of the copy.
    fn copy_of_the_repository(&self) -> PathBuf {
        let copy = self.0.join("repository");
        for folder in ["scripts", "web"] {
            std::fs::create_dir_all(copy.join(folder)).expect("a folder");
            for entry in std::fs::read_dir(Path::new(ROOT).join(folder)).expect("the folder") {
                let entry = entry.expect("an entry");
                std::fs::copy(entry.path(), copy.join(folder).join(entry.file_name()))
                    .expect("a copy");
            }
        }
        copy
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Run the script of the repository at `root` to put the page in `out`, from the stand-ins in `prebuilt`.
fn run_script(root: &Path, out: &Path, prebuilt: &Path) -> Output {
    Command::new("sh")
        .arg(root.join("scripts/build-web.sh"))
        .arg(out)
        .env("HELIXLOOM_WEB_PREBUILT", prebuilt)
        .output()
        .expect("sh runs the script")
}

fn said(run: &Output) -> String {
    format!(
        "{}{}",
        String::from_utf8_lossy(&run.stdout),
        String::from_utf8_lossy(&run.stderr)
    )
}

/// Package a page whose program is `wasm` into the folder `dist` of `scratch`; that folder.
fn packaged(scratch: &Scratch, wasm: &[u8]) -> PathBuf {
    let out = scratch.0.join("dist");
    let run = run_script(Path::new(ROOT), &out, &scratch.prebuilt(wasm));
    assert!(run.status.success(), "{}", said(&run));
    out
}

fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .expect("the folder")
        .map(|entry| {
            entry
                .expect("an entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    names.sort();
    names
}

fn text(dir: &Path, name: &str) -> String {
    std::fs::read_to_string(dir.join(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
}

/// The one file of `dir` that the script named after a checksum: `begins`, sixteen letters of the checksum, `ends`.
fn the_one(dir: &Path, begins: &str, ends: &str) -> String {
    let found: Vec<String> = names(dir)
        .into_iter()
        .filter(|name| {
            name.strip_prefix(begins)
                .and_then(|rest| rest.strip_suffix(ends))
                .is_some_and(|sum| sum.len() == 16 && sum.bytes().all(|b| b.is_ascii_hexdigit()))
        })
        .collect();
    assert_eq!(found.len(), 1, "{begins}…{ends} among {:?}", names(dir));
    found[0].clone()
}

/// The files the worker says it keeps, and the version it says it is.
fn kept_by(worker: &str) -> (Vec<String>, String) {
    let list = worker
        .split_once("const FILES = [")
        .and_then(|(_, rest)| rest.split_once(']'))
        .expect("a list of files")
        .0;
    let mut files: Vec<String> = list
        .split(',')
        .map(|file| file.trim().trim_matches('"').to_string())
        .filter(|file| !file.is_empty())
        .collect();
    files.sort();
    let version = worker
        .split_once("const VERSION = \"")
        .and_then(|(_, rest)| rest.split_once('"'))
        .expect("a version")
        .0;
    (files, version.to_string())
}

/// What a packed copy holds once it is unpacked.
fn unpacked(packed: &Path) -> Vec<u8> {
    let run = Command::new("gzip")
        .arg("-dc")
        .arg(packed)
        .output()
        .expect("gzip unpacks");
    assert!(run.status.success(), "{}", packed.display());
    run.stdout
}

#[test]
fn the_page_is_packaged_with_a_manifest_its_icons_and_a_worker_that_names_every_file_of_the_build()
{
    let scratch = Scratch::new("package");
    let out = packaged(&scratch, b"\0asm one");
    let script = the_one(&out, "helixloom-", ".js");
    let program = the_one(&out, "helixloom-", "_bg.wasm");
    let style = the_one(&out, "style-", ".css");

    // the page names its files, the manifest and an icon
    let page = text(&out, "index.html");
    for needed in [
        script.as_str(),
        program.as_str(),
        style.as_str(),
        "rel=\"manifest\" href=\"./manifest.webmanifest\"",
        "href=\"./icon-192.png\"",
    ] {
        assert!(page.contains(needed), "{needed} in\n{page}");
    }
    // the script names the program by its new name
    assert!(text(&out, &script).contains(&program));

    // the worker keeps the page and every file the page needs, and nothing else: no packed copy, not itself
    let worker = text(&out, "sw.js");
    assert!(
        !worker.contains("__"),
        "nothing is left to be filled in:\n{worker}"
    );
    let (kept, version) = kept_by(&worker);
    let mut wanted = vec![
        "./".to_string(),
        format!("./{script}"),
        format!("./{program}"),
        format!("./{style}"),
        "./manifest.webmanifest".to_string(),
        "./icon-192.png".to_string(),
        "./icon-512.png".to_string(),
    ];
    wanted.sort();
    assert_eq!(kept, wanted);
    assert!(
        version.len() == 16 && version.bytes().all(|b| b.is_ascii_hexdigit()),
        "{version}"
    );
    for file in kept.iter().filter(|file| *file != "./") {
        assert!(out.join(file).is_file(), "{file} is in the folder");
    }

    // the manifest and the icons are the repository's, as they are
    for same in ["manifest.webmanifest", "icon-192.png", "icon-512.png"] {
        assert_eq!(
            std::fs::read(out.join(same)).expect("packaged"),
            std::fs::read(Path::new(ROOT).join("web").join(same)).expect("in the repository"),
            "{same}"
        );
    }
    // a packed copy beside each text and the program, and it is the copy of that file; a picture is packed already
    let all = names(&out);
    for packed in [
        "index.html",
        "sw.js",
        "manifest.webmanifest",
        script.as_str(),
        program.as_str(),
        style.as_str(),
    ] {
        assert_eq!(
            unpacked(&out.join(format!("{packed}.gz"))),
            std::fs::read(out.join(packed)).expect("the file"),
            "{packed}.gz"
        );
    }
    assert!(!all.iter().any(|name| name.ends_with(".png.gz")), "{all:?}");
    // and nothing is in the folder that the script did not mean to put there
    assert_eq!(
        all.len(),
        8 + 6,
        "eight files and six packed copies: {all:?}"
    );
}

#[test]
fn another_build_has_other_names_and_another_version_and_takes_away_only_what_the_script_wrote() {
    let scratch = Scratch::new("again");
    let out = packaged(&scratch, b"\0asm one");
    let first = the_one(&out, "helixloom-", "_bg.wasm");
    let first_script = the_one(&out, "helixloom-", ".js");
    let (_, version_one) = kept_by(&text(&out, "sw.js"));
    // files of the person's own, some of them named much as the script names its files
    let own = [
        "notes.txt",
        "helixloom-notes.js",
        "helixloom-old_bg.wasm",
        "style-mine.css",
        "helixloom-notes.js.gz",
    ];
    for name in own {
        std::fs::write(out.join(name), "someone's own file").expect("a file");
    }
    // and a packed copy that is no longer the copy of anything
    std::fs::write(out.join("sw.js.gz"), "left over").expect("a file");

    packaged(&scratch, b"\0asm two");
    let second = the_one(&out, "helixloom-", "_bg.wasm");
    let (_, version_two) = kept_by(&text(&out, "sw.js"));
    assert_ne!(first, second, "another program has another name");
    assert_ne!(
        version_one, version_two,
        "and the worker is another version"
    );
    assert!(!out.join(&first).exists() && !out.join(format!("{first}.gz")).exists());
    // the script that loads the program names the program, so it is another file too, and the old one is gone
    assert!(!out.join(&first_script).exists() && !out.join(format!("{first_script}.gz")).exists());
    for name in own {
        assert_eq!(text(&out, name), "someone's own file", "{name}");
    }
    assert_eq!(
        unpacked(&out.join("sw.js.gz")),
        std::fs::read(out.join("sw.js")).expect("the worker")
    );

    // the same program again: the same names, the same version, so a browser sees that nothing is new
    packaged(&scratch, b"\0asm two");
    assert_eq!(the_one(&out, "helixloom-", "_bg.wasm"), second);
    assert_eq!(kept_by(&text(&out, "sw.js")).1, version_two);
}

#[test]
fn the_script_says_when_it_builds_nothing_and_takes_a_folder_named_from_where_it_was_called() {
    let scratch = Scratch::new("said");
    scratch.prebuilt(b"\0asm one");
    // called from the scratch folder, with both folders named from there
    let run = Command::new("sh")
        .arg(Path::new(ROOT).join("scripts/build-web.sh"))
        .arg("out")
        .env("HELIXLOOM_WEB_PREBUILT", "prebuilt")
        .current_dir(&scratch.0)
        .output()
        .expect("sh runs the script");
    assert!(run.status.success(), "{}", said(&run));
    assert!(scratch.0.join("out/index.html").is_file());
    let told = String::from_utf8_lossy(&run.stderr).into_owned();
    assert!(
        told.contains("HELIXLOOM_WEB_PREBUILT is set: nothing is built"),
        "a page of stand-in files is not passed off as a build: {told}"
    );
}

#[test]
fn the_script_refuses_the_folders_its_own_sources_are_in_and_leaves_them_as_they_are() {
    let scratch = Scratch::new("refuses");
    let copy = scratch.copy_of_the_repository();
    let prebuilt = scratch.prebuilt(b"\0asm one");
    let before = names(&copy.join("web"));
    for own in [copy.join("web"), copy.clone()] {
        let run = run_script(&copy, &own, &prebuilt);
        assert!(!run.status.success(), "{}: {}", own.display(), said(&run));
        assert!(
            said(&run).contains("the script's own sources are there"),
            "{}",
            said(&run)
        );
        assert_eq!(names(&copy.join("web")), before, "{}", own.display());
        assert!(!copy.join("index.html").exists() && !copy.join("sw.js").exists());
    }
    // and any other folder of the copy is taken
    let run = run_script(&copy, &copy.join("dist"), &prebuilt);
    assert!(run.status.success(), "{}", said(&run));
}

#[test]
fn a_page_that_cannot_be_put_together_leaves_the_page_that_was_there_whole() {
    let scratch = Scratch::new("whole");
    let copy = scratch.copy_of_the_repository();
    let out = copy.join("dist");
    let run = run_script(&copy, &out, &scratch.prebuilt(b"\0asm one"));
    assert!(run.status.success(), "{}", said(&run));
    let there = names(&out);
    let page = text(&out, "index.html");
    // something the page needs is gone from the sources, and the program is another
    std::fs::remove_file(copy.join("web/icon-512.png")).expect("removed");
    let run = run_script(&copy, &out, &scratch.prebuilt(b"\0asm two"));
    assert!(!run.status.success(), "{}", said(&run));
    assert_eq!(
        names(&out),
        there,
        "nothing was taken away and nothing was added"
    );
    assert_eq!(text(&out, "index.html"), page);
}

#[test]
fn the_worker_is_another_version_when_only_a_file_that_keeps_its_name_is_another() {
    let scratch = Scratch::new("version");
    let copy = scratch.copy_of_the_repository();
    let out = copy.join("dist");
    let prebuilt = scratch.prebuilt(b"\0asm one");
    let version = |what: &str| {
        let run = run_script(&copy, &out, &prebuilt);
        assert!(run.status.success(), "{what}: {}", said(&run));
        kept_by(&text(&out, "sw.js")).1
    };
    let mut seen = vec![version("as it is")];
    // the manifest, an icon and the worker itself keep their names from build to build: a browser learns that one of
    // them is new only from the worker being another version
    let manifest = copy.join("web/manifest.webmanifest");
    std::fs::write(
        &manifest,
        text(&copy.join("web"), "manifest.webmanifest").replace("standalone", "minimal-ui"),
    )
    .expect("written");
    seen.push(version("another manifest"));
    let mut icon = std::fs::read(copy.join("web/icon-192.png")).expect("the icon");
    icon.push(0);
    std::fs::write(copy.join("web/icon-192.png"), icon).expect("written");
    seen.push(version("another icon"));
    let worker = copy.join("web/sw.js");
    std::fs::write(
        &worker,
        format!("{}// one more line\n", text(&copy.join("web"), "sw.js")),
    )
    .expect("written");
    seen.push(version("another worker"));
    for (at, one) in seen.iter().enumerate() {
        assert!(
            !seen[..at].contains(one),
            "the versions are all different: {seen:?}"
        );
    }
}

#[test]
fn the_icons_are_pictures_of_the_sizes_their_names_say() {
    for (name, size) in [("icon-192.png", 192u32), ("icon-512.png", 512)] {
        let bytes = std::fs::read(Path::new(ROOT).join("web").join(name)).expect("the icon");
        assert_eq!(&bytes[..8], b"\x89PNG\r\n\x1a\n", "{name} is a PNG");
        let number =
            |at: usize| u32::from_be_bytes(bytes[at..at + 4].try_into().expect("four bytes"));
        assert_eq!((number(16), number(20)), (size, size), "{name}");
    }
}

/// What the worker and the page must hold, as `tests/web/sw_harness.mjs` names them: each is looked for by name, so that
/// one that is taken out of the harness is missed here.
const HELD: [&str; 12] = [
    "it keeps every file of its build, asked of the server itself, in a cache of its page and its build",
    "a build of which one file cannot be had is not kept at all",
    "a build that the browser has no room to keep whole leaves nothing kept",
    "a page that the server sends by way of a redirect is kept as an answer of its own",
    "with a network the page is the newest; without one it is the kept one",
    "without a network and with nothing kept the page is what a browser shows when it reaches nothing",
    "a file of the build comes from what is kept, and the network is not asked",
    "a file that is kept by no cache yet is fetched",
    "everything else is left to the browser: above all a call to a workstation",
    "when it takes over it removes the caches of this page's older builds, and no other cache",
    "the manifest is what a browser needs to install the page",
    "the page asks for the worker once it has loaded, and a browser without workers is left in peace",
];

#[test]
fn the_worker_keeps_the_page_for_use_without_a_network_and_leaves_every_other_request_alone() {
    let node = Command::new("node").arg("--version").output();
    if !node.is_ok_and(|out| out.status.success()) {
        // where the tests are run for everyone (the variable CI is set there), a worker that was not run is a failure;
        // on a person's own computer it is said, and the other tests of this file still hold the packaging
        assert!(
            std::env::var_os("CI").is_none(),
            "node was not found: the service worker cannot be tested (tests/web/sw_harness.mjs)"
        );
        eprintln!("node was not found: the service worker was not run (tests/web/sw_harness.mjs)");
        return;
    }
    let scratch = Scratch::new("worker");
    let out = packaged(&scratch, b"\0asm one");
    let run = Command::new("node")
        .arg(Path::new(ROOT).join("tests/web/sw_harness.mjs"))
        .arg(out.join("sw.js"))
        .arg(out.join("manifest.webmanifest"))
        .arg(out.join("index.html"))
        .output()
        .expect("node runs the harness");
    let told = said(&run);
    assert!(run.status.success(), "{told}");
    for held in HELD {
        assert!(
            told.contains(&format!("held: {held}\n")),
            "{held}\nin\n{told}"
        );
    }
    assert!(told.contains("ALL HELD"), "{told}");
}
