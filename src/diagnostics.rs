//! What a person copies into a report of a problem (`docs/ROADMAP.md` 1.8, "copy diagnostics"): the version of the program, the
//! computer it runs on, the programs and model files it found, and the note of the last fault. It is made of things that
//! cannot hold the person's work: no sequence, no name of a file that was opened, no address or token of a workstation, not
//! what the program said when it stopped (only where in its code that was), and the home folder is written `~`, so that the
//! text can be pasted where others read it. What it does show of the computer: where the other programs are (their paths). Nothing is sent anywhere: the text goes to
//! the clipboard, and the person decides where it goes from there.

use std::path::{Path, PathBuf};

/// The file in the settings folder that keeps the note of the last fault.
pub const CRASH_FILE: &str = "last-crash.txt";
/// The most of that file that is read, the most letters of one thing in a note, and the longest a version is.
const MOST_CRASH_BYTES: usize = 4_000;
const MOST_FIELD: usize = 200;
const MOST_VERSION: usize = 24;

/// A program that was looked for: its name, and where it is and in which version when it was found.
pub struct Looked {
    pub label: String,
    pub found: Option<(String, String)>,
}

/// What the text is made from.
pub struct Facts {
    pub version: String,
    pub os: String,
    pub arch: String,
    /// What the program was built with (`gpu`, `cuda`).
    pub features: Vec<String>,
    /// The programs of the registry; none while they are still being looked for.
    pub programs: Option<Vec<Looked>>,
    /// The package managers that were found, by their names.
    pub managers: Vec<String>,
    /// The groups of model files, and whether each is installed.
    pub models: Vec<(String, bool)>,
    /// Which folder the settings are kept in, in words: its path can hold a name.
    pub settings: String,
    /// Whether a workstation is connected, and what it said of itself; never its address.
    pub workstation: String,
    pub appearance: String,
    pub last_crash: Option<String>,
    /// The home folder, which is written `~` wherever it is in a path.
    pub home: Option<String>,
}

/// `text` with the home folder written `~` wherever it stands, alone or at the start of a path, with either kind of
/// slash after it: a path says who the person is otherwise. A folder whose name only begins like the home folder's is left.
pub fn tilde(text: &str, home: Option<&str>) -> String {
    let Some(home) = home
        .map(|h| h.trim_end_matches(['/', '\\']))
        .filter(|h| !h.is_empty())
    else {
        return text.to_string();
    };
    let mut out = String::new();
    let mut rest = text;
    while let Some(at) = rest.find(home) {
        let after = &rest[at + home.len()..];
        // the folder itself or something in it, not `/Users/adam` for `/Users/ada`
        let whole = after
            .chars()
            .next()
            .is_none_or(|c| !(c.is_alphanumeric() || matches!(c, '_' | '-' | '.')));
        out.push_str(&rest[..at]);
        out.push_str(if whole { "~" } else { home });
        rest = after;
    }
    out.push_str(rest);
    out
}

/// A version as a program said it, if it is one: letters, digits and the few signs a version has. A program whose first line
/// is something else (a path, a sentence) has no version here rather than that line in a report.
pub fn plain_version(version: &str) -> String {
    let plain = !version.is_empty()
        && version.len() <= MOST_VERSION
        && version
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '.' | '+' | '-' | '_'));
    if plain {
        version.to_string()
    } else {
        String::new()
    }
}

/// Where in the code a fault was, without the folders of the computer the program was built on: from `src/` on for the
/// program's own files, from the name of the package on for a file of another package.
pub fn place_of(file: &str, line: u32) -> String {
    let file = file.replace('\\', "/");
    let short = match file.find("/registry/src/") {
        // …/registry/src/<the index>/<package>/src/x.rs
        Some(at) => file[at + "/registry/src/".len()..]
            .split_once('/')
            .map_or(file.as_str(), |(_, package)| package),
        None if file.starts_with('/') || file.contains(":/") => file
            .rfind("/src/")
            .map_or("a file of the program", |at| &file[at + 1..]),
        None => file.as_str(),
    };
    format!("{short}:{line}")
}

/// The diagnostics as text, a line for each thing.
pub fn text(f: &Facts) -> String {
    let home = f.home.as_deref();
    let mut out = format!(
        "HelixLoom {}\nSystem: {} on {}\nBuilt with: {}\n",
        f.version,
        f.os,
        f.arch,
        if f.features.is_empty() {
            "no optional features".to_string()
        } else {
            f.features.join(", ")
        }
    );
    out.push_str(&format!("Appearance: {}\n", f.appearance));
    out.push_str(&format!("Settings folder: {}\n", f.settings));
    out.push_str(&format!("Workstation: {}\n", f.workstation));
    out.push_str("Programs:\n");
    match &f.programs {
        None => out.push_str("  still being looked for\n"),
        Some(programs) => {
            for p in programs {
                match &p.found {
                    Some((path, version)) => out.push_str(&format!(
                        "  {}: {}{}\n",
                        p.label,
                        tilde(path, home),
                        if version.is_empty() {
                            String::new()
                        } else {
                            format!(" ({version})")
                        }
                    )),
                    None => out.push_str(&format!("  {}: not found\n", p.label)),
                }
            }
        }
    }
    out.push_str(&format!(
        "Package managers: {}\n",
        if f.managers.is_empty() {
            "none found".to_string()
        } else {
            f.managers.join(", ")
        }
    ));
    out.push_str("Model files:\n");
    for (group, installed) in &f.models {
        out.push_str(&format!(
            "  {group}: {}\n",
            if *installed {
                "installed"
            } else {
                "not installed"
            }
        ));
    }
    match &f.last_crash {
        Some(note) => {
            out.push_str("Last fault:\n");
            for line in note.lines() {
                out.push_str(&format!("  {}\n", tilde(line, home)));
            }
        }
        None => out.push_str("Last fault: none noted\n"),
    }
    out
}

/// The note of a fault: when, in which version, on which thread, and where in the code. What the program said as it
/// stopped is not kept: it can hold a name from the person's work, and the place is what finds the fault.
pub fn crash_note(date: &str, version: &str, thread: &str, place: &str) -> String {
    let one_line = |text: &str| -> String {
        text.chars()
            .map(|c| if c.is_control() { ' ' } else { c })
            .take(MOST_FIELD)
            .collect()
    };
    format!(
        "{}, HelixLoom {}, thread {}\nat {}\n",
        one_line(date),
        one_line(version),
        one_line(thread),
        one_line(place)
    )
}

/// The note of the last fault that is kept in `dir`, if there is one.
pub fn last_crash(dir: &Path) -> Option<String> {
    let bytes = std::fs::read(dir.join(CRASH_FILE)).ok()?;
    let text = String::from_utf8_lossy(&bytes[..bytes.len().min(MOST_CRASH_BYTES)]).into_owned();
    (!text.trim().is_empty()).then_some(text)
}

/// The note of a fault of this version of the program, at `location` (a file as the compiler names it, and a line) on
/// `thread`, `seconds` after the start of 1970.
fn fault_note(location: Option<(&str, u32)>, thread: &str, seconds: u64) -> String {
    let place = location.map_or_else(
        || "an unknown place".to_string(),
        |(file, line)| place_of(file, line),
    );
    crash_note(
        &crate::connect::utc_date(seconds),
        env!("CARGO_PKG_VERSION"),
        thread,
        &place,
    )
}

/// The settings folder in words: its path can hold a name, and that it is not the usual one is what a report needs.
pub fn settings_in_words(kept: bool, named: bool) -> &'static str {
    match (kept, named) {
        (false, _) => "none",
        (true, true) => "a folder named by HELIXLOOM_CONFIG_DIR",
        (true, false) => "the usual folder",
    }
}

/// A program that was looked for, as the text holds it: where it is, and its version when what it printed looks like one.
pub fn looked(label: &str, found: Option<(&Path, &str)>) -> Looked {
    Looked {
        label: label.to_string(),
        found: found.map(|(path, version)| (path.display().to_string(), plain_version(version))),
    }
}

/// From now on a fault (a panic, on any thread) leaves a note in `dir` before the program says what it always said. A fault
/// that only stops one job (the window goes on) leaves its note too. The note is of the last fault only.
pub fn note_faults_in(dir: PathBuf) {
    let before = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let seconds = web_time::SystemTime::now()
            .duration_since(web_time::UNIX_EPOCH)
            .map_or(0, |d| d.as_secs());
        let note = fault_note(
            info.location().map(|l| (l.file(), l.line())),
            std::thread::current().name().unwrap_or("unnamed"),
            seconds,
        );
        // a note that cannot be written is not worth a second fault
        let _ = std::fs::write(dir.join(CRASH_FILE), note);
        before(info);
    }));
}

#[cfg(test)]
mod tests;
