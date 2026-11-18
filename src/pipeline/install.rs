//! Installing a program that HelixLoom calls, when the person asks for it (`docs/ROADMAP.md` 2.1). The registry knows the name of
//! each program's package (`tools::Program`); this finds the package managers that are on the computer and runs one of them with
//! that name and nothing else. Nothing is installed unless the person says so: the command is shown first, and both the window
//! (Tools → Programs HelixLoom Calls) and `helixloom tools install NAME` ask. A workstation never installs anything on a request:
//! this is for the person at the keyboard, or for the owner of a server at its command line.
//!
//! Homebrew installs for the whole computer. micromamba, mamba and conda install into the environment that is active, so they
//! are used only when one is (`CONDA_PREFIX`): HelixLoom makes no environment and changes none that was not chosen.

use super::tools::{self, Program};
use std::io::{BufRead, BufReader, Read};
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::sync::mpsc;
use std::time::Duration;

/// How many of the last lines a package manager printed are given when it fails.
const TAIL: usize = 3;
/// The most bytes of one line that are passed on in one piece.
const MOST_LINE: u64 = 4_096;
/// How often a package manager that prints nothing is asked whether it has ended, and how long its last lines are waited for
/// once it has.
const BEAT: Duration = Duration::from_millis(100);
const LINGER: Duration = Duration::from_millis(300);
/// The longest that a process which the package manager left behind, and which goes on printing, is listened to.
const DRAIN: Duration = Duration::from_secs(1);

/// The package managers that are looked for, in the order they are preferred.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Manager {
    Homebrew,
    Micromamba,
    Mamba,
    Conda,
}

impl Manager {
    pub const ALL: [Manager; 4] = [
        Manager::Homebrew,
        Manager::Micromamba,
        Manager::Mamba,
        Manager::Conda,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Manager::Homebrew => "Homebrew",
            Manager::Micromamba => "micromamba",
            Manager::Mamba => "mamba",
            Manager::Conda => "conda",
        }
    }

    /// The name of its program.
    pub fn command(self) -> &'static str {
        match self {
            Manager::Homebrew => "brew",
            Manager::Micromamba => "micromamba",
            Manager::Mamba => "mamba",
            Manager::Conda => "conda",
        }
    }

    /// The manager a word names: its label or the name of its program, in any case.
    pub fn parse(word: &str) -> Option<Manager> {
        let word = word.to_ascii_lowercase();
        Manager::ALL
            .into_iter()
            .find(|m| word == m.command() || word == m.label().to_ascii_lowercase())
    }

    /// The name of `program`'s package with this manager, if it has one.
    fn package(self, program: &Program) -> Option<&'static str> {
        let name = match self {
            Manager::Homebrew => program.homebrew,
            Manager::Micromamba | Manager::Mamba | Manager::Conda => program.conda,
        };
        (!name.is_empty()).then_some(name)
    }

    /// Whether it installs into an environment that must be active.
    fn needs_environment(self) -> bool {
        self != Manager::Homebrew
    }
}

/// A package manager that is on this computer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Found {
    pub manager: Manager,
    pub path: PathBuf,
}

/// The package managers on this computer, the preferred first.
pub fn managers() -> Vec<Found> {
    Manager::ALL
        .into_iter()
        .filter_map(|manager| {
            tools::locate(&[manager.command()]).map(|path| Found { manager, path })
        })
        .collect()
}

/// What `manager` is told to install `program`: what follows its own name on the command line. None when it has no package of
/// the program.
pub fn arguments(program: &Program, manager: Manager) -> Option<Vec<String>> {
    let package = manager.package(program)?;
    let words: &[&str] = match manager {
        Manager::Homebrew => &["install"],
        // the packages of these programs are in Bioconda, which stands on conda-forge
        Manager::Micromamba | Manager::Mamba | Manager::Conda => &[
            "install",
            "--yes",
            "--channel",
            "conda-forge",
            "--channel",
            "bioconda",
        ],
    };
    Some(
        words
            .iter()
            .copied()
            .chain([package])
            .map(String::from)
            .collect(),
    )
}

/// What installing a program will run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Plan {
    /// The program's name in the registry.
    pub program: &'static str,
    pub manager: Manager,
    pub path: PathBuf,
    pub args: Vec<String>,
    /// The environment that conda and its kin install into: the one that is active. None for Homebrew.
    pub environment: Option<String>,
}

impl Plan {
    /// The command as a person would type it: `brew install hmmer`.
    pub fn line(&self) -> String {
        format!("{} {}", self.manager.command(), self.args.join(" "))
    }

    /// What the person is told before anything runs: the command, where the program comes from, its licence.
    pub fn notice(&self, program: &Program) -> String {
        let place = match self.manager {
            Manager::Homebrew => "for this computer, from Homebrew's servers".to_string(),
            _ => format!(
                "into the environment that is active ({}), from the conda-forge and Bioconda channels",
                self.environment.as_deref().unwrap_or("none")
            ),
        };
        // the very program that will be started is named: `brew` could be another one than the person thinks of
        format!(
            "To install {} ({}), HelixLoom will run:\n\n    {}\n\nwith the {} at {}. This downloads the program {place}. Its \
             licence is {} as Homebrew lists it; it is another project's program, and HelixLoom only starts the package manager.",
            program.label,
            program.role,
            self.line(),
            self.manager.command(),
            self.path.display(),
            program.licence
        )
    }
}

/// How to install `program` with one of the managers `found`: the one that is `wanted`, or the first that can. `environment` is
/// the environment that is active for conda and its kin (`CONDA_PREFIX`), if one is. When none can, why not and what to do.
pub fn plan(
    program: &Program,
    found: &[Found],
    wanted: Option<Manager>,
    environment: Option<&str>,
) -> Result<Plan, String> {
    let by_hand = tools::install_ways(program);
    let by_hand = if by_hand.is_empty() {
        String::new()
    } else {
        format!(" By hand: {}.", by_hand.join("   or   "))
    };
    if found.is_empty() {
        return Err(format!(
            "No package manager was found on this computer (Homebrew, micromamba, mamba or conda), so HelixLoom cannot install {} \
             for you.{by_hand}",
            program.label
        ));
    }
    if let Some(wanted) = wanted {
        if !found.iter().any(|f| f.manager == wanted) {
            let have: Vec<&str> = found.iter().map(|f| f.manager.label()).collect();
            return Err(format!(
                "{} was not found on this computer; it has {}.",
                wanted.label(),
                have.join(", ")
            ));
        }
    }
    let mut why_not = Vec::new();
    for f in found
        .iter()
        .filter(|f| wanted.is_none_or(|w| w == f.manager))
    {
        let Some(args) = arguments(program, f.manager) else {
            why_not.push(format!(
                "{} has no package of {}",
                f.manager.label(),
                program.label
            ));
            continue;
        };
        if f.manager.needs_environment() && environment.is_none_or(|e| e.trim().is_empty()) {
            why_not.push(format!(
                "{0} installs into the environment that is active, and none is: activate one (`{1} activate NAME`) and start \
                 HelixLoom from it",
                f.manager.label(),
                f.manager.command()
            ));
            continue;
        }
        return Ok(Plan {
            program: program.label,
            manager: f.manager,
            path: f.path.clone(),
            args,
            environment: f
                .manager
                .needs_environment()
                .then(|| environment.map(|e| e.trim().to_string()))
                .flatten(),
        });
    }
    Err(format!(
        "HelixLoom cannot install {} here: {}.{by_hand}",
        program.label,
        why_not.join("; ")
    ))
}

/// Whether an answer to "Run it?" is a yes.
pub fn agreed(answer: &str) -> bool {
    matches!(answer.trim().to_ascii_lowercase().as_str(), "y" | "yes")
}

/// A line that a program printed, without what a terminal would obey: the codes that colour text, move the cursor, set a
/// title or write to the clipboard are taken out, and so is every other control character but the tab.
pub fn clean(line: &str) -> String {
    let mut out = String::new();
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        if c != '\u{1b}' {
            // a tab stays; no other control character does, and no mark that turns the order of the text after it
            let turns = matches!(c, '\u{202a}'..='\u{202e}' | '\u{2066}'..='\u{2069}');
            if c == '\t' || !(c.is_control() || turns) {
                out.push(c);
            }
            continue;
        }
        match chars.peek() {
            // ESC [ … and a letter: colours, the cursor
            Some('[') => {
                chars.next();
                for c in chars.by_ref() {
                    if ('\u{40}'..='\u{7e}').contains(&c) {
                        break;
                    }
                }
            }
            // ESC ] (and P, _, ^, X) … to a bell or to ESC \: a title, the clipboard, strings for the terminal itself
            Some(']' | 'P' | '_' | '^' | 'X') => {
                chars.next();
                while let Some(c) = chars.next() {
                    if c == '\u{7}' {
                        break;
                    }
                    if c == '\u{1b}' {
                        chars.next();
                        break;
                    }
                }
            }
            // any other code is the escape and one character
            Some(_) => {
                chars.next();
            }
            None => {}
        }
    }
    out
}

/// A piece of a stream as it is shown: without the end of its line, cleaned, and of a line that was drawn over and over (a
/// bar of progress, which goes back to the start of its line each time) only what was drawn last.
fn shown(piece: &[u8]) -> String {
    let text = String::from_utf8_lossy(piece);
    let text = text.trim_end_matches(['\n', '\r']);
    clean(text.rsplit('\r').next().unwrap_or(text))
}

/// What one of the package manager's streams prints, a line at a time; a line of any length comes in pieces of about
/// `MOST_LINE` bytes, so that a program that never ends its line cannot fill the memory. A character that the cut of a long
/// line falls in the middle of goes whole into the next piece.
fn pass_on(stream: impl Read, send: &mpsc::Sender<String>) {
    let mut reader = BufReader::new(stream);
    let mut carried: Vec<u8> = Vec::new();
    loop {
        let mut piece = std::mem::take(&mut carried);
        match reader
            .by_ref()
            .take(MOST_LINE)
            .read_until(b'\n', &mut piece)
        {
            Ok(0) | Err(_) => {
                if !piece.is_empty() {
                    let _ = send.send(shown(&piece));
                }
                return;
            }
            Ok(_) => {
                if piece.last() != Some(&b'\n') {
                    if let Err(cut) = std::str::from_utf8(&piece) {
                        if cut.error_len().is_none() {
                            carried = piece.split_off(cut.valid_up_to());
                        }
                    }
                }
                if send.send(shown(&piece)).is_err() {
                    return;
                }
            }
        }
    }
}

/// The package manager while it runs: when the run is left before it has ended (a fault on the way), it is stopped and
/// waited for, so that none is left behind half alive.
struct Running(std::process::Child);

impl Drop for Running {
    fn drop(&mut self) {
        // one that has ended already takes no harm from either
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

/// Run the plan. Every line the package manager prints goes to `say` as it comes (both of its streams, cleaned of what a
/// terminal would obey); when it does not end well, its last lines are the error. When it has ended, what it printed last is
/// read and the run is over: a process it left behind that still holds its output is not waited for.
pub fn run(plan: &Plan, mut say: impl FnMut(&str)) -> Result<(), String> {
    // a browser starts no program
    if cfg!(target_arch = "wasm32") {
        return Err("A page in a browser cannot install a program.".into());
    }
    let mut child = Running(
        Command::new(&plan.path)
            .args(&plan.args)
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| format!("Could not start {}: {e}", plan.path.display()))?,
    );
    let (send, lines) = mpsc::channel::<String>();
    if let Some(stream) = child.0.stdout.take() {
        let send = send.clone();
        std::thread::spawn(move || pass_on(stream, &send));
    }
    if let Some(stream) = child.0.stderr.take() {
        let send = send.clone();
        std::thread::spawn(move || pass_on(stream, &send));
    }
    drop(send);
    let mut last: Vec<String> = Vec::new();
    let mut heard = |line: String, say: &mut dyn FnMut(&str)| {
        say(&line);
        if !line.trim().is_empty() {
            last.push(line);
            if last.len() > TAIL {
                last.remove(0);
            }
        }
    };
    let waited =
        |e: std::io::Error| format!("{} could not be waited for: {e}", plan.manager.command());
    let status = loop {
        match lines.recv_timeout(BEAT) {
            Ok(line) => heard(line, &mut say),
            Err(mpsc::RecvTimeoutError::Timeout) => {}
            // both streams are closed: it has ended, or is about to
            Err(mpsc::RecvTimeoutError::Disconnected) => break child.0.wait().map_err(waited)?,
        }
        if let Some(status) = child.0.try_wait().map_err(waited)? {
            // what it printed last; a process it left behind that goes on printing is listened to for a moment only
            let until = std::time::Instant::now() + DRAIN;
            while std::time::Instant::now() < until {
                match lines.recv_timeout(LINGER) {
                    Ok(line) => heard(line, &mut say),
                    Err(_) => break,
                }
            }
            break status;
        }
    };
    if status.success() {
        return Ok(());
    }
    Err(format!(
        "{} stopped ({status}){}",
        plan.manager.command(),
        if last.is_empty() {
            ".".to_string()
        } else {
            format!(": {}", last.join(" | "))
        }
    ))
}

#[cfg(test)]
mod tests;
