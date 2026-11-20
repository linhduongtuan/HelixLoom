//! Projects: the work on the pages as one file that can be kept, sent to another person and opened on another computer or in a
//! browser. **File ▸ Save Project As…** writes the open sequences (each as GenBank, with its features) and what the Alignment, Trees,
//! Structure and Pipelines pages show, in the text those pages keep for the session (`ui::results`); **File ▸ Open Project…** (or
//! dropping a `.hlproj` file on the window, or naming it on the command line) puts them back beside what is open.
//!
//! A project is made of the same lines as the session (`kept.rs`): one line to a sequence and one to a result, a field that holds a
//! tab or a line end written with a backslash, and a last line that says how many of each there are and the checksum of all the
//! lines before it, so that a file that was cut short or changed is refused and not half opened. Unlike the session it holds every
//! text itself: a sequence that is also in a file of the computer is not left to be read from there, which another computer cannot.
//! The sequences come back tied to no file and with unsaved changes, so that saving one asks where, and closing it asks whether.
//!
//! Not in a project: what the tool windows show, how a page is scrolled or zoomed, and the places of the windows.

use crate::kept::{fields, hash, line};
use crate::seq::io;
use crate::ui::app::App;
use crate::ui::platform::route::Source;
use crate::ui::session::{SavedDoc, SavedResult, Session};

/// What a project file's name ends in.
pub const EXTENSION: &str = "hlproj";
/// The name a project is offered under.
pub const FILE: &str = "project.hlproj";
const NAME: &str = "helixloom-project";
const VERSION: &str = "1";
/// The most that is read of a project file. The work of one person is far under it, and a browser's page far under what it would
/// take to hold a file of this size and its texts at once.
const MOST_BYTES: usize = 256 * 1024 * 1024;

/// The checksum as the end line writes it.
fn seal_of(body: &str) -> String {
    format!("{:016x}", hash(body))
}

/// A body (all the lines before the end) with its end line.
fn seal(body: &str, docs: usize, results: usize) -> String {
    format!(
        "{body}{}",
        line(&[
            "end",
            &docs.to_string(),
            &results.to_string(),
            &seal_of(body)
        ])
    )
}

/// The text of a project made of `session`, whose sequences and results must all hold their texts.
pub fn write(session: &Session) -> Result<String, String> {
    if session.docs.is_empty() && session.results.is_empty() {
        return Err("Nothing to save in a project yet.".into());
    }
    let mut body = line(&[NAME, VERSION]);
    body += &line(&["workspace", &session.workspace]);
    body += &line(&["active", &session.active.to_string()]);
    for doc in &session.docs {
        let text = doc
            .text
            .as_deref()
            .ok_or_else(|| format!("The text of {} is not there to keep.", doc.name))?;
        body += &line(&["doc", &doc.name, text]);
    }
    for result in &session.results {
        let text = result
            .text
            .as_deref()
            .ok_or_else(|| format!("The text of the {} is not there to keep.", result.kind))?;
        body += &line(&["result", &result.kind, text]);
    }
    Ok(seal(&body, session.docs.len(), session.results.len()))
}

/// The session that a project holds; why not, when it is no project, is of a later version, was cut short or changed, or holds nothing.
pub fn read(text: &str) -> Result<Session, String> {
    read_at_most(text, MOST_BYTES)
}

fn read_at_most(text: &str, most: usize) -> Result<Session, String> {
    if text.len() > most {
        return Err(format!(
            "This file is too large for a project ({} megabytes at the most).",
            most / (1024 * 1024)
        ));
    }
    let first = text.lines().next().map(fields).unwrap_or_default();
    if first.first().map(String::as_str) != Some(NAME) {
        return Err("This is not a HelixLoom project file.".into());
    }
    if first.get(1).map(String::as_str) != Some(VERSION) {
        return Err(
            "This project was made by a later version of HelixLoom than this one: it cannot be read here.".into(),
        );
    }
    // the end line is the last one: what is before it is what its checksum is of
    let end_at = text.rfind("\nend\t").map(|at| at + 1);
    let Some(end_at) = end_at else {
        return Err("This project file was cut short: its last line is not there.".into());
    };
    let (body, end) = text.split_at(end_at);
    let end = fields(end.trim_end_matches(['\n', '\r']));
    let [_, docs, results, sum] = end.as_slice() else {
        return Err("This project file was cut short: its last line is not whole.".into());
    };
    if *sum != seal_of(body) {
        return Err("This project file has changed since it was written (its checksum is not the one it carries).".into());
    }
    let mut session = Session::default();
    for line in body.lines() {
        let f = fields(line);
        match f.iter().map(String::as_str).collect::<Vec<_>>().as_slice() {
            ["workspace", workspace] => session.workspace = workspace.to_string(),
            ["active", n] => session.active = n.parse().unwrap_or(0),
            ["doc", name, text] => session.docs.push(SavedDoc {
                name: name.to_string(),
                path: None,
                dirty: true,
                text: Some(text.to_string()),
                keep: None,
            }),
            ["result", kind, text] => session.results.push(SavedResult {
                kind: kind.to_string(),
                text: Some(text.to_string()),
                keep: None,
            }),
            // the name line, and a line that a later version adds
            _ => {}
        }
    }
    if (
        session.docs.len().to_string(),
        session.results.len().to_string(),
    ) != (docs.clone(), results.clone())
    {
        return Err(
            "This project file was cut short: it holds fewer sequences or results than it says."
                .into(),
        );
    }
    if session.docs.is_empty() && session.results.is_empty() {
        return Err("This project file holds nothing.".into());
    }
    Ok(session)
}

impl App {
    /// The work as a session that holds every text: the sequences that are worth keeping (not the example plasmid, not an empty new
    /// one), and what the pages show.
    fn project_session(&self) -> Session {
        let kept = self.kept_docs();
        let front = match self.tab {
            crate::ui::app::Tab::Doc(i) => i,
            crate::ui::app::Tab::Assembly => self.last_doc,
        };
        let docs = kept
            .iter()
            .map(|(_, view)| SavedDoc {
                name: view.doc.name.clone(),
                path: None,
                dirty: true,
                text: Some(io::to_genbank(&view.doc)),
                keep: None,
            })
            .collect();
        let results = self
            .result_stamps()
            .into_iter()
            .filter_map(|(kind, _)| {
                Some(SavedResult {
                    kind: kind.to_string(),
                    text: Some(self.result_text(kind)?),
                    keep: None,
                })
            })
            .collect();
        Session {
            workspace: self.workspace.view_name().to_string(),
            active: kept.iter().position(|(i, _)| *i == front).unwrap_or(0),
            docs,
            lost: 0,
            results,
        }
    }

    /// File ▸ Save Project As…: the work as one file, to be saved where the person says (downloaded, in a browser).
    pub(in crate::ui) fn save_project(&mut self) {
        let session = self.project_session();
        match write(&session) {
            Ok(text) => {
                let said = crate::ui::figure::keep(FILE, text.as_bytes());
                if !said.is_empty() {
                    self.say(said);
                }
            }
            Err(why) => self.say(format!(
                "{why} Open a sequence, or make an alignment, a tree or a structure, and try again."
            )),
        }
    }

    /// Put a project back beside what is open (see `restore_session`); why not, when it cannot be read, and then nothing is changed.
    pub fn open_project(&mut self, file: Source) {
        let shown = file.shown();
        let read_it = file
            .bytes()
            .map_err(|e| e.to_string())
            .and_then(|bytes| {
                if bytes.len() > MOST_BYTES {
                    return Err("This file is too large for a project.".to_string());
                }
                String::from_utf8(bytes.into_owned())
                    .map_err(|_| "This is not a HelixLoom project file.".to_string())
            })
            .and_then(|text| read(&text));
        let session = match read_it {
            Ok(session) => session,
            Err(why) => {
                self.say(format!("Could not open {shown}: {why}"));
                return;
            }
        };
        let done = self.restore_session(&session);
        let results = crate::ui::results::named(done.results.iter().map(String::as_str));
        let mut said = format!(
            "Opened the project {shown}: {} sequence{}",
            done.opened,
            if done.opened == 1 { "" } else { "s" }
        );
        if !results.is_empty() {
            said += &format!(", with {results}");
        }
        if done.lost > 0 || !done.results_lost.is_empty() {
            said += &format!(
                "; {} sequence{} and {} result{} could not be put back",
                done.lost,
                if done.lost == 1 { "" } else { "s" },
                done.results_lost.len(),
                if done.results_lost.len() == 1 {
                    ""
                } else {
                    "s"
                }
            );
        }
        self.say(said);
    }
}

#[cfg(test)]
mod tests;
