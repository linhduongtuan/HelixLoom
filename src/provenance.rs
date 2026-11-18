//! What made a result (roadmap 1.3): the inputs by name, size and checksum, the settings, the programs with their versions, the day
//! and where it ran. A result that carries this can say what has changed since it was made, and be made again. Nothing here draws:
//! `ui/provenance.rs` shows a record, and each page makes its own.

use crate::pipeline::flow::digest;
use web_time::{SystemTime, UNIX_EPOCH};

/// The most changes that are listed; the rest are counted.
pub const MOST_CHANGES: usize = 8;
/// How much of a checksum is shown in a window (the text has all of it).
const SHORT: usize = 12;

/// One input of an analysis: its name, how many bytes it holds and a checksum of them (not of the name).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Input {
    pub name: String,
    pub digest: String,
    pub size: usize,
}

impl Input {
    pub fn of(name: &str, bytes: &[u8]) -> Input {
        Input {
            name: name.to_string(),
            digest: digest(bytes),
            size: bytes.len(),
        }
    }

    /// The start of the checksum, for a window.
    pub fn short(&self) -> &str {
        &self.digest[..SHORT.min(self.digest.len())]
    }
}

/// What a page knows when a run starts: enough to make the record when the answer comes.
#[derive(Clone, Debug, PartialEq)]
pub struct Asked {
    pub analysis: &'static str,
    pub inputs: Vec<Input>,
    pub parameters: Vec<(String, String)>,
    pub ran_on: &'static str,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Provenance {
    pub analysis: String,
    pub inputs: Vec<Input>,
    pub parameters: Vec<(String, String)>,
    /// The programs that made the answer, with their versions as they said them; this program last.
    pub programs: Vec<String>,
    pub ran_on: String,
    /// The day it was made, YYYY-MM-DD in UTC.
    pub date: String,
    pub seconds: f64,
}

/// Today, YYYY-MM-DD in UTC.
pub fn today() -> String {
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    crate::connect::utc_date(seconds)
}

/// Text on one line: tabs and line ends become spaces.
fn flat(text: &str) -> String {
    text.chars()
        .map(|c| if c.is_control() { ' ' } else { c })
        .collect()
}

impl Asked {
    /// The record of the answer: what was asked, the programs that answered (this program is added, once), today, and the time.
    pub fn made(self, programs: Vec<String>, seconds: f64) -> Provenance {
        let own = format!("HelixLoom {}", env!("CARGO_PKG_VERSION"));
        let mut programs: Vec<String> = programs.into_iter().filter(|p| *p != own).collect();
        programs.push(own);
        Provenance {
            analysis: self.analysis.to_string(),
            inputs: self.inputs,
            parameters: self.parameters,
            programs,
            ran_on: self.ran_on.to_string(),
            date: today(),
            seconds,
        }
    }
}

impl Provenance {
    /// The record as text, for copying beside a result: a line for each thing, the inputs a line each with the whole checksum.
    pub fn to_text(&self) -> String {
        let settings = if self.parameters.is_empty() {
            "none".to_string()
        } else {
            self.parameters
                .iter()
                .map(|(k, v)| format!("{}={}", flat(k), flat(v)))
                .collect::<Vec<_>>()
                .join("; ")
        };
        let mut text =
            format!(
            "{}\nMade on {} on {} in {:.2} s\nPrograms: {}\nSettings: {settings}\nInputs ({}):\n",
            flat(&self.analysis),
            self.date,
            flat(&self.ran_on),
            self.seconds,
            self.programs.iter().map(|p| flat(p)).collect::<Vec<_>>().join("; "),
            self.inputs.len()
        );
        for i in &self.inputs {
            text.push_str(&format!("  {}\t{}\t{}\n", flat(&i.name), i.size, i.digest));
        }
        text
    }

    /// The record as the program keeps it with a result, to read back itself (`from_kept`): a line for each thing it
    /// holds, the fields of a line with a tab between them. `to_text` is for a person and cannot be read back.
    pub fn kept(&self) -> String {
        use crate::kept::line;
        let mut text = line(&["analysis", &self.analysis]);
        text.push_str(&line(&["ran_on", &self.ran_on]));
        text.push_str(&line(&["date", &self.date]));
        // every digit: the record that comes back is the record that was kept
        text.push_str(&line(&["seconds", &format!("{:?}", self.seconds)]));
        for program in &self.programs {
            text.push_str(&line(&["program", program]));
        }
        for (key, value) in &self.parameters {
            text.push_str(&line(&["parameter", key, value]));
        }
        for input in &self.inputs {
            text.push_str(&line(&[
                "input",
                &input.name,
                &input.size.to_string(),
                &input.digest,
            ]));
        }
        text
    }

    /// The record that `kept` wrote. None for a text that is not one: it has no name of an analysis, or a time or a size
    /// in it is no number (half a record would say things about the result that nobody recorded). A line that this
    /// version does not know is passed over, so that a later version may say more.
    pub fn from_kept(text: &str) -> Option<Provenance> {
        let mut record = Provenance {
            analysis: String::new(),
            inputs: Vec::new(),
            parameters: Vec::new(),
            programs: Vec::new(),
            ran_on: String::new(),
            date: String::new(),
            seconds: 0.0,
        };
        let mut named = false;
        for line in text.lines() {
            let fields = crate::kept::fields(line);
            let fields: Vec<&str> = fields.iter().map(String::as_str).collect();
            match fields.as_slice() {
                ["analysis", analysis] => {
                    record.analysis = analysis.to_string();
                    named = true;
                }
                ["ran_on", ran_on] => record.ran_on = ran_on.to_string(),
                ["date", date] => record.date = date.to_string(),
                ["seconds", seconds] => {
                    record.seconds = seconds
                        .parse()
                        .ok()
                        .filter(|s: &f64| s.is_finite() && *s >= 0.0)?;
                }
                ["program", program] => record.programs.push(program.to_string()),
                ["parameter", key, value] => {
                    record.parameters.push((key.to_string(), value.to_string()));
                }
                ["input", name, size, digest] => record.inputs.push(Input {
                    name: name.to_string(),
                    size: size.parse().ok()?,
                    digest: digest.to_string(),
                }),
                _ => {}
            }
        }
        named.then_some(record)
    }

    /// What differs between the record and the inputs and settings as they are now, a line each; nothing when the result would be
    /// made from the same things. At most `MOST_CHANGES` are listed, and the rest counted.
    pub fn changed_since(&self, inputs: &[Input], parameters: &[(String, String)]) -> Vec<String> {
        let mut changes = Vec::new();
        // inputs that are the same (by name and checksum) pair off; of the rest, one of the same name has changed
        let mut now: Vec<&Input> = inputs.iter().collect();
        let mut unmatched: Vec<&Input> = Vec::new();
        for old in &self.inputs {
            match now.iter().position(|n| *n == old) {
                Some(k) => {
                    now.remove(k);
                }
                None => unmatched.push(old),
            }
        }
        for old in unmatched {
            match now.iter().position(|n| n.name == old.name) {
                Some(k) => {
                    now.remove(k);
                    changes.push(format!("{} has changed", old.name));
                }
                None => changes.push(format!("{} is no longer there", old.name)),
            }
        }
        changes.extend(now.iter().map(|n| format!("{} is new", n.name)));
        if changes.is_empty() && self.inputs.as_slice() != inputs {
            changes.push("the order of the inputs has changed".to_string());
        }
        for (key, was) in &self.parameters {
            match parameters.iter().find(|(k, _)| k == key) {
                Some((_, is)) if is == was => {}
                Some((_, is)) => changes.push(format!("{key} was {was}, is now {is}")),
                None => changes.push(format!("{key} is no longer a setting")),
            }
        }
        for (key, is) in parameters {
            if !self.parameters.iter().any(|(k, _)| k == key) {
                changes.push(format!("{key} is new ({is})"));
            }
        }
        if changes.len() > MOST_CHANGES {
            let more = changes.len() - MOST_CHANGES;
            changes.truncate(MOST_CHANGES);
            changes.push(format!("and {more} more"));
        }
        changes
    }
}

#[cfg(test)]
mod tests;
