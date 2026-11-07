//! Computing on another machine. `helixloom serve` runs on a workstation: it takes jobs sent over the network (an
//! alignment, a tree estimate, a pipeline, PhyloCNN estimates, a structure prediction), runs each in a process of its
//! own, on a GPU of its own when it needs one, keeps the results on disk and hands them to any computer that asks with
//! the right token. The app on a Mac, a Windows PC or a Linux machine is the interface: it sends the work and shows
//! what comes back (`ui/platform/workstation.rs`).
//!
//! Nothing here is in the original program; it is all new work for the person's workstation.

pub mod cli;
pub mod client;
pub mod gpus;
pub mod http;
pub mod jobs;
pub mod run;
pub mod server;
pub mod settings;
pub mod species;
pub mod tasks;
pub mod weight;
pub mod wire;

/// The port the server listens on unless told otherwise.
pub const DEFAULT_PORT: u16 = 7878;
/// Bumped when the messages change in a way that an older program cannot read.
pub const PROTOCOL: u32 = 1;
/// A token shorter than this is refused: the token is the only thing that keeps other people off the workstation, and the
/// server does not slow down a wrong guess. (`helixloom token` makes one of 64 characters.)
pub const MIN_TOKEN: usize = 32;
/// The longest label or error kept with a job, in characters.
pub const MAX_TEXT: usize = 200;

/// Text on one line: control characters (tabs, line breaks, terminal escapes) become spaces and the ends are trimmed.
/// Nothing is cut short, so that a long name or a long token comes back whole.
pub fn clean_line(text: &str) -> String {
    let flat: String = text
        .chars()
        .map(|c| if c.is_control() { ' ' } else { c })
        .collect();
    flat.trim().to_string()
}

/// The same, no longer than `MAX_TEXT` characters: for a label, an error, or what a server says about itself.
pub fn one_line(text: &str) -> String {
    clean_line(text).chars().take(MAX_TEXT).collect()
}

#[cfg(test)]
mod tests {
    use super::{clean_line, one_line};

    #[test]
    fn text_is_made_one_line_and_short() {
        assert_eq!(one_line("a\tb\nc\r\nd"), "a b c  d");
        assert_eq!(one_line("  padded  "), "padded");
        assert_eq!(one_line(&"x".repeat(500)).chars().count(), 200);
        assert_eq!(one_line("é ü ñ"), "é ü ñ");
    }

    #[test]
    fn a_long_name_or_token_is_cleaned_but_not_cut() {
        assert_eq!(clean_line("a\tb\u{1b}[2J"), "a b [2J");
        let long = "n".repeat(500);
        assert_eq!(clean_line(&format!("  {long}\n")), long);
    }
}
