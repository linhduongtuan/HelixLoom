//! A notice outside the window, for a job that took long and ended while the person was looking at something else. It is the operating
//! system's own: on a Mac the notification centre (through `osascript`), on Linux `notify-send` where it is installed; elsewhere, and
//! in a browser, nothing. A notice that cannot be made is not an error: the toast in the window is there too.

/// The AppleScript that shows `body` with `title`: both are put in as text, with the characters that would end it escaped.
pub fn applescript(title: &str, body: &str) -> String {
    let quote = |text: &str| {
        let cleaned: String = text
            .chars()
            .map(|c| if c.is_control() { ' ' } else { c })
            .collect();
        format!("\"{}\"", cleaned.replace('\\', "\\\\").replace('"', "\\\""))
    };
    format!(
        "display notification {} with title {}",
        quote(body),
        quote(title)
    )
}

/// `text` with the characters of markup written as text.
pub fn escape_markup(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Show a notice, and do not wait for it.
#[cfg(not(target_arch = "wasm32"))]
pub fn desktop(title: &str, body: &str) {
    use std::process::{Command, Stdio};
    let mut command = if cfg!(target_os = "macos") {
        // by its place: a program that is first on the PATH is not to be trusted with what it is told
        let mut c = Command::new("/usr/bin/osascript");
        c.arg("-e").arg(applescript(title, body));
        c
    } else if cfg!(target_os = "linux") {
        let mut c = Command::new("notify-send");
        // many notification daemons read the body as markup
        c.arg("--app-name=HelixLoom")
            .arg("--")
            .arg(escape_markup(title))
            .arg(escape_markup(body));
        c
    } else {
        return;
    };
    command
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null());
    // a thread waits for it, so that a finished program is not left as a zombie; no program that is missing matters
    std::thread::spawn(move || {
        let _ = command.status();
    });
}

#[cfg(target_arch = "wasm32")]
pub fn desktop(_title: &str, _body: &str) {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_script_holds_the_text_as_text_whatever_it_contains() {
        assert_eq!(
            applescript("HelixLoom", "Phylogeny finished in 42 s"),
            "display notification \"Phylogeny finished in 42 s\" with title \"HelixLoom\""
        );
        // quotes and backslashes cannot end the string and run what follows
        let script = applescript("a\"b", "say \"x\" & do shell script \"rm -rf ~\"\\");
        assert_eq!(
            script,
            "display notification \"say \\\"x\\\" & do shell script \\\"rm -rf ~\\\"\\\\\" with title \"a\\\"b\""
        );
        assert_eq!(escape_markup("a & <b>"), "a &amp; &lt;b&gt;");
        // line breaks cannot start another statement
        assert!(!applescript("t", "one\ntwo\rthree").contains(['\n', '\r']));
    }
}
