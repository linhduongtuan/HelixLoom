//! What the app and the workstation send each other: a container of named sections. A section is a line
//! `@name length`, then exactly `length` bytes, then a line break. Sequences, trees, structures and tables are text or
//! bytes of any kind, so nothing is escaped or guessed: a reader knows how long a section is before it reads it, and a
//! container that does not add up (names that are not plain, a length that runs past the end, bytes left over) is
//! refused.

use std::collections::HashMap;

/// A container holds at most this many sections.
pub const MAX_SECTIONS: usize = 64;
const MAX_NAME: usize = 32;
/// The line that starts a section is never longer than this.
const MAX_HEAD: usize = 64;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Container {
    sections: Vec<(String, Vec<u8>)>,
}

fn valid_name(name: &str) -> bool {
    (1..=MAX_NAME).contains(&name.len())
        && name
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'_')
}

impl Container {
    pub fn new() -> Container {
        Container::default()
    }

    /// This container with section `name` added (`name` is lower case letters, digits and `_`, at most 32).
    pub fn with(mut self, name: &str, data: impl Into<Vec<u8>>) -> Container {
        debug_assert!(valid_name(name), "{name}");
        self.sections.push((name.to_string(), data.into()));
        self
    }

    pub fn bytes(&self, name: &str) -> Option<&[u8]> {
        self.sections
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, d)| d.as_slice())
    }

    /// The section as text; `None` when it is missing or is not UTF-8.
    pub fn text(&self, name: &str) -> Option<&str> {
        std::str::from_utf8(self.bytes(name)?).ok()
    }

    /// The section as text, or the words that say what is wrong.
    pub fn require(&self, name: &str) -> Result<&str, String> {
        match self.bytes(name) {
            None => Err(format!("The message has no \"{name}\" section.")),
            Some(b) => {
                std::str::from_utf8(b).map_err(|_| format!("The \"{name}\" section is not text."))
            }
        }
    }

    pub fn encode(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.sections.iter().map(|(_, d)| d.len() + 48).sum());
        for (name, data) in &self.sections {
            out.extend(format!("@{name} {}\n", data.len()).as_bytes());
            out.extend(data);
            out.push(b'\n');
        }
        out
    }

    pub fn decode(bytes: &[u8]) -> Result<Container, String> {
        const BAD: &str = "The message is not a container of sections";
        let mut rest = bytes;
        let mut out = Container::new();
        while !rest.is_empty() {
            if out.sections.len() >= MAX_SECTIONS {
                return Err(format!("{BAD}: more than {MAX_SECTIONS} sections."));
            }
            let head_end = rest
                .iter()
                .take(MAX_HEAD)
                .position(|&b| b == b'\n')
                .ok_or_else(|| format!("{BAD}: a section starts with no header line."))?;
            let head = std::str::from_utf8(&rest[..head_end])
                .map_err(|_| format!("{BAD}: a header line is not text."))?;
            let (name, length) = head
                .strip_prefix('@')
                .and_then(|h| h.split_once(' '))
                .ok_or_else(|| format!("{BAD}: a header line is not \"@name length\"."))?;
            if !valid_name(name) {
                return Err(format!("{BAD}: \"{name}\" is not a section name."));
            }
            if out.sections.iter().any(|(n, _)| n == name) {
                return Err(format!("{BAD}: section \"{name}\" appears twice."));
            }
            let length: usize = length
                .parse()
                .map_err(|_| format!("{BAD}: the length of \"{name}\" is not a number."))?;
            let start = head_end + 1;
            // the data, then its line break: both must be there
            let stop = start
                .checked_add(length)
                .filter(|&s| s < rest.len())
                .ok_or_else(|| format!("{BAD}: section \"{name}\" runs past the end."))?;
            if rest[stop] != b'\n' {
                return Err(format!(
                    "{BAD}: section \"{name}\" is not followed by a line break."
                ));
            }
            out.sections
                .push((name.to_string(), rest[start..stop].to_vec()));
            rest = &rest[stop + 1..];
        }
        Ok(out)
    }
}

/// The lines `key=value` of a text, as a map (a line with no `=` is ignored; the first `=` splits).
pub fn meta(text: &str) -> HashMap<&str, &str> {
    text.lines().filter_map(|l| l.split_once('=')).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_container_comes_back_as_it_went_in_whatever_its_bytes() {
        let binary: Vec<u8> = (0..=255u8).collect();
        let c = Container::new()
            .with("meta", "engine=mafft\nbootstrap=100")
            .with("data", binary.clone())
            .with("empty", Vec::new())
            .with("fasta", ">a\nACGT\n>b\nAC-T\n");
        let bytes = c.encode();
        let back = Container::decode(&bytes).expect("a container");
        assert_eq!(back, c);
        assert_eq!(back.bytes("data"), Some(binary.as_slice()));
        assert_eq!(back.text("empty"), Some(""));
        assert_eq!(back.text("missing"), None);
        // a section that holds the container's own syntax is just bytes
        let tricky = Container::new().with("a", "@b 3\nxyz\n").with("c", "@@");
        assert_eq!(
            Container::decode(&tricky.encode()).expect("a container"),
            tricky
        );
        assert_eq!(Container::decode(b"").expect("empty"), Container::new());
    }

    #[test]
    fn sections_that_do_not_add_up_are_refused() {
        let good = Container::new()
            .with("a", "hello")
            .with("b", "world")
            .encode();
        assert!(Container::decode(&good).is_ok());
        // cut anywhere inside a section: refused, never a panic (the end of the first section is a whole container)
        let first = Container::new().with("a", "hello").encode().len();
        for k in (1..good.len()).filter(|&k| k != first) {
            let cut = Container::decode(&good[..k]);
            assert!(cut.is_err(), "cut at {k} read as {cut:?}");
        }
        assert_eq!(
            Container::decode(&good[..first])
                .expect("one section")
                .text("a"),
            Some("hello")
        );
        for bad in [
            &b"a 3\nabc\n"[..],                   // no @
            b"@a\nabc\n",                         // no length
            b"@A 3\nabc\n",                       // upper case name
            b"@a-b 3\nabc\n",                     // character not allowed
            b"@ 3\nabc\n",                        // no name
            b"@a x\nabc\n",                       // length is not a number
            b"@a -1\nabc\n",                      // negative length
            b"@a 99999999999999999999999\nabc\n", // length overflows
            b"@a 3\nabcd\n",                      // data longer than its length
            b"@a 3\nabc",                         // no line break after the data
            b"@a 1\nx\n@a 1\ny\n",                // twice
            b"@a 1\nx\njunk",                     // bytes left over
        ] {
            assert!(
                Container::decode(bad).is_err(),
                "{:?}",
                String::from_utf8_lossy(bad)
            );
        }
        // a name that is too long, a header line that never ends, too many sections
        assert!(Container::decode(format!("@{} 0\n\n", "a".repeat(33)).as_bytes()).is_err());
        assert!(Container::decode(&vec![b'@'; 1000]).is_err());
        let many: Vec<u8> = (0..=MAX_SECTIONS)
            .flat_map(|i| format!("@s{i} 0\n\n").into_bytes())
            .collect();
        assert!(Container::decode(&many).is_err());
    }

    #[test]
    fn a_claimed_length_far_beyond_the_bytes_does_not_allocate_it() {
        let huge = format!("@a {}\nshort\n", usize::MAX / 2);
        assert!(Container::decode(huge.as_bytes()).is_err());
    }

    #[test]
    fn required_text_says_what_is_missing_or_wrong() {
        let c = Container::new()
            .with("t", "text")
            .with("b", vec![0xff, 0xfe]);
        assert_eq!(c.require("t"), Ok("text"));
        assert!(c.require("nope").unwrap_err().contains("no \"nope\""));
        assert!(c.require("b").unwrap_err().contains("not text"));
    }

    #[test]
    fn key_value_lines_make_a_map() {
        let m = meta("a=1\nb=two=2\nno equals\n\nc=");
        assert_eq!(m["a"], "1");
        assert_eq!(m["b"], "two=2");
        assert_eq!(m["c"], "");
        assert_eq!(m.len(), 3);
    }
}
