//! Text that the program keeps and reads back itself (the session's list, a result of a page, the record of what made
//! it): lines, each a few fields with a tab between them. A field may hold anything the person wrote, a tab or a line end
//! too, so those are written with a backslash and read back as they were.

/// A field as it is written on a line: no tab and no line end is left in it.
pub fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            c => out.push(c),
        }
    }
    out
}

/// A field as it was before it was written.
pub fn unescape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('r') => out.push('\r'),
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

/// A line of `fields`, ended.
pub fn line(fields: &[&str]) -> String {
    let written: Vec<String> = fields.iter().map(|field| escape(field)).collect();
    format!("{}\n", written.join("\t"))
}

/// The fields of a line, as they were before they were written.
pub fn fields(line: &str) -> Vec<String> {
    line.split('\t').map(unescape).collect()
}

/// A checksum of `text` (FNV-1a, 64 bits): to tell that a text that was read back is the one that was written, not to keep anything
/// secret. It reads each byte once, so a text of megabytes costs a moment, and is made only for a text that is written or read.
pub fn hash(text: &str) -> u64 {
    text.bytes().fold(0xcbf2_9ce4_8422_2325, |sum, byte| {
        (sum ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// The word that `table` writes `value` with (none, for a value that the table does not hold).
pub fn word<T: PartialEq + Copy>(table: &[(T, &'static str)], value: T) -> &'static str {
    table
        .iter()
        .find(|(held, _)| *held == value)
        .map_or("", |(_, word)| word)
}

/// The value that `word` stands for in `table`.
pub fn value<T: Copy>(table: &[(T, &'static str)], word: &str) -> Option<T> {
    table
        .iter()
        .find(|(_, held)| *held == word)
        .map(|(value, _)| *value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_checksum_is_the_one_of_the_published_fnv_1a_and_tells_texts_apart() {
        // the published test values of FNV-1a (64 bits)
        assert_eq!(hash(""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(hash("a"), 0xaf63_dc4c_8601_ec8c);
        assert_eq!(hash("foobar"), 0x8594_4171_f739_67e8);
        // the same length and letters in another order, and one letter changed in a long text
        assert_ne!(hash("ACGT"), hash("TGCA"));
        let long = "ACGT".repeat(10_000);
        let changed = format!("{}C{}", &long[..20_000], &long[20_001..]);
        assert_ne!(hash(&long), hash(&changed));
    }

    #[test]
    fn a_line_of_fields_is_read_back_as_it_was_written_whatever_the_fields_hold() {
        for written in [
            vec!["plain", "fields"],
            vec![
                "a\ttab",
                "a\nline end",
                "a\r\nWindows one",
                "a \\ backslash",
                "\\t that is no tab",
            ],
            vec!["", "", ""],
            vec!["one"],
            vec!["ends with a backslash\\"],
        ] {
            let kept = line(&written);
            assert!(
                kept.ends_with('\n') && !kept[..kept.len() - 1].contains(['\n', '\r']),
                "{kept:?}"
            );
            assert_eq!(kept.matches('\t').count(), written.len() - 1, "{kept:?}");
            assert_eq!(fields(kept.trim_end_matches('\n')), written);
        }
        // a line that ends in half an escape (cut short) keeps its backslash
        assert_eq!(unescape("cut\\"), "cut\\");
    }
}
