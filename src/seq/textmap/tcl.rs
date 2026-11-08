//! The few string commands of Tcl that the original's Text Map is written with, as they behave in Tcl 8.5 (the version the
//! original runs on), on rows of characters: a count that is not positive repeats nothing, a range or a replacement that falls
//! outside the string is cut to it, and one that does not fall inside it at all is empty or changes nothing.

/// `string repeat s n`: nothing for an `n` that is not positive.
pub(super) fn repeat(s: &str, n: i64) -> Vec<char> {
    let count = usize::try_from(n).unwrap_or(0);
    s.chars().cycle().take(s.chars().count() * count).collect()
}

/// `string repeat " " n`.
pub(super) fn spaces(n: i64) -> Vec<char> {
    repeat(" ", n)
}

/// `string range s first last`, with `last` counted from the first character (`end-N` is `len - 1 - N`).
pub(super) fn range(s: &[char], first: i64, last: i64) -> Vec<char> {
    let len = s.len() as i64;
    let first = first.max(0);
    let last = last.min(len - 1);
    if first > last {
        return Vec::new();
    }
    s[first as usize..=last as usize].to_vec()
}

/// `string replace s first last new`.
pub(super) fn replace(s: &[char], first: i64, last: i64, new: &[char]) -> Vec<char> {
    let len = s.len() as i64;
    let first = first.max(0);
    let last = last.min(len - 1);
    if first > last || first >= len {
        return s.to_vec();
    }
    let mut out = s[..first as usize].to_vec();
    out.extend_from_slice(new);
    out.extend_from_slice(&s[last as usize + 1..]);
    out
}

/// A row of characters from text.
pub(super) fn chars(s: &str) -> Vec<char> {
    s.chars().collect()
}

pub(super) fn text(s: &[char]) -> String {
    s.iter().collect()
}

/// Tcl's `%`: the result has the sign of the divisor.
pub(super) fn modulo(a: i64, b: i64) -> i64 {
    a.rem_euclid(b)
}

/// Tcl's `/` on integers: the quotient is rounded down.
pub(super) fn divide(a: i64, b: i64) -> i64 {
    a.div_euclid(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(s: &str) -> Vec<char> {
        chars(s)
    }

    // each line is what `tclsh` 8.5.9 answers for the same command
    #[test]
    fn the_string_commands_behave_as_tcl_8_5_does() {
        assert!(repeat("x", -3).is_empty() && repeat("x", 0).is_empty());
        assert_eq!(text(&repeat("x", 2)), "xx");
        assert_eq!(text(&range(&c("abc"), 1, 3 - 1)), "bc");
        assert!(range(&c("abc"), 5, 2).is_empty());
        assert!(range(&c("abc"), 0, 3 - 1 - 5).is_empty());
        assert_eq!(text(&range(&c("abc"), -2, 1)), "ab");
        assert!(range(&c("abc"), 2, 1).is_empty());
        let r = |first, last| text(&replace(&c("abcdef"), first, last, &c("XY")));
        assert_eq!(r(1, 2), "aXYdef");
        assert_eq!(r(4, 9), "abcdXY");
        assert_eq!(r(7, 9), "abcdef");
        assert_eq!(r(3, 1), "abcdef");
        assert_eq!(r(-2, 1), "XYcdef");
        assert_eq!(r(2, -1), "abcdef");
        assert_eq!(text(&replace(&c("abc"), 3, 3, &c("Z"))), "abc");
        assert_eq!((modulo(-1, 5), modulo(7, 3)), (4, 1));
        assert_eq!((divide(-1, 3), divide(-7, 2)), (-1, -4));
    }
}
