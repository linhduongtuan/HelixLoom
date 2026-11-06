//! The original's Jump To and Select From-To dialogs (`jump_to`, AppMain.tcl line 13169, and
//! `select_from_to`, line 13097): the rules that turn what was typed into a caret position or a selection.
//! Tested against the original's own lines run with Tcl.

/// What a number box holds once its leading zeros are dropped (`string trimleft … 0`): a box with nothing
/// or only zeros in it counts as empty.
fn number(text: &str) -> Option<u64> {
    let digits = text.trim_start_matches('0');
    (!digits.is_empty()).then(|| digits.parse().unwrap_or(u64::MAX))
}

/// Jump To: the caret goes in front of base `n`, which is position `n - 1` counting from 0. An empty box, a
/// zero or a 1 all go to the very start.
pub fn jump_to(text: &str) -> usize {
    number(text).map_or(0, |n| usize::try_from(n - 1).unwrap_or(usize::MAX))
}

/// Select From-To: bases `from` to `to`, counted from 1 and both included, as the 0-based half-open range to
/// select. Numbers the wrong way round are swapped (only when both are filled in), an empty `from` means 1 and
/// an empty `to` the end of the `len` bases. A range past the end is left for the caller to clip.
pub fn select_from_to(from: &str, to: &str, len: usize) -> (usize, usize) {
    let (mut from, mut to) = (number(from), number(to));
    if let (Some(a), Some(b)) = (from, to) {
        if a > b {
            (from, to) = (Some(b), Some(a));
        }
    }
    let first = from.unwrap_or(1) - 1;
    let last = to.unwrap_or(len as u64);
    (
        usize::try_from(first).unwrap_or(usize::MAX),
        usize::try_from(last).unwrap_or(usize::MAX),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Lines from the original's own code: kind | from (or the jump box) | to | length | what it selected.
    const GOLDEN: &str = "\
range|5|10|100|4 10
range|10|5|100|4 10
range|||100|0 100
range|0|20|100|0 20
range|5||100|4 100
range||7|100|0 7
range|007|010|50|6 10
range|20|0|100|19 100
range|300|20|100|19 300
range|1|1|100|0 1
range|100|100|100|99 100
range|9|10|100|8 10
range|10|9|100|8 10
range|00|00|100|0 100
range|1|100|100|0 100
range|42||42|41 42
range|0|0|100|0 100
range||1|100|0 1
range|2||1|1 1
jump|5||100|4 4
jump|1||100|0 0
jump|0||100|0 0
jump|||100|0 0
jump|007||100|6 6
jump|100||100|99 99
jump|250||100|249 249
jump|10||100|9 9
jump|00||100|0 0
";

    #[test]
    fn both_dialogs_do_what_the_originals_own_code_does() {
        let mut checked = 0;
        for line in GOLDEN.lines() {
            let f: Vec<&str> = line.split('|').collect();
            let want: Vec<usize> = f[4].split(' ').map(|n| n.parse().unwrap()).collect();
            let got = match f[0] {
                "range" => {
                    let (a, b) = select_from_to(f[1], f[2], f[3].parse().unwrap());
                    vec![a, b]
                }
                _ => {
                    let at = jump_to(f[1]);
                    vec![at, at]
                }
            };
            assert_eq!(got, want, "{line}");
            checked += 1;
        }
        assert_eq!(checked, 28);
    }
}
