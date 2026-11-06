//! What Silent Sites and Diagnostic Sites share: the heading of the report, the letters a pattern letter stands for, and
//! the recognition sequence without the `N`s at its ends. The original writes both reports with the same lines
//! (`silent_sites` L35791, `do_diagnostic_site` L36062) and the same string maps (L35946, L36102).

use super::enzymes::Enzyme;

#[cfg(test)]
pub mod golden;

/// The first lines of a result window: the date, the sequence's name, the stretch searched, the report's name, a blank line.
/// The original leaves a space after the name and after the stretch.
pub fn header(date: &str, name: &str, from: usize, to: usize, heading: &str) -> String {
    format!("{date}\n{name} \nFrom {from} to {to}. \n{heading}\n\n")
}

/// Whether sequence letter `c` counts as matching pattern letter `p`. The original's `string map` makes an ambiguity code
/// stand for itself, for every code that is part of it and for the plain bases it covers (`B` is `BYSKCGT`); a plain base
/// only matches itself. Both arguments are upper case.
pub fn pattern_covers(p: u8, c: u8) -> bool {
    let covered: &[u8] = match p {
        b'B' => b"BYSKCGT",
        b'D' => b"DRWKAGT",
        b'H' => b"HMWYACT",
        b'K' => b"KGT",
        b'M' => b"MAC",
        b'N' => b"NBDHVKMRYWSACGT",
        b'R' => b"RAG",
        b'S' => b"SGC",
        b'V' => b"VMSRACG",
        b'W' => b"WAT",
        b'Y' => b"YCT",
        _ => return p == c,
    };
    covered.contains(&c)
}

/// The recognition sequence the way `read_enzymes` keeps it (`flatpat`, L64461): the `N`s at both ends are cut off.
pub fn flat_pattern(e: &Enzyme) -> Vec<u8> {
    // any other letter is dropped first (`regsub -all {[^ACGTNWSYRKMBDHV]}`: a `U` goes), then the `N`s at the ends
    let kept: Vec<u8> = e
        .site
        .iter()
        .copied()
        .filter(|c| b"ACGTNWSYRKMBDHV".contains(c))
        .collect();
    let first = kept.iter().position(|&c| c != b'N');
    let last = kept.iter().rposition(|&c| c != b'N');
    match (first, last) {
        (Some(a), Some(b)) => kept[a..=b].to_vec(),
        _ => Vec::new(),
    }
}

/// The text of one sequence letter in a report: the letter, upper case.
pub fn upper(c: u8) -> char {
    char::from(c.to_ascii_uppercase())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn heading_has_the_originals_spaces_and_blank_line() {
        assert_eq!(
            header("DATE", "pUC19", 1, 300, "Silent Sites"),
            "DATE\npUC19 \nFrom 1 to 300. \nSilent Sites\n\n"
        );
    }

    #[test]
    fn a_code_covers_itself_its_parts_and_the_bases_but_a_base_only_itself() {
        assert!(pattern_covers(b'R', b'A'));
        assert!(pattern_covers(b'R', b'R'));
        assert!(!pattern_covers(b'R', b'C'));
        assert!(pattern_covers(b'B', b'Y'));
        assert!(pattern_covers(b'N', b'N'));
        assert!(!pattern_covers(b'A', b'N'));
        assert!(pattern_covers(b'T', b'T'));
    }

    #[test]
    fn ends_of_n_are_cut_off_but_inner_ones_stay() {
        let e = Enzyme::parse("X", "NNGAATTCNN").expect("a valid pattern");
        assert_eq!(flat_pattern(&e), b"GAATTC".to_vec());
        let sfi = Enzyme::parse("SfiI", "GGCCNNNN^NGGCC").expect("a valid pattern");
        assert_eq!(flat_pattern(&sfi), b"GGCCNNNNNGGCC".to_vec());
    }

    #[test]
    fn a_letter_the_original_does_not_keep_is_dropped_before_the_ends_are_cut() {
        // `read_enzymes` keeps only A C G T N W S Y R K M B D H V of the pattern, so a U never reaches `flatpat`
        let e = Enzyme::parse("X", "NGAUCN").expect("a valid pattern");
        assert_eq!(flat_pattern(&e), b"GAC".to_vec());
    }
}
