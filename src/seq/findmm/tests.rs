use super::*;

/// 300 cases made by running the original's `do_find_mm` on random sequences and patterns, tab separated:
/// text, circular, pattern, case sensitive, rev-com, as literal, mismatches, minimum 3' bases, wrap, action,
/// caret, selection start, and what it did.
const GOLDEN: &str = include_str!("golden.txt");

fn shown(outcome: &Outcome) -> String {
    match outcome {
        Outcome::Selected {
            start,
            end,
            message,
        } => format!("{start} {end};{message}"),
        Outcome::NotFound { message } => format!(";{message}"),
        Outcome::Highlighted { tags, message } => {
            let parts: Vec<String> = tags
                .iter()
                .map(|t| {
                    format!(
                        "{}:{}-{}",
                        if t.reverse { "foundr" } else { "foundf" },
                        t.start,
                        t.end
                    )
                })
                .collect();
            format!("{};{message}", parts.join(","))
        }
    }
}

#[test]
fn every_search_does_what_the_originals_own_code_does() {
    let mut cases = 0;
    for line in GOLDEN.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let q = Query {
            pattern: f[2].to_string(),
            case_sensitive: f[3] == "1",
            rev_com: f[4] == "1",
            literal: f[5] == "1",
            mismatches: f[6].parse().unwrap(),
            min_3_bases: f[7].parse().ok(),
            wrap: f[8] == "1",
        };
        let action = match f[9] {
            "next" => Action::Next,
            "prev" => Action::Prev,
            _ => Action::Highlight,
        };
        let got = find(
            f[0].as_bytes(),
            f[1] == "1",
            &q,
            action,
            f[10].parse().unwrap(),
            f[11].parse().ok(),
        );
        assert_eq!(shown(&got), f[12], "case {cases}: {line}");
        cases += 1;
    }
    assert_eq!(cases, 300);
}

#[test]
fn the_defaults_are_the_originals() {
    let q = Query::default();
    assert!(!q.case_sensitive && q.rev_com && !q.literal && q.wrap);
    assert_eq!((q.mismatches, q.min_3_bases), (0, None));
}

#[test]
fn several_patterns_and_the_reverse_complement_are_tried_in_order() {
    let q = Query {
        pattern: "GGATCC|GAATTC".into(),
        ..Query::default()
    };
    // forward GAATTC at 9, and GGATCC's reverse complement (itself) at 0: Find Next from 1 reaches 9
    let text = b"GGATCCAAAGAATTCAAA";
    match find(text, false, &q, Action::Next, 1, None) {
        Outcome::Selected { start, end, .. } => assert_eq!((start, end), (9, 15)),
        other => panic!("{other:?}"),
    }
    // Find Prev skips a match that starts one base before the caret when the sequence is not circular
    let aa = Query {
        pattern: "AA".into(),
        rev_com: false,
        ..Query::default()
    };
    match find(b"CAAAC", false, &aa, Action::Prev, 3, None) {
        Outcome::Selected { start, .. } => assert_eq!(start, 1),
        other => panic!("{other:?}"),
    }
}
