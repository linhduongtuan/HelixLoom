use super::*;

/// 240 cases made by running the original's `next_ORF` and `orf_list` on random sequences:
/// text | circular | strand | minimum | start rule | caret | selection start | direction | what it selected.
const GOLDEN: &str = include_str!("golden.txt");

fn settings(circular: bool, strand: Strand) -> Settings {
    Settings {
        circular,
        strand,
        min: 0,
        rule: StartRule::Met,
    }
}

#[test]
fn find_next_selects_what_the_originals_own_code_selects() {
    let mut cases = 0;
    for line in GOLDEN.lines() {
        let f: Vec<&str> = line.split('|').collect();
        let s = Settings {
            circular: f[1] == "1",
            strand: match f[2] {
                "1" => Strand::Top,
                "-1" => Strand::Bottom,
                _ => Strand::Both,
            },
            min: f[3].parse().unwrap(),
            rule: if f[4] == "Met" {
                StartRule::Met
            } else {
                StartRule::AfterStop
            },
        };
        let got = find_next(
            f[0].as_bytes(),
            &s,
            f[5].parse().unwrap(),
            f[6].parse().ok(),
            f[7] == "forward",
        );
        let shown = match got {
            Err(NoStopCodon) => "error".to_string(),
            Ok(None) => "bell".to_string(),
            Ok(Some(found)) => format!("{} {}", found.select.0, found.select.1),
        };
        assert_eq!(shown, f[8], "case {cases}: {line}");
        cases += 1;
    }
    assert_eq!(cases, 240);
}

#[test]
fn the_list_shows_the_originals_habits() {
    let top = settings(false, Strand::Top);
    // two ATGs share one stop; both are listed, in order of start
    let both = orf_list(b"ATGAAAATGAAATAA", &top).unwrap();
    assert_eq!(
        both.iter().map(|e| (e.start, e.end)).collect::<Vec<_>>(),
        [(0, 14), (6, 14)]
    );
    // a circular sequence is read twice round, so an ORF can run across the origin
    let wrap = orf_list(b"AAATAAATGAAA", &settings(true, Strand::Top)).unwrap();
    assert_eq!((wrap[0].start, wrap[0].end), (6, 5));
    // with no stop codon anywhere the original fails
    assert_eq!(orf_list(b"ATGCCCCCC", &top), Err(NoStopCodon));
    // a stop in another frame is taken for the ORF's own: it ends at the end of the sequence
    let stale = orf_list(b"ATGAAAAAA", &top).unwrap();
    assert_eq!((stale[0].start, stale[0].end), (0, 9));
    // a bottom-strand ORF is listed with the top-strand coordinates of its stop and start
    let bottom = orf_list(b"TTATTTCAT", &settings(false, Strand::Bottom)).unwrap();
    assert_eq!(
        bottom,
        [Entry {
            start: 0,
            end: 8,
            reverse: true
        }]
    );
}
