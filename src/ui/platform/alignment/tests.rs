use super::*;
use crate::pipeline::steps::Aligned;

/// Twelve proteins of nine columns: the first five are alike, the next four vary, and the last is a gap in two rows.
fn family() -> Aligned {
    let rows = [
        "MKTAYIAK-",
        "MKTAYIAKQ",
        "MKTAYLAKQ",
        "MKTAYIVKQ",
        "MKTAYIAKQ",
        "MKTAYVAKQ",
        "MKTAYIAKQ",
        "MKTAYIGKQ",
        "MKTAYIAKQ",
        "MKTAYIAKQ",
        "MKTAYIAK-",
        "MKTAYIAKQ",
    ];
    Aligned {
        rows: rows
            .iter()
            .enumerate()
            .map(|(k, r)| (format!("row {k}"), r.as_bytes().to_vec()))
            .collect(),
        engine: "test".into(),
        seconds: 0.0,
    }
}

#[test]
fn the_bars_follow_the_chosen_score_and_the_blocks_follow_the_bars() {
    let mut page = AlignmentPage::default();
    page.set_result(family());
    // identity to begin with: a column is conserved when every row is alike (no gap counts for the commonest symbol)
    let (rows, values, score) = page.conservation().expect("an alignment");
    assert_eq!((rows.len(), values.len(), score), (12, 9, Score::Identity));
    assert_eq!(values[0], Some(1.0));
    assert!(values[5].is_some_and(|v| v < 1.0));
    assert_eq!(
        page.stats
            .as_ref()
            .expect("stats")
            .bars
            .blocks
            .iter()
            .map(|b| (b.from, b.to))
            .collect::<Vec<_>>(),
        vec![(0, 4)]
    );
    // a new choice is taken up by the next refresh, and not before
    page.score = Score::Jsd;
    assert_eq!(
        page.conservation().expect("an alignment").2,
        Score::Identity
    );
    page.refresh_scores();
    let (_, jsd, score) = page.conservation().expect("an alignment");
    assert_eq!(score, Score::Jsd);
    assert_eq!(jsd.len(), 9);
    // the columns that are alike score high, the last one (gaps) is lowered, and none is the identity share
    assert!(jsd[0].expect("a score") > jsd[5].expect("a score"));
    assert_ne!(jsd[0], Some(1.0));
    let bars = &page.stats.as_ref().expect("stats").bars;
    assert!(bars.highlight > 0.0 && bars.highlight < 1.0);
    assert!(bars.blocks.iter().all(|b| b.to - b.from + 1 >= MIN_BLOCK));
    // going back gives the identity again
    page.score = Score::Identity;
    page.refresh_scores();
    assert_eq!(
        page.conservation().expect("an alignment").2,
        Score::Identity
    );
}

#[test]
fn the_revision_moves_when_the_scores_do_and_not_before() {
    let mut page = AlignmentPage::default();
    assert_eq!(page.revision(), 0);
    page.set_result(family());
    let first = page.revision();
    assert!(first > 0);
    // nothing has changed
    page.refresh_scores();
    assert_eq!(page.revision(), first);
    // a choice does not count until the bars are made again from it
    page.score = Score::Jsd;
    assert_eq!(page.revision(), first);
    page.refresh_scores();
    let second = page.revision();
    assert!(second > first);
    page.refresh_scores();
    assert_eq!(page.revision(), second);
    // a new alignment
    page.set_result(family());
    assert!(page.revision() > second);
}

#[test]
fn an_alignment_that_arrives_when_the_score_is_chosen_already_is_scored_by_it() {
    let mut page = AlignmentPage {
        score: Score::Jsd,
        ..AlignmentPage::default()
    };
    page.set_result(family());
    assert_eq!(page.conservation().expect("an alignment").2, Score::Jsd);
}

#[test]
fn a_page_with_no_alignment_has_no_conservation() {
    let mut page = AlignmentPage::default();
    assert!(page.conservation().is_none());
    page.refresh_scores();
    assert!(page.conservation().is_none());
}

#[test]
fn a_sequence_with_the_same_name_gets_a_number() {
    let mut page = AlignmentPage::default();
    page.add("a", b"ACGT");
    page.add("a", b"ACGT"); // the very same sequence is not added twice
    page.add("a", b"TTTT");
    page.add("a", b"GGGG");
    let names: Vec<&str> = page.inputs.iter().map(|i| i.name.as_str()).collect();
    assert_eq!(names, ["a", "a #2", "a #3"]);
    assert!(page.inputs.iter().all(|i| i.keep));
}

#[test]
fn pasted_text_is_fasta_or_one_sequence() {
    let fasta = sequences_in(">one\nACGT\nAC\n>two\nGGGG\n", "x");
    assert_eq!(
        fasta,
        vec![
            ("one".into(), b"ACGTAC".to_vec()),
            ("two".into(), b"GGGG".to_vec())
        ]
    );
    let bare = sequences_in("acgt acgt\n123\n", "Pasted");
    assert_eq!(bare, vec![("Pasted".to_string(), b"acgtacgt".to_vec())]);
    assert!(sequences_in("  \n 12 \n", "x").is_empty());
}
