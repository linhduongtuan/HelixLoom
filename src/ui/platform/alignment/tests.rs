use super::*;
use crate::pipeline::steps::Aligned;
use crate::ui::app::{App, Startup};
use std::time::{Duration, Instant};

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
        fasta.expect("fasta"),
        vec![
            ("one".into(), b"ACGTAC".to_vec()),
            ("two".into(), b"GGGG".to_vec())
        ]
    );
    let bare = sequences_in("acgt acgt\n123\n", "Pasted");
    assert_eq!(
        bare.expect("bare"),
        vec![("Pasted".to_string(), b"acgtacgt".to_vec())]
    );
    assert!(sequences_in("  \n 12 \n", "x").expect("nothing").is_empty());
}

#[test]
fn a_clustal_text_gives_its_sequences_without_the_gaps_and_one_that_cannot_be_read_gives_its_reason(
) {
    let text = "CLUSTAL W (1.83) multiple sequence alignment\n\nfirst  AC-GT-A\nsecond ACAG.TA\n\n";
    let found = sequences_in(text, "x").expect("clustal");
    assert_eq!(
        found,
        vec![
            ("first".to_string(), b"ACGTA".to_vec()),
            ("second".to_string(), b"ACAGTA".to_vec())
        ]
    );
    // a page takes them as it takes FASTA
    let mut page = AlignmentPage::default();
    found.iter().for_each(|(n, b)| page.add(n, b));
    assert_eq!(page.inputs.len(), 2);
    // rows of different lengths are no alignment: the reason, and nothing is added
    let ragged = "CLUSTAL W\n\nx ACGT\ny AC\n";
    let why = sequences_in(ragged, "x").expect_err("ragged");
    assert!(why.contains("not of one length"), "{why}");
    // a pasted protein that spells a header is a protein
    let protein = sequences_in("MUSCLEMUSCLE", "Pasted").expect("a sequence");
    assert_eq!(protein.len(), 1);
}

#[test]
fn the_progressive_method_aligns_with_the_settings_the_page_shows_and_the_record_names_them() {
    let mut app = App::new(Startup::default());
    let ctx = egui::Context::default();
    let options = clustal::Options {
        guide: clustal::Guide::Upgma,
        free_end_gaps: false,
        refine: 1,
        fast: false,
    };
    let page = &mut app.platform.alignment;
    for d in crate::sample::tree_demo() {
        page.add(&d.name, &d.seq);
    }
    // the engine that was chosen carries other settings than the ones the page shows now: the page's are the ones that count
    page.engine = Engine::Progressive(clustal::Options::default());
    page.progressive = options;
    app.start_alignment(&ctx);
    let end = Instant::now() + Duration::from_secs(60);
    while app.platform.alignment.result.is_none() {
        app.alignment_tick(&ctx);
        assert!(Instant::now() < end, "gave up waiting");
        std::thread::sleep(Duration::from_millis(15));
    }
    let page = &app.platform.alignment;
    let result = page.result.as_ref().expect("an alignment");
    for part in ["UPGMA guide tree", "end gaps cost", "refined, up to 1 pass"] {
        assert!(result.engine.contains(part), "{}", result.engine);
    }
    let (record, engine) = page.made.as_ref().expect("the record");
    assert!(*engine == Engine::Progressive(options));
    let text = record.to_text();
    for part in ["progressive (ClustalW-style)", "UPGMA", "refinement passes"] {
        assert!(text.contains(part), "{text}");
    }
    // and kept, then put back on another page, it can be run again with them
    let kept = page.kept();
    let mut back = AlignmentPage::default();
    back.put_back(&kept).expect("put back");
    assert_eq!(back.progressive, options);
    assert!(back
        .made
        .as_ref()
        .is_some_and(|(_, e)| *e == Engine::Progressive(options)));
}

#[test]
fn a_file_named_as_clustal_that_is_not_one_is_said_so_and_adds_nothing() {
    let folder = std::env::temp_dir().join(format!("helixloom-aln-{}", std::process::id()));
    std::fs::create_dir_all(&folder).expect("a folder");
    let path = folder.join("x.aln");
    std::fs::write(&path, ">a\nACGT\n>b\nACGA\n").expect("a file");
    let mut page = AlignmentPage::default();
    page.add_source(Source::Path(&path));
    assert!(page.inputs.is_empty(), "no one long junk sequence");
    assert!(page.message.contains("Clustal"), "{}", page.message);
    // and a real one is read, whatever the case of its extension
    let real = folder.join("y.ALN");
    std::fs::write(&real, "CLUSTAL W\n\nu  AC-GT\nv  ACAGT\n").expect("a file");
    page.add_source(Source::Path(&real));
    assert_eq!(page.inputs.len(), 2);
    std::fs::remove_dir_all(&folder).ok();
}
