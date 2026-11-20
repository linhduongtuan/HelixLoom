use super::*;
use crate::pipeline::steps::Aligned;
use crate::seq::codonalign::Codons;
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

/// One codon for each amino acid (a second one for alanine when `alt`), as the letters a coding sequence is written with.
fn coding(protein: &str, alt: bool) -> Vec<u8> {
    protein
        .bytes()
        .flat_map(|a| {
            match a {
                b'M' => "ATG",
                b'K' => "AAA",
                b'T' => "ACC",
                b'A' if alt => "GCT",
                b'A' => "GCC",
                b'Y' => "TAC",
                b'I' => "ATC",
                b'V' => "GTT",
                b'Q' => "CAG",
                b'R' => "CGT",
                b'S' => "TCT",
                b'F' => "TTC",
                b'H' => "CAC",
                _ => "NNN",
            }
            .bytes()
        })
        .collect()
}

/// Run the alignment the page is asked for to its end.
fn wait_for_alignment(app: &mut App, ctx: &egui::Context) {
    let end = Instant::now() + Duration::from_secs(60);
    while app.platform.alignment.is_running() {
        app.alignment_tick(ctx);
        assert!(Instant::now() < end, "gave up waiting");
        std::thread::sleep(Duration::from_millis(15));
    }
}

fn coding_page(app: &mut App) -> Vec<(String, Vec<u8>)> {
    let proteins = [
        ("whole", "MKTAYIAKQRQISFVKSHFSRQ", false),
        ("missing two", "MKTIAKQRQISFVKSHFSRQ", true),
        ("one changed", "MKTAYVAKQRQISFVKSHFSRQ", false),
        ("two more", "MKTAYIAKQKKRQISFVKSHFSRQ", true),
    ];
    let page = &mut app.platform.alignment;
    page.clear_inputs();
    let seqs: Vec<(String, Vec<u8>)> = proteins
        .iter()
        .map(|(n, p, alt)| (n.to_string(), coding(p, *alt)))
        .collect();
    seqs.iter().for_each(|(n, b)| page.add(n, b));
    seqs
}

#[test]
fn aligned_by_translation_the_dna_keeps_its_codons_whole_and_the_record_says_how_they_were_read() {
    let mut app = App::new(Startup::default());
    let ctx = egui::Context::default();
    let seqs = coding_page(&mut app);
    // the built-in method is chosen, which has no scores for amino acids: the proteins are aligned another way, and the result says so
    app.platform.alignment.codons = Some(Codons::default());
    app.start_alignment(&ctx);
    wait_for_alignment(&mut app, &ctx);
    let page = &app.platform.alignment;
    let result = page.result.as_ref().expect("an alignment");
    assert_eq!(result.rows.len(), 4);
    for ((name, original), (got, row)) in seqs.iter().zip(&result.rows) {
        assert_eq!(name, got);
        let bases: Vec<u8> = row.iter().copied().filter(|&b| b != b'-').collect();
        assert_eq!(&bases, original, "{name} loses or moves no base");
        assert_eq!(row.len() % 3, 0);
        for codon in row.chunks(3) {
            let gaps = codon.iter().filter(|&&b| b == b'-').count();
            assert!(gaps == 0 || gaps == 3, "{name}: {codon:?}");
        }
    }
    // the two missing residues are one run of six gaps, and the two extra ones are six gaps in the other rows
    let gaps = |name: &str| {
        let row = &result.rows.iter().find(|r| r.0 == name).expect("a row").1;
        row.iter().filter(|&&b| b == b'-').count()
    };
    assert_eq!((gaps("whole"), gaps("missing two")), (6, 12));
    assert!(
        result.engine.contains("carried back to the codons"),
        "{}",
        result.engine
    );
    assert!(
        result.engine.contains("reading frame 1"),
        "{}",
        result.engine
    );
    let (record, _) = page.made.as_ref().expect("the record");
    assert!(
        record
            .to_text()
            .contains("reading frame 1, genetic code 1 (Standard)"),
        "{}",
        record.to_text()
    );
    // the viewer takes it for DNA
    assert!(!crate::seq::alnstats::is_protein(
        &result.rows.iter().map(|r| r.1.clone()).collect::<Vec<_>>()
    ));
}

#[test]
fn the_translation_is_kept_with_the_page_and_run_again_from_the_record() {
    let mut app = App::new(Startup::default());
    let ctx = egui::Context::default();
    coding_page(&mut app);
    app.platform.alignment.codons = Some(Codons { frame: 0, code: 1 });
    app.start_alignment(&ctx);
    wait_for_alignment(&mut app, &ctx);
    // kept, and put back, it is the same choice; a line that says another thing than a choice leaves the page as it was
    let kept = app.platform.alignment.kept();
    assert!(kept.contains("codons\t0\t1"), "{kept}");
    let mut back = AlignmentPage::default();
    back.put_back(&kept).expect("put back");
    assert_eq!(back.codons, Some(Codons { frame: 0, code: 1 }));
    let odd = kept.replace("codons\t0\t1", "codons\t7\t1");
    let mut other = AlignmentPage::default();
    other.put_back(&odd).expect("put back");
    assert_eq!(other.codons, None);
    let before = app.platform.alignment.stamp();
    app.platform.alignment.codons = Some(Codons { frame: 1, code: 1 });
    assert_ne!(app.platform.alignment.stamp(), before);
    // run again takes the choice from the record of what is shown, not from what the page shows now
    app.platform.alignment.codons = None;
    let mut out = Vec::new();
    app.apply_alignment_clicks(
        Clicks {
            made: made_ui::Asks {
                again: true,
                ..made_ui::Asks::default()
            },
            ..Clicks::default()
        },
        &ctx,
        &mut out,
    );
    assert_eq!(
        app.platform.alignment.codons,
        Some(Codons { frame: 0, code: 1 })
    );
    wait_for_alignment(&mut app, &ctx);
    let engine = &app
        .platform
        .alignment
        .result
        .as_ref()
        .expect("an alignment")
        .engine;
    assert!(
        engine.contains("genetic code 2 (Vertebrate mitochondrial)"),
        "{engine}"
    );
}

#[test]
fn a_protein_or_a_sequence_without_a_codon_is_refused_on_the_page_with_the_reason() {
    let mut app = App::new(Startup::default());
    let ctx = egui::Context::default();
    let page = &mut app.platform.alignment;
    page.add("dna", &coding("MKTAYIAKQ", false));
    page.add("protein", b"MKVLAAGIVGLLLAQPDEFHKRSTW");
    page.codons = Some(Codons::default());
    app.start_alignment(&ctx);
    wait_for_alignment(&mut app, &ctx);
    let page = &app.platform.alignment;
    assert!(page.result.is_none());
    assert!(
        page.message.contains("protein does not look like DNA"),
        "{}",
        page.message
    );
    // and the same when it is one codon short
    let page = &mut app.platform.alignment;
    page.clear_inputs();
    page.add("dna", &coding("MKTAYIAKQ", false));
    page.add("short", b"AT");
    app.start_alignment(&ctx);
    wait_for_alignment(&mut app, &ctx);
    assert!(
        app.platform
            .alignment
            .message
            .contains("short has no whole codon"),
        "{}",
        app.platform.alignment.message
    );
}

#[test]
fn the_choice_is_on_the_page_with_the_frame_and_the_code_once_it_is_ticked() {
    let mut app = App::new(Startup::default());
    app.workspace = crate::ui::platform::Workspace::Alignment;
    let mut screen = crate::ui::testing::Screen::new();
    screen.settle(&mut app);
    assert!(screen.text().contains("Align by the translation"));
    assert!(!screen.text().contains("Reading frame"));
    app.platform.alignment.codons = Some(Codons::default());
    screen.settle(&mut app);
    let text = screen.text();
    assert!(
        text.contains("Reading frame") && text.contains("Genetic code"),
        "{text}"
    );
}
