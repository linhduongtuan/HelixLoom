use super::*;
use crate::connect::{Reply, Table};
use crate::ui::app::{Action, App, Startup};
use crate::ui::testing::Screen;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Mutex;

static NEXT: AtomicU32 = AtomicU32::new(0);

struct Folder(PathBuf);

impl Folder {
    fn new() -> Folder {
        let path = std::env::temp_dir().join(format!(
            "helixloom-fetchwin-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

const GENBANK: &str = "LOCUS       AB123456                  24 bp    DNA     linear   BCT 01-JAN-2020\nDEFINITION  A test.\nACCESSION   AB123456\nFEATURES             Location/Qualifiers\n     CDS             1..12\nORIGIN\n        1 atgaaatttg gccctaaaat ttaa\n//\n";
const FASTA: &str = ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF\n";
const PDB: &str = "HEADER    UBIQUITIN\nATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N\nATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00 10.38           C\nATOM      3  C   MET A   1      26.913  26.639   3.531  1.00  9.62           C\nATOM      4  O   MET A   1      27.886  26.463   4.263  1.00  9.62           O\nATOM      5  CB  MET A   1      25.112  24.880   3.649  1.00 13.77           C\nATOM      6  CG  MET A   1      25.353  24.860   5.134  1.00 16.29           C\nATOM      7  SD  MET A   1      23.930  23.959   5.904  1.00 17.17           S\nATOM      8  CE  MET A   1      24.447  23.984   7.620  1.00 16.11           C\nATOM      9  N   GLN A   2      26.335  27.770   3.258  1.00  9.27           N\nATOM     10  CA  GLN A   2      26.850  29.021   3.898  1.00  9.07           C\nATOM     11  C   GLN A   2      26.100  29.253   5.202  1.00  8.72           C\nATOM     12  O   GLN A   2      24.865  29.024   5.330  1.00  8.22           O\nEND\n";

fn table() -> Arc<Table> {
    let reply = |source: Source, id: &str, body: &str| {
        (
            Request::check(source, id).unwrap().url(),
            Reply {
                status: 200,
                body: body.to_string(),
            },
        )
    };
    Arc::new(Table {
        replies: vec![
            reply(Source::NcbiNucleotide, "AB123456.1", GENBANK),
            reply(Source::UniProt, "P69905", FASTA),
            reply(Source::Pdb, "1A3N", PDB),
        ],
        asked: Mutex::new(Vec::new()),
    })
}

/// An app with the window open, answered by `net`.
fn app(net: &Arc<Table>) -> (App, Screen) {
    let mut app = App::new(Startup::default());
    app.fetch_win.fetcher = Some(net.clone());
    app.run(Action::Fetch, &egui::Context::default());
    let mut screen = Screen::new();
    screen.settle(&mut app);
    (app, screen)
}

fn fetch(app: &mut App, screen: &mut Screen, source: Source, text: &str) {
    app.fetch_win.source = source;
    app.fetch_win.text = text.to_string();
    screen.click(app, "Fetch");
    let until = std::time::Instant::now() + Duration::from_secs(30);
    while app.fetch_win.is_running() {
        assert!(std::time::Instant::now() < until, "no answer");
        std::thread::sleep(Duration::from_millis(10));
        screen.frame(app, vec![]);
    }
    screen.settle(app);
}

#[test]
fn the_action_opens_the_window_with_the_database_and_the_accession_box() {
    let (_, screen) = app(&table());
    let text = screen.text();
    for want in [
        "Fetch from a database",
        "Database",
        "NCBI nucleotide",
        "Fetch",
        "Nothing is kept: there is no settings folder.",
    ] {
        assert!(text.contains(want), "no \"{want}\":\n{text}");
    }
}

#[test]
fn a_nucleotide_record_opens_in_sequences_with_where_it_came_from() {
    let net = table();
    let (mut app, mut screen) = app(&net);
    let before = app.docs.len();
    fetch(
        &mut app,
        &mut screen,
        Source::NcbiNucleotide,
        " ab123456.1 ",
    );
    assert_eq!(app.docs.len(), before + 1);
    let doc = &app.docs[app.docs.len() - 1].doc;
    assert_eq!(doc.len(), 24);
    assert!(
        doc.comment
            .contains("Fetched from NCBI nucleotide AB123456.1"),
        "{}",
        doc.comment
    );
    assert_eq!(app.workspace, Workspace::Sequences);
    assert_eq!(
        app.fetch_win.message,
        "Opened AB123456.1 (24 bases), downloaded."
    );
    let text = screen.text();
    assert!(
        text.contains("WHERE IT CAME FROM") && text.contains("Licence: NCBI data are public"),
        "{text}"
    );
}

#[test]
fn a_protein_goes_to_the_alignment_page_and_a_structure_to_the_structure_page() {
    let net = table();
    let (mut app, mut screen) = app(&net);
    let before = app.platform.alignment.inputs.len();
    fetch(&mut app, &mut screen, Source::UniProt, "P69905");
    assert_eq!(app.platform.alignment.inputs.len(), before + 1);
    assert!(
        app.fetch_win
            .message
            .contains("Added P69905 (47 residues) to the list of the Alignment page"),
        "{}",
        app.fetch_win.message
    );

    fetch(&mut app, &mut screen, Source::Pdb, "1a3n");
    assert!(app.platform.structure.loaded.is_some());
    assert_eq!(app.workspace, Workspace::Structure);
    assert!(
        app.fetch_win
            .message
            .starts_with("Opened the structure 1A3N"),
        "{}",
        app.fetch_win.message
    );
}

#[test]
fn a_bad_accession_is_refused_before_anything_is_asked() {
    let net = table();
    let (mut app, mut screen) = app(&net);
    fetch(&mut app, &mut screen, Source::UniProt, "P69905&x=1");
    assert!(
        app.fetch_win
            .message
            .contains("is not an accession of UniProt"),
        "{}",
        app.fetch_win.message
    );
    assert!(net.asked.lock().unwrap().is_empty());
}

#[test]
fn a_record_that_is_not_there_is_said_so() {
    let net = table();
    let (mut app, mut screen) = app(&net);
    fetch(&mut app, &mut screen, Source::NcbiNucleotide, "ZZ999999.1");
    assert_eq!(
        app.fetch_win.message,
        "NCBI nucleotide has no record called ZZ999999.1."
    );
    assert!(app.fetch_win.last.is_none());
}

#[test]
fn a_record_kept_in_the_settings_folder_is_not_asked_for_again() {
    let f = Folder::new();
    let net = table();
    let (mut app, mut screen) = app(&net);
    app.fetch_win.cache_dir = Some(f.0.join("downloads"));
    fetch(&mut app, &mut screen, Source::UniProt, "P69905");
    assert!(app.fetch_win.message.contains("downloaded"));
    fetch(&mut app, &mut screen, Source::UniProt, "P69905");
    assert!(
        app.fetch_win.message.contains("from the copy kept here"),
        "{}",
        app.fetch_win.message
    );
    assert_eq!(net.asked.lock().unwrap().len(), 1);
    assert!(screen.text().contains("kept since"));
    app.fetch_win.refresh = true;
    fetch(&mut app, &mut screen, Source::UniProt, "P69905");
    assert_eq!(
        net.asked.lock().unwrap().len(),
        2,
        "a fresh copy was asked for"
    );
}

#[test]
fn a_startup_with_a_settings_folder_keeps_the_records_in_it() {
    let f = Folder::new();
    let app = App::new(Startup {
        settings_dir: Some(f.0.clone()),
        ..Startup::default()
    });
    assert_eq!(app.fetch_win.cache_dir, Some(f.0.join("downloads")));
    assert_eq!(App::new(Startup::default()).fetch_win.cache_dir, None);
}
