//! What a browser brings (a name and bytes, no path) is opened by the same loaders as a file on this computer.

use super::*;
use crate::sample;
use crate::seq::{io, scf};
use crate::ui::app::Startup;
use crate::ui::platform::route::Source;

const NEWICK: &str = "((A:1,B:2):1,(C:1,D:1):2);";

fn app() -> App {
    App::new(Startup::default())
}

/// A file the browser brought, as `take_arrived` gets it.
fn arrived(wanted: Wanted, name: &str, bytes: &[u8]) -> Arrived {
    Arrived {
        wanted,
        files: vec![(name.to_string(), Ok(bytes.to_vec()))],
    }
}

#[test]
fn a_file_without_a_path_goes_where_its_name_says() {
    let mut app = app();
    let before = app.docs.len();
    let genbank = io::to_genbank(&sample::demo_plasmid());
    app.open_source(Source::Bytes {
        name: "pDemo.gb",
        bytes: genbank.as_bytes(),
    });
    assert_eq!(app.workspace, Workspace::Sequences);
    assert_eq!(app.docs.len(), before + 1);
    assert!(app.docs[before].doc.path.is_none());
    assert!(!app.docs[before].doc.dirty);
    assert_eq!(app.status, "Opened 1 sequence from pDemo.gb");

    app.open_source(Source::Bytes {
        name: "1UBQ.PDB",
        bytes: include_bytes!("../../../assets/structure/1ubq.pdb"),
    });
    assert_eq!(app.workspace, Workspace::Structure);
    let loaded = app.platform.structure.loaded.as_ref();
    assert_eq!(loaded.map(|l| l.name.as_str()), Some("1UBQ.PDB"));

    app.open_source(Source::Bytes {
        name: "strains.nwk",
        bytes: NEWICK.as_bytes(),
    });
    assert_eq!(app.workspace, Workspace::Trees);
    assert!(app.tree_win.tree.is_some());

    let trace = scf::to_scf(&sample::demo_trace());
    app.open_source(Source::Bytes {
        name: "run1.scf",
        bytes: &trace,
    });
    assert_eq!(app.workspace, Workspace::Signals);
    let loaded = app.platform.signals.loaded.as_ref();
    assert_eq!(loaded.map(|l| l.name.as_str()), Some("run1.scf"));
}

#[test]
fn bytes_and_a_path_open_the_same_sequence() {
    let genbank = io::to_genbank(&sample::demo_plasmid());
    let folder = crate::pipeline::tools::work_folder("incoming").expect("a folder");
    let path = folder.join("pDemo.gb");
    std::fs::write(&path, &genbank).expect("a file");

    let (mut by_path, mut by_bytes) = (app(), app());
    by_path.open_file(&path);
    by_bytes.open_source(Source::Bytes {
        name: "pDemo.gb",
        bytes: genbank.as_bytes(),
    });
    let (a, b) = (
        &by_path.docs.last().expect("a sequence").doc,
        &by_bytes.docs.last().expect("a sequence").doc,
    );
    assert_eq!((&a.name, &a.seq, a.circular), (&b.name, &b.seq, b.circular));
    assert_eq!(a.features.len(), b.features.len());
    assert_eq!(a.path.as_deref(), Some(path.as_path()));
    assert!(b.path.is_none());
}

#[test]
fn a_file_asked_for_by_a_page_goes_to_that_page_whatever_its_name() {
    let mut app = app();
    app.take_arrived(arrived(Wanted::Tree, "tree.txt", NEWICK.as_bytes()));
    assert!(app.tree_win.tree.is_some());
    assert!(app.platform.trees.message.starts_with("Opened tree.txt"));

    app.take_arrived(arrived(Wanted::Fastq, "reads.txt", b"@r1\nACGT\n+\nIIII\n"));
    assert_eq!(app.reads_win.reads.len(), 1);

    app.take_arrived(arrived(Wanted::TreeToolsTree, "t.txt", NEWICK.as_bytes()));
    assert_eq!(app.treetools_win.newick, NEWICK);

    app.take_arrived(arrived(Wanted::PhyloCnnTree, "t.txt", NEWICK.as_bytes()));
    assert_eq!(app.phylocnn_win.newick, NEWICK);

    app.take_arrived(arrived(
        Wanted::Alignment,
        "two.fa",
        b">a\nACGT\n>b\nACGA\n",
    ));
    assert_eq!(app.platform.alignment.inputs.len(), 2);
}

#[test]
fn files_for_the_assembler_arrive_as_one_batch() {
    let mut app = app();
    app.take_arrived(Arrived {
        wanted: Wanted::AssemblyReads,
        files: vec![
            ("a.fa".into(), Ok(b">a\nACGTACGTAC\n".to_vec())),
            ("b.fa".into(), Ok(b">b\nTTGGCCAATT\n".to_vec())),
        ],
    });
    assert_eq!(app.assembly.reads.len(), 2);
    assert_eq!(app.assembly.message, "Added 2 sequences from 2 files.");
}

#[test]
fn a_file_the_browser_could_not_read_is_said_so() {
    let mut app = app();
    app.take_arrived(Arrived {
        wanted: Wanted::Anything,
        files: vec![("big.gb".into(), Err("File is too large".into()))],
    });
    assert_eq!(app.status, "Couldn't read big.gb: File is too large");
}

#[test]
fn a_file_that_is_not_text_is_refused_where_text_is_needed() {
    let mut app = app();
    app.take_arrived(arrived(Wanted::Tree, "tree.nwk", &[0xff, 0xfe, 0x00]));
    assert!(app.tree_win.tree.is_none());
    let said = &app.platform.trees.message;
    assert!(said.starts_with("Could not read tree.nwk"), "{said}");
}
