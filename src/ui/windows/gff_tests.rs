use crate::seq::doc::Doc;
use crate::ui::app::{App, Startup};
use crate::ui::platform::route::Source;

const GFF: &str = "##gff-version 3
pGM\tProkka\tCDS\t101\t400\t.\t+\t0\tID=c1;Name=repA
pGM\tProkka\tCDS\t501\t700\t.\t-\t0\tID=c2
";

fn app_with(name: &str, length: usize) -> App {
    let mut app = App::new(Startup::default());
    app.docs.clear();
    app.add_doc(Doc::new(name, vec![b'A'; length]));
    app
}

#[test]
fn an_annotation_opened_while_a_sequence_is_open_puts_its_genes_on_it_and_undo_takes_them_off() {
    let mut app = app_with("pGM", 1000);
    app.open_source(Source::Bytes {
        name: "p.gff3",
        bytes: GFF.as_bytes(),
    });
    let names: Vec<&str> = app.docs[0]
        .doc
        .features
        .iter()
        .map(|f| f.name.as_str())
        .collect();
    assert_eq!(names, ["repA", "c2"]);
    assert!(app.status.starts_with("Added 2 features"), "{}", app.status);
    app.docs[0].undo();
    assert!(app.docs[0].doc.features.is_empty());
}

#[test]
fn an_annotation_that_cannot_be_used_says_why_and_changes_nothing() {
    // no sequence open
    let mut none = App::new(Startup::default());
    none.docs.clear();
    none.load_gff(Source::Bytes {
        name: "p.gff3",
        bytes: GFF.as_bytes(),
    });
    assert!(none.status.contains("Open the sequence"), "{}", none.status);
    // a file that is no annotation, and one of another sequence
    let mut app = app_with("pGM", 1000);
    app.load_gff(Source::Bytes {
        name: "p.gff3",
        bytes: b"not an annotation\n",
    });
    assert!(
        app.status.contains("Could not use p.gff3"),
        "{}",
        app.status
    );
    let two = format!("{GFF}other\t.\tCDS\t1\t30\t.\t+\t0\tID=x\n");
    let mut other = app_with("something else", 1000);
    other.load_gff(Source::Bytes {
        name: "p.gff3",
        bytes: two.as_bytes(),
    });
    assert!(
        other.status.contains("pGM") && other.status.contains("something else"),
        "{}",
        other.status
    );
    assert!(app.docs[0].doc.features.is_empty() && other.docs[0].doc.features.is_empty());
}
