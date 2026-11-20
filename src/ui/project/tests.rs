use super::*;
use crate::pipeline::steps::Aligned;
use crate::seq::doc::Doc;
use crate::seq::features::Feature;
use crate::ui::app::{App, Startup};
use crate::ui::figure::take_saved;
use crate::ui::platform::route::{kind_of, Kind, Source};
use crate::ui::platform::Workspace;
use std::path::Path;

fn doc(name: &str, bases: &str, circular: bool) -> Doc {
    let mut d = Doc::new(name, bases.as_bytes().to_vec());
    d.circular = circular;
    d.features
        .push(Feature::new("a\tfeature", "gene", 1, 4, true));
    d
}

/// A session that holds everything a project holds: names that have a tab and a line end in them, texts of several lines.
fn session() -> Session {
    Session {
        workspace: "alignment".into(),
        active: 1,
        docs: vec![
            SavedDoc {
                name: "pET28a\tGFP\nv2".into(),
                path: None,
                dirty: true,
                text: Some("LOCUS       x\nORIGIN\n//\n".into()),
                keep: None,
            },
            SavedDoc {
                name: "second".into(),
                path: None,
                dirty: true,
                text: Some("LOCUS       y\r\n".into()),
                keep: None,
            },
        ],
        lost: 0,
        results: vec![SavedResult {
            kind: "alignment".into(),
            text: Some("engine\tbuiltin\ninput\t1\ta\tACGT\n".into()),
            keep: None,
        }],
    }
}

#[test]
fn a_project_is_read_back_as_it_was_written_with_every_text_whole() {
    let text = write(&session()).expect("written");
    assert!(text.starts_with("helixloom-project\t1\n"), "{text}");
    assert_eq!(read(&text), Ok(session()));
    // nothing but lines: no line end and no tab is left inside a field
    assert!(text.lines().all(|l| !l.is_empty()));
    // a line that a later version adds is passed over
    let (body, end) = text.rsplit_once("end\t").expect("an end line");
    let _ = end;
    let with_extra = format!("{body}later\twhatever\n");
    let sealed = seal(&with_extra, 2, 1);
    assert_eq!(read(&sealed), Ok(session()));
}

#[test]
fn what_cannot_be_a_project_is_not_written() {
    let mut s = session();
    s.docs[0].text = None;
    assert!(write(&s).expect_err("a text is missing").contains("text"));
    let mut s = session();
    s.results[0].text = None;
    assert!(write(&s).is_err());
    assert!(write(&Session::default())
        .expect_err("nothing")
        .contains("Nothing"));
}

#[test]
fn a_file_that_is_cut_short_changed_of_a_later_version_or_not_a_project_is_refused_in_words() {
    let text = write(&session()).expect("written");
    let lines: Vec<&str> = text.lines().collect();
    // cut short, at the end line and before it
    let cut = lines[..lines.len() - 1].join("\n");
    assert!(read(&cut).expect_err("cut").contains("cut short"));
    let half = &text[..text.len() / 2];
    assert!(read(half).expect_err("cut").contains("cut short"));
    // one letter changed
    let changed = text.replacen("second", "secone", 1);
    assert!(read(&changed).expect_err("changed").contains("changed"));
    // another version, and something else altogether
    let later = text.replacen("helixloom-project\t1", "helixloom-project\t2", 1);
    assert!(read(&later).expect_err("later").contains("later version"));
    for other in ["", ">a\nACGT\n", "LOCUS x\n", "helixloom-session\t1\n"] {
        assert!(read(other)
            .expect_err("not a project")
            .contains("not a HelixLoom project"));
    }
    // a project with nothing in it
    let empty = seal(
        "helixloom-project\t1\nworkspace\tsequences\nactive\t0\n",
        0,
        0,
    );
    assert!(read(&empty).expect_err("empty").contains("nothing"));
    // a size that is over the limit is refused before it is looked at
    assert!(read_at_most(&text, 10)
        .expect_err("too large")
        .contains("too large"));
}

#[test]
fn the_work_on_the_pages_is_saved_as_a_project_and_opened_again_in_another_window() {
    let mut first = App::new(Startup::default());
    first.docs.clear();
    first.add_doc(doc("circular one", "ACGTACGTACGT", true));
    first.add_doc(doc("linear two", "GGGGCCCCAAAA", false));
    first.platform.alignment.add("one", b"ACGTA");
    first.platform.alignment.add("two", b"ACTTA");
    first.platform.alignment.set_result(Aligned {
        rows: vec![
            ("one".to_string(), b"ACGTA".to_vec()),
            ("two".to_string(), b"ACTTA".to_vec()),
        ],
        engine: "built in".to_string(),
        seconds: 0.5,
    });
    first.workspace = Workspace::Alignment;
    let _ = take_saved();
    first.save_project();
    let saved = take_saved();
    assert_eq!(saved.len(), 1);
    assert_eq!(saved[0].0, FILE);
    assert!(saved[0].0.ends_with(".hlproj"));

    let mut second = App::new(Startup::default());
    second.docs.clear();
    second.open_source(Source::Bytes {
        name: "mine.hlproj",
        bytes: &saved[0].1,
    });
    let names: Vec<&str> = second.docs.iter().map(|d| d.doc.name.as_str()).collect();
    assert_eq!(names, ["circular one", "linear two"]);
    let one = &second.docs[0].doc;
    assert_eq!(one.seq, b"ACGTACGTACGT");
    assert!(one.circular && one.features.iter().any(|f| f.name == "a\tfeature"));
    assert!(!second.docs[1].doc.circular);
    // tied to no file: saving asks where, and closing asks whether
    assert!(second
        .docs
        .iter()
        .all(|d| d.doc.path.is_none() && d.doc.dirty));
    // the alignment page is as it was, and the workspace
    let rows = &second
        .platform
        .alignment
        .result
        .as_ref()
        .expect("an alignment")
        .rows;
    assert_eq!(rows[1].1, b"ACTTA");
    assert_eq!(second.platform.alignment.inputs.len(), 2);
    assert_eq!(second.workspace, Workspace::Alignment);
    // the tab that was in front is in front again (the second sequence: it was the last that was opened)
    assert!(matches!(second.tab, crate::ui::app::Tab::Doc(1)));
}

#[test]
fn a_project_that_is_opened_says_what_came_back_and_one_that_cannot_be_read_changes_nothing() {
    let mut app = App::new(Startup::default());
    app.docs.clear();
    app.add_doc(doc("kept", "ACGTACGT", false));
    let text = write(&session_of(&["from the project"])).expect("written");
    app.open_source(Source::Bytes {
        name: "p.hlproj",
        bytes: text.as_bytes(),
    });
    assert!(app.docs.iter().any(|d| d.doc.name == "from the project"));
    assert!(app.docs.iter().any(|d| d.doc.name == "kept"));
    let said = app.status.clone();
    assert!(
        said.contains("p.hlproj") && said.contains("1 sequence"),
        "{said}"
    );
    // a project that was cut short opens nothing, and says why
    let count = app.docs.len();
    app.open_source(Source::Bytes {
        name: "cut.hlproj",
        bytes: &text.as_bytes()[..text.len() / 2],
    });
    assert_eq!(app.docs.len(), count);
    let said = app.status.clone();
    assert!(
        said.contains("cut.hlproj") && said.contains("cut short"),
        "{said}"
    );
}

#[test]
fn there_is_nothing_to_save_when_nothing_is_open_and_the_example_is_not_work() {
    let mut app = App::new(Startup::default());
    let _ = take_saved();
    app.save_project();
    assert!(take_saved().is_empty());
    assert!(app.status.contains("Nothing"), "{}", app.status);
}

#[test]
fn the_ending_of_the_name_makes_a_file_a_project_whatever_its_case() {
    assert_eq!(kind_of(Path::new("work/my.hlproj")), Kind::Project);
    assert_eq!(kind_of(Path::new("MY.HLPROJ")), Kind::Project);
    assert_eq!(kind_of(Path::new("my.gb")), Kind::Sequence);
}

/// A session of sequences with these names.
fn session_of(names: &[&str]) -> Session {
    Session {
        workspace: "sequences".into(),
        active: 0,
        docs: names
            .iter()
            .map(|n| SavedDoc {
                name: n.to_string(),
                path: None,
                dirty: true,
                text: Some(crate::seq::io::to_genbank(&doc(n, "ACGTACGTAC", false))),
                keep: None,
            })
            .collect(),
        lost: 0,
        results: Vec::new(),
    }
}
