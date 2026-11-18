//! What the pages show is kept with the session: written with it, offered with it, and put back on its pages.

use super::*;
use crate::pipeline::steps::Aligned;
use crate::seq::newick::NwTree;
use crate::ui::results::{self, ALIGNMENT, STRUCTURE, TREE};
use crate::ui::session::Memory;
use crate::ui::testing::Screen;
use std::rc::Rc;

fn tree(text: &str) -> crate::seq::tree::Node {
    NwTree::parse_all(text).expect("a tree")[0].to_node()
}

/// Something on each of the three pages.
fn fill(app: &mut App) {
    let page = &mut app.platform.alignment;
    page.add("one", b"ACGTA");
    page.add("two", b"ACTTA");
    page.set_result(Aligned {
        rows: vec![
            ("one".to_string(), b"ACGTA".to_vec()),
            ("two".to_string(), b"ACTTA".to_vec()),
        ],
        engine: "built in".to_string(),
        seconds: 0.5,
    });
    app.show_tree(tree("((a:1,b:2)90:1,c:3);"), "shown");
    app.platform.structure.load_example();
}

fn texts(app: &App) -> [Option<String>; 3] {
    [ALIGNMENT, TREE, STRUCTURE].map(|kind| app.result_text(kind))
}

fn kinds(session: &Session) -> Vec<&str> {
    session.results.iter().map(|r| r.kind.as_str()).collect()
}

#[test]
fn what_the_pages_show_is_offered_with_the_last_session_and_put_back_on_its_pages() {
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(dirty("kept", "ACGTACGT"));
    fill(&mut app);
    app.workspace = Workspace::Trees;
    save_now(&mut app);
    assert_eq!(
        folder.names("session"),
        [
            "1-0.gb",
            "1-alignment.result",
            "1-structure.result",
            "1-tree.result",
            "manifest.txt"
        ]
    );
    let shown = texts(&app);
    assert!(shown.iter().all(Option::is_some));
    drop(app);
    let mut next = folder.app();
    let offer = next.keeper.offer.clone().expect("an offer");
    assert_eq!(kinds(&offer), [ALIGNMENT, TREE, STRUCTURE]);
    // nothing is on a page before the person has said so
    assert!(next.result_stamps().is_empty());
    next.answer_offer(true);
    assert_eq!(texts(&next), shown);
    assert_eq!(next.workspace, Workspace::Trees);
    assert_eq!(
        next.status,
        "Reopened 1 sequence from the last session, with the Alignment page, the tree and the structure"
    );
    // what came back is the session now: written once more, it is the same
    save_now(&mut next);
    drop(next);
    let again = offer_of(&folder).expect("an offer");
    assert_eq!(kinds(&again), [ALIGNMENT, TREE, STRUCTURE]);
}

#[test]
fn with_no_sequence_open_what_a_page_shows_is_still_a_session_and_nothing_at_all_is_none() {
    let folder = Folder::new();
    let mut app = folder.app();
    assert!(!app.keeper.pending(&app.session_signature()));
    app.flush_session();
    assert!(folder.names("session").is_empty());
    app.platform.alignment.add("pasted", b"ACGT");
    app.workspace = Workspace::Alignment;
    assert!(app.keeper.pending(&app.session_signature()));
    app.flush_session();
    assert!(!app.keeper.pending(&app.session_signature()));
    assert_eq!(
        folder.names("session"),
        ["1-alignment.result", "manifest.txt"]
    );
    drop(app);
    let mut next = folder.app();
    let offer = next.keeper.offer.clone().expect("an offer");
    assert_eq!((offer.docs.len(), kinds(&offer)), (0, vec![ALIGNMENT]));
    next.answer_offer(true);
    assert_eq!(next.platform.alignment.inputs.len(), 1);
    assert_eq!(next.workspace, Workspace::Alignment);
    assert_eq!(
        next.status,
        "Put back the Alignment page from the last session"
    );
    // "Start fresh" takes it all away
    drop(next);
    let mut fresh = folder.app();
    fresh.answer_offer(false);
    assert!(fresh.result_stamps().is_empty());
    assert!(folder.names("session").is_empty());
    drop(fresh);
    assert_eq!(offer_of(&folder), None);
}

#[test]
fn a_result_that_has_not_changed_is_not_written_again_one_that_has_is_and_one_that_is_gone_goes() {
    let folder = Folder::new();
    let mut app = folder.app();
    fill(&mut app);
    save_now(&mut app);
    let results = |folder: &Folder| -> Vec<String> {
        let names = folder.names("session").into_iter();
        names.filter(|n| n.ends_with(".result")).collect()
    };
    assert_eq!(
        results(&folder),
        ["1-alignment.result", "1-structure.result", "1-tree.result"]
    );
    // a sequence is opened: the results stay in the files they are in
    app.add_doc(dirty("new", "ACGT"));
    save_now(&mut app);
    assert_eq!(
        results(&folder),
        ["1-alignment.result", "1-structure.result", "1-tree.result"]
    );
    // another tree: it alone is written
    app.show_tree(tree("(x:1,y:1);"), "another");
    save_now(&mut app);
    assert_eq!(
        results(&folder),
        ["1-alignment.result", "1-structure.result", "3-tree.result"]
    );
    // the Tree window is asked for a tree of one sequence: there is no tree now, and none is kept
    app.build_tree();
    assert!(app.tree_win.tree.is_none());
    save_now(&mut app);
    assert_eq!(
        results(&folder),
        ["1-alignment.result", "1-structure.result"]
    );
    // a save that failed is not trusted: the next writes everything
    app.keeper.fail("failed".to_string(), Instant::now());
    let _ = app.keeper.take_problem();
    save_now(&mut app);
    assert_eq!(
        results(&folder),
        ["5-alignment.result", "5-structure.result"]
    );
    drop(app);
    let back = offer_of(&folder).expect("an offer");
    assert_eq!(kinds(&back), [ALIGNMENT, STRUCTURE]);
}

#[test]
fn in_a_browser_a_result_with_no_room_is_said_once_the_sequences_are_kept_and_it_is_kept_when_there_is_room(
) {
    let table = Rc::new(Memory::default());
    let mut app = page_on(&table);
    app.add_doc(dirty("typed", "ACGTACGT"));
    app.platform.structure.load_example();
    app.show_tree(tree("(a:1,b:1);"), "shown");
    // room for the sequence and the tree, and not for a structure
    table.room.set(Some(4_000));
    app.flush_session();
    let said = app.keeper.take_problem().expect("something to say");
    assert_eq!(
        said,
        "Could not keep the structure for next time: there is no room for it here, or it could not be written. \
         The open sequences are kept."
    );
    // it is shown as a failure is
    assert_eq!(
        crate::ui::toasts::tone_of(&said),
        crate::ui::design::Tone::Bad
    );
    let kept = Store::on(table.clone()).load_session().expect("a session");
    assert_eq!((kept.docs.len(), kinds(&kept)), (1, vec![TREE]));
    // the next save says nothing more of it
    app.add_doc(dirty("more", "GGGG"));
    app.flush_session();
    assert_eq!(app.keeper.take_problem(), None);
    assert!(!app.keeper.pending(&app.session_signature()));
    // room is made: it is tried at every save although it has not changed, so the next one keeps it, with nothing said
    table.room.set(None);
    app.add_doc(dirty("and more", "TTTT"));
    app.flush_session();
    assert_eq!(app.keeper.take_problem(), None);
    let kept = Store::on(table.clone()).load_session().expect("a session");
    assert_eq!((kept.docs.len(), kinds(&kept)), (3, vec![TREE, STRUCTURE]));
    // and when it is left out another time, that is said another time
    table.room.set(Some(4_000));
    app.platform.structure.load_example();
    app.flush_session();
    assert!(app
        .keeper
        .take_problem()
        .is_some_and(|said| said.contains("the structure")));
}

#[test]
fn turning_the_remembering_off_removes_the_results_and_on_again_they_are_written_anew() {
    let folder = Folder::new();
    let mut app = folder.app();
    fill(&mut app);
    save_now(&mut app);
    assert_eq!(app.keeper.set_keeping(false).as_deref(), Some(NOT_KEPT_NOW));
    assert!(folder.names("session").is_empty());
    assert!(!app.keeper.pending(&app.session_signature()));
    assert_eq!(app.keeper.set_keeping(true).as_deref(), Some(KEPT_AGAIN));
    app.flush_session();
    assert_eq!(
        folder.names("session"),
        [
            "1-alignment.result",
            "1-structure.result",
            "1-tree.result",
            "manifest.txt"
        ]
    );
}

#[test]
fn a_kept_result_that_cannot_be_read_is_said_and_the_rest_comes_back() {
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(dirty("kept", "ACGTACGT"));
    fill(&mut app);
    save_now(&mut app);
    drop(app);
    // a list of an earlier version has no length and checksum to tell the texts by: these two are read as they are
    let list = folder.0.join("session/manifest.txt");
    let earlier: String = std::fs::read_to_string(&list)
        .expect("a list")
        .lines()
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let keep = if fields.len() == 5 && fields[0] == "result" {
                3
            } else {
                fields.len()
            };
            format!("{}\n", fields[..keep].join("\t"))
        })
        .collect();
    std::fs::write(&list, earlier).expect("a list");
    folder.file("session/1-tree.result", "newick\t((broken;\n");
    folder.file(
        "session/1-structure.result",
        "name\tx\npdb\tnot a structure\n",
    );
    let mut next = folder.app();
    next.answer_offer(true);
    assert!(next.platform.alignment.result.is_some());
    assert!(next.tree_win.tree.is_none() && next.platform.structure.loaded.is_none());
    assert_eq!(
        next.status,
        "Reopened 1 sequence from the last session, with the Alignment page; \
         the tree and the structure could not be put back"
    );
}

#[test]
fn a_result_that_is_not_the_text_the_list_gives_the_length_and_checksum_of_is_not_offered_and_the_rest_is(
) {
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(dirty("kept", "ACGTACGT"));
    fill(&mut app);
    save_now(&mut app);
    drop(app);
    // a text of the same length under the name of the tree's file (what a refused write leaves behind, or another page wrote)
    let tree = folder.0.join("session/1-tree.result");
    let old = std::fs::read_to_string(&tree).expect("the tree");
    let other: String = old.chars().rev().collect();
    assert_eq!(other.len(), old.len());
    std::fs::write(&tree, other).expect("another text");
    let mut next = folder.app();
    let offer = next.keeper.offer.clone().expect("an offer");
    assert_eq!(kinds(&offer), [ALIGNMENT, STRUCTURE]);
    next.answer_offer(true);
    assert!(next.tree_win.tree.is_none());
    assert_eq!(
        next.status,
        "Reopened 1 sequence from the last session, with the Alignment page and the structure"
    );
}

#[test]
fn the_question_names_what_the_pages_showed_and_asks_about_them_alone_when_no_sequence_was_open() {
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(dirty("kept", "ACGTACGT"));
    fill(&mut app);
    save_now(&mut app);
    drop(app);
    let mut next = folder.app();
    let mut screen = Screen::new();
    screen.settle(&mut next);
    let text = screen.text();
    assert!(text.contains("Reopen your last session?"), "{text}");
    assert!(
        text.contains("Last time 1 sequence was open, 1 with changes that were not saved."),
        "{text}"
    );
    assert!(
        text.contains("Also kept: the Alignment page, the tree and the structure."),
        "{text}"
    );
    drop(next);
    // results alone
    let alone = Folder::new();
    let mut app = alone.app();
    app.show_tree(tree("(a:1,b:1);"), "shown");
    save_now(&mut app);
    drop(app);
    let mut next = alone.app();
    let mut screen = Screen::new();
    screen.settle(&mut next);
    let text = screen.text();
    assert!(text.contains("Reopen your last session?"), "{text}");
    assert!(
        text.contains("No sequence was open last time. Kept from it: the tree."),
        "{text}"
    );
    assert!(!text.contains("cannot be read any more"), "{text}");
    screen.click(&mut next, "Reopen");
    assert!(next.keeper.offer.is_none());
    assert!(next.tree_win.tree.is_some());
    assert_eq!(results::named([TREE]), "the tree");
}

#[test]
fn what_is_said_after_reopening_is_true_for_one_and_for_several_and_when_nothing_came_back() {
    let docs = |n: usize| Session {
        docs: vec![SavedDoc::default(); n],
        ..Session::default()
    };
    let kinds = |kinds: &[&str]| kinds.iter().map(|k| k.to_string()).collect::<Vec<_>>();
    let done = Restored::default;
    for (session, restored, said) in [
        (
            docs(2),
            Restored { opened: 2, ..done() },
            "Reopened 2 sequences from the last session",
        ),
        (
            docs(1),
            Restored { opened: 1, results: kinds(&[TREE]), ..done() },
            "Reopened 1 sequence from the last session, with the tree",
        ),
        (
            docs(0),
            Restored { results: kinds(&[TREE, STRUCTURE]), ..done() },
            "Put back the tree and the structure from the last session",
        ),
        (
            docs(0),
            Restored { results_lost: kinds(&[TREE]), ..done() },
            "Nothing was put back from the last session: the tree could not be put back",
        ),
        (
            docs(3),
            Restored { opened: 1, already_open: 1, lost: 1, ..done() },
            "Reopened 1 sequence from the last session; 1 was open already; \
             1 could not be reopened (its file has gone or cannot be read)",
        ),
        (
            docs(5),
            Restored { opened: 1, already_open: 2, lost: 2, ..done() },
            "Reopened 1 sequence from the last session; 2 were open already; \
             2 could not be reopened (their files have gone or cannot be read)",
        ),
        (
            docs(1),
            Restored { opened: 1, results_lost: kinds(&[STRUCTURE]), ..done() },
            "Reopened 1 sequence from the last session; the structure could not be put back",
        ),
        (
            // the only sequence of the session is one that is gone: it was a session of sequences all the same
            Session {
                lost: 1,
                ..docs(0)
            },
            Restored {
                lost: 1,
                results: kinds(&[TREE]),
                ..done()
            },
            "Reopened 0 sequences from the last session, with the tree; \
             1 could not be reopened (its file has gone or cannot be read)",
        ),
        (
            // a kind that a later version kept is not one this version can show: it is said, not passed over
            docs(1),
            Restored { opened: 1, results_lost: kinds(&[TREE, "pipeline"]), ..done() },
            "Reopened 1 sequence from the last session; the tree and a result of another version could not be put back",
        ),
        (
            docs(1),
            Restored { opened: 1, results_lost: kinds(&["pipeline"]), ..done() },
            "Reopened 1 sequence from the last session; a result of another version could not be put back",
        ),
    ] {
        assert_eq!(reopened(&session, &restored), said);
    }
}

#[test]
fn a_result_with_no_room_when_no_sequence_is_open_is_said_without_speaking_of_sequences() {
    let table = Rc::new(Memory::default());
    let mut app = page_on(&table);
    app.platform.structure.load_example();
    app.show_tree(tree("(a:1,b:1);"), "shown");
    // room for the tree or for the list of what is kept, and not for both: the list is what must be written
    let held: usize = table
        .held
        .borrow()
        .iter()
        .map(|(k, v)| k.len() + v.len())
        .sum();
    table.room.set(Some(held + 100));
    app.flush_session();
    assert_eq!(
        crate::ui::session::Shelf::names(&*table),
        ["session.owner", "session/manifest.txt"]
    );
    assert_eq!(
        app.keeper.take_problem().as_deref(),
        Some("Could not keep the tree and the structure for next time: there is no room for them here, or they could not be written.")
    );
}

#[test]
fn the_question_says_what_cannot_be_reopened_in_words_that_fit_one_and_several() {
    for (gone, said) in [
        (1, "1 more cannot be reopened: it cannot be read."),
        (2, "2 more cannot be reopened: they cannot be read."),
    ] {
        let folder = Folder::new();
        let mut app = folder.app();
        for name in ["one", "two", "three"] {
            app.add_doc(dirty(name, "ACGT"));
        }
        save_now(&mut app);
        drop(app);
        for n in 0..gone {
            let text = folder.0.join(format!("session/1-{n}.gb"));
            std::fs::remove_file(text).expect("removed");
        }
        let mut next = folder.app();
        let mut screen = Screen::new();
        screen.settle(&mut next);
        let text = screen.text();
        assert!(text.contains(said), "{text}");
    }
}

#[test]
fn a_run_of_the_pipelines_page_is_kept_offered_and_put_back_like_the_other_pages_results() {
    let folder = Folder::new();
    let mut app = folder.app();
    let run = crate::pipeline::Outcome {
        log: vec!["Align: 2 sequences.".to_string()],
        ..Default::default()
    };
    let page = &mut app.platform.pipelines;
    page.show_fetched(run, "fetched".to_string());
    app.workspace = Workspace::Pipelines;
    save_now(&mut app);
    assert_eq!(
        folder.names("session"),
        ["1-pipelines.result", "manifest.txt"]
    );
    drop(app);
    let mut next = folder.app();
    let offer = next.keeper.offer.clone().expect("an offer");
    assert_eq!(kinds(&offer), [results::PIPELINES]);
    assert!(next.platform.pipelines.outcome.is_none());
    next.answer_offer(true);
    assert_eq!(
        next.status,
        "Put back the Pipelines page from the last session"
    );
    let run = next.platform.pipelines.outcome.as_ref().expect("a run");
    assert_eq!(run.log, ["Align: 2 sequences."]);
    assert_eq!(next.workspace, Workspace::Pipelines);
}

#[test]
fn a_run_of_a_pipeline_file_is_kept_and_put_back_and_one_that_a_run_stands_in_the_way_of_is_said_truly(
) {
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(dirty("kept", "ACGTACGT"));
    app.show_example_flow_run();
    let report = app.platform.pipelines.flow.report.clone();
    assert!(report.is_some());
    save_now(&mut app);
    assert_eq!(
        folder.names("session"),
        ["1-0.gb", "1-pipelines.result", "manifest.txt"]
    );
    drop(app);
    let mut next = folder.app();
    next.answer_offer(true);
    assert_eq!(next.platform.pipelines.flow.report, report);
    assert_eq!(
        next.status,
        "Reopened 1 sequence from the last session, with the Pipelines page"
    );
    // a run has been asked for on the page (an example at the start): the kept run is not put over it, and what is said
    // is that it could not be put back, which is so; not that it could not be read, which it could
    save_now(&mut next);
    drop(next);
    let mut busy = folder.app();
    busy.platform.pipelines.run_requested = true;
    busy.answer_offer(true);
    assert!(busy.platform.pipelines.flow.report.is_none());
    assert_eq!(
        busy.status,
        "Reopened 1 sequence from the last session; the Pipelines page could not be put back"
    );
}
