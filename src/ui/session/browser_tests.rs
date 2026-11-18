//! The store in a browser: the same things kept in a table of names and texts (`medium::Memory` stands where the browser's
//! storage is), by the same code that keeps them in the settings folder.

use super::medium::{Medium, Memory};
use super::*;
use std::rc::Rc;

/// A store on a table, and the table to look into.
fn page() -> (Store, Rc<Memory>) {
    let table = Rc::new(Memory::default());
    (Store::on(table.clone()), table)
}

fn names(table: &Memory) -> Vec<String> {
    table.names()
}

fn unsaved(name: &str, text: &str) -> SavedDoc {
    SavedDoc {
        name: name.into(),
        path: None,
        dirty: true,
        text: Some(text.into()),
        keep: None,
    }
}

fn session(docs: Vec<SavedDoc>) -> Session {
    Session {
        workspace: "alignment".into(),
        active: 1,
        docs,
        lost: 0,
        results: Vec::new(),
    }
}

#[test]
fn a_session_kept_in_a_browser_comes_back_as_it_was() {
    let (store, table) = page();
    assert!(store.in_a_browser());
    assert_eq!(store.load_session(), None);
    let kept = session(vec![
        unsaved("pDemo", "LOCUS one\n"),
        unsaved("a\tname", "LOCUS two\n"),
    ]);
    let files = store.save_session(&kept).expect("saved");
    assert_eq!(
        files,
        [Some("1-0.gb".to_string()), Some("1-1.gb".to_string())]
    );
    assert_eq!(
        names(&table),
        ["session/1-0.gb", "session/1-1.gb", "session/manifest.txt"]
    );
    // another visit: a store on the same table
    let again = Store::on(table.clone());
    assert_eq!(again.load_session(), Some(kept));
}

#[test]
fn a_new_save_takes_the_old_texts_away_and_a_save_that_fails_leaves_the_session_before_it_whole() {
    let (store, table) = page();
    let first = session(vec![unsaved("d", "first")]);
    store.save_session(&first).expect("saved");
    // something of another page of the site, and something of this program's that is no sequence
    table.set("session/notes", "mine").expect("set");
    table.set("places.txt", "").expect("set");

    // the list of what is open cannot be written: the texts of the new save are there, and the session is still the first
    *table.refused.borrow_mut() = Some("session/manifest.txt".into());
    let second = session(vec![unsaved("d", "second"), unsaved("e", "more")]);
    assert!(store.save_session(&second).is_err());
    assert_eq!(store.load_session(), Some(first.clone()));

    // a text cannot be written: the same
    *table.refused.borrow_mut() = Some("session/2-1.gb".into());
    assert!(store.save_session(&second).is_err());
    assert_eq!(store.load_session(), Some(first));

    // and when it goes through, the texts of the saves before it are gone, and nothing else is
    *table.refused.borrow_mut() = None;
    store.save_session(&second).expect("saved");
    assert_eq!(store.load_session(), Some(second));
    assert_eq!(
        names(&table),
        [
            "places.txt",
            "session/2-0.gb",
            "session/2-1.gb",
            "session/manifest.txt",
            "session/notes"
        ]
    );
}

#[test]
fn a_sequence_that_has_not_changed_keeps_its_text_in_a_browser_too() {
    let (store, table) = page();
    store
        .save_session(&session(vec![unsaved("d", "same"), unsaved("e", "old")]))
        .expect("saved");
    // the first has not changed since: its text is kept as it is, under the name it has
    let unchanged = SavedDoc {
        text: None,
        keep: Some("1-0.gb".into()),
        ..unsaved("d", "")
    };
    let files = store
        .save_session(&session(vec![unchanged.clone(), unsaved("e", "new")]))
        .expect("saved");
    assert_eq!(
        files,
        [Some("1-0.gb".to_string()), Some("2-1.gb".to_string())]
    );
    assert_eq!(
        names(&table),
        ["session/1-0.gb", "session/2-1.gb", "session/manifest.txt"]
    );
    let back = store.load_session().expect("a session");
    assert_eq!(back.docs[0].text.as_deref(), Some("same"));
    assert_eq!(back.docs[1].text.as_deref(), Some("new"));
    // a kept text that is gone is an error, not a sequence lost in silence
    table.remove("session/1-0.gb");
    let gone = store.save_session(&session(vec![unchanged]));
    assert!(gone.is_err_and(|why| why.contains("1-0.gb")));
}

#[test]
fn a_full_storage_says_so_and_what_was_kept_before_stays() {
    let (store, table) = page();
    let small = session(vec![unsaved("d", "ACGT")]);
    store.save_session(&small).expect("saved");
    let used: usize = table
        .held
        .borrow()
        .iter()
        .map(|(k, v)| k.len() + v.len())
        .sum();
    // room for what is there and a little more: a long sequence does not fit
    table.room.set(Some(used + 200));
    let long = session(vec![unsaved("d", &"ACGT".repeat(500))]);
    let why = store.save_session(&long).expect_err("it does not fit");
    assert_eq!(why, browser::FULL);
    assert_eq!(store.load_session(), Some(small));
    // the settings are told the same way
    table.room.set(Some(used));
    let appearance = crate::ui::appearance::Appearance::default();
    assert_eq!(
        store.save_appearance(ThemeMode::Dark, &appearance),
        Err(browser::FULL.to_string())
    );
}

#[test]
fn nothing_open_clears_the_session_and_a_session_can_be_put_aside_and_cleared() {
    let (store, table) = page();
    store
        .save_session(&session(vec![unsaved("d", "one")]))
        .expect("saved");
    table.set("appearance.txt", "theme\tdark\n").expect("set");
    // put aside: it is no longer the session, and it is all there under the other name
    store.archive_session().expect("put aside");
    assert_eq!(store.load_session(), None);
    assert_eq!(
        names(&table),
        [
            "appearance.txt",
            "session.previous/1-0.gb",
            "session.previous/manifest.txt"
        ]
    );
    // putting aside again when there is no session is nothing; a second one takes the place of the first
    store.archive_session().expect("nothing to put aside");
    store
        .save_session(&session(vec![unsaved("e", "two"), unsaved("f", "three")]))
        .expect("saved");
    store.archive_session().expect("put aside");
    assert_eq!(
        names(&table),
        [
            "appearance.txt",
            "session.previous/1-0.gb",
            "session.previous/1-1.gb",
            "session.previous/manifest.txt"
        ]
    );
    assert_eq!(table.get("session.previous/1-0.gb").as_deref(), Some("two"));
    store.clear_previous().expect("cleared");
    assert_eq!(names(&table), ["appearance.txt"]);
    // with nothing open there is nothing to come back to
    store
        .save_session(&session(vec![unsaved("d", "one")]))
        .expect("saved");
    store.save_session(&session(Vec::new())).expect("cleared");
    assert_eq!(names(&table), ["appearance.txt"]);
    store.clear_session().expect("nothing to clear");
}

#[test]
fn a_session_of_another_version_is_put_aside_not_written_over_and_a_text_outside_the_session_is_not_read(
) {
    let (store, table) = page();
    table
        .set("session/manifest.txt", "helixloom-session\t99\nsomething\n")
        .expect("set");
    table.set("session/9-0.gb", "theirs").expect("set");
    assert_eq!(store.load_session(), None);
    store
        .save_session(&session(vec![unsaved("d", "mine")]))
        .expect("saved");
    assert_eq!(
        table.get("session.unreadable/9-0.gb").as_deref(),
        Some("theirs")
    );
    assert_eq!(
        store.load_session().expect("a session").docs[0]
            .text
            .as_deref(),
        Some("mine")
    );
    // a list that names a text of another folder: it is not this program's, and is counted as lost
    table.set("appearance.txt", "secret").expect("set");
    table
        .set(
            "session/manifest.txt",
            "helixloom-session\t1\ngeneration\t5\ndoc\t1\td\t\t../appearance.txt\ndoc\t1\te\t\t1-0.gb\n",
        )
        .expect("set");
    let back = store.load_session().expect("a session");
    assert_eq!(back.lost, 1);
    assert_eq!(back.docs.len(), 1);
    assert_eq!(back.docs[0].text.as_deref(), Some("mine"));
}

#[test]
fn the_appearance_the_places_and_a_mark_are_kept_in_a_browser() {
    use crate::ui::appearance::Appearance;
    let (store, table) = page();
    assert_eq!(store.appearance(), None);
    assert_eq!(store.places(), Places::default());
    assert!(!store.has_mark("tour-seen"));

    let appearance = Appearance::default().with_spec("scale=150,contrast=high");
    store
        .save_appearance(ThemeMode::Dark, &appearance)
        .expect("kept");
    let places = Places::parse("12\t34\tORFs\n");
    store.save_places(&places).expect("kept");
    store.keep_mark("tour-seen").expect("kept");
    assert_eq!(names(&table), ["appearance.txt", "places.txt", "tour-seen"]);

    // another visit
    let again = Store::on(table.clone());
    assert_eq!(again.appearance(), Some((ThemeMode::Dark, appearance)));
    assert_eq!(again.places(), places);
    assert!(again.has_mark("tour-seen"));
    // what is not a text of this program's reads as the defaults, and only the first few thousand letters are looked at
    table
        .set(
            "appearance.txt",
            &format!("{}theme\tdark\n", "é".repeat(5000)),
        )
        .expect("set");
    assert_eq!(
        again.appearance(),
        Some((ThemeMode::default(), Appearance::default()))
    );
}

#[test]
fn a_page_has_no_lock_no_recent_files_and_no_place_for_its_window() {
    let (store, table) = page();
    // no file to lock: the page takes the session over with a mark (the test below), and keeps nothing else for it
    assert!(matches!(store.lock_session(), Lock::Taken(_)));
    assert_eq!(names(&table), ["session.owner"]);
    // a file that was opened in a page has no path to open it by again
    store
        .add_recent(Path::new("/no/such/file.gb"))
        .expect("nothing fails");
    assert!(store.recent().is_empty());
    assert_eq!(store.window(), None);
}

#[test]
fn what_a_page_keeps_is_under_names_of_its_own() {
    assert_eq!(browser::entry("places.txt"), "helixloom/places.txt");
    assert_eq!(
        browser::kept_name("helixloom/session/1-0.gb"),
        Some("session/1-0.gb")
    );
    // another program's entry, and the bare prefix, are not this program's
    assert_eq!(browser::kept_name("other/places.txt"), None);
    assert_eq!(browser::kept_name("helixloom"), None);
    assert_eq!(browser::kept_name("helixloom/"), None);
    assert_eq!(
        browser::kept_name(&browser::entry("tour-seen")),
        Some("tour-seen")
    );
}

fn mark_of(store: &Store) -> String {
    match store.lock_session() {
        Lock::Taken(mark) => mark,
        other => panic!("a page takes the session over: {other:?}"),
    }
}

#[test]
fn the_page_opened_last_owns_the_session_and_a_page_before_it_no_longer_does() {
    let (first, table) = page();
    let mine = mark_of(&first);
    assert!(first.owns(Some(&mine)));
    // another page of the site, opened after it (in the same thousandth of a second, as a test is)
    let second = Store::on(table.clone());
    let theirs = mark_of(&second);
    assert_ne!(mine, theirs);
    assert!(
        theirs.parse::<u64>().expect("a number") > mine.parse::<u64>().expect("a number"),
        "{theirs} after {mine}"
    );
    assert!(second.owns(Some(&theirs)));
    assert!(
        !first.owns(Some(&mine)),
        "the first page's mark is no longer the one kept"
    );
    // the mark is no part of the session: clearing the session leaves it
    second.clear_session().expect("cleared");
    assert!(second.owns(Some(&theirs)));
    // a mark that is not a number (someone else's) is taken over all the same
    table.set("session.owner", "not a number").expect("set");
    assert!(!second.owns(Some(&theirs)));
    let third = mark_of(&Store::on(table.clone()));
    assert!(third.parse::<u64>().is_ok());
    // a page that could write no mark (its storage is full) has none to lose: it goes on, and is told when a save fails
    table.room.set(Some(0));
    let full = Store::on(table.clone());
    assert!(matches!(full.lock_session(), Lock::Unavailable));
    assert!(full.owns(None));
}

#[test]
fn a_folder_of_the_computer_is_always_its_copys_to_write() {
    let dir = std::env::temp_dir().join(format!("helixloom-owns-{}", std::process::id()));
    let store = Store::open(Some(&dir)).expect("a store");
    assert!(!store.in_a_browser());
    assert!(store.owns(None) && store.owns(Some("anything")));
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_text_of_a_table_is_read_whole_or_up_to_so_many_bytes_cut_where_a_letter_ends() {
    let table = Rc::new(Memory::default());
    table.set("note", "aé").expect("set");
    let medium = Medium::Shelf(table.clone());
    assert_eq!(medium.read("note", None).as_deref(), Some("aé"));
    assert_eq!(medium.read("note", Some(99)).as_deref(), Some("aé"));
    assert_eq!(medium.read("note", Some(3)).as_deref(), Some("aé"));
    // the é is two bytes: two bytes of the text end inside it, and are cut back to the a
    assert_eq!(medium.read("note", Some(2)).as_deref(), Some("a"));
    assert_eq!(medium.read("note", Some(1)).as_deref(), Some("a"));
    assert_eq!(medium.read("note", Some(0)).as_deref(), Some(""));
    assert_eq!(medium.read("no such name", None), None);
    assert_eq!(medium.read("no such name", Some(4)), None);
    assert!(medium.holds("note") && !medium.holds("no such name"));
}

#[test]
fn a_folder_of_a_table_is_the_names_that_begin_with_it_and_a_stroke_and_no_others() {
    let (store, table) = page();
    // only a session that was put aside is there: it is not the session
    table
        .set("session.previous/manifest.txt", "x")
        .expect("set");
    table
        .set("session.previous/1-0.gb", "kept aside")
        .expect("set");
    let medium = Medium::Shelf(table.clone());
    assert!(!medium.has_folder("session"));
    assert!(medium.has_folder("session.previous"));
    assert!(medium.names_in("session").is_empty());
    assert_eq!(
        medium.names_in("session.previous"),
        ["1-0.gb", "manifest.txt"]
    );
    // so there is nothing to put aside, and clearing the session leaves what was put aside
    store.archive_session().expect("nothing to put aside");
    store.clear_session().expect("nothing to clear");
    medium.remove_folder("session").expect("nothing to remove");
    assert_eq!(
        names(&table),
        ["session.previous/1-0.gb", "session.previous/manifest.txt"]
    );
    assert_eq!(
        table.get("session.previous/1-0.gb").as_deref(),
        Some("kept aside")
    );
}

#[test]
fn putting_a_session_aside_takes_the_place_of_all_that_was_put_aside_before() {
    let (store, table) = page();
    // put aside before: three texts
    for (name, text) in [
        ("1-0.gb", "old"),
        ("1-1.gb", "older"),
        ("manifest.txt", "old list"),
    ] {
        table
            .set(&format!("session.previous/{name}"), text)
            .expect("set");
    }
    store
        .save_session(&session(vec![unsaved("d", "new")]))
        .expect("saved");
    store.archive_session().expect("put aside");
    // exactly the session that was just put aside: the second text of the one before is gone
    assert_eq!(
        names(&table),
        ["session.previous/1-0.gb", "session.previous/manifest.txt"]
    );
    assert_eq!(table.get("session.previous/1-0.gb").as_deref(), Some("new"));
}

#[test]
fn a_session_that_cannot_be_put_aside_whole_stays_where_it_is_whole() {
    let (store, table) = page();
    let kept = session(vec![unsaved("d", "one"), unsaved("e", "two")]);
    store.save_session(&kept).expect("saved");
    let before = names(&table);
    // the last of its texts cannot be written under the other name (the table is full at that point)
    *table.refused.borrow_mut() = Some("session.previous/manifest.txt".into());
    assert!(store.archive_session().is_err());
    // nothing has moved: the session is the session, and nothing is half put aside
    assert_eq!(names(&table), before);
    assert_eq!(store.load_session(), Some(kept.clone()));
    // and when it can be, it is all there
    *table.refused.borrow_mut() = None;
    store.archive_session().expect("put aside");
    assert_eq!(store.load_session(), None);
    assert_eq!(
        names(&table),
        [
            "session.previous/1-0.gb",
            "session.previous/1-1.gb",
            "session.previous/manifest.txt"
        ]
    );
}

#[test]
fn a_session_of_another_version_that_cannot_be_put_aside_is_not_saved_into() {
    let (store, table) = page();
    table
        .set("session/manifest.txt", "helixloom-session\t99\n")
        .expect("set");
    // named as this program names its texts
    table.set("session/1-0.gb", "theirs").expect("set");
    *table.refused.borrow_mut() = Some("session.unreadable/manifest.txt".into());
    let why = store
        .save_session(&session(vec![unsaved("d", "mine")]))
        .expect_err("not saved");
    assert!(why.contains("another version"), "{why}");
    assert_eq!(table.get("session/1-0.gb").as_deref(), Some("theirs"));
    assert_eq!(
        table.get("session/manifest.txt").as_deref(),
        Some("helixloom-session\t99\n")
    );
    // put aside twice, what is aside is the second, not both
    *table.refused.borrow_mut() = None;
    store
        .save_session(&session(vec![unsaved("d", "mine")]))
        .expect("saved");
    assert_eq!(
        table.get("session.unreadable/1-0.gb").as_deref(),
        Some("theirs")
    );
    table
        .set("session/manifest.txt", "something else\n")
        .expect("set");
    store
        .save_session(&session(vec![unsaved("d", "mine again")]))
        .expect("saved");
    assert_eq!(
        Medium::Shelf(table.clone()).names_in("session.unreadable"),
        ["1-0.gb", "manifest.txt"]
    );
    assert_eq!(
        table.get("session.unreadable/1-0.gb").as_deref(),
        Some("mine")
    );
    assert_eq!(
        table.get("session.unreadable/manifest.txt").as_deref(),
        Some("something else\n")
    );
    // the numbering goes on from nothing after a list that was not this program's
    assert!(table.get("session/1-0.gb").is_some());
}

#[test]
fn a_text_of_the_session_that_is_gone_is_counted_and_the_rest_come_back() {
    let (store, table) = page();
    store
        .save_session(&session(vec![unsaved("d", "one"), unsaved("e", "two")]))
        .expect("saved");
    table.remove("session/1-0.gb");
    let back = store.load_session().expect("a session");
    assert_eq!(back.lost, 1);
    assert_eq!(back.docs.len(), 1);
    assert_eq!(back.docs[0].name, "e");
    assert_eq!(back.docs[0].text.as_deref(), Some("two"));
}

#[test]
fn the_order_a_table_gives_its_names_in_does_not_matter() {
    // a browser gives the names of its table in no order to count on: here, last first
    let (store, table) = page();
    table.backwards.set(true);
    let first = session(vec![unsaved("d", "one"), unsaved("e", "two")]);
    store.save_session(&first).expect("saved");
    assert_eq!(store.load_session(), Some(first));
    let second = session(vec![unsaved("d", "three")]);
    store.save_session(&second).expect("saved");
    assert_eq!(store.load_session(), Some(second.clone()));
    let mut left = names(&table);
    left.sort();
    assert_eq!(left, ["session/2-0.gb", "session/manifest.txt"]);
    store.archive_session().expect("put aside");
    let mut aside = names(&table);
    aside.sort();
    assert_eq!(
        aside,
        ["session.previous/2-0.gb", "session.previous/manifest.txt"]
    );
    store.clear_previous().expect("cleared");
    assert!(names(&table).is_empty());
}

#[test]
fn a_place_for_the_window_and_an_empty_list_of_recent_files_are_kept_like_anything_else() {
    let (store, table) = page();
    let place = WindowPlace {
        x: 10.0,
        y: 20.0,
        width: 900.0,
        height: 700.0,
    };
    store.save_window(&place).expect("kept");
    assert_eq!(Store::on(table.clone()).window(), Some(place));
    store.clear_recent().expect("kept");
    assert_eq!(table.get("recent.txt").as_deref(), Some(""));
}

#[test]
fn a_mark_makes_the_settings_folder_it_is_kept_in_and_a_folder_of_its_name_is_no_mark() {
    let dir = std::env::temp_dir().join(format!("helixloom-mark-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    let store = Store::open(Some(&dir)).expect("a store");
    // nothing has made the folder yet: keeping a mark does
    assert!(!dir.exists());
    store.keep_mark("tour-seen").expect("kept");
    assert!(store.has_mark("tour-seen"));
    // a folder that happens to have a mark's name is not that mark
    fs::create_dir_all(dir.join("other-mark")).expect("a folder");
    assert!(!store.has_mark("other-mark"));
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_mark_is_the_time_it_was_taken_at_so_two_pages_that_read_the_same_last_mark_take_different_ones(
) {
    let thousandths = || {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("a time after 1970")
            .as_millis() as u64
    };
    let (store, _table) = page();
    let before = thousandths();
    let mark: u64 = mark_of(&store).parse().expect("a number");
    // not one more than nothing: the time, which two pages opened a moment apart do not share
    assert!(
        (before..=thousandths()).contains(&mark),
        "{mark} is not the time"
    );
}

fn result(kind: &str, text: &str) -> SavedResult {
    SavedResult {
        kind: kind.into(),
        text: Some(text.into()),
        keep: None,
    }
}

fn unchanged(kind: &str, file: &str) -> SavedResult {
    SavedResult {
        kind: kind.into(),
        text: None,
        keep: Some(file.into()),
    }
}

fn with_results(docs: Vec<SavedDoc>, results: Vec<SavedResult>) -> Session {
    Session {
        results,
        ..session(docs)
    }
}

fn letters_held(table: &Memory) -> usize {
    table
        .held
        .borrow()
        .iter()
        .map(|(k, v)| k.len() + v.len())
        .sum()
}

#[test]
fn what_the_pages_show_is_kept_beside_the_sequences_and_alone_it_is_a_session_too() {
    let (store, table) = page();
    let kept = with_results(
        vec![unsaved("d", "LOCUS one\n")],
        vec![
            result("alignment", "rows\twith tabs\nand lines\n"),
            result("tree", "(a,b);"),
        ],
    );
    let files = store.keep_session(&kept).expect("saved");
    assert_eq!(files.docs, [Some("1-0.gb".to_string())]);
    assert_eq!(
        files.results,
        [
            Some("1-alignment.result".to_string()),
            Some("1-tree.result".to_string())
        ]
    );
    assert_eq!(
        names(&table),
        [
            "session/1-0.gb",
            "session/1-alignment.result",
            "session/1-tree.result",
            "session/manifest.txt"
        ]
    );
    assert_eq!(Store::on(table.clone()).load_session(), Some(kept));
    // with no sequence open, what a page shows is still something to come back to
    let alone = with_results(Vec::new(), vec![result("structure", "ATOM\n")]);
    store.keep_session(&alone).expect("saved");
    assert_eq!(
        names(&table),
        ["session/2-structure.result", "session/manifest.txt"]
    );
    assert_eq!(store.load_session(), Some(alone));
    // and with neither there is nothing
    store.keep_session(&session(Vec::new())).expect("cleared");
    assert!(names(&table).is_empty());
}

#[test]
fn a_result_that_has_not_changed_keeps_its_file_and_one_that_its_page_no_longer_shows_takes_its_file_with_it(
) {
    let (store, table) = page();
    let docs = || vec![unsaved("d", "one")];
    store
        .keep_session(&with_results(
            docs(),
            vec![result("alignment", "A"), result("tree", "T")],
        ))
        .expect("saved");
    // a file of that ending that this program did not write is nobody's to remove
    table.set("session/my-notes.result", "mine").expect("set");
    let second = with_results(docs(), vec![unchanged("alignment", "1-alignment.result")]);
    let files = store.keep_session(&second).expect("saved");
    assert_eq!(files.results, [Some("1-alignment.result".to_string())]);
    assert_eq!(
        names(&table),
        [
            "session/1-alignment.result",
            "session/2-0.gb",
            "session/manifest.txt",
            "session/my-notes.result"
        ]
    );
    let back = store.load_session().expect("a session");
    assert_eq!(back.results, [result("alignment", "A")]);
    // a kept result whose file is gone is left out, which the save says by giving no file for it; the sequences are kept
    table.remove("session/1-alignment.result");
    let files = store.keep_session(&second).expect("saved");
    assert_eq!(files.results, [None]);
    let back = store.load_session().expect("a session");
    assert_eq!((back.docs.len(), back.results.len()), (1, 0));
}

#[test]
fn a_result_that_does_not_fit_is_left_out_and_a_new_one_takes_the_room_of_the_one_it_replaces() {
    let (store, table) = page();
    let tree = |text: &str| with_results(vec![unsaved("d", "ACGT")], vec![result("tree", text)]);
    store.keep_session(&tree("(a,b);")).expect("saved");
    // room for a tree of 600 letters, and not for two of them
    table.room.set(Some(letters_held(&table) + 700));
    let files = store.keep_session(&tree(&"a".repeat(600))).expect("saved");
    assert_eq!(files.results, [Some("2-tree.result".to_string())]);
    // the tree has changed: the one of the last save, which the page no longer shows, makes room for it
    let files = store.keep_session(&tree(&"b".repeat(600))).expect("saved");
    assert_eq!(files.results, [Some("3-tree.result".to_string())]);
    let back = store.load_session().expect("a session");
    assert_eq!(back.results, [result("tree", &"b".repeat(600))]);
    // one that fits nowhere is left out, and the sequences are kept
    let huge = with_results(
        vec![unsaved("d", "TTTT")],
        vec![result("tree", &"c".repeat(5000))],
    );
    let files = store.keep_session(&huge).expect("saved");
    assert_eq!(
        (files.docs, files.results),
        (vec![Some("4-0.gb".to_string())], vec![None])
    );
    let back = store.load_session().expect("a session");
    assert_eq!(back.docs[0].text.as_deref(), Some("TTTT"));
    assert!(back.results.is_empty());
    assert_eq!(names(&table), ["session/4-0.gb", "session/manifest.txt"]);
}

#[test]
fn the_sequences_come_first_when_the_results_that_are_kept_leave_them_no_room() {
    let (store, table) = page();
    store
        .keep_session(&with_results(
            vec![unsaved("d", "ACGT")],
            vec![result("structure", &"x".repeat(2000))],
        ))
        .expect("saved");
    table.room.set(Some(letters_held(&table) + 100));
    // the sequence has grown, and beside the structure it does not fit: the structure is what gives way
    let grown = with_results(
        vec![unsaved("d", &"ACGT".repeat(300))],
        vec![unchanged("structure", "1-structure.result")],
    );
    let files = store.keep_session(&grown).expect("saved");
    assert_eq!(
        (files.docs, files.results),
        (vec![Some("2-0.gb".to_string())], vec![None])
    );
    let back = store.load_session().expect("a session");
    assert_eq!(
        back.docs[0].text.as_deref(),
        Some("ACGT".repeat(300).as_str())
    );
    assert!(back.results.is_empty());
    // with no result to give way, a sequence that does not fit is still an error, and what was kept stays
    let more = session(vec![unsaved("d", &"ACGT".repeat(2000))]);
    assert_eq!(
        store.keep_session(&more).map(|_| ()),
        Err(browser::FULL.to_string())
    );
    assert_eq!(
        store.load_session().expect("a session").docs[0].text,
        back.docs[0].text
    );
}

#[test]
fn a_result_under_a_name_that_is_not_a_page_or_in_a_file_that_is_elsewhere_or_gone_is_neither_written_nor_read(
) {
    let (store, table) = page();
    let files = store
        .keep_session(&with_results(
            vec![unsaved("d", "one")],
            vec![
                result("tree", "T"),
                result("../tree", "out"),
                result("", "none"),
            ],
        ))
        .expect("saved");
    assert_eq!(
        files.results,
        [Some("1-tree.result".to_string()), None, None]
    );
    assert_eq!(
        names(&table),
        [
            "session/1-0.gb",
            "session/1-tree.result",
            "session/manifest.txt"
        ]
    );
    table.set("elsewhere.txt", "not a result").expect("set");
    let manifest = table.get("session/manifest.txt").expect("a manifest");
    // a file that is elsewhere, a kind that is no page's, a file of another kind, and a file that is gone
    let more = "result\ttree\t../elsewhere.txt\nresult\tTree Of Life\t1-tree.result\n\
                result\talignment\t1-tree.result\nresult\talignment\t1-alignment.result\n";
    table
        .set("session/manifest.txt", &format!("{manifest}{more}"))
        .expect("set");
    let back = store.load_session().expect("a session");
    assert_eq!(back.results, [result("tree", "T")]);
}

#[test]
fn a_save_that_fails_for_another_reason_than_room_leaves_the_results_that_are_kept_alone() {
    let (store, table) = page();
    store
        .keep_session(&with_results(
            vec![unsaved("d", "one")],
            vec![result("tree", "T")],
        ))
        .expect("saved");
    // the text of a sequence, kept from the last save, is gone: that is an error, and making room would not mend it
    table.remove("session/1-0.gb");
    let same = SavedDoc {
        text: None,
        keep: Some("1-0.gb".into()),
        ..unsaved("d", "")
    };
    let again = with_results(vec![same], vec![unchanged("tree", "1-tree.result")]);
    let why = store.keep_session(&again).expect_err("a text is gone");
    assert!(why.contains("1-0.gb"), "{why}");
    assert_eq!(table.get("session/1-tree.result").as_deref(), Some("T"));
    // a session of another version that cannot be put aside is not saved into, and its files are not removed either
    let (store, table) = page();
    table
        .set("session/manifest.txt", "helixloom-session\t99\n")
        .expect("set");
    table.set("session/7-tree.result", "theirs").expect("set");
    *table.refused.borrow_mut() = Some("session.unreadable/manifest.txt".into());
    let new = with_results(vec![unsaved("d", "one")], vec![result("tree", "T")]);
    assert!(store.keep_session(&new).is_err());
    assert_eq!(
        table.get("session/7-tree.result").as_deref(),
        Some("theirs")
    );
}

fn list_of(table: &Memory) -> String {
    table.get("session/manifest.txt").expect("a list")
}

fn result_line(table: &Memory) -> Vec<String> {
    let list = list_of(table);
    let line = list
        .lines()
        .find(|l| l.starts_with("result\t"))
        .expect("a result");
    line.split('\t').map(str::to_string).collect()
}

#[test]
fn a_result_is_listed_with_its_length_and_checksum_and_a_text_that_is_not_that_one_is_left_out() {
    let (store, table) = page();
    let tree = "(one,(two,three));";
    store
        .keep_session(&with_results(
            vec![unsaved("d", "LOCUS")],
            vec![result("tree", tree)],
        ))
        .expect("saved");
    let line = result_line(&table);
    assert_eq!(line[..3], ["result", "tree", "1-tree.result"]);
    assert_eq!(line[3], tree.len().to_string());
    assert_eq!(line[4], format!("{:016x}", crate::kept::hash(tree)));
    assert_eq!(store.load_session().expect("a session").results.len(), 1);
    // a text of the same length (what a refused write leaves behind) is not the one that was listed: left out, the sequence stays
    table
        .set(
            "session/1-tree.result",
            "(one,(two,three));".replace("two", "TWO").as_str(),
        )
        .expect("set");
    let back = store.load_session().expect("a session");
    assert_eq!((back.docs.len(), back.results.len()), (1, 0));
    // a text of another length is left out too
    table
        .set("session/1-tree.result", "(one,two);")
        .expect("set");
    assert!(store.load_session().expect("a session").results.is_empty());
}

#[test]
fn a_result_that_is_kept_as_it_is_is_listed_again_with_the_same_length_and_checksum_without_being_read(
) {
    let (store, table) = page();
    store
        .keep_session(&with_results(
            vec![unsaved("d", "one")],
            vec![result("tree", "(a,b);")],
        ))
        .expect("saved");
    let first = result_line(&table);
    // the next save keeps the file as it is, and the list says of it what it said
    let again = with_results(
        vec![unsaved("d", "two")],
        vec![unchanged("tree", "1-tree.result")],
    );
    store.keep_session(&again).expect("saved");
    assert_eq!(result_line(&table), first);
    // and it is still checked: a text put in its place is not read back
    table.set("session/1-tree.result", "(b,a);").expect("set");
    assert!(store.load_session().expect("a session").results.is_empty());
}

#[test]
fn a_list_of_an_earlier_version_that_gives_no_checksum_is_read_as_it_was_and_a_checksum_of_one_field_is_not(
) {
    let (store, table) = page();
    store
        .keep_session(&with_results(
            vec![unsaved("d", "one")],
            vec![result("tree", "(a,b);")],
        ))
        .expect("saved");
    let line = result_line(&table);
    let list = list_of(&table);
    let without = list.replace(&line.join("\t"), &line[..3].join("\t"));
    table.set("session/manifest.txt", &without).expect("set");
    assert_eq!(store.load_session().expect("a session").results.len(), 1);
    // a line with a length and no checksum is damaged, not older: it is not read
    let damaged = list.replace(&line.join("\t"), &line[..4].join("\t"));
    table.set("session/manifest.txt", &damaged).expect("set");
    assert!(store.load_session().expect("a session").results.is_empty());
    // a result that is kept as it is, from a list that gave no checksum, is listed without one, and is read as before
    table.set("session/manifest.txt", &without).expect("set");
    let again = with_results(
        vec![unsaved("d", "two")],
        vec![unchanged("tree", "1-tree.result")],
    );
    store.keep_session(&again).expect("saved");
    assert_eq!(result_line(&table).len(), 3);
    assert_eq!(store.load_session().expect("a session").results.len(), 1);
}
