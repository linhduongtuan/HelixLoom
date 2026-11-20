use super::stand_in::{Idb, Mode};
use super::*;
use crate::ui::session::medium::Memory;
use crate::ui::session::{SavedDoc, SavedResult, Session, Store};
use std::cell::Cell;

const ROOM_FOR_ALL: usize = 1 << 20;

fn small() -> Rc<Memory> {
    Rc::new(Memory::default())
}

fn rc(text: &str) -> Rc<str> {
    Rc::from(text)
}

#[test]
fn a_text_is_there_at_once_and_the_database_gets_it_behind_it() {
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    large.set("session/1-tree.result", "(a,b);").expect("kept");
    // asked for and not yet written: it is read, listed and held already
    assert_eq!(
        large.get("session/1-tree.result").as_deref(),
        Some("(a,b);")
    );
    assert!(large.holds("session/1-tree.result"));
    assert_eq!(large.names(), ["session/1-tree.result"]);
    assert!(idb.stored.borrow().is_empty());
    assert_eq!(idb.waiting(), 1);
    idb.finish();
    assert_eq!(
        idb.entries(),
        [("session/1-tree.result".to_string(), "(a,b);".to_string())]
    );
    // a page that is opened then finds it
    let next = idb.open(ROOM_FOR_ALL);
    assert_eq!(next.get("session/1-tree.result").as_deref(), Some("(a,b);"));
    assert_eq!(next.take_failure(), None);
    // nothing was written for the opening: what is loaded is what the database holds already
    assert_eq!(idb.ops.borrow().len(), 1);
}

#[test]
fn changes_go_to_the_database_one_at_a_time_each_when_the_one_before_is_over_and_in_the_order_made()
{
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    large.set("session/1-a.result", "one").expect("kept");
    // the first is handed over at once; the others wait behind it
    large.set("session/2-a.result", "two").expect("kept");
    large.remove("session/1-a.result");
    large.set("session/2-b.result", "three").expect("kept");
    assert_eq!((idb.waiting(), idb.ops.borrow().len()), (1, 1));
    idb.finish();
    assert_eq!((idb.waiting(), idb.ops.borrow().len()), (1, 2));
    idb.finish();
    idb.finish();
    assert_eq!(idb.waiting(), 1);
    assert_eq!(
        *idb.ops.borrow(),
        [
            Op::Put("session/1-a.result".into(), rc("one")),
            Op::Put("session/2-a.result".into(), rc("two")),
            Op::Delete("session/1-a.result".into()),
            Op::Put("session/2-b.result".into(), rc("three")),
        ]
    );
    idb.finish();
    assert_eq!(idb.waiting(), 0);
    assert_eq!(idb.names(), ["session/2-a.result", "session/2-b.result"]);
    // with nothing waiting, a change is handed over at once again
    large.remove("session/2-a.result");
    assert_eq!(idb.waiting(), 1);
}

#[test]
fn a_change_to_a_name_takes_the_place_of_the_ones_that_still_wait_for_that_name() {
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    large.set("session/0-blocker.result", "b").expect("kept");
    // two texts for one name, a text and then its removal, a removal and then a text
    large.set("session/1-a.result", "first").expect("kept");
    large.set("session/1-a.result", "second").expect("kept");
    large.set("session/1-b.result", "doomed").expect("kept");
    large.remove("session/1-b.result");
    large.remove("session/1-c.result");
    large.set("session/1-c.result", "back").expect("kept");
    for _ in 0..4 {
        idb.finish();
    }
    assert_eq!(idb.waiting(), 0);
    let asked: Vec<Op> = idb.ops.borrow().clone();
    assert_eq!(
        asked,
        [
            Op::Put("session/0-blocker.result".into(), rc("b")),
            Op::Put("session/1-a.result".into(), rc("second")),
            Op::Delete("session/1-b.result".into()),
            Op::Put("session/1-c.result".into(), rc("back")),
        ]
    );
    // what was never asked for was not held, and a text that is not there is not asked to be taken away
    assert_eq!(
        idb.names(),
        [
            "session/0-blocker.result",
            "session/1-a.result",
            "session/1-c.result"
        ]
    );
}

#[test]
fn a_text_that_would_not_fit_in_memory_is_refused_and_what_a_change_frees_is_there_to_use() {
    let idb = Idb::new(Mode::Now);
    let name = "session/1-a.result";
    // room for 100 bytes of names and texts
    let large = idb.open(100);
    large.set(name, &"x".repeat(60)).expect("fits");
    let why = large
        .set("session/1-b.result", &"y".repeat(60))
        .expect_err("does not fit");
    assert_eq!(why, FULL_LARGE);
    assert!(!large.holds("session/1-b.result"));
    // the same name again is counted once: a text that takes the place of one is measured against what is left of it
    large.set(name, &"z".repeat(70)).expect("replaces");
    assert_eq!(idb.entries()[0].1, "z".repeat(70));
    // a refused text asked nothing of the database, and one that is removed frees its room
    assert_eq!(idb.ops.borrow().len(), 2);
    large.remove(name);
    large
        .set("session/1-b.result", &"y".repeat(60))
        .expect("fits now");
    // what was loaded is held whatever its size, and takes room
    let full = Large::new(
        idb.clone(),
        vec![("session/9-c.result".into(), "c".repeat(500))],
        100,
    );
    assert!(full.holds("session/9-c.result"));
    assert!(full.set("session/9-d.result", "d").is_err());
    full.remove("session/9-c.result");
    full.set("session/9-d.result", "d").expect("room again");
}

#[test]
fn a_change_the_database_refuses_leaves_no_trace_in_memory_and_what_it_may_still_hold_is_taken_away_with_the_next_change(
) {
    let idb = Idb::new(Mode::Later);
    let name = "session/1-tree.result";
    idb.stored.borrow_mut().insert(name.into(), "A".into());
    let large = idb.open(ROOM_FOR_ALL);
    // B was to take the place of A, and is refused: the table holds neither, and the database still holds A under the name
    large.set(name, "B").expect("kept in memory");
    idb.fail("QuotaExceededError");
    assert!(!large.holds(name));
    assert_eq!(large.take_failure().as_deref(), Some("QuotaExceededError"));
    assert_eq!(idb.stored.borrow().get(name).map(String::as_str), Some("A"));
    // nothing is asked of the database at once (a database that refuses everything would be asked for ever)
    assert_eq!((idb.waiting(), idb.ops.borrow().len()), (0, 1));
    // the next change asked takes A away first
    idb.mode.set(Mode::Now);
    large.set("session/2-x.result", "x").expect("kept");
    assert_eq!(idb.names(), ["session/2-x.result"]);
    assert_eq!(
        idb.ops.borrow()[1..],
        [
            Op::Delete(name.into()),
            Op::Put("session/2-x.result".into(), rc("x"))
        ]
    );
    // and it is asked for once, not again
    large.set("session/3-x.result", "y").expect("kept");
    assert_eq!(idb.ops.borrow().len(), 4);
}

#[test]
fn a_removal_the_database_refuses_is_asked_again_with_the_next_change_even_though_memory_has_forgotten_the_name(
) {
    let idb = Idb::new(Mode::Later);
    let name = "session/1-tree.result";
    idb.stored.borrow_mut().insert(name.into(), "A".into());
    let large = idb.open(ROOM_FOR_ALL);
    large.remove(name);
    idb.fail("InvalidStateError");
    assert!(!large.holds(name) && large.names().is_empty());
    assert_eq!(idb.names(), [name]);
    // memory does not hold the name, so asking to remove it asks nothing: the database is told by the next change
    large.remove(name);
    assert_eq!(idb.ops.borrow().len(), 1);
    idb.mode.set(Mode::Now);
    large.set("session/2-x.result", "x").expect("kept");
    assert_eq!(idb.names(), ["session/2-x.result"]);
}

#[test]
fn a_newer_text_that_waits_takes_a_refused_ones_place_and_nothing_is_taken_away_for_it() {
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    let name = "session/1-a.result";
    large.set(name, "older").expect("kept");
    large.set(name, "newer").expect("kept");
    idb.fail("QuotaExceededError");
    // the refused one is not what the table holds, so it is not dropped; the newer is next, and the database is not asked to
    // take the name away
    assert_eq!(large.get(name).as_deref(), Some("newer"));
    idb.finish();
    assert_eq!(idb.entries(), [(name.to_string(), "newer".to_string())]);
    assert!(!idb
        .ops
        .borrow()
        .iter()
        .any(|op| matches!(op, Op::Delete(_))));
}

#[test]
fn a_clear_the_database_refuses_is_asked_again_with_the_next_change() {
    let idb = Idb::new(Mode::Later);
    for name in ["session/1-a.result", "session/1-b.result"] {
        idb.stored.borrow_mut().insert(name.into(), "t".into());
    }
    let large = idb.open(ROOM_FOR_ALL);
    large.clear();
    assert!(large.names().is_empty());
    idb.fail("InvalidStateError");
    assert_eq!(idb.stored.borrow().len(), 2);
    idb.mode.set(Mode::Now);
    large.set("session/2-c.result", "c").expect("kept");
    assert_eq!(idb.names(), ["session/2-c.result"]);
}

#[test]
fn one_text_the_database_has_no_room_for_does_not_take_the_others_with_it() {
    let idb = Idb::new(Mode::Now);
    idb.quota.set(Some(100));
    let large = idb.open(ROOM_FOR_ALL);
    large
        .set("session/1-a.result", &"a".repeat(50))
        .expect("kept");
    large
        .set("session/1-b.result", &"b".repeat(500))
        .expect("kept in memory");
    large
        .set("session/1-c.result", &"c".repeat(60))
        .expect("kept");
    assert_eq!(idb.names(), ["session/1-a.result", "session/1-c.result"]);
    // the table holds what the database does
    assert_eq!(large.names(), idb.names());
    assert_eq!(large.take_failure().as_deref(), Some("QuotaExceededError"));
}

#[test]
fn the_first_reason_is_kept_until_it_is_told_and_the_page_is_woken_at_every_refusal() {
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    let woken = Rc::new(Cell::new(0));
    let counted = woken.clone();
    large.set_waker(Rc::new(move || counted.set(counted.get() + 1)));
    large.set("session/1-a.result", "a").expect("kept");
    large.set("session/1-b.result", "b").expect("kept");
    idb.fail("first reason");
    idb.fail("second reason");
    assert_eq!(woken.get(), 2);
    assert_eq!(large.take_failure().as_deref(), Some("first reason"));
    assert_eq!(large.take_failure(), None, "told once");
    // the next reason is new
    large.set("session/1-c.result", "c").expect("kept");
    idb.fail("third reason");
    assert_eq!(large.take_failure().as_deref(), Some("third reason"));
    // a change that went well wakes nobody
    large.set("session/1-d.result", "d").expect("kept");
    idb.finish();
    assert_eq!(woken.get(), 3);
}

#[test]
fn a_result_is_renamed_inside_the_table_with_the_room_of_one_text_and_nothing_is_asked_twice() {
    let idb = Idb::new(Mode::Now);
    // room for the text once under either name (81 and 90 bytes of name and text), and not twice
    let large = idb.open(100);
    let (from, to) = ("session/1-tree.result", "session.previous/1-tree.result");
    large.set(from, &"t".repeat(60)).expect("fits");
    large.rename(from, to).expect("renamed");
    assert!(!large.holds(from));
    assert_eq!(large.get(to), Some("t".repeat(60)));
    assert_eq!(idb.names(), [to]);
    assert_eq!(
        idb.ops.borrow()[1..],
        [
            Op::Put(to.into(), rc(&"t".repeat(60))),
            Op::Delete(from.into())
        ]
    );
    // the room that is left is the room of one text: 90 bytes are held, and 22 more do not fit
    assert!(large.set("session/9-z.result", "z").is_err());
    large.remove(to);
    large.set("session/9-z.result", "z").expect("fits");
    // a name that is not held is nothing to move, and a rename that would not fit is refused with the table as it was
    assert_eq!(large.rename("session/none.result", to), Ok(()));
    let tight = Large::new(idb.clone(), vec![(from.to_string(), "t".repeat(60))], 85);
    assert!(tight.rename(from, to).is_err());
    assert!(tight.holds(from) && !tight.holds(to));
    // onto a name that is held, it replaces it, and counts once
    let two = Large::new(
        idb.clone(),
        vec![
            (from.to_string(), "1".repeat(20)),
            (to.to_string(), "2".repeat(20)),
        ],
        90,
    );
    two.rename(from, to).expect("renamed over");
    assert_eq!(two.get(to), Some("1".repeat(20)));
    // 50 bytes are held (the text counted once, under its new name): 31 more fit, and 41 do not
    assert!(two.set("session/9-z.result", &"z".repeat(30)).is_err());
    assert!(two.set("session/9-z.result", &"z".repeat(10)).is_ok());
}

#[test]
fn clearing_takes_everything_away_in_memory_at_once_and_in_the_database_in_its_turn_whoever_put_it_there(
) {
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    large.set("session/0-blocker.result", "b").expect("kept");
    large.set("session/1-a.result", "a").expect("kept");
    // what another page of the site put in the database after this one was opened is not known to this one
    idb.stored
        .borrow_mut()
        .insert("session/7-other.result".into(), "o".into());
    large.clear();
    assert!(large.names().is_empty());
    // what waited is of no use now: the change that is out, and the clear
    assert_eq!(idb.waiting(), 1);
    idb.finish();
    idb.finish();
    assert_eq!(idb.waiting(), 0);
    assert!(idb.entries().is_empty());
    assert_eq!(idb.ops.borrow().len(), 2);
    // and the table is used again after it
    large.set("session/2-c.result", "c").expect("kept");
    idb.finish();
    assert_eq!(idb.names(), ["session/2-c.result"]);
}

#[test]
fn only_results_go_to_the_larger_table_and_every_other_name_stays_in_the_small_one() {
    let table = small();
    let idb = Idb::new(Mode::Now);
    let roomy = Roomy::new(table.clone(), Some(idb.open(ROOM_FOR_ALL)));
    roomy.set("session/manifest.txt", "list").expect("kept");
    roomy.set("session/1-0.gb", "LOCUS").expect("kept");
    roomy.set("appearance.txt", "dark").expect("kept");
    roomy.set("session/1-tree.result", "(a,b);").expect("kept");
    roomy
        .set("session.previous/1-tree.result", "(c,d);")
        .expect("kept");
    assert_eq!(
        table.names(),
        ["appearance.txt", "session/1-0.gb", "session/manifest.txt"]
    );
    assert_eq!(
        idb.names(),
        ["session.previous/1-tree.result", "session/1-tree.result"]
    );
    // all of them are read, held and listed through the one
    assert_eq!(
        roomy.get("session/1-tree.result").as_deref(),
        Some("(a,b);")
    );
    assert_eq!(roomy.get("session/1-0.gb").as_deref(), Some("LOCUS"));
    assert!(roomy.holds("session/1-tree.result") && roomy.holds("appearance.txt"));
    assert!(!roomy.holds("session/2-tree.result"));
    assert_eq!(roomy.names().len(), 5);
    roomy.remove("session/1-tree.result");
    roomy.remove("session/1-0.gb");
    assert!(!roomy.holds("session/1-tree.result") && !roomy.holds("session/1-0.gb"));
    assert_eq!(idb.entries().len(), 1);
    assert_eq!(table.names(), ["appearance.txt", "session/manifest.txt"]);
}

#[test]
fn a_result_beyond_the_small_tables_room_is_kept_in_the_larger_one() {
    let table = small();
    table.room.set(Some(500));
    let idb = Idb::new(Mode::Now);
    let roomy = Roomy::new(table.clone(), Some(idb.open(ROOM_FOR_ALL)));
    let big = "x".repeat(10_000);
    roomy
        .set("session/1-alignment.result", &big)
        .expect("the larger table takes it");
    assert_eq!(roomy.get("session/1-alignment.result"), Some(big.clone()));
    // a name that is not a result is still the small table's, with its room
    assert!(roomy.set("session/1-0.gb", &big).is_err());
    // with no larger table the result is refused as the small one refuses
    let alone = Roomy::new(table, None);
    assert!(alone.set("session/2-alignment.result", &big).is_err());
}

#[test]
fn a_result_that_an_earlier_version_kept_in_the_small_table_is_read_listed_and_taken_away_from_both(
) {
    let table = small();
    table
        .set("session/1-tree.result", "(old,tree);")
        .expect("set");
    let idb = Idb::new(Mode::Now);
    let roomy = Roomy::new(table.clone(), Some(idb.open(ROOM_FOR_ALL)));
    assert_eq!(
        roomy.get("session/1-tree.result").as_deref(),
        Some("(old,tree);")
    );
    assert!(roomy.holds("session/1-tree.result"));
    assert_eq!(roomy.names(), ["session/1-tree.result"]);
    // the same name in both is listed once, and the larger one's text is the one read
    roomy
        .set("session/1-tree.result", "(new,tree);")
        .expect("kept");
    assert_eq!(roomy.names(), ["session/1-tree.result"]);
    assert_eq!(
        roomy.get("session/1-tree.result").as_deref(),
        Some("(new,tree);")
    );
    roomy.remove("session/1-tree.result");
    assert!(roomy.names().is_empty() && table.names().is_empty() && idb.entries().is_empty());
}

#[test]
fn a_result_is_moved_through_the_two_tables_to_where_a_result_goes_whichever_table_it_is_in() {
    let table = small();
    let idb = Idb::new(Mode::Now);
    let roomy = Roomy::new(table.clone(), Some(idb.open(ROOM_FOR_ALL)));
    // in the larger table: moved there, with nothing written twice
    roomy.set("session/1-a.result", "A").expect("kept");
    roomy
        .rename("session/1-a.result", "session.previous/1-a.result")
        .expect("moved");
    assert_eq!(idb.names(), ["session.previous/1-a.result"]);
    // in the small table (an earlier version): it is moved into the larger one
    table.set("session/1-b.result", "B").expect("set");
    roomy
        .rename("session/1-b.result", "session.previous/1-b.result")
        .expect("moved");
    assert!(table.names().is_empty());
    assert_eq!(
        idb.names(),
        ["session.previous/1-a.result", "session.previous/1-b.result"]
    );
    // a text of the session goes the small table's way, and a name that is not there is nothing to move
    roomy.set("session/1-0.gb", "LOCUS").expect("kept");
    roomy
        .rename("session/1-0.gb", "session.previous/1-0.gb")
        .expect("moved");
    assert_eq!(table.names(), ["session.previous/1-0.gb"]);
    assert_eq!(roomy.rename("session/none.result", "x/none.result"), Ok(()));
}

#[test]
fn a_result_takes_the_room_of_the_sequences_only_where_it_is_in_the_small_table() {
    let table = small();
    table.set("session/1-old.result", "old").expect("set");
    let idb = Idb::new(Mode::Now);
    let roomy = Roomy::new(table.clone(), Some(idb.open(ROOM_FOR_ALL)));
    roomy.set("session/2-new.result", "new").expect("kept");
    assert!(roomy.shares_room("session/1-0.gb"));
    assert!(roomy.shares_room("session/manifest.txt"));
    assert!(
        roomy.shares_room("session/1-old.result"),
        "it is in the small table"
    );
    assert!(
        !roomy.shares_room("session/2-new.result"),
        "it is in the larger one"
    );
    assert!(!roomy.shares_room("session/3-not-yet.result"));
    // with no larger table every result does
    assert!(Roomy::new(table, None).shares_room("session/3-not-yet.result"));
}

#[test]
fn clearing_what_is_kept_apart_takes_the_whole_database_and_leaves_the_small_table_alone() {
    let table = small();
    table.set("appearance.txt", "dark").expect("set");
    let idb = Idb::new(Mode::Now);
    let roomy = Roomy::new(table.clone(), Some(idb.open(ROOM_FOR_ALL)));
    roomy.set("session/1-a.result", "A").expect("kept");
    idb.stored
        .borrow_mut()
        .insert("session/9-another-tab.result".into(), "X".into());
    roomy.clear_apart();
    assert!(idb.entries().is_empty() && roomy.names() == ["appearance.txt"]);
    // a page with no larger table has nothing apart
    Roomy::new(table, None).clear_apart();
}

#[test]
fn with_no_larger_table_the_small_one_keeps_the_results_as_it_did_before() {
    let table = small();
    let roomy = Roomy::new(table.clone(), None);
    roomy.set("session/1-tree.result", "(a,b);").expect("kept");
    assert_eq!(table.names(), ["session/1-tree.result"]);
    assert_eq!(roomy.take_problem(), None);
    roomy.remove("session/1-tree.result");
    assert!(table.names().is_empty());
}

#[test]
fn a_late_failure_of_the_database_is_told_through_the_two_tables_once() {
    let idb = Idb::new(Mode::Refuse);
    let roomy = Roomy::new(small(), Some(idb.open(ROOM_FOR_ALL)));
    assert_eq!(roomy.take_problem(), None);
    // the write is refused as it is handed over; the store has been told that it was kept
    roomy
        .set("session/1-tree.result", "(a,b);")
        .expect("kept in memory");
    assert!(!roomy.holds("session/1-tree.result"));
    assert_eq!(roomy.take_problem().as_deref(), Some("QuotaExceededError"));
    assert_eq!(roomy.take_problem(), None);
}

#[test]
fn which_names_are_results_is_what_the_session_writes() {
    // the ending the session gives a result's file (`session.rs`) is the one that sends it to the larger table
    let table = small();
    let idb = Idb::new(Mode::Now);
    let store = Store::on(Rc::new(Roomy::new(
        table.clone(),
        Some(idb.open(ROOM_FOR_ALL)),
    )));
    store
        .keep_session(&with_a_result("LOCUS", "(a,b);"))
        .expect("kept");
    assert_eq!(idb.names(), ["session/1-tree.result"]);
    assert!(table.names().iter().all(|name| !name.ends_with(".result")));
    for name in [
        "session/12-alignment.result",
        "session.previous/3-tree.result",
        "x.result",
    ] {
        assert!(is_large(name), "{name}");
    }
    for name in [
        "session/1-0.gb",
        "session/manifest.txt",
        "appearance.txt",
        "result",
        "a.results",
    ] {
        assert!(!is_large(name), "{name}");
    }
}

fn store_on(table: &Rc<Memory>, large: Large) -> Store {
    Store::on(Rc::new(Roomy::new(table.clone(), Some(large))))
}

fn with_a_result(doc: &str, text: &str) -> Session {
    Session {
        workspace: "trees".into(),
        active: 0,
        docs: vec![SavedDoc {
            name: "d".into(),
            path: None,
            dirty: true,
            text: Some(doc.into()),
            keep: None,
        }],
        lost: 0,
        results: vec![SavedResult {
            kind: "tree".into(),
            text: Some(text.into()),
            keep: None,
        }],
    }
}

#[test]
fn putting_the_session_aside_moves_the_results_in_the_database_with_the_rest_and_clearing_takes_them_away(
) {
    let (table, idb) = (small(), Idb::new(Mode::Now));
    let store = store_on(&table, idb.open(ROOM_FOR_ALL));
    store
        .keep_session(&with_a_result("LOCUS", "(a,b);"))
        .expect("kept");
    assert_eq!(idb.names(), ["session/1-tree.result"]);
    store.archive_session().expect("put aside");
    // the text itself moved, whole
    assert_eq!(
        idb.entries(),
        [(
            "session.previous/1-tree.result".to_string(),
            "(a,b);".to_string()
        )]
    );
    assert_eq!(
        table.names(),
        ["session.previous/1-0.gb", "session.previous/manifest.txt"]
    );
    assert_eq!(store.load_session(), None, "it is no longer the session");
    store.clear_previous().expect("cleared");
    assert!(idb.names().is_empty() && table.names().is_empty());
}

#[test]
fn a_session_whose_results_are_more_than_half_the_room_can_still_be_put_aside_because_a_move_is_not_a_copy(
) {
    let (table, idb) = (small(), Idb::new(Mode::Now));
    // room for a result once (21 bytes of name and 100 of text) and its move (30 bytes of name), and not for it twice
    let store = store_on(&table, idb.open(200));
    let kept = with_a_result("LOCUS", &"t".repeat(100));
    store.keep_session(&kept).expect("kept");
    store
        .archive_session()
        .expect("put aside, with no room for a second copy");
    assert_eq!(idb.entries()[0].1, "t".repeat(100));
    // and the session that was put aside is not held up the sequences that come after it: a new session is kept beside it
    let next = with_a_result("LOCUS", "(c,d);");
    store
        .keep_session(&next)
        .expect("a new session beside the one that was put aside");
    assert_eq!(store.load_session().expect("a session").results.len(), 1);
}

#[test]
fn a_save_that_fails_because_the_small_table_is_full_leaves_the_results_in_the_database_alone() {
    let (table, idb) = (small(), Idb::new(Mode::Now));
    let store = store_on(&table, idb.open(ROOM_FOR_ALL));
    store
        .keep_session(&with_a_result("LOCUS", "(a,b);"))
        .expect("kept");
    // the small table is now full: a longer sequence does not fit
    let held: usize = table
        .held
        .borrow()
        .iter()
        .map(|(k, v)| k.len() + v.len())
        .sum();
    table.room.set(Some(held + 10));
    let why = store
        .keep_session(&with_a_result(&"ACGT".repeat(500), "(a,b);"))
        .expect_err("no room");
    assert_eq!(why, crate::ui::session::browser::FULL);
    // making room would not have mended it: the result is not taken away for it
    assert_eq!(idb.entries().len(), 1);
    assert_eq!(
        store
            .load_session()
            .expect("the session is as it was")
            .results
            .len(),
        1
    );
    // and a result that is in the small table (no database) is taken away to make room, as it was
    let alone = small();
    let before = Store::on(alone.clone());
    before
        .keep_session(&with_a_result("LOCUS", &"t".repeat(300)))
        .expect("kept");
    let held: usize = alone
        .held
        .borrow()
        .iter()
        .map(|(k, v)| k.len() + v.len())
        .sum();
    alone.room.set(Some(held + 100));
    before
        .keep_session(&with_a_result(&"ACGT".repeat(20), "(a,b);"))
        .expect("kept after making room");
    assert_eq!(
        before.load_session().expect("a session").docs[0]
            .text
            .as_deref()
            .map(str::len),
        Some(80)
    );
}

#[test]
fn a_result_that_an_earlier_version_left_in_the_full_small_table_is_taken_away_to_make_room_for_the_sequences(
) {
    // the small table is full of yesterday's result; the page has a database now
    let table = small();
    Store::on(table.clone())
        .keep_session(&with_a_result("LOCUS", &"t".repeat(2_000)))
        .expect("kept by the earlier version");
    let held: usize = table
        .held
        .borrow()
        .iter()
        .map(|(k, v)| k.len() + v.len())
        .sum();
    table.room.set(Some(held + 100));
    let idb = Idb::new(Mode::Now);
    let store = store_on(&table, idb.open(ROOM_FOR_ALL));
    // a longer sequence does not fit beside it: the result in the small table goes, the sequence is kept
    let longer = with_a_result(&"ACGT".repeat(40), "(a,b);");
    store.keep_session(&longer).expect("kept after making room");
    let back = store.load_session().expect("a session");
    assert_eq!(back.docs[0].text.as_deref().map(str::len), Some(160));
    assert_eq!(back.results.len(), 1);
    // the new result is in the database, and the old one is gone from the small table
    assert!(table.names().iter().all(|name| !name.ends_with(".result")));
    assert_eq!(idb.names().len(), 1);
}

#[test]
fn forgetting_the_sessions_takes_the_database_empty_including_what_another_page_put_there() {
    let (table, idb) = (small(), Idb::new(Mode::Now));
    let store = store_on(&table, idb.open(ROOM_FOR_ALL));
    store
        .keep_session(&with_a_result("LOCUS", "(a,b);"))
        .expect("kept");
    idb.stored
        .borrow_mut()
        .insert("session/9-another-tab.result".into(), "X".into());
    store.forget_sessions().expect("forgotten");
    assert!(idb.entries().is_empty());
    assert!(table.names().is_empty());
}

/// The names the database was asked for, in the order it was, as text: `put X`, `delete X`, `clear`.
fn asked(idb: &Idb) -> Vec<String> {
    idb.ops
        .borrow()
        .iter()
        .map(|op| match op {
            Op::Put(name, _) => format!("put {name}"),
            Op::Delete(name) => format!("delete {name}"),
            Op::Clear => "clear".to_string(),
        })
        .collect()
}

#[test]
fn a_newer_text_that_was_written_after_an_older_one_was_refused_is_not_taken_away_by_a_later_change(
) {
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    let name = "session/1-a.result";
    large.set(name, "older").expect("kept");
    large.set(name, "newer").expect("kept");
    idb.fail("QuotaExceededError");
    // the newer one is written, and then another change comes: the database must still hold the newer one
    idb.finish();
    idb.mode.set(Mode::Now);
    large.set("session/2-b.result", "b").expect("kept");
    assert_eq!(
        idb.stored.borrow().get(name).map(String::as_str),
        Some("newer")
    );
    assert_eq!(
        asked(&idb),
        [
            format!("put {name}"),
            format!("put {name}"),
            "put session/2-b.result".to_string()
        ]
    );
}

#[test]
fn what_the_database_may_hold_is_not_asked_to_be_taken_away_when_the_same_name_is_put_again_or_already_waits_to_go(
) {
    let (x, y, z) = (
        "session/1-x.result",
        "session/2-y.result",
        "session/1-z.result",
    );
    // a refused text, and the same name put again: the put takes the place of what was there
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    large.set(x, "one").expect("kept");
    idb.fail("QuotaExceededError");
    idb.mode.set(Mode::Now);
    large.set(x, "two").expect("kept");
    assert_eq!(asked(&idb), [format!("put {x}"), format!("put {x}")]);
    // a refused text whose name then has a removal waiting behind another change: it is not asked for twice
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    large.set(x, "one").expect("kept");
    large.set(z, "z").expect("kept");
    large.remove(x);
    idb.fail("QuotaExceededError");
    // z is out and the removal of x waits behind it: a change for another name comes now, and finds x to be taken away
    large.set(y, "y").expect("kept");
    idb.mode.set(Mode::Now);
    idb.finish();
    assert_eq!(
        asked(&idb),
        [
            format!("put {x}"),
            format!("put {z}"),
            format!("delete {x}"),
            format!("put {y}")
        ]
    );
}

#[test]
fn what_an_earlier_refusal_left_is_taken_away_before_the_changes_that_wait_and_not_after_them() {
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    large.set("session/1-a.result", "a").expect("kept");
    large.set("session/1-b.result", "b").expect("kept");
    large.set("session/1-c.result", "c").expect("kept");
    idb.fail("QuotaExceededError");
    // b is out, c waits; the next change finds a to be taken away, and it goes before c
    large.set("session/1-d.result", "d").expect("kept");
    while idb.waiting() > 0 {
        idb.finish();
    }
    assert_eq!(
        asked(&idb),
        [
            "put session/1-a.result",
            "put session/1-b.result",
            "delete session/1-a.result",
            "put session/1-c.result",
            "put session/1-d.result"
        ]
    );
}

#[test]
fn a_clear_that_the_database_took_leaves_nothing_to_take_away_and_one_after_a_refused_one_is_not_asked_again(
) {
    // what an earlier refusal left is of no use once everything is taken away
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    large.set("session/1-a.result", "a").expect("kept");
    idb.fail("QuotaExceededError");
    idb.mode.set(Mode::Now);
    large.clear();
    large.set("session/2-b.result", "b").expect("kept");
    assert_eq!(
        asked(&idb),
        ["put session/1-a.result", "clear", "put session/2-b.result"]
    );
    // a clear that is refused is asked again with the next change; a clear that is asked and taken in the meantime is not asked a third time
    let idb = Idb::new(Mode::Later);
    let large = idb.open(ROOM_FOR_ALL);
    large.clear();
    idb.fail("InvalidStateError");
    idb.mode.set(Mode::Now);
    large.clear();
    large.set("session/2-b.result", "b").expect("kept");
    assert_eq!(asked(&idb), ["clear", "clear", "put session/2-b.result"]);
    assert_eq!(idb.names(), ["session/2-b.result"]);
}

#[test]
fn the_room_is_used_up_to_the_last_byte_and_what_is_taken_out_gives_its_room_back() {
    let idb = Idb::new(Mode::Now);
    let name = "session/1-a.result";
    // 100 bytes of room: a name of 18 and a text of 82 fill it exactly
    let large = idb.open(100);
    assert!(large.set(name, &"x".repeat(82)).is_ok(), "fills it exactly");
    assert!(large.set(name, &"x".repeat(83)).is_err(), "a byte over");
    // a rename that comes to the room exactly (the new name is 9 bytes longer, the text 9 shorter)
    let tight = Large::new(idb.clone(), vec![(name.to_string(), "t".repeat(60))], 87);
    assert!(tight.rename(name, "session.previous/1-a.result").is_ok());
    // clearing gives all the room back
    large.clear();
    assert!(large.set(name, &"y".repeat(82)).is_ok());
    // a text that the database refused gives its room back
    let refusing = Idb::new(Mode::Refuse);
    let large = refusing.open(100);
    large
        .set(name, &"x".repeat(60))
        .expect("kept in memory, refused at once");
    assert!(!large.holds(name));
    refusing.mode.set(Mode::Now);
    assert!(
        large.set("session/1-b.result", &"x".repeat(60)).is_ok(),
        "the refused one took no room"
    );
}
