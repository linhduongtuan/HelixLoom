//! The larger table of a page (`session/large.rs`): results are kept in the browser's database behind the small table, which keeps
//! the open sequences as it did.

use super::*;
use crate::seq::newick::NwTree;
use crate::ui::results::{self, STRUCTURE, TREE};
use crate::ui::session::large::stand_in::{Idb, Mode};
use crate::ui::session::large::Op;
use crate::ui::session::large::{Roomy, ROOM};
use crate::ui::session::Memory;
use std::rc::Rc;

/// A page that keeps the open sequences in `table` and the results in `idb`.
fn page_with(table: &Rc<Memory>, idb: &Rc<Idb>) -> App {
    let mut app = App::new(Startup::default());
    let roomy = Roomy::new(table.clone(), Some(idb.open(ROOM)));
    app.keeper = Keeper::with(Some(Store::on(Rc::new(roomy))));
    app.docs.clear();
    app
}

fn tree(text: &str) -> crate::seq::tree::Node {
    NwTree::parse_all(text).expect("a tree")[0].to_node()
}

fn used(table: &Memory) -> usize {
    table
        .held
        .borrow()
        .iter()
        .map(|(k, v)| k.len() + v.len())
        .sum()
}

fn result_names(names: Vec<String>) -> Vec<String> {
    names
        .into_iter()
        .filter(|n| n.ends_with(".result"))
        .collect()
}

#[test]
fn a_result_that_the_small_table_has_no_room_for_is_kept_in_the_database_and_comes_back() {
    let table = Rc::new(Memory::default());
    let idb = Idb::new(Mode::Now);
    let mut app = page_with(&table, &idb);
    app.add_doc(dirty("typed", "ACGTACGT"));
    // the small table has room for the sequence and the list, and not for a structure of 320 thousand letters
    table.room.set(Some(used(&table) + 20_000));
    app.platform.structure.load_antibody_example();
    let structure = app.result_text(STRUCTURE).expect("a structure");
    assert!(structure.len() > 300_000);
    app.flush_session();
    assert_eq!(app.keeper.take_problem(), None);
    let stored = result_names(idb.entries().into_iter().map(|(name, _)| name).collect());
    assert_eq!(stored.len(), 1, "{stored:?}");
    assert!(stored[0].ends_with("-structure.result"));
    assert!(result_names(crate::ui::session::Shelf::names(&*table)).is_empty());
    drop(app);
    // the next visit: a page opened on the same tables
    let mut next = page_with(&table, &idb);
    let offer = next.keeper.offer.clone().expect("an offer");
    assert_eq!(offer.docs.len(), 1);
    assert_eq!(offer.results.len(), 1);
    next.answer_offer(true);
    assert_eq!(next.result_text(STRUCTURE), Some(structure));
    assert_eq!(
        next.status,
        "Reopened 1 sequence from the last session, with the structure"
    );
}

#[test]
fn the_sequences_are_in_the_small_table_when_the_save_is_over_and_a_page_closed_before_the_database_answers_loses_the_result_only(
) {
    let table = Rc::new(Memory::default());
    let idb = Idb::new(Mode::Later);
    let mut app = page_with(&table, &idb);
    app.add_doc(dirty("typed", "ACGTACGT"));
    app.show_tree(tree("(a:1,b:1);"), "shown");
    app.flush_session();
    // the sequence and the list are there, as a page that is closed this moment finds them; the result is on its way
    let names = crate::ui::session::Shelf::names(&*table);
    assert!(names.iter().any(|n| n.ends_with(".gb")), "{names:?}");
    assert!(
        names.iter().any(|n| n.ends_with("manifest.txt")),
        "{names:?}"
    );
    assert_eq!(idb.waiting(), 1);
    assert!(idb.entries().is_empty());
    drop(app);
    let mut next = page_with(&table, &idb);
    let offer = next.keeper.offer.clone().expect("an offer");
    assert_eq!((offer.docs.len(), offer.results.len()), (1, 0));
    next.answer_offer(true);
    assert_eq!(next.status, "Reopened 1 sequence from the last session");
    assert_eq!(next.docs[0].doc.name, "typed");
}

#[test]
fn a_database_that_refuses_a_write_after_the_save_is_said_once_in_red_and_the_result_is_written_again_at_the_next_save(
) {
    let table = Rc::new(Memory::default());
    let idb = Idb::new(Mode::Refuse);
    let mut app = page_with(&table, &idb);
    app.add_doc(dirty("typed", "ACGTACGT"));
    app.show_tree(tree("(a:1,b:1);"), "shown");
    app.flush_session();
    // the save itself went well; the database said no a moment later, and the next frame says so
    let ctx = egui::Context::default();
    app.keep_work(&ctx);
    let said = app.status.clone();
    assert_eq!(
        said,
        "Could not keep what the pages show for next time: the browser's database refused it (QuotaExceededError)."
    );
    assert_eq!(
        crate::ui::toasts::tone_of(&said),
        crate::ui::design::Tone::Bad
    );
    // the sequences are kept all the same
    let offer = Store::on(table.clone()).load_session().expect("a session");
    assert_eq!(offer.docs.len(), 1);
    // nothing more is said while nothing is asked of it
    app.status.clear();
    app.keep_work(&ctx);
    assert_eq!(app.status, "");
    // the work changes and the result is written again, and refused again: it is the same reason, and it is not said again
    app.add_doc(dirty("one more", "CCCC"));
    app.flush_session();
    app.keep_work(&ctx);
    assert_eq!(
        idb.ops
            .borrow()
            .iter()
            .filter(|op| matches!(op, Op::Put(..)))
            .count(),
        2
    );
    assert_eq!(app.status, "");
    // the database allows writes again, and the work changes: the result is written, with nothing said
    idb.mode.set(Mode::Now);
    app.add_doc(dirty("more", "GGGG"));
    app.flush_session();
    app.keep_work(&ctx);
    assert_eq!(app.status, "");
    let stored = result_names(idb.entries().into_iter().map(|(name, _)| name).collect());
    assert_eq!(stored.len(), 1, "{stored:?}");
    assert!(stored[0].ends_with("-tree.result"));
}

#[test]
fn turning_the_remembering_off_takes_the_results_out_of_the_database_too_and_on_again_writes_them()
{
    let table = Rc::new(Memory::default());
    let idb = Idb::new(Mode::Now);
    let mut app = page_with(&table, &idb);
    app.add_doc(dirty("typed", "ACGTACGT"));
    app.show_tree(tree("(a:1,b:1);"), "shown");
    app.flush_session();
    assert_eq!(idb.entries().len(), 1);
    assert_eq!(app.keeper.set_keeping(false).as_deref(), Some(NOT_KEPT_NOW));
    assert!(idb.entries().is_empty());
    assert!(result_names(crate::ui::session::Shelf::names(&*table)).is_empty());
    assert_eq!(app.keeper.set_keeping(true).as_deref(), Some(KEPT_AGAIN));
    app.flush_session();
    assert_eq!(idb.entries().len(), 1);
}

#[test]
fn results_that_the_small_table_kept_before_are_offered_put_back_and_then_moved_to_the_database() {
    let table = Rc::new(Memory::default());
    // an earlier version of the page: no database, the results in the small table
    let mut before = page_on(&table);
    before.add_doc(dirty("typed", "ACGTACGT"));
    before.show_tree(tree("(a:1,b:1);"), "shown");
    before.flush_session();
    let old = result_names(crate::ui::session::Shelf::names(&*table));
    assert_eq!(old.len(), 1, "{old:?}");
    drop(before);
    // this version, with a database
    let idb = Idb::new(Mode::Now);
    let mut app = page_with(&table, &idb);
    let offer = app.keeper.offer.clone().expect("an offer");
    assert_eq!(offer.results.len(), 1);
    app.answer_offer(true);
    assert!(app.tree_win.tree.is_some());
    assert_eq!(
        app.status,
        "Reopened 1 sequence from the last session, with the tree"
    );
    // the next save writes what came back into the database, and takes the old copy out of the small table
    app.add_doc(dirty("more", "GGGG"));
    app.flush_session();
    assert!(result_names(crate::ui::session::Shelf::names(&*table)).is_empty());
    assert_eq!(idb.entries().len(), 1);
    assert!(idb.entries()[0].0.ends_with("-tree.result"));
    assert_eq!(results::called(TREE), Some("the tree"));
}

#[test]
fn a_page_without_a_database_keeps_the_results_in_the_small_table_as_it_did() {
    let table = Rc::new(Memory::default());
    let mut app = App::new(Startup::default());
    let roomy = Roomy::new(table.clone(), None);
    app.keeper = Keeper::with(Some(Store::on(Rc::new(roomy))));
    app.docs.clear();
    app.add_doc(dirty("typed", "ACGTACGT"));
    app.show_tree(tree("(a:1,b:1);"), "shown");
    app.flush_session();
    assert_eq!(
        result_names(crate::ui::session::Shelf::names(&*table)).len(),
        1
    );
    assert_eq!(app.keeper.take_problem(), None);
}

#[test]
fn a_refusal_that_came_between_two_frames_is_known_before_the_next_save_works_out_what_to_write() {
    let table = Rc::new(Memory::default());
    let idb = Idb::new(Mode::Later);
    let mut app = page_with(&table, &idb);
    app.add_doc(dirty("typed", "ACGTACGT"));
    app.show_tree(tree("(a:1,b:1);"), "shown");
    app.flush_session();
    assert_eq!(idb.waiting(), 1);
    // the database refuses while no frame is drawn; the sequence changes and is written on the way out of the page
    idb.fail("QuotaExceededError");
    idb.mode.set(Mode::Now);
    app.add_doc(dirty("more", "GGGG"));
    app.flush_session();
    // the result was not taken to be written: it is written again with this save, and the manifest names it
    let stored = result_names(idb.entries().into_iter().map(|(name, _)| name).collect());
    assert_eq!(stored.len(), 1, "{stored:?}");
    let both = Store::on(Rc::new(Roomy::new(table.clone(), Some(idb.open(ROOM)))));
    let back = both.load_session().expect("a session");
    assert_eq!((back.docs.len(), back.results.len()), (2, 1));
    // and the reason is told, once
    assert!(app
        .keeper
        .take_problem()
        .is_some_and(|said| said.contains("QuotaExceededError")));
    assert_eq!(app.keeper.take_problem(), None);
}

#[test]
fn a_refusal_is_not_told_as_a_failure_to_keep_when_nothing_is_being_kept() {
    let table = Rc::new(Memory::default());
    let idb = Idb::new(Mode::Later);
    let mut app = page_with(&table, &idb);
    app.add_doc(dirty("typed", "ACGTACGT"));
    app.show_tree(tree("(a:1,b:1);"), "shown");
    app.flush_session();
    // the person says not to remember: the removal of what was kept is queued behind the write that is out
    assert_eq!(app.keeper.set_keeping(false).as_deref(), Some(NOT_KEPT_NOW));
    idb.fail("InvalidStateError");
    let ctx = egui::Context::default();
    app.status.clear();
    app.keep_work(&ctx);
    assert_eq!(app.status, "");
    assert_eq!(app.keeper.take_problem(), None);
}
