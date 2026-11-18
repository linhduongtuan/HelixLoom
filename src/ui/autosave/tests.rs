use super::*;
use crate::ui::app::Startup;
use std::sync::atomic::{AtomicU32, Ordering};

static NEXT: AtomicU32 = AtomicU32::new(0);

/// How often, and how long apart, a folder that is found locked is tried again (about two seconds in all).
const LOCK_TRIES: usize = 200;
const LOCK_PAUSE: Duration = Duration::from_millis(10);

/// A settings folder of its own for one test, removed when the test ends.
struct Folder(PathBuf);

impl Folder {
    fn new() -> Folder {
        let path = std::env::temp_dir().join(format!(
            "helixloom-keeper-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }

    /// A copy of the program that keeps its settings here, with nothing open (the example plasmid is gone), for a test in which
    /// no other copy is alive: a folder that is found locked is tried again for a while before it is believed.
    ///
    /// The lock of the session folder is let go when the last holder of the open file lets go, and a program that another test is
    /// starting at that moment (they run side by side) has the file too until it has started: the lock is then refused for a few
    /// tens of milliseconds although no copy of the program has it. (Measured with a small program that locks, drops and locks a
    /// file again while other threads start programs: refused about once in 150 rounds, never for longer than 0.1 s; with no
    /// program started it was never refused in 278,000 rounds.)
    fn app(&self) -> App {
        for _ in 0..LOCK_TRIES {
            let app = self.app_while_another_may_be_alive();
            if !app.keeper.other_copy {
                return app;
            }
            drop(app);
            std::thread::sleep(LOCK_PAUSE);
        }
        self.app_while_another_may_be_alive()
    }

    /// A copy of the program as `app` makes it, without trying again: for a test in which another copy is alive on purpose.
    fn app_while_another_may_be_alive(&self) -> App {
        let mut app = App::new(Startup {
            settings_dir: Some(self.0.clone()),
            ..Startup::default()
        });
        app.docs.clear();
        app
    }

    fn file(&self, name: &str, text: &str) -> PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, text).expect("a file");
        path
    }

    fn names(&self, sub: &str) -> Vec<String> {
        let mut v: Vec<String> = std::fs::read_dir(self.0.join(sub))
            .map(|d| {
                d.flatten()
                    .map(|e| e.file_name().to_string_lossy().into_owned())
                    .collect()
            })
            .unwrap_or_default();
        v.sort();
        v
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn dirty(name: &str, seq: &str) -> Doc {
    let mut d = Doc::new(name, seq.as_bytes().to_vec());
    d.dirty = true;
    d
}

/// Save the session of `app` now, as the keeper does after a change has waited.
fn save_now(app: &mut App) {
    let signature = app.session_signature();
    let session = app.session_now(&signature);
    app.keeper.write(session, signature, Instant::now());
    assert!(app.keeper.problem.is_none(), "{:?}", app.keeper.problem);
}

/// Another copy's keeper must not be alive while a test makes the next: it holds the lock.
fn offer_of(folder: &Folder) -> Option<Session> {
    folder.app().keeper.offer.clone()
}

#[test]
fn the_last_session_is_offered_and_reopened_as_it_was_names_unsaved_changes_tab_and_workspace() {
    let folder = Folder::new();
    let on_disk = folder.file(
        "clean.gb",
        &io::to_genbank(&Doc::new("clean", b"AAAACCCCGGGG".to_vec())),
    );
    let mut app = folder.app();
    assert!(
        app.keeper.offer.is_none(),
        "a first run has nothing to offer"
    );
    let mut edited = Doc::new("edited", b"GGGGTTTT".to_vec());
    edited.path = Some(folder.file(
        "edited.gb",
        &io::to_genbank(&Doc::new("edited", b"GGGG".to_vec())),
    ));
    edited.dirty = true;
    let mut clean = Doc::new("clean", b"AAAACCCCGGGG".to_vec());
    clean.path = Some(on_disk);
    // a name that a GenBank locus would change: it must come back as it was
    for d in [
        dirty("pET28a-GFP construct v2", "ACGTACGTACGT"),
        edited,
        clean,
    ] {
        app.add_doc(d);
    }
    app.tab = Tab::Doc(1);
    app.workspace = Workspace::Alignment;
    save_now(&mut app);
    drop(app);

    let mut next = folder.app();
    assert_eq!(next.keeper.offer.as_ref().map(|s| s.docs.len()), Some(3));
    // the program opened with its example plasmid, which is not wanted beside the sequences that were open
    next.add_doc(crate::sample::demo_plasmid());
    next.docs.rotate_right(1);
    next.answer_offer(true);
    assert!(next.keeper.offer.is_none());
    assert_eq!(next.docs.len(), 3, "the example is gone");
    let back = &next.docs;
    assert_eq!(
        (
            back[0].doc.name.as_str(),
            back[0].doc.dirty,
            &back[0].doc.seq[..]
        ),
        ("pET28a-GFP construct v2", true, &b"ACGTACGTACGT"[..]),
        "the name as it was, not as a locus writes it"
    );
    assert_eq!(
        (back[1].doc.dirty, &back[1].doc.seq[..]),
        (true, &b"GGGGTTTT"[..])
    );
    assert!(back[1].doc.path.is_some());
    assert_eq!(
        (back[2].doc.dirty, &back[2].doc.seq[..]),
        (false, &b"AAAACCCCGGGG"[..])
    );
    assert_eq!(next.tab, Tab::Doc(1));
    assert_eq!(next.workspace, Workspace::Alignment);
}

#[test]
fn work_done_before_the_offer_is_answered_is_kept_and_the_offered_session_is_put_aside_not_lost() {
    let folder = Folder::new();
    let mut first = folder.app();
    first.add_doc(dirty("old one", "AAAAAAAA"));
    first.add_doc(dirty("old two", "CCCCCCCC"));
    save_now(&mut first);
    drop(first);

    let mut next = folder.app();
    assert!(next.keeper.offer.is_some());
    // the person goes on without answering: what they do is kept
    next.add_doc(dirty("new", "GGGGGGGG"));
    next.flush_session();
    assert_eq!(
        folder.names("session.previous").len(),
        3,
        "the old session, whole"
    );
    assert!(folder
        .names("session")
        .contains(&"manifest.txt".to_string()));
    // the offer is still there to be answered, and reopening it brings the old sequences back beside the new
    assert!(next.keeper.offer.is_some());
    next.answer_offer(true);
    let names: Vec<&str> = next.docs.iter().map(|d| d.doc.name.as_str()).collect();
    assert_eq!(names, ["new", "old one", "old two"]);
    // and once the session holds them all, what was put aside is removed
    next.flush_session();
    assert!(folder.names("session.previous").is_empty());
    let kept = next.keeper.saved.clone().expect("saved");
    assert_eq!(kept.docs.len(), 3);
    drop(next);
    assert_eq!(offer_of(&folder).map(|s| s.docs.len()), Some(3));
}

#[test]
fn starting_fresh_removes_what_was_offered_and_keeps_the_work_done_since() {
    let folder = Folder::new();
    let mut first = folder.app();
    first.add_doc(dirty("old", "AAAAAAAA"));
    save_now(&mut first);
    drop(first);

    // not answered, and the person works: the old one is put aside; then they start fresh
    let mut next = folder.app();
    next.add_doc(dirty("new", "GGGGGGGG"));
    next.flush_session();
    next.answer_offer(false);
    assert!(folder.names("session.previous").is_empty());
    drop(next);
    let kept = offer_of(&folder).expect("the new work is kept");
    assert_eq!(kept.docs.len(), 1);
    assert_eq!(kept.docs[0].name, "new");

    // answered at once, with nothing done since: the session itself goes
    let mut again = folder.app();
    again.answer_offer(false);
    drop(again);
    assert!(offer_of(&folder).is_none());
}

#[test]
fn a_change_is_written_when_things_go_quiet_or_when_it_has_waited_long_and_not_before() {
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(dirty("mine", "ACGTACGTAC"));
    let start = Instant::now();
    let at = |ms: u64| start + Duration::from_millis(ms);
    let sig = app.session_signature();
    assert_eq!(app.keeper.due(&sig, at(0)), Due::Wait);
    assert_eq!(app.keeper.due(&sig, at(1_900)), Due::Wait);
    assert_eq!(app.keeper.due(&sig, at(2_000)), Due::Now);
    // typing goes on: each change starts the quiet again...
    let mut busy = app.session_signature();
    for ms in (2_100..14_000).step_by(1_000) {
        busy.docs[0].len += 1;
        assert_eq!(app.keeper.due(&busy, at(ms)), Due::Wait, "{ms}");
    }
    // ...but not for ever
    busy.docs[0].len += 1;
    assert_eq!(app.keeper.due(&busy, at(17_000)), Due::Now);
    // once written, the same state is not written again
    save_now(&mut app);
    assert_eq!(
        app.keeper.due(&app.session_signature(), at(100_000)),
        Due::Never
    );
    // an edit, a change of the circle, a rename, a path and a move to another workspace are each a change
    let mut seen = vec![app.session_signature()];
    app.docs[0].replace_selection(b"AC");
    seen.push(app.session_signature());
    app.docs[0].doc.circular = !app.docs[0].doc.circular;
    seen.push(app.session_signature());
    app.docs[0].doc.name = "renamed".into();
    seen.push(app.session_signature());
    app.docs[0].doc.path = Some(folder.0.join("somewhere.gb"));
    seen.push(app.session_signature());
    app.workspace = Workspace::Trees;
    seen.push(app.session_signature());
    for pair in seen.windows(2) {
        assert_ne!(pair[0], pair[1]);
    }
    // closing every sequence is a change too: the folder is then emptied
    app.docs.clear();
    assert_eq!(
        app.keeper.due(&app.session_signature(), at(200_000)),
        Due::Wait
    );
    assert_eq!(
        app.keeper.due(&app.session_signature(), at(203_000)),
        Due::Now
    );
    save_now(&mut app);
    drop(app);
    assert!(offer_of(&folder).is_none());
}

#[test]
fn a_sequence_that_has_not_changed_is_not_written_again_while_another_is_edited() {
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(dirty("big", "ACGTACGTACGTACGT"));
    app.add_doc(dirty("small", "TTTT"));
    save_now(&mut app);
    assert_eq!(
        folder.names("session"),
        ["1-0.gb", "1-1.gb", "manifest.txt"]
    );
    let big_text = std::fs::read(folder.0.join("session/1-0.gb")).expect("a file");
    // the second is edited: only it is written, and the first keeps its file
    app.docs[1].replace_selection(b"GG");
    save_now(&mut app);
    assert_eq!(
        folder.names("session"),
        ["1-0.gb", "2-1.gb", "manifest.txt"]
    );
    assert_eq!(
        std::fs::read(folder.0.join("session/1-0.gb")).expect("a file"),
        big_text
    );
    drop(app);
    let back = offer_of(&folder).expect("a session");
    assert_eq!(back.docs.len(), 2);
    // what is in the files is the sequences (GenBank writes them in blocks and in lower case: read them back)
    let sequence = |i: usize| {
        let text = back.docs[i].text.as_deref().expect("a text");
        io::parse_any(text, "x").expect("a record")[0].seq.clone()
    };
    assert_eq!(sequence(0), b"ACGTACGTACGTACGT");
    assert_eq!(sequence(1).len(), 6);
}

#[test]
fn a_save_that_fails_is_said_once_and_tried_again_later_and_then_works() {
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(dirty("mine", "ACGTACGT"));
    // a file where the session folder should be: the save cannot be made
    let blocker = folder.file("session", "not a folder");
    let start = Instant::now();
    let sig = app.session_signature();
    let session = app.session_now(&sig);
    app.keeper.write(session.clone(), sig.clone(), start);
    let said = app.keeper.take_problem().expect("a problem");
    assert!(said.starts_with("Could not save the session"), "{said}");
    assert_eq!(
        crate::ui::toasts::tone_of(&said),
        crate::ui::design::Tone::Bad
    );
    // the same failure again is not said again, and nothing is tried before the wait is over
    app.keeper
        .write(session.clone(), sig.clone(), start + Duration::from_secs(1));
    assert_eq!(app.keeper.take_problem(), None);
    assert_eq!(
        app.keeper.due(&sig, start + Duration::from_secs(2)),
        Due::Wait
    );
    assert_eq!(
        app.keeper
            .due(&sig, start + RETRY_AFTER + Duration::from_secs(3)),
        Due::Now
    );
    // the cause goes away: the next try works, and a new failure would be said again
    std::fs::remove_file(blocker).expect("removed");
    app.keeper
        .write(session, sig, start + RETRY_AFTER + Duration::from_secs(3));
    assert!(app.keeper.saved.is_some());
    assert_eq!(app.keeper.said_problem, None);
}

#[test]
fn two_copies_of_the_program_do_not_share_a_session_and_the_second_says_so() {
    let folder = Folder::new();
    let mut first = folder.app();
    first.add_doc(dirty("one", "ACGT"));
    save_now(&mut first);
    // the second copy starts while the first runs: it is offered nothing, keeps nothing, and tells the person
    let mut second = folder.app_while_another_may_be_alive();
    assert!(second.keeper.other_copy);
    assert!(second.keeper.offer.is_none());
    second.add_doc(dirty("two", "TTTT"));
    assert_eq!(
        second.keeper.due(
            &second.session_signature(),
            Instant::now() + Duration::from_secs(100)
        ),
        Due::Never
    );
    second.flush_session();
    assert!(second
        .keeper
        .take_problem()
        .is_some_and(|t| t.contains("Another copy")));
    // and the first one's session is as it was
    assert_eq!(folder.names("session"), ["1-0.gb", "manifest.txt"]);
    drop(first);
    drop(second);
    // when the first has ended, a copy that starts keeps the session
    let third = folder.app();
    assert!(!third.keeper.other_copy);
    assert_eq!(third.keeper.offer.as_ref().map(|s| s.docs.len()), Some(1));
}

#[test]
fn a_sequence_whose_file_is_open_already_is_not_opened_twice_and_unsaved_work_on_it_comes_back_as_a_copy(
) {
    let folder = Folder::new();
    let path = folder.file(
        "shared.gb",
        &io::to_genbank(&Doc::new("shared", b"ACGTACGT".to_vec())),
    );
    let mut app = folder.app();
    let mut clean = Doc::new("shared", b"ACGTACGT".to_vec());
    clean.path = Some(path.clone());
    app.add_doc(clean);
    save_now(&mut app);
    drop(app);
    // the next run has the file open already (it was given on the command line, say)
    let mut next = folder.app();
    let mut open = Doc::new("shared", b"ACGTACGT".to_vec());
    open.path = Some(path.clone());
    next.add_doc(open);
    next.answer_offer(true);
    assert_eq!(next.docs.len(), 1, "not opened twice");
    // and if the saved one had unsaved changes, it is a copy that cannot be saved over the file
    let session = Session {
        workspace: "sequences".into(),
        active: 0,
        docs: vec![SavedDoc {
            name: "shared".into(),
            path: Some(path),
            dirty: true,
            text: Some(io::to_genbank(&Doc::new("shared", b"GGGGCCCC".to_vec()))),
            keep: None,
        }],
        lost: 0,
        results: Vec::new(),
    };
    let done = next.restore_session(&session);
    assert_eq!((done.opened, done.already_open, done.lost), (1, 0, 0));
    let copy = next.docs.last().expect("the copy");
    assert_eq!(copy.doc.name, "shared (recovered)");
    assert!(copy.doc.path.is_none() && copy.doc.dirty);
    assert_eq!(&copy.doc.seq[..], b"GGGGCCCC");
}

#[test]
fn a_sequence_that_cannot_be_read_back_is_counted_and_the_tab_in_front_is_the_right_one() {
    let folder = Folder::new();
    let gone = folder.file(
        "gone.gb",
        &io::to_genbank(&Doc::new("gone", b"AAAAAAAA".to_vec())),
    );
    let mut app = folder.app();
    let mut lost = Doc::new("gone", b"AAAAAAAA".to_vec());
    lost.path = Some(gone.clone());
    for d in [
        lost,
        dirty("second", "CCCCCCCC"),
        dirty("third", "GGGGGGGG"),
    ] {
        app.add_doc(d);
    }
    app.tab = Tab::Doc(1);
    save_now(&mut app);
    drop(app);
    std::fs::remove_file(&gone).expect("removed");
    let mut next = folder.app();
    let offer = next.keeper.offer.clone().expect("an offer");
    let done = next.restore_session(&offer);
    assert_eq!((done.opened, done.lost), (2, 1));
    // the tab that was in front was "second": the first that came back
    assert_eq!(next.tab, Tab::Doc(0));
    assert_eq!(next.docs[0].doc.name, "second");
    // a record that is not in a file of several is lost, not replaced by another
    let many = folder.file(
        "many.gb",
        &(io::to_genbank(&Doc::new("a", b"ACGT".to_vec()))
            + &io::to_genbank(&Doc::new("b", b"TTTT".to_vec()))),
    );
    let saved = |name: &str| SavedDoc {
        name: name.into(),
        path: Some(many.clone()),
        dirty: false,
        text: None,
        keep: None,
    };
    assert_eq!(
        read_saved(&saved("b")).map(|d| d.seq),
        Some(b"TTTT".to_vec())
    );
    assert!(read_saved(&saved("c")).is_none());
}

#[test]
fn the_example_plasmid_and_an_empty_new_sequence_are_not_kept_until_somebody_has_touched_them() {
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(crate::sample::demo_plasmid());
    app.add_doc(Doc::new("Untitled", Vec::new()));
    assert!(app.session_signature().docs.is_empty());
    assert_eq!(
        app.keeper.due(&app.session_signature(), Instant::now()),
        Due::Never
    );
    // the person's own sequence: only that is kept, and it is the one in front
    app.add_doc(dirty("mine", "ACGTACGT"));
    let signature = app.session_signature();
    assert_eq!((signature.docs.len(), signature.active), (1, 0));
    assert_eq!(app.session_now(&signature).docs[0].name, "mine");
    // an example that was edited is theirs now; so is an untitled one that has bases in it
    app.docs[0].replace_selection(b"AC");
    app.docs[1].replace_selection(b"ACGT");
    assert_eq!(app.session_signature().docs.len(), 3);
    // reopening twice does not pile up examples
    save_now(&mut app);
    drop(app);
    for _ in 0..2 {
        let mut next = folder.app();
        next.add_doc(crate::sample::demo_plasmid());
        next.answer_offer(true);
        assert_eq!(
            next.docs
                .iter()
                .filter(|d| d.doc.name == crate::sample::DEMO_NAME)
                .count(),
            1,
            "the edited one that was kept, and no other example"
        );
        next.flush_session();
    }
}

#[test]
fn recent_files_are_remembered_across_runs_and_one_that_has_gone_is_dropped_when_asked_for() {
    let folder = Folder::new();
    let (a, b) = (
        folder.file("a.gb", &io::to_genbank(&Doc::new("a", b"AC".to_vec()))),
        folder.file("b.gb", &io::to_genbank(&Doc::new("b", b"GT".to_vec()))),
    );
    let mut app = folder.app();
    app.keeper.remember(&a);
    app.keeper.remember(&b);
    app.keeper.remember(&a);
    assert_eq!(app.keeper.recent, [a.clone(), b.clone()]);
    drop(app);
    let mut next = folder.app();
    assert_eq!(next.keeper.recent, [a.clone(), b.clone()]);
    let labels: Vec<String> = next.recent_items().into_iter().map(|(l, _)| l).collect();
    assert!(labels[0].starts_with("a.gb"), "{labels:?}");
    // opening one makes a tab and puts it first
    let before = next.docs.len();
    next.open_recent(1);
    assert_eq!(next.docs.len(), before + 1);
    assert_eq!(next.keeper.recent[0], b);
    // one that has gone is said so and is not offered again
    std::fs::remove_file(&a).expect("removed");
    let index = next
        .keeper
        .recent
        .iter()
        .position(|p| *p == a)
        .expect("listed");
    next.open_recent(index);
    assert!(!next.keeper.recent.contains(&a));
    assert!(next.status.contains("not there"), "{}", next.status);
    next.keeper.forget_recent();
    drop(next);
    assert!(folder.app().keeper.recent.is_empty());
}

#[test]
fn a_recent_file_is_labelled_with_its_name_and_its_folder() {
    assert_eq!(
        recent_label(Path::new("/a/plasmids/pUC19.gb")),
        "pUC19.gb   —   plasmids"
    );
    assert_eq!(recent_label(Path::new("pUC19.gb")), "pUC19.gb");
}

#[test]
fn without_a_settings_folder_nothing_is_offered_or_kept() {
    let mut app = App::new(Startup::default());
    assert!(app.keeper.offer.is_none());
    assert_eq!(
        app.keeper.due(
            &app.session_signature(),
            Instant::now() + Duration::from_secs(60)
        ),
        Due::Never
    );
    app.keeper.remember(Path::new("/no/such/file.gb"));
    assert!(app.keeper.recent.is_empty());
    assert!(Keeper::window_place(None).is_none());
}

/// Run one frame of `ctx` in which the window is at `place` (outer corner, inner size), on a screen of `screen` if that is known.
fn look_at_window(
    app: &mut App,
    ctx: &egui::Context,
    place: (f32, f32, f32, f32),
    maximized: bool,
    screen: Option<(f32, f32)>,
    now: Instant,
) -> egui::FullOutput {
    let info = egui::ViewportInfo {
        outer_rect: Some(egui::Rect::from_min_size(
            egui::pos2(place.0, place.1),
            egui::vec2(place.2, place.3 + 28.0),
        )),
        inner_rect: Some(egui::Rect::from_min_size(
            egui::pos2(place.0, place.1 + 28.0),
            egui::vec2(place.2, place.3),
        )),
        maximized: Some(maximized),
        monitor_size: screen.map(|(w, h)| egui::vec2(w, h)),
        ..Default::default()
    };
    let input = egui::RawInput {
        viewports: std::iter::once((egui::ViewportId::ROOT, info)).collect(),
        ..crate::ui::headless::frame_input(0)
    };
    let mut output = ctx.run_ui(input, |_| app.keep_window(ctx, now));
    // egui insists that what a frame hands over (here, the textures it made) is dealt with before it is dropped
    output.textures_delta.clear();
    output
}

#[test]
fn the_place_of_the_window_is_written_a_moment_after_it_stops_moving_and_not_while_it_is_maximised()
{
    let folder = Folder::new();
    let mut app = folder.app();
    let ctx = egui::Context::default();
    let start = Instant::now();
    let store = Store::open(Some(&folder.0)).expect("a store");
    assert_eq!(store.window(), None);
    let look = |app: &mut App, place, maximized, now| {
        look_at_window(app, &ctx, place, maximized, None, now)
    };

    look(&mut app, (100.0, 50.0, 1000.0, 700.0), false, start);
    // not yet: the window may still be moving
    assert_eq!(store.window(), None);
    look(
        &mut app,
        (100.0, 50.0, 1000.0, 700.0),
        false,
        start + WINDOW_SAVE_AFTER,
    );
    assert_eq!(
        store.window(),
        Some(WindowPlace {
            x: 100.0,
            y: 50.0,
            width: 1000.0,
            height: 700.0
        })
    );
    // moved, then moved again before the time was up: the place that is kept is the last one
    let later = start + Duration::from_secs(10);
    look(&mut app, (300.0, 80.0, 1000.0, 700.0), false, later);
    look(
        &mut app,
        (320.0, 90.0, 1100.0, 720.0),
        false,
        later + Duration::from_secs(1),
    );
    assert_eq!(store.window().map(|w| w.x), Some(100.0), "still waiting");
    look(
        &mut app,
        (320.0, 90.0, 1100.0, 720.0),
        false,
        later + Duration::from_secs(3),
    );
    assert_eq!(
        store.window(),
        Some(WindowPlace {
            x: 320.0,
            y: 90.0,
            width: 1100.0,
            height: 720.0
        })
    );
    // a maximised window is not the place to come back to
    let last = later + Duration::from_secs(60);
    look(&mut app, (0.0, 0.0, 1920.0, 1080.0), true, last);
    look(
        &mut app,
        (0.0, 0.0, 1920.0, 1080.0),
        true,
        last + Duration::from_secs(10),
    );
    assert_eq!(store.window().map(|w| w.width), Some(1100.0));
    // and the next run starts there
    assert_eq!(
        Keeper::window_place(Some(&folder.0)).map(|w| w.x),
        Some(320.0)
    );
}

#[test]
fn the_place_of_a_zoomed_window_is_kept_in_the_windows_own_units_so_that_a_start_does_not_change_it(
) {
    // egui reports the window in points, which are the window's units divided by the zoom: a window of 1280 by 800 that is zoomed to
    // 125 per cent is reported as 1024 by 640. What is kept has to be 1280 by 800, or every start would make it smaller (and at 80
    // per cent larger).
    for (zoom, reported, kept) in [
        (
            1.0,
            (100.0, 50.0, 1280.0, 800.0),
            (100.0, 50.0, 1280.0, 800.0),
        ),
        (
            1.25,
            (80.0, 40.0, 1024.0, 640.0),
            (100.0, 50.0, 1280.0, 800.0),
        ),
        (
            0.8,
            (125.0, 62.5, 1600.0, 1000.0),
            (100.0, 50.0, 1280.0, 800.0),
        ),
        (
            2.0,
            (50.0, 25.0, 640.0, 400.0),
            (100.0, 50.0, 1280.0, 800.0),
        ),
    ] {
        let folder = Folder::new();
        let mut app = folder.app();
        let ctx = egui::Context::default();
        let start = Instant::now();
        // the zoom takes hold in the frame after it is set
        ctx.set_zoom_factor(zoom);
        look_at_window(&mut app, &ctx, reported, false, None, start);
        look_at_window(&mut app, &ctx, reported, false, None, start);
        assert_eq!(ctx.zoom_factor(), zoom);
        look_at_window(
            &mut app,
            &ctx,
            reported,
            false,
            None,
            start + WINDOW_SAVE_AFTER,
        );
        let store = Store::open(Some(&folder.0)).expect("a store");
        assert_eq!(
            store.window(),
            Some(WindowPlace {
                x: kept.0,
                y: kept.1,
                width: kept.2,
                height: kept.3
            }),
            "zoom {zoom}"
        );
    }
}

#[test]
fn a_window_that_starts_where_no_screen_shows_it_is_moved_back_once_and_one_that_is_seen_is_left() {
    let moved = |output: &egui::FullOutput| {
        output
            .viewport_output
            .get(&egui::ViewportId::ROOT)
            .is_some_and(|v| {
                v.commands
                    .iter()
                    .any(|c| matches!(c, egui::ViewportCommand::OuterPosition(_)))
            })
    };
    let folder = Folder::new();
    let ctx = egui::Context::default();
    let now = Instant::now();
    // left on a monitor that has gone: nothing of it is on a 1920 by 1080 screen
    let mut app = folder.app();
    let output = look_at_window(
        &mut app,
        &ctx,
        (5000.0, 3000.0, 1000.0, 700.0),
        false,
        Some((1920.0, 1080.0)),
        now,
    );
    assert!(moved(&output));
    // it is asked once
    let output = look_at_window(
        &mut app,
        &ctx,
        (5000.0, 3000.0, 1000.0, 700.0),
        false,
        Some((1920.0, 1080.0)),
        now,
    );
    assert!(!moved(&output));
    // partly off the screen but plenty of it shows: left where it is
    let mut app = folder.app();
    let output = look_at_window(
        &mut app,
        &ctx,
        (1500.0, 200.0, 1000.0, 700.0),
        false,
        Some((1920.0, 1080.0)),
        now,
    );
    assert!(!moved(&output));
    // the size of the screen is not known yet: nothing is decided
    let mut app = folder.app();
    let output = look_at_window(
        &mut app,
        &ctx,
        (5000.0, 3000.0, 1000.0, 700.0),
        false,
        None,
        now,
    );
    assert!(!moved(&output));
    assert!(!app.keeper.on_screen_checked);
}

/// One whole frame, in which the window is asked to close if `closing`.
fn frame(app: &mut App, ctx: &egui::Context, closing: bool) {
    let info = egui::ViewportInfo {
        events: if closing {
            vec![egui::ViewportEvent::Close]
        } else {
            Vec::new()
        },
        ..Default::default()
    };
    let input = egui::RawInput {
        viewports: std::iter::once((egui::ViewportId::ROOT, info)).collect(),
        ..crate::ui::headless::frame_input(0)
    };
    let mut output = ctx.run_ui(input, |ui| app.draw(ui));
    output.textures_delta.clear();
}

#[test]
fn work_that_is_waiting_to_be_written_is_written_when_the_window_is_closed_and_not_before() {
    let folder = Folder::new();
    let mut app = folder.app();
    let ctx = egui::Context::default();
    app.add_doc(dirty("typed", "ACGTTGCA"));
    // an ordinary frame: the change has only just been seen, so it waits
    frame(&mut app, &ctx, false);
    assert!(folder.names("session").is_empty(), "written too soon");
    // the window is closed: it is written at once
    frame(&mut app, &ctx, true);
    drop(app);
    let offer = offer_of(&folder).expect("the session was written");
    assert_eq!(offer.docs.len(), 1);
    assert_eq!(offer.docs[0].name, "typed");
    assert!(offer.docs[0].dirty);
}

#[test]
fn the_program_ending_writes_what_is_waiting() {
    use eframe::App as _;
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(dirty("typed", "ACGTTGCA"));
    // Cmd+Q on a Mac ends the program without a close request: the exit hook writes
    app.on_exit();
    drop(app);
    assert_eq!(offer_of(&folder).map(|s| s.docs.len()), Some(1));
}

#[test]
fn the_offer_of_the_last_session_is_drawn_over_the_window_and_stays_until_it_is_answered() {
    let folder = Folder::new();
    let mut app = folder.app();
    app.add_doc(dirty("waiting", "ACGT"));
    app.add_doc(dirty("another", "TTTT"));
    save_now(&mut app);
    drop(app);
    let mut next = folder.app();
    let ctx = egui::Context::default();
    for _ in 0..3 {
        frame(&mut next, &ctx, false);
    }
    assert!(next.keeper.offer.is_some(), "an unanswered offer stays");
    // forgetting unsaved work is asked about once more (the second question is drawn too)
    next.keeper.confirm_fresh = true;
    for _ in 0..3 {
        frame(&mut next, &ctx, false);
    }
    assert!(next.keeper.offer.is_some());
    next.keeper.confirm_fresh = false;
    // nothing was lost by looking: the session is still what it was
    drop(next);
    assert_eq!(offer_of(&folder).map(|s| s.docs.len()), Some(2));
}

/// Drag from `from` to `to` with the primary button, in frames of their own as a hand would.
fn drag(screen: &mut crate::ui::testing::Screen, app: &mut App, from: egui::Pos2, to: egui::Pos2) {
    let button = |pos, pressed| egui::Event::PointerButton {
        pos,
        button: egui::PointerButton::Primary,
        pressed,
        modifiers: egui::Modifiers::NONE,
    };
    screen.frame(app, vec![egui::Event::PointerMoved(from)]);
    screen.frame(app, vec![button(from, true)]);
    let step = (to - from) / 4.0;
    for k in 1..=4 {
        screen.frame(app, vec![egui::Event::PointerMoved(from + step * k as f32)]);
    }
    screen.frame(app, vec![button(to, false)]);
    screen.frame(app, vec![]);
}

/// Where the tool window `title` is, as egui has it.
fn corner(ctx: &egui::Context, title: &str) -> egui::Pos2 {
    ctx.memory(|m| m.area_rect(crate::ui::design::window_id(title)))
        .unwrap_or_else(|| panic!("{title} was not drawn"))
        .min
}

/// Let the keeper look at where the tool windows are, at the moment `now` (egui keeps where they are between frames).
fn note_places(app: &mut App, ctx: &egui::Context, now: Instant) {
    app.keep_places(ctx, now);
}

#[test]
fn a_tool_window_that_was_dragged_is_remembered_and_the_next_run_opens_it_there() {
    let folder = Folder::new();
    let mut app = folder.app();
    app.find_win.open = true;
    app.motif_win.open = true;
    let mut screen = crate::ui::testing::Screen::new();
    screen.settle(&mut app);
    let ctx = screen.ctx().clone();
    let (find_start, motif_start) = (corner(&ctx, "Find"), corner(&ctx, "Motif Finder"));
    // the title bar of Find, moved by a hand
    let grab = find_start + egui::vec2(30.0, 12.0);
    drag(&mut screen, &mut app, grab, grab + egui::vec2(120.0, 70.0));
    let moved = corner(&ctx, "Find");
    assert!(
        (moved - find_start - egui::vec2(120.0, 70.0)).length() < 3.0,
        "{find_start:?} to {moved:?}"
    );
    assert_eq!(
        corner(&ctx, "Motif Finder"),
        motif_start,
        "the other window did not move"
    );
    // a moment later the place is written, and not before; a window that was not moved is not kept
    let store = Store::open(Some(&folder.0)).expect("a store");
    let start = Instant::now();
    note_places(&mut app, &ctx, start);
    assert!(store.places().is_empty(), "the window may still be moving");
    note_places(&mut app, &ctx, start + WINDOW_SAVE_AFTER);
    let kept = store.places();
    let (x, y) = kept.get("Find").expect("Find was kept");
    assert!(
        (x - moved.x).abs() < 1.0 && (y - moved.y).abs() < 1.0,
        "{x} {y} against {moved:?}"
    );
    assert_eq!(
        kept.len(),
        1,
        "Motif Finder was left where it opens: {:?}",
        kept.to_text()
    );
    drop(app);
    // the next run opens Find where it was left, and the other where it always opens
    let mut next = folder.app();
    next.find_win.open = true;
    next.motif_win.open = true;
    let mut again = crate::ui::testing::Screen::new();
    again.settle(&mut next);
    let ctx = again.ctx().clone();
    assert!(
        (corner(&ctx, "Find") - moved).length() < 1.0,
        "{:?}",
        corner(&ctx, "Find")
    );
    assert_eq!(corner(&ctx, "Motif Finder"), motif_start);
    // put back at its own place, it is forgotten again
    let there = corner(&ctx, "Find");
    let grab = there + egui::vec2(30.0, 12.0);
    drag(&mut again, &mut next, grab, grab - egui::vec2(120.0, 70.0));
    let start = Instant::now();
    note_places(&mut next, &ctx, start);
    note_places(&mut next, &ctx, start + WINDOW_SAVE_AFTER);
    assert!(
        Store::open(Some(&folder.0))
            .expect("a store")
            .places()
            .is_empty(),
        "back where it opens"
    );
}

#[test]
fn a_place_that_is_far_off_the_screen_is_brought_back_and_without_a_settings_folder_nothing_is_kept(
) {
    let folder = Folder::new();
    let store = Store::open(Some(&folder.0)).expect("a store");
    let mut far = crate::ui::places::Places::default();
    far.set("Find", 15_000.0, 15_000.0);
    store.save_places(&far).expect("saved");
    let mut app = folder.app();
    app.find_win.open = true;
    let mut screen = crate::ui::testing::Screen::new();
    screen.settle(&mut app);
    let rect = screen
        .ctx()
        .memory(|m| m.area_rect(crate::ui::design::window_id("Find")))
        .expect("drawn");
    let screen_rect = screen.ctx().content_rect();
    assert!(
        screen_rect.contains(rect.min) && rect.min.x < screen_rect.right() - 50.0,
        "{rect:?} is on the screen {screen_rect:?}"
    );
    // no settings folder: the places are neither read nor written
    let mut bare = App::new(Startup::default());
    bare.find_win.open = true;
    let mut bare_screen = crate::ui::testing::Screen::new();
    bare_screen.settle(&mut bare);
    let ctx = bare_screen.ctx().clone();
    let grab = corner(&ctx, "Find") + egui::vec2(30.0, 12.0);
    drag(
        &mut bare_screen,
        &mut bare,
        grab,
        grab + egui::vec2(60.0, 60.0),
    );
    let start = Instant::now();
    note_places(&mut bare, &ctx, start);
    note_places(&mut bare, &ctx, start + WINDOW_SAVE_AFTER);
    assert!(bare.keeper.places_for_windows().is_empty());
}

/// A page of a browser: an app that keeps things in `table`, which stands where the browser's storage for the site is, with
/// nothing open.
fn page_on(table: &std::rc::Rc<crate::ui::session::Memory>) -> App {
    let mut app = App::new(Startup::default());
    app.keeper = Keeper::with(Some(Store::on(table.clone())));
    app.docs.clear();
    app
}

#[test]
fn in_a_browser_the_open_sequences_are_offered_again_at_the_next_visit_and_reopened_as_they_were() {
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let mut app = page_on(&table);
    assert!(app.keeper.has_settings() && app.keeper.kept_in_a_browser());
    assert!(
        app.keeper.offer.is_none(),
        "a first visit has nothing to offer"
    );
    // a page has no files: a sequence that was dropped on it or pasted has no path, and is kept as text
    app.add_doc(dirty("pasted in the page", "ACGTACGTACGT"));
    app.add_doc(dirty("a second one", "GGGGCCCC"));
    app.tab = Tab::Doc(1);
    app.workspace = Workspace::Alignment;
    save_now(&mut app);
    assert!(
        crate::ui::session::Shelf::names(table.as_ref())
            .iter()
            .any(|name| name == "session/manifest.txt"),
        "kept in the table"
    );
    drop(app);

    // the next visit
    let mut next = page_on(&table);
    let offer = next
        .keeper
        .offer
        .clone()
        .expect("the last session is offered");
    assert_eq!(offer.docs.len(), 2);
    assert!(offer
        .docs
        .iter()
        .all(|d| d.text.is_some() && d.path.is_none()));
    next.answer_offer(true);
    let names: Vec<&str> = next.docs.iter().map(|d| d.doc.name.as_str()).collect();
    assert_eq!(names, ["pasted in the page", "a second one"]);
    assert_eq!(next.docs[0].doc.seq, b"ACGTACGTACGT");
    assert_eq!(next.workspace, Workspace::Alignment);
    assert_eq!(next.tab, Tab::Doc(1));
}

#[test]
fn in_a_browser_the_appearance_is_kept_for_the_next_visit_and_a_full_storage_is_told_once() {
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let mut app = page_on(&table);
    let chosen = Appearance::default().with_spec("scale=125");
    let now = Instant::now();
    app.keeper
        .appearance_starts_as(ThemeMode::default(), Appearance::default());
    assert_eq!(
        app.keeper.keep_appearance(ThemeMode::Dark, &chosen, now),
        None
    );
    assert_eq!(
        page_on(&table).keeper.saved_appearance(),
        Some((ThemeMode::Dark, chosen))
    );
    // no room for another text: the keeper says why, once, and tries again later
    table.room.set(Some(0));
    let other = Appearance::default().with_spec("scale=150");
    assert!(app
        .keeper
        .keep_appearance(ThemeMode::Dark, &other, now)
        .is_some());
    assert!(app.keeper.appearance_failed);
    let told = app.keeper.take_problem().expect("told");
    assert!(told.contains("browser's storage"), "{told}");
    assert_eq!(app.keeper.take_problem(), None);
}

#[test]
fn a_page_opened_later_takes_the_session_over_and_the_earlier_page_stops_keeping_and_says_so_once()
{
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let look = |name: &str| crate::ui::session::Shelf::get(table.as_ref(), name);
    let mut first = page_on(&table);
    first.add_doc(dirty("made in the first page", "ACGTACGT"));
    save_now(&mut first);
    let list = look("session/manifest.txt").expect("kept");

    // the same site in another tab: it is offered what the first page kept
    let mut second = page_on(&table);
    let offered = second.keeper.offer.clone().expect("offered");
    assert_eq!(offered.docs[0].name, "made in the first page");

    // the first page goes on working: it finds that the session is no longer its own, writes nothing, and says why
    first.add_doc(dirty("made after the other page opened", "GGGG"));
    let signature = first.session_signature();
    let session = first.session_now(&signature);
    first.keeper.write(session, signature, Instant::now());
    assert_eq!(
        look("session/manifest.txt"),
        Some(list.clone()),
        "nothing was written"
    );
    assert_eq!(first.keeper.take_problem().as_deref(), Some(TAKEN_OVER));
    // and it does not try again, or say it again
    let signature = first.session_signature();
    assert!(!first.keeper.pending(&signature));
    assert_eq!(first.keeper.due(&signature, Instant::now()), Due::Never);
    first.flush_session();
    assert_eq!(first.keeper.take_problem(), None);

    // the second page works on without answering: what it was offered is put aside, and its own work is the session
    second.add_doc(dirty("made in the second page", "TTTT"));
    save_now(&mut second);
    assert!(look("session.previous/manifest.txt").is_some());
    let mine = look("session/manifest.txt").expect("kept");
    assert!(mine.contains("made in the second page"), "{mine}");
    // nothing the first page does now touches it
    first.flush_session();
    assert_eq!(look("session/manifest.txt"), Some(mine));
}

#[test]
fn a_page_that_was_taken_over_does_not_remove_the_other_pages_session_by_starting_fresh() {
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let look = |name: &str| crate::ui::session::Shelf::get(table.as_ref(), name);
    let mut before = page_on(&table);
    before.add_doc(dirty("kept from an earlier visit", "ACGT"));
    save_now(&mut before);
    drop(before);
    // a page with that session on offer, not answered yet
    let mut first = page_on(&table);
    assert!(first.keeper.offer.is_some());
    // another page opens, reopens the session and works on
    let mut second = page_on(&table);
    second.answer_offer(true);
    second.add_doc(dirty("more", "GG"));
    save_now(&mut second);
    let theirs = look("session/manifest.txt").expect("kept");
    // "Start fresh" in the first page would remove the session it was offered: it is the second page's now, and stays
    first.answer_offer(false);
    assert_eq!(look("session/manifest.txt"), Some(theirs));
    assert_eq!(first.keeper.take_problem().as_deref(), Some(TAKEN_OVER));
}

#[test]
fn only_a_keeper_on_a_browsers_table_says_that_it_keeps_things_in_a_browser() {
    let folder = Folder::new();
    let app = folder.app();
    assert!(app.keeper.has_settings());
    assert!(!app.keeper.kept_in_a_browser());
    let none = App::new(Startup::default());
    assert!(!none.keeper.has_settings() && !none.keeper.kept_in_a_browser());
}

/// Whether `table` holds anything of a session: the one that is kept, one that was put aside, or one of another version
/// that could not be read.
fn holds_a_session(table: &crate::ui::session::Memory) -> bool {
    crate::ui::session::Shelf::names(table).iter().any(|name| {
        ["session/", "session.previous/", "session.unreadable/"]
            .iter()
            .any(|folder| name.starts_with(folder))
    })
}

/// Whether the site of `table` is set to remember the work.
fn remembers(table: &std::rc::Rc<crate::ui::session::Memory>) -> bool {
    Store::on(table.clone()).sessions_wanted()
}

#[test]
fn remembering_turned_off_removes_what_was_kept_and_keeps_nothing_until_it_is_turned_on_again() {
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let mut app = page_on(&table);
    assert_eq!(
        app.keeper.keeps_work(),
        Some(true),
        "the work is remembered until the person says not to"
    );
    assert_eq!(app.keeper.set_keeping(true), None, "it is on already");
    app.add_doc(dirty("a sequence nobody else should read", "ACGTACGT"));
    save_now(&mut app);
    assert!(holds_a_session(&table));

    assert_eq!(app.keeper.set_keeping(false).as_deref(), Some(NOT_KEPT_NOW));
    assert_eq!(app.keeper.keeps_work(), Some(false));
    assert!(!remembers(&table));
    assert!(!holds_a_session(&table), "what was kept is removed");
    assert_eq!(app.keeper.set_keeping(false), None, "it is off already");
    // nothing is written while it is off, whatever is done and however the page is left
    app.add_doc(dirty("another", "GGGG"));
    let signature = app.session_signature();
    assert!(!app.keeper.pending(&signature));
    assert_eq!(app.keeper.due(&signature, Instant::now()), Due::Never);
    app.flush_session();
    assert!(!holds_a_session(&table));

    // on again in the same page, with nothing changed since: what is open is written, though it was written before
    assert_eq!(app.keeper.set_keeping(true).as_deref(), Some(KEPT_AGAIN));
    assert_eq!(app.keeper.keeps_work(), Some(true));
    assert!(remembers(&table));
    app.flush_session();
    assert!(app.keeper.problem.is_none(), "{:?}", app.keeper.problem);
    let kept = Store::on(table.clone()).load_session().expect("a session");
    let names: Vec<&str> = kept.docs.iter().map(|d| d.name.as_str()).collect();
    assert_eq!(names, ["a sequence nobody else should read", "another"]);
    assert!(
        kept.docs.iter().all(|d| d.text.is_some()),
        "with their texts"
    );

    // off and on once more with nothing done in between: the same work is written again, though this page wrote it before
    app.keeper.set_keeping(false);
    assert!(!holds_a_session(&table));
    app.keeper.set_keeping(true);
    app.flush_session();
    let again = Store::on(table.clone()).load_session().expect("a session");
    assert_eq!(again.docs.len(), 2);
    assert!(again.docs.iter().all(|d| d.text.is_some()));

    // off again, and the page is left: the next visit is offered nothing, the switch is as it was left, nothing is kept
    app.keeper.set_keeping(false);
    drop(app);
    let mut next = page_on(&table);
    assert!(next.keeper.offer.is_none());
    assert_eq!(next.keeper.keeps_work(), Some(false));
    assert_eq!(next.keeper.take_problem(), None);
    next.add_doc(dirty("typed at the next visit", "TTTT"));
    next.flush_session();
    assert!(!holds_a_session(&table));

    // turned on at that visit: kept from then on, and offered at the visit after
    assert_eq!(next.keeper.set_keeping(true).as_deref(), Some(KEPT_AGAIN));
    next.flush_session();
    assert!(holds_a_session(&table));
    drop(next);
    let later = page_on(&table);
    assert_eq!(later.keeper.keeps_work(), Some(true));
    let offer = later.keeper.offer.clone().expect("offered again");
    assert_eq!(offer.docs[0].name, "typed at the next visit");
}

#[test]
fn a_session_on_offer_one_put_aside_and_one_that_could_not_be_read_go_too() {
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let mut before = page_on(&table);
    before.add_doc(dirty("kept from an earlier visit", "ACGT"));
    save_now(&mut before);
    drop(before);
    // on offer, not answered: turning it off answers it, and what was offered is gone
    let mut app = page_on(&table);
    assert!(app.keeper.offer.is_some());
    app.keeper.confirm_fresh = true;
    app.keeper.set_keeping(false);
    assert!(app.keeper.offer.is_none() && !app.keeper.confirm_fresh);
    assert!(!holds_a_session(&table));

    // on offer, and the person worked on, so that the offered session was put aside
    app.keeper.set_keeping(true);
    app.add_doc(dirty("first", "ACGT"));
    app.flush_session();
    drop(app);
    let mut app = page_on(&table);
    assert!(app.keeper.offer.is_some());
    app.add_doc(dirty("made without answering", "GGGG"));
    save_now(&mut app);
    // and a session of another version of the program, which was moved out of the way when it could not be read
    crate::ui::session::Shelf::set(table.as_ref(), "session.unreadable/manifest.txt", "older")
        .expect("set");
    crate::ui::session::Shelf::set(table.as_ref(), "session.unreadable/0-0.gb", "ACGT")
        .expect("set");
    let names = crate::ui::session::Shelf::names(table.as_ref());
    for folder in ["session/", "session.previous/", "session.unreadable/"] {
        assert!(
            names.iter().any(|n| n.starts_with(folder)),
            "{folder} in {names:?}"
        );
    }
    app.keeper.set_keeping(false);
    let names = crate::ui::session::Shelf::names(table.as_ref());
    assert!(
        !names.iter().any(|n| n.starts_with("session/")
            || n.starts_with("session.previous/")
            || n.starts_with("session.unreadable/")),
        "{names:?}"
    );
    assert!(!app.keeper.archived);
}

#[test]
fn remembering_turned_off_in_one_tab_stops_the_tab_that_keeps_the_session_and_it_says_so_once() {
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let mut before = page_on(&table);
    before.add_doc(dirty("kept from an earlier visit", "ACGT"));
    save_now(&mut before);
    drop(before);
    let mut first = page_on(&table);
    // the tab that was opened last keeps the session; it has the offer of the earlier visit on its screen
    let mut second = page_on(&table);
    assert!(second.keeper.offer.is_some());
    // the person turns the remembering off in the other tab: what was kept is removed
    assert_eq!(
        first.keeper.set_keeping(false).as_deref(),
        Some(NOT_KEPT_NOW)
    );
    assert!(!holds_a_session(&table));
    assert_eq!(
        first.keeper.take_problem(),
        None,
        "that tab was told already"
    );
    // the second tab's switch shows it at once, though it has not tried to write yet
    assert_eq!(second.keeper.keeps_work(), Some(false));
    // and when it has something to write it writes nothing, takes its offer away, and says why
    second.add_doc(dirty("made after that", "GGGG"));
    let signature = second.session_signature();
    let session = second.session_now(&signature);
    second.keeper.write(session, signature, Instant::now());
    assert!(!holds_a_session(&table));
    assert!(second.keeper.offer.is_none());
    assert_eq!(
        second.keeper.take_problem().as_deref(),
        Some(TURNED_OFF_ELSEWHERE)
    );
    // it has stopped: it does not look again, or say it again
    let signature = second.session_signature();
    assert!(!second.keeper.pending(&signature));
    assert_eq!(second.keeper.due(&signature, Instant::now()), Due::Never);
    second.flush_session();
    assert_eq!(second.keeper.take_problem(), None);
    assert!(!holds_a_session(&table));

    // turned on again in the second tab: it keeps the session, and the first, which turned it off, keeps nothing
    assert_eq!(second.keeper.set_keeping(true).as_deref(), Some(KEPT_AGAIN));
    second.flush_session();
    let kept = Store::on(table.clone()).load_session().expect("a session");
    assert_eq!(kept.docs[0].name, "made after that");
    first.add_doc(dirty("made in the tab that turned it off", "TT"));
    first.flush_session();
    let kept = Store::on(table.clone()).load_session().expect("a session");
    assert_eq!(kept.docs.len(), 1);
}

#[test]
fn a_page_that_is_to_remember_nothing_takes_the_session_from_nobody_and_removes_one_it_finds() {
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let mut app = page_on(&table);
    app.add_doc(dirty("written by a page that had not heard yet", "ACGT"));
    save_now(&mut app);
    drop(app);
    // the mark alone: as when another page was told to remember nothing just as this one wrote
    Store::on(table.clone())
        .want_sessions(false)
        .expect("the mark");
    assert!(holds_a_session(&table));
    let whose = crate::ui::session::Shelf::get(table.as_ref(), "session.owner");
    let mut next = page_on(&table);
    assert!(next.keeper.offer.is_none());
    assert!(!holds_a_session(&table));
    assert_eq!(next.keeper.take_problem(), None);
    assert_eq!(
        crate::ui::session::Shelf::get(table.as_ref(), "session.owner"),
        whose,
        "the page did not take the session over"
    );
}

#[test]
fn a_storage_that_takes_no_note_of_the_choice_still_loses_what_was_kept_and_keeps_nothing_in_this_run(
) {
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let mut app = page_on(&table);
    app.add_doc(dirty("a sequence nobody else should read", "ACGTACGT"));
    save_now(&mut app);
    *table.refused.borrow_mut() = Some("session.off".to_string());
    let said = app.keeper.set_keeping(false).expect("something to say");
    assert!(
        said.contains("could not be noted") && said.contains("session.off may not be written"),
        "{said}"
    );
    assert!(
        !holds_a_session(&table),
        "what was kept is removed all the same"
    );
    assert_eq!(app.keeper.keeps_work(), Some(false), "off for this run");
    app.add_doc(dirty("another", "GGGG"));
    app.flush_session();
    assert!(!holds_a_session(&table));
    // the site was never told, so another tab of it goes on remembering
    *table.refused.borrow_mut() = None;
    let mut other = page_on(&table);
    other.add_doc(dirty("the other tab's", "TTTT"));
    save_now(&mut other);
    // the person turns it on again in the first tab: what the other tab kept is not a left-over to be removed
    assert_eq!(app.keeper.set_keeping(true).as_deref(), Some(KEPT_AGAIN));
    let kept = Store::on(table.clone()).load_session().expect("a session");
    assert_eq!(kept.docs[0].name, "the other tab's");
    // and off, now that the storage takes the note
    assert_eq!(app.keeper.set_keeping(false).as_deref(), Some(NOT_KEPT_NOW));
    assert!(!remembers(&table) && !holds_a_session(&table));
}

/// The settings folder as a store, to look at what the program left in it.
fn store_of(folder: &Folder) -> Store {
    Store::open(Some(&folder.0)).expect("a store")
}

/// Turn the remembering on in a copy that had it off, as a person would. The lock of the folder is now and then refused
/// for a moment although nobody has it (`Folder::app` says why): the copy is then started again, which waits for it.
fn turn_on(folder: &Folder, mut app: App) -> App {
    let said = app.keeper.set_keeping(true);
    if app.keeper.other_copy {
        drop(app);
        return folder.app();
    }
    assert_eq!(said.as_deref(), Some(KEPT_AGAIN));
    app
}

#[test]
fn in_a_settings_folder_turning_it_off_removes_the_sessions_and_the_recent_files_and_lasts() {
    let folder = Folder::new();
    let opened = folder.file("opened.gb", "LOCUS");
    let mut app = folder.app();
    app.keeper.remember(&opened);
    assert_eq!(app.keeper.recent, std::slice::from_ref(&opened));
    app.add_doc(dirty("on a shared computer", "ACGTACGT"));
    save_now(&mut app);
    assert!(folder.names("session").iter().any(|n| n == "manifest.txt"));
    // a session that was put aside, and one of another version that could not be read
    for aside in ["session.previous", "session.unreadable"] {
        std::fs::create_dir_all(folder.0.join(aside)).expect("a folder");
        folder.file(&format!("{aside}/0-0.gb"), "ACGT");
    }

    assert_eq!(app.keeper.set_keeping(false).as_deref(), Some(NOT_KEPT_NOW));
    for gone in ["session", "session.previous", "session.unreadable"] {
        assert!(
            folder
                .names(gone)
                .iter()
                .all(|n| !n.ends_with(".gb") && n != "manifest.txt"),
            "{gone}: {:?}",
            folder.names(gone)
        );
    }
    assert!(app.keeper.recent.is_empty());
    assert!(
        store_of(&folder).recent().is_empty(),
        "the list of recent files is emptied"
    );
    // while it is off nothing is noted: not the work, not the files that are opened
    app.keeper.remember(&opened);
    assert!(app.keeper.recent.is_empty() && store_of(&folder).recent().is_empty());
    app.add_doc(dirty("more", "GG"));
    app.flush_session();
    assert!(folder.names("session").iter().all(|n| n != "manifest.txt"));
    drop(app);

    // the next run: as it was left; and a session that got there all the same is removed, not offered
    std::fs::create_dir_all(folder.0.join("session.previous")).expect("a folder");
    folder.file("session.previous/0-0.gb", "ACGT");
    let mut next = folder.app();
    assert!(!next.keeper.other_copy);
    assert!(next.keeper.offer.is_none());
    assert_eq!(next.keeper.keeps_work(), Some(false));
    assert!(folder.names("session.previous").is_empty());
    next.keeper.remember(&opened);
    assert!(next.keeper.recent.is_empty());
    // a copy that remembers nothing holds no lock: a second copy is not told that another has the session
    let second = folder.app_while_another_may_be_alive();
    assert!(!second.keeper.other_copy);
    drop(second);

    // and on again: this copy keeps the session and the recent files as before
    let mut next = turn_on(&folder, next);
    assert_eq!(next.keeper.keeps_work(), Some(true));
    next.keeper.remember(&opened);
    assert_eq!(store_of(&folder).recent(), [opened]);
    next.add_doc(dirty("kept again", "TTTT"));
    save_now(&mut next);
    drop(next);
    assert!(offer_of(&folder).is_some());
}

#[test]
fn turning_it_on_while_another_copy_has_taken_the_session_says_so_and_leaves_the_session_to_it() {
    let folder = Folder::new();
    store_of(&folder).want_sessions(false).expect("the mark");
    // two copies started while nothing was remembered: neither holds the session
    let first = folder.app();
    let mut second = folder.app_while_another_may_be_alive();
    assert_eq!(first.keeper.keeps_work(), Some(false));
    assert_eq!(second.keeper.keeps_work(), Some(false));
    // the first is told to remember again, and takes the session
    let mut first = turn_on(&folder, first);
    first.add_doc(dirty("the first copy's", "ACGT"));
    save_now(&mut first);
    // the second is told the same: the session is the first's, and it says so
    assert_eq!(second.keeper.set_keeping(true).as_deref(), Some(OTHER_COPY));
    assert_eq!(
        second.keeper.keeps_work(),
        None,
        "nothing to choose in this copy"
    );
    second.add_doc(dirty("the second copy's", "GG"));
    second.flush_session();
    drop(second);
    drop(first);
    let offer = offer_of(&folder).expect("the first copy's session");
    assert_eq!(offer.docs.len(), 1);
    assert_eq!(offer.docs[0].name, "the first copy's");
}

#[test]
fn the_switch_has_nothing_to_say_where_nothing_can_be_kept_or_another_copy_keeps_the_session() {
    let mut none = App::new(Startup::default());
    assert_eq!(none.keeper.keeps_work(), None);
    assert_eq!(none.keeper.set_keeping(false), None);
    assert_eq!(none.keeper.keeps_work(), None);

    let folder = Folder::new();
    let mut first = folder.app();
    first.add_doc(dirty("the first copy's", "ACGT"));
    save_now(&mut first);
    let mut second = folder.app_while_another_may_be_alive();
    assert!(second.keeper.other_copy);
    assert_eq!(second.keeper.keeps_work(), None);
    // the session is the first copy's: the second may not remove it, or turn the first's remembering off
    assert_eq!(second.keeper.set_keeping(false), None);
    assert!(folder.names("session").iter().any(|n| n == "manifest.txt"));
    assert!(store_of(&folder).sessions_wanted());
    assert_eq!(first.keeper.keeps_work(), Some(true));
}

/// The File menu of `app`, open on the screen.
fn file_menu(app: &mut App) -> crate::ui::testing::Screen {
    let mut screen = crate::ui::testing::Screen::new();
    screen.settle(app);
    screen.click(app, "File");
    screen.settle(app);
    screen
}

#[test]
fn the_file_menu_has_the_switch_and_a_click_turns_the_remembering_off_and_tells_the_person() {
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let mut app = page_on(&table);
    app.add_doc(dirty("kept so far", "ACGT"));
    save_now(&mut app);
    let mut screen = file_menu(&mut app);
    screen.click(&mut app, KEEP_SWITCH);
    screen.settle(&mut app);
    assert_eq!(app.keeper.keeps_work(), Some(false));
    assert!(!holds_a_session(&table));
    assert_eq!(app.status, NOT_KEPT_NOW);
    // the command palette has the same command, and it turns it on again
    let command = crate::ui::commands::all()
        .iter()
        .find(|c| matches!(c.action, Action::ToggleKeeping))
        .expect("a command for the palette");
    assert!(
        command.title.starts_with("Remember my work"),
        "{}",
        command.title
    );
    app.run(command.action.clone(), &screen.ctx().clone());
    assert_eq!(app.keeper.keeps_work(), Some(true));
    assert_eq!(app.status, KEPT_AGAIN);
}

#[test]
fn the_switch_is_not_in_the_menu_where_there_is_nothing_to_choose_and_the_command_says_why() {
    // nowhere to keep anything
    let mut none = App::new(Startup::default());
    none.docs.clear();
    let screen = file_menu(&mut none);
    assert!(screen.find("Save As…").is_some(), "the menu is open");
    assert!(screen.find(KEEP_SWITCH).is_none());
    none.run(Action::ToggleKeeping, &screen.ctx().clone());
    assert_eq!(none.status, NOWHERE_TO_KEEP);
    // another copy keeps the session
    let folder = Folder::new();
    let _first = folder.app();
    let mut second = folder.app_while_another_may_be_alive();
    let screen = file_menu(&mut second);
    assert!(screen.find("Save As…").is_some() && screen.find(KEEP_SWITCH).is_none());
    second.run(Action::ToggleKeeping, &screen.ctx().clone());
    assert_eq!(second.status, OTHER_COPY);
}

#[test]
fn while_the_last_session_is_on_offer_the_command_asks_for_an_answer_first_and_removes_nothing() {
    let table = std::rc::Rc::new(crate::ui::session::Memory::default());
    let mut before = page_on(&table);
    before.add_doc(dirty("work that was not saved", "ACGT"));
    save_now(&mut before);
    drop(before);
    let mut app = page_on(&table);
    assert!(app.keeper.offer.is_some());
    let ctx = egui::Context::default();
    app.run(Action::ToggleKeeping, &ctx);
    assert_eq!(app.status, ANSWER_FIRST);
    assert!(app.keeper.offer.is_some());
    assert!(holds_a_session(&table));
    assert_eq!(app.keeper.keeps_work(), Some(true));
}

#[cfg(unix)]
#[test]
fn a_list_of_recent_files_that_cannot_be_removed_is_not_shown_while_nothing_is_to_be_remembered() {
    use std::os::unix::fs::PermissionsExt;
    let folder = Folder::new();
    let opened = folder.file("opened.gb", "LOCUS");
    let mut app = folder.app();
    app.keeper.remember(&opened);
    drop(app);
    store_of(&folder).want_sessions(false).expect("the mark");
    // a settings folder that takes no more writes: the list is there, and cannot be emptied
    let read_only = |on: bool| {
        let mode = if on { 0o555 } else { 0o755 };
        std::fs::set_permissions(&folder.0, std::fs::Permissions::from_mode(mode)).expect("set");
    };
    read_only(true);
    if std::fs::write(folder.0.join("probe"), "x").is_ok() {
        // whoever runs the tests may write anywhere (root): there is no folder that refuses, and nothing to try
        read_only(false);
        return;
    }
    let mut next = folder.app();
    read_only(false);
    assert_eq!(
        store_of(&folder).recent(),
        [opened],
        "the list is still in the folder"
    );
    assert!(next.keeper.recent.is_empty(), "and is not shown");
    let said = next.keeper.take_problem().expect("the person is told");
    assert!(
        said.starts_with("What was not to be remembered could not be removed"),
        "{said}"
    );
}

mod results;

mod large;
