use super::*;
use crate::seq::io;
use std::sync::atomic::{AtomicU32, Ordering};

static NEXT: AtomicU32 = AtomicU32::new(0);

/// A folder of its own for one test, removed when the test ends.
struct Folder(PathBuf);

impl Folder {
    fn new() -> Folder {
        let path = std::env::temp_dir().join(format!(
            "helixloom-session-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }

    fn store(&self) -> Store {
        Store::open(Some(&self.0)).expect("a store")
    }

    fn file(&self, name: &str, text: &str) -> PathBuf {
        let path = self.0.join(name);
        fs::write(&path, text).expect("a file");
        path
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn doc(name: &str, seq: &str) -> Doc {
    let mut d = Doc::new(name, seq.as_bytes().to_vec());
    d.circular = true;
    d
}

#[test]
fn text_with_tabs_newlines_and_backslashes_reads_back_as_it_was_written_into_a_line_of_fields() {
    for text in [
        "",
        "plain",
        "a\tb",
        "two\nlines\r\nand a return",
        "back\\slash",
        "ends with a backslash\\",
        "looks like an escape: \\t \\n \\\\",
        "ünïcode — ✓",
    ] {
        let line = escape(text);
        assert!(!line.contains(['\t', '\n', '\r']), "{line:?}");
        assert_eq!(unescape(&line), text);
    }
}

#[test]
fn a_recent_file_goes_to_the_top_once_and_the_list_has_a_limit() {
    let list: Vec<PathBuf> = ["a", "b", "c"].iter().map(PathBuf::from).collect();
    assert_eq!(
        push_recent(&list, Path::new("b"), 5),
        ["b", "a", "c"].map(PathBuf::from)
    );
    assert_eq!(
        push_recent(&list, Path::new("d"), 3),
        ["d", "a", "b"].map(PathBuf::from)
    );
    assert_eq!(push_recent(&[], Path::new("x"), 3), [PathBuf::from("x")]);
}

#[test]
fn recent_files_are_kept_across_runs_latest_first_and_a_file_that_has_gone_is_left_out() {
    let folder = Folder::new();
    let (a, b, c) = (
        folder.file("a.gb", "x"),
        folder.file("b.gb", "x"),
        folder.file("c with\ttab.gb", "x"),
    );
    let store = folder.store();
    assert!(store.recent().is_empty());
    for p in [&a, &b, &c, &a] {
        store.add_recent(p).expect("added");
    }
    // a new run reads the same list
    let again = folder.store();
    assert_eq!(again.recent(), [a.clone(), c.clone(), b.clone()]);
    fs::remove_file(&c).expect("removed");
    assert_eq!(again.recent(), [a.clone(), b]);
    // more than the limit: the oldest fall off
    for k in 0..MOST_RECENT + 3 {
        let p = folder.file(&format!("n{k}.gb"), "x");
        store.add_recent(&p).expect("added");
    }
    assert_eq!(store.recent().len(), MOST_RECENT);
    store.clear_recent().expect("cleared");
    assert!(store.recent().is_empty());
}

#[test]
fn the_place_of_the_window_is_kept_and_nonsense_is_not_taken() {
    let folder = Folder::new();
    let store = folder.store();
    assert_eq!(store.window(), None);
    let place = WindowPlace {
        x: 120.5,
        y: 60.0,
        width: 1280.0,
        height: 800.0,
    };
    store.save_window(&place).expect("saved");
    assert_eq!(store.window(), Some(place));
    for bad in [
        "1 2 3",
        "1 2 3 4 5",
        "a b c d",
        "0 0 10 10",
        "0 0 100000 800",
        "NaN 0 800 600",
        "inf 0 800 600",
        "99999 0 800 600",
    ] {
        fs::write(folder.0.join(WINDOW), bad).expect("written");
        assert_eq!(store.window(), None, "{bad}");
    }
}

#[test]
fn a_session_comes_back_as_it_was_with_text_only_for_what_its_file_does_not_hold() {
    let folder = Folder::new();
    let store = folder.store();
    let mut untitled = doc("pasted\tplasmid", "ACGTACGTAC");
    untitled.dirty = true;
    let mut edited = doc("edited", "GGGGCCCCAA");
    edited.path = Some(folder.0.join("edited.gb"));
    edited.dirty = true;
    let mut clean = doc("clean", "TTTTAAAA");
    clean.path = Some(folder.0.join("clean.gb"));
    let saved: Vec<SavedDoc> = [&untitled, &edited, &clean]
        .iter()
        .map(|d| SavedDoc::of(d))
        .collect();
    assert!(saved[0].text.is_some() && saved[1].text.is_some());
    assert!(
        saved[2].text.is_none(),
        "a clean sequence is read back from its file"
    );
    let session = Session {
        workspace: "alignment".into(),
        active: 2,
        docs: saved,
        lost: 0,
        results: Vec::new(),
    };
    store.save_session(&session).expect("saved");
    let back = folder.store().load_session().expect("a session");
    assert_eq!(back, session);
    // the text is the sequence: it reads back as one
    let docs = io::parse_any(back.docs[0].text.as_deref().expect("text"), "x").expect("a record");
    assert_eq!(docs[0].seq, b"ACGTACGTAC");
    assert!(docs[0].circular);
}

#[test]
fn a_new_save_replaces_the_old_files_and_a_save_that_stopped_halfway_leaves_the_old_session_whole()
{
    let folder = Folder::new();
    let store = folder.store();
    let session = |text: &str| Session {
        workspace: "sequences".into(),
        active: 0,
        docs: vec![SavedDoc {
            name: "d".into(),
            path: None,
            dirty: true,
            text: Some(text.into()),
            keep: None,
        }],
        lost: 0,
        results: Vec::new(),
    };
    store.save_session(&session("first")).expect("saved");
    let names = |dir: &Path| -> Vec<String> {
        let mut v: Vec<String> = fs::read_dir(dir)
            .expect("a folder")
            .flatten()
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .collect();
        v.sort();
        v
    };
    let dir = folder.0.join(SESSION);
    assert_eq!(names(&dir), ["1-0.gb", "manifest.txt"]);
    // the program stops after it wrote the next generation's file and before it replaced the manifest
    fs::write(dir.join("2-0.gb"), "half of the second").expect("written");
    fs::write(dir.join(".manifest.txt.tmp"), "half a manifest").expect("written");
    let back = store.load_session().expect("a session");
    assert_eq!(back.docs[0].text.as_deref(), Some("first"));
    // the next save goes on from the manifest, and leaves only what it names
    store.save_session(&session("second")).expect("saved");
    assert_eq!(names(&dir), ["2-0.gb", "manifest.txt"]);
    assert_eq!(
        store.load_session().expect("a session").docs[0]
            .text
            .as_deref(),
        Some("second")
    );
}

#[test]
fn a_file_of_the_session_that_is_gone_is_counted_and_the_rest_come_back() {
    let folder = Folder::new();
    let store = folder.store();
    let docs: Vec<SavedDoc> = ["one", "two", "three"]
        .iter()
        .map(|n| SavedDoc {
            name: n.to_string(),
            path: None,
            dirty: true,
            text: Some(format!("text of {n}")),
            keep: None,
        })
        .collect();
    store
        .save_session(&Session {
            workspace: "home".into(),
            active: 0,
            docs,
            lost: 0,
            results: Vec::new(),
        })
        .expect("saved");
    fs::remove_file(folder.0.join(SESSION).join("1-1.gb")).expect("removed");
    let back = store.load_session().expect("a session");
    assert_eq!(back.lost, 1);
    assert_eq!(
        back.docs
            .iter()
            .map(|d| d.name.as_str())
            .collect::<Vec<_>>(),
        ["one", "three"]
    );
}

#[test]
fn nothing_open_clears_the_session_and_a_manifest_that_is_not_one_is_not_read() {
    let folder = Folder::new();
    let store = folder.store();
    assert!(store.load_session().is_none());
    store
        .save_session(&Session {
            workspace: "home".into(),
            active: 0,
            docs: vec![SavedDoc {
                name: "d".into(),
                path: None,
                dirty: true,
                text: Some("t".into()),
                keep: None,
            }],
            lost: 0,
            results: Vec::new(),
        })
        .expect("saved");
    assert!(store.load_session().is_some());
    store.save_session(&Session::default()).expect("cleared");
    assert!(!folder.0.join(SESSION).exists());
    assert!(store.load_session().is_none());
    // clearing what is not there is not an error
    store.clear_session().expect("nothing to clear");
    // a manifest of another program, or of another version, or garbage
    fs::create_dir_all(folder.0.join(SESSION)).expect("a folder");
    for text in [
        "",
        "hello",
        "helixloom-session\t2\ndoc\t1\ta\t\t\n",
        "\u{0}\u{1}\u{2}",
    ] {
        fs::write(folder.0.join(SESSION).join(MANIFEST), text).expect("written");
        assert!(store.load_session().is_none(), "{text:?}");
    }
}

#[test]
fn without_a_settings_folder_nothing_is_kept() {
    assert!(Store::open(None).is_none());
}

#[test]
fn a_manifest_cannot_make_the_program_read_a_file_outside_the_session_folder() {
    let folder = Folder::new();
    let store = folder.store();
    let secret = folder.file("secret.txt", "LOCUS secret");
    let dir = folder.0.join(SESSION);
    fs::create_dir_all(&dir).expect("a folder");
    fs::write(dir.join("1-0.gb"), "inside").expect("written");
    let manifest = format!(
        "{MAGIC}\ngeneration\t1\nworkspace\thome\nactive\t0\n\
         doc\t1\tinside\t\t1-0.gb\n\
         doc\t1\tup\t\t../secret.txt\n\
         doc\t1\tabsolute\t\t{}\n\
         doc\t1\thidden\t\t.manifest.txt.tmp\n",
        secret.display()
    );
    fs::write(dir.join(MANIFEST), manifest).expect("written");
    let back = store.load_session().expect("a session");
    assert_eq!(back.docs.len(), 1);
    assert_eq!(back.docs[0].text.as_deref(), Some("inside"));
    assert_eq!(back.lost, 3);
}

fn one(name: &str, text: &str) -> Session {
    Session {
        workspace: "sequences".into(),
        active: 0,
        docs: vec![SavedDoc {
            name: name.into(),
            path: None,
            dirty: true,
            text: Some(text.into()),
            keep: None,
        }],
        lost: 0,
        results: Vec::new(),
    }
}

#[cfg(unix)]
#[test]
fn the_session_folder_and_its_files_are_private_to_the_person() {
    use std::os::unix::fs::PermissionsExt;
    let folder = Folder::new();
    folder
        .store()
        .save_session(&one("d", "secret"))
        .expect("saved");
    let mode = |p: PathBuf| fs::metadata(p).expect("metadata").permissions().mode() & 0o777;
    let dir = folder.0.join(SESSION);
    assert_eq!(mode(dir.clone()), 0o700);
    assert_eq!(mode(dir.join(MANIFEST)), 0o600);
    assert_eq!(mode(dir.join("1-0.gb")), 0o600);
}

#[test]
fn only_files_the_program_wrote_are_cleaned_away_and_a_failed_write_leaves_no_temporary_file() {
    let folder = Folder::new();
    let store = folder.store();
    store.save_session(&one("d", "first")).expect("saved");
    let dir = folder.0.join(SESSION);
    // a person's own file in the folder, and a leftover of an interrupted save
    fs::write(dir.join("notes.txt"), "mine").expect("written");
    fs::write(dir.join("report.gb"), "mine too").expect("written");
    fs::write(dir.join(".1-0.gb.tmp"), "half").expect("written");
    store.save_session(&one("d", "second")).expect("saved");
    let mut names: Vec<String> = fs::read_dir(&dir)
        .expect("a folder")
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    assert_eq!(names, ["2-0.gb", "manifest.txt", "notes.txt", "report.gb"]);
    // a write that fails removes its temporary file: a folder where the file should be
    let target = folder.0.join("blocked");
    fs::create_dir_all(&target).expect("a folder");
    assert!(write_atomic(&target, b"x", false).is_err());
    assert!(!folder.0.join(".blocked.tmp").exists());
}

#[test]
fn a_session_can_be_put_aside_and_a_session_of_another_version_is_put_aside_not_written_over() {
    let folder = Folder::new();
    let store = folder.store();
    store.archive_session().expect("nothing to put aside");
    store.save_session(&one("old", "old text")).expect("saved");
    store.archive_session().expect("put aside");
    assert!(!folder.0.join(SESSION).exists());
    assert!(folder.0.join(PREVIOUS).join(MANIFEST).exists());
    store.save_session(&one("new", "new text")).expect("saved");
    assert_eq!(store.load_session().expect("a session").docs[0].name, "new");
    store.clear_previous().expect("removed");
    assert!(!folder.0.join(PREVIOUS).exists());
    // nothing to remove is not a failure
    store.clear_previous().expect("nothing there");
    // a manifest that a newer version wrote: it is moved aside, whole, before this version writes its own
    fs::write(
        folder.0.join(SESSION).join(MANIFEST),
        "helixloom-session\t2\nsomething new\n",
    )
    .expect("written");
    store.save_session(&one("mine", "my text")).expect("saved");
    assert_eq!(
        fs::read_to_string(folder.0.join(UNREADABLE).join(MANIFEST)).expect("kept"),
        "helixloom-session\t2\nsomething new\n"
    );
    assert_eq!(
        store.load_session().expect("a session").docs[0].name,
        "mine"
    );
}

#[test]
fn a_sequence_that_has_not_changed_keeps_its_file_and_a_kept_file_that_is_gone_is_an_error_not_a_loss(
) {
    let folder = Folder::new();
    let store = folder.store();
    let files = store.save_session(&one("d", "first")).expect("saved");
    assert_eq!(files, [Some("1-0.gb".to_string())]);
    let mut again = one("d", "unused");
    again.docs[0].text = None;
    again.docs[0].keep = files[0].clone();
    let files = store.save_session(&again).expect("saved");
    assert_eq!(
        files,
        [Some("1-0.gb".to_string())],
        "the same file, not a new one"
    );
    assert_eq!(
        store.load_session().expect("a session").docs[0]
            .text
            .as_deref(),
        Some("first")
    );
    // the file is gone: the save says so (and writes nothing), so that the program writes the text itself next time
    fs::remove_file(folder.0.join(SESSION).join("1-0.gb")).expect("removed");
    assert!(store.save_session(&again).is_err());
}

#[test]
fn only_one_copy_can_hold_the_session_lock_and_it_is_free_again_when_that_copy_ends() {
    let folder = Folder::new();
    let store = folder.store();
    let held = store.lock_session();
    assert!(matches!(held, Lock::Held(_)), "{held:?}");
    assert!(matches!(store.lock_session(), Lock::Busy));
    drop(held);
    // free again. A program that another test is starting at this very moment has the file too until it has started, and
    // the lock is then refused for some tens of milliseconds though no copy holds it (the tests of the keeper say how
    // that was measured): it is asked for a moment, not once
    let again = (0..200).any(|_| {
        let free = matches!(store.lock_session(), Lock::Held(_));
        if !free {
            std::thread::sleep(std::time::Duration::from_millis(10));
        }
        free
    });
    assert!(
        again,
        "the lock is free once the copy that held it has ended"
    );
}

#[test]
fn names_with_a_colon_or_a_separator_are_not_names_of_files_in_a_folder() {
    for bad in ["", ".hidden", "a/b", "a\\b", "C:x", "..", "a\0b"] {
        assert!(!plain_name(bad), "{bad:?}");
    }
    for good in ["1-0.gb", "12-345.gb", "manifest.txt"] {
        assert!(plain_name(good), "{good:?}");
    }
    assert!(is_ours("1-0.gb") && is_ours("12-345.gb") && is_ours(".manifest.txt.tmp"));
    for not in [
        "notes.txt",
        "report.gb",
        "a-b.gb",
        "1-.gb",
        "-1.gb",
        "manifest.txt",
        "1-0.gb.bak",
    ] {
        assert!(!is_ours(not), "{not}");
    }
}

#[test]
fn the_appearance_file_is_read_as_a_short_plain_file_and_nothing_else() {
    let folder = Folder::new();
    let store = folder.store();
    assert_eq!(store.appearance(), None, "none kept");
    folder.file("appearance.txt", "theme dark\nscale 125\n");
    assert_eq!(
        store.appearance().map(|(t, a)| (t, a.scale)),
        Some((ThemeMode::Dark, 125))
    );
    // no more than the first few thousand bytes are read, whatever the size: a line after them is not seen
    let big = format!(
        "scale 90\n{}density compact\n",
        "#\n".repeat(MAX_APPEARANCE_BYTES as usize)
    );
    folder.file("appearance.txt", &big);
    let read = store.appearance().expect("read");
    assert_eq!(read.1.scale, 90);
    assert_eq!(read.1.density, crate::ui::appearance::Density::Comfortable);
    // a byte that is not text, or a byte-order mark, does not stop the rest being read
    std::fs::write(
        folder.0.join("appearance.txt"),
        b"\xef\xbb\xbftheme light\nscale 110\n\xff\xfe\n",
    )
    .expect("a file");
    assert_eq!(
        store.appearance().map(|(t, a)| (t, a.scale)),
        Some((ThemeMode::Light, 110))
    );
    // a folder of that name is not a file: nothing is read
    std::fs::remove_file(folder.0.join("appearance.txt")).expect("removed");
    std::fs::create_dir(folder.0.join("appearance.txt")).expect("a folder");
    assert_eq!(store.appearance(), None);
}

#[cfg(unix)]
#[test]
fn a_device_in_the_place_of_the_appearance_file_is_not_read() {
    // a link to /dev/zero would go on for ever, and a pipe would wait for ever
    let folder = Folder::new();
    std::os::unix::fs::symlink("/dev/zero", folder.0.join("appearance.txt")).expect("a link");
    assert_eq!(folder.store().appearance(), None);
}

#[test]
fn the_places_of_the_tool_windows_are_kept_read_with_a_limit_and_a_damaged_file_gives_none() {
    let folder = Folder::new();
    let store = folder.store();
    assert!(store.places().is_empty(), "nothing kept yet");
    let mut places = crate::ui::places::Places::default();
    places.set("Find", 300.0, 200.0);
    places.set("Motif Finder", 150.5, 90.0);
    store.save_places(&places).expect("saved");
    assert_eq!(store.places(), places);
    assert_eq!(
        fs::read_to_string(folder.0.join(PLACES)).expect("a file"),
        places.to_text()
    );
    // lines that are not places do not matter
    fs::write(folder.0.join(PLACES), "junk\n300 200 Find\nx y z\n").expect("written");
    assert_eq!(store.places().len(), 1);
    // a file far larger than a list can be is read only so far, and a folder where the file should be is no file
    fs::write(folder.0.join(PLACES), "1 2 Window\n".repeat(100_000)).expect("written");
    assert!(store.places().len() <= crate::ui::places::MOST_PLACES);
    // what lies past the limit is not read at all: a place before it is found, one after it is not
    let padding = "junk\n".repeat(MAX_PLACES_BYTES as usize / 5 + 10);
    fs::write(
        folder.0.join(PLACES),
        format!("300 200 Early\n{padding}400 300 Late\n"),
    )
    .expect("written");
    let read = store.places();
    assert_eq!(read.get("Early"), Some((300.0, 200.0)));
    assert_eq!(read.get("Late"), None);
    fs::remove_file(folder.0.join(PLACES)).expect("removed");
    fs::create_dir(folder.0.join(PLACES)).expect("a folder");
    assert!(store.places().is_empty());
}
