use super::*;
use crate::ui::app::{App, Startup};
use crate::ui::platform::Workspace;
use crate::ui::testing::Screen;
use eframe::egui::{pos2, Event, Modifiers, PointerButton, Pos2};

fn shown(reverse: bool, kept: Option<(usize, usize)>) -> Shown {
    Shown {
        bases: 480,
        reverse,
        kept,
    }
}

fn peak() -> Peak {
    Peak {
        at: 156,
        call: b'G',
        quality: Some(42),
        scan: Some(1902),
        signal: vec![(b'G', 1480), (b'A', 120), (b'C', 30), (b'T', -25)],
    }
}

#[test]
fn one_base_reads_as_its_number_its_call_how_sure_it_is_and_its_signals() {
    let line = words(&Picked::One(peak()), shown(false, Some((60, 400))));
    assert_eq!(
        line,
        "Base 157 of 480 · called G · quality 42 (1 wrong in 15,849) · scan 1,902 · signal G 1,480, A 120, C 30, T -25 · kept"
    );
    // outside what trimming keeps (its end is not kept), and with nothing trimmed
    let off = words(&Picked::One(peak()), shown(false, Some((60, 156))));
    assert!(off.ends_with(" · trimmed off"), "{off}");
    let first = Peak { at: 60, ..peak() };
    let kept = words(&Picked::One(first), shown(false, Some((60, 156))));
    assert!(
        kept.starts_with("Base 61 of 480") && kept.ends_with(" · kept"),
        "{kept}"
    );
    let plain = words(&Picked::One(peak()), shown(false, None));
    assert!(plain.ends_with("T -25"), "{plain}");
}

#[test]
fn a_mirrored_read_is_numbered_from_the_other_end_and_its_letters_are_the_complements() {
    let line = words(&Picked::One(peak()), shown(true, Some((60, 400))));
    assert_eq!(
        line,
        "Base 324 of 480 · called C · quality 42 (1 wrong in 15,849) · scan 1,902 · signal C 1,480, T 120, G 30, A -25 · kept"
    );
    // the first base of the file is the last of the picture
    let first = Peak { at: 0, ..peak() };
    assert!(words(&Picked::One(first), shown(true, None)).starts_with("Base 480 of 480 · "));
}

#[test]
fn a_base_without_quality_or_peak_says_so_and_a_stronger_channel_than_the_call_is_named() {
    let bare = Peak {
        quality: None,
        scan: None,
        signal: Vec::new(),
        ..peak()
    };
    assert_eq!(
        words(&Picked::One(bare), shown(false, None)),
        "Base 157 of 480 · called G · no quality in the file"
    );
    // called C, though A is stronger at the peak
    let mixed = Peak {
        call: b'c',
        signal: vec![(b'A', 700), (b'C', 600), (b'G', 0), (b'T', 0)],
        ..peak()
    };
    let line = words(&Picked::One(mixed.clone()), shown(false, None));
    assert!(
        line.ends_with(" · signal A 700, C 600, G 0, T 0 · the strongest signal is A"),
        "{line}"
    );
    // mirrored, it is named as the picture names it
    let mirrored = words(&Picked::One(mixed.clone()), shown(true, None));
    assert!(
        mirrored.ends_with(" · the strongest signal is T"),
        "{mirrored}"
    );
    // as strong as the call is not stronger, and a call that is no base has no channel of its own
    let level = Peak {
        signal: vec![(b'A', 600), (b'C', 600), (b'G', 0), (b'T', 0)],
        ..mixed.clone()
    };
    assert!(!words(&Picked::One(level), shown(false, None)).contains("strongest"));
    let unsure = Peak {
        call: b'N',
        ..mixed
    };
    assert!(!words(&Picked::One(unsure), shown(false, None)).contains("strongest"));
}

fn stretch() -> Stretch {
    Stretch {
        from: 149,
        to: 180,
        mean: Some(38.3),
        lowest: Some((12, 162)),
        poor: Some(2),
        unsure: 1,
    }
}

#[test]
fn several_bases_read_as_their_range_their_qualities_and_what_trimming_keeps_of_them() {
    let line = words(&Picked::Several(stretch()), shown(false, Some((60, 170))));
    assert_eq!(
        line,
        "Bases 150–180 of 480 (31) · mean quality 38.3 · lowest 12 at base 163 · 2 below quality 20 · 1 not A, C, G or T · 21 kept by the trimming"
    );
    // mirrored: the same bases, numbered from the other end
    let mirrored = words(&Picked::Several(stretch()), shown(true, None));
    assert_eq!(
        mirrored,
        "Bases 301–331 of 480 (31) · mean quality 38.3 · lowest 12 at base 318 · 2 below quality 20 · 1 not A, C, G or T"
    );
    // what trimming keeps is counted in the bases of the file, however the picture numbers them
    let both = words(&Picked::Several(stretch()), shown(true, Some((60, 170))));
    assert!(
        both.starts_with("Bases 301–331 of 480 (31) · ")
            && both.ends_with(" · 21 kept by the trimming"),
        "{both}"
    );
    // all kept, none kept, and one base of them kept at either end
    let kept = |range| words(&Picked::Several(stretch()), shown(false, Some(range)));
    assert!(kept((149, 180)).ends_with(" · all kept"));
    assert!(kept((0, 480)).ends_with(" · all kept"));
    assert!(kept((180, 300)).ends_with(" · all trimmed off"));
    assert!(kept((0, 149)).ends_with(" · all trimmed off"));
    assert!(kept((179, 300)).ends_with(" · 1 kept by the trimming"));
    assert!(kept((0, 150)).ends_with(" · 1 kept by the trimming"));
}

#[test]
fn good_bases_say_that_none_is_poor_and_a_file_without_qualities_says_that() {
    let good = Stretch {
        poor: Some(0),
        unsure: 0,
        ..stretch()
    };
    assert_eq!(
        words(&Picked::Several(good), shown(false, None)),
        "Bases 150–180 of 480 (31) · mean quality 38.3 · lowest 12 at base 163 · none below quality 20"
    );
    let bare = Stretch {
        mean: None,
        lowest: None,
        poor: None,
        unsure: 0,
        ..stretch()
    };
    assert_eq!(
        words(&Picked::Several(bare), shown(false, None)),
        "Bases 150–180 of 480 (31) · no quality in the file"
    );
}

/// The Signals page with the example trace.
fn window() -> (App, Screen) {
    let mut app = App::new(Startup::default());
    app.load_example(Workspace::Signals);
    app.workspace = Workspace::Signals;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    (app, screen)
}

/// Where the ruler of the trace writes `number`: over the peak of that base. The trace is right of the options, and
/// the ruler is the highest place such a number can be.
fn ruler(screen: &Screen, number: usize) -> Pos2 {
    screen
        .find_all(&number.to_string())
        .into_iter()
        .filter(|at| at.x > 420.0)
        .min_by(|a, b| a.y.total_cmp(&b.y))
        .unwrap_or_else(|| panic!("no {number} on the ruler:\n{}", screen.text()))
}

fn press(at: Pos2, pressed: bool) -> Event {
    Event::PointerButton {
        pos: at,
        button: PointerButton::Primary,
        pressed,
        modifiers: Modifiers::NONE,
    }
}

#[test]
fn a_click_on_a_base_says_what_it_is_and_deselect_lets_it_go() {
    let (mut app, mut screen) = window();
    assert!(screen.text().contains(HINT), "{}", screen.text());
    assert_eq!(app.platform.signals.picked_words(), None);

    // in the traces, under the ruler's 20
    let at = ruler(&screen, 20) + eframe::egui::vec2(0.0, 120.0);
    screen.click_at(&mut app, at);
    screen.settle(&mut app);
    let words = app
        .platform
        .signals
        .picked_words()
        .expect("a base is picked");
    let trace = &app.platform.signals.loaded.as_ref().expect("a trace").trace;
    let call = char::from(trace.bases[19]);
    assert!(
        words.starts_with(&format!(
            "Base 20 of 480 · called {call} · quality {} (1 wrong in ",
            trace.quality[19]
        )),
        "{words}"
    );
    assert!(
        words.contains(&format!(" · scan {} · signal ", trace.peaks[19])),
        "{words}"
    );
    // trimming keeps bases 10 to 464 of the example: this one is kept
    assert!(words.ends_with(" · kept"), "{words}");
    // the page writes those words, and no longer how to pick
    let text = screen.text();
    assert!(text.contains(&words), "{text}");
    assert!(!text.contains(HINT));

    screen.click(&mut app, CLEAR);
    screen.settle(&mut app);
    assert_eq!(app.platform.signals.picked_words(), None);
    assert!(screen.text().contains(HINT));
    assert!(screen.find(CLEAR).is_none());
}

#[test]
fn a_drag_across_bases_says_how_good_they_are_together_and_mirroring_renumbers_them() {
    let (mut app, mut screen) = window();
    let from = ruler(&screen, 20) + eframe::egui::vec2(0.0, 120.0);
    let to = pos2(ruler(&screen, 40).x, from.y);
    screen.frame(&mut app, vec![Event::PointerMoved(from)]);
    screen.frame(&mut app, vec![press(from, true)]);
    screen.frame(&mut app, vec![Event::PointerMoved(to)]);
    // the stretch is picked while the button is still down, so that it is seen growing under the pointer
    assert_eq!(app.platform.signals.selection, Some((19, 40)));
    screen.frame(&mut app, vec![press(to, false)]);
    screen.settle(&mut app);
    let words = app
        .platform
        .signals
        .picked_words()
        .expect("bases are picked");
    assert!(
        words.starts_with("Bases 20–40 of 480 (21) · mean quality "),
        "{words}"
    );
    assert!(screen.text().contains(&words), "{}", screen.text());
    // the same bases in a read shown reverse-complemented: numbered from the other end
    app.platform.signals.reverse = true;
    screen.settle(&mut app);
    let mirrored = app.platform.signals.picked_words().expect("still picked");
    assert!(
        mirrored.starts_with("Bases 441–461 of 480 (21) · mean quality "),
        "{mirrored}"
    );
    assert!(screen.text().contains(&mirrored), "{}", screen.text());
}

#[test]
fn another_trace_has_nothing_picked() {
    let (mut app, mut screen) = window();
    let at = ruler(&screen, 30) + eframe::egui::vec2(0.0, 120.0);
    screen.click_at(&mut app, at);
    screen.settle(&mut app);
    assert!(app.platform.signals.picked_words().is_some());
    app.platform
        .signals
        .load("again.ab1", crate::sample::demo_trace());
    assert_eq!(app.platform.signals.picked_words(), None);
    // and a page with no trace has no words, whatever was picked
    let mut bare = crate::ui::platform::signals::SignalsPage::default();
    bare.pick(Some((3, 9)));
    assert_eq!(bare.picked_words(), None);
}

#[test]
fn the_pictures_of_a_picked_base_and_of_picked_bases_show_their_words() {
    for (view, begins) in [
        ("peak", "Base 157 of 480 · called "),
        ("peaks", "Bases 150–180 of 480 (31) · mean quality "),
    ] {
        let mut app = App::new(Startup {
            view: view.to_string(),
            ..Startup::default()
        });
        assert_eq!(app.workspace, Workspace::Signals);
        let mut screen = Screen::new();
        screen.settle(&mut app);
        let words = app.platform.signals.picked_words().expect("picked");
        assert!(words.starts_with(begins), "{view}: {words}");
        assert!(screen.text().contains(&words), "{}", screen.text());
    }
}

#[test]
fn the_page_says_whether_trimming_keeps_the_base_and_says_nothing_of_it_when_nothing_is_trimmed() {
    let (mut app, _) = window();
    let page = &mut app.platform.signals;
    assert_eq!(page.kept(), Some((9, 464)), "what the example keeps");
    page.pick(Some((0, 1)));
    let first = page.picked_words().expect("picked");
    assert!(
        first.starts_with("Base 1 of 480 · ") && first.ends_with(" · trimmed off"),
        "{first}"
    );
    page.pick(Some((9, 10)));
    assert!(page.picked_words().expect("picked").ends_with(" · kept"));
    page.pick(Some((0, 12)));
    let some = page.picked_words().expect("picked");
    assert!(some.ends_with(" · 3 kept by the trimming"), "{some}");
    // with trimming off no base is kept or trimmed off
    page.trim = false;
    let plain = page.picked_words().expect("picked");
    assert!(
        !plain.contains("kept") && !plain.contains("trimmed"),
        "{plain}"
    );
}

#[test]
fn a_page_without_a_trace_does_not_say_how_to_pick_bases() {
    let mut app = App::new(Startup::default());
    app.workspace = Workspace::Signals;
    let mut screen = Screen::new();
    screen.settle(&mut app);
    let text = screen.text();
    assert!(text.contains("No trace yet"), "{text}");
    assert!(!text.contains(HINT), "{text}");
}

/// Drag from `from` to `to`, the last move and the letting go in one frame, as a quick hand does.
fn drag(screen: &mut Screen, app: &mut App, from: Pos2, to: Pos2) {
    let half = pos2((from.x + to.x) / 2.0, from.y);
    screen.frame(app, vec![Event::PointerMoved(from)]);
    screen.frame(app, vec![press(from, true)]);
    screen.frame(app, vec![Event::PointerMoved(half)]);
    screen.frame(app, vec![Event::PointerMoved(to), press(to, false)]);
    screen.settle(app);
}

#[test]
fn a_drag_to_the_left_picks_the_same_bases_and_the_last_move_counts_though_the_button_is_let_go_in_it(
) {
    let (mut app, mut screen) = window();
    let y = ruler(&screen, 20).y + 120.0;
    let (twenty, forty) = (ruler(&screen, 20).x, ruler(&screen, 40).x);
    drag(&mut screen, &mut app, pos2(forty, y), pos2(twenty, y));
    let words = app.platform.signals.picked_words().expect("picked");
    assert!(words.starts_with("Bases 20–40 of 480 (21) · "), "{words}");
    // and to the right, with the last move in the frame that ends it
    drag(&mut screen, &mut app, pos2(twenty, y), pos2(forty, y));
    let words = app.platform.signals.picked_words().expect("picked");
    assert!(words.starts_with("Bases 20–40 of 480 (21) · "), "{words}");
    assert_eq!(app.platform.signals.selection, Some((19, 40)));
}

#[test]
fn in_a_mirrored_picture_a_drag_picks_the_bases_under_it_numbered_as_the_ruler_numbers_them() {
    let (mut app, mut screen) = window();
    app.platform.signals.reverse = true;
    screen.settle(&mut app);
    let y = ruler(&screen, 20).y + 120.0;
    let (twenty, forty) = (ruler(&screen, 20).x, ruler(&screen, 40).x);
    drag(&mut screen, &mut app, pos2(twenty, y), pos2(forty, y));
    let words = app.platform.signals.picked_words().expect("picked");
    assert!(words.starts_with("Bases 20–40 of 480 (21) · "), "{words}");
    // which are the last bases of the file but nineteen
    assert_eq!(app.platform.signals.selection, Some((440, 461)));
}

#[test]
fn a_click_beside_the_bases_picks_the_base_nearest_to_it() {
    let (mut app, mut screen) = window();
    // left of the first base: the tenth is some 130 points into the trace
    let at = ruler(&screen, 10) + eframe::egui::vec2(-140.0, 120.0);
    screen.click_at(&mut app, at);
    screen.settle(&mut app);
    let words = app.platform.signals.picked_words().expect("the first base");
    assert!(words.starts_with("Base 1 of 480 · "), "{words}");
}

#[test]
fn the_line_is_one_row_whether_it_says_little_or_much() {
    let (mut app, mut screen) = window();
    app.platform.signals.pick(Some((149, 180)));
    screen.settle(&mut app);
    let short = screen.find(CLEAR).expect("the link");
    let ruler_then = ruler(&screen, 20);
    // one base says more than fits beside the link at this width
    app.platform.signals.pick(Some((156, 157)));
    screen.settle(&mut app);
    let long = screen.find(CLEAR).expect("the link");
    assert_eq!(long, short, "the link keeps its place");
    assert!(long.x < 1280.0 - 18.0, "{long:?}");
    assert_eq!(ruler(&screen, 20), ruler_then, "and the trace under it too");
}
