use super::*;
use crate::checks::{self, CHECKS};
use crate::ui::app::{Action, App, Startup};
use crate::ui::commands;
use crate::ui::testing::Screen;

/// The app with the page open, as the command of the Help menu opens it.
fn opened() -> (App, Screen) {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    app.run(Action::Checks, screen.ctx());
    screen.settle(&mut app);
    (app, screen)
}

/// The page with `words` in its filter.
fn filtered(words: &str) -> (App, Screen) {
    let (mut app, mut screen) = opened();
    app.checks_win.query = words.into();
    screen.settle(&mut app);
    (app, screen)
}

#[test]
fn the_command_opens_the_page_with_a_count_for_each_grade_and_the_first_group() {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert!(!screen.text().contains(TITLE));
    app.run(Action::Checks, screen.ctx());
    screen.settle(&mut app);
    let text = screen.text();
    assert!(text.contains(TITLE), "{text}");
    let all: Vec<&checks::Check> = CHECKS.iter().collect();
    for (standard, n) in checks::tally(&all) {
        let chip = format!("{n} {}", standard.label());
        assert!(text.contains(&chip), "{chip} in {text}");
    }
    // the groups are headings in capitals, and the first analyses are listed under the first
    assert!(text.contains(&CHECKS[0].group.to_uppercase()), "{text}");
    assert!(text.contains(CHECKS[0].analysis), "{text}");
    // nothing is open: what an analysis was compared with is read after a click
    assert!(!text.contains(CHECKS[0].against), "{text}");
}

#[test]
fn the_command_is_in_the_palette_and_does_not_take_the_person_to_the_sequence() {
    let command = commands::all()
        .iter()
        .find(|c| matches!(c.action, Action::Checks))
        .expect("a command for the page");
    assert_eq!(command.title, "How each analysis was checked");
    assert_eq!(command.group, "Help");
    assert!(!commands::acts_on_the_sequence(&Action::Checks));
}

#[test]
fn a_filter_leaves_the_checks_that_hold_its_words_and_says_when_there_is_none() {
    let (mut app, mut screen) = filtered("primer3");
    let text = screen.text();
    assert!(text.contains("Primer design (Primer3)"), "{text}");
    assert!(!text.contains("Species tree from gene trees"), "{text}");
    assert!(!text.contains(NONE), "{text}");
    // the counts are of what is shown
    let shown = checks::matching("primer3");
    for (standard, n) in checks::tally(&shown) {
        let chip = format!("{n} {}", standard.label());
        assert!(text.contains(&chip), "{chip} in {text}");
    }
    app.checks_win.query = "nosuchword".into();
    screen.settle(&mut app);
    let text = screen.text();
    assert!(text.contains(NONE), "{text}");
    assert!(!text.contains("Primer design (Primer3)"), "{text}");
}

#[test]
fn a_click_on_an_analysis_shows_what_it_was_compared_with_and_its_tests_and_another_hides_them() {
    // a check whose texts are short, so that all of them are on the screen at once
    let words = "another program's bootstrap";
    let (mut app, mut screen) = filtered(words);
    let found = checks::matching(words);
    assert_eq!(found.len(), 1);
    let check = found[0];
    assert!(!check.tests.is_empty());
    assert!(!screen.text().contains(check.against));
    screen.click(&mut app, check.analysis);
    screen.settle(&mut app);
    let text = screen.text();
    for part in [check.against, check.result, check.limits] {
        assert!(text.contains(part), "{part} in {text}");
    }
    for (file, test) in check.tests {
        assert!(
            text.contains(test) && text.contains(file),
            "{test} in {text}"
        );
    }
    screen.click(&mut app, check.analysis);
    screen.settle(&mut app);
    assert!(!screen.text().contains(check.against));
}

#[test]
fn copy_puts_the_checks_that_are_shown_on_the_clipboard_as_the_command_line_prints_them() {
    let (mut app, mut screen) = filtered("primer3");
    assert_eq!(screen.copied(), None);
    screen.click(&mut app, COPY);
    screen.settle(&mut app);
    let shown = checks::matching("primer3");
    assert!(!shown.is_empty() && shown.len() < CHECKS.len());
    assert_eq!(screen.copied(), Some(checks::to_text(&shown).as_str()));
}

#[test]
fn the_view_named_checks_starts_with_the_page_open() {
    let mut app = App::new(Startup {
        view: "checks".into(),
        ..Startup::default()
    });
    let mut screen = Screen::new();
    screen.settle(&mut app);
    assert!(app.checks_win.open);
    assert!(screen.text().contains(TITLE));
}
