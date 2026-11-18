//! Where the tool windows were left, to open them there again (`places.txt` in the settings folder): a window title and the corner of
//! the window, in points, a line each. What is read is not trusted: a line that is not understood is left out and spoils no other,
//! a corner that is not a number or is far away is dropped, and the list is bounded.

use std::collections::BTreeMap;

/// The most places kept, the longest title, and the furthest a corner may be from the origin (points).
pub const MOST_PLACES: usize = 200;
pub const MOST_TITLE: usize = 100;
const FAR: f32 = 20_000.0;
/// The nearest to the top left a window opens: clear of the rail and of the top bar, so that its title can be grabbed.
pub const FLOOR: (f32, f32) = (96.0, 50.0);

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Places(BTreeMap<String, (f32, f32)>);

/// Whether a title and a corner can be kept.
fn keepable(title: &str, x: f32, y: f32) -> bool {
    !title.is_empty()
        && title.chars().count() <= MOST_TITLE
        && !title.chars().any(char::is_control)
        && x.is_finite()
        && y.is_finite()
        && x.abs() <= FAR
        && y.abs() <= FAR
}

impl Places {
    #[cfg(test)]
    pub fn get(&self, title: &str) -> Option<(f32, f32)> {
        self.0.get(title).copied()
    }

    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&String, &(f32, f32))> {
        self.0.iter()
    }

    /// Keep the corner of the window `title`; whether that changed anything by more than a point. A place that cannot be kept, or one
    /// more than there is room for, is not.
    pub fn set(&mut self, title: &str, x: f32, y: f32) -> bool {
        if !keepable(title, x, y) {
            return false;
        }
        match self.0.get(title) {
            Some(&(a, b)) if (a - x).abs() <= 1.0 && (b - y).abs() <= 1.0 => false,
            Some(_) => {
                self.0.insert(title.to_string(), (x, y));
                true
            }
            None if self.0.len() >= MOST_PLACES => false,
            None => {
                self.0.insert(title.to_string(), (x, y));
                true
            }
        }
    }

    /// Forget the place of `title`; whether there was one.
    pub fn remove(&mut self, title: &str) -> bool {
        self.0.remove(title).is_some()
    }

    /// A line each: the two numbers, then the title (which may have spaces).
    pub fn to_text(&self) -> String {
        self.0
            .iter()
            .map(|(title, (x, y))| format!("{x} {y} {title}\n"))
            .collect()
    }

    pub fn parse(text: &str) -> Places {
        let mut places = Places::default();
        for line in text.lines() {
            let mut parts = line.splitn(3, ' ');
            let (Some(x), Some(y), Some(title)) = (parts.next(), parts.next(), parts.next()) else {
                continue;
            };
            if let (Ok(x), Ok(y)) = (x.parse::<f32>(), y.parse::<f32>()) {
                places.set(title, x, y);
            }
        }
        places
    }
}

/// A place brought to where a window can be grabbed: not above or left of the floor.
pub fn floored(x: f32, y: f32) -> (f32, f32) {
    (x.max(FLOOR.0), y.max(FLOOR.1))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn places_go_to_text_and_back_and_a_line_that_is_not_understood_spoils_no_other() {
        let mut places = Places::default();
        assert!(places.set("Find", 300.5, 200.0));
        assert!(places.set("Species Tree from Gene Trees", 120.0, 80.25));
        let text = places.to_text();
        assert_eq!(
            text,
            "300.5 200 Find\n120 80.25 Species Tree from Gene Trees\n"
        );
        assert_eq!(Places::parse(&text), places);
        let messy = "nonsense\n\n1 2\nx y Find\n300.5 200 Find\n  \n5 6 \nNaN 3 Bad\n7 inf Bad\n120 80.25 Species Tree from Gene Trees\n";
        assert_eq!(Places::parse(messy), places);
    }

    #[test]
    fn a_place_changes_only_when_the_window_moved_and_a_place_that_cannot_be_kept_is_refused() {
        let mut places = Places::default();
        assert!(places.set("Find", 100.0, 100.0));
        assert!(
            !places.set("Find", 100.5, 99.5),
            "within a point is the same place"
        );
        assert!(places.set("Find", 103.0, 100.0));
        assert_eq!(places.get("Find"), Some((103.0, 100.0)));
        for (title, x, y) in [
            ("", 1.0, 1.0),
            ("Find", f32::NAN, 1.0),
            ("Find", 1.0, f32::INFINITY),
            ("Find", 30_000.0, 1.0),
            ("Find", 1.0, -30_000.0),
            ("Tab\there", 1.0, 1.0),
        ] {
            assert!(!places.set(title, x, y), "{title:?} {x} {y}");
        }
        assert!(!places.set(&"x".repeat(MOST_TITLE + 1), 1.0, 1.0));
        assert!(places.set(&"x".repeat(MOST_TITLE), 1.0, 1.0));
        assert!(places.remove("Find") && !places.remove("Find"));
    }

    #[test]
    fn the_list_is_bounded_whatever_a_file_holds() {
        let text: String = (0..MOST_PLACES + 50)
            .map(|k| format!("{k} {k} Window {k}\n"))
            .collect();
        assert_eq!(Places::parse(&text).len(), MOST_PLACES);
        // a full list still takes a change to a window it has
        let mut places = Places::parse(&text);
        assert!(places.set("Window 3", 500.0, 500.0));
        assert!(!places.set("One more", 5.0, 5.0));
    }

    #[test]
    fn a_window_is_never_put_above_or_left_of_the_floor() {
        assert_eq!(floored(10.0, 5.0), FLOOR);
        assert_eq!(floored(300.0, 200.0), (300.0, 200.0));
        assert_eq!(floored(-500.0, 900.0), (FLOOR.0, 900.0));
        // the floor is clear of the rail (92 wide) and the top bar (46 high)
        assert!(FLOOR.0 > 92.0 && FLOOR.1 > 46.0);
    }
}
