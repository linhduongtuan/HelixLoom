//! The line of the Signals page that says what the picked base of the trace is, or the picked bases (`docs/ROADMAP.md` 1.4).
//! It speaks as the picture does: in a read that is shown reverse-complemented the bases are numbered from the other end
//! and the letters are the complements.

use super::SignalsPage;
use crate::seq::peak::{one_in, picked, Peak, Picked, Stretch, POOR};
use crate::ui::map::commify;
use crate::ui::platform::chromatogram::shown_base;

/// What the line says while nothing is picked.
pub const HINT: &str = "Click a base, or drag across some, to see what was called there.";
/// The link that lets the pick go.
pub const CLEAR: &str = "Deselect";

/// How the picture shows the trace: how many bases it has, whether it is mirrored, and the bases that trimming keeps
/// (counted as the file counts them; none when nothing is trimmed).
#[derive(Clone, Copy, Debug)]
pub struct Shown {
    pub bases: usize,
    pub reverse: bool,
    pub kept: Option<(usize, usize)>,
}

impl Shown {
    /// The number the ruler writes over base `at` (counted from 0 in the file).
    fn number(&self, at: usize) -> usize {
        if self.reverse {
            self.bases.saturating_sub(at)
        } else {
            at + 1
        }
    }

    fn letter(&self, base: u8) -> char {
        char::from(shown_base(base, self.reverse))
    }
}

/// A signal with a comma at every thousand; one below nothing keeps its sign.
fn signal_text(signal: i16) -> String {
    let size = commify(usize::from(signal.unsigned_abs()));
    if signal < 0 {
        format!("-{size}")
    } else {
        size
    }
}

fn peak_words(p: &Peak, s: Shown) -> String {
    let mut parts = vec![
        format!("Base {} of {}", commify(s.number(p.at)), commify(s.bases)),
        format!("called {}", s.letter(p.call)),
        match p.quality {
            Some(q) => format!("quality {q} (1 wrong in {})", commify(one_in(q) as usize)),
            None => "no quality in the file".to_string(),
        },
    ];
    if let Some(scan) = p.scan {
        parts.push(format!("scan {}", commify(scan as usize)));
        let signals: Vec<String> = p
            .signal
            .iter()
            .map(|&(base, signal)| format!("{} {}", s.letter(base), signal_text(signal)))
            .collect();
        parts.push(format!("signal {}", signals.join(", ")));
        // another channel is stronger than the call's own at the call's peak: a mixed peak, or a call to look at again
        let call = p.call.to_ascii_uppercase();
        let own = p.signal.iter().find(|&&(base, _)| base == call);
        if let (Some(&(top, most)), Some(&(_, own))) = (p.signal.first(), own) {
            if most > own {
                parts.push(format!("the strongest signal is {}", s.letter(top)));
            }
        }
    }
    if let Some((from, to)) = s.kept {
        let kept = (from..to).contains(&p.at);
        parts.push(if kept { "kept" } else { "trimmed off" }.to_string());
    }
    parts.join(" · ")
}

fn stretch_words(st: &Stretch, s: Shown) -> String {
    let count = st.to - st.from;
    // mirrored, the last base of the file is the first of the picture
    let (a, b) = (s.number(st.from), s.number(st.to - 1));
    let mut parts = vec![format!(
        "Bases {}–{} of {} ({})",
        commify(a.min(b)),
        commify(a.max(b)),
        commify(s.bases),
        commify(count)
    )];
    match st.mean {
        Some(mean) => parts.push(format!("mean quality {mean:.1}")),
        None => parts.push("no quality in the file".to_string()),
    }
    if let Some((quality, at)) = st.lowest {
        parts.push(format!(
            "lowest {quality} at base {}",
            commify(s.number(at))
        ));
    }
    match st.poor {
        Some(0) => parts.push(format!("none below quality {POOR}")),
        Some(poor) => parts.push(format!("{} below quality {POOR}", commify(poor))),
        None => {}
    }
    if st.unsure > 0 {
        parts.push(format!("{} not A, C, G or T", commify(st.unsure)));
    }
    if let Some((from, to)) = s.kept {
        let kept = st.to.min(to).saturating_sub(st.from.max(from));
        parts.push(if kept == count {
            "all kept".to_string()
        } else if kept == 0 {
            "all trimmed off".to_string()
        } else {
            format!("{} kept by the trimming", commify(kept))
        });
    }
    parts.join(" · ")
}

/// What is picked, in a line.
pub fn words(picked: &Picked, shown: Shown) -> String {
    match picked {
        Picked::One(peak) => peak_words(peak, shown),
        Picked::Several(stretch) => stretch_words(stretch, shown),
    }
}

impl SignalsPage {
    /// Pick bases `from` up to `to` (not counted), as a drag across them does; none lets the pick go.
    pub(in crate::ui) fn pick(&mut self, bases: Option<(usize, usize)>) {
        self.selection = bases;
    }

    /// What the line says of the bases that are picked; nothing when none are.
    pub(in crate::ui) fn picked_words(&self) -> Option<String> {
        let trace = &self.loaded.as_ref()?.trace;
        let (from, to) = self.selection?;
        let shown = Shown {
            bases: trace.bases.len(),
            reverse: self.reverse,
            kept: self.kept(),
        };
        picked(trace, from, to).map(|picked| words(&picked, shown))
    }
}

#[cfg(test)]
mod tests;
