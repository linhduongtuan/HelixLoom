//! The Alignment page as the session keeps it: the list of sequences with their ticks, the choices of the page, and the alignment
//! that is shown with the record of what made it. What is kept is put back at the next start, when the person says so.

use super::{AlignmentPage, Engine, Input, Score};
use crate::kept::{fields, line, value, word};
use crate::pipeline::steps::Aligned;
use crate::provenance::Provenance;
use crate::seq::clustal::{Guide, Options, MOST_REFINE};
use crate::seq::codonalign::Codons;
use crate::seq::dna::GENETIC_CODES;
use crate::ui::platform::alnview::Colouring;
use std::hash::{Hash, Hasher};

const COLOURINGS: [(Colouring, &str); 3] = [
    (Colouring::ByType, "type"),
    (Colouring::Differences, "differences"),
    (Colouring::Plain, "plain"),
];
const SCORES: [(Score, &str); 3] = [
    (Score::Identity, "identity"),
    (Score::Entropy, "entropy"),
    (Score::Jsd, "jsd"),
];

/// Said on the page when the list of sequences came back and the alignment that was kept with it did not.
const ALIGNMENT_NOT_READ: &str =
    "The alignment that was kept cannot be read (its rows are empty or not of one length); the sequences are back.";

/// The settings of the progressive method as the kept text writes them: `GUIDE:PASSES:FAST:FREE`, the guide tree by its word, the
/// passes of refinement, and 1 or 0 for fast distances and for free end gaps.
fn options_code(o: &Options) -> String {
    format!(
        "{}:{}:{}:{}",
        o.guide.name(),
        o.refine,
        u8::from(o.fast),
        u8::from(o.free_end_gaps)
    )
}

/// The settings from their code; none when the code is not one (a word that is not a guide tree, more passes than the program
/// allows, anything but 0 or 1 for the two switches).
fn options_of(code: &str) -> Option<Options> {
    let parts: Vec<&str> = code.split(':').collect();
    let [guide, refine, fast, free] = parts.as_slice() else {
        return None;
    };
    let switch = |s: &str| match s {
        "0" => Some(false),
        "1" => Some(true),
        _ => None,
    };
    Some(Options {
        guide: Guide::from_name(guide)?,
        refine: refine.parse().ok().filter(|r| *r <= MOST_REFINE)?,
        fast: switch(fast)?,
        free_end_gaps: switch(free)?,
    })
}

/// An engine as the kept text writes it: `builtin`, `mafft`, or `progressive:` and the settings (so that a record that is kept can
/// be run again with them).
fn engine_code(engine: Engine) -> String {
    match engine {
        Engine::BuiltIn => "builtin".to_string(),
        Engine::Mafft => "mafft".to_string(),
        Engine::Progressive(o) => format!("progressive:{}", options_code(&o)),
    }
}

fn engine_of(code: &str) -> Option<Engine> {
    match code {
        "builtin" => Some(Engine::BuiltIn),
        "mafft" => Some(Engine::Mafft),
        other => other
            .strip_prefix("progressive:")
            .and_then(options_of)
            .map(Engine::Progressive),
    }
}

fn letters(bases: &[u8]) -> String {
    String::from_utf8_lossy(bases).into_owned()
}

impl AlignmentPage {
    /// A number that is another one whenever what `kept` gives has changed; none when the page holds nothing to keep. It is asked
    /// at every frame, so it reads no sequence: the list counts its own changes (`listed`), and an alignment its own (`revision`).
    pub(in crate::ui) fn stamp(&self) -> Option<u64> {
        if self.inputs.is_empty() && self.result.is_none() {
            return None;
        }
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        (self.listed, self.revision).hash(&mut hasher);
        (self.result.is_some(), self.made.is_some()).hash(&mut hasher);
        (self.engine_now(), self.progressive, self.codons).hash(&mut hasher);
        word(&COLOURINGS, self.colouring).hash(&mut hasher);
        word(&SCORES, self.score).hash(&mut hasher);
        for input in &self.inputs {
            input.keep.hash(&mut hasher);
        }
        Some(hasher.finish())
    }

    /// The page as text that `put_back` reads.
    pub(in crate::ui) fn kept(&self) -> String {
        let mut text = line(&["engine", &engine_code(self.engine_now())]);
        text += &line(&["progressive", &options_code(&self.progressive)]);
        if let Some(c) = self.codons {
            text += &line(&["codons", &c.frame.to_string(), &c.code.to_string()]);
        }
        text += &line(&["colouring", word(&COLOURINGS, self.colouring)]);
        text += &line(&["score", word(&SCORES, self.score)]);
        for input in &self.inputs {
            let ticked = if input.keep { "1" } else { "0" };
            text += &line(&["input", ticked, &input.name, &letters(&input.bases)]);
        }
        let Some(aligned) = &self.result else {
            return text;
        };
        text += &line(&[
            "aligned",
            &aligned.engine,
            &format!("{:?}", aligned.seconds),
        ]);
        for (name, row) in &aligned.rows {
            text += &line(&["row", name, &letters(row)]);
        }
        if let Some((made, engine)) = &self.made {
            text += &line(&["made", &engine_code(*engine), &made.kept()]);
        }
        text
    }

    /// Show what `text` (of `kept`) holds, in place of what the page shows (an alignment that is being made is given up with it). A
    /// line that this version does not know is passed over. An alignment whose rows are empty or not of one length is not shown:
    /// the list of sequences comes back without it, and the page says so. A text with nothing of the page in it changes nothing.
    pub(in crate::ui) fn put_back(&mut self, text: &str) -> Result<(), String> {
        let mut back = AlignmentPage::default();
        let (mut aligned, mut made) = (None::<Aligned>, None);
        for kept in text.lines() {
            let kept = fields(kept);
            let kept: Vec<&str> = kept.iter().map(String::as_str).collect();
            match kept.as_slice() {
                ["engine", name] => {
                    back.engine = engine_of(name).unwrap_or_default();
                    if let Engine::Progressive(o) = back.engine {
                        back.progressive = o;
                    }
                }
                ["progressive", code] => {
                    // a line that means nothing leaves what the engine line gave
                    if let Some(options) = options_of(code) {
                        back.progressive = options;
                    }
                }
                ["codons", frame, code] => {
                    // a choice that is none of the page's leaves it as it was
                    back.codons = frame
                        .parse::<usize>()
                        .ok()
                        .zip(code.parse::<usize>().ok())
                        .filter(|&(f, c)| f <= 2 && c < GENETIC_CODES.len())
                        .map(|(frame, code)| Codons { frame, code });
                }
                ["colouring", name] => {
                    back.colouring = value(&COLOURINGS, name).unwrap_or(back.colouring)
                }
                ["score", name] => back.score = value(&SCORES, name).unwrap_or_default(),
                ["input", ticked, name, bases] => back.inputs.push(Input {
                    name: name.to_string(),
                    bases: bases.as_bytes().to_vec(),
                    keep: *ticked == "1",
                }),
                ["aligned", engine, seconds] => {
                    let seconds = seconds.parse().ok();
                    aligned = Some(Aligned {
                        rows: Vec::new(),
                        engine: engine.to_string(),
                        seconds: seconds
                            .filter(|s: &f64| s.is_finite() && *s >= 0.0)
                            .unwrap_or(0.0),
                    });
                }
                ["row", name, row] => {
                    if let Some(aligned) = aligned.as_mut() {
                        aligned
                            .rows
                            .push((name.to_string(), row.as_bytes().to_vec()));
                    }
                }
                ["made", engine, record] => {
                    made = Provenance::from_kept(record).zip(engine_of(engine))
                }
                _ => {}
            }
        }
        let aligned = aligned.filter(|a| !a.rows.is_empty());
        let columns = aligned.as_ref().map_or(0, Aligned::columns);
        let unreadable = aligned
            .as_ref()
            .is_some_and(|a| columns == 0 || a.rows.iter().any(|r| r.1.len() != columns));
        if unreadable && back.inputs.is_empty() {
            return Err(
                "The alignment that was kept cannot be read: its rows are empty or not of one length."
                    .into(),
            );
        }
        if back.inputs.is_empty() && aligned.is_none() {
            return Err("Nothing of the alignment that was kept can be read.".into());
        }
        // the list does not wait for its alignment
        let aligned = aligned.filter(|_| !unreadable);
        if unreadable {
            back.message = ALIGNMENT_NOT_READ.to_string();
        }
        if let Some(aligned) = aligned {
            back.set_result(aligned);
            // the record is of the alignment: without the one there is not the other
            back.made = made;
        }
        // the counts go on from where they were: what was made from the page before this is out of date, not up to date again
        (back.listed, back.revision) = (self.listed + 1, self.revision + 1);
        *self = back;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::provenance::{Asked, Input as Given};

    fn record() -> Provenance {
        Asked {
            analysis: "alignment",
            inputs: vec![Given::of("one\twith a tab", b"ACGT")],
            parameters: vec![("engine".to_string(), "MAFFT".to_string())],
            ran_on: "this computer",
        }
        .made(vec!["MAFFT v7.526".to_string()], 1.25)
    }

    fn aligned() -> Aligned {
        Aligned {
            rows: vec![
                ("one\twith a tab".to_string(), b"ACGT-A".to_vec()),
                ("two".to_string(), b"AC-TTA".to_vec()),
            ],
            engine: "MAFFT v7.526".to_string(),
            seconds: 1.25,
        }
    }

    /// A page with every choice away from where it starts.
    fn page() -> AlignmentPage {
        let mut page = AlignmentPage::default();
        page.add("one\twith a tab", b"ACGTA");
        page.add("two", b"ACTTA");
        page.add("left\nout", b"GGGG");
        page.inputs[2].keep = false;
        (page.engine, page.colouring, page.score) =
            (Engine::Mafft, Colouring::Differences, Score::Entropy);
        page.set_result(aligned());
        page.made = Some((record(), Engine::Mafft));
        page
    }

    #[test]
    fn the_page_comes_back_as_it_was_kept_with_its_ticks_its_choices_its_alignment_and_the_record()
    {
        let text = page().kept();
        let mut back = AlignmentPage::default();
        back.put_back(&text).expect("put back");
        let listed: Vec<(&str, &[u8], bool)> = back
            .inputs
            .iter()
            .map(|i| (i.name.as_str(), &i.bases[..], i.keep))
            .collect();
        assert_eq!(
            listed,
            [
                ("one\twith a tab", &b"ACGTA"[..], true),
                ("two", &b"ACTTA"[..], true),
                ("left\nout", &b"GGGG"[..], false)
            ]
        );
        assert!(back.engine == Engine::Mafft && back.colouring == Colouring::Differences);
        assert_eq!(back.score, Score::Entropy);
        let shown = back.result.as_ref().expect("an alignment");
        assert_eq!(shown.rows, aligned().rows);
        assert_eq!(
            (shown.engine.as_str(), shown.seconds),
            ("MAFFT v7.526", 1.25)
        );
        // what the viewer draws is worked out anew, with the score that was chosen
        let stats = back.stats.as_ref().expect("the columns");
        assert_eq!((stats.columns.len(), stats.bars.score), (6, Score::Entropy));
        let (made, engine) = back.made.as_ref().expect("the record");
        assert!(*made == record() && *engine == Engine::Mafft);
        // and kept again it is the same text
        assert_eq!(back.kept(), text);
    }

    #[test]
    fn a_list_with_no_alignment_is_kept_too_and_a_record_is_not_kept_without_its_alignment() {
        let mut page = AlignmentPage::default();
        page.add("only", b"ACGT");
        page.made = Some((record(), Engine::BuiltIn));
        let text = page.kept();
        assert!(
            !text.contains("made") && !text.contains("aligned"),
            "{text}"
        );
        let mut back = AlignmentPage::default();
        back.put_back(&text).expect("put back");
        assert_eq!(back.inputs.len(), 1);
        assert!(back.result.is_none() && back.made.is_none());
        // a record whose alignment is not in the text is not put on the page
        let stray = format!("{text}{}", line(&["made", "mafft", &record().kept()]));
        back.put_back(&stray).expect("put back");
        assert!(back.made.is_none());
    }

    #[test]
    fn the_stamp_is_another_one_after_every_change_that_is_kept_and_the_same_one_otherwise() {
        assert_eq!(AlignmentPage::default().stamp(), None);
        let mut page = page();
        let mut seen = vec![page.stamp().expect("a stamp")];
        assert_eq!(page.stamp(), Some(seen[0]), "nothing has changed");
        let mut changed = |page: &AlignmentPage, what: &str| {
            let stamp = page.stamp().expect("a stamp");
            assert!(!seen.contains(&stamp), "{what}");
            seen.push(stamp);
        };
        page.inputs[0].keep = false;
        changed(&page, "a tick");
        page.engine = Engine::BuiltIn;
        changed(&page, "the engine");
        page.colouring = Colouring::Plain;
        changed(&page, "the colouring");
        page.score = Score::Jsd;
        changed(&page, "the score");
        page.made = None;
        changed(&page, "the record");
        page.set_result(aligned());
        changed(&page, "an alignment");
        page.add("more", b"TTTT");
        changed(&page, "one more sequence");
        // the list emptied and filled again with sequences of the same names and lengths is not the list it was
        let again: Vec<Input> = page.inputs.clone();
        page.clear_inputs();
        for input in &again {
            page.add(&input.name, &input.bases);
            page.inputs.last_mut().expect("it").keep = input.keep;
        }
        changed(&page, "the same names over again");
        // what the page only shows for the moment is not kept, and does not move the stamp
        let before = page.stamp();
        (page.message, page.paste) = ("said".to_string(), ">typed".to_string());
        (page.selected, page.reveal) = (Some(2), true);
        assert_eq!(page.stamp(), before);
    }

    #[test]
    fn a_text_that_is_not_an_alignment_changes_nothing_and_a_line_from_a_later_version_is_passed_over(
    ) {
        let mut shown = page();
        let before = shown.kept();
        for (text, why) in [
            ("", "Nothing of the alignment"),
            ("engine\tmafft\nscore\tjsd\n", "Nothing of the alignment"),
            (
                "aligned\tMAFFT\t1\nrow\ta\tACGT\nrow\tb\tAC\n",
                "empty or not of one length",
            ),
            (
                "aligned\tMAFFT\t1\nrow\ta\t\nrow\tb\t\n",
                "empty or not of one length",
            ),
        ] {
            let said = shown.put_back(text).expect_err(text);
            assert!(said.contains(why), "{said}");
            assert_eq!(shown.kept(), before);
        }
        // words and numbers that mean nothing are the ones the page starts with; rows before their alignment are no rows
        let odd = "engine\tclustal\ncolouring\tneon\nscore\tbest\nfrom a later version\twith\tfields\nrow\tstray\tACGT\n\
                   input\t1\ta\tACGT\naligned\tX\tnot a number\nrow\ta\tACGT\nmade\tmafft\tnot a record\n";
        shown.put_back(odd).expect("put back");
        assert!(shown.engine == Engine::BuiltIn && shown.colouring == Colouring::ByType);
        assert_eq!(shown.score, Score::Identity);
        let aligned = shown.result.as_ref().expect("an alignment");
        assert_eq!((aligned.rows.len(), aligned.seconds), (1, 0.0));
        assert!(shown.made.is_none());
    }

    #[test]
    fn a_list_comes_back_although_its_alignment_cannot_be_read_and_the_page_says_so() {
        let mut shown = AlignmentPage::default();
        let kept = "engine\tmafft\ninput\t1\ta\tACGT\ninput\t0\tb\tAC\naligned\tMAFFT\t1\nrow\ta\tACGT\nrow\tb\tAC\n";
        shown.put_back(kept).expect("the list is read");
        assert_eq!(shown.inputs.len(), 2);
        assert!(shown.engine == Engine::Mafft && !shown.inputs[1].keep);
        assert!(shown.result.is_none() && shown.made.is_none());
        assert_eq!(shown.message, ALIGNMENT_NOT_READ);
        // one that is read whole says nothing
        let whole = "input\t1\ta\tACGT\naligned\tMAFFT\t1\nrow\ta\tACGT\n";
        shown.put_back(whole).expect("put back");
        assert!(shown.result.is_some() && shown.message.is_empty());
    }

    #[test]
    fn a_time_that_is_no_time_is_none() {
        for seconds in ["-5", "inf", "NaN", "soon"] {
            let mut shown = AlignmentPage::default();
            let kept = format!("aligned\tMAFFT\t{seconds}\nrow\ta\tACGT\n");
            shown.put_back(&kept).expect("put back");
            let aligned = shown.result.as_ref().expect("an alignment");
            assert_eq!(aligned.seconds, 0.0, "{seconds}");
        }
    }

    #[test]
    fn putting_back_moves_the_counts_on_so_that_what_was_made_from_the_page_before_is_out_of_date()
    {
        let mut shown = page();
        let (stamp, revision) = (shown.stamp(), shown.revision);
        let text = shown.kept();
        shown.put_back(&text).expect("put back");
        assert!(shown.revision > revision);
        assert_ne!(shown.stamp(), stamp);
    }

    /// Settings of the progressive method that are none of the defaults.
    fn settled() -> Options {
        Options {
            guide: Guide::Upgma,
            free_end_gaps: false,
            refine: 3,
            fast: true,
        }
    }

    #[test]
    fn the_progressive_method_and_its_settings_come_back_as_they_were_kept_with_the_record_that_names_them(
    ) {
        let mut page = page();
        (page.engine, page.progressive) = (Engine::Progressive(Options::default()), settled());
        page.made = Some((record(), Engine::Progressive(settled())));
        let text = page.kept();
        assert!(text.contains("engine\tprogressive:upgma:3:1:0\n"), "{text}");
        let mut back = AlignmentPage::default();
        back.put_back(&text).expect("put back");
        // the page shows the settings that were kept, and the engine is asked for with them
        assert_eq!(back.progressive, settled());
        assert!(back.engine_now() == Engine::Progressive(settled()));
        let (_, engine) = back.made.as_ref().expect("the record");
        assert!(*engine == Engine::Progressive(settled()));
        assert_eq!(back.kept(), text);
    }

    #[test]
    fn the_settings_of_the_progressive_method_are_kept_while_another_engine_is_the_choice() {
        let mut page = page();
        page.progressive = settled();
        assert!(page.engine == Engine::Mafft);
        let text = page.kept();
        let mut back = AlignmentPage::default();
        back.put_back(&text).expect("put back");
        assert!(back.engine == Engine::Mafft);
        assert_eq!(back.progressive, settled());
        // so that choosing the method again shows what was set
        back.engine = Engine::Progressive(Options::default());
        assert!(back.engine_now() == Engine::Progressive(settled()));
    }

    #[test]
    fn settings_that_mean_nothing_are_the_ones_the_page_starts_with() {
        for code in [
            "progressive:nj:11:0:1", // more passes than the program makes
            "progressive:nj:-1:0:1",
            "progressive:nj:x:0:1",
            "progressive:tree:0:0:1", // not a guide tree
            "progressive:nj:0:2:1",   // a switch is 0 or 1
            "progressive:nj:0:0:yes",
            "progressive:nj:0:0",     // one field short
            "progressive:nj:0:0:1:9", // one too many
            "progressive:",
            "progressive",
        ] {
            let mut shown = AlignmentPage::default();
            let kept = format!("engine\t{code}\nprogressive\t{code}\ninput\t1\ta\tACGT\n");
            shown.put_back(&kept).expect("put back");
            assert!(shown.engine == Engine::BuiltIn, "{code}");
            assert_eq!(shown.progressive, Options::default(), "{code}");
        }
        // the largest number of passes the program makes is read as it is
        let mut shown = AlignmentPage::default();
        let most = format!("engine\tprogressive:nj:{MOST_REFINE}:0:1\ninput\t1\ta\tACGT\n");
        shown.put_back(&most).expect("put back");
        assert_eq!(shown.progressive.refine, MOST_REFINE);
    }

    #[test]
    fn a_progressive_line_that_means_nothing_leaves_the_settings_the_engine_line_gave() {
        let mut shown = AlignmentPage::default();
        let kept = "engine\tprogressive:upgma:3:1:0\nprogressive\tgarbage\ninput\t1\ta\tACGT\n";
        shown.put_back(kept).expect("put back");
        assert_eq!(shown.progressive, settled());
        assert!(shown.engine_now() == Engine::Progressive(settled()));
    }

    #[test]
    fn the_stamp_follows_the_settings_of_the_progressive_method() {
        let mut page = page();
        page.engine = Engine::Progressive(Options::default());
        let mut seen = vec![page.stamp().expect("a stamp")];
        for (what, set) in [
            (
                "the guide tree",
                Options {
                    guide: Guide::Upgma,
                    ..Options::default()
                },
            ),
            (
                "the passes",
                Options {
                    refine: 2,
                    ..Options::default()
                },
            ),
            (
                "fast distances",
                Options {
                    fast: true,
                    ..Options::default()
                },
            ),
            (
                "the end gaps",
                Options {
                    free_end_gaps: false,
                    ..Options::default()
                },
            ),
        ] {
            page.progressive = set;
            let stamp = page.stamp().expect("a stamp");
            assert!(!seen.contains(&stamp), "{what}");
            seen.push(stamp);
        }
        // settings that are kept are kept whichever engine is the choice
        page.engine = Engine::Mafft;
        page.progressive = Options {
            refine: 7,
            ..Options::default()
        };
        let stamp = page.stamp().expect("a stamp");
        assert!(!seen.contains(&stamp));
    }
}
