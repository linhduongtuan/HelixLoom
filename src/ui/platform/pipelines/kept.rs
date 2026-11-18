//! The Pipelines page as the session keeps it: the run that is shown in each of its two modes (the fixed route, and a pipeline
//! file run on a folder), with what was ticked and the record of what made it. A run is kept in the text it already has for
//! the workstation (`remote::tasks`, `remote::flow`), which is read back with the same care as an answer from the network; so a
//! name in it comes back as that text writes it (on one line, without spaces at its ends), and a run that the text cannot hold
//! (a file with no name, more files than a job may have) is not put back.
//!
//! The folder of a file run and the folder of its cache are not kept: a path is not taken from a kept text, and nothing is
//! listed when a session is put back. They are chosen again to run again.

use super::flow::Mode;
use super::{Kept, PipelinesPage, Source};
use crate::kept::{fields, line, value, word};
use crate::pipeline::Outcome;
use crate::provenance::Provenance;
use crate::remote::flow::FlowResult;
use crate::remote::tasks::{decode_pipeline_response, encode_pipeline_response};
use crate::remote::wire::Container;
use std::hash::{Hash, Hasher};

const MODES: [(Mode, &str); 2] = [(Mode::Route, "route"), (Mode::File, "file")];
const SOURCES: [(Source, &str); 2] = [(Source::Open, "open"), (Source::Example, "example")];
/// Where a run of a pipeline file was done, as the page says it.
const PLACES: [&str; 2] = ["this computer", "the workstation"];

/// Said on the page when the run of the last session is put back.
pub const RUN_IS_BACK: &str = "The run of your last session.";
/// Said when a run is going or has been asked for: what was kept is not put over it (the run would take what was put back away
/// when it starts, or its answer would come over it).
const RUNNING: &str =
    "A run is going, or about to start, on the Pipelines page: the run that was kept is not put over it.";
const NOTHING: &str = "Nothing of the pipeline run that was kept can be read.";

fn ticks_text(ticked: [bool; 5]) -> String {
    ticked.iter().map(|t| if *t { '1' } else { '0' }).collect()
}

/// What a run was ticked with, as it was kept: a source that is one, and five boxes each `0` or `1`.
fn ticks_of(source: &str, ticked: &str) -> Option<Kept> {
    let boxes: Option<Vec<bool>> = ticked
        .chars()
        .map(|c| match c {
            '1' => Some(true),
            '0' => Some(false),
            _ => None,
        })
        .collect();
    Some((value(&SOURCES, source)?, boxes?.try_into().ok()?))
}

/// A container as text; none when a section of it is not text (a run's are).
fn text_of(container: &Container) -> Option<String> {
    String::from_utf8(container.encode()).ok()
}

/// A run of the fixed route as it was kept. Its alignment is read by the pages that take it, so its rows are of one length.
fn route_of(run: &str) -> Result<Outcome, String> {
    let outcome = decode_pipeline_response(&Container::decode(run.as_bytes())?)?;
    let ragged = outcome
        .aligned
        .as_ref()
        .is_some_and(|a| a.rows.iter().any(|row| row.1.len() != a.columns()));
    if ragged {
        return Err("the rows of its alignment are not of one length.".into());
    }
    Ok(outcome)
}

fn file_of(run: &str) -> Result<FlowResult, String> {
    FlowResult::decode(&Container::decode(run.as_bytes())?)
}

impl PipelinesPage {
    /// A number that is another one whenever what `kept_text` gives has changed; none when the page shows no run. Asked at every
    /// frame: it reads no run (each mode counts its own changes, `shown`).
    pub(in crate::ui) fn stamp(&self) -> Option<u64> {
        if self.outcome.is_none() && self.flow.report.is_none() {
            return None;
        }
        let (source, ticked) = self.kept();
        let flow = &self.flow;
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        (word(&MODES, self.mode), word(&SOURCES, source), ticked).hash(&mut hasher);
        // of each mode, only what is kept of it: nothing while it shows no run
        let route = self.outcome.as_ref();
        route
            .map(|_| (self.shown, self.made.is_some()))
            .hash(&mut hasher);
        let file = flow.report.as_ref();
        file.map(|_| (flow.shown, flow.made.is_some(), flow.ran_on))
            .hash(&mut hasher);
        Some(hasher.finish())
    }

    /// The runs that the page shows as text that `put_back` reads; none when it shows no run.
    pub(in crate::ui) fn kept_text(&self) -> Option<String> {
        self.stamp()?;
        let (source, ticked) = self.kept();
        let mut text = line(&["mode", word(&MODES, self.mode)]);
        text += &line(&["ticks", word(&SOURCES, source), &ticks_text(ticked)]);
        let route = self.outcome.as_ref().map(encode_pipeline_response);
        if let Some(run) = route.as_ref().and_then(text_of) {
            text += &line(&["route", &run]);
            if let Some((made, (source, ticked))) = &self.made {
                let (source, ticked) = (word(&SOURCES, *source), ticks_text(*ticked));
                text += &line(&["route_made", source, &ticked, &made.kept()]);
            }
        }
        let flow = &self.flow;
        let file = flow.report.as_ref().map(|report| FlowResult {
            pipeline: flow.pipeline.clone(),
            report: report.clone(),
            // the time is in the record; the text only has to have one
            seconds: 0.0,
        });
        if let Some(run) = file.and_then(|result| text_of(&result.encode())) {
            text += &line(&["file", flow.ran_on, &run]);
            if let Some((made, _)) = &flow.made {
                text += &line(&["file_made", &made.kept()]);
            }
        }
        Some(text)
    }

    /// Show the runs that `text` (of `kept_text`) holds, in place of the ones the page shows. A run that cannot be read is not
    /// shown, and its mode of the page says why; when none can be read, or a run is going, nothing changes and why is said.
    pub(in crate::ui) fn put_back(&mut self, text: &str) -> Result<(), String> {
        if self.is_running() || self.run_requested {
            return Err(RUNNING.into());
        }
        let (mut mode, mut ticks) = (None, None);
        let (mut route, mut route_made) = (None, None);
        let (mut file, mut file_made) = (None, None);
        for kept in text.lines() {
            let kept = fields(kept);
            let kept: Vec<&str> = kept.iter().map(String::as_str).collect();
            match kept.as_slice() {
                ["mode", name] => mode = value(&MODES, name),
                ["ticks", source, ticked] => ticks = ticks_of(source, ticked),
                ["route", run] => route = Some(route_of(run)),
                ["route_made", source, ticked, record] => {
                    route_made = Provenance::from_kept(record).zip(ticks_of(source, ticked))
                }
                ["file", ran_on, run] => {
                    let place = PLACES.iter().find(|place| *place == ran_on);
                    file = Some(file_of(run).map(|result| (place.copied().unwrap_or(""), result)));
                }
                ["file_made", record] => file_made = Provenance::from_kept(record),
                _ => {}
            }
        }
        if !matches!(route, Some(Ok(_))) && !matches!(file, Some(Ok(_))) {
            let why = [route.and_then(Result::err), file.and_then(Result::err)];
            let why = why.into_iter().flatten().next();
            return Err(why.map_or(NOTHING.to_string(), |why| format!("{NOTHING} {why}")));
        }
        let not_read = |why: String| format!("The run that was kept cannot be read: {why}");
        self.mode = mode.unwrap_or(self.mode);
        if let Some(ticks) = ticks {
            self.restore(ticks);
        }
        match route {
            Some(Ok(outcome)) => {
                self.set_outcome(Some(outcome));
                (self.made, self.message) = (route_made, RUN_IS_BACK.to_string());
            }
            unread => {
                self.set_outcome(None);
                self.made = None;
                self.message = unread.and_then(Result::err).map_or(String::new(), not_read);
            }
        }
        let page = &mut self.flow;
        match file {
            Some(Ok((ran_on, result))) => {
                (page.pipeline, page.ran_on) = (result.pipeline, ran_on);
                page.set_report(Some(result.report));
                page.made = file_made.map(|made| (made, page.pipeline.clone()));
                page.message = RUN_IS_BACK.to_string();
                // the folder that is chosen, if one is, holds other files for another pipeline
                page.count_files();
            }
            unread => {
                page.forget_run();
                page.message = unread.and_then(Result::err).map_or(String::new(), not_read);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::steps::{Aligned, PrimerPair, Primers, Tree};
    use crate::provenance::{Asked, Input};
    use crate::ui::app::{App, Startup};
    use crate::ui::job::Job;
    use eframe::egui;
    use std::path::PathBuf;

    fn record(analysis: &'static str) -> Provenance {
        Asked {
            analysis,
            inputs: vec![Input::of("Isolate 1", b"ACGTACGT")],
            parameters: vec![("steps".to_string(), "align, tree, primers".to_string())],
            ran_on: "this computer",
        }
        .made(vec!["MAFFT v7.526".to_string()], 3.5)
    }

    /// A run of the fixed route with every part.
    fn outcome() -> Outcome {
        Outcome {
            log: vec![
                "Align: 2 sequences.".to_string(),
                "Tree: built in.".to_string(),
            ],
            aligned: Some(Aligned {
                rows: vec![
                    ("Isolate 1".to_string(), b"ACGT-ACGT".to_vec()),
                    ("Isolate 2".to_string(), b"ACGTTACGT".to_vec()),
                ],
                engine: "MAFFT v7.526".to_string(),
                seconds: 0.25,
            }),
            tree: Some(Tree {
                newick: "(Isolate_1:0.1,Isolate_2:0.2);".to_string(),
                engine: "IQ-TREE 2.3".to_string(),
                model: "HKY+F".to_string(),
                ln_likelihood: Some(-1234.5678),
                seconds: 1.5,
            }),
            primers: Some(Primers {
                pairs: vec![PrimerPair {
                    left: "ACGTACGTACGTACGTACGT".to_string(),
                    right: "TTGCATTGCATTGCATTGCA".to_string(),
                    left_tm: 59.87,
                    right_tm: 60.12,
                    left_start: 12,
                    right_start: 240,
                    product_size: 248,
                    penalty: 0.3125,
                }],
                engine: "primer3 2.6.1".to_string(),
                seconds: 0.03,
            }),
        }
    }

    /// The parts of a run as a test can compare them (the types of a run are not compared whole).
    fn parts(o: &Outcome) -> (Vec<String>, String, String, String) {
        let aligned = o.aligned.as_ref().map(|a| (&a.rows, &a.engine, a.seconds));
        let tree = o.tree.as_ref();
        let tree = tree.map(|t| (&t.newick, &t.engine, &t.model, t.ln_likelihood, t.seconds));
        let primers = o.primers.as_ref().map(|p| (&p.pairs, &p.engine, p.seconds));
        (
            o.log.clone(),
            format!("{aligned:?}"),
            format!("{tree:?}"),
            format!("{primers:?}"),
        )
    }

    /// A page that shows a run of the route, with every box away from where it starts.
    fn route_page() -> PipelinesPage {
        let mut page = PipelinesPage::default();
        let ticked = (Source::Example, [false, true, false, false, true]);
        page.restore(ticked);
        page.set_outcome(Some(outcome()));
        page.made = Some((
            record("Pipeline"),
            (Source::Open, [true, true, false, true, false]),
        ));
        page
    }

    /// An app whose Pipelines page shows the run of the example pipeline file, with its record and its folders.
    fn file_app() -> App {
        let mut app = App::new(Startup::default());
        app.show_example_flow_run();
        let page = &mut app.platform.pipelines;
        page.mode = Mode::File;
        let flow = &mut page.flow;
        assert!(flow.report.is_some(), "{}", flow.message);
        flow.cache = Some(PathBuf::from("a cache\twith a tab"));
        flow.made = Some((record("Pipeline file"), flow.pipeline.clone()));
        app
    }

    #[test]
    fn a_run_of_the_route_comes_back_with_its_parts_what_was_ticked_and_its_record() {
        let shown = route_page();
        let text = shown.kept_text().expect("a run");
        let mut back = PipelinesPage::default();
        back.put_back(&text).expect("put back");
        assert_eq!(back.mode, Mode::Route);
        assert_eq!(
            back.kept(),
            (Source::Example, [false, true, false, false, true])
        );
        let run = back.outcome.as_ref().expect("a run");
        assert_eq!(parts(run), parts(&outcome()));
        let (made, ticked) = back.made.as_ref().expect("the record");
        assert_eq!(*made, record("Pipeline"));
        assert_eq!(*ticked, (Source::Open, [true, true, false, true, false]));
        assert_eq!(back.message, RUN_IS_BACK);
        assert!(back.flow.report.is_none() && back.flow.made.is_none());
        assert_eq!(back.kept_text(), Some(text));
        // a run with one part only, and no record (it was fetched from the workstation's list)
        let mut fetched = PipelinesPage::default();
        let only = Outcome {
            primers: outcome().primers,
            ..Outcome::default()
        };
        fetched.show_fetched(only, "fetched".to_string());
        back.put_back(&fetched.kept_text().expect("a run"))
            .expect("put back");
        let run = back.outcome.as_ref().expect("a run");
        assert!(run.aligned.is_none() && run.tree.is_none() && run.primers.is_some());
        assert!(back.made.is_none());
    }

    #[test]
    fn a_run_of_a_pipeline_file_comes_back_with_its_pipeline_its_report_where_it_ran_and_its_record(
    ) {
        let app = file_app();
        let shown = &app.platform.pipelines;
        let text = shown.kept_text().expect("a run");
        let mut back = PipelinesPage::default();
        // another pipeline is on the page: the one of the run takes its place
        back.flow.pipeline.name = "another".to_string();
        back.put_back(&text).expect("put back");
        assert_eq!(back.mode, Mode::File);
        let (flow, was) = (&back.flow, &shown.flow);
        assert_eq!(flow.pipeline, was.pipeline);
        assert_eq!(flow.report, was.report);
        assert_eq!(flow.ran_on, "this computer");
        // the folders are not kept: a path is not taken from a kept text
        assert!(was.folder.is_some() && was.cache.is_some());
        assert!(flow.folder.is_none() && flow.cache.is_none());
        assert!(
            !text.contains("example-loci") && !text.contains("a cache"),
            "{text}"
        );
        let (made, pipeline) = flow.made.as_ref().expect("the record");
        assert!(*made == record("Pipeline file") && *pipeline == was.pipeline);
        assert_eq!(flow.message, RUN_IS_BACK);
        assert!(back.outcome.is_none());
        assert_eq!(back.kept_text(), Some(text));
        // the folders that are chosen on the page stay as they are, whatever a kept text says of folders
        let chosen = (
            Some(PathBuf::from("chosen")),
            Some(PathBuf::from("its cache")),
        );
        (back.flow.folder, back.flow.cache) = chosen.clone();
        let named = format!(
            "{}folder\t/elsewhere\ncache\t/elsewhere/cache\n",
            back.kept_text().expect("a run")
        );
        back.put_back(&named).expect("put back");
        assert_eq!((back.flow.folder.clone(), back.flow.cache.clone()), chosen);
        // a place that is no place of a run is none, and a run with no record comes back so
        let mut bare = file_app();
        let flow = &mut bare.platform.pipelines.flow;
        (flow.made, flow.ran_on) = (None, "the workstation");
        let text = bare.platform.pipelines.kept_text().expect("a run");
        back.put_back(&text).expect("put back");
        assert!(back.flow.made.is_none());
        assert_eq!(back.flow.ran_on, "the workstation");
        back.put_back(&text.replace("the workstation", "the moon"))
            .expect("put back");
        assert_eq!(back.flow.ran_on, "");
    }

    #[test]
    fn both_modes_are_kept_together_and_a_page_with_no_run_keeps_nothing() {
        let mut page = PipelinesPage::default();
        page.restore((Source::Example, [true; 5]));
        assert_eq!(
            (page.stamp(), page.kept_text()),
            (None, None),
            "boxes alone are no run"
        );
        let mut app = file_app();
        app.platform.pipelines.set_outcome(Some(outcome()));
        let text = app.platform.pipelines.kept_text().expect("two runs");
        let mut back = PipelinesPage::default();
        back.put_back(&text).expect("put back");
        assert!(back.outcome.is_some() && back.flow.report.is_some());
        assert_eq!(back.kept_text(), Some(text));
    }

    #[test]
    fn the_stamp_is_another_one_after_every_change_that_is_kept_and_the_same_one_otherwise() {
        let mut app = file_app();
        let page = &mut app.platform.pipelines;
        let mut seen = vec![page.stamp().expect("a stamp")];
        assert_eq!(page.stamp(), Some(seen[0]), "nothing has changed");
        let mut changed = |page: &PipelinesPage, what: &str| {
            let stamp = page
                .stamp()
                .unwrap_or_else(|| panic!("no stamp after {what}"));
            assert!(!seen.contains(&stamp), "{what}");
            seen.push(stamp);
        };
        page.set_outcome(Some(outcome()));
        changed(page, "a run of the route");
        page.set_outcome(Some(outcome()));
        changed(page, "the same run again");
        page.made = Some((record("Pipeline"), page.kept()));
        changed(page, "its record");
        page.show_fetched(outcome(), "fetched".to_string());
        changed(page, "a run fetched from the workstation");
        page.mode = Mode::Route;
        changed(page, "the mode");
        page.source = Source::Example;
        changed(page, "the source");
        for tick in 0..5 {
            let (source, mut ticked) = page.kept();
            ticked[tick] = !ticked[tick];
            page.restore((source, ticked));
            changed(page, "a box");
        }
        page.flow.made = None;
        changed(page, "the record of the file's run");
        page.flow.ran_on = "the workstation";
        changed(page, "where the file's run was done");
        let report = page.flow.report.clone();
        page.flow.set_report(report);
        changed(page, "the same report again");
        page.flow.forget_run();
        changed(page, "the file's run forgotten");
        // what the page only says for the moment is not kept, nor are the folders; and of a mode that shows no run, nothing is
        let before = (page.stamp(), page.kept_text());
        (page.message, page.flow.message) = ("said".to_string(), "said".to_string());
        (page.flow.kept, page.flow.size) = (Some(3), Some((4, 5)));
        (page.flow.folder, page.flow.cache) = (Some(PathBuf::from("other")), None);
        assert!(page.flow.report.is_none(), "only the route shows a run now");
        page.flow.ran_on = "this computer";
        page.flow.made = Some((record("Pipeline file"), page.flow.pipeline.clone()));
        page.flow.forget_run();
        assert_eq!((page.stamp(), page.kept_text()), before);
        // the same the other way round: the route's own count does not move the stamp while only the file shows a run
        let report = file_app().platform.pipelines.flow.report.clone();
        page.flow.set_report(report);
        page.set_outcome(None);
        let before = (page.stamp().expect("the file's run"), page.kept_text());
        page.set_outcome(None);
        page.made = None;
        assert_eq!(
            (page.stamp().expect("the file's run"), page.kept_text()),
            before
        );
        page.flow.forget_run();
        assert_eq!(page.stamp(), None, "no run in either mode");
    }

    #[test]
    fn a_run_that_cannot_be_read_is_not_shown_and_the_other_mode_comes_back_without_it() {
        let app = file_app();
        let file = app.platform.pipelines.kept_text().expect("a run");
        let route = route_page().kept_text().expect("a run");
        let mut shown = route_page();
        let before = (shown.stamp(), shown.kept_text());
        // nothing of a run: the page stays as it is, and why is said
        // an alignment whose rows are of two lengths, kept as a run is (a text cut by hand is refused before its rows are read)
        let mut ragged = route_page();
        let mut run = outcome();
        run.aligned.as_mut().expect("an alignment").rows[1].1 = b"ACGTT".to_vec();
        ragged.set_outcome(Some(run));
        let ragged = ragged.kept_text().expect("a run");
        let cut = &route[..route.len() / 2];
        for (text, why) in [
            ("", ""),
            ("mode\tfile\nticks\texample\t11111\n", ""),
            ("route\tnot a container\n", "not a container of sections"),
            (ragged.as_str(), "not of one length"),
            (cut, "runs past the end"),
            ("file\tthis computer\t@meta 3\\nxyz\\n\n", "seconds"),
        ] {
            let said = shown.put_back(text).expect_err(text);
            // with no run in the text there is no more to say; else the reason follows
            match why {
                "" => assert_eq!(said, NOTHING),
                why => assert!(said.starts_with(NOTHING) && said.contains(why), "{said}"),
            }
            assert_eq!((shown.stamp(), shown.kept_text()), before);
        }
        // one mode's run is damaged: the other comes back, and the damaged one's side of the page says so
        let both = format!("{file}route\tnot a container\n");
        shown.put_back(&both).expect("the file's run is read");
        assert!(shown.outcome.is_none() && shown.made.is_none());
        assert!(
            shown
                .message
                .starts_with("The run that was kept cannot be read:"),
            "{}",
            shown.message
        );
        assert!(shown.flow.report.is_some());
        let both = format!("{route}file\tthis computer\tnot a container\n");
        shown.put_back(&both).expect("the route's run is read");
        assert!(
            shown.outcome.is_some() && shown.flow.report.is_none() && shown.flow.made.is_none()
        );
        assert!(shown
            .flow
            .message
            .starts_with("The run that was kept cannot be read:"));
        // words that mean nothing are passed over, and what was ticked stays as it is; a record without its run is none
        let mut plain = PipelinesPage::default();
        let odd = route
            .replace("ticks\texample\t01001", "ticks\tnowhere\t01001")
            .replace("mode\troute", "mode\tsideways\nlater\tline");
        plain.put_back(&odd).expect("put back");
        assert_eq!(plain.kept(), PipelinesPage::default().kept());
        assert_eq!(plain.mode, Mode::Route);
        for ticks in ["0100", "010011", "01a01"] {
            assert_eq!(ticks_of("open", ticks), None, "{ticks}");
        }
        let stray = format!(
            "{file}{}",
            line(&["route_made", "open", "11111", &record("x").kept()])
        );
        plain.put_back(&stray).expect("put back");
        assert!(plain.outcome.is_none() && plain.made.is_none());
    }

    #[test]
    fn a_run_that_was_kept_is_not_put_over_a_run_that_is_going_or_has_been_asked_for() {
        let text = route_page().kept_text().expect("a run");
        let going = Job::spawn(&egui::Context::default(), || {
            Err("never looked at".to_string())
        });
        let mut busy = PipelinesPage {
            job: Some(going),
            ..PipelinesPage::default()
        };
        assert_eq!(busy.put_back(&text), Err(RUNNING.to_string()));
        assert!(busy.outcome.is_none());
        // asked for and not started yet (the programs are still being looked for): the run would take what was put back away,
        // and would be made with the boxes of the kept run and not with the ones that were ticked for it
        let mut asked = PipelinesPage {
            run_requested: true,
            source: Source::Example,
            ..PipelinesPage::default()
        };
        let before = asked.kept();
        assert_eq!(asked.put_back(&text), Err(RUNNING.to_string()));
        assert!(asked.outcome.is_none() && asked.run_requested);
        assert_eq!(asked.kept(), before);
    }

    #[test]
    fn putting_a_run_back_counts_as_a_change_of_the_run_in_each_mode() {
        let mut app = file_app();
        let page = &mut app.platform.pipelines;
        page.set_outcome(Some(outcome()));
        let text = page.kept_text().expect("two runs");
        let without = |word: &str| -> String {
            let lines = text.lines().filter(|l| !l.starts_with(word));
            lines.map(|l| format!("{l}\n")).collect()
        };
        // the same text put back is still another showing: each mode's count goes on
        let (route, file) = (page.shown, page.flow.shown);
        let before = page.stamp();
        page.put_back(&text).expect("put back");
        assert!(page.shown > route && page.flow.shown > file);
        assert_ne!(page.stamp(), before);
        // and a mode that the text has no run for is counted as changed too: its run is gone
        let (route, file) = (page.shown, page.flow.shown);
        page.put_back(&without("file")).expect("put back");
        assert!(page.flow.report.is_none() && page.flow.shown > file && page.shown > route);
        page.put_back(&text).expect("put back");
        let route = page.shown;
        page.put_back(&without("route")).expect("put back");
        assert!(page.outcome.is_none() && page.shown > route);
    }
}
