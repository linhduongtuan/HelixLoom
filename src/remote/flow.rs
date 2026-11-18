//! The pipeline-file job: a pipeline file and the files of a folder go in, what the run did and what it made come out
//! (`pipeline::flow`). Like the other kinds it runs the same function on the workstation and on this computer, and travels
//! as a container of sections: the pipeline as its own text, the files one after another with their names and sizes as
//! lines, and the report as tab-separated lines.
//!
//! On the workstation the answers of the steps are kept in a folder of the server's own (`--data-dir`), so that a run
//! sent again, with one locus changed, does only the steps of that locus and the gathering. The request names no folder.

use super::clean_line;
use super::tasks::{number, MAX_RESIDUES, MAX_SEQUENCES, MAX_TREES};
use super::wire::{meta, Container};
use crate::pipeline::flow::{
    self as engine, Cache, ItemReport, Kind, Pipeline, Report, Status, StepReport, Value,
};
use std::collections::HashSet;
use std::path::Path;
use web_time::Instant;

/// The most files a request may hold: as many as gene trees a species-tree job may hold.
pub const MAX_LOCI: usize = MAX_TREES;
/// The longest pipeline file a request may carry: a pipeline is a few lines.
const MAX_PIPELINE_BYTES: usize = 64 << 10;
/// The most bytes of files a request may hold (twice the letters a job may hold: names and line ends come with them).
const MAX_LOCI_BYTES: usize = 2 * MAX_RESIDUES;

/// What to run: the pipeline and the files it is run on, a name and what it holds each.
#[derive(Clone, Debug, PartialEq)]
pub struct FlowJob {
    pub pipeline: Pipeline,
    pub loci: Vec<(String, Value)>,
}

/// Texts one after another: their names and sizes, one to a line, and the texts as one.
fn pack(items: &[(&str, String)]) -> (String, String, String) {
    let names = items.iter().map(|(n, _)| clean_line(n)).collect::<Vec<_>>();
    let sizes = items
        .iter()
        .map(|(_, t)| t.len().to_string())
        .collect::<Vec<_>>();
    let body = items.iter().map(|(_, t)| t.as_str()).collect::<String>();
    (names.join("\n"), sizes.join("\n"), body)
}

/// The texts that `pack` made, with their names; refused when the sizes do not cut `body` exactly into the names' texts.
fn unpack<'a>(
    names: &'a str,
    sizes: &str,
    body: &'a str,
    most: usize,
) -> Result<Vec<(&'a str, &'a str)>, String> {
    let bad = |why: &str| format!("The files in the message do not add up: {why}.");
    let (names, sizes): (Vec<&str>, Vec<&str>) = (names.lines().collect(), sizes.lines().collect());
    if names.len() != sizes.len() {
        return Err(bad("there are not as many sizes as names"));
    }
    if names.len() > most {
        return Err(format!(
            "A job may hold at most {most} files; this one has {}.",
            names.len()
        ));
    }
    let mut start = 0usize;
    let mut out = Vec::with_capacity(names.len());
    for (name, size) in names.into_iter().zip(sizes) {
        let size: usize = size.parse().map_err(|_| bad("a size is not a number"))?;
        let end = start
            .checked_add(size)
            .ok_or_else(|| bad("a size is too large"))?;
        // `get` is `None` past the end and inside a character
        let text = body
            .get(start..end)
            .ok_or_else(|| bad("a size runs past the end or cuts a character"))?;
        out.push((name, text));
        start = end;
    }
    if start != body.len() {
        return Err(bad("bytes are left over"));
    }
    Ok(out)
}

/// The rows of what a file holds, which a request may not exceed.
fn rows_within_limit(name: &str, value: &Value) -> Result<(), String> {
    if let Value::Sequences(rows) | Value::Alignment(rows) = value {
        if rows.len() > MAX_SEQUENCES {
            return Err(format!(
                "{name} has {} sequences, and a job may hold at most {MAX_SEQUENCES} in a file.",
                rows.len()
            ));
        }
    }
    Ok(())
}

impl FlowJob {
    pub fn encode(&self) -> Container {
        let texts: Vec<(&str, String)> = self
            .loci
            .iter()
            .map(|(name, value)| (name.as_str(), value.to_text()))
            .collect();
        let (names, sizes, body) = pack(&texts);
        Container::new()
            .with("meta", format!("files={}\n", self.loci.len()))
            .with("pipeline", self.pipeline.to_text())
            .with("names", names)
            .with("sizes", sizes)
            .with("loci", body)
    }

    /// A request comes from the network: the pipeline is read and checked, the files are cut out by their sizes and each is
    /// read as what the pipeline takes, all within what a job may hold.
    pub fn decode(c: &Container) -> Result<FlowJob, String> {
        let text = c.require("pipeline")?;
        if text.len() > MAX_PIPELINE_BYTES {
            return Err(
                "The pipeline file in the message is far longer than a pipeline is.".into(),
            );
        }
        let pipeline = Pipeline::parse(text)?;
        let body = c.require("loci")?;
        if body.len() > MAX_LOCI_BYTES {
            return Err(format!(
                "A job may hold at most {MAX_LOCI_BYTES} bytes of files; this one has {}.",
                body.len()
            ));
        }
        let files = unpack(c.require("names")?, c.require("sizes")?, body, MAX_LOCI)?;
        if files.is_empty() {
            return Err("The message holds no files to run the pipeline on.".into());
        }
        let mut seen = HashSet::new();
        let loci = files
            .into_iter()
            .map(|(name, text)| {
                if name.is_empty() || !seen.insert(name) {
                    return Err(format!(
                        "The files must have names that differ; \"{name}\" is empty or comes twice."
                    ));
                }
                let value = Value::from_text(pipeline.input, text)
                    .map_err(|why| format!("{name}: {why}"))?;
                rows_within_limit(name, &value)?;
                Ok((name.to_string(), value))
            })
            .collect::<Result<Vec<_>, String>>()?;
        Ok(FlowJob { pipeline, loci })
    }
}

/// What a run gave: the report of the engine, whole, and the pipeline it is a report of (the steps in it are named by
/// their places in the pipeline's list).
#[derive(Clone, Debug, PartialEq)]
pub struct FlowResult {
    pub pipeline: Pipeline,
    pub report: Report,
    pub seconds: f64,
}

/// Run the job here, keeping the answers of its steps in `cache` when there is a folder for them. A folder that cannot be
/// used does not stop the run: every step is then done, and the report says so.
pub fn run(job: FlowJob, cache: Option<&Path>) -> Result<FlowResult, String> {
    let started = Instant::now();
    let (cache, trouble) = match cache.map(Cache::open).transpose() {
        Ok(cache) => (cache, None),
        Err(why) => (
            None,
            Some(format!(
                "The workstation's cache could not be used, so every step was run: {why}"
            )),
        ),
    };
    let FlowJob { pipeline, loci } = job;
    let available = engine::available_for(&pipeline);
    let mut report = engine::run(&pipeline, loci, cache.as_ref(), &available)?;
    report.warnings.extend(trouble);
    Ok(FlowResult {
        pipeline,
        report,
        seconds: started.elapsed().as_secs_f64(),
    })
}

// --- the result on the wire ---------------------------------------------------------------------------------------------------

/// What each file of the run has when its own steps are done: the kind of the last of them.
fn each_gives(p: &Pipeline) -> Kind {
    match p.per_file() {
        0 => p.input,
        n => p.steps[n - 1].gives(),
    }
}

fn step_row(index: usize, r: &StepReport) -> String {
    let (word, why) = match &r.status {
        Status::Ran => ("ran", String::new()),
        Status::Cached => ("cached", String::new()),
        Status::Failed(why) => ("failed", clean_line(why)),
    };
    format!(
        "step\t{index}\t{word}\t{}\t{}\t{why}\n",
        r.seconds,
        clean_line(&r.note)
    )
}

impl FlowResult {
    pub fn encode(&self) -> Container {
        let p = &self.pipeline;
        let mut rows = String::new();
        let mut outputs: Vec<(&str, String)> = Vec::new();
        for item in &self.report.items {
            let name = clean_line(&item.name);
            rows.push_str(&format!(
                "item\t{name}\t{}\n",
                u8::from(item.output.is_some())
            ));
            // a file's reports are those of its first steps, in order (a step that fails ends them)
            for (k, r) in item.steps.iter().enumerate() {
                rows.push_str(&step_row(k, r));
            }
            if let Some(value) = &item.output {
                outputs.push((item.name.as_str(), value.to_text()));
            }
        }
        rows.push_str("gathered\n");
        // the steps done once follow the steps done for each file
        for (k, r) in self.report.gathered.iter().enumerate() {
            rows.push_str(&step_row(p.per_file() + k, r));
        }
        for name in &self.report.left_out {
            rows.push_str(&format!("left\t{}\n", clean_line(name)));
        }
        for text in &self.report.warnings {
            rows.push_str(&format!("warn\t{}\n", clean_line(text)));
        }
        let (names, sizes, body) = pack(&outputs);
        let mut c = Container::new()
            .with("meta", format!("seconds={}\n", self.seconds))
            .with("pipeline", p.to_text())
            .with("report", rows)
            .with("out_names", names)
            .with("out_sizes", sizes)
            .with("outputs", body);
        if let Some(value) = &self.report.output {
            c = c.with("output", value.to_text());
        }
        c
    }

    /// A reply comes from the network: every row is checked against the pipeline it says it is for.
    pub fn decode(c: &Container) -> Result<FlowResult, String> {
        let bad = |why: &str| format!("The pipeline result in the message is damaged: {why}.");
        let m = meta(c.require("meta")?);
        let seconds: f64 = number(&m, "seconds")?;
        if !(seconds.is_finite() && seconds >= 0.0) {
            return Err(bad("the time is not a number of seconds"));
        }
        let text = c.require("pipeline")?;
        if text.len() > MAX_PIPELINE_BYTES {
            return Err(bad("the pipeline is far longer than a pipeline is"));
        }
        let pipeline = Pipeline::parse(text)?;
        let (split, steps) = (pipeline.per_file(), pipeline.steps.len());
        let last = pipeline.check()?;
        let outputs = unpack(
            c.require("out_names")?,
            c.require("out_sizes")?,
            c.require("outputs")?,
            MAX_LOCI,
        )?;
        let mut outputs = outputs.into_iter();

        let mut report = Report::default();
        // 0: the files, 1: the steps done once, 2: the left out and the warnings
        let mut stage = 0u8;
        let mut open: Option<ItemReport> = None;
        let close = |open: &mut Option<ItemReport>, report: &mut Report| {
            if let Some(item) = open.take() {
                report.items.push(item);
            }
        };
        for row in c.require("report")?.lines() {
            let fields: Vec<&str> = row.split('\t').collect();
            match fields.as_slice() {
                ["item", name, has_output] if stage == 0 => {
                    close(&mut open, &mut report);
                    if report.items.len() >= MAX_LOCI {
                        return Err(bad("there are too many files"));
                    }
                    let output = match *has_output {
                        "0" => None,
                        "1" => {
                            let (out_name, text) =
                                outputs.next().ok_or_else(|| bad("an answer is missing"))?;
                            if out_name != *name {
                                return Err(bad("the answers are not in the order of the files"));
                            }
                            Some(
                                Value::from_text(each_gives(&pipeline), text)
                                    .map_err(|why| bad(&format!("{name}: {why}")))?,
                            )
                        }
                        _ => return Err(bad("a file does not say whether it has an answer")),
                    };
                    open = Some(ItemReport {
                        name: (*name).to_string(),
                        steps: Vec::new(),
                        output,
                    });
                }
                ["step", index, word, seconds, note, why] if stage <= 1 => {
                    let index: usize = index.parse().map_err(|_| bad("a step has no place"))?;
                    // the steps done for each file are the first ones, the others are done once
                    let allowed = if stage == 0 { 0..split } else { split..steps };
                    if !allowed.contains(&index) {
                        return Err(bad("a step is not one of those done there"));
                    }
                    let seconds: f64 = seconds.parse().map_err(|_| bad("a step has no time"))?;
                    if !(seconds.is_finite() && seconds >= 0.0) {
                        return Err(bad("a step has no time"));
                    }
                    let status = match *word {
                        "ran" => Status::Ran,
                        "cached" => Status::Cached,
                        "failed" => Status::Failed((*why).to_string()),
                        _ => return Err(bad("a step is neither run, found nor failed")),
                    };
                    let one = StepReport {
                        step: pipeline.steps[index].clone(),
                        status,
                        seconds,
                        note: (*note).to_string(),
                    };
                    match (stage, open.as_mut()) {
                        (0, Some(item)) if item.steps.len() < split => item.steps.push(one),
                        (0, _) => return Err(bad("a step belongs to no file")),
                        (_, _) if report.gathered.len() < steps => report.gathered.push(one),
                        _ => return Err(bad("too many steps")),
                    }
                }
                ["gathered"] if stage == 0 => {
                    close(&mut open, &mut report);
                    stage = 1;
                }
                ["left", name] if stage >= 1 && report.left_out.len() < MAX_LOCI => {
                    stage = 2;
                    report.left_out.push((*name).to_string());
                }
                ["warn", text] if stage >= 1 && report.warnings.len() < MAX_LOCI + 16 => {
                    stage = 2;
                    report.warnings.push((*text).to_string());
                }
                _ => return Err(bad("a row is not one that a report has there")),
            }
        }
        if stage == 0 {
            return Err(bad("the report stops before its steps done once"));
        }
        if outputs.next().is_some() {
            return Err(bad("there are answers that no file has"));
        }
        report.output = c
            .bytes("output")
            .map(|_| Value::from_text(last, c.require("output")?))
            .transpose()
            .map_err(|why| bad(&why))?;
        Ok(FlowResult {
            pipeline,
            report,
            seconds,
        })
    }
}

#[cfg(test)]
mod tests;
