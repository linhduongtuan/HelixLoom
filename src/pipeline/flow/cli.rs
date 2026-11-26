//! `helixloom pipeline`: run a pipeline file on a folder with no window, and with no server.

use super::export::{source_of, write_folder, Programs};
use super::{available_for, read_folder, run, Cache, Pipeline, Report, Status, Value};
use crate::pipeline::asm;
use std::path::{Path, PathBuf};

pub const USAGE: &str = "\
helixloom pipeline: pipelines as files of steps (`docs/PIPELINES.md`)

  helixloom pipeline run PIPELINE INPUT [--cache DIR] [--out FILE] [--out-dir DIR]
      run the pipeline on the files of the folder INPUT (or on one file); with --cache, a step whose input and
      settings have not changed is not run again, so a run that stopped is resumed by running it again
  helixloom pipeline export PIPELINE INPUT [--out-dir DIR] [--mafft NAME] [--iqtree NAME] [--helixloom NAME]
      write the pipeline as a Snakefile and a main.nf (Nextflow) into DIR (default: workflow) that run it on the files of the
      folder INPUT using the command-line programs for its steps
  helixloom pipeline check PIPELINE
      read the file and say what it does
  helixloom pipeline example [unaligned|gfa|plasmid|plasmid-short|plasmid-long|nanopore-plasmid|minion-reads]
      print a pipeline file for aligned loci, unaligned loci, a folder of GFA graphs, or one of the assembly routes (`docs/ROUTES.md`)
  helixloom pipeline clear-cache DIR
      forget what a cache holds

--out FILE       write a single result (such as a species tree or one graph's summary) to FILE
--out-dir DIR    write each input's result (tree, alignment or graph summary) to DIR, named after the input
";

struct Run {
    pipeline: PathBuf,
    input: PathBuf,
    cache: Option<PathBuf>,
    out: Option<PathBuf>,
    out_dir: Option<PathBuf>,
}

fn parse_run(args: &[String]) -> Result<Run, String> {
    let mut positional = Vec::new();
    let (mut cache, mut out, mut out_dir) = (None, None, None);
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .map(PathBuf::from)
                .ok_or(format!("{name} needs a value"))
        };
        match arg.as_str() {
            "--cache" => cache = Some(value("--cache")?),
            "--out" => out = Some(value("--out")?),
            "--out-dir" => out_dir = Some(value("--out-dir")?),
            flag if flag.starts_with("--") => return Err(format!("There is no option {flag}.")),
            _ => positional.push(PathBuf::from(arg)),
        }
    }
    match positional.as_slice() {
        [pipeline, input] => Ok(Run {
            pipeline: pipeline.clone(),
            input: input.clone(),
            cache,
            out,
            out_dir,
        }),
        _ => Err("run needs a pipeline file and an input folder (or file).".into()),
    }
}

/// What `export` is asked: the pipeline, the folder it will read, where to write the workflows and what the programs are called.
struct Export {
    pipeline: PathBuf,
    input: PathBuf,
    out_dir: PathBuf,
    programs: Programs,
}

fn parse_export(args: &[String]) -> Result<Export, String> {
    let mut positional = Vec::new();
    let mut out_dir = PathBuf::from("workflow");
    let mut programs = Programs::default();
    let mut it = args.iter();
    while let Some(arg) = it.next() {
        let mut value = |name: &str| it.next().cloned().ok_or(format!("{name} needs a value"));
        match arg.as_str() {
            "--out-dir" => out_dir = PathBuf::from(value("--out-dir")?),
            "--mafft" => programs.mafft = value("--mafft")?,
            "--iqtree" => programs.iqtree = value("--iqtree")?,
            "--helixloom" => programs.helixloom = value("--helixloom")?,
            flag if flag.starts_with("--") => return Err(format!("There is no option {flag}.")),
            _ => positional.push(PathBuf::from(arg)),
        }
    }
    match positional.as_slice() {
        [pipeline, input] => Ok(Export {
            pipeline: pipeline.clone(),
            input: input.clone(),
            out_dir,
            programs,
        }),
        _ => Err("export needs a pipeline file and the folder of the files it will run on.".into()),
    }
}

fn load(path: &Path) -> Result<Pipeline, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("Could not read {}: {e}", path.display()))?;
    Pipeline::parse(&text).map_err(|why| format!("{}: {why}", path.display()))
}

fn describe(p: &Pipeline) -> String {
    let mut text = format!("{}\ninput: {}\n", p.name, p.input.word());
    for (k, step) in p.steps.iter().enumerate() {
        let settings: Vec<String> = step
            .params()
            .iter()
            .map(|(k, v)| format!("{k}={v}"))
            .collect();
        text.push_str(&format!(
            "  {}. {} {}  ({} to {})\n",
            k + 1,
            step.name(),
            settings.join(" "),
            step.takes().word(),
            step.gives().word()
        ));
    }
    if let Some(split) = p.split_text() {
        text.push_str(&split);
        text.push('\n');
    }
    text
}

/// The extension a value is written with.
fn extension(v: &Value) -> &'static str {
    match v {
        Value::Reads(_) => "fastq",
        Value::Sequences(_) | Value::Alignment(_) => "fasta",
        Value::Tree(_) | Value::Trees(_) => "tre",
        Value::Gfa(_) => "gfa",
        Value::GraphSummary(_) => "tsv",
    }
}

fn write_to(path: &Path, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("Could not write {}: {e}", path.display()))
}

fn outputs(r: &Run, report: &Report) -> Result<Vec<String>, String> {
    let mut said = Vec::new();
    if let (Some(path), Some(value)) = (&r.out, &report.output) {
        write_to(path, &value.to_text())?;
        said.push(format!("Wrote {} to {}.", value.describe(), path.display()));
    } else if let Some(path) = &r.out {
        said.push(format!(
            "Nothing to write to {}: the run gave no single answer.",
            path.display()
        ));
    }
    if let Some(dir) = &r.out_dir {
        std::fs::create_dir_all(dir)
            .map_err(|e| format!("Could not make {}: {e}", dir.display()))?;
        let mut n = 0;
        for item in &report.items {
            if let Some(value) = &item.output {
                write_to(
                    &dir.join(format!("{}.{}", item.name, extension(value))),
                    &value.to_text(),
                )?;
                n += 1;
            }
        }
        said.push(format!("Wrote {n} files to {}.", dir.display()));
    }
    Ok(said)
}

/// What the command says, and whether everything worked.
fn work(args: &[String]) -> Result<(String, bool), String> {
    let Some(command) = args.first() else {
        return Ok((USAGE.to_string(), true));
    };
    let rest = &args[1..];
    match command.as_str() {
        "check" => {
            let [path] = rest else {
                return Err("check needs a pipeline file.".into());
            };
            let p = load(Path::new(path))?;
            let ends = p.check()?;
            Ok((
                format!("{}ends with {}.\n", describe(&p), ends.word()),
                true,
            ))
        }
        "example" => {
            let p = match rest {
                [] => Pipeline::loci_to_species(),
                [word] if word == "unaligned" => Pipeline::unaligned_loci_to_species(),
                [word] if word == "gfa" => Pipeline::inspect_graphs(),
                [word] => asm::routes::named(word).ok_or_else(|| {
                    let routes: Vec<&str> = asm::routes::NAMES.iter().map(|n| n.0).collect();
                    format!(
                        "Use pipeline example, pipeline example unaligned, pipeline example gfa, or pipeline example with one of {}.",
                        routes.join(", ")
                    )
                })?,
                _ => return Err("pipeline example takes one name, or none.".into()),
            };
            Ok((p.to_text(), true))
        }
        "export" => {
            let e = parse_export(rest)?;
            let p = load(&e.pipeline)?;
            p.check()?;
            if !e.input.is_dir() {
                return Err(format!(
                    "{} is not a folder: a workflow reads the files of a folder.",
                    e.input.display()
                ));
            }
            let source = source_of(&e.input, p.input)?;
            let written = write_folder(&e.out_dir, &p, &source, &e.programs)?;
            let names: Vec<String> = written.iter().map(|w| w.display().to_string()).collect();
            Ok((
                format!(
                    "Wrote {} for the files *.{} of {}.\nRun it with `snakemake --cores 4` or `nextflow run main.nf` in {}.\n",
                    names.join(" and "),
                    source.ext,
                    source.folder,
                    e.out_dir.display()
                ),
                true,
            ))
        }
        "clear-cache" => {
            let [dir] = rest else {
                return Err("clear-cache needs the cache folder.".into());
            };
            let n = Cache::open(Path::new(dir))?.clear()?;
            Ok((format!("Forgot {n} answers.\n"), true))
        }
        "run" => {
            let r = parse_run(rest)?;
            let p = load(&r.pipeline)?;
            let items = if r.input.is_dir() {
                read_folder(&r.input, p.input)?
            } else {
                let text = std::fs::read_to_string(&r.input)
                    .map_err(|e| format!("Could not read {}: {e}", r.input.display()))?;
                let name = r
                    .input
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .unwrap_or("input")
                    .to_string();
                vec![(
                    name,
                    Value::from_text(p.input, &text)
                        .map_err(|w| format!("{}: {w}", r.input.display()))?,
                )]
            };
            let cache = r.cache.as_deref().map(Cache::open).transpose()?;
            let report = run(&p, items, cache.as_ref(), &available_for(&p))?;
            let mut text = report.lines().join("\n");
            text.push('\n');
            for line in outputs(&r, &report)? {
                text.push_str(&line);
                text.push('\n');
            }
            text.push_str(&report.summary());
            text.push('\n');
            let ok = !report
                .items
                .iter()
                .flat_map(|i| i.steps.iter())
                .chain(report.gathered.iter())
                .any(|s| matches!(s.status, Status::Failed(_)));
            Ok((text, ok))
        }
        other => Err(format!(
            "There is no pipeline command called {other}. `helixloom pipeline` lists them."
        )),
    }
}

/// Run the command; the exit code: 0 when every step worked, 1 when one failed or the command could not run.
pub fn main(args: &[String]) -> i32 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{USAGE}");
        return 0;
    }
    match work(args) {
        Ok((text, ok)) => {
            print!("{text}");
            i32::from(!ok)
        }
        Err(why) => {
            eprintln!("{why}");
            1
        }
    }
}

#[cfg(test)]
mod tests;
