//! Running a job on this machine: the same functions the app calls for itself, whoever asked. The server starts one
//! process for each job (`helixloom run-job`), which reads the request from the job's folder, runs it here and writes the
//! response beside it.

use super::species;
use super::tasks::{CnnResult, Folded, Request, Response};
use crate::pipeline::steps::{Aligned, Named};
use crate::pipeline::tools::{self, Installed};
use crate::pipeline::{self, predict, steps, Available};
use crate::seq::estimate::{estimate, EstimateOptions};
use crate::seq::msa::Msa;
use crate::seq::newick::NwTree;
use crate::seq::phylocnn::{estimate_many, Engine, Priors, Task};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

/// What the server decides for the jobs it runs.
#[derive(Clone, Debug, Default)]
pub struct Env {
    /// The program that folds a protein on a GPU: its words; the input and output files are added after them.
    pub fold_command: Option<Vec<String>>,
    /// How long that program may take, in seconds.
    pub fold_seconds: u64,
}

/// How often a running folding program is looked at.
const FOLD_POLL: Duration = Duration::from_millis(200);

/// Estimate every tree in `text` with the PhyloCNN networks. One tree always runs on the CPU: sending a single tree to
/// the GPU costs more than it saves.
pub fn run_phylocnn(
    text: &str,
    sampling: f64,
    task: Task,
    priors: Priors,
    gpu: bool,
) -> Result<CnnResult, String> {
    let trees = NwTree::parse_all(text)?;
    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    let (engine, name) = if gpu && trees.len() > 1 {
        (Engine::Gpu, "GPU")
    } else {
        (Engine::Cpu, "CPU")
    };
    #[cfg(not(all(feature = "gpu", not(target_arch = "wasm32"))))]
    let (engine, name) = {
        let _ = gpu;
        (Engine::Cpu, "CPU")
    };
    let batch = estimate_many(&trees, sampling, task, priors, engine)?;
    Ok(CnnResult {
        task,
        per_tree: batch.per_tree,
        skipped: batch.skipped,
        seconds: batch.seconds,
        engine: name.to_string(),
        trees: trees.len(),
        tips: trees.first().map_or(0, NwTree::leaf_count),
    })
}

fn installed(use_installed: bool) -> Available {
    if !use_installed {
        return Available::default();
    }
    let find = |p: &tools::Program| -> Option<Installed> { tools::detect(p) };
    Available {
        mafft: find(&tools::MAFFT),
        iqtree: find(&tools::IQTREE),
        primer3: find(&tools::PRIMER3),
    }
}

/// Fold with the workstation's own program: `command input.fasta output.pdb`, which must write a PDB file and stop
/// within `limit`. The program is named by the person who runs the server, never by a request.
fn fold_local(command: &[String], sequence: &[u8], limit: Duration) -> Result<String, String> {
    let (program, args) = command
        .split_first()
        .ok_or("The folding command is empty.")?;
    let dir = tools::work_folder("fold")?;
    let (input, output) = (dir.join("input.fasta"), dir.join("output.pdb"));
    std::fs::write(
        &input,
        format!(">query\n{}\n", String::from_utf8_lossy(sequence)),
    )
    .map_err(|e| format!("Could not write the sequence for the folding program: {e}"))?;
    let mut child = Command::new(program)
        .args(args)
        .arg(&input)
        .arg(&output)
        .stdin(Stdio::null())
        .spawn()
        .map_err(|e| format!("Could not start {program}: {e}"))?;
    let start = Instant::now();
    let outcome = loop {
        match child.try_wait() {
            Ok(Some(status)) if status.success() => break Ok(()),
            Ok(Some(status)) => break Err(format!("{program} stopped with {status}.")),
            Ok(None) if start.elapsed() > limit => {
                let _ = child.kill();
                let _ = child.wait();
                break Err(format!(
                    "{program} did not finish within {} seconds and was stopped.",
                    limit.as_secs()
                ));
            }
            Ok(None) => std::thread::sleep(FOLD_POLL),
            Err(e) => break Err(format!("Could not wait for {program}: {e}")),
        }
    };
    let pdb = outcome.and_then(|()| {
        std::fs::read_to_string(&output)
            .map_err(|e| format!("{program} wrote no structure file: {e}"))
    });
    let _ = std::fs::remove_dir_all(&dir);
    let pdb = pdb?;
    if !pdb.lines().any(|l| l.starts_with("ATOM")) {
        return Err(format!("{program} wrote a file with no atoms in it."));
    }
    Ok(pdb)
}

/// Align `seqs` with MAFFT when asked for and found here, else with the built-in aligner.
fn align(seqs: &[Named], mafft: bool) -> Result<Aligned, String> {
    let dir = tools::work_folder("align")?;
    let found = if mafft {
        tools::detect(&tools::MAFFT)
    } else {
        None
    };
    let aligned = steps::align(seqs, found.as_ref(), &dir);
    let _ = std::fs::remove_dir_all(&dir);
    aligned
}

/// Fold `sequence` with the workstation's own program (`local`), else through the ESM Atlas web service.
fn fold(sequence: &[u8], local: bool, env: &Env) -> Result<Folded, String> {
    let start = Instant::now();
    let (pdb, engine) = if local {
        let command = env
            .fold_command
            .as_ref()
            .ok_or("This workstation has no folding program configured.")?;
        let limit = Duration::from_secs(env.fold_seconds.max(1));
        (
            fold_local(command, sequence, limit)?,
            "the workstation's folding program",
        )
    } else {
        let curl = tools::detect(&tools::CURL).ok_or(
            "curl was not found on the workstation, so it cannot ask the ESM Atlas service.",
        )?;
        (
            predict::predict(&curl.path, predict::ESM_ATLAS_URL, sequence)?,
            "the ESM Atlas web service",
        )
    };
    Ok(Folded {
        pdb,
        engine: engine.to_string(),
        seconds: start.elapsed().as_secs_f64(),
    })
}

/// Run `request` here and give the response.
pub fn execute(request: Request, env: &Env) -> Result<Response, String> {
    match request {
        Request::Align { seqs, mafft } => align(&seqs, mafft).map(Response::Align),
        Request::Tree {
            alignment,
            analysis,
            model,
            bootstrap,
            generations,
        } => {
            let msa = Msa {
                names: alignment.iter().map(|r| r.0.clone()).collect(),
                rows: alignment.into_iter().map(|r| r.1).collect(),
            };
            let options = EstimateOptions {
                analysis,
                model,
                bootstrap,
                generations,
            };
            estimate(&msa, &options).map(Response::Tree)
        }
        Request::Pipeline {
            seqs,
            target,
            settings,
            choice,
        } => {
            let available = installed(choice.use_installed);
            pipeline::run(&seqs, target, &settings, &choice, &available).map(Response::Pipeline)
        }
        Request::PhyloCnn {
            newick,
            sampling,
            task,
            priors,
            gpu,
        } => run_phylocnn(&newick, sampling, task, priors, gpu).map(Response::PhyloCnn),
        Request::Fold { sequence, local } => fold(&sequence, local, env).map(Response::Fold),
        Request::Species(job) => species::run(&job).map(Response::Species),
    }
}

#[cfg(test)]
mod tests;
