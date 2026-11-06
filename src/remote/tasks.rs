//! The jobs the workstation runs and how each travels. The app builds a `Request` and the server decodes it; the server
//! builds a `Response` and the app decodes it. Both are containers of sections (`wire.rs`): sequences as FASTA, trees as
//! Newick, settings and numbers as `key=value` lines, tables as tab-separated lines. Numbers are written the way Rust
//! prints them, which reads back as exactly the same number.

use super::wire::{meta, Container};
use super::{clean_line, one_line};
use crate::pipeline::steps::{fasta_text, read_fasta, Aligned, Named, PrimerPair, Primers, Tree};
use crate::pipeline::{Choice, Outcome};
use crate::seq::estimate::{Analysis, Estimate};
use crate::seq::newick::NwTree;
use crate::seq::phylocnn::{output_names, Estimate as CnnEstimate, Priors, Task};
use crate::seq::primer3::Settings;
use crate::seq::tree::DistanceModel;
use std::collections::HashMap;
use std::str::FromStr;

/// The most bootstrap resamples and Bayesian generations a request may ask for.
pub const MAX_BOOTSTRAP: usize = 10_000;
pub const MAX_GENERATIONS: usize = 10_000_000;
/// What a request may hold. Each is far more than a real study needs, and a job's memory grows with them (a tree of N
/// sequences builds a table of N by N numbers, the PhyloCNN networks keep a large input for every tree, the built-in
/// aligner indexes every letter), so together they keep what one request can ask of the workstation to a few gigabytes.
pub const MAX_TREES: usize = 20_000;
pub const MAX_SEQUENCES: usize = 5_000;
/// Letters of sequence in all, in a request.
pub const MAX_RESIDUES: usize = 64 * 1024 * 1024;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Align,
    Tree,
    Pipeline,
    PhyloCnn,
    Fold,
}

impl Kind {
    pub const ALL: [Kind; 5] = [
        Kind::Align,
        Kind::Tree,
        Kind::Pipeline,
        Kind::PhyloCnn,
        Kind::Fold,
    ];

    /// The name used in messages and in the list of jobs.
    pub fn name(self) -> &'static str {
        match self {
            Kind::Align => "align",
            Kind::Tree => "tree",
            Kind::Pipeline => "pipeline",
            Kind::PhyloCnn => "phylocnn",
            Kind::Fold => "fold",
        }
    }

    pub fn parse(name: &str) -> Option<Kind> {
        Kind::ALL.into_iter().find(|k| k.name() == name)
    }

    /// What the page calls it.
    pub fn label(self) -> &'static str {
        match self {
            Kind::Align => "Alignment",
            Kind::Tree => "Tree estimation",
            Kind::Pipeline => "Pipeline",
            Kind::PhyloCnn => "PhyloCNN estimates",
            Kind::Fold => "Structure prediction",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Request {
    Align {
        seqs: Vec<Named>,
        mafft: bool,
    },
    Tree {
        alignment: Vec<Named>,
        analysis: Analysis,
        model: DistanceModel,
        bootstrap: usize,
        generations: usize,
    },
    Pipeline {
        seqs: Vec<Named>,
        /// Where primers go on the first sequence: start and length.
        target: Option<(usize, usize)>,
        settings: Settings,
        choice: Choice,
    },
    PhyloCnn {
        newick: String,
        sampling: f64,
        task: Task,
        priors: Priors,
        gpu: bool,
    },
    Fold {
        sequence: Vec<u8>,
        /// Fold with the workstation's own program on a GPU, instead of asking the ESM Atlas service.
        local: bool,
    },
}

fn model_name(m: DistanceModel) -> &'static str {
    match m {
        DistanceModel::P => "p",
        DistanceModel::JukesCantor => "jc",
        DistanceModel::Kimura2 => "k2p",
    }
}

fn task_name(t: Task) -> &'static str {
    match t {
        Task::ChooseModel => "choose",
        Task::Bd => "bd",
        Task::Bdei => "bdei",
        Task::Bdss => "bdss",
    }
}

fn priors_name(p: Priors) -> &'static str {
    match p {
        Priors::Simulated => "simulated",
        Priors::Hiv => "hiv",
    }
}

macro_rules! settings_codec {
    ($($field:ident),* $(,)?) => {
        /// The primer settings as `key=value` lines.
        fn settings_text(s: &Settings) -> String {
            let mut out = String::new();
            $( out.push_str(&format!("{}={}\n", stringify!($field), s.$field)); )*
            out
        }

        /// The settings the lines give; one they leave out keeps primer3's default.
        fn settings_parse(text: &str) -> Result<Settings, String> {
            let m = meta(text);
            let mut s = Settings::default();
            $(
                if let Some(v) = m.get(stringify!($field)) {
                    s.$field = v.parse().map_err(|_| {
                        format!("The primer setting {} is not a number.", stringify!($field))
                    })?;
                }
            )*
            Ok(s)
        }
    };
}

settings_codec!(
    opt_size,
    min_size,
    max_size,
    opt_tm,
    min_tm,
    max_tm,
    min_gc,
    max_gc,
    salt_mm,
    divalent_mm,
    dntp_mm,
    dna_nm,
    max_poly_x,
    gc_clamp,
    max_end_gc,
    max_self_any,
    max_self_end,
    pair_compl_any,
    pair_compl_end,
    max_diff_tm,
    product_min,
    product_max,
    num_return,
);

fn lines(items: &[(&str, String)]) -> String {
    items.iter().map(|(k, v)| format!("{k}={v}\n")).collect()
}

/// A number from the message's `key=value` lines.
fn number<T: FromStr>(m: &HashMap<&str, &str>, key: &str) -> Result<T, String> {
    m.get(key)
        .and_then(|v| v.parse().ok())
        .ok_or_else(|| format!("The message's \"{key}\" is missing or is not a number."))
}

fn flag(m: &HashMap<&str, &str>, key: &str) -> bool {
    m.get(key) == Some(&"1")
}

fn fasta(seqs: &[Named]) -> String {
    let clean: Vec<Named> = seqs
        .iter()
        .map(|(name, bases)| (clean_line(name), bases.clone()))
        .collect();
    fasta_text(&clean)
}

fn sequences(c: &Container, section: &str) -> Result<Vec<Named>, String> {
    let seqs = read_fasta(c.require(section)?);
    if seqs.is_empty() || seqs.iter().any(|(_, b)| b.is_empty()) {
        return Err(format!(
            "The \"{section}\" section has no sequences, or an empty one."
        ));
    }
    Ok(seqs)
}

/// The sequences of a request, if there are not more of them, or more letters in them, than a job may be asked to hold.
fn within_limits(rows: Vec<Named>) -> Result<Vec<Named>, String> {
    if rows.len() > MAX_SEQUENCES {
        return Err(format!(
            "A job may hold at most {MAX_SEQUENCES} sequences; this one has {}.",
            rows.len()
        ));
    }
    let letters: usize = rows.iter().map(|r| r.1.len()).sum();
    if letters > MAX_RESIDUES {
        return Err(format!(
            "A job may hold at most {MAX_RESIDUES} letters of sequence; this one has {letters}."
        ));
    }
    Ok(rows)
}

fn bounded(what: &str, value: usize, most: usize) -> Result<usize, String> {
    if value > most {
        return Err(format!("{what} may be at most {most}."));
    }
    Ok(value)
}

fn bit(b: bool) -> String {
    u8::from(b).to_string()
}

fn encode_pipeline_request(
    seqs: &[Named],
    target: Option<(usize, usize)>,
    settings: &Settings,
    choice: &Choice,
) -> Container {
    Container::new()
        .with(
            "meta",
            lines(&[
                ("align", bit(choice.align)),
                ("tree", bit(choice.tree)),
                ("primers", bit(choice.primers)),
                ("use_installed", bit(choice.use_installed)),
                ("bootstrap", choice.bootstrap.to_string()),
                (
                    "target",
                    target.map_or(String::new(), |(a, n)| format!("{a},{n}")),
                ),
            ]),
        )
        .with("settings", settings_text(settings))
        .with("sequences", fasta(seqs))
}

fn encode_phylocnn_request(
    newick: &str,
    sampling: f64,
    task: Task,
    priors: Priors,
    gpu: bool,
) -> Container {
    Container::new()
        .with(
            "meta",
            lines(&[
                ("sampling", sampling.to_string()),
                ("task", task_name(task).into()),
                ("priors", priors_name(priors).into()),
                ("gpu", bit(gpu)),
            ]),
        )
        .with("trees", newick.to_string())
}

/// A tree request: the alignment (every row the same length) and how to build the tree from it.
fn decode_tree_request(c: &Container, m: &HashMap<&str, &str>) -> Result<Request, String> {
    let rows = within_limits(sequences(c, "alignment")?)?;
    if rows.iter().any(|r| r.1.len() != rows[0].1.len()) {
        return Err("The rows of an alignment must all be the same length.".into());
    }
    Ok(Request::Tree {
        alignment: rows,
        analysis: m
            .get("analysis")
            .and_then(|a| Analysis::from_name(a))
            .ok_or("The message's \"analysis\" is not a method this server knows.")?,
        model: match m.get("model").copied() {
            Some("p") => DistanceModel::P,
            Some("jc") => DistanceModel::JukesCantor,
            Some("k2p") => DistanceModel::Kimura2,
            _ => return Err("The message's \"model\" is not a distance model.".into()),
        },
        bootstrap: bounded("The bootstrap", number(m, "bootstrap")?, MAX_BOOTSTRAP)?,
        generations: bounded(
            "The number of generations",
            number(m, "generations")?,
            MAX_GENERATIONS,
        )?,
    })
}

/// A pipeline request: the sequences, where primers go, the primer settings and which steps to run.
fn decode_pipeline_request(c: &Container, m: &HashMap<&str, &str>) -> Result<Request, String> {
    let seqs = within_limits(sequences(c, "sequences")?)?;
    let target = match m.get("target").copied().unwrap_or("") {
        "" => None,
        t => {
            let (a, n) = t
                .split_once(',')
                .and_then(|(a, n)| Some((a.parse().ok()?, n.parse().ok()?)))
                .ok_or("The message's \"target\" is not \"start,length\".")?;
            Some((a, n))
        }
    };
    Ok(Request::Pipeline {
        seqs,
        target,
        settings: settings_parse(c.require("settings")?)?,
        choice: Choice {
            align: flag(m, "align"),
            tree: flag(m, "tree"),
            primers: flag(m, "primers"),
            use_installed: flag(m, "use_installed"),
            bootstrap: bounded("The bootstrap", number(m, "bootstrap")?, MAX_BOOTSTRAP)?,
        },
    })
}

/// A PhyloCNN request: the trees (not more than a job may hold), which question to ask of them and with which priors.
fn decode_phylocnn_request(c: &Container, m: &HashMap<&str, &str>) -> Result<Request, String> {
    let text = c.require("trees")?;
    // every tree ends with a semicolon: counting them costs no memory, parsing them all would
    let trees = text.matches(';').count();
    if trees > MAX_TREES {
        return Err(format!(
            "A job may hold at most {MAX_TREES} trees; this one has {trees}."
        ));
    }
    Ok(Request::PhyloCnn {
        newick: text.to_string(),
        sampling: number(m, "sampling")?,
        task: Task::ALL
            .into_iter()
            .find(|t| Some(&task_name(*t)) == m.get("task"))
            .ok_or("The message's \"task\" is not a PhyloCNN question.")?,
        priors: Priors::ALL
            .into_iter()
            .find(|p| Some(&priors_name(*p)) == m.get("priors"))
            .ok_or("The message's \"priors\" is not a PhyloCNN prior.")?,
        gpu: flag(m, "gpu"),
    })
}

impl Request {
    pub fn kind(&self) -> Kind {
        match self {
            Request::Align { .. } => Kind::Align,
            Request::Tree { .. } => Kind::Tree,
            Request::Pipeline { .. } => Kind::Pipeline,
            Request::PhyloCnn { .. } => Kind::PhyloCnn,
            Request::Fold { .. } => Kind::Fold,
        }
    }

    /// Whether the job needs a GPU of its own. A PhyloCNN job does only when it has many trees and this build can use a
    /// GPU: one tree always runs on the processor (`run_phylocnn`), and holding a GPU for it would keep real GPU jobs
    /// waiting.
    pub fn needs_gpu(&self) -> bool {
        match self {
            Request::PhyloCnn {
                gpu: true, newick, ..
            } => cfg!(feature = "gpu") && newick.matches(';').count() > 1,
            Request::Fold { local, .. } => *local,
            _ => false,
        }
    }

    pub fn encode(&self) -> Container {
        match self {
            Request::Align { seqs, mafft } => Container::new()
                .with(
                    "meta",
                    lines(&[("engine", (if *mafft { "mafft" } else { "builtin" }).into())]),
                )
                .with("sequences", fasta(seqs)),
            Request::Tree {
                alignment,
                analysis,
                model,
                bootstrap,
                generations,
            } => Container::new()
                .with(
                    "meta",
                    lines(&[
                        ("analysis", analysis.name().into()),
                        ("model", model_name(*model).into()),
                        ("bootstrap", bootstrap.to_string()),
                        ("generations", generations.to_string()),
                    ]),
                )
                .with("alignment", fasta(alignment)),
            Request::Pipeline {
                seqs,
                target,
                settings,
                choice,
            } => encode_pipeline_request(seqs, *target, settings, choice),
            Request::PhyloCnn {
                newick,
                sampling,
                task,
                priors,
                gpu,
            } => encode_phylocnn_request(newick, *sampling, *task, *priors, *gpu),
            Request::Fold { sequence, local } => Container::new()
                .with("meta", lines(&[("local", bit(*local))]))
                .with("sequence", sequence.clone()),
        }
    }

    pub fn decode(kind: Kind, c: &Container) -> Result<Request, String> {
        let text = c.require("meta")?;
        let m = meta(text);
        match kind {
            Kind::Align => Ok(Request::Align {
                seqs: within_limits(sequences(c, "sequences")?)?,
                mafft: m.get("engine") == Some(&"mafft"),
            }),
            Kind::Tree => decode_tree_request(c, &m),
            Kind::Pipeline => decode_pipeline_request(c, &m),
            Kind::PhyloCnn => decode_phylocnn_request(c, &m),
            Kind::Fold => Ok(Request::Fold {
                sequence: crate::pipeline::predict::clean_sequence(c.require("sequence")?)?,
                local: flag(&m, "local"),
            }),
        }
    }
}

/// What a run of the PhyloCNN networks gave: the estimates of every tree that could be used.
#[derive(Clone, Debug, PartialEq)]
pub struct CnnResult {
    pub task: Task,
    pub per_tree: Vec<Vec<CnnEstimate>>,
    pub skipped: Vec<(usize, String)>,
    pub seconds: f64,
    /// "CPU" or "GPU".
    pub engine: String,
    pub trees: usize,
    /// The tips of the first tree.
    pub tips: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Folded {
    pub pdb: String,
    pub engine: String,
    pub seconds: f64,
}

pub enum Response {
    Align(Aligned),
    Tree(Estimate),
    Pipeline(Outcome),
    PhyloCnn(CnnResult),
    Fold(Folded),
}

fn primer_rows(p: &Primers) -> String {
    p.pairs
        .iter()
        .map(|q| {
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                q.left,
                q.right,
                q.left_tm,
                q.right_tm,
                q.left_start,
                q.right_start,
                q.product_size,
                q.penalty
            )
        })
        .collect()
}

fn parse_primer_rows(text: &str) -> Result<Vec<PrimerPair>, String> {
    let bad = || "A primer pair in the message is damaged.".to_string();
    text.lines()
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let [left, right, ltm, rtm, ls, rs, size, penalty] = f.as_slice() else {
                return Err(bad());
            };
            Ok(PrimerPair {
                left: (*left).to_string(),
                right: (*right).to_string(),
                left_tm: ltm.parse().map_err(|_| bad())?,
                right_tm: rtm.parse().map_err(|_| bad())?,
                left_start: ls.parse().map_err(|_| bad())?,
                right_start: rs.parse().map_err(|_| bad())?,
                product_size: size.parse().map_err(|_| bad())?,
                penalty: penalty.parse().map_err(|_| bad())?,
            })
        })
        .collect()
}

fn newick_node(text: &str) -> Result<crate::seq::tree::Node, String> {
    NwTree::parse_all(text)?
        .first()
        .map(NwTree::to_node)
        .ok_or_else(|| "The message holds no tree.".into())
}

/// The lines that say which engine did a piece of work and how long it took.
fn timing(engine: &str, seconds: f64) -> String {
    lines(&[
        ("engine", one_line(engine)),
        ("seconds", seconds.to_string()),
    ])
}

/// The engine a reply's `key=value` lines name; none if they do not.
fn engine_of(m: &HashMap<&str, &str>) -> String {
    m.get("engine").unwrap_or(&"").to_string()
}

fn encode_pipeline_response(o: &Outcome) -> Container {
    let mut c = Container::new().with("log", o.log.join("\n"));
    if let Some(a) = &o.aligned {
        c = c
            .with("aligned_meta", timing(&a.engine, a.seconds))
            .with("aligned", fasta(&a.rows));
    }
    if let Some(t) = &o.tree {
        c = c
            .with(
                "tree_meta",
                lines(&[
                    ("engine", one_line(&t.engine)),
                    ("model", one_line(&t.model)),
                    (
                        "ln_likelihood",
                        t.ln_likelihood.map_or(String::new(), |l| l.to_string()),
                    ),
                    ("seconds", t.seconds.to_string()),
                ]),
            )
            .with("tree", t.newick.clone());
    }
    if let Some(p) = &o.primers {
        c = c
            .with("primers_meta", timing(&p.engine, p.seconds))
            .with("primers", primer_rows(p));
    }
    c
}

fn decode_pipeline_response(c: &Container) -> Result<Outcome, String> {
    let mut out = Outcome {
        log: c.require("log")?.lines().map(str::to_string).collect(),
        ..Outcome::default()
    };
    if c.bytes("aligned").is_some() {
        let m = meta(c.require("aligned_meta")?);
        out.aligned = Some(Aligned {
            rows: sequences(c, "aligned")?,
            engine: engine_of(&m),
            seconds: number(&m, "seconds")?,
        });
    }
    if c.bytes("tree").is_some() {
        let m = meta(c.require("tree_meta")?);
        out.tree = Some(Tree {
            newick: c.require("tree")?.to_string(),
            engine: engine_of(&m),
            model: m.get("model").unwrap_or(&"").to_string(),
            ln_likelihood: m.get("ln_likelihood").and_then(|l| l.parse().ok()),
            seconds: number(&m, "seconds")?,
        });
    }
    if c.bytes("primers").is_some() {
        let m = meta(c.require("primers_meta")?);
        out.primers = Some(Primers {
            pairs: parse_primer_rows(c.require("primers")?)?,
            engine: engine_of(&m),
            seconds: number(&m, "seconds")?,
        });
    }
    Ok(out)
}

fn encode_phylocnn_response(r: &CnnResult) -> Container {
    let rows: String = r
        .per_tree
        .iter()
        .enumerate()
        .map(|(i, estimates)| {
            let values: Vec<String> = estimates.iter().map(|e| e.value.to_string()).collect();
            format!("{i}\t{}\n", values.join(","))
        })
        .collect();
    let skipped: String = r
        .skipped
        .iter()
        .map(|(i, why)| format!("{i}\t{}\n", one_line(why)))
        .collect();
    Container::new()
        .with(
            "meta",
            lines(&[
                ("task", task_name(r.task).into()),
                ("seconds", r.seconds.to_string()),
                ("engine", one_line(&r.engine)),
                ("trees", r.trees.to_string()),
                ("tips", r.tips.to_string()),
            ]),
        )
        .with("estimates", rows)
        .with("skipped", skipped)
}

fn decode_phylocnn_response(c: &Container) -> Result<CnnResult, String> {
    let m = meta(c.require("meta")?);
    let task = Task::ALL
        .into_iter()
        .find(|t| Some(&task_name(*t)) == m.get("task"))
        .ok_or("The reply's \"task\" is not a PhyloCNN question.")?;
    let names = output_names(task);
    let bad = || "A PhyloCNN estimate in the reply is damaged.".to_string();
    let per_tree = c
        .require("estimates")?
        .lines()
        .map(|l| {
            let (_, values) = l.split_once('\t').ok_or_else(bad)?;
            let values: Vec<f64> = values
                .split(',')
                .map(|v| v.parse().map_err(|_| bad()))
                .collect::<Result<_, _>>()?;
            if values.len() != names.len() {
                return Err(bad());
            }
            Ok(names
                .iter()
                .zip(values)
                .map(|(&name, value)| CnnEstimate { name, value })
                .collect())
        })
        .collect::<Result<Vec<Vec<CnnEstimate>>, String>>()?;
    let skipped = c
        .require("skipped")?
        .lines()
        .filter_map(|l| {
            let (i, why) = l.split_once('\t')?;
            Some((i.parse().ok()?, why.to_string()))
        })
        .collect();
    Ok(CnnResult {
        task,
        per_tree,
        skipped,
        seconds: number(&m, "seconds")?,
        engine: engine_of(&m),
        trees: number(&m, "trees")?,
        tips: number(&m, "tips")?,
    })
}

/// `into_…` methods that take the answer to a question apart: the answer of the kind asked for, or why not.
macro_rules! into_variant {
    ($($method:ident, $variant:ident, $type:ty, $what:literal;)*) => {$(
        pub fn $method(self) -> Result<$type, String> {
            match self {
                Response::$variant(answer) => Ok(answer),
                _ => Err(concat!("The workstation answered with something other than ", $what, ".").into()),
            }
        }
    )*};
}

impl Response {
    into_variant! {
        into_align, Align, Aligned, "an alignment";
        into_tree, Tree, Estimate, "a tree";
        into_pipeline, Pipeline, Outcome, "a pipeline result";
        into_phylocnn, PhyloCnn, CnnResult, "PhyloCNN estimates";
        into_fold, Fold, Folded, "a structure";
    }

    pub fn encode(&self) -> Container {
        match self {
            Response::Align(a) => Container::new()
                .with("meta", timing(&a.engine, a.seconds))
                .with("rows", fasta(&a.rows)),
            // the exact writer: a tree made on the workstation is the very tree the app would have made
            Response::Tree(e) => Container::new()
                .with("newick", e.tree.newick_exact())
                .with("details", e.details.join("\n")),
            Response::Pipeline(o) => encode_pipeline_response(o),
            Response::PhyloCnn(r) => encode_phylocnn_response(r),
            Response::Fold(f) => Container::new()
                .with("meta", timing(&f.engine, f.seconds))
                .with("pdb", f.pdb.clone()),
        }
    }

    pub fn decode(kind: Kind, c: &Container) -> Result<Response, String> {
        match kind {
            Kind::Align => {
                let m = meta(c.require("meta")?);
                Ok(Response::Align(Aligned {
                    rows: sequences(c, "rows")?,
                    engine: engine_of(&m),
                    seconds: number(&m, "seconds")?,
                }))
            }
            Kind::Tree => Ok(Response::Tree(Estimate {
                tree: newick_node(c.require("newick")?)?,
                details: c.require("details")?.lines().map(str::to_string).collect(),
            })),
            Kind::Pipeline => decode_pipeline_response(c).map(Response::Pipeline),
            Kind::PhyloCnn => decode_phylocnn_response(c).map(Response::PhyloCnn),
            Kind::Fold => {
                let m = meta(c.require("meta")?);
                Ok(Response::Fold(Folded {
                    pdb: c.require("pdb")?.to_string(),
                    engine: engine_of(&m),
                    seconds: number(&m, "seconds")?,
                }))
            }
        }
    }
}

#[cfg(test)]
mod tests;
