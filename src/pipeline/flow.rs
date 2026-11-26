//! Pipelines as data: a pipeline is a short text file of typed steps, run on one input or on every file of a folder, and a step
//! whose input and settings have not changed is not run again. A run that stopped is therefore resumed by running it again: what
//! was done is found in the cache, and only what was not done is done.
//!
//! The file looks like this (`#` starts a comment):
//!
//! ```text
//! helixloom-pipeline 1
//! name: Gene trees, then a species tree
//! input: alignment
//! step tree engine=builtin bootstrap=0
//! step species
//! ```
//!
//! The steps up to the first one that takes many trees are done for each file of the folder (each locus); the trees that come
//! out are gathered, and the rest of the steps are done once on all of them.

use super::asm::AsmStep;
use super::steps::{self, fasta_text, read_fasta, Aligned, Named};
use super::tools::{self, Installed};
use super::Available;
use crate::seq::fastq::{parse_fastq, to_fastq, FastqRead};
use crate::seq::newick::NwTree;
use crate::seq::par::par_queue;
use crate::seq::{astral, gfa};
use std::borrow::Cow;
use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use web_time::Instant;

/// The version of the file format; a file made by a later one is refused.
pub const VERSION: u32 = 1;
/// Raised when a step would answer differently for the same input: it makes every cached answer of the earlier kind unused.
const ENGINE: u32 = 1;

// --- what flows ---------------------------------------------------------------------------------------------------------------

/// The kind of what flows from one step to the next.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// Short or long reads with their qualities, as FASTQ.
    Reads,
    Sequences,
    Alignment,
    Tree,
    Trees,
    Gfa,
    GraphSummary,
}

impl Kind {
    pub fn word(self) -> &'static str {
        match self {
            Kind::Reads => "reads",
            Kind::Sequences => "sequences",
            Kind::Alignment => "alignment",
            Kind::Tree => "tree",
            Kind::Trees => "trees",
            Kind::Gfa => "gfa",
            Kind::GraphSummary => "graph-summary",
        }
    }

    pub fn parse(word: &str) -> Option<Kind> {
        [
            Kind::Reads,
            Kind::Sequences,
            Kind::Alignment,
            Kind::Tree,
            Kind::Trees,
            Kind::Gfa,
            Kind::GraphSummary,
        ]
        .into_iter()
        .find(|k| k.word() == word)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Reads(Vec<FastqRead>),
    Sequences(Vec<Named>),
    /// Rows of one length.
    Alignment(Vec<Named>),
    /// Newick.
    Tree(String),
    Trees(Vec<String>),
    /// Original graph text. The inspection step validates it and reports failures per input.
    Gfa(String),
    GraphSummary(gfa::Summary),
}

impl Value {
    pub fn kind(&self) -> Kind {
        match self {
            Value::Reads(_) => Kind::Reads,
            Value::Sequences(_) => Kind::Sequences,
            Value::Alignment(_) => Kind::Alignment,
            Value::Tree(_) => Kind::Tree,
            Value::Trees(_) => Kind::Trees,
            Value::Gfa(_) => Kind::Gfa,
            Value::GraphSummary(_) => Kind::GraphSummary,
        }
    }

    /// FASTA for sequences, Newick for trees, GFA for graphs and tab-separated metrics for summaries.
    pub fn to_text(&self) -> String {
        self.as_text().into_owned()
    }

    /// Borrow raw GFA so hashing a graph does not need a second copy of it.
    pub fn as_text(&self) -> Cow<'_, str> {
        match self {
            Value::Reads(reads) => Cow::Owned(to_fastq(reads)),
            Value::Sequences(rows) | Value::Alignment(rows) => Cow::Owned(fasta_text(rows)),
            Value::Tree(newick) => Cow::Owned(format!("{}\n", newick.trim())),
            Value::Trees(all) => {
                Cow::Owned(all.iter().map(|t| format!("{}\n", t.trim())).collect())
            }
            Value::Gfa(text) => Cow::Borrowed(text),
            Value::GraphSummary(summary) => Cow::Owned(summary.to_tsv()),
        }
    }

    pub fn from_text(kind: Kind, text: &str) -> Result<Value, String> {
        // a file made by some editors starts with a byte order mark, which is not part of the first name
        let without_bom = text.trim_start_matches('\u{feff}');
        let trees = || -> Vec<String> {
            without_bom
                .lines()
                .map(str::trim)
                .filter(|l| !l.is_empty())
                .map(str::to_string)
                .collect()
        };
        Ok(match kind {
            Kind::Reads => Value::Reads(parse_fastq(without_bom)?),
            Kind::Sequences => Value::Sequences(read_fasta(without_bom)),
            Kind::Alignment => {
                let rows = read_fasta(without_bom);
                let width = rows.first().map_or(0, |r| r.1.len());
                if let Some(bad) = rows.iter().find(|r| r.1.len() != width) {
                    return Err(format!(
                        "{} has {} columns and the first row {width}: an alignment has rows of one length.",
                        bad.0,
                        bad.1.len()
                    ));
                }
                Value::Alignment(rows)
            }
            Kind::Tree => match trees().as_slice() {
                [one] => Value::Tree(one.clone()),
                other => return Err(format!("Expected one tree, found {}.", other.len())),
            },
            Kind::Trees => Value::Trees(trees()),
            Kind::Gfa => Value::Gfa(text.to_owned()),
            Kind::GraphSummary => Value::GraphSummary(gfa::Summary::from_tsv(text)?),
        })
    }

    /// A few words for a report: "6 sequences, 600 columns".
    pub fn describe(&self) -> String {
        match self {
            Value::Reads(reads) => format!(
                "{} reads, {} bases",
                reads.len(),
                reads.iter().map(|r| r.seq.len()).sum::<usize>()
            ),
            Value::Sequences(rows) => format!("{} sequences", rows.len()),
            Value::Alignment(rows) => format!(
                "{} sequences, {} columns",
                rows.len(),
                rows.first().map_or(0, |r| r.1.len())
            ),
            Value::Tree(_) => "a tree".to_string(),
            Value::Trees(all) => format!("{} trees", all.len()),
            Value::Gfa(_) => "an assembly graph".to_owned(),
            Value::GraphSummary(summary) => format!(
                "{} segments, {} links, {} components",
                summary.segments, summary.links, summary.components
            ),
        }
    }
}

// --- steps and pipelines ------------------------------------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Step {
    /// Sequences into an alignment: the built-in centre-star method, or the installed MAFFT.
    Align { mafft: bool },
    /// An alignment into a maximum-likelihood tree: the built-in search, or the installed IQ-TREE (with bootstrap replicates).
    Tree { iqtree: bool, bootstrap: usize },
    /// Gene trees into the species tree that agrees with most of their quartets (the ASTRAL method).
    Species,
    /// Validate one assembly graph and count its segments, links, components and known bases.
    GfaInspect,
    /// A step of the assembly routes (`pipeline::asm`, `docs/ROUTES.md`).
    Asm(AsmStep),
}

impl Step {
    pub fn name(&self) -> &'static str {
        match self {
            Step::Align { .. } => "align",
            Step::Tree { .. } => "tree",
            Step::Species => "species",
            Step::GfaInspect => "gfa-inspect",
            Step::Asm(a) => a.name(),
        }
    }

    pub fn takes(&self) -> Kind {
        match self {
            Step::Align { .. } => Kind::Sequences,
            Step::Tree { .. } => Kind::Alignment,
            Step::Species => Kind::Trees,
            Step::GfaInspect => Kind::Gfa,
            Step::Asm(a) => a.takes(),
        }
    }

    pub fn gives(&self) -> Kind {
        match self {
            Step::Align { .. } => Kind::Alignment,
            Step::Tree { .. } | Step::Species => Kind::Tree,
            Step::GfaInspect => Kind::GraphSummary,
            Step::Asm(a) => a.gives(),
        }
    }

    /// The settings, all of them, in the order they are written.
    pub fn params(&self) -> Vec<(&'static str, String)> {
        let engine = |external: bool, name: &str| {
            (
                "engine",
                if external { name } else { "builtin" }.to_string(),
            )
        };
        match self {
            Step::Align { mafft } => vec![engine(*mafft, "mafft")],
            Step::Tree { iqtree, bootstrap } => vec![
                engine(*iqtree, "iqtree"),
                ("bootstrap", bootstrap.to_string()),
            ],
            Step::Species | Step::GfaInspect => Vec::new(),
            Step::Asm(a) => a.params(),
        }
    }

    fn parse(name: &str, params: &[(String, String)]) -> Result<Step, String> {
        let mut step = match name {
            "align" => Step::Align { mafft: false },
            "tree" => Step::Tree {
                iqtree: false,
                bootstrap: 0,
            },
            "species" => Step::Species,
            "gfa-inspect" => Step::GfaInspect,
            other => {
                return match AsmStep::parse(other, params)? {
                    Some(step) => Ok(Step::Asm(step)),
                    None => Err(format!(
                        "there is no step called {other} (the steps are align, tree, species, gfa-inspect, trim, filter, assemble, hybrid, flye, concatemer, circularize, rotate and medaka)"
                    )),
                }
            }
        };
        for (key, value) in params {
            match (&mut step, key.as_str(), value.as_str()) {
                (Step::Align { mafft }, "engine", "builtin") => *mafft = false,
                (Step::Align { mafft }, "engine", "mafft") => *mafft = true,
                (Step::Tree { iqtree, .. }, "engine", "builtin") => *iqtree = false,
                (Step::Tree { iqtree, .. }, "engine", "iqtree") => *iqtree = true,
                (Step::Tree { bootstrap, .. }, "bootstrap", n) => {
                    *bootstrap = n.parse().map_err(|_| {
                        format!("bootstrap must be a whole number of replicates, not {n}")
                    })?;
                }
                (_, "engine", other) => {
                    return Err(format!("{name} has no engine called {other}"));
                }
                (_, other, _) => return Err(format!("{name} has no setting called {other}")),
            }
        }
        Ok(step)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Pipeline {
    pub name: String,
    pub input: Kind,
    pub steps: Vec<Step>,
}

impl Pipeline {
    /// A separate graph summary for every GFA file in the input folder.
    pub fn inspect_graphs() -> Pipeline {
        Pipeline {
            name: "Inspect assembly graphs".to_owned(),
            input: Kind::Gfa,
            steps: vec![Step::GfaInspect],
        }
    }

    /// The first target of the roadmap: a folder of loci, each an alignment, into gene trees and then the species tree.
    pub fn loci_to_species() -> Pipeline {
        Pipeline {
            name: "Gene trees, then a species tree".to_string(),
            input: Kind::Alignment,
            steps: vec![
                Step::Tree {
                    iqtree: false,
                    bootstrap: 0,
                },
                Step::Species,
            ],
        }
    }

    /// The same from loci that are not aligned yet.
    pub fn unaligned_loci_to_species() -> Pipeline {
        let mut p = Pipeline::loci_to_species();
        p.name = "Align, gene trees, then a species tree".to_string();
        p.input = Kind::Sequences;
        p.steps.insert(0, Step::Align { mafft: false });
        p
    }

    pub fn parse(text: &str) -> Result<Pipeline, String> {
        let mut lines = text
            .lines()
            .enumerate()
            .map(|(k, l)| (k + 1, l.trim()))
            .filter(|(_, l)| !l.is_empty() && !l.starts_with('#'));
        let (n, head) = lines.next().ok_or("The pipeline file is empty.")?;
        let version = head
            .strip_prefix("helixloom-pipeline")
            .map(str::trim)
            .and_then(|v| v.parse::<u32>().ok())
            .ok_or(format!(
                "line {n}: a pipeline file starts with \"helixloom-pipeline {VERSION}\""
            ))?;
        if version > VERSION {
            return Err(format!(
                "This pipeline file is version {version}; this program reads up to version {VERSION}. Update the program."
            ));
        }
        let (mut name, mut input, mut steps) = (String::new(), None, Vec::new());
        for (n, line) in lines {
            let at = |why: String| format!("line {n}: {why}");
            if let Some(rest) = line.strip_prefix("name:") {
                name = rest.trim().to_string();
            } else if let Some(rest) = line.strip_prefix("input:") {
                let word = rest.trim();
                input = Some(Kind::parse(word).ok_or(at(format!(
                    "unknown input kind {word} (expected reads, sequences, alignment, tree, trees or gfa)"
                )))?);
            } else if let Some(rest) = line.strip_prefix("step") {
                let mut words = rest.split_whitespace();
                let which = words.next().ok_or(at("a step needs a name".into()))?;
                let mut params = Vec::new();
                for word in words {
                    let (k, v) = word
                        .split_once('=')
                        .ok_or(at(format!("{word}: settings are written key=value")))?;
                    params.push((k.to_string(), v.to_string()));
                }
                steps.push(Step::parse(which, &params).map_err(at)?);
            } else {
                return Err(at(format!("{line}: expected name:, input: or step")));
            }
        }
        let p = Pipeline {
            name,
            input: input.ok_or(
                "The pipeline file says no input: (reads, sequences, alignment, tree, trees or gfa).",
            )?,
            steps,
        };
        p.check()?;
        Ok(p)
    }

    pub fn to_text(&self) -> String {
        let mut out = format!("helixloom-pipeline {VERSION}\n");
        let _ = writeln!(out, "name: {}", self.name);
        let _ = writeln!(out, "input: {}", self.input.word());
        for step in &self.steps {
            let _ = write!(out, "step {}", step.name());
            for (k, v) in step.params() {
                let _ = write!(out, " {k}={v}");
            }
            out.push('\n');
        }
        out
    }

    /// Whether each step takes what the one before gives, and what the pipeline ends with. Trees are gathered for a step that
    /// takes many: a tree made for each file becomes the list of them.
    pub fn check(&self) -> Result<Kind, String> {
        if self.steps.is_empty() {
            return Err("A pipeline needs at least one step.".into());
        }
        let mut now = self.input;
        for (k, step) in self.steps.iter().enumerate() {
            if let Step::Tree {
                iqtree: true,
                bootstrap: 1..=999,
            } = step
            {
                return Err(format!(
                    "Step {} (tree) asks for fewer than 1000 bootstrap replicates; IQ-TREE's ultrafast bootstrap makes at least 1000, so the number here would not be what is run.",
                    k + 1
                ));
            }
            if let Step::Tree {
                iqtree: false,
                bootstrap: 1..,
            } = step
            {
                return Err(format!(
                    "Step {} (tree) asks for bootstrap replicates, which only engine=iqtree makes: the built-in search does not, and would leave them out without saying so.",
                    k + 1
                ));
            }
            let gathers = now == Kind::Tree && step.takes() == Kind::Trees;
            if step.takes() != now && !gathers {
                return Err(format!(
                    "Step {} ({}) takes {}, but {} is what comes to it.",
                    k + 1,
                    step.name(),
                    step.takes().word(),
                    now.word()
                ));
            }
            now = step.gives();
        }
        Ok(now)
    }

    /// Which steps are done for each file, and what happens to the trees after them; `None` when nothing is gathered.
    pub fn split_text(&self) -> Option<String> {
        let each = self.per_file();
        if each >= self.steps.len() {
            return None;
        }
        if each == 0 {
            return Some("Each file is a gene tree; they are gathered for the steps.".to_string());
        }
        let which = if each == 1 {
            "Step 1 is".to_string()
        } else {
            format!("Steps 1 to {each} are")
        };
        Some(format!(
            "{which} done for each file; the trees are then gathered for the rest."
        ))
    }

    /// How many steps are done for each file: those before the first that takes many trees. When the input is trees already,
    /// every step is done for each file.
    pub fn per_file(&self) -> usize {
        if self.input == Kind::Trees {
            return self.steps.len();
        }
        self.steps
            .iter()
            .position(|s| s.takes() == Kind::Trees)
            .unwrap_or(self.steps.len())
    }
}

// --- the cache ----------------------------------------------------------------------------------------------------------------

/// Two 64-bit FNV-1a sums of `bytes`, in hexadecimal: enough to tell two inputs apart, with no other crate and the same
/// answer on every machine and in every release (unlike the standard library's hashers).
pub fn digest(bytes: &[u8]) -> String {
    let (mut a, mut b) = (0xcbf2_9ce4_8422_2325u64, 0x8422_2325_cbf2_9ce4u64);
    for &byte in bytes {
        a = (a ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
        b = (b ^ u64::from(byte))
            .wrapping_mul(0x0100_0000_01b3)
            .rotate_left(5);
    }
    format!("{a:016x}{b:016x}")
}

/// What an answer is filed under: the step, its settings, the version of the program that does it (when it is an installed
/// one) and the input.
pub fn key(step: &Step, tool_version: &str, input: &Value) -> String {
    let settings: String = step
        .params()
        .iter()
        .map(|(k, v)| format!("{k}={v};"))
        .collect();
    // the version of this program is in the key: a release that changes what a built-in step does does not read the answers of
    // an earlier one (`ENGINE` is for a change made between releases)
    let text = format!(
        "{ENGINE}|{}|{}|{settings}|{tool_version}|{}|{}",
        env!("CARGO_PKG_VERSION"),
        step.name(),
        input.kind().word(),
        digest(input.as_text().as_bytes())
    );
    format!("{}-{}", step.name(), digest(text.as_bytes()))
}

/// Whether a file starts the way this program starts the answers it keeps: so that forgetting the answers of a cache never
/// removes another program's `.value` file from a folder that was chosen by mistake.
fn is_ours(path: &Path) -> bool {
    use std::io::BufRead;
    std::fs::File::open(path)
        .ok()
        .and_then(|f| std::io::BufReader::new(f).lines().next()?.ok())
        .is_some_and(|first| first == CACHE_HEAD)
}

/// A folder of answers, one file for each, named by `key`.
pub struct Cache {
    dir: PathBuf,
}

const CACHE_HEAD: &str = "helixloom-cache 2";
/// What a step's note goes on with when its answer could not be kept in the cache (then a colon and why).
const NOT_KEPT: &str = "; not kept in the cache";

impl Cache {
    pub fn open(dir: &Path) -> Result<Cache, String> {
        std::fs::create_dir_all(dir)
            .map_err(|e| format!("Could not make the cache folder {}: {e}", dir.display()))?;
        Ok(Cache {
            dir: dir.to_path_buf(),
        })
    }

    fn file(&self, key: &str) -> PathBuf {
        self.dir.join(format!("{key}.value"))
    }

    /// The answer and the note that went with it; `None` when there is none, or when the file is not one this program wrote.
    pub fn get(&self, key: &str, kind: Kind) -> Option<(Value, String)> {
        let text = std::fs::read_to_string(self.file(key)).ok()?;
        let mut parts = text.splitn(5, '\n');
        let (head, kind_line, note_line, digest_line, body) = (
            parts.next()?,
            parts.next()?,
            parts.next()?,
            parts.next()?,
            parts.next()?,
        );
        if head != CACHE_HEAD || kind_line != format!("kind: {}", kind.word()) {
            return None;
        }
        // a file cut short, or changed, is not the answer that was kept
        if digest_line != format!("digest: {}", digest(body.as_bytes())) {
            return None;
        }
        let note = note_line.strip_prefix("note: ")?.to_string();
        Value::from_text(kind, body).ok().map(|v| (v, note))
    }

    /// File an answer: written beside, then moved into place, so that no one reads half of it.
    pub fn put(&self, key: &str, value: &Value, note: &str) -> Result<(), String> {
        use std::sync::atomic::{AtomicU64, Ordering};
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let note = note.replace('\n', " ");
        let body = value.to_text();
        let text = format!(
            "{CACHE_HEAD}\nkind: {}\nnote: {note}\ndigest: {}\n{body}",
            value.kind().word(),
            digest(body.as_bytes())
        );
        let temporary = self.dir.join(format!(
            ".{key}.{}-{}.tmp",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&temporary, text).map_err(|e| e.to_string())?;
        std::fs::rename(&temporary, self.file(key)).map_err(|e| {
            let _ = std::fs::remove_file(&temporary);
            e.to_string()
        })
    }

    /// How many answers a folder keeps, without making the folder if it is not there.
    pub fn kept_in(dir: &Path) -> usize {
        Cache {
            dir: dir.to_path_buf(),
        }
        .len()
    }

    /// How many answers are kept.
    pub fn len(&self) -> usize {
        self.answers().count()
    }

    fn answers(&self) -> impl Iterator<Item = PathBuf> {
        std::fs::read_dir(&self.dir)
            .into_iter()
            .flatten()
            .flatten()
            .map(|e| e.path())
            .filter(|p| p.extension().is_some_and(|x| x == "value") && is_ours(p))
    }

    /// Forget every answer; how many there were.
    pub fn clear(&self) -> Result<usize, String> {
        let mut n = 0;
        for path in self.answers() {
            std::fs::remove_file(&path)
                .map_err(|e| format!("Could not remove {}: {e}", path.display()))?;
            n += 1;
        }
        Ok(n)
    }
}

// --- running ------------------------------------------------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq)]
pub enum Status {
    Ran,
    /// The answer of an earlier run was found.
    Cached,
    Failed(String),
}

#[derive(Clone, Debug, PartialEq)]
pub struct StepReport {
    pub step: Step,
    pub status: Status,
    pub seconds: f64,
    /// What made the answer ("centre-star alignment (built in)", "IQ-TREE v3.1.4"), or why it was not kept.
    pub note: String,
}

impl StepReport {
    /// The program that made the answer of an align or a tree step, as its note names it (the note of a tree step goes on
    /// with the model, and either can end with why the answer was not kept). None for a step that failed, and for the
    /// species and graph inspection steps, whose notes describe their results.
    pub fn program(&self) -> Option<&str> {
        if matches!(self.status, Status::Failed(_))
            || matches!(self.step, Step::Species | Step::GfaInspect)
        {
            return None;
        }
        let named = self.note.split(NOT_KEPT).next()?;
        let named = named.split(", model ").next()?.trim();
        (!named.is_empty()).then_some(named)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ItemReport {
    pub name: String,
    pub steps: Vec<StepReport>,
    /// What the last step gave, when none failed.
    pub output: Option<Value>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Report {
    pub items: Vec<ItemReport>,
    /// The steps done once, on everything that was gathered.
    pub gathered: Vec<StepReport>,
    /// The names of the files whose chain failed, and so are not in what was gathered.
    pub left_out: Vec<String>,
    /// What the person should know about the answer: a file whose tree was too small to be gathered, and so on.
    pub warnings: Vec<String>,
    /// What the whole run gives: the last gathered step's answer, or the only file's.
    pub output: Option<Value>,
}

impl Report {
    fn reports(&self) -> impl Iterator<Item = &StepReport> {
        self.items
            .iter()
            .flat_map(|i| i.steps.iter())
            .chain(self.gathered.iter())
    }

    pub fn ran(&self) -> usize {
        self.reports().filter(|r| r.status == Status::Ran).count()
    }

    pub fn cached(&self) -> usize {
        self.reports()
            .filter(|r| r.status == Status::Cached)
            .count()
    }

    pub fn failed(&self) -> usize {
        self.reports()
            .filter(|r| matches!(r.status, Status::Failed(_)))
            .count()
    }

    /// One line: how many steps were done, how many were found, how many failed.
    pub fn summary(&self) -> String {
        let mut text = format!(
            "{} steps run, {} found in the cache",
            self.ran(),
            self.cached()
        );
        if self.failed() > 0 {
            let _ = write!(text, ", {} failed", self.failed());
        }
        if !self.left_out.is_empty() {
            let _ = write!(text, "; left out: {}", self.left_out.join(", "));
        }
        text.push('.');
        text
    }

    /// Every step, one to a line, for a log.
    pub fn lines(&self) -> Vec<String> {
        let one = |who: &str, r: &StepReport| {
            let what = match &r.status {
                Status::Ran => format!("ran in {:.2} s", r.seconds),
                Status::Cached => "found in the cache".to_string(),
                Status::Failed(why) => format!("FAILED: {why}"),
            };
            let note = if r.note.is_empty() || matches!(r.status, Status::Failed(_)) {
                String::new()
            } else {
                format!(" ({})", r.note)
            };
            format!("{who}{}: {what}{note}", r.step.name())
        };
        let mut lines = Vec::new();
        for item in &self.items {
            let who = format!("{}: ", item.name);
            lines.extend(item.steps.iter().map(|r| one(&who, r)));
        }
        lines.extend(self.gathered.iter().map(|r| one("all: ", r)));
        lines
    }
}

/// Do one step on `input`, or take its answer from the cache.
fn do_step(
    step: &Step,
    name: &str,
    input: &Value,
    cache: Option<&Cache>,
    available: &Available,
) -> (StepReport, Option<Value>) {
    if let Step::Asm(asm) = step {
        return do_asm_step(step, asm, name, input, cache);
    }
    let installed: Option<&Installed> = match step {
        Step::Align { mafft: true } => available.mafft.as_ref(),
        Step::Tree { iqtree: true, .. } => available.iqtree.as_ref(),
        _ => None,
    };
    let report = |status, seconds, note: String| StepReport {
        step: step.clone(),
        status,
        seconds,
        note,
    };
    let wants_program = matches!(
        step,
        Step::Align { mafft: true } | Step::Tree { iqtree: true, .. }
    );
    if wants_program && installed.is_none() {
        let which = if matches!(step, Step::Align { .. }) {
            "MAFFT"
        } else {
            "IQ-TREE"
        };
        return (
            report(
                Status::Failed(format!(
                    "{which} is not installed here; the step asks for it by name, so the built-in engine is not used instead"
                )),
                0.0,
                String::new(),
            ),
            None,
        );
    }
    let version = installed.map_or("", |i| i.version.as_str());
    let cache = cache.map(|c| (c, key(step, version, input)));
    if let Some((value, note)) = cache.as_ref().and_then(|(c, key)| c.get(key, step.gives())) {
        return (report(Status::Cached, 0.0, note), Some(value));
    }
    let start = Instant::now();
    match work(step, input, installed) {
        Ok((value, note)) => {
            let saved = cache.map(|(c, key)| c.put(&key, &value, &note));
            let note = match saved {
                Some(Err(why)) => format!("{note}{NOT_KEPT}: {why}"),
                _ => note,
            };
            (
                report(Status::Ran, start.elapsed().as_secs_f64(), note),
                Some(value),
            )
        }
        Err(why) => (
            report(
                Status::Failed(why),
                start.elapsed().as_secs_f64(),
                String::new(),
            ),
            None,
        ),
    }
}

/// A step of the assembly routes: like `do_step`, with the program it calls (if one) checked first, and what it reads besides its input
/// (a file of long reads, of genes) in the key of its answer.
fn do_asm_step(
    step: &Step,
    asm: &AsmStep,
    name: &str,
    input: &Value,
    cache: Option<&Cache>,
) -> (StepReport, Option<Value>) {
    let report = |status, seconds, note: String| StepReport {
        step: step.clone(),
        status,
        seconds,
        note,
    };
    if let Err(why) = asm.ready() {
        return (report(Status::Failed(why), 0.0, String::new()), None);
    }
    let cache = cache.map(|c| (c, key(step, &asm.context(name), input)));
    if let Some((value, note)) = cache.as_ref().and_then(|(c, key)| c.get(key, step.gives())) {
        return (report(Status::Cached, 0.0, note), Some(value));
    }
    let start = Instant::now();
    match asm.work(input, name) {
        Ok((value, note)) => {
            let saved = cache.map(|(c, key)| c.put(&key, &value, &note));
            let note = match saved {
                Some(Err(why)) => format!("{note}{NOT_KEPT}: {why}"),
                _ => note,
            };
            (
                report(Status::Ran, start.elapsed().as_secs_f64(), note),
                Some(value),
            )
        }
        Err(why) => (
            report(
                Status::Failed(why),
                start.elapsed().as_secs_f64(),
                String::new(),
            ),
            None,
        ),
    }
}

/// The step itself.
fn work(
    step: &Step,
    input: &Value,
    installed: Option<&Installed>,
) -> Result<(Value, String), String> {
    // a folder for the files of an installed program; the built-in engines write nothing
    let folder = || {
        if installed.is_some() {
            tools::WorkFolder::new("flow")
        } else {
            Ok(tools::WorkFolder::none())
        }
    };
    match (step, input) {
        (Step::GfaInspect, Value::Gfa(text)) => {
            let value = Value::GraphSummary(gfa::inspect(text.as_bytes())?);
            let note = value.describe();
            Ok((value, note))
        }
        (Step::Align { .. }, Value::Sequences(rows)) => {
            let done = steps::align(rows, installed, &folder()?)?;
            Ok((Value::Alignment(done.rows), done.engine))
        }
        (Step::Tree { bootstrap, .. }, Value::Alignment(rows)) => {
            let aligned = Aligned {
                rows: rows.clone(),
                engine: String::new(),
                seconds: 0.0,
            };
            let done = steps::tree(&aligned, installed, *bootstrap, &folder()?)?;
            let note = format!("{}, model {}", done.engine, done.model);
            Ok((Value::Tree(done.newick), note))
        }
        (Step::Species, Value::Trees(all)) => {
            let found = astral::infer(&all.join("\n"), &astral::Options::default())?;
            let mut note = format!(
                "{} gene trees, {} species, {:.0} % of the quartets agree",
                found.genes,
                found.taxa.len(),
                found.normalized * 100.0
            );
            if found.skipped > 0 {
                let _ = write!(note, "; {} too small gene trees left out", found.skipped);
            }
            Ok((Value::Tree(found.newick), note))
        }
        (step, input) => Err(format!(
            "{} cannot take {}",
            step.name(),
            input.kind().word()
        )),
    }
}

/// Do `steps` one after the other on `start`; the reports, and the answer when none failed.
fn chain(
    steps: &[Step],
    name: &str,
    start: &Value,
    cache: Option<&Cache>,
    available: &Available,
) -> (Vec<StepReport>, Option<Value>) {
    let mut reports = Vec::new();
    let mut now = Cow::Borrowed(start);
    for step in steps {
        let (report, value) = do_step(step, name, &now, cache, available);
        reports.push(report);
        match value {
            Some(v) => now = Cow::Owned(v),
            None => return (reports, None),
        }
    }
    (reports, Some(now.into_owned()))
}

/// The fewest leaves a gene tree may have for the species-tree step to use it (the method's own limit).
const MOST_LEAVES_NEEDED: usize = 4;

/// How many leaves a Newick tree has; none when it cannot be read.
fn leaf_count(newick: &str) -> usize {
    NwTree::parse_all(newick)
        .ok()
        .and_then(|trees| {
            trees
                .first()
                .map(|t| t.nodes.iter().filter(|n| n.children.is_empty()).count())
        })
        .unwrap_or(0)
}

/// `work`, with a fault in it (a panic) turned into a failed step of the first of `steps`, so that one file's fault does not end a
/// run of many.
fn guarded(
    steps: &[Step],
    first: Option<&Step>,
    work: impl FnOnce() -> (Vec<StepReport>, Option<Value>),
) -> (Vec<StepReport>, Option<Value>) {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(work)) {
        Ok(done) => done,
        Err(payload) => {
            let why = payload
                .downcast_ref::<&str>()
                .map(|s| (*s).to_string())
                .or_else(|| payload.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "no message".to_string());
            let step = steps.first().or(first).cloned().unwrap_or(Step::Species);
            (
                vec![StepReport {
                    step,
                    status: Status::Failed(format!("the program faulted: {why}")),
                    seconds: 0.0,
                    note: String::new(),
                }],
                None,
            )
        }
    }
}

/// Run `pipeline` on `items` (a name and an input each): the first steps for every item at once on all the cores, then the
/// rest once on the trees that came out. A file that fails is left out of what is gathered and named in the report; the others
/// go on.
pub fn run(
    pipeline: &Pipeline,
    items: Vec<(String, Value)>,
    cache: Option<&Cache>,
    available: &Available,
) -> Result<Report, String> {
    pipeline.check()?;
    if items.is_empty() {
        return Err("There is nothing to run on: no input.".into());
    }
    if let Some((name, v)) = items.iter().find(|(_, v)| v.kind() != pipeline.input) {
        return Err(format!(
            "{name} is {}, and the pipeline takes {}.",
            v.kind().word(),
            pipeline.input.word()
        ));
    }
    let split = pipeline.per_file();
    let each = par_queue(&items, |(name, value)| {
        let (steps, output) = guarded(&pipeline.steps[..split], pipeline.steps.first(), || {
            chain(&pipeline.steps[..split], name, value, cache, available)
        });
        ItemReport {
            name: name.clone(),
            steps,
            output,
        }
    });
    let mut report = Report {
        left_out: each
            .iter()
            .filter(|i| i.output.is_none())
            .map(|i| i.name.clone())
            .collect(),
        ..Report::default()
    };
    if split < pipeline.steps.len() {
        let mut trees: Vec<String> = Vec::new();
        for item in &each {
            let Some(Value::Tree(t)) = &item.output else {
                continue;
            };
            // a gene tree of three leaves says nothing about how species are related: it is left out, and said to be
            let leaves = leaf_count(t);
            if leaves < MOST_LEAVES_NEEDED {
                report.left_out.push(item.name.clone());
                report.warnings.push(format!(
                    "{}: its gene tree has {leaves} leaves, and a species tree needs at least {MOST_LEAVES_NEEDED} in each gene tree it uses, so it was left out",
                    item.name
                ));
            } else {
                trees.push(t.clone());
            }
        }
        if trees.is_empty() {
            report.gathered.push(StepReport {
                step: pipeline.steps[split].clone(),
                status: Status::Failed(
                    "no gene tree was made, so there is nothing to gather".into(),
                ),
                seconds: 0.0,
                note: String::new(),
            });
        } else {
            let gathered = Value::Trees(trees);
            let (reports, output) =
                guarded(&pipeline.steps[split..], pipeline.steps.get(split), || {
                    chain(&pipeline.steps[split..], "", &gathered, cache, available)
                });
            report.gathered = reports;
            report.output = output;
        }
    } else if let [only] = each.as_slice() {
        report.output = only.output.clone();
    }
    report.items = each;
    Ok(report)
}

/// The installed programs a pipeline asks for by name, and no other: looking for them runs them.
pub fn available_for(p: &Pipeline) -> Available {
    let wants = |f: fn(&Step) -> bool| p.steps.iter().any(f);
    Available {
        mafft: wants(|s| matches!(s, Step::Align { mafft: true }))
            .then(|| tools::detect(&tools::MAFFT))
            .flatten(),
        iqtree: wants(|s| matches!(s, Step::Tree { iqtree: true, .. }))
            .then(|| tools::detect(&tools::IQTREE))
            .flatten(),
        primer3: None,
    }
}

// --- reading a folder -----------------------------------------------------------------------------------------------------------

/// The endings of the names of files that hold `kind`.
pub fn extensions(kind: Kind) -> &'static [&'static str] {
    match kind {
        Kind::Reads => &["fastq", "fq"],
        Kind::Sequences | Kind::Alignment => &["fasta", "fa", "fas", "fna", "faa", "aln"],
        Kind::Tree | Kind::Trees => &["tre", "tree", "nwk", "newick", "nex"],
        Kind::Gfa => &["gfa"],
        Kind::GraphSummary => &["tsv"],
    }
}

fn short_name(p: &Path) -> String {
    p.file_name().map_or_else(
        || p.display().to_string(),
        |n| n.to_string_lossy().into_owned(),
    )
}

/// The files of `dir` that a pipeline taking `kind` reads (by the ending of the name), in no order.
fn accepted_files(dir: &Path, kind: Kind) -> std::io::Result<Vec<PathBuf>> {
    let extensions = extensions(kind);
    Ok(std::fs::read_dir(dir)?
        .flatten()
        .map(|e| e.path())
        .filter(|p| {
            p.is_file()
                && p.extension()
                    .and_then(|x| x.to_str())
                    .is_some_and(|x| extensions.contains(&x.to_ascii_lowercase().as_str()))
        })
        .collect())
}

/// How many files `read_folder` would read, and how many bytes they hold: a listing of the folder, not a reading of the files.
/// `None` when the folder cannot be listed.
pub fn folder_size(dir: &Path, kind: Kind) -> Option<(usize, u64)> {
    folder_stamp(dir, kind).map(|stamp| (stamp.files, stamp.bytes))
}

/// A folder as a listing shows it, without reading its files: what says whether it is still what a run was made from.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Stamp {
    pub files: usize,
    pub bytes: u64,
    /// When the newest of the files was last changed, in seconds since 1970; 0 where that is not known.
    pub newest: u64,
    /// A checksum of the listing: every file's name, size and time of change, in the order of the names. A file that is
    /// renamed, or put back from a copy that is older, changes it.
    pub listing: String,
}

pub fn folder_stamp(dir: &Path, kind: Kind) -> Option<Stamp> {
    let mut files = accepted_files(dir, kind).ok()?;
    files.sort();
    let (mut bytes, mut newest, mut listing) = (0, 0, String::new());
    for path in &files {
        let known = std::fs::metadata(path).ok();
        let size = known.as_ref().map_or(0, std::fs::Metadata::len);
        let changed = known
            .and_then(|m| m.modified().ok())
            .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
            .map_or(0, |d| d.as_secs());
        bytes += size;
        newest = newest.max(changed);
        let _ = writeln!(listing, "{}\t{size}\t{changed}", short_name(path));
    }
    Some(Stamp {
        files: files.len(),
        bytes,
        newest,
        listing: digest(listing.as_bytes()),
    })
}

/// The named inputs in a folder, sorted by filename. GFA files retain their original text for the inspection step.
pub fn read_folder(dir: &Path, kind: Kind) -> Result<Vec<(String, Value)>, String> {
    let mut files =
        accepted_files(dir, kind).map_err(|e| format!("Could not read {}: {e}", dir.display()))?;
    files.sort();
    if files.is_empty() {
        return Err(format!(
            "{} holds no file of {} (the names end in {}).",
            dir.display(),
            kind.word(),
            extensions(kind).join(", ")
        ));
    }
    let mut seen: std::collections::HashMap<String, &PathBuf> = std::collections::HashMap::new();
    for path in &files {
        let stem = path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_default();
        if let Some(first) = seen.insert(stem.clone(), path) {
            return Err(format!(
                "{} and {} would both be called {stem}: give them names that differ before the dot.",
                short_name(first),
                short_name(path)
            ));
        }
    }
    files
        .iter()
        .map(|path| {
            let name = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("input")
                .to_string();
            let text = std::fs::read_to_string(path)
                .map_err(|e| format!("Could not read {}: {e}", path.display()))?;
            let value = if kind == Kind::Gfa {
                Value::Gfa(text)
            } else {
                Value::from_text(kind, &text).map_err(|why| format!("{}: {why}", path.display()))?
            };
            Ok((name, value))
        })
        .collect()
}

pub mod cli;
pub mod export;

#[cfg(test)]
pub(crate) mod tests;
