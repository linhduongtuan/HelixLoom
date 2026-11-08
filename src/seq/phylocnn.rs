//! PhyloCNN (Perez & Gascuel 2025, github.com/manolofperez/phyloCNN): estimate the parameters of
//! birth-death epidemiological models (BD, BDEI, BDSS) and choose between them from a time-scaled
//! phylogeny, with a convolutional network. This is the paper's tree encoding (`Encoding/
//! PhyloCNN_Encoding_PhyloDyn.py`) and its published 2-generation networks, whose weights were
//! converted to a small binary file (kept in the models folder, `models.rs`: the authors' repository has no licence, so
//! the files are not part of this program) and are run here in plain Rust: grouped
//! convolution, batch normalisation, global average pooling and dense layers, exactly as Keras
//! computes them. The preprocessing (leaves and inner nodes as two sorted, zero-padded channels, the
//! sampling probability as a 19th column) follows the paper's training notebooks. Checked against
//! the output the authors' own notebook prints for their Zurich HIV tree.

use super::newick::NwTree;
use super::par::par_map;
use crate::models;

/// Trees larger than this do not fit the trained networks.
pub const MAX_TIPS: usize = 500;
/// Trees below this are smaller than anything the networks were trained on (the paper simulates 200 to 500).
pub const MIN_TRAINED_TIPS: usize = 200;
const FEATURES: usize = 18;
/// The mean branch length every tree is rescaled to (`TARGET_AVG_BL`).
const TARGET_MEAN_BRANCH: f64 = 1.0;

/// One encoded tree: a row of 18 features for every node, root first, level by level.
pub struct Encoding {
    pub rows: Vec<[f64; FEATURES]>,
    /// The factor the branch lengths were divided by.
    pub rescale: f64,
}

/// The paper's encoding of a rooted binary tree.
pub fn encode(tree: &NwTree) -> Result<Encoding, String> {
    let mut t = tree.clone();
    if t.nodes[t.root].children.len() == 1 {
        let child = t.nodes[t.root].children[0];
        t.nodes[child].parent = None;
        t.root = child;
    }
    // inner nodes with a branch of length zero become polytomies
    for i in t.postorder() {
        let (leaf, dist, parent) = (t.is_leaf(i), t.nodes[i].dist, t.nodes[i].parent);
        let Some(parent) = parent.filter(|_| !leaf && dist == 0.0) else {
            continue;
        };
        let kids = std::mem::take(&mut t.nodes[i].children);
        let at = t.nodes[parent].children.iter().position(|&c| c == i);
        if let Some(at) = at {
            t.nodes[parent]
                .children
                .splice(at..=at, kids.iter().copied());
        }
        kids.iter().for_each(|&k| t.nodes[k].parent = Some(parent));
        t.nodes[i].parent = None;
    }
    let order = t.levelorder();
    let mean = order.iter().map(|&i| t.nodes[i].dist).sum::<f64>() / order.len() as f64;
    let rescale = mean / TARGET_MEAN_BRANCH;
    if !(rescale > 0.0 && rescale.is_finite()) {
        return Err("The branch lengths of this tree average zero.".into());
    }
    for &i in &order {
        t.nodes[i].dist /= rescale;
    }
    if let Some(&bad) = order
        .iter()
        .find(|&&i| !t.is_leaf(i) && t.nodes[i].children.len() != 2)
    {
        return Err(format!(
            "PhyloCNN needs a binary tree, but a node has {} children.",
            t.nodes[bad].children.len()
        ));
    }
    let n = t.nodes.len();
    let mut leaves = vec![0.0f64; n];
    for i in t.postorder() {
        leaves[i] = if t.is_leaf(i) {
            1.0
        } else {
            t.nodes[i].children.iter().map(|&c| leaves[c]).sum()
        };
    }
    let mut from_root = vec![0.0f64; n];
    for &i in &order {
        from_root[i] = t.nodes[i]
            .parent
            .map_or(0.0, |p| from_root[p] + t.nodes[i].dist);
    }
    let deepest = order
        .iter()
        .filter(|&&i| t.is_leaf(i))
        .map(|&i| from_root[i])
        .fold(0.0, f64::max);
    let up = |i: usize| t.nodes[i].parent;
    // the distance down to a descendant's age, or -1 when there is none
    let mut rows = Vec::with_capacity(order.len());
    for &i in &order {
        let grand = up(i).and_then(up);
        let great = grand.and_then(up);
        let mut row = [-1.0f64; FEATURES];
        row[0] = deepest - from_root[i];
        row[1] = from_root[i];
        row[2] = if i == t.root { 0.0 } else { t.nodes[i].dist };
        row[3] = leaves[i];
        row[4] = grand.map_or(-1.0, |g| from_root[i] - from_root[g]);
        row[9] = great.map_or(-1.0, |g| from_root[i] - from_root[g]);
        if let [a, b] = t.nodes[i].children[..] {
            // the child with fewer leaves is the minor one; on a tie the second child is
            let (minor, major) = if leaves[a] >= leaves[b] {
                (b, a)
            } else {
                (a, b)
            };
            row[5] = t.nodes[minor].dist;
            row[6] = t.nodes[major].dist;
            row[7] = leaves[minor];
            row[8] = leaves[major];
            for (child, first) in [(minor, 10), (major, 12)] {
                if let [x, y] = t.nodes[child].children[..] {
                    let (small, large) = if leaves[x] >= leaves[y] {
                        (y, x)
                    } else {
                        (x, y)
                    };
                    row[first] = from_root[small] - from_root[i];
                    row[first + 1] = from_root[large] - from_root[i];
                    row[first + 4] = leaves[small];
                    row[first + 5] = leaves[large];
                }
            }
        }
        rows.push(row);
    }
    Ok(Encoding { rows, rescale })
}

/// The network's input, 500 rows x 19 features x 2 channels (leaves, then inner nodes), as the
/// training notebooks build it: each channel sorted by distance from the root, the inner nodes
/// followed by a copy of the last one, both padded with zeros; the sampling probability is the
/// 19th feature. Values are single precision, as the files the notebooks read.
pub fn model_input(enc: &Encoding, sampling: f64) -> Result<Vec<f32>, String> {
    let mut leaves: Vec<[f64; FEATURES + 1]> = Vec::new();
    let mut inner: Vec<[f64; FEATURES + 1]> = Vec::new();
    for row in &enc.rows {
        let mut full = [sampling; FEATURES + 1];
        full[..FEATURES].copy_from_slice(row);
        if row[3] == 1.0 {
            leaves.push(full);
        } else if row[3] > 1.0 {
            inner.push(full);
        }
    }
    if leaves.len() > MAX_TIPS {
        return Err(format!(
            "{} tips: the networks take at most {MAX_TIPS}.",
            leaves.len()
        ));
    }
    let by_age = |a: &[f64; FEATURES + 1], b: &[f64; FEATURES + 1]| a[1].total_cmp(&b[1]);
    leaves.sort_by(by_age);
    inner.sort_by(by_age);
    let Some(&last) = inner.last() else {
        return Err("The tree needs at least two tips.".into());
    };
    inner.push(last);
    let mut x = vec![0.0f32; MAX_TIPS * (FEATURES + 1) * 2];
    for (channel, rows) in [&leaves, &inner].into_iter().enumerate() {
        for (r, row) in rows.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                x[(r * (FEATURES + 1) + j) * 2 + channel] = v as f32;
            }
        }
    }
    Ok(x)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Act {
    Elu,
    Relu,
    Softmax,
    Linear,
}

enum Layer {
    Conv {
        groups: usize,
        kw: usize,
        cin: usize,
        filters: usize,
        act: Act,
        w: Vec<f32>,
    },
    Norm {
        eps: f32,
        gamma: Vec<f32>,
        beta: Vec<f32>,
        mean: Vec<f32>,
        var: Vec<f32>,
    },
    Pool,
    Dense {
        inputs: usize,
        outputs: usize,
        act: Act,
        w: Vec<f32>,
        b: Vec<f32>,
    },
}

/// A trained network: the layers of one of the paper's Keras models.
pub struct Network {
    input: [usize; 3],
    layers: Vec<Layer>,
}

enum State {
    /// rows x width x channels
    Grid(Vec<f32>),
    /// rows x channels
    Rows(Vec<f32>, usize),
    Flat(Vec<f32>),
}

fn activate(act: Act, x: &mut [f32]) {
    match act {
        Act::Elu => x
            .iter_mut()
            .for_each(|v| *v = if *v > 0.0 { *v } else { v.exp() - 1.0 }),
        Act::Relu => x.iter_mut().for_each(|v| *v = v.max(0.0)),
        Act::Linear => {}
        Act::Softmax => {
            let top = x.iter().copied().fold(f32::MIN, f32::max);
            x.iter_mut().for_each(|v| *v = (*v - top).exp());
            let sum: f32 = x.iter().sum();
            x.iter_mut().for_each(|v| *v /= sum);
        }
    }
}

impl Network {
    /// Read a network file written by the converter (header lines, a blank line, f32 weights).
    pub fn parse(bytes: &[u8]) -> Result<Network, String> {
        let cut = bytes
            .windows(2)
            .position(|w| w == b"\n\n")
            .ok_or("The network file has no header.")?;
        let header =
            std::str::from_utf8(&bytes[..cut]).map_err(|_| "The network header isn't text.")?;
        let (words, _) = bytes[cut + 2..].as_chunks::<4>();
        let mut weights = words.iter().map(|w| f32::from_le_bytes(*w));
        let mut take = |n: usize| -> Result<Vec<f32>, String> {
            let v: Vec<f32> = weights.by_ref().take(n).collect();
            (v.len() == n)
                .then_some(v)
                .ok_or_else(|| "The network file is cut short.".to_string())
        };
        let mut lines = header.lines();
        if lines.next() != Some("PCNN1") {
            return Err("This is not a PhyloCNN network file.".into());
        }
        let shape: Vec<usize> = lines
            .next()
            .and_then(|l| l.strip_prefix("input "))
            .map(|l| {
                l.split_whitespace()
                    .filter_map(|n| n.parse().ok())
                    .collect()
            })
            .unwrap_or_default();
        let input: [usize; 3] = shape
            .try_into()
            .map_err(|_| "The network input shape is not rows, width, channels.")?;
        let mut layers = Vec::new();
        for line in lines {
            let mut words = line.split_whitespace();
            if words.next() != Some("layer") {
                return Err(format!("Unexpected line in the network header: {line}"));
            }
            let kind = words.next().unwrap_or_default();
            let pairs: Vec<(&str, &str)> = words.filter_map(|w| w.split_once('=')).collect();
            let get = |key: &str| -> Result<&str, String> {
                pairs
                    .iter()
                    .find(|(k, _)| *k == key)
                    .map(|(_, v)| *v)
                    .ok_or_else(|| format!("Layer '{kind}' lacks {key}."))
            };
            let num = |key: &str| -> Result<usize, String> {
                get(key)?
                    .parse()
                    .map_err(|_| format!("Layer '{kind}': {key} is not a number."))
            };
            let act = |key: &str| -> Result<Act, String> {
                match get(key)? {
                    "elu" => Ok(Act::Elu),
                    "relu" => Ok(Act::Relu),
                    "softmax" => Ok(Act::Softmax),
                    "linear" => Ok(Act::Linear),
                    other => Err(format!("Unknown activation '{other}'.")),
                }
            };
            layers.push(match kind {
                "conv" => {
                    let (groups, kh, kw) = (num("groups")?, num("kh")?, num("kw")?);
                    let (cin, filters) = (num("cin")?, num("filters")?);
                    if kh != 1 || groups == 0 || cin % groups != 0 || filters % groups != 0 {
                        return Err(
                            "Only one-row convolutions with whole groups are supported.".into()
                        );
                    }
                    Layer::Conv {
                        groups,
                        kw,
                        cin,
                        filters,
                        act: act("act")?,
                        w: take(kw * (cin / groups) * filters)?,
                    }
                }
                "bn" => {
                    let n = num("n")?;
                    let eps = get("eps")?
                        .parse()
                        .map_err(|_| "Bad batch-normalisation epsilon.")?;
                    Layer::Norm {
                        eps,
                        gamma: take(n)?,
                        beta: take(n)?,
                        mean: take(n)?,
                        var: take(n)?,
                    }
                }
                "gap" => Layer::Pool,
                "dense" => {
                    let (inputs, outputs) = (num("in")?, num("out")?);
                    Layer::Dense {
                        inputs,
                        outputs,
                        act: act("act")?,
                        w: take(inputs * outputs)?,
                        b: take(outputs)?,
                    }
                }
                other => return Err(format!("Unknown layer '{other}'.")),
            });
        }
        Ok(Network { input, layers })
    }

    /// Run the network on one input of the shape it was trained on.
    pub fn run(&self, input: &[f32]) -> Result<Vec<f32>, String> {
        let [rows, width, channels] = self.input;
        if input.len() != rows * width * channels {
            return Err(format!(
                "The network takes {rows} x {width} x {channels} numbers."
            ));
        }
        let mut state = State::Grid(input.to_vec());
        for layer in &self.layers {
            state = match (layer, state) {
                (
                    Layer::Conv {
                        groups,
                        kw,
                        cin,
                        filters,
                        act,
                        w,
                    },
                    State::Grid(x),
                ) if *kw == width && *cin == channels => {
                    let per_group = cin / groups;
                    let mut out = vec![0.0f32; rows * filters];
                    for r in 0..rows {
                        for f in 0..*filters {
                            let base = f / (filters / groups) * per_group;
                            let mut sum = 0.0f32;
                            for j in 0..*kw {
                                for c in 0..per_group {
                                    sum += x[(r * width + j) * channels + base + c]
                                        * w[(j * per_group + c) * filters + f];
                                }
                            }
                            out[r * filters + f] = sum;
                        }
                    }
                    activate(*act, &mut out);
                    State::Rows(out, *filters)
                }
                (
                    Layer::Conv {
                        groups,
                        kw: 1,
                        cin,
                        filters,
                        act,
                        w,
                    },
                    State::Rows(x, have),
                ) if have == *cin => {
                    let per_group = cin / groups;
                    let mut out = vec![0.0f32; rows * filters];
                    for r in 0..rows {
                        for f in 0..*filters {
                            let base = f / (filters / groups) * per_group;
                            out[r * filters + f] = (0..per_group)
                                .map(|c| x[r * have + base + c] * w[c * filters + f])
                                .sum();
                        }
                    }
                    activate(*act, &mut out);
                    State::Rows(out, *filters)
                }
                (
                    Layer::Norm {
                        eps,
                        gamma,
                        beta,
                        mean,
                        var,
                    },
                    State::Rows(mut x, have),
                ) if gamma.len() == have => {
                    for (i, v) in x.iter_mut().enumerate() {
                        let c = i % have;
                        *v = gamma[c] * (*v - mean[c]) / (var[c] + eps).sqrt() + beta[c];
                    }
                    State::Rows(x, have)
                }
                (Layer::Pool, State::Rows(x, have)) => {
                    let mut mean = vec![0.0f32; have];
                    for (i, v) in x.iter().enumerate() {
                        mean[i % have] += v;
                    }
                    mean.iter_mut().for_each(|m| *m /= rows as f32);
                    State::Flat(mean)
                }
                (
                    Layer::Dense {
                        inputs,
                        outputs,
                        act,
                        w,
                        b,
                    },
                    State::Flat(x),
                ) if x.len() == *inputs => {
                    let mut out: Vec<f32> = (0..*outputs)
                        .map(|o| {
                            b[o] + (0..*inputs).map(|i| x[i] * w[i * outputs + o]).sum::<f32>()
                        })
                        .collect();
                    activate(*act, &mut out);
                    State::Flat(out)
                }
                _ => return Err("The network's layers do not fit together.".into()),
            };
        }
        match state {
            State::Flat(v) => Ok(v),
            _ => Err("The network does not end in a vector.".into()),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Task {
    /// Which of BD, BDEI and BDSS fits the tree.
    ChooseModel,
    Bd,
    Bdei,
    Bdss,
}

impl Task {
    pub const ALL: [Task; 4] = [Task::ChooseModel, Task::Bd, Task::Bdei, Task::Bdss];

    /// Whether a network exists for this question with these priors: the HIV priors were only
    /// trained for model choice and BDSS.
    pub fn trained_with(self, priors: Priors) -> bool {
        network_file(self, priors).is_some()
    }

    pub fn label(self) -> &'static str {
        match self {
            Task::ChooseModel => "Which model fits (BD, BDEI or BDSS)",
            Task::Bd => "BD: birth-death",
            Task::Bdei => "BDEI: with an exposed (incubation) stage",
            Task::Bdss => "BDSS: with superspreaders",
        }
    }
}

/// Which simulations the network was trained on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Priors {
    /// The paper's general epidemiological priors.
    Simulated,
    /// Priors chosen for the Zurich HIV data; only model choice and BDSS were trained.
    Hiv,
}

impl Priors {
    pub const ALL: [Priors; 2] = [Priors::Simulated, Priors::Hiv];

    pub fn label(self) -> &'static str {
        match self {
            Priors::Simulated => "general priors",
            Priors::Hiv => "HIV priors",
        }
    }

    /// What the training outbreaks looked like (README of the PhyloCNN repository).
    pub fn note(self) -> &'static str {
        match self {
            Priors::Simulated => {
                "Trained on outbreaks with R0 1–5, infectious period 1–10, 200–500 tips, sampling 0.01–1."
            }
            Priors::Hiv => "Trained on outbreaks simulated with priors chosen for HIV.",
        }
    }
}

/// The file of the trained network for a question with these priors; none for the questions the HIV priors were not trained
/// for.
fn network_file(task: Task, priors: Priors) -> Option<&'static str> {
    Some(match (priors, task) {
        (Priors::Simulated, Task::ChooseModel) => "Trained_2Generation_PhyDyn.bin",
        (Priors::Simulated, Task::Bd) => "Trained_2Generation_BD.bin",
        (Priors::Simulated, Task::Bdei) => "Trained_2Generation_BDEI.bin",
        (Priors::Simulated, Task::Bdss) => "Trained_2Generation_BDSS.bin",
        (Priors::Hiv, Task::ChooseModel) => "Trained_2Generation_PhyDyn_HIV.bin",
        (Priors::Hiv, Task::Bdss) => "Trained_2Generation_BDSS_HIV.bin",
        (Priors::Hiv, _) => return None,
    })
}

/// The trained network for a question, from the models folder: the files belong to the authors of PhyloCNN and are not part
/// of this program (see `models`).
fn network_bytes(task: Task, priors: Priors) -> Result<Vec<u8>, String> {
    let file = network_file(task, priors).ok_or(NO_HIV_NETWORK)?;
    models::read(&models::PHYLOCNN, file)
}

/// The Zurich HIV tree from the paper's repository (Rasmussen et al. 2017): 200 sequences. Also from the models folder.
pub fn zurich_tree() -> Result<String, String> {
    models::read_text(&models::PHYLOCNN, "ZurichHIV-tree.nwk")
}

/// The tree the PhyloCNN window offers as its example, and whether it is the Zurich HIV tree: where that file is not
/// installed it is a made-up tree of as many tips.
pub fn example_tree() -> (String, bool) {
    match zurich_tree() {
        Ok(text) => (text, true),
        Err(_) => (crate::sample::made_up_tree(EXAMPLE_TIPS), false),
    }
}

/// The number of tips of the Zurich HIV tree.
const EXAMPLE_TIPS: usize = 200;

#[derive(Clone, Debug, PartialEq)]
pub struct Estimate {
    pub name: &'static str,
    pub value: f64,
}

/// What the network's outputs are called for a question, in order.
pub fn output_names(task: Task) -> &'static [&'static str] {
    match task {
        Task::ChooseModel => &["BD", "BDEI", "BDSS"],
        Task::Bd => &["R0", "infectious period"],
        Task::Bdei => &["R0", "infectious period", "incubation period"],
        Task::Bdss => &[
            "R0",
            "infectious period",
            "X_SS (superspreader transmission ratio)",
            "f_SS (share of superspreaders)",
        ],
    }
}

/// Name the network's raw outputs and bring the times back to the tree's own units: the networks work in
/// rescaled time, so these outputs are multiplied by the rescale factor.
fn finish(task: Task, raw: &[f32], rescale: f64) -> Result<Vec<Estimate>, String> {
    let names = output_names(task);
    if raw.len() != names.len() {
        return Err("The network returned an unexpected number of values.".into());
    }
    let in_time = |i: usize| match task {
        Task::Bd | Task::Bdss => i == 1,
        Task::Bdei => i == 1 || i == 2,
        Task::ChooseModel => false,
    };
    Ok(names
        .iter()
        .zip(raw)
        .enumerate()
        .map(|(i, (&name, &v))| Estimate {
            name,
            value: f64::from(v) * if in_time(i) { rescale } else { 1.0 },
        })
        .collect())
}

fn check_sampling(sampling: f64) -> Result<(), String> {
    if sampling > 0.0 && sampling <= 1.0 {
        Ok(())
    } else {
        Err("The sampling probability must be above 0 and at most 1.".into())
    }
}

const NO_HIV_NETWORK: &str = "Only model choice and BDSS were trained with the HIV priors.";

/// Run one of the trained networks on a time-scaled tree. `sampling` is the share of infected
/// people the tree contains (0 to 1). Times come back in the units of the tree's branch lengths.
// ponytail: only the tests use the one-tree call, the window goes through `estimate_many` (a batch of one).
#[cfg(test)]
pub fn estimate(
    tree: &NwTree,
    sampling: f64,
    task: Task,
    priors: Priors,
) -> Result<Vec<Estimate>, String> {
    check_sampling(sampling)?;
    let bytes = network_bytes(task, priors)?;
    let encoding = encode(tree)?;
    let out = Network::parse(&bytes)?.run(&model_input(&encoding, sampling)?)?;
    finish(task, &out, encoding.rescale)
}

/// Where the networks run for many trees.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Engine {
    /// The plain Rust code on every CPU core.
    Cpu,
    /// candle on the GPU (not in the browser, where candle does not build).
    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    Gpu,
}

/// The answer for a set of trees: each tree's estimates, the trees that could not be used, and how long
/// the networks took (the encoding is not counted).
pub struct Batch {
    pub per_tree: Vec<Vec<Estimate>>,
    pub skipped: Vec<(usize, String)>,
    pub seconds: f64,
    /// How many GPUs worked on the trees (none when the processor did).
    pub gpus: usize,
}

/// Run a network on every tree. A tree that cannot be encoded (not binary, too many tips) is skipped and
/// reported, as long as at least one tree can be used.
pub fn estimate_many(
    trees: &[NwTree],
    sampling: f64,
    task: Task,
    priors: Priors,
    engine: Engine,
) -> Result<Batch, String> {
    check_sampling(sampling)?;
    let net = Network::parse(&network_bytes(task, priors)?)?;
    let encoded: Vec<Result<(f64, Vec<f32>), String>> = par_map(trees, |t| {
        let e = encode(t)?;
        Ok((e.rescale, model_input(&e, sampling)?))
    });
    let (mut rescales, mut inputs, mut skipped) = (Vec::new(), Vec::new(), Vec::new());
    for (i, item) in encoded.into_iter().enumerate() {
        match item {
            Ok((rescale, input)) => {
                rescales.push(rescale);
                inputs.push(input);
            }
            Err(why) => skipped.push((i, why)),
        }
    }
    if inputs.is_empty() {
        let why = skipped
            .first()
            .map_or("There is no tree.", |s| s.1.as_str());
        return Err(format!("No tree can be used. {why}"));
    }
    let start = web_time::Instant::now();
    let (raw, gpus): (Vec<Vec<f32>>, usize) = match engine {
        Engine::Cpu => (
            par_map(&inputs, |x| net.run(x))
                .into_iter()
                .collect::<Result<_, _>>()?,
            0,
        ),
        #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
        Engine::Gpu => gpu::run_batch(&net, &inputs)?,
    };
    let seconds = start.elapsed().as_secs_f64();
    let per_tree = raw
        .iter()
        .zip(&rescales)
        .map(|(r, &rescale)| finish(task, r, rescale))
        .collect::<Result<_, _>>()?;
    Ok(Batch {
        per_tree,
        skipped,
        seconds,
        gpus,
    })
}

/// One estimate across many trees: the mean and the middle 95 % (what a posterior sample is summarised by).
#[derive(Clone, Debug, PartialEq)]
pub struct Spread {
    pub name: &'static str,
    pub mean: f64,
    pub median: f64,
    pub low: f64,
    pub high: f64,
}

/// A quantile with the values sorted, between neighbours (`p` from 0 to 1).
fn quantile(sorted: &[f64], p: f64) -> f64 {
    let at = p * (sorted.len() - 1) as f64;
    let (lo, hi) = (at.floor() as usize, at.ceil() as usize);
    sorted[lo] + (sorted[hi] - sorted[lo]) * (at - lo as f64)
}

/// Mean, median and the 2.5 % and 97.5 % points of every estimate over the trees.
pub fn summarize(per_tree: &[Vec<Estimate>]) -> Vec<Spread> {
    let Some(first) = per_tree.first() else {
        return Vec::new();
    };
    (0..first.len())
        .map(|k| {
            let mut values: Vec<f64> = per_tree.iter().map(|t| t[k].value).collect();
            values.sort_by(f64::total_cmp);
            Spread {
                name: first[k].name,
                mean: values.iter().sum::<f64>() / values.len() as f64,
                median: quantile(&values, 0.5),
                low: quantile(&values, 0.025),
                high: quantile(&values, 0.975),
            }
        })
        .collect()
}

#[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
pub mod gpu;

#[cfg(test)]
mod tests;
