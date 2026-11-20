//! How the annotations of a tree show when it is drawn, after the choices FigTree gives (colour and width of branches by an annotation,
//! labels, shapes and bars at nodes, a legend), worked out as plain numbers and colours that the viewer and the SVG picture both use.
//!
//! The rules are this program's own reading of what FigTree's documentation says; no source was read. A number gets a colour from a
//! gradient between the least and the greatest value there is (a dark violet through green to yellow, a scale that stays readable in
//! grey and to the colour-blind), anything else a colour for each different value.

use super::newick::NwTree;
use super::treelayout::Laid;

pub type Rgb = [u8; 3];

/// Colours for the different values of a text annotation, and for groups.
pub const PALETTE: [Rgb; 8] = [
    [0x2f, 0x6f, 0xed],
    [0xd6, 0x45, 0x45],
    [0x1e, 0x9e, 0x5a],
    [0xc7, 0x7c, 0x0e],
    [0x8a, 0x4f, 0xd0],
    [0x0e, 0x9a, 0xa7],
    [0xc2, 0x3c, 0x8f],
    [0x6b, 0x7a, 0x2f],
];

/// The ends and the middle of the gradient for a number.
const GRADIENT: [Rgb; 3] = [[68, 1, 84], [33, 145, 140], [253, 231, 37]];

/// A branch is drawn this many times as wide as usual at the greatest value of the annotation that sets its width, and as it is at the least.
const MOST_WIDTH: f32 = 4.0;
/// The radius, in points, of a node's shape at the least and at the greatest value.
pub const SHAPE_RADII: (f32, f32) = (2.5, 7.0);
/// The radius of a node's shape when the annotation is not a number.
const SHAPE_PLAIN: f32 = 3.5;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StyleOptions {
    /// Colour the branches by this annotation.
    pub colour_by: Option<String>,
    /// Make the branches wider where this (numeric) annotation is greater.
    pub width_by: Option<String>,
    /// Write these annotations beside the inner nodes, one after the other with a slash between; an empty key is the node's own label.
    pub labels: Vec<String>,
    /// Put a shape on each node that has this annotation; its size follows the annotation when it is a number.
    pub shape_by: Option<String>,
    /// Draw a bar across each node for this annotation, which is a range of heights written `{least,greatest}` (a 95 % HPD interval).
    pub bars_by: Option<String>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Legend {
    Gradient {
        key: String,
        low: f64,
        high: f64,
    },
    Categories {
        key: String,
        items: Vec<(String, Rgb)>,
    },
}

/// What is drawn on each node of a layout, in the order of `Laid::nodes`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Decor {
    pub colour: Vec<Option<Rgb>>,
    /// How many times as wide as usual the branch to the node is drawn.
    pub width: Vec<f32>,
    pub label: Vec<String>,
    pub shape: Vec<Option<f32>>,
    /// The ends of the bar across the node, in the layout's distance from the root.
    pub bar: Vec<Option<(f32, f32)>>,
    pub legend: Option<Legend>,
}

fn mix(a: Rgb, b: Rgb, t: f64) -> Rgb {
    let t = t.clamp(0.0, 1.0);
    [0, 1, 2].map(|k| (f64::from(a[k]) + (f64::from(b[k]) - f64::from(a[k])) * t).round() as u8)
}

/// The colour of a value `t` of the way (0 to 1) along the gradient.
pub fn gradient(t: f64) -> Rgb {
    if t < 0.5 {
        mix(GRADIENT[0], GRADIENT[1], t * 2.0)
    } else {
        mix(GRADIENT[1], GRADIENT[2], (t - 0.5) * 2.0)
    }
}

/// What an annotation holds over a tree: numbers, or words.
enum Values {
    Numbers { low: f64, high: f64 },
    Words(Vec<String>),
    None,
}

/// The value of annotation `key` on node `i` of the tree as written; the node's own label when `key` is empty and the node is not a tip.
fn raw(t: &NwTree, i: usize, key: &str) -> Option<String> {
    if key.is_empty() {
        let n = &t.nodes[i];
        return (!n.children.is_empty() && !n.name.is_empty()).then(|| n.name.clone());
    }
    t.attr(i, key).map(str::to_string)
}

fn survey(t: &NwTree, key: &str) -> Values {
    let found: Vec<String> = (0..t.nodes.len()).filter_map(|i| raw(t, i, key)).collect();
    if found.is_empty() {
        return Values::None;
    }
    let numbers: Vec<f64> = found
        .iter()
        .filter_map(|v| v.parse::<f64>().ok())
        .filter(|v| v.is_finite())
        .collect();
    if numbers.len() == found.len() {
        let low = numbers.iter().copied().fold(f64::INFINITY, f64::min);
        let high = numbers.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        return Values::Numbers { low, high };
    }
    let mut words: Vec<String> = Vec::new();
    for v in found {
        if !words.contains(&v) {
            words.push(v);
        }
    }
    Values::Words(words)
}

/// How far (0 to 1) `v` is between the least and the greatest value; a tree whose values are all one number has every value in the middle.
fn along(v: f64, low: f64, high: f64) -> f64 {
    if high > low {
        (v - low) / (high - low)
    } else {
        0.5
    }
}

/// A number as a label: whole numbers as they are, others to two decimals (0.98), very small or large ones in scientific form.
pub fn format_number(v: f64) -> String {
    if v != 0.0 && (v.abs() < 1e-3 || v.abs() >= 1e6) {
        format!("{v:.2e}")
    } else if (v - v.round()).abs() < 1e-9 {
        format!("{v:.0}")
    } else {
        format!("{v:.2}")
    }
}

fn label_text(t: &NwTree, i: usize, keys: &[String]) -> String {
    let parts: Vec<String> = keys
        .iter()
        .filter_map(|k| {
            let v = raw(t, i, k)?;
            Some(
                v.parse::<f64>()
                    .ok()
                    .filter(|x| x.is_finite())
                    .map_or(v, format_number),
            )
        })
        .collect();
    parts.join("/")
}

/// The decoration of each node of `laid`, which was laid out from `t` (or from `t.to_node()`: the same tree, in the same order).
pub fn decorate(t: &NwTree, laid: &Laid, o: &StyleOptions) -> Decor {
    let order = t.preorder();
    let of = |i: usize| laid.nodes[i].origin.min(order.len().saturating_sub(1));
    let at = |i: usize| order[of(i)];
    let n = laid.nodes.len();
    let mut d = Decor {
        colour: vec![None; n],
        width: vec![1.0; n],
        label: vec![String::new(); n],
        shape: vec![None; n],
        bar: vec![None; n],
        legend: None,
    };
    if t.nodes.is_empty() || n == 0 {
        return d;
    }
    if let Some(key) = &o.colour_by {
        colour(t, key, &mut d, &at);
    }
    if let Some(key) = &o.width_by {
        if let Values::Numbers { low, high } = survey(t, key) {
            for (i, width) in d.width.iter_mut().enumerate() {
                if let Some(v) = t.attr_number(at(i), key) {
                    *width = 1.0 + (MOST_WIDTH - 1.0) * along(v, low, high) as f32;
                }
            }
        }
    }
    if !o.labels.is_empty() {
        for i in 0..n {
            if !laid.nodes[i].is_tip() || laid.nodes[i].cartoon.is_some() {
                d.label[i] = label_text(t, at(i), &o.labels);
            }
        }
    }
    if let Some(key) = &o.shape_by {
        shapes(t, key, &mut d, &at);
    }
    if let Some(key) = &o.bars_by {
        let deepest = laid.nodes.iter().map(|p| p.radius).fold(0.0, f32::max);
        for (i, bar) in d.bar.iter_mut().enumerate() {
            *bar = t
                .attr_range(at(i), key)
                .map(|(lo, hi)| (deepest - hi as f32, deepest - lo as f32));
        }
    }
    d
}

fn colour(t: &NwTree, key: &str, d: &mut Decor, at: &dyn Fn(usize) -> usize) {
    match survey(t, key) {
        Values::Numbers { low, high } => {
            for (i, slot) in d.colour.iter_mut().enumerate() {
                if let Some(v) = t.attr_number(at(i), key) {
                    *slot = Some(gradient(along(v, low, high)));
                }
            }
            d.legend = Some(Legend::Gradient {
                key: key.to_string(),
                low,
                high,
            });
        }
        Values::Words(words) => {
            let tone =
                |w: &str| PALETTE[words.iter().position(|x| x == w).unwrap_or(0) % PALETTE.len()];
            for (i, slot) in d.colour.iter_mut().enumerate() {
                if let Some(v) = t.attr(at(i), key) {
                    *slot = Some(tone(v));
                }
            }
            let items = words.iter().map(|w| (w.clone(), tone(w))).collect();
            d.legend = Some(Legend::Categories {
                key: key.to_string(),
                items,
            });
        }
        Values::None => {}
    }
}

fn shapes(t: &NwTree, key: &str, d: &mut Decor, at: &dyn Fn(usize) -> usize) {
    let kind = survey(t, key);
    for (i, slot) in d.shape.iter_mut().enumerate() {
        *slot = match (&kind, t.attr(at(i), key)) {
            (Values::Numbers { low, high }, Some(v)) => v.parse::<f64>().ok().map(|v| {
                SHAPE_RADII.0 + (SHAPE_RADII.1 - SHAPE_RADII.0) * along(v, *low, *high) as f32
            }),
            (Values::Words(_), Some(_)) => Some(SHAPE_PLAIN),
            _ => None,
        };
    }
}

/// What the numbers on the axis under a tree are.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum AxisMode {
    /// The distance from the root: 0 at the root.
    #[default]
    FromRoot,
    /// The time before the tips (the height of a node): 0 at the tips, growing towards the root.
    BeforeTips,
    /// A date: the tips are at the date given, and a node is as long before it as its distance from the tips says.
    Dates,
}

/// The axis under a tree (FigTree's time scale and scale axis): which numbers it shows, the scale of its units to the branch lengths, and
/// whether lines run up from its ticks.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Axis {
    pub mode: AxisMode,
    /// How many units of the axis a branch length of 1 is (a rate of substitutions per year gives years per length).
    pub scale: f64,
    /// The date of the tips, for `AxisMode::Dates`.
    pub tip_date: f64,
    pub grid: bool,
}

impl Default for Axis {
    fn default() -> Self {
        Axis {
            mode: AxisMode::FromRoot,
            scale: 1.0,
            tip_date: 0.0,
            grid: false,
        }
    }
}

/// The step between ticks: 1, 2 or 5 times a power of ten, the smallest that gives at most `target` ticks over `span`.
fn tick_step(span: f64, target: usize) -> f64 {
    if !(span.is_finite() && span > 0.0) || target == 0 {
        return 1.0;
    }
    let raw = span / target as f64;
    let power = 10f64.powf(raw.log10().floor());
    [1.0, 2.0, 5.0, 10.0]
        .into_iter()
        .map(|m| m * power)
        .find(|&s| s >= raw * (1.0 - 1e-9))
        .unwrap_or(power * 10.0)
}

impl Axis {
    /// The number on the axis at distance `x` from the root, in a tree whose tips are farthest at `deepest`.
    pub fn value(&self, x: f32, deepest: f32) -> f64 {
        let (x, deepest) = (f64::from(x), f64::from(deepest));
        match self.mode {
            AxisMode::FromRoot => x * self.scale,
            AxisMode::BeforeTips => (deepest - x) * self.scale,
            AxisMode::Dates => self.tip_date - (deepest - x) * self.scale,
        }
    }

    /// The distance from the root at which the axis shows `value`: the other way round from `value`.
    pub fn distance_of(&self, value: f64, deepest: f32) -> f32 {
        let deepest = f64::from(deepest);
        let scale = if self.scale == 0.0 { 1.0 } else { self.scale };
        (match self.mode {
            AxisMode::FromRoot => value / scale,
            AxisMode::BeforeTips => deepest - value / scale,
            AxisMode::Dates => deepest - (self.tip_date - value) / scale,
        }) as f32
    }

    /// The ticks along the axis: where each is (distance from the root) and the number written at it, about `target` of them, at round numbers.
    pub fn ticks(&self, deepest: f32, target: usize) -> Vec<(f32, String)> {
        if !(deepest.is_finite() && deepest > 0.0) {
            return Vec::new();
        }
        let (a, b) = (self.value(0.0, deepest), self.value(deepest, deepest));
        let (low, high) = (a.min(b), a.max(b));
        let step = tick_step(high - low, target);
        let decimals = (-step.log10().floor()).max(0.0) as usize;
        let first = (low / step - 1e-9).ceil() as i64;
        let last = (high / step + 1e-9).floor() as i64;
        (first..=last)
            .take(200)
            .map(|k| {
                let v = k as f64 * step;
                (self.distance_of(v, deepest), format!("{v:.decimals$}"))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests;
