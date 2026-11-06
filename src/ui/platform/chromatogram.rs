//! The chromatogram, drawn by the rules of the original's trace window (`new_abi_window`, AppMain.tcl:30863, 0x16b63a):
//! the bases the machine called along the top with a number above every tenth one, under them the four colour traces
//! scaled so that the tallest peak of all four just reaches the top of the trace area, and behind them the graph of
//! how sure each call was (`abi_draw_phred`, 0x173324). Reverse-complementing mirrors the whole picture and
//! complements the letters (`abi_revcom`, 0x1754b2). Only the scans on screen are drawn; bases can be picked with the
//! pointer, and the weak ends that trimming would drop are shaded.

use crate::seq::abif::Trace;
use crate::seq::dna::complement;
use crate::ui::design::{self, Look};
use crate::ui::theme::mix;
use eframe::egui::{
    pos2, vec2, Align2, Color32, FontId, Rect, ScrollArea, Sense, Shape, Stroke, Ui,
};

/// The trace area starts this far under the top of the window: the original moves the traces to y = 35 and scales them to
/// the height of the canvas less 45 (`new_abi_window`, AppMain.tcl:31110-31111, 0x16b63a); the base numbers and the base
/// letters share the space above.
const TRACE_TOP: f32 = 35.0;
const RULER_HEIGHT: f32 = 14.0;
const LETTERS_HEIGHT: f32 = TRACE_TOP - RULER_HEIGHT;
/// Room left under the baseline: the 45 less the 35 above (the original leaves ten points).
const BOTTOM_MARGIN: f32 = 10.0;
const MIN_TRACE_HEIGHT: f32 = 120.0;
const SIDE: f32 = 24.0;
/// A base number on the ruler for every this many bases.
const RULER_STEP: usize = 10;
/// The quality that fills the whole height of the quality graph: the original divides each quality by 100.
const QUALITY_FULL: f32 = 100.0;
/// The quality graph is #8080FF mixed into the page, 20 parts in a hundred for its fill and 35 for its edge.
const QUALITY_COLOUR: Color32 = Color32::from_rgb(0x80, 0x80, 0xff);
const QUALITY_FILL: f32 = 20.0;
const QUALITY_EDGE: f32 = 35.0;

pub struct View<'a> {
    pub trace: &'a Trace,
    pub look: &'a Look,
    /// Points on the screen for one scan: the original's horizontal scale, 0.5 to 5 and 1 when a trace opens.
    pub x_scale: f32,
    /// The original's vertical scale, 1 to 5: at 1 the tallest peak of the four traces reaches the top of the trace
    /// area, at 5 everything stands five times as high.
    pub y_scale: f32,
    /// Which of the A, C, G, T channels are drawn (by the letter the file gives the channel).
    pub shown: [bool; 4],
    /// The bases kept when trimming (`[from, to)`), or none to shade nothing.
    pub kept: Option<(usize, usize)>,
    /// The quality graph behind the traces (the original's "Toggle Phred").
    pub quality_graph: bool,
    /// Mirror the picture and complement the letters.
    pub reverse: bool,
}

/// What the pointer did over the chromatogram this frame.
#[derive(Default)]
pub struct Touched {
    /// The bases picked by dragging across them, as `[from, to)`.
    pub picked: Option<(usize, usize)>,
    /// A click that hit no base: the pick is cleared.
    pub cleared: bool,
}

/// The base whose peak is nearest to scan `x` (the peaks run in order).
pub fn base_at(peaks: &[u32], x: f32) -> Option<usize> {
    if peaks.is_empty() {
        return None;
    }
    let after = peaks.partition_point(|&p| (p as f32) < x);
    let before = after.saturating_sub(1);
    [before, after.min(peaks.len() - 1)]
        .into_iter()
        .min_by(|&a, &b| {
            (peaks[a] as f32 - x)
                .abs()
                .total_cmp(&(peaks[b] as f32 - x).abs())
        })
}

/// Where a scan is drawn along the trace: where it is, or mirrored about the middle of the trace when the view is
/// reverse-complemented.
pub fn mirrored(scan: f32, scans: usize, reverse: bool) -> f32 {
    if reverse {
        scans.saturating_sub(1) as f32 - scan
    } else {
        scan
    }
}

/// A base as it is shown: itself, or its complement in a reverse-complemented view (with the original's IUPAC pairs).
pub fn shown_base(base: u8, reverse: bool) -> u8 {
    if reverse {
        complement(base)
    } else {
        base
    }
}

/// How much one unit of signal is drawn: the tallest signal of any channel fills the trace area at scale 1.
pub fn trace_scale(trace: &Trace, trace_height: f32, y_scale: f32) -> f32 {
    y_scale * trace_height / f32::from(trace.tallest())
}

/// One step of the quality graph: the stretch of scans a base stands for and how high it reaches, 0 to 1.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Step {
    pub from: f32,
    pub to: f32,
    pub level: f32,
}

/// The steps of `abi_draw_phred`: each base's step runs from halfway to the base before to halfway to the base after,
/// the first from the start of the trace and the last to its end, and stands as high as the base's quality out of 100
/// (a quality above 100 is drawn as 100).
pub fn quality_steps(peaks: &[u32], quality: &[u8], scans: usize) -> Vec<Step> {
    let n = peaks.len().min(quality.len());
    let end = scans.saturating_sub(1) as f32;
    let middle = |a: u32, b: u32| (a as f32 + b as f32) / 2.0;
    (0..n)
        .map(|i| {
            let from = if i == 0 {
                0.0
            } else {
                middle(peaks[i - 1], peaks[i])
            };
            let to = if i + 1 == n {
                end
            } else {
                middle(peaks[i], peaks[i + 1])
            };
            Step {
                from,
                to: to.max(from),
                level: (f32::from(quality[i]) / QUALITY_FULL).min(1.0),
            }
        })
        .collect()
}

pub fn show(ui: &mut Ui, v: &View, selection: Option<(usize, usize)>) -> Touched {
    let t = v.trace;
    let look = v.look;
    let scans = t.scans();
    let height = ui
        .available_size()
        .y
        .max(RULER_HEIGHT + LETTERS_HEIGHT + MIN_TRACE_HEIGHT);
    let letters_bottom = RULER_HEIGHT + LETTERS_HEIGHT;
    let trace_height = height - letters_bottom - BOTTOM_MARGIN;
    let baseline = letters_bottom + trace_height;
    let scale = trace_scale(t, trace_height, v.y_scale);
    let width = SIDE * 2.0 + scans as f32 * v.x_scale;
    let mut touched = Touched::default();
    ScrollArea::horizontal()
        .auto_shrink([false, false])
        .show_viewport(ui, |ui, view| {
            let (rect, response) =
                ui.allocate_exact_size(vec2(width, height), Sense::click_and_drag());
            let p = ui.painter_at(rect);
            // a scan, as it is drawn, on the screen
            let place = |scan: f32, y: f32| {
                let drawn = mirrored(scan, scans, v.reverse);
                pos2(rect.left() + SIDE + drawn * v.x_scale, rect.top() + y)
            };
            let scan_under =
                |x: f32| mirrored((x - rect.left() - SIDE) / v.x_scale, scans, v.reverse);
            // the scans that are on screen, with one to spare
            let (left, right) = (
                (view.min.x - SIDE) / v.x_scale,
                (view.max.x - SIDE) / v.x_scale,
            );
            let (low, high) = if v.reverse {
                let end = scans.saturating_sub(1) as f32;
                (end - right, end - left)
            } else {
                (left, right)
            };
            let (first, last) = (
                low.floor().max(0.0) as usize,
                (high.ceil().max(0.0) as usize + 1).min(scans),
            );

            // the bases that are not kept, and the bases picked
            let band = |from_scan: f32, to_scan: f32, colour: Color32| {
                let (a, b) = (place(from_scan, 0.0).x, place(to_scan, 0.0).x);
                p.rect_filled(
                    Rect::from_min_max(
                        pos2(a.min(b), rect.top() + RULER_HEIGHT),
                        pos2(a.max(b), rect.top() + baseline),
                    ),
                    0.0,
                    colour,
                );
            };
            let edge = |base: usize, right: bool| {
                let at_peak = t.peaks[base] as f32;
                let step = if base + 1 < t.peaks.len() {
                    t.peaks[base + 1] as f32 - at_peak
                } else {
                    8.0
                };
                at_peak + if right { step / 2.0 } else { -step / 2.0 }
            };
            if let Some((from, to)) = v.kept.filter(|&(f, e)| f < e && e <= t.peaks.len()) {
                let shade = look.dim.gamma_multiply(0.18);
                if from > 0 {
                    band(0.0, edge(from, false), shade);
                }
                if to < t.peaks.len() {
                    band(edge(to - 1, true), scans as f32, shade);
                }
            }
            if let Some((from, to)) = selection.filter(|&(f, e)| f < e && e <= t.peaks.len()) {
                band(
                    edge(from, false),
                    edge(to - 1, true),
                    look.accent.gamma_multiply(0.22),
                );
            }

            // the quality graph, behind the traces
            if v.quality_graph && t.has_quality() {
                let range = baseline - letters_bottom;
                let steps = quality_steps(&t.peaks, &t.quality, scans);
                let fill = mix(QUALITY_COLOUR, look.surface, QUALITY_FILL);
                for s in &steps {
                    let top = baseline - s.level * range;
                    let (a, b) = (place(s.from, 0.0).x, place(s.to, 0.0).x);
                    p.rect_filled(
                        Rect::from_min_max(
                            pos2(a.min(b), rect.top() + top),
                            pos2(a.max(b), rect.top() + baseline),
                        ),
                        0.0,
                        fill,
                    );
                }
                let mut outline = vec![place(0.0, baseline)];
                for s in &steps {
                    let top = baseline - s.level * range;
                    outline.push(place(s.from, top));
                    outline.push(place(s.to, top));
                }
                outline.push(place(scans.saturating_sub(1) as f32, baseline));
                p.add(Shape::line(
                    outline,
                    Stroke::new(1.0, mix(QUALITY_COLOUR, look.surface, QUALITY_EDGE)),
                ));
            }

            // the traces, each in the colour of the base it stands for as shown
            let step = ((1.0 / v.x_scale).floor() as usize).max(1);
            for (k, channel) in t.channels.iter().enumerate() {
                let base = t.order[k];
                let slot = b"ACGT".iter().position(|&b| b == base);
                if !slot.is_some_and(|s| v.shown[s]) {
                    continue;
                }
                let points: Vec<_> = (first..last.min(channel.len()))
                    .step_by(step)
                    .map(|s| place(s as f32, baseline - f32::from(channel[s].max(0)) * scale))
                    .collect();
                p.add(Shape::line(
                    points,
                    Stroke::new(
                        1.5,
                        design::trace_color(shown_base(base, v.reverse), look.dark),
                    ),
                ));
            }
            p.line_segment(
                [place(0.0, baseline), place(scans as f32, baseline)],
                Stroke::new(1.0, look.border),
            );

            // the bases, numbered along the picture from its left
            let letters = v.x_scale * 8.0 >= 7.0;
            let n = t.bases.len().min(t.peaks.len());
            for (i, (&base, &scan)) in t.bases.iter().zip(&t.peaks).enumerate() {
                let x = place(scan as f32, 0.0).x;
                if x < rect.left() + view.min.x - 20.0 || x > rect.left() + view.max.x + 20.0 {
                    continue;
                }
                let shown = shown_base(base, v.reverse);
                let colour = design::trace_color(shown, look.dark);
                let centre = pos2(x, rect.top() + RULER_HEIGHT + LETTERS_HEIGHT / 2.0);
                if letters {
                    p.text(
                        centre,
                        Align2::CENTER_CENTER,
                        char::from(shown).to_string(),
                        FontId::monospace(14.0),
                        colour,
                    );
                } else {
                    p.circle_filled(centre, 1.6, colour);
                }
                let number = if v.reverse { n - i } else { i + 1 };
                if number % RULER_STEP == 0 {
                    p.line_segment(
                        [
                            pos2(x, rect.top() + RULER_HEIGHT - 5.0),
                            pos2(x, rect.top() + RULER_HEIGHT - 1.0),
                        ],
                        Stroke::new(1.0, look.dim),
                    );
                    p.text(
                        pos2(x, rect.top() + RULER_HEIGHT / 2.0 - 2.0),
                        Align2::CENTER_CENTER,
                        number.to_string(),
                        FontId::proportional(10.0),
                        look.dim,
                    );
                }
            }

            // picking bases by dragging across them
            if let Some(pointer) = response.interact_pointer_pos() {
                if let Some(here) = base_at(&t.peaks, scan_under(pointer.x)) {
                    let origin = ui
                        .ctx()
                        .input(|i| i.pointer.press_origin())
                        .and_then(|o| base_at(&t.peaks, scan_under(o.x)))
                        .unwrap_or(here);
                    if response.dragged() || response.drag_stopped() {
                        touched.picked = Some((origin.min(here), origin.max(here) + 1));
                    }
                }
            }
            if response.clicked() {
                if let Some(pointer) = response.interact_pointer_pos() {
                    touched.picked = base_at(&t.peaks, scan_under(pointer.x)).map(|b| (b, b + 1));
                    touched.cleared = touched.picked.is_none();
                }
            }
            if let Some(pointer) = response.hover_pos() {
                if let Some(b) = base_at(&t.peaks, scan_under(pointer.x)) {
                    let quality = if t.has_quality() {
                        format!(", quality {}", t.quality.get(b).copied().unwrap_or(0))
                    } else {
                        String::new()
                    };
                    let number = if v.reverse { n - b } else { b + 1 };
                    response.on_hover_text(format!(
                        "base {number} {}{quality}, scan {}",
                        char::from(shown_base(t.bases[b], v.reverse)),
                        t.peaks[b]
                    ));
                }
            }
        });
    touched
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_base_nearest_to_a_scan_is_found_at_either_end_and_between() {
        let peaks = [10u32, 22, 35, 47];
        assert_eq!(base_at(&peaks, -5.0), Some(0));
        assert_eq!(base_at(&peaks, 10.0), Some(0));
        assert_eq!(base_at(&peaks, 15.0), Some(0));
        assert_eq!(base_at(&peaks, 17.0), Some(1));
        assert_eq!(base_at(&peaks, 41.0), Some(2));
        assert_eq!(base_at(&peaks, 44.0), Some(3));
        assert_eq!(base_at(&peaks, 900.0), Some(3));
        assert_eq!(base_at(&[], 3.0), None);
    }

    #[test]
    fn the_quality_graph_steps_at_the_middle_between_bases_as_the_original_draws_it() {
        let steps = quality_steps(&[10, 22, 35], &[40, 100, 250], 50);
        assert_eq!(
            steps,
            vec![
                Step {
                    from: 0.0,
                    to: 16.0,
                    level: 0.4
                },
                Step {
                    from: 16.0,
                    to: 28.5,
                    level: 1.0
                },
                // the last runs to the end of the trace, and a quality above 100 stands as 100
                Step {
                    from: 28.5,
                    to: 49.0,
                    level: 1.0
                },
            ]
        );
        // the steps meet, so the graph has no gap
        assert!(steps.windows(2).all(|w| w[0].to == w[1].from));
        assert_eq!(quality_steps(&[], &[], 50), vec![]);
        // one base stands for the whole trace
        assert_eq!(
            quality_steps(&[5], &[50], 20),
            vec![Step {
                from: 0.0,
                to: 19.0,
                level: 0.5
            }]
        );
    }

    #[test]
    fn reversing_mirrors_a_scan_about_the_middle_of_the_trace_and_complements_the_letters() {
        assert_eq!(mirrored(0.0, 101, false), 0.0);
        assert_eq!(mirrored(0.0, 101, true), 100.0);
        assert_eq!(mirrored(30.0, 101, true), 70.0);
        // twice is the way it was
        assert_eq!(mirrored(mirrored(37.0, 101, true), 101, true), 37.0);
        assert_eq!(shown_base(b'A', true), b'T');
        assert_eq!(shown_base(b'g', true), b'c');
        assert_eq!(shown_base(b'N', true), b'N');
        assert_eq!(shown_base(b'A', false), b'A');
    }

    #[test]
    fn the_tallest_peak_of_all_four_traces_fills_the_trace_area_and_the_scale_multiplies_it() {
        let trace = Trace {
            channels: [vec![10, 400, 20], vec![0, 100, 0], vec![0; 3], vec![0; 3]],
            ..Trace::default()
        };
        assert_eq!(trace_scale(&trace, 200.0, 1.0), 0.5);
        assert_eq!(trace_scale(&trace, 200.0, 3.0), 1.5);
        // a trace with no signal does not divide by zero
        assert!(trace_scale(&Trace::default(), 200.0, 1.0).is_finite());
    }
}
