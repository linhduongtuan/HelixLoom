//! Contig viewer: position ruler, coverage plot, colour-coded consensus and the stacked reads.

use crate::seq::assemble::{Contig, Placed};
use crate::seq::dna::{base_matches, revcom};
use crate::ui::map::commify;
use crate::ui::theme::{self, Palette};
use eframe::egui::{
    self, pos2, vec2, Align2, Color32, FontId, Rect, ScrollArea, Sense, Shape, Stroke, Ui,
};

const RULER_H: f32 = 26.0;
const COVERAGE_H: f32 = 48.0;
const CONSENSUS_H: f32 = 22.0;
/// Space between the coverage plot and the consensus row (holds the row's label).
const CONSENSUS_GAP: f32 = 14.0;
const LEFT_PAD: f32 = 12.0;
/// Zoom (pixels per base) at which letters replace coloured blocks.
const LETTER_ZOOM: f32 = 7.0;
/// Below this zoom single bases can't be told apart, so the consensus is drawn as one bar.
const BLOCK_ZOOM: f32 = 2.0;
pub const MIN_ZOOM: f32 = 0.05;
pub const MAX_ZOOM: f32 = 24.0;

pub struct ContigView {
    pub zoom: f32,
    /// Index into `contig.placed`.
    pub selected: Option<usize>,
    /// Fit the whole contig into the window the next time it is drawn.
    pub fit_pending: bool,
}

impl Default for ContigView {
    fn default() -> Self {
        ContigView {
            zoom: 1.0,
            selected: None,
            fit_pending: true,
        }
    }
}

/// Row for each read when reads are stacked without touching (two bases apart); also the row count.
pub fn pack_rows(spans: &[(usize, usize)]) -> (Vec<usize>, usize) {
    let mut ends: Vec<usize> = Vec::new();
    let rows = spans
        .iter()
        .map(|&(s, e)| match ends.iter().position(|&end| end + 2 <= s) {
            Some(r) => {
                ends[r] = e;
                r
            }
            None => {
                ends.push(e);
                ends.len() - 1
            }
        })
        .collect();
    (rows, ends.len())
}

pub fn base_color(b: u8) -> Color32 {
    match b.to_ascii_uppercase() {
        b'A' => Color32::from_rgb(52, 160, 64),
        b'C' => Color32::from_rgb(58, 100, 220),
        b'G' => Color32::from_rgb(214, 140, 20),
        b'T' => Color32::from_rgb(210, 58, 58),
        _ => Color32::from_gray(130),
    }
}

/// The read as it lies in the contig (reverse-complemented if it was joined the other way round).
pub fn oriented(seq: &[u8], placed: &Placed) -> Vec<u8> {
    if placed.reverse {
        revcom(seq)
    } else {
        seq.to_vec()
    }
}

/// `(mismatches, bases compared)` of a read against the consensus; N bases and ambiguous
/// consensus calls that allow the base are not counted as mismatches.
pub fn compare_to_consensus(contig: &Contig, placed: &Placed, seq: &[u8]) -> (usize, usize) {
    let bases = oriented(seq, placed);
    let mut compared = 0;
    let mut mismatches = 0;
    for (k, &b) in bases.iter().enumerate() {
        if b == b'N' {
            continue;
        }
        compared += 1;
        if !base_matches(contig.consensus_at(placed.offset + k), b, false) {
            mismatches += 1;
        }
    }
    (mismatches, compared)
}

/// Ruler spacing: the smallest round step that keeps labels at least 70 px apart.
fn ruler_step(zoom: f32) -> usize {
    [
        1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10_000, 20_000, 50_000,
    ]
    .into_iter()
    .find(|&s| s as f32 * zoom >= 70.0)
    .unwrap_or(100_000)
}

struct Geometry {
    zoom: f32,
    row_h: f32,
    font: FontId,
    letters: bool,
}

pub fn show(
    ui: &mut Ui,
    view: &mut ContigView,
    contig: &Contig,
    inputs: &[(String, Vec<u8>)],
    pal: &Palette,
) {
    let ratio = ui
        .ctx()
        .fonts_mut(|f| f.glyph_width(&FontId::monospace(10.0), 'A'))
        / 10.0;
    if view.fit_pending {
        view.zoom = ((ui.available_width() - 2.0 * LEFT_PAD - 40.0) / contig.span.max(1) as f32)
            .clamp(MIN_ZOOM, MAX_ZOOM);
        view.fit_pending = false;
    }
    let zoom = view.zoom.clamp(MIN_ZOOM, MAX_ZOOM);
    let letters = zoom >= LETTER_ZOOM;
    let geo = Geometry {
        zoom,
        row_h: if letters {
            (zoom / ratio * 1.25).clamp(17.0, 34.0)
        } else {
            17.0
        },
        font: FontId::monospace(zoom / ratio),
        letters,
    };
    let spans: Vec<(usize, usize)> = contig
        .placed
        .iter()
        .map(|p| (p.offset, p.offset + inputs[p.read].1.len()))
        .collect();
    let (rows, row_count) = pack_rows(&spans);
    let width = contig.span as f32 * zoom + LEFT_PAD * 2.0 + 40.0;
    let reads_top = RULER_H + COVERAGE_H + CONSENSUS_GAP + CONSENSUS_H + 10.0;
    let height = reads_top + row_count as f32 * geo.row_h + 24.0;

    ScrollArea::both()
        .auto_shrink([false, false])
        .show_viewport(ui, |ui, vp| {
            let (rect, _) = ui.allocate_exact_size(
                vec2(width.max(vp.width()), height.max(vp.height())),
                Sense::hover(),
            );
            let p = ui.painter().with_clip_rect(rect.intersect(ui.clip_rect()));
            p.rect_filled(rect, 0.0, pal.text_bg);
            let x_of = |col: f32| rect.left() + LEFT_PAD + col * zoom;
            let first =
                (((vp.min.x - LEFT_PAD) / zoom).floor().max(0.0) as usize).saturating_sub(1);
            let last =
                ((((vp.max.x - LEFT_PAD) / zoom).ceil().max(0.0)) as usize + 1).min(contig.span);
            let ctx = Ctx {
                p: &p,
                rect,
                pal,
                geo: &geo,
                first,
                last,
                x_of: &x_of,
            };
            ruler(&ctx);
            coverage(&ctx, contig);
            consensus_row(&ctx, contig);
            let mut clicked = None;
            for (k, pl) in contig.placed.iter().enumerate() {
                let seq = &inputs[pl.read].1;
                let (x0, x1) = (x_of(pl.offset as f32), x_of((pl.offset + seq.len()) as f32));
                if x1 < rect.left() + vp.min.x || x0 > rect.left() + vp.max.x {
                    continue;
                }
                let bar = Rect::from_min_size(
                    pos2(x0, rect.top() + reads_top + rows[k] as f32 * geo.row_h),
                    vec2(x1 - x0, geo.row_h - 2.0),
                );
                read_bar(
                    &ctx,
                    contig,
                    pl,
                    &inputs[pl.read],
                    bar,
                    view.selected == Some(k),
                );
                let resp = ui.interact(bar, ui.id().with(("read", k)), Sense::click());
                if resp.clicked() {
                    clicked = Some(k);
                }
                resp.on_hover_ui(|ui| {
                    let (mism, compared) = compare_to_consensus(contig, pl, seq);
                    ui.label(egui::RichText::new(&inputs[pl.read].0).strong());
                    ui.label(format!(
                        "{} • {}–{} • {} bp{}",
                        if pl.reverse {
                            "reverse complement"
                        } else {
                            "forward"
                        },
                        commify(pl.offset + 1),
                        commify(pl.offset + seq.len()),
                        commify(seq.len()),
                        if pl.contained {
                            " • inside a longer read"
                        } else {
                            ""
                        }
                    ));
                    ui.label(format!(
                        "{mism} mismatch{} in {} compared bases",
                        if mism == 1 { "" } else { "es" },
                        commify(compared)
                    ));
                });
            }
            if contig.circular {
                origin_marker(&ctx, contig, height);
            }
            if let Some(k) = clicked {
                view.selected = Some(k);
            }
        });
}

struct Ctx<'a> {
    p: &'a egui::Painter,
    rect: Rect,
    pal: &'a Palette,
    geo: &'a Geometry,
    first: usize,
    last: usize,
    x_of: &'a dyn Fn(f32) -> f32,
}

fn ruler(c: &Ctx) {
    let step = ruler_step(c.geo.zoom);
    let y = c.rect.top() + RULER_H;
    let mut pos = (c.first / step) * step;
    while pos <= c.last {
        let x = (c.x_of)(pos as f32);
        c.p.line_segment([pos2(x, y - 7.0), pos2(x, y)], Stroke::new(1.0, c.pal.fg));
        c.p.text(
            pos2(x + 2.0, y - 9.0),
            Align2::LEFT_BOTTOM,
            commify(pos.max(1)),
            FontId::proportional(10.0),
            c.pal.dim,
        );
        pos += step;
    }
    c.p.line_segment(
        [pos2(c.rect.left(), y), pos2(c.rect.right(), y)],
        Stroke::new(1.0, c.pal.line),
    );
}

fn coverage(c: &Ctx, contig: &Contig) {
    let top = c.rect.top() + RULER_H + 2.0;
    let bottom = top + COVERAGE_H - 6.0;
    let max = contig.depth.iter().copied().max().unwrap_or(1).max(1) as f32;
    let fill = c.pal.accent.gamma_multiply(0.45);
    for col in c.first..c.last.min(contig.depth.len()) {
        let h = contig.depth[col] as f32 / max * (bottom - top);
        let x = (c.x_of)(col as f32);
        c.p.rect_filled(
            Rect::from_min_max(pos2(x, bottom - h), pos2(x + c.geo.zoom.max(1.0), bottom)),
            0.0,
            fill,
        );
    }
    c.p.line_segment(
        [pos2(c.rect.left(), bottom), pos2(c.rect.right(), bottom)],
        Stroke::new(1.0, c.pal.line),
    );
    c.p.text(
        pos2(c.rect.left() + 4.0, top),
        Align2::LEFT_TOP,
        format!("coverage (max {}×)", max as u32),
        FontId::proportional(10.0),
        c.pal.dim,
    );
}

fn consensus_row(c: &Ctx, contig: &Contig) {
    let top = c.rect.top() + RULER_H + COVERAGE_H + CONSENSUS_GAP;
    c.p.text(
        pos2(c.rect.left() + 4.0, top - 2.0),
        Align2::LEFT_BOTTOM,
        "consensus (red marks: reads disagree)",
        FontId::proportional(10.0),
        c.pal.dim,
    );
    if c.geo.zoom < BLOCK_ZOOM {
        let (x0, x1) = ((c.x_of)(c.first as f32), (c.x_of)(c.last as f32));
        c.p.rect_filled(
            Rect::from_min_max(pos2(x0, top), pos2(x1, top + CONSENSUS_H - 2.0)),
            0.0,
            c.pal.dim.gamma_multiply(0.55),
        );
    }
    for col in c.first..c.last {
        let b = contig.consensus_at(col);
        let x = (c.x_of)(col as f32);
        let cell = Rect::from_min_size(pos2(x, top), vec2(c.geo.zoom, CONSENSUS_H - 2.0));
        if c.geo.letters {
            c.p.text(
                cell.center(),
                Align2::CENTER_CENTER,
                (b as char).to_string(),
                c.geo.font.clone(),
                base_color(b),
            );
        } else if c.geo.zoom >= BLOCK_ZOOM {
            c.p.rect_filled(cell, 0.0, base_color(b).gamma_multiply(0.85));
        }
        if contig
            .disagree
            .get(col % contig.consensus.len().max(1))
            .copied()
            .unwrap_or(false)
        {
            c.p.rect_filled(
                Rect::from_min_size(pos2(x, top - 1.0), vec2(c.geo.zoom.max(2.0), CONSENSUS_H)),
                0.0,
                Color32::from_rgb(220, 40, 40).gamma_multiply(if c.geo.zoom < BLOCK_ZOOM {
                    1.0
                } else {
                    0.35
                }),
            );
        }
    }
}

fn read_bar(
    c: &Ctx,
    contig: &Contig,
    pl: &Placed,
    input: &(String, Vec<u8>),
    bar: Rect,
    selected: bool,
) {
    let dark = c.pal.dark;
    let tint = match (pl.reverse, dark) {
        (false, false) => Color32::from_rgb(196, 218, 255),
        (true, false) => Color32::from_rgb(255, 219, 184),
        (false, true) => Color32::from_rgb(48, 78, 128),
        (true, true) => Color32::from_rgb(128, 84, 44),
    };
    let fill = if pl.contained {
        theme::mix(tint, c.pal.text_bg, 60.0)
    } else {
        tint
    };
    let head = (bar.height() * 0.5).min(bar.width());
    let mid = bar.center().y;
    let pts = if pl.reverse {
        vec![
            pos2(bar.left() + head, bar.top()),
            pos2(bar.right(), bar.top()),
            pos2(bar.right(), bar.bottom()),
            pos2(bar.left() + head, bar.bottom()),
            pos2(bar.left(), mid),
        ]
    } else {
        vec![
            pos2(bar.left(), bar.top()),
            pos2(bar.right() - head, bar.top()),
            pos2(bar.right(), mid),
            pos2(bar.right() - head, bar.bottom()),
            pos2(bar.left(), bar.bottom()),
        ]
    };
    let edge = if selected {
        Stroke::new(2.0, c.pal.accent)
    } else {
        Stroke::new(1.0, theme::mix(Color32::BLACK, fill, 25.0))
    };
    c.p.add(Shape::convex_polygon(pts, fill, edge));
    let seq = oriented(&input.1, pl);
    if c.geo.letters {
        let a = c.first.saturating_sub(pl.offset);
        let b = c.last.saturating_sub(pl.offset).min(seq.len());
        let disagrees = |k: usize| {
            seq[k] != b'N' && !base_matches(contig.consensus_at(pl.offset + k), seq[k], false)
        };
        for (k, &base) in seq
            .iter()
            .enumerate()
            .take(b)
            .skip(a)
            .filter(|(k, _)| disagrees(*k))
        {
            let x = bar.left() + k as f32 * c.geo.zoom;
            c.p.rect_filled(
                Rect::from_min_size(pos2(x, bar.top()), vec2(c.geo.zoom, bar.height())),
                0.0,
                Color32::from_rgb(214, 48, 48),
            );
            c.p.text(
                pos2(x + c.geo.zoom / 2.0, mid),
                Align2::CENTER_CENTER,
                (base as char).to_string(),
                c.geo.font.clone(),
                Color32::WHITE,
            );
        }
        // matching bases in one piece of text; disagreeing columns are blank (drawn above)
        let text: String = (a..b)
            .map(|k| if disagrees(k) { ' ' } else { seq[k] as char })
            .collect();
        c.p.text(
            pos2(bar.left() + a as f32 * c.geo.zoom, mid),
            Align2::LEFT_CENTER,
            text,
            c.geo.font.clone(),
            theme::text_on(fill).gamma_multiply(0.8),
        );
    } else {
        let a = c.first.saturating_sub(pl.offset);
        let b = c.last.saturating_sub(pl.offset).min(seq.len());
        for (k, &base) in seq.iter().enumerate().take(b).skip(a) {
            if base != b'N' && !base_matches(contig.consensus_at(pl.offset + k), base, false) {
                let x = bar.left() + k as f32 * c.geo.zoom;
                c.p.rect_filled(
                    Rect::from_min_size(
                        pos2(x, bar.top()),
                        vec2(c.geo.zoom.max(2.0), bar.height()),
                    ),
                    0.0,
                    Color32::from_rgb(214, 48, 48),
                );
            }
        }
        if bar.width() > 60.0 {
            c.p.with_clip_rect(bar.shrink(1.0)).text(
                pos2(bar.left() + 4.0, mid),
                Align2::LEFT_CENTER,
                &input.0,
                FontId::proportional(10.0),
                theme::text_on(fill),
            );
        }
    }
}

/// Where the circle closes: everything to the right repeats the start of the consensus.
fn origin_marker(c: &Ctx, contig: &Contig, height: f32) {
    let x = (c.x_of)(contig.consensus.len() as f32);
    let top = c.rect.top() + RULER_H;
    c.p.line_segment(
        [pos2(x, top), pos2(x, c.rect.top() + height)],
        Stroke::new(1.5, Color32::from_rgb(200, 60, 60)),
    );
    c.p.rect_filled(
        Rect::from_min_max(pos2(x, top), pos2(c.rect.right(), c.rect.top() + height)),
        0.0,
        c.pal.dim.gamma_multiply(0.12),
    );
    c.p.text(
        pos2(x + 4.0, top + 2.0),
        Align2::LEFT_TOP,
        "circle closes here: bases to the right repeat the start",
        FontId::proportional(10.0),
        Color32::from_rgb(200, 60, 60),
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::assemble::{assemble, Params, Read};

    #[test]
    fn reads_are_stacked_without_touching() {
        let spans = [(0, 100), (50, 150), (110, 200), (160, 260), (300, 400)];
        let (rows, count) = pack_rows(&spans);
        assert_eq!(rows, vec![0, 1, 0, 1, 0]);
        assert_eq!(count, 2);
        assert_eq!(pack_rows(&[]), (vec![], 0));
        let (rows, _) = pack_rows(&[(0, 100), (101, 200)]);
        assert_eq!(
            rows,
            vec![0, 1],
            "reads one base apart still need separate rows"
        );
    }

    #[test]
    fn ruler_labels_stay_at_least_70_pixels_apart() {
        assert_eq!(ruler_step(24.0), 5);
        assert_eq!(ruler_step(1.0), 100);
        assert!(ruler_step(0.5) as f32 * 0.5 >= 70.0);
    }

    #[test]
    fn mismatches_are_counted_against_the_consensus() {
        let genome: Vec<u8> = b"ACGTTGCAAGCTTGACCATGGCATCGATCGGATCCATTAGC".repeat(4);
        let mut bad = genome[10..130].to_vec();
        bad[20] = if bad[20] == b'A' { b'C' } else { b'A' };
        let reads = vec![
            Read {
                name: "a".into(),
                seq: genome[..100].to_vec(),
            },
            Read {
                name: "b".into(),
                seq: genome[10..130].to_vec(),
            },
            Read {
                name: "c".into(),
                seq: genome[10..130].to_vec(),
            },
            Read {
                name: "bad".into(),
                seq: bad,
            },
        ];
        let asm = assemble(
            &reads,
            &Params {
                min_overlap: 20,
                ..Params::reads()
            },
        );
        assert_eq!(asm.contigs.len(), 1);
        let c = &asm.contigs[0];
        let inputs: Vec<(String, Vec<u8>)> = reads
            .iter()
            .map(|r| (r.name.clone(), crate::seq::assemble::normalize(&r.seq)))
            .collect();
        let by_name = |n: &str| c.placed.iter().find(|p| inputs[p.read].0 == n).unwrap();
        assert_eq!(compare_to_consensus(c, by_name("b"), &inputs[1].1).0, 0);
        let (mism, compared) = compare_to_consensus(c, by_name("bad"), &inputs[3].1);
        assert_eq!((mism, compared), (1, 120));
        assert_eq!(c.conflicts, 1);
        assert!(c.disagree.iter().filter(|&&d| d).count() == 1);
    }
}
