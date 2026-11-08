//! Digestion: fragment table and a gel picture next to a size ladder.

use crate::seq::enzymes::{collect_cuts, digest, Enzyme, Fragment, Site};
use crate::ui::app::App;
use crate::ui::map::commify;
use crate::ui::state::DocView;
use crate::ui::theme::{self, Palette};
use eframe::egui::{self, pos2, vec2, Align2, FontId, RichText, ScrollArea, Sense, Stroke, Ui};
use std::collections::BTreeSet;

/// Size ladders shipped with the original (DNA_ladders.txt and the NEB folder).
pub const LADDERS: &[(&str, &[usize])] = &[
    (
        "1Kb plus",
        &[
            12000, 11000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000, 2000, 1650, 1650, 1000,
            850, 650, 500, 400, 300, 200, 100,
        ],
    ),
    (
        "lambda Hind",
        &[23130, 9416, 6557, 4361, 2322, 2027, 564, 125],
    ),
    (
        "1 kb DNA Ladder (N3232)",
        &[
            10002, 8001, 6001, 5001, 4001, 3001, 3001, 3001, 2000, 1500, 1000, 517, 500,
        ],
    ),
    (
        "100 bp DNA Ladder (N3231)",
        &[
            1517, 1200, 1000, 1000, 1000, 900, 800, 700, 600, 517, 517, 517, 500, 500, 500, 400,
            300, 200, 100,
        ],
    ),
    (
        "50 bp DNA Ladder (N3236)",
        &[
            1350, 916, 766, 700, 650, 600, 550, 500, 500, 500, 450, 400, 350, 300, 250, 200, 200,
            200, 150, 100, 50,
        ],
    ),
];

/// The original's lowest band darkness for a faint fragment (gel_min_gray).
const GEL_MIN_GRAY: f32 = 10.0;
const GEL_HEIGHT: f32 = 210.0;
const BAND_WIDTH: f32 = 50.0;
const BAND_SPACING: f32 = 5.0;

#[derive(Default)]
pub struct DigestWin {
    pub open: bool,
    pub selection_only: bool,
    pub ladder: usize,
}

pub struct DigestResult {
    pub fragments: Vec<Fragment>,
    pub total: usize,
    /// Offset of the digested region (for a selection-only digest).
    pub origin: usize,
    pub enzymes: Vec<String>,
}

/// Digest the whole sequence, or just the selection (treated as a linear piece).
pub fn compute(
    dv: &mut DocView,
    enzymes: &[Enzyme],
    chosen: &BTreeSet<String>,
    selection_only: bool,
) -> DigestResult {
    let region = dv.selection().filter(|_| selection_only);
    let circular = dv.doc.circular && region.is_none();
    let (a, b) = region.unwrap_or((0, dv.doc.len()));
    let all_sites = dv.sites(enzymes);
    let picked: Vec<(&Enzyme, Vec<Site>)> = enzymes
        .iter()
        .zip(all_sites)
        .filter(|(e, _)| chosen.contains(&e.name))
        .map(|(e, sites)| {
            let inside: Vec<Site> = sites
                .iter()
                .filter(|s| {
                    s.cuts_inside
                        && (region.is_none()
                            || (s.top_cut >= a
                                && s.top_cut <= b
                                && s.bottom_cut >= a
                                && s.bottom_cut <= b))
                })
                .map(|s| Site {
                    start: s.start.saturating_sub(a),
                    top_cut: s.top_cut.saturating_sub(a),
                    bottom_cut: s.bottom_cut.saturating_sub(a),
                    ..s.clone()
                })
                .collect();
            (e, inside)
        })
        .collect();
    let refs: Vec<(&Enzyme, &[Site])> = picked.iter().map(|(e, s)| (*e, s.as_slice())).collect();
    DigestResult {
        fragments: digest(b - a, circular, &collect_cuts(&refs)),
        total: b - a,
        origin: a,
        enzymes: picked.iter().map(|(e, _)| e.name.clone()).collect(),
    }
}

/// Vertical position of a band: the original's log scale (larger fragments sit higher).
pub struct GelScale {
    scale: f32,
    origin: f32,
}

impl GelScale {
    pub fn new(shortest: usize, longest: usize) -> GelScale {
        let (s, l) = ((shortest.max(10)) as f32, longest.max(11) as f32);
        let scale = if l > s {
            GEL_HEIGHT / (l.log10() - s.log10())
        } else {
            20.0
        };
        GelScale {
            scale,
            origin: 30.0 + scale * l.log10(),
        }
    }

    pub fn y(&self, size: usize) -> f32 {
        self.origin - self.scale * (size.max(1) as f32).log10()
    }
}

fn draw_gel(ui: &mut Ui, pal: &Palette, ladder: &[usize], frags: &[Fragment]) {
    let (rect, resp) = ui.allocate_exact_size(
        vec2(2.0 * (BAND_WIDTH + BAND_SPACING) + 110.0, GEL_HEIGHT + 56.0),
        Sense::hover(),
    );
    let p = ui.painter_at(rect);
    p.rect_filled(rect, 2.0, pal.text_bg);
    p.rect_stroke(
        rect,
        2.0,
        Stroke::new(1.0, pal.line),
        egui::StrokeKind::Inside,
    );
    let smallest = frags
        .iter()
        .map(|f| f.length)
        .chain(ladder.last().copied())
        .min()
        .unwrap_or(10);
    let largest = frags
        .iter()
        .map(|f| f.length)
        .chain(ladder.first().copied())
        .max()
        .unwrap_or(100);
    let scale = GelScale::new(smallest, largest);
    let lane1 = rect.left() + 46.0;
    let lane2 = lane1 + BAND_WIDTH + BAND_SPACING;
    let font = FontId::proportional(10.0);
    p.text(
        pos2(lane1 + BAND_WIDTH / 2.0, rect.top() + 24.0),
        Align2::CENTER_BOTTOM,
        "1",
        FontId::proportional(12.0),
        pal.fg,
    );
    p.text(
        pos2(lane2 + BAND_WIDTH / 2.0, rect.top() + 24.0),
        Align2::CENTER_BOTTOM,
        "2",
        FontId::proportional(12.0),
        pal.fg,
    );
    let y0 = rect.top();
    let mut last: Option<(usize, f32)> = None;
    let mut last_label_y = f32::MIN;
    for &size in ladder {
        let y = y0 + scale.y(size);
        let thick = match last {
            Some((s, _)) if s == size => 5.0,
            _ => 3.0,
        };
        // name each ladder band once, skipping labels that would print on top of the previous one
        if last.is_none_or(|(s, _)| s != size) && y - last_label_y >= 9.0 {
            p.text(
                pos2(lane1 - 4.0, y),
                Align2::RIGHT_CENTER,
                size.to_string(),
                font.clone(),
                pal.dim,
            );
            last_label_y = y;
        }
        p.line_segment(
            [pos2(lane1, y), pos2(lane1 + BAND_WIDTH, y)],
            Stroke::new(thick, pal.fg),
        );
        last = Some((size, y));
    }
    let max_mass = frags.iter().map(|f| f.mass_pct).max().unwrap_or(1).max(1) as f32;
    for f in frags {
        let y = y0 + scale.y(f.length);
        let gray = (100.0 * f.mass_pct as f32 / max_mass).max(GEL_MIN_GRAY);
        let color = theme::mix(pal.fg, pal.text_bg, gray);
        p.line_segment(
            [pos2(lane2, y), pos2(lane2 + BAND_WIDTH, y)],
            Stroke::new(3.0, color),
        );
        p.text(
            pos2(lane2 + BAND_WIDTH + 5.0, y),
            Align2::LEFT_CENTER,
            commify(f.length),
            font.clone(),
            pal.fg,
        );
    }
    if let Some(pos) = resp.hover_pos() {
        let near = frags.iter().min_by(|a, b| {
            (y0 + scale.y(a.length) - pos.y)
                .abs()
                .total_cmp(&(y0 + scale.y(b.length) - pos.y).abs())
        });
        if let Some(f) =
            near.filter(|f| pos.x > lane2 && (y0 + scale.y(f.length) - pos.y).abs() < 6.0)
        {
            resp.on_hover_text(format!("{} bp ({}–{})", commify(f.length), f.left, f.right));
        }
    }
}

impl App {
    pub fn digest_window(&mut self, ctx: &egui::Context) {
        if !self.digest_win.open {
            return;
        }
        let Some(i) = self.current_doc() else {
            self.digest_win.open = false;
            return;
        };
        let pal = theme::palette(self.prefs.dark);
        let mut open = true;
        crate::ui::design::window("Digest")
            .open(&mut open)
            .default_size([560.0, 420.0])
            .show(ctx, |ui| {
                let chosen = self.selected.clone();
                let has_selection = self.docs[i].selection().is_some();
                let result = compute(
                    &mut self.docs[i],
                    &self.enzymes,
                    &chosen,
                    self.digest_win.selection_only && has_selection,
                );
                ui.horizontal(|ui| {
                    ui.label(
                        RichText::new(if result.enzymes.is_empty() {
                            "No enzymes selected".to_string()
                        } else {
                            result.enzymes.join(" + ")
                        })
                        .strong(),
                    );
                    if ui.button("Select enzymes…").clicked() {
                        self.enzyme_win.open = true;
                    }
                });
                ui.horizontal(|ui| {
                    ui.add_enabled(
                        has_selection,
                        egui::Checkbox::new(&mut self.digest_win.selection_only, "Selection only"),
                    );
                    ui.label("Ladder");
                    egui::ComboBox::from_id_salt("ladder")
                        .selected_text(LADDERS[self.digest_win.ladder].0)
                        .show_ui(ui, |ui| {
                            for (k, (name, _)) in LADDERS.iter().enumerate() {
                                ui.selectable_value(&mut self.digest_win.ladder, k, *name);
                            }
                        });
                });
                ui.separator();
                ui.horizontal_top(|ui| {
                    draw_gel(
                        ui,
                        &pal,
                        LADDERS[self.digest_win.ladder].1,
                        &result.fragments,
                    );
                    ui.vertical(|ui| {
                        ScrollArea::vertical().max_height(300.0).show(ui, |ui| {
                            egui::Grid::new("fragments")
                                .striped(true)
                                .spacing([12.0, 2.0])
                                .show(ui, |ui| {
                                    for h in ["Size", "Site 1", "Site 2", "Mass %"] {
                                        ui.label(RichText::new(h).strong());
                                    }
                                    ui.end_row();
                                    for f in &result.fragments {
                                        ui.label(RichText::new(commify(f.length)).monospace());
                                        ui.label(site_text(&f.left, f.start + result.origin));
                                        ui.label(site_text(&f.right, f.end + result.origin));
                                        ui.label(RichText::new(f.mass_pct.to_string()).monospace());
                                        ui.end_row();
                                    }
                                });
                        });
                        ui.label(
                            RichText::new(format!(
                                "{} fragments, {} bp total",
                                result.fragments.len(),
                                commify(result.total)
                            ))
                            .small(),
                        );
                        if ui.button("Copy table").clicked() {
                            let text = std::iter::once("Size\tSite 1\tSite 2\tMass %".to_string())
                                .chain(result.fragments.iter().map(|f| {
                                    format!(
                                        "{}\t{} {}\t{} {}\t{}",
                                        f.length,
                                        f.left,
                                        f.start + result.origin,
                                        f.right,
                                        f.end + result.origin,
                                        f.mass_pct
                                    )
                                }))
                                .collect::<Vec<_>>()
                                .join("\n");
                            ui.ctx().copy_text(text);
                        }
                    });
                });
            });
        self.digest_win.open = open;
    }
}

fn site_text(label: &str, pos: usize) -> String {
    if label.is_empty() || label == "uncut" {
        label.to_string()
    } else {
        format!("{label} {pos}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::doc::Doc;
    use crate::seq::enzymes::default_enzymes;

    fn plasmid() -> DocView {
        let mut seq = vec![b'T'; 600];
        seq[100..106].copy_from_slice(b"GAATTC");
        seq[420..426].copy_from_slice(b"AAGCTT");
        let mut d = Doc::new("p", seq);
        d.circular = true;
        DocView::new(d)
    }

    fn names(v: &[&str]) -> BTreeSet<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn digesting_with_two_enzymes_splits_the_circle() {
        let r = compute(
            &mut plasmid(),
            &default_enzymes(),
            &names(&["EcoRI", "HindIII"]),
            false,
        );
        // EcoRI cuts after base 101 and HindIII after base 421
        assert_eq!(
            r.fragments.iter().map(|f| f.length).collect::<Vec<_>>(),
            vec![320, 280]
        );
        assert_eq!(r.fragments.iter().map(|f| f.length).sum::<usize>(), 600);
        assert_eq!(r.enzymes, vec!["EcoRI", "HindIII"]);
    }

    #[test]
    fn no_enzymes_means_one_uncut_fragment() {
        let r = compute(&mut plasmid(), &default_enzymes(), &names(&[]), false);
        assert_eq!(r.fragments.len(), 1);
        assert_eq!(r.fragments[0].left, "uncut");
    }

    #[test]
    fn selection_only_digests_the_selected_piece_as_linear() {
        let mut dv = plasmid();
        dv.select(50, 500);
        let r = compute(
            &mut dv,
            &default_enzymes(),
            &names(&["EcoRI", "HindIII"]),
            true,
        );
        assert_eq!(r.total, 450);
        assert_eq!(r.fragments.iter().map(|f| f.length).sum::<usize>(), 450);
        assert_eq!(r.fragments.len(), 3);
        assert_eq!(r.origin, 50);
    }

    #[test]
    fn gel_positions_follow_the_log_scale() {
        let g = GelScale::new(100, 10_000);
        assert!((g.y(10_000) - 30.0).abs() < 1e-3);
        assert!((g.y(100) - 240.0).abs() < 1e-3);
        assert!(g.y(1000) > g.y(5000));
        assert!(LADDERS
            .iter()
            .all(|(_, l)| l.windows(2).all(|w| w[0] >= w[1])));
    }
}
