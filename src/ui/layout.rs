//! Row layout of the sequence editor: how many bases fit per row, which lane each feature
//! sits in, and the height of every row.

use crate::seq::doc::Doc;

pub const ROW_GAP: f32 = 10.0;
const LANE_PAD: usize = 1;

pub struct Layout {
    pub per_row: usize,
    pub n: usize,
    pub rows: usize,
    /// Lane of each feature (index aligned with `doc.features`).
    pub lanes: Vec<usize>,
    /// Number of feature lanes shown in each row.
    pub row_lanes: Vec<usize>,
    /// Top of each row; one extra entry holds the total height.
    pub row_y: Vec<f32>,
    pub line_h: f32,
    pub lane_h: f32,
    pub second_strand: bool,
}

impl Layout {
    pub fn new(
        doc: &Doc,
        per_row: usize,
        line_h: f32,
        lane_h: f32,
        second_strand: bool,
        trans_rows: usize,
    ) -> Layout {
        let per_row = per_row.max(1);
        let n = doc.len();
        let rows = n.div_ceil(per_row).max(1);
        let lanes = assign_lanes(doc);
        let mut row_lanes = vec![0usize; rows];
        for (f, &lane) in doc.features.iter().zip(&lanes) {
            for &(s, e) in &f.spans {
                if e <= s {
                    continue;
                }
                for slot in &mut row_lanes[s / per_row..=((e - 1) / per_row).min(rows - 1)] {
                    *slot = (*slot).max(lane + 1);
                }
            }
        }
        let mut row_y = Vec::with_capacity(rows + 1);
        let mut y = 0.0;
        for &l in &row_lanes {
            row_y.push(y);
            let strands = if second_strand { 2.0 } else { 1.0 };
            y += line_h * strands + l as f32 * lane_h + trans_rows as f32 * line_h * 0.9 + ROW_GAP;
        }
        row_y.push(y);
        Layout {
            per_row,
            n,
            rows,
            lanes,
            row_lanes,
            row_y,
            line_h,
            lane_h,
            second_strand,
        }
    }

    pub fn total_height(&self) -> f32 {
        *self.row_y.last().unwrap_or(&0.0)
    }

    /// Base range `[start, end)` shown in row `r`.
    pub fn row_range(&self, r: usize) -> (usize, usize) {
        let s = r * self.per_row;
        (s.min(self.n), (s + self.per_row).min(self.n))
    }

    pub fn row_of_base(&self, base: usize) -> usize {
        (base / self.per_row).min(self.rows - 1)
    }

    /// Row containing vertical offset `y` (clamped).
    pub fn row_at_y(&self, y: f32) -> usize {
        let i = self.row_y.partition_point(|&top| top <= y);
        i.saturating_sub(1).min(self.rows - 1)
    }

    pub fn lanes_y(&self, r: usize) -> f32 {
        self.row_y[r] + self.line_h * if self.second_strand { 2.0 } else { 1.0 }
    }

    pub fn trans_y(&self, r: usize) -> f32 {
        self.lanes_y(r) + self.row_lanes[r] as f32 * self.lane_h
    }
}

/// Greedy interval colouring: the first lane whose existing features don't touch this one.
fn assign_lanes(doc: &Doc) -> Vec<usize> {
    let mut order: Vec<usize> = (0..doc.features.len()).collect();
    order.sort_by_key(|&i| {
        (
            doc.features[i].start(),
            std::cmp::Reverse(doc.features[i].len()),
        )
    });
    let mut lanes: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut out = vec![0usize; doc.features.len()];
    for i in order {
        let spans = &doc.features[i].spans;
        let free = |taken: &Vec<(usize, usize)>| {
            spans.iter().all(|&(s, e)| {
                taken
                    .iter()
                    .all(|&(ts, te)| e + LANE_PAD <= ts || te + LANE_PAD <= s)
            })
        };
        let lane = match lanes.iter().position(free) {
            Some(l) => l,
            None => {
                lanes.push(Vec::new());
                lanes.len() - 1
            }
        };
        lanes[lane].extend(spans.iter().copied());
        out[i] = lane;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::features::Feature;

    fn doc_with(features: Vec<Feature>) -> Doc {
        let mut d = Doc::new("t", vec![b'A'; 250]);
        d.features = features;
        d
    }

    #[test]
    fn overlapping_features_get_different_lanes() {
        let d = doc_with(vec![
            Feature::new("a", "gene", 10, 60, false),
            Feature::new("b", "gene", 40, 90, false),
            Feature::new("c", "gene", 100, 120, false),
        ]);
        let l = Layout::new(&d, 100, 18.0, 17.0, true, 0);
        assert_eq!(l.lanes, vec![0, 1, 0]);
        assert_eq!(l.rows, 3);
        assert_eq!(l.row_lanes, vec![2, 1, 0]);
    }

    #[test]
    fn rows_grow_with_their_lanes_and_translation() {
        let d = doc_with(vec![Feature::new("a", "gene", 10, 60, false)]);
        let plain = Layout::new(&d, 100, 18.0, 17.0, false, 0);
        let rich = Layout::new(&d, 100, 18.0, 17.0, true, 3);
        assert!(rich.total_height() > plain.total_height());
        assert_eq!(plain.row_range(2), (200, 250));
        assert_eq!(plain.row_at_y(plain.row_y[1] + 1.0), 1);
        assert_eq!(plain.row_at_y(1e9), plain.rows - 1);
        assert_eq!(plain.row_of_base(99), 0);
        assert_eq!(plain.row_of_base(100), 1);
    }

    #[test]
    fn wrapping_features_reserve_lanes_on_both_ends() {
        let mut f = Feature::new("w", "gene", 240, 250, false);
        f.spans.push((0, 20));
        let d = doc_with(vec![f, Feature::new("x", "gene", 5, 15, false)]);
        let l = Layout::new(&d, 100, 18.0, 17.0, true, 0);
        assert_ne!(l.lanes[0], l.lanes[1]);
        assert_eq!(l.row_lanes[2], 1);
    }

    #[test]
    fn an_empty_sequence_still_has_one_row() {
        let l = Layout::new(&Doc::new("e", Vec::new()), 60, 18.0, 17.0, true, 0);
        assert_eq!(l.rows, 1);
        assert_eq!(l.row_range(0), (0, 0));
    }
}
