use super::*;
use crate::seq::minimap::sim::{genome, noisy, Rng};

fn anchors(from: (usize, usize), step: usize, n: usize) -> Vec<Anchor> {
    (0..n)
        .map(|i| Anchor {
            x: from.0 + i * step,
            y: from.1 + i * step,
        })
        .collect()
}

/// What a CIGAR consumes of the read and of the target.
fn spans(cigar: &str) -> (usize, usize) {
    let (mut q, mut t, mut n) = (0, 0, 0usize);
    for c in cigar.chars() {
        match c {
            '0'..='9' => n = n * 10 + c.to_digit(10).unwrap() as usize,
            'M' => (q, t, n) = (q + n, t + n, 0),
            'I' => (q, n) = (q + n, 0),
            'D' => (t, n) = (t + n, 0),
            other => panic!("{other}"),
        }
    }
    (q, t)
}

/// The matches and the columns that a CIGAR gives when it is walked along the two sequences.
fn walked(cigar: &str, q: &[u8], t: &[u8], (qs, ts): (usize, usize)) -> (usize, usize) {
    let (mut qi, mut ti, mut n, mut matches, mut columns) = (qs, ts, 0usize, 0, 0);
    for c in cigar.chars() {
        match c {
            '0'..='9' => n = n * 10 + c.to_digit(10).unwrap() as usize,
            'M' => {
                matches += (0..n).filter(|&i| q[qi + i] == t[ti + i]).count();
                (qi, ti, columns, n) = (qi + n, ti + n, columns + n, 0);
            }
            'I' => (qi, columns, n) = (qi + n, columns + n, 0),
            _ => (ti, columns, n) = (ti + n, columns + n, 0),
        }
    }
    (matches, columns)
}

#[test]
fn a_read_that_is_the_target_is_one_run_of_matches() {
    let t = genome(500, 1);
    let b = along(&anchors((100, 100), 10, 20), &t, &t, 15).unwrap();
    // from the first anchor to the end of the last: 100 to 295 plus 15
    assert_eq!(
        b,
        Base {
            cigar: "205M".into(),
            matches: 205,
            block: 205,
            differences: 0
        }
    );
}

#[test]
fn a_difference_between_two_anchors_is_aligned_and_counted() {
    let t = genome(300, 2);
    // a substitution
    let mut q = t.clone();
    q[60] = if q[60] == b'A' { b'C' } else { b'A' };
    let a = [Anchor { x: 20, y: 20 }, Anchor { x: 100, y: 100 }];
    let b = along(&a, &q, &t, 15).unwrap();
    assert_eq!(
        (b.cigar.as_str(), b.matches, b.block, b.differences),
        ("95M", 94, 95, 1)
    );
    // three extra bases in the read
    let mut q = t[..60].to_vec();
    q.extend_from_slice(b"TTT");
    q.extend_from_slice(&t[60..]);
    let a = [Anchor { x: 20, y: 20 }, Anchor { x: 100, y: 103 }];
    let b = along(&a, &q, &t, 15).unwrap();
    assert_eq!(spans(&b.cigar), (98, 95), "{}", b.cigar);
    assert!(b.cigar.contains("3I"), "{}", b.cigar);
    assert_eq!((b.matches, b.differences), (95, 3));
    // two bases of the target missing in the read
    let mut q = t[..60].to_vec();
    q.extend_from_slice(&t[62..]);
    let a = [Anchor { x: 20, y: 20 }, Anchor { x: 100, y: 98 }];
    let b = along(&a, &q, &t, 15).unwrap();
    assert_eq!(spans(&b.cigar), (93, 95), "{}", b.cigar);
    assert!(b.cigar.contains("2D"), "{}", b.cigar);
}

#[test]
fn anchors_on_one_line_that_overlap_make_one_block_and_one_that_does_not_fit_is_passed_over() {
    let t = genome(200, 3);
    // k-mers of 15 three bases apart: they overlap, the block is 12 + 15
    let b = along(&anchors((10, 10), 3, 5), &t, &t, 15).unwrap();
    assert_eq!(b.cigar, "27M");
    // an anchor that comes before the end of the block on another line has no place in it
    let a = [
        Anchor { x: 10, y: 10 },
        Anchor { x: 20, y: 28 },
        Anchor { x: 60, y: 60 },
    ];
    let b = along(&a, &t, &t, 15).unwrap();
    assert_eq!(spans(&b.cigar), (65, 65), "{}", b.cigar);
    assert_eq!(along(&[], &t, &t, 15), None);
}

#[test]
fn along_a_noisy_read_the_cigar_is_consistent_with_the_sequences() {
    let g = genome(4000, 4);
    let mut r = Rng::new(9);
    let q = noisy(&g[500..3500], 0.03, 0.03, 0.03, &mut r);
    // anchors where the k-mers are exact: find them the plain way along the diagonal of each position, resyncing as an index would
    let mut found = Vec::new();
    let mut shift: i64 = 0;
    let mut y = 0;
    while y + 15 <= q.len() {
        let mut hit = None;
        for d in -8i64..=8 {
            let x = (500 + y) as i64 + shift + d;
            if x >= 0
                && (x as usize) + 15 <= g.len()
                && g[x as usize..x as usize + 15] == q[y..y + 15]
            {
                hit = Some((x as usize, d));
                break;
            }
        }
        if let Some((x, d)) = hit {
            found.push(Anchor { x, y });
            shift += d;
        }
        y += 12;
    }
    assert!(found.len() > 40, "{}", found.len());
    let b = along(&found, &q, &g, 15).unwrap();
    let (qs, ts) = (found[0].y, found[0].x);
    let (qspan, tspan) = spans(&b.cigar);
    let last = found[found.len() - 1];
    assert_eq!((qspan, tspan), (last.y + 15 - qs, last.x + 15 - ts));
    assert_eq!(walked(&b.cigar, &q, &g, (qs, ts)), (b.matches, b.block));
    assert_eq!(b.block, b.matches + b.differences);
    assert!(
        b.matches as f64 / b.block as f64 > 0.8,
        "{}",
        b.matches as f64 / b.block as f64
    );
}

#[test]
fn a_gap_too_large_to_align_in_one_piece_is_still_written_and_consumes_both() {
    let t = genome(16_000, 5);
    let q = genome(16_000, 6);
    // 7,000 against 7,000 is 49 million cells: past the limit of one table
    let a = [Anchor { x: 100, y: 100 }, Anchor { x: 7200, y: 7200 }];
    let b = along(&a, &q, &t, 15).unwrap();
    assert_eq!(spans(&b.cigar), (7115, 7115));
    // and an uneven one: 8,000 of the read against 7,000 of the target
    let a = [Anchor { x: 100, y: 100 }, Anchor { x: 7200, y: 8200 }];
    let b = along(&a, &q, &t, 15).unwrap();
    assert_eq!(spans(&b.cigar), (8115, 7115));
    assert!(b.cigar.contains('I'));
}

#[test]
fn an_anchor_that_follows_the_last_one_on_the_read_with_a_jump_on_the_target_is_a_deletion() {
    let t = genome(100, 8);
    // the first k-mer ends at 15 on both, the next begins at 15 on the read and at 20 on the target
    let mut q = t[..15].to_vec();
    q.extend_from_slice(&t[20..35]);
    let a = [Anchor { x: 0, y: 0 }, Anchor { x: 20, y: 15 }];
    let b = along(&a, &q, &t, 15).unwrap();
    assert_eq!(
        (b.cigar.as_str(), b.matches, b.block, b.differences),
        ("15M5D15M", 30, 35, 5)
    );
    // and the same with the target joined and the read jumping: an insertion
    let mut q = t[..15].to_vec();
    q.extend_from_slice(b"TTTTT");
    q.extend_from_slice(&t[15..30]);
    let a = [Anchor { x: 0, y: 0 }, Anchor { x: 15, y: 20 }];
    let b = along(&a, &q, &t, 15).unwrap();
    assert_eq!((b.cigar.as_str(), b.differences), ("15M5I15M", 5));
}
