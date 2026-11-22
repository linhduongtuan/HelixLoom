use super::*;
use crate::seq::minimap::sim::{genome, read, Rng};
use crate::seq::minimap::{map_read, Options};

fn mapped(cigar: bool) -> (Index, Vec<u8>, Hit) {
    let g = genome(20_000, 31);
    let o = Options {
        cigar,
        ..Options::default()
    };
    let index = Index::build(
        vec![
            ("chr1".into(), genome(5000, 32)),
            ("chr2".into(), g.clone()),
        ],
        &o,
    )
    .unwrap();
    let mut r = Rng::new(2);
    let q = read(&g, 3000, 2500, true, 0.02, &mut r);
    let hit = map_read(&index, &q, &o, &|_| false).remove(0);
    (index, q, hit)
}

#[test]
fn a_line_has_the_twelve_columns_of_paf_and_the_tags_of_a_chain() {
    let (index, q, h) = mapped(false);
    let text = line("read1", q.len(), &index, &h, false);
    let f: Vec<&str> = text.split('\t').collect();
    assert_eq!(&f[..2], ["read1", q.len().to_string().as_str()]);
    assert_eq!(
        (
            f[2].parse::<usize>().unwrap(),
            f[3].parse::<usize>().unwrap()
        ),
        (h.qstart, h.qend)
    );
    assert_eq!((f[4], f[5], f[6]), ("-", "chr2", "20000"));
    assert_eq!(
        (
            f[7].parse::<usize>().unwrap(),
            f[8].parse::<usize>().unwrap()
        ),
        (h.tstart, h.tend)
    );
    // without the alignment: the bases the anchors cover, and the longer of the two spans
    assert_eq!(
        (
            f[9].parse::<usize>().unwrap(),
            f[10].parse::<usize>().unwrap()
        ),
        (h.covered, (h.qend - h.qstart).max(h.tend - h.tstart))
    );
    assert_eq!(f[11], h.mapq.to_string());
    assert_eq!(f[12], "tp:A:P");
    assert_eq!(f[13], format!("cm:i:{}", h.anchors.len()));
    assert_eq!(f[14], format!("s1:i:{}", h.score.round() as i64));
    assert!(f.iter().any(|t| t.starts_with("dv:f:")));
    assert!(f
        .iter()
        .all(|t| !t.starts_with("cg:") && !t.starts_with("NM:")));
    assert!(!text.contains('\n'));
}

#[test]
fn with_the_alignment_the_matches_and_the_block_come_from_it_and_the_cigar_and_nm_are_there() {
    let (index, q, h) = mapped(true);
    let text = line("r", q.len(), &index, &h, false);
    let f: Vec<&str> = text.split('\t').collect();
    let b = h.base.as_ref().expect("an alignment");
    assert_eq!(
        (
            f[9].parse::<usize>().unwrap(),
            f[10].parse::<usize>().unwrap()
        ),
        (b.matches, b.block)
    );
    assert!(f.contains(&format!("NM:i:{}", b.differences).as_str()));
    assert!(f.contains(&format!("cg:Z:{}", b.cigar).as_str()));
    // the matches are fewer than the block and most of it
    assert!(b.matches < b.block && b.matches as f64 > 0.85 * b.block as f64);
}

#[test]
fn a_secondary_says_so_a_primary_with_a_rival_gives_its_score_and_overlaps_have_no_quality() {
    let (index, q, mut h) = mapped(false);
    h.best_secondary = 123.4;
    let t = line("r", q.len(), &index, &h, false);
    assert!(
        t.contains("\ts2:i:123\t") || t.contains("\ts2:i:123"),
        "{t}"
    );
    h.primary = false;
    let t = line("r", q.len(), &index, &h, false);
    assert!(t.contains("tp:A:S") && !t.contains("s2:i:"), "{t}");
    // an overlap: no quality (255), no secondary score, every chain a primary
    h.primary = true;
    let t = line("r", q.len(), &index, &h, true);
    let f: Vec<&str> = t.split('\t').collect();
    assert_eq!(f[11], "255");
    assert!(!t.contains("s2:i:"));
    // the divergence is left out when it could not be told
    h.divergence = None;
    assert!(!line("r", q.len(), &index, &h, false).contains("dv:f:"));
}
