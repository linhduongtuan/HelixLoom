use super::*;
use crate::seq::doc::Doc;

/// 140 cases made by running the original's own PCR code (`search_with_mm3`, `pcr_dialog_precompute`,
/// `pcr_dialog_do_pcr`, `Tm`, cut out of AppMain.tcl unchanged) on random templates and primers.
/// Tab separated, a `T` line for each case (template, circular, min match, max mismatches, min 3' match)
/// followed by `S` lines (a primer and every site the scan lists) and `P` lines (a forward and a reverse
/// site, and what `do_pcr` made of them).
const GOLDEN: &str = include_str!("golden.txt");

/// 240 cases made by running the original's own `get_tags_get_coords` (which `do_pcr` reaches through
/// `get_tags2` to copy the template's features into the product) on random features: the template's
/// length, the stretch the product copies `[first, last)`, where it starts in the product, the
/// feature's spans `a-b,a-b`, and the coordinates the original gives (empty when nothing is copied).
const COORDS: &str = include_str!("coords.txt");

const TEMPLATE_NAME: &str = "pDemo";

fn primer(name: &str, seq: &str) -> PcrPrimer {
    PcrPrimer {
        name: name.to_string(),
        seq: seq.as_bytes().to_vec(),
    }
}

fn site(forward: bool, location: &str, matched: &str) -> Site {
    Site {
        forward,
        location: location.parse().expect("a location"),
        matched: matched.parse().expect("a match length"),
        mismatches: 0,
    }
}

/// The sites as the dialog's table lists them, `dir:location:matchlength:mismatches`, put in the order
/// Tcl's `lsort` gives the same rows written with spaces.
fn sites_as_listed(sites: &[Site]) -> String {
    let mut rows: Vec<[String; 4]> = sites
        .iter()
        .map(|s| {
            [
                if s.forward { ">>>" } else { "<<<" }.to_string(),
                s.location.to_string(),
                s.matched.to_string(),
                s.mismatches.to_string(),
            ]
        })
        .collect();
    rows.sort_by_key(|r| r.join(" "));
    let rows: Vec<String> = rows.iter().map(|r| r.join(":")).collect();
    rows.join(" ")
}

fn hash(seq: &[u8]) -> u64 {
    seq.iter()
        .fold(0, |h, &c| (h * 131 + u64::from(c)) % 1_000_000_007)
}

fn hex(color: [u8; 3]) -> String {
    format!("#{:02x}{:02x}{:02x}", color[0], color[1], color[2])
}

/// The features as the golden file writes them: `start-end,colour,name,reverse,qualifiers`.
fn features_as_listed(features: &[Feature]) -> String {
    let rows: Vec<String> = features
        .iter()
        .map(|f| {
            let (start, end) = f.spans[0];
            let color = if f.reverse { f.rev_color } else { f.fwd_color };
            let (key, value) = &f.qualifiers[0];
            format!(
                "{start}-{end},{},{},{},{key} {{{value}}}",
                hex(color),
                f.name,
                u8::from(f.reverse)
            )
        })
        .collect();
    rows.join(";")
}

/// What the golden file holds for `ok`: the fields after the messages.
fn product_as_listed(p: &Product) -> String {
    let text = String::from_utf8_lossy(&p.seq);
    let n = p.seq.len();
    [
        n.to_string(),
        text[..n.min(30)].to_string(),
        text[n.saturating_sub(30)..].to_string(),
        hash(&p.seq).to_string(),
        p.comment.join("\\n"),
        p.name.clone(),
        features_as_listed(&p.features),
    ]
    .join("|")
}

#[test]
fn every_scan_finds_the_sites_the_originals_own_code_finds() {
    let (mut template, mut circular, mut scan) = (Vec::new(), false, PcrScan::default());
    let mut scans = 0;
    for line in GOLDEN.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "T" => {
                template = f[1].as_bytes().to_vec();
                circular = f[2] == "1";
                scan = PcrScan {
                    min_match: f[3].parse().unwrap(),
                    max_mismatch: f[4].parse().unwrap(),
                    min_3_match: f[5].parse().unwrap(),
                };
            }
            "S" => {
                let found = find_sites(&template, circular, f[1].as_bytes(), &scan);
                assert_eq!(
                    sites_as_listed(&found),
                    f[2],
                    "primer {} on {} template, {scan:?}",
                    f[1],
                    if circular { "a circular" } else { "a linear" }
                );
                scans += 1;
            }
            _ => {}
        }
    }
    assert_eq!(scans, 280);
}

#[test]
fn every_product_is_what_the_originals_own_code_makes() {
    let (mut template, mut circular) = (Vec::new(), false);
    let (mut made, mut refused) = (0, 0);
    for line in GOLDEN.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "T" => {
                template = f[1].as_bytes().to_vec();
                circular = f[2] == "1";
            }
            "P" => {
                let got = make_product(
                    &template,
                    circular,
                    TEMPLATE_NAME,
                    &primer("F1", f[1]),
                    &site(true, f[2], f[3]),
                    &primer("R1", f[4]),
                    &site(false, f[5], f[6]),
                );
                let what = format!("{} {}…, circular {circular}: {line}", f[2], f[5]);
                if let Some(message) = f[7].strip_prefix("rejected:") {
                    assert_eq!(got.unwrap_err(), message, "{what}");
                    refused += 1;
                    continue;
                }
                let (messages, rest) = f[7]
                    .strip_prefix("ok|")
                    .and_then(|r| r.split_once('|'))
                    .expect("an ok result");
                let p = got.unwrap_or_else(|e| panic!("refused ({e}): {what}"));
                // the first message holds the mismatch warnings, one to a line; each overlap warning is
                // another, and the messages are joined with a space
                let expected: Vec<String> = messages
                    .replace("\\n", "\n")
                    .replace(" Warning: ", "\nWarning: ")
                    .lines()
                    .filter(|m| !m.is_empty())
                    .map(str::to_string)
                    .collect();
                assert_eq!(p.warnings, expected, "{what}");
                assert_eq!(product_as_listed(&p), rest, "{what}");
                made += 1;
            }
            _ => {}
        }
    }
    assert_eq!((made, refused), (65, 42));
}

#[test]
fn the_template_features_a_product_copies_are_those_the_originals_own_code_copies() {
    let mut copied = 0;
    for line in COORDS.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let number = |text: &str| text.parse::<usize>().expect("a number");
        let (len, first, last, offset) = (number(f[0]), number(f[1]), number(f[2]), number(f[3]));
        let spans: Vec<(usize, usize)> = f[4]
            .split(',')
            .map(|s| {
                let (a, b) = s.split_once('-').expect("a span");
                (number(a), number(b))
            })
            .collect();
        let mut feature = Feature::new("f", "misc_feature", spans[0].0, spans[0].1, false);
        feature.spans = spans;
        let mut template = Doc::new("t", vec![b'A'; len]);
        template.features.push(feature);
        let mut got: Vec<usize> = template
            .features_in(first, last, false)
            .into_iter()
            .flat_map(|f| f.shifted(offset).spans)
            .flat_map(|(a, b)| [a, b])
            .collect();
        got.sort_unstable();
        let listed: Vec<usize> = f
            .get(5)
            .map_or(Vec::new(), |c| c.split_whitespace().map(number).collect());
        // a feature that only touches the stretch is listed as an empty range at its start, which tags
        // no base, so nothing is drawn
        let expected: Vec<usize> = listed
            .chunks(2)
            .filter(|pair| pair[0] != pair[1])
            .flatten()
            .copied()
            .collect();
        assert_eq!(got, expected, "{line}");
        copied += usize::from(!expected.is_empty());
    }
    assert!(copied > 80, "only {copied} features were copied");
}

fn synthetic() -> Vec<u8> {
    let mut seed = 99u64;
    (0..400)
        .map(|_| {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[((seed >> 33) % 4) as usize]
        })
        .collect()
}

fn products(
    template: &[u8],
    circular: bool,
    fwd: &[u8],
    rev: &[u8],
    scan: &PcrScan,
) -> Vec<Product> {
    let name = |n: &str, s: &[u8]| PcrPrimer {
        name: n.to_string(),
        seq: s.to_vec(),
    };
    pcr_products(
        template,
        circular,
        TEMPLATE_NAME,
        &name("Forward", fwd),
        &name("Reverse", rev),
        scan,
    )
}

fn scan(min_match: usize) -> PcrScan {
    PcrScan {
        min_match,
        ..PcrScan::default()
    }
}

#[test]
fn the_scan_defaults_are_the_original_dialogs() {
    let d = PcrScan::default();
    assert_eq!((d.min_match, d.max_mismatch, d.min_3_match), (17, 0, 4));
    assert_eq!(PCR_WARN_BELOW, 15);
}

#[test]
fn products_carry_both_primers_with_their_tails_and_a_comment() {
    let mut t = vec![b'T'; 400];
    let f_site = b"GATTACAGGCATCCAGTCAGTC";
    let r_site = revcom(b"CCTAGGTCAAGTCCGATCGTAA");
    t[100..122].copy_from_slice(f_site);
    t[300..322].copy_from_slice(&r_site);
    let mut fwd = b"GGATCC".to_vec();
    fwd.extend_from_slice(f_site);
    let mut rev = b"CTCGAG".to_vec();
    rev.extend_from_slice(b"CCTAGGTCAAGTCCGATCGTAA");
    let found = products(&t, false, &fwd, &rev, &scan(15));
    assert_eq!(found.len(), 1);
    let p = &found[0];
    assert_eq!(p.length(), 6 + (322 - 100) + 6);
    assert!(p.seq.starts_with(b"GGATCCGATTACAGG"));
    assert!(p.seq.ends_with(&revcom(&rev)));
    assert_eq!((p.fwd, p.rev), ((100, 122), (300, 322)));
    assert!(p.warnings.is_empty());
    assert_eq!(p.name, "PCR pDemo Forward Reverse");
    assert_eq!(p.comment[..2], ["PCR:", "template: pDemo"]);
    assert_eq!(p.comment.last().unwrap(), "Product length: 234");
    // the two primer features: the whole forward primer, then the whole reverse primer
    assert_eq!(p.features[0].spans, vec![(0, 28)]);
    assert_eq!(p.features[1].spans, vec![(206, 234)]);
    assert!(!p.features[0].reverse && p.features[1].reverse);
    assert_eq!(p.template_span, Some((100, 322)));
    assert_eq!(p.template_offset, 6);
}

#[test]
fn mismatches_are_tolerated_inside_a_primer_but_never_at_its_3_prime_end() {
    let f_site = b"GATTACAGGCATCCAGTCAGTC";
    let r_site = b"CCTAGGTCAAGTCCGATCGTAA";
    let mut t = vec![b'T'; 400];
    t[100..122].copy_from_slice(f_site);
    t[300..322].copy_from_slice(&revcom(r_site));
    let mut mutated = f_site.to_vec();
    mutated[8] = b'T'; // one mismatch, 13 bases from the 3' end
    assert!(products(&t, false, &mutated, r_site, &scan(20)).is_empty());
    let loose = PcrScan {
        min_match: 20,
        max_mismatch: 1,
        min_3_match: 4,
    };
    let found = products(&t, false, &mutated, r_site, &loose);
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].fwd, (100, 122));
    assert_eq!(
        found[0].warnings,
        vec!["Warning: First forward primer mismatch 13 bases from 3' end.".to_string()]
    );
    // a mismatch among the last four bases is never accepted
    let mut bad_end = f_site.to_vec();
    *bad_end.last_mut().unwrap() = b'A';
    assert!(products(&t, false, &bad_end, r_site, &loose).is_empty());
}

#[test]
fn a_short_exact_3_prime_match_warns_like_the_original() {
    let f_site = b"GATTACAGGCATCCAGTCAGTC";
    let r_site = b"CCTAGGTCAAGTCCGATCGTAA";
    let mut t = vec![b'T'; 400];
    t[100..122].copy_from_slice(f_site);
    t[300..322].copy_from_slice(&revcom(r_site));
    let mut tailed = b"GGATCA".to_vec(); // its last base differs from the template's C
    tailed.extend_from_slice(&f_site[10..]); // so only 12 bases match the template
    let found = products(&t, false, &tailed, r_site, &scan(10));
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].warnings.len(), 1);
    assert!(found[0].warnings[0].contains("mismatch 12 bases"));
    // the comment gives the Tm of the matching 12 bases and of the whole primer
    assert!(found[0]
        .comment
        .iter()
        .any(|l| l.contains("Tm1:") && l.contains("Tm2:")));
}

#[test]
fn products_can_span_the_origin_of_a_circle() {
    let mut t = vec![b'T'; 200];
    let f_site = b"GATTACAGGCATCCAGTCAGTC";
    let r_site = b"CCTAGGTCAAGTCCGATCGTAA";
    t[180..200].copy_from_slice(&f_site[..20]);
    t[..2].copy_from_slice(&f_site[20..]);
    t[40..62].copy_from_slice(&revcom(r_site));
    let found = products(&t, true, f_site, r_site, &scan(15));
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].length(), 62 + 20);
    assert_eq!(found[0].fwd, (180, 2));
    assert_eq!(found[0].template_span, None);
    assert!(products(&t, false, f_site, r_site, &scan(15)).is_empty());
}

#[test]
fn primers_that_touch_or_overlap_on_a_line_do_not_face_each_other() {
    // the original wants the forward 3' end strictly left of the reverse start: touching is refused
    let t = synthetic();
    let (f, r) = (primer("F", "ACGT"), primer("R", "ACGT"));
    for (fi, ri) in [(100, 100), (120, 100)] {
        let refused = make_product(
            &t,
            false,
            "t",
            &f,
            &site(true, &fi.to_string(), "4"),
            &r,
            &site(false, &ri.to_string(), "4"),
        );
        assert_eq!(refused.unwrap_err(), NOT_FACING);
    }
    // on a circle the same pair simply goes round the other way
    let round = make_product(
        &t,
        true,
        "t",
        &f,
        &site(true, "120", "4"),
        &r,
        &site(false, "100", "4"),
    );
    assert!(round.is_ok());
}

#[test]
fn primers_that_overlap_on_a_circle_warn_in_the_originals_words() {
    let t = synthetic();
    let f = primer("F", &String::from_utf8_lossy(&t[100..124]));
    let r = primer("R", &String::from_utf8_lossy(&revcom(&t[110..134])));
    let found = products(&t, true, &f.seq, &r.seq, &scan(17));
    assert_eq!(found.len(), 1);
    // the reverse primer sits inside the forward one's reach: the product is the whole circle and more
    assert!(found[0].length() > t.len());
    assert!(found[0]
        .warnings
        .iter()
        .any(|w| w.starts_with("Warning: Reverse primer amplifies part of forward primer")));
    assert!(found[0]
        .warnings
        .iter()
        .any(|w| w.starts_with("Warning: Forwrd primer amplifies part of reverse primer")));
}

#[test]
fn short_or_missing_matches_produce_nothing() {
    let t = synthetic();
    assert!(products(
        &t,
        false,
        b"ACGTACGTACGTAAAAAAAA",
        b"CCCCCCCCCCCCCCCCCCCCCCGG",
        &scan(15)
    )
    .is_empty());
    assert!(find_sites(&t, false, b"ACGT", &scan(17)).is_empty());
    assert!(find_sites(&[], true, b"ACGTACGTACGTACGTACGT", &scan(17)).is_empty());
}
