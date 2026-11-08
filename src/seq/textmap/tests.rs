use super::features::sequence_tags;
use super::genes::{translate, uppercase_cdna, STANDARD_CODE};
use super::out::Pos;
use super::*;
use crate::seq::enzymes::{all_sites, default_enzymes};
use crate::seq::features::Rgb;

/// 70 cases made by running the original's `do_text_map` (with `get_tags_exact`, `get_cds`, `copy_tags`, `linear_enz_sort`,
/// `translate`, `uppercase_cdna` and `findenzymes` behind it) against real Tk text widgets: a random sequence with features,
/// the part looked at, the enzymes, the dialog's settings and the order of its frames, and what the original wrote: the text of
/// its window and the ranges and colour of every tag in it (tests/tcl/text_map_driver.tcl).
const GOLDEN: &str = include_str!("golden.txt");

/// A tag of the page: the colour behind it and its ranges.
type Tag = (Option<Rgb>, Vec<(Pos, Pos)>);

enum Expect {
    Error(String),
    Page {
        text: String,
        tags: BTreeMap<String, Tag>,
    },
}

struct Case {
    number: usize,
    circular: bool,
    seq: Vec<u8>,
    title: String,
    region: Region,
    selection: Option<(usize, usize)>,
    row_width: usize,
    methylated: bool,
    enzymes: Vec<String>,
    settings: Settings,
    features: Vec<MapFeature>,
    /// The `fn` range of each feature the original made one for, as it listed it.
    joins: Vec<(usize, Vec<usize>)>,
    expect: Expect,
}

/// The frame the original keeps under a name.
fn frame_named(name: &str) -> Frame {
    let names = [
        "index",
        "enzymes",
        "dna",
        "translation",
        "sec_strand",
        "genes",
        "features",
        "space",
    ];
    let at = names.iter().position(|n| *n == name).expect("a frame");
    Frame::ALL[at]
}

fn color(text: &str) -> Rgb {
    let n = u32::from_str_radix(text.trim_start_matches('#'), 16).expect("a colour");
    [(n >> 16) as u8, (n >> 8) as u8, n as u8]
}

fn numbers(words: &[&str]) -> Vec<usize> {
    words.iter().map(|w| w.parse().expect("a number")).collect()
}

fn pairs(list: &[usize]) -> Vec<(usize, usize)> {
    list.chunks(2).map(|p| (p[0], p[1])).collect()
}

fn setting(settings: &mut Settings, key: &str, value: &str) {
    let on = value == "1";
    let mark = |v: &str| {
        if v == "_" {
            ' '
        } else {
            v.chars().next().expect("a mark")
        }
    };
    let count = || value.parse::<usize>().expect("a number");
    match key {
        "width" => settings.width = count(),
        "space" => settings.space = count(),
        "show_space" => settings.show_space = on,
        "numbers" => {
            settings.numbers = LineNumbers::ALL
                .into_iter()
                .find(|n| n.label() == value)
                .expect("a numbers choice")
        }
        "copy" => settings.copy_highlight = on,
        "second" => settings.second_strand = on,
        "show_trans" => settings.show_translation = on,
        "trans" => {
            let (frames, letters) = value.split_once('_').expect("frames_letters");
            settings.reading = Reading::new(frames.parse().expect("frames"), letters == "3");
        }
        "show_enz" => settings.show_enzymes = on,
        "genes" => settings.genes = on,
        "show_index" => settings.show_index = on,
        "char1" => settings.index_char1 = mark(value),
        "char2" => settings.index_char2 = mark(value),
        "sp1" => settings.index_spacing1 = count(),
        "sp2" => settings.index_spacing2 = count(),
        "nums" => settings.index_numbers = on,
        "show_features" => settings.show_features = on,
        "hidden" => settings.show_hidden_features = on,
        "crop" => settings.crop_features = on,
        other => panic!("a setting that is not known: {other}"),
    }
}

/// `feature i dir hidden fwd rev type | ranges… | fn… | name`
fn feature_of(line: &str) -> (MapFeature, Option<(usize, Vec<usize>)>) {
    let parts: Vec<&str> = line.splitn(4, " | ").collect();
    let words: Vec<&str> = parts[0].split(' ').collect();
    let number: usize = words[1].parse().expect("a number");
    let numbers_of = |text: &str| numbers(&text.split_whitespace().collect::<Vec<_>>());
    let ranges = numbers_of(parts[1]);
    let join = numbers_of(parts[2]);
    let feature = MapFeature {
        name: parts[3].to_string(),
        kind: words[6].to_string(),
        spans: pairs(&ranges),
        reverse: words[2] == "1",
        hidden: words[3] == "1",
        fwd: color(words[4]),
        rev: color(words[5]),
    };
    let joined = !join.is_empty();
    (feature, joined.then_some((number, join)))
}

fn position(word: &str) -> Pos {
    let (line, col) = word.split_once('.').expect("line.column");
    (
        line.parse().expect("a line"),
        col.parse().expect("a column"),
    )
}

fn tag_of(words: &[&str]) -> (String, Tag) {
    let background = (words[2] != "-").then(|| color(words[2]));
    let places: Vec<Pos> = words[3..].iter().map(|w| position(w)).collect();
    let ranges = places.chunks(2).map(|p| (p[0], p[1])).collect();
    (words[1].to_string(), (background, ranges))
}

fn parse(file: &str) -> Vec<Case> {
    let mut cases = Vec::new();
    let mut lines = file.lines();
    while let Some(head) = lines.next() {
        let Some(number) = head.strip_prefix("CASE textmap ") else {
            continue;
        };
        let mut case = Case {
            number: number.parse().expect("a case number"),
            circular: false,
            seq: Vec::new(),
            title: String::new(),
            region: Region::All,
            selection: None,
            row_width: 0,
            methylated: false,
            enzymes: Vec::new(),
            settings: Settings::default(),
            features: Vec::new(),
            joins: Vec::new(),
            expect: Expect::Error(String::new()),
        };
        let (mut text, mut tags) = (Vec::new(), BTreeMap::new());
        let mut error = None;
        let mut in_text = false;
        for line in lines.by_ref() {
            if line == "END" {
                break;
            }
            if in_text && line.starts_with('|') {
                text.push(line[1..].to_string());
                continue;
            }
            let words: Vec<&str> = line.split(' ').collect();
            match words[0] {
                "circular" => case.circular = words[1] == "1",
                "seq" => case.seq = words[1].as_bytes().to_vec(),
                "title" => case.title = line["title ".len()..].to_string(),
                "region" if words[1] == "all" => case.region = Region::All,
                "region" => {
                    let (a, b) = (words[1].parse().expect("a"), words[2].parse().expect("b"));
                    case.region = Region::Part(a, b);
                }
                "rowwidth" => case.row_width = words[1].parse().expect("a width"),
                "selection" if words[1] == "-" => case.selection = None,
                "selection" => {
                    case.selection =
                        Some((words[1].parse().expect("a"), words[2].parse().expect("b")))
                }
                "methylated" => case.methylated = words[1] == "1",
                "enz" => case.enzymes = words[1..].iter().map(|w| w.to_string()).collect(),
                "settings" => {
                    for pair in &words[1..] {
                        let (key, value) = pair.split_once('=').expect("key=value");
                        setting(&mut case.settings, key, value);
                    }
                }
                "order" => {
                    case.settings.order = words[1..].iter().map(|w| frame_named(w)).collect()
                }
                "feature" => {
                    let (feature, join) = feature_of(line);
                    case.features.push(feature);
                    case.joins.extend(join);
                }
                "ERROR" => error = Some(line["ERROR ".len()..].to_string()),
                "out" => in_text = true,
                "tags" => in_text = false,
                "tag" => {
                    let (name, value) = tag_of(&words);
                    tags.insert(name, value);
                }
                _ => {}
            }
        }
        case.expect = match error {
            Some(message) => Expect::Error(message),
            None => Expect::Page {
                text: text.join("\n"),
                tags,
            },
        };
        cases.push(case);
    }
    cases
}

/// The enzymes of a case with their sites as the build finds them, numbered from 1, in the order the case lists them.
fn enzymes_of(case: &Case) -> Vec<(String, Vec<usize>)> {
    let all = default_enzymes();
    let sites = all_sites(&all, &case.seq, case.circular, case.methylated);
    let found: BTreeMap<String, Vec<usize>> = all
        .iter()
        .zip(sites)
        .map(|(e, s)| {
            (
                e.name.clone(),
                s.iter().map(|site| site.start + 1).collect(),
            )
        })
        .collect();
    case.enzymes
        .iter()
        .map(|name| (name.clone(), found[name].clone()))
        .collect()
}

fn run(case: &Case) -> Result<Option<TextMap>, Error> {
    let sites = enzymes_of(case);
    let enzymes: Vec<EnzymeSites> = sites
        .iter()
        .map(|(name, sites)| EnzymeSites { name, sites })
        .collect();
    text_map(&Request {
        seq: &case.seq,
        circular: case.circular,
        title: &case.title,
        date: "DATE",
        region: case.region,
        selection: case.selection,
        row_width: case.row_width,
        enzymes: &enzymes,
        features: &case.features,
        code: STANDARD_CODE,
        settings: &case.settings,
    })
}

fn check(cases: &[Case]) {
    for case in cases {
        let made = run(case);
        match (&case.expect, made) {
            (Expect::Error(want), Err(got)) => assert_eq!(&got.0, want, "case {}", case.number),
            (Expect::Page { text, tags }, Ok(Some(page))) => {
                assert_eq!(&page.text, text, "the text of case {}", case.number);
                let got: BTreeMap<String, Tag> = page
                    .tags
                    .iter()
                    .filter(|t| !t.ranges.is_empty())
                    .map(|t| (t.name.clone(), (t.background, t.ranges.clone())))
                    .collect();
                assert_eq!(&got, tags, "the tags of case {}", case.number);
            }
            (_, Ok(_)) => panic!("case {} should have gone the other way", case.number),
            (_, Err(e)) => panic!("case {}: {e}", case.number),
        }
    }
}

#[test]
fn the_page_is_the_one_the_originals_code_makes_in_every_case() {
    let cases = parse(GOLDEN);
    assert_eq!(cases.len(), 70);
    check(&cases);
}

#[test]
fn the_join_tags_are_the_ones_the_original_puts_on_a_sequence() {
    for case in parse(GOLDEN) {
        let tags = sequence_tags(&case.features, case.seq.len(), case.circular);
        let made: Vec<(usize, Vec<(usize, usize)>)> = tags
            .iter()
            .filter(|t| t.is_fn)
            .map(|t| (t.feature, t.ranges.clone()))
            .collect();
        // Tk lists the ranges of a tag in order, and the original adds `first..len` and `0..last` for one through the origin
        let want: Vec<(usize, Vec<(usize, usize)>)> = case
            .joins
            .iter()
            .map(|(number, list)| {
                let mut ranges = pairs(list);
                ranges.sort_unstable();
                (*number, ranges)
            })
            .collect();
        assert_eq!(made, want, "case {}", case.number);
    }
}

#[test]
fn translating_a_row_reads_the_codons_as_the_original_does() {
    let t = |dna: &str, three, forward, spacing| {
        translate(
            dna.as_bytes(),
            STANDARD_CODE.as_bytes(),
            three,
            forward,
            spacing,
        )
    };
    assert_eq!(t("ATGGCCTAA", false, true, 2), "M  A  *");
    assert_eq!(t("ATGGCCTAA", true, true, 0), "MetAlaEnd");
    // a codon with a letter that is not a base is X, a trailing part of a codon is left out, and u is t
    assert_eq!(t("ATGNCCAUGA", false, true, 2), "M  X  M");
    assert_eq!(t("ATGNCCAUGA", true, true, 0), "MetXxxMet");
    // the other strand reads each codon as its complement (the order does not turn round)
    assert_eq!(t("ATGGCCTAA", false, false, 2), "H  G  L");
    assert_eq!(t("", false, true, 2), "");
}

#[test]
fn uppercase_cdna_takes_the_small_letters_out_and_says_where() {
    // a lead of small letters, an intron between capitals (gt...ag), a tail of small letters after a stop
    let c = uppercase_cdna(b"ccATGAAAgtaaggagCCCTGAtttt");
    assert_eq!(String::from_utf8_lossy(&c.dna), "ATGAAACCCTGA");
    assert_eq!(c.introns, vec![(0, 2), (6, 8), (12, 4)]);
}

#[test]
fn a_selection_through_the_origin_is_the_error_the_original_stops_with() {
    let settings = Settings::default();
    let made = text_map(&Request {
        seq: b"ACGTACGT",
        circular: true,
        title: "x",
        date: "DATE",
        region: Region::ThroughOrigin,
        selection: None,
        row_width: 0,
        enzymes: &[],
        features: &[],
        code: STANDARD_CODE,
        settings: &settings,
    });
    assert_eq!(
        made.err().map(|e| e.0),
        Some("can't read \"seq_len\": no such variable".to_string())
    );
}

/// Checks any file the driver makes, as the committed one is checked: `HELIXLOOM_TEXTMAP_GOLDEN=file cargo test textmap --
/// --ignored` (the driver's first line says how to make one).
#[test]
#[ignore = "needs a file made with tests/tcl/text_map_driver.tcl"]
fn a_larger_file_from_the_driver_gives_the_same_pages() {
    let Ok(path) = std::env::var("HELIXLOOM_TEXTMAP_GOLDEN") else {
        return;
    };
    let cases = parse(&std::fs::read_to_string(path).expect("a file"));
    assert!(!cases.is_empty());
    check(&cases);
}

#[test]
fn a_long_sequence_with_many_features_is_made_into_a_page_without_trouble() {
    // 300,000 bases, 2,000 features of 300 bases (every 150th base on), on a circle, with the colours put on the DNA lines
    let seq: Vec<u8> = (0..300_000).map(|i| b"ACGT"[(i * 7 + i / 5) % 4]).collect();
    let features: Vec<MapFeature> = (0..2_000)
        .map(|i| MapFeature {
            name: format!("feature {i}"),
            kind: if i % 3 == 0 { "CDS" } else { "gene" }.to_string(),
            spans: vec![(i * 150, i * 150 + 300)],
            reverse: i % 2 == 1,
            hidden: false,
            fwd: [200, 100, 50],
            rev: [50, 100, 200],
        })
        .collect();
    let settings = Settings::default();
    let started = std::time::Instant::now();
    let page = text_map(&Request {
        seq: &seq,
        circular: true,
        title: "long",
        date: "DATE",
        region: Region::All,
        selection: None,
        row_width: 60,
        enzymes: &[],
        features: &features,
        code: STANDARD_CODE,
        settings: &settings,
    })
    .expect("a page")
    .expect("something to map");
    assert_eq!(page.dna_rows.len(), 3_000);
    assert!(page.tags.iter().any(|t| t.name == "f1999#"));
    // a guard against a page that takes minutes: it takes about a second in a build with the tests' settings
    assert!(started.elapsed().as_secs() < 120, "{:?}", started.elapsed());
}
