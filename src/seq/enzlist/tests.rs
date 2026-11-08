use super::*;
use crate::seq::enzymes::{all_sites, default_enzymes};
use crate::seq::sitesearch::golden::{self, Case};

/// 72 cases made by running the original's `list_enzymes_window` (with `list_enzymes_text`, `enz_site_usage_filter`,
/// `enz_count` and `findenzymes` behind it) on random sequences: the sequence, the selection, the enzymes chosen in the
/// Enzyme Selector, whether Dam/Dcm methylation was on, the choices of the dialog, the sites the original found, and the text
/// of its window.
const GOLDEN: &str = include_str!("golden.txt");

/// 48 more, made with the driver's dense mode (seed 4242): runs of a site repeated at its own period, so that sites overlap
/// as far as they can, and half of the selections starting and ending on the first bases of sites or one before them.
const DENSE: &str = include_str!("golden_dense.txt");

fn entry(name: &str, pattern: &str, sites: &[usize]) -> Entry {
    Entry {
        name: name.into(),
        pattern: pattern.into(),
        sites: sites.to_vec(),
    }
}

/// The entries of the default enzymes for a case: their sites as the build finds them, numbered from 1.
fn entries_of(case: &Case) -> Vec<Entry> {
    let methylated = case.say("methylated") == "1";
    let enzymes = default_enzymes();
    let sites = all_sites(&enzymes, &case.seq, case.circular, methylated);
    enzymes
        .iter()
        .zip(sites)
        .map(|(e, s)| Entry {
            name: e.name.clone(),
            pattern: e.pattern_text(),
            sites: s.iter().map(|site| site.start + 1).collect(),
        })
        .collect()
}

fn format_of(case: &Case) -> Format {
    Format::ALL[case
        .say("format")
        .parse::<usize>()
        .expect("a format number")]
}

fn which_of(case: &Case) -> Which {
    let name = case.say("which");
    let found = Which::ALL.into_iter().find(|w| {
        let word = match w {
            Which::All => "all",
            Which::Selected => "selected",
            Which::Absent => "absent",
            Which::Unique => "unique",
        };
        word == name
    });
    found.expect("a list name")
}

/// The bases looked at: the selection when the case asks for it and there is one.
fn span_of(case: &Case) -> Option<Span> {
    let asked = case.say("region") == "selection";
    case.selection
        .filter(|_| asked)
        .map(|(start, end)| Span { start, end })
}

/// The text of every case as HelixLoom writes it, against the original's window text.
fn check_lists(cases: &[Case]) {
    for (n, case) in cases.iter().enumerate() {
        let entries = entries_of(case);
        let (format, which, span) = (format_of(case), which_of(case), span_of(case));
        let columns: usize = case.say("columns").parse().expect("a column count");
        let selected = case.chosen.iter().cloned().collect();
        let chosen = choose(&entries, which, &selected, span);
        let body = text(&chosen, format, columns, span);
        let report = report("DATE", "TITLE", which, span, case.seq.len(), &body);
        assert_eq!(report, case.text, "case {}", n + 1);
    }
}

/// The sites of all 116 enzymes of every case, against the lists the original's `findenzymes` made: what differs, one line each.
fn site_differences(cases: &[Case]) -> Vec<String> {
    let mut found = Vec::new();
    for (n, case) in cases.iter().enumerate() {
        for e in entries_of(case) {
            let theirs = case.sites.iter().find(|(name, _)| *name == e.name);
            let theirs = theirs.map(|(_, at)| at.clone()).unwrap_or_default();
            if e.sites != theirs {
                let line = format!(
                    "case {}: sites of {}: {:?}, the original's {theirs:?}",
                    n + 1,
                    e.name,
                    e.sites
                );
                found.push(line);
            }
        }
    }
    found
}

fn check_sites(cases: &[Case]) {
    let differences = site_differences(cases);
    assert!(differences.is_empty(), "{}", differences.join("\n"));
}

#[test]
fn lists_equal_the_originals_on_every_case() {
    let cases = golden::parse(GOLDEN);
    assert_eq!(cases.len(), 72);
    check_lists(&cases);
}

#[test]
fn the_sites_behind_the_lists_equal_the_originals_in_number_and_order() {
    check_sites(&golden::parse(GOLDEN));
}

#[test]
fn lists_and_sites_equal_the_originals_on_the_dense_cases() {
    let cases = golden::parse(DENSE);
    assert_eq!(cases.len(), 48);
    check_sites(&cases);
    check_lists(&cases);
}

#[test]
fn the_dense_cases_put_sites_on_each_edge_of_the_selection() {
    // the number of cases in which some site is numbered `edge(start, end)` with a selection that was asked for
    let cases = golden::parse(DENSE);
    let having = |edge: fn(usize, usize) -> usize| {
        let on_edge = |case: &&Case| {
            let Some(Span { start, end }) = span_of(case) else {
                return false;
            };
            let place = edge(start, end);
            case.sites.iter().any(|(_, at)| at.contains(&place))
        };
        cases.iter().filter(on_edge).count()
    };
    assert!(
        having(|a, _| a) >= 3,
        "a site on the base before the selection"
    );
    assert!(having(|a, _| a + 1) >= 3, "a site on its first base");
    assert!(having(|_, b| b) >= 3, "a site on its last base");
    assert!(having(|_, b| b + 1) >= 3, "a site on the base after it");
}

/// The same two checks on a file that the driver made with other seeds and more cases, for a wider comparison:
/// `HELIXLOOM_ENZLIST_GOLDEN=/tmp/golden.txt cargo test enzlist -- --ignored`
#[test]
#[ignore = "needs HELIXLOOM_ENZLIST_GOLDEN, a file made by tests/tcl/enzyme_list_driver.tcl"]
fn a_larger_file_made_by_the_driver_gives_the_same_lists_and_sites() {
    let path =
        std::env::var("HELIXLOOM_ENZLIST_GOLDEN").expect("HELIXLOOM_ENZLIST_GOLDEN names a file");
    let file = std::fs::read_to_string(path).expect("the file reads");
    let cases = golden::parse(&file);
    assert!(!cases.is_empty(), "no cases in the file");
    check_sites(&cases);
    check_lists(&cases);
}

#[test]
fn the_cases_cover_every_format_list_and_both_regions() {
    let cases = golden::parse(GOLDEN);
    for f in 0..4 {
        let has = |c: &&Case| c.say("format") == f.to_string();
        assert!(cases.iter().filter(has).count() >= 10, "format {f}");
    }
    for which in ["all", "selected", "absent", "unique"] {
        assert!(cases.iter().any(|c| c.say("which") == which), "{which}");
    }
    let (with, without) = cases
        .iter()
        .partition::<Vec<_>, _>(|c| span_of(c).is_some());
    assert!(with.len() >= 20 && without.len() >= 20);
    assert!(cases.iter().any(|c| c.say("methylated") == "1"));
    assert!(cases.iter().any(|c| c.say("columns") == "1"));
    assert!(cases.iter().any(|c| c.chosen.is_empty()));
}

#[test]
fn the_enzyme_names_and_patterns_equal_the_enzyme_file_the_original_reads() {
    // Format 2 writes each pattern as the file gives it; one case lists them all
    let cases = golden::parse(GOLDEN);
    let case = cases
        .iter()
        .find(|c| {
            format_of(c) == Format::Positions && which_of(c) == Which::All && span_of(c).is_none()
        })
        .expect("a case that lists every enzyme with its pattern");
    let all = entries_of(case);
    let listed = case.text.lines().filter(|l| l.len() > 30).count();
    assert_eq!(listed, all.len());
    for e in &all {
        let line = format!("{:<15}{:<15}", e.name, e.pattern);
        assert!(case.text.contains(&line), "{} {}", e.name, e.pattern);
    }
}

fn three() -> Vec<Entry> {
    vec![
        entry("SmaI", "CCC^GGG", &[10, 40]),
        entry("EcoRI", "G^AATTC", &[20]),
        entry("PstI", "CTGCA^G", &[]),
    ]
}

#[test]
fn names_fill_the_columns_and_a_row_that_is_not_full_has_no_newline() {
    let entries = three();
    let all: Vec<&Entry> = entries.iter().collect();

    let two = text(&all, Format::Names, 2, None);
    let one = text(&all, Format::Names, 1, None);

    assert_eq!(two, "EcoRI          PstI           \nSmaI           ");
    assert_eq!(one, "EcoRI          \nPstI           \nSmaI           \n");
}

#[test]
fn names_and_counts_are_padded_to_twelve_and_eight() {
    let entries = three();
    let all: Vec<&Entry> = entries.iter().collect();

    let list = text(&all, Format::Counts, 3, None);

    assert_eq!(
        list,
        "EcoRI       1       PstI        0       SmaI        2       \n"
    );
}

#[test]
fn positions_go_the_fewest_sites_first_and_equal_counts_keep_their_order() {
    let mut entries = three();
    entries.push(entry("ApaI", "GGGCC^C", &[5]));
    let all: Vec<&Entry> = entries.iter().collect();

    let list = text(&all, Format::Positions, 6, None);

    // no sites: the count is padded to seven and so is the empty list of places
    let expected = "PstI           CTGCA^G        0             \n\
                    EcoRI          G^AATTC        1      20     \n\
                    ApaI           GGGCC^C        1      5      \n\
                    SmaI           CCC^GGG        2      10 40  \n";
    assert_eq!(list, expected);
}

#[test]
fn by_site_sorts_by_the_pattern_as_written() {
    let entries = three();
    let all: Vec<&Entry> = entries.iter().collect();

    let list = text(&all, Format::BySite, 6, None);

    let expected = "CCC^GGG        SmaI           \n\
                    CTGCA^G        PstI           \n\
                    G^AATTC        EcoRI          \n";
    assert_eq!(list, expected);
}

#[test]
fn a_field_longer_than_its_width_is_not_cut() {
    let entries = [entry("XcmI", "CCANNNNN^NNNNTGG", &[3])];
    let all: Vec<&Entry> = entries.iter().collect();

    let list = text(&all, Format::BySite, 6, None);

    assert_eq!(list, "CCANNNNN^NNNNTGGXcmI           \n");
}

#[test]
fn the_selection_counts_sites_from_its_first_base_but_the_list_shows_one_base_earlier() {
    // the selection is bases 11 to 20: ten bases come before it
    let span = Some(Span { start: 10, end: 20 });
    let entries = [
        entry("Before", "A^AAA", &[10]),
        entry("First", "A^AAA", &[11]),
        entry("Last", "A^AAA", &[20]),
        entry("After", "A^AAA", &[21]),
    ];
    let counts: Vec<usize> = entries.iter().map(|e| e.count_in(span)).collect();
    let all: Vec<&Entry> = entries.iter().collect();

    let list = text(&all, Format::Counts, 1, span);

    assert_eq!(counts, [0, 1, 1, 0]);
    let expected =
        "After       0       \nBefore      1       \nFirst       1       \nLast        1       \n";
    assert_eq!(list, expected);
}

#[test]
fn absent_and_unique_look_at_the_selected_bases_only() {
    let span = Some(Span { start: 10, end: 20 });
    let entries = vec![
        entry("Inside", "A^AAA", &[15]),
        entry("Twice", "A^AAA", &[12, 18]),
        entry("Outside", "A^AAA", &[30]),
        entry("None", "A^AAA", &[]),
    ];
    let none = BTreeSet::new();

    let absent: Vec<&str> = choose(&entries, Which::Absent, &none, span)
        .iter()
        .map(|e| e.name.as_str())
        .collect();
    let unique: Vec<&str> = choose(&entries, Which::Unique, &none, span)
        .iter()
        .map(|e| e.name.as_str())
        .collect();
    let absent_all: Vec<&str> = choose(&entries, Which::Absent, &none, None)
        .iter()
        .map(|e| e.name.as_str())
        .collect();

    assert_eq!(absent, ["Outside", "None"]);
    assert_eq!(unique, ["Inside"]);
    assert_eq!(absent_all, ["None"]);
}

#[test]
fn selected_enzymes_come_in_sorted_order_and_none_selected_lists_nothing() {
    let entries = three();
    let picked: BTreeSet<String> = ["SmaI", "EcoRI", "Gone"]
        .iter()
        .map(|s| s.to_string())
        .collect();

    let some: Vec<&str> = choose(&entries, Which::Selected, &picked, None)
        .iter()
        .map(|e| e.name.as_str())
        .collect();
    let none = choose(&entries, Which::Selected, &BTreeSet::new(), None);

    assert_eq!(some, ["EcoRI", "SmaI"]);
    assert!(none.is_empty());
    assert_eq!(text(&none, Format::Positions, 6, None), "");
}

#[test]
fn the_heading_leaves_a_space_after_the_name_and_the_bases() {
    let all = report("DATE", "pUC19", Which::Unique, None, 2686, "BODY");
    let part = report(
        "DATE",
        "pUC19",
        Which::Selected,
        Some(Span { start: 10, end: 20 }),
        2686,
        "",
    );

    assert_eq!(
        all,
        "DATE\npUC19 \nFrom 0 to 2686. \n\nUnique Sites\n\nBODY"
    );
    assert_eq!(part, "DATE\npUC19 \nFrom 10 to 20. \n\nSelected Sites\n\n");
}

#[test]
fn the_dialog_choices_are_named_as_in_the_original() {
    let formats: Vec<&str> = Format::ALL.iter().map(|f| f.label()).collect();
    let lists: Vec<&str> = Which::ALL.iter().map(|w| w.label()).collect();
    let titles: Vec<&str> = Which::ALL.iter().map(|w| w.title()).collect();

    assert_eq!(
        formats,
        [
            "Name",
            "Name, #Sites",
            "Name, Site, #Sites, Positions",
            "Site, Name"
        ]
    );
    assert_eq!(
        lists,
        [
            "All",
            "Selected Enzymes Only",
            "Absent Enzymes",
            "Unique Enzymes"
        ]
    );
    assert_eq!(
        titles,
        [
            "All Enzymes",
            "Selected Sites",
            "Absent Sites",
            "Unique Sites"
        ]
    );
    let with_columns: Vec<bool> = Format::ALL.iter().map(|f| f.uses_columns()).collect();
    assert_eq!(with_columns, [true, true, false, false]);
}

#[test]
fn a_column_count_outside_one_to_ten_is_kept_inside_it() {
    let entries: Vec<Entry> = (0..12)
        .map(|k| entry(&format!("E{k:02}"), "A^A", &[]))
        .collect();
    let all: Vec<&Entry> = entries.iter().collect();

    let none = text(&all, Format::Names, 0, None);
    let many = text(&all, Format::Names, 99, None);

    assert_eq!(none.lines().count(), 12);
    assert_eq!(many.lines().count(), 2);
}
