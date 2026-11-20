use super::*;

/// The source of a file of the program, by its path from the top of the repository.
fn source(file: &str) -> String {
    let path = format!("{}/{file}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// Why `test` is no test that runs with the others in this source, or nothing when it is one. It reads the text, so there is one
/// thing it cannot see: a file whose module its parent compiles only in some builds (the tests of the GPUs are so, and the table
/// says it of them).
fn fault(text: &str, test: &str) -> Option<&'static str> {
    if text.lines().any(|l| l.trim_start().starts_with("#![cfg")) {
        return Some("is in a file that only some builds compile");
    }
    let Some(at) = text.find(&format!("fn {test}(")) else {
        return Some("is not in the file");
    };
    // the attributes of a function are the lines just above it that start with `#[` (comments may lie between them)
    let before = text[..at].rsplit_once('\n').map_or("", |(lines, _)| lines);
    let attributes: Vec<&str> = before
        .lines()
        .rev()
        .map(str::trim)
        .take_while(|l| l.starts_with("#[") || l.starts_with("//"))
        .filter(|l| l.starts_with("#["))
        .collect();
    if !attributes.contains(&"#[test]") {
        return Some("is not marked as a test");
    }
    // `#[ignore]`, and `#[cfg_attr(…, ignore)]` too
    if attributes.iter().any(|a| a.contains("ignore")) {
        return Some("is left out of the usual run");
    }
    if attributes.iter().any(|a| a.contains("cfg")) {
        return Some("is compiled only in some builds");
    }
    None
}

#[test]
fn every_test_the_page_names_is_in_the_file_it_names_and_runs_with_the_others() {
    let mut wrong = Vec::new();
    for check in CHECKS {
        for (file, test) in check.tests {
            if let Some(why) = fault(&source(file), test) {
                wrong.push(format!("{}: {test} {why} ({file})", check.analysis));
            }
        }
    }
    assert!(wrong.is_empty(), "{}", wrong.join("\n"));
}

#[test]
fn a_test_that_is_missing_unmarked_or_left_out_is_told_from_one_that_runs() {
    let text = "#[test]\nfn runs() {}\n\nfn helper() {}\n\n#[test]\n#[ignore = \"needs files\"]\nfn left_out() {}\n\n    \
                #[test]\n    // a comment between\n    #[should_panic]\n    fn indented() {}\n\n#[test]\n\
                #[cfg(feature = \"gpu\")]\nfn gated() {}\n\n#[test]\n#[cfg_attr(not(unix), ignore)]\nfn sometimes() {}\n\n\
                #[test]\nfn runs_all() {}\n";
    assert_eq!(fault(text, "runs"), None);
    assert_eq!(fault(text, "indented"), None);
    assert_eq!(fault(text, "runs_all"), None);
    // a name that is only the start or the end of another's is not that test
    assert_eq!(fault(text, "run"), Some("is not in the file"));
    assert_eq!(fault(text, "all"), Some("is not in the file"));
    assert_eq!(fault(text, "nowhere"), Some("is not in the file"));
    assert_eq!(fault(text, "helper"), Some("is not marked as a test"));
    assert_eq!(
        fault(text, "sometimes"),
        Some("is left out of the usual run")
    );
    assert_eq!(
        fault(text, "gated"),
        Some("is compiled only in some builds")
    );
    assert_eq!(
        fault("#![cfg(unix)]\n#[test]\nfn runs() {}\n", "runs"),
        Some("is in a file that only some builds compile")
    );
    assert_eq!(
        fault(text, "left_out"),
        Some("is left out of the usual run")
    );
}

#[test]
fn a_check_that_claims_a_comparison_names_a_test_and_every_check_says_what_was_not_checked() {
    for check in CHECKS {
        let name = check.analysis;
        assert!(
            check.standard == Standard::Untested || !check.tests.is_empty(),
            "{name} claims a comparison and names no test"
        );
        for (what, text) in [
            ("an analysis", check.analysis),
            ("a group", check.group),
            ("what it was compared with", check.against),
            ("how close it came", check.result),
            ("what was not checked", check.limits),
        ] {
            assert!(!text.trim().is_empty(), "{name} has no text for {what}");
        }
    }
    let mut names: Vec<&str> = CHECKS.iter().map(|c| c.analysis).collect();
    names.sort_unstable();
    let count = names.len();
    names.dedup();
    assert_eq!(names.len(), count, "two checks have one name");
}

#[test]
fn the_checks_come_in_their_groups_in_the_order_of_the_table_and_each_once() {
    let all: Vec<&Check> = CHECKS.iter().collect();
    let groups = grouped(&all);
    assert_eq!(groups[0].0, CHECKS[0].group);
    let total: usize = groups.iter().map(|(_, list)| list.len()).sum();
    assert_eq!(total, CHECKS.len());
    for (group, list) in &groups {
        assert!(list.iter().all(|c| c.group == *group));
    }
    // a group is listed once, however its checks lie in the table
    let mut names: Vec<&str> = groups.iter().map(|(g, _)| *g).collect();
    let count = names.len();
    names.sort_unstable();
    names.dedup();
    assert_eq!(names.len(), count);
    assert!(grouped(&[]).is_empty());
}

#[test]
fn a_search_finds_the_checks_that_hold_every_word_whatever_the_case() {
    assert_eq!(matching("").len(), CHECKS.len());
    assert_eq!(matching("   ").len(), CHECKS.len());
    let astral: Vec<&str> = matching("ASTRAL 5.7.8")
        .iter()
        .map(|c| c.analysis)
        .collect();
    assert_eq!(astral, ["Species tree from gene trees"]);
    // the words of a grade find every check that has it, and a word of what was not checked finds a check too
    for standard in Standard::ALL {
        let found = matching(standard.label());
        for check in CHECKS.iter().filter(|c| c.standard == standard) {
            assert!(found.contains(&check), "{} by its grade", check.analysis);
        }
    }
    assert!(matching("pfam")
        .iter()
        .any(|c| c.analysis == "Protein domains (HMMER, Pfam)"));
    // every word must be there, not one of them
    assert!(matching("primer3 nosuchword").is_empty());
    assert!(matching("nosuchword").is_empty());
}

#[test]
fn a_check_is_found_by_a_word_of_any_of_its_texts_of_its_grade_and_of_its_tests() {
    let check = Check {
        group: "Groupword",
        analysis: "Nameword",
        standard: Standard::Independent,
        against: "Againstword.",
        result: "Resultword.",
        limits: "Limitsword.",
        tests: &[("src/fileword.rs", "testword_one")],
    };
    for word in [
        "groupword",
        "nameword",
        "second implementation",
        "againstword",
        "resultword",
        "limitsword",
        "fileword",
        "testword",
    ] {
        assert!(check.holds(word), "{word}");
    }
    assert!(!check.holds("otherword"));
    assert!(!check.holds("reference"));
}

#[test]
fn the_tally_counts_every_check_once_under_its_grade_strongest_first() {
    let all: Vec<&Check> = CHECKS.iter().collect();
    let counts = tally(&all);
    assert_eq!(
        counts.map(|(s, _)| s),
        [
            Standard::Reference,
            Standard::Independent,
            Standard::Known,
            Standard::Untested
        ]
    );
    assert_eq!(counts.iter().map(|(_, n)| n).sum::<usize>(), CHECKS.len());
    for (standard, n) in counts {
        assert_eq!(n, CHECKS.iter().filter(|c| c.standard == standard).count());
        assert!(n > 0, "no check is graded {}", standard.label());
    }
    assert_eq!(tally(&[]).map(|(_, n)| n), [0; 4]);
}

#[test]
fn the_text_gives_each_check_its_grade_what_it_was_compared_with_and_what_was_not_checked() {
    let a = Check {
        group: "Trees",
        analysis: "A made-up method",
        standard: Standard::Known,
        against: "Three cases worked by hand.",
        result: "Equal.",
        limits: "Nothing else.",
        tests: &[("src/a.rs", "first_case"), ("src/b.rs", "second_case")],
    };
    let b = Check {
        group: "Trees",
        analysis: "Another",
        standard: Standard::Untested,
        against: "Nothing.",
        result: "Unknown.",
        limits: "All of it.",
        tests: &[],
    };
    assert_eq!(
        to_text(&[&a, &b]),
        "How each analysis was checked: 2 analyses \
         (0 reference program, 0 second implementation, 1 known answers, 1 not checked)\n\
         \n\
         Trees\n\
         \n\
         A made-up method [known answers]\n\
         \x20 compared with: Three cases worked by hand.\n\
         \x20 result: Equal.\n\
         \x20 not checked: Nothing else.\n\
         \x20 tests: first_case (src/a.rs); second_case (src/b.rs)\n\
         \n\
         Another [not checked]\n\
         \x20 compared with: Nothing.\n\
         \x20 result: Unknown.\n\
         \x20 not checked: All of it.\n\
         \x20 tests: none\n"
    );
    // one analysis is said in the singular
    assert!(to_text(&[&a]).starts_with("How each analysis was checked: 1 analysis ("));
}

#[test]
fn every_grade_has_a_label_and_a_line_that_says_what_it_means() {
    for standard in Standard::ALL {
        assert!(!standard.label().is_empty() && standard.meaning().ends_with('.'));
    }
    assert!(Standard::Reference < Standard::Untested);
}

#[test]
fn a_result_is_led_to_the_check_of_the_method_that_made_it_by_its_kind_and_its_settings() {
    let of = |analysis: &str, settings: &[(&str, &str)]| {
        let settings: Vec<(String, String)> = settings
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
        of_result(analysis, &settings).map(|c| c.analysis)
    };
    for (analysis, settings, check) in [
        (
            "Alignment",
            vec![("engine", "built in")],
            "Alignment (built in)",
        ),
        (
            "Alignment",
            vec![("engine", "progressive (ClustalW-style)")],
            "Alignment (progressive, after ClustalW)",
        ),
        (
            "Alignment",
            vec![("engine", "MAFFT")],
            "MAFFT, IQ-TREE and the other programs HelixLoom runs",
        ),
        (
            "Tree estimation",
            vec![
                ("method", "Neighbour-joining"),
                ("distances", "Jukes-Cantor"),
            ],
            "Distance trees (neighbour-joining, UPGMA)",
        ),
        (
            "Tree estimation",
            vec![("method", "UPGMA")],
            "Distance trees (neighbour-joining, UPGMA)",
        ),
        (
            "Tree estimation",
            vec![("method", "Maximum likelihood")],
            "Maximum-likelihood trees",
        ),
        (
            "Tree estimation",
            vec![("method", "Bayesian")],
            "Bayesian trees",
        ),
        (
            "Species tree",
            vec![("weights", "none"), ("gene families", "no")],
            "Species tree from gene trees",
        ),
        (
            "Species tree",
            vec![("weights", "hybrid"), ("gene families", "no")],
            "Weighted species trees (by support and branch length)",
        ),
        (
            "Species tree",
            vec![("weights", "none"), ("gene families", "yes")],
            "Gene families (duplications and losses)",
        ),
        (
            "Motif search",
            vec![("width", "10")],
            "Motif Finder (scan, discovery, comparison)",
        ),
        (
            "Protein domain search",
            vec![],
            "Protein domains (HMMER, Pfam)",
        ),
        (
            "Structure prediction",
            vec![],
            "Structure prediction (a folding program)",
        ),
        (
            "PhyloCNN estimates",
            vec![],
            "PhyloCNN (epidemic parameters from a tree)",
        ),
        (
            "Pipeline file",
            vec![],
            "Pipeline files (run, cache, export)",
        ),
    ] {
        // the name is that of a check in the table: a check that is renamed fails here
        assert_eq!(
            of(analysis, &settings),
            Some(check),
            "{analysis} {settings:?}"
        );
    }
    assert_eq!(
        of("Tree estimation", &[("method", "A method of tomorrow")]),
        None
    );
    assert_eq!(of("Something else", &[]), None);
    // a run of the fixed route is several methods: no one check is its own
    assert_eq!(of("Pipeline", &[("steps", "align, tree, primers")]), None);
    // an engine that the table does not know is not taken for the built-in one
    assert_eq!(
        of("Alignment", &[("engine", "an aligner of tomorrow")]),
        None
    );
    assert_eq!(of("Alignment", &[]), None);
}
