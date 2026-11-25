use super::*;

fn error(input: &str) -> String {
    inspect(input.as_bytes()).expect_err("the graph should be rejected")
}

#[test]
fn summary_counts_link_records_and_components_with_unknown_lengths() {
    let input = "# example\r\nH\tVN:Z:1.0\r\n\
                 L\ta\t+\tb\t-\t2M\n\
                 S\tc\t*\tLN:i:10\n\
                 L\tb\t+\ta\t-\t2M\n\
                 S\ta\tACG\n\
                 L\tc\t-\tc\t+\t*\n\
                 S\td\t*\n\
                 S\tb\tAC";
    assert_eq!(
        inspect(input.as_bytes()).unwrap(),
        Summary {
            segments: 4,
            links: 3,
            components: 3,
            isolated: 1,
            known_bases: 15,
            unknown_lengths: 1,
        }
    );
}

#[test]
fn link_between_two_connected_groups_merges_the_whole_component() {
    let input = "S\ta\tA\nS\tb\tA\nS\tc\tA\nS\td\tA\n\
                 L\ta\t+\tb\t+\t1M\nL\tc\t+\td\t+\t1M\n\
                 L\ta\t+\tc\t+\t1M\n";
    let summary = inspect(input.as_bytes()).unwrap();
    assert_eq!(
        (summary.components, summary.isolated, summary.links),
        (1, 0, 3)
    );
    let with_cycle = format!("{input}L\tb\t+\td\t+\t1M\n");
    let summary = inspect(with_cycle.as_bytes()).unwrap();
    assert_eq!(
        (summary.components, summary.isolated, summary.links),
        (1, 0, 4)
    );
}

#[test]
fn header_only_is_an_empty_graph_but_comments_alone_are_not_a_graph() {
    assert_eq!(
        inspect("H\tVN:Z:1.0\n".as_bytes()).unwrap(),
        Summary::default()
    );
    for input in ["", "\n", "# H\tVN:Z:1.0\n\n"] {
        assert!(error(input).contains("no GFA"));
    }
}

#[test]
fn duplicate_segment_is_reported_even_if_the_first_has_no_length() {
    let why = error("# comment\nS\tx\t*\n\nS\tx\tA");
    assert!(
        why.contains("line 4") && why.contains("line 2") && why.contains("'x'"),
        "{why}"
    );
}

#[test]
fn undefined_endpoints_are_reported_at_the_referencing_line() {
    for input in [
        "H\tVN:Z:1.0\nL\tx\t+\ta\t-\t*\nS\ta\tA",
        "S\ta\tA\nL\ta\t+\tx\t-\t0M",
        "H\tVN:Z:1.0\nL\tx\t+\tx\t+\t1M",
    ] {
        let why = error(input);
        assert!(
            why.contains("line 2") && why.contains("undefined") && why.contains("'x'"),
            "{why}"
        );
    }
}

#[test]
fn bounds_use_the_largest_consumption_at_either_end_across_all_links() {
    // These links precede both definitions. The self-link consumes one from x and three from x.
    let input = "L\tx\t+\tx\t-\t1M2I\nL\tx\t+\ty\t+\t1M\nL\tx\t-\ty\t+\t*\nS\tx\tAA\nS\ty\tA";
    let why = error(input);
    assert!(
        why.contains("line 1") && why.contains("3 bases") && why.contains("'x'"),
        "{why}"
    );
    let input = "S\tx\tAA\nS\ty\tAAA\nL\tx\t-\ty\t+\t1M2D";
    let why = error(input);
    assert!(why.contains("line 3") && why.contains("'x'"), "{why}");
    // The same CIGAR fits when the endpoint consuming three bases has no known length.
    let summary = inspect("S\tx\t*\nS\ty\tA\nL\tx\t-\ty\t+\t1M2D".as_bytes()).unwrap();
    assert_eq!((summary.known_bases, summary.unknown_lengths), (1, 1));
}

#[test]
fn zero_is_a_known_length_and_total_length_cannot_wrap() {
    let empty = inspect("S\tx\t*\tLN:i:0\nL\tx\t+\tx\t+\t0M".as_bytes()).unwrap();
    assert_eq!(
        (empty.known_bases, empty.unknown_lengths, empty.isolated),
        (0, 0, 0)
    );
    let input = format!("S\tx\t*\tLN:i:{}\nS\ty\tA", u64::MAX);
    let why = error(&input);
    assert!(
        why.contains("line 2") && why.contains("sum of segment lengths"),
        "{why}"
    );
}

#[test]
fn malformed_record_and_read_errors_include_the_input_line() {
    assert!(error("# comment\nS\tx").contains("line 2"));
    assert!(inspect(&b"S\tx\tA\n\xff"[..])
        .unwrap_err()
        .contains("line 2"));
}

#[test]
fn summary_tsv_is_the_cli_report_and_reads_workstation_size_counts() {
    let text = "segments\t3\nlinks\t2\ncomponents\t2\nisolated_segments\t1\nknown_segment_bases\t0\nsegments_without_length\t2\n";
    let graph = b"S\ta\t*\nS\tb\t*\tLN:i:0\nS\tc\t*\nL\ta\t+\tb\t-\t*\nL\tb\t-\tb\t+\t0M\n";
    let expected = inspect(&graph[..]).unwrap();
    assert_eq!(expected.to_tsv(), text);
    assert_eq!(Summary::from_tsv(text).unwrap(), expected);
    let reordered = text.lines().rev().collect::<Vec<_>>().join("\r\n");
    assert_eq!(Summary::from_tsv(&reordered).unwrap(), expected);

    let large = Summary {
        segments: u64::MAX,
        links: u64::MAX - 1,
        components: 1,
        isolated: 0,
        known_bases: u64::MAX,
        unknown_lengths: u64::MAX - 1,
    };
    for s in [large, Summary::default()] {
        assert_eq!(Summary::from_tsv(&s.to_tsv()).unwrap(), s);
    }
}

#[test]
fn summary_tsv_rejects_missing_duplicate_unknown_and_corrupt_counts() {
    let good = Summary::default().to_tsv();
    for (text, diagnostic) in [
        (String::new(), "missing"),
        (good.replace("links\t0\n", ""), "links"),
        (good.replace("links\t0", "segments\t0"), "duplicate"),
        (good.replace("components\t0", "component\t0"), "unknown"),
        (good.replace("segments\t0\n", "segments 0\n"), "line"),
        (good.replace("links\t0", "links\t1\t2"), "links"),
        (good.replace("links\t0", "links\t-1"), "links"),
        (good.replace("links\t0", "links\t1.0"), "links"),
        (
            good.replace("links\t0", "links\t18446744073709551616"),
            "links",
        ),
        (format!("{good}extra\t1\n"), "unknown"),
    ] {
        let why = Summary::from_tsv(&text).unwrap_err();
        assert!(why.contains(diagnostic), "{text:?}: {why}");
    }
}

#[test]
fn summary_tsv_checks_counts_for_graphs_with_self_links_and_repeated_links() {
    // (segments, links, components, isolated): a linked singleton needs a self-link.
    for (segments, links, components, isolated, possible) in [
        (0, 0, 0, 0, true),
        (0, 1, 0, 0, false),
        (1, 2, 1, 0, true),
        (1, 0, 1, 0, false),
        (1, 2, 1, 1, false),
        (2, 0, 2, 2, true),
        (2, 1, 3, 0, false),
        (2, 1, 1, 2, false),
        (2, 1, 0, 0, false),
        (2, 1, 1, 1, false),
        (3, 1, 1, 0, false),
        (3, 2, 1, 0, true),
        (3, 2, 3, 0, false),
        (3, 3, 3, 0, true),
    ] {
        let summary = Summary {
            segments,
            links,
            components,
            isolated,
            ..Summary::default()
        };
        assert_eq!(
            Summary::from_tsv(&summary.to_tsv()).is_ok(),
            possible,
            "{summary:?}"
        );
    }
    let isolated = Summary {
        segments: 1,
        isolated: 1,
        components: 1,
        ..Summary::default()
    };
    for (unknown_lengths, known_bases, possible) in
        [(2, 0, false), (1, 1, false), (1, 0, true), (0, 0, true)]
    {
        let summary = Summary {
            unknown_lengths,
            known_bases,
            ..isolated
        };
        assert_eq!(
            Summary::from_tsv(&summary.to_tsv()).is_ok(),
            possible,
            "{summary:?}"
        );
    }
}

const TABLE_HEADER: &str = "segment\tlength_bp\tdepth\tkmer_count\tlink_records\n";

fn table(input: &str) -> String {
    let mut out = Vec::new();
    segment_table(input.as_bytes(), &mut out).unwrap();
    String::from_utf8(out).unwrap()
}

#[test]
fn segment_rows_follow_declarations_and_count_each_incident_link_record_once() {
    let input = "L\ta\t+\tb\t-\t*\n\
                 L\ta\t+\tb\t-\t*\n\
                 L\tb\t+\tb\t-\t0M\n\
                 S\tb\t*\tLN:i:0\tKC:i:0\tDP:f:-0.000\n\
                 S\tx\t*\n\
                 S\ta\tAAA\tKC:i:18446744073709551615\tDP:f:1.000000000000009\n\
                 L\ta\t-\ta\t-\t1M\n";
    assert_eq!(
        table(input),
        format!(
            "{TABLE_HEADER}\
        b\t0\t0\t0\t3\n\
        x\t.\t.\t.\t0\n\
        a\t3\t1.000000000000009\t18446744073709551615\t3\n"
        )
    );
}

#[test]
fn header_only_writes_column_names_and_depths_keep_small_and_large_values() {
    assert_eq!(table("H\tVN:Z:1.0\n"), TABLE_HEADER);
    let input = "S\tx\tA\tDP:f:1.23e-20\nS\ty\tA\tDP:f:5.67e30\n";
    assert_eq!(
        table(input),
        format!("{TABLE_HEADER}x\t1\t1.23e-20\t.\t0\ny\t1\t5.67e30\t.\t0\n")
    );
}

#[test]
fn segment_table_reports_read_and_cross_record_errors_before_writing_the_header() {
    for invalid in [
        &b"S\ta\tAAA\nL\ta\t-\tb\t+\t0M\n"[..],
        &b"S\ta\tAAA\nL\ta\t-\ta\t+\t4M\n"[..],
        &b"S\ta\tAAA\n\xff"[..],
    ] {
        let mut out = Vec::new();
        let why = segment_table(invalid, &mut out).unwrap_err();
        assert!(why.contains("line 2"), "{why}");
        assert!(out.is_empty());
    }
}

#[test]
fn segment_table_propagates_errors_writing_rows() {
    let mut bytes = vec![0; TABLE_HEADER.len()];
    let why = segment_table("S\ta\tA\n".as_bytes(), bytes.as_mut_slice()).unwrap_err();
    assert_eq!(bytes, TABLE_HEADER.as_bytes());
    assert!(why.contains("standard output"), "{why}");
}
