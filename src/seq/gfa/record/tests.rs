use super::{parse, Record};

#[test]
fn reported_lengths_distinguish_sequence_declared_zero_and_unknown() {
    for (fields, length) in [
        ("aB.=z", Some(5)),
        ("aB.=z\tLN:i:5", Some(5)),
        ("*\tLN:i:+5", Some(5)),
        ("*\tLN:i:0", Some(0)),
        ("*\tLN:i:18446744073709551615", Some(u64::MAX)),
        ("*", None),
    ] {
        assert_eq!(
            parse(&format!("S\tread=*,1\t{fields}")),
            Ok(Record::Segment {
                name: "read=*,1",
                length,
                depth: None,
                kmer_count: None,
            })
        );
    }
    for fields in [
        "",
        "A*",
        "A C",
        "AB7",
        "A\tLN:i:2",
        "*\tLN:i:18446744073709551616",
    ] {
        assert!(parse(&format!("S\ts\t{fields}")).is_err(), "{fields}");
    }
}

#[test]
fn checks_the_name_on_segments_and_both_link_ends() {
    for id in [
        "",
        "has space",
        "=a",
        "*",
        "b+,c",
        "b-,c",
        "x\u{7f}",
        "é",
        "a\n",
    ] {
        for line in [
            format!("S\t{id}\tAC"),
            format!("L\t{id}\t+\ts\t-\t*"),
            format!("L\ts\t+\t{id}\t-\t*"),
        ] {
            assert!(parse(&line).is_err(), "{line:?}");
        }
    }
}

#[test]
fn checks_known_tags_and_accepts_printable_annotations_with_colons() {
    for v in ["1.0", "1.1", "1.2"] {
        assert_eq!(
            parse(&format!("H\tVN:Z:{v}\tx1:Z:my graph")),
            Ok(Record::Header)
        );
    }
    for bad in ["VN:i:1", "VN:Z:2.0", "VN:Z:1", "VN:Z:1.3"] {
        assert!(parse(&format!("H\t{bad}")).is_err(), "{bad}");
    }
    assert_eq!(parse("H"), Ok(Record::Header));
    assert_eq!(
        parse("S\ts\tAA\tLN:i:2\tKC:i:+4\tDP:f:+.25E-1\tUR:Z:file:/assembly.fa"),
        Ok(Record::Segment {
            name: "s",
            length: Some(2),
            depth: Some(0.025),
            kmer_count: Some(4),
        })
    );
    for field in [
        "LN:Z:2",
        "LN:i:-1",
        "LN:i:2.0",
        "KC:i:1.5",
        "KC:i:-1",
        "KC:f:1.0",
        "DP:i:3",
        "DP:f:-1",
        "DP:f:-1e-100",
        "DP:f:nan",
        "DP:f:inf",
        "DP:f:7e100",
        "DP:f:1.e3",
        "kc:i:1\tkc:i:2",
        "LN:i:2\tLN:i:2",
        "DP:f:0.25\tDP:f:1",
        "KC:i:0\tKC:i:2",
        "2X:Z:abc",
        "x?:i:2",
        "abc:i:3",
        "id::a",
        "id:Q:3",
        "id:Z:",
        "id:Z:a\u{1b}",
        "annotation",
    ] {
        assert!(parse(&format!("S\ts\tAA\t{field}")).is_err(), "{field}");
    }
    assert!(parse("L\ta\t+\tb\t-\t*\tID:Z:x\tID:Z:y").is_err());
}

#[test]
fn segment_metadata_preserves_independent_counts_and_decimal_depth() {
    for (fields, depth, kmer_count) in [
        ("", None, None),
        ("\tDP:f:0", Some(0.0), None),
        ("\tKC:i:0", None, Some(0)),
        ("\tKC:i:0\tDP:f:-0.0", Some(-0.0), Some(0)),
        ("\tDP:f:1e-100\tLN:i:2\tKC:i:1", Some(1e-100), Some(1)),
        (
            "\tKC:i:18446744073709551615\tDP:f:12.3456789012345",
            Some(12.3456789012345),
            Some(u64::MAX),
        ),
    ] {
        assert_eq!(
            parse(&format!("S\ts\tAC{fields}")),
            Ok(Record::Segment {
                name: "s",
                length: Some(2),
                depth,
                kmer_count
            }),
        );
    }
}

#[test]
fn cigar_counts_the_from_segment_as_reference_and_the_to_segment_as_query() {
    for (cigar, from_span, to_span) in [
        ("0M", Some(0), Some(0)),
        ("*", None, None),
        ("31M", Some(31), Some(31)),
        ("1M2I4D8N16S32H64P128X256=", Some(397), Some(403)),
        ("18446744073709551615D", Some(u64::MAX), Some(0)),
        ("18446744073709551615I", Some(0), Some(u64::MAX)),
    ] {
        for (ao, bo) in [("+", "+"), ("-", "+"), ("+", "-"), ("-", "-")] {
            assert_eq!(
                parse(&format!("L\tA\t{ao}\tB\t{bo}\t{cigar}")),
                Ok(Record::Link {
                    from: "A",
                    to: "B",
                    from_span,
                    to_span
                })
            );
        }
    }
}

#[test]
fn corrupt_links_never_become_zero_length_or_unknown_overlaps() {
    for cigar in [
        "",
        "3",
        "M",
        "-2M",
        "2m",
        "1M+",
        "1M1Z",
        "1M2",
        "1M*",
        "18446744073709551616M",
        "18446744073709551615D1M",
        "18446744073709551615I1M",
    ] {
        assert!(
            parse(&format!("L\ta\t-\tb\t+\t{cigar}")).is_err(),
            "{cigar}"
        );
    }
    for line in [
        "L\ta\tx\tb\t+\t*",
        "L\ta\t-\tb\tz\t*",
        "L\ta\t?\tb\t+\t*",
        "L\ta\t+\tb\t?\t*",
        "L\ta\t-\tb\t-",
    ] {
        assert!(parse(line).is_err(), "{line}");
    }
}

#[test]
fn a_record_outside_the_supported_subset_is_an_error() {
    assert_eq!(parse(""), Ok(Record::Skip));
    assert_eq!(parse("# comment\tS\tbad"), Ok(Record::Skip));
    for record in ["P", "W", "C", "J", "X", "S id A"] {
        let err = parse(record).unwrap_err();
        assert!(err.contains(record), "{err}");
    }
    for line in ["S", "S\ta", "L", "L\ta", "L\ta\t+", "L\ta\t+\tb"] {
        assert!(parse(line).is_err(), "{line}");
    }
}
