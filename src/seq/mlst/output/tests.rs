use super::*;
use crate::seq::mlst::call::Status;

fn typing() -> Typing {
    Typing {
        scheme: "neis".into(),
        calls: vec![
            ("abcZ".into(), Call::Exact(vec!["3".into()])),
            (
                "adk".into(),
                Call::Novel(vec![("7".into(), b"ACGT".to_vec())]),
            ),
            ("fumC".into(), Call::Partial("4".into())),
            ("gdh".into(), Call::Missing),
            ("pdhC".into(), Call::Exact(vec!["1".into(), "2".into()])),
        ],
        sequence_type: None,
        score: 71.5,
        status: Status::Ok,
    }
}

#[test]
fn a_line_has_the_file_the_scheme_the_type_and_each_locus_with_its_call() {
    let t = typing();
    assert_eq!(
        line("g.fa", Some(&t), '\t'),
        "g.fa\tneis\t-\tabcZ(3)\tadk(~7)\tfumC(4?)\tgdh(-)\tpdhC(1,2)"
    );
    assert_eq!(
        line("g.fa", Some(&t), ','),
        "g.fa,neis,-,abcZ(3),adk(~7),fumC(4?),gdh(-),pdhC(1,2)"
    );
    let known = Typing {
        sequence_type: Some("11".into()),
        ..t
    };
    assert!(line("g.fa", Some(&known), '\t').starts_with("g.fa\tneis\t11\t"));
    // no scheme chosen
    assert_eq!(line("g.fa", None, '\t'), "g.fa\t-\t-");
}

#[test]
fn the_full_table_has_the_status_and_the_score_and_the_calls_joined_by_semicolons() {
    assert_eq!(
        FULL_HEADER.join("\t"),
        "FILE\tSCHEME\tST\tSTATUS\tSCORE\tALLELES"
    );
    assert_eq!(
        full_line("g.fa", Some(&typing()), '\t'),
        "g.fa\tneis\t-\tOK\t72\tabcZ(3);adk(~7);fumC(4?);gdh(-);pdhC(1,2)"
    );
    assert_eq!(full_line("g.fa", None, ','), "g.fa,-,-,NONE,0,");
}

#[test]
fn json_is_an_array_with_an_object_for_each_genome_and_is_escaped() {
    let rows = vec![
        ("one".to_string(), "dir/one.fa".to_string(), Some(typing())),
        ("we\"ird\\name".to_string(), "x.fa".to_string(), None),
    ];
    let text = json(&rows);
    assert!(text.starts_with("[\n") && text.ends_with("\n]\n"));
    assert!(text.contains("\"id\": \"one\", \"filename\": \"dir/one.fa\", \"scheme\": \"neis\", \"sequence_type\": \"-\""));
    assert!(text.contains("\"alleles\": {\"abcZ\": \"3\", \"adk\": \"~7\", \"fumC\": \"4?\", \"gdh\": \"-\", \"pdhC\": \"1,2\"}"));
    assert!(text.contains("\"id\": \"we\\\"ird\\\\name\""));
    assert!(text.contains("\"scheme\": \"-\", \"sequence_type\": \"-\", \"alleles\": {}"));
    assert_eq!(text.matches("\"filename\"").count(), 2);
    assert_eq!(json_text("a\nb\t\u{1}"), "\"a\\nb\\t\\u0001\"");
}

#[test]
fn a_novel_allele_is_written_as_fasta_named_for_the_scheme_the_locus_and_the_md5_of_its_bases() {
    let t = typing();
    let want = format!(">neis.adk-{} g.fa\nACGT\n", md5::hex(b"ACGT"));
    assert_eq!(novel_fasta("g.fa", &t), want);
    assert_eq!(md5::hex(b"ACGT").len(), 32);
    // two novel places at one locus are two records, and none that are not novel
    let two = Typing {
        calls: vec![(
            "adk".into(),
            Call::Novel(vec![
                ("7".into(), b"AC".to_vec()),
                ("8".into(), b"GT".to_vec()),
            ]),
        )],
        ..typing()
    };
    assert_eq!(novel_fasta("g.fa", &two).matches('>').count(), 2);
    let none = Typing {
        calls: vec![("adk".into(), Call::Exact(vec!["7".into()]))],
        ..typing()
    };
    assert_eq!(novel_fasta("g.fa", &none), "");
}
