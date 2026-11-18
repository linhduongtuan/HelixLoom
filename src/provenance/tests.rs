use super::*;

fn record() -> Provenance {
    Provenance {
        analysis: "Alignment".into(),
        inputs: vec![
            Input::of("seqA", b"ACGTACGT"),
            Input::of("seqB", b"ACGAACGT"),
        ],
        parameters: vec![
            ("engine".into(), "MAFFT".into()),
            ("bootstrap".into(), "0".into()),
        ],
        programs: vec!["MAFFT v7.526".into(), "HelixLoom 0.2.0".into()],
        ran_on: "this computer".into(),
        date: "2026-11-13".into(),
        seconds: 0.4217,
    }
}

#[test]
fn an_input_is_its_name_its_size_and_a_checksum_of_its_bytes() {
    let a = Input::of("seqA", b"ACGTACGT");
    assert_eq!((a.name.as_str(), a.size), ("seqA", 8));
    assert_eq!(a.digest.len(), 32);
    assert_eq!(a, Input::of("seqA", b"ACGTACGT"));
    assert_ne!(
        a.digest,
        Input::of("seqA", b"ACGTACGA").digest,
        "one letter changes it"
    );
    assert_eq!(
        a.digest,
        Input::of("other name", b"ACGTACGT").digest,
        "the name is not in it"
    );
    assert_eq!(a.short(), &a.digest[..12]);
}

#[test]
fn a_record_is_written_as_text_that_says_everything_it_holds() {
    let text = record().to_text();
    let a = Input::of("seqA", b"ACGTACGT");
    let b = Input::of("seqB", b"ACGAACGT");
    assert_eq!(
        text,
        format!(
            "Alignment\n\
             Made on 2026-11-13 on this computer in 0.42 s\n\
             Programs: MAFFT v7.526; HelixLoom 0.2.0\n\
             Settings: engine=MAFFT; bootstrap=0\n\
             Inputs (2):\n  seqA\t8\t{}\n  seqB\t8\t{}\n",
            a.digest, b.digest
        )
    );
    // nothing to say is said, not left out
    let bare = Provenance {
        parameters: Vec::new(),
        inputs: Vec::new(),
        ..record()
    };
    let text = bare.to_text();
    assert!(
        text.contains("Settings: none\n") && text.contains("Inputs (0):\n"),
        "{text}"
    );
    // a name that would break the lines is flattened
    let odd = Provenance {
        inputs: vec![Input::of("a\tb\nc", b"AC")],
        ..record()
    };
    assert!(odd.to_text().contains("  a b c\t2\t"), "{}", odd.to_text());
}

#[test]
fn nothing_has_changed_when_the_inputs_and_settings_are_the_ones_recorded() {
    let r = record();
    assert!(r.changed_since(&r.inputs, &r.parameters).is_empty());
}

#[test]
fn what_changed_is_said_for_inputs_that_changed_went_or_came_and_for_their_order() {
    let r = record();
    let same = r.parameters.clone();
    let (a, b) = (r.inputs[0].clone(), r.inputs[1].clone());
    let edited = Input::of("seqB", b"ACGAACGG");
    assert_eq!(
        r.changed_since(&[a.clone(), edited], &same),
        ["seqB has changed"]
    );
    assert_eq!(
        r.changed_since(std::slice::from_ref(&a), &same),
        ["seqB is no longer there"]
    );
    let extra = Input::of("seqC", b"TTTT");
    assert_eq!(
        r.changed_since(&[a.clone(), b.clone(), extra], &same),
        ["seqC is new"]
    );
    assert_eq!(
        r.changed_since(&[b.clone(), a.clone()], &same),
        ["the order of the inputs has changed"]
    );
    // a renamed input is one that went and one that came
    let renamed = Input {
        name: "seqZ".into(),
        ..b.clone()
    };
    assert_eq!(
        r.changed_since(&[a.clone(), renamed], &same),
        ["seqB is no longer there", "seqZ is new"]
    );
    // two inputs of one name are told apart by what they hold
    let twice = Provenance {
        inputs: vec![Input::of("x", b"AAAA"), Input::of("x", b"CCCC")],
        ..record()
    };
    assert!(twice.changed_since(&twice.inputs, &same).is_empty());
    assert_eq!(
        twice.changed_since(&[Input::of("x", b"AAAA"), Input::of("x", b"GGGG")], &same),
        ["x has changed"]
    );
}

#[test]
fn what_changed_is_said_for_settings_that_changed_went_or_came() {
    let r = record();
    let p = |pairs: &[(&str, &str)]| -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|(k, v)| ((*k).to_string(), (*v).to_string()))
            .collect()
    };
    assert_eq!(
        r.changed_since(&r.inputs, &p(&[("engine", "built in"), ("bootstrap", "0")])),
        ["engine was MAFFT, is now built in"]
    );
    assert_eq!(
        r.changed_since(&r.inputs, &p(&[("engine", "MAFFT")])),
        ["bootstrap is no longer a setting"]
    );
    assert_eq!(
        r.changed_since(
            &r.inputs,
            &p(&[("engine", "MAFFT"), ("bootstrap", "0"), ("seed", "7")])
        ),
        ["seed is new (7)"]
    );
    // the order of the settings does not matter
    assert!(r
        .changed_since(&r.inputs, &p(&[("bootstrap", "0"), ("engine", "MAFFT")]))
        .is_empty());
}

#[test]
fn a_long_list_of_changes_is_cut_and_says_how_many_more_there_are() {
    let many: Vec<Input> = (0..30)
        .map(|k| Input::of(&format!("s{k}"), b"AC"))
        .collect();
    let r = Provenance {
        inputs: many,
        ..record()
    };
    let changes = r.changed_since(&[], &r.parameters);
    assert_eq!(changes.len(), MOST_CHANGES + 1);
    assert_eq!(changes[0], "s0 is no longer there");
    assert_eq!(
        changes[MOST_CHANGES],
        format!("and {} more", 30 - MOST_CHANGES)
    );
}

#[test]
fn what_was_asked_becomes_the_record_when_the_answer_comes_with_this_program_and_the_day() {
    let asked = Asked {
        analysis: "Tree estimation",
        inputs: vec![Input::of("row", b"AC-T")],
        parameters: vec![("method".into(), "neighbour-joining".into())],
        ran_on: "the workstation",
    };
    let made = asked.clone().made(vec!["IQ-TREE 3.1.4".into()], 2.5);
    assert_eq!(made.analysis, "Tree estimation");
    assert_eq!(made.inputs, asked.inputs);
    assert_eq!(made.parameters, asked.parameters);
    assert_eq!(made.ran_on, "the workstation");
    assert_eq!(made.seconds, 2.5);
    assert_eq!(
        made.programs,
        [
            "IQ-TREE 3.1.4".to_string(),
            format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))
        ]
    );
    // the day is a date of the Gregorian calendar, of today
    let date = &made.date;
    assert_eq!(date.len(), 10);
    assert!(
        date.as_bytes()[4] == b'-' && date.as_bytes()[7] == b'-',
        "{date}"
    );
    assert!(date.as_str() >= "2026-11-13", "{date}");
    // an engine that is this program itself is not named twice
    let own = asked.made(
        vec![format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))],
        0.1,
    );
    assert_eq!(own.programs.len(), 1);
}

#[test]
fn a_record_that_is_kept_comes_back_as_it_was_whatever_its_words_hold() {
    let plain = record();
    assert_eq!(Provenance::from_kept(&plain.kept()), Some(plain));
    // names and settings are the person's: a tab, a line end, a backslash or nothing at all in them comes back as it was
    let awkward = Provenance {
        analysis: "An analysis\twith a tab".into(),
        inputs: vec![
            Input::of("a name\nover two lines", b"ACGT"),
            Input::of("", b""),
            Input::of("back\\slash", b"TTTT"),
        ],
        parameters: vec![
            ("a key".into(), "a value\twith a tab".into()),
            ("empty".into(), String::new()),
        ],
        programs: vec!["a program 1.0\r".into(), "HelixLoom 0.2.0".into()],
        ran_on: "the workstation".into(),
        date: "2026-11-18".into(),
        seconds: 0.123_456_789_012_345,
    };
    assert_eq!(Provenance::from_kept(&awkward.kept()), Some(awkward));
    // a record of nothing but its name
    let bare = Provenance {
        analysis: "Bare".into(),
        inputs: Vec::new(),
        parameters: Vec::new(),
        programs: Vec::new(),
        ran_on: String::new(),
        date: String::new(),
        seconds: 0.0,
    };
    assert_eq!(Provenance::from_kept(&bare.kept()), Some(bare));
}

#[test]
fn a_kept_record_that_cannot_be_read_is_no_record_and_a_line_it_does_not_know_is_passed_over() {
    let kept = record().kept();
    // a line from a later version of the program
    let later = format!("{kept}something new\tto come\n");
    assert_eq!(Provenance::from_kept(&later), Some(record()));
    // nothing, or no name of the analysis, or a time or a size that is no number: no record, not half of one
    assert_eq!(Provenance::from_kept(""), None);
    let nameless: String = kept
        .lines()
        .filter(|l| !l.starts_with("analysis\t"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(Provenance::from_kept(&nameless), None);
    let seconds = kept
        .lines()
        .find(|l| l.starts_with("seconds\t"))
        .expect("the time");
    assert_eq!(
        Provenance::from_kept(&kept.replace(seconds, "seconds\tsoon")),
        None
    );
    let input = kept
        .lines()
        .find(|l| l.starts_with("input\t"))
        .expect("an input");
    let mut fields: Vec<&str> = input.split('\t').collect();
    fields[2] = "many";
    assert_eq!(
        Provenance::from_kept(&kept.replace(input, &fields.join("\t"))),
        None
    );
    // a time that is not a time (it would be shown, and added up)
    assert_eq!(
        Provenance::from_kept(&kept.replace(seconds, "seconds\tNaN")),
        None
    );
    assert_eq!(
        Provenance::from_kept(&kept.replace(seconds, "seconds\t-1")),
        None
    );
}
