use super::*;

fn feature(name: &str, kind: &str, from: usize, to: usize) -> FeatureInfo {
    FeatureInfo {
        name: name.into(),
        kind: kind.into(),
        from,
        to,
    }
}

struct Held {
    sequences: Vec<String>,
    features: Vec<FeatureInfo>,
    recent: Vec<String>,
}

impl Held {
    fn new() -> Held {
        Held {
            sequences: vec!["pUC19".into(), "pET28a".into(), "reads.fastq".into()],
            features: vec![
                feature("lacZ alpha", "CDS", 146, 469),
                feature("AmpR", "CDS", 1626, 2486),
                feature("ori", "rep_origin", 867, 1455),
                feature("T7 promoter", "promoter", 20, 38),
            ],
            recent: vec![
                "pBR322.gb".into(),
                "ladder.fasta".into(),
                "pUC19-copy.gb".into(),
            ],
        }
    }

    fn work(&self, current: Option<usize>, bases: usize) -> Work<'_> {
        Work {
            sequences: &self.sequences,
            current,
            bases,
            features: &self.features,
            recent: &self.recent,
        }
    }
}

fn titles(entries: &[Entry]) -> Vec<&str> {
    entries.iter().map(|e| e.title.as_str()).collect()
}

fn jump(query: &str, held: &Held) -> Option<usize> {
    answers(query, &held.work(Some(0), 2686))
        .lead
        .into_iter()
        .find_map(|e| match e.action {
            Action::JumpTo(n) => Some(n),
            _ => None,
        })
}

#[test]
fn a_position_typed_in_any_of_its_usual_forms_is_a_jump_to_that_base() {
    let held = Held::new();
    for query in [
        "go to 1200",
        "goto 1200",
        "jump 1200",
        "jump to 1200",
        "base 1200",
        "pos 1200",
        "position 1200",
        "  Go To   1200  ",
        "1200",
        "GOTO1200",
    ] {
        assert_eq!(jump(query, &held), Some(1200), "{query:?}");
    }
    let a = answers("go to 1200", &held.work(Some(0), 2686));
    assert_eq!(titles(&a.lead), ["Jump to base 1200"]);
    assert_eq!(a.lead[0].group, "Sequence");
    assert!(a.trail.is_empty(), "a position is not a name");
    // the last base is a base, not the end
    let last = answers("go to 2686", &held.work(Some(0), 2686));
    assert_eq!(titles(&last.lead), ["Jump to base 2686"]);
    // past the end it goes to the end, and says so
    let past = answers("go to 99999", &held.work(Some(0), 2686));
    assert_eq!(
        titles(&past.lead),
        ["Jump to the end of the sequence (2686 bases)"]
    );
    assert!(matches!(past.lead[0].action, Action::JumpTo(2686)));
    // not a position: no number, a zero, or something after it
    for query in [
        "go to",
        "go to abc",
        "go to 0",
        "go to 12x",
        "jump to 1 2",
        "base",
        "",
    ] {
        assert_eq!(jump(query, &held), None, "{query:?}");
    }
    // with no sequence open, or an empty one, there is nowhere to go
    assert!(answers("go to 5", &held.work(None, 0)).lead.is_empty());
    assert!(answers("go to 5", &held.work(Some(0), 0)).lead.is_empty());
}

#[test]
fn find_and_a_pattern_is_a_search_of_the_sequence_in_view_with_the_pattern_as_typed() {
    let held = Held::new();
    let a = answers("find GAATTC", &held.work(Some(0), 2686));
    assert_eq!(titles(&a.lead), ["Find GAATTC in the sequence"]);
    assert!(matches!(&a.lead[0].action, Action::FindText(p) if p == "GAATTC"));
    assert_eq!(a.lead[0].group, "Sequence");
    // the case of the pattern is kept (the Find window decides whether it matters), spaces round it are not
    let b = answers("  Search   ggatCC ", &held.work(Some(0), 2686));
    assert!(
        matches!(&b.lead[0].action, Action::FindText(p) if p == "ggatCC"),
        "{:?}",
        b.lead
    );
    // a name is not searched for when a pattern is
    assert!(a.trail.is_empty());
    assert!(answers("find", &held.work(Some(0), 10)).lead.is_empty());
    assert!(answers("find   ", &held.work(Some(0), 10)).lead.is_empty());
    assert!(answers("find GAATTC", &held.work(None, 0)).lead.is_empty());
    let long = format!("find {}", "A".repeat(MOST_PATTERN + 1));
    assert!(answers(&long, &held.work(Some(0), 10)).lead.is_empty());
}

#[test]
fn the_names_of_the_open_sequences_the_features_and_the_recent_files_answer_a_few_letters() {
    let held = Held::new();
    let a = answers("puc", &held.work(Some(1), 5000));
    assert_eq!(titles(&a.trail), ["Switch to pUC19", "Open pUC19-copy.gb"]);
    assert!(matches!(a.trail[0].action, Action::SwitchDoc(0)));
    assert!(matches!(a.trail[1].action, Action::OpenRecent(2)));
    assert_eq!(a.trail[0].group, "Open sequences");
    assert_eq!(a.trail[1].group, "Recent files");
    assert!(a.lead.is_empty());
    // the sequence in view is not offered: it is there already
    let a = answers("puc", &held.work(Some(0), 5000));
    assert_eq!(titles(&a.trail), ["Open pUC19-copy.gb"]);
    // a feature by its name or its kind
    let a = answers("lac", &held.work(Some(0), 2686));
    assert_eq!(titles(&a.trail), ["lacZ alpha — CDS 146–469"]);
    assert!(matches!(a.trail[0].action, Action::SelectFeature(0)));
    assert_eq!(a.trail[0].group, "Features");
    let a = answers("cds", &held.work(Some(0), 2686));
    assert_eq!(a.trail.len(), 2, "{:?}", titles(&a.trail));
    // the case does not matter
    assert_eq!(answers("AMPR", &held.work(Some(0), 2686)).trail.len(), 1);
    // one letter says too little; nothing typed says nothing
    assert!(answers("p", &held.work(Some(0), 2686)).trail.is_empty());
    assert!(answers("", &held.work(Some(0), 2686)).trail.is_empty());
    assert!(answers("   ", &held.work(Some(0), 2686)).trail.is_empty());
    // with no sequence in view there are no features to offer, but the open ones and the recent files still are
    let a = answers("pet", &held.work(None, 0));
    assert_eq!(titles(&a.trail), ["Switch to pET28a"]);
    assert!(
        answers("lac", &held.work(None, 0)).trail.is_empty(),
        "features belong to a sequence in view"
    );
}

#[test]
fn each_kind_of_answer_is_limited_so_that_a_common_word_does_not_fill_the_list() {
    let mut held = Held::new();
    held.features = (0..40)
        .map(|k| feature(&format!("gene {k}"), "CDS", k, k + 1))
        .collect();
    held.recent = (0..20).map(|k| format!("gene-{k}.gb")).collect();
    held.sequences = (0..20).map(|k| format!("gene sequence {k}")).collect();
    let a = answers("gene", &held.work(Some(0), 100));
    let count = |group: &str| a.trail.iter().filter(|e| e.group == group).count();
    assert_eq!(count("Features"), MOST_FEATURES);
    assert_eq!(count("Recent files"), MOST_RECENT);
    assert_eq!(count("Open sequences"), MOST_SEQUENCES);
}
