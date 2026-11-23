use super::*;
use crate::seq::mlst::search::Found;
use crate::seq::mlst::tests::scheme;

/// A find of allele `id`: `covered` of `len` bases, `matches` of `columns` the same.
fn found(id: &str, exact: bool, covered: usize, matches: usize, columns: usize) -> Found {
    Found {
        allele: id.to_string(),
        exact: exact.then(|| id.to_string()),
        covered,
        allele_len: 100,
        matches,
        columns,
        sequence: vec![b'A'; covered],
    }
}

#[test]
fn a_locus_is_called_exact_then_novel_then_partial_then_missing() {
    let p = Params::default();
    // exact wins over everything else
    let all = [
        found("7", true, 100, 100, 100),
        found("3", false, 100, 99, 100),
        found("5", false, 50, 50, 50),
    ];
    assert_eq!(call_locus(&all, &p), Call::Exact(vec!["7".into()]));
    // novel: all of the allele, alike enough; the closest is named
    assert_eq!(
        call_locus(&[found("3", false, 100, 99, 100)], &p).text(),
        "~3"
    );
    // partial: some of it, alike enough
    assert_eq!(
        call_locus(&[found("5", false, 40, 40, 40)], &p).text(),
        "5?"
    );
    // nothing alike enough, or nothing
    assert_eq!(
        call_locus(&[found("5", false, 100, 80, 100)], &p),
        Call::Missing
    );
    assert_eq!(call_locus(&[], &p), Call::Missing);
}

#[test]
fn the_limits_of_identity_and_coverage_are_themselves_allowed() {
    let p = Params::default();
    // 95 % alike: novel; 94.9 % is not
    assert_eq!(
        call_locus(&[found("1", false, 100, 95, 100)], &p).text(),
        "~1"
    );
    assert_eq!(
        call_locus(&[found("1", false, 100, 949, 1000)], &p),
        Call::Missing
    );
    // 10 % of the allele: partial; 9 % is not
    assert_eq!(
        call_locus(&[found("2", false, 10, 10, 10)], &p).text(),
        "2?"
    );
    assert_eq!(call_locus(&[found("2", false, 9, 9, 9)], &p), Call::Missing);
    // a find that is all of the allele and under the identity does not become a partial one
    assert_eq!(
        call_locus(&[found("2", false, 100, 90, 100)], &p),
        Call::Missing
    );
    // the limits are the ones given
    let lenient = Params {
        min_identity: 80.0,
        min_coverage: 5.0,
        ..p
    };
    assert_eq!(
        call_locus(&[found("2", false, 100, 90, 100)], &lenient).text(),
        "~2"
    );
    assert_eq!(
        call_locus(&[found("2", false, 5, 5, 5)], &lenient).text(),
        "2?"
    );
}

#[test]
fn several_places_are_called_together_in_order_and_the_best_partial_is_named() {
    let p = Params::default();
    let two = [
        found("10", true, 100, 100, 100),
        found("9", true, 100, 100, 100),
        found("10", true, 100, 100, 100),
    ];
    // numbers in numeric order (9 before 10), and each once
    assert_eq!(call_locus(&two, &p).text(), "9,10");
    let letters = [
        found("B2", true, 100, 100, 100),
        found("A1", true, 100, 100, 100),
    ];
    assert_eq!(call_locus(&letters, &p).text(), "A1,B2");
    let novel = [
        found("3", false, 100, 98, 100),
        found("4", false, 100, 97, 100),
    ];
    assert_eq!(call_locus(&novel, &p).text(), "~3,~4");
    // of several partial finds the one that matches most
    let partial = [
        found("1", false, 30, 30, 30),
        found("2", false, 60, 60, 60),
        found("3", false, 50, 50, 50),
    ];
    assert_eq!(call_locus(&partial, &p).text(), "2?");
}

fn typed(calls: Vec<Call>) -> Typing {
    let s = scheme("synth", calls.len(), 6);
    summarize(
        &s,
        calls
            .into_iter()
            .enumerate()
            .map(|(l, c)| (format!("L{l}"), c))
            .collect(),
    )
}

fn ex(id: usize) -> Call {
    Call::Exact(vec![id.to_string()])
}

#[test]
fn the_score_is_ninety_sixty_three_and_eighteen_over_the_loci_and_ten_for_a_known_type() {
    // type 1 of this scheme: allele (1 + l) % 6 + 1: 3, 4, 5, 6, 1, 2, 3
    let known = typed((0..7).map(|l| ex((1 + l) % 6 + 1)).collect());
    assert_eq!(
        (known.sequence_type.as_deref(), known.status),
        (Some("1"), Status::Perfect)
    );
    assert!((known.score - 100.0).abs() < 1e-9);
    let unknown = typed(std::iter::repeat_n(ex(1), 7).collect());
    assert!((unknown.score - 90.0).abs() < 1e-9);
    assert_eq!(
        (unknown.sequence_type, unknown.status),
        (None, Status::Novel)
    );
    let mut one_novel: Vec<Call> = std::iter::repeat_n(ex(1), 6).collect();
    one_novel.push(Call::Novel(vec![("2".into(), b"ACGT".to_vec())]));
    let t = typed(one_novel);
    assert!((t.score - (6.0 * 90.0 + 63.0) / 7.0).abs() < 1e-9);
    assert_eq!(t.status, Status::Ok);
    let mut partial: Vec<Call> = std::iter::repeat_n(ex(1), 4).collect();
    partial.extend(std::iter::repeat_n(Call::Partial("1".into()), 3));
    let t = typed(partial);
    assert!((t.score - (4.0 * 90.0 + 3.0 * 18.0) / 7.0).abs() < 1e-9);
    // under 70, and none of the other statuses: bad
    assert_eq!(t.status, Status::Bad);
}

#[test]
fn the_status_follows_the_calls_in_the_order_the_readme_gives() {
    assert_eq!(
        typed(std::iter::repeat_n(Call::Missing, 7).collect()).status,
        Status::None
    );
    assert_eq!(
        typed(std::iter::repeat_n(Call::Missing, 7).collect()).score,
        0.0
    );
    let mut mixed: Vec<Call> = std::iter::repeat_n(ex(1), 6).collect();
    mixed.push(Call::Exact(vec!["1".into(), "2".into()]));
    let t = typed(mixed);
    assert_eq!((t.status, t.sequence_type.clone()), (Status::Mixed, None));
    // several alleles still score as exact
    assert!((t.score - 90.0).abs() < 1e-9);
    let mut missing: Vec<Call> = std::iter::repeat_n(ex(1), 6).collect();
    missing.push(Call::Missing);
    assert_eq!(typed(missing).status, Status::Missing);
    // the words
    let words: Vec<&str> = [
        Status::Perfect,
        Status::Novel,
        Status::None,
        Status::Mixed,
        Status::Missing,
        Status::Bad,
        Status::Ok,
    ]
    .iter()
    .map(|s| s.word())
    .collect();
    assert_eq!(
        words,
        ["PERFECT", "NOVEL", "NONE", "MIXED", "MISSING", "BAD", "OK"]
    );
}

#[test]
fn a_call_is_written_as_the_readme_writes_it() {
    assert_eq!(Call::Exact(vec!["12".into()]).text(), "12");
    assert_eq!(Call::Exact(vec!["1".into(), "2".into()]).text(), "1,2");
    assert_eq!(Call::Novel(vec![("5".into(), Vec::new())]).text(), "~5");
    assert_eq!(Call::Partial("5".into()).text(), "5?");
    assert_eq!(Call::Missing.text(), "-");
    assert_eq!(id_order("9", "10"), std::cmp::Ordering::Less);
    assert_eq!(id_order("A", "B"), std::cmp::Ordering::Less);
}

#[test]
fn a_novel_allele_needs_every_base_of_the_allele_one_short_is_partial() {
    let p = Params::default();
    // all 100 bases, all the same but one: novel; 99 bases, every one the same: partial
    assert_eq!(
        call_locus(&[found("3", false, 100, 99, 100)], &p).text(),
        "~3"
    );
    assert_eq!(
        call_locus(&[found("3", false, 99, 99, 99)], &p).text(),
        "3?"
    );
}

#[test]
fn a_locus_with_two_alleles_gives_no_type_even_when_the_first_of_them_would() {
    // type 1 of this scheme has allele (1 + l) % 6 + 1 at locus l; here locus 0 has that allele and another
    let mut calls: Vec<Call> = (0..7).map(|l| ex((1 + l) % 6 + 1)).collect();
    calls[0] = Call::Exact(vec!["2".into(), "5".into()]);
    let t = typed(calls);
    assert_eq!((t.sequence_type, t.status), (None, Status::Mixed));
    // each of the seven loci scores as exact, and nothing is added for a type
    assert_eq!(t.score, 90.0);
}

#[test]
fn a_score_of_exactly_seventy_is_not_bad_and_a_whole_number_is_a_whole_number() {
    // 4 exact, 4 novel and 1 partial of 9 loci: (360 + 252 + 18) / 9 = 70
    let novel = || Call::Novel(vec![("2".into(), b"ACGT".to_vec())]);
    let mut calls: Vec<Call> = std::iter::repeat_n(ex(1), 4).collect();
    calls.extend(std::iter::repeat_with(novel).take(4));
    calls.push(Call::Partial("1".into()));
    let t = typed(calls);
    assert_eq!(t.score, 70.0);
    assert_eq!(t.status, Status::Ok);
    // one exact less and one novel more: 67
    let mut calls: Vec<Call> = std::iter::repeat_n(ex(1), 3).collect();
    calls.extend(std::iter::repeat_with(novel).take(5));
    calls.push(Call::Partial("1".into()));
    let t = typed(calls);
    assert_eq!(t.score, 67.0);
    assert_eq!(t.status, Status::Bad);
    // a known type of seven loci is 100, not a hair over or under
    let known = typed((0..7).map(|l| ex((1 + l) % 6 + 1)).collect());
    assert_eq!(known.score, 100.0);
}
