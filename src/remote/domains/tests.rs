use super::*;
use crate::pipeline::domains::testing::{pressed_db, stand_in, Folder, TABLE};

fn job() -> DomainsJob {
    DomainsJob {
        queries: vec![
            ("first protein".to_string(), vec![b'M'; 60]),
            ("second".to_string(), vec![b'G'; 30]),
        ],
        cutoff: Cutoff::Evalue(0.01),
    }
}

fn decoded(job: &DomainsJob) -> Result<DomainsJob, String> {
    let c = Container::decode(&job.encode().encode()).expect("a container");
    let m = meta(c.require("meta").expect("meta"));
    DomainsJob::decode(&c, &m)
}

/// The same sections with the text of `name` changed by `edit`.
fn edited(c: &Container, name: &str, edit: impl Fn(&str) -> String) -> Container {
    ["meta", "names", "sequences", "domains"]
        .into_iter()
        .filter_map(|n| c.text(n).map(|t| (n, t.to_string())))
        .fold(Container::new(), |c, (n, t)| {
            c.with(n, if n == name { edit(&t) } else { t })
        })
}

#[test]
fn the_job_goes_through_bytes_and_comes_back_as_it_was() {
    let mut job = job();
    assert_eq!(decoded(&job).expect("a job"), job);
    job.cutoff = Cutoff::Curated;
    assert_eq!(decoded(&job).expect("a job"), job);
}

#[test]
fn a_job_with_a_cutoff_or_sequences_that_make_no_sense_is_refused() {
    let c = job().encode();
    let refused = |name: &str, edit: &dyn Fn(&str) -> String, why: &str| {
        let c = edited(&c, name, edit);
        let m = meta(c.require("meta").expect("meta"));
        let said = DomainsJob::decode(&c, &m).expect_err(why);
        assert!(!said.is_empty(), "{why}");
    };
    refused(
        "meta",
        &|t| t.replace("evalue=0.01", "evalue=0"),
        "no E-value",
    );
    refused(
        "meta",
        &|t| t.replace("evalue=0.01", "evalue=1000"),
        "too large",
    );
    refused(
        "meta",
        &|t| t.replace("evalue=0.01", "evalue=NaN"),
        "not a number",
    );
    refused("meta", &|t| t.replace("evalue=0.01", "evalue=x"), "text");
    refused(
        "names",
        &|t| t.lines().next().unwrap_or("").to_string() + "\n",
        "one name short",
    );
    refused("sequences", &|_| String::new(), "no sequences");
    refused("sequences", &|t| t.replace('M', "m"), "lower case");
    refused("sequences", &|t| t.replace('M', "1"), "a digit");
    refused(
        "sequences",
        &|t| t.replace("G", "GG").replace("\n", "\n\n"),
        "a blank sequence",
    );
    refused(
        "sequences",
        &|_| format!("{}\n", "A".repeat(MAX_LENGTH + 1)),
        "too long",
    );
}

#[cfg(unix)]
#[test]
fn the_run_gives_the_domains_with_the_names_and_lengths_and_what_made_them() {
    let folder = Folder::new("remote-domains");
    let db = pressed_db(&folder);
    let program = stand_in(&folder, TABLE, "ok");
    let result = run_with(&job(), &program, &db).expect("a result");
    assert_eq!(result.names, ["first protein", "second"]);
    assert_eq!(result.lengths, [60, 30]);
    assert_eq!(result.hits.len(), 2);
    assert_eq!(result.program, "HMMER 3.4");
    assert_eq!(result.database, "Pfam-A.hmm");
    // and it goes through bytes unchanged
    let back = DomainsResult::decode(&result.encode()).expect("a result that decodes");
    assert_eq!(back, result);
}

#[cfg(unix)]
#[test]
fn a_damaged_result_is_refused_with_the_reason() {
    let folder = Folder::new("remote-domains-bad");
    let db = pressed_db(&folder);
    let result = run_with(&job(), &stand_in(&folder, TABLE, "ok"), &db).expect("a result");
    let c = result.encode();
    let damage = |name: &str, edit: &dyn Fn(&str) -> String| {
        DomainsResult::decode(&edited(&c, name, edit)).expect_err("damaged")
    };
    assert!(damage("names", &|t| t.replacen("\t60", "", 1)).contains("no length"));
    assert!(damage("names", &|t| t.replacen("\t60", "\t0", 1)).contains("out of range"));
    assert!(damage("names", &|t| t.replacen("\t60", "\tx", 1)).contains("out of range"));
    assert!(
        damage("domains", &|t| t.replacen("0\tPkinase", "9\tPkinase", 1)).contains("not there")
    );
    assert!(
        damage("domains", &|t| t.replacen("\t5\t55\t", "\t5\t61\t", 1)).contains("sequence of 60")
    );
    assert!(damage("domains", &|t| t.replacen("\t264\t", "\tx\t", 1)).contains("no count"));
    assert!(damage("domains", &|t| t.replacen("\t60\t", "\ty\t", 1)).contains("no number"));
    assert!(damage("domains", &|t| t
        .lines()
        .next()
        .unwrap_or("")
        .replacen('\t', "", 1)
        + "\n")
    .contains("wrong number"));
    assert!(damage("meta", &|t| t.replace("seconds=", "secs=")).contains("seconds"));
}

#[test]
fn a_computer_without_the_program_or_the_database_says_what_to_install() {
    // (the search finds hmmscan on the PATH and Pfam in the models folder; where either is missing the reason says what to do)
    match run(&job()) {
        Ok(result) => assert_eq!(result.names.len(), 2),
        Err(why) => assert!(
            why.contains("hmmscan was not found")
                || why.contains("Pfam database is not installed")
                || why.contains("hmmscan stopped"),
            "{why}"
        ),
    }
}
