use super::*;

const WORD: &[u8] = b"ACCGTTAGCA";

/// `count` sequences of 120 random bases with WORD in each (forward strand), the same on every run.
fn planted(count: usize) -> Vec<Vec<u8>> {
    let mut state: u64 = 0x1234_5678_9ABC_DEF1;
    let mut next = move |n: usize| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 33) % n as u64) as usize
    };
    (0..count)
        .map(|_| {
            let mut seq: Vec<u8> = (0..120).map(|_| b"ACGT"[next(4)]).collect();
            let at = next(120 - WORD.len() + 1);
            seq[at..at + WORD.len()].copy_from_slice(WORD);
            seq
        })
        .collect()
}

fn job() -> MotifJob {
    MotifJob {
        seqs: planted(12),
        settings: Settings {
            width: 10,
            motifs: 2,
            both_strands: false,
            seeds: 7,
            iterations: 30,
            effort: 2.5,
        },
    }
}

fn decoded(job: &MotifJob) -> Result<MotifJob, String> {
    let c = Container::decode(&job.encode().encode()).expect("a container");
    let m = meta(c.require("meta").expect("meta"));
    MotifJob::decode(&c, &m)
}

/// The same sections with the text of `name` changed by `edit`.
fn edited(c: &Container, name: &str, edit: impl Fn(&str) -> String) -> Container {
    ["meta", "sequences", "motifs"]
        .into_iter()
        .filter_map(|n| c.text(n).map(|t| (n, t.to_string())))
        .fold(Container::new(), |c, (n, t)| {
            c.with(n, if n == name { edit(&t) } else { t })
        })
}

#[test]
fn the_job_goes_through_bytes_and_comes_back_as_it_was() {
    let job = job();
    assert_eq!(decoded(&job).expect("a job"), job);
    let default = MotifJob {
        settings: Settings::default(),
        ..job
    };
    assert_eq!(decoded(&default).expect("a job"), default);
}

#[test]
fn a_job_with_settings_out_of_range_or_what_are_no_sequences_is_refused() {
    let c = job().encode();
    let refused = |name: &str, edit: &dyn Fn(&str) -> String, why: &str| {
        let c = edited(&c, name, edit);
        let m = meta(c.require("meta").expect("meta"));
        let said = MotifJob::decode(&c, &m).expect_err(why);
        assert!(!said.is_empty(), "{why}");
        said
    };
    let set = |key: &'static str, value: &'static str| {
        move |t: &str| {
            t.lines()
                .map(|l| {
                    if l.starts_with(&format!("{key}=")) {
                        format!("{key}={value}\n")
                    } else {
                        format!("{l}\n")
                    }
                })
                .collect::<String>()
        }
    };
    for (key, value) in [
        ("width", "2"),
        ("width", "41"),
        ("motifs", "0"),
        ("motifs", "51"),
        ("seeds", "0"),
        ("seeds", "1001"),
        ("iterations", "0"),
        ("iterations", "10001"),
        ("effort", "0"),
        ("effort", "5000"),
        ("effort", "NaN"),
        ("width", "wide"),
    ] {
        refused("meta", &set(key, value), &format!("{key}={value}"));
    }
    // a character that is no letter in a sequence, and too many sequences
    let said = refused("sequences", &|t| t.replace('A', "A1"), "digits");
    assert!(said.contains("no letter"), "{said}");
    refused(
        "sequences",
        &|_| "ACGTACGTAC\n".repeat(MAX_SEQUENCES + 1),
        "too many",
    );
    // windows that would not fit in memory are refused before any is made
    let big = edited(&c, "sequences", |_| format!("{}\n", "A".repeat(30_000_000))).with(
        "meta",
        c.text("meta")
            .expect("meta")
            .replace("width=10", "width=40"),
    );
    let m = meta(big.require("meta").expect("meta"));
    assert!(MotifJob::decode(&big, &m)
        .unwrap_err()
        .contains("use fewer bases"));
    // the settings are all needed
    let without = edited(&c, "meta", |t| {
        t.lines()
            .filter(|l| !l.starts_with("effort="))
            .map(|l| format!("{l}\n"))
            .collect()
    });
    let m = meta(without.require("meta").expect("meta"));
    assert!(MotifJob::decode(&without, &m).is_err());
}

#[test]
fn the_run_gives_what_the_search_gives_and_the_result_goes_through_bytes_unchanged() {
    let job = job();
    let result = run(&job).expect("a result");
    let direct = discover(&job.seqs, &job.settings).expect("a search");
    assert_eq!(result.found, direct);
    assert_eq!(result.found.len(), 2);
    // the planted word is the first motif, found in nearly every sequence
    assert!(result.found[0].sites.len() >= 10);
    let back = MotifResult::decode(&result.encode()).expect("a result that decodes");
    assert_eq!(back, result, "the floats come back to the last bit");
    result.fits(&job.seqs).expect("it is about these sequences");
}

#[test]
fn a_result_is_about_the_sequences_that_were_sent_and_no_others() {
    let job = job();
    let result = run(&job).expect("a result");
    // fewer sequences: a site is in one that is not there
    assert!(result.fits(&job.seqs[..1]).is_err());
    // shorter sequences: the window runs off the end
    let short: Vec<Vec<u8>> = job.seqs.iter().map(|s| s[..30].to_vec()).collect();
    assert!(result.fits(&short).is_err());
    // a window must end inside its sequence: cut a sequence exactly at the end of its last window, and one base before it
    let k = result.found[0].sites[0].sequence;
    let end = result
        .found
        .iter()
        .flat_map(|f| f.sites.iter().map(|s| (s, f.motif.width())))
        .filter(|(s, _)| s.sequence == k)
        .map(|(s, width)| s.position + width)
        .max()
        .expect("a window");
    let cut = |len: usize| {
        let mut seqs = job.seqs.clone();
        seqs[k].truncate(len);
        seqs
    };
    assert!(result.fits(&cut(end)).is_ok());
    assert!(result.fits(&cut(end - 1)).is_err());
    assert!(MotifResult {
        found: Vec::new(),
        seconds: 0.0
    }
    .fits(&[])
    .is_ok());
}

#[test]
fn a_damaged_result_is_refused_with_the_reason() {
    let result = run(&job()).expect("a result");
    let c = result.encode();
    assert!(MotifResult::decode(&c).is_ok());
    let damage = |name: &str, edit: &dyn Fn(&str) -> String| {
        MotifResult::decode(&edited(&c, name, edit)).expect_err("damaged")
    };
    // the number of motifs it says, against what it has
    let said = damage("meta", &|t| t.replace("motifs=2", "motifs=1"));
    assert!(said.contains("more than it says"), "{said}");
    let said = damage("meta", &|t| t.replace("motifs=2", "motifs=3"));
    assert!(said.contains("fewer than it says"), "{said}");
    assert!(damage("meta", &|t| t.replace("motifs=2", "motifs=999")).contains("too many"));
    // a heading, a column, a site
    assert!(damage("motifs", &|t| t.replacen("motif\t", "motive\t", 1)).contains("no heading"));
    assert!(damage("motifs", &|t| t.replacen("\t10\t", "\t99\t", 1)).contains("out of range"));
    let lines = |t: &str| t.lines().map(str::to_string).collect::<Vec<String>>();
    let with_line = |n: usize, new: &'static str| {
        move |t: &str| {
            lines(t)
                .into_iter()
                .enumerate()
                .map(|(k, l)| format!("{}\n", if k == n { new.to_string() } else { l }))
                .collect::<String>()
        }
    };
    // line 0 is the heading; lines 1 to 10 are the columns; then the sites
    assert!(damage("motifs", &with_line(1, "0.5\t0.5\t0.5\t0.5")).contains("no distribution"));
    assert!(damage("motifs", &with_line(1, "0.5\t0.5\t0")).contains("not four numbers"));
    assert!(damage("motifs", &with_line(1, "a\tb\tc\td")).contains("probability"));
    assert!(damage("motifs", &with_line(1, "-0.1\t0.6\t0.3\t0.2")).contains("no distribution"));
    assert!(damage("motifs", &with_line(11, "0\t1\t2")).contains("wrong fields"));
    assert!(damage("motifs", &with_line(11, "0\t1\t2\t3.0\t0.5")).contains("strand"));
    assert!(damage("motifs", &with_line(11, "0\t1\t0\tNaN\t0.5")).contains("out of range"));
    assert!(damage("motifs", &with_line(11, "0\t1\t0\t3.0\t1.5")).contains("out of range"));
    assert!(damage("motifs", &with_line(11, "x\t1\t0\t3.0\t0.5")).contains("sequence number"));
    // cut short
    let whole = c.text("motifs").expect("motifs").to_string();
    for keep in [1, 5, 12] {
        let short = whole.lines().take(keep).collect::<Vec<_>>().join("\n");
        assert!(
            MotifResult::decode(&edited(&c, "motifs", |_| format!("{short}\n"))).is_err(),
            "{keep} lines"
        );
    }
}
