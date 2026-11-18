use super::*;
use crate::pipeline::flow::tests::{loci, locus};
use std::path::PathBuf;

/// A folder of its own for one test, removed when the test ends.
struct Folder(PathBuf);

impl Folder {
    fn new() -> Folder {
        Folder(crate::pipeline::tools::work_folder("flowjob").expect("a folder"))
    }

    fn cache(&self) -> PathBuf {
        self.0.join("cache")
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn job(n: u64) -> FlowJob {
    FlowJob {
        pipeline: Pipeline::loci_to_species(),
        loci: loci(n),
    }
}

fn through_bytes(c: &Container) -> Container {
    Container::decode(&c.encode()).expect("a container")
}

/// The sections of a request, to be put together again with one of them changed.
fn request_with(good: &Container, name: &str, new: &str) -> Container {
    ["meta", "pipeline", "names", "sizes", "loci"]
        .into_iter()
        .fold(Container::new(), |c, section| {
            let text = if section == name {
                new
            } else {
                good.text(section).expect("a section")
            };
            c.with(section, text)
        })
}

fn result_with(good: &Container, name: &str, new: &str) -> Container {
    [
        "meta",
        "pipeline",
        "report",
        "out_names",
        "out_sizes",
        "outputs",
        "output",
    ]
    .into_iter()
    .fold(Container::new(), |c, section| {
        match (section == name, good.text(section)) {
            (true, _) => c.with(section, new),
            (false, Some(text)) => c.with(section, text),
            (false, None) => c,
        }
    })
}

// --- the request ------------------------------------------------------------------------------------------------------------

#[test]
fn the_job_goes_through_bytes_and_comes_back_as_it_was() {
    for pipeline in [
        Pipeline::loci_to_species(),
        Pipeline::unaligned_loci_to_species(),
    ] {
        let loci = match pipeline.input {
            Kind::Sequences => (1..=3)
                .map(|k| (format!("locus{k}"), Value::Sequences(locus(k * 5, 90))))
                .collect(),
            _ => loci(3),
        };
        let sent = FlowJob { pipeline, loci };
        let back = FlowJob::decode(&through_bytes(&sent.encode())).expect("a job");
        assert_eq!(back, sent);
    }
}

#[test]
fn a_request_that_does_not_add_up_is_refused_and_says_what_is_wrong() {
    let good = job(2).encode();
    let (names, sizes, body) = (
        good.text("names").unwrap().to_string(),
        good.text("sizes").unwrap().to_string(),
        good.text("loci").unwrap().to_string(),
    );
    let first_size = sizes.lines().next().unwrap().to_string();
    let refused = |c: Container, wanted: &str| {
        let why = FlowJob::decode(&c).expect_err(wanted);
        assert!(why.contains(wanted), "{wanted:?} is not in {why:?}");
    };
    // sizes that do not cut the files out exactly
    refused(
        request_with(&good, "sizes", &first_size),
        "not as many sizes",
    );
    refused(
        request_with(&good, "sizes", &format!("{first_size}\n99999999")),
        "past the end",
    );
    refused(
        request_with(&good, "sizes", &format!("{first_size}\n1")),
        "left over",
    );
    refused(
        request_with(&good, "sizes", &format!("{first_size}\nx")),
        "not a number",
    );
    refused(
        request_with(&good, "sizes", &format!("{first_size}\n{}", usize::MAX)),
        "too large",
    );
    // a size that ends inside a character
    let wide = request_with(&good, "loci", &format!("{}é", &body[..body.len() - 1]));
    let cut = format!(
        "{first_size}\n{}",
        body.len() - first_size.parse::<usize>().unwrap()
    );
    refused(request_with(&wide, "sizes", &cut), "cuts a character");
    // names
    refused(request_with(&good, "names", "locus1\nlocus1"), "differ");
    refused(request_with(&good, "names", "locus1\n"), "not as many");
    refused(
        request_with(&good, "names", &format!("locus1\n\n{names}")),
        "not as many",
    );
    // nothing to run on
    let empty = request_with(&good, "names", "");
    let empty = request_with(&empty, "sizes", "");
    refused(request_with(&empty, "loci", ""), "no files");
    // a pipeline that is not one, or is far too long
    refused(
        request_with(
            &good,
            "pipeline",
            "helixloom-pipeline 1\ninput: alignment\nstep blast\n",
        ),
        "no step called blast",
    );
    refused(
        request_with(&good, "pipeline", &"#".repeat(MAX_PIPELINE_BYTES + 1)),
        "far longer",
    );
    refused(request_with(&good, "pipeline", ""), "empty");
    // a file that is not what the pipeline takes: rows of different lengths
    let ragged = ">a\nACGT\n>b\nAC\n";
    let one = request_with(&good, "names", "ragged");
    let one = request_with(&one, "sizes", &ragged.len().to_string());
    refused(request_with(&one, "loci", ragged), "ragged");
    // too many files, and too many sequences in one
    let many = "x\n".repeat(MAX_LOCI + 1);
    let many = request_with(&good, "names", &many);
    refused(
        request_with(&many, "sizes", &"1\n".repeat(MAX_LOCI + 1)),
        "at most",
    );
    let rows: String = (0..=MAX_SEQUENCES).map(|k| format!(">s{k}\nA\n")).collect();
    let big = request_with(&good, "names", "wide");
    let big = request_with(&big, "sizes", &rows.len().to_string());
    refused(request_with(&big, "loci", &rows), "at most");
    // a section that is not there
    let missing = Container::new().with("meta", "files=1\n");
    assert!(FlowJob::decode(&missing).unwrap_err().contains("pipeline"));
}

// --- running ------------------------------------------------------------------------------------------------------------------

#[test]
fn the_run_gives_what_the_engine_gives_and_a_run_sent_again_finds_its_answers() {
    let folder = Folder::new();
    let asked = job(5);
    let direct = engine::run(
        &asked.pipeline,
        asked.loci.clone(),
        None,
        &Default::default(),
    )
    .expect("a report");
    let first = run(asked.clone(), Some(&folder.cache())).expect("a result");
    assert_eq!(first.report.output, direct.output);
    // five gene trees and the species tree
    assert_eq!((first.report.ran(), first.report.cached()), (6, 0));
    let Some(Value::Tree(tree)) = &first.report.output else {
        panic!("a species tree");
    };
    assert!(tree.starts_with('(') && tree.contains("taxon1"));
    // sent again: nothing is run
    let again = run(asked.clone(), Some(&folder.cache())).expect("a result");
    assert_eq!((again.report.ran(), again.report.cached()), (0, 6));
    assert_eq!(again.report.output, first.report.output);
    // one locus changed: its tree and the gathering are run again, and the four others are found
    let mut changed = asked;
    changed.loci[2].1 = Value::Alignment(locus(999, 150));
    let after = run(changed, Some(&folder.cache())).expect("a result");
    assert_eq!((after.report.ran(), after.report.cached()), (2, 4));
}

#[test]
fn a_cache_that_cannot_be_used_is_said_and_the_run_goes_on() {
    let folder = Folder::new();
    // a regular file where the folder should be
    let blocked = folder.0.join("blocked");
    std::fs::write(&blocked, "not a folder").expect("a file");
    let done = run(job(3), Some(&blocked)).expect("a result");
    assert_eq!((done.report.ran(), done.report.cached()), (4, 0));
    assert!(done.report.output.is_some());
    assert!(
        done.report
            .warnings
            .iter()
            .any(|w| w.contains("cache could not be used")),
        "{:?}",
        done.report.warnings
    );
    // with no folder at all there is nothing to say
    let plain = run(job(3), None).expect("a result");
    assert!(plain.report.warnings.is_empty());
}

#[test]
fn a_run_with_nothing_to_run_on_is_refused_by_the_engine_not_the_wire() {
    let none = FlowJob {
        pipeline: Pipeline::loci_to_species(),
        loci: Vec::new(),
    };
    assert!(run(none, None).unwrap_err().contains("nothing to run on"));
}

// --- the result -----------------------------------------------------------------------------------------------------------------

/// A result that has a locus that failed, one that was left out for being small, and the steps of both kinds.
fn mixed() -> FlowResult {
    let mut asked = job(4);
    // one sequence has no tree
    asked.loci.push((
        "lonely".into(),
        Value::Alignment(vec![("a".into(), b"ACGTACGTAC".to_vec())]),
    ));
    run(asked, None).expect("a result")
}

#[test]
fn a_result_goes_through_bytes_and_comes_back_whole() {
    let done = mixed();
    assert_eq!(done.report.failed(), 1);
    assert_eq!(done.report.left_out, ["lonely"]);
    let back = FlowResult::decode(&through_bytes(&done.encode())).expect("a result");
    assert_eq!(back, done);
    assert_eq!(back.encode(), done.encode());
    // the pipelines of other shapes: the steps of each file are the first ones, the others are done once
    let unaligned = FlowJob {
        pipeline: Pipeline::unaligned_loci_to_species(),
        loci: (1..=3)
            .map(|k| (format!("locus{k}"), Value::Sequences(locus(k * 5, 90))))
            .collect(),
    };
    let done = run(unaligned, None).expect("a result");
    assert_eq!(done.report.items[0].steps.len(), 2);
    let back = FlowResult::decode(&through_bytes(&done.encode())).expect("a result");
    assert_eq!(back, done);
    // a pipeline that ends with a tree for each file gives no whole answer for many files, and the file's for one
    let trees =
        Pipeline::parse("helixloom-pipeline 1\ninput: alignment\nstep tree\n").expect("a pipeline");
    for n in [1, 3] {
        let done = run(
            FlowJob {
                pipeline: trees.clone(),
                loci: loci(n),
            },
            None,
        )
        .expect("a result");
        assert_eq!(done.report.output.is_some(), n == 1);
        assert!(done.report.items.iter().all(|i| i.output.is_some()));
        let back = FlowResult::decode(&through_bytes(&done.encode())).expect("a result");
        assert_eq!(back, done, "{n} files");
    }
}

#[test]
fn a_damaged_result_is_refused() {
    let good = mixed().encode();
    let report = good.text("report").unwrap().to_string();
    let refused = |c: Container, wanted: &str| {
        let why = FlowResult::decode(&c).expect_err(wanted);
        assert!(why.contains(wanted), "{wanted:?} is not in {why:?}");
    };
    assert!(FlowResult::decode(&good).is_ok());
    // the time
    for bad in ["NaN", "-1", "inf"] {
        refused(
            result_with(&good, "meta", &format!("seconds={bad}\n")),
            "time",
        );
    }
    refused(result_with(&good, "meta", "seconds=x\n"), "seconds");
    // the steps: one that is not of the pipeline, one that is not of the place it is in, a status that is none
    refused(
        result_with(
            &good,
            "report",
            &report.replacen("step\t0\t", "step\t7\t", 1),
        ),
        "not one of those done there",
    );
    refused(
        result_with(
            &good,
            "report",
            &report.replacen("step\t0\t", "step\t1\t", 1),
        ),
        "not one of those done there",
    );
    refused(
        result_with(&good, "report", &report.replacen("\tran\t", "\tmaybe\t", 1)),
        "neither run, found nor failed",
    );
    refused(
        result_with(
            &good,
            "report",
            &report.replacen("\tran\t", "\tran\tx\t", 1),
        ),
        "a row is not one",
    );
    // rows in the wrong places
    refused(
        result_with(&good, "report", "step\t0\tran\t0\t\t\n"),
        "belongs to no file",
    );
    refused(
        result_with(&good, "report", &format!("gathered\n{report}")),
        "not one that a report has there",
    );
    refused(result_with(&good, "report", "item\ta\t0\n"), "stops before");
    refused(
        result_with(&good, "report", &format!("{report}nonsense\n")),
        "not one that a report has there",
    );
    refused(
        result_with(&good, "report", &report.replacen("\t1\n", "\t2\n", 1)),
        "whether it has an answer",
    );
    // the answers: fewer or more than the files that have one, out of order, not what the pipeline gives
    let (names, sizes) = (
        good.text("out_names").unwrap().to_string(),
        good.text("out_sizes").unwrap().to_string(),
    );
    let mut dropped = names.lines().collect::<Vec<_>>();
    dropped.pop();
    let dropped_sizes = sizes
        .lines()
        .take(dropped.len())
        .collect::<Vec<_>>()
        .join("\n");
    let body = good.text("outputs").unwrap();
    let kept = &body[..dropped_sizes
        .lines()
        .map(|s| s.parse::<usize>().unwrap())
        .sum::<usize>()];
    let missing = result_with(&good, "out_names", &dropped.join("\n"));
    let missing = result_with(&missing, "out_sizes", &dropped_sizes);
    refused(
        result_with(&missing, "outputs", kept),
        "an answer is missing",
    );
    let swapped = names.lines().rev().collect::<Vec<_>>().join("\n");
    refused(
        result_with(&good, "out_names", &swapped),
        "not in the order",
    );
    refused(
        result_with(&good, "out_names", &format!("{names}\nextra")),
        "not as many sizes",
    );
    let extra_sizes = format!("{sizes}\n1");
    let extra = result_with(&good, "out_names", &format!("{names}\nextra"));
    let extra = result_with(&extra, "out_sizes", &extra_sizes);
    refused(
        result_with(&extra, "outputs", &format!("{body}x")),
        "no file has",
    );
    refused(
        result_with(&good, "output", "(A,B,C);\n(A,B,D);\n"),
        "one tree",
    );
    // a pipeline that is not one
    refused(
        result_with(&good, "pipeline", "helixloom-pipeline 1\n"),
        "input",
    );
    // sections that are not there
    assert!(FlowResult::decode(&Container::new().with("meta", "seconds=1\n")).is_err());
}
