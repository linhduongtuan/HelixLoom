use super::*;
use std::sync::atomic::{AtomicU32, Ordering};

static NEXT: AtomicU32 = AtomicU32::new(0);

/// A folder of its own for one test, removed when the test ends.
struct Folder(PathBuf);

impl Folder {
    fn new() -> Folder {
        let path = std::env::temp_dir().join(format!(
            "helixloom-flow-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }

    fn file(&self, name: &str, text: &str) -> PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, text).expect("a file");
        path
    }

    fn cache(&self) -> Cache {
        Cache::open(&self.0.join("cache")).expect("a cache")
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn next(seed: &mut u64) -> u64 {
    *seed = seed
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    *seed >> 33
}

/// A locus of six taxa in two groups of three: the groups differ in many places, and each taxon in a few of its own.
pub(crate) fn locus(seed: u64, len: usize) -> Vec<Named> {
    let mut s = seed;
    let base = |s: &mut u64| b"ACGT"[(next(s) % 4) as usize];
    let root: Vec<u8> = (0..len).map(|_| base(&mut s)).collect();
    let change = |seq: &mut Vec<u8>, times: usize, s: &mut u64| {
        for _ in 0..times {
            let at = next(s) as usize % len;
            seq[at] = base(s);
        }
    };
    let mut groups = [root.clone(), root];
    for g in &mut groups {
        change(g, len / 8, &mut s);
    }
    (1..=6)
        .map(|i| {
            let mut seq = groups[(i - 1) / 3].clone();
            change(&mut seq, len / 30, &mut s);
            (format!("taxon{i}"), seq)
        })
        .collect()
}

pub(crate) fn loci(n: u64) -> Vec<(String, Value)> {
    (1..=n)
        .map(|k| (format!("locus{k}"), Value::Alignment(locus(k * 17, 150))))
        .collect()
}

fn none() -> Available {
    Available::default()
}

// --- the file -------------------------------------------------------------------------------------------------------------------

#[test]
fn a_pipeline_is_written_and_read_back_the_same() {
    for p in [
        Pipeline::loci_to_species(),
        Pipeline::unaligned_loci_to_species(),
        Pipeline::inspect_graphs(),
    ] {
        assert_eq!(Pipeline::parse(&p.to_text()), Ok(p));
    }
}

#[test]
fn the_file_may_have_comments_blank_lines_and_leave_out_the_defaults() {
    let p = Pipeline::parse(
        "# my pipeline\n\nhelixloom-pipeline 1\nname: Loci\ninput: sequences\n\nstep align\nstep tree\nstep species\n",
    )
    .unwrap();
    assert_eq!(p.name, "Loci");
    assert_eq!(p.input, Kind::Sequences);
    assert_eq!(
        p.steps,
        [
            Step::Align { mafft: false },
            Step::Tree {
                iqtree: false,
                bootstrap: 0
            },
            Step::Species
        ]
    );
}

#[test]
fn a_file_that_cannot_be_read_says_which_line_and_why() {
    let head = "helixloom-pipeline 1\ninput: alignment\n";
    let err = |body: &str| Pipeline::parse(&format!("{head}{body}")).unwrap_err();
    assert!(err("step blast\n").contains("line 3") && err("step blast\n").contains("blast"));
    assert!(err("step tree colour=red\n").contains("no setting called colour"));
    assert!(err("step tree bootstrap=many\n").contains("bootstrap must be a whole number"));
    assert!(err("step tree engine=fasttree\n").contains("no engine called fasttree"));
    assert!(err("step tree bootstrap\n").contains("key=value"));
    assert!(err("hello\n").contains("expected name:, input: or step"));
    assert!(Pipeline::parse("input: tree\n")
        .unwrap_err()
        .contains("starts with"));
    assert!(Pipeline::parse("helixloom-pipeline 1\nstep species\n")
        .unwrap_err()
        .contains("no input"));
    assert!(Pipeline::parse("helixloom-pipeline 1\ninput: dna\n")
        .unwrap_err()
        .contains("sequences, alignment, tree, trees or gfa"));
    assert!(Pipeline::parse("").unwrap_err().contains("empty"));
}

#[test]
fn a_file_of_a_later_version_is_refused_with_what_to_do() {
    let why = Pipeline::parse("helixloom-pipeline 9\ninput: trees\nstep species\n").unwrap_err();
    assert!(why.contains("version 9") && why.contains("Update"), "{why}");
}

#[test]
fn each_step_must_take_what_the_one_before_gives() {
    let mut p = Pipeline::loci_to_species();
    assert_eq!(p.check(), Ok(Kind::Tree));
    p.steps.insert(0, Step::Align { mafft: false });
    let why = p.check().unwrap_err();
    assert!(
        why.contains("Step 1 (align) takes sequences") && why.contains("alignment"),
        "{why}"
    );
    p.input = Kind::Sequences;
    assert_eq!(p.check(), Ok(Kind::Tree));
    p.steps.clear();
    assert!(p.check().unwrap_err().contains("at least one step"));
}

#[test]
fn trees_are_gathered_for_a_step_that_takes_many_but_never_made_up() {
    // a tree for each file becomes the list; an alignment cannot become trees
    let p = Pipeline {
        name: String::new(),
        input: Kind::Alignment,
        steps: vec![Step::Species],
    };
    assert!(p.check().unwrap_err().contains("takes trees"));
}

#[test]
fn gfa_inspection_takes_only_graphs_and_its_result_is_a_summary() {
    let p = Pipeline::parse("helixloom-pipeline 1\ninput: gfa\nstep gfa-inspect").unwrap();
    assert_eq!(p.steps, [Step::GfaInspect]);
    assert_eq!(p.check(), Ok(Kind::GraphSummary));
    assert_eq!(p.per_file(), 1);
    assert_eq!(p.split_text(), None);
    assert_eq!(Kind::parse("graph-summary"), Some(Kind::GraphSummary));
    let rejected =
        |body: &str| Pipeline::parse(&format!("helixloom-pipeline 1\n{body}")).unwrap_err();
    assert!(rejected("input: alignment\nstep gfa-inspect").contains("takes gfa"));
    assert!(rejected("input: gfa\nstep gfa-inspect\nstep tree").contains("graph-summary"));
    assert!(rejected("input: gfa\nstep gfa-inspect option=1").contains("no setting called option"));
}

#[test]
fn bootstrap_with_the_built_in_search_is_refused_and_not_left_out_quietly() {
    let mut p = Pipeline::loci_to_species();
    p.steps[0] = Step::Tree {
        iqtree: false,
        bootstrap: 100,
    };
    assert!(p.check().unwrap_err().contains("engine=iqtree"));
    p.steps[0] = Step::Tree {
        iqtree: true,
        bootstrap: 1000,
    };
    assert!(p.check().is_ok());
}

#[test]
fn the_steps_before_the_first_that_takes_many_are_done_for_each_file() {
    assert_eq!(Pipeline::loci_to_species().per_file(), 1);
    assert_eq!(Pipeline::unaligned_loci_to_species().per_file(), 2);
    let on_trees = Pipeline {
        name: String::new(),
        input: Kind::Trees,
        steps: vec![Step::Species],
    };
    assert_eq!(on_trees.per_file(), 1);
    let no_gather = Pipeline {
        name: String::new(),
        input: Kind::Sequences,
        steps: vec![Step::Align { mafft: false }],
    };
    assert_eq!(no_gather.per_file(), 1);
}

#[test]
fn the_split_is_said_in_words_that_fit_one_step_or_many() {
    assert_eq!(
        Pipeline::loci_to_species().split_text().unwrap(),
        "Step 1 is done for each file; the trees are then gathered for the rest."
    );
    assert!(Pipeline::unaligned_loci_to_species()
        .split_text()
        .unwrap()
        .starts_with("Steps 1 to 2 are done"));
    let nothing = Pipeline {
        name: String::new(),
        input: Kind::Trees,
        steps: vec![Step::Species],
    };
    assert_eq!(nothing.split_text(), None);
}

// --- values and keys --------------------------------------------------------------------------------------------------------------

#[test]
fn values_go_to_text_and_back() {
    let rows = locus(1, 40);
    for v in [
        Value::Sequences(rows.clone()),
        Value::Alignment(rows),
        Value::Sequences(vec![
            ("empty".into(), Vec::new()),
            ("odd name | with, marks (1)".into(), b"AC-*".to_vec()),
        ]),
        Value::Tree("((a:1,b:1):1,c:2);".into()),
        Value::Trees(vec!["((a,b),c);".into(), "((a,c),b);".into()]),
    ] {
        assert_eq!(Value::from_text(v.kind(), &v.to_text()), Ok(v));
    }
}

#[test]
fn an_alignment_whose_rows_differ_in_length_is_not_one() {
    let why = Value::from_text(Kind::Alignment, ">a\nACGT\n>b\nACG\n").unwrap_err();
    assert!(
        why.contains("rows of one length") && why.starts_with('b'),
        "{why}"
    );
    assert!(Value::from_text(Kind::Tree, "(a,b);\n(a,c);\n").is_err());
}

#[test]
fn the_digest_is_the_same_every_time_and_tells_inputs_apart() {
    assert_eq!(digest(b"ACGT"), digest(b"ACGT"));
    assert_ne!(digest(b"ACGT"), digest(b"ACGA"));
    assert_ne!(digest(b""), digest(b"\0"));
    assert_eq!(digest(b"ACGT").len(), 32);
    // a fixed number, worked out apart from this code: the cache of one release is read by the next
    assert_eq!(digest(b"ACGT"), "9a90178ba8feda4e2fb0da422cf1fff0");
}

#[test]
fn what_an_answer_is_filed_under_depends_on_the_step_its_settings_its_program_and_its_input() {
    let input = Value::Alignment(locus(1, 40));
    let other = Value::Alignment(locus(2, 40));
    let tree = |bootstrap| Step::Tree {
        iqtree: true,
        bootstrap,
    };
    let base = key(&tree(100), "v3.1", &input);
    assert_eq!(base, key(&tree(100), "v3.1", &input));
    assert_ne!(base, key(&tree(200), "v3.1", &input), "a setting");
    assert_ne!(
        base,
        key(&tree(100), "v3.2", &input),
        "the program's version"
    );
    assert_ne!(base, key(&tree(100), "v3.1", &other), "the input");
    assert_ne!(
        key(&Step::Species, "", &Value::Trees(vec![])),
        key(
            &Step::Tree {
                iqtree: false,
                bootstrap: 0
            },
            "",
            &Value::Trees(vec![])
        ),
        "the step"
    );
}

// --- the cache --------------------------------------------------------------------------------------------------------------------

#[test]
fn an_answer_kept_is_found_again_with_its_note() {
    let f = Folder::new();
    let cache = f.cache();
    let v = Value::Tree("((a,b),c);".into());
    assert!(cache.get("k", Kind::Tree).is_none());
    cache.put("k", &v, "built in\nsecond line").unwrap();
    assert_eq!(
        cache.get("k", Kind::Tree),
        Some((v, "built in second line".to_string()))
    );
    assert_eq!(cache.len(), 1);
    assert!(
        cache.get("k", Kind::Alignment).is_none(),
        "another kind is not this answer"
    );
}

#[test]
fn a_damaged_or_foreign_file_is_not_an_answer_and_is_replaced() {
    let f = Folder::new();
    let cache = f.cache();
    let v = Value::Tree("((a,b),c);".into());
    cache.put("k", &v, "n").unwrap();
    let path = f.0.join("cache").join("k.value");
    for damage in [
        "",
        "helixloom-cache 1\nkind: tree\n",
        "something else\n",
        "helixloom-cache 1\nkind: tree\nnote: x\n(a,b;\n\n(c;\n",
    ] {
        std::fs::write(&path, damage).unwrap();
        assert!(cache.get("k", Kind::Tree).is_none(), "{damage:?}");
    }
    cache.put("k", &v, "n").unwrap();
    assert!(cache.get("k", Kind::Tree).is_some());
}

#[test]
fn keeping_an_answer_leaves_no_half_written_file_and_clearing_forgets_all() {
    let f = Folder::new();
    let cache = f.cache();
    for k in 0..5 {
        cache
            .put(&format!("k{k}"), &Value::Tree("(a,b);".into()), "")
            .unwrap();
    }
    let names: Vec<String> = std::fs::read_dir(f.0.join("cache"))
        .unwrap()
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    assert!(names.iter().all(|n| n.ends_with(".value")), "{names:?}");
    assert_eq!(cache.clear(), Ok(5));
    assert_eq!(cache.len(), 0);
}

// --- running ----------------------------------------------------------------------------------------------------------------------

#[test]
fn loci_become_gene_trees_and_a_species_tree() {
    let report = run(&Pipeline::loci_to_species(), loci(4), None, &none()).unwrap();
    assert_eq!(report.items.len(), 4);
    assert!(report
        .items
        .iter()
        .all(|i| matches!(i.output, Some(Value::Tree(_)))));
    assert_eq!((report.ran(), report.cached(), report.failed()), (5, 0, 0));
    let Some(Value::Tree(species)) = &report.output else {
        panic!("a species tree: {:?}", report.output);
    };
    for k in 1..=6 {
        assert!(species.contains(&format!("taxon{k}")), "{species}");
    }
    assert!(report.left_out.is_empty());
    assert_eq!(report.lines().len(), 5);
}

#[test]
fn a_second_run_finds_every_step_in_the_cache_and_gives_the_same_answer() {
    let f = Folder::new();
    let cache = f.cache();
    let p = Pipeline::unaligned_loci_to_species();
    let unaligned: Vec<(String, Value)> = loci(3)
        .into_iter()
        .map(|(n, v)| match v {
            Value::Alignment(rows) => (n, Value::Sequences(rows)),
            other => (n, other),
        })
        .collect();
    let first = run(&p, unaligned.clone(), Some(&cache), &none()).unwrap();
    assert_eq!((first.ran(), first.cached()), (7, 0));
    let second = run(&p, unaligned, Some(&cache), &none()).unwrap();
    assert_eq!((second.ran(), second.cached(), second.failed()), (0, 7, 0));
    assert_eq!(first.output, second.output);
    assert!(second
        .lines()
        .iter()
        .all(|l| l.contains("found in the cache")));
    assert!(second
        .summary()
        .starts_with("0 steps run, 7 found in the cache"));
}

#[test]
fn changing_one_locus_runs_again_only_what_depends_on_it() {
    let f = Folder::new();
    let cache = f.cache();
    let p = Pipeline::loci_to_species();
    run(&p, loci(4), Some(&cache), &none()).unwrap();
    let mut changed = loci(4);
    changed[2].1 = Value::Alignment(locus(999, 150));
    let report = run(&p, changed, Some(&cache), &none()).unwrap();
    // its gene tree is new, so the species tree, which is made of all of them, is new; the other three trees are found
    assert_eq!((report.ran(), report.cached()), (2, 3));
    let status = |name: &str| {
        report.items.iter().find(|i| i.name == name).unwrap().steps[0]
            .status
            .clone()
    };
    assert_eq!(status("locus3"), Status::Ran);
    assert_eq!(status("locus1"), Status::Cached);
}

#[test]
fn a_run_that_failed_in_one_locus_goes_on_with_the_rest_and_is_resumed_by_running_it_again() {
    let f = Folder::new();
    let cache = f.cache();
    let p = Pipeline::loci_to_species();
    let mut items = loci(3);
    items.push((
        "broken".into(),
        Value::Alignment(locus(5, 150).into_iter().take(2).collect()),
    ));
    let first = run(&p, items.clone(), Some(&cache), &none()).unwrap();
    assert_eq!((first.ran(), first.failed()), (4, 1));
    assert_eq!(first.left_out, ["broken"]);
    assert!(
        first.output.is_some(),
        "the species tree of the three that worked"
    );
    assert!(first.summary().contains("1 failed") && first.summary().contains("broken"));
    assert!(first
        .lines()
        .iter()
        .any(|l| l.starts_with("broken: tree: FAILED: A tree needs at least three")));

    // the broken locus is mended; running again does only what was not done
    items[3].1 = Value::Alignment(locus(6, 150));
    let again = run(&p, items, Some(&cache), &none()).unwrap();
    assert_eq!((again.ran(), again.cached(), again.failed()), (2, 3, 0));
    assert!(again.left_out.is_empty());
}

#[test]
fn a_damaged_cache_file_costs_one_step_not_the_run() {
    let f = Folder::new();
    let cache = f.cache();
    let p = Pipeline::loci_to_species();
    run(&p, loci(3), Some(&cache), &none()).unwrap();
    let one = std::fs::read_dir(f.0.join("cache"))
        .unwrap()
        .flatten()
        .map(|e| e.path())
        .find(|p| {
            p.file_name()
                .unwrap()
                .to_string_lossy()
                .starts_with("tree-")
        })
        .unwrap();
    std::fs::write(&one, "garbage").unwrap();
    let report = run(&p, loci(3), Some(&cache), &none()).unwrap();
    assert_eq!((report.ran(), report.cached(), report.failed()), (1, 3, 0));
}

#[test]
fn without_a_cache_everything_runs_every_time() {
    let p = Pipeline::loci_to_species();
    let a = run(&p, loci(3), None, &none()).unwrap();
    let b = run(&p, loci(3), None, &none()).unwrap();
    assert_eq!((a.ran(), b.ran(), b.cached()), (4, 4, 0));
    assert_eq!(
        a.output, b.output,
        "the same answer, whichever core did which locus"
    );
}

#[test]
fn a_step_that_asks_for_a_program_that_is_not_there_fails_and_does_not_use_another() {
    let mut p = Pipeline::unaligned_loci_to_species();
    p.steps[0] = Step::Align { mafft: true };
    let items = vec![("l".to_string(), Value::Sequences(locus(1, 60)))];
    let report = run(&p, items, None, &none()).unwrap();
    let Status::Failed(why) = &report.items[0].steps[0].status else {
        panic!("{:?}", report.items[0].steps);
    };
    assert!(
        why.contains("MAFFT is not installed") && why.contains("not used instead"),
        "{why}"
    );
}

#[test]
fn the_input_is_checked_before_anything_runs() {
    let p = Pipeline::loci_to_species();
    assert!(run(&p, Vec::new(), None, &none())
        .unwrap_err()
        .contains("no input"));
    let wrong = vec![("a".to_string(), Value::Sequences(locus(1, 30)))];
    let why = run(&p, wrong, None, &none()).unwrap_err();
    assert!(
        why.contains("a is sequences") && why.contains("takes alignment"),
        "{why}"
    );
}

#[test]
fn a_single_file_gives_its_own_answer_when_nothing_is_gathered() {
    let p = Pipeline {
        name: String::new(),
        input: Kind::Alignment,
        steps: vec![Step::Tree {
            iqtree: false,
            bootstrap: 0,
        }],
    };
    let one = run(&p, loci(1), None, &none()).unwrap();
    assert!(matches!(one.output, Some(Value::Tree(_))));
    let two = run(&p, loci(2), None, &none()).unwrap();
    assert!(
        two.output.is_none(),
        "two files, two answers: they are in the items"
    );
    assert!(two.items.iter().all(|i| i.output.is_some()));
}

#[test]
fn when_no_gene_tree_could_be_made_the_species_step_says_so() {
    let p = Pipeline::loci_to_species();
    let items = vec![(
        "only".to_string(),
        Value::Alignment(locus(1, 60).into_iter().take(2).collect()),
    )];
    let report = run(&p, items, None, &none()).unwrap();
    assert_eq!(report.failed(), 2);
    assert!(report.output.is_none());
    assert!(report
        .lines()
        .iter()
        .any(|l| l.contains("nothing to gather")));
}

// --- a folder ---------------------------------------------------------------------------------------------------------------------

#[test]
fn a_folder_is_read_by_name_and_only_the_files_the_pipeline_takes() {
    let f = Folder::new();
    f.file("b.fasta", ">x\nACGT\n>y\nACGA\n");
    f.file("a.FA", ">x\nACGT\n>y\nACGA\n");
    f.file("notes.txt", "not a sequence file");
    f.file("c.tre", "(a,b);");
    let got = read_folder(&f.0, Kind::Alignment).unwrap();
    let names: Vec<&str> = got.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(names, ["a", "b"]);
    assert!(matches!(got[0].1, Value::Alignment(_)));
    let trees = read_folder(&f.0, Kind::Trees).unwrap();
    assert_eq!(trees.len(), 1);
    assert_eq!(trees[0].1, Value::Trees(vec!["(a,b);".into()]));
}

#[test]
fn a_folder_without_the_files_or_with_a_file_that_is_not_an_alignment_says_so() {
    let f = Folder::new();
    f.file("notes.txt", "x");
    let why = read_folder(&f.0, Kind::Alignment).unwrap_err();
    assert!(
        why.contains("holds no file of alignment") && why.contains("fasta"),
        "{why}"
    );
    f.file("ragged.fasta", ">x\nACGT\n>y\nAC\n");
    let why = read_folder(&f.0, Kind::Alignment).unwrap_err();
    assert!(
        why.contains("ragged.fasta") && why.contains("rows of one length"),
        "{why}"
    );
    assert!(
        read_folder(&f.0, Kind::Sequences).is_ok(),
        "unaligned sequences may differ in length"
    );
    assert!(read_folder(&f.0.join("nowhere"), Kind::Sequences).is_err());
}

#[test]
fn a_folder_of_loci_runs_end_to_end_from_files() {
    let f = Folder::new();
    for (k, (_, v)) in loci(3).into_iter().enumerate() {
        f.file(&format!("gene{}.fasta", k + 1), &v.to_text());
    }
    let items = read_folder(&f.0, Kind::Alignment).unwrap();
    let cache = f.cache();
    let report = run(&Pipeline::loci_to_species(), items, Some(&cache), &none()).unwrap();
    assert_eq!((report.ran(), report.failed()), (4, 0));
    assert_eq!(report.items[0].name, "gene1");
}

#[test]
fn a_folder_of_graphs_preserves_input_and_resumes_one_failed_inspection() {
    let f = Folder::new();
    let cache = f.cache();
    let p = Pipeline::inspect_graphs();
    let good = "L\ta\t+\tb\t-\t1M\r\nS\ta\tAC\r\nS\tb\t*\r\nS\tc\t*\tLN:i:4";
    let invalid = "H\tVN:Z:1.0\nL\tx\t+\tx\t+\t0M";
    f.file("good.GFA", good);
    f.file("bad.gfa", invalid);
    f.file("unrelated.tsv", "old report");
    let inputs = read_folder(&f.0, Kind::Gfa).unwrap();
    assert_eq!(inputs[1], ("good".into(), Value::Gfa(good.into())));
    assert_eq!(inputs.len(), 2);
    let run_folder = || {
        run(
            &p,
            read_folder(&f.0, p.input).unwrap(),
            Some(&cache),
            &none(),
        )
        .unwrap()
    };
    let first = run_folder();
    assert_eq!((first.ran(), first.failed(), first.cached()), (1, 1, 0));
    assert_eq!(first.left_out, ["bad"]);
    assert!(first.items[0].output.is_none());
    assert!(first
        .lines()
        .iter()
        .any(|l| l.contains("bad: gfa-inspect: FAILED: line 2")));
    let expected = "segments\t3\nlinks\t1\ncomponents\t2\nisolated_segments\t1\nknown_segment_bases\t6\nsegments_without_length\t1\n";
    let checked = first.items[1].output.as_ref().unwrap();
    assert_eq!(checked.to_text(), expected);
    assert_eq!(checked.kind(), Kind::GraphSummary);
    assert_eq!(first.items[1].steps[0].program(), None);
    assert_eq!(cache.len(), 1);

    f.file("bad.gfa", "H\tVN:Z:1.0\n");
    let second = run_folder();
    assert_eq!((second.ran(), second.failed(), second.cached()), (1, 0, 1));
    assert_eq!(second.items[1].output.as_ref(), Some(checked));
    assert_eq!(
        second.items[0].output,
        Some(Value::GraphSummary(gfa::Summary::default()))
    );
    f.file("good.GFA", &format!("# changed\n{good}"));
    let third = run_folder();
    assert_eq!((third.ran(), third.failed(), third.cached()), (1, 0, 1));
    assert_eq!(third.items[1].output.as_ref(), Some(checked));
}

#[test]
fn gfa_and_graph_summary_keep_their_formats_in_the_cache_and_wire_values() {
    let marked = "\u{feff}H\tVN:Z:1.0\n";
    let graph = Value::from_text(Kind::Gfa, marked).unwrap();
    assert_eq!(graph.as_text(), marked);
    assert!(work(&Step::GfaInspect, &graph, None)
        .unwrap_err()
        .contains("line 1"));
    assert!(work(&Step::GfaInspect, &Value::Gfa(String::new()), None).is_err());
    let f = Folder::new();
    let cache = f.cache();
    let input = Value::Gfa("S\ts\tAC".into());
    let cache_key = key(&Step::GfaInspect, "", &input);
    let saved = run(
        &Pipeline::inspect_graphs(),
        vec![("graph".into(), input.clone())],
        Some(&cache),
        &none(),
    )
    .unwrap();
    let item = saved.output.unwrap();
    assert_eq!(
        Value::from_text(Kind::GraphSummary, &item.to_text()),
        Ok(item.clone())
    );
    assert_eq!(cache.get(&cache_key, Kind::GraphSummary).unwrap().0, item);
    let malformed = "segments\t0\n";
    std::fs::write(
        cache.file(&cache_key),
        format!(
            "{CACHE_HEAD}\nkind: graph-summary\nnote: old\ndigest: {}\n{malformed}",
            digest(malformed.as_bytes())
        ),
    )
    .unwrap();
    assert!(cache.get(&cache_key, Kind::GraphSummary).is_none());
    let again = run(
        &Pipeline::inspect_graphs(),
        vec![("graph".into(), input)],
        Some(&cache),
        &none(),
    )
    .unwrap();
    assert_eq!((again.ran(), again.cached()), (1, 0));
    assert_eq!(again.output, Some(item));
}

// --- what the review found ------------------------------------------------------------------------------------------------------

#[test]
fn a_cache_file_that_was_cut_short_or_changed_is_not_an_answer() {
    let f = Folder::new();
    let cache = f.cache();
    let v = Value::Tree("((a:1,b:1):1,(c:1,d:1):1);".into());
    cache.put("k", &v, "n").unwrap();
    let path = f.0.join("cache").join("k.value");
    let whole = std::fs::read_to_string(&path).unwrap();
    // the file as a crash can leave it: the head, then part of the tree
    let cut = &whole[..whole.len() - 8];
    std::fs::write(&path, cut).unwrap();
    assert!(
        cache.get("k", Kind::Tree).is_none(),
        "a cut tree is one line and would read as a tree"
    );
    // the tree changed, the digest not
    std::fs::write(&path, whole.replace("(a:1", "(z:1")).unwrap();
    assert!(cache.get("k", Kind::Tree).is_none());
    std::fs::write(&path, &whole).unwrap();
    assert_eq!(cache.get("k", Kind::Tree).map(|(v, _)| v), Some(v));
}

#[test]
fn many_threads_writing_and_reading_one_answer_never_see_half_of_it() {
    let f = Folder::new();
    let cache = f.cache();
    let long = |c: char| {
        Value::Tree(format!(
            "({},{});",
            c.to_string().repeat(4000),
            c.to_string().repeat(4000)
        ))
    };
    let (a, b) = (long('a'), long('b'));
    cache.put("k", &a, "").unwrap();
    std::thread::scope(|s| {
        for k in 0..4 {
            let (cache, value) = (&cache, if k % 2 == 0 { &a } else { &b });
            s.spawn(move || {
                for _ in 0..60 {
                    cache.put("k", value, "").unwrap();
                }
            });
        }
        s.spawn(|| {
            for _ in 0..400 {
                let (got, _) = cache
                    .get("k", Kind::Tree)
                    .expect("a whole answer, whichever was last");
                assert!(got == a || got == b);
            }
        });
    });
    let names: Vec<String> = std::fs::read_dir(f.0.join("cache"))
        .unwrap()
        .flatten()
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect();
    assert_eq!(names, ["k.value"], "no temporary file is left");
}

#[test]
fn forgetting_the_answers_leaves_every_file_that_is_not_the_programs() {
    let f = Folder::new();
    let cache = f.cache();
    cache.put("a", &Value::Tree("(a,b);".into()), "").unwrap();
    cache.put("b", &Value::Tree("(a,b);".into()), "").unwrap();
    let foreign = f.0.join("cache").join("notes.value");
    std::fs::write(&foreign, "somebody else's file").unwrap();
    assert_eq!(cache.len(), 2);
    assert_eq!(cache.clear(), Ok(2));
    assert!(foreign.exists());
}

#[test]
fn counting_a_folder_that_is_not_there_does_not_make_it() {
    let f = Folder::new();
    let missing = f.0.join("not-yet");
    assert_eq!(Cache::kept_in(&missing), 0);
    assert!(!missing.exists());
}

#[test]
fn a_gene_tree_too_small_to_use_is_left_out_and_said_to_be() {
    let mut items = loci(3);
    // four taxa would do; three leaves cannot say how species are related
    let small: Vec<Named> = locus(77, 150).into_iter().take(3).collect();
    items.push(("small".into(), Value::Alignment(small)));
    let report = run(&Pipeline::loci_to_species(), items, None, &none()).unwrap();
    assert_eq!(report.failed(), 0, "its tree was made");
    assert_eq!(report.left_out, ["small"]);
    assert_eq!(report.warnings.len(), 1);
    assert!(
        report.warnings[0].starts_with("small: its gene tree has 3 leaves"),
        "{}",
        report.warnings[0]
    );
    assert!(report.summary().contains("left out: small"));
    assert!(matches!(report.output, Some(Value::Tree(_))));
    let species = report.gathered[0].note.clone();
    assert!(species.starts_with("3 gene trees"), "{species}");
}

#[test]
fn bootstrap_with_iqtree_needs_the_thousand_that_iqtree_makes_at_least() {
    let mut p = Pipeline::loci_to_species();
    for (n, ok) in [
        (0, true),
        (100, false),
        (999, false),
        (1000, true),
        (5000, true),
    ] {
        p.steps[0] = Step::Tree {
            iqtree: true,
            bootstrap: n,
        };
        assert_eq!(p.check().is_ok(), ok, "bootstrap {n}");
    }
    p.steps[0] = Step::Tree {
        iqtree: true,
        bootstrap: 100,
    };
    assert!(p.check().unwrap_err().contains("at least 1000"));
}

#[test]
fn a_byte_order_mark_is_not_part_of_the_first_name() {
    let text = "\u{feff}>first\nACGT\n>second\nACGA\n";
    let Ok(Value::Sequences(rows)) = Value::from_text(Kind::Sequences, text) else {
        panic!("sequences");
    };
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].0, "first");
}

#[test]
fn a_fault_in_one_file_is_a_failed_step_and_not_the_end_of_the_run() {
    let step = Step::Tree {
        iqtree: false,
        bootstrap: 0,
    };
    let (reports, output) = guarded(std::slice::from_ref(&step), None, || {
        panic!("a made-up fault")
    });
    assert!(output.is_none());
    assert_eq!(reports.len(), 1);
    assert_eq!(reports[0].step, step);
    let Status::Failed(why) = &reports[0].status else {
        panic!("{:?}", reports[0].status);
    };
    assert!(
        why.contains("the program faulted: a made-up fault"),
        "{why}"
    );
    let (reports, output) = guarded(&[], Some(&Step::Species), || {
        (Vec::new(), Some(Value::Tree("(a,b);".into())))
    });
    assert!(
        reports.is_empty() && output.is_some(),
        "no fault, no change"
    );
}

#[test]
fn two_files_that_would_be_called_the_same_are_refused_even_when_other_files_sort_between() {
    let f = Folder::new();
    f.file("a.fa", ">x\nAC\n");
    f.file("a.fasta", ">x\nAC\n");
    let why = read_folder(&f.0, Kind::Sequences).unwrap_err();
    assert!(
        why.contains("a.fa and a.fasta would both be called a"),
        "{why}"
    );
    std::fs::remove_file(f.0.join("a.fa")).unwrap();
    // "a.fas.fa" sorts between "a.fa" and "a.fasta", and its name is a.fas
    f.file("a.fa", ">x\nAC\n");
    f.file("a.fas.fa", ">x\nAC\n");
    assert!(read_folder(&f.0, Kind::Sequences).is_err());
    std::fs::remove_file(f.0.join("a.fa")).unwrap();
    assert!(
        read_folder(&f.0, Kind::Sequences).is_ok(),
        "a.fasta and a.fas.fa differ"
    );
}

#[test]
fn gene_trees_can_be_the_input_and_each_file_is_one_gene_tree() {
    let p = Pipeline {
        name: String::new(),
        input: Kind::Tree,
        steps: vec![Step::Species],
    };
    assert_eq!(p.check(), Ok(Kind::Tree));
    assert_eq!(p.per_file(), 0);
    assert_eq!(
        p.split_text().unwrap(),
        "Each file is a gene tree; they are gathered for the steps."
    );
    let first = run(&Pipeline::loci_to_species(), loci(4), None, &none()).unwrap();
    let trees: Vec<(String, Value)> = first
        .items
        .iter()
        .map(|i| (i.name.clone(), i.output.clone().unwrap()))
        .collect();
    let again = run(&p, trees, None, &none()).unwrap();
    assert!(again.items.iter().all(|i| i.steps.is_empty()));
    assert_eq!(
        again.output, first.output,
        "the same gene trees give the same species tree"
    );
}

/// What the built-in steps give for a fixed input, as digests. If this fails, a built-in step has begun to answer differently
/// for the same input: raise `ENGINE`, so that answers kept by earlier runs are not served as if they were new, and put the new
/// digests here. (They were taken from the program as it was when the cache was built; they say that an answer has changed, not
/// that it is right.)
#[test]
fn the_built_in_steps_answer_as_they_did_when_the_engine_number_was_set() {
    let rows: Vec<Named> = locus(5, 90);
    let sequences = Value::Sequences(rows.clone());
    let (aligned, _) = work(&Step::Align { mafft: false }, &sequences, None).unwrap();
    let (tree, _) = work(
        &Step::Tree {
            iqtree: false,
            bootstrap: 0,
        },
        &aligned,
        None,
    )
    .unwrap();
    let Value::Tree(newick) = &tree else {
        panic!("a tree")
    };
    let trees = Value::Trees(
        (1..=4)
            .map(|k| {
                let (t, _) = work(
                    &Step::Tree {
                        iqtree: false,
                        bootstrap: 0,
                    },
                    &Value::Alignment(locus(k, 90)),
                    None,
                )
                .unwrap();
                t.to_text().trim().to_string()
            })
            .collect(),
    );
    let (species, _) = work(&Step::Species, &trees, None).unwrap();
    let got = [
        digest(aligned.to_text().as_bytes()),
        digest(newick.as_bytes()),
        digest(species.to_text().as_bytes()),
    ];
    assert_eq!(
        got, GOLDEN,
        "ENGINE is {ENGINE}: a built-in step answers differently now"
    );
}

const GOLDEN: [&str; 3] = [
    "051644ca0eb4f6aeeb2dbc31fa2d61e6",
    "0e9111598a620c18e817255c4e7649c6",
    "9cab2bcd39c5246534216dfd21069706",
];
