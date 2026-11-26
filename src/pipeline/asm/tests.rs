//! The steps of the assembly routes as a pipeline: their settings, the files they read, the routes as pipeline files, and a route run whole.

use super::*;
use crate::pipeline::flow::{key, read_folder, run, Cache, Pipeline, Status, Step};
use crate::pipeline::Available;
use crate::seq::debruijn::sim::{random, reads as short_reads};
use crate::seq::dna::revcom;
use crate::seq::fastq::FastqRead;
use crate::seq::hybrid::sim::{long_reads, ring};

pub(crate) fn folder(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("helixloom-asm-{name}-{}", std::process::id()));
    std::fs::remove_dir_all(&d).ok();
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Reads of quality 40 from bases.
pub(crate) fn fastq_reads(seqs: &[Vec<u8>]) -> Vec<FastqRead> {
    seqs.iter()
        .enumerate()
        .map(|(i, s)| FastqRead {
            name: format!("r{i}"),
            seq: s.clone(),
            qual: vec![40; s.len()],
        })
        .collect()
}

fn step(line: &str) -> Result<Step, String> {
    let text = format!("helixloom-pipeline 1\ninput: reads\n{line}\n");
    let text = if line.contains("hybrid") {
        text.replace("input: reads", "input: gfa")
    } else if line.starts_with("step concatemer")
        || line.starts_with("step circularize")
        || line.starts_with("step rotate")
        || line.starts_with("step medaka")
    {
        text.replace("input: reads", "input: sequences")
    } else {
        text
    };
    Pipeline::parse(&text).map(|p| p.steps[0].clone())
}

#[test]
fn every_step_has_the_defaults_that_are_written_in_its_file() {
    let p = |line: &str| step(line).unwrap().params();
    assert_eq!(
        p("step trim"),
        [
            ("adapter", "AGATCGGAAGAGC".to_string()),
            ("quality", "20".into()),
            ("min-length", "20".into()),
            ("cut", "tail".into())
        ]
    );
    assert_eq!(
        p("step filter"),
        [
            ("min-length", "2000".to_string()),
            ("target-bases", "500000000".into()),
            ("window", "100".into())
        ]
    );
    assert_eq!(
        p("step assemble"),
        [
            ("ks", "auto".to_string()),
            ("cov-cutoff", "auto".into()),
            ("min-qual", "10".into())
        ]
    );
    assert_eq!(
        p("step hybrid long=x"),
        [
            ("long", "x".to_string()),
            ("mode", "normal".into()),
            ("min-quality", "auto".into()),
            ("single-depth", "auto".into()),
            ("overlap", "0".into())
        ]
    );
    assert_eq!(
        p("step flye"),
        [
            ("mode", "nano-raw".to_string()),
            ("genome-size", "auto".into())
        ]
    );
    assert_eq!(
        p("step concatemer"),
        [
            ("min-identity", "95".to_string()),
            ("min-unit", "200".into())
        ]
    );
    assert_eq!(
        p("step circularize"),
        [
            ("min-overlap", "30".to_string()),
            ("min-identity", "95".into())
        ]
    );
    assert_eq!(p("step rotate genes=g.fa"), [("genes", "g.fa".to_string())]);
    assert_eq!(p("step medaka reads=r"), [("reads", "r".to_string())]);
}

#[test]
fn what_a_step_takes_and_gives_is_what_the_routes_need() {
    let kinds = |line: &str| {
        let s = step(line).unwrap();
        (s.takes().word(), s.gives().word())
    };
    assert_eq!(kinds("step trim"), ("reads", "reads"));
    assert_eq!(kinds("step filter"), ("reads", "reads"));
    assert_eq!(kinds("step assemble"), ("reads", "gfa"));
    assert_eq!(kinds("step hybrid long=x"), ("gfa", "sequences"));
    assert_eq!(kinds("step flye"), ("reads", "sequences"));
    for s in [
        "step concatemer",
        "step circularize",
        "step rotate genes=g",
        "step medaka reads=r",
    ] {
        assert_eq!(kinds(s), ("sequences", "sequences"), "{s}");
    }
}

#[test]
fn every_setting_that_is_not_a_good_one_is_refused_with_the_step_and_the_setting() {
    let why = |line: &str| step(line).err().unwrap();
    assert!(why("step trim colour=red").contains("trim has no setting called colour"));
    assert!(why("step trim quality=high").contains("quality must be a number"));
    assert!(why("step trim adapter=AC").contains("adapter is none or at least four bases"));
    assert!(why("step trim adapter=ACG").contains("adapter is none or at least four bases"));
    assert!(why("step trim adapter=ACGX").contains("adapter is none"));
    assert!(why("step trim cut=middle").contains("cut is tail, right, off"));
    assert!(why("step filter window=wide").contains("window must be a number"));
    assert!(why("step assemble ks=22").contains("odd"));
    assert!(why("step assemble ks=33,21").contains("must rise"));
    assert!(why("step assemble ks=x").contains("list of numbers"));
    assert!(why("step assemble cov-cutoff=-1")
        .contains("cov-cutoff is auto or off or a number of 0 or more"));
    // a number that is not one, and a number that is not finite
    assert!(why("step assemble cov-cutoff=inf")
        .contains("cov-cutoff is auto or off or a number of 0 or more"));
    assert!(why("step assemble cov-cutoff=nan")
        .contains("cov-cutoff is auto or off or a number of 0 or more"));
    assert!(why("step hybrid long=x mode=fast").contains("mode is conservative, normal, bold"));
    assert!(why("step hybrid long=x min-quality=101")
        .contains("min-quality is auto or a number from 0 to 100"));
    assert!(why("step hybrid long=x single-depth=0")
        .contains("single-depth is auto or a number above 0"));
    assert!(why("step hybrid long=x overlap=big").contains("overlap must be a number"));
    assert!(why("step flye mode=fast").contains("mode is nano-raw, nano-corr, nano-hq"));
    assert!(why("step flye genome-size=big").contains("genome-size is auto or a size"));
    assert!(why("step flye genome-size=5mm").contains("genome-size is auto or a size"));
    assert!(why("step concatemer min-identity=101").contains("a share from 0 to 100"));
    assert!(why("step circularize min-overlap=x").contains("min-overlap must be a number"));
}

#[test]
fn a_step_that_needs_a_file_says_which_when_it_is_not_given() {
    assert!(step("step hybrid")
        .err()
        .unwrap()
        .contains("hybrid needs long="));
    assert!(step("step rotate")
        .err()
        .unwrap()
        .contains("rotate needs genes="));
    assert!(step("step medaka")
        .err()
        .unwrap()
        .contains("medaka needs reads="));
}

#[test]
fn the_settings_that_are_good_are_taken() {
    assert!(step("step assemble ks=21,33,55 cov-cutoff=off min-qual=0").is_ok());
    assert!(step("step assemble cov-cutoff=2.5").is_ok());
    assert!(step("step flye genome-size=60k mode=nano-hq").is_ok());
    assert!(step("step flye genome-size=2.6m").is_ok());
    assert!(step("step hybrid long=a mode=bold min-quality=0 single-depth=12.5").is_ok());
    assert!(step("step trim adapter=none cut=right").is_ok());
    assert!(step("step concatemer min-identity=100 min-unit=1").is_ok());
}

#[test]
fn a_step_of_another_name_is_not_one_of_these() {
    assert_eq!(AsmStep::parse("align", &[]).unwrap(), None);
    let e = Pipeline::parse("helixloom-pipeline 1\ninput: reads\nstep sort\n")
        .err()
        .unwrap();
    assert!(
        e.contains("there is no step called sort") && e.contains("medaka"),
        "{e}"
    );
}

#[test]
fn a_word_is_added_to_a_name_once_and_looked_for_whole() {
    assert_eq!(tag("contig_1", "circular=true"), "contig_1 circular=true");
    assert_eq!(
        tag("contig_1 circular=true", "circular=true"),
        "contig_1 circular=true"
    );
    assert!(has_tag("1 length=5 circular=true", "circular=true"));
    assert!(!has_tag("1 length=5 notcircular=true", "circular=true"));
    assert!(!has_tag("circular=truely", "circular=true"));
}

#[test]
fn the_file_that_a_step_reads_is_the_one_named_or_the_one_in_the_folder_with_the_input_s_name() {
    let d = folder("side");
    std::fs::write(d.join("a.fastq"), "").unwrap();
    std::fs::write(d.join("b.fa"), "").unwrap();
    let hybrid = |long: &str| step(&format!("step hybrid long={long}")).unwrap();
    let Step::Asm(by_file) = hybrid(d.join("a.fastq").to_str().unwrap()) else {
        panic!()
    };
    assert_eq!(
        by_file.side_file("whatever").unwrap(),
        Some(d.join("a.fastq"))
    );
    let Step::Asm(by_folder) = hybrid(d.to_str().unwrap()) else {
        panic!()
    };
    assert_eq!(by_folder.side_file("a").unwrap(), Some(d.join("a.fastq")));
    assert_eq!(by_folder.side_file("b").unwrap(), Some(d.join("b.fa")));
    let e = by_folder.side_file("c").unwrap_err();
    assert!(
        e.contains("holds no file called c") && e.contains("fastq, fq, fasta, fa"),
        "{e}"
    );
    let Step::Asm(nowhere) = hybrid("/no/such/place") else {
        panic!()
    };
    assert!(nowhere
        .side_file("a")
        .unwrap_err()
        .contains("is not a file or a folder"));
    // a step with no side file has none
    let Step::Asm(trim) = step("step trim").unwrap() else {
        panic!()
    };
    assert_eq!(trim.side_file("a").unwrap(), None);
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn what_the_key_holds_changes_when_the_side_file_does() {
    let d = folder("context");
    let long = d.join("long.fasta");
    std::fs::write(&long, ">r\nACGT\n").unwrap();
    let Step::Asm(h) = step(&format!("step hybrid long={}", long.display())).unwrap() else {
        panic!()
    };
    let before = h.context("x");
    assert_eq!(before, h.context("x"));
    std::fs::write(&long, ">r\nACGA\n").unwrap();
    assert_ne!(before, h.context("x"));
    // a file that is not there is not an error here: the step says so when it runs
    std::fs::remove_file(&long).unwrap();
    assert!(h.context("x").contains("missing"));
    let Step::Asm(t) = step("step trim").unwrap() else {
        panic!()
    };
    assert_eq!(t.context("x"), "|");
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn each_route_is_a_pipeline_file_that_reads_checks_and_writes_back_the_same() {
    for (name, _) in routes::NAMES {
        let p = routes::named(name).unwrap_or_else(|| panic!("{name}"));
        let text = p.to_text();
        let back = Pipeline::parse(&text).unwrap_or_else(|e| panic!("{name}: {e}\n{text}"));
        assert_eq!(back, p, "{name}");
        back.check().unwrap_or_else(|e| panic!("{name}: {e}"));
    }
    assert!(routes::named("nothing").is_none());
}

#[test]
fn the_routes_have_the_steps_of_their_papers_in_order() {
    let words = |name: &str| -> Vec<&'static str> {
        routes::named(name)
            .unwrap()
            .steps
            .iter()
            .map(Step::name)
            .collect()
    };
    assert_eq!(words("plasmid"), ["trim", "assemble", "hybrid"]);
    assert_eq!(words("plasmid-short"), ["trim", "assemble"]);
    assert_eq!(words("plasmid-long"), ["hybrid"]);
    assert_eq!(
        words("nanopore-plasmid"),
        ["flye", "concatemer", "circularize", "rotate", "medaka"]
    );
    assert_eq!(words("minion-reads"), ["filter"]);
    assert_eq!(routes::named("plasmid-long").unwrap().input.word(), "gfa");
    // the numbers of the paper
    let p = routes::named("minion-reads").unwrap().to_text();
    assert!(
        p.contains("min-length=2000 target-bases=500000000 window=100"),
        "{p}"
    );
    let n = routes::named("nanopore-plasmid").unwrap().to_text();
    assert!(n.contains("min-overlap=30 min-identity=95"), "{n}");
}

#[test]
fn a_step_that_does_not_fit_what_comes_to_it_is_refused_when_the_file_is_read() {
    let e = Pipeline::parse("helixloom-pipeline 1\ninput: reads\nstep hybrid long=x\n")
        .err()
        .unwrap();
    assert!(
        e.contains("Step 1 (hybrid) takes gfa, but reads is what comes to it"),
        "{e}"
    );
    let e = Pipeline::parse("helixloom-pipeline 1\ninput: sequences\nstep trim\n")
        .err()
        .unwrap();
    assert!(e.contains("trim") && e.contains("takes reads"), "{e}");
}

/// A folder of the short reads of a plasmid with a repeat that they cannot span, and a folder of long reads that can.
fn plasmid_folders(name: &str) -> (PathBuf, PathBuf, Vec<u8>) {
    let d = folder(name);
    let (short, long) = (d.join("short"), d.join("long"));
    std::fs::create_dir_all(&short).unwrap();
    std::fs::create_dir_all(&long).unwrap();
    let r = ring(3000, 1200, 1);
    let pieces = short_reads(&r.genome, true, 100, 40.0, 0.0, 2);
    std::fs::write(
        short.join("sample1.fastq"),
        crate::seq::fastq::to_fastq(&fastq_reads(&pieces)),
    )
    .unwrap();
    let reads = long_reads(&r.genome, true, 5000, 40.0, 0.02, 3);
    let fasta: String = reads
        .iter()
        .map(|(n, s)| format!(">{n}\n{}\n", String::from_utf8_lossy(s)))
        .collect();
    std::fs::write(long.join("sample1.fasta"), fasta).unwrap();
    (short, long, r.genome)
}

#[test]
fn the_galaxy_route_runs_whole_and_gives_the_plasmid_as_one_circle_and_the_cache_keeps_it() {
    let (short, long, genome) = plasmid_folders("route1");
    let text = routes::plasmid()
        .to_text()
        .replace("long=long-reads", &format!("long={}", long.display()));
    let p = Pipeline::parse(&text).unwrap();
    let cache = Cache::open(&short.parent().unwrap().join("cache")).unwrap();
    let items = read_folder(&short, Kind::Reads).unwrap();
    let first = run(&p, items.clone(), Some(&cache), &Available::default()).unwrap();
    let item = &first.items[0];
    assert!(
        item.steps.iter().all(|s| matches!(s.status, Status::Ran)),
        "{:?}",
        item.steps
    );
    let Some(Value::Sequences(rows)) = &item.output else {
        panic!("{:?}", item.output)
    };
    assert_eq!(rows.len(), 1);
    assert!(has_tag(&rows[0].0, "circular=true"), "{}", rows[0].0);
    assert!(rows[0].0.contains("length=11400"), "{}", rows[0].0);
    let twice: Vec<u8> = [&genome[..], &genome].concat();
    let s = &rows[0].1;
    assert!(
        twice.windows(s.len()).any(|w| w == &s[..])
            || revcom(&twice).windows(s.len()).any(|w| w == &s[..])
    );
    // a second run finds all three steps in the cache
    let again = run(&p, items, Some(&cache), &Available::default()).unwrap();
    assert!(again.items[0]
        .steps
        .iter()
        .all(|s| matches!(s.status, Status::Cached)));
    // other long reads: the graph is kept, the joining is done again
    std::fs::write(long.join("sample1.fasta"), ">r\nACGTACGTACGT\n").unwrap();
    let other = run(
        &p,
        read_folder(&short, Kind::Reads).unwrap(),
        Some(&cache),
        &Available::default(),
    )
    .unwrap();
    let statuses: Vec<bool> = other.items[0]
        .steps
        .iter()
        .map(|s| matches!(s.status, Status::Cached))
        .collect();
    assert_eq!(statuses[..2], [true, true]);
    assert!(!statuses[2]);
    std::fs::remove_dir_all(short.parent().unwrap()).ok();
}

#[test]
fn the_nanopore_checks_collapse_a_concatemer_close_the_circle_and_turn_it_to_the_gene() {
    let d = folder("route2");
    let circle = random(3000, 5);
    // the assembler made two copies and a bit: the end of the contig is the start again. Cut in two, the unit of 3,030 bases ends with the
    // 30 bases that begin it
    let contig = [&circle[..], &circle, &circle[..60]].concat();
    let gene = circle[1000..1400].to_vec();
    let genes = d.join("start-genes.fasta");
    std::fs::write(
        &genes,
        format!(">dnaA\n{}\n", String::from_utf8_lossy(&gene)),
    )
    .unwrap();
    let text = format!(
        "helixloom-pipeline 1\ninput: sequences\nstep concatemer\nstep circularize\nstep rotate genes={}\n",
        genes.display()
    );
    let p = Pipeline::parse(&text).unwrap();
    let report = run(
        &p,
        vec![(
            "plasmid1".to_string(),
            Value::Sequences(vec![("contig_1".to_string(), contig)]),
        )],
        None,
        &Available::default(),
    )
    .unwrap();
    let Some(Value::Sequences(rows)) = &report.items[0].output else {
        panic!("{:?}", report.items[0].steps)
    };
    assert_eq!(
        rows[0].0,
        "contig_1 concatemer=2 circular=true rotated=true"
    );
    assert_eq!(rows[0].1.len(), 3000);
    assert_eq!(&rows[0].1[..400], &gene[..]);
    std::fs::remove_dir_all(&d).ok();
}

#[test]
fn a_step_that_calls_a_program_that_is_not_installed_fails_and_says_so() {
    if crate::pipeline::tools::detect(&crate::pipeline::tools::FLYE).is_some() {
        return;
    }
    let p = routes::nanopore_plasmid();
    let reads = Value::Reads(fastq_reads(&[random(3000, 6)]));
    let report = run(&p, vec![("a".into(), reads)], None, &Available::default()).unwrap();
    let first = &report.items[0].steps[0];
    let Status::Failed(why) = &first.status else {
        panic!("{:?}", first.status)
    };
    assert!(
        why.contains("Flye is not installed here") && why.contains("nothing else is used instead"),
        "{why}"
    );
}

#[test]
fn the_keys_of_two_steps_that_differ_only_in_a_setting_differ() {
    let a = step("step trim quality=20").unwrap();
    let b = step("step trim quality=21").unwrap();
    let input = Value::Reads(fastq_reads(&[random(100, 7)]));
    assert_ne!(key(&a, "", &input), key(&b, "", &input));
    assert_eq!(key(&a, "", &input), key(&a, "", &input));
    assert_ne!(key(&a, "|x", &input), key(&a, "|y", &input));
}
