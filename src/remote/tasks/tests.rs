use super::*;
use crate::remote::run::{execute, Env};

fn isolates() -> Vec<Named> {
    crate::sample::tree_demo()
        .into_iter()
        .map(|d| (d.name, d.seq))
        .collect()
}

/// A request of every kind, with settings that are not the defaults.
fn requests() -> Vec<Request> {
    // every field of the primer settings changed, so one the codec forgets shows
    let settings = Settings {
        opt_size: 21,
        min_size: 19,
        max_size: 26,
        opt_tm: 61.5,
        min_tm: 58.25,
        max_tm: 64.0,
        min_gc: 25.0,
        max_gc: 75.5,
        salt_mm: 55.0,
        divalent_mm: 2.0,
        dntp_mm: 0.7,
        dna_nm: 40.0,
        max_poly_x: 4,
        gc_clamp: 1,
        max_end_gc: 3,
        max_self_any: 7,
        max_self_end: 2,
        pair_compl_any: 9,
        pair_compl_end: 4,
        max_diff_tm: 3.5,
        product_min: 90,
        product_max: 410,
        num_return: 7,
    };
    vec![
        Request::Align {
            seqs: isolates(),
            mafft: true,
        },
        Request::Tree {
            alignment: vec![
                ("a".into(), b"AC-GT".to_vec()),
                ("b b".into(), b"ACCGT".to_vec()),
            ],
            analysis: Analysis::Bayesian,
            model: DistanceModel::JukesCantor,
            bootstrap: 250,
            generations: 20_000,
        },
        Request::Pipeline {
            seqs: isolates(),
            target: Some((200, 60)),
            settings,
            choice: Choice {
                align: true,
                tree: false,
                primers: true,
                use_installed: false,
                bootstrap: 1000,
            },
        },
        Request::PhyloCnn {
            newick: "((A:1,B:1):2,C:3);\n((A:1,B:2):1,C:3);".into(),
            sampling: 0.125,
            task: Task::Bdei,
            priors: Priors::Simulated,
            gpu: true,
        },
        Request::Fold {
            sequence: b"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQ".to_vec(),
            local: true,
        },
    ]
}

#[test]
fn every_kind_of_request_comes_back_as_it_went() {
    for request in requests() {
        let kind = request.kind();
        let bytes = request.encode().encode();
        let back = Request::decode(kind, &Container::decode(&bytes).expect("a container"))
            .expect("a request");
        assert_eq!(back, request, "{kind:?}");
    }
    assert_eq!(
        Kind::ALL.map(Kind::name),
        ["align", "tree", "pipeline", "phylocnn", "fold"]
    );
    for k in Kind::ALL {
        assert_eq!(Kind::parse(k.name()), Some(k));
    }
    assert_eq!(Kind::parse("rm -rf"), None);
}

#[test]
fn only_a_gpu_job_says_it_needs_a_gpu() {
    let needs: Vec<bool> = requests().iter().map(Request::needs_gpu).collect();
    // the PhyloCNN request has two trees, and a build without the GPU code has no use for one
    assert_eq!(needs, [false, false, false, cfg!(feature = "gpu"), true]);
    let cnn = |newick: &str, gpu: bool| Request::PhyloCnn {
        newick: newick.into(),
        sampling: 0.5,
        task: Task::Bd,
        priors: Priors::Simulated,
        gpu,
    };
    assert!(!cnn("(A:1,B:1);\n(A:1,B:2);", false).needs_gpu());
    // one tree runs on the processor whatever was asked, so it does not hold a GPU
    assert!(!cnn("(A:1,B:1);", true).needs_gpu());
}

#[test]
fn names_that_would_break_a_fasta_line_are_made_one_line() {
    let request = Request::Align {
        seqs: vec![
            ("two\nlines\t>x".into(), b"ACGT".to_vec()),
            ("plain".into(), b"ACGA".to_vec()),
        ],
        mafft: false,
    };
    let Ok(Request::Align { seqs, .. }) = Request::decode(Kind::Align, &request.encode()) else {
        panic!("a request");
    };
    assert_eq!(seqs.len(), 2);
    assert_eq!(seqs[0].0, "two lines >x");
    assert_eq!(seqs[1], ("plain".to_string(), b"ACGA".to_vec()));
    // a long name is not cut, so that two long names that begin alike do not become one
    let long = |tail: char| format!("{}{tail}", "a long FASTA header ".repeat(20));
    let request = Request::Align {
        seqs: vec![(long('1'), b"ACGT".to_vec()), (long('2'), b"ACGA".to_vec())],
        mafft: false,
    };
    let Ok(Request::Align { seqs, .. }) = Request::decode(Kind::Align, &request.encode()) else {
        panic!("a request");
    };
    assert_eq!(
        (seqs[0].0.clone(), seqs[1].0.clone()),
        (long('1'), long('2'))
    );
}

#[test]
fn requests_that_are_wrong_are_refused_in_words() {
    let good = Request::Tree {
        alignment: vec![
            ("a".into(), b"ACGT".to_vec()),
            ("b".into(), b"ACGA".to_vec()),
        ],
        analysis: Analysis::NeighborJoining,
        model: DistanceModel::P,
        bootstrap: 10,
        generations: 1000,
    }
    .encode();
    assert!(Request::decode(Kind::Tree, &good).is_ok());
    let with = |meta: &str| {
        Container::new()
            .with("meta", meta)
            .with("alignment", ">a\nACGT\n>b\nACGA\n")
    };
    let ok = "analysis=nj\nmodel=p\nbootstrap=10\ngenerations=1000\n";
    assert!(Request::decode(Kind::Tree, &with(ok)).is_ok());
    for (what, meta) in [
        (
            "method",
            "analysis=parsimony\nmodel=p\nbootstrap=10\ngenerations=1000\n",
        ),
        (
            "model",
            "analysis=nj\nmodel=x\nbootstrap=10\ngenerations=1000\n",
        ),
        (
            "bootstrap",
            "analysis=nj\nmodel=p\nbootstrap=ten\ngenerations=1000\n",
        ),
        (
            "too many resamples",
            "analysis=nj\nmodel=p\nbootstrap=99999\ngenerations=1000\n",
        ),
        (
            "too many generations",
            "analysis=nj\nmodel=p\nbootstrap=1\ngenerations=99999999999\n",
        ),
        ("no generations", "analysis=nj\nmodel=p\nbootstrap=1\n"),
    ] {
        assert!(Request::decode(Kind::Tree, &with(meta)).is_err(), "{what}");
    }
    // sections that are missing, or hold no sequences
    assert!(Request::decode(Kind::Tree, &Container::new().with("meta", ok)).is_err());
    assert!(Request::decode(
        Kind::Tree,
        &Container::new().with("meta", ok).with("alignment", "")
    )
    .is_err());
    assert!(Request::decode(
        Kind::Tree,
        &Container::new()
            .with("meta", ok)
            .with("alignment", ">a\n>b\nAC\n")
    )
    .is_err());
    // rows of different lengths are not an alignment
    let ragged = Container::new()
        .with("meta", ok)
        .with("alignment", ">a\nACGT\n>b\nAC\n");
    assert!(Request::decode(Kind::Tree, &ragged)
        .unwrap_err()
        .contains("same length"));
    assert!(Request::decode(Kind::Align, &Container::new().with("sequences", ">a\nAC\n")).is_err());
    // protein that is not protein, and a sequence too short or too long for the service
    let fold = |s: &str| {
        Container::new()
            .with("meta", "local=0\n")
            .with("sequence", s)
    };
    assert!(Request::decode(Kind::Fold, &fold("MQIFVKTLTGKTITLEVEP")).is_ok());
    assert!(Request::decode(Kind::Fold, &fold("MQIFVKTLTGKTITLEVEP1X")).is_err());
    assert!(Request::decode(Kind::Fold, &fold("MQI")).is_err());
    assert!(Request::decode(Kind::Fold, &fold(&"A".repeat(401))).is_err());
    // a PhyloCNN question that does not exist
    let cnn = |task: &str, priors: &str| {
        Container::new()
            .with(
                "meta",
                format!("sampling=0.5\ntask={task}\npriors={priors}\ngpu=0\n"),
            )
            .with("trees", "(A:1,B:1);")
    };
    assert!(Request::decode(Kind::PhyloCnn, &cnn("bd", "simulated")).is_ok());
    assert!(Request::decode(Kind::PhyloCnn, &cnn("xx", "simulated")).is_err());
    assert!(Request::decode(Kind::PhyloCnn, &cnn("bd", "none")).is_err());
    // a pipeline whose primer target or settings are not numbers
    let sent = requests()[2].encode();
    let section = |name: &str| sent.text(name).expect("a section").to_string();
    let changed = |meta: String, settings: String| {
        Container::new()
            .with("meta", meta)
            .with("settings", settings)
            .with("sequences", section("sequences"))
    };
    assert!(Request::decode(
        Kind::Pipeline,
        &changed(section("meta"), section("settings"))
    )
    .is_ok());
    // the edits below only mean something if the text they change is really there
    assert!(
        section("meta").contains("target=200,60") && section("settings").contains("opt_size=21")
    );
    let target = section("meta").replace("target=200,60", "target=200");
    assert!(Request::decode(Kind::Pipeline, &changed(target, section("settings"))).is_err());
    let size = section("settings").replace("opt_size=21", "opt_size=big");
    assert!(Request::decode(Kind::Pipeline, &changed(section("meta"), size)).is_err());
}

#[test]
fn settings_left_out_keep_primer3s_defaults() {
    assert_eq!(settings_parse("").expect("defaults"), Settings::default());
    let s = settings_parse("opt_size=22\nnot_a_setting=5\n").expect("settings");
    assert_eq!(s.opt_size, 22);
    assert_eq!(s.min_size, Settings::default().min_size);
    assert_eq!(
        settings_parse(&settings_text(&Settings::default())).expect("settings"),
        Settings::default()
    );
}

/// A response goes through bytes and comes out as the same container.
fn through_bytes(response: &Response) -> Response {
    let kind = match response {
        Response::Align(_) => Kind::Align,
        Response::Tree(_) => Kind::Tree,
        Response::Pipeline(_) => Kind::Pipeline,
        Response::PhyloCnn(_) => Kind::PhyloCnn,
        Response::Fold(_) => Kind::Fold,
    };
    let bytes = response.encode().encode();
    let back = Response::decode(kind, &Container::decode(&bytes).expect("a container"))
        .expect("a response");
    assert_eq!(
        back.encode(),
        response.encode(),
        "{kind:?} encodes the same after a round trip"
    );
    back
}

#[test]
fn an_alignment_a_tree_and_a_pipeline_come_back_as_the_server_made_them() {
    let env = Env::default();
    let Response::Align(aligned) = execute(
        Request::Align {
            seqs: isolates(),
            mafft: false,
        },
        &env,
    )
    .expect("aligned") else {
        panic!("an alignment");
    };
    let Response::Align(back) = through_bytes(&Response::Align(aligned.clone())) else {
        panic!("an alignment");
    };
    assert_eq!(
        (&back.rows, &back.engine, back.seconds),
        (&aligned.rows, &aligned.engine, aligned.seconds)
    );

    let tree = execute(
        Request::Tree {
            alignment: aligned.rows.clone(),
            analysis: Analysis::NeighborJoining,
            model: DistanceModel::Kimura2,
            bootstrap: 20,
            generations: 1000,
        },
        &env,
    )
    .expect("a tree");
    let Response::Tree(back) = through_bytes(&tree) else {
        panic!("a tree");
    };
    let Response::Tree(made) = tree else {
        panic!("a tree")
    };
    // the very same tree, not only the same picture: names, lengths and supports to the last digit
    assert_eq!(back.tree, made.tree);
    assert_eq!(back.details, made.details);
    assert!(back.details.iter().any(|d| d.contains("bootstrap")));

    let choice = Choice {
        align: true,
        tree: true,
        primers: true,
        use_installed: false,
        bootstrap: 0,
    };
    let run = execute(
        Request::Pipeline {
            seqs: isolates(),
            target: Some((200, 60)),
            settings: Settings::default(),
            choice,
        },
        &env,
    )
    .expect("a pipeline");
    let Response::Pipeline(back) = through_bytes(&run) else {
        panic!("a pipeline");
    };
    let Response::Pipeline(made) = run else {
        panic!("a pipeline")
    };
    assert_eq!(back.log, made.log);
    assert_eq!(
        back.aligned.as_ref().map(|a| &a.rows),
        made.aligned.as_ref().map(|a| &a.rows)
    );
    assert_eq!(
        back.tree.as_ref().map(|t| &t.newick),
        made.tree.as_ref().map(|t| &t.newick)
    );
    assert_eq!(
        back.tree.as_ref().map(|t| t.ln_likelihood),
        made.tree.as_ref().map(|t| t.ln_likelihood)
    );
    assert_eq!(
        back.primers.as_ref().map(|p| &p.pairs),
        made.primers.as_ref().map(|p| &p.pairs)
    );
    assert!(!back.primers.expect("primers").pairs.is_empty());
}

#[test]
fn phylocnn_estimates_come_back_whole() {
    if crate::models::phylocnn_missing() {
        return;
    }
    let Response::PhyloCnn(made) = execute(
        Request::PhyloCnn {
            newick: crate::seq::phylocnn::example_tree().0.trim().to_string(),
            sampling: 0.25,
            task: Task::Bdss,
            priors: Priors::Hiv,
            gpu: false,
        },
        &Env::default(),
    )
    .expect("estimates") else {
        panic!("estimates");
    };
    assert_eq!((made.trees, made.tips, made.per_tree.len()), (1, 200, 1));
    assert_eq!(made.per_tree[0].len(), 4);
    let Response::PhyloCnn(back) = through_bytes(&Response::PhyloCnn(made.clone())) else {
        panic!("estimates");
    };
    assert_eq!(back, made);
    // trees that were left out are reported with their reasons
    let skipped = CnnResult {
        skipped: vec![(2, "not binary".into()), (5, "too many tips".into())],
        ..made
    };
    let Response::PhyloCnn(back) = through_bytes(&Response::PhyloCnn(skipped.clone())) else {
        panic!("estimates");
    };
    assert_eq!(back.skipped, skipped.skipped);
}

#[test]
fn a_structure_comes_back_whole() {
    let fold = Folded {
        pdb: "ATOM      1  N   MET A   1\nEND\n".into(),
        engine: "a program".into(),
        seconds: 12.5,
    };
    let Response::Fold(back) = through_bytes(&Response::Fold(fold.clone())) else {
        panic!("a structure");
    };
    assert_eq!(back, fold);
}

#[test]
fn an_answer_is_taken_only_as_the_kind_it_is() {
    let fold = Folded {
        pdb: "ATOM\n".into(),
        engine: "x".into(),
        seconds: 1.0,
    };
    assert_eq!(Response::Fold(fold.clone()).into_fold(), Ok(fold));
    let aligned = Aligned {
        rows: vec![("a".into(), b"AC".to_vec())],
        engine: "e".into(),
        seconds: 0.5,
    };
    assert_eq!(
        Response::Align(aligned).into_align().map(|a| a.rows.len()),
        Ok(1)
    );
    let why = Response::Fold(Folded {
        pdb: String::new(),
        engine: String::new(),
        seconds: 0.0,
    })
    .into_align()
    .err()
    .expect("a different answer");
    assert!(why.contains("an alignment"), "{why}");
    assert!(Response::Align(Aligned {
        rows: Vec::new(),
        engine: String::new(),
        seconds: 0.0
    })
    .into_tree()
    .is_err());
}

#[test]
fn damaged_responses_are_refused_not_trusted() {
    let t = Container::new().with("meta", "task=bd\nseconds=1\nengine=CPU\ntrees=1\ntips=3\n");
    // wrong number of values for the question, a value that is not a number, a row with no tab
    for rows in ["0\t1.5\n", "0\t1.5,x\n", "0 1.5,2\n"] {
        let c = t.clone().with("estimates", rows).with("skipped", "");
        assert!(Response::decode(Kind::PhyloCnn, &c).is_err(), "{rows:?}");
    }
    let ok = t
        .clone()
        .with("estimates", "0\t1.5,2\n")
        .with("skipped", "");
    assert!(Response::decode(Kind::PhyloCnn, &ok).is_ok());
    assert!(Response::decode(
        Kind::Tree,
        &Container::new()
            .with("newick", "not a tree")
            .with("details", "")
    )
    .is_err());
    assert!(Response::decode(
        Kind::Align,
        &Container::new()
            .with("meta", "seconds=1\n")
            .with("rows", "")
    )
    .is_err());
    let pairs = Container::new()
        .with("log", "")
        .with("primers_meta", "seconds=1\n")
        .with("primers", "A\tB\t1\t2\t3\n");
    assert!(Response::decode(Kind::Pipeline, &pairs).is_err());
}

#[test]
fn a_request_may_not_ask_a_job_to_hold_more_than_a_few_gigabytes() {
    let rows = |n: usize| -> String { (0..n).map(|i| format!(">s{i}\nACGT\n")).collect() };
    let align = |text: String| {
        Container::new()
            .with("meta", "engine=builtin\n")
            .with("sequences", text)
    };
    assert!(Request::decode(Kind::Align, &align(rows(MAX_SEQUENCES))).is_ok());
    let why = Request::decode(Kind::Align, &align(rows(MAX_SEQUENCES + 1))).expect_err("too many");
    assert!(why.contains("at most 5000 sequences"), "{why}");
    // the same for the alignment of a tree and for a pipeline
    let tree = |text: String| {
        Container::new()
            .with(
                "meta",
                "analysis=nj\nmodel=p\nbootstrap=0\ngenerations=1000\n",
            )
            .with("alignment", text)
    };
    assert!(Request::decode(Kind::Tree, &tree(rows(MAX_SEQUENCES))).is_ok());
    assert!(Request::decode(Kind::Tree, &tree(rows(MAX_SEQUENCES + 1))).is_err());
    let sent = requests()[2].encode();
    let pipeline = |text: String| {
        Container::new()
            .with("meta", sent.text("meta").expect("meta"))
            .with("settings", sent.text("settings").expect("settings"))
            .with("sequences", text)
    };
    assert!(Request::decode(Kind::Pipeline, &pipeline(rows(MAX_SEQUENCES))).is_ok());
    assert!(Request::decode(Kind::Pipeline, &pipeline(rows(MAX_SEQUENCES + 1))).is_err());

    // too many letters in all, however few the sequences
    let long = format!(">a\n{}\n", "A".repeat(MAX_RESIDUES + 1));
    let why = Request::decode(Kind::Align, &align(long)).expect_err("too long");
    assert!(why.contains("letters of sequence"), "{why}");

    // trees are counted by their semicolons, without being parsed
    let trees = |n: usize| {
        Container::new()
            .with("meta", "sampling=0.5\ntask=bd\npriors=simulated\ngpu=0\n")
            .with("trees", "(A:1,B:1);".repeat(n))
    };
    assert!(Request::decode(Kind::PhyloCnn, &trees(MAX_TREES)).is_ok());
    let why = Request::decode(Kind::PhyloCnn, &trees(MAX_TREES + 1)).expect_err("too many trees");
    assert!(why.contains("at most 20000 trees"), "{why}");
    // a posterior sample of BEAST's own default size (10,001 trees) is not turned away
    assert!(Request::decode(Kind::PhyloCnn, &trees(10_001)).is_ok());
}

#[test]
fn a_tree_comes_back_exactly_whatever_its_names_lengths_and_supports() {
    use crate::seq::tree::Node;
    let leaf = |name: &str, length: f64| Node {
        name: name.into(),
        length,
        support: None,
        children: Vec::new(),
    };
    let inner = |support: Option<f64>, length: f64, children: Vec<Node>| Node {
        name: String::new(),
        length,
        support,
        children,
    };
    let tree = inner(
        None,
        0.0,
        vec![
            inner(
                Some(87.5),
                0.000_000_1,
                vec![leaf("Isolate 1", 0.1 + 0.2), leaf("a'b (c): d,e;[f]", 1e-9)],
            ),
            inner(
                Some(100.0),
                12_345.678_901_234_5,
                vec![leaf("plain", 0.0), leaf("é ü ñ", 2.5)],
            ),
            leaf("123", 3.0),
        ],
    );
    let back = Response::Tree(crate::seq::estimate::Estimate {
        tree: tree.clone(),
        details: vec!["made here".into()],
    });
    let Response::Tree(back) = through_bytes(&back) else {
        panic!("a tree");
    };
    assert_eq!(back.tree, tree);
    // the text for a file is still plain and short
    assert_eq!(leaf("a b", 1.0).newick(), "a_b;");
}
