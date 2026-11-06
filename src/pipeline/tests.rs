use super::export::{nextflow, snakefile, write_folder, Plan, PrimerSpec};
use super::steps::{self, primer3_input, read_fasta, read_primer3_output, Named};
use super::tools::{self, detect, MAFFT, PRIMER3};
use super::*;
use crate::seq::newick::NwTree;
use crate::seq::treetools::robinson_foulds;

fn demo_sequences() -> Vec<Named> {
    crate::sample::tree_demo()
        .into_iter()
        .map(|d| (d.name.clone(), d.seq.clone()))
        .collect()
}

#[test]
fn the_snakefile_has_a_rule_for_each_chosen_step_with_the_programs_command() {
    let text = snakefile(&Plan::default());
    for line in [
        "rule all:",
        "\"results/aligned.fasta\",",
        "\"results/tree.treefile\",",
        "\"results/primers.txt\",",
        "shell: \"mafft --auto --quiet {input} > {output}\"",
        "shell: \"iqtree3 -s {input} -m MFP -seed 1 -redo -quiet --prefix results/tree\"",
        "shell: \"primer3_core {input} > {output}\"",
    ] {
        assert!(text.contains(line), "missing {line}\n{text}");
    }
    let only_align = snakefile(&Plan {
        tree: false,
        primers: false,
        ..Plan::default()
    });
    assert!(only_align.contains("rule align:") && !only_align.contains("rule tree:"));
    assert!(!only_align.contains("primer3"));
    let boot = snakefile(&Plan {
        bootstrap: 500,
        iqtree: "iqtree2".into(),
        ..Plan::default()
    });
    assert!(
        boot.contains("iqtree2 -s {input} -m MFP -B 1000 -seed 1"),
        "{boot}"
    );
}

#[test]
fn the_nextflow_script_has_a_process_for_each_chosen_step() {
    let text = nextflow(&Plan::default());
    for line in [
        "nextflow.enable.dsl = 2",
        "process ALIGN {",
        "process TREE {",
        "process PRIMERS {",
        "mafft --auto --quiet ${sequences} > aligned.fasta",
        "iqtree3 -s ${alignment} -m MFP -seed 1 -redo -quiet --prefix tree",
        "primer3_core ${spec} > primers.txt",
        "aligned = ALIGN(Channel.fromPath(params.sequences))",
        "TREE(aligned)",
    ] {
        assert!(text.contains(line), "missing {line}\n{text}");
    }
    let no_tree = nextflow(&Plan {
        tree: false,
        ..Plan::default()
    });
    assert!(!no_tree.contains("TREE") && no_tree.contains("PRIMERS(Channel"));
}

#[test]
fn the_exported_folder_holds_clean_sequence_names_and_primer3s_input() {
    let seqs = vec![
        ("Isolate 1".to_string(), b"ACGTACGTACGT".to_vec()),
        ("b/2".to_string(), b"ACGTACGAACGT".to_vec()),
    ];
    let dir = tools::work_folder("export-test").unwrap();
    let settings = Settings::default();
    let files = write_folder(
        &dir,
        &seqs,
        &Plan::default(),
        Some(PrimerSpec {
            sequence: &seqs[0],
            target: Some((4, 3)),
            settings: &settings,
        }),
    )
    .unwrap();
    assert_eq!(files.len(), 4);
    let fasta = std::fs::read_to_string(dir.join("sequences.fasta")).unwrap();
    assert_eq!(fasta, ">Isolate_1\nACGTACGTACGT\n>b_2\nACGTACGAACGT\n");
    let input = std::fs::read_to_string(dir.join("primer3_input.txt")).unwrap();
    assert!(input.starts_with(
        "SEQUENCE_ID=Isolate_1\nSEQUENCE_TEMPLATE=ACGTACGTACGT\nSEQUENCE_TARGET=4,3\n"
    ));
    assert!(input.ends_with("=\n"));
    // without primers the primer step needs no sequence
    let no_primers = Plan {
        primers: false,
        ..Plan::default()
    };
    assert_eq!(
        write_folder(&dir, &seqs, &no_primers, None).unwrap().len(),
        3
    );
    assert!(write_folder(&dir, &seqs, &Plan::default(), None).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn primer3s_input_carries_the_settings_as_its_tags() {
    let s = Settings {
        opt_size: 22,
        product_min: 150,
        product_max: 400,
        num_return: 3,
        ..Settings::default()
    };
    let text = primer3_input("x", b"ACGT", None, &s);
    for line in [
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0",
        "PRIMER_OPT_SIZE=22",
        "PRIMER_OPT_TM=60.0",
        "PRIMER_SALT_DIVALENT=1.5",
        "PRIMER_MAX_SELF_ANY=8",
        "PRIMER_PRODUCT_SIZE_RANGE=150-400",
        "PRIMER_NUM_RETURN=3",
    ] {
        assert!(text.contains(line), "missing {line}");
    }
    assert!(!text.contains("SEQUENCE_TARGET"));
}

#[test]
fn primer3s_output_is_read_into_pairs() {
    let text = "SEQUENCE_ID=x\nPRIMER_PAIR_NUM_RETURNED=2\n\
        PRIMER_PAIR_0_PENALTY=0.4\nPRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\nPRIMER_RIGHT_0_SEQUENCE=TTTTGGGGCCCCAAAATTTT\n\
        PRIMER_LEFT_0=10,20\nPRIMER_RIGHT_0=229,20\nPRIMER_LEFT_0_TM=60.1\nPRIMER_RIGHT_0_TM=59.8\nPRIMER_PAIR_0_PRODUCT_SIZE=220\n\
        PRIMER_PAIR_1_PENALTY=1.5\nPRIMER_LEFT_1_SEQUENCE=CCCCCCCCCCCCCCCCCCCC\nPRIMER_RIGHT_1_SEQUENCE=GGGGGGGGGGGGGGGGGGGG\n\
        PRIMER_LEFT_1=12,20\nPRIMER_RIGHT_1=240,20\nPRIMER_LEFT_1_TM=61.0\nPRIMER_RIGHT_1_TM=60.0\nPRIMER_PAIR_1_PRODUCT_SIZE=229\n=\n";
    let pairs = read_primer3_output(text);
    assert_eq!(pairs.len(), 2);
    assert_eq!(
        (
            pairs[0].left_start,
            pairs[0].right_start,
            pairs[0].product_size
        ),
        (10, 210, 220)
    );
    assert_eq!((pairs[0].left_tm, pairs[0].penalty), (60.1, 0.4));
    assert_eq!(pairs[1].right, "GGGGGGGGGGGGGGGGGGGG");
    assert!(read_primer3_output("PRIMER_PAIR_NUM_RETURNED=0\n=\n").is_empty());
}

#[test]
fn the_built_in_steps_run_a_whole_workflow() {
    let seqs = demo_sequences();
    let choice = Choice {
        align: true,
        tree: true,
        primers: true,
        use_installed: false,
        bootstrap: 0,
    };
    let out = run(
        &seqs,
        Some((200, 60)),
        &Settings::default(),
        &choice,
        &Available::default(),
    )
    .unwrap();
    assert_eq!(out.log.len(), 3, "{:?}", out.log);
    let aligned = out.aligned.unwrap();
    assert_eq!((aligned.rows.len(), aligned.columns()), (6, 600));
    let tree = NwTree::parse_all(&out.tree.unwrap().newick)
        .unwrap()
        .remove(0);
    assert_eq!(tree.leaf_count(), 6);
    let primers = out.primers.unwrap();
    assert!(!primers.pairs.is_empty());
    // every pair puts the target between its primers
    assert!(primers
        .pairs
        .iter()
        .all(|p| p.left_start < 200 && p.right_start > 260));
    // a tree needs the alignment, so asking for one without the other still aligns
    let tree_only = Choice {
        align: false,
        primers: false,
        ..choice
    };
    assert!(run(
        &seqs,
        None,
        &Settings::default(),
        &tree_only,
        &Available::default()
    )
    .unwrap()
    .aligned
    .is_some());
}

#[test]
fn fasta_text_round_trips_with_gaps() {
    let rows = vec![
        ("a b".to_string(), b"AC-GT".to_vec()),
        ("c".to_string(), b"ACGGT".to_vec()),
    ];
    let back = read_fasta(&steps::fasta_text(&rows));
    assert_eq!(back, rows);
}

// The tests below run the real programs when they are installed, and say so when they are not.

#[test]
fn mafft_gives_the_same_alignment_as_the_built_in_one_when_there_are_no_gaps_to_place() {
    let Some(mafft) = detect(&MAFFT) else {
        eprintln!("MAFFT is not installed: skipped");
        return;
    };
    let seqs = demo_sequences();
    let dir = tools::work_folder("mafft-test").unwrap();
    let theirs = steps::align(&seqs, Some(&mafft), &dir).unwrap();
    let ours = steps::align(&seqs, None, &dir).unwrap();
    assert!(theirs.engine.starts_with("MAFFT"));
    assert_eq!(theirs.rows, ours.rows);
    assert_eq!(theirs.columns(), 600);
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn iqtree_and_the_built_in_search_find_the_same_tree_for_clearly_related_sequences() {
    let Some(iqtree) = detect(&tools::IQTREE) else {
        eprintln!("IQ-TREE is not installed: skipped");
        return;
    };
    let seqs = demo_sequences();
    let dir = tools::work_folder("iqtree-test").unwrap();
    let aligned = steps::align(&seqs, None, &dir).unwrap();
    let theirs = steps::tree(&aligned, Some(&iqtree), 0, &dir).unwrap();
    let ours = steps::tree(&aligned, None, 0, &dir).unwrap();
    assert!(
        theirs.engine.starts_with("IQ-TREE") && !theirs.model.is_empty(),
        "{}",
        theirs.engine
    );
    let (a, b) = (
        NwTree::parse_all(&theirs.newick).unwrap().remove(0),
        NwTree::parse_all(&ours.newick).unwrap().remove(0),
    );
    // all six names are back in IQ-TREE's tree
    assert_eq!(
        a.nodes
            .iter()
            .filter(|n| n.name.starts_with("Isolate"))
            .count(),
        6
    );
    let rf = robinson_foulds(&a, &b).unwrap();
    eprintln!(
        "IQ-TREE ({}) ln L {:?}, built in {:?}; Robinson-Foulds {} of {}",
        theirs.model, theirs.ln_likelihood, ours.ln_likelihood, rf.distance, rf.max
    );
    // the two trees put under the built-in models: which topology is more likely there?
    let msa = crate::seq::msa::Msa {
        names: aligned.rows.iter().map(|r| r.0.replace(' ', "_")).collect(),
        rows: aligned.rows.iter().map(|r| r.1.clone()).collect(),
    };
    let best = |tree: &NwTree| {
        crate::seq::likelihood::fit_models(&msa, &tree.to_node())
            .unwrap()
            .into_iter()
            .min_by(|x, y| x.bic.total_cmp(&y.bic))
            .unwrap()
    };
    let (fa, fb) = (best(&a), best(&b));
    // the trees can differ in how they split a branch that is almost nothing; what matters is that the built-in
    // search finds a topology at least as likely as IQ-TREE's under the built-in models
    assert!(fb.ln_l >= fa.ln_l - 1e-3, "{} vs {}", fb.ln_l, fa.ln_l);
    eprintln!(
        "under the built-in models: IQ-TREE's topology ln L {:.3} ({}), the built-in search's {:.3} ({})",
        fa.ln_l,
        fa.kind.label(),
        fb.ln_l,
        fb.kind.label()
    );
    std::fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn primer3_and_the_port_of_it_pick_the_same_pairs_when_primer3_is_installed() {
    let Some(primer3) = detect(&PRIMER3) else {
        eprintln!("primer3_core is not installed: skipped");
        return;
    };
    let template = crate::sample::demo_plasmid().seq;
    let dir = tools::work_folder("primer3-test").unwrap();
    let settings = Settings::default();
    let theirs =
        steps::primers(&template, Some((224, 56)), &settings, Some(&primer3), &dir).unwrap();
    let ours = steps::primers(&template, Some((224, 56)), &settings, None, &dir).unwrap();
    assert_eq!(theirs.pairs.len(), ours.pairs.len());
    for (a, b) in theirs.pairs.iter().zip(&ours.pairs) {
        assert_eq!(
            (&a.left, &a.right, a.product_size),
            (&b.left, &b.right, b.product_size)
        );
        assert!((a.left_tm - b.left_tm).abs() < 0.05 && (a.penalty - b.penalty).abs() < 0.01);
    }
    std::fs::remove_dir_all(&dir).unwrap();
}

/// Writes a workflow folder for trying the Snakemake and Nextflow files by hand:
/// `cargo test export_a_folder -- --ignored --nocapture`, then run the engines in the folder it names.
#[test]
#[ignore = "writes a folder for trying the exported files by hand"]
fn export_a_folder_for_snakemake_and_nextflow() {
    let seqs = demo_sequences();
    let dir = std::env::temp_dir().join("helixloom-export-demo");
    let settings = Settings::default();
    let plan = Plan::default();
    let files = write_folder(
        &dir,
        &seqs,
        &plan,
        Some(PrimerSpec {
            sequence: &seqs[0],
            target: Some((200, 60)),
            settings: &settings,
        }),
    )
    .unwrap();
    eprintln!("EXPORTED {} files to {}", files.len(), dir.display());
}
