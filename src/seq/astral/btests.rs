//! Bootstrap support: what the replicates are, and that the support is the share of them that have the branch.

use super::bootstrap::{Boot, Resample};
use super::{infer, Options};

fn options(replicates: usize, resample: Resample, texts: Vec<String>) -> Options {
    Options {
        bootstrap: Boot {
            replicates,
            seed: 7,
            resample,
        },
        site_trees: texts,
        ..Options::default()
    }
}

/// The support on the branch of the one inner branch of a tree of four species (that of C and D, the side without A).
fn only_support(text: &str, o: &Options) -> f64 {
    let out = infer(text, o).expect("a species tree");
    assert_eq!(out.branches.len(), 1, "{}", out.newick);
    out.branches[0].bootstrap.expect("a support")
}

const AB: &str = "((A,B),(C,D));\n";
const AC: &str = "((A,C),(B,D));\n";

#[test]
fn genes_that_all_say_the_same_give_every_branch_full_support() {
    let text = AB.repeat(30);
    let out = infer(&text, &options(40, Resample::Genes, Vec::new())).expect("a tree");
    assert_eq!(out.branches[0].bootstrap, Some(100.0));
    let info = out.bootstrap.expect("what the bootstrap did");
    assert_eq!((info.asked, info.done, info.failed), (40, 40, 0));
    // the label on the branch is the support, and the log says what it is
    assert!(out.newick.contains(")100:"), "{}", out.newick);
    assert!(
        out.log
            .iter()
            .any(|l| l.contains("Bootstrap: 40 replicates resampling the genes")),
        "{:?}",
        out.log
    );
}

#[test]
fn genes_that_are_split_give_a_support_between_none_and_full_that_is_the_same_for_the_same_seed() {
    // forty genes join A with B and thirty-five A with C: the main tree has A and B together, and about half of the
    // resamplings of the genes have more of the others
    let text = format!("{}{}", AB.repeat(40), AC.repeat(35));
    let first = only_support(&text, &options(200, Resample::Genes, Vec::new()));
    assert!(first > 20.0 && first < 95.0, "{first}");
    assert_eq!(
        first,
        only_support(&text, &options(200, Resample::Genes, Vec::new()))
    );
    // another seed draws other genes
    let mut other = options(200, Resample::Genes, Vec::new());
    other.bootstrap.seed = 8;
    let second = only_support(&text, &other);
    assert!(second > 20.0 && second < 95.0, "{second}");
    // the support is a share of the replicates that gave a tree, as many places as 1/200 of them
    assert!(
        ((first * 2.0).round() - first * 2.0).abs() < 1e-9,
        "{first}"
    );
}

/// A gene's bootstrap trees: `first` of them say A with B, the rest A with C.
fn bootstrap_trees(total: usize, first: usize) -> String {
    format!("{}{}", AB.repeat(first), AC.repeat(total - first))
}

#[test]
fn resampling_the_sites_takes_the_same_bootstrap_tree_of_every_gene_for_a_replicate() {
    // five genes, each with twenty bootstrap trees, the first ten saying A with B and the last ten A with C: the replicates
    // 0 to 9 give the tree with A and B together and 10 to 19 the other, so half of them have the branch of the main tree
    let genes = AB.repeat(5);
    let sites: Vec<String> = (0..5).map(|_| bootstrap_trees(20, 10)).collect();
    let o = options(20, Resample::Sites, sites);
    assert_eq!(only_support(&genes, &o), 50.0);
    // fewer replicates, all in the first half: every one has it
    let o = options(
        10,
        Resample::Sites,
        (0..5).map(|_| bootstrap_trees(20, 10)).collect(),
    );
    assert_eq!(only_support(&genes, &o), 100.0);
}

#[test]
fn resampling_genes_and_sites_draws_both_and_gives_a_support_in_between() {
    let genes = AB.repeat(6);
    let sites: Vec<String> = (0..6).map(|_| bootstrap_trees(20, 12)).collect();
    let both = only_support(&genes, &options(200, Resample::Both, sites.clone()));
    assert!(both > 20.0 && both < 100.0, "{both}");
    assert_eq!(
        both,
        only_support(&genes, &options(200, Resample::Both, sites))
    );
}

#[test]
fn what_the_bootstrap_cannot_use_is_refused_with_the_reason() {
    let genes = AB.repeat(4);
    let why = |o: Options| infer(&genes, &o).expect_err("refused");
    // sites need bootstrap trees, one text for each gene
    assert!(why(options(5, Resample::Sites, Vec::new())).contains("4 gene trees"));
    let three: Vec<String> = (0..3).map(|_| bootstrap_trees(10, 5)).collect();
    assert!(why(options(5, Resample::Both, three)).contains("3 genes"));
    // a gene with fewer bootstrap trees than replicates, for the sites alone
    let few: Vec<String> = (0..4).map(|_| bootstrap_trees(3, 2)).collect();
    let said = why(options(5, Resample::Sites, few.clone()));
    assert!(said.contains("Gene 1 has 3 bootstrap trees"), "{said}");
    // genes and sites draw from them, so that few is enough
    assert!(infer(&genes, &options(5, Resample::Both, few)).is_ok());
    // a species the gene trees do not have, and a gene with no bootstrap tree
    let odd: Vec<String> = (0..4).map(|_| "((A,B),(C,E));\n".to_string()).collect();
    assert!(why(options(2, Resample::Sites, odd)).contains("species E"));
    let none: Vec<String> = vec![String::new(); 4];
    assert!(why(options(2, Resample::Sites, none)).contains("No tree"));
    // the number of replicates, a given tree, and what is not built together with it
    assert!(why(options(1_001, Resample::Genes, Vec::new())).contains("from 1 to 1000"));
    // weights go with resampling the genes, whose supports and lengths are in the gene trees, but not with bootstrap trees, which
    // have none to weigh by
    let mut weighted = options(5, Resample::Both, bootstrap_trees_for(4));
    weighted.weights = super::Weights::Support;
    let said =
        infer(&"((A:1,B:1)0.9:1,(C:1,D:1)0.9:1);\n".repeat(4), &weighted).expect_err("refused");
    assert!(
        said.contains("cannot be used together") && said.contains("bootstrap trees"),
        "{said}"
    );
    let mut families = options(5, Resample::Genes, Vec::new());
    families.families = true;
    assert!(infer(&genes, &families)
        .expect_err("refused")
        .contains("not built"));
    assert!(super::score_tree(
        &genes,
        "((A,B),(C,D));",
        &options(5, Resample::Genes, Vec::new())
    )
    .expect_err("refused")
    .contains("only scored"));
}

#[test]
fn genes_that_were_left_out_as_too_small_take_their_bootstrap_trees_with_them() {
    // the second gene has three species and is left out; the bootstrap trees are given for all three genes, in order
    let genes = format!("{AB}((A,B),C);\n{AB}");
    let sites = vec![
        bootstrap_trees(6, 6),
        "((A,B),C);\n".repeat(6),
        bootstrap_trees(6, 6),
    ];
    let o = options(6, Resample::Sites, sites);
    assert_eq!(only_support(&genes, &o), 100.0);
}

#[test]
fn a_run_without_replicates_has_no_support() {
    let out = infer(&AB.repeat(5), &Options::default()).expect("a tree");
    assert!(out.bootstrap.is_none() && out.branches.iter().all(|b| b.bootstrap.is_none()));
    assert!(!out.log.iter().any(|l| l.contains("Bootstrap")));
}

#[test]
fn the_names_of_the_ways_to_resample_are_read_back() {
    for r in Resample::ALL {
        assert_eq!(Resample::parse(r.name()), Some(r));
    }
    assert_eq!(Resample::parse(" Both "), Some(Resample::Both));
    assert_eq!(Resample::parse("sample"), None);
    assert_eq!(Resample::default(), Resample::Genes);
}

#[test]
fn a_list_of_bootstrap_files_is_read_with_paths_taken_from_its_folder() {
    use super::bootstrap::read_list;
    let dir = crate::pipeline::tools::work_folder("boot-list").expect("a folder");
    let write = |name: &str, text: &str| std::fs::write(dir.join(name), text).expect("a file");
    write("g1.tre", "((A,B),(C,D));\n");
    write("g2.tre", "((A,C),(B,D));\n");
    write("list.txt", "# the genes, in order\n\ng1.tre\n  g2.tre  \n");
    let texts = read_list(&dir.join("list.txt")).expect("the texts");
    assert_eq!(texts, ["((A,B),(C,D));\n", "((A,C),(B,D));\n"]);
    // a list with nothing in it, a file that is not there, and a list that is not there
    write("empty.txt", "# nothing\n\n");
    assert!(read_list(&dir.join("empty.txt"))
        .unwrap_err()
        .contains("names no file"));
    write("missing.txt", "g1.tre\nnot-there.tre\n");
    assert!(read_list(&dir.join("missing.txt"))
        .unwrap_err()
        .contains("not-there.tre"));
    assert!(read_list(&dir.join("nope.txt"))
        .unwrap_err()
        .contains("Could not read the list"));
    let _ = std::fs::remove_dir_all(&dir);
}

/// Bootstrap trees for `genes` genes, as the sites need them.
fn bootstrap_trees_for(genes: usize) -> Vec<String> {
    (0..genes).map(|_| bootstrap_trees(6, 6)).collect()
}

/// A gene tree of four species that joins A with B (`ab`) or A with C, every branch with this support and length 1.
fn weighted_gene(ab: bool, support: &str) -> String {
    if ab {
        format!("((A:1,B:1){support}:1,(C:1,D:1){support}:1);\n")
    } else {
        format!("((A:1,C:1){support}:1,(B:1,D:1){support}:1);\n")
    }
}

#[test]
fn a_bootstrap_with_weights_resamples_the_genes_with_their_supports_and_lengths() {
    // forty genes that join A with B, well supported, and thirty-five that join A with C, with a support of 0: they weigh nothing,
    // so every resampling of the genes still has A with B (the unweighted bootstrap of the same genes does not: its support is
    // the share of resamplings in which forty of the genes outnumber thirty-five)
    let text = format!(
        "{}{}",
        weighted_gene(true, "0.9").repeat(40),
        weighted_gene(false, "0").repeat(35)
    );
    let mut weighted = options(200, Resample::Genes, Vec::new());
    weighted.weights = super::Weights::Support;
    let out = infer(&text, &weighted).expect("a species tree");
    assert_eq!(out.branches.len(), 1);
    assert_eq!(out.branches[0].bootstrap, Some(100.0), "{}", out.newick);
    assert_eq!(out.weights, super::Weights::Support);
    let plain = only_support(&text, &options(200, Resample::Genes, Vec::new()));
    assert!(
        plain < 99.0,
        "the unweighted support of the same genes is {plain}"
    );
    // the same seed gives the same support, and the genes are drawn as they are without weights
    let again = infer(&text, &weighted).expect("a species tree");
    assert_eq!(again.branches[0].bootstrap, out.branches[0].bootstrap);
    // the other way round: genes that weigh the same give a support between none and all, as the unweighted ones do
    let even = format!(
        "{}{}",
        weighted_gene(true, "0.9").repeat(40),
        weighted_gene(false, "0.9").repeat(35)
    );
    let support = infer(&even, &weighted).expect("a species tree").branches[0]
        .bootstrap
        .expect("a support");
    assert!(support > 20.0 && support < 95.0, "{support}");
}

#[test]
fn the_replicate_trees_and_their_greedy_consensus_are_kept_only_when_asked_for() {
    let text = AB.repeat(30);
    let mut keep = options(40, Resample::Genes, Vec::new());
    keep.keep_replicates = true;
    let out = infer(&text, &keep).expect("a species tree");
    assert_eq!(out.replicate_trees.len(), 40);
    assert!(
        out.replicate_trees
            .iter()
            .all(|t| t == &out.replicate_trees[0]),
        "{:?}",
        out.replicate_trees
    );
    let first = &out.replicate_trees[0];
    assert!(
        first.ends_with(");") && (first.contains("(C,D)") || first.contains("(D,C)")),
        "{first}"
    );
    // all the trees agree, so the consensus is that tree with every branch at 100
    assert_eq!(out.consensus.as_deref(), Some("(A,(B,(C,D)100));"));
    // not asked for: nothing is kept
    let plain = infer(&text, &options(40, Resample::Genes, Vec::new())).expect("a species tree");
    assert!(plain.replicate_trees.is_empty() && plain.consensus.is_none());
    // the trees are of the bootstrap: without replicates there are none
    let none = Options {
        keep_replicates: true,
        ..Options::default()
    };
    let said = infer(&text, &none).expect_err("refused");
    assert!(
        said.contains("replicate trees") && said.contains("bootstrap"),
        "{said}"
    );
}

#[test]
fn the_consensus_of_split_genes_has_the_branch_most_replicates_have_with_its_share() {
    let text = format!("{}{}", AB.repeat(40), AC.repeat(35));
    let mut keep = options(100, Resample::Genes, Vec::new());
    keep.keep_replicates = true;
    let out = infer(&text, &keep).expect("a species tree");
    assert_eq!(out.replicate_trees.len(), 100);
    // with four species the three splits cannot go together, so the consensus has the one that most replicates have: A with B
    // (the clade C, D) or A with C (the clade B, D)
    let having = |a: &str, b: &str| {
        out.replicate_trees
            .iter()
            .filter(|t| t.contains(&format!("({a},{b})")) || t.contains(&format!("({b},{a})")))
            .count()
    };
    let (cd, bd) = (having("C", "D"), having("B", "D"));
    assert!(
        cd + bd >= 99,
        "{cd} + {bd} of 100 (a tie of the three taxa is the only other tree)"
    );
    assert!(
        cd > 0 && bd > 0,
        "{cd} {bd}: the genes are split, so both appear"
    );
    let consensus = out.consensus.expect("a consensus");
    let (clade, best) = if cd >= bd { ("C,D", cd) } else { ("B,D", bd) };
    let label = format!("{:.1}", best as f64)
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string();
    assert!(
        consensus.contains(&format!("({clade}){label}")),
        "{consensus} should have ({clade}) at {label}"
    );
    // and the main tree's support is the share of the replicates that have its branch
    let main = out.branches[0].bootstrap.expect("a support");
    assert!(
        (main - cd as f64).abs() < 1e-9 || (main - bd as f64).abs() < 1e-9,
        "{main}"
    );
}

#[test]
fn the_bootstrap_itself_refuses_bootstrap_trees_of_the_sites_when_the_genes_are_weighted() {
    use super::bootstrap::{run, SiteTrees};
    use super::genes::read_weighted;
    use super::search;
    let genes = read_weighted(
        &weighted_gene(true, "0.9").repeat(4),
        0,
        super::Weights::Support,
    )
    .expect("weighted genes");
    let so = search::Options {
        weights: super::Weights::Support,
        exact: false,
        extra: 1,
        polylimit: 8,
        polish: true,
        root: 0,
    };
    let found = search::search(&genes, &so).expect("a species tree");
    let sites = SiteTrees::read(&genes, &bootstrap_trees_for(4)).expect("bootstrap trees");
    for resample in [Resample::Sites, Resample::Both] {
        let boot = Boot {
            replicates: 3,
            seed: 1,
            resample,
        };
        let said = run(&genes, Some(&sites), &so, &found.sp, &boot, false)
            .err()
            .expect("refused");
        assert!(said.contains("cannot be used together"), "{said}");
    }
    // the genes themselves are resampled with their weights
    let boot = Boot {
        replicates: 3,
        seed: 1,
        resample: Resample::Genes,
    };
    assert!(run(&genes, None, &so, &found.sp, &boot, false).is_ok());
}
