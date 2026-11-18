//! `helixloom astral` as a pipeline runs it: the real program, started with files and options, with none of the program's own
//! code in between. The gene trees are the example of the species-tree window (40 genes on 8 species simulated from the tree
//! in `assets/species/truth.tre`); the score 1733 is what this engine and ASTRAL 5.7.8 both give for them.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const GENES: &str = include_str!("../assets/species/genes.tre");
const TRUTH: &str = include_str!("../assets/species/truth.tre");
const SCORE: &str = "Final quartet score is: 1733";

/// A folder for one test, removed when the test is over.
struct Folder(PathBuf);

impl Folder {
    fn new(label: &str) -> Folder {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let path = std::env::temp_dir().join(format!(
            "helixloom-astral-{label}-{}-{nanos}",
            std::process::id()
        ));
        std::fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }

    fn file(&self, name: &str, text: &str) -> PathBuf {
        let path = self.0.join(name);
        std::fs::write(&path, text).expect("a file");
        path
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn astral(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_helixloom"))
        .arg("astral")
        .args(args)
        .output()
        .expect("the program starts")
}

fn path(p: &Path) -> &str {
    p.to_str().expect("a path that is text")
}

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn stdout(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).into_owned()
}

#[test]
fn it_writes_the_species_tree_to_a_file_and_the_log_to_the_standard_error() {
    let dir = Folder::new("search");
    let genes = dir.file("genes.tre", GENES);
    let found = dir.0.join("species.tre");
    let out = astral(&["-i", path(&genes), "-o", path(&found)]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stdout(&out).is_empty(), "the tree goes to the file");
    let log = stderr(&out);
    assert!(
        log.contains("40 gene trees on 8 taxa") && log.contains(SCORE),
        "{log}"
    );
    let tree = std::fs::read_to_string(&found).expect("the tree");
    assert!(tree.ends_with(";\n"), "{tree}");
    for species in "ABCDEFGH".chars() {
        assert!(tree.contains(species), "{species} in {tree}");
    }
    // the local posterior probability and a length on every inner branch, as ASTRAL writes them
    assert!(tree.contains(")0.6:") || tree.contains(")1:"), "{tree}");
}

#[test]
fn the_same_input_gives_the_same_bytes_every_time_whatever_the_option_spelling() {
    let dir = Folder::new("same");
    let genes = dir.file("genes.tre", GENES);
    let first = astral(&["-i", path(&genes)]);
    assert!(first.status.success(), "{}", stderr(&first));
    assert!(stdout(&first).trim_end().ends_with(';'));
    for _ in 0..4 {
        // a last argument that is no option is the input, as in ASTER
        let again = astral(&[path(&genes)]);
        assert_eq!(again.stdout, first.stdout);
    }
    let long = astral(&["--input", path(&genes), "--branch-annotate", "3"]);
    assert_eq!(long.stdout, first.stdout);
}

#[test]
fn a_given_tree_is_scored_and_annotated_with_what_each_option_asks_for() {
    let dir = Folder::new("score");
    let genes = dir.file("genes.tre", GENES);
    let truth = dir.file("truth.tre", TRUTH);
    let out = astral(&["-i", path(&genes), "-q", path(&truth), "-t", "2"]);
    assert!(out.status.success(), "{}", stderr(&out));
    let log = stderr(&out);
    assert!(
        log.contains(SCORE),
        "the tree the genes came from has the best score: {log}"
    );
    assert!(!log.contains("Searched"), "nothing was searched: {log}");
    let tree = stdout(&out);
    for part in ["pp1=", "q1=", "f1=", "QC=", "EN="] {
        assert!(tree.contains(part), "{part} in {tree}");
    }
    let quartets = astral(&["-i", path(&genes), "-q", path(&truth), "-t", "8"]);
    assert!(stdout(&quartets).contains("'[q1="), "{}", stdout(&quartets));
    let none = astral(&["-i", path(&genes), "-q", path(&truth), "-t", "0"]);
    assert!(!stdout(&none).contains("'["), "{}", stdout(&none));
}

#[test]
fn the_exact_search_and_the_search_of_clusters_find_the_same_best_score_on_eight_species() {
    let dir = Folder::new("exact");
    let genes = dir.file("genes.tre", GENES);
    let exact = astral(&["-i", path(&genes), "-x"]);
    assert!(exact.status.success(), "{}", stderr(&exact));
    assert!(stderr(&exact).contains(SCORE), "{}", stderr(&exact));
    let plain = astral(&["-i", path(&genes), "-p", "0", "--no-polish"]);
    assert!(stderr(&plain).contains(SCORE), "{}", stderr(&plain));
}

#[test]
fn what_cannot_be_done_is_said_on_the_standard_error_with_a_failing_exit_code() {
    let dir = Folder::new("refuse");
    let genes = dir.file("genes.tre", GENES);
    let family = dir.file("family.tre", "((A,B),(A,C),D);\n");
    let few = dir.file("few.tre", "((A,B),C);\n((A,C),B);\n");
    let cases: [(Vec<&str>, &str); 6] = [
        (vec![], "No gene trees given"),
        (
            vec!["-i", "/nonexistent/genes.tre"],
            "/nonexistent/genes.tre",
        ),
        (vec!["-i", path(&family)], "ASTRAL-Pro"),
        (vec!["-i", path(&few)], "at least four taxa"),
        (vec!["-i", path(&genes), "--nonsense", "1"], "--nonsense"),
        (vec!["-i", path(&genes), "-t", "5"], "-t"),
    ];
    for (args, said) in cases {
        let out = astral(&args);
        assert_eq!(out.status.code(), Some(1), "{args:?}");
        assert!(stderr(&out).contains(said), "{args:?}: {}", stderr(&out));
        assert!(stdout(&out).is_empty(), "{args:?}");
    }
}

#[test]
fn help_lists_the_options_and_succeeds_and_the_other_name_of_the_command_works() {
    let help = astral(&["--help"]);
    assert!(help.status.success());
    let text = stdout(&help);
    for part in [
        "Usage: helixloom astral",
        "--score-tree",
        "--branch-annotate",
        "--exact",
        "--lambda",
    ] {
        assert!(text.contains(part), "{part}");
    }
    let dir = Folder::new("alias");
    let genes = dir.file("genes.tre", GENES);
    let out = Command::new(env!("CARGO_BIN_EXE_helixloom"))
        .args(["species-tree", "-i", path(&genes)])
        .output()
        .expect("the program starts");
    assert!(out.status.success() && String::from_utf8_lossy(&out.stderr).contains(SCORE));
}

#[test]
fn names_with_spaces_and_brackets_are_written_in_quotes_and_the_tree_reads_back() {
    let dir = Folder::new("quotes");
    let genes = dir.file(
        "genes.tre",
        "(('Homo sapiens','Pan troglodytes'),('Nomascus (gibbon)','Pongo pygmaeus'));\n\
         (('Homo sapiens','Nomascus (gibbon)'),('Pan troglodytes','Pongo pygmaeus'));\n\
         (('Homo sapiens','Pan troglodytes'),('Nomascus (gibbon)','Pongo pygmaeus'));\n",
    );
    let out = astral(&["-i", path(&genes)]);
    assert!(out.status.success(), "{}", stderr(&out));
    let tree = stdout(&out);
    for name in [
        "'Homo sapiens'",
        "'Nomascus (gibbon)'",
        "'Pan troglodytes'",
        "'Pongo pygmaeus'",
    ] {
        assert!(tree.contains(name), "{name} in {tree}");
    }
    // the text that was written is a tree that can be given back: scored against the same genes it is the same tree
    let species = dir.file("species.tre", &tree);
    let again = astral(&["-i", path(&genes), "-q", path(&species)]);
    assert!(again.status.success(), "{}", stderr(&again));
    let score = |out: &Output| {
        stderr(out)
            .lines()
            .find(|l| l.starts_with("Final quartet score"))
            .map(str::to_string)
    };
    assert!(score(&out).is_some() && score(&out) == score(&again));
}

#[test]
fn a_given_tree_is_hung_from_the_species_asked_for_and_one_that_is_not_there_is_refused() {
    let dir = Folder::new("root");
    let genes = dir.file("genes.tre", GENES);
    let truth = dir.file("truth.tre", TRUTH);
    let out = astral(&["-i", path(&genes), "-q", path(&truth), "--root", "H"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stdout(&out).starts_with("(H,"), "{}", stdout(&out));
    let nobody = astral(&["-i", path(&genes), "-q", path(&truth), "--root", "NOPE"]);
    assert_eq!(nobody.status.code(), Some(1));
    assert!(
        stderr(&nobody).contains("root taxon NOPE"),
        "{}",
        stderr(&nobody)
    );
}

/// Gene trees that disagree: three, each sure that A goes with B, and four that join A with C and are hardly sure of it.
fn disagreeing_genes() -> String {
    let ab = "(A:0.1,B:0.1,(C:0.1,D:0.1)0.95:0.1);\n".repeat(3);
    let ac = "(A:0.1,C:0.1,(B:0.1,D:0.1)0.05:0.1);\n".repeat(4);
    format!("{ab}{ac}")
}

#[test]
fn weights_make_the_sure_gene_trees_count_more_and_the_log_says_so() {
    let folder = Folder::new("weights");
    let genes = folder.file("genes.tre", &disagreeing_genes());
    let tree = |args: &[&str]| {
        let mut all = vec!["-i", path(&genes), "-t", "0"];
        all.extend_from_slice(args);
        let out = astral(&all);
        assert!(out.status.success(), "{}", stderr(&out));
        (stdout(&out).trim().to_string(), stderr(&out))
    };
    // counted one each, the four that join A with C win; weighted, the three that are sure of A with B
    let (plain, plain_log) = tree(&[]);
    let (weighted, log) = tree(&["--weights", "support"]);
    assert_eq!(plain, "(A,(C,(B,D)));", "{plain_log}");
    assert_eq!(weighted, "(A,(B,(C,D)));", "{log}");
    assert!(
        log.contains("Weights: quartets weighted by the support"),
        "{log}"
    );
    assert!(
        log.contains("Final weighted quartet score is: 2.850000"),
        "{log}"
    );
    assert!(log.contains("Final quartet score is: 3"), "{log}");
    assert!(!plain_log.contains("Weights:"), "{plain_log}");
    // the same through the short option, and a tree that is scored rather than searched
    assert_eq!(tree(&["-w", "hybrid"]).0, weighted);
    let given = folder.file("species.tre", "((A,B),(C,D));");
    let (_, scored) = tree(&["-w", "support", "-q", path(&given)]);
    assert!(
        scored.contains("Final weighted quartet score is: 2.850000"),
        "{scored}"
    );
}

#[test]
fn weights_that_the_gene_trees_cannot_give_stop_the_run_with_the_reason_and_a_failing_exit_code() {
    let folder = Folder::new("weights-missing");
    let bare = folder.file("bare.tre", "((A,B),(C,D));\n((A,B),(C,D));");
    for (weights, why) in [
        ("support", "no branch support"),
        ("hybrid", "no branch support"),
        ("length", "no branch lengths"),
    ] {
        let out = astral(&["-i", path(&bare), "--weights", weights]);
        assert!(!out.status.success(), "{weights}");
        assert!(stderr(&out).contains(why), "{weights}: {}", stderr(&out));
    }
    let out = astral(&["-i", path(&bare), "--weights", "heavy"]);
    assert!(!out.status.success());
    assert!(
        stderr(&out).contains("none, support, length or hybrid"),
        "{}",
        stderr(&out)
    );
}

const FAMILIES: &str = "(((A,B),(C,(D,E))),((A,B),(C,(D,E))));\n\
                        (((A,B),(C,(D,E))),((A,B),(D,E)));\n\
                        ((A,(B,(C,(D,E)))),((A,B),(D,E)));\n";

#[test]
fn gene_families_are_searched_with_the_families_option_and_refused_without_it() {
    let folder = Folder::new("families");
    let genes = folder.file("families.tre", FAMILIES);
    // as single-copy trees a repeated species is refused, with what to do
    let out = astral(&["-i", path(&genes)]);
    assert!(!out.status.success());
    assert!(stderr(&out).contains("more than once"), "{}", stderr(&out));
    // as gene families they give the species tree, and the log says what was read
    let out = astral(&["-i", path(&genes), "--families"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert_eq!(
        stdout(&out).trim(),
        "(A,(B,(C,(D,E))));",
        "{}",
        stderr(&out)
    );
    let log = stderr(&out);
    assert!(log.contains("3 gene family trees on 5 species"), "{log}");
    assert!(log.contains("ASTRAL-Pro"), "{log}");
    assert!(!log.contains("Weights:"), "{log}");
    // the short option, a given tree scored, and the combination with weights refused
    assert!(astral(&["-i", path(&genes), "-f"]).status.success());
    let given = folder.file("species.tre", "((A,B),(C,(D,E)));");
    let out = astral(&["-i", path(&genes), "-f", "-q", path(&given)]);
    assert!(out.status.success(), "{}", stderr(&out));
    let out = astral(&["-i", path(&genes), "-f", "-w", "support"]);
    assert!(!out.status.success());
    assert!(
        stderr(&out).contains("cannot be used together"),
        "{}",
        stderr(&out)
    );
}

#[test]
fn bootstrap_support_is_written_on_the_branches_from_genes_alone_or_from_files_of_bootstrap_trees()
{
    let folder = Folder::new("bootstrap");
    let genes = folder.file("genes.tre", &"((A,B),(C,D));\n".repeat(5));
    // genes alone: all agree, so every replicate has the branch
    let out = astral(&["-i", path(&genes), "--gene-only", "-r", "20"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stdout(&out).contains(")100:"), "{}", stdout(&out));
    let log = stderr(&out);
    assert!(
        log.contains(
            "Bootstrap: 20 replicates resampling the genes, seed 692; 20 gave a species tree"
        ),
        "{log}"
    );
    // the seed is the option that changes what is drawn, and the same one gives the same text
    let again = astral(&["-i", path(&genes), "--gene-only", "-r", "20", "--seed", "3"]);
    assert!(stderr(&again).contains("seed 3"), "{}", stderr(&again));
    // the sites: every gene has twenty bootstrap trees, ten saying A with B and ten saying A with C, and replicate r takes the
    // r-th of each: half the replicates have the branch
    let both_ways = format!(
        "{}{}",
        "((A,B),(C,D));\n".repeat(10),
        "((A,C),(B,D));\n".repeat(10)
    );
    let mut list = String::from("# five genes\n");
    for g in 1..=5 {
        folder.file(&format!("boot{g}.tre"), &both_ways);
        list.push_str(&format!("boot{g}.tre\n"));
    }
    let list = folder.file("list.txt", &list);
    let out = astral(&["-i", path(&genes), "-b", path(&list), "-r", "20"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stdout(&out).contains(")50:"), "{}", stdout(&out));
    assert!(
        stderr(&out).contains("resampling the sites"),
        "{}",
        stderr(&out)
    );
    let both = astral(&["-i", path(&genes), "-b", path(&list), "-g", "-r", "20"]);
    assert!(both.status.success(), "{}", stderr(&both));
    assert!(
        stderr(&both).contains("genes and their sites"),
        "{}",
        stderr(&both)
    );
    // what is asked that cannot be done stops the run with the reason
    for (args, why) in [
        (vec!["-i", path(&genes), "-r", "5"], "needs -b LIST"),
        (
            vec!["-i", path(&genes), "-b", "no-such-list.txt"],
            "Could not read the list",
        ),
        (
            vec!["-i", path(&genes), "-b", path(&list), "-w", "support"],
            "cannot be used together",
        ),
        (
            vec!["-i", path(&genes), "--replicate-trees", "r.tre"],
            "are for a bootstrap",
        ),
        (vec!["-i", path(&genes), "--gene-only", "-f"], "not built"),
    ] {
        let out = astral(&args);
        assert!(!out.status.success(), "{args:?}");
        assert!(stderr(&out).contains(why), "{args:?}: {}", stderr(&out));
    }
}

/// The labels of the inner branches of a Newick tree: the text between a closing bracket and the colon that follows it (quoted
/// labels in single quotes are taken whole).
fn labels(newick: &str) -> Vec<String> {
    let text: Vec<char> = newick.trim().chars().collect();
    let mut found = Vec::new();
    let mut i = 0;
    while i < text.len() {
        if text[i] == ')' {
            let mut label = String::new();
            i += 1;
            let mut quoted = false;
            while i < text.len() && (quoted || !matches!(text[i], ':' | ',' | ')' | ';')) {
                quoted ^= text[i] == '\'';
                label.push(text[i]);
                i += 1;
            }
            found.push(label);
        } else {
            i += 1;
        }
    }
    found
}

/// The three frequencies `f1 f2 f3` of every branch in a tree annotated with `-t 2`.
fn frequencies(newick: &str) -> Vec<[f64; 3]> {
    labels(newick)
        .iter()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let value = |key: &str| -> f64 {
                let start = l.find(&format!("{key}=")).expect("the key") + key.len() + 1;
                let end = l[start..].find([';', ']']).map_or(l.len(), |e| start + e);
                l[start..end].parse().expect("a number")
            };
            [value("f1"), value("f2"), value("f3")]
        })
        .collect()
}

#[test]
fn the_polytomy_test_labels_each_branch_with_the_p_value_the_chi_squared_of_its_frequencies_gives()
{
    let folder = Folder::new("polytomy");
    let genes = folder.file("genes.tre", GENES);
    let truth = folder.file("truth.tre", TRUTH);
    // the example's 40 gene trees are complete and binary, so the counts of the test are the frequencies of `-t 2`
    let by_t10 = astral(&["-i", path(&genes), "-q", path(&truth), "-t", "10"]);
    assert!(by_t10.status.success(), "{}", stderr(&by_t10));
    let by_name = astral(&["-i", path(&genes), "-q", path(&truth), "--polytomy"]);
    assert_eq!(stdout(&by_name), stdout(&by_t10), "-t 10 is --polytomy");
    let log = stderr(&by_t10);
    assert!(
        log.contains("Polytomy test (Sayyari and Mirarab 2018): 5 of 5 branches tested"),
        "{log}"
    );
    let all = astral(&["-i", path(&genes), "-q", path(&truth), "-t", "2"]);
    let counts = frequencies(&stdout(&all));
    // (the root's brackets carry no branch, so their labels are empty)
    let ps: Vec<String> = labels(&stdout(&by_t10))
        .into_iter()
        .filter(|l| !l.is_empty())
        .collect();
    assert_eq!((counts.len(), ps.len()), (5, 5), "{}", stdout(&by_t10));
    for (f, p) in counts.iter().zip(&ps) {
        let n: f64 = f.iter().sum();
        let chi2: f64 = f.iter().map(|x| (x - n / 3.0).powi(2)).sum::<f64>() / (n / 3.0);
        let expected = (-chi2 / 2.0).exp();
        let written: f64 = p.parse().unwrap_or_else(|_| panic!("a number, not {p:?}"));
        assert!(
            (written - expected).abs() <= 1e-3 + 0.06 * expected,
            "{f:?}: {written} against {expected}"
        );
    }
    // the same program also finds the tree itself, and labels it the same way
    let searched = astral(&["-i", path(&genes), "--polytomy"]);
    assert!(searched.status.success(), "{}", stderr(&searched));
    assert_eq!(
        labels(&stdout(&searched))
            .iter()
            .filter(|l| !l.is_empty())
            .count(),
        5
    );
    // labels off stay off; the log still says what the test found
    let off = astral(&[
        "-i",
        path(&genes),
        "-q",
        path(&truth),
        "-t",
        "0",
        "--polytomy",
    ]);
    assert!(!stdout(&off).contains(')') || labels(&stdout(&off)).iter().all(String::is_empty));
    assert!(
        stderr(&off).contains("No branch is labelled"),
        "{}",
        stderr(&off)
    );
}

#[test]
fn the_polytomy_test_leaves_out_branches_with_too_few_gene_trees_and_gives_way_to_a_bootstrap_label(
) {
    let folder = Folder::new("polytomy-few");
    // four gene trees: fewer than the ten the test is made for
    let few = folder.file("few.tre", &"((A,B),(C,D));\n".repeat(4));
    let out = astral(&["-i", path(&few), "--polytomy"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(
        stderr(&out).contains("0 of 1 branches tested"),
        "{}",
        stderr(&out)
    );
    assert!(
        labels(&stdout(&out)).iter().all(String::is_empty),
        "an untested branch has only its length: {}",
        stdout(&out)
    );
    // with a bootstrap, the support is the label and the log says where the p-values are
    let genes = folder.file("genes.tre", GENES);
    let out = astral(&["-i", path(&genes), "--gene-only", "-r", "5", "--polytomy"]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stdout(&out).contains(")100:"), "{}", stdout(&out));
    assert!(
        stderr(&out).contains("the p-values are in the table of branches"),
        "{}",
        stderr(&out)
    );
    // gene families have no table of quartet counts to test
    let families = folder.file("families.tre", "((A,A),(B,C));\n");
    let out = astral(&["-i", path(&families), "--families", "--polytomy"]);
    assert!(!out.status.success());
    assert!(
        stderr(&out).contains("polytomy test") && stderr(&out).contains("cannot be used together"),
        "{}",
        stderr(&out)
    );
}

#[test]
fn the_replicate_trees_and_their_consensus_are_written_when_asked_for_and_weights_go_with_the_genes(
) {
    let folder = Folder::new("consensus");
    let text = format!(
        "{}{}",
        "((A,B),(C,D));\n".repeat(30),
        "((A,C),(B,D));\n".repeat(20)
    );
    let genes = folder.file("genes.tre", &text);
    let (replicates, consensus) = (
        folder.0.join("replicates.tre"),
        folder.0.join("consensus.tre"),
    );
    let out = astral(&[
        "-i",
        path(&genes),
        "--gene-only",
        "-r",
        "25",
        "--replicate-trees",
        path(&replicates),
        "--consensus",
        path(&consensus),
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    let trees = std::fs::read_to_string(&replicates).expect("the replicate trees");
    assert_eq!(trees.lines().count(), 25, "{trees}");
    assert!(
        trees
            .lines()
            .all(|l| l.starts_with("(A,") && l.ends_with(");")),
        "{trees}"
    );
    let joined = std::fs::read_to_string(&consensus).expect("the consensus");
    assert_eq!(joined.lines().count(), 1);
    // four species: one branch, with the share of the replicates that have the most frequent grouping of the two
    assert!(
        joined.starts_with("(A,(")
            && joined.trim_end().ends_with(");")
            && joined.matches(')').count() == 3,
        "{joined}"
    );
    // the main tree is still what the standard output has, with the support as its label
    assert!(stdout(&out).contains(")"), "{}", stdout(&out));
    // weights go with resampling the genes: genes with supports and lengths
    let weighted = folder.file(
        "weighted.tre",
        &format!(
            "{}{}",
            "((A:1,B:1)0.9:1,(C:1,D:1)0.9:1);\n".repeat(20),
            "((A:1,C:1)0:1,(B:1,D:1)0:1);\n".repeat(15)
        ),
    );
    let out = astral(&[
        "-i",
        path(&weighted),
        "--gene-only",
        "-r",
        "20",
        "-w",
        "support",
    ]);
    assert!(out.status.success(), "{}", stderr(&out));
    assert!(stdout(&out).contains(")100:"), "{}", stdout(&out));
}
