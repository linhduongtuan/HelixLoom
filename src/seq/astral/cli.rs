//! `helixloom astral`: the species-tree search from the command line, with the options of ASTRAL's that this engine has,
//! so that a pipeline (a Snakemake or Nextflow file, a job on a cluster) can call it as it would call ASTRAL.

use super::{bootstrap, infer, score_tree, Options, Resample, Weights};
use std::time::Instant;

pub const USAGE: &str = "\
Usage: helixloom astral -i GENE_TREES [-o SPECIES_TREE] [options]    (also: helixloom species-tree)

  Estimates the unrooted species tree that shares the most quartets with a set of unrooted gene trees (the method of
  ASTRAL), with a length in coalescent units and a local posterior probability on every inner branch. The gene trees are
  Newick trees, one after the other, with one leaf per species in each (they may lack species and have polytomies).

  -i, --input FILE        the gene trees (a last argument that is no option is the input too)
  -o, --output FILE       where to write the species tree (the standard output when left out)
  -q, --score-tree FILE   do not search: score and annotate the binary species tree in FILE
  -t, --branch-annotate N what each branch carries: 0 nothing, 1 quartet support in percent, 2 everything,
                          3 local posterior probability (the default), 4 the three posteriors, 8 the three quartet supports,
                          10 the p-value of the test for a polytomy (also --polytomy; see below)
      --polytomy          test every branch for being a polytomy (Sayyari and Mirarab 2018): a chi-squared test of whether
                          the gene trees show its three resolutions equally often, for a branch that at least 10 gene trees
                          resolve a quartet around. The p-values are the labels of the branches (a branch not tested has
                          only its length; a bootstrap support, if asked for, is the label instead). Not with --families
  -x, --exact             look at every group of taxa (at most 12 taxa)
  -f, --families          the gene trees are gene family trees (ASTRAL-Pro, Zhang et al. 2020): a species may occur in a tree
                          more than once, the leaf labels being the species names; a tree with a root of two branches is
                          rooted as it is, one with three is rooted where duplications and losses cost least. Binary trees
                          only; no branch support or length is given on the species tree; not with --weights
  -p, --extraLevel N      0: only the gene trees' own splits; 1 (the default): also the ways to resolve polytomies
  -c, --lambda X          the Yule prior of the branch lengths (default 0.5; 0 gives maximum likelihood lengths)
  -m, --minleaves N       leave out gene trees with fewer leaves
  -r, --replicates N      bootstrap replicates (100 when -b or --gene-only is given without it): the search is run again on
                          resampled genes, and the number on each branch is the share of those trees (in percent) that have it
  -b, --bootstrap LIST    a file that lists, one path to a line, the files of bootstrap trees of the genes, in the order of the
                          gene trees (a relative path is taken from the folder of the list): the sites are resampled, replicate
                          r taking the r-th bootstrap tree of every gene, so each file needs at least N trees
  -g, --genes-and-sites   with -b: draw the genes with replacement too, and for each a bootstrap tree of its own
      --gene-only         resample the genes only, from the gene trees alone (no -b)
  -s, --seed N            the seed of the resampling (default 692); the same seed gives the same support
                          (not with --families or -q; with --weights only when the genes are resampled, not the sites)
      --replicate-trees FILE  also write the species tree of every replicate, one to a line, to FILE
      --consensus FILE    also write the greedy consensus of the replicate trees to FILE: the splits most replicates have,
                          each kept if it fits the ones before, labelled with the share of the replicates that have it
  -w, --weights NAME      weight the quartets of the gene trees (weighted ASTRAL, Zhang and Mirarab 2022): none (the default),
                          support (by the support of the branches between the two pairs of a quartet), length (by the
                          lengths of the paths of the pairs) or hybrid (both). Supports are the labels of the inner
                          branches (0 to 1, or percentages), lengths the branch lengths; the branch support on the output
                          is still ASTRAL's, counted without the weights
      --root TAXON        hang the tree from this taxon when writing it (also --outgroup)
      --no-polish         do not swap neighbouring groups after the search
      --polylimit N       the largest polytomy whose pairs of branches are joined (default 8)
  -h, --help              this text

The log goes to the standard error. See docs/SPECIES.md.
";

struct Args {
    /// Where to write the replicate trees and their consensus.
    replicate_out: Option<String>,
    consensus_out: Option<String>,
    /// The list of the files of bootstrap trees.
    list: Option<String>,
    gene_only: bool,
    with_genes: bool,
    /// The replicates asked for with -r, when it was.
    replicates: Option<usize>,
    input: Option<String>,
    output: Option<String>,
    score: Option<String>,
    options: Options,
}

fn number<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} needs a number, not {value}."))
}

/// Set what an option that takes a value asks for.
fn set_value(a: &mut Args, flag: &str, value: &str) -> Result<(), String> {
    let o = &mut a.options;
    match flag {
        "-i" | "--input" => a.input = Some(value.into()),
        "-o" | "--output" => a.output = Some(value.into()),
        "-q" | "--score-tree" => a.score = Some(value.into()),
        "-t" | "--branch-annotate" => {
            o.annotate = number(flag, value)?;
            if o.annotate == 10 {
                // ASTRAL's -t 10: the test, whose p-values the branches then carry
                o.annotate = Options::default().annotate;
                o.polytomy = true;
            } else if ![0, 1, 2, 3, 4, 8].contains(&o.annotate) {
                return Err("-t takes 0, 1, 2, 3, 4, 8 or 10.".into());
            }
        }
        "-p" | "--extraLevel" => o.extra = number::<u8>(flag, value)?.min(1),
        "-c" | "--lambda" => {
            o.lambda = number(flag, value)?;
            if !(o.lambda >= 0.0 && o.lambda.is_finite()) {
                return Err("-c takes a number that is not negative.".into());
            }
        }
        "-m" | "--minleaves" => o.min_leaves = number(flag, value)?,
        "-r" | "--replicates" => a.replicates = Some(number(flag, value)?),
        "-s" | "--seed" => o.bootstrap.seed = number(flag, value)?,
        "-b" | "--bootstrap" => a.list = Some(value.into()),
        "--replicate-trees" => a.replicate_out = Some(value.into()),
        "--consensus" => a.consensus_out = Some(value.into()),
        "-w" | "--weights" => {
            o.weights = Weights::parse(value)
                .ok_or("-w takes none, support, length or hybrid.".to_string())?;
        }
        "--polylimit" => o.polylimit = number(flag, value)?,
        "--root" | "--outgroup" => o.root = Some(value.into()),
        _ => {
            return Err(format!(
                "The option {flag} is not known; `helixloom astral --help` lists them."
            ))
        }
    }
    Ok(())
}

/// What the bootstrap options together ask for: how to resample and how many replicates, or no bootstrap.
fn settle_bootstrap(a: &mut Args) -> Result<(), String> {
    let resample = match (a.list.is_some(), a.gene_only, a.with_genes) {
        (false, false, false) => {
            if a.replicates.is_some() {
                return Err(
                    "-r needs -b LIST (the bootstrap trees of the genes) or --gene-only.".into(),
                );
            }
            if a.replicate_out.is_some() || a.consensus_out.is_some() {
                return Err("--replicate-trees and --consensus are for a bootstrap: ask for one with --gene-only or -b LIST.".into());
            }
            return Ok(());
        }
        (true, true, _) => {
            return Err("--gene-only uses the gene trees alone: leave out -b.".into());
        }
        (false, true, true) => {
            return Err("-g is for -b: with --gene-only the genes are resampled anyway.".into());
        }
        (false, false, true) => return Err("-g needs -b LIST.".into()),
        (false, true, false) => Resample::Genes,
        (true, false, false) => Resample::Sites,
        (true, false, true) => Resample::Both,
    };
    a.options.bootstrap.resample = resample;
    a.options.keep_replicates = a.replicate_out.is_some() || a.consensus_out.is_some();
    a.options.bootstrap.replicates = a.replicates.unwrap_or(100);
    if a.options.bootstrap.replicates == 0 {
        return Err("-r needs at least 1 replicate.".into());
    }
    Ok(())
}

fn parse(args: &[String]) -> Result<Args, String> {
    let mut a = Args {
        replicate_out: None,
        consensus_out: None,
        list: None,
        gene_only: false,
        with_genes: false,
        replicates: None,
        input: None,
        output: None,
        score: None,
        options: Options::default(),
    };
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        match flag.as_str() {
            "-x" | "--exact" => a.options.exact = true,
            "-f" | "--families" => a.options.families = true,
            "-g" | "--genes-and-sites" => a.with_genes = true,
            "--gene-only" => a.gene_only = true,
            "--no-polish" => a.options.polish = false,
            "--polytomy" => a.options.polytomy = true,
            f if f.starts_with('-') => {
                let value = it.next().ok_or_else(|| format!("{f} needs a value."))?;
                set_value(&mut a, f, value)?;
            }
            // a last argument that is no option is the input, as in ASTER
            bare if a.input.is_none() => a.input = Some(bare.into()),
            bare => return Err(format!("{bare} is one argument too many.")),
        }
    }
    settle_bootstrap(&mut a)?;
    Ok(a)
}

fn read(path: &str) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("Could not read {path}: {e}"))
}

fn work(args: &[String]) -> Result<(), String> {
    let mut a = parse(args)?;
    if let Some(list) = &a.list {
        a.options.site_trees = bootstrap::read_list(std::path::Path::new(list))?;
    }
    let input = a.input.ok_or("No gene trees given: use -i FILE.")?;
    let genes = read(&input)?;
    let started = Instant::now();
    let found = match &a.score {
        Some(path) => score_tree(&genes, &read(path)?, &a.options)?,
        None => infer(&genes, &a.options)?,
    };
    for line in &found.log {
        eprintln!("{line}");
    }
    if let Some(n) = found.clusters {
        eprintln!(
            "Searched {n} groups of taxa in {:.3} s",
            started.elapsed().as_secs_f64()
        );
    }
    if let Some(path) = &a.replicate_out {
        std::fs::write(
            path,
            found
                .replicate_trees
                .iter()
                .map(|t| format!("{t}\n"))
                .collect::<String>(),
        )
        .map_err(|e| format!("Could not write {path}: {e}"))?;
    }
    if let (Some(path), Some(tree)) = (&a.consensus_out, &found.consensus) {
        std::fs::write(path, format!("{tree}\n"))
            .map_err(|e| format!("Could not write {path}: {e}"))?;
    }
    match a.output {
        Some(path) => std::fs::write(&path, format!("{}\n", found.newick))
            .map_err(|e| format!("Could not write {path}: {e}")),
        None => {
            println!("{}", found.newick);
            Ok(())
        }
    }
}

/// Run the command; the exit code.
pub fn run(args: &[String]) -> i32 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{USAGE}");
        return 0;
    }
    match work(args) {
        Ok(()) => 0,
        Err(why) => {
            eprintln!("{why}");
            1
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(s: &str) -> Vec<String> {
        s.split_whitespace().map(str::to_string).collect()
    }

    #[test]
    fn the_options_are_read_like_astrals() {
        let a = parse(&words(
            "-i g.tre -o s.tre -t 2 -x -p 0 -c 0 -m 5 --root t3 --no-polish",
        ))
        .expect("options");
        assert_eq!(
            (a.input.as_deref(), a.output.as_deref()),
            (Some("g.tre"), Some("s.tre"))
        );
        assert_eq!(
            (a.options.annotate, a.options.exact, a.options.extra),
            (2, true, 0)
        );
        assert_eq!((a.options.lambda, a.options.min_leaves), (0.0, 5));
        assert_eq!(a.options.root.as_deref(), Some("t3"));
        assert!(!a.options.polish);
        // a last argument that is no option is the input, as in ASTER
        assert_eq!(
            parse(&words("-t 1 genes.tre"))
                .expect("options")
                .input
                .as_deref(),
            Some("genes.tre")
        );
        let q = parse(&words("-i g -q sp.tre")).expect("options");
        assert_eq!(q.score.as_deref(), Some("sp.tre"));
        assert_eq!(a.options.weights, Weights::None);
        assert!(!a.options.families);
        assert!(parse(&words("-i g -f")).expect("options").options.families);
        // no bootstrap unless asked for
        assert_eq!(a.options.bootstrap.replicates, 0);
        assert!(
            parse(&words("--families -i g"))
                .expect("options")
                .options
                .families
        );
        for (word, weights) in [
            ("hybrid", Weights::Hybrid),
            ("support", Weights::Support),
            ("length", Weights::Length),
            ("none", Weights::None),
        ] {
            let w = parse(&words(&format!("-i g --weights {word}"))).expect("options");
            assert_eq!(w.options.weights, weights);
        }
        assert_eq!(
            parse(&words("-i g -w Hybrid"))
                .expect("options")
                .options
                .weights,
            Weights::Hybrid
        );
    }

    #[test]
    fn the_bootstrap_options_say_how_to_resample_and_how_many_times() {
        let o = |text: &str| parse(&words(text)).expect("options").options.bootstrap;
        let b = o("-i g --gene-only");
        assert_eq!(
            (b.resample, b.replicates, b.seed),
            (Resample::Genes, 100, 692)
        );
        let b = o("-i g --gene-only -r 30 -s 5");
        assert_eq!((b.replicates, b.seed), (30, 5));
        let b = o("-i g -b list.txt");
        assert_eq!((b.resample, b.replicates), (Resample::Sites, 100));
        let b = o("-i g -b list.txt -g -r 7");
        assert_eq!((b.resample, b.replicates), (Resample::Both, 7));
        // what makes no sense together is said
        let refused = |text: &str, why: &str| {
            let said = parse(&words(text)).err().expect("refused");
            assert!(said.contains(why), "{text}: {said}");
        };
        refused("-i g -r 10", "needs -b LIST");
        refused("-i g -b l --gene-only", "leave out -b");
        refused("-i g -g", "needs -b");
        refused("-i g --gene-only -g", "-g is for -b");
        refused("-i g --gene-only -r 0", "at least 1");
        refused("-i g -r many --gene-only", "needs a number");
    }

    #[test]
    fn the_polytomy_test_is_asked_for_with_t_10_as_in_astral_or_with_its_own_name() {
        let o = |text: &str| parse(&words(text)).expect("options").options;
        assert!(!o("-i g").polytomy);
        for text in ["-i g -t 10", "-i g --polytomy", "-i g --branch-annotate 10"] {
            let o = o(text);
            assert!(o.polytomy, "{text}");
            assert_eq!(
                o.annotate, 3,
                "the labels are the p-values, not another level"
            );
        }
        // labels switched off stay off whatever else is asked
        let off = o("-i g -t 0 --polytomy");
        assert!(off.polytomy && off.annotate == 0);
        assert!(parse(&words("-i g -t 5"))
            .err()
            .is_some_and(|e| e.contains("10")));
    }

    #[test]
    fn a_wrong_option_is_said_and_refused() {
        assert!(parse(&words("-i")).is_err());
        assert!(parse(&words("-i g -t 5")).is_err());
        assert!(parse(&words("-i g -c -1")).is_err());
        assert!(parse(&words("-i g --nonsense")).is_err());
        assert!(parse(&words("-i g -m many")).is_err());
        assert!(parse(&words("-i g -w heavy"))
            .err()
            .is_some_and(|e| e.contains("hybrid")));
        assert_eq!(run(&words("-o only.tre")), 1);
        assert_eq!(run(&words("--help")), 0);
    }
}
