//! `helixloom astral`: the species-tree search from the command line, with the options of ASTRAL's that this engine has,
//! so that a pipeline (a Snakemake or Nextflow file, a job on a cluster) can call it as it would call ASTRAL.

use super::{infer, score_tree, Options};
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
                          3 local posterior probability (the default), 4 the three posteriors, 8 the three quartet supports
  -x, --exact             look at every group of taxa (at most 12 taxa)
  -p, --extraLevel N      0: only the gene trees' own splits; 1 (the default): also the ways to resolve polytomies
  -c, --lambda X          the Yule prior of the branch lengths (default 0.5; 0 gives maximum likelihood lengths)
  -m, --minleaves N       leave out gene trees with fewer leaves
      --root TAXON        hang the tree from this taxon when writing it (also --outgroup)
      --no-polish         do not swap neighbouring groups after the search
      --polylimit N       the largest polytomy whose pairs of branches are joined (default 8)
  -h, --help              this text

The log goes to the standard error. See docs/SPECIES.md.
";

struct Args {
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
            if ![0, 1, 2, 3, 4, 8].contains(&o.annotate) {
                return Err("-t takes 0, 1, 2, 3, 4 or 8.".into());
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

fn parse(args: &[String]) -> Result<Args, String> {
    let mut a = Args {
        input: None,
        output: None,
        score: None,
        options: Options::default(),
    };
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        match flag.as_str() {
            "-x" | "--exact" => a.options.exact = true,
            "--no-polish" => a.options.polish = false,
            f if f.starts_with('-') => {
                let value = it.next().ok_or_else(|| format!("{f} needs a value."))?;
                set_value(&mut a, f, value)?;
            }
            // a last argument that is no option is the input, as in ASTER
            bare if a.input.is_none() => a.input = Some(bare.into()),
            bare => return Err(format!("{bare} is one argument too many.")),
        }
    }
    Ok(a)
}

fn read(path: &str) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("Could not read {path}: {e}"))
}

fn work(args: &[String]) -> Result<(), String> {
    let a = parse(args)?;
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
    }

    #[test]
    fn a_wrong_option_is_said_and_refused() {
        assert!(parse(&words("-i")).is_err());
        assert!(parse(&words("-i g -t 5")).is_err());
        assert!(parse(&words("-i g -c -1")).is_err());
        assert!(parse(&words("-i g --nonsense")).is_err());
        assert!(parse(&words("-i g -m many")).is_err());
        assert_eq!(run(&words("-o only.tre")), 1);
        assert_eq!(run(&words("--help")), 0);
    }
}
