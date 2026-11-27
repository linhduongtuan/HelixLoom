//! `helixloom count`: how many reads fall on each gene, from SAM files and an annotation.

use super::annotation;
use super::assign::{count, Index, Rules, Splice, Tally};
use super::output::{summary, table, Column};
use crate::seq::variants::sam;

pub const USAGE: &str = "\
Usage: helixloom count [options] -a ANNOTATION ALIGNMENTS.sam...

  The number of reads on each gene (or exon, with -f), after featureCounts (Liao, Smyth and Shi, Bioinformatics 2014), written from its
  users' guide and paper (docs/COUNT.md). ANNOTATION is GTF, GFF3 or a SAF table (GeneID Chr Start End Strand); the alignments are SAM (not
  BAM), for example from `helixloom map -a`. Pairs are not counted as fragments yet.

  -a FILE                 the annotation
  -o FILE                 the count table; the summary goes to FILE.summary (without -o, the table to the standard output and the summary
                          to the standard error)
  -t TYPE[,TYPE]          the lines of the annotation to count (default exon)
  -g ATTRIBUTE            the attribute that names the gene a line belongs to (default gene_id)
  -f                      count each feature (exon) and not each gene
  -O                      count a read that overlaps more than one gene for each of them (otherwise it is ambiguous)
  -M                      count a multi-mapping read (the NH tag) at each place (otherwise not at all)
  --fraction              with -M and/or -O: share the count between the places and the genes
  --primary               count primary alignments only (then -M is not needed and secondaries are left out)
  -Q N                    the lowest mapping quality that counts (default 0)
  -s 0|1|2                strand: 0 either (default), 1 the read's strand must be the feature's, 2 it must be the other
  --minOverlap N          the fewest overlapping bases (default 1)
  --fracOverlap F         the fewest overlapping bases as a share of the read (soft clips counted in its length)
  --fracOverlapFeature F  the fewest overlapping bases as a share of the feature or gene
  --nonOverlap N          the most bases of the read that may lie outside
  --nonOverlapFeature N   the most bases of the feature or gene that may be left uncovered
  --largestOverlap        give a read that overlaps several targets to the one it overlaps most
  --ignoreDup             leave out reads marked as duplicates
  --splitOnly             count only alignments that skip over a junction (N in the CIGAR)
  --nonSplitOnly          count only alignments that do not
  --each-read             count each read of a pair on its own
  -h, --help              this text
";

struct Args {
    annotation: String,
    output: Option<String>,
    types: Vec<String>,
    attribute: String,
    rules: Rules,
    files: Vec<String>,
}

fn number<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} needs a whole number, not {value}."))
}

fn share(flag: &str, value: &str) -> Result<f64, String> {
    value
        .parse::<f64>()
        .ok()
        .filter(|v| (0.0..=1.0).contains(v))
        .ok_or_else(|| format!("{flag} needs a share between 0 and 1, not {value}."))
}

fn parse(args: &[String]) -> Result<Args, String> {
    let mut a = Args {
        annotation: String::new(),
        output: None,
        types: vec!["exon".into()],
        attribute: "gene_id".into(),
        rules: Rules::default(),
        files: Vec::new(),
    };
    let (mut split, mut whole) = (false, false);
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value."))
        };
        match flag.as_str() {
            "-a" => a.annotation = value("-a")?,
            "-o" => a.output = Some(value("-o")?),
            "-t" => {
                a.types = value("-t")?
                    .split(',')
                    .filter(|t| !t.is_empty())
                    .map(String::from)
                    .collect()
            }
            "-g" => a.attribute = value("-g")?,
            "-f" => a.rules.feature_level = true,
            "-O" => a.rules.multi_overlap = true,
            "-M" => a.rules.multi_mapping = true,
            "--fraction" => a.rules.fraction = true,
            "--primary" => a.rules.primary_only = true,
            "-Q" => a.rules.min_mapq = number("-Q", &value("-Q")?)?,
            "-s" => a.rules.strandedness = number("-s", &value("-s")?)?,
            "--minOverlap" => {
                a.rules.min_overlap = number("--minOverlap", &value("--minOverlap")?)?
            }
            "--fracOverlap" => {
                a.rules.frac_overlap = share("--fracOverlap", &value("--fracOverlap")?)?
            }
            "--fracOverlapFeature" => {
                a.rules.frac_overlap_feature =
                    share("--fracOverlapFeature", &value("--fracOverlapFeature")?)?
            }
            "--nonOverlap" => {
                a.rules.non_overlap = Some(number("--nonOverlap", &value("--nonOverlap")?)?)
            }
            "--nonOverlapFeature" => {
                a.rules.non_overlap_feature = Some(number(
                    "--nonOverlapFeature",
                    &value("--nonOverlapFeature")?,
                )?)
            }
            "--largestOverlap" => a.rules.largest_overlap = true,
            "--ignoreDup" => a.rules.ignore_duplicates = true,
            "--splitOnly" => split = true,
            "--nonSplitOnly" => whole = true,
            "--each-read" => a.rules.each_read = true,
            other if other.starts_with('-') && other.len() > 1 => {
                return Err(format!(
                    "{other} is not an option (`helixloom count --help`)."
                ))
            }
            path => a.files.push(path.to_string()),
        }
    }
    a.rules.splice = match (split, whole) {
        (true, true) => return Err("--splitOnly and --nonSplitOnly leave nothing to count.".into()),
        (true, false) => Splice::Split,
        (false, true) => Splice::Whole,
        (false, false) => Splice::Any,
    };
    if a.annotation.is_empty() {
        return Err("Name the annotation with -a (GTF, GFF3 or SAF).".into());
    }
    if a.files.is_empty() {
        return Err("Name at least one SAM file to count.".into());
    }
    if a.types.is_empty() {
        return Err("-t needs at least one type.".into());
    }
    if a.rules.strandedness > 2 {
        return Err("-s must be 0, 1 or 2.".into());
    }
    if a.rules.fraction && !a.rules.multi_mapping && !a.rules.multi_overlap {
        return Err("--fraction needs -M, -O or both: otherwise there is nothing to share.".into());
    }
    Ok(a)
}

/// The counts of each file, and whether any read lay on a sequence that the annotation names.
fn counts(a: &Args, index: &Index<'_>) -> Result<(Vec<Tally>, bool), String> {
    let mut tallies = Vec::new();
    let mut seen = false;
    for path in &a.files {
        let text =
            std::fs::read_to_string(path).map_err(|e| format!("Could not read {path}: {e}"))?;
        let sam = sam::parse(&text).map_err(|e| format!("{path}: {e}"))?;
        seen |= sam
            .records
            .iter()
            .any(|r| r.target.as_ref().is_some_and(|t| index.names(t)));
        tallies.push(count(index, &sam.records, &a.rules).map_err(|e| format!("{path}: {e}"))?);
    }
    Ok((tallies, seen))
}

fn write(path: &str, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("Could not write {path}: {e}"))
}

fn run_inner(args: &[String]) -> Result<(), String> {
    let a = parse(args)?;
    let text = std::fs::read_to_string(&a.annotation)
        .map_err(|e| format!("Could not read {}: {e}", a.annotation))?;
    let annotation = annotation::parse(&text, &a.types, &a.attribute)
        .map_err(|e| format!("{}: {e}", a.annotation))?;
    let index = Index::new(&annotation, a.rules.feature_level);
    let (tallies, seen) = counts(&a, &index)?;
    if !seen {
        eprintln!("No read lies on a sequence that the annotation names: check that the two name the sequences alike (chr1 and 1, for example).");
    }
    let columns: Vec<Column<'_>> = a
        .files
        .iter()
        .zip(&tallies)
        .map(|(name, tally)| Column { name, tally })
        .collect();
    let counted = table(&index, &columns, a.rules.feature_level, a.rules.fraction);
    match &a.output {
        Some(path) => {
            write(path, &counted)?;
            write(&format!("{path}.summary"), &summary(&columns))?;
        }
        None => {
            print!("{counted}");
            eprint!("{}", summary(&columns));
        }
    }
    Ok(())
}

pub fn run(args: &[String]) -> i32 {
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{USAGE}");
        return 0;
    }
    match run_inner(args) {
        Ok(()) => 0,
        Err(why) => {
            eprintln!("{why}");
            1
        }
    }
}

#[cfg(test)]
mod tests;
