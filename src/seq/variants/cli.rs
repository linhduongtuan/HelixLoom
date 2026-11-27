//! `helixloom variants`: map reads to a reference (or read their SAM), pile them up, and write the differences as VCF, the consensus and the
//! depth.

use super::call::{call_target, Calls, Params};
use super::mapped::map_to_sam;
use super::output::{consensus, depth_table, summary, vcf, Target};
use super::pileup::{pile, Filter};
use super::reads::{self, Read};
use super::sam::{self, Sam};
use crate::seq::minimap::Options;

pub const USAGE: &str = "\
Usage: helixloom variants [options] REFERENCE READS

  The differences between reads and a reference, and the consensus, for a sample with one copy of the sequence (a plasmid, a bacterial
  genome). REFERENCE is FASTA. READS is FASTA or FASTQ (not compressed; mapped here with the mapper of `helixloom map`) or a SAM file (not
  BAM). The answer is a VCF with one haploid sample (docs/VARIANTS.md): the model is this program's own, and bcftools was not run.

  -o FILE            where to write the VCF (the standard output when left out)
  --consensus FILE   the consensus as FASTA: the reference with the calls made, N where too few reads
  --depth FILE       the reads over each position, as a table
  --sam FILE         also write the mapping as SAM (when READS is FASTA or FASTQ)
  --sample NAME      the name of the sample in the VCF (default: sample)
  -q N               lowest mapping quality of a read that counts (default 0)
  -Q N               lowest quality of a base that counts (default 13)
  -d N               most reads that count at one position (default 250)
  -D N               fewest reads with a base to call a site, or to write a base of the consensus (default 5)
  --min-qual N       lowest quality of a call that is made (default 20)
  --prior F          the chance that a base differs from the reference (default 0.001)
  --indel-prior F    the chance that a place has an insertion or a deletion (default 0.0001)
  --indel-error F    the chance that a read shows an indel that is not there, or misses one that is (default 0.1)
  --default-quality N  the quality of every base of a read that has none (FASTA, or * in a SAM; default 15)
  -k N, -w N         the minimizers of the mapping (default 15 and 10, as `helixloom map`)
  -h, --help         this text
";

struct Args {
    reference: String,
    reads: String,
    output: Option<String>,
    consensus: Option<String>,
    depth: Option<String>,
    sam: Option<String>,
    sample: String,
    filter: Filter,
    params: Params,
    mapping: Options,
}

fn number<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} needs a number, not {value}."))
}

fn parse(args: &[String]) -> Result<Args, String> {
    let mut a = Args {
        reference: String::new(),
        reads: String::new(),
        output: None,
        consensus: None,
        depth: None,
        sam: None,
        sample: "sample".into(),
        filter: Filter::default(),
        params: Params::default(),
        mapping: Options {
            cigar: true,
            ..Options::default()
        },
    };
    let mut files = Vec::new();
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value."))
        };
        match flag.as_str() {
            "-o" | "--output" => a.output = Some(value("-o")?),
            "--consensus" => a.consensus = Some(value("--consensus")?),
            "--depth" => a.depth = Some(value("--depth")?),
            "--sam" => a.sam = Some(value("--sam")?),
            "--sample" => a.sample = value("--sample")?,
            "-q" => a.filter.min_mapq = number("-q", &value("-q")?)?,
            "-Q" => a.filter.min_base_quality = number("-Q", &value("-Q")?)?,
            "-d" => a.filter.max_depth = number("-d", &value("-d")?)?,
            "-D" => a.params.min_depth = number("-D", &value("-D")?)?,
            "--min-qual" => a.params.min_quality = number("--min-qual", &value("--min-qual")?)?,
            "--prior" => a.params.snp_prior = number("--prior", &value("--prior")?)?,
            "--indel-prior" => {
                a.params.indel_prior = number("--indel-prior", &value("--indel-prior")?)?
            }
            "--indel-error" => {
                a.params.indel_error = number("--indel-error", &value("--indel-error")?)?
            }
            "--default-quality" => {
                a.filter.default_quality =
                    number("--default-quality", &value("--default-quality")?)?
            }
            "-k" => a.mapping.k = number("-k", &value("-k")?)?,
            "-w" => a.mapping.w = number("-w", &value("-w")?)?,
            other if other.starts_with('-') && other.len() > 1 => {
                return Err(format!(
                    "{other} is not an option (`helixloom variants --help`)."
                ))
            }
            path => files.push(path.to_string()),
        }
    }
    let [reference, reads] = <[String; 2]>::try_from(files).map_err(|_| {
        "Name a reference and the reads: helixloom variants REFERENCE READS.".to_string()
    })?;
    (a.reference, a.reads) = (reference, reads);
    if !(1..=32).contains(&a.mapping.k) {
        return Err(format!("k must be from 1 to 32, not {}.", a.mapping.k));
    }
    if a.mapping.w == 0 {
        return Err("w must be 1 or more.".into());
    }
    if a.filter.max_depth == 0 || a.filter.max_depth > usize::from(u16::MAX) {
        return Err(format!(
            "-d must be from 1 to {}, not {}.",
            u16::MAX,
            a.filter.max_depth
        ));
    }
    if a.filter.default_quality > 93 {
        return Err("The default quality must be from 0 to 93.".into());
    }
    a.params.check()?;
    Ok(a)
}

/// A SAM file begins with a header line `@HD`, `@SQ`, `@RG`, `@PG` or `@CO` and a tab; a FASTQ read's name could begin with `@HD` but not
/// continue with a tab.
fn is_sam(path: &str, text: &str) -> bool {
    path.to_ascii_lowercase().ends_with(".sam")
        || ["@HD\t", "@SQ\t", "@RG\t", "@PG\t", "@CO\t"]
            .iter()
            .any(|h| text.starts_with(h))
}

/// The reference, as targets with a name and the bases.
fn reference(path: &str) -> Result<Vec<(String, Vec<u8>)>, String> {
    let targets: Vec<(String, Vec<u8>)> = reads::load(path)?
        .into_iter()
        .map(|r| (r.name, r.seq))
        .collect();
    for (i, (name, seq)) in targets.iter().enumerate() {
        if seq.is_empty() {
            return Err(format!("The reference sequence {name} has no bases."));
        }
        if targets[..i].iter().any(|(other, _)| other == name) {
            return Err(format!(
                "Two sequences of the reference are named {name}; their names must differ."
            ));
        }
    }
    Ok(targets)
}

/// The alignments of the reads: mapped here, or read from a SAM file whose targets agree with the reference.
fn alignments(a: &Args, targets: &[(String, Vec<u8>)]) -> Result<(Sam, bool), String> {
    let text = std::fs::read_to_string(&a.reads)
        .map_err(|e| format!("Could not read {}: {e}", a.reads))?;
    if is_sam(&a.reads, &text) {
        if a.sam.is_some() {
            return Err("The reads are a SAM file already: --sam writes the mapping that this program makes.".into());
        }
        let sam = sam::parse(&text).map_err(|e| format!("{}: {e}", a.reads))?;
        for (name, length) in &sam.targets {
            match targets.iter().find(|t| &t.0 == name) {
                None => {
                    return Err(format!(
                        "The SAM file has the target {name}, which the reference does not."
                    ))
                }
                Some(t) if t.1.len() != *length => {
                    return Err(format!(
                        "The SAM file gives {name} a length of {length}, and the reference has {}.",
                        t.1.len()
                    ))
                }
                Some(_) => {}
            }
        }
        if let Some(r) = sam.records.iter().find(|r| {
            r.target
                .as_ref()
                .is_some_and(|t| !targets.iter().any(|(n, _)| n == t))
        }) {
            return Err(format!(
                "The read {} lies on {}, which the reference does not have.",
                r.name,
                r.target.as_deref().unwrap_or("*")
            ));
        }
        return Ok((sam, true));
    }
    let loaded: Vec<Read> = reads::from_text(&text).map_err(|e| format!("{}: {e}", a.reads))?;
    Ok((map_to_sam(targets, &loaded, &a.mapping)?, false))
}

fn write(path: &Option<String>, text: &str) -> Result<(), String> {
    match path {
        Some(p) => std::fs::write(p, text).map_err(|e| format!("Could not write {p}: {e}")),
        None => {
            print!("{text}");
            Ok(())
        }
    }
}

fn run_inner(args: &[String]) -> Result<(), String> {
    let a = parse(args)?;
    let targets = reference(&a.reference)?;
    let (sam, was_sam) = alignments(&a, &targets)?;
    if let (Some(path), false) = (&a.sam, was_sam) {
        write(&Some(path.clone()), &sam::write(&sam))?;
    }
    let lengths: Vec<(String, usize)> = targets.iter().map(|(n, s)| (n.clone(), s.len())).collect();
    let piles = pile(&lengths, &sam.records, &a.filter);
    let calls: Vec<Calls> = targets
        .iter()
        .zip(&piles)
        .map(|((_, seq), p)| call_target(seq, p, &a.params))
        .collect();
    let results: Vec<Target<'_>> = lengths
        .iter()
        .zip(&calls)
        .map(|((name, length), calls)| Target {
            name,
            length: *length,
            calls,
        })
        .collect();
    write(&a.output, &vcf(&results, &a.sample, &a.params))?;
    if let Some(path) = &a.consensus {
        write(&Some(path.clone()), &consensus(&results))?;
    }
    if let Some(path) = &a.depth {
        write(&Some(path.clone()), &depth_table(&results))?;
    }
    let used: usize = piles.iter().map(|p| p.reads).sum();
    let seen = sam.records.len();
    eprint!(
        "{used} alignments used of {seen} lines.\n{}",
        summary(&results)
    );
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
