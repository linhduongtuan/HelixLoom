//! `helixloom kat`: the k-mer analyses of `seq::kmer` from the command line (`hist`, `gcp`, `comp`, `sect`).

use super::{check_k, comparison, gc_matrix, section, spectrum, summarize, Counts};
use crate::pipeline::steps::read_fasta;
use crate::seq::fastq::parse_fastq;

pub const DEFAULT_K: usize = 27;

pub const USAGE: &str = "\
Usage: helixloom kat COMMAND [options] INPUT...

  K-mer analysis in the manner of KAT (Mapleson et al., Bioinformatics 2017), written from its paper and manual (docs/KAT.md). An
  INPUT is a FASTA or FASTQ file (not compressed); several files that are to be one input are joined with commas: a.fq,b.fq.
  k-mers are canonical (a k-mer and its reverse complement are one) and hold only A, C, G, T (U counts as T).

  hist INPUT...                 the spectrum: how many distinct k-mers occur once, twice ...; with the valley, the peak and the genome
                                size that they say
  gcp INPUT...                  the same by the number of G and C in the k-mer: a matrix
  comp INPUT1 INPUT2            how often each k-mer occurs in one against the other: a matrix, with row and column 0 for what only
                                one has
  sect -c READS SEQUENCES       how often each k-mer of each sequence occurs in READS: statistics for each sequence, and with
                                --coverage the numbers along it

  -k N          the length of the k-mer, 1 to 32 (default 27)
  -o FILE       where to write the answer (the standard output when left out)
      --low N   hist: leave out k-mers that occur fewer than N times (default 1)
      --high N  hist, gcp: the last bin holds all that occur N times or more (default 10000)
      --high1 N, --high2 N   comp: the same for the first and the second input (default 1000)
  -c INPUT      sect: the reads whose k-mers are counted
      --coverage FILE   sect: the number for each k-mer along each sequence, one sequence to two lines (- where a k-mer has a letter
                        that is no base)
  -h, --help    this text
";

struct Args {
    command: String,
    k: usize,
    output: Option<String>,
    low: u32,
    high: u32,
    high1: u32,
    high2: u32,
    counts: Option<String>,
    coverage: Option<String>,
    inputs: Vec<String>,
}

fn number<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} needs a number, not {value}."))
}

fn parse(args: &[String]) -> Result<Args, String> {
    let mut it = args.iter();
    let command = it
        .next()
        .ok_or("Name a command: hist, gcp, comp or sect (`helixloom kat --help`).")?
        .clone();
    let mut a = Args {
        command,
        k: DEFAULT_K,
        output: None,
        low: 1,
        high: 10_000,
        high1: 1000,
        high2: 1000,
        counts: None,
        coverage: None,
        inputs: Vec::new(),
    };
    while let Some(flag) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value."))
        };
        match flag.as_str() {
            "-k" => a.k = number("-k", &value("-k")?)?,
            "-o" | "--output" => a.output = Some(value("-o")?),
            "--low" => a.low = number("--low", &value("--low")?)?,
            "--high" => a.high = number("--high", &value("--high")?)?,
            "--high1" => a.high1 = number("--high1", &value("--high1")?)?,
            "--high2" => a.high2 = number("--high2", &value("--high2")?)?,
            "-c" => a.counts = Some(value("-c")?),
            "--coverage" => a.coverage = Some(value("--coverage")?),
            f if f.starts_with('-') && f.len() > 1 => {
                return Err(format!(
                    "The option {f} is not known; `helixloom kat --help` lists them."
                ))
            }
            bare => a.inputs.push(bare.to_string()),
        }
    }
    check_k(a.k)?;
    if a.high == 0 || a.high1 == 0 || a.high2 == 0 {
        return Err("--high, --high1 and --high2 need at least 1.".into());
    }
    Ok(a)
}

/// The sequences of an input: FASTA or FASTQ files, joined with commas.
fn read_input(spec: &str) -> Result<Vec<(String, Vec<u8>)>, String> {
    let mut out = Vec::new();
    for path in spec.split(',').filter(|p| !p.is_empty()) {
        let text =
            std::fs::read_to_string(path).map_err(|e| format!("Could not read {path}: {e}"))?;
        let first = text.trim_start().chars().next();
        match first {
            Some('>') => out.extend(read_fasta(&text)),
            Some('@') => out.extend(
                parse_fastq(&text)
                    .map_err(|e| format!("{path}: {e}"))?
                    .into_iter()
                    .map(|r| (r.name, r.seq)),
            ),
            _ => {
                return Err(format!(
                    "{path} is neither FASTA nor FASTQ: it should begin with > or @."
                ))
            }
        }
    }
    if out.is_empty() {
        return Err(format!("{spec} holds no sequences."));
    }
    Ok(out)
}

fn count_of(spec: &str, k: usize) -> Result<Counts, String> {
    let seqs = read_input(spec)?;
    let slices: Vec<&[u8]> = seqs.iter().map(|(_, s)| s.as_slice()).collect();
    Counts::count(&slices, k)
}

/// A row or column of a matrix as far as the last that is not 0.
fn last_used(counts: impl Iterator<Item = u64>) -> usize {
    counts
        .enumerate()
        .filter(|(_, n)| *n > 0)
        .map(|(i, _)| i)
        .last()
        .unwrap_or(0)
}

fn hist(a: &Args) -> Result<String, String> {
    if a.inputs.is_empty() {
        return Err("hist needs an input.".into());
    }
    let c = count_of(&a.inputs.join(","), a.k)?;
    let spec = spectrum(&c, a.low, a.high);
    let s = summarize(&spec);
    let mut out = format!(
        "# k: {}\n# distinct k-mers: {}\n# total k-mers: {}\n",
        c.k,
        c.distinct(),
        c.total()
    );
    match (s.valley, s.peak, s.genome_size) {
        (Some(v), Some(p), Some(g)) => {
            out += &format!("# valley: {v}\n# peak: {p}\n# estimated genome size: {g}\n# share of distinct k-mers below the valley: {:.3}\n", s.error_share);
            if let Some(h) = s.half_peak {
                out += &format!("# half peak: {h}\n");
            }
        }
        _ => out += "# no peak after the errors: too few reads for a size\n",
    }
    out += "###\n";
    let first = a.low.max(1) as usize;
    let last = last_used(spec.iter().copied())
        .max(first)
        .min(a.high as usize);
    for (m, n) in spec.iter().enumerate().take(last + 1).skip(first) {
        out += &format!("{m}\t{n}\n");
    }
    Ok(out)
}

fn gcp(a: &Args) -> Result<String, String> {
    if a.inputs.is_empty() {
        return Err("gcp needs an input.".into());
    }
    let c = count_of(&a.inputs.join(","), a.k)?;
    let rows = gc_matrix(&c, a.high);
    let width = rows
        .iter()
        .map(|r| last_used(r.iter().copied()))
        .max()
        .unwrap_or(1)
        .max(1);
    let mut out = format!("# k: {}\n# rows: the number of G and C in the k-mer; columns: how often it occurs\n###\ngc", c.k);
    for m in 1..=width {
        out += &format!("\t{m}");
    }
    out.push('\n');
    for (g, row) in rows.iter().enumerate() {
        out += &g.to_string();
        for n in &row[1..=width] {
            out += &format!("\t{n}");
        }
        out.push('\n');
    }
    Ok(out)
}

fn comp(a: &Args) -> Result<String, String> {
    let [first, second] = a.inputs.as_slice() else {
        return Err("comp needs two inputs.".into());
    };
    let (x, y) = (count_of(first, a.k)?, count_of(second, a.k)?);
    let rows = comparison(&x, &y, a.high1, a.high2)?;
    let shared: u64 = rows.iter().skip(1).flat_map(|r| r.iter().skip(1)).sum();
    let only_first: u64 = rows.iter().map(|r| r[0]).sum();
    let only_second: u64 = rows[0].iter().sum();
    let mut out = format!(
        "# k: {}\n# distinct k-mers in both: {shared}\n# only in the first: {only_first}\n# only in the second: {only_second}\n# rows: occurrences in the first; columns: in the second\n###\n",
        x.k
    );
    let width = rows
        .iter()
        .map(|r| last_used(r.iter().copied()))
        .max()
        .unwrap_or(0);
    let height = last_used(rows.iter().map(|r| r.iter().sum()));
    for row in &rows[..=height] {
        out += &row[..=width]
            .iter()
            .map(u64::to_string)
            .collect::<Vec<_>>()
            .join("\t");
        out.push('\n');
    }
    Ok(out)
}

fn sect(a: &Args) -> Result<(String, String), String> {
    let reads = a
        .counts
        .as_deref()
        .ok_or("sect needs the reads whose k-mers are counted: -c READS.")?;
    let [file] = a.inputs.as_slice() else {
        return Err("sect needs one file of sequences.".into());
    };
    let c = count_of(reads, a.k)?;
    let mut stats = String::from("seq_name\tmedian\tmean\tgc%\tlength\tkmers\tinvalid_kmers\t%invalid\tnon_zero_kmers\t%non_zero\n");
    let mut along = String::new();
    for (name, seq) in read_input(file)? {
        let s = section(&c, &seq);
        let kmers = s.coverage.len();
        let pct = |n: usize, of: usize| {
            if of == 0 {
                0.0
            } else {
                100.0 * n as f64 / of as f64
            }
        };
        stats += &format!(
            "{name}\t{:.1}\t{:.2}\t{:.2}\t{}\t{kmers}\t{}\t{:.2}\t{}\t{:.2}\n",
            s.median,
            s.mean,
            s.gc,
            seq.len(),
            s.invalid,
            pct(s.invalid, kmers),
            s.non_zero,
            pct(s.non_zero, kmers - s.invalid)
        );
        along += &format!(">{name}\n");
        along += &s
            .coverage
            .iter()
            .map(|v| v.map_or("-".to_string(), |n| n.to_string()))
            .collect::<Vec<_>>()
            .join(" ");
        along.push('\n');
    }
    Ok((stats, along))
}

fn run_inner(args: &[String]) -> Result<(), String> {
    let a = parse(args)?;
    let write = |file: &Option<String>, text: &str| match file {
        Some(f) => std::fs::write(f, text).map_err(|e| format!("Could not write {f}: {e}")),
        None => {
            print!("{text}");
            Ok(())
        }
    };
    match a.command.as_str() {
        "hist" => write(&a.output, &hist(&a)?),
        "gcp" => write(&a.output, &gcp(&a)?),
        "comp" => write(&a.output, &comp(&a)?),
        "sect" => {
            let (stats, along) = sect(&a)?;
            if let Some(file) = &a.coverage {
                std::fs::write(file, along).map_err(|e| format!("Could not write {file}: {e}"))?;
            }
            write(&a.output, &stats)
        }
        other => Err(format!(
            "{other} is not a command here: hist, gcp, comp or sect (`helixloom kat --help`)."
        )),
    }
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
