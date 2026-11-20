//! `helixloom gff`: what `gffread` does for a GFF3 or GTF file, from the command line: read an annotation, filter it, write it as
//! GFF3 or GTF, and cut the spliced transcripts, coding sequences and proteins out of a genome.

use super::{
    cds_sequence, parse, protein, spliced, write_gff3, write_gtf, Annotation, Filters, Genome,
    Range, Strand, Transcript,
};
use std::collections::HashSet;

pub const USAGE: &str = "\
Usage: helixloom gff ANNOTATION [-g GENOME.fa] [options]

  Reads a GFF3 or GTF file (the format is told from the file), keeps the transcripts that pass the filters, and writes them as
  GFF3 (or GTF), and the sequences they cut out of a genome. After gffread's options of the same names; see docs/GFF.md.

  -g FILE        the genome, FASTA (needed for -w -x -y -N -V -J)
  -o FILE        write the annotation here (the standard output when none of -w -x -y --table is asked for)
  -T             write GTF in place of GFF3
  -F             keep all the attributes of the transcripts
  -w FILE        the spliced exons of each transcript, FASTA
  -x FILE        the coding sequence of each transcript that has one, FASTA
  -y FILE        the protein of each, FASTA (a stop inside is written as . ; -S writes it as *)
  -W             put the coordinates in the names of -w and -x (loc: and exons: or segs:)
  -S             write stops as * in -y
      --table T  one line for each transcript with the fields T (comma-separated): @id @geneid @genename @chr @start @end @strand
                 @numexons @exons @cds @covlen @cdslen, or the name of an attribute
  Filters:
  -C             only transcripts that have CDS            --nc         only those that have none
  -U             no transcript of a single exon            -i N         none with an intron longer than N
  -l N           none with an intron shorter than N        -N           none with an intron that is not GT-AG, GC-AG or AT-AC
  -V             none with a stop inside its CDS           -J           only complete CDS: ATG ... stop, none inside
  -r RANGE       only transcripts that overlap [strand]seq:start-end   -R   only those wholly inside it (with -r)
      --ids FILE   only the ids in the file (one to a line)             --nids FILE   not those
  -h, --help     this text

Not done: merging transcripts into loci (-M -K -Q -Y), BED and TLF, -u -j -O -t -L -Z and the other options of gffread.
";

#[derive(Default)]
struct Args {
    annotation: Option<String>,
    genome: Option<String>,
    output: Option<String>,
    gtf: bool,
    keep_all: bool,
    spliced: Option<String>,
    cds: Option<String>,
    protein: Option<String>,
    coordinates: bool,
    star: bool,
    table: Option<String>,
    filters: Filters,
    range: Option<Range>,
    inside: bool,
    ids: Option<String>,
    not_ids: Option<String>,
}

fn parse_args(args: &[String]) -> Result<Args, String> {
    let mut a = Args::default();
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        let mut value = |name: &str| {
            it.next()
                .cloned()
                .ok_or_else(|| format!("{name} needs a value."))
        };
        match flag.as_str() {
            "-g" => a.genome = Some(value("-g")?),
            "-o" => a.output = Some(value("-o")?),
            "-w" => a.spliced = Some(value("-w")?),
            "-x" => a.cds = Some(value("-x")?),
            "-y" => a.protein = Some(value("-y")?),
            "--table" => a.table = Some(value("--table")?),
            "--ids" => a.ids = Some(value("--ids")?),
            "--nids" => a.not_ids = Some(value("--nids")?),
            "-r" => a.range = Some(Range::parse(&value("-r")?)?),
            "-i" => a.filters.max_intron = Some(number("-i", &value("-i")?)?),
            "-l" => a.filters.min_intron = Some(number("-l", &value("-l")?)?),
            "-T" => a.gtf = true,
            "-F" => a.keep_all = true,
            "-W" => a.coordinates = true,
            "-S" => a.star = true,
            "-C" => a.filters.coding_only = true,
            "--nc" => a.filters.noncoding_only = true,
            "-U" => a.filters.no_single_exon = true,
            "-N" => a.filters.canonical_splice = true,
            "-V" => a.filters.no_inframe_stop = true,
            "-J" => a.filters.complete_cds = true,
            "-R" => a.inside = true,
            f if f.starts_with('-') && f.len() > 1 => {
                return Err(format!(
                    "The option {f} is not known here; `helixloom gff --help` lists what is done."
                ))
            }
            bare if a.annotation.is_none() => a.annotation = Some(bare.to_string()),
            bare => return Err(format!("{bare} is one argument too many.")),
        }
    }
    if a.inside && a.range.is_none() {
        return Err("-R goes with -r.".into());
    }
    Ok(a)
}

fn number(flag: &str, value: &str) -> Result<usize, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} needs a number, not {value}."))
}

fn read(path: &str) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("Could not read {path}: {e}"))
}

fn ids(path: &str) -> Result<HashSet<String>, String> {
    Ok(read(path)?
        .lines()
        .map(|l| l.split_whitespace().next().unwrap_or("").to_string())
        .filter(|l| !l.is_empty())
        .collect())
}

fn wrapped(name: &str, bases: &[u8]) -> String {
    let mut out = format!(">{name}\n");
    for chunk in bases.chunks(70) {
        out += &String::from_utf8_lossy(chunk);
        out.push('\n');
    }
    out
}

/// The coordinates for `-W`: the exons as genomic ranges.
fn exons_text(t: &Transcript) -> String {
    t.exons
        .iter()
        .map(|(a, b)| format!("{a}-{b}"))
        .collect::<Vec<_>>()
        .join(",")
}

/// The CDS segments as ranges within the coding sequence, counted from its 5' end.
fn segments_text(t: &Transcript, length: usize) -> String {
    let mut lengths: Vec<usize> = t.cds.iter().map(|c| c.end - c.start + 1).collect();
    if t.strand == Strand::Minus {
        lengths.reverse();
    }
    let skip = lengths.iter().sum::<usize>().saturating_sub(length);
    let mut at = 0;
    let mut parts = Vec::new();
    for (k, n) in lengths.into_iter().enumerate() {
        let n = if k == 0 { n.saturating_sub(skip) } else { n };
        parts.push(format!("{}-{}", at + 1, at + n));
        at += n;
    }
    parts.join(",")
}

fn field(t: &Transcript, name: &str) -> String {
    let cds_len: usize = t.cds.iter().map(|c| c.end - c.start + 1).sum();
    let exon_len: usize = t.exons.iter().map(|(a, b)| b - a + 1).sum();
    match name {
        "@id" => t.id.clone(),
        "@geneid" => t.gene_id.clone().unwrap_or_else(|| ".".into()),
        "@genename" => t.gene_name.clone().unwrap_or_else(|| ".".into()),
        "@chr" => t.seqname.clone(),
        "@start" => t.start.to_string(),
        "@end" => t.end.to_string(),
        "@strand" => t.strand.symbol().to_string(),
        "@numexons" => t.exons.len().to_string(),
        "@exons" => exons_text(t),
        "@cds" => t
            .cds
            .iter()
            .map(|c| format!("{}-{}", c.start, c.end))
            .collect::<Vec<_>>()
            .join(","),
        "@covlen" => exon_len.to_string(),
        "@cdslen" => cds_len.to_string(),
        attribute => t
            .attributes
            .iter()
            .find(|(k, _)| k == attribute)
            .map_or_else(|| ".".to_string(), |(_, v)| v.clone()),
    }
}

fn run_inner(args: &[String]) -> Result<(), String> {
    let a = parse_args(args)?;
    let path = a
        .annotation
        .as_deref()
        .ok_or("Name the annotation file (GFF3 or GTF): `helixloom gff --help`.")?;
    let genome = a
        .genome
        .as_deref()
        .map(|p| read(p).map(|t| Genome::from_fasta(&t)))
        .transpose()?;
    let needs_genome = a.spliced.is_some() || a.cds.is_some() || a.protein.is_some();
    if needs_genome && genome.is_none() {
        return Err("-w, -x and -y cut sequences out of a genome: give it with -g.".into());
    }
    let mut filters = a.filters.clone();
    filters.range = a.range.clone().map(|r| (r, a.inside));
    filters.only_ids = a.ids.as_deref().map(ids).transpose()?;
    filters.not_ids = a.not_ids.as_deref().map(ids).transpose()?;
    let mut annotation = parse(&read(path)?).map_err(|e| format!("{path}: {e}"))?;
    let before = annotation.transcripts.len();
    let mut kept = Vec::new();
    for t in annotation.transcripts {
        if filters.keeps(&t, genome.as_ref())? {
            kept.push(t);
        }
    }
    annotation = Annotation { transcripts: kept };
    annotation.sort_by_location();
    eprintln!(
        "{} of {before} transcripts kept.",
        annotation.transcripts.len()
    );
    let write = |file: &str, text: String| {
        std::fs::write(file, text).map_err(|e| format!("Could not write {file}: {e}"))
    };
    if let (Some(g), true) = (&genome, needs_genome) {
        let (mut w, mut x, mut y) = (String::new(), String::new(), String::new());
        let stop = if a.star { b'*' } else { b'.' };
        for t in &annotation.transcripts {
            let name = |extra: String| {
                if a.coordinates && !extra.is_empty() {
                    format!("{} {extra}", t.id)
                } else {
                    t.id.clone()
                }
            };
            let loc = format!(
                "loc:{}({}){}-{}",
                t.seqname,
                t.strand.symbol(),
                t.start,
                t.end
            );
            if a.spliced.is_some() {
                w += &wrapped(
                    &name(format!("{loc} exons:{}", exons_text(t))),
                    &spliced(t, g)?,
                );
            }
            if a.cds.is_some() {
                if let Some(cds) = cds_sequence(t, g)? {
                    x += &wrapped(
                        &name(format!("{loc} segs:{}", segments_text(t, cds.len()))),
                        &cds,
                    );
                }
            }
            if a.protein.is_some() {
                if let Some(p) = protein(t, g, stop)? {
                    y += &wrapped(&t.id, &p);
                }
            }
        }
        for (file, text) in [(&a.spliced, w), (&a.cds, x), (&a.protein, y)] {
            if let Some(file) = file {
                write(file, text)?;
            }
        }
    }
    if let Some(table) = &a.table {
        let fields: Vec<&str> = table.split(',').collect();
        let text: String = annotation
            .transcripts
            .iter()
            .map(|t| {
                fields
                    .iter()
                    .map(|f| field(t, f))
                    .collect::<Vec<_>>()
                    .join("\t")
                    + "\n"
            })
            .collect();
        match &a.output {
            Some(file) => write(file, text)?,
            None => print!("{text}"),
        }
        return Ok(());
    }
    let text = if a.gtf {
        write_gtf(&annotation, a.keep_all)
    } else {
        write_gff3(&annotation, a.keep_all)
    };
    match (&a.output, needs_genome) {
        (Some(file), _) => write(file, text)?,
        (None, false) => print!("{text}"),
        // the sequences were what was asked for: no annotation in the way, unless -o says where it goes
        (None, true) => {}
    }
    Ok(())
}

pub fn run(args: &[String]) -> i32 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
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
