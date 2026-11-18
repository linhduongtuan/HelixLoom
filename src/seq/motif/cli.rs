//! `helixloom motif`: the motif search of the Motif Finder from the command line, so that a pipeline (a Snakemake or Nextflow file,
//! a job on a cluster) can call it: sequences in, motifs in MEME's text format and the places they were found out.

use super::{background, discover, meme_text, read_sequences, Found, Settings, EFFORT_RANGE};
use web_time::Instant;

pub const USAGE: &str = "\
Usage: helixloom motif -i SEQUENCES [-o MOTIFS.meme] [options]

  Looks for motifs that DNA sequences share, by expectation-maximisation (MEME's method in its simplest form: at most one site in
  a sequence, on either strand), and writes them in MEME's text format, which FIMO, Tomtom and other programs read.

  -i, --input FILE      the sequences: FASTA, or one sequence on each line (a last argument that is no option is the input too)
  -o, --output FILE     where to write the motifs (the standard output when left out)
      --sites FILE      also write the places found, one to a line with a header: motif, sequence, position (counted from 1),
                        strand, score in bits, probability that it is a site
  -w, --width N         columns of a motif, 3 to 40 (default 8)
  -n, --motifs N        how many motifs, one after the other: the places of one are hidden from the next search (default 1)
      --one-strand      look on the strand as written only (the default is both)
      --seeds N         starting words followed to the end (default 12)
      --iterations N    rounds of the search from each start (default 60)
      --effort X        how much work the search may spend over a big set, 1 being the budget of the Motif Finder window
                        (about five seconds); 0.01 to 1000 (default 1). A small set is searched in full whatever this is
  -h, --help            this text

The log goes to the standard error. See docs/MOTIFS.md.
";

struct Args {
    input: Option<String>,
    output: Option<String>,
    sites: Option<String>,
    settings: Settings,
}

fn number<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} needs a number, not {value}."))
}

fn parse(args: &[String]) -> Result<Args, String> {
    let mut a = Args {
        input: None,
        output: None,
        sites: None,
        settings: Settings::default(),
    };
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        match flag.as_str() {
            "--one-strand" => a.settings.both_strands = false,
            f if f.starts_with('-') => {
                let value = it.next().ok_or_else(|| format!("{f} needs a value."))?;
                match f {
                    "-i" | "--input" => a.input = Some(value.clone()),
                    "-o" | "--output" => a.output = Some(value.clone()),
                    "--sites" => a.sites = Some(value.clone()),
                    "-w" | "--width" => a.settings.width = number(f, value)?,
                    "-n" | "--motifs" => a.settings.motifs = number(f, value)?,
                    "--seeds" => a.settings.seeds = number(f, value)?,
                    "--iterations" => a.settings.iterations = number(f, value)?,
                    "--effort" => a.settings.effort = number(f, value)?,
                    _ => {
                        return Err(format!(
                            "The option {f} is not known; `helixloom motif --help` lists them."
                        ))
                    }
                }
            }
            bare if a.input.is_none() => a.input = Some(bare.into()),
            bare => return Err(format!("{bare} is one argument too many.")),
        }
    }
    let s = &a.settings;
    if s.motifs == 0 || s.seeds == 0 || s.iterations == 0 {
        return Err("-n, --seeds and --iterations need at least 1.".into());
    }
    if !(EFFORT_RANGE.0..=EFFORT_RANGE.1).contains(&s.effort) {
        return Err(format!(
            "--effort is from {} to {}.",
            EFFORT_RANGE.0, EFFORT_RANGE.1
        ));
    }
    Ok(a)
}

/// The places the motifs were found, a line each, with a header (the position counted from 1, as the sequence is written).
pub fn sites_table(found: &[Found], names: &[String]) -> String {
    let mut text = String::from("motif\tsequence\tposition\tstrand\tscore\tprobability\n");
    for f in found {
        for s in &f.sites {
            text.push_str(&format!(
                "{}\t{}\t{}\t{}\t{:.3}\t{:.3}\n",
                f.motif.name.replace(' ', "_"),
                names[s.sequence].replace(['\t', '\n'], " "),
                s.position + 1,
                if s.reverse { '-' } else { '+' },
                s.score,
                s.probability
            ));
        }
    }
    text
}

fn read(path: &str) -> Result<String, String> {
    std::fs::read_to_string(path).map_err(|e| format!("Could not read {path}: {e}"))
}

fn write(path: &str, text: &str) -> Result<(), String> {
    std::fs::write(path, text).map_err(|e| format!("Could not write {path}: {e}"))
}

fn work(args: &[String]) -> Result<(), String> {
    let a = parse(args)?;
    let input = a.input.ok_or("No sequences given: use -i FILE.")?;
    let (names, seqs) = read_sequences(&read(&input)?)?;
    let started = Instant::now();
    let found = discover(&seqs, &a.settings)?;
    eprintln!(
        "{} sequences, {} bases; {} of {} motifs found in {:.2} s",
        seqs.len(),
        seqs.iter().map(Vec::len).sum::<usize>(),
        found.len(),
        a.settings.motifs,
        started.elapsed().as_secs_f64()
    );
    for f in &found {
        eprintln!(
            "{}: {} ({} columns, {:.1} bits), {} of {} sequences have it, {:.1} bits better than the background",
            f.motif.name,
            f.motif.consensus(),
            f.motif.width(),
            f.motif.information(),
            f.sites.len(),
            seqs.len(),
            f.gain
        );
    }
    let meme = meme_text(&found, &background(&seqs));
    if let Some(path) = &a.sites {
        write(path, &sites_table(&found, &names))?;
    }
    match a.output {
        Some(path) => write(&path, &meme),
        None => {
            print!("{meme}");
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
    use super::super::Site;
    use super::*;
    use crate::seq::motif::Motif;

    fn words(s: &str) -> Vec<String> {
        s.split_whitespace().map(str::to_string).collect()
    }

    #[test]
    fn the_options_are_read_and_the_defaults_are_the_windows() {
        let a = parse(&words("-i seqs.fa -o m.meme --sites s.tsv -w 12 -n 3 --one-strand --seeds 5 --iterations 20 --effort 4"))
            .expect("options");
        assert_eq!(
            (a.input.as_deref(), a.output.as_deref(), a.sites.as_deref()),
            (Some("seqs.fa"), Some("m.meme"), Some("s.tsv"))
        );
        let s = a.settings;
        assert_eq!((s.width, s.motifs, s.seeds, s.iterations), (12, 3, 5, 20));
        assert!(!s.both_strands && s.effort == 4.0);
        // the defaults are those of the Motif Finder, and the last bare argument is the input
        let d = parse(&words("seqs.fa")).expect("options");
        let w = Settings::default();
        assert_eq!(d.input.as_deref(), Some("seqs.fa"));
        assert_eq!(
            (d.settings.width, d.settings.motifs, d.settings.effort),
            (w.width, w.motifs, w.effort)
        );
        assert!(d.settings.both_strands);
    }

    #[test]
    fn a_wrong_option_is_said_and_refused() {
        let said = |text: &str| parse(&words(text)).err().expect("refused");
        assert!(said("-i").contains("needs a value"));
        assert!(said("-i s -w wide").contains("needs a number"));
        assert!(said("-i s --nonsense 1").contains("not known"));
        assert!(said("-i s -n 0").contains("at least 1"));
        assert!(said("-i s --effort 0").contains("--effort"));
        assert!(said("-i s --effort 5000").contains("--effort"));
        assert!(said("a.fa b.fa").contains("one argument too many"));
        assert_eq!(run(&words("--sites only.tsv")), 1);
        assert_eq!(run(&words("--help")), 0);
    }

    #[test]
    fn the_places_found_are_a_table_counted_from_one_with_the_strand() {
        let motif = Motif::from_sites("my motif", &[&b"GATTACA"[..]; 2]).expect("a motif");
        let found = [Found {
            motif,
            sites: vec![
                Site {
                    sequence: 1,
                    position: 0,
                    reverse: false,
                    score: 9.5,
                    probability: 0.99,
                },
                Site {
                    sequence: 0,
                    position: 41,
                    reverse: true,
                    score: 7.25,
                    probability: 0.6,
                },
            ],
            gain: 3.0,
        }];
        let names = ["first\tone".to_string(), "second".to_string()];
        assert_eq!(
            sites_table(&found, &names),
            "motif\tsequence\tposition\tstrand\tscore\tprobability\n\
             my_motif\tsecond\t1\t+\t9.500\t0.990\n\
             my_motif\tfirst one\t42\t-\t7.250\t0.600\n"
        );
        assert_eq!(sites_table(&[], &names).lines().count(), 1);
    }
}
