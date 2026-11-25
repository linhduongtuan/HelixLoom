//! The file and standard-input entry points for GFA inspection.

use super::{inspect, segment_table};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

const USAGE: &str = "\
Usage: helixloom gfa inspect FILE
       helixloom gfa inspect -
       helixloom gfa segments FILE
       helixloom gfa segments -

Read an uncompressed GFA 1 graph made of header, segment and link records.
Use - to read standard input. Both commands write tab-separated output.

inspect    Print six graph totals. Segment lengths include overlaps. A component is a
           group connected by links in either direction. A segment with no links is isolated.
segments   List each segment in file order with length, recorded depth (DP), k-mer count
           (KC) and incident link count. A dot means unknown. Each self-link counts once.
";

fn run_report(command: &str, input: &str) -> Result<(), String> {
    let reader: Box<dyn BufRead> = match input {
        "-" => Box::new(io::stdin().lock()),
        path => Box::new(BufReader::new(File::open(path).map_err(|e| e.to_string())?)),
    };
    let mut output = BufWriter::new(io::stdout().lock());
    let output_error = |e: io::Error| format!("standard output: {e}");
    match command {
        "segments" => segment_table(reader, &mut output)?,
        _ => {
            let summary = inspect(reader)?;
            output
                .write_all(summary.to_tsv().as_bytes())
                .map_err(output_error)?;
        }
    }
    output.flush().map_err(output_error)
}

pub fn run(args: &[String]) -> i32 {
    let help = |arg: &str| matches!(arg, "-h" | "--help");
    let (command, input) = match args {
        [flag] if help(flag) => {
            print!("{USAGE}");
            return 0;
        }
        [cmd, input] if matches!(cmd.as_str(), "inspect" | "segments") => (cmd, input),
        _ => {
            eprint!("{USAGE}");
            return 1;
        }
    };
    if help(input) {
        print!("{USAGE}");
        return 0;
    }
    if let Err(why) = run_report(command, input) {
        let source = if input == "-" {
            "standard input"
        } else {
            input
        };
        eprintln!("gfa {command}: {source}: {why}");
        return 1;
    }
    0
}
