//! The routes as pipeline files: `helixloom pipeline example NAME` prints each (`docs/ROUTES.md`).

use super::{AsmStep, Cut};
use crate::pipeline::flow::{Kind, Pipeline, Step};

fn trim() -> Step {
    Step::Asm(AsmStep::Trim {
        adapter: "AGATCGGAAGAGC".into(),
        quality: 20,
        min_length: 20,
        cut: Cut::Tail,
    })
}

fn assemble() -> Step {
    Step::Asm(AsmStep::Assemble {
        ks: "auto".into(),
        cutoff: "auto".into(),
        min_qual: 10,
    })
}

fn hybrid() -> Step {
    Step::Asm(AsmStep::Hybrid {
        long: "long-reads".into(),
        mode: "normal".into(),
        min_quality: "auto".into(),
        single_depth: "auto".into(),
        overlap: 0,
    })
}

/// The Galaxy route for plasmids from short and long reads, whole: the short reads cleaned, assembled into a graph, and the graph joined with the
/// long reads of the folder `long-reads` (a file of the same name before the dot as each short-read file).
pub fn plasmid() -> Pipeline {
    Pipeline {
        name: "Plasmids from short and long reads: trim, assemble, join with the long reads".into(),
        input: Kind::Reads,
        steps: vec![trim(), assemble(), hybrid()],
    }
}

/// The first half of the Galaxy route, which ends at the graph: the place where a person looks at it (`helixloom gfa segments`) and decides what
/// is a repeat, before the second half.
pub fn plasmid_short() -> Pipeline {
    Pipeline {
        name: "Plasmids, part 1: trim the short reads and assemble the graph".into(),
        input: Kind::Reads,
        steps: vec![trim(), assemble()],
    }
}

/// The second half, from a folder of graphs.
pub fn plasmid_long() -> Pipeline {
    Pipeline {
        name: "Plasmids, part 2: join the graph with the long reads".into(),
        input: Kind::Gfa,
        steps: vec![hybrid()],
    }
}

/// The long-read route for plasmids (NanoSeq's order): Flye, the concatemers collapsed, the circles closed, turned to a start gene, polished.
pub fn nanopore_plasmid() -> Pipeline {
    Pipeline {
        name: "Plasmids from nanopore reads: Flye, concatemers, circles, start gene, Medaka".into(),
        input: Kind::Reads,
        steps: vec![
            Step::Asm(AsmStep::Flye {
                mode: "nano-raw".into(),
                genome_size: "auto".into(),
            }),
            Step::Asm(AsmStep::Concatemer {
                min_identity: 95,
                min_unit: 200,
            }),
            Step::Asm(AsmStep::Circularize {
                min_overlap: 30,
                min_identity: 95,
            }),
            Step::Asm(AsmStep::Rotate {
                genes: "start-genes.fasta".into(),
            }),
            Step::Asm(AsmStep::Medaka {
                reads: "reads".into(),
            }),
        ],
    }
}

/// The long-read selection of the multiplex-MinION paper, for a folder of reads that are already sorted by barcode: the reads that go to
/// `hybrid` in `plasmid_long`.
pub fn minion_reads() -> Pipeline {
    Pipeline {
        name: "Multiplex MinION: select the long reads of each barcode (2 kb, 500 Mbp, a 100 base window)".into(),
        input: Kind::Reads,
        steps: vec![Step::Asm(AsmStep::Filter {
            min_length: 2000,
            target_bases: 500_000_000,
            window: 100,
        })],
    }
}

/// The names `pipeline example` takes, with what each is.
pub const NAMES: [(&str, &str); 5] = [
    (
        "plasmid",
        "short and long reads to plasmids in one go (the Galaxy route without its human step)",
    ),
    (
        "plasmid-short",
        "part 1 of the Galaxy route: trim the short reads, assemble the graph",
    ),
    (
        "plasmid-long",
        "part 2 of the Galaxy route: join a folder of graphs with the long reads",
    ),
    (
        "nanopore-plasmid",
        "plasmids from nanopore reads alone (the NanoSeq route)",
    ),
    (
        "minion-reads",
        "select the long reads of each barcode (the multiplex-MinION route)",
    ),
];

/// The route with this name.
pub fn named(word: &str) -> Option<Pipeline> {
    match word {
        "plasmid" => Some(plasmid()),
        "plasmid-short" => Some(plasmid_short()),
        "plasmid-long" => Some(plasmid_long()),
        "nanopore-plasmid" => Some(nanopore_plasmid()),
        "minion-reads" => Some(minion_reads()),
        _ => None,
    }
}

#[cfg(test)]
mod tests;
