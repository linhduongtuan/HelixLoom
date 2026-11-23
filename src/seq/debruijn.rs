//! A de Bruijn assembler for short reads in the manner of SPAdes (Bankevich et al., J Comput Biol 2012, 19:455; Nurk et al., J Comput Biol 2013,
//! 20:714; the SPAdes manual), for **small data, plasmids and small genomes, from reads that are not paired**. Written from those papers and the
//! manual; SPAdes is GPL licensed, and none of its source was read or copied. What is here, and what is the program's own choice, is in
//! `docs/DEBRUIJN.md`.
//!
//! * `pack`: k-mers packed two bits a base into 128 bits (k up to 63, odd), canonical forms, and the counting of the k-mers of the reads;
//! * `unitig`: the graph of the k-mers with every unbranched path made one segment (the condensed graph of the paper);
//! * `graph`: the bidirected graph of segments, links of k-1 bases, and the merging of paths that the removal of an edge left unbranched;
//! * `clean`: tip clipping, bulge removal and the depth filter;
//! * `assemble`: the steps in order, for one k or for several, each k starting from the contigs of the one before;
//! * `reads`: FASTA or FASTQ reads, cut at unknown or poor bases;
//! * `output`: contigs as FASTA (named as the SPAdes manual names them) and the graph as GFA;
//! * `cli`: `helixloom dbg`.

pub mod assemble;
mod clean;
pub mod cli;
pub mod graph;
mod output;
pub mod pack;
pub mod reads;
#[cfg(test)]
mod sim;
#[cfg(test)]
mod tests;
mod unitig;
