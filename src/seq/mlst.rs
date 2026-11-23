//! Multi-locus sequence typing in the manner of mlst (Torsten Seemann): the alleles of a scheme's loci are looked for in the contigs of a
//! genome, each locus is called, and the calls are a sequence type. Written from mlst's README (mlst is GPL licensed, and none of its source
//! was read or copied), with the program's own aligner and mapper in place of BLAST.
//!
//! * `scheme`: a scheme folder (profile table and one FASTA of alleles for each locus), which the person downloaded: PubMLST's data may not
//!   be bundled;
//! * `search`: where each locus lies in the contigs, on either strand, and how well;
//! * `call`: the call for each locus (`n`, `~n`, `n?`, `-`, `n,m`), the score, the sequence type and its status;
//! * `output`: TSV, CSV, the full table, JSON, and the novel alleles as FASTA named by their MD5.

mod call;
pub mod cli;
mod md5;
mod output;
mod scheme;
mod search;
#[cfg(test)]
mod tests;
