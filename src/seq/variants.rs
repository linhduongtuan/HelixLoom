//! Reads against a reference: the alignments as SAM, the pileup, the consensus, and the differences as VCF, for one copy of the sequence (a
//! plasmid, a bacterial chromosome, a haploid genome). `helixloom variants`, and `-a` of `helixloom map`.
//!
//! * **SAM** (`sam`): the mapper's hits written as lines (soft clips, flags, the strand), and a SAM file read, from the specification;
//! * **the pileup** (`pileup`): for each position, the reads' bases with the evidence of their qualities, the deletions, the insertions;
//! * **the calls and the consensus** (`call`): a substitution by the posterior of the four bases from the qualities, an insertion or a
//!   deletion by the share of the reads that show it, in the manner of the approach of Li (*Bioinformatics* 2011); **the model is this
//!   program's own and says so** (its formulas could not be read), and `bcftools` was not run;
//! * **output** (`output`): VCF 4.3 with one haploid sample, the consensus as FASTA, the depth as a table.
//!
//! Written from the SAM and VCF specifications, the paper and the manuals of `bcftools mpileup` and `bcftools call` (the defaults that the
//! manual gives are used where they are named: a base of quality 13, mapping quality 0, depth 250). No source of samtools or bcftools was read.

pub mod call;
pub mod cli;
pub mod mapped;
pub mod output;
pub mod pileup;
pub mod reads;
pub mod sam;
#[cfg(test)]
pub(crate) mod sim;
