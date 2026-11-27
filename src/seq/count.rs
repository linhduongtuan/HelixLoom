//! Counting reads over genes, after featureCounts of the Subread package (Liao, Smyth and Shi, *Bioinformatics* 2014, 30:923). `helixloom
//! count`.
//!
//! * **the annotation** (`annotation`): the exons (or other types) of a GTF or GFF3 file grouped into genes by an attribute, or a SAF table;
//! * **the rules** (`assign`): which alignments are counted, for which feature or gene, and why the others are not;
//! * **the output** (`output`): the count table and the summary.
//!
//! Subread is under the GPL (version 3 or later, as its R package states), and HelixLoom is MIT: **none of its source was read or copied**.
//! This is written from the users' guide of Subread (the chapter on read summarization) and the paper. featureCounts was not run, so nothing says
//! the numbers are the same. Not built: counting pairs as fragments, BAM, junction counts, read shifting and extension, `--maxMOp` (every
//! aligned block is looked at).

pub mod annotation;
pub mod assign;
pub mod cli;
pub mod output;
