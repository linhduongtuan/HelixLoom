# Counting reads over genes, after featureCounts

`helixloom count` says how many reads fall on each gene (or exon), from SAM files and an annotation. It is the first stage of the plan "Read
counting and long-read assembly, after Subread and NECAT" (`docs/ROADMAP.md` 2.10).

```
helixloom count -a ANNOTATION [-o counts.tsv] [options] reads.sam...
```

`ANNOTATION` is GTF, GFF3 or a SAF table (`GeneID Chr Start End Strand`, positions from 1, both ends counted). The alignments are SAM (not BAM),
for example from `helixloom map -a` (`docs/VARIANTS.md`) or from another aligner. With `-o` the table goes to the file and the summary to the
file with `.summary` after it; without `-o` the table is written out and the summary goes to the standard error.

## Where it comes from, and the licence

featureCounts is part of the Subread package (Liao, Smyth and Shi, *Bioinformatics* 2014, 30:923). Subread is **under the GPL (version 3 or later,
as the Bioconductor page of its R version states)** and HelixLoom is MIT, so **none of its source was read or copied**. This was written from the
chapter on read summarization of the Subread users' guide (read in full) and the paper (read through its methods and results). featureCounts was
**not run**, so nothing says the numbers are the same.

## The rules

* **Features and genes.** The lines of the annotation whose type is `-t` (default `exon`) are features; those with the same value of the attribute
  `-g` (default `gene_id`) are one gene (a meta-feature). A SAF table has both in its columns. A line of the type without the attribute is an
  error that names it. With `-f` each feature is counted on its own.
* **Overlap.** An alignment overlaps a feature when at least one aligned base (`M`, `=`, `X`) lies in it: bases that are deleted, skipped over (a
  junction), inserted or clipped do not count. The overlap with a gene is the number of aligned bases in any of its features, each base once. A read
  that spans two exons of one gene counts once for it.
* **Filters**, in the guide's order; the first that stops an alignment is the one it is counted under: unmapped; mapping quality under `-Q`;
  duplicate (`--ignoreDup`); multi-mapping (`NH` above 1, unless `-M` or `--primary`); secondary (`--primary`); split (`--nonSplitOnly`) or not
  split (`--splitOnly`); no feature; overlap too short (`--minOverlap`, `--fracOverlap` of the read with its soft clips,
  `--fracOverlapFeature`, `--nonOverlap`, `--nonOverlapFeature`: all those asked must be met); ambiguous (more than one target, without `-O`).
  `--largestOverlap` gives a read to the target it overlaps most (two equal ones stay ambiguous).
* **Strand** (`-s`): 0 either (default), 1 the read's strand must be the feature's, 2 the other. A feature with no strand fits either (this
  program's choice).
* **Counts.** One for each target of an assigned alignment. `--fraction` (with `-M`, `-O` or both) divides it by `NH` and by the number of
  targets. The summary counts **alignments**, not reads (as the guide says): a read mapped in several places, or split (a primary and
  supplementary lines), has a line for each. `helixloom map -a` writes the `NH` tag that the multi-mapping rule reads.
* **The output** is the count table (`Geneid Chr Start End Strand Length` and a column for each file; a gene lists the sequence, start, end and
  strand of each feature, joined by `;`, and its length counts each base once) and the summary (`Status` and the categories of the guide that are
  built here).

## What was checked

39 tests: GTF, GFF3 (escapes, a `##FASTA` section) and SAF read, with each fault named; the overlap at its first and last base; aligned blocks
only; a read over two exons; a base in two features of one gene counted once; a long feature found past a short one that starts after it; the
ambiguity, `-O`, `--largestOverlap` and the three strand settings; each filter, their order and the flags that turn them off; every overlap rule
at the exact limit; the fractions; the refusal of paired reads; the table and summary to the letter; and a run on reads simulated inside three genes
and mapped by `helixloom map`, whose counts are the numbers simulated (a read from between genes is not counted).

**The code was broken on purpose** in 61 places, one at a time (57 in the counting code and the line reader it uses, 4 in the SAM flags of split
reads): the tests caught 58 at the first run. One breakage did not compile and was rewritten; two got through and got a test each (a feature
that ends exactly where a read begins, looked at past a longer feature on the other strand; and a share that is exactly met when its product in
binary is a little over the whole number, where the first test had used numbers with no such rounding). All 61 are caught.

A later independent reading found that the warning about sequence names compared every read with every feature (slow exactly when the names do
not match); it now looks the name up.

## Not done, and not known

* **Paired reads are not counted as fragments** (`--countReadPairs`, `-B`, `-C`, `-P`, `-d`, `-D`): reads of pairs are refused unless
  `--each-read` says to count each read on its own, so that nothing counts a pair twice without being asked. This is the biggest gap for RNA-seq.
* **No BAM**, no junction counts (`-J`), no read shifting or extension (`--readShiftType`, `--readExtension5/3`, `--read2pos`), no chromosome alias
  file (`-A`), no groups by read group, no in-built annotations, no negative `--minOverlap`.
* **Every aligned block is looked at** (the program's `--maxMOp` limit is not applied).
* **Never compared with featureCounts**, and never run on real data. Whether it counts supplementary lines the same way is not known: here each
  is an alignment, because the guide lists no filter for them.
