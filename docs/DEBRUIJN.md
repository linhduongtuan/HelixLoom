# A de Bruijn assembler for short reads, after SPAdes

`helixloom dbg` assembles short reads with a de Bruijn graph, in the manner of SPAdes, for **small data: plasmids and small genomes, from reads that
are not paired**. It is the fifth stage of the plan of `docs/ROADMAP.md` 2.7. It was written from the SPAdes papers and manual; SPAdes is GPL
licensed, so none of its source was read or copied. It was not compared with SPAdes, and it was never run on real reads.

```
helixloom dbg [options] READS...           contigs as FASTA (standard output), a summary on the error output
helixloom dbg -o DIR READS...              DIR/contigs.fasta and DIR/assembly_graph.gfa
```

`helixloom debruijn` and `helixloom spades` are other names for the same command; the options are this program's, not SPAdes's (a few are named alike:
`-k`, `-o`, `-1`, `-2`, `-s`, `--cov-cutoff`).

## What it does

1. **Reads** (FASTA or FASTQ, not compressed; `-1`, `-2` and `-s` all name reads, and a pair is used as two sets of reads that are not paired) are cut
   into the runs of bases that are A, C, G or T (U as T) and, for FASTQ, of Phred quality 10 or more (`--min-qual`).
2. **For each k** (default: those of 21, 33 and 55 that are not longer than the longest read; `-k 21,33,55` to choose, odd, rising, 3 to 63):
   the k-mers are counted, canonical (a k-mer and its reverse complement are one), on both strands; the graph of the k-mers is made and every
   unbranched path becomes one segment (the *condensed* graph of the paper: its edges are segments, written out as bases); the graph is cleaned:
   * **tips** are clipped, the shortest first: a segment that goes nowhere at one end, is joined at the other to places that have a way on without
     it, is not longer than the longest read (in k-mers, `--tip-length`), and has a coverage that is not over 2.5 times that of the best way on
     beside it (`--tip-ratio`);
   * **bulges** are removed: of paths that lead from one place to the same other place (each path a single segment that nothing else leads into), a
     path of at most 150 k-mers (`--bulge-length`) that spells nearly the bases of a path with more coverage (the edit distance at most a tenth
     of the longer: `--bulge-diff`) is removed, and the path that stays gains its coverage; as paths are removed, a path that was in the middle of
     others is one segment, and the bulges round it are found in the next round (up to 20 rounds of tips and bulges);
   * the **depth filter** removes the segments whose coverage is under a cutoff (`--cov-cutoff`: a number, `off`, or `auto`, the default), and
     the tips and bulges that this leaves are cleaned again.
3. **The next k** starts from the reads, and from the segments of the graph before, put in as k-mers that have the coverage the segment had
   wherever the reads gave the k-mer no count (the iteration of the paper): a repeat that is longer than the smaller k and shorter than the larger is
   resolved by the larger, and where the reads are too few for the larger k the smaller one's work is kept.
4. **At the last k**, a segment that is joined to nothing and is shorter than 200 bases is removed (`--isolated-length`).

The contigs are the segments of the last graph, longest first, named as the SPAdes manual names them, `NODE_<number>_length_<bases>_cov_<coverage>`,
the coverage being the average number of times a k-mer of the segment was counted at the last k (the manual does not say how it is averaged; this is a
mean over the k-mers). Only one strand of each is written. A segment that is a whole circle (its end leads to its own start, and to nothing else) is
written once round, with the word `circular` after its name. `--min-length` (default 200) keeps shorter contigs out of the FASTA; the numbers do not
change. The GFA (1.0) has every segment (`LN`, `KC`, the k-mer count, and `DP`, the depth, as the manual names them) and a link for each place two
segments meet, overlapping by k-1 bases; there are no paths.

## What is from the papers, and what is this program's choice

Only part of the papers could be read in full: of Bankevich et al. 2012 the abstract, sections 1 to 4.4 and figures 1 to 3; of Nurk et al. 2013 the
abstract and sections 1, 2.1 to 2.5 and the start of 3.1. The rest was read through a summary, so what comes from there is second-hand and is marked
(second-hand).

* **From the text read**: the de Bruijn graph of the reads, the k-mers as edges (§3.2); the condensed graph, its edges being the unbranched paths; the
  coverage of an edge as how often its k-mers occur in the reads; reads and their reverse complements used together; every change made to both strands at
  once (Nurk §2.2); small k merges repeats, large k breaks low coverage (Bankevich §3.3); tips and bulges as the false structures that errors make
  (Nurk §2.2, §3.1). From the manual (output page): the contig names, the k-mer coverage of the largest k, the GFA tags DP and KC, the default
  lists of k (21, 33, 55 for reads of 100 bases; more for longer reads).
* **Second-hand**: that a tip is at most as long as a read, with a relative coverage test of 2.5 for ordinary data and none for single cells; that a
  bulge is at most 150 long; that the graph for the next k is made from the reads and the old edges; that short isolated paths go last (about 200).
  These are used as stated, in k-mers, and may be wrong in detail.
* **This program's choice** (the papers do not say, or were not read to it): the **similarity test of a bulge** (edit distance at most a tenth of the
  longer path, any path with more coverage that is alike, not only the best); the **automatic cutoff** (a tenth of the coverage that half the bases are at
  or above, among the segments of two k or more, at least 1.5 and at most half of that coverage); the quality limit 10; the 200 of `--min-length`; the
  coverage a carried k-mer is given (the coverage of its segment at the smaller k); k to 63 and not 127 (a k-mer is held in 128 bits); the order of
  contigs; that a circle is written once round and marked; the names of the options. The papers' erroneous-connection removal with its
  increasing cutoff, and the flow method that finds the chimeric edges of single-cell data (Nurk §2.2 to §2.4), are **not done**.

## What was measured

On sequences of random bases and reads simulated from them (substitution errors only, both strands, no pairs, no quality values):

* perfect reads of a line give the line back, at one k, at three, and with the default; a circle gives one circle of its own length; with 0.5 % errors
  it is still one circle; with 1 % errors at 60 times, a line of 5,000 bases comes back whole and nothing else is left;
* two copies that differ by one base (a bulge) are one contig, which has the coverage of both; a chromosome of 3 kb at 30 times and a plasmid of
  1.5 kb at 90 times are two circles, the plasmid's coverage three times the other's;
* a repeat of 28 bases is not resolved at k = 21 and is at k = 33 (and at 21, 33); a repeat of 300 bases, longer than the reads, is **not** resolved, and
  the graph shows the repeat linked to what comes before and after it;
* a stretch seen at 1.6 times beside a genome at 40 is taken away by the automatic cutoff and stays with `--cov-cutoff off`; the same reads, or the same
  reads turned to the other strand, give the same contigs, and so do two runs.

A throwaway test (not kept in the suite), on random genomes with no repeats, one run each, on a Mac that other work was loading: 100 kb at 50 times with
0.3 % errors, at 50 times with 1 % errors and at 100 times with 0.1 % errors, each a single contig equal to the genome in 0.5 to 2.6 seconds; at 10 times,
two contigs, 99.98 % of the genome's 31-mers, no wrong 31-mer; and 1 Mb at 40 times with 0.3 % errors, a single contig equal to the genome, in 8.6
seconds. The built command was tried once on two FASTQ files of simulated pairs (read as unpaired) from a circular chromosome of 20 kb at 40 times and a circular
plasmid of 3 kb at 120 times, with 0.3 % errors: two circular contigs, each equal to its sequence up to where the circle is cut and the strand, with
k-mer coverage 15.6 and 47.1, in half a second. These are the easiest cases (no repeat, no bias in the reads, no indel errors) and say nothing about real data.

**Deliberate breakages.** 94 one-line changes were made to the code on purpose, one at a time (a bound off by one, a sign, a strand, a comparison, a
branch removed), to see whether the tests notice. 80 were caught and 14 were not. Of the 14: nine were closed with tests (a tip that another
segment can only reach through; a hairpin that is not a tip; the share of bases two paths may differ by, rounded down; the k-mer minimum itself;
one segment cut away and the path it left made whole again; the rounding of the k-mer count in the GFA; N50 at the exact half; a cutoff of 0),
**one was a fault found by it** (a bulge whose two paths both lead back to the segment they left, which is every bulge in a circle with one
place of difference, was never removed: the guard against it is gone, and a circle of two copies that differ by a base is one circle), one was
dead code (the walk that stops at the start of a circle also stops at a k-mer passed before, so the check for coming round was removed), and
**four are left, which change nothing that was found**: a guard in the merging of runs against a segment joined to itself (the merging stops at a
segment it has placed anyway), the comparison of lengths before the edit distance (a speed-up: the distance is never under the difference in length),
the skipping of places with fewer than two ways on in the bulge search (also a speed-up), and a guard against comparing a path with itself in a group
(reachable only by a segment that is in the group in both directions, which no test makes). The 13 changes that touch what was changed were run again and all
caught.

## Not done, and not known

* **Never run on real reads, never compared with SPAdes** (or any assembler). What it does with Illumina errors that are not uniform, with indels, with
  uneven coverage, with real repeats and with real plasmids is not known.
* **No error correction** (BayesHammer's, which SPAdes runs first), so the graph of reads with errors is far larger than the genome before it is cleaned
  (at k = 21, 1 % errors, 60 times: about twenty times); the cleaning gets the answer, but the work and the memory are those of the larger graph.
* **No paired reads** (a pair is two reads), so a repeat longer than a read is not resolved; **no scaffolds**; no paths in the GFA; no FASTG; no
  `contigs.paths`; no read length over what the k allows (k to 63).
* **The whole table of k-mers is in memory**, about 32 bytes for each distinct k-mer and the graph on top: bacterial genomes, plasmids and small
  data, not a human genome. The counting is shared among the cores; the cleaning is one thread.
* The depth filter removes any segment under the cutoff, as the manual's `--cov-cutoff` does, and not only the erroneous connections; a stretch of a
  true genome with locally low coverage is cut out. A short plasmid (under 200 bases) is taken as an isolated fragment and removed at the last k
  (`--isolated-length 0` keeps it).
* Compressed reads, the Windows and Linux builds (not tried), and any interface other than the command line.
