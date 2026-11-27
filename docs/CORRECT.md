# Correcting long reads, after NECAT

`helixloom correct` corrects noisy long reads (nanopore) with the reads that overlap them. It is the second stage of the plan "Read counting and
long-read assembly, after Subread and NECAT" (`docs/ROADMAP.md` 2.10), made for plasmids and small genomes.

```
helixloom correct [options] READS...     FASTA out
```

`READS` are FASTA or FASTQ files (not compressed); the corrected reads come out as FASTA, a summary on the standard error.

## Where it comes from, and the licence

NECAT (Chen et al., *Nature Communications* 2021, 12:60) states **no licence** (no licence file, nothing in its README), so **nothing of its code
or text was used**: this follows the method as the paper describes it, read in its open copy (PMC7782737; the supplementary notes were not). NECAT
was **not run**, and the accuracy it reports is not claimed here.

## What it does

For each read (the template):

1. **Candidates** are the other reads that the program's mapper chains to it, best chain score first (the paper chains pairs of k-mers that it
   scores by distance, DDF, which it does not define).
2. **A threshold of identity** for the template: of the 50 best candidates, those whose alignment covers over 60 % of either read are measured, and
   an alignment must have at least their **mean identity less five standard deviations**, and never under 0.5. The paper writes the threshold as
   `d0 − 5D` and combines it with a global error of 0.5 so that the stricter wins; what `d` measures and the formula of `D` could not be read, and
   identity with the standard deviation is the reading under which the two behave as the paper says.
3. **Supports**: the 200 best candidates, in order; one is passed over when every position of its span has more than 12 supports already, and only an
   alignment that meets the threshold counts.
4. **The consensus** (the program's own: the paper uses FALCON's, which it does not describe): at each position that at least 4 supports cover, a
   deletion when more reads delete it than have a base there (the template counted), else the base with the most evidence from the qualities (a tie
   keeps the template's), and an insertion after it when over half of those with a base there insert the same bases.
5. **Two rounds** (`--rounds`): every round but the last writes each read whole, correcting the supported positions, counting each alignment only
   over its longest run of 500-base blocks with at most half errors (the paper aligns blocks of 500 out from a seed and stops at a block over
   half errors); the last corrects the corrected reads again with the whole alignments and writes only the supported stretches of at least 500
   bases (the paper writes only corrected pieces in its second step; the minimum is the program's). A read split into pieces has `_1`, `_2` …
   after its name; a read with no supported stretch is dropped.

## Options

`-o`, `--rounds` (1 to 5, default 2), `--coverage` (12), `--candidates` (200), `--min-support` (4), `--min-length` (500), `-k` (15), `-w` (5).

## What was measured

By the tests (29): the threshold on hand-worked numbers; candidates in order, one per read; an alignment's span, identity and whether it is long;
the coverage limit (13 identical reads taken of 30), the number looked at and the threshold's effect; the part of a record cut out of an alignment,
the run of good blocks at exactly half errors, the vote (support, ties, deletions and insertions at their exact majorities), pieces; and the whole
correction on simulated reads: errors cut from over 8 % to under 1 %, one round keeping every read, a read alone kept by one round and dropped by
two, a read with an unsupported middle written in two pieces.

**A sweep with a script that was not kept**: reads of 2,000 to 4,000 bases simulated from a random genome of 20,000 bases with the same chance of a
substitution, an insertion and a deletion, either strand, each corrected read measured against the bases it came from (the share of aligned columns
that match):

| errors in all | depth | rounds | genome | identity before | after | bases kept |
|---|---|---|---|---|---|---|
| 12 % | 30× | 1 | random | 89.6 % | 98.09 % | 99 % |
| 12 % | 30× | 2 | random | 89.6 % | 99.89 % | 97 % |
| 12 % | 30× | 3 | random | 89.6 % | 99.98 % | 97 % |
| 12 % | 15× | 2 | random | 89.6 % | 99.65 % | 95 % |
| 12 % | 8× | 2 | random | 89.6 % | 99.16 % | 82 % |
| 18 % | 30× | 2 | random | 85.1 % | 98.21 % | 92 % |
| 6 % | 30× | 2 | random | 94.5 % | 100.00 % | 99 % |
| 12 % | 30× | 2 | a repeat | 89.6 % | 99.85 % | 97 % |
| 6 % | 30× | 2 | a repeat | 94.5 % | 99.92 % | 99 % |

"A repeat" is a stretch of 4,000 bases copied to a second place with 1 base in 100 changed. The correction is lower there (99.85 % against 99.89 %,
99.92 % against 100.00 %): **reads of one copy support the other, and some of the differences between the copies are lost.** A run of a 20,000-base
genome at 30× took 1 to 5 seconds on a Mac under load, all rows with the default k of 15. The simulated errors fall at random places and independently; **real nanopore errors do
not** (they gather in homopolymers and repeat), so this is the best case.

**The code was broken on purpose** in 57 places, one at a time: 49 were caught at the first run. One did not compile and was rewritten, one
changed nothing (the mapper gives the best place first, so the code that looked for it was removed), and six got through and got tests (an
alignment long only for the template; the threshold measured on long alignments only; a candidate passed over only when its whole span is
covered; a tie of the template with another base whatever the template's; an insertion in exactly half of the reads; a piece of exactly the
shortest length). All 56 that remain are caught. An independent reading of the code (two reviewers and a skeptic, none of them running it)
found nothing in the correction but the U letters of RNA reads, which are now read as T, and four small faults in the variant caller and the
counter, now fixed and tested.

## Not done, and not known

* **Never run on real reads**, and never compared with NECAT.
* **Small genomes only**: every read and an index of all of them are held in memory, and each read is aligned with up to 200 others.
* **Homopolymers and systematic errors** are not modelled; a vote cannot correct an error that most reads share.
* **Repeats and mixed samples** are blurred, as above: a candidate from the other copy is a support when it passes the threshold.
* The paper's **DDF seeds, DALIGN realignment of hard blocks, trimming of corrected reads before assembly** and its FALCON consensus are not
  built; the corrected reads are not yet assembled (stage 3 of the plan).
