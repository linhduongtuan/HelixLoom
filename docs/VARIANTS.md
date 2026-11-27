# Reads against a reference: SAM, consensus and variants

`helixloom variants REFERENCE READS` maps reads to a reference (or reads their SAM), piles them up, and writes the differences as a VCF, the
consensus as FASTA and the depth as a table, for a sample that has **one copy** of the sequence: a plasmid, a bacterial chromosome, a haploid
genome. `helixloom map -a` writes the mapping as SAM. This is the first stage of the plan "Reads and variants" of `docs/ROADMAP.md` 2.9, and
the second use of the mapper of `docs/MAP.md` (stage 3 of 2.7).

```
helixloom variants [options] REFERENCE READS     VCF out; --consensus FILE, --depth FILE, --sam FILE
helixloom map -a [options] TARGET READS          SAM in place of PAF
```

`REFERENCE` is FASTA. `READS` is FASTA or FASTQ (not compressed; mapped here with the mapper of `helixloom map`, with its base-level alignment), or
a SAM file (text, not BAM). A sequence is named by the first word of its header.

## What it does, and what it was written from

**Written from** the SAM and VCF specifications of the SAM/BAM working group, the paper that the approach follows (Li, *Bioinformatics* 2011,
27:2987, read in full for the parts below, with the formulas as images that **could not be read**), and the manual pages of `bcftools mpileup` and
`bcftools call` (read for the defaults and the tags). **No source of samtools or bcftools was read or copied**, and neither program was run:
nothing here claims the same calls as theirs.

1. **SAM (`variants::sam`).** A mapping becomes a line: flag 16 for the reverse strand, 256 for a secondary (which has no bases, as the
   specification allows), 2048 for each further place of a read that lies in several (a chimera, or a read across the start of a circle), `NH`
   for the number of lines of a place, the position from 1, the mapping quality of the mapper, a CIGAR of the aligned bases with the **ends that were not aligned
   as soft clips** (the mapper aligns from the first anchor to the last), the read as it lies against the reference (the reverse complement for the
   other strand, qualities turned round), and the tags `NM`, `tp`, `cm`, `s1` and `dv`. A read that does not map has a line of its own (flag 4).
   The reader takes the `@SQ` lines, refuses a line with fewer than 11 columns, a flag or position that is not a number, a CIGAR that is not runs
   of a number and an operation, a CIGAR whose bases are not the read's, and qualities that do not match the bases.
2. **The pileup (`variants::pileup`).** For each position: how many reads have each of A, C, G and T and the sum of the evidence that they
   give, the reads that have a deletion there, and the insertions that reads have after it. An alignment counts when it is mapped and not a
   secondary, not a duplicate and not marked as failing the sequencer's checks (so a primary, or a supplementary, the other part of a split read,
   as `bcftools mpileup` counts it by default), with a mapping quality of at least `-q` (0); a base counts when its quality is at least `-Q` (13);
   at most `-d` (250) reads count at a position. These defaults are those that the `bcftools mpileup` manual gives. **A base's quality is held at
   30**, because the errors of one read are not independent of the others' and a very high quality would make a site look surer than it is. Base
   alignment quality (BAQ) is not done.
3. **A substitution.** If the sample has the base `a`, a read base of error chance `e = 10^(−q/10)` is seen with chance `1 − e` when it is `a` and
   `e/3` for each other base, and the reads are independent. The log-likelihood of `a` is then a constant plus, for each read that shows `a`,
   `ln(3(1 − e)/e)` (never below 0). The prior of the reference base is `1 − π` and each other base has `π/3` (`--prior`, 0.001). The call is the
   base of the best posterior, and the quality is `−10 log10` of the posterior of the reference base, held at 999. **The model is this program's
   own**: the paper's formulas (which add a correction for the dependence between errors, and a prior that is the allele-frequency spectrum of a
   population) could not be read, and a haploid sample has no population.
4. **An insertion or a deletion.** The reads give no quality for these. With `n` of `N` reads over the position showing it and a chance
   `ε` (`--indel-error`, 0.1) that a read shows an indel that is not there or misses one that is, the log-likelihood ratio of "the sample has it"
   to "it has not" is `(2n − N) ln((1 − ε)/ε)`, and its prior is `ρ` (`--indel-prior`, 0.0001). A run of deleted positions is one deletion, and its
   quality, depth and counts are those of its weakest position. Of the insertions after a position the commonest is taken. A deletion is written
   as VCF writes it, with the base before it (or after it, at the very start); an insertion after the base before it.
5. **The consensus.** The reference with the calls made applied (quality at least `--min-qual`, 20): a base that is called, a deleted run left out,
   an insertion put in, and **`N` where fewer than `-D` (5) reads have a base** and no deletion was called. Where the reference is `N`, the base
   with the most evidence is written only when it is clearly ahead of the next (by the lowest quality), otherwise `N`. U in a read or a reference is
   read as T; a SAM line with position 0 lies on no sequence; an insertion that follows a skip (`N` in a CIGAR) is not counted.
6. **The VCF** (version 4.3) has one haploid sample: `GT:DP:AD` with `1`, the reads at the site and the reads with the reference and with the
   alternate. FILTER is `PASS`: a call under the lowest quality is not made at all. `--depth` writes the reads over each position.

## What was measured

By the tests (83 for this stage: 82 in the module and one for `map -a`): the weight of a base against a hand-worked value; the pileup on records made by hand (the filters, the cap, a
deletion, an insertion that follows a base and does not follow a deletion or the start of a read, soft clips, skips); the calls on hand-made
pileups with numbers worked out apart from the code (a substitution of ten reads at quality 20 is 212.5; two reads at quality 13 are 3.18; of 20
reads, 13 showing a deletion is 17.3 and 14 is 36.3; five reads over a deletion are 8.4); the mapper's hits written as SAM and read back (both
strands, the clips, the secondary, the supplementary of a split read, the read that did not map); and a run of the command on reads simulated from a sample with a substitution, a
deletion of two bases and an insertion of one, which found those three, nothing else, and the consensus of the sample (the `N` at the two ends
where few reads reach aside), and gave the same VCF from the SAM that it had written.

**An accuracy sweep, with a script that was not kept** (reads simulated from a 3,000-base random reference with two substitutions, a deletion of
two bases and an insertion of one, 1,500 bases long, with the same chance of a substitution, an insertion and a deletion in each base, all bases of
one quality; ten references to a setting; counts are over the 40 differences planted in each setting):

| each kind of error | read quality | depth | found of 40 | false calls | the consensus |
|---|---|---|---|---|---|
| 1 % (3 % in all) | 20 | 15× | 39 | 0 | right but for one run that missed a difference |
| 1 % | 20 | 30× and 60× | 40 | 0 | right |
| 2 % (6 % in all) | 15 | 15× | 34 | 0 | wrong in the runs that missed a difference |
| 2 % | 15 | 30× and 60× | 40 | 0 | right |
| 3 % (9 % in all) | 14 | 15×, 30×, 60× | 35, 39, 40 | 0 | right at 60× |
| 4 % (12 % in all) | 14 | 15×, 30×, 60× | 32, 37, 38 | 0 | wrong in the runs that missed one |
| 4 % | 10, `-Q 10` | 15×, 30×, 60× | 31, 36, 38 | 0 | the same |
| 4 % | 10, default `-Q 13` | any | 14 to 18 | 0 | **no base is counted**: they are under the lowest quality, so no substitution is called |

No false call in any of the 150 runs. What is missed at the higher error rates is nearly all insertions and deletions (a missed one makes the
consensus the wrong length, which is what "wrong" means here). The simulated errors are independent and at random places, which real nanopore
errors are not (they gather in homopolymers), so this is the best case.

**The code was broken on purpose** in 103 places, one at a time (a boundary moved, a sign turned, a branch left out, a constant changed, a sort
key, a clamp, a flag): the tests caught 90 at the first run. The 13 that got through were each looked at, and **none was a change that changes
nothing**: each got a test and was caught on a second run. They were: the cap on the quality of an indel (999, not 255); which of two other
bases wins a tie (the earlier of A, C, G, T); the quality at which an indel and a substitution are called, to the last digit; an insertion with one
base that is not A, C, G or T among good ones; the depth of an insertion (the reads that have a base before it, not those that delete it); a
reference longer or shorter than the pileup; a deletion that covers a whole target; the depth table counting every read over a position; the cap
of reads at a position for deletions too; the flag of a secondary on the other strand (it must keep flag 16); the `NM` tag (the bases that differ,
not those that match); and `-k 32` being accepted. One round was started under a heavy load on the machine and was stopped and restarted from the
28th with a lower optimisation level for the program itself, which does not change what a test sees. Four more breakages of the flags of split
reads (added later, `docs/COUNT.md`) were all caught.

## Options

`-o` (the VCF; the standard output when left out), `--consensus FILE`, `--depth FILE`, `--sam FILE` (also write the mapping), `--sample NAME`,
`-q`, `-Q`, `-d`, `-D`, `--min-qual`, `--prior`, `--indel-prior`, `--indel-error`, `--default-quality` (for a read with no qualities: 15), `-k`, `-w`.
`helixloom variants --help` lists them. A summary line for each target (length, mean depth, positions written as `N`, the numbers of each kind
of difference) goes to the standard error.

## Not done, and not known

* **Never run against samtools or bcftools**, and never on real reads: the numbers above are on simulated ones, whose errors are independent.
* **One copy only.** No diploid or mixed samples, no ploidy option, no calling of variants at several frequencies in one sample.
* **No BAQ and no model of dependent errors.** The quality is held at 30 and a base below 13 is not counted. Qualities at high depth are too sure.
* **Insertions are found by the exact bases that reads write.** Reads that write a long insertion in different ways (errors inside it) split it
  among them, and it is found less often. Neighbouring differences are not joined into one event (a substitution beside an insertion is two lines),
  and an insertion is skipped after a position that is deleted.
* **Mapping.** The mapper does not extend the ends of an alignment, so the first and last bases of each read (45 on average here) are soft clips
  and do not count; the depth near the ends of a sequence is lower for it. A repeat is mapped by its best chain only (a secondary is not counted).
* **No BAM, CRAM, VCF compression or index; no region option; no genotype quality or likelihoods (PL) in the VCF; no filters beyond the lowest
  quality, depth and mapping quality.**
* **No pipeline step yet**, and no read-QC report beyond `seq::fastq`'s: those are the next stages of the plan.
