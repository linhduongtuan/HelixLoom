# The assembly routes as pipeline files

Three published routes to a finished plasmid or small genome, written as **pipeline files** (`docs/PIPELINES.md`): a short text file of typed steps,
run on every file of a folder, that does not repeat work it has done. It is the seventh stage of the plan of `docs/ROADMAP.md` 2.7, and it joins the
built-in tools of stages 2 to 6 (reads cleaned, `helixloom dbg`, `helixloom hybrid`) with the programs that stay outside (Flye, Medaka).

```
helixloom pipeline example plasmid > plasmid.txt            the Galaxy route, whole
helixloom pipeline example plasmid-short > part1.txt        its first half: trim, assemble the graph
helixloom pipeline example plasmid-long > part2.txt         its second half: join the graph with the long reads
helixloom pipeline example nanopore-plasmid > nanopore.txt  the NanoSeq route
helixloom pipeline example minion-reads > minion.txt        the long-read selection of the multiplex-MinION route
helixloom pipeline check plasmid.txt                        what it does, step by step
helixloom pipeline run plasmid.txt ./short-reads --cache ./cache --out-dir ./plasmids
```

A **folder holds one sample in each file**: one FASTQ of the short reads (join the two files of a pair with `cat`: pairs are not understood), or one FASTQ
of long reads. A step that needs more than what flows to it names a file or a folder in a setting (`long=`, `reads=`, `genes=`): the folder holds, for
each input file, a file of the same name before the dot (`sample1.fastq` goes with `long/sample1.fastq`). A path has no space in it. What each step
found is written as a word after the name in the FASTA: `circular=true`, `concatemer=3`, `rotated=true`. The answers are kept in the cache, and the key of one
holds the checksum of the file that the step read, so changing the long reads makes the joining run again and not the assembly before it.

## The steps

| Step | Takes | Gives | Settings (the defaults) | What does it |
|------|-------|-------|-------------------------|--------------|
| `trim` | reads | reads | `adapter=AGATCGGAAGAGC` (or `none`), `quality=20`, `min-length=20`, `cut=tail` (`right`, `off`) | fastp's rules (`seq::fastq`) |
| `filter` | reads | reads | `min-length=2000`, `target-bases=500000000`, `window=100` | the selection of the multiplex-MinION paper (`seq::readfilter`) |
| `assemble` | reads | gfa | `ks=auto` (or `21,33,55`), `cov-cutoff=auto` (or `off`, a number), `min-qual=10` | `helixloom dbg` (`docs/DEBRUIJN.md`) |
| `hybrid` | gfa | sequences | `long=PATH` (required), `mode=normal` (`conservative`, `bold`), `min-quality=auto`, `single-depth=auto`, `overlap=0` | `helixloom hybrid` (`docs/HYBRID.md`) |
| `flye` | reads | sequences | `mode=nano-raw` (`nano-corr`, `nano-hq`), `genome-size=auto` (or `5m`, `60k`) | **the installed Flye** |
| `concatemer` | sequences | sequences | `min-identity=95`, `min-unit=200` | `seq::circular::find_concatemer` |
| `circularize` | sequences | sequences | `min-overlap=30`, `min-identity=95` | `seq::circular::end_overlap` |
| `rotate` | sequences | sequences | `genes=PATH` (required) | `seq::hybrid::rotate` (`docs/HYBRID.md`) |
| `medaka` | sequences | sequences | `reads=PATH` (required) | **the installed Medaka** |

`trim` cuts the adapter at the end of a read, then the end of low quality (a window of four bases whose mean is under `quality`), and drops the reads that are
then too short or hold too many poor or unknown bases; it is **fastp's** way and not Trim Galore's. `flye` and `medaka` are called by name like MAFFT and
IQ-TREE are: if the program is not installed the step fails and says so, and nothing else is used instead (`helixloom tools` says how to install them).
**Neither Flye nor Medaka is installed on the computer this was built on, so neither was run**: the command lines follow what each program's documentation
gives (`flye --nano-raw reads.fastq --out-dir flye-out --threads N [--genome-size S]`, then the contigs of `assembly.fasta`, those that `assembly_info.txt`
marks circular getting `circular=true`; `medaka_consensus -i READS -d draft.fasta -o medaka-out -t N`, then `consensus.fasta`), and the glue is tried on
stand-in programs that write the files the real ones write. Canu is not wrapped: a Canu assembly can be the input of a pipeline that starts at
`input: sequences`.

`concatemer` cuts a contig that is several head-to-tail copies of one sequence to the first of its equal parts. It tries 10 copies down to 2 and takes the
most for which each other part maps to the first (the program's mapper) over 90 % of its length, at least 95 % alike (the mapper's estimate of the
divergence), the unit at least 200 bases. An assembly of 2 copies and 60 bases more has a unit of 3,030 bases for a plasmid of 3,000: the 30 bases that
repeat are what `circularize` then finds. **`circularize`** compares the last bases of a contig with its first, base to base, for every length from 30 to 5,000
(never more than half the contig), and where they are 95 % alike (the most alike, the longer if they tie) cuts the overlap off and says `circular=true`;
one that says it already is left alone. **`rotate`** turns each circle to start at the best match of the nucleotide genes in the file (`docs/HYBRID.md`).

## Route 1: short and long reads, after the Galaxy workflow for plasmid collapse

The workflow is the one of Schäfer et al. (*Microb Genom* 2024, "A practical guide and Galaxy workflow to avoid inter-plasmidic repeat collapse and false gene
loss in Unicycler's hybrid assemblies"); its repository is GPL licensed, so only its README was read (by someone else, on 25 November 2026), **and the paper
itself could not be fetched** (it answered with a challenge page). What the README says: the reads are checked and filtered (FastQC and Trim Galore for the
short reads, NanoStat and NanoFilt for the long), a short-read assembly is made as a first look at the graph (and BUSCO), then Unicycler makes the hybrid
assembly; Unicycler can collapse a repeat that is on several plasmids into one sequence, so **the multiplicities of the segments often have to be corrected
by a person**, who looks at the SPAdes graph, and the hybrid assembly is run again from the corrected graph (which, the README says, Galaxy cannot yet do);
the result is evaluated with BUSCO, a length computation, BWA-MEM and samtools. The README gives no cut-off, no length, no Unicycler mode and no depth.

* **`plasmid`**: `trim`, `assemble`, `hybrid long=…`: the whole route with no human step.
* **`plasmid-short`** and **`plasmid-long`**: the same cut in two where the person is. Run the first on the short reads (`--out-dir graphs`), look at the
  graphs (`helixloom gfa segments graphs/sample1.gfa` lists the depth, the k-mers and the links of each segment), and run the second on the folder of graphs.
  **What a person can change here is what the program lets them: `single-depth=` (the depth of one copy, when the plasmid is not at the chromosome's copy number),
  `mode=` and `min-quality=` of `hybrid`; the multiplicity of one segment cannot be edited, and the graph is not given back to Unicycler.** The program estimates the
  multiplicities (`docs/HYBRID.md`) and does not collapse repeats of different plasmids by itself: that is what a wrong estimate does, and it is what the person
  is there to catch.
* Not in it: FastQC, NanoStat, NanoFilt (`filter` is the other paper's selection, below), BUSCO, BWA-MEM, samtools. `gfa-inspect` (`docs/PIPELINES.md`) summarises a graph.

## Route 2: nanopore reads alone, the order of NanoSeq

**The README of NanoSeq could not be found again to be read for this**: web searches found other pipelines (a Nextflow pipeline for synthetic-biology constructs,
a clone-validation workflow from Oxford Nanopore, a Snakemake plasmid workflow) and none that matched. The order here is that of the roadmap's note of 19 November,
which says it was taken from the README: Canu or Flye, a concatemer check, a circularisation check (a 30 base overlap at 95 % identity), rotation, Medaka. Those
numbers, the order and what each check means are **not verified against that README**; **how each check is done is this program's own** (above), not NanoSeq's.

* **`nanopore-plasmid`**: `flye`, `concatemer`, `circularize min-overlap=30 min-identity=95`, `rotate genes=start-genes.fasta`, `medaka reads=reads`. Give the genes
  (the *dnaA* or *repA* of the species, as nucleotides) and the folder of the same reads for Medaka (the example names `reads`; change it). To start from a Canu assembly,
  take the steps from `concatemer` on, with `input: sequences`.

## Route 3: multiplex MinION, after Wick et al. 2017

The paper (*Microb Genom* 3:e000132, read in full as text; its supplement was not) completes bacterial genomes from 12 barcoded isolates on one MinION flow cell
with Illumina reads of each. Its method: Albacore base-calls and sorts by barcode; Porechop trims and sorts by barcode too, and **a read is kept only where the
two agree**; each barcode's reads are brought to 500 Mbp (88 times) by dropping the reads under 2 kb, and then, if there are more, the reads with the worst stretch
of quality (the mean Phred over a window of 100 bases) until 500 Mbp are left; Unicycler (v0.4.0) assembles each isolate from its Illumina and its MinION reads and
polishes with Pilon; a replicon is complete when it is one circular contig, and one that is not is finished by hand with Bandage and BLAST. The paper does not say which
Unicycler mode or k it used and gives no barcode threshold.

* **`minion-reads`** is the selection: `filter min-length=2000 target-bases=500000000 window=100`, for a folder of reads that are **already sorted by barcode**
  (basecalling and demultiplexing, and the agreement of two demultiplexers, are not done: `docs/ROADMAP.md` does not plan them). How a read is scored by its worst
  window, and that of two as bad the later goes first, are this program's own.
* Then **`plasmid-short`** (the Illumina reads of each barcode, the same name before the dot) and **`plasmid-long`** with `long=` the folder that `minion-reads` wrote.
  Pilon's polishing and the finishing by hand are not done.

## What was measured

* The steps one by one, on reads and contigs made of random bases: the adapter, the quality tail, the minimum length, the selection of long reads by length and by
  worst window (with its ties), a concatemer of two to five copies with a few changes, one that is not, a circle with an overlap of 30 bases and of 25, with one base
  changed and with two, a start gene on either strand; the settings (every one that is refused, with the step and the setting named), the files a step reads
  (a file, a folder, none, a path that is not there), and the key of a cached answer changing with the file.
* **Route 1 whole, in the test suite**: reads of a circle `X R Y R Z` of 11.4 kb (a repeat of 1.2 kb twice) cut as 100-base reads, 40 times, and long reads of 5 kb with
  6 % errors; `trim`, `assemble` and `hybrid` give one circular replicon of 11,400 bases, equal to the circle up to where it is cut; a second run finds all three
  answers in the cache; other long reads make only the joining run again.
* **From the command line** (not in the suite, one run): an 8.6 kb plasmid with two copies of an 0.8 kb repeat, 3,440 short reads at 40 times and 35 long reads of 2 to 12 kb at 30 times
  with 8 % errors: `helixloom pipeline run` gave one circular replicon of 8,600 bases, equal to the plasmid, in 0.07 seconds, and a second run found all three steps in the cache.
* The checks of the NanoSeq order (a concatemer, a circle, a start gene) on one made-up contig of two copies and 60 bases more: one contig of 3,000 bases that starts at the gene.
* Flye and Medaka on **stand-in programs** that record their arguments and write the files the real ones write: the arguments, the circles read from `assembly_info.txt`,
  the names kept through Medaka, a program that stops, writes nothing or gives the wrong number of contigs.

**Deliberate breakages.** 70 one-line changes were made to the new code on purpose, one at a time (a bound, a sign, a default, a branch removed, a name swapped), to see whether
the tests notice. 57 were caught and 13 were not. Ten of the 13 were closed with tests (a tie between overlaps going to the longer, a unit of exactly the shortest, an adapter of
three bases, a cutoff that is not finite, a window whose mean is exactly the quality, the cutoff and the mode reaching the assembler and the joining, a Medaka that gives more
contigs, the `.fq` ending, the file names of the long-read route) and caught on a second run. **Three are left, which change nothing that was found**: a map of a unit to the
first that covers exactly the share asked for (an alignment's length cannot be set to the base), the choice of the primary alignment of a unit to the first (a unit is
mapped to one copy of itself), and a target that the kept reads equal exactly (nothing is dropped either way).

## Not done, and not known

* **Never run on real reads**; none of the three routes was compared with the program it follows. **Flye and Medaka were not run.** The numbers of route 2 are the roadmap's.
* **Not exported**: a pipeline with one of these steps is refused by `pipeline export` (no Snakefile or main.nf is written for it), and the workstation reads the files named
  in the settings from its own disk, not from this computer's: run these pipelines here.
* **No evaluation steps** (BUSCO, QUAST, read mapping) and no quality control of reads (FastQC, NanoStat); the names of the replicons say their length, depth and whether
  they are circles, and `gfa-inspect` counts the segments of a graph.
* A sample is one file: paired files are joined first; the same name before the dot is how the long reads, the genes and the reads for Medaka are found.
* The checks of a plasmid (`concatemer`, `circularize`) work on bases that match one for one, so a circle whose overlap has an insertion or a deletion in it is cut at a shorter one
  or not at all, and a contig of copies that differ by more than 5 % is not seen as one. A unit that is cut at `length ÷ copies` carries what the contig had over, divided by the copies.
