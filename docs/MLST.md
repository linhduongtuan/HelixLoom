# Multi-locus sequence typing, after mlst

`helixloom mlst` types a genome against a typing scheme: it finds the alleles of the scheme's loci in the contigs, calls each locus, and gives the
sequence type. It is the fourth stage of the plan of `docs/ROADMAP.md` 2.7, written from the README of mlst (Torsten Seemann; mlst is GPL licensed, so
none of its source was read or copied) with the program's own mapper and aligner in place of BLAST. It was not compared with mlst: that needs mlst
and a scheme, and a scheme is a download.

```
helixloom mlst --scheme FOLDER [options] GENOME.fa...        type against one scheme
helixloom mlst --schemes FOLDER [options] GENOME.fa...       a folder of schemes: the best for each genome
```

## The data is yours to download

**PubMLST's data is not included**, and may not be bundled in software (its terms of 13 May 2026: data submitted since 1 January 2025 is for
non-commercial use and may not be redistributed in software). A scheme is a folder that you downloaded, in the layout that mlst's README gives for
adding a scheme by hand:

* `NAME.txt`: the profile table, tab separated, with a header. The first column is `ST`, then one column for each locus, then optionally other
  columns; `clonal_complex`, `mlst_clade`, `species`, `CC` and `Lineage` are ignored. (If the folder has one `.txt` file and it is not named for the
  folder, it is taken.)
* `LOCUS.tfa` (or `.fa`, `.fasta`, `.fna`) for each locus: nucleotide FASTA, one record for each allele, named `LOCUS_ID` or `LOCUS-ID`. An id may be
  a number or letters.

The scheme's name is the folder's name.

## What it does

1. **Find the places.** The genome's contigs are indexed with the program's mapper (`helixloom map`, 11-mers in windows of 3); every allele of a locus is
   mapped to them on both strands. Alleles that land on one place on a contig are the candidates for that place.
2. **Align the best candidates.** The 8 that share most minimizers with the place are aligned over it, with a margin of contig on each side, by the
   program's aligner, the ends that overhang being free. The stretch of contig that the alignment covers is *what is found*.
3. **Call the place.** If what is found is exactly one of the alleles of the locus, it is that allele (a look-up, so it does not matter which candidate
   found it). Otherwise the closest allele is named, and how much of it was found and how alike it is decide the call.
4. **Call the locus** from its places.

## The calls

As mlst's README writes them:

| Call | Meaning |
| --- | --- |
| `n` | the allele n, intact: all of it, every base the same |
| `~n` | a novel allele: all of the allele present, alike enough (`--minid`); n is the closest known allele |
| `n?` | a partial match to n: some of it (`--mincov` or more), alike enough |
| `-` | missing: too little of it, or too unlike |
| `n,m` | more than one allele at the locus |

The defaults, **which the README does not state and are this program's choice**: `--minid` 95 (percent), `--mincov` 10 (percent of the allele).
Where a locus has several places, it is called by the best kind among them (exact, then novel, then partial) and reports all the places of that kind.

## The score and the sequence type

For a scheme of N loci: 90/N for each exact locus, 63/N for each novel, 18/N for each partial, nothing for a missing one, and 10 more when the alleles
are a known sequence type: a score out of 100. A locus with several alleles scores as exact. The sequence type is given when every locus has one exact
allele and the alleles are a row of the profile table; otherwise it is `-`.

With `--full` the status is PERFECT (a known type), NOVEL (exact alleles that are in no profile), NONE (no allele matches), MIXED (a locus has several
alleles), MISSING (an allele is missing), BAD (none of those and a score under 70) or OK.

With `--schemes`, each genome is typed against every scheme and the one that scores best is used, if it scores at least `--minscore` (50). **The README
does not say how a tie goes**: here the first scheme in alphabetical order of the folders wins. A genome that no scheme types is written as
`FILE`, `-`, `-`.

## Output

One line for each genome, tab separated: the file, the scheme, the sequence type, and `locus(call)` for each locus. `--csv` for commas; `--full` for a
header and the columns FILE, SCHEME, ST, STATUS, SCORE, ALLELES; `--json FILE` for an array of objects (`id`, `filename`, `scheme`, `sequence_type`,
`alleles`); `--novel FILE` for the novel alleles as FASTA, named `scheme.locus-MD5 file` with the MD5 of the sequence; `--nopath`, `--label NAME`,
`-o FILE`.

## What was measured

Everything on made-up schemes of random alleles that differ from each other by a few bases, and genomes of three contigs of random bases with alleles put in
them on both strands (no real scheme is in the repository): the alleles of a known type are called exactly and give the type with a score of 100; exact
alleles that are in no profile give 90 and no type; an allele with three bases changed is `~n`, with the closest named and its bases written; an allele
cut short by the end of a contig is `n?`; one that is not there is `-`; two different alleles at a locus are `n,m` in numeric order; an allele 10 % unlike is
missing until `--minid` is lowered; one unknown base makes a novel allele; an allele split over two contigs is partial; the best of two schemes is chosen
and a tie goes to the first; the scores are worked out by hand (a score that is a whole number, 70 or 100, comes out as exactly that); an allele
that is the middle of another is not a place of its own (the longer one is called); two alleles at places that overlap in coordinates but are on
other strands or contigs stay two places; the closest allele is named among more than the 8 that are aligned; the MD5 equals the RFC's test suite.

**Deliberate breakages.** 61 one-line changes were made to the code on purpose, one at a time (a comparison off by one, a constant, a sort order,
a branch removed), to see whether the tests notice. 45 were caught, 14 were not and 2 did not build. The 14 were closed with tests (a score of
exactly 70, a locus with two alleles whose first would make a type, an allele one base short of novel, places on other strands and contigs, an
allele inside another, the closest allele among more than 8, a blank line in the profile table, the limits 0 and 100, the order of a folder of
schemes), one was dead code (the contig names are never written) and was removed, and the 16 that were run again, the 2 that did not build
among them, are all caught but one: **a hit that begins exactly where the previous place ends is made part of it by `<=` and the tests do not
notice**, because two hits touch exactly only by chance and no genome here makes them.

## Not done, and not known

* **Not compared with mlst** on a real scheme and genome. The thresholds above are guesses, and mlst's choices at the edges (an allele that runs off a
  contig, a base `N`, a scheme tied with another) are not documented in what was read.
* An end base that does not match may be left out of the stretch by the free overhang, so an allele with a changed base at its very end can be called
  partial where mlst might call it novel. A base `N` in a contig counts as a mismatch.
* Only the first 8 candidates at a place are aligned. A genome with a locus several times over in different places is searched for each.
* BLAST is not used, so the hits differ from mlst's in detail. No compressed input, GenBank or EMBL input, `--legacy`, `--info`, `--list`, `--fofn`,
  threads, or a scheme download (that is yours to do).
