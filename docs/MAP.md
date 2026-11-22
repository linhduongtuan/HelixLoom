# Mapping reads and finding overlaps, after minimap2

`helixloom map` maps long reads to a reference, or finds the overlaps between reads, and writes PAF. It is the third stage of the plan of
`docs/ROADMAP.md` 2.7, written from the paper (Li, *Bioinformatics* 2018, 34:3094), the manual page and the README of minimap2 (which is MIT
licensed; none of its source was read or copied). It was not compared with minimap2: that needs a download, and nothing here claims the same
output.

```
helixloom map [options] TARGET READS        map the reads to the target, PAF out
helixloom map -x ava-ont [options] READS    overlaps between the reads, each pair once
```

`TARGET` and `READS` are FASTA or FASTQ files (not compressed; several files that are one input are joined with commas). A sequence is named by
the first word of its header.

## What it does

1. **Minimizers.** Every k-mer (k 15) is hashed as the smaller of itself and its reverse complement (a k-mer that is its own reverse
   complement is not taken), and of every window of w (10) k-mers the one with the smallest hash is kept, with its strand. Letters other than
   A, C, G and T break a sequence into runs. The hash is the finalizer of SplitMix64 (the paper does not name its own).
2. **Index.** All the minimizers of all targets, found by hash. The most frequent are ignored when a read is looked up: the share `-f` (0.0002)
   of the distinct minimizers, never fewer than the 10 occurrences below which nothing is.
3. **Chaining.** The places where a read and a target share a minimizer (anchors) are chained as the paper describes: a dynamic programme in
   the order of the target, an anchor gaining the bases it adds (at most k) and paying `0.01·k·|l| + 0.5·log2|l|` for a difference `l` between
   the gap on the read and the gap on the target, with the gap on either at most `-g` (10,000) and the difference at most `-r` (500), looking
   back until 50 anchors have not improved the score. Chains are taken from the best end back, each anchor in one chain; a chain has at least
   `-n` (3) anchors and `-m` (40) score.
4. **Primary and secondary.** From the best score down, a chain that overlaps a primary on the read by half (`-M`) of the shorter is its
   secondary. A secondary scoring under `-p` (0.8) of its primary, or beyond `-N` (5), is left out.
5. **Mapping quality**: the paper's `40 (1 − f2/f1) min(1, m/10) ln f1`, held between 0 and 60 (`f1` the score of the best chain, `f2` of the best
   secondary, `m` the anchors of the best).
6. **PAF.** Twelve columns, then `tp` (`P` or `S`), `cm` (anchors), `s1` (score), `s2` (the best secondary's score), `dv` (the divergence that the
   chain suggests, from the paper's `ln(n/m)/k`). Without `-c`, column 10 is the bases the anchors cover and 11 the longer of the two spans.
7. **`-c`: base-level alignment.** The anchors are exact matches; what lies between two of them is aligned by the program's own aligner (affine
   gaps), and written as `cg` (the CIGAR) and `NM`; columns 10 and 11 are then counted from it. The alignment runs from the first anchor to the
   last: **the ends are not extended** (minimap2 extends them with a drop-off rule), and a gap too large for one table is written as a run of
   matches against the shorter side with the rest a gap.
8. **Overlaps (`-x ava-ont`)**: `-k15 -w5 -m100 -r2000`; the reads are indexed and mapped to themselves; a read is not mapped to itself, a pair is
   reported once (the earlier read is the query), every chain is reported as a primary, and the quality is 255 (none).

## Options

`-x map-ont|ava-ont`, `-k`, `-w`, `-f`, `-g`, `-r`, `-n`, `-m`, `-p`, `-N`, `-M`, `-c`, `-o FILE`. A flag after `-x` changes what the preset sets.
The values are minimap2's defaults as its manual page gives them (`-r` as its first number); `ava-ont` is `-k15 -Xw5 -e0 -m100 -r2k`.

## What was measured

* The sketch equals a plain brute-force one for every k and w tried, with breaks, lower case and RNA; a sequence and its reverse complement give the
  same minimizers mirrored with the other strand; the chain scores are worked out by hand from the paper's formula; mapping quality is the formula.
* Reads simulated from a known sequence with 3 % each of substitutions, insertions and deletions (about 10 % error): all 1,500 reads of 1 to 10
  kb from a random million-base genome were mapped to where they came from, on the right strand, none unmapped and none wrong with a quality of
  30 or more; 8 Mb of reads mapped in 0.06 s on all cores of this Mac (a release build; `cargo test --release minimap_accuracy -- --ignored
  --nocapture`). **A random genome has no repeats and no low-complexity regions: this says the algorithm works, not that it is as accurate as
  minimap2 on a real genome.**
* A read inside a repeat of two nearly identical copies has a second mapping and a quality of 10 or less; a read from a copy that runs into the
  sequence beside it is certain; a read of two distant places has two primaries that do not overlap on the read; reads from no place map nowhere.
* With `-c` the CIGAR consumes exactly the spans of the read and the target and its matches equal those counted by walking it along the two
  sequences.
* The code was broken on purpose in 47 places (the strand of a minimizer, ties in a window, the gap cost, the limits of a gap and of the bandwidth, the
  score of a chain that runs into one taken before it, the overlap that makes a secondary, the quality, the PAF strand, the pair that is skipped in
  overlaps, and more): all 47 were caught, 8 of them only after tests were written for them.
* Overlaps: reads simulated to overlap in known places are found once each, with the overlap to within 250 bases, on the right strand, and a
  read that overlaps nobody has no line.

## Not done

The SIMD kernel and the drop-off extension of the ends; splice mode; short reads, in pairs or not (`sr`); homopolymer-compressed minimizers
(`-H`, for PacBio CLR); the other presets (`map-hifi`, `map-pb`, `asm5`, `asm10`, `asm20`, `lr:hq`); minimizers sampled from the frequent ones
(`-e`); the `cs`, `MD`, `AS`, `ms`, `nn` and `rl` tags; supplementary alignments (`SA`); SAM output; compressed input; an index kept on disk;
and any claim of the same output as minimap2.
