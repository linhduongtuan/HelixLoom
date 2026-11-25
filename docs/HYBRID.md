# Joining an assembly graph with long reads, after Unicycler

`helixloom hybrid` takes the assembly graph that short reads give (a GFA, from `helixloom dbg` or from SPAdes) and long reads, and joins the
segments that the long reads show to follow each other, through the repeats that the short reads could not span, into longer replicons, circles
among them. It is the sixth stage of the plan of `docs/ROADMAP.md` 2.7. It was written from the paper of Unicycler (Wick, Judd, Gorrie and Holt,
*PLoS Comput Biol* 2017, 13:e1005595) alone: Unicycler is GPL licensed, so none of its source was read or copied, and its README could not be
fetched. It was not compared with Unicycler and was never run on real reads.

```
helixloom hybrid [options] GRAPH.gfa LONGREADS...          replicons as FASTA (standard output), a summary on the error output
helixloom hybrid -o DIR GRAPH.gfa LONGREADS...             DIR/assembly.fasta and DIR/assembly.gfa
```

`helixloom unicycler` is another name for the same command. The options are this program's. `helixloom gfa inspect GRAPH.gfa` reports on the same graphs.

## What it does

1. **Reads the graph**: the segments with their bases and depth (`DP`, or `KC` over the k-mers of the segment, or 1) and the links, which must all
   overlap by the same k−1 bases and agree in their bases (a de Bruijn graph does). A graph with no link takes `--overlap`.
2. **Gives each segment a number of copies** (the paper's *multiplicity*). A single copy has the depth of the genome as a whole, taken as the depth
   that half of the graph's bases are at or below (`--single-depth` to give it). *Seeds*: a segment within a quarter of that depth (`--depth-tolerance`)
   with at most one way on at each end is one copy. Then a greedy spreading: where the ways on from an end are each reached from that end only, the
   copies that leave are the copies that arrive, so a segment whose neighbours at an end are all known is their sum, and one neighbour that is not
   known is what is left; a number is taken only if the depth agrees (within 0.4 of it). When nothing more follows, the largest segment that could
   be one copy (within half of the depth) is set to one, and so on. What is left gets its depth rounded, and the report counts these.
3. **Finds the joins that the long reads show.** Each read is mapped (the program's mapper, `helixloom map`) to the segments that are single copy.
   Two alignments of at least 300 bases (`--min-align`) that follow each other along a read, the first stopping within 200 bases (`--end-slack`)
   of the end of its segment and the second starting within 200 bases of the start of its, are a join: the read leaves the first segment, crosses a gap,
   and enters the second, each in the direction its alignment goes. A join and the same join read on the other strand are one; the reads that make
   one join are one candidate, with the gap each shows.
4. **Makes each candidate a bridge**: the routes through the repeat segments of the graph (up to 8 of them, none that is single copy) from the end of the
   first to the start of the second are looked for; the one whose length is nearest the reads' median gap, within 100 bases and a fifth of the gap, is the
   bridge, and it is the graph's own bases, as the paper says, which have a short read's accuracy. A route of none is a link that is already there. With no
   such route, the bridge is the bases of the one read whose gap is the median (no consensus of reads is made), which the paper says is more error-prone,
   and is rated lower; it needs a gap of more than nothing.
5. **Scores it** from 0 to 100: the product of a factor each for the number of reads (n/(n+1)), how well the reads' gaps and the route agree, the lengths
   of the two segments (each to 1 kb), how much of each alignment the minimizers cover, and how alike the two depths are, halved for a bridge of the
   reads' own bases. **The formula is this program's** (the paper lists these as what quality depends on and gives no formula). A bridge must reach
   25 (`--mode conservative`), 10 (normal, the default) or 1 (bold), the paper's numbers on this other scale, or `--min-quality`; `--min-reads` asks for
   more reads.
6. **Puts the bridges in, best first.** A bridge is a segment of its own between the two (the route's bases, or the reads'), and every other link at the
   two ends is taken away: a single copy has one way on at each end, and each end takes one bridge. A repeat that bridges have used as many times as it
   has copies is removed.
7. **Joins what is unbranched.** Bold joins every run that has no branch; conservative and normal join only the single-copy segments and the bridges. A
   segment whose end leads to its own start, and to nothing else, is a circle, written once round and marked `circular=true`.
8. **Rotates the circles** (with `--start-genes`): a FASTA of nucleotide genes that you give, such as *dnaA* or *repA* of the species (none is bundled), is
   looked for in each circle with the mapper, and the circle is turned and, if the gene lies on the other strand, flipped, so that the gene is first on the
   forward strand. The paper does this with TBLASTN, which is not done: a distant gene is not found.

The FASTA names are `>NUMBER length=BASES depth=DEPTHx`, with `circular=true` for a circle (the form is this program's choice; the README that shows Unicycler's
could not be read). The numbers are by length, longest first; the depth is the segment's average k-mer depth. The GFA is `helixloom dbg`'s.

## What is from the paper, and what is this program's choice

Only part of the paper could be read as text: of about 106,000 characters, the first 40,000 (the introduction, the design: short-read assembly,
multiplicity, bridges, semi-global alignment, bridging with long reads, bridge application, the final steps, the modes, and the start of the results);
the rest (limitations, summary, the legends of the supplement) was read through a summary. The README and the supplementary figures were not read.

* **From the paper**: that a segment of one copy has the depth of the genome and one connection at each end; the seeds, the greedy spreading and the
  largest segment set to one when it stalls; that long reads are aligned to the single-copy segments; that a bridge's bases are a route through the
  graph that the reads' gap picks, and the reads' own bases only when there is none; what quality depends on; that bridges are applied from the best down and
  a repeat's copies used are taken off; the three modes, with the thresholds 25, 10, 1 and the differences in what is joined; that a circle is one segment
  whose end leads to its start; the start gene (*dnaA* or *repA*) and the turning and flipping.
* **Second-hand** (the summary of the unread part): that the paper scales quality to 0 to 100 as a product of functions.
* **This program's choice**: every number (the tolerances 0.25, 0.4 and 0.5 of the copy numbers, the 8 repeats of a route, the 100 bases and fifth of a gap,
  the 300 and 200 bases of the alignments, the 1 kb and 0.4 of the score); the mapper in place of the paper's semi-global alignment; the median read in
  place of a consensus; the formula; that conservative and normal join alike; the nucleotide search for the start gene; the form of the FASTA names.
* **Not done**: the bridges that SPAdes' own paths give from read pairs (this program has no pairs); the branch-and-bound search of the graph; the removal
  of contigs under half the median depth; polishing (Pilon, Racon and Medaka stay external); TBLASTN.

## What was measured

On sequences of random bases, graphs made by `helixloom dbg` from error-free short reads, and long reads simulated with 6 % errors:

* a circle `X R Y R Z` whose repeat `R` (1,200 bases) is twice in it: the short-read graph has three segments and two ways through the repeat; the long reads
  make two bridges and the circle is one circular segment equal to the genome up to where it is cut and the strand, in each of the three modes; with
  `--min-quality 100` no bridge is made and the graph is as it was; reads of another genome change nothing;
* two segments with 300 bases between them that the graph does not have: joined with the reads' bases, the ends the graph's own, the length within 150 bases;
  a segment lacking 200 bases of its circle: closed into a circle; two segments that follow each other with no gap and no link: not joined;
* the bridges, one against another (the better wins an end), a bridge that is a link (the other links go), the depth of a bridge, bold against normal;
  a gene at either strand, across the join, with one base in twenty changed, the best of several, and none;
* a throwaway run, not kept in the suite: a circular chromosome of 50.5 kb with three copies of a 1.5 kb repeat, a circular plasmid of 8.6 kb with two copies of
  an 0.8 kb one, short reads at 40 times with 0.2 % errors through `helixloom dbg` (eight segments, none circular), and 229 long reads of 2 to 20 kb at 30 times
  with 8 % errors (3 % changes, 2.5 % put in, 2.5 % taken out): five bridges, and two circular replicons, each of the length of its sequence and equal to it up
  to where the circle is cut and the strand. A single run on random sequences, the easiest case; real repeats, real plasmids at other copy numbers and real
  reads are not tested.

**Deliberate breakages.** 87 one-line changes were made to the code on purpose, one at a time (a bound, a sign, a strand, a sort order, a branch removed), to see
whether the tests notice. 60 were caught and 27 were not. Of the 27, fifteen were closed with tests (a link with no overlap field, a segment of exactly k
bases, the default tolerance, a segment that more than one thing leads into, the rounding of a depth, a secondary alignment that is not a join, a bridge of exactly
k bases, the better of two bridges for one end, the fewest reads, how the reads' gaps and the route's distance are rated, a gap that is a read's and not the
graph's, the default share of a gene, a line that holds the start gene, a quality of exactly 100); one was a guard that could never be reached (a bridge
between a segment that had been removed, which nothing in the loop removes), and was deleted; and **eleven are left, which change nothing that was found**:
the boundaries where an alignment is exactly as long as `--min-align`, stops exactly `--end-slack` from the end, or has a share of the gene exactly at the minimum,
or a bridge's quality is exactly the threshold (an alignment's length cannot be set to the base); a copy count that takes a neighbour that has none left, and
one that sets a value of one (the largest-segment rule gets there as well); the tie that picks the smaller number among segments of one length; a gap of
exactly nothing (no bases either way); a repeat of one copy being "used up" (a single-copy segment is never on a route); and the check that an alignment
is the primary one when the gene is found (the same gene twice in the doubled circle gives the same start). The 15 that were closed were run again and all
caught.

## Not done, and not known

* **Never run on real reads, never compared with Unicycler** (or any assembler).
* **Plasmids that do not have the chromosome's copy number**: their depth does not fit the single copy, so their segments come out as repeats and are not
  bridged (the paper says the same of its method); `--single-depth` helps only for one. A whole plasmid in one piece is not changed.
* A gap of nothing or less that has no link, no route and no overlap in the graph is not bridged. A read that runs across three segments makes two joins,
  not a path of three at once. Reads that are chimeras make joins that are wrong, and only the quality (the number of reads) stands against them.
* The long reads are not corrected or polished; a bridge from the reads' own bases has their errors (about the rate of the reads). No Pilon or Racon step.
* The graph must be a de Bruijn graph (links of k−1 bases that agree); one from an assembler with other overlaps is refused.
* Windows and Linux builds not tried.
