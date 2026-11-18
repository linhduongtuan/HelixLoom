# Motifs in DNA (the Motif Finder and `helixloom motif`)

A *motif* is a short pattern that many sequences share, written as a position weight matrix: for each column, how likely each
base is. This program finds motifs by the method of MEME in its simplest form, scans sequences for more places that match one
(with exact p-values, the way of FIMO), and compares motifs. It is written from the papers, not from the programs' code
(`CREDITS.md`).

## Using it

* **In the program**: Tools → *Motif Finder*. Paste DNA sequences (FASTA, or one on each line) or press *Use an example*, set the
  width, the number of motifs, the strands and the effort, and press *Find motifs*. The window shows each motif as a logo with its
  consensus, information and sites; *Scan* lists more windows that match the chosen motif with a p-value at most the one you set;
  *Copy as MEME text* gives the motif to other programs.
* **Where it runs is your choice**, as for every heavy job: between the settings and the button, *This computer*, *Workstation*
  and *Automatic* (a set of more than 300,000 bases goes to the workstation, a smaller one stays here). This computer takes at
  most 300,000 bases (60,000 in a browser) and 2,000 sequences; a job on the workstation takes up to 100,000 sequences and as
  many bases as fit in 4 GB of memory for their windows (about 16 million for a motif 8 columns wide on both strands). The
  workstation's answer is shown in the window if it still holds the sequences the answer is about; a result opened later from the
  Workstation page is not shown for other sequences (the sequences are not kept with it).
* **From a script or a cluster**: `helixloom motif -i sequences.fa -o motifs.meme --sites places.tsv`:

  | Option | Meaning |
  |--------|---------|
  | `-i`, `--input FILE` | the sequences: FASTA, or one on each line (a last argument that is no option is the input too) |
  | `-o`, `--output FILE` | where to write the motifs, as one MEME file (the standard output when left out) |
  | `--sites FILE` | also write the places found: motif, sequence, position (from 1), strand, score in bits, probability |
  | `-w`, `--width N` | columns of a motif, 3 to 40 (default 8) |
  | `-n`, `--motifs N` | how many motifs, one after the other: the places of one are hidden from the next search (default 1) |
  | `--one-strand` | look on the strand as written only (default: both) |
  | `--seeds N`, `--iterations N` | starting words followed to the end (12), rounds from each start (60) |
  | `--effort X` | how much work the search may spend over a big set; 1 is the window's budget of about five seconds; 0.01 to 1000 |

  The log goes to the standard error and the exit code is 0 when it worked.
* **On the workstation**: the job kind is `motif`. A request holds the sequences (one to a line) and the settings (width, motifs,
  strands, seeds, iterations, effort); the answer holds each motif's matrix, gain and sites. The app checks that the answer is
  about the sequences it sent (every site is a window inside a sequence) before it shows it.

## What it does

1. **Discovery**: every window of every sequence (both strands unless told otherwise) is a candidate site. The words held by the
   most sequences start the search; each start is followed by expectation-maximisation, with at most one site in a sequence
   (Bailey and Elkan 1994). The motif with the largest *gain* (how much likelier the sequences are with it than with the
   background alone, in bits) is kept, its sites (the best window of each sequence where it is more likely than not a site) are
   hidden from the next search, and the search runs again for the next motif.
2. **The work is bounded**: the search can go over about 5 × 10⁸ window-columns in all, so a big set is searched less thoroughly
   (fewer starting words, fewer seeds and rounds) rather than for ever; a small set is searched in full. `--effort` and the
   *Effort* setting scale that budget, which is the one way to search a big set more thoroughly.
3. **The scan** scores a window against the motif as log-odds in bits (against the sequences' base frequencies) and gives the
   p-value of the score exactly, from the distribution of the scores of all random windows, rounded to hundredths of a bit.
4. **Comparison** of two motifs is the best mean Pearson correlation of their columns over the ways they can overlap, either strand.

## Limits

* DNA only (A, C, G, T; U is read as T). Protein motifs are not built.
* At most one site in a sequence and one motif width for a search: this is MEME's simplest model (its "zero or one occurrence per
  sequence" with a fixed width), not its "any number of repeats" or its width choice. STREME-style discovery of many short
  motifs is not built either.
* **Checked on planted motifs only**: a word with one base in ten changed, put in random sequences on either strand, is found
  where it was put and as the consensus; two motifs are found one after the other; the exact p-value against scoring all 256
  words of a four-column motif; the command and the job against the same search done in the program. **Not compared with MEME**
  on real promoter sets, so there is no claim to MEME's accuracy.
* Known motifs to compare with (JASPAR is free to use but is a download, not made without your yes) and scanning on the GPUs are
  not built.
