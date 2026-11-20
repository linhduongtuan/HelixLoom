# Aligning DNA by its translation

For coding DNA, aligning the letters can put a gap inside a codon or match a codon's third base with another's first, which no protein would
do. On the **Alignment** page, tick **Align by the translation (keep codons whole)**, choose the reading frame (1, 2 or 3) and the genetic code,
and press **Align**:

1. Every sequence is read in the reading frame: the bases before the first codon (0, 1 or 2) are a *head*, the whole codons are translated, and
   the 0 to 2 bases left at the end are a *tail*. Gaps and spaces in what you gave are not part of a sequence; `U` is read as `T`.
2. The proteins are aligned: by MAFFT when it is chosen (and found where the job runs), else by the built-in progressive method (with its
   settings on the page, with BLOSUM62 scores). The built-in centre-star method has no scores for amino acids, so choosing it means the
   progressive method with its usual settings for the proteins. On a workstation, the proteins are what is sent.
3. The alignment is carried back: every residue becomes its codon, every gap becomes three gaps; the heads come first, and the tails, padded
   with gaps, last. So the columns of the DNA alignment come in threes, and a column of three holds the codons of one amino acid.

A stop codon, or a codon with an ambiguous base, is written `X` in the protein that is aligned (the codon itself is kept and comes back).

The result is a DNA alignment like any other: the conservation, the consensus, a tree, and **Pick primers and probes…** work on it. The record
(*How this was made*) says the frame and the code, and **Run again** reads them from it. The choice is kept with the page.

## What it refuses

A sequence that does not look like DNA (the answer names it), one with no whole codon in the frame, a choice of frame or code that does not
exist, and an aligner that returns a changed or shortened sequence, are all refused with the reason, and nothing is shown in place of the
alignment that was there.

## Limits

* One reading frame for all the sequences. Sequences that are in different frames, or have frame shifts, stop codons in the middle, or
  introns, are not repaired: translate them in the right frame first, or align the letters.
* The frame is not found for you.
* The protein alignment is not shown (the DNA alignment is what the page keeps).
* Not compared with pal2nal, TranslatorX or MACSE on the same sequences: none was run. The cases are worked out by hand, and for random
  sequences with random gaps it is checked that no base is lost or moved and that every gap is a whole codon (`Help ▸ Checks`).
