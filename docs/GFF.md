# GFF3 and GTF: `helixloom gff` and importing genes onto a sequence

What `gffread` does for an annotation, built into HelixLoom (stage 1 of the plan in `docs/ROADMAP.md`, 2.7). Written from the GFF3
specification (Sequence Ontology 1.26), the GTF pages of Ensembl and BioPerl, and `gffread`'s usage text and paper. `gffread` is MIT
licensed; none of its source was copied.

## On the command line

```sh
helixloom gff annotation.gff3 -g genome.fa -w transcripts.fa -x cds.fa -y proteins.fa   # the sequences
helixloom gff annotation.gtf -C -U -o filtered.gff3                                     # a filtered annotation
helixloom gff annotation.gff3 --table @id,@chr,@start,@end,@strand,@numexons,@cdslen    # one line for each transcript
```

* **Reads** GFF3 and GTF (told from the `##gff-version 3` line or the attributes). A transcript is what exons and CDS name as their
  `Parent` (GFF3) or `transcript_id` (GTF); a CDS with no parent (Prokka, Bakta and Prodigal write them so) is a transcript of its own; a
  gene with CDS directly under it stands for its transcript. When exons are given, UTR lines add nothing; when CDS is given, start and
  stop codon lines add nothing; with CDS and UTR only, the exons are made from them.
* **Writes** GFF3 (a line for each transcript, then its exons and CDS with `Parent`) or, with `-T`, GTF. `-F` keeps all the attributes.
* **Cuts sequences** with `-g` (FASTA, found by the first word of each header): `-w` the spliced exons, `-x` the coding sequence (the
  phase of its first segment taken off), `-y` the protein (standard code; the stop codon at the end is left out, one inside is `.`, or `*`
  with `-S`). On the minus strand the transcript is the reverse complement. `-W` puts `loc:` and `exons:` or `segs:` in the names.
* **Filters** (as in `gffread`): `-C`, `--nc`, `-U`, `-i N`, `-l N`, `-r [strand]seq:start-end` with `-R`, `-N` (introns that are not GT-AG,
  GC-AG or AT-AC; needs `-g`), `-V` (a stop inside the CDS; needs `-g`), `-J` (a CDS that begins ATG, ends with a stop, and has none inside),
  `--ids` and `--nids`.

## On a sequence that is open

**Features ▸ Import Features from GFF/GTF…**, a `.gff`, `.gff3` or `.gtf` file opened or dropped on the window, or the command palette, puts
the genes on the open sequence as features: a `CDS` for each coding transcript (with its segments, on its strand), an `mRNA` too when the
exons say more than the CDS, and a feature of its own kind (`tRNA`, `rRNA`, `ncRNA`) for the rest, with the gene's name and the product.
The annotation's name for the sequence must be the sequence's name (or its first word), or the annotation must name only one sequence.
Features that reach beyond the sequence are left out and counted. Undo takes them off.

## What was checked, and what was not

Checked against cases worked out by hand (a gene of two exons on each strand, phases 1 and 2, a stop inside, GFF3 and GTF of the same gene,
what the writers write read back, every filter, every error with its line): see *Help ▸ How Each Analysis Was Checked*. **Not compared with
`gffread` itself**: its repository's `examples/` folder (an annotation, a genome and the outputs its test script expects) would be the
test, and fetching it is a download that has not been asked for. Rules taken on trust from `gffread`'s documentation and not verified:
that the stop codon is left out of the protein (inferred from its test output), that phase trims the first segment of the CDS, and the
default `.` for a stop inside.

**Not done:** merging transcripts into loci (`-M -K -Q -Y`), BED and TLF, streaming, `-L` (Ensembl GTF), `-u -j -O -t -Z`, a genetic code
other than the standard, and a window for it (the command line and the import are what there is).
