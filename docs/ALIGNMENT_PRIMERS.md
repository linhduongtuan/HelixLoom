# Primers and probes from an alignment

On the **Alignment** page, *Pick primers and probes…* (also in the command palette) opens a window that designs primer pairs, and a TaqMan
probe between the primers if asked, from the columns that the aligned DNA sequences share, so that one order serves a family of sequences
(strains of a virus, copies of a gene). The method is this program's arrangement of parts it already had (primer3's rules in `seq::primer3`,
the probe rules in `seq::probe`); it is not a port of another program.

## What it does

1. Every column has a commonest base. A column is **usable** when at least the chosen share of all the sequences (a gap counts against)
   carry that base and the gaps are not the majority. With 100 % every sequence carries it.
2. The consensus of the usable columns, with `N` everywhere else, is the template. primer3's rules refuse any oligo that is not made of A, C,
   G and T, so no primer or probe can lie on an unusable column; the amplicon between the primers may.
3. Pairs (and probes) are picked from that template with primer3's rules. Positions are alignment columns, and the product size limits are
   asked in columns; the table gives the true sizes, row by row.
4. Each oligo is laid against every sequence. It **fits exactly** where it matches in every column, and **fits** where it has at most the
   mismatches allowed (a gap counts as one), and, for a primer, none in its last five bases, where a mismatch stops the polymerase most. A
   design **covers** the sequences that all its oligos fit. Designs are ranked by the sequences they cover, then by penalty.

## The window

* Size, product size, melting temperature and number of pairs, as in the Primer3 window; *Also pick a TaqMan probe* (and *Product 70–150*).
* **Columns must agree in at least … % of the sequences**; **mismatches** allowed; **only columns … to …** (a region).
* The table: left and right primer (5′ to 3′ as ordered) with their Tm, the probe, the product (fewest to most bases over the sequences
  covered), how many sequences the design covers (hover: how many match perfectly), and the penalty.
* **Show** marks the oligos on the alignment (a bar under the header and an outline over the rows, in the colours of the key) and scrolls to
  them. **Consensus** opens the consensus, without its gaps, as a sequence with the primers and the probe as features.
* Under the table, for the chosen design: the columns of each oligo, where it lies in the first sequence's own numbering, and the sequences
  it does not cover by name. **Copy table** and **Save table…** give one line to an oligo (tab-separated) for an order sheet or a
  spreadsheet.

## Limits

* The oligo is the consensus: no degenerate (IUPAC) oligos, which would cover more sequences with one order.
* The melting temperature and the penalties are primer3's for the consensus oligo, not for the oligo with the mismatches of one sequence. How
  many sequences an oligo fits is a count of mismatches, not a prediction that it primes.
* Nothing is checked against a genome for other places the oligos could bind: that needs one.
* An alignment of more than 50,000 columns is narrowed to a region first; the search runs while the window waits, as the Primer3 window's
  does.
* Not compared with another program that designs primers from an alignment. The cases are worked out by hand (`Help ▸ Checks`).
