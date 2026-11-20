# K-mer analysis: `helixloom kat`

K-mer counting and the readings that follow from it, in the manner of KAT, the K-mer Analysis Toolkit (Mapleson et al., *Bioinformatics*
2017, 33:574). Written from the paper and the manual; KAT is GPL licensed and none of its source was read or copied (stage 2 of the plan in
`docs/ROADMAP.md`, 2.7).

```sh
helixloom kat hist -k 27 reads.fq                      # the spectrum, with the valley, the peak and the genome size it says
helixloom kat gcp -k 27 reads.fq                       # the spectrum split by the number of G and C: a matrix
helixloom kat comp -k 27 reads.fq assembly.fa          # how often each k-mer is in one against the other
helixloom kat sect -k 27 -c reads.fq assembly.fa       # how well the reads cover each sequence of the assembly
```

* **Input**: FASTA or FASTQ files, not compressed; files that are to be one input are joined with commas (`a.fq,b.fq`).
* **k-mers** are canonical (a k-mer and its reverse complement are one, whichever is smaller as a number), made of A, C, G and T (`U`
  counts as T); a window that holds any other letter is not counted. `-k` is 1 to 32 (default 27, the length the manual's examples use:
  the default of KAT itself was not in what was read).
* **`hist`** writes `#` lines (k, distinct and total k-mers, and what the spectrum says), a line `###`, and `multiplicity<TAB>count`
  for each bin up to the last that is used. `--high N` makes the last bin hold all that occur N times or more (default 10,000); `--low N`
  leaves out those that occur fewer than N times.
* **`gcp`** is the same as a matrix: a row for each number of G and C (0 to k), a column for each multiplicity.
* **`comp`** is a matrix of the multiplicity in the first input (rows) against the second (columns), with row 0 and column 0 for what only
  one of them has; the `#` lines say how many distinct k-mers are in both and in each alone.
* **`sect`** takes the reads with `-c` and a file of sequences, and writes, for each sequence: the median and mean multiplicity of its
  k-mers, its G+C, the number of k-mers, of those that hold a letter that is no base, and of those that occur at all. `--coverage FILE`
  also writes the multiplicity of each k-mer along each sequence (`-` for one with a letter that is no base). KAT's own table has a last
  column (`%non_zero_corrected`) whose meaning was not in what was read: it is not written.

## Reading the spectrum

A spectrum of reads has a mound of errors at low multiplicity (a k-mer with a mistake in it occurs about once), a valley, and a peak at
the coverage of the genome's k-mers (a genome with two different copies has a second peak at half of that), with repeats above. The plain
rules here (`seq::kmer::summarize`): smooth the counts over three bins; the **valley** is where they stop falling after multiplicity 2;
the **peak** is the highest bin from there on; the **genome size** is the sum of multiplicity times count from the valley on, over the
peak (the textbook estimate; it counts the last bin at its own multiplicity, which underestimates a genome with long repeats); a **half
peak** is a local maximum within one of half the peak that is at least a tenth as high. A spectrum that never rises again has no peak and
no size (too few reads). **This is not KAT's own fit** (`kat_distanalysis.py` fits distributions to the peaks; its method is not in the
manual pages that were read), so a size here and one from KAT may differ. Where the first rise is only a small bump on a long error tail,
the valley may be put after it.

## What was checked, and what was not

Checked: every count against a plain string-based counter (all the windows of random sequences with other letters, lower case and U, for
k = 1, 2, 3, 5, 11, 21, 27, 31, 32); a case worked out by hand; the matrices; reads simulated from a random genome of 20,000 bases
(30-fold, both strands, one base in 200 wrong): the peak came where the k-mer coverage is, and the size within 10 % of the genome; and
synthetic spectra with one peak and two. **Not compared with KAT** (its program and its example data were not run: that needs a download),
so no claim of the same numbers is made. Counting is in memory: a bacterial genome's reads, not a human genome's. **Not done:** `kat
filter`, `blob`, `cvg` and the plots; compressed input; counting on disk.
