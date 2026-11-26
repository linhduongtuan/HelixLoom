# Pipelines as files, and connectors to databases

Two things that turn the tools into a platform: a **pipeline** is a short text file of typed steps that can be run on a whole folder
and that does not repeat work it has done, and a **connector** fetches a record from a public database by its accession and
remembers where it came from. Both were built in November 2026 and tried on a Mac only.

## A pipeline file

```text
helixloom-pipeline 1
name: Gene trees, then a species tree
input: alignment
step tree engine=builtin bootstrap=0
step species
```

* The first line is the format and its version (a file of a later version is refused, with what to do).
* `input:` is what each file of the folder holds: `sequences` (FASTA, not aligned), `alignment` (FASTA, rows of one length) or
  `trees` (Newick), or `gfa` (an uncompressed GFA assembly graph), or `reads` (FASTQ, names ending in `.fastq` or `.fq`).
* `step NAME key=value …` are the steps, in order. A line that starts with `#` is a comment (a `#` after something is not).

| Step | Takes | Gives | Settings |
|------|-------|-------|----------|
| `align` | sequences | alignment | `engine=builtin` (centre-star) or `engine=mafft` |
| `tree` | alignment | tree | `engine=builtin` (likelihood search) or `engine=iqtree`; `bootstrap=N` (only with `iqtree`) |
| `species` | trees | tree | none (the ASTRAL method: `helixloom astral --help`) |
| `gfa-inspect` | gfa | graph-summary | none. Produces the same six metrics as `helixloom gfa inspect`. |
| `trim`, `filter`, `assemble`, `hybrid`, `flye`, `concatemer`, `circularize`, `rotate`, `medaka` | reads, gfa or sequences | reads, gfa or sequences | the assembly routes: see `docs/ROUTES.md`. |

Each step must take what the one before gives; the file is refused otherwise, saying which step. **A tree made for each file is
gathered into the list of trees that `species` takes**, so the steps before the first one that takes trees are done for each file
of the folder (on all the cores at once) and the rest once. `engine=mafft` and `engine=iqtree` mean *that program*: if it is not
installed the step fails and says so, rather than quietly using the built-in engine and giving other answers. `bootstrap` with the
built-in search is refused for the same reason, and with `iqtree` it must be 0 or at least 1000 (IQ-TREE's ultrafast bootstrap makes at
least 1000, so a smaller number would not be what is run).

### The cache, and resuming

With a cache folder, the answer of every step is kept under a key made from the step, its settings, the version of the installed
program that made it (when there is one), and the input. A step whose key is there is **not run again**. That is also how a run
that stopped is resumed: run it again, and what was done is found, and only the rest is done. Change one locus, and only its steps
and the gathering step run again. Answers are written beside, then moved into place, and each carries a checksum of its body, so a
crash leaves no half file, and a file that is cut short or changed counts as not there. The key holds the version of this program
as well as the step and the input: **a new release does not read the answers of an earlier one** (and a test fails when a built-in step
begins to answer differently between releases, to make the engine number in the code be raised). `clear-cache` removes only files that
this program wrote.

A file whose chain fails (too few sequences for a tree; a fault in a step is caught and counted as a failure too) is **left out of what is
gathered and named in the report**; the others go on, and the species tree is made from the rest. So is a locus whose gene tree has
**fewer than four leaves**, which says nothing about how species are related (the report says so in a line of its own). The run says
"N failed; left out: …" and the command exits with 1 when a step failed. Look at that line: a species tree from three loci of ten
is still a species tree. A folder is read before anything runs: a file of the wrong kind (an alignment whose rows differ in length), or
two files that would have the same name (`a.fa` and `a.fasta`), stop the run at once, naming the file. Files whose names do not end in
one of the accepted endings are not read, and not mentioned.

### Using it

```bash
helixloom pipeline example > loci.txt          # a pipeline file to start from (`example unaligned` aligns first)
helixloom pipeline check loci.txt              # what it does, and that its steps fit
helixloom pipeline run loci.txt ./loci --cache ./cache --out species.tre --out-dir ./gene-trees
helixloom pipeline clear-cache ./cache
helixloom pipeline export loci.txt ./loci --out-dir ./workflow   # a Snakefile and a main.nf for the same pipeline (below)
```

In the window: **Pipelines ▸ Pipeline file** (next to *Fixed route*, at the right of the page's header): open or save a file, start from a built-in example, choose the
folder and the cache folder, and run. The page lists what each step did (run now, or from the cache), and the species tree can be
copied or drawn in Trees. **The run that the page shows is kept with your session** (18 November 2026), in either mode of the page, with its
record: the next start offers it again with the open sequences (`docs/WEB.md`, *What the page keeps*). The folder and the cache folder are not
kept with it: choose them again to run again.

### Evaluating a folder of assembly graphs

```sh
helixloom pipeline example gfa > inspect-graphs.txt
helixloom pipeline run inspect-graphs.txt ./graphs --cache ./graph-cache --out-dir ./summaries
```

The pipeline takes `input: gfa` and runs `step gfa-inspect` once per file. For example, `graphs/sample.1.gfa` produces
`summaries/sample.1.tsv`. Each TSV gives the segment, link, component and isolated-segment counts, the sum of known segment lengths,
and the number of segments of unknown length. A malformed graph record is reported under its input file, while valid files still produce
summaries. The command returns a failing exit status if any graph failed. As with other pipelines, an unreadable file or invalid UTF-8 stops
folder loading. When you reuse an output folder, an older TSV for an input that now fails or has been removed can still be on disk.
To check a single input and choose its output filename, pass that input file and `--out summary.tsv`.
On a rerun the cache reuses summaries for unchanged inputs.
The input folder's graph files are read into memory for the batch. [Supported GFA records and metric definitions](DEBRUIJN.md#inspecting-a-graph)

On the desktop, choose the assembly graph preset under **Pipelines ▸ Pipeline file**, choose the input folder, and press **Run pipeline**.
The **Graph summaries** table has one row per successful file. **Copy TSV** copies the six metrics for that file.

### Exporting to Snakemake and Nextflow

`helixloom pipeline export PIPELINE INPUT [--out-dir DIR] [--mafft NAME] [--iqtree NAME] [--helixloom NAME]` writes the pipeline as a
**`Snakefile`** and a **`main.nf`** (Nextflow, DSL 2) into `DIR` (default `workflow`), which do for every file of the folder `INPUT` what
`pipeline run` does: `mafft --auto --quiet`, `iqtree3 -s … -m MFP [-B n] -seed 1 -redo -quiet --prefix …` (the flags the engine itself uses), the trees put
together, and `helixloom astral -i … -m 4` (a gene tree of fewer than four leaves is left out, as the engine leaves it out). The programs are
called by the names given (`iqtree2` for an older IQ-TREE, a path for `helixloom`); a name may hold letters, digits and `/ . _ + -` only.

* **Export requires a command for every step.** The graph pipeline uses `helixloom gfa inspect` and publishes a TSV for each input
  in `results/summaries/`. An alignment or tree step with `engine=builtin` (the default) is done by this program's own code,
  which no other workflow engine can run, so the export is refused and says which step and what to change (`engine=mafft`,
  `engine=iqtree`) or to use `pipeline run`.
* The workflow finds its files by **one ending** (`*.fasta`): a folder that has `.fa` and `.fasta` files is refused. The folder is written into the
  files as an **absolute path** that may hold letters, digits and `/ . _ + - ~ @ , : =` only (no space, quote, brace or star), so the files
  are tied to this computer's folder; edit the one line (`glob_wildcards(...)`, `params.input`) to move them.
* A pipeline that ends in trees (no `species` step) leaves one tree for each file; one whose input is `trees` gets a species tree from each file.
* Sequence names must be plain (letters, digits, `- _ .`): IQ-TREE refuses others, and the engine here, which renames sequences itself, does not
  have that limit.
* **Tried:** the text of both files is compared with the text expected, for the pipelines above and the errors; the `main.nf` of a pipeline
  (align with MAFFT, a tree with IQ-TREE 3.1.4, `species`) was read by Nextflow (`-preview`) and **run for real** on three loci of five
  sequences with the MAFFT and IQ-TREE installed here and this program's `astral`: 7 processes ran and the species tree it wrote is the very tree
  that `helixloom pipeline run` gives for the same folder. **Snakemake is not installed on the computer this was built on, so the Snakefile has not been run or
  parsed by it**; it follows Snakemake's documented `glob_wildcards` and `expand` idioms.

### On the workstation

With a workstation connected (the Workstation page) a pipeline file goes there like the other jobs: the choice *This computer / Workstation /
Automatic* is drawn on the page, and with *Automatic* a run is sent when its folder holds more than 100 files, or more than a million bytes, or when a step asks for
more than 200 bootstrap replicates (the same limits as the other jobs; the bytes of the files stand for the letters). The files are read on this computer and sent with the
pipeline in one request (at most 20,000 files, 128 MiB of files, 5,000 sequences in a file, and the server's `--max-upload-mb`); the report that comes back is the
same report a run here makes, and opens on the same page, saying it ran on the workstation. A job is in the Workstation page's list like any other, and its result can be
fetched from there later.

**The answers are kept on the workstation, in its own folder** (`flow-cache`, beside `jobs` in the server's `--data-dir`): a request names no folder, so the cache folder chosen on this
page is not used for a run there. A run sent again with one locus changed does only that locus's steps and the gathering; the cache is not emptied by the server's
`--keep-days`, which removes results only, and has no size limit: on the workstation `helixloom pipeline clear-cache DATA-DIR/flow-cache` forgets it all. If the folder
cannot be written the run still goes, every step is done, and the report says so. `engine=mafft` and `engine=iqtree` mean the programs installed on the
workstation, as they mean those installed here for a run here.

**Tried** on a Mac only, with the server and the app on the same computer (the app's own tests, and `tests/remote.rs` with the real program, in which a second job finds all five
answers the first kept). Not tried: a workstation across a network, a Linux or Windows server, a run of thousands of files, or IQ-TREE or MAFFT started by a job of the server.

## Connectors: File ▸ Fetch from a Database…

| Database | Accession looks like | What comes | Goes to |
|----------|---------------------|------------|---------|
| NCBI nucleotide (E-utilities `efetch`) | `NC_000913.3`, `U00096.3`, `NZ_CP009072.1` | GenBank file | a new tab in Sequences |
| NCBI protein | `NP_414542.1` | FASTA | the list of the Alignment page |
| UniProt (REST) | `P69905`, `A0A024R161`, an isoform `P69905-2` | FASTA | the list of the Alignment page |
| PDB (RCSB) | `1A3N` | PDB file | the Structure page |
| AlphaFold DB | a UniProt accession | PDB file (found through the entry's `pdbUrl`, which must be a file of `alphafold.ebi.ac.uk`) | the Structure page |

* The accession is checked against the form of the database's own accessions before it goes into an address, so nothing typed or
  pasted can add to a request.
* The web is reached with **`curl`** (the same way the structure prediction is), restricted to `https` (also for redirects, of which
  there are at most five), without the person's `.curlrc`, one minute at most, and a limit of 200 MB that curl applies only when the
  server says the size beforehand. NCBI requests are spaced 350 ms apart (its limit without a key is three a second). Records are kept in the `downloads`
  folder of the settings folder: asking again for a kept record asks the database nothing (*Get a fresh copy* does).
* **Where it came from is kept**: beside each kept file a `.provenance` file says the database, the accession, the address, the day,
  the licence and a checksum (a file changed on the disk is fetched again), and a sequence opened from a record carries the same line
  in its comment, which a GenBank file saved from it keeps. A record that cannot be kept (the folder cannot be written) is still opened,
  with a line saying it was not kept. A PDB entry that has no file in the legacy PDB format (large structures come only as mmCIF) says so.
* The licence lines are written from what the databases say on their policy pages (NCBI: public data, but a submitter may claim
  rights; UniProt and AlphaFold DB: CC BY 4.0; PDB: CC0). They were not fetched by the program: **check them before relying on one
  for a publication.**

**Tried live** on 11 November 2026 from one Mac: UniProt, the PDB (including the complex 1VFB) and the AlphaFold DB's file server answered; NCBI's `efetch` redirected
that computer to its abuse page (so NCBI's address is unconfirmed), and the AlphaFold DB's entry service answered 403, which is why the program then looks for the file by its name
(`AF-<accession>-F1-model_v6.pdb`, then v5, v4). The live check is an ignored test that must be asked for: `cargo test --bin helixloom -- --ignored live --nocapture`.

Not built: searching a database (only an accession is fetched), NCBI Datasets, ENA/SRA and Ensembl, an NCBI API key or e-mail (NCBI asks
for them; without them the limit is three requests a second), downloads of many records at once, and any of it in the browser build
(a browser page has no `curl`, and the window is hidden there).
