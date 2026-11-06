# HelixLoom roadmap: what to add for a better experience and better pipelines

This is the answer to "which aspects can I add and implement to improve the end-user experience and the bioinformatics
pipelines". It is written for the person who decides what comes next, in plain words, and every item says **why** it
matters, **what exactly** to build, **what it builds on** in the program today, and a rough **size**
(S: days, M: 1 to 3 weeks, L: 1 to 2 months of focused work). Items marked ★ change the daily experience the most for
the effort.

## Where HelixLoom stands

Working today, on the computer:

* **Sequences**: plasmid editor (maps, restriction enzymes, ORFs, features, primers, PCR, contig assembler), calculators.
* **Signals**: Sanger traces (.ab1, .scf) with quality, trimming and export; sequencing reads (FASTQ) summary and cleaning.
* **Alignment** (built in or MAFFT), **Phylogeny** (distance, maximum likelihood, Bayesian, IQ-TREE), **Trees** (drawing,
  Newick tools), **Structure** (ESMFold prediction, PDB viewer), **Pipelines** (several steps in one run, export to
  Snakemake and Nextflow).
* **Workstation**: `helixloom serve` computes the heavy jobs on a machine with GPUs; the app chooses where each job runs
  (this computer, the workstation, or automatically).
* **New in this release**: the name and interface of HelixLoom (top bar, command palette, toasts, task-first Home, a
  design system), a models folder for files that belong to other projects, and a build for the browser.

Not yet: anything that needs a licence decision, anything that has only been tried on one Mac (see `DEPLOY.md`), and the
items below.

---

## Part 1. The experience

### 1.1 Make everything findable
* ★ **Command palette, second step (S).** It exists (Cmd+K). Make it search *inside the work* too: feature names, enzyme
  names, sequence names, open results, recent files, and commands with arguments ("go to 1200", "find GAATTC", "align
  these 6"). Why: people should never have to know which menu something is in. Builds on `ui/commands.rs`.
* **Recent files and a library (M).** A list of recent documents on Home, and a *library* panel for a project's files with
  search and tags. Why: today every session starts from an empty window.

### 1.2 Never lose work
* ★ **Autosave and session restore (M).** Save open documents and results to a project folder as you go, atomically, and
  reopen everything where it was. Why: nothing frightens a researcher more than a lost afternoon.
* **Undo history you can see (S).** A list of edits with names ("Insert 12 bases", "Reverse-complement") to jump back to.
  Builds on the undo stack of `DocView`.
* **Project folder format (M).** One folder per project: `project.json` (what is in it), the files in open formats, a
  `results/` folder with each run's inputs and settings. Why: shareable, versionable with git, readable in ten years.

### 1.3 See what is happening
* ★ **Job tray (M).** The top bar already shows how many analyses are running; make it a panel: each job with progress,
  where it runs, how long, a cancel button, and the result when done (click to open). Desktop notification when a long job
  ends. Builds on `Job`, the Workstation job list and the toasts.
* **Provenance on every result (M).** Store what produced a result: input checksums, parameters, program versions, date.
  Offer "run again" and "what changed since". Why: this is what makes an analysis reproducible and publishable.

### 1.4 One selection, shown everywhere
* ★ **Linked views and an inspector (M).** Select a feature, an ORF, a peak or a column in one view and see it highlighted
  in all the others; a side panel shows details of what is selected (length, Tm, translation, annotations, links).
* **Compare mode (M).** Two documents side by side (plasmid versions, two traces, two trees) with differences marked.

### 1.5 Guide the newcomer
* **First-run tour and sample project (S).** Three screens that point at the search box, the rail and the job tray; a
  sample project that opens with a plasmid, a trace and a tree.
* **Empty states that say what to do (S).** An empty alignment page should offer "Add open sequences / Open FASTA / Try the
  example" as buttons (several pages already do; make it the rule).
* **Explain in place (M).** A "what is this?" next to every unusual term (Tm, Q score, bootstrap, pLDDT), and error
  messages that say what to do next.

### 1.6 Polish and comfort
* **Typography (S).** Bundle a modern font pair (Inter for text, JetBrains Mono for sequences, both under the SIL Open Font
  Licence) so that headings have real weights; today the built-in font has one.
* **Appearance settings (S).** UI scale, density (compact, comfortable), accent colour, high-contrast theme,
  colour-blind-safe sequence palettes. The theme already follows the system.
* **Accessibility (M).** Keyboard-only use of every page, focus rings, screen-reader labels (egui supports AccessKit),
  minimum sizes of hit targets.
* **Languages (M).** English and Vietnamese first, with all text in one table.
* **Window placement (S).** Tool windows open under the top bar and remember where they were; docking panels later (L).

### 1.7 Getting results out
* **Publication figures (M).** SVG, PDF and high-resolution PNG of maps, alignments, trees, traces and structures, with
  fonts and colours as on screen.
* **Report builder (M).** One click makes an HTML or PDF report of a project: methods text written from the provenance,
  figures, tables, and the programs cited with versions.

### 1.8 Feeling like a real application
* **Installers (M).** macOS app bundle with signing and notarization, Windows installer, Linux AppImage or Flatpak;
  file associations (.gb, .fasta, .ab1, .pdb, .nwk) and "open with HelixLoom".
* **Updates and crash reports (M).** Check for a new version; opt-in crash reports; a "copy diagnostics" button.

---

## Part 2. The bioinformatics pipelines

### 2.1 The pipeline engine
* ★ **Pipelines as data (L).** Today the Pipelines page runs a fixed route. Make a pipeline a list of steps with typed
  inputs and outputs and parameters (a graph underneath, a simple list or a node canvas on top), saved with a project,
  versioned, and runnable on one sample or a whole folder (batch). Why: this turns HelixLoom from a set of tools into a
  platform. Builds on `pipeline/` and the Snakemake/Nextflow export.
* **Caching and resume (M).** A step whose inputs and parameters have not changed is not run again; a failed run resumes
  where it stopped.
* **Tool registry (M).** One place that knows each external program (MAFFT, IQ-TREE, primer3, fastp, minimap2, samtools,
  HMMER, MMseqs2, …): how to detect it, its version, how to install it (micromamba or a container), and its licence.
* **Run anywhere (L).** The same pipeline on this computer, the workstation, a Slurm cluster, or a cloud queue; export to
  Nextflow, Snakemake and CWL; import them. Builds on the Workstation and the `Compute` choice.
* **Reproducibility package (M).** Export a project as an RO-Crate: data, workflow, parameters, versions, checksums.

### 2.2 Data in and out
* ★ **Connectors (M).** Fetch from NCBI (Entrez, Datasets), UniProt, PDB, AlphaFold DB, ENA/SRA, Ensembl; keep a local
  cache; remember where each file came from (and its licence).
* **Formats (M).** GenBank, EMBL, GFF3, FASTA, FASTQ, BAM/CRAM with index, VCF, SnapGene `.dna`, SBOL; large files read by
  index, not loaded whole.
* **Genome viewer (L).** A track-based viewer (like IGV) for reads, variants and annotations on a reference.

### 2.3 Analyses to add, by question
* ★ **Conserved motifs and domains** (you asked for these): motif discovery (MEME/STREME), domain search against Pfam,
  InterPro and CDD with HMMER, per-column conservation scores, mapping conservation onto a 3D structure, protein
  language-model embeddings (ESM-2, ProtT5) for similarity and for the effect of a variant (ESM-1v), clustering of families
  (MMseqs2). Each is a *model group* in the models folder and a job kind on the Workstation (see `DEPLOY.md` section 10).
* ★ **Epitope and paratope** (you asked for these): B-cell epitope prediction (BepiPred-3, DiscoTope-3), antibody
  numbering and region marking (ANARCI), paratope prediction (Parapred, Paragraph), antibody-antigen complex prediction
  (AlphaFold-Multimer, Boltz, Chai; check each licence), docking (HADDOCK, ClusPro), interface analysis (contacts,
  buried area), binding-affinity change for mutations, T-cell epitopes (NetMHCpan). Validate against the SAbDab and IEDB
  benchmark sets and show the accuracy on the page before the page claims anything.
* **Cloning design (M).** Gibson and Golden Gate planning with overhang fidelity, primer design with off-target checks,
  restriction-cloning planner, codon optimization with organism tables, in-silico gel and PCR, CRISPR guide design with
  off-targets.
* **Sequence verification (M).** Map Sanger or nanopore reads to the reference, call the consensus, list differences,
  and give a pass or fail with the evidence (the kind of report sequencing services return).
* **Phylogenetics (M).** IQ-TREE 3 with model finding and ultrafast bootstrap, partitioned models, dating (TreeTime,
  BEAST2), ancestral reconstruction, tree comparison and tanglegrams.
* **Reads and variants (L).** fastp and MultiQC reports, minimap2 or bwa-mem2 mapping, samtools, bcftools or DeepVariant
  calling, assembly of long reads (Flye) and polishing (Medaka), taxonomic profiling (Kraken2).
* **Structure (M).** AlphaFold DB lookup, ColabFold, structure alignment (TM-align, Foldseek), pocket detection.

### 2.4 Trust
* ★ **A benchmark page (M).** For every analysis, say how it was checked and how accurate it is (the golden tests already
  compare Tm, PCR, tree building, DSSP, ABIF reading and the calculators with the reference programs; show that to the user).
* **Model cards and licences (S).** Every model in the models folder lists its source, its licence, the data it was trained
  on, and what it must not be used for.
* **Reference checks in CI (M).** Run the golden comparisons on every change, on macOS, Linux and Windows.

### 2.5 Compute
* **From one workstation to a server for a lab (L).** Accounts and sign-in (OIDC), TLS, quotas, job queues with
  priorities, GPU scheduling, storage on S3 or NFS, audit log, container images, a REST API with an OpenAPI description and a
  Python package so that notebooks can drive the same pipelines. Builds on `src/remote/`.

---

## Part 3. A web version alongside the computer version

One codebase, two targets: the same Rust program compiled to WebAssembly runs in the browser (`docs/WEB.md`). Phases:

1. **W1 (this release)**: the program builds for the web; open files by dropping or choosing them, download results;
   small jobs run in the page.
2. **W2 (M)**: the browser connects to a Workstation over `fetch` with the same token, and the server also serves the web
   files, so heavy jobs run on the GPUs and results come back to the browser.
3. **W3 (M)**: projects kept in the browser (IndexedDB) and shareable links to read-only results.
4. **W4 (L)**: accounts and shared projects on the server; background jobs in web workers; installable app (PWA) that
   works offline for the editor and viewers; WebGPU for rendering large trees and alignments.

## Part 4. The new interface, step by step

Phase 1 is in this release (`docs/DESIGN.md`): colours and typography, top bar, command palette, toasts, Home.

2. **Phase 2 (M)**: bring every workspace to the same look: one header, one controls panel, one results area, real empty
   states, linked views, the inspector.
3. **Phase 3 (L)**: docking panels (drag a panel to a new place), tabs for results, saved layouts.
4. **Phase 4 (M)**: motion and polish: smooth transitions, loading skeletons, the typography pair, an illustration set.

## Part 5. Suggested order

| Milestone | Contents | Why first |
|-----------|----------|-----------|
| **A (2 to 3 weeks)** | Autosave and session restore; recent files; job tray; window placement | Removes the fear of losing work and the feeling of not knowing what runs |
| **B (3 to 4 weeks)** | Linked views and inspector; empty states; first-run tour; typography | Makes the daily screens feel finished |
| **C (4 to 6 weeks)** | Pipelines as data with caching; connectors (NCBI, UniProt, PDB) | Turns the tools into a platform |
| **D (4 to 6 weeks)** | Motifs and domains (HMMER, Pfam, conservation on structure) | First of your two research directions; light compute |
| **E (6 to 8 weeks)** | Epitope and paratope models on the Workstation; benchmark page | Second direction; needs GPUs, licences and validation |
| **F (parallel)** | Web W2 and W3; installers and updates; CI on three systems | Reach |

## Questions for the owner

1. **Licence** for the code (MIT, Apache-2.0, GPL-3.0, or closed)? The repository has none yet, so nobody else may reuse it.
2. **Who is it for first**: wet-lab biologists who clone and sequence, or computational people who run pipelines?
3. **Which models first**, and may some be built into the program or must all run as separate programs on the workstation?
4. **Hosting** for the web version: your own server, a university server, or a cloud?
5. **PhyloCNN**: ask its authors for permission to share the trained networks, or leave them as a local add-on?
