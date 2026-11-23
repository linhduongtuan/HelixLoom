# HelixLoom roadmap: what to add for a better experience and better pipelines

This is the answer to "which aspects can I add and implement to improve the end-user experience and the bioinformatics
pipelines". It is written for the person who decides what comes next, in plain words, and every item says **why** it
matters, **what exactly** to build, **what it builds on** in the program today, and a rough **size**
(S: days, M: 1 to 3 weeks, L: 1 to 2 months of focused work). Items marked ★ change the daily experience the most for
the effort.

## In short (November 2026)

* **Where a heavy job runs is already your choice**, for every kind of job: once a workstation is connected (Workstation
  page), every page and window that starts a heavy job shows the choice *This computer*, *Workstation* (the machine with
  several GPUs, over `helixloom serve`) or *Automatic* (big jobs go to the workstation, small ones stay), with a line that says
  where this run will go; the choice is kept between sessions. The new species-tree job uses it like the others
  (`docs/SPECIES.md`). **How many GPUs one job gets is also your choice** now (Workstation page, *GPUs for one job*: 1, 2, …
  up to what the workstation has, or all): a job that can use several GPUs (PhyloCNN with many trees shares them out) is
  given that many at once, the others wait their turn in the order they were sent, and `--fold-gpus` lets a folding program
  use more than one. What is not there yet: choosing *which* GPUs by number, the dynamic programming of an *alignment* on GPUs by a fused kernel (see 2.5; an eager PyTorch version exists, opt-in, and loses 60 to 1,000 times to the processor, whose alignments are 7 times faster since November 2026), the CUDA path for
  several GPUs run on real hardware (see the next point), other tensor libraries than candle (burn, tch-rs and dfdx are not
  built), and showing a cost or time estimate before a job is sent. **The web version has the same choice** (W2, below): the browser
  reaches the workstation through `fetch`, and a job there no longer stalls the page.
* **Trees on several GPUs with candle or PyTorch: built, and not yet run on a CUDA card.** The same work (distances, bootstrap, the
  likelihood of a tree) can run through candle, built into the program, or through **PyTorch in Python worker processes, one per GPU**
  (`HELIXLOOM_TORCH_PYTHON`; no CUDA toolkit needed to build the program). On a Mac's GPU, PyTorch did the likelihood of a 64-tip tree over
  100,000 patterns 16 times faster than the processor where candle was slower than it, and the distances 95 times faster against candle's 61
  (`DEPLOY.md` 1.7); nobody knows yet how either does on an RTX 3090 with an RTX 5060 Ti. Everything below about candle holds for both.
* **Trees on several GPUs with candle: built, and not yet run on a CUDA card.** Distances between all pairs of sequences,
  the distances of every bootstrap resample, the likelihood of a tree (so maximum-likelihood search, the fit of the four models
  and the Bayesian sampler) all run on every GPU the job is given, cut by each card's speed (measured, or `HELIXLOOM_GPU_WEIGHTS`)
  for a card pair such as an RTX 3090 with an RTX 5060 Ti. A job is sent to the GPUs only when it is big enough, and a tree job on a workstation without GPUs runs on
  its processor. Distances equal the processor's to the last digit; the likelihood agrees to about seven digits. It was tested on
  candle's CPU device in unequal pools, with many threads, and on a Mac's GPU, where it showed a real fault (several threads on one Metal device gave silently
  wrong numbers) that is now guarded against, and an independent review of the code, which read candle's CUDA backend, found problems that are now fixed (a server without GPUs
  refusing a big tree, a card that cannot run the work disabling the other, silent fallbacks to the processor, a noise level in the likelihood that
  could pass for a better tree). **The speed-up on your cards is unknown** until `helixloom gpu-check --bench` has run on them
  (`DEPLOY.md` 1.7): on the Mac the distances were 5 to 38 times faster and a likelihood was slower than the processor below about 6 million cells, because candle's small operations are
  limited by launching. Not on the candle GPUs: alignment (fast on the processor now: see the table below), neighbour-joining (PyTorch does it) and UPGMA themselves, species trees, pipelines.
* **What runs on the GPUs, and what does not** (November 2026; "PyTorch" needs `HELIXLOOM_TORCH_PYTHON`, see `DEPLOY.md` 1.7):

  | Work | On the GPUs? |
  | --- | --- |
  | Distances between sequences, and every bootstrap resample's | yes (candle or PyTorch): equal to the processor's to the last digit |
  | Likelihood of a tree: maximum-likelihood search, the fit of the four models, the Bayesian sampler | yes (candle or PyTorch): agrees to about seven digits |
  | Choosing the centre of an alignment (the built-in aligner's k-mer counts) | not GPU work: made one pass on the processor (each 6-mer's number of holders), the same sums as the pairwise way, which took half a minute for 600 sequences |
  | Neighbour-joining of thousands of sequences | yes (PyTorch, on cards with 64-bit floats): the processor's tree where no pairs tie |
  | PhyloCNN | yes (candle) |
  | The pairwise alignments of the built-in aligner (dynamic programming) | **no, and no longer needed**: they are on the processor, eight pairs to a vector instruction, 1.3 ns a cell on one core (10.4 before); **a PyTorch version exists, opt-in** (`HELIXLOOM_ACCEL_ALIGN_CELLS`): the same alignments to the last bit, 60 to 1,000 times slower than the processor on a CPU worker and a Mac's GPU; only a fused kernel could beat the processor, and it needs a CUDA machine to write and check (2.5) |
  | MAFFT, IQ-TREE, folding programs | other programs: they use what they use (a folding program is given GPUs) |
  | UPGMA, species trees (ASTRAL), conservation scores, pipelines, drawing | processor |
* **ASTRAL is rebuilt and checked against the real ASTRAL; the rest of ASTER is not.** Tools → *Species Tree from Gene Trees*
  and `helixloom astral` give the same scores, supports and branch lengths as ASTRAL 5.7.8 on every comparison that was run,
  several times faster. Weighted ASTRAL, gene families (ASTRAL-Pro3), CASTER and WASTER are the next steps (2.3).
* **Conserved motifs and domains, then epitopes and paratopes: steps 1 and 3 of the plan are built (3 in its native form), the rest is planned** (2.6).
  Built: per-column conservation (Jensen-Shannon divergence after Capra and Singh 2007, and Shannon entropy, with sequence
  weights, a gap penalty and a window), conserved blocks, and the same scores as colours on a 3D structure. The formulas are
  checked against hand-worked cases, a closed form, and a separate implementation of the same formulas (which an independent
  review showed can check the arithmetic but not the data: it had passed a program with three background frequencies swapped);
  **the accuracy that the paper reports for predicting functional residues was not reproduced here** (no benchmark data, and
  the authors' program was not run), so treat the scores as a guide to where to look, not as a verdict. Known limits: lining a
  structure up with the alignment runs on the window's own thread (a few alignments of a chain against the likeliest rows: a
  chain of thousands of residues takes seconds); the turquoise-maroon colours are poor for red-green colour blindness (a
  residue's score is written when it is clicked, and a colour-blind-safe ramp belongs with 1.6). Next, in the plan's order and each needing a decision of yours (which models, which licences):
  protein language models, the antibody tools and hold-out benchmark sets that have to pass before a page claims anything
  (the HMMER domain search against Pfam is built and was run for real, see step 2 of 2.6; the page that says how each analysis
  was checked is built, 2.4).
  **Motifs (step 3), native part built: Tools → Motif Finder.** `seq/motif.rs` scans DNA with a position weight matrix on both strands, with the
  p-value of a score worked out exactly (checked against scoring all 256 words of a four-column motif one by one: it caught a rounding fault that was
  fixed), finds motifs by expectation-maximisation from the words of the sequences (MEME's method in its simplest form: at most one site in a sequence,
  either strand, several motifs one after the other), compares motifs by correlation, and writes MEME's text format for FIMO and Tomtom. The window shows the
  logo, the sites and a scan, and runs the search on its own thread. **Checked on planted motifs only** (a word with one base in ten changed, on either strand, found where it was put
  and as the consensus; two motifs found one after the other): *not* against MEME on real promoter sets, so it has no claim to MEME's accuracy, and it
  is bounded in time on large sets (about five seconds of work whatever the size: the starting words, ranked by how many sequences hold them, and the iterations are cut down, so a big set is searched less thoroughly; the window takes at most 300,000 bases, 60,000 in a browser). **The command `helixloom motif` and the workstation job (kind `motif`, with the choice of where it runs in the window) are built, November 2026: `docs/MOTIFS.md`**
  (an *effort* setting scales the budget of work, and a search refuses windows that would not fit in 4 GB of memory; checked on planted motifs, not against MEME).
  Not built: known motifs to compare with (JASPAR is free to
  use but is a download), protein motifs, `STREME`-style discovery of many short motifs, and scanning on the GPUs (a scan is a
  convolution of one-hot sequences, which PyTorch does well: a step for when sets are large).
* **The ten changes that would improve the daily experience and the pipelines the most**, in the order to do them:
  1. autosave and session restore (1.2); 2. a job tray with progress, cancel and notifications (1.3); 3. pipelines as data,
  with caching and a batch mode (2.1), starting with *a folder of locus alignments → gene trees → species tree*;
  4. connectors to NCBI, UniProt, PDB and AlphaFold DB, remembering where each file came from (2.2); 5. provenance on every
  result, with "run again" (1.3; built for nine kinds of result: alignments, trees, species trees, motifs, protein domains, PhyloCNN estimates, predicted structures, and runs of the Pipelines page's fixed route and of pipeline files); 6. linked views and an inspector (1.4; built for the Sequences page, for a column of the alignment with the residue of the structure that is lined up with it, for the ORF list and for the bases of a trace); 7. a benchmark page that shows how each analysis
  was checked (2.4; built: Help → *How Each Analysis Was Checked* and `helixloom checks`); 8. publication figures and a report builder (1.7; built in part: SVG and PNG figures of five views, and an HTML report of the nine kinds of result that carry a record); 9. the tool registry that detects, versions and
  installs the external programs (2.1; built: `helixloom tools`, `helixloom tools install NAME` and Tools → *Programs HelixLoom Calls…*, which install only on a yes); 10. installers and updates, so that others can run it at all (1.8; built in part: a Mac app and its disk image, signed ad hoc, and Help → *Copy Diagnostics*; a version check and a developer signature wait for questions 11 and 12).

## Where HelixLoom stands

Working today, on the computer:

* **Sequences**: plasmid editor (maps, restriction enzymes, ORFs, features, primers, PCR, contig assembler), calculators.
* **Signals**: Sanger traces (.ab1, .scf) with quality, trimming and export; sequencing reads (FASTQ) summary and cleaning.
* **Alignment** (built in, which aligns eight pairs at once on the processor, or MAFFT), **Phylogeny** (distance, maximum likelihood, Bayesian, IQ-TREE, and the species tree
  from many gene trees, the ASTRAL method), **Trees** (drawing, Newick tools), **Structure** (ESMFold prediction, PDB
  viewer), **Pipelines** (several steps in one run, export to Snakemake and Nextflow).
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
  **Built in part (November 2026, `ui/answers.rs`):** the palette answers `go to 1200` (also `goto`, `jump`, `base`, `pos`, or the number alone; past the end it goes to the end),
  `find GAATTC` (also `search`, `seek`; the Find window's settings decide the match), and two letters or more of the name of an open
  sequence (switch to it), of a feature of the sequence in view or of its kind (select it, as a click on it does) and of a recent file (open it).
  **Not built:** enzyme names, ORFs, open results (an alignment, a tree), the peaks of a trace, and commands such as "align these 6" that take an argument
  of another kind.
* **Recent files and a library (M).** A list of recent documents on Home, and a *library* panel for a project's files with
  search and tags. Why: today every session starts from an empty window.
  **Built: File ▸ Open Recent** (the twelve latest files opened or saved, the latest first, a file that has gone left out and dropped when chosen,
  Clear Menu). **Also built: the list on Home** (November 2026): the five latest files, each a button that opens it, under the first card of Home ("Pick up where you left off"; nothing when there are none, and none are kept without a settings folder). **Not built:** the library panel.

### 1.2 Never lose work
* ★ **Autosave and session restore (M).** Save open documents and results to a project folder as you go, atomically, and
  reopen everything where it was. Why: nothing frightens a researcher more than a lost afternoon.
  **Built for the open sequences (November 2026):** when the open sequences have been quiet for two seconds (or a change has waited fifteen), and
  when the window is closed or the program ends, they are written to the `session` folder of the settings folder (`~/Library/Application Support/HelixLoom`
  on a Mac, `~/.config/helixloom` on Linux, `%APPDATA%\HelixLoom` on Windows). A sequence that has its own file and no unsaved change is kept as a path
  only; any other (a new or pasted one, or one with changes not saved) as GenBank text, under its own name, so nothing that was typed is lost to a
  crash; one that has not changed since the last save is not written again. A save writes new files, puts them on disk (`fsync`), replaces the manifest in one
  step, then removes the old files, so that whatever moment the program stops at, a power cut too, the folder holds a whole session (tested by stopping a
  save halfway); the folder is private to the person (0700, files 0600 on Unix), only files the program wrote are ever removed, and a session of another
  version is put aside (`session.unreadable`) instead of written over. At the next start a dialog that blocks the window asks "Reopen your last session?", names the
  sequences and which have unsaved changes, and reopens them as they were, with the tab and workspace; "Start fresh" for unsaved work asks once more.
  If the person goes on working without answering, the offered session is put aside (`session.previous`) and the new work is kept as usual, so closing at any moment
  loses nothing. A file that is open already is not opened again (unsaved changes to it come back as a copy that is not tied to the file); one that cannot be
  read is counted and said. The example plasmid and an empty new sequence are not kept until touched. Two copies of the program running at once do not share a session
  (a lock file): the second keeps none and says so. A save that fails is said once, in red, and tried again after five seconds.
  **Built, 18 November 2026: what the Alignment, Trees, Structure and Pipelines pages show is kept with the session** (`src/ui/results.rs`): the Alignment page's list of
  sequences with their ticks, its choices and its alignment, the tree, the structure, and the runs of the Pipelines page in both of its modes, each with the record of what made it. They are written beside the sequences
  (a file each, `generation-kind.result`; one that has not changed is not written again), named in the question at the next start (with no sequence open they are a
  session by themselves) and put back on their pages by *Reopen*; one that cannot be read is said, and the rest comes back. **The sequences never wait for a result:**
  a result that cannot be written is left out and said once, and tried again at every save; when a sequence cannot be written while results are kept, the results
  make room. A run of the Pipelines page is kept in the text it has for the workstation; the folders of a pipeline file's run are not kept. **Not built:** what the tool windows show (PhyloCNN, species trees, domains, motifs), a project folder
  with `project.json` and a `results/` folder, and a sequence's undo history. **In a browser the same is kept since 16 November 2026**, in the browser's storage for the site instead of a folder (`docs/WEB.md`, *What the page keeps*). **Known limits:** the text of a very large unsaved
  sequence is written on the window's own thread, which takes a moment for tens of megabytes; and the list of recent files looks at each file once at the start, which can wait on a
  disk that has gone away.
* **Undo history you can see (S).** A list of edits with names ("Insert 12 bases", "Reverse-complement") to jump back to.
  Builds on the undo stack of `DocView`.
* **Project folder format (M).** One folder per project: `project.json` (what is in it), the files in open formats, a
  `results/` folder with each run's inputs and settings. Why: shareable, versionable with git, readable in ten years.

### 1.3 See what is happening
* ★ **Job tray (M).** The top bar already shows how many analyses are running; make it a panel: each job with progress,
  where it runs, how long, a cancel button, and the result when done (click to open). Desktop notification when a long job
  ends. Builds on `Job`, the Workstation job list and the toasts.
  **Built, without progress bars and without a cancel button:** the top bar's count is a button that opens a panel of what is running (how long
  each has run, and a Show button that goes to its page) and of the last eight that ended (how long each took, how long ago, and an Open button);
  an end is a toast that names the analysis ("Phylogeny ended after 2 min 34 s."; "ended", because the tray cannot tell an answer from a failure: the page says which) and, for one that ran twenty seconds or more while the
  window was not in front, a notice from the operating system (the notification centre on a Mac through `/usr/bin/osascript`, `notify-send` on Linux where it
  is installed; **neither was tried on a Linux machine**, and nothing is shown on Windows or in a browser). **The tray and the notice work only while the window is being drawn:** a minimised or fully covered window
  runs no frame, so an analysis that ends then is noticed when the window comes back (the session is kept meanwhile). The pill counts a motif search too; what ended is forgotten after fifteen minutes. There is **no cancel**: the analyses run on
  threads or wait for a workstation and none can be stopped half way (the Workstation page stops a job there), and **no progress**: the programs they
  call (MAFFT, IQ-TREE) report none. A job that the person starts from another computer is on the Workstation page, not in the tray.
* **Provenance on every result (M).** Store what produced a result: input checksums, parameters, program versions, date.
  Offer "run again" and "what changed since". Why: this is what makes an analysis reproducible and publishable.
  **Built for three results (November 2026; `src/provenance.rs`, `src/ui/provenance.rs`):** an alignment (Alignment page), a tree (Phylogeny page) and a species tree
  (its window) each carry a record of what made them: every input by name, size and checksum (the ticked sequences; the rows of the alignment; the gene trees, the tree to score and
  the bootstrap trees), the settings the method uses, the programs with the versions they gave (MAFFT's, and this program's), the day (UTC), where it ran and how long it took. Under the
  result, *How this was made* (closed until clicked) shows the record, says **what has changed since** (an input that changed, went or came, their order, a setting that changed) or that
  nothing has, **Copy the record** gives it as text, and **Run again with these settings** puts the settings of the result back and runs on the inputs as they are now. A result that
  came from elsewhere (an alignment from a pipeline, a tree opened from a file, a job fetched from the workstation's list) carries no record, and shows none; nor does a tree built in the
  Phylogenetic Tree window, or one rooted at its midpoint or ladderized on the Trees page (running the estimation again would not give that tree). When MAFFT is chosen but is not installed
  where the alignment runs, the built-in aligner runs, the choice on the page goes back to it, and the record says what ran. An independent review of this (two reviewers and a skeptic)
  found those two faults and three tests that proved too little; all five are mended.
  **Built for six more results (14 November 2026)**, each with the same section (*How this was made*, what has changed since, *Copy the record*, *Run again with these
  settings*, and the link to how its kind of analysis was checked):
  * **the motifs** of the Motif Finder: every sequence by its name, and the six settings of the search;
  * **the domains** of Protein Domains: every protein (by the first word of its line, as HMMER names it), the threshold (Pfam's own, or the E-value), and as programs the HMMER
    and the Pfam file that the answer names (an answer that names none, as an older workstation's might, leaves no empty name);
  * **the estimates** of PhyloCNN: the trees as they were given, what was estimated, the priors and the sampling probability, and as program the authors' networks with whether
    they ran on the CPU or a GPU. The networks are not installed on the Mac this was written on, so the record of a real run was not seen: the tests give the window a made-up
    answer, and one test takes the real path where the networks are installed;
  * **a predicted structure** (Structure page): the sequence that was sent, who was asked to fold it (ESMFold through the ESM Atlas service, or the workstation's own folding
    program) and what the answer says made it; the time is from the asking to the answer. *Run again* is offered only while who would fold now is who made the one that is
    shown. Opening any other structure takes the record away;
  * **a run of the fixed route** (Pipelines): every sequence, the steps, whether the installed programs were used, the bootstrap, where the primers go, and the primer settings
    (Tools → Primer3) as an input with a checksum; as programs what each of the three parts says made it. A result of the fixed route that is fetched from the
    workstation's list takes the record of the page's own last run away;
  * **a run of a pipeline file**: every file of the folder by its name, with the size and the checksum of what was read from it (its sequences or trees as text, not the bytes
    on the disk), the pipeline step by step, and the folder by its name, the number of its files, their bytes, when the newest was changed, and the start of a checksum of
    its listing (every file's name, size and time of change). As programs it names the aligner and the tree builder of the steps (without the model, which a tree step's
    note goes on with); what the species step found is not a program. *What has changed since* does not read the files again (a frame must not): while the record is
    **open** the folder is listed every two seconds (on the thread that draws: a folder on a slow drive would make the page stutter then), and a different listing says that
    the folder changed: a file added, taken away, renamed, resized, or changed at another time, even put back from an older copy. A file written again within the same
    second to the same size is not seen as changed.
  A result fetched from the workstation's list carries no record in any of the six (what it was made from is not known here), and a result that goes takes its record with it.
  An independent review of the six (three reviewers, each on a part) found ten things, all mended with a test each: where the example of the Motif Finder and of Protein
  Domains ran was said from the search before it; a pipeline result fetched from the workstation was shown beside the record of the page's last run; what the species step
  found was listed among the programs; two lines under figures of the report said more than is always true (see 1.7); the folder was listed every two seconds on every page
  whether the record was open or not, and a renamed file or one put back from an older copy was not seen; and an answer that names no program left an empty name.
  The code of the six records, their parts of the report and the lookup of their checks was then broken on purpose in 94 places, one at a time: the tests caught 93, and
  the one they let through (a tree with an infinite length must still be drawn in a report) got its test.
  **Not built:** the other results (reads, primers picked in the Primer3 window, digests, contig assemblies, a tree built in the Phylogenetic Tree window); the record saved beside a file that is saved, or kept in a
  project (since 18 November 2026 the records of the alignment, the tree, the structure and the runs of the Pipelines page are kept with the session, with the result itself; the record of any other result lives as long as the result is on the page); the inputs themselves (only their checksums are kept, so *Run again* uses the inputs as they are now, and the list says
  if they differ); for a run on the workstation the version named is this app's, since the workstation's is not reported with a result; and the checksums are made anew on every frame
  while the section is open, which a very large set would feel.

### 1.4 One selection, shown everywhere
* ★ **Linked views and an inspector (M).** Select a feature, an ORF, a peak or a column in one view and see it highlighted
  in all the others; a side panel shows details of what is selected (length, Tm, translation, annotations, links).
  **Built for the Sequences page (November 2026):** the sequence view, the feature table and both maps (circular and linear) select the same
  feature: a click on its arrow on a map selects it (it was only a place before), and the selected feature is outlined heavily on the maps. The
  **Inspector** (View ▸ Show Inspector, under the map) says what is selected: for a feature its name, type, location as GenBank writes it, length,
  strand, GC, Tm (only up to 100 bases, the original's limit), the original's ORF summary, the translation of a CDS on its own strand with the
  genetic code in use (the first 60 residues), its notes (qualifiers) and the other features that share a base with it (click one to select it); for
  bases selected without a feature, the same numbers.
  **Built for a column of the alignment and the residue of the structure (16 November 2026; `ui/platform/linked.rs`, `ui/platform/alignment/column.rs`,
  `ui/platform/structure/link.rs`):** on the Alignment page a click on a letter, or on the ruler, the bar or the consensus letter above it, selects its
  **column**: it is outlined over all the rows and tinted at the top, and one line above the alignment says what it holds (*Column 44 of 76 · consensus I ·
  I 10, A 2 · no gaps · Identity 0.83*, with *among the conserved* and the conserved block it is in where that is so; a column that is mostly gaps says
  that it is not scored). A click on it again, or *Deselect*, lets it go. **Once a structure has been coloured by the alignment** (*Colour by
  conservation* on the Structure page, which lines each chain up with a row), each residue knows its column, and the column and the residue are one
  selection on both pages: a column that is selected picks its residue in 3D and on the sequence strip, a residue that is picked (in the picture or on
  the strip) selects its column, and the Alignment page scrolls that column into the middle when it is next shown; letting go of one lets go of the
  other (*Reset view* too). The line of the column has a link to the residue (*Show I 44 of chain A in 3D*), the Structure page has one back
  (*I 44 is column 44 of the alignment*), and the label of the picked residue in the picture names its column. With copies of one chain every copy is
  in the column: the first is picked, the others are counted on the link, and a copy that is picked by hand stays picked. A residue that no row has
  (a tag, a chain that no row matches) selects no column; a column that no residue is in picks none.
  **What ends the link:** another alignment, or another score chosen on the Alignment page, takes the colours off the structure and the columns with
  them (*Colour by conservation* makes them again); another structure has no columns until it is coloured. A structure that is not lined up leaves the
  column alone, and a column leaves such a structure alone. Two views for pictures show it: `--view linked-column` and `--view linked-residue` (an
  alignment of ubiquitin with variants that are **made up for the picture**, and say so in their names).
  **Built for an ORF of the ORF list and for the bases of a trace (16 November 2026; `ui/windows/orfs.rs`, `seq/peak.rs`,
  `ui/platform/signals/pick.rs`):**
  * **The ORF list** marks the ORF whose bases are exactly the selected ones, wherever they were selected (a row or a bar of the list, Find Next
    ORF, a feature on a map, a drag): its row is marked and scrolled to, its bar on the map of the six frames has a heavy line round it, and the
    heading names it (*2 ORFs in 2,724 bp · the selection is the forward ORF at 288–902*), so that colour is not the only sign. The list scrolls to
    the ORF once, when it becomes the selected one or the list is opened on it, and can be scrolled away from it. Of an ORF that runs past the end of
    a circular sequence a click selects the part before the end, as it did; the heading says that it is that part, and the line is round that piece of
    the bar only. What the selected bases are is in the inspector, which reads the same selection.
  * **A trace** (Signals page): one line above the trace says what the picked base is (*Base 157 of 480 · called G · quality 42 (1 wrong in 15,849) ·
    scan 1,902 · signal G 1,480, A 120, C 30, T 25 · kept*): its number as the ruler writes it, the call, how sure the machine was and what that is as
    a chance, the scan of its peak, the signal of the four channels there with the strongest first (and *the strongest signal is A* when another
    channel than the call's own is the strongest: a mixed peak, or a call to look at again), and whether trimming keeps it. Of bases picked by a drag
    it says how good they are together (*Bases 150–180 of 480 (31) · mean quality 56.6 · lowest 52 at base 165 · none below quality 20 · all kept*,
    and how many calls are not A, C, G or T). In a read that is shown reverse-complemented the line speaks as the picture does: the bases are numbered
    from the other end and the letters are the complements. A file with no qualities says so.
  * **Three older faults that the new tests, the pictures and the review showed, all mended with a test each.** Bases picked by dragging across a
    trace collapsed to the one base under the pointer when the button was let go, so *Drag across bases to pick them* had not worked (the first base
    of a drag is kept from the frame the drag begins now; a test drags as a hand does, to the right, to the left and in a mirrored picture). The ORF
    window grew a button's height in every frame until it was as tall as the screen, its title over the top bar and its button off the bottom (it
    keeps the size it is given, and a test draws ninety frames of it). And the map of the six frames wrote base 1 as the end of an ORF that ends on
    the last base of the sequence, where the table wrote the last base (both ask one function now).
  * The line that says what is selected is one function (`design::selection_line`), which the column of the alignment uses too.
  * **An independent review** (two readers: the code, the tests) found, beside the third fault above: a hint that promised what does not happen (a
    click beside the bases does not let the pick go: it picks the nearest base, and a test says so now); a pick that could be a few bases short when
    the last move and the letting go came in one frame, and that began at another base if the trace scrolled under a held button; a list that did not
    scroll to the selected ORF when it was closed and opened again; and the heading and the line for an ORF that runs past the end (above). All are
    mended. The reader of the tests named twelve gaps, each closed with a test: whether trimming keeps a base as the page itself says it; a drag to
    the left and in a mirrored picture; the list's scrolling (to the ORF, not pinned to it, again when opened anew); the line round the bar; a click on
    a bar; the mark of a row by its colour and not by a count; the size of the window to fifteen points; the mean of bases of which some have no
    quality; and a click on each of the two links of the alignment's line when both are there. The code of this part was then broken on purpose in 73
    places, one at a time: the tests caught 71, and the two they let through (the stretch must be picked while the button is still down, not only when
    it is let go; an ORF that ends on the last base is not a part of an ORF) got their tests and are caught.
  **Not built:** a base of a trace tied to the base of the sequence that was made from it (*Open in Sequences* keeps no tie to the trace), a row of
  the alignment with the tip of a tree, moving
  the column with the arrow keys, a column in the palette's answers (`column 120`), links out of the inspector, and Compare mode. **Known limits of
  the trace's line:** the mean is the plain mean of the qualities (not of the chances they stand for); *the strongest signal* is read at the one scan
  the file names as the peak; a click anywhere in the trace picks the nearest base, so the pick is let go with *Deselect* or *Clear the pick*. **Known limits of the
  column:** the link rests on the lining up that *Colour by conservation* makes, so it has that step's limits (a chain is matched to a row that is at
  least 70 % the same; the lining up runs on the window's own thread); choosing another score ends the link although the columns are the same columns;
  facts too long for the line are cut short there and read whole under the pointer; the line and the outline are not in a report (a report is of the
  alignment), but a figure saved from the page shows the outline, as the window does; the column is scrolled to about the middle (a point to eight left
  of it, by the edge of the frame and a scroll bar). Tried on a Mac only, and both pages were looked at as pictures.
  **An independent review** (three readers, each on a part: the link and what it keeps, the viewer and the facts of a column, the tests) found one fault
  and three small things, all mended with a test each: *Colour by conservation* pressed while a column that no residue is in was selected let go of that
  column (the column stays now, whatever residue was picked); a column that is a gap in every row was written with an empty place between two dots; a
  click in the spare room under the last row selected the column above it; and the facts of the column copied the whole alignment in every frame (they
  are read where they are now). The reader of the tests named twelve places where a test proved less than it should, and each has its test or its
  mend now: that another score ends the link; the copies of a chain (a test of the name had no copies in it); a click after the view has scrolled; the
  last column of a block; that the page does not scroll to a column that is not there; the tint above the column; a click on a residue in the picture
  itself, and beside the structure; which of the two wins when a column is selected and another residue is picked before they are lined up; that a click
  on a column goes through the link; a pixel distance that a test had by heart; and the two places that said what a gap is, which are one now. **Left as
  they are:** the place a test clicks is worked out from the numbers of the ruler, so a column drawn a few points off would not be noticed; and one
  thing a reader thought possible cannot happen (colouring a second time and matching nothing: an alignment that has changed has taken the colours off
  before the button can be pressed). The code of the column, the link and the viewer's part in them was then broken on purpose in 79 places, one at a
  time, and the tests caught every one. **Known limits:** a feature is "the selected one" only while the selection is
  exactly its bases (drag or arrow off it and the panel reads the bases instead), but Edit and Delete Feature, and the highlighted row of the feature table,
  still use the last feature clicked, as they did before; the translation ignores a `/codon_start` qualifier, as the rest of the program does; the inspector takes
  the lower part of the map's column (about 45 % of it at the start, and it can be dragged), so on a small window the map is shorter. Tried on a Mac only.
* **Compare mode (M).** Two documents side by side (plasmid versions, two traces, two trees) with differences marked.

### 1.5 Guide the newcomer
* **First-run tour and sample project (S).** Three screens that point at the search box, the rail and the job tray; a
  sample project that opens with a plasmid, a trace and a tree.
  **Built:** a card on the first run that points, with a ring round it, at the search box, then the pages on the left, then where jobs run (the
  Jobs button of the job tray is not there until something has run, so the card points at the pill beside it and says what will appear); Next, Skip
  the tour, Esc. It is shown once for a settings folder (the file `tour-seen` in it), not when the program was started with files or a picture, and
  not while a saved session is being offered. Help ▸ Show the Tour brings it back. The last card offers **Help ▸ Open the Sample Project** (the demo plasmid,
  the demo trace in Signals, the demo tree in Trees; the plasmid is not opened twice, **but the trace and the tree that were there are replaced**, as the pages' own *Example* buttons do). **Built, 16 November 2026: the tour in the browser** (that it was seen is kept in the browser's storage). **Not built:** a project that is saved as one thing (it is three pages' examples opened together).
* **Empty states that say what to do (S).** An empty alignment page should offer "Add open sequences / Open FASTA / Try the
  example" as buttons (several pages already do; make it the rule).
  **Built for the six pages that wait for a result:** Alignment (Add open sequences, Open FASTA…, Try the example), Signals, Trees and Structure (Open…,
  Try the example), Phylogeny (Try the example, Add sequences in Alignment) and Pipelines (Try the example, Open a sequence) now put the buttons in the empty
  card itself (`design::empty_state`); tests find each button on the drawn screen and click it. **Not done:** the Workstation job list, the tool windows
  and the dialogs (they say what is missing in words), and a check that a new page uses the shared card.
* **Explain in place (M).** A "what is this?" next to every unusual term (Tm, Q score, bootstrap, pLDDT), and error
  messages that say what to do next.

### 1.6 Polish and comfort
* **Typography (S).** Bundle a modern font pair (Inter for text, JetBrains Mono for sequences, both under the SIL Open Font
  Licence) so that headings have real weights; today the built-in font has one.
  **Built (November 2026):** Inter for the text and JetBrains Mono for the sequences, protein and code, each one variable font file in
  `assets/fonts` (together about 1.2 MB, which is added to the program and to the web build); text that stands out (headings, names,
  table heads) is drawn at weight 600 of the same file (`RichText::semibold`, which every former `strong()` now uses). egui's own fonts
  stay behind them for signs that these two lack. **Not done:** a bold weight for the monospace sequence text, and a check on Linux and Windows (a variable font is drawn by egui itself, so it should look the same, but
  only a Mac was tried). Inter is wider than the old font, so some text wraps to more lines than before (the cards on the Home page).
* **Appearance settings (S).** UI scale, density (compact, comfortable), accent colour, high-contrast theme,
  colour-blind-safe sequence palettes. The theme already follows the system.
  **Built (November 2026):** View ▸ Appearance… (and the command palette; `--view appearance` draws it) holds the theme, which is now kept
  between runs (it was not), and: eight sizes of the whole window, 80 % to 200 % (Cmd +, Cmd − and Cmd 0 from anywhere, in the steps of that
  list and not egui's tenths, which are turned off); two densities; six accents (the indigo as designed, blue, teal, green, purple, and a
  graphite for those who find colour distracting); a high-contrast setting; and colour-blind-safe colours for the bases and residues of
  the alignment, the trace letters and the icons that show them, the consensus of the contig view, the logo of the Motif Finder and the editor's
  palette. Everything is used at once and kept in `appearance.txt` of the settings folder, a line to a setting, so that a line that is not
  understood spoils no other. The window says whether the choices are kept (there may be no settings folder; a write may fail; or the run was started
  with settings on the command line).
  * *What each one guarantees, and what a test measures.* **Accents:** each reads on the page (4.5 to 1 or more) and carries its own text,
    in both themes, and the text of a chosen button reads on its ground (4.5; the designed indigo is left as it was, and its chosen buttons read at 3.9
    to 1 on a light page and **2.7 on a dark one**, below 4.5: any other accent, or high contrast, mends that, and the indigo itself may be changed
    if you say so). **High contrast:** text, dimmed text, links and the status colours at 7 to 1 or more, and the edges of boxes at 3 to 1,
    against the page, the cards, the grounds that controls sit on (a box, a button under the pointer) and a chosen button's tint, in both themes with every
    accent; lines a half thicker; letters in the alignment black or white by their cell. **Density:** compact keeps the ordinary buttons and rows of the
    window at 24 points or more, the least that WCAG 2.2 (2.5.8) allows for a target (it is about 19 CSS pixels at 80 %); **small buttons
    (18 points, measured; 18 places make one), the rows of the inspector (15) and those of the alignment (20) are not changed by it.** **Defaults:** they leave the
    colours, the spacing and the editor's palette exactly as they were (a test says so); the one change at default settings is the left rail's, which
    did not count the space between its buttons and ran off a window that was short, and now squeezes its buttons correctly and scrolls when even the
    shortest do not fit. **Safe colours:** see below.
  * *How the safe colours were checked.* The bases and trace letters are Okabe and Ito's palette, softened for cells that hold letters; the
    eight groups of amino acids need more colours than that palette has, so they were found by a search for the largest smallest distance
    between any two of them, as seen with normal vision and with simulated protan, deutan and tritan vision (Machado, Oliveira and
    Fernandes 2009, in CIE Lab), kept light enough for black letters and not loud. The smallest distance between two colours of a set, over
    the four kinds of vision, is now: bases 14.2 on a light page and 15.3 on a dark one (8.7 and 13.9 with the standard colours); residues
    19.3 and 13.4 (4.2 and 4.6: the standard dark ones are see-through, and are measured as they look over the page); trace letters 17.0 and 16.1
    (4.3 for both). The tests carry their own copy of the maths, and the numbers they pin came from a separate script in Python.
    **This is a measurement on a model of colour vision, not a trial with people, and the matrices were written from memory of the paper and not
    compared with it** (they pass the checks that a grey stays grey and that red and green come close for a person who cannot tell them). Colour is
    never the only sign: the letters are in the cells too. With the standard colours the letters of a cell keep the text colour as they always did,
    though in the standard dark theme three of their cells (the green A at 4.3 to 1, the orange G at 3.5 and the yellow proline at 3.8, as they look over
    the page) are below 4.5 to 1 with it, which was not changed (high contrast and the safe colours mend it).
  * *Mended on the way (the pictures at a zoom showed them).* The window's size and place, which egui reports at the zoom, were kept at the zoom and so
    shrank (or grew) at every start; they are kept in the window's own units now. The cells of the alignment were as wide as one letter rounded up
    and drifted from the letters of a long row (40 points in 200 letters at 80 %); the width is now shared out over a run. The rail (above). The screenshot at
    a zoom drew a picture the size of the zoom. A write of `appearance.txt` that fails is told once and tried again every five seconds; settings from the
    command line (`--appearance`, `--dark`) are put on the kept ones and a run given them keeps none (it used to write them over what was chosen); a value
    that is another option (`--appearance --dark`) is no longer taken as the value; the file is read only if it is a plain file, and no more than 4 KiB of
    it, and a byte-order mark does not hide its first line. **Two copies of the program** each write all six lines from their own memory, so that
    the last to change a setting wins the others' changes made meanwhile; that was not changed.
  * **Not done:** the safe setting does not reach the colours of the 3D structure, the map's features, the trees and the conservation bars; the web build
    keeps the choices in the browser's storage since 16 November 2026 (`docs/WEB.md`); a colour of the accent of your own choosing; following the system's own contrast and size settings; a size
    for the sequence text apart from the whole window's; **at 175 % and 200 % the top bar of a window 1,280 points wide is too narrow for its menus, which
    overlap the search box and the job lights** (the bar does not make room as it does for a narrower window; the window has to be about 1,400 points wide at
    175 %, 1,600 at 200 %); pages were drawn at every size without an error but only the sequence page, the alignment and the Appearance window were looked at;
    and a check on Linux and Windows (the keys there are Ctrl +, Ctrl − and Ctrl 0). The swatches of the accent have a ring for the pointer and for the keyboard
    (Tab, then Enter) and a test for that; the other controls of the window are egui's own and were not given a keyboard test.
* **Accessibility (M).** Keyboard-only use of every page, focus rings, screen-reader labels (egui supports AccessKit),
  minimum sizes of hit targets.
* **Languages (M).** English and Vietnamese first, with all text in one table.
* **Window placement (S).** Tool windows open under the top bar and remember where they were; docking panels later (L).
  **Built for the main window:** its size and position are kept (written a moment after it stops moving; not while it is maximised, minimised or
  full screen) and it starts there. **Built for the tool windows too (November 2026):** a tool window that was moved opens there the next run (`places.txt` in the settings folder, a line to a window:
  its corner in points and its title; written a moment after the last of them stops moving, as the main window's is; a window that is back where it opens by itself is forgotten; at most
  200 are kept and a line that is not understood spoils no other). A place is never put above or left of the rail and the top bar, and egui keeps a window on the screen, so one left on a
  monitor that has gone comes back onto the screen that is there (tested with a place 15,000 points away). **Not built:** the tool windows' *sizes*; the places are in points, so a
  run at another zoom opens them at the same points, which are another share of the screen; and, for the main window, only the numbers are checked to be sane (deleting `window.txt` in the settings folder resets it).

### 1.7 Getting results out
* **Publication figures (M).** SVG, PDF and high-resolution PNG of maps, alignments, trees, traces and structures, with
  fonts and colours as on screen.
  **Built in part, November 2026: *Save figure* on the map, the alignment, the trace, the tree (Trees and Phylogeny pages) and
  the structure** (`src/ui/figure.rs`, `src/ui/figure/svg.rs`). The view is drawn once more, alone and off the screen, by the
  very function that draws it in the window, and what it draws is saved as **SVG** (a line stays a line and text stays text, in
  stretches of one font and colour with the width they have on screen; the fonts are named, Inter and JetBrains Mono, with the
  usual ones after them) or as a **PNG** with four pixels to a point (fewer for a view so large that it would pass 24 million
  pixels). A figure is what the window shows (the same zoom, the same
  place it is scrolled to, the same turn of a structure, the colours of the light or the dark theme) without the buttons around
  it, and without the hint of what the pointer does; the map is given at least 900 by 700 points, so that the names of all its
  enzymes fit, which they do not in its small panel. Checked: each kind of shape against the SVG element expected, the button
  clicked on each page in tests that draw the window, and all five figures looked at, the SVGs in a browser and the PNGs as
  pictures. **Not built:** PDF (open the SVG in a browser or a drawing program and export it); the fonts are named, not carried
  in the file, so a computer without them draws the text in its nearest font, at the same width; pictures, soft shadows, what a view paints through a callback and a colour that is added to what lies under it are left
  out; a number that is no number is written as 0 and a character that a file may not hold is dropped, so that a figure always
  opens; a mesh whose corners differ in colour becomes triangles of one colour each; a
  figure is the whole room of its view, not cut to what is drawn in it; the figure of an alignment is the part in view, not the
  whole alignment; there is no choice of size or of pixels to a point; and the other views (Text Map, the digest's gel, the
  contig view, the domains and the motif logos) have no button yet.
* **Report builder (M).** One click makes an HTML or PDF report of a project: methods text written from the provenance,
  figures, tables, and the programs cited with versions.
  **Built in part, November 2026: File → *Save a Report of the Results…*** (`src/report.rs`, `src/ui/platform/report.rs`). One
  HTML file that stands by itself, with a part for each result that carries the record of how it was made (the alignment, the
  tree the Phylogeny page estimated, the species tree, and since 14 November 2026 the motifs with their logos, the protein
  domains drawn along their sequences, a predicted structure as it is turned on its page, the PhyloCNN estimates, and a run of
  the Pipelines page, by its fixed route or by a pipeline file, with what each step did and its tree): the **methods in sentences written from the record** (what made it, the
  day, where, how long, the settings, how many inputs), **how its method was checked** (the grade, what it was compared with
  and what was not checked, from the table of 2.4), **the figure** drawn by the view of the window (in the file itself, as SVG,
  in the light colours whatever the window wears), the settings as a table, the tree in Newick format and the scores of a
  species tree, **every input with its size and its checksum**, and at the end the programs with their versions. An alignment
  is drawn by its first 120 columns and 60 sequences, and the line under it says so; a tree of more than 296 tips has its
  tips closer than their names are tall, and the line under it says that too; a species tree whose text cannot be read is not
  drawn, and its part says why. A name from somebody's file that is written like markup (`<script>…`) is text in the page and
  in its figure, which a test takes through the whole way. Checked: the text and the page against what
  is expected (with text that tries to add to the page), the command on the examples of the three results in tests that draw the
  window, and a report of all three read in a browser. The six newer parts are checked by a test that reads the page for each
  part's heading, check, drawing and table (a table of more than 200 rows is cut, and its part says so); a run of the fixed
  route has no one check, and its part says where the checks of its steps are; the part of a pipeline file says that its check
  is of the running of the file, not of the methods of its steps. The line under a tree speaks of the scale bar only when one is
  drawn (a tree whose branches have no length has none), and the line under the domains says that there are eight colours, so
  that two families can share one. **A report of the six newer parts was read in a browser (15 November 2026), which showed a
  fault that no test had:** every figure named its first cut `cut0`, and a browser looks a name up in the whole page, so every
  figure after the first was cut to the room of the first figure's cut (of the domains only a sliver was left). The cuts of a
  figure are now named after the figure (a checksum of what it draws), a test reads a page of five figures for it, and the page
  was read again: the logo, the domains, the structure and both trees are as they are in the windows. **Not built:** PDF (print the page from a browser); a report of the
  results that have no record (reads, primers of the Primer3 window, digests, assemblies); the references to cite for each program and method
  (the report names the programs and their versions, not their papers); a choice of what goes in; the map, the traces and the
  structures of a project (a report is of results with a record, and there is no project yet that holds the rest).

### 1.8 Feeling like a real application
* **Installers (M).** macOS app bundle with signing and notarization, Windows installer, Linux AppImage or Flatpak;
  file associations (.gb, .fasta, .ab1, .pdb, .nwk) and "open with HelixLoom".
  **Built in part, November 2026: a Mac app** (`scripts/bundle-macos.sh`). The script makes a release build into
  `HelixLoom.app` (the program, its `Info.plist`, an icon) and a disk image of it to hand on. The icon is the program's own
  mark on a tile, drawn by the program at each size an icon holds (`helixloom icon FILE.png [SIZE]`). Made and tried here: the
  app is 25 MB and its disk image 10 MB; its description passes `plutil -lint`, its signature verifies, the program inside it
  draws its window to a picture, and opened as an app it starts, keeps its window's place in its settings folder, and closes.
  **Not built, and why:** it is signed *ad hoc*, not with a developer identity, and not notarized (both need an Apple Developer
  account: question 12), so on another Mac the first start needs a right-click and Open, **which was not tried on another
  Mac**; it declares no kinds of file and has no "Open With" (a file double-clicked in the Finder is handed to an app in a way
  the program does not take yet, so the app does not promise it; files are opened from the app, dropped on its window or named
  on the command line); there is no Windows installer and no Linux AppImage or Flatpak (neither system is here to try one on).
* **Updates and crash reports (M).** Check for a new version; opt-in crash reports; a "copy diagnostics" button.
  **Built in part, 20 November 2026: Help ▸ Check for Updates… and `helixloom update`** (`src/update.rs`, `docs/USING.md`). Question 11 (where releases are
  published) is answered by the releases page of the repository. **Nothing is asked of any server until the button is pressed** (no check at start,
  none in the background); then one request goes to `api.github.com` (made by `curl`, sending the program's name and nothing of the work). The window
  says whether a later release is out, shows its notes, and, if the person says so, downloads the file for this computer (the disk image of a Mac with
  Apple silicon, the zip of Windows) into the Downloads folder, checks it against the release's `SHA256SUMS` (a file that does not match is removed),
  and says what to do with it. **The program does not replace itself**: an ad hoc signed Mac app cannot be replaced behind Gatekeeper's back, and
  rewriting its own files cannot be tried on every system. The web version updates when its page is loaded again (the window says so and has a Reload
  button). Tested: SHA-256 against the published test values, the reading of GitHub's answer, the comparison of versions, the choice of file, the
  checksum and every refusal of a file or address, with a table in place of curl. **Tried against GitHub itself on the Mac** (`helixloom update` answered *up to date, the latest release is v0.2.0*); the
  download was not tried against a real file (no release later than `v0.2.0` exists), nor was the window in a browser or anything on Windows. **Not built:** a check at start (it would be the program asking about itself without
  being asked), an updater that installs, crash reports.
  **Built in part, November 2026: Help → *Copy Diagnostics*** (`src/diagnostics.rs`, `src/ui/diagnostics.rs`). It puts on the
  clipboard what a report of a problem needs: the version, the system, what the program was built with, the appearance, which
  settings folder is used (in words: the usual one, or one that was named), whether a workstation is connected and how many
  GPUs it has (**never its address or its token**), each program of the registry with where it is and its version, the package
  managers found, which model files are installed, and the note of the last fault. It is **made of things that cannot hold
  the work**: no sequence, no name of a file that was opened, and not what the program said when a fault stopped it (a
  message can carry a name from the work), only where in the code that was. The home folder is written `~` wherever it stands,
  a version is shown only when it looks like one, and the place of a fault is told without the folders of the computer the
  program was built on. **What it does show of the computer** is the paths of the other programs, which a folder you named
  yourself can be part of: read it before you paste it. **Nothing is sent anywhere.**
  **A fault leaves a note**: when the program, or one job of it, stops by a fault of its own (a panic, on any thread), the day,
  the version, the thread and the place in the code are written to `last-crash.txt` in the settings folder, for the
  diagnostics. The note is of the last fault only, with its day, and it is kept until the next one; the command line and the
  workstation server leave none.
  After a read for privacy, which found that the first version kept the fault's message (a sequence name could ride in it),
  hid the home folder only before a slash, and showed the settings folder's path: all three are changed as above.
  **Not built:** checking for a new version (it needs a place where releases are published: question 11; until then the
  program asks nothing of any server about itself); sending a report of a fault anywhere (it needs somewhere to send it and
  your yes; today nothing leaves the computer); an updater.

---

## Part 2. The bioinformatics pipelines

### 2.1 The pipeline engine
* ★ **Pipelines as data (L).** Today the Pipelines page runs a fixed route. Make a pipeline a list of steps with typed
  inputs and outputs and parameters (a graph underneath, a simple list or a node canvas on top), saved with a project,
  versioned, and runnable on one sample or a whole folder (batch). Why: this turns HelixLoom from a set of tools into a
  platform. Builds on `pipeline/` and the Snakemake/Nextflow export.
  **Built (November 2026), for the first target of the list above, *a folder of locus alignments → gene trees → species tree*:** a pipeline is a
  text file of typed steps (`align`, `tree`, `species`; each must take what the one before gives, checked when the file is read), run on one file or on every file of a
  folder (the per-file steps on all the cores, the trees then gathered for the rest), from the window (**Pipelines ▸ Pipeline file**) or with
  `helixloom pipeline run FILE FOLDER`. The format, the steps and the limits are in `docs/PIPELINES.md`. **Not built:** a graph with branches or joins
  (a pipeline is a list; trees are gathered by one fixed rule), steps for the other tools (primers, reads, motifs, structure), a node canvas, and saving a pipeline
  with a project. **Built, November 2026: the Snakemake and Nextflow export of a pipeline file** (`helixloom pipeline export`, `docs/PIPELINES.md`): for steps that other programs do (`engine=mafft`,
  `engine=iqtree`, `species`); the Nextflow file was run for real and gave the very species tree the engine gives; the Snakefile has not been run (Snakemake is not installed here).
  **Built (November 2026): a pipeline file runs on the workstation** as a job of its own kind (`flow`, `src/remote/flow.rs`): the files are read here and sent with the pipeline,
  the workstation runs the same function and keeps the answers of the steps in its own folder (`flow-cache` in `--data-dir`, so that a run sent again does only what changed),
  and the same report comes back and opens on the Pipelines page. With the choice on Automatic a run goes there when it passes the limits of the other jobs. Tried on a Mac with
  the server on the same computer, not across a network (`docs/PIPELINES.md`). **Not built:** a limit on the size of that cache, a cache shared between the workstation and this computer,
  `helixloom pipeline run` sending a run to a server from the command line, and a folder named on the workstation instead of sent.
* **Caching and resume (M).** A step whose inputs and parameters have not changed is not run again; a failed run resumes
  where it stopped.
  **Built:** with a cache folder every step's answer is kept under a key of the step, its settings, the installed program's version and the input,
  and a step whose key is there is not run again; a run that stopped is resumed by running it again (tests: change one locus and only its steps and
  the gathering step run again; a damaged cache file costs one step). A file that fails is left out of what is gathered and named in the report (exit code 1).
  **Not built:** a cache size limit or expiry (`clear-cache` forgets all), a key that follows the content of a step's program beyond its version line,
  and sharing a cache between computers.
* **Tool registry (M).** One place that knows each external program (MAFFT, IQ-TREE, primer3, fastp, minimap2, samtools,
  HMMER, MMseqs2, …): how to detect it, its version, how to install it (micromamba or a container), and its licence.
  **Built in part, November 2026: `helixloom tools`** lists the eight programs HelixLoom calls (MAFFT, IQ-TREE, HMMER, primer3, fastp, curl,
  Snakemake, Nextflow) with what each is used for, whether and where it is found and its version, its licence (as Homebrew's formulae list
  them that month; check a program's own before passing it on) and the commands that install it (`brew install …`, `conda install -c bioconda …`);
  the Workstation page shows the same programs for the server.
  **Built, November 2026: installing one on request** (`src/pipeline/install.rs`), from the window **Tools → *Programs HelixLoom Calls…*** (each program
  with what it is for, where it was found and in which version, its licence; *Install…* on one that is not found) and with
  `helixloom tools install NAME [--with brew|micromamba|mamba|conda] [--yes]` for the owner of a server. HelixLoom finds the package managers that are on
  the computer (Homebrew, then micromamba, mamba, conda) and runs the first that has the program's package with **the name the registry has for it and
  nothing a person typed**; no shell is started. **The command is shown first** (`brew install hmmer`), with where the program comes from and its
  licence, and nothing runs without a yes: a button in the window (which can also only copy the command), a typed `y` at the command line. What the
  package manager prints is shown as it comes, and when it ends well the programs are looked for again. conda and its kin install into the environment
  that is active and are used only when one is: HelixLoom makes no environment. A workstation never installs anything on a request.
  After a read for security (which found no way to change what is run and none to install without a yes): the notice names **the very program that
  will be started** (`with the brew at /opt/homebrew/bin/brew`) and the environment conda would install into; what the package manager prints is
  cleaned of the codes a terminal obeys (colours, a title, a write to the clipboard) and of the marks that turn the order of text; a line that never
  ends comes in pieces, with no letter cut in two; a bar of progress is shown by its last state; a process that the package manager leaves behind
  holding its output is not waited for, and one that goes on printing is listened to for a second only; the package manager is stopped and waited
  for if the run is left by a fault; and a search for the programs that was going when an install ended is not believed (another is made).
  An install cannot be stopped from the window once it runs. A code that a very long line is cut in the middle of shows as a few stray letters.
  Checked with stand-ins for the package managers (a script that prints what it was given, one that fails), with the window clicked in tests, and with
  the real command on a program that is installed (it says so) and on one that is not, answered no (nothing ran). **No real install was made through
  it**: HMMER was installed with `brew install hmmer` by hand before this was built. **Not built:** a container (Docker, Apptainer) or an environment
  of HelixLoom's own; pinning a version, or a pipeline that says which version it needs; removing or updating a program; minimap2, samtools, MMseqs2
  and the others, which no step calls yet; licence text checked against each program's own files; an install on the workstation from the app.
* **Run anywhere (L).** The same pipeline on this computer, the workstation, a Slurm cluster, or a cloud queue; export to
  Nextflow, Snakemake and CWL; import them. Builds on the Workstation and the `Compute` choice.
* **Reproducibility package (M).** Export a project as an RO-Crate: data, workflow, parameters, versions, checksums.

### 2.2 Data in and out
* ★ **Connectors (M).** Fetch from NCBI (Entrez, Datasets), UniProt, PDB, AlphaFold DB, ENA/SRA, Ensembl; keep a local
  cache; remember where each file came from (and its licence).
  **Built (November 2026) for NCBI nucleotide and protein, UniProt, the PDB and the AlphaFold DB, by accession:** File ▸ Fetch from a Database… keeps each record in
  the settings folder with a provenance file (database, accession, address, day, licence, checksum), opens it where it belongs (a nucleotide record as a
  document whose comment carries the provenance, proteins in the Alignment page's list, structures in the Structure page), and asks the database nothing for a record that
  is kept. The accessions are checked against each database's form before they are put in an address; the web is reached with `curl` over `https` only.
  **Not built:** searching (only accessions), NCBI Datasets, ENA/SRA and Ensembl, an NCBI key or e-mail, many records at once, and the browser build. **The licence lines were written
  from the databases' policy pages, not fetched: check them.** The tests use a table, and `curl` against a server on this computer. **Tried live once, on 11 November 2026 from one Mac** (the opt-in test
  `live_one_record_of_each_database_and_a_complex`, run with `cargo test --bin helixloom -- --ignored live`; about 490 KB in all): **UniProt** (P69905), **the PDB** (1UBQ and the
  complex 1VFB) and **the AlphaFold DB's file server** (P69905, model v6) answered and the files parsed; **NCBI's `efetch` redirected this computer to its abuse page, so its
  address is not confirmed** (NCBI asks for a `tool` and an `email` in a request; the program sends `tool` only, and never an e-mail address unless it is given one), and **the AlphaFold DB's entry service answered 403
  here**, so when it refuses, the files are tried by their names (versions 6, 5, 4: only 6 existed that day). NCBI protein was not tried.
* **Formats (M).** GenBank, EMBL, GFF3, FASTA, FASTQ, BAM/CRAM with index, VCF, SnapGene `.dna`, SBOL; large files read by
  index, not loaded whole.
* **Genome viewer (L).** A track-based viewer (like IGV) for reads, variants and annotations on a reference.

### 2.3 Analyses to add, by question
* ★ **Conserved motifs and domains** (you asked for these; planned in stages in 2.6): motif discovery (MEME/STREME), domain
  search against Pfam, InterPro and CDD with HMMER, per-column conservation scores, mapping conservation onto a 3D
  structure, protein language-model embeddings (ESM-2, ProtT5) for similarity and for the effect of a variant (ESM-1v),
  clustering of families (MMseqs2). Each is a *model group* in the models folder and a job kind on the Workstation (see
  `DEPLOY.md` section 10).
* ★ **Epitope and paratope** (you asked for these; planned in stages in 2.6): B-cell epitope prediction (BepiPred-3,
  DiscoTope-3), antibody numbering and region marking (ANARCI), paratope prediction (Parapred, Paragraph),
  antibody-antigen complex prediction (AlphaFold-Multimer, Boltz, Chai; check each licence), docking (HADDOCK, ClusPro),
  interface analysis (contacts, buried area), binding-affinity change for mutations, T-cell epitopes (NetMHCpan). Validate
  against the SAbDab and IEDB benchmark sets and show the accuracy on the page before the page claims anything.
* ★ **More ways to align, and probes for real-time PCR** (asked on 19 November 2026: "techniques like ClustalW, ClustalX, etc.", and probe
  design added to the primer tool for PCR and TaqMan real-time PCR). **Built, 19 November 2026, tried on a Mac only:**
  (a) a **progressive aligner after ClustalW** (`src/seq/clustal.rs`): pairwise distances (alignment, or shared words for big groups),
  a neighbour-joining or UPGMA guide tree rooted at its midpoint, sequence weights from the tree, profile-to-profile alignment with
  position-specific gap opening and free or costly end gaps, BLOSUM62 for proteins, and optional refinement by cutting the tree's edges
  and realigning; chosen on the Alignment page beside the built-in and MAFFT, sent to the workstation as a job, and kept with the
  session; (b) the **Clustal format** (`.aln`) read and written, with ClustalW's `* : .` marks of conservation (the groups of its
  documentation); (c) **TaqMan probe design** between the two primers of each pair (`src/seq/probe.rs`, in the Primer3 window), after
  Applied Biosystems' published rules, with probes of both strands and the probe's Tm set above the primers'. Checked: see *Help ▸ How
  Each Analysis Was Checked*. On simulated sequences the progressive aligner is about as accurate as the centre-star method (a little behind at
  low divergence, a little ahead at some higher ones), so what it adds is ClustalW's vocabulary (guide tree, weights, profile alignment) and
  refinement, not a gain in accuracy that has been shown. Three independent readers read the code (the aligner, the probe design, the format and
  its wiring); what they found was mended and given tests (the distances of the guide tree had counted a one-letter overlap as no distance; RNA
  was left out of them; a refinement that could not be aligned in one piece failed the whole alignment; a Clustal header of other programs was
  not known; the count at the end of a line was not looked at; features could be put on another sequence than the pairs were picked for).
  **The code was then broken on purpose in 99 places, one at a time: the tests caught 92 at the first run; the 7 that got through were looked at, 5
  of them got tests (the joining of the stretches between primers, a count that is not a number, the residues that have no partner in the
  distance, whether refinement can return a worse alignment, whether more than one pass of it helps) and were then caught; 2 stand,
  and are the same in everything the tests can see** (which way round the 3' end of a probe is scored against a primer, which the score
  of the whole stretch settles for probes of the sizes used; and the limit of ten passes, which the loop reaches only when a pass still
  improves). A first round that was stopped half way left one breakage in the source: found by checking every anchor, and put back by hand. **Not built, and said so:** a port of ClustalW's or ClustalX's code (this follows ClustalW's description,
  so it does not give ClustalW's alignments); ClustalX's colouring of conserved columns and its quality scores; T-Coffee, MUSCLE
  and other methods; the progressive method in pipeline files; minor-groove-binder probes; a probe-only mode for primers you already
  have; checks that a probe is unique in a genome or avoids a known variant; and no comparison with Primer Express or primer3's
  internal oligo. The gap costs were tuned on simulated sequences, which have more insertions and deletions than most real data.
* **Cloning design (M).** Gibson and Golden Gate planning with overhang fidelity, primer design with off-target checks,
  restriction-cloning planner, codon optimization with organism tables, in-silico gel and PCR, CRISPR guide design with
  off-targets.
* **Sequence verification (M).** Map Sanger or nanopore reads to the reference, call the consensus, list differences,
  and give a pass or fail with the evidence (the kind of report sequencing services return).
* **Phylogenetics (M).** IQ-TREE 3 with model finding and ultrafast bootstrap, partitioned models, dating (TreeTime,
  BEAST2), ancestral reconstruction, tree comparison and tanglegrams.
* ★ **Species trees, the rest of the ASTER family** (the ASTRAL core is built, `docs/SPECIES.md`). In the order that pays
  most: **a batch pipeline** from a folder of locus alignments to one gene tree each (IQ-TREE or the built-in methods, on the
  workstation) to the species tree (M, needs 2.1; **built, November 2026**: see 2.1, on this computer and as a workstation job, tried on a Mac only); **weighted ASTRAL** (gene trees with branch lengths or supports as
  weights, which is more accurate when the gene trees are poorly resolved; M; **built, November 2026, as the paper defines the
  weights, not yet compared with wastral, and without the weighted branch support**: `docs/SPECIES.md`); **multi-copy gene families** (ASTRAL-Pro3,
  with the orthology and paralogy that families bring; L; **built, November 2026, as the score of ASTRAL-Pro (2020), not
  ASTRAL-Pro 2's or 3's algorithms, without branch supports, and not compared with its program**: `docs/SPECIES.md`); **bootstrap over genes and sites** (**built, November 2026**: `docs/SPECIES.md`; not compared with ASTRAL's; **the replicate trees, their greedy consensus and the bootstrap with the weights are built too** (November 2026)) and a *polytomy test* (M; **built, November 2026**, as the paper of Sayyari and Mirarab (2018) states it, `-t 10` / `--polytomy`: `docs/SPECIES.md`; not compared with ASTRAL's);
  **CASTER and WASTER**, species trees straight from whole-genome alignments and from reads (L, and they need the
  alignment and read pages of 2.3 to be solid first); **branch lengths in substitutions per site** (CASTLES-II; M);
  **D\*** and the test for hybridization (M). The paper for each is in `docs/SPECIES.md`. ASTER is under the AGPL-3.0 and
  HelixLoom is MIT (question 6), so each is written from its paper, as the ASTRAL core was.
* **Reads and variants (L).** fastp and MultiQC reports, minimap2 or bwa-mem2 mapping, samtools, bcftools or DeepVariant
  calling, assembly of long reads (Flye) and polishing (Medaka), taxonomic profiling (Kraken2).
* **Structure (M).** AlphaFold DB lookup, ColabFold, structure alignment (TM-align, Foldseek), pocket detection.

### 2.7 Genome assembly, mapping, typing and annotation tools (asked on 19 November 2026)

The owner asked for the functions of ten repositories to be rebuilt in HelixLoom: Unicycler, the Galaxy workflow for plasmid collapse
with Unicycler, `hybridassembly`, the multiplex-MinION protocol of the Unicycler author, SPAdes, NanoSeq, minimap2, KAT, mlst and gffread.
Three readers read their pages (README, wiki, papers, format specifications; nothing was cloned or downloaded) and reported.

**The rule for each, from its licence** (the same as for ASTER, 2.3): **minimap2 and gffread are MIT**: reading is allowed, and what is
written is still from the paper and the documentation. **Unicycler, SPAdes, KAT, mlst, NanoSeq and the Galaxy workflow are GPL (2 or 3)**:
HelixLoom is MIT, so nothing of their source is read for copying, and each is rebuilt from its paper and manual only. **The MinION protocol
repository has no licence**: nothing of its code or text is used, only the method as its published paper describes it. **`hybridassembly`
is MIT** but is only a Nextflow pipeline of other programs. **PubMLST's data may not be bundled** (its terms of 13 May 2026: data submitted
since 1 January 2025 is for non-commercial use and may not be redistributed in software): mlst typing takes a scheme folder that the person
downloaded.

**The order, by what each needs and what each is worth to a plasmid and small-genome user** (each stage ends with its tests, a round of
deliberate breakages, an honest list of what is not done, and no claim of equivalence with the original until it has been run against it,
which needs the original's test data or its program: a download that waits for a yes):

1. **gffread: GFF3 and GTF** (MIT): **built, 19 November 2026** (`docs/GFF.md`, `helixloom gff`, Features ▸ Import Features from GFF/GTF…). The code was broken on purpose in 25 places: all 25 were caught, 7 of them only after tests were written for them (the phase of a minus strand
   CDS, a line of eight columns, the kind of a gene with CDS under it, an intron of one base, the AT-AC pair, a partial last codon, the width of a FASTA line).
2. **KAT: k-mer counting, spectra, comparison** (GPL, from the paper) — **built, 20 November 2026** (`docs/KAT.md`, `helixloom kat hist|gcp|comp|sect`;
   a count checked against a plain string counter, a simulated genome's spectrum and size; not compared with KAT; `filter`, `blob`, the plots and
   counting on disk not done; with the update check, the code was broken on purpose in 30 places and all 30 were caught, 5 of them only after tests were
   written for them: G counted among the G and C, the plateau at the bottom of the valley, the start of the valley, the exact size and the height a half peak must have). The plan as written: a canonical k-mer counter, `hist`, `gcp`, `comp`, `sect`, a peak finder
   and a genome-size estimate (the textbook one: k-mers above the error valley over the homozygous peak, which is not KAT's own fit, whose
   method the pages do not give). The k-mer counter is the base of stage 5.
3. **minimap2: minimizers, chaining, PAF** (MIT, from the paper) — **built, 21 November 2026** (`docs/MAP.md`, `helixloom map`; checked by reads simulated from a known
   sequence and by hand-worked cases; the code was broken on purpose in 47 places and all 47 were caught, 8 of them only after tests were written for
   them (how the share of frequent minimizers rounds, the exact limits of a gap and of the bandwidth, a chain that runs into one taken before it, an overlap of
   exactly half, the union of the bases the anchors cover, the divergence, and an anchor that follows the last on the read with a jump on the target);
   **not compared with minimap2**; the ends are not extended, and splice, short reads, SAM and the other presets are not done). The plan as written: (w,k)-minimizers (k 15, w 10), an index, anchors, the chaining of the paper
   (gap cost 0.01·w̄·|l| + 0.5·log2|l|, 50 steps back), primary and secondary chains, mapping quality, PAF; the overlap mode (`ava-ont`);
   base-level alignment from HelixLoom's banded aligner. Out of reach for now: the SSE kernel, splice mode, paired short reads, parity with
   minimap2's output. Tested by reads simulated from a known sequence (the original's own examples need a download).
4. **mlst typing** (GPL, from the README) — **built, 22 November 2026** (`docs/MLST.md`, `helixloom mlst`; checked on made-up schemes and genomes, **not compared with mlst**; 61 deliberate breakages, all but one boundary caught after tests were added; the thresholds and the tie rule are stated guesses; no BLAST, no compressed input). The plan as written: a scheme folder (profile table and one FASTA of alleles per locus), contigs scanned on both
   strands with HelixLoom's aligner, calls `n`, `~n`, `n?`, `-` and `n,m`, the score (+90/N exact, +63/N novel, +18/N partial, +10 for a known
   sequence type), the sequence type, TSV, CSV and JSON, novel alleles as FASTA with MD5 names. Which scheme is chosen automatically, and
   how ties go, is not documented: a stated guess.
5. **A de Bruijn assembler, SPAdes-lite** (GPL, from the papers) — **built, 23 November 2026** (`docs/DEBRUIJN.md`, `helixloom dbg`; only parts of the two papers could be read in full, the rest second-hand; checked on simulated reads of random sequences; **never run on real reads or compared with SPAdes**; 94 deliberate breakages, 80 caught at first and the rest closed with tests or a fix except four that change nothing; no error correction, no paired reads, k to 63). The plan as written: k-mer counting, a graph at one k or several (21, 33, 55), tip clipping,
   bubble removal, a depth filter, unitigs as GFA and FASTA. Useful for small, plasmid-scale reads only; BayesHammer-quality error
   correction, repeat resolution with read pairs and genome-scale speed are out of reach.
6. **Hybrid bridging, Unicycler-lite** (GPL, from the paper): read a GFA, give each contig a multiplicity from its depth and links, bridge
   single-copy contigs with long reads (stage 3), merge, detect and close circular replicons, rotate them to a start gene. The paper gives
   no bridge-scoring formula, so the scoring is this program's own and says so. No polishing (Racon, Pilon and Medaka stay external).
7. **Pipelines** that say the routes of the Galaxy plasmid-collapse workflow (reads → trim → filter → short-read assembly → hybrid assembly
   → evaluate; its multiplicity correction is a human step), of NanoSeq (Canu or Flye → concatemer check → circularisation check → rotate →
   Medaka; the 30 bp overlap and 95 % identity rules are in its README), and of the multiplex-MinION protocol (the paper's method), as pipeline
   files that call the external programs where they are installed and, as stages 2 to 6 are built, the built-in ones.

**Not planned**: basecalling, Canu, Flye, Medaka, Racon, Pilon, SPAdes' own error correction and scaffolding, and anything that needs the
GPL source. **Stage 1 (gffread) has not been run against `gffread`**: its examples folder is small, and fetching it needs a yes.

### 2.8 Tree drawing and editing after FigTree and TreeGraph 2, and primers from an alignment (asked on 20 November 2026)

The owner asked to improve the tree visualisation with what FigTree and TreeGraph 2 do, and for a button that picks primers and probes from
an alignment. Both tree programs are GPL, so **no source was read**: two agents read their documentation, tutorials, change logs and the
TreeGraph 2 paper and returned feature lists in their own words; the code was written from those lists (`docs/TREES.md`).

* **Built, November 2026 (`docs/TREES.md`):** NEXUS files with `translate` tables and BEAST-style `[&key=value]` annotations (`seq::nexus`, and
  the Newick reader keeps the annotations); the radial (equal-angle) layout; collapsing a clade to a triangle; colour, width, labels, shapes
  and bars by annotations with a legend (`seq::treestyle`); an axis in distance, time before the tips or dates; highlight, find, rotate,
  order, root on a branch, take out weak branches, rename tips, cut out a clade (`seq::treeedit`); support values from other trees: the
  share of replicates with each branch, its strongest rival, or the value on the same branch of another tree.
* **Built, 20 November 2026: align DNA by its translation (`docs/CODON_ALIGNMENT.md`)**, asked by the owner after the first look at the
  primers window: the Alignment page's *Align by the translation* translates every sequence (reading frame 1 to 3, any genetic code the
  program has), has the proteins aligned (MAFFT when chosen, else the progressive method) and puts the codons back (`seq::codonalign`). Not
  done: choosing a frame for each sequence, finding the frame, frame shifts and introns, a protein view of the result, the command line.
* **Built: primers and probes from an alignment (`docs/ALIGNMENT_PRIMERS.md`):** the button *Pick primers and probes…* on the Alignment
  page; oligos only on columns that enough sequences share, each laid against every sequence, a ranking by the sequences covered, the oligos
  marked on the alignment, the consensus opened with them as features, a table to copy or save.
* **Not run against the originals:** FigTree and TreeGraph 2 were not run, so nothing says the drawing or the numbers are the same; the
  rules the documentation leaves open (how a rooting shares a branch's length, how a missing branch is valued, how the radial layout
  is drawn) are this program's choices and are said so in the text. **Not looked at on a screen by a person**: the tests draw the page
  without a window and read its text and shapes.
* **Not done:** see the end of `docs/TREES.md` (fish-eye zoom, polar controls, taxon colours in NEXUS, tip-annotation import with
  ancestral states, pie charts, separate labels with their own formats, PhyloXML, NeXML, XTG, cut and paste of nodes, a command line) and
  the limits of `docs/ALIGNMENT_PRIMERS.md` (no degenerate oligos, no genome check). The `helixloom` command line has nothing for either.

### 2.4 Trust
* ★ **A benchmark page (M).** For every analysis, say how it was checked and how accurate it is (the golden tests already
  compare Tm, PCR, tree building, DSSP, ABIF reading and the calculators with the reference programs; show that to the user).
  **Built, November 2026: Help → *How Each Analysis Was Checked*, and `helixloom checks [WORDS…]`** (`src/checks.rs`, the table in
  `src/checks/table.rs`, the page in `src/ui/checks.rs`). Forty analyses, each with a grade (*reference program*: compared with the
  program it follows; *second implementation*: a separate, plain implementation or the definition worked out the slow way; *known
  answers*: hand-worked or planted cases; *not checked*), what it was compared with, how close it came with the tolerance the test
  uses, **what was not checked**, and the tests by name. Today: 16 reference program, 7 second implementation, 15 known answers,
  2 not checked (a folding program's structure and MAFFT's or IQ-TREE's answer, which are those programs' own); the domain
  search, which had never met a real HMMER when the page was made, has since been run against HMMER 3.4 and Pfam 39.0 and
  is among the known answers. The page has a filter, counts of what is shown and *Copy as text*; the *How this was
  made* section of a result links to it. **What keeps it honest:** a test reads the source and fails when the table names a test
  that is not there, is not marked as a test, is left out of the usual run (`#[ignore]`, also through `cfg_attr`), or is compiled
  only in some builds (`#[cfg(…)]` on the test or on its file), so the page cannot claim a check that nobody makes. It cannot see
  a module that its *parent* compiles only in some builds: the tests of the GPUs are so, and the table says it of them. What the
  table *says* about a test (the tolerance, the number of cases) is written by hand from reading it and no test checks that.
  **It was read once against the tests it names** by two independent readers, who found twelve statements to correct, all
  corrected: a bootstrap test described as asking for full support asks for more than 60 %; the support of split genes is only
  checked to lie in a range; the exact search of ASTRAL was compared on six of the twelve data sets; four lists of what is not
  ported or not compared were incomplete; a count was 3,800 lines and not 3,900; and two tests had names that said more than they
  assert, which are renamed (`bootstrap_supports_a_clear_split…`, `…within_three_parts_in_a_thousand`). The stricter test then
  found one more by itself: the tests that run the stand-in for HMMER are compiled on Unix only, so they are no longer listed. **Not built:** the numbers are not measured when the page opens (it shows what the tests assert,
  not a run on your computer); no hold-out benchmark sets (2.6, step 7: SAbDab, IEDB, SKEMPI are downloads); the link from a
  result opens the whole page, not the one check of that result's method; and the three comparisons that are run by hand (real
  trace files against Biopython and the original, the Nextflow file run for real) are named in the text only.
* **Model cards and licences (S).** Every model in the models folder lists its source, its licence, the data it was trained
  on, and what it must not be used for.
* **Reference checks in CI (M).** Run the golden comparisons on every change, on macOS, Linux and Windows.

### 2.5 Compute

* **GPU compute with candle: what is done and what comes next (M each unless said).** Done: the pool of cards with speed weights,
  exact distances and bootstrap on the cards, the likelihood on the cards, neighbour-joining (PyTorch), routing of tree jobs by size, `helixloom gpu-check [--bench]`.
  Next, in this order: (1) run `gpu-check --bench` on the Ubuntu machine and set the sizes it suggests (the program's own are from a Mac); (2) **keep partial
  likelihoods and re-root at the branch being fitted**: **done for the lengths of branches, on the processor and with PyTorch**. The search tried the whole tree about 20 times for each branch length; now what lies beyond both ends of the branch is pruned once and each length costs 16 products for a pattern (`EdgeScan` in `seq/likelihood.rs`, and
  `edge_begin`, `edge_eval` and `edge_end` in `accel/torch/worker.py`, with `EdgeScan` in `seq/likelihood/torch.rs`, where the two ends of a branch stay on the card). The likelihood is that of the whole tree to 1e-9 on the processor, for every branch of random trees with gaps and ambiguity codes and every length from 1e-6 to 5,
  and to 5e-6 on PyTorch workers (1, 2 and 3 of them, and from `par_map` threads); the search for a topology (`improve_topology`) is **4.7, 7.1 and 10.1 times faster at 10, 20 and 40 tips** (600 to 1,500 columns, on a loaded M1 Max), and setting 12 branches on a PyTorch worker (on the CPU device, whose cost for each operation is like a card's launch) **3.8, 9.7 and 14.4 times**, which is the 10 to 20 times that was expected. The search can take another path, though: a greedy search that takes the first change that helps follows
  differences of 1e-4 in a gain, so from a bad start (a chain of tips) the final likelihood was the same in 21 of 40 random cases, higher in 11 and lower in 8 (on average 12.8 higher), and from a tree one swap off the truth both ways end on the same tree. What is left (M): candle keeps no partial likelihoods (its engine says so and the whole tree is tried, as before); the partial vectors of a branch could serve its neighbours (a real re-rooting: one pruning for every branch instead of two); the Bayesian sampler's branch moves still try the whole tree; (3) with PyTorch on CUDA, replay each tree's whole evaluation as a CUDA graph (the tree is the same for the ~20 trial lengths of a branch, only the matrices change), which removes
  nearly all the launching, and `torch.compile` for the fused tail (M, only worth doing once `gpu-check --bench` has been read on the CUDA cards); with candle, a
  [4, patterns] layout for the four-wide reductions and, where TF32 is off (candle's default), a matrix product for the messages (S);
  (4) thresholds in *patterns* rather than cells once CUDA is measured (S); (5) **retry the share of a card that fails on another card: done for PyTorch** (a worker that dies has the patterns it was given done by the next worker, which holds every column, and the result's report says "a GPU failed (why) and another GPU did its share"; tested by killing one of two workers while four threads of `par_map` ask for the likelihood, and the other answers for all of them; the same retry for a likelihood cut among the workers, which happens from 200,000 patterns for each worker, is in the code and not tried by a test, whose alignments are smaller); not done: **spread the model fits by card speed instead of by turns** (S: with four fits and two cards of unequal speed the gain is about a tenth, and more for the chains and resamples of the sampler and the bootstrap), the same retry for candle, and for the shares of distances (S); (6) **alignment**: done: the choice of a centre is one pass on the processor now (an independent review showed the sum over all pairs of shared 6-mers is, for each 6-mer, its number of holders minus one: no GPU
  was needed, and an earlier version of this plan had put it on the GPUs and called it 61 times faster than a slow baseline), and neighbour-joining is on the GPUs with PyTorch. The pairwise dynamic programming (`seq/align.rs`) was 10.4 ns for each cell on one core and the speed-up was on the processor: (a) **done**. The three choices at
  each cell are branchless now (M against X, then the winner against Y: two conditional moves each, the same ties and f32 sums in the same order), and `align_many` runs eight pairs in lock step, a lane for each pair,
  which the compiler turns into 128-bit vector instructions in a release build (no `unsafe`, no `std::simd`, no new dependency; pairs of unequal length are sorted into groups and padded, and nothing in a padded column is ever read).
  It is bit-identical to the old `align`, which the test module keeps, unchanged, as its oracle: about 850 pairs (lengths 0, 1, 2 and up, DNA, two letters, protein letters, mixed case, related pairs with changes and indels,
  repeats and poly-A, free ends or not), the whole alignment compared and the score to the last bit, and a whole centre-star alignment against the pair-by-pair one (changing a tie, an end-cell rule or a trace code fails it).
  Measured on an M1 Max, one core, for each cell, 500 to 4,000 bases, flat in the size and the same for unrelated and related pairs: **before 10.3 ns (10.9 in a release build); one pair now 2.9 ns (3.5 times); eight in lock step 1.33 ns in a
  release build (8 times)**, and 3.6 ns in the dev and test profile (`opt-level = 1`, where the compiler makes no vector instructions at all, so there the lanes are scalar code and the one-pair path is the faster by a fifth;
  the dev profile still gets 2.9 times from `align_many`; `opt-level = 2` for this crate would make it 1.9 ns, `opt-level = 3` the release figure). A centre star of 200 sequences of 1,500 bases (`center_star`, all cores): **0.13 s in a release build, 0.70 s before; 0.32 s in the test profile, 0.65 s before**,
  of which the alignments are now 0.10 s (0.28 s) and the choice of a centre 0.02 s (back to back with the old code on a machine loaded by other work: the 199 alignments 0.87 s before and 0.12 s after in a release build, 0.74 s and 0.29 s in the test profile, with identical alignments). **What is left of (a), each small (S) and none urgent:** (a1) the one byte of trace for each cell and pair costs about half of the lock step (0.62 ns without it); a flag-bit coding of it was 12% faster in a release build and
  much slower for one pair, so it was not taken; (a2) a group's trace must fit in 64 MB, so sequences over about 2,800 bases go one pair at a time at 2.9 ns (a limit of 128 MB would take them to 4,000 bases, at 1 GB of memory over 10 threads: your choice);
  (a3) `par_map` shares the 25 groups of a job among the cores in equal parts, so the slowest core (an efficiency core on this Mac) sets the time (0.10 s measured against about 0.06 s), which a queue that each core takes the next group from would fix, and which a workstation with equal cores does not need;
  (a4) **done (questions 9 and 10 were answered yes)**: the dev and test profile now has `opt-level = 2` for this crate too (it made no vector instructions at 1; 1.9 ns a cell instead of 3.6, at the price of slower compiles), and the scoring is in whole numbers (+10, −1, −10, −2), which makes every implementation equal by construction and changes alignments only where rounding broke a tie. (b) **the eager version is done, opt-in; the fused one is not**: `accel/torch/align.rs` and `align_batch` in `worker.py` fill the tables of many pairs along their anti-diagonals on the PyTorch workers (shared among the GPUs of the job by speed), with the processor's 32-bit sums and ties, so every alignment is **the same to the last bit** (the old aligner is the oracle: about 12 sets of up to 9 pairs on pools of 1, 2 and 3 workers, free ends or not, and a line-by-line numpy port of the kernel on 88 batches, on the CPU device and on the Mac's GPU, which has no integer shifts or `index_copy_` in this old PyTorch, so the code goes round them), a failed worker leaves the pairs to the processor, and what the workers do not take (an empty sequence, a pair too large) the processor does and the results come back in order. It is used **only on request** (`HELIXLOOM_ACCEL=force`, or `HELIXLOOM_ACCEL_ALIGN_CELLS` set; infinite by default) because it is **slower than the processor wherever it was measured**: 52 to 74 ns a cell and pair on a CPU-device worker against about 1 on the processor's cores (61 to 72 times), and 292 to 949 ns on this Mac's GPU (413 to 1,008 times; a loaded machine, 32 to 64 pairs of 300 to 500), each diagonal being about sixty small operations. `helixloom gpu-check --bench` now times it against the processor and suggests the variable only if the cards win by 1.5 times. What is left (L, needs a CUDA machine): one diagonal step under `torch.compile` or as a CUDA graph, which fuses the adds, compares and selects, and is the one thing that could beat the cores (only if a job really reaches 10^11 cells: 10,000 sequences of 3,000 bases), checked bit for bit by the same tests on the cards. (c) **done** (`seq/tree/joins.rs`): the processor's neighbour-joining makes the same joins to the last bit as before (a copy of the old code is the test's oracle: about 190 matrices of 0 to 130 clusters, and of 450 and 801, with ties, zeros of either sign, matrices that are not symmetric, infinities and not-a-numbers, the work shared in 2, 3, 4 and 7 parts; changing a tie, the order of a sum or a step of the packing fails it), and faster: the table is one flat array that a join does not shift (the clusters that are left are a list of places, and when a quarter of the places are empty the table is packed into a smaller one), the row sums are made four rows at a time (each still added in the order of the list), and from 800 clusters on the sums and the search are shared among the cores through a work queue (`par_queue` in `seq/par.rs`: each thread takes the next run of rows, so a slower core does fewer). Measured on an M1 Max, on a machine that other work was loading, for 1,200, 2,400 and 4,000 clusters: 1.27 s, 10.9 s and 50 s before; 0.71 s, 5.9 s and 27 s on one core (1.8 times); 0.46 s, 2.3 s and 6.5 to 8.8 s on all cores. What is left: each join starts the threads twice (about 0.4 ms), so the cores give 4 to 8 times at 4,000 and little below 1,200, which threads that stay and wait at a barrier would remove (M); the work is still n cubed, and RapidNJ's sorted rows would cut it and can be made to take the same pair, at n squared more memory (M); `JOIN_CELLS` in `accel/mod.rs`, from which neighbour-joining goes to the GPUs, was a guess against the old processor time (twenty seconds of one core at about 2,700 sequences, which is now 8 s on one core and about 3 s on all), and `gpu-check --bench` decides it; (7) other libraries than candle (burn, tch-rs, dfdx) only if candle's CUDA numbers disappoint (L).
* **From one workstation to a server for a lab (L).** Accounts and sign-in (OIDC), TLS, quotas, job queues with
  priorities, GPU scheduling, storage on S3 or NFS, audit log, container images, a REST API with an OpenAPI description and a
  Python package so that notebooks can drive the same pipelines. Builds on `src/remote/`.

### 2.6 A plan for conserved motifs and domains, then epitopes and paratopes

Not built yet. The order is by how little each step needs (no model, then a program, then a trained model on a GPU) and by
how soon it can be checked against something known. Each step is a job kind of the workstation (`DEPLOY.md` section 10), so
each gets the choice of where it runs, the queue, the GPU and the way back for free, and each paints its result on the
pages that exist (the alignment viewer, the 3D structure viewer, the tables).

| Step | What | How | Licence note (check before building) | Size |
|------|------|-----|------|------|
| 1 | **Conservation** per alignment column (Shannon entropy, Jensen-Shannon divergence of Capra and Singh 2007), conserved blocks, and the same colours on a 3D structure through the residue numbering. **Built** (Alignment page: *Conservation*; Structure page: *Colour by conservation*): sequence weights of Henikoff and Henikoff 1994, the BLOSUM62 background, a gap penalty, a window of 3 columns on each side, a score only for columns that are at most 30 % gaps; blocks are runs of 4 or more columns in the top quarter of the scores; a structure's chains are lined up with the best-matching row by the program's pairwise aligner, so missing residues and tags are no obstacle | native Rust; `seq/conservation.rs`, `seq/fold.rs`, the alignment and structure pages | none: published formulas | S to M, done |
| 2 | **Domains**: search a protein against Pfam, with the domain architecture drawn along the sequence. **Built, November 2026, and run against a real HMMER 3.4 and Pfam 39.0** (Tools → Protein Domains, the `domains` job of the workstation, `docs/DOMAINS.md`): `hmmscan` is run as a program with `--domtblout`, its table is read and drawn along each sequence, here or on the workstation. In every run it is checked with tables made to the documented format and a stand-in program; by hand, where HMMER and Pfam are installed, the protein of HMMER's tutorial gives the ten domains `hmmscan` itself writes, and the four proteins of the bundled structures are found in the families they are known for, in about two seconds. Not built: Pfam's clan rule for overlapping families, InterProScan, CDD | HMMER 3 (`hmmscan`) run as a program, Pfam as a model group in the models folder; InterProScan and CDD later | HMMER is BSD-3, Pfam is CC0, CDD is public domain; InterProScan is Apache-2.0 but some of the databases it calls are not free to pass on | M |
| 3 | **Motifs**: scan with position weight matrices, discover motifs by expectation-maximization, compare them to known ones. **The scan, the discovery and the comparison with each other are built** (Tools → Motif Finder), and so are the command `helixloom motif` and the workstation job (`docs/MOTIFS.md`); known motifs and protein motifs are not | native Rust for scanning and a simple discovery; MEME, STREME, FIMO and Tomtom when installed | the MEME Suite licence restricts commercial use, so call it, never include it | M |
| 4 | **Protein language models** on the workstation: per-residue embeddings (ESM-2), the effect of a variant (ESM-1v style masked marginals), similarity search and family clustering (MMseqs2) | a program of the workstation's own, named by the server owner and never by a request, exactly as `--fold-command` works for folding; weights in the models folder | ESM-2 and ESM-1v are MIT; ProtT5 and the newer ESM models have their own terms; MMseqs2 and Foldseek are GPL-3.0, so only ever run as separate programs | M to L |
| 5 | **Antibody regions and the paratope**: number an antibody (ANARCI), mark the framework and the CDRs, predict which residues touch the antigen (Parapred, Paragraph), predict B-cell epitopes on an antigen (BepiPred-3, DiscoTope-3) and T-cell epitopes (NetMHCpan, MHCflurry) | numbering is native or ANARCI as a program; the predictors as programs of the workstation; scores per residue painted on the structure | ANARCI is BSD-3, MHCflurry is Apache-2.0; BepiPred, DiscoTope and NetMHCpan are free for academic use only (DTU); Parapred and Paragraph: read the repository | L |
| 6 | **The interaction itself**: predict the antibody-antigen complex on the GPUs (Boltz-2, Chai-1, AlphaFold-Multimer; AlphaFold 3 only under its non-commercial terms), then **analyse the interface natively** (residue contacts, buried surface by the Shrake-Rupley method, hydrogen bonds), dock with HADDOCK3 when a model is poor, and score a prediction against a known complex with DockQ | job kind on the workstation with the GPU; the analysis is native Rust on the structure page. **The native analysis is built** (Structure page: *Interface*; `seq/interface.rs`): choose which chains are on each side, and it measures the solvent-accessible surface (Shrake-Rupley, 1.4 Å probe, Bondi radii, hydrogens left out) of each group alone and together, so every residue's buried area, the interface area (half of what is buried), the residues in contact (heavy atoms within 4.5 Å) and the **candidate** hydrogen bonds (a donor and an acceptor by atom name within 3.5 Å, no angle: some are not real bonds); the structure and the sequence strip are coloured by buried area, a table can be clicked to pick a residue in 3D and copied for a spreadsheet. Checked against an independent implementation on ubiquitin cut in two, and against the exact areas of one and two atoms. **Known limits:** an ion has the radius of carbon (wrong for it), a ligand or sugar is measured and listed but cannot be picked in 3D (the viewer draws protein chains), the hydrogen-bond rule is by atom names within 3.5 Å with no angle, the first model and first conformer only are read, and in the browser build the measuring runs on the page's own thread and stalls it on a large complex. **Not built:** the predicting of the complex, docking, DockQ, a hydrogen-bond angle, salt bridges, side chains' relative accessibility. **Tried on a real complex** (PDB 1VFB: the antibody D1.3 Fv, chains A and B, with hen egg-white lysozyme, chain C; fetched live; test
  `an_antibody_and_its_antigen_from_the_pdb`, which runs with the others on the file `assets/structure/1vfb.pdb`; the *Antibody example* button of the Structure page loads it): 1372 Å² buried in all (686 Å² on each side), 39 residue contacts and 18 hydrogen-bond candidates; **all 19 residues of the paratope are in CDRs**
  (light His30, Tyr32, Tyr49, Tyr50, Thr53, Trp92, Ser93; heavy Thr30, Gly31, Tyr32, Trp52, Asp54, Arg99, Asp100, Tyr101, Arg102 and others), the epitope is the two lysozyme segments 18-27 and 116-129 that
  are known for this complex, and Gln121 is the antigen residue that buries most. That the interface falls where it is known to fall is a check on the method, not a measure of how exact the areas are; the usual size of such interfaces is
  some hundreds of square angstroms on each side (the literature figure for this complex was not looked up) | Boltz is MIT, AlphaFold-Multimer is Apache-2.0 with CC BY 4.0 parameters, HADDOCK3 is Apache-2.0, DockQ is MIT; Chai-1, AlphaFold 3 and Rosetta/FoldX have terms to read | L |
| 7 | **The benchmark page** (2.4) for all of the above | hold-out sets: SAbDab (antibody-antigen structures), IEDB (epitopes), SKEMPI (binding changes on mutation) | SAbDab is CC BY 4.0; IEDB and SKEMPI are free to use | M |

Rules that go with the plan:

* **A model does not ship until its published numbers are reproduced here**, on data it was not trained on (split by
  sequence identity and by date, so that near copies of the training set cannot flatter it), with the metrics the field uses:
  AUC-ROC, area under the precision-recall curve and MCC for residues, DockQ and the share of correct complexes for docking.
  The page shows these numbers next to every result.
* **Say what is still hard.** Predicting which antibody binds which antigen, blind, is right in only a minority of cases even
  with the best current models, and their confidence scores are not calibrated for it: show the interface confidence (ipTM,
  pLDDT), the benchmark rate for that kind of case, and a plain "treat as a hypothesis to test at the bench".
* **Licences first**: before a model is added, its card in the models folder says where it comes from, under what terms, what
  it was trained on and what it must not be used for (2.4). A model that is free for academics only is run by the person
  who has the licence, on their own workstation, and is never part of the repository or of a download.
* **Nothing from a request is a program**: a model is named by the server owner's options, as folding is (`DEPLOY.md`).

---

## Part 3. A web version alongside the computer version

One codebase, two targets: the same Rust program compiled to WebAssembly runs in the browser (`docs/WEB.md`). Phases:

1. **W1 (this release)**: the program builds for the web; open files by dropping or choosing them, download results;
   small jobs run in the page.
2. **W2 (M), done**: the browser connects to a Workstation over `fetch` with the same token, and the server also serves the
   web files (`helixloom serve --web-dir dist`), so heavy jobs run on the GPUs and results come back to the browser. A page
   hosted elsewhere is allowed by name (`--allow-origin`). Tried in a real browser against a real server: connect, an
   alignment, a Bayesian tree and a whole pipeline (with the workstation's MAFFT and IQ-TREE) on the workstation, the
   wrong-token and server-down messages, a page of another origin (allowed and not allowed) and removing a job; `docs/WEB.md` lists what was not tried. Still to do: a job tray with progress
   and cancel (1.3). The address and the token are kept for the tab since 17 November 2026 (W3).
3. **W3 (M)**: projects kept in the browser (IndexedDB) and shareable links to read-only results.
   **Built in part, 16 November 2026: the page keeps what the desktop app keeps in its settings folder** (the open sequences,
   offered again at the next visit; the appearance; the places of the tool windows; that the tour was seen), in the browser's
   storage for the site, through the same code, and writes what is waiting when the page is left; **a script builds the web
   version without trunk** (`scripts/build-web.sh`) with a packed copy of each file, and **the server sends the packed copy**
   to browsers that take gzip (5.4 MB for the 16.4 MB program). `docs/WEB.md` says what was tried in a browser and what was
   not. **Two tabs never write into one session:** the tab opened last keeps it, and an older tab stops and says so. An
   independent review (two readers: the code, the tests) found that first of all (two tabs would have written over each
   other's work), and also that a session which could not be put aside whole was left half moved and then lost, that a
   packed copy of the page left over from another build would have been sent (it is as long in every build: the server now
   compares the checksum too), and that the build script named its files after nothing where `shasum` is missing; all are
   mended with a test each, and the twelve gaps the reader of the tests named are closed. The code of this step was then
   broken on purpose in 74 places, one at a time: the tests caught 69; four more are caught by tests written for them (a
   mark makes the settings folder it is kept in; a folder of a mark's name is no mark; a tab's mark is the time, so two tabs
   that read the same last mark do not take the same one; a packed copy with the right checksum and the wrong length is
   refused); the last is the line that gives a page the browser's storage, which no test on a computer can run and only the
   trial in a browser covers.
   **17 November 2026: a switch that keeps nothing of the work, for a computer that others use** (File ▸ *Remember My
   Work for Next Time*, and the command palette; `docs/WEB.md` says what it does and what it does not cover). Off removes
   what was kept (the session, one on offer, one put aside, one of an older version that could not be read, and in the
   desktop app the list of recent files) and nothing more of them is kept, in that browser or settings folder, until it is
   turned on again; another tab of the page stops too. Three independent readers went through it (the code, the tests, the
   privacy of what is left behind) and what they found is mended, each with a test: the choice was noted before anything was
   removed, so a storage that was full kept everything and said the opposite (what was kept goes first now); the list of
   recent files and a session of an older version were not covered; a copy that was told to remember again while another
   copy had the session said that it did; the command removed a session that was on offer without the question about work
   that was not saved (it asks for the answer first); and the menu promised more than the switch does (its hint now names
   what is not covered). **The reader of the privacy also found an older fault that has nothing to do with the switch: every
   alignment with MAFFT, every pipeline and every pipeline file left the person's sequences in the system's temporary
   folder and never removed them (on Linux, where that folder is everyone's, in a folder that every user of the computer
   could read; on a Mac the temporary folder is the person's own), and a structure prediction on the workstation left its
   folder whenever it failed.** On the Mac this was written on there were 4,363 such folders, most of them from the tests.
   A run's folder now removes itself
   when the run is over, however it ends, and only the person who runs the program may look into it meanwhile
   (`pipeline::tools::WorkFolder`; the bare path that let a run forget is kept for tests only). Folders that earlier
   versions left behind are still there: they are named `helixloom-align-…`, `helixloom-workflow-…`,
   `helixloom-flow-…` and `helixloom-fold-…` in the temporary folder, and may be deleted. The switch and the folder of a
   run were then broken on purpose in 46 places, one at a time (two more were left out: one would have deleted the whole
   temporary folder of the computer, one changes nothing): the tests caught 43; two more are caught by tests written for
   them (the same work is written again when the remembering is turned off and on with nothing done in between; a list of
   recent files that cannot be removed is not shown); the last was a line that nothing needed, and is gone. **17 November 2026: the connection to the workstation is kept for the tab**, by the owner's decision of that day: the
   address, the token, where jobs run and the GPUs of a job, in the table a browser gives one tab, so that loading the page
   again does not ask for the token, and nowhere that stays. Two independent readers (security and code; the tests) went
   through it, and what they found is mended with tests: the page's words promised more than a tab's table gives (a
   duplicated tab has a copy, a browser can bring a closed tab back: the words say so now, and name *Forget* as the sure
   way); when the browser would not take a new connection the one before it stayed, and the page would have come back
   connected to where it had been (the old entry is removed); an entry that makes no connection left its token in the
   table and on the page with no button to remove it (it is removed and not shown); and the page said *Connected* before
   the workstation had answered. `docs/WEB.md` says what a tab's table is and is not, and that HelixLoom should be hosted
   at an address of its own. The code was then broken on purpose in 15 places, and the tests caught all 15; the lines that
   reach the browser's own tables cannot be run by a test on a computer, so they were tried in a real browser against a
   local server (`docs/WEB.md` has what was seen: the token in the tab's table and not in the site's, the page connected
   again after it was loaded again, a second tab without it, *Forget* emptying the table). The formatting and the lints of
   the three builds are clean and every test passes: 1,606 in the app and 1,534 without its optional parts. **Not mended,
   from the readers:** the connection is kept when *Connect* is pressed, before the workstation has answered, so a
   mistyped address is tried again with the token when the page is loaded again (the desktop app does the same with its
   file). **18 November 2026: the results of three pages are kept, in the page and in the desktop app** (the owner's
   choice of that day for the next step): the Alignment page with its sequences and alignment, the tree, and the
   structure, each with its record, through the code that keeps the sequences (`src/ui/results.rs`,
   `Store::keep_session`; `docs/WEB.md`, *What the page keeps*). A result is kept when there is room for it and never
   stands in the way of the sequences. Three independent readers went through it, and what they found is mended with
   tests: a save that failed for another reason than room removed the results that were kept (now only a text that could
   not be written makes them give way); a tree that was opened or pasted kept, and showed on the Phylogeny page, the
   lines about the tree before it; a kept alignment that could not be read took a good list of sequences with it; and the
   words said after reopening were wrong for one and for none. It was then broken on purpose in 79 places: the tests
   caught 77, a test was written for one more, and the last changes nothing that can be seen. Tried in two real browsers,
   one of them with its storage really full (`docs/WEB.md`). The formatting and the lints of the three builds are clean
   and every test passes: 1,641 in the app and 1,569 without its optional parts. **Known limits:** the storage is still
   the few megabytes of `localStorage`; how a page is scrolled, turned or coloured is not kept; "Start fresh" removes
   kept results without asking once more (it asks only for sequences with changes that were not saved); a result kept by
   a later version is said and not shown. **Later the same day, the runs of the Pipelines page** (the owner's choice again, "Yes, that next"): the run that the
   page shows in each of its two modes (the fixed route, and a pipeline file run on a folder), with what was ticked and
   its record, kept in the text that a run already has for the workstation (`remote::tasks`, `remote::flow`) and read back
   with the care an answer from the network is read with (`src/ui/platform/pipelines/kept.rs`). Three independent readers
   went through it, and what they found is mended with tests: a run that had been asked for and had not started was
   overwritten by the run that was put back, and the boxes of the kept run took the place of the ones just ticked (a kept
   run is now not put over a run that is going or has been asked for); the page said a run "could not be read" when it had
   only been refused for that reason (the words are now "could not be put back"); the stamp changed when a folder was
   chosen while only the route's run was shown, which wrote the whole session again for nothing; and the folders of a file
   run, which came from the kept text and were listed at once on the thread that draws, are no longer kept at all (a path
   is not taken from a kept text: they are chosen again to run again). It was then broken on purpose in 55 places: the
   tests caught 52 at once. The other three each take away one of two counts that follow each other (the start of a run
   and its end both count that it has changed), so that either alone changes nothing the keeper can see; with both of a
   pair taken away the tests caught it, in both modes, and the three are left as they are (my reading, not a proof).
   Tried in the built-in browser (the example run kept, the page loaded again, the question *No sequence was open last
   time. Kept from it: the Pipelines page.*, and the run with its record back after *Reopen*) and as a picture of the
   desktop window's question, drawn from a throwaway settings folder (`docs/WEB.md`). The formatting and the lints of the
   three builds are clean and every test passes: 1,652 in the app and 1,580 without its optional parts. **Known limits of
   this step**, from the readers and not tried: a name in a run comes back on one line, without spaces at its ends (that is
   how the workstation's text writes it); a file run with a file that has no name, or with more than 20,000 files, is not
   put back, and the page says so; a run is written whole at each change of it. **19 November 2026: the results go to the
   browser's database (IndexedDB)** (the owner's choice, "Yes, that next"; `src/ui/session/large.rs`, `browser.rs`, `docs/WEB.md`,
   *Where the results are kept*): the open sequences stay in `localStorage` (written at once, so a tab that is closed that
   instant keeps them) and only the result files go to the database, which is read whole before the page starts into a table in
   memory (160 MB at most), the changes being handed to it one at a time in the order made. Three independent readers went through
   the first version, and what they found is mended with tests: a refused change left the database and the table disagreeing (a
   refused removal, an older text left under a name: it could come back as the result of another session; now what may be left is
   taken away with the next change, and the list of the session gives the length and a checksum of each result, so that a text that
   is not the listed one is left out); a refused result took the others of its batch with it (one change is now one transaction);
   a result of an earlier version that was still in `localStorage` could not be removed to make room for the sequences; putting a
   session aside copied every result, so a session over half the memory could never be put aside and held up the sequences (it is
   now a move inside the table); the queue kept superseded texts alive (a change for a name takes the place of those waiting); the
   refusal was told at every save, after the next write instead of before it, and even after the remembering was turned off; turning
   it off did not clear another tab's results; and a page without `localStorage` read the database for nothing. Tried in Brave without
   a window (a real Chromium and its IndexedDB): the example structure is in the database and not in `localStorage`; with
   `localStorage` full to its refusal the antibody (320,596 letters) is kept all the same, which the page could not do before; a
   6,024,435-letter alignment put in the database as an earlier version writes it is offered, read and shown (1,500 sequences of 4,000
   letters) and written again after a change, `localStorage` holding only the list; 6 MB in the database added about 0.1 s to the start
   (154 to 165 ms against 60 to 131 ms without, three starts each); and with a database that never answers the page started after
   about 4 s with the results in `localStorage`. The formatting and the lints of the three builds are clean and every test passes: 1,691
   in the app and 1,619 without its optional parts. **The deliberate-breakage round was run afterwards (19 November 2026): 72 breakages,
   all 72 caught (61 on the first run; 11 got new tests and were then caught).** **Not done:** the mended version was **not read again**
   by independent readers; a database that really refuses (quota) was only tried with a stand-in, not in a browser; Firefox, Safari, two tabs and a
   tab closed in the middle of a write were not tried. **Known limits:** a database that stops answering without an error holds the
   results back until the page is loaded again (the sequences are not affected); results that another tab wrote after this page opened
   are not offered by it until it is loaded again. **Not built:** the tool windows' results,
   projects with names kept in the browser, shareable links.
   **19 November 2026: project files** (`src/ui/project.rs`; the next step of W3 that needs no decision of the owner's): **File ▸ Save Project As…**
   writes the open sequences (each as GenBank with its features) and what the Alignment, Trees, Structure and Pipelines pages show, in the texts
   those pages keep for the session, as one `.hlproj` file; **File ▸ Open Project…**, the Open dialog, a drop on the window and the command line
   open it (the sequences come back beside those that are open, tied to no file, with unsaved changes). The file ends in a line with the numbers of
   sequences and results and the checksum of all the lines before it: a file that was cut short or changed, one of a later version, one that is no
   project and one over 256 MB are refused in words and open nothing. Tested: a write and a read back of texts with tabs and line ends in them,
   every kind of refusal, and a save in one window and an open in another (names, circular or linear, features, the alignment's rows and inputs, the
   workspace, the tab in front). The code was then broken on purpose in 14 places (the checksum, the counts, the version, the end of a
   name, the size limit, the front tab, the workspace, the unsaved flag): all 14 were caught, one only after the test of the front tab was added.
   **Not built, and said so:** projects with names kept in the browser's storage and listed in the app (a file is the unit, kept where
   the person puts it); links to results that another person can open (they need a place to hold the result and a decision about who may read it:
   the workstation keeps results for 7 days behind its token); what the tool windows show. Not tried in a browser or on another computer; the limit of 256 MB was
   not tried against a browser's memory. **Not done for the
   switch:** it does not cover the records fetched from a database or the workstation's connection (the hint says so);
   the files of the settings folder other than the session are written readable by others on the computer where the
   system's defaults say so, and a file half written when the program is killed is not looked for; in a browser the tabs'
   check and write are not one step, so two tabs can cross within a thousandth of a second.
4. **W4 (L)**: accounts and shared projects on the server; background jobs in web workers; installable app (PWA) that
   works offline for the editor and viewers; WebGPU for rendering large trees and alignments.
   **Built in part, 17 November 2026: the page can be installed and opens without a network** (`docs/WEB.md`, *Installed,
   and without a network*): a manifest, two icons drawn by the program itself, and a service worker that keeps the files
   of its build and answers nothing else, so that a call to a workstation is never looked at or kept. Three independent
   readers (the worker as a browser runs it; the build script; the tests) found, and it is mended with tests: two
   HelixLoom pages in two folders of one site would have deleted each other's kept files (a cache is named after the
   page's folder too); a page that a server sends by way of a redirect would not have opened without a network; a server
   that says nothing about keeping could have handed the worker the page of the build before; the script took away files
   by a pattern wide enough to take a person's own, would have deleted its own sources if it was pointed at their folder,
   took the old page away before it knew the new one could be made, and did not say when it built nothing. The worker and
   the script were then broken on purpose in 38 places and the tests caught all 38; the formatting and the lints of the
   three builds are clean and every test passes (1,606 in the app, 1,534 without its optional parts, 8 of the packaged page
   and its worker). Tried in a real browser: with the server stopped the page opened from what the worker had kept, and a
   call to the workstation was not answered from it (`docs/WEB.md` has what was seen and what was not tried). Not done: a
   word in the page when a newer build is waiting to take over; a page that waits on a network that never answers is not
   given the kept one (it is given it when the network says no); results are not kept for use without a network.

## Part 4. The new interface, step by step

Phase 1 is in this release (`docs/DESIGN.md`): colours and typography, top bar, command palette, toasts, Home.

2. **Phase 2 (M)**: bring every workspace to the same look: one header, one controls panel, one results area, real empty
   states, linked views, the inspector.
   **Built in part.** The header was on every page already; the empty states, the inspector and two linked pairs are under
   1.4 and 1.5. **16 November 2026: one controls panel and one results area.** The seven pages that have controls and
   results (Signals, Alignment, Phylogeny, Trees, Structure, Pipelines, Workstation) each built their own column, 270 to 340
   points wide, so the results began at seven different places and jumped sideways from page to page. They are now drawn by
   one piece of code (`design::controls_panel`, `design::results_area`): one card of controls, 320 points wide, with the
   results beside it from the same place on. What the pictures of the seven pages, and then an independent reader of the
   pages, showed to differ went with it:
   * the switch between the modes of Pipelines was above its controls and is in the header, where that of Signals is;
   * *Save figure* of Phylogeny was above the tree and is among the header's actions, as on the other pages, and that of
     Signals is there before a trace is opened, pale, so that the header does not move;
   * the controls of Structure began with buttons and begin with a heading, as the others do;
   * **where an analysis runs and the button that runs it are one piece** (`run_button`): the line was several controls
     above the button on Alignment and Pipelines, and under it for a pipeline file; the room above the button was 0, 6 or 8
     points; a pipeline file had no wheel while it ran and called its button *Run the pipeline* where the route says *Run
     pipeline*. The button was a plain one on three pages (*Align*, *Estimate tree*, *Predict structure*) and is the filled
     one on all (**my choice of rule, not the owner's: say so if the next step, such as *Estimate a tree* in the header of
     Alignment, should be the only filled button**);
   * the frame round a view (an alignment, a tree, a structure, the strip of its residues, a trace) was written out six
     times and is one function (`view_frame`), which also gives the views the thicker outline of the high-contrast setting,
     as cards have; the titled cards of Pipelines and the Workstation were two copies of one function and the cards of a
     pipeline file had smaller titles (`titled_card`); the room between cards was 8 or 10 points (`CARD_GAP`); *Estimate a
     tree instead* on Trees was a small button among ordinary ones.

   `docs/DESIGN.md` (*The frame of a page*) has the rules and `ui/platform/same_look.rs` holds them for the seven pages,
   with and without their examples. Looked at in pictures of each page before and after (the desktop build, light theme,
   1280 by 800 points). **Reviewed by three independent readers** (the code, the tests, what still differs): the reader of
   the code found no fault; the reader of the tests named eleven gaps, of which the ones that could let a fault through are
   closed (the tests waited for no search for programs and so depended on the computer: they now run with none installed;
   where the results begin was only checked on an empty page; a second card in the column, a second switch left behind, a
   heading that is only capital letters, a button that can never be pressed all passed) and two are stated here (one size
   of window is tried; the waits are real time). Three older tests of the Structure page failed on the way: the new heading
   had pushed two buttons under the edge of the column on the test screen. A person would scroll, so the click of the tests
   now turns the wheel over the column when what it is to click is out of sight, which mends them and every later one.
   The code of this step was then broken on purpose in 29 places, one at a time (the width, the scrolling and the card of
   the column, the margin of the results, the outline of a view, the size of a heading, the run button and its line in each
   of its states, the switch, the figure button, the results of three pages, the click that scrolls): the tests caught all
   29, four of them only because of tests that were written while the list was made (a view has the outline of a card; the
   figure button is pale until there is something to save; where an analysis runs is said right above its button; the
   button cannot be pressed while its analysis runs, and a wheel turns beside it). With all of it in, the formatting and
   the lints of the three builds (the app, the app without its optional parts, the browser) are clean and every test
   passes: 1,578 in the app and 1,506 without the optional parts, and the tests of the command line and the server.
   **Not done, from the reader of the pages:** a message after a run is red on Signals when no trace is open and grey
   elsewhere, and under the button, at the end of the card or under the hint depending on the page; in a short window the
   results of Alignment and Phylogeny do not scroll (those of Pipelines and the Workstation do); the Reads tab of Signals is
   one card across the page, with no column; the buttons that send a result on are in the header on three pages, in the
   controls on Signals and in the result cards on Pipelines; hints are 11.5, 12 or 12.5 points; the sections of Structure
   change order when a structure is opened; Trees alone has a card at the right, which narrows the tree when one is opened.
   **17 November 2026: the pages in a small window.** A tab of a browser on a small screen has far fewer points than the
   window the pages were drawn for, and so has a window with the appearance's larger sizes. Pictures of the seven pages in
   a window of 853 by 533 points showed three faults, all mended: **the buttons of a header were drawn over its line of
   purpose** (the title and that line now have the room the buttons leave and are cut short with the rest under the
   pointer; the buttons go on a row of their own when even the title has none); **the 3D structure was not drawn at all**
   (two choices of the row above it stood in rooms of their own, which do not go to the next row with the others: they
   pushed the rest of the row and the whole viewer out of the window; every choice is one widget now); and **a tree was a
   sliver 109 points wide** between its controls and the card that says what it is (that card goes under the tree when the
   tree would have less than 360 points). The rail was right already: it squeezes its buttons and then scrolls.
   `ui/platform/same_look.rs` holds the frame and the header for all seven pages in windows of 853 by 533 and 640 by 400
   points too, and the header is tested by itself in rooms of 1,100, 620 and 260 points. The pictures were needed twice
   over: the first mend of the header passed its tests and still left Trees with nothing but its header, because the row
   the buttons were given took the whole height of the page and the tests of the pages only looked for the controls, not
   for where they were (they now say that the page begins right under the header). The code of this step was then broken
   on purpose in 15 places: the tests caught 10, and the five that got through each showed something: the test of the
   header compared words with words, so a button's frame over the line was not seen; the buttons never had to go on to a
   second row with the letters the tests draw; and two of the structure's choices were made another way than the ones the
   test clicked. The title now decides when the buttons go under it (it is never cut), every choice of the structure's row
   is the same widget, and with the tests that were added the round was run again: 14 of 14 caught. **Two independent
   readers** (the code, the tests) then found, and it is mended with tests: the card under the tree took whatever height
   its contents asked for, so a picked clade with many tips left the tree no room and made it jump under the pointer that
   picked it (the card has one height now and scrolls); the wrapping fault that hid the structure was waiting in two more
   places, the header of Signals and the PCR window; the least width of a tree was held by two windows only, so most
   values of it passed (four windows now, two on each side of it); and the tests of the header did not measure how high it
   is, or that its line is whole where it has room. A test of the session's lock that is older than this work failed once
   in the full run: it asked for the lock once, at a moment when another test was starting a program, which is known to
   refuse it for some tens of milliseconds; it asks for a moment now. Not mended, from the readers: the width of a header's
   actions is remembered from the frame before, so the first frame of a page, and the frame after its actions change,
   draws them over the end of the line; the tests draw with egui's own letters, which are narrower than the program's, so
   a width that just fits in a test may not in the window (the pictures are the check for that). With all of this in, the
   formatting and the lints of the three builds are clean and every test passes: 1,600 in the app and 1,528 without its
   optional parts, and the tests of the command line and the server. **Still not done for small windows:** below about 700 points of width the column of controls leaves the
   results too little to work in (it does not fold away or go above the results); in a short window the results of
   Alignment and Phylogeny do not scroll; a phone was not tried. **Looked at in a browser on 17 November 2026**, once
   `wasm-bindgen` was on this Mac again and the page could be built: the frame in the dark theme, Trees and Structure in
   a window of 853 by 533 points, and the switch of W3 turned off and on with the browser's storage read after each
   (`docs/WEB.md` has what was seen). The dark theme of the other pages, and Firefox and Safari, were not looked at.
3. **Phase 3 (L)**: docking panels (drag a panel to a new place), tabs for results, saved layouts.
4. **Phase 4 (M)**: motion and polish: smooth transitions, loading skeletons, the typography pair, an illustration set.

## Part 5. Suggested order

| Milestone | Contents | Why first |
|-----------|----------|-----------|
| **A (2 to 3 weeks)** | Autosave and session restore; recent files; job tray; window placement. **Built in part, November 2026: see 1.2, 1.3 and 1.6 for what is and is not** | Removes the fear of losing work and the feeling of not knowing what runs |
| **B (3 to 4 weeks)** | Linked views and inspector; empty states; first-run tour; typography. **Built in part, November 2026: the inspector, the column of an alignment linked with the residue of a structure, the selected ORF in the ORF list, what the picked bases of a trace are, the empty states of six pages and the tour (1.4, 1.5); typography is built (1.6)** | Makes the daily screens feel finished |
| **C (4 to 6 weeks)** | Pipelines as data with caching; connectors (NCBI, UniProt, PDB). **Built in part, November 2026: pipeline files with a cache and a batch run (the loci → gene trees → species tree route), here and on the workstation, and the connectors by accession; see 2.1 and 2.2 for what is not** | Turns the tools into a platform |
| **D (4 to 6 weeks)** | Motifs and domains (HMMER, Pfam, conservation on structure) | First of your two research directions; light compute |
| **E (6 to 8 weeks)** | Epitope and paratope models on the Workstation; benchmark page | Second direction; needs GPUs, licences and validation |
| **F (parallel)** | Web W2 and W3; installers and updates; CI on three systems | Reach |
| **G (after C)** | The species-tree family: a batch pipeline from locus alignments to the species tree, weighted ASTRAL, gene families, bootstrap. **Built in part, November 2026: the batch pipeline (here and on the workstation), the weights of weighted ASTRAL (not compared with wastral, without its weighted branch support) and the gene family score of ASTRAL-Pro (2020; not compared with its program, without branch supports), and the bootstrap over genes and sites (not compared with ASTRAL's)** | Builds on the ASTRAL core that is done and on the pipeline engine |

## Questions for the owner

**Answered by the owner, November 2026** (questions 1 to 5 and 7 to 10; 6 follows from 1):

1. **Licence** for the code: **MIT**, open to everyone (`LICENSE`; the dependencies were checked, see below). *Was: MIT, Apache-2.0, GPL-3.0, or closed; the repository had none, so nobody else could reuse it.*
2. **Who is it for first**: **computational people who run pipelines**, not (first) wet-lab biologists. So the order of work favours the pipelines, the
   workstation and the batch methods over the editor's polish.
3. **Which models first**, and may some be built into the program: **all run as separate programs on the workstation**; none is built in.
4. **Hosting** for the web version: **the owner's own server**.
5. **PhyloCNN**: the trained networks **stay a local add-on**; they are not shared or bundled, and no permission is asked of the authors.
6. **The AGPL-3.0 question.** **Settled by 1: with MIT, no code of ASTER (AGPL-3.0) may be ported; the species-tree methods are written from the papers, as the
   ASTRAL core was.** ASTER is under the AGPL-3.0. If you chose that licence for HelixLoom, parts of ASTER could be
   ported directly instead of being written again from the papers, which is faster; the price is that everyone who uses or
   hosts HelixLoom must offer the source of their version. The ASTRAL core was written from the papers, so it is free of
   the question either way (ASTRAL itself is Apache-2.0).
7. **Which of the species-tree extras first**: **the batch pipeline from locus alignments, then the weights, then gene families**; the genome-alignment
   methods (see 2.3) are not first.
8. **Which step of the motif, domain, epitope and paratope plan next (2.6)**: HMMER domain search against Pfam (a database to
   download), motif discovery, or the antibody and epitope predictors, and which licences you hold: BepiPred, DiscoTope and
   NetMHCpan are free for academic use only, and the protein language models and Boltz-2 need weights downloaded to the workstation.
   **Answered: all three, HMMER domain search, motif discovery and the antibody and epitope predictors** (the order is not given; the plan
   below goes by how little each needs: motif discovery, then HMMER, then the predictors). **The licences: the answer ticked "academic
   licences", "language-model weights" and "Boltz-2 weights" and also "none yet", which cannot all be true, so none is taken as held
   until the owner says which**: nothing is downloaded without a yes to the filename, source and size, and a model that is free for
   academics only is built to run as the licence holder's own program on their own workstation, never bundled.
9. **Optimised development builds?** The dev and test profile is `opt-level = 1` (and 2 for the libraries), which makes no vector
   instructions, so `cargo run` and `cargo test` align 2.9 times faster than before where a release build is 7 times; `opt-level = 2` for
   `helixloom` itself in `[profile.dev.package.helixloom]` gives about 1.9 ns a cell, at the price of slower compiles (2.5, item 6a4).
   **Answered: yes, and done** (`Cargo.toml`, `[profile.dev.package.helixloom] opt-level = 2`).
10. **The aligner's scoring in whole numbers** (+10, −1, −10, −2, instead of +1, −0.1, −1, −0.2)? Every implementation, the cards' too,
    would then equal the processor's by construction; alignments change only where rounding broke a tie today (2.5, item 6a).
    **Answered: yes, and done** (November 2026): match +10, mismatch −1, first gap −10, further gaps −2 (`MATCH` and `Scoring` in
    `seq/align.rs`, the PyTorch worker, the Align window, `MECHANICS.md`). Scores shown are ten times the original's. Measured with a
    throwaway script (not kept) that mirrors the aligner in 32-bit floats with the old and the new numbers, on 600 random pairs of
    DNA (related and unrelated) aligned end to end: about a quarter of the alignments differ, and in **every** one both alignments
    score the optimum exactly, so the difference is the choice among equally good alignments, which the old rounding made
    arbitrarily and the whole numbers make by the rule (a match, then a gap in `a`, then a gap in `b`).

**Said by the owner on 16 November 2026:**

* **The web app first.** HelixLoom is to reach macOS, Windows and Linux, and the way there is a web app that is built and
  deployed, before an installer for each system. So the work goes to what the web version lacks (part 3; `docs/WEB.md`),
  and the installers of 1.8 wait.
* **No Apple Developer ID** (question 12): the Mac app stays signed ad hoc, and signing and notarization are not planned.
* **The code is not pushed to the GitHub repository for now**, so the checks there (`.github/workflows/ci.yml`) do not run:
  they are run here by hand. Question 11 (where releases are published) is still open, and matters less while the web app
  is the way in.
* Next after the roadmap's order: the interface phases 2, 3 and 4 (part 4).

**Said by the owner on 19 November 2026:** asked to update the README, write how to deploy on Ubuntu and how to use the app and the web
app, and push the code to GitHub (done: commit `4c55985`). This lifts "not pushed for now" above; a push still waits for the owner's
word each time, as does a commit.

**Said by the owner on 17 November 2026:**

* **A page may keep the workstation's token until its tab is closed** (asked: "May the web page keep the workstation's
  token until its tab is closed?"; answered: "Yes, until closed"). So the page keeps the connection in the table the
  browser gives one tab, and never in the one that stays (part 3, W3; `docs/WEB.md`, *The connection, kept for the tab*).

**Open, November 2026** (both decide what the last item of the list of ten, installers and updates, can still do):

11. **Where will releases be published?** The releases page of the repository on GitHub, your own server (as for the web
    version), or nowhere yet. The app can check for a new version only against such a place; until you say which, it checks
    nothing and contacts nothing about itself.
12. **An Apple Developer ID for the Mac app?** With one (a yearly fee to Apple) the app can be signed and notarized, and opens
    on any Mac without a warning. Without one it stays signed ad hoc: it works, and other people must right-click and choose
    Open the first time. The same question will come for Windows (a code-signing certificate).

**The licence check behind question 1 (November 2026).** Of the 481 packages in `Cargo.lock`, 474 were read from the local registry: every one is MIT,
Apache-2.0, BSD, ISC, Zlib, Unicode, BSL-1.0 or Unlicense, or offers one of those among its choices (`r-efi` twice, `MIT OR Apache-2.0 OR LGPL-2.1-or-later`,
and `self_cell`, `Apache-2.0 OR GPL-2.0-only`: used under the permissive choice). **Not read:** seven packages that are only built with the optional
CUDA feature and are not in the local registry (`cudarc` twice, `candle-kernels`, `bindgen_cuda`, `ug-cuda`, `float8`, `libloading`), so check them
when the CUDA build is first made. The bundled fonts and the 1VFB structure keep their own licences (SIL OFL 1.1, CC0), listed in `CREDITS.md`. What
the program re-implements from ApE, ASTRAL, PhyloCNN and the papers is written from their behaviour or their papers and uses none of their code
(`CREDITS.md` says so for each); the one place that depends on a licence that does not exist is PhyloCNN's repository, which has no licence file, so its
trained weights are neither included nor shared. *This is the author's reading of the files, not legal advice.*
