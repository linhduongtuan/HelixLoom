# Using HelixLoom: as an app and as a web app

HelixLoom is one program in two forms. The **app** is a window on your computer (macOS, Windows, Linux). The **web app** is the
same program compiled for the browser: you open an address and it runs in a tab. They look and behave alike; what differs is
listed in *App or web app?* below.

## Start it

**The app.** Build it once (you need [Rust](https://rustup.rs); on Ubuntu see [UBUNTU.md](UBUNTU.md)), then run it:

```sh
git clone https://github.com/linhduongtuan/HelixLoom.git && cd HelixLoom
cargo run --release                     # the window opens
cargo run --release -- my.gb trace.ab1  # open files at the start
```

On a Mac, `scripts/bundle-macos.sh` makes `target/bundle/HelixLoom.app` and a disk image. It is signed ad hoc: on another Mac the
first start needs a right-click and **Open** ([DEPLOY.md](../DEPLOY.md), *A Mac app*).

**The web app.** Someone has to host it once: a folder of static files, made by `scripts/build-web.sh` and handed out by
`helixloom serve --web-dir dist` or by any web server ([WEB.md](WEB.md); on Ubuntu, [UBUNTU.md](UBUNTU.md)). You then open its
address, for example `http://127.0.0.1:7878/` through an SSH tunnel. Your browser needs WebGL 2 (a current Chrome, Edge, Firefox
or Safari should have it; **only Chromium-based browsers (Brave) have been tried**). The first load downloads about 5 MB when the
server sends the packed copy (`helixloom serve --web-dir` does), 16 MB when it does not.
**Install it like a program** with the browser's *Install* button, when it offers one (it does so only for a page opened from
`https://` or from the computer itself): it then opens in its own window and also without a network.

Address options: `?view=trees` opens a page with its example (`?view=phylogeny`, `?view=pipelines`, …), and `?dark` starts in
the dark theme.

## The first five minutes

1. **Home** shows a card for each kind of analysis, each with **Try an example**: press one and the page fills with a made-up
   alignment, trace, tree or structure, so you can see what it does before you use your own data. The first run also shows a
   three-card tour (**Help ▸ Show the Tour** brings it back), and **Help ▸ Open the Sample Project** opens a plasmid, a trace and a
   tree together.
2. **Press Cmd+K** (Ctrl+K on Windows and Linux): one box that finds every command and tool.
3. **Open your own files**: **File ▸ Open…**, or drop a file on the window. The end of its name decides where it opens: a trace
   (`.ab1`, `.abi`, `.scf`) in *Signals*, a structure (`.pdb`, `.ent`) in *Structure*, a Newick tree (`.nwk`, `.newick`, `.tre`,
   `.tree`, `.treefile`) in *Trees*, and anything else (GenBank, FASTA, plain letters) in *Sequences* — or, while *Alignment* or
   *Phylogeny* is showing, added to the alignment's list. Sequencing reads (FASTQ) have their own tools on the *Signals* page.
   **File ▸ Fetch from a Database…** gets a record of NCBI, UniProt, the PDB or the AlphaFold DB by its accession (in the app; not
   in the browser).
4. **Work page by page** (the rail at the left):

   | Page | What you do there |
   | --- | --- |
   | **Sequences** | edit a plasmid or a sequence; circular and linear maps, restriction sites and digests, ORFs, features, primers and PCR (**Tools ▸ Primer3 Pair Design…** picks primer pairs, and a TaqMan probe between them when you tick the box); assemble reads; calculators |
   | **Signals** | look at a Sanger trace and the quality of each base, trim the weak ends, export; summarise and clean FASTQ reads |
   | **Alignment** | gather sequences (*Add open sequences*, *Open FASTA…* (it reads GenBank, FASTA and Clustal `.aln`), or paste), choose the method (built-in centre-star, built-in progressive after ClustalW, or MAFFT) and press **Align**; a colour-coded viewer with consensus, identity and conservation; click a letter to select its column; save as FASTA or Clustal; tick **Align by the translation** to align coding DNA through its proteins and keep the codons whole (`docs/CODON_ALIGNMENT.md`); **Pick primers and probes…** designs primers and a TaqMan probe from the columns the sequences share and marks them on the alignment (`docs/ALIGNMENT_PRIMERS.md`) |
   | **Phylogeny** | build a tree from the alignment: neighbour-joining, UPGMA, maximum likelihood, Bayesian; bootstrap support; the species tree from many gene trees |
   | **Trees** | draw a tree as a phylogram, cladogram, circle or radial (unrooted) fan; read Newick and NEXUS with their annotations; colour, label and mark by them; collapse, highlight, root, rotate, order; add support from other trees; an axis; save as SVG, NEXUS or a figure (`docs/TREES.md`) |
   | **Structure** | open a PDB file, or predict a structure from a sequence; turn it around in 3D, colour it by confidence, secondary structure, chain, or the alignment's conservation; measure an interface between chains |
   | **Pipelines** | align, build a tree and pick primers in one run (*Fixed route*), or run a pipeline file on a folder of loci (*Pipeline file*) and keep the answer of every step |
   | **Workstation** | send the heavy jobs to a computer with GPUs and read the results here (below) |

5. **Save what you made.** **File ▸ Save**, **Save As…** and **Export as FASTA…** for sequences; **Save figure** on a map, an alignment, a trace, a tree or a structure
   (SVG, or a PNG with four pixels to a point); **File ▸ Save a Report of the Results…** writes one HTML page of the results with
   the methods in sentences and every input with its checksum. In the browser, a save is a download.
   **File ▸ Save Project As…** writes the open sequences and what the Alignment, Trees, Structure and Pipelines pages show as one
   `.hlproj` file, to keep, to send to someone, or to move to another computer or a browser; **File ▸ Open Project…** (or dropping the file on the
   window) puts it back beside what is open. The sequences come back tied to no file and with unsaved changes, so that saving one asks where. A
   file that was cut short or changed is refused and nothing is opened.

Under a result, **How this was made** says what made it (every input with its checksum, the settings, the programs and their
versions, the day, where it ran), what has changed since, and runs it again with the same settings.
**Help ▸ How Each Analysis Was Checked** says, for each analysis, what it was compared with and what was *not* checked.
HelixLoom has not been validated for clinical or diagnostic use.

## Your work is remembered

The open sequences, and what the Alignment, Trees, Structure and Pipelines pages show, are saved a few seconds after each change
and offered the next time you start: *Reopen your last session?* The app keeps them in its settings folder (`~/Library/Application
Support/HelixLoom` on a Mac, `~/.config/helixloom` on Linux, `%APPDATA%\HelixLoom` on Windows); the web app keeps them in the
browser, for that site, on that computer ([WEB.md](WEB.md), *What the page keeps*).

**On a computer that others use, turn it off:** **File ▸ Remember My Work for Next Time** removes what was kept and keeps nothing
more until you turn it on again. *Start fresh*, at the question, forgets the last session. In the browser, clearing the site's
data does too. Nothing is ever sent anywhere by this.

## Updates

**Help ▸ Check for Updates…** (or `helixloom update`) asks GitHub whether a later release exists. Nothing is asked until you press the button, and
the program never checks by itself. If one is out, the window shows what is new and can download the file for your computer (Mac with Apple
silicon, Windows) into your Downloads folder, checked against the release's checksums; you then open the disk image or unzip the file yourself, since
the program does not replace itself. The web version updates when you reload the page. See [DEPLOY.md](../DEPLOY.md) for updating the workstation.

## Heavy jobs on a workstation (optional)

Alignments, trees, species trees, pipelines, PhyloCNN estimates and structure predictions can run on a computer with GPUs instead of
yours ([DEPLOY.md](../DEPLOY.md); on Ubuntu, [UBUNTU.md](UBUNTU.md)). Open **Workstation**, type the address and the token, press
**Connect**, then choose **Where jobs run**: *This computer*, *Workstation*, or *Automatic* (big jobs go there, small ones stay).
The same choice sits beside the button of every page that starts a job. You can close the app while a job runs: the result waits
on the workstation for 7 days, and **Open** beside the job brings it back. If the workstation cannot be reached the run says so; it
never quietly falls back to your computer. When the page is served by the workstation itself, the address is already filled in.

## App or web app?

| | App | Web app |
| --- | --- | --- |
| Editor, alignment, trees, structure viewer, traces, calculators | yes | yes, on the page's own thread: a long job on this computer stalls the page while it runs (sent to a workstation, it does not) |
| Open and save files | file dialogs, drag and drop | file dialogs and drag and drop; saving is a download |
| MAFFT, IQ-TREE, HMMER, primer3, Snakemake, Nextflow, fastp, curl on this computer | used when installed (**Tools ▸ Programs HelixLoom Calls…**) | not available: the built-in methods are used; sent to a workstation, they run there |
| Fetch from a database; online structure prediction | yes | no (needs a program the page cannot start) |
| Workstation | settings kept in a file | the address and token are kept for the **tab** only; *Forget* removes them |
| Work remembered | settings folder | the browser's storage, and its database for big results |
| Works without a network | yes | when installed |
| GPU engine on this computer | yes (`--features cuda` or Metal on a Mac) | no |

## If something does not work

* **Nothing opens a window / the page stays on "Loading…"**: see [UBUNTU.md](UBUNTU.md) section 9, or [WEB.md](WEB.md); a browser
  without WebGL 2 cannot draw the page.
* **"The workstation cannot be reached"**: the address, the token, or the tunnel; [DEPLOY.md](../DEPLOY.md) section 6.
* **A report of a problem**: **Help ▸ Copy Diagnostics** copies the version, the system, the programs found and the place of the last
  fault, and nothing of your work (no sequence, no file name); nothing is sent anywhere.
