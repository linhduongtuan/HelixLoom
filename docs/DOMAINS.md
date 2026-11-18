# Protein domains (Tools → Protein Domains, with HMMER and Pfam)

A *domain* is a part of a protein that folds and works on its own and recurs in many proteins. **Pfam** is a database of protein
families, each described by a profile hidden Markov model; **HMMER 3**'s `hmmscan` searches a protein against all of them and
says where in the protein each family matches. This window runs that search and draws the answer along the sequence.

HMMER and Pfam are other projects' program and data, not part of HelixLoom: the program runs `hmmscan` as it runs MAFFT and IQ-TREE
and reads the table it writes (`--domtblout`). Nothing of HMMER's code is here.

## What you need on the computer that runs the search

1. **HMMER 3** with `hmmscan` on the PATH (HMMER is BSD-3-Clause): for example `brew install hmmer` on a Mac, `apt install hmmer` on
   Ubuntu, or `conda install -c bioconda hmmer`. The Workstation page shows `HMMER 3.4` among the programs it found.
2. **Pfam-A**, pressed: download `Pfam-A.hmm.gz` from `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/` (about 658 MB when
   this was written), unpack it, run `hmmpress Pfam-A.hmm` (it makes `Pfam-A.hmm.h3f`, `.h3i`, `.h3m` and `.h3p` beside it), and put
   the five files in the `pfam` folder of the models folder (`docs/MODELS.md`; on a workstation, the workstation's).

The search **never downloads anything itself**, and a request cannot name a program or a database: where either is missing the
window says what to install. The program is found on the computer that runs the job, so with the choice set to *Workstation* it
is the workstation's HMMER and Pfam that are used.

## Using it

Paste protein sequences (FASTA, or one on each line) or press *Use an example* (ubiquitin), choose how a hit is judged, and press
*Find domains*:

* **Pfam's own thresholds** (`--cut_ga`, the default): each Pfam family carries a curated gathering threshold, the score above which
  Pfam itself says that a sequence belongs to the family. This is the setting that makes the domains the ones Pfam would report.
* **E-value at most X**: one cutoff for the sequence and the domain alike (`-E X --domE X`; 0 < X ≤ 100).

Where it runs is the usual three-way choice (*This computer*, *Workstation*, *Automatic*: more than 100 sequences go to the
workstation). Each sequence is drawn as a bar with its domains as coloured boxes under it, the same family always in the same
colour (a colour-blind-safe set of eight, so unrelated families can share one); a box says its family when it is wide enough,
and the hover text gives the accession, description, residues, score, E-value and the part of the profile that matched.
Domains that overlap are drawn one under the other. *Copy as a table* gives every domain, tab-separated, for a spreadsheet.

At most 10,000 sequences, 20 million residues in all and 100,000 residues in a sequence are searched at once.

## What the numbers are

HMMER's own, from the domain table (HMMER 3.4 User's Guide): the **score** is the bit score of the domain; the **E-value** shown
is the *independent* E-value (what the comparison would have received if this were the only domain found; it counts all the
profiles in the database, so it is the stringent measure; the conditional one is in the copied table); the residues drawn are the
**alignment** (the envelope, which may reach beyond it, is in the hover text); **Profile** is the part of the family's model that
matched, out of its length.

## Limits

* **Run against a real HMMER and the real Pfam by hand, not in the usual tests.** In November 2026, with HMMER 3.4 (Homebrew)
  and Pfam 39.0 (`Pfam-A.hmm.gz`, 658,389,433 bytes, checked against the EBI's MD5 and pressed: about 8 GB in the `pfam` folder
  with the text file, 4.2 GB without it and without the download):
  * the protein of HMMER's own tutorial (`7LESS_DROME`, 2,554 residues) against its two tutorial families gives, through the
    program, **the ten domains that `hmmscan` writes when it is run directly**, with the same places on the profile and on the
    sequence, the same scores and the same accuracies, and the seven of them that pass Pfam's own thresholds;
  * the four proteins of the structures that come with the program (the example of the window) against the whole of Pfam, in
    about two seconds on the Mac it was run on: ubiquitin is in the ubiquitin family (PF00240), both variable domains of the antibody D1.3 in the V-set
    (PF07686), and hen egg-white lysozyme in the C-type lysozymes (PF00062), each over more than 85 % of its length.

  Both are tests that are left out of the usual run, because HMMER and Pfam are not in the repository
  (`src/remote/domains/real.rs`: `HELIXLOOM_HMMER_TUTORIAL=/opt/homebrew/share/doc/hmmer/tutorial cargo test --bin helixloom --
  real_hmmer real_pfam --ignored`). What runs on every change is the reading of tables made by hand to the documented format
  and a program that stands in for `hmmscan`. **Not tried:** a large set of sequences, a search on the workstation with a real
  HMMER (the job was tried with the stand-in), Linux and Windows. The first real run showed nothing wrong in the reading; it did
  show that a FASTA name is its first word (two sequences named "antibody …" looked alike) and that the names did not keep one
  column, which is put right.
* **Pfam's clans are not resolved.** HMMER reports every family that matches; Pfam's rule that two families of one clan may not
  overlap is applied by its own `pfam_scan.pl` and is not here. Overlaps are drawn and said.
* Protein only. InterProScan, CDD, SMART and the other databases are not searched; one database (`Pfam-A.hmm`) is used.
* The sequence names go to HMMER as `q1`, `q2`, … and are matched back by number, so that no name can matter to the program.
* There is no `helixloom` command for this: `hmmscan` is the command, and this window and the job add the drawing, the choice of
  where to run it and the matching of the table to the sequences.
