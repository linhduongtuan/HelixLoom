# Running the heavy work on a workstation with GPUs

HelixLoom can send its heavy jobs to another computer, a workstation with several GPUs, and bring the results back to the
Mac, Windows or Linux computer you work on. The workstation runs one program, `helixloom serve`. Your own computers run
the HelixLoom app as usual and have one extra page, **Workstation**, where you connect to it.

```
  your Mac / Windows PC / Linux PC                      the workstation (several GPUs)
  ┌─────────────────────────────┐   jobs go over the    ┌────────────────────────────────┐
  │ HelixLoom app                │   network (a private  │ helixloom serve                 │
  │  Alignment, Phylogeny,      │ ────────────────────▶ │  queue ─▶ one process per job  │
  │  Pipelines, PhyloCNN,       │   network or an SSH   │  each job sees only its GPUs   │
  │  Structure                  │   tunnel)             │  results kept on disk          │
  │  Workstation page: connect, │ ◀──────────────────── │  (7 days, and across restarts) │
  │  jobs list, Open results    │   results come back   └────────────────────────────────┘
  └─────────────────────────────┘
```

What is sent, and what is not:

| In the app | Where it runs when the job goes to the workstation |
| --- | --- |
| Alignment (built-in or MAFFT) | the workstation's processor |
| Phylogeny: neighbour-joining, UPGMA, maximum likelihood, Bayesian, bootstrap | the workstation's processor, and its GPUs when the job is big enough to be worth them (maximum likelihood and Bayesian over many alignment columns; the distances and the bootstrap of many sequences): see step 1.7 |
| Pipelines (align, tree, primers) | the workstation's processor, with its MAFFT, IQ-TREE and primer3 |
| PhyloCNN epidemic estimates | the workstation, on one GPU or on several at once (see "GPUs for one job" below) when "many trees run on the GPU" is chosen |
| Structure prediction | the workstation's GPU with your folding program (one GPU, or more if the server is told the program uses them: `--fold-gpus`), or the ESM Atlas web service asked from the workstation |
| Plasmid editing, the contig assembler, sequencing reads, chromatograms, tree drawing, primer3 window, tree tools | stays on your computer: they are quick, or they are about what is on your screen |

**You choose where the jobs run**, and can change it at any time, on the Workstation page or on any page that starts a job
(Alignment, Phylogeny, Pipelines, PhyloCNN, Structure), where the same choice sits beside the button:

| Choice | What happens |
| --- | --- |
| **This computer** (the start) | Every job runs on the computer you are using. Nothing is sent anywhere. |
| **Workstation** | Every job is uploaded, computed on the workstation and brought back. |
| **Automatic** | Big jobs go to the workstation and small ones stay here, where there is nothing to upload. A job is big when it has more than 100 sequences or 1,000,000 letters, asks for more than 200 bootstrap resamples or 100,000 Bayesian generations, has more than 100 trees, or can use a GPU there (PhyloCNN with the GPU ticked, a structure the workstation folds on its own GPU). The page says which way the job as you have set it up would go. |

Workstation and Automatic need a connected workstation; with none, everything runs on your computer. These sizes are
estimates, not measurements; if a job should go the other way, pick This computer or Workstation for it.

**GPUs for one job.** When the connected workstation has two or more GPUs, the Workstation page (under *Where jobs run*)
also asks how many of them one job is given: **1 GPU** (the start: every job that needs a GPU takes one, so several jobs run
side by side), **2 GPUs**, **3 GPUs** and so on up to what the workstation has, or **All GPUs**. A job that can spread its
work over several GPUs (PhyloCNN with many trees: it shares the trees out over them) is given that many at once, or fewer if
it has fewer trees to share or the workstation has fewer GPUs. A big tree job (step 1.7) uses every GPU it is given, so choose **All GPUs** to have it use both of two. A job that
cannot use more than one GPU (a structure, unless the server was told its program does) is given one whatever you choose. The
choice is remembered with the address and token.

## 1. On the workstation

You need a Linux workstation (Ubuntu is assumed below), with the NVIDIA driver installed: `nvidia-smi` should list your
GPUs. Nothing else is needed for jobs that run on the processor. For the GPU parts you also need the CUDA toolkit (to build
the GPU version, step 1.2) and Python (for the folding program, step 5).

### 1.1 Get the program there

Get the code onto the workstation (clone it there; copying a folder from another computer must leave out its `target` folder,
which holds that computer's build), then install Rust once and build. `docs/UBUNTU.md` has the whole path on Ubuntu, with the web
version and the service:

```sh
ssh you@workstation
sudo apt install build-essential pkg-config git curl libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev libxkbcommon-dev libssl-dev
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh     # follow the questions, then open a new shell
git clone https://github.com/linhduongtuan/HelixLoom.git && cd HelixLoom
```

### 1.2 Build it

```sh
cargo build --release                       # processor only
cargo build --release --features cuda       # also lets PhyloCNN and big trees use the NVIDIA GPUs through candle (needs the CUDA toolkit)
```

**Or use PyTorch for the big tree jobs, and build with no CUDA toolkit at all.** The ordinary build runs the tensor work of tree jobs in
Python worker processes, one for each GPU, when it is told which Python to use. You need a Python with PyTorch (and NumPy) on the
workstation; the program installs nothing:

```sh
cargo build --release                       # no --features cuda needed for this
sudo python3 -m venv /opt/helixloom/torch-env       # where the service's user can read it (the unit hides /home)
sudo /opt/helixloom/torch-env/bin/pip install torch numpy
# it must list both cards, and sm_86 and sm_120 (RTX 3090 and RTX 5060 Ti); if sm_120 is missing, install the build for CUDA 12.8 or later from pytorch.org
/opt/helixloom/torch-env/bin/python -c 'import torch; print(torch.cuda.device_count(), torch.cuda.get_arch_list())'
export HELIXLOOM_TORCH_PYTHON=/opt/helixloom/torch-env/bin/python   # for the server: the `Environment=` line in the service file (step 3)
./target/release/helixloom gpu-check --bench
```

This only covers the distances, the bootstrap and the likelihood of a tree (step 1.7); PhyloCNN still needs `--features cuda`. The path
of the Python is the server owner's, in the server's environment, as `--fold-command` is: nothing sent to the server ever names a program.

The program is then `target/release/helixloom`. If the `cuda` build stops with a message about `nvcc` or CUDA, the toolkit
is missing or not on the path: `nvcc --version` must work first.

**Two different cards (an RTX 3090 with an RTX 5060 Ti, say).** candle compiles its kernels once, for one compute capability
(`CUDA_COMPUTE_CAP`; if it is not set, the build asks `nvidia-smi`, which gives whichever card it lists first). Set it to the
*oldest* card's, so that the newer one compiles the same kernels when it first uses them:

```sh
# the CUDA toolkit from NVIDIA's repository, version 12.8 or later (the RTX 5060 Ti is Blackwell and needs it, for the program's
# build and for the libcublas it loads while running; Ubuntu's own `nvidia-cuda-toolkit` is 12.0, too old), and a driver of 570
# or later with the open kernel modules, which Blackwell cards need
CUDA_COMPUTE_CAP=86 cargo build --release --features cuda   # 86 is the RTX 3090; use 8.9 -> 89 for an RTX 4090, and so on
export CUDA_CACHE_MAXSIZE=4294967296        # the driver keeps the compiled kernels in ~/.nv/ComputeCache: room, and a home folder the server's user can write
./target/release/helixloom gpu-check        # "Card 0 alone", "Card 1 alone" and "All cards together" must all pass
```

**This recipe was written from reading the code and candle's sources; it has not been run on a CUDA card.** It is what to
try first. If `gpu-check` fails on the newer card only, the message says which check and the reason: it is most likely a
`libcublas` older than 12.8 (the 5060 Ti is then left out of the pool, with a line on the standard error saying so, and the
3090 works alone), or a driver too old to compile the 3090's kernels for it.

### 1.3 The programs that the jobs call (optional, recommended)

```sh
sudo apt install mafft primer3 curl
```

IQ-TREE comes from <https://iqtree.org>, from conda (`conda install -c bioconda iqtree`) or from your package manager: the
server looks for `iqtree3`, `iqtree2` or `iqtree` on its path. A program that is missing is replaced by the built-in
method, and the Workstation page shows what the server found.

`helixloom tools` lists the programs the server can call, which are found and in which version. Where Homebrew, micromamba,
mamba or conda is installed, `helixloom tools install NAME` installs one with it: it prints the command first
(`micromamba install --yes --channel conda-forge --channel bioconda mafft`) and asks before it runs it (`--yes` answers).
conda and its kin install into the environment that is active, so activate the one the server runs in first. Nothing is ever
installed because a job asked for it.

### 1.4 Make a token and start the server

The token is the password that your computers send. Make one, keep it in a file only you can read, and start the server:

```sh
./target/release/helixloom token > token.txt
chmod 600 token.txt
./target/release/helixloom serve --token-file token.txt --gpus auto
```

It prints where it listens, the name of the workstation, and each GPU it will hand out. In a second shell on the
workstation, ask it how it is:

```sh
./target/release/helixloom status --server 127.0.0.1:7878 --token-file token.txt
```

Keep the token to yourself: anyone who has it can use the workstation. Use the one that `helixloom token` makes (it comes
from the system's random numbers and cannot be guessed): a token that a person invents can be guessed by a computer that
can reach the port, since the server does not slow down a wrong guess. The token is never put on the command line (a
command line can be read by other users with `ps`): it comes from a file, or from the environment variable
`HELIXLOOM_TOKEN` (which the server takes out of its environment as soon as it has read it, so that the programs it runs
never see it). A token shorter than 32 characters is refused.

### 1.5 The options of `serve`

| Option | Meaning |
| --- | --- |
| `--listen ADDRESS` | where to listen; default `127.0.0.1:7878`, this computer only (see step 2) |
| `--token-file FILE` | the file that holds the token (or set `HELIXLOOM_TOKEN`) |
| `--data-dir DIR` | where jobs and results are kept (`jobs`), and the answers of the steps of pipeline files (`flow-cache`, which the server never empties: `helixloom pipeline clear-cache DIR/flow-cache` does); default `./helixloom-data` |
| `--gpus auto`, `none`, or a list | which GPUs to hand out: `auto` asks `nvidia-smi`; `0,1,3` or `0-3` names them |
| `--cpu-jobs N` | jobs that need no GPU at once; default 2. A job may itself use many cores: choose with that in mind |
| `--max-upload-mb N` | the biggest request accepted; default 512 |
| `--keep-days N` | days a finished job's result is kept; default 7 |
| `--job-hours N` | hours one job may run before it is stopped; default 24 |
| `--fold-command "PROGRAM words"` | a program that folds a protein on a GPU (step 5) |
| `--fold-gpus N` | how many GPUs one folding job is given; default 1. Raise it only if your folding program uses several GPUs (it is shown all of them) |
| `--name NAME` | the name your computers see; default the workstation's own |
| `--web-dir DIR` | also hand out the web version found in `DIR` (what `scripts/build-web.sh` or `trunk build --release` makes, `docs/WEB.md`): a browser that opens `http://HOST:7878/` gets the app, and its Workstation page is already pointed at this server. A file with a packed copy beside it (`name.gz`, as the script makes) is sent packed to a browser that takes gzip: the program is a third of its size that way. The folder must hold only the build, because it is handed to anyone: the server refuses to start if a file in it holds the token, or if the token's file or the `--data-dir` is inside it. It is read once, when the server starts: restart it after a new build |
| `--allow-origin ORIGIN` | let the web version hosted at `ORIGIN` (`https://host` or `http://host:port`, nothing after it) call this server from there. Not needed for the page `--web-dir` hands out. May be given more than once |

### 1.6 How the GPUs are shared

Each GPU does one job at a time. A job that needs GPUs takes the lowest-numbered free ones (one, unless the person who sent
it chose more under *GPUs for one job* and the job can use them) and is shown those GPUs only (the server sets
`CUDA_VISIBLE_DEVICES` for that job, a list such as `0,1`, and numbers the GPUs as `nvidia-smi` does), so two jobs never meet
on a GPU. A job that wants several waits until that many are free. If the GPUs are busy, the next GPU job waits in the queue
in the order the jobs were sent, and none behind it overtakes it, so a job that wants four GPUs is not kept waiting for ever
by a stream of jobs that want one. Jobs that need no GPU go ahead on the processor in the meantime. What a job may be given is
cut to what it can use: PhyloCNN no more GPUs than it has trees, a structure no more than `--fold-gpus`, and no job more
than the workstation has. Every job runs in a process of its own: a job that crashes, runs too long (`--job-hours`) or is stopped from the
app cannot take the server down, and the server cleans up after it, including any program the job started. What one request
may hold is limited (5,000 sequences, 64 million letters of sequence in all, 20,000 trees), which keeps the memory of one job
to a few gigabytes; the service file (step 3) also caps the memory of the server and all its jobs together, so that a job
that is asked to do something enormous is stopped instead of slowing the whole machine.

What uses a GPU: PhyloCNN estimates of many trees (with the CUDA build; it divides the trees among the GPUs it is given, and
the estimates come back in the order of the trees), big tree jobs (step 1.7), and structure prediction with `--fold-command`.
The PhyloCNN networks are files of their authors' and are not part of the program: install them on the workstation too (a
`models` folder next to the server's settings, or the folder named by `HELIXLOOM_MODELS`; see `docs/MODELS.md`).
Alignments and pipelines run on the processor, and so do tree jobs too small to be worth a GPU. A tree job is given GPUs by its
size alone (the server never opens a GPU itself) and waits for them in the queue like any GPU job; on a workstation with no GPU
(`--gpus none`) it is not refused but runs on the processor.

### 1.7 Trees on the GPUs

A tree job runs in a process that is shown the GPUs it was given, and uses them when it is big enough, through **candle** (built into the
program) or **PyTorch** (Python workers, one per GPU; chosen when `HELIXLOOM_TORCH_PYTHON` names a Python, or with
`HELIXLOOM_ACCEL_BACKEND=torch|candle`). Both do the same work and give the same answers:

| What | How it uses the GPUs | How exact |
| --- | --- | --- |
| Distances between all pairs of sequences, and the resampled distances of a bootstrap | the columns are written as 0 and 1 and the three counts of each pair come from three matrix products; every card does its share of the rows, and a bootstrap sends the columns to the cards once and then only a weight for each column per resample | the counts are whole numbers a 32-bit float holds exactly, and the distance is worked out from them by the processor's own function: **equal to the processor's to the last digit** |
| The likelihood of a tree (maximum-likelihood search and the fit of its four models; the Bayesian sampler) | Felsenstein's pruning written as tensor operations over all the distinct columns; one big likelihood is cut into ranges of columns, one per card in proportion to its speed, and the four model fits and the sampler's chains are spread over the cards | 32-bit logarithms summed in 64 bits: **agrees to about seven digits, not to the last**; on GPUs a change of topology must raise the log-likelihood by at least a millionth of its size (0.001 at the least) to be taken as better, so that noise is not mistaken for a better tree |

One more thing goes to the GPUs with PyTorch (not with candle):

| What | How it uses the GPUs | How exact |
| --- | --- | --- |
| Neighbour-joining of thousands of sequences | the n by n distance matrix stays on one card in 64 bits, and every join is a few passes over it (the matrix is cut down as clusters leave) | **the processor's joins and tree**: tied pairs are taken earliest in the list first, as the processor takes them, but the row sums are added in another order, so a near-tie, one in the last digit, can be decided the other way (about one matrix in ten with many identical sequences and fractional distances, in the checks made): the tree has the same total length and differs only in how a polytomy is resolved. A card with no 64-bit floats (a Mac's GPU) leaves it to the processor. Not above 15,000 sequences. The processor's own joins are 1.8 times faster on one core since November 2026 and shared among its cores from 800 sequences on, to the last bit the same joins (4,000 sequences: 50 s before, 7 to 9 s on a 10-core Mac), so the cards have less to win here |

**What is not on the GPUs, and why.** MAFFT and IQ-TREE are other programs. UPGMA, species trees (ASTRAL), the conservation scores and
pipelines run on the processor; a pipeline's job is given no GPU. **Alignment is not GPU work in this program, and the part that was slow needed no GPU.**
The built-in aligner chooses a centre sequence, by how many 6-mers each shares with the others, and aligns every other sequence to it. The choice used to intersect
the sets of 6-mers of every pair of sequences (half a minute for 600 sequences of 1,500 bases); it now counts, for each 6-mer, how many sequences hold it, which
gives the same sums in one pass (`seq/msa.rs`, checked against the pairwise way on random, related, protein and degenerate sets). The alignments themselves (the dynamic
programming of the centre against each other sequence) took 10.4 ns for each cell on one core, shared among all cores; they take **1.3 ns now in a release build** (2.9 ns for one pair alone, 3.6 ns in a
debug or test build), because the program does eight pairs at once in vector instructions on the processor, with results identical to the old code to the last bit (`seq/align.rs`, checked against a copy of it on about 850 pairs):
a centre star of 200 sequences of 1,500 bases takes 0.13 s on a Mac with 10 cores (it took 0.70 s). **They are not on the GPUs unless you ask, on purpose.** Eager PyTorch,
along anti-diagonals with about 45 small tensor operations for each, took 20 ns for each cell at 64 pairs of 1,000 by 1,000 and 10 ns at 256 pairs on a Mac's GPU, where one core took 10.4 ns then (all eight, 1.3 ns) and takes 1.3 ns
now (all eight, about 0.17 ns); that is mostly the launching of operations, so a card would do better, perhaps twice the cores (an estimate from the memory traffic, not a measurement), and
only a fused kernel would beat the cores by much: it cannot be written and checked without a CUDA machine, and the processor's answer is now eight times nearer. **They can go to the PyTorch workers on request** (`HELIXLOOM_ACCEL=force`, or `HELIXLOOM_ACCEL_ALIGN_CELLS` set to the cells of all the tables together from which they go: `accel/torch/align.rs`): the tables are filled along their anti-diagonals for many pairs at once, with the processor's 32-bit sums and ties, the pairs shared among the GPUs of the job by speed, and the alignments come out **the same to the last bit** (checked against the old aligner on 1, 2 and 3 workers of the CPU device, and, in the worker's own check, on a Mac's GPU), with the processor taking over if a worker fails; but it is **slower than the processor wherever it has been measured**: 61 to 72 times on a worker of the CPU device and 413 to 1,008 times on this Mac's GPU (a busy machine, 32 to 64 pairs of 300 to 500 bases), because every diagonal is about sixty small operations. So nothing sends alignments there unless asked, and `helixloom gpu-check --bench` times them on your cards and says whether to set the variable (a fused kernel is the one thing that could win, and needs a CUDA machine). Sequences over about 2,800 bases, whose eight tables would take over 64 MB, are aligned one pair at a time (2.9 ns a cell). PhyloCNN, the one other user of the
GPUs, has its own path. When a card fails in the middle of a job, the processor finishes it and the job's details and log say which card and why.

**How big is big enough.** Below a size the processor is faster, because a GPU has to be fed and answered. The server decides
by size alone, and the job decides again from its own alignment (patterns times tips for a likelihood; sequences squared times
columns, for all resamples together, for distances). The starting values (10 million cells for a likelihood, 100 million for
distances) were measured on **one machine, a Mac's GPU**, so on your cards they are a place to start from. Find the real ones with
`helixloom gpu-check --bench` (a few minutes), which tests each card and then all together against the processor's answers, times both on growing
sizes and prints the variables to set.

| Variable | Meaning |
| --- | --- |
| `HELIXLOOM_TORCH_PYTHON=/path/to/python` | the Python (with PyTorch and NumPy) whose workers do the tensor work; with it set, PyTorch is the backend unless `HELIXLOOM_ACCEL_BACKEND=candle` |
| `HELIXLOOM_ACCEL_BACKEND=auto\|torch\|candle` | which library; `auto` (the start) is PyTorch where a Python is named, else candle |
| `HELIXLOOM_TORCH_DEVICES=cuda:0,cuda:1` | the devices PyTorch uses, instead of every card it finds (`cpu` tries it with no GPU; `mps` is a Mac's GPU) |
| `HELIXLOOM_ACCEL=off\|auto\|force` | `off`: processor only; `auto` (the start): GPUs for big work; `force`: GPUs whatever the size (to try them) |
| `HELIXLOOM_GPU_WEIGHTS=0.6,0.4` | the share of the work of each GPU, in `nvidia-smi` order; not set, each card's speed is measured when the job starts. Measured speed rates matrix products, and the likelihood is limited by memory speed, so for an RTX 3090 with an RTX 5060 Ti `2,1` (their memory bandwidths are about 936 and 448 GB/s) may suit it better: try both |
| `HELIXLOOM_ACCEL_ALIGN_CELLS` | the cells of all the tables of a set of pairwise alignments from which they go to the PyTorch workers: **infinite unless you set it**, because the processor was 60 to 1,000 times faster wherever this was measured (`gpu-check --bench` says for your cards) |
| `HELIXLOOM_ACCEL_LIKELIHOOD_CELLS`, `HELIXLOOM_ACCEL_DISTANCE_CELLS`, `HELIXLOOM_ACCEL_JOIN_CELLS` | the sizes above which the GPUs are used: patterns times tips; sequences squared times columns; sequences cubed (what `gpu-check --bench` suggests) |

Set them in the service file (step 3) for the server, or in the shell for `gpu-check`. The server numbers the GPUs as `nvidia-smi`
does (`CUDA_DEVICE_ORDER=PCI_BUS_ID`), and so does `gpu-check`; inside a job, the cards it was given are named `CUDA 0`, `CUDA 1`... counting from the first of
those, whatever their numbers on the machine.

**PyTorch or candle.** On the one GPU that has been measured, a Mac's (`helixloom gpu-check --bench --quick`, warm, 256 sequences by 20,000 columns
for distances and 64 tips by 100,000 patterns for the likelihood), against one core of the processor:

| | distances | likelihood of one tree |
| --- | --- | --- |
| candle (Metal) | 61 times faster | 0.8 times: *slower* than the processor |
| PyTorch (Metal, `mps`) | 95 times faster | 16 times faster |

candle's small reductions and its broadcast products are slow, and the likelihood launches about nine of them for every inner node; PyTorch does each
node with a matrix product and a fused reduction. That is one machine and one GPU: **how either does on an RTX 3090 and an RTX 5060 Ti is not known until
`gpu-check --bench` has run on them.** PyTorch costs seconds to start in each job process that uses it (importing it and opening the cards, in parallel for all
the cards), so its start for the size of distances that go to the GPUs is higher (5·10^9 cells, ten seconds or so of one core; the likelihood's, 10^7 cells, is the
same as for candle, since a search makes hundreds of evaluations). Each worker is a process of its own, one per GPU, and ends when the job's process does. What PyTorch prints, and the line that says which cards
were found and how the work is cut (`HelixLoom: PyTorch 2.x, CUDA 12.x (built for …): share of the likelihood …`), are in the job's log. The workers
run with `python -I`, which ignores the user's own packages (`pip install --user`): install PyTorch and NumPy in a virtual environment, as above. A worker
that takes more than five minutes to start, or more than half an hour over one request (a GPU that has stopped answering), is killed and the processor does the
work. A likelihood too short to keep two cards busy (under about 200,000 patterns for each) goes whole to the faster card: with a tree of a few hundred tips, the time is
the launching of kernels, not the card, and a second card would add nothing (an estimate, not measured on your cards). **Setting the length of a branch no longer tries the whole tree for every length** (since November 2026): the workers keep what lies beyond both ends of the branch and try each length with one small product, so a search for a tree makes about a tenth to a twentieth of the launches (setting 12 branches took 3.8, 9.7 and 14.4 times less on a worker of the CPU device, for trees of 10, 20 and 40 tips; the processor's search is 4.7 to 10 times faster by the same change); the lengths found agree with the whole tree's to the noise of 32 bits (a few hundredths of the logarithm on a flat peak). If one of two workers dies in the middle of a job, the other does the work it was given, and the job's report says so. Without `HELIXLOOM_TORCH_PYTHON`, `HELIXLOOM_ACCEL_BACKEND=torch` runs `python3`
from the server's `PATH`.

**What to expect, and what is not known (candle).** Tried only with candle's CPU and Metal devices on a Mac. There, the distances were
5 times faster on a GPU at 64 sequences by 20,000 columns and 33 to 38 times at 256 sequences, and the likelihood was *slower*
than the processor below about 6 million cells and 1.2 times faster at 6.4 million: candle launches about nine small operations
for every inner node of a tree, so it is limited by launching, not by arithmetic, and a CUDA card will win by more or less than
this; nobody knows how much until `gpu-check --bench` has run on it. The likelihood search also evaluates the whole tree about 20 times
for every branch length it adjusts; keeping the partial results (re-rooting at the branch) would cut that 10 to 20 times on the
GPUs *and* the processor, and is the next step for speed (`docs/ROADMAP.md`). Whether CUDA needs the one-thread-at-a-time
rule per card that Metal needed (several threads on one Metal device gave silently wrong numbers, found by the tests) has not been tried;
the program takes turns on every card whatever the card, which costs a little overlap at most. Each GPU holds a full copy of the alignment's columns (the
tips as 32-bit numbers, and for distances 28 bytes a column and sequence, 4 GB at most per card; a bigger one is refused and the processor does it, with a
line in the log), and the partial results of a likelihood take about 20 bytes per pattern and inner node alive at once, which for
a few million patterns of a few hundred tips can exceed the 16 GB card: its failure is then reported and the processor finishes.

## 2. Reaching the workstation from your computers

The server does not encrypt anything: the token and your sequences travel as they are. That is fine inside a network you
trust, and it is the reason the server listens on `127.0.0.1` (this computer only) until you say otherwise. **Never put the
port on the open internet.** There are three ways to connect safely, from the best to the simplest.

1. **An SSH tunnel** (nothing to open on the workstation; the server stays on `127.0.0.1`). On each of your computers, keep
   this running (Windows 10 and 11 have `ssh` in PowerShell too):

   ```sh
   ssh -N -L 7878:127.0.0.1:7878 you@workstation
   ```

   In the app, the address is then `127.0.0.1:7878`.

2. **A private network such as Tailscale or WireGuard.** Start the server on the workstation's address on that network
   (`--listen 100.64.0.5:7878`, say) and use that address in the app.

3. **Your own trusted network.** Start the server with `--listen 0.0.0.0:7878` and let only your computers in:

   ```sh
   sudo ufw allow from 192.168.1.0/24 to any port 7878 proto tcp     # your network's addresses
   ```

### 2.1 From a browser

The web version (`docs/WEB.md`) reaches the server through the browser's `fetch`, with the same token. The simplest way is to
let the server hand out the page itself, so that the page and the server have the same address and the browser has nothing
to allow:

```sh
scripts/build-web.sh                       # on the computer where you build (docs/WEB.md): makes dist/, with packed copies
                                           # (or: trunk build --release)
helixloom serve --listen 127.0.0.1:7878 --token-file token.txt --web-dir dist
```

Then open `http://127.0.0.1:7878/` (through the SSH tunnel of way 1, on the computer you work at), go to Workstation, type
the token and press Connect: the address is already there. The browser keeps the address and the token for that tab, so
loading the page again does not ask for them; *Forget* on the Workstation page removes them, and a new tab starts without
them. On a computer that others use, press *Forget* before you leave (`docs/WEB.md` says why closing the tab is not enough).

The page can also be **installed** like a program and opened without a network (`docs/WEB.md`, *Installed, and without a
network*). Browsers allow that only to a page that is opened over `https://` or from the computer itself, so through the
SSH tunnel (`http://127.0.0.1:7878/`) it is offered, and at `http://workstation.local:7878/` it is not.

Two things the browser insists on:

- **A page opened over `https://` cannot call an `http://` address** (browsers make an exception for `localhost` and
  `127.0.0.1`). The server does not speak TLS. Either open the page from the server (`http://`, as above), or put a proxy
  that has a certificate in front of the server (Caddy, nginx) and type its `https://` address into the page. Tell the proxy
  that jobs are large: nginx refuses a body over 1 MB unless `client_max_body_size` is at least `--max-upload-mb`, and its
  `proxy_read_timeout` and `proxy_send_timeout` must be longer than the longest upload or result takes to move.
- **A page that the server does not hand out is another origin**, and the browser asks the server first whether it may call
  it. Name that page when you start the server, exactly as the browser writes it, and the server answers: `--allow-origin
  https://lab.example.org`. Without it the browser refuses and the page says it cannot reach the workstation. The server never
  answers "anyone"; the token is still needed for every job. The address that such a page fills in is its own, not the
  workstation's: replace it before you type the token, or the token goes to the computer that hosts the page.

## 3. As a service that starts at boot (Linux, systemd)

`deploy/helixloom-server.service` is a ready-made unit; its first lines are the commands that make a user for it, a place
for the token and the service itself:

```sh
sudo useradd --system --create-home --home-dir /var/lib/helixloom --shell /usr/sbin/nologin helixloom
sudo mkdir -p /opt/helixloom && sudo cp target/release/helixloom /opt/helixloom/
sudo install -d -o helixloom -g helixloom -m 700 /etc/helixloom
sudo -u helixloom sh -c '/opt/helixloom/helixloom token > /etc/helixloom/token' && sudo chmod 600 /etc/helixloom/token
sudo cp deploy/helixloom-server.service /etc/systemd/system/
sudo systemctl daemon-reload && sudo systemctl enable --now helixloom-server
journalctl -u helixloom-server -f           # what the server says; each job's own log is in its folder
```

`sudo cat /etc/helixloom/token` shows the token to paste into the app. After an update of the program, copy the new
`helixloom` to `/opt/helixloom/` and `sudo systemctl restart helixloom-server`: results and waiting jobs are kept, and a job
that was running at that moment is marked as failed (run it again).

## 4. On each of your computers (Mac, Windows, Linux)

Build and run the app as usual (`cargo run --release` in the `app` folder). The computer needs no GPU, no MAFFT and no
IQ-TREE for this. Rust comes from <https://rustup.rs>; on a Mac it asks for Apple's command line tools, on Windows for the
Visual Studio C++ build tools, and on Linux the usual development packages for graphical programs may be needed.

1. Open **Workstation** in the rail (or Go ▸ Workstation, ⌥⌘9 on a Mac).
2. Type the workstation's **address** (`127.0.0.1:7878` with the tunnel above, or the address of step 2) and paste the
   **token**; press **Connect**. The page shows the workstation's name, its GPUs and which of them are busy, and the programs it
   found.
3. Under **Where jobs run** choose **Workstation** (every job goes there) or **Automatic** (the big ones do). The pages
   beside their buttons say where the job as you have set it up would run, and a job that goes there offers the
   workstation's programs (MAFFT, IQ-TREE, a GPU for PhyloCNN) instead of your computer's.
4. Work as usual. A job that goes there shows "…on the workstation…" while it runs, and its result appears where it always
   does. Choose **This computer** to compute on your own computer again. If the workstation cannot be reached, the run
   says so: it never quietly falls back to your computer, because a ten-hour job should not start on a laptop by mistake.

You can close the app while a job runs. The result waits on the workstation (7 days unless `--keep-days` says otherwise):
open the Workstation page later and press **Open** beside the job. **Remove** deletes a finished job and its result from the
workstation; **Stop** stops a running one.

The address and token are remembered in a small file that only you can read (on a Mac and Linux; on Windows it lies in
your own profile):

| System | File |
| --- | --- |
| macOS | `~/Library/Application Support/HelixLoom/workstation.txt` |
| Linux | `~/.config/helixloom/workstation.txt` (or under `$XDG_CONFIG_HOME`) |
| Windows | `%APPDATA%\HelixLoom\workstation.txt` |

**Forget** on the Workstation page clears it, and your choice of where jobs run with it. The file keeps the choice as
`compute=local`, `compute=workstation` or `compute=automatic` (a file from before the choice had three ways says `use=1`,
which is read as Workstation). Setting `HELIXLOOM_CONFIG_DIR` puts the file in another folder.

### A Mac app

`scripts/bundle-macos.sh` makes `target/bundle/HelixLoom.app` and a disk image of it (`HelixLoom-0.2.0.dmg`) from a release
build, with the program's own mark as its icon. It is signed ad hoc: it starts on the Mac that made it, and on another Mac
the first start needs a right-click and **Open**, because it is not signed with a developer identity and not notarized by
Apple. With an Apple Developer account, sign it with your identity instead (`codesign --force --options runtime --timestamp
--sign "Developer ID Application: …" HelixLoom.app`), have the disk image notarized (`xcrun notarytool submit … --wait`) and
staple the answer to the app (`xcrun stapler staple HelixLoom.app`); **these three were not tried here**, for want of an
account. A file double-clicked in the Finder does not open in the app yet: open it from the app or drop it on its window.

### Releases for end users

A release is made by pushing a tag: `git tag v0.2.1 && git push origin v0.2.1` (the version in `Cargo.toml` should be the tag's).
`.github/workflows/release.yml` then builds on GitHub's machines and puts on the Release of that tag
`HelixLoom-v0.2.1-macos-arm64.dmg` (Apple Silicon, signed ad hoc as above), `HelixLoom-v0.2.1-windows-x64-NOT-TRIED.zip` and
`SHA256SUMS`, and publishes the server's container image (below). **The Windows file is marked in its name and in the notes because
nobody has started HelixLoom on Windows**: its build may fail (the Mac disk image and the release go on without it), and a build that
comes out may still not start. It is not signed: Windows shows SmartScreen's warning, and a console window may open behind the
program, because the same program is also the command line. Signing it needs a certificate (SignPath Foundation signs open-source
programs for free; Azure Trusted Signing is cheap), and Apple's notarization an Apple Developer account: none is used. The programs
are built without the GPU engine. **The workflow ran once, for `v0.2.0` on 19 November 2026, and all four jobs came out**: the Mac disk image (9.7 MB), the Windows zip (9.7 MB) and the container image were built, and the Release holds the first two with `SHA256SUMS`. **None of the files has been started on a computer other than the one that made the Mac app, or on Windows at all**; the container image was built, not run. By hand, for the Mac file only:
`scripts/bundle-macos.sh`, then `gh release create v0.2.1 target/bundle/HelixLoom-0.2.1.dmg --title "HelixLoom 0.2.1"`.

### Checks on every change

`.github/workflows/ci.yml` runs on every push to `main`, every pull request and by hand: formatting, lints and the tests on macOS (the system
the program has been tried on), the same on Linux and Windows (**not tried by hand: those two report but do not turn a change red**; a job that is
green there for a while is the reason to take `continue-on-error` out), and the lints of the browser build. The GPU build is not tried there (a
CI machine has no GPU).

### In a container

The `Dockerfile` builds the server and the web page into one image: see its first lines for the commands. It computes on the
processor, listens on `0.0.0.0` inside the container (publish the port on `127.0.0.1` only, as the example does) and reads the
token from `/run/secrets/token`. The release workflow pushes it to `ghcr.io/linhduongtuan/helixloom-server`; a package of a personal
account starts private, so make it public in the package's settings if others are to pull it. **Built by the release workflow for `v0.2.0`, not run**: it
follows the steps of `docs/UBUNTU.md`, and the libraries it installs to run the program are a guess from what the build needs.

## 5. Folding proteins on the GPUs (optional)

Without a folding program, the Structure page's prediction goes to the ESM Atlas web service, asked from the workstation
with its `curl` (the sequence leaves your network, up to 400 residues). To fold on your own GPUs instead, give the server
a program with `--fold-command`. The server runs it as `PROGRAM input.fasta output.pdb` on the GPU the job was given,
waits for the PDB file and stops the program after `--job-hours`. The program is named by you, in the server's own options:
a request from the network can never name a program.

`deploy/esmfold_predict.py` is an example that uses ESMFold from the `transformers` library:

```sh
sudo mkdir -p /opt/helixloom && sudo python3 -m venv /opt/helixloom/esmfold-env
sudo /opt/helixloom/esmfold-env/bin/pip install torch transformers accelerate
sudo cp deploy/esmfold_predict.py /opt/helixloom/
./target/release/helixloom serve --token-file token.txt --gpus auto \
    --fold-command "/opt/helixloom/esmfold-env/bin/python /opt/helixloom/esmfold_predict.py"
```

The first run downloads the model's weights (several gigabytes) and needs a GPU with plenty of memory; long proteins need
more (lower `CHUNK_SIZE` in the script if memory runs out). Try the script by hand on one sequence first (its top lines show
how). The Workstation page then says "folds proteins on its GPUs", and the Structure page says that the sequence is folded
on a GPU of the workstation.

## 6. If something does not work

| What you see | What to do |
| --- | --- |
| Anything you want to report | **Help ▸ Copy Diagnostics** in the app copies the version, the system, the programs it found (with their paths) and the day and place in the code of the last fault: paste that into the report. It holds nothing of your work, no address and no token. |
| "The workstation refused the token" | The token differs from the one in the file the server was started with. Paste it again, with nothing before or after. |
| "Cannot reach 127.0.0.1:7878" | The SSH tunnel is not running, the server is not running (`systemctl status helixloom-server`), or a firewall is in the way. |
| "This workstation speaks version 2 of the protocol…" | The app and the server are from different versions: build both from the same folder. |
| The page lists no GPU | `nvidia-smi` must work for the user that runs the server; or name the GPUs yourself with `--gpus 0,1`. |
| A GPU job says "no GPU set up" | The server runs with `--gpus none`, or `--gpus auto` found none. (A big tree job is not refused for this: it runs on the processor.) |
| A big tree job did not use the GPUs | Look at the job's details and log: they say "the processor (the GPUs could not be used: …)" or which card failed. Otherwise the job was too small (`gpu-check --bench` finds the size), or the server's GPUs setting gave it none, or the program was built without `--features cuda`. |
| `gpu-check` says "PyTorch could not be started" or "No GPU can be used" with PyTorch | The Python named by `HELIXLOOM_TORCH_PYTHON` cannot `import torch`, or PyTorch sees no CUDA device: run `/path/to/python -c 'import torch; print(torch.cuda.device_count(), torch.cuda.get_arch_list())'`. The RTX 5060 Ti needs a PyTorch built for CUDA 12.8 or later (`sm_120` in that list), and a driver of 570 or later. PyTorch's own messages are in the job's log (`jobs/<id>/log.txt`). |
| `gpu-check` leaves a card out, or prints "…is left out of the pool" | The card could not run the speed test: the CUDA toolkit's libraries are older than the card needs (step 1.2), or the driver is. The other cards go on working. |
| A job was given fewer GPUs than "GPUs for one job" says | It can use no more (a structure needs `--fold-gpus` on the server; PhyloCNN needs at least as many trees as GPUs), or the workstation has fewer. The job's line on the Workstation page says which GPUs it got. |
| A job waits although a GPU is free | The job in front of it wants more GPUs than are free; jobs wait in the order they were sent. Choose fewer GPUs for one job, or wait for the job in front. |
| A PhyloCNN run on the GPU fails with "This computer has no GPU the program can use" | The program on the workstation was built without `--features cuda`, or the CUDA driver is not usable by it: build again as in step 1.2. |
| "…has no folding program configured" | Start the server with `--fold-command` (step 5), or let the app use the web service. |
| A job failed | The reason is on the job's line; the job's own log is `jobs/<id>/log.txt` in the server's data folder. |
| An upload is "larger than this server accepts" | Start the server with a bigger `--max-upload-mb`. |
| "A job may hold at most 5000 sequences", "…letters of sequence", "…20000 trees" | Split the work into several jobs: these limits keep one job's memory in bounds. |
| "The tree is nested more than 2000 levels deep" | The tree the workstation made (or a tree file) is a ladder of more than 2000 tips, which the app does not read: drawing it would take more room than the program has. Make the tree from fewer sequences (cluster the near-identical ones first). |
| "The queue is full", "The workstation keeps 2000 jobs already" | 1000 jobs wait, or 2000 are kept: Remove some on the Workstation page. |
| "Too many connections" | One address may have 8 connections at once, and the server 64: wait a moment, or close other programs that talk to it. |
| "The workstation is taking in other uploads" | Two of the largest uploads are being taken in: try again in a moment. |

`helixloom help` prints the commands of the program: `serve`, `token`, `status`.

## 7. What the server does and does not allow

- A request is one of nine kinds of job (alignment, tree, pipeline, PhyloCNN, structure, species tree, pipeline file, motif search, domain search) with its data. It is
  checked when it arrives. Nothing in it is used as a file name, a command or a program.
- The token is checked from the header before any body is read, so a stranger without the token cannot make the server
  take in a large upload.
- With `--web-dir`, the files of the page (and nothing else) are sent to anyone who asks, without a token: it is the program,
  not your data. The folder is read into memory when the server starts and a request can only look names up in that list,
  so no request can name a file of the workstation (`/../token.txt`, `%2e%2e` and the like are simply not in the list).
  Links, hidden files and names with unusual characters are left out. Every address under `/v1/` still needs the token (a
  folder called `v1` in the page's folder is refused). What is in the folder is public, so the server refuses to start if the
  token is in a file there, or the token's file or the jobs' folder is inside it: keep the web build in a folder of its own.
- The page is sent to at most 16 computers at a time (the others are told to try again), and each file within a minute or so
  (the 13 MB program), so that a computer that only pretends to read cannot hold the places that the person with the token
  needs for long. That is a limit and not a lock: a determined computer that can reach the port can still keep it busy (the
  next point).
- Anyone who can reach the port can learn whether an origin is in the `--allow-origin` list, by sending it as `Origin`: the
  answer to a refusal has to say so for the page to be able to read it.
- With `--allow-origin`, a page of exactly that origin may read the server's answers (a refusal too, so that it can say why).
  No other origin gets any permission, and there is no setting that allows every origin.
- Requests, headers and uploads have size limits; at most 64 connections are served at once, and at most 8 from one
  address; a connection must send its whole header within 5 seconds, however slowly, or it is dropped; an upload is given
  the time its size needs at 50 KB a second (the largest, about three hours), however slowly it comes; and only as many
  bytes of uploads are taken in at the same moment as two of the largest fit in memory (a further upload is told to try
  again). A job's request is deleted from the disk when the job is over, and at most 2,000 jobs are kept at once.
- A computer that can reach the port can still make it busy without the token, by opening connections (from many addresses
  at once, or as fast as it can). That is why the port is for computers you trust only: reach it through an SSH tunnel, or
  let only your own computers in with a firewall (step 2).
- The app is just as careful with what a server sends it: a reply may be only so large and take only so long, text from a
  server is cut to one short line, and a tree nested more than 2,000 levels deep is refused. Even so, an unencrypted
  connection can be tampered with by anyone on the way: another reason for step 2.
- Each job runs in a process of its own, in its own folder, with its own temporary folder, as the user that runs the server.
  Run the server as a user of its own, with no other rights (the systemd unit does).
- The data folder is readable by that user only. Results are removed after `--keep-days`.
- Nothing is encrypted. The token and the data are as private as the network they cross: see step 2.

## 8. What has been tried, and what has not

Tried, on a Mac, with the server and the app on the same computer: the whole path (the app sends an alignment, a tree, a
pipeline, a pipeline file over a folder of loci, PhyloCNN estimates, a species tree, a motif search, a domain search or a structure; the server queues it and runs it in a process of its own; the result comes
back and opens in its page), the token and the limits, results kept across a restart, stopping a running job and ending its
process, and the sharing of (pretend) GPUs, with a stand-in folding program that reports which GPUs it was shown: one GPU to
a job, then several GPUs to one job, and the order in which jobs wait for them. The sharing of one batch of PhyloCNN trees
over several devices was tried with the processor standing in for the devices (the estimates equal those of a single
device). These are automatic tests (`cargo test`, including `tests/remote.rs`, which starts the real program).

**Not tried:** a real NVIDIA GPU or the `cuda` build (this includes the multi-GPU alignment, tree and bootstrap code of step 1.7, which was
tried with candle's CPU device in every pool shape, with several threads at once, and on a Mac's GPU; whether a 3090 and a 5060 Ti in one process work, and
what the CUDA build recipe of step 1.2 does, is for `helixloom gpu-check` to say on your machine; the same for the PyTorch workers: they were run on a Mac's CPU and GPU, never on CUDA), PhyloCNN on several real GPUs at once (the code that opens GPU number 0, 1,
2… of the CUDA build has never been compiled or run, and how much faster it is than one GPU is not measured), the systemd unit
and the ESMFold example (no Linux GPU machine was at hand), running the server on Windows, building and running the app on Windows and Linux, an SSH
tunnel or Tailscale (they are standard tools, but they were not run), and uploads of hundreds of megabytes. Expect to try each
of these once on your own machines; if one fails, the message on the Workstation page or in `journalctl` says where.

## 9. For scripts: the messages

The server speaks plain HTTP/1.1, one request per connection, with `Authorization: Bearer <token>`. A body is a
sequence of sections, each `@name length`, a line break, exactly `length` bytes and a line break.

| Request | Meaning |
| --- | --- |
| `GET /v1/status` | name, version, protocol, running and waiting jobs, GPUs and what they do, the programs found |
| `GET /v1/jobs` | every job, one line each: id, kind, status, times, the GPUs it was given (`0,1`), how many GPUs it waits for, label, error |
| `POST /v1/jobs` | send a job: the sections of its kind and a `job` section (`kind=align`, `label=…`, and optionally `gpus=N`: how many GPUs the sender asks for, 0 for every GPU, one when it says nothing); answers `201` with its `id` |
| `GET /v1/jobs/ID` | one job's line |
| `GET /v1/jobs/ID/result` | the result of a finished job (`409` while it is not ready or if it failed) |
| `DELETE /v1/jobs/ID` | remove a job; stop it if it is running |

`src/remote/tasks.rs` lists the sections of each kind of job and of each result.

## 10. Adding another analysis later (motifs and domains, epitopes and paratopes)

A new analysis gets the choice of where it runs, the upload, the queue, the GPU and the way back for the price of four
small pieces. The server does not change for it, and the Workstation page needs one line: the arm of `open_remote_result`
that shows a fetched result in the page it belongs to (the compiler asks for it).

1. **A kind of job** in `src/remote/tasks.rs`: a `Kind`, a `Request` that holds its input (sequences, a structure, the
   settings) and a `Response` that holds what comes back, each with its sections in `encode` and `decode`. Check every
   number and every length on `decode`, as the other kinds do: a request comes from the network.
2. **What runs** in `src/remote/run.rs`, in the arm of `execute` for the new request. The same function runs the job on
   the workstation and, when you choose This computer, on your own computer, so the two cannot differ. A model that is a
   program of its own (a Python model on a GPU, say) is called like the folding program is: its command comes from the
   server's own options, never from a request.
3. **How big it is**: say what makes a job of this kind slow, either in a `Weight` constructor next to the others in
   `src/remote/weight.rs` (sequences, letters, structures, a GPU) or by adding a field and a limit to `Weight`. This is
   all Automatic needs to decide; `needs_gpu` on the request tells the server to give the job a GPU, and `most_gpus` (in
   `src/remote/tasks.rs`) says how many it can use at once: one, unless the model spreads its work over several GPUs the way
   PhyloCNN does (`src/seq/phylocnn/gpu.rs`, `run_batch_across`).
4. **A page** that builds the request and asks `platform.place(&weight)` where it runs: `Place::Here` runs `execute` on a
   worker thread, `Place::There(remote)` calls `remote.run(&request, label)`. Drawing `compute_picker` beside the button
   gives the page the same choice as the others, and `platform.programs_at(&place)` gives the programs of the computer it
   goes to. `Alignment` is the shortest page to copy.

The species-tree job (`src/remote/species.rs`, the `species` kind, `src/ui/windows/species.rs`) is the newest worked example
of exactly these four pieces: its `Weight` is the leaves of its gene trees, its function is `species::run`, and the same
function is `helixloom astral` on the command line. The motif search (`src/remote/motif.rs`, the `motif` kind,
`src/ui/windows/motifs.rs`, `helixloom motif`) is the same shape with a smaller `Weight` (the bases of its sequences). The domain
search (`src/remote/domains.rs`, the `domains` kind, `src/ui/windows/domains.rs`) runs another program, HMMER's `hmmscan`, which the
server finds on its own PATH and reads Pfam from its own models folder (`docs/DOMAINS.md`); the Workstation page lists HMMER among
the programs it found.

The analyses you have in mind fit this directly: conserved motifs and domains are a function of an alignment or a set of
sequences (a few sections in, a table or annotations out); an epitope and paratope interaction is a function of two
sequences or two structures (an antigen and an antibody) with a score per residue pair as the result, heavy enough that it
will usually want the GPU. Heavy models belong on the workstation, and the choice of Automatic keeps the light ones on the
laptop.
