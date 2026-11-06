# Running the heavy work on a workstation with GPUs

HelixLoom can send its heavy jobs to another computer, a workstation with several GPUs, and bring the results back to the
Mac, Windows or Linux computer you work on. The workstation runs one program, `helixloom serve`. Your own computers run
the HelixLoom app as usual and have one extra page, **Workstation**, where you connect to it.

```
  your Mac / Windows PC / Linux PC                      the workstation (several GPUs)
  ┌─────────────────────────────┐   jobs go over the    ┌────────────────────────────────┐
  │ HelixLoom app                │   network (a private  │ helixloom serve                 │
  │  Alignment, Phylogeny,      │ ────────────────────▶ │  queue ─▶ one process per job  │
  │  Pipelines, PhyloCNN,       │   network or an SSH   │  each GPU job sees one GPU     │
  │  Structure                  │   tunnel)             │  results kept on disk          │
  │  Workstation page: connect, │ ◀──────────────────── │  (7 days, and across restarts) │
  │  jobs list, Open results    │   results come back   └────────────────────────────────┘
  └─────────────────────────────┘
```

What is sent, and what is not:

| In the app | Where it runs when the job goes to the workstation |
| --- | --- |
| Alignment (built-in or MAFFT) | the workstation's processor |
| Phylogeny: neighbour-joining, UPGMA, maximum likelihood, Bayesian, bootstrap | the workstation's processor |
| Pipelines (align, tree, primers) | the workstation's processor, with its MAFFT, IQ-TREE and primer3 |
| PhyloCNN epidemic estimates | the workstation, on one GPU when "many trees run on the GPU" is chosen |
| Structure prediction | the workstation's GPU with your folding program, or the ESM Atlas web service asked from the workstation |
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

## 1. On the workstation

You need a Linux workstation (Ubuntu is assumed below), with the NVIDIA driver installed: `nvidia-smi` should list your
GPUs. Nothing else is needed for jobs that run on the processor. For the GPU parts you also need the CUDA toolkit (to build
the GPU version, step 1.2) and Python (for the folding program, step 5).

### 1.1 Get the program there

Copy the `app` folder of this project to the workstation (not its `target` folder, which holds this computer's build), then
install Rust once and build:

```sh
rsync -a --exclude target app/ you@workstation:HelixLoom-app/       # from your Mac or Linux computer
ssh you@workstation
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh     # follow the questions, then open a new shell
sudo apt install build-essential pkg-config        # more libraries if the build says one is missing
cd HelixLoom-app
```

### 1.2 Build it

```sh
cargo build --release                       # processor only
cargo build --release --features cuda       # also lets PhyloCNN use the NVIDIA GPUs (needs the CUDA toolkit)
```

The program is then `target/release/helixloom`. If the `cuda` build stops with a message about `nvcc` or CUDA, the toolkit
is missing or not on the path: `nvcc --version` must work first.

### 1.3 The programs that the jobs call (optional, recommended)

```sh
sudo apt install mafft primer3 curl
```

IQ-TREE comes from <https://iqtree.org>, from conda (`conda install -c bioconda iqtree`) or from your package manager: the
server looks for `iqtree3`, `iqtree2` or `iqtree` on its path. A program that is missing is replaced by the built-in
method, and the Workstation page shows what the server found.

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
| `--data-dir DIR` | where jobs and results are kept; default `./helixloom-data` |
| `--gpus auto`, `none`, or a list | which GPUs to hand out: `auto` asks `nvidia-smi`; `0,1,3` or `0-3` names them |
| `--cpu-jobs N` | jobs that need no GPU at once; default 2. A job may itself use many cores: choose with that in mind |
| `--max-upload-mb N` | the biggest request accepted; default 512 |
| `--keep-days N` | days a finished job's result is kept; default 7 |
| `--job-hours N` | hours one job may run before it is stopped; default 24 |
| `--fold-command "PROGRAM words"` | a program that folds a protein on a GPU (step 5) |
| `--name NAME` | the name your computers see; default the workstation's own |

### 1.6 How the GPUs are shared

Each GPU does one job at a time. A job that needs a GPU takes the lowest-numbered free one and is shown that GPU only (the
server sets `CUDA_VISIBLE_DEVICES` for that job, and numbers the GPUs as `nvidia-smi` does), so two jobs never meet on a GPU.
If every GPU is busy, the next GPU job waits in the queue, and jobs that need no GPU go ahead on the processor in the
meantime. Every job runs in a process of its own: a job that crashes, runs too long (`--job-hours`) or is stopped from the
app cannot take the server down, and the server cleans up after it, including any program the job started. What one request
may hold is limited (5,000 sequences, 64 million letters of sequence in all, 20,000 trees), which keeps the memory of one job
to a few gigabytes; the service file (step 3) also caps the memory of the server and all its jobs together, so that a job
that is asked to do something enormous is stopped instead of slowing the whole machine.

What uses a GPU: PhyloCNN estimates of many trees (with the CUDA build), and structure prediction with `--fold-command`.
The PhyloCNN networks are files of their authors' and are not part of the program: install them on the workstation too (a
`models` folder next to the server's settings, or the folder named by `HELIXLOOM_MODELS`; see `docs/MODELS.md`).
Alignments, trees and pipelines run on the processor.

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
| "The workstation refused the token" | The token differs from the one in the file the server was started with. Paste it again, with nothing before or after. |
| "Cannot reach 127.0.0.1:7878" | The SSH tunnel is not running, the server is not running (`systemctl status helixloom-server`), or a firewall is in the way. |
| "This workstation speaks version 2 of the protocol…" | The app and the server are from different versions: build both from the same folder. |
| The page lists no GPU | `nvidia-smi` must work for the user that runs the server; or name the GPUs yourself with `--gpus 0,1`. |
| A GPU job says "no GPU set up" | The server runs with `--gpus none`, or `--gpus auto` found none. |
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

- A request is one of five kinds of job (alignment, tree, pipeline, PhyloCNN, structure) with its data. It is checked when
  it arrives. Nothing in it is used as a file name, a command or a program.
- The token is checked from the header before any body is read, so a stranger without the token cannot make the server
  take in a large upload.
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
pipeline, PhyloCNN estimates or a structure; the server queues it and runs it in a process of its own; the result comes
back and opens in its page), the token and the limits, results kept across a restart, stopping a running job and ending its
process, and the sharing of two (pretend) GPUs, with a stand-in folding program that reports which GPU it was shown. These
are automatic tests (`cargo test`, including `tests/remote.rs`, which starts the real program).

**Not tried:** a real NVIDIA GPU or the `cuda` build, several real GPUs at once, the systemd unit and the ESMFold example (no
Linux GPU machine was at hand), running the server on Windows, building and running the app on Windows and Linux, an SSH
tunnel or Tailscale (they are standard tools, but they were not run), and uploads of hundreds of megabytes. Expect to try each
of these once on your own machines; if one fails, the message on the Workstation page or in `journalctl` says where.

## 9. For scripts: the messages

The server speaks plain HTTP/1.1, one request per connection, with `Authorization: Bearer <token>`. A body is a
sequence of sections, each `@name length`, a line break, exactly `length` bytes and a line break.

| Request | Meaning |
| --- | --- |
| `GET /v1/status` | name, version, protocol, running and waiting jobs, GPUs and what they do, the programs found |
| `GET /v1/jobs` | every job, one line each: id, kind, status, times, GPU, label, error |
| `POST /v1/jobs` | send a job: the sections of its kind and a `job` section (`kind=align`, `label=…`); answers `201` with its `id` |
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
   all Automatic needs to decide; `needs_gpu` on the request tells the server to give the job one GPU.
4. **A page** that builds the request and asks `platform.place(&weight)` where it runs: `Place::Here` runs `execute` on a
   worker thread, `Place::There(remote)` calls `remote.run(&request, label)`. Drawing `compute_picker` beside the button
   gives the page the same choice as the others, and `platform.programs_at(&place)` gives the programs of the computer it
   goes to. `Alignment` is the shortest page to copy.

The analyses you have in mind fit this directly: conserved motifs and domains are a function of an alignment or a set of
sequences (a few sections in, a table or annotations out); an epitope and paratope interaction is a function of two
sequences or two structures (an antigen and an antibody) with a score per residue pair as the result, heavy enough that it
will usually want the GPU. Heavy models belong on the workstation, and the choice of Automatic keeps the light ones on the
laptop.
