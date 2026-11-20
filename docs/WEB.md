# HelixLoom in a browser

The same code runs as a desktop app and, compiled to WebAssembly (`wasm32-unknown-unknown`), in a browser tab. The desktop
build is unchanged. The page is one canvas that eframe draws on, so it needs a browser with WebGL 2 (or WebGPU).

## Build and run

You need the Rust `wasm32-unknown-unknown` target and [trunk](https://trunkrs.dev) (the Rust web bundler, `cargo install --locked trunk`;
it downloads `wasm-bindgen` and `wasm-opt` on first use). Take care that the `trunk` on your PATH is that one: the
`trunk` of trunk.io is a different program (`trunk --version` of the right one starts with 0.). Take care of the Rust on
your PATH too: a Rust from Homebrew has no `wasm32` standard library, whatever `rustup target list` says; put rustup's
toolchain first (`export PATH="$HOME/.rustup/toolchains/stable-aarch64-apple-darwin/bin:$PATH"`, with your own toolchain's
name) before `trunk build`, or the build stops with "can't find crate for `core`".

```sh
rustup target add wasm32-unknown-unknown
trunk build --release          # makes dist/ (index.html is the page; Trunk.toml sets the folder and relative addresses)
cd dist && python3 -m http.server 8080    # any static file server will do; then open http://localhost:8080/
```

**Or, without trunk: `scripts/build-web.sh [FOLDER]`** (November 2026). It makes the same folder (`dist/` when none is named)
with what is on the computer and downloads nothing: it puts rustup's toolchain first by itself, builds the program for the
browser with the folders of the building computer kept out of it, lets `wasm-bindgen` write the script (the one on the PATH, in cargo's bin folder,
or the one trunk keeps in its cache, in the very version `Cargo.lock` names; it says what to install when there is none.
**That cache is not a safe place:** on 17 November 2026 the system had cleaned it away on the Mac this was written on, and
the page could not be rebuilt until the tool was fetched again (the release's own file for the Mac, 8.0 MB, checked against
the checksum published beside it, and put in cargo's bin folder, where it stays; `cargo install --locked wasm-bindgen-cli
--version …` puts one there too),
names each file after what is in it, and puts a **packed copy** (`.gz`) beside each file. `helixloom serve --web-dir` sends the
packed copy to a browser that takes gzip (nginx does with `gzip_static on`): the program is 16.4 MB and its packed copy
5.4 MB. The server sends a copy only when it is the copy of the file beside it: it compares the checksum and the length that a
gzip file ends with against the file's own (it has nothing to unpack with), so a copy left over from another build is not
sent in the place of a new file (the page, `index.html`, is as long in every build, so the length alone would not tell). It does not run `wasm-opt`, which trunk does: that made the file 12.3 MB here, and 5.1 MB packed, a gain of a
sixteenth over the wire for half a minute of building. Made and tried here on a Mac: the folder it makes was opened in a
browser from a plain file server and from `helixloom serve --web-dir`, which sent the packed program (checked with `curl`:
unpacked, it is the file byte for byte). **Not tried:** the script on Linux (it uses `shasum`, or `sha256sum` where only that
is there, and stops when neither is) or in the Windows shells.

`index.html` asks trunk for `--no-default-features`: the `gpu` feature (candle) does not build for the browser, and the code
that uses it is switched off for `wasm32` whatever the features say.

Without trunk: `cargo build --release --target wasm32-unknown-unknown --no-default-features`, then `wasm-bindgen --target web`
on `target/wasm32-unknown-unknown/release/helixloom.wasm`.

**Host the files in `dist/`**: `index.html`, one `helixloom-<hash>.js`, one `helixloom-<hash>_bg.wasm` (16.4 MB in November
2026; serve it compressed, 5.4 MB over the wire: the script above makes the packed copies), and one `style-<hash>.css`. The addresses are relative, so the folder can sit under
any path. The server must send `.wasm` as `application/wasm`; no special headers are needed (no threads, no
`SharedArrayBuffer`).

Without trunk (`cargo build --release --target wasm32-unknown-unknown --no-default-features`, then `wasm-bindgen --target web
--no-typescript --out-dir DIR --out-name helixloom` and an `index.html` that imports `./helixloom.js` and calls `init()`), the
page also runs: that is how it was tried here, with the `wasm-bindgen` that trunk keeps in its cache. `wasm-opt` (which trunk
runs) makes the `.wasm` file smaller.

Address options (like the command line's): `?view=trees` opens a workspace with its example, or any screen `--view` knows
(`?view=phylogeny`, `?view=pipelines`, ...); `?dark` starts in the dark theme (the browser's own setting is followed otherwise).

**Before hosting it publicly**: the wasm file contains the paths of the machine that built it (in panic messages). Build with
`RUSTFLAGS="--remap-path-prefix=$HOME=~"` so that none of them is in it.

## What works in the browser, and what does not yet

Checked in a real browser engine (headless Brave, WebGL 2): the page starts without a panic and draws Home and the Alignment,
Phylogeny, Pipelines (the example ran align, tree and primer design inside the page), Workstation, Structure, Trees and
Signals pages, in light and dark. **Written, compiled and tested natively, but not yet clicked through in a browser**: choosing a
file, dropping a file, and the downloads. Please try those first.

| | Browser |
|---|---|
| Plasmid editor: maps, enzymes, digest, ORFs, features, find, PCR, primer3 port, contig assembler | works |
| Alignment, Phylogeny (NJ, UPGMA, ML, Bayesian), Trees, Structure viewer, Signals (traces), calculators, read tools | work, on the page's own thread |
| File → Open, drag and drop, and the Open buttons of the pages (names and bytes go to the same routing as a path does) | works |
| Saving: sequences (GenBank, FASTA), alignment FASTA, Newick, tree SVG, PDB, SCF, cleaned FASTQ, contigs | works, as a browser download |
| Long jobs (likelihood, Bayesian, many PhyloCNN trees) | run inline on this computer: the page stalls while they run (a Bayesian tree of three short sequences, 400,000 generations, froze it for 0.6 s). Sent to the workstation they do not: the page stays responsive (its timer never slipped by more than 64 ms) |
| Workstation (jobs on a GPU server): connect, status, job list, send alignments, trees, species trees, pipelines, PhyloCNN and structure predictions, fetch and open results, stop or remove a job | works, through `fetch` (see below) |
| MAFFT, IQ-TREE, HMMER, primer3, fastp, curl, Snakemake, Nextflow | not available on this computer (the built-in methods are used); online structure prediction needs curl, so it is not available. **Sent to the workstation they run there**: a pipeline from the page used the workstation's MAFFT and IQ-TREE |
| Save figure (SVG, PNG), the report of the results (HTML), Copy Diagnostics, How Each Analysis Was Checked | compiled for the browser and tested natively; **not yet clicked through in a browser** (a figure and a report are saved as a download; a large PNG is drawn on the page's own thread) |
| Programs HelixLoom Calls, installing a program | not available: a page starts no program, and the window says so |
| Save workflow files (a folder), install PhyloCNN networks (a folder), PhyloCNN networks (model files) | not available: the page says so |
| The appearance, the places of the tool windows, that the tour was seen, **the open sequences** and **what the Alignment, Trees, Structure and Pipelines pages show** kept between visits | kept, in the browser's own storage for the site (`localStorage`) since November 2026, the pages' results since the 18th in its database (IndexedDB) since the 19th, each when there is room for it: see *What the page keeps* below |
| The workstation's address and token kept between visits | kept **for the tab**, since 17 November 2026: loading the page again does not ask for them, and they are not kept beyond the tab (the owner's decision; see *The connection, kept for the tab* below) |
| Recent files | not kept: a file that was opened in a page has no path to open it by again |
| GPU engine (candle) | not in the browser |

## What the page keeps

A page has no settings folder, so until November 2026 the web version kept nothing. It now keeps, in the storage the browser
gives the page's site (`localStorage`, under names that begin with `helixloom/`), what the desktop app keeps in its settings
folder, through the same code (`src/ui/session.rs` writes to a folder or to that storage: `src/ui/session/medium.rs`,
`src/ui/session/browser.rs`):

* **the open sequences**, as they are a moment after each change (two seconds of quiet, or fifteen of steady typing), and
  at once when the page is left or hidden (a tab that is closed or reloaded draws no more frames). At the next visit the page
  asks *Reopen your last session?*, as the desktop app does. The texts of the sequences are written first and the list of
  what is open last, so a page that is closed in the middle of a save comes back with the session before it;
* **what the pages show** (18 November 2026): the Alignment page (its list of sequences with their ticks, its choices and
  its alignment), the tree that the Trees and Phylogeny pages show, the structure, and the runs of the Pipelines page
  (the fixed route's alignment, tree and primers with its log and what was ticked; and the report of a pipeline file
  with its pipeline), each with the record of what made it (*How this was made*). They are written with the open sequences and offered with them: the question names them
  (*Also kept: the Alignment page, the tree and the structure*), and with no sequence open they are a session by
  themselves (*No sequence was open last time. Kept from it: …*). *Reopen* puts each back on its page; one that cannot be
  read is said and the rest comes back. A result that has not changed is not written again. **A result is kept when there
  is room for it, and the sequences never wait for one:** a result that does not fit is left out, and the page says so
  once, in red (*Could not keep the structure for next time: there is no room for it here, or it could not be written*);
  it is tried again at every later save, so it is kept as soon as there is room; a new result takes the room of the one it
  replaces; and when a sequence or the list of what is open cannot be written while results are kept, the results are
  removed to make room and the save is tried once more without them. **A run of the Pipelines page** is kept in the text
  it already has for the workstation, and read back with the care an answer from the network is read with: a name in it
  comes back on one line, without spaces at its ends, and a run that this text cannot hold (a file with no name, more
  than 20,000 files) is not put back. The folder of a pipeline file's run and the folder of its cache are not kept (the
  desktop app has them, a page has none): a path is not taken from what was kept and nothing is listed when a session
  is put back, so they are chosen again to run again. A kept run is not put over a run that is going or has been asked
  for. **Not kept:** what the tool windows show (PhyloCNN, species trees, domains, motifs), a trace of the Signals
  page, and how a page is scrolled, turned, coloured or zoomed. **Where the results are kept (19 November 2026).** The open sequences, the list of what is open, the appearance and the rest stay in
  `localStorage` (about five megabytes for a site in the Chromium that was tried: the example alignment takes 8 thousand letters of it,
  ubiquitin 80 thousand and an antibody with its antigen 320 thousand). That table takes a write at the moment it is asked for, so a
  page that is closed that instant keeps them. **The results go to the browser's database (IndexedDB) instead**, which holds hundreds
  of megabytes. A database answers later than a write is asked for, and the code that keeps things reads and writes at once, so: the
  page opens the database and reads all of it **before it starts** (6 megabytes added about a tenth of a second to the start, in
  Brave), holds it in memory (at most 160 megabytes of results; one that would not fit is refused as a full `localStorage` refuses),
  and hands every change to the database one at a time, in the order it was made, a change for a name taking the place of those that
  still wait for it. A result whose write has not finished when the page is closed is simply not there at the next visit; the
  sequences are, because they are in `localStorage`. More: a change that the database **refuses** is not kept, the page says so once
  for each reason (*Could not keep what the pages show for next time: the browser's database refused it (QuotaExceededError)*), and
  what the database may still hold under that name is taken away with the next change; the list of the session gives the **length and a
  checksum of each result**, and a text that is not the one listed (what a refused write leaves behind, or another tab wrote) is left
  out, so that it can never come back as the result of another session (a list of an earlier version, which has none, is read as
  before); results that an earlier version kept in `localStorage` are still offered and put back, and move to the database at the next
  save; putting the last session aside moves each result inside the table, with no second copy; turning the remembering off removes
  the whole table of results, also what another tab put there; and a browser with no database, or one that does not open it in 4
  seconds or read it in 20 (the page then starts without it), keeps the results in `localStorage` as before. **Limits:** a database that
  stops answering without an error holds the results back until the page is loaded again (the sequences are not affected); results
  that another tab wrote after this page opened are not offered by it until it is loaded again; the start takes longer the more
  results there are; projects with names and links to results are not built;
* **the appearance** (theme, size, density, accent, contrast, colours): a choice is kept *in this browser*, and the window
  says so. Until a theme is chosen the page follows the browser's own light or dark setting (it used to be given that setting
  as if it had been asked for in the address, which would have kept every dark-mode visitor from keeping anything);
* **where the tool windows were left**, and **that the first-run tour was seen** (the tour is shown at a first visit now).

**What it means for the person:** this is kept by the browser, on that computer, for that site: another browser, a private
window or another address of the same server starts empty, and clearing the site's data removes it. The sequences and the
results of the pages are there
as plain text, as they are in the settings folder of the desktop app: anyone who uses that browser after you can read them,
and is offered them (*Reopen your last session?*; *Start fresh* removes them). **On a shared computer, turn the remembering
off: File ▸ Remember My Work for Next Time** (17 November 2026; the command
palette finds it as *Remember my work*). Turning it off removes what was kept of the work, first, and then notes the choice:
the session (the sequences and the results of the pages), one that is on offer, one that was put aside, one of an older
version that could not be read, and (in the desktop app) the list of recent files. Nothing more of the work is kept from then on, in that browser, until it is turned on again;
another tab of the page shows the switch off at once and stops keeping the next time it has something to keep, and says so.
While the last session is on offer the command asks for an answer to that first (*Start fresh* asks before it removes work
that was not saved, and the switch must not go round that question). When the storage is so full that the choice cannot be
noted, what was kept is removed all the same and the page says that it will remember again at the next visit.
**What the switch does not cover**, and its hint in the menu says so: files you saved yourself; in the desktop app the records
fetched from a database (kept in `downloads` in the settings folder, to save fetching them twice) and the workstation's
address and token (the *Forget* button of the Workstation page removes those); the appearance, the places of the windows
and the mark of the tour, which say nothing of the work. Two tabs that act in the same thousandth of a second can still
cross (one writes its session just as the other removes it): what is left then is removed, not offered, the next time the
page is opened with the remembering off.
The browser gives a site a few
megabytes: plasmids fit many times over; a set of sequences that does not fit is not half kept (the list of what is open is
written last), and the page says once that it is not keeping the work. A result of a page that does not fit is left out by
itself (above): the sequences are kept all the same.

**Two tabs.** The desktop app lets only one copy keep the session, by a lock file. A page has no file to lock and cannot tell
whether another tab is still open, so the rule is the other way round: **the tab that was opened last keeps the session**
(it writes a mark when it starts), and a tab that was open before it finds another mark than its own the next time it has
something to keep, writes nothing more, and says *Another tab of HelixLoom keeps the session now, so what you change in this
tab is not kept: save what you need*. The new tab is offered what the older one had kept. So two tabs never write into one
session; what the older tab changes after the newer one opened is in that tab only. (The mark is looked at just before a
save, not held through it: two tabs that save in the same few thousandths of a second as one of them is opened could still
cross. That was not seen, and cannot be ruled out.) **Tried in a real browser** (the built-in one of the tool this was written with, on a Mac): the tour at a first visit
and not at the second; a theme chosen, the page reloaded, and the theme still there though the browser is set to dark; three
bases typed into the example plasmid, the page reloaded, the session offered and reopened with the bases in it; and, 0.6
seconds after another edit, the storage still held the old sequence until the page was told it was being left, and held the
new one at once after. With **two tabs** open on the page: an edit in the first was kept; the second, opened after it, was offered that work; and
an edit made in the first after that was not written, and the first tab said why (that message was then made shorter: it did
not fit the notice at the bottom, and a sign in it that the font lacks was taken out; the shorter one was not looked at again
in a browser). **Not tried:** Firefox and Safari; a storage that is really full or switched off (tested with a stand-in
only); a phone.

**Tried again on 17 November 2026, with the page built from the code as it then was** (the same browser, the page served
on this Mac): the frame of the pages in the dark theme (Alignment); in a window of 853 by 533 points, Trees with its
controls, its tree across the whole width and the card that says what it is under the tree, and Structure with its viewer
drawn under three rows of choices; and the switch: three bases typed, the storage held the session; *Remember My Work for
Next Time* clicked in the File menu, and the storage held the note of the choice and no session, the page said so in a line
that fits its notice; the page loaded again, nothing offered, the switch still off, three more bases typed and the page
told it was being left, and still no session; the switch clicked again, and the session was back. With a second tab opened
on the page, an edit in the first showed *Another tab of HelixLoom keeps the session now…* in one line of the notice, which
is the shorter message that had not been looked at. **Not tried:** the same in Firefox and Safari; turning the switch off
in one tab while another keeps the session (tested with the stand-in only).

**Tried on 18 November 2026: the results of the pages** (the page built from the code of that day, served on this Mac).
In the built-in browser: the example alignment (6 sequences, 600 columns), a neighbour-joining tree estimated from it, and
the example structure. The storage then held `1-alignment.result` (8,030 letters), `3-tree.result` (896),
`5-structure.result` (79,580) and the list that names them, and no sequence (the example plasmid, untouched, is not kept);
the alignment was written once, and not again when the tree and the structure came. The page was loaded again and asked
*Reopen your last session? No sequence was open last time. Kept from it: the Alignment page, the tree and the structure.*
After *Reopen* the Structure page was in front, as it had been, with ubiquitin drawn; the Alignment page had its six
sequences, their ticks, the alignment and *How this was made*; the Phylogeny page had the tree with its bootstrap labels
and its record; and the notice said *Put back the Alignment page, the tree and the structure from the last session*.
**A storage that is really full** could not be had in that browser (it took 43 megabytes without refusing), so it was
tried in Brave without a window (a Chromium, driven over its debugging protocol, with a profile that was thrown away):
the storage was filled until the browser refused (`QuotaExceededError`, after 4 megabytes and 933 kilobytes) and 30
kilobytes were given back; the antibody example (320,596 letters as it is kept) was then shown in place of ubiquitin: it
was left out, the list was written naming no result, and the page said so in its notice; room was made and another page
opened, and the structure was kept at that save, as the rule above says. **Not tried:** Firefox, Safari and a phone;
results of real data near the size of the storage; two tabs with results.

**Tried on 19 November 2026: the larger storage**, in Brave without a window (a real Chromium with its IndexedDB, driven over its debugging
protocol; the page built from the code of that day, served on this Mac). The example structure was in the database as
`session/1-structure.result` (79,580 letters) and not in `localStorage`, which held the list (121 bytes) and the mark of the tab;
after *Reopen* it was written again under the next generation. With `localStorage` filled until the browser refused a write
(`QuotaExceededError`, after 4 megabytes and 1,010 kilobytes) and 30 kilobytes given back, the antibody (320,596 letters) was kept all
the same, in the database: the page could not do that before. A result of 6,024,435 letters, put in the database by a script with a list
that gives no checksum (as an earlier version writes), was offered, read and shown (1,500 sequences of 4,000 letters), and after a change
written again (`102-alignment.result`, with the list, 127 bytes, in `localStorage`). With that 6 MB in the database the page was
started in 154 to 165 ms, against 60 to 131 ms with none (three starts each; the very first start, which also compiles the program,
took 551 ms). With a database whose `open` never calls back the page was still on *Loading HelixLoom…* at 1.5 s and had started at
6.1 s, and kept its results in `localStorage`. **Not tried:** Firefox and Safari; a database that really refuses a write (the
browser's never did: that is tested with a stand-in only); two tabs; a tab closed in the middle of a write; results over 6 megabytes.
The deliberate-breakage round was run for this step on 19 November 2026, after the page was tried: 72 breakages of the code were each
put in by a script and the tests run, and **all 72 were caught** (61 on the first run; the 11 that were not got new tests, and then were).
The mended version was not read again by independent readers.

**Tried on 18 November 2026: the runs of the Pipelines page** (the page built from the code of that day). In the built-in
browser: *Try the example* on Pipelines (align, tree and primers of six isolates); the storage then held
`1-pipelines.result` (5,833 letters) and the list that names it, and nothing else. The page was loaded again and asked
*Reopen your last session? No sequence was open last time. Kept from it: the Pipelines page.* After *Reopen* the page
showed what ran, the alignment, the tree and the primer pairs as before, the notice said *Put back the Pipelines page from
the last session*, and *How this was made* said *Nothing has changed since: running it again gives this result*. The
desktop window was drawn (a picture, no window) from a throwaway settings folder after a run of a pipeline file had been
drawn into it: the folder held `1-pipelines.result` and the question was the same one; *Reopen* was not clicked there.
**Not tried:** a pipeline file's run in a browser (a page has no folders to choose), a run from the workstation kept and put
back, and the full storage with a run in it (the rule is the same as for the other results, and is tested with the stand-in).

**How far the tests reach, for the runs of the Pipelines page.** The keeping of them was broken on purpose in 55 places,
one at a time. The tests caught 52 at once. The three left each take away one of two counts that follow each other (a run's
start and its end both count that it has changed): either alone changes nothing the keeper can see, and with both of a
pair taken away the tests caught it, for the route and for a pipeline file.

**How far the tests reach, for the results.** The keeping of results was broken on purpose in 79 places, one at a time. The
tests caught 77 at once and, with a test written for one that got through, 78. The one left takes the count away from the
emptying of the Alignment page's list, which changes nothing that can be seen: the list also counts every sequence that is
added to it.

**How far the tests reach.** The code that keeps things in a browser and the server's packed copies were broken on purpose
in 74 places, one at a time. The tests caught 69 at once and, with four tests written for what got through, 73. The one left
is the line that hands a page the browser's storage (`Store::of_the_page`): it exists only in the browser build, so no test
that runs on a computer reaches it, and the trial above is what shows that it works.

## Installed, and without a network

Since 17 November 2026 the page that `scripts/build-web.sh` makes can be **installed** like a program (the browser offers
it: *Install HelixLoom* in the address bar or the menu of Chrome and Edge, *Add to Dock* in Safari), and **opens without a
network** once it has been opened with one. The program runs in the browser, so everything goes on without a network but
the work that is sent to a workstation, which says that it cannot be reached, as it does when the server is down.

What makes it so is three things the script puts beside the page (`web/` has their sources):

* **a manifest** (`manifest.webmanifest`): the name, the colours, that it opens in a window of its own, and
* **two icons** (`icon-192.png`, `icon-512.png`), which are the program's own mark, drawn by the program
  (`helixloom icon web/icon-192.png 192`; run it again when the mark changes);
* **a service worker** (`sw.js`), which keeps the files of the build it came with. It is written to answer as little as
  it can: the page itself (the newest when the server answers, so that a new build is seen at once; the kept one when it
  does not) and the files of its own build, by their exact names (they are named after what is in them, so a kept one is
  never an old one). **Every other request is left to the browser as if there were no worker**, above all the calls to a
  workstation, which carry a token and bring back results: they are never looked at and never kept. A build is kept
  whole or not at all, in a cache named after the page's folder and the build, and the caches of that page's older builds
  go when a new worker takes over (when the last tab of the old one is closed; until then a tab without a network shows
  the build it was opened with). Another HelixLoom page in another folder of the same site has caches of its own, which
  this one does not touch. The files are asked of the server itself when they are kept, not of what the browser holds by
  its own rules, and a page that the server sends by way of a redirect is kept as an answer of its own, since a browser
  refuses to open a page from a redirected one.

**It needs what browsers ask of a page that wants this: `https://`, or the computer itself** (`http://localhost`,
`http://127.0.0.1`). The page of a workstation that is reached as `http://workstation.local:7878` works as before but is
not installed and not kept: put a proxy with a certificate in front (`DEPLOY.md`, the two things the browser insists on) or
open it through an SSH tunnel as `http://127.0.0.1:7878`. A page built with trunk has no worker (only the script writes
one, because the worker has to name the files of the build).

**The script is careful with the folder it is given.** The page is put together in a folder of its own and goes into the
named one only when all of it is there, so a build that fails leaves the page that was there whole. It takes away only
what it wrote before: the files it names after a checksum, by that very shape of name (a `helixloom-notes.js` of your own
stays), and the five that keep their names. It refuses the repository's root and its `web` folder, where its own sources
have those names.

**How it is tested.** `tests/web.rs` runs the script on two small files that stand for the program
(`HELIXLOOM_WEB_PREBUILT`; the script says on its error stream when it builds nothing) and checks what it puts in the
folder: the page names its files, the worker lists exactly the page and what it needs, every packed copy unpacks to its
file, another build has other names and another version, so has a build in which only the manifest, an icon or the worker
is another, files that are not the script's stay, a build that cannot be finished changes nothing, and the two folders of
sources are refused (tried on a copy of the script and its sources, not on the repository). It then runs the worker itself
in node with a cache and a network that can be switched off (`tests/web/sw_harness.mjs`): twelve things it must hold, among
them that fifteen kinds of request it has no business with are left alone and the network is asked for none of them, and
that the page's own short script asks for the worker only where a browser has workers. On a computer without node that
test says so and tests nothing; where the tests are run for everyone (`CI` is set) it fails instead. The worker and the
script were then broken on purpose in 38 places, one at a time, and the tests caught all 38. **What no test here can
show** is a real browser: the stand-in keeps plain objects where a browser keeps answers that can be read once, and knows
nothing of a browser's rules for requests. That is what the trial below is for.

**Tried in a real browser on 17 November 2026** (Brave 157, which is Chromium, run without a window and with a profile
made for the trial; the page served by `helixloom serve` on this Mac): the page loaded, the worker became active and its
cache, named after the page's address and the build, held exactly the seven files; a call to the workstation without a
token was answered by the server (401) and nothing of it was in the cache. Then the server was stopped, so that nothing
answered at its address, and the page was asked for again with another view: it opened, drawn from what was kept (the
Trees page with its example), under the worker's control, and the call to the workstation had no answer. The headers the
server sends with the worker, the manifest and an icon were read too (`text/javascript`, `application/manifest+json`,
`image/png`, each to be asked about every time). **Not tried:** installing the page from a browser's menu, and the window
of its own that it then opens in; Safari and Firefox; a new build taking over from an old one in an open tab; a server
that redirects; the page under a folder of a site (the tests cover that with the stand-in only). The browser that is built
into the tool this was written with could not be used: it fails to fetch any worker's script, from this server and from a
plain one alike, and the page then says in the console that it is not kept and goes on as before, which is what it should
do in a browser that refuses.

## Reaching a workstation from the browser

A page cannot open the plain TCP connection that the program uses, but the workstation server speaks HTTP, so the page sends
the same requests through the browser's `fetch`, with the same bearer token (`src/remote/fetch.rs`; the client is written
once, as `async` functions, and the program waits for them on a thread). The jobs are started with `Job::spawn_async`, so a
job on the workstation goes on between frames and the page does not stall.

**The simplest way** is to let the server hand the page out, so that the page and the server have one address:

```sh
helixloom serve --listen 127.0.0.1:7878 --token-file token.txt --web-dir dist
```

Open `http://127.0.0.1:7878/` (through an SSH tunnel if the server is on another computer), go to Workstation (the address is
already filled in), type the token and press Connect. The server reads `dist/` once, when it starts, into memory, and a request
can only look names up in that list; the page is sent to anyone, without a token, and every address under `/v1/` still needs
one. Keep the build in a folder of its own: the server refuses to start if the token, its file or its jobs' folder is in the
one it hands out. **A page hosted somewhere else** can call the server if the server names it: `--allow-origin https://lab.example.org`
(exactly the origin the browser shows; there is no "any origin"). Two limits of the browser: a page opened over `https://`
cannot call an `http://` address (open it from the server, or put a proxy with a certificate in front and type its `https://`
address), and the server does not encrypt anything, so use an SSH tunnel or a private network (`DEPLOY.md`, 2 and 2.1).

**The connection, kept for the tab.** Until 17 November 2026 a page kept nothing of the connection, and the token had to
be typed again every time the page was loaded. The owner of the program decided that a page may keep the token until its
tab is closed, and no longer. So the page keeps the address, the token, where jobs run and how many GPUs a job gets in
the table the browser gives **one tab** (`sessionStorage`, under `helixloom/workstation.txt`), and never in the site's
table, which stays (`ui/session/browser.rs` has both tables; `Keep` in `ui/platform/workstation.rs` is where the
connection goes). Loading the page again in the tab connects as before; *Forget* removes the entry at once; another tab
starts without it. **What that table is and is not**, which the page's own words under *Keep it private* say too: a tab
that is duplicated has a copy of its own, and a browser that brings a closed tab back brings the entry back, so on a
computer that others use the sure way is *Forget*, not closing the tab. Any page of the same site that is loaded in that
tab can read and write the table: **host HelixLoom at an address of its own** (the workstation's own server does that),
not in a folder of a site that has other pages. What the page does about it: it keeps only what made a connection, it
removes an entry that makes none and shows nothing of it, and when the browser will not take a new connection it removes
the one before, so that the page never comes back connected to where it was instead of where the person last connected.
The token goes nowhere else: not into the site's table, an address, a report, the diagnostics or a message.
**Tried in a real browser on 17 November 2026** (the page served by `helixloom serve` on this Mac, with a token made up
for the trial and removed after it): the token typed and *Connect* pressed, and the browser's two tables read: the tab's
held the address, the token, where jobs run and the GPUs, the site's held nothing of them, and neither did the address
of the page; *Workstation* chosen under *Where jobs run* and the page loaded again: connected without asking, the choice
as it was left; a second tab of the page: its table empty, not connected; *Forget*: the tab's table empty at once. This
trial is the only thing that covers the lines that reach the browser's tables, which no test on a computer can run: a
mistake there (the site's table in the place of the tab's) would pass every test. **Not tried:** Firefox and Safari; a
duplicated tab and a tab brought back by the browser, which the page's words warn of; a browser that gives the page no
storage.

**Tried** (a real browser, the real server on the same Mac, no GPU): the page served by the server; connect; status and the job
list; an alignment and a Bayesian tree sent to the workstation and shown in their pages; a whole pipeline (the six-isolate
example: align with the workstation's MAFFT, a tree with its IQ-TREE, primer pairs) run there and shown; a malformed job refused by the server
and its message shown; the wrong-token message; the server-down message; a page of another origin that the server allows (the
preflight and the calls) and one that it does not (refused by the browser); removing a job (a `DELETE` from another origin).
**Not tried**: Firefox and Safari; an `https://` address; structure, PhyloCNN and species-tree jobs sent from the page (they use
the same call as the alignment, the tree and the pipeline, but were not run); a server with the page cap (16 at once) or the
slow-reader deadline met by a real browser (both are tested with the server's own code only); results of many megabytes; many jobs at once; a
workstation with real GPUs.

## Next steps for the web platform

1. **Browser storage**: the open sequences, what the Alignment, Trees, Structure and Pipelines pages show, and the
   preferences are kept, and the remembering can be turned off (above); the results are in the browser's database
   (IndexedDB), since 19 November 2026. Still to do: what the tool windows show, projects with names, and links to a result
   that others can read. The connection to the workstation is kept for the tab (above).
2. **Web workers** for the jobs that run on this computer: `Job::spawn` and `par_map` are the two places to replace (a worker
   pool with `wasm-bindgen-rayon` needs cross-origin isolation); a long local job then no longer stalls the page either.
3. **Installable app**: built (*Installed, and without a network*, above). Still to do there: a word in the page when a
   newer build is waiting to take over. The results that are kept (1) are in the browser's database, so they are there
   without a network too; large ones need the larger storage of 1.
4. Put the workstation server on a cargo feature so the web build does not compile it (the command line is already left out of
   the web build; `#![cfg_attr(target_arch = "wasm32", allow(dead_code))]` in `main.rs` silences what only the server uses), and
   tune the release profile for size.
