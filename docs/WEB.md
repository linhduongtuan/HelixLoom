# HelixLoom in a browser

The same code runs as a desktop app and, compiled to WebAssembly (`wasm32-unknown-unknown`), in a browser tab. The desktop
build is unchanged. The page is one canvas that eframe draws on, so it needs a browser with WebGL 2 (or WebGPU).

## Build and run

You need the Rust `wasm32-unknown-unknown` target and [trunk](https://trunkrs.dev) (the Rust web bundler, `cargo install --locked trunk`;
it downloads `wasm-bindgen` and `wasm-opt` on first use). Take care that the `trunk` on your PATH is that one: the
`trunk` of trunk.io is a different program.

```sh
rustup target add wasm32-unknown-unknown
trunk build --release          # makes dist/ (index.html is the page; Trunk.toml sets the folder and relative addresses)
cd dist && python3 -m http.server 8080    # any static file server will do; then open http://localhost:8080/
```

`index.html` asks trunk for `--no-default-features`: the `gpu` feature (candle) does not build for the browser, and the code
that uses it is switched off for `wasm32` whatever the features say.

Without trunk: `cargo build --release --target wasm32-unknown-unknown --no-default-features`, then `wasm-bindgen --target web`
on `target/wasm32-unknown-unknown/release/helixloom.wasm`.

**Host the files in `dist/`**: `index.html`, one `helixloom-<hash>.js`, one `helixloom-<hash>_bg.wasm` (about 9 MB; serve it
compressed, about 2.5 MB over the wire), and one `style-<hash>.css`. The addresses are relative, so the folder can sit under
any path. The server must send `.wasm` as `application/wasm`; no special headers are needed (no threads, no
`SharedArrayBuffer`).

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
| Long jobs (likelihood, Bayesian, many PhyloCNN trees) | run inline: the page stalls while they run |
| Workstation (jobs on a GPU server) | not available: the page says so |
| MAFFT, IQ-TREE, primer3, fastp, curl, Snakemake, Nextflow | not available (the built-in methods are used); online structure prediction needs curl, so it is not available |
| Save workflow files (a folder), install PhyloCNN networks (a folder), PhyloCNN networks (model files) | not available: the page says so |
| Settings and the workstation's token kept between visits | not kept |
| GPU engine (candle) | not in the browser |

## Next steps for the web platform

1. **Workstation over fetch**: a small HTTP gateway on the workstation server that speaks the same protocol with the same bearer
   token over `fetch`, serves the web files itself (so there is no cross-origin problem), and answers CORS if hosted elsewhere.
2. **Browser storage**: keep the open sequences, preferences and (with care) the connection in `localStorage`/IndexedDB.
3. **Web workers** for jobs: `Job::spawn` and `par_map` are the two places to replace (a worker pool with `wasm-bindgen-rayon` needs
   cross-origin isolation); the page then no longer stalls.
4. **Installable app**: a manifest and a service worker, so it works offline once loaded.
5. Put the workstation server and command line behind a cargo feature so the web build does not compile them
   (`#![cfg_attr(target_arch = "wasm32", allow(dead_code))]` in `main.rs` silences what only they use), and tune the release profile for size.
