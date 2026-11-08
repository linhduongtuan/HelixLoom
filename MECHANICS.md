# Mechanics sheet — ApE (A plasmid Editor) 3.1.9 → HelixLoom

Source of truth: `/Applications/ApE.app/Contents/Resources/Scripts/AppMain.tcl` (67,911 lines of Tcl/Tk 8.6) and the
data files in `Accessory Files/`. Nothing was pasted from them: every rule was re-written in Rust, and the numeric ones
are checked by unit tests against values printed by the original's own procs (run with the system `tclsh`).

## How the original was read (deep mode)

* **There is no machine code to decompile.** `Contents/MacOS` is empty, the Tcl/Tk framework links dangle, and the bundle
  holds only Tcl text. Tried with the Craft Anything tools: the target opened as an unknown file and as a raw binary
  (base 0 and base 0xa72d3) gave no functions and no strings; the Ghidra import failed with "No load spec found";
  `native_function_pseudocode` at 0xa72d3 returned only the stub `void entry(){goto label_0;}` because those bytes are
  ASCII text, not instructions.
* **What did work.** `native_bytes_find` located the text of `proc Tm` at 0xa72d3 and `proc extinction` at 0xa7bc5, and
  `native_bytes_show` returned the readable text at those addresses. So **"address" below = byte offset of the `proc`
  keyword in `AppMain.tcl`** (hex; the same number is the file offset). The other offsets were computed the same way
  (`grep -b`, byte counting) and the first two were confirmed with the tool. `Lnnn` = line number in the file.
* **Cross-checks.** Tm, GC %, extinction and the primer self-complementarity align were run through the original's own
  procs with `tclsh`; the printed values are the test vectors in `dna.rs` and `primers.rs`.
* **Run again from scratch in a later session, with the tools, in this order** (evidence ids from the session log):
  `core_target_open` on the bundle (a folder, `ev_fa764408…`) → the real launcher is `Wish` per `Info.plist` but
  `Contents/MacOS` is empty, and `file` finds no Mach-O/ELF anywhere in the bundle → the main program is therefore
  `AppMain.tcl` (3,245,692 bytes, `core_target_triage` `ev_aa8ca6fd…`: "not a format we recognise") → reopened as raw
  arm64 with the Ghidra language pinned: `native_strings_find` (two names, `ev_02b308b7…`, `ev_53f39fae…`) found 0
  strings (it only searches data sections), `native_functions_list` found 1 "function" (`entry`, 0 instructions,
  `ev_69538851…`), and `native_function_pseudocode` answered **function_not_found** ("Ghidra could not disassemble code
  at this address") for `Tm` 0xa72d3, `search_with_mm3` 0x2a43a7, `pcr_dialog_do_pcr` 0x2944ee, `pcr_dialog_scan`
  0x297522 and `pcr_dialog` 0x28656f. The bytes there are ASCII, not instructions. `native_bytes_find`
  (`ev_e80801d5…`: 29 `proc pcr…` hits; `ev_b4abd2d9…`: no `proc fasta`) and `native_bytes_show` (`ev_64921f7c…`: the
  text of `search_with_mm3` at 0x2a43a7) are the tools that do read this program; the text was then read in full.
* **Run once more, in the order asked for, in the last session** (the person asked for it again): `core_target_open`
  on `AppMain.tcl` → "a file we don't recognise yet" (`t_4467f482ee05`); `native_strings_find` on that target →
  `not_supported: the file is malformed (Unknown file magic)`; reopened as raw arm64 with the Ghidra language
  `AARCH64:LE:64:v8A`: `native_strings_find` for `translation_window` → 0 strings (`ev_9a866f0a…`),
  `native_functions_list` → 1 function, `entry`, 0 instructions (`ev_01903f84…`), and `native_function_pseudocode`
  answered **function_not_found** one by one for `translation_window` 0x1a8803, its molecular-weight line 0x1aaa55,
  `translate` 0xa51b0 and `clip_copy_special` 0x90351. `native_bytes_find` found each procedure's text: `translation_window`
  0x1a8803 (`ev_06829752…`), the molecular-weight line 0x1aaa57 (`ev_95dc55ca…`), `clip_copy_special` 0x90351
  (`ev_d68edd69…`), `bankit_text` 0x8e450 (`ev_2bb7c9d2…`), `feature_coords_text` 0x1d05cf (`ev_e5d471b6…`), `com`
  0xa4b5a (`ev_e91c6397…`) and the Tcl 8.5 `revcom` 0xa49b0 (`ev_a8304599…`). Those texts were then read and **run**:
  the original's own lines for `revcom`, `com`, `rev`, `gen_translation_map`, `translate`, the whole line-building part
  of `translation_window` (lines 36549–36555, 36567–36598, 36605–36681) and the conversions of `clip_copy_special`
  (lines 12368–12391) were cut out of the file unchanged and run with the system Tcl 8.5 against a stand-in for the
  text widget and the clipboard. What they printed for 21 Translate-window cases and 45 Copy Special cases is the
  expected output of the unit tests in `src/seq/translation/golden.rs` and `src/seq/copyspecial/golden.rs`.
* **And once more** (a sixth request), in the same order: `core_target_open` → unknown file; `native_strings_find` →
  "Unknown file magic"; reopened as raw arm64 (Ghidra `AARCH64:LE:64:v8A`): `native_strings_find` for `find_orf_min` →
  0 strings (`ev_12f96c95…`), `native_functions_list` → the single `entry` (`ev_01903f84…`), and
  `native_function_pseudocode` → **function_not_found** for `orf_list` 0xaa4a3 and `next_ORF` 0xa9b56.
  `native_bytes_find` located `next_ORF` 0xa9b56 (`ev_b63de354…`), `orf_list` 0xaa4a3 (`ev_206881a7…`), `do_find_mm`
  0x9fdda (`ev_8dfe4e13…`), `jump_to` 0x976e3 (`ev_263f1f8a…`), `select_from_to` 0x96c8f (`ev_7db28337…`) and
  `check_char` 0x99474 (`ev_a299bbcf…`). Their text was read and the original's own lines were run with Tcl on random
  inputs: 240 Find-Next-ORF cases, 300 Find cases (`do_find_mm` itself, whole procedure) and 28 Jump / Select
  From-To cases are the expected output of the tests in `src/seq/nextorf/`, `src/seq/findmm/` and `src/seq/goto.rs`.
  Of the 240 ORF cases the earlier build, which listed ORFs its own way, chose a different one in **57**.
* **A seventh request**, same order: `core_target_open` → unknown file; `native_strings_find` → "Unknown file magic"; raw
  arm64: 0 strings for `annotate` (`ev_58b0c664…`), the single `entry` (`ev_01903f84…`), **function_not_found** for
  `pcr_dialog_do_pcr`, `translation_window` and `pcr_dialog_filter`; `native_bytes_find` found the text of the seven PCR
  procedures (ids in section 8, "PCR"), which were read and run: 280 scans, 107 primer pairs and 240 feature-copy cases.
  Of 280 primers the earlier build found different annealing sites for **50**, and of 84 primer pairs it predicted
  different products for **13**.
* **An eighth request** (again in the same order, this time for what the Signals workspace stands on). `core_target_open`
  on the bundle → a folder (`t_f903046d7ecc`; `core_target_triage` `ev_fa764408…`); `core_archive_show` → 7,976 entries,
  7,898 files, `Contents/MacOS` empty and no `Tcl` or `Tk` library inside the two framework folders either
  (`ev_2f05c8f0…`): **the bundle holds no machine code at all**. `core_target_open` on `AppMain.tcl` → unknown file
  (`t_4467f482ee05`); `native_strings_find` → "Unknown file magic"; reopened as raw arm64 (`AARCH64:LE:64:v8A`):
  `native_strings_find` for `abif` → 0 strings (`ev_b044664c…`), `native_functions_list` → the single `entry`
  (`ev_01903f84…`), `native_function_pseudocode` → **function_not_found** for `read_abi_info`. `native_bytes_find`
  located the text of eleven procedures (table below). Searching `^proc abi_` found **the original's whole Sanger trace
  window, which section 7 had wrongly listed as not built**: it was read in full, `read_abi_info` and `abi_save_scf` were
  cut out unchanged and run with the system `tclsh` 8.5.9, and the build now follows them where it had followed the file
  format's description or its own ideas (section 8, "Sanger trace window", lists every difference).
* **A ninth request** (the person said nothing had been decompiled and named the same tools in the same order; this time
  for what the sequence view, the trace window and the Workstation stand on). `core_target_open` on the bundle → a
  folder (`t_f903046d7ecc`); on `AppMain.tcl` → "a file we don't recognise yet" (`t_4467f482ee05`); `native_strings_find`
  → `not_supported: … Unknown file magic`; reopened as raw arm64 (`AARCH64:LE:64:v8A`): `native_strings_find` for `proc `
  → 0 strings (`ev_7e80127d…`), `native_functions_list` → the single `entry` (`ev_01903f84…`, the id it always gives),
  `native_function_pseudocode` → **function_not_found** for `read_abi_info` 0x1684b9. The bytes are Tcl text, not
  instructions, so what reads this program is `native_bytes_find` and the file itself. `native_bytes_find` located
  `proc textarea_width_changed` twice, 0x2cad6e (L61796) and 0x2d0413 (L62344; the second definition replaces the first
  when the file loads) (`ev_cfb4875f…`), the `<Configure>` binding that calls it, 0x1d0b4 (L2067, `ev_3a1eed4d…`),
  `set info(default_textarea_width) 75`, 0x3031da (L66907, `ev_710c3e1f…`), and `set info(abi_default_canvas_height) 400`,
  0x310f12 (L67219, `ev_dee6598a…`). The text was read in full and **changed the build in two places**:
  1. **The row width of the sequence view.** The earlier build stopped at 75 bases a row (the sheet said "up to 75 bases,
     fewer when the window is narrower", marked estimated). The original uses 75 only until its window has a size; from
     then on it sets the width again at every resize to `int((window width − 4) / character width)` (`textarea_width_changed`,
     bound to `<Configure>`), so a wider window shows more bases in a row, with no upper limit. Its ruler is written in
     tens (`    *` and the number, right aligned) while the string is shorter than the row, so the last ten can be cut
     short: at 75 columns the star of the eighth ten shows and its number does not. HelixLoom now does both. Its one
     added limit is a lower one, 10 bases a row, where the original has none (estimated).
  2. **The trace window's vertical layout.** `new_abi_window` moves the four traces to y = 35 and scales them to the
     height of the canvas less 45 (L31110–31111), which leaves 10 points under them; the base letters (at one line
     height) and the base numbers (at 0) sit above. The earlier build started the traces 42 points down, a figure of its
     own; it now starts them at 35 and leaves 10. The first window opens 650 × 400 points (`abi_default_canvas_width`,
     `abi_default_canvas_height`), which HelixLoom does not copy: its Signals page fills the space it is given.
  **For the Workstation the original has nothing to read.** The text holds no `nvidia`, `cuda`, `gpu`, `workstation`,
  `bearer`, `ssh` or `socket -server` anywhere (searched, 0 hits). Its only network code is a handful of `::http::geturl`
  calls (the web-update check with `-timeout 1000`, L9588; the REBASE enzyme download; a LabArchives upload with
  `-timeout 4000`, L33906; BLAST submission with `-timeout 20000` and `-blocksize 4096`, L49603, and its status poll with
  `-timeout 4000`, L49685). Those are small requests to public sites and do not carry over to jobs that take minutes to
  hours, so the Workstation keeps its own deadlines, which stay **estimated** (section 9).
* **A tenth request** (the same words, again): `core_target_open` on the bundle (`t_f903046d7ecc`) and on `AppMain.tcl`
  (`t_4467f482ee05`) → unknown file; `native_strings_find` → "Unknown file magic"; raw arm64 (`AARCH64:LE:64:v8A`):
  `native_strings_find` for `mol_calc` → 0 strings (`ev_6923511e…`), `native_functions_list` → the single `entry`
  (`ev_01903f84…`), `native_function_pseudocode` → **function_not_found** at 0x18dde0, where `native_bytes_find` finds the
  text of `proc mol_calc` (`ev_ecc8fe64…`). Setting the original's **Tools menu** (L1293–1294, 0x1245a) beside the build's
  found two tools the build did not have: **Calculator** (⌘1), an arithmetic window (the build's "Calculator" was a window of
  its own for Tm and GC, now "Tm and GC Calculator"), and the **Molecular Reaction Calculator** (⌘2). Both are now built from
  the original's procedures (section 8, "Calculator and Molecular Reaction Calculator"), which were also run unchanged with the
  system `tclsh` 8.5.9 on about 3,900 random cases to check the Rust code against.
* **An eleventh request** (once more the same words; the session had no target open, so each was opened again):
  `core_target_open` on the bundle (`t_f903046d7ecc`) and on `AppMain.tcl` (`t_4467f482ee05`); `core_target_triage` →
  "isn't a format we recognise yet" (`ev_aa8ca6fd…`); `Contents/MacOS` is an empty folder, and the Tcl and Tk frameworks hold
  only `Resources` and `_CodeSignature` (no `Tcl` or `Tk` library file); `native_bytes_find` for the Mach-O magics `CF FA ED FE`
  and `CA FE BA BE` in `AppMain.tcl` → 0 matches (`ev_44667c45…`, `ev_5f8ce4b5…`). Opened as raw arm64
  (`AARCH64:LE:64:v8A`): `native_strings_find` for `proc ` → 0 strings (`ev_5168e123…`), `native_functions_list` → the single
  `entry` of size 0 (`ev_01903f84…`), `native_function_pseudocode` for `entry` → `goto label_0` (`ev_3fd231b9…`) and for 0x18dde0
  → **function_not_found**, where `native_bytes_find` finds the text of `proc mol_calc` (`ev_b09717d2…`). Following the eighth
  request's lesson the original's **File, Edit and Enzymes menus** (L1072–1198) were set beside the build's, which the tenth
  request had done for Tools only: the Enzymes menu has two entries that nobody had listed as missing, **Silent Sites** (L1194
  @0x101e1, `ev_06f4b51f…`) and **Find Diagnostic Site** (L1195), now built from `silent_sites` (@0x1a0b2c, `ev_76470122…`),
  `coordinate_translation` (@0x1a24be, `ev_84c22195…`), `do_diagnostic_site` (@0x1a3e5b, `ev_f091f3f4…`),
  `add_daignostic_site_dialog` (@0x1a32b3, `ev_e4a6fdca…`), `degen_transl` (@0x2ea304, `ev_87dfe624…`), `triplet_transl`
  (@0x2e9fae, `ev_f7a46f03…`) and `seq_len_mod` (@0x2aa1ab, `ev_c596debe…`) (section 8, "Silent Sites and Find Diagnostic
  Site"). Still missing from the same menus, in order of size: partial digests (`enz_digest_partial_dialog` L17587 @0xc5edc,
  `ev_948071a4…`), List Enzymes and Quick Lists (`list_enzymes_dialog` L16305, `list_enzymes_window` L16378), Text Map
  (`enz_text_map` L19401, `do_text_map` L19730 @0xdec9f, 450 lines, `ev_7bbf5f47…`), Print and Export as PDF, Insert Piece,
  LabArchives. The File menu's **Duplicate Selection** and **Duplicate Full Sequence** (L1077–1078; `duplicate_window`
  L7750 @0x5c69b, `ev_7b52804b…`), which the same comparison found, are built too (section 8, table, "Duplicate Selection").
* **A thirteenth request** (the person said again that nothing had been decompiled; the session had no target open, so each was
  opened again). `core_target_open` on the bundle (`t_f903046d7ecc`) and on `AppMain.tcl` (`t_4467f482ee05`, 3,245,692 bytes,
  sha256 4467f482…) → unknown file; `native_strings_find` → "Unknown file magic"; reopened as raw arm64 (`AARCH64:LE:64:v8A`):
  `native_strings_find` for `species` → 0 strings (`ev_f4773af1…`) and for `list_enzymes_window` → 0 (`ev_1c54c4d6…`),
  `native_functions_list` → the single `entry` (`ev_69538851…`), `native_function_pseudocode` for `entry` → `goto label_0`
  (`ev_a2941ded…`) and **function_not_found**, one call after the other and each waited for, for `list_enzymes_window`,
  `list_enzymes_text` 0xb8d47, `list_enzymes_dialog` 0xb7316, `enz_site_usage_filter` 0xb4584, `enz_count` 0xb4790 and 0x3070f7.
  The bytes there are Tcl text, so what reads this program is `native_bytes_find` (`proc list_enzymes_` → 0xb7316, 0xb82e9,
  0xb8d47, `ev_4469f099…`; `proc enz_count` 0xb4790, `ev_e7f5b09d…`; `proc enz_site_usage_filter` 0xb4584, `ev_88f514cd…`;
  `Quick Lists` → 0x1005f, 0x2b68c7, 0x2da470, `ev_c3658ef7…`) and `native_bytes_show` (`list_enzymes_window` 0xb82e9,
  `ev_ab176ea7…`; `list_enzymes_text` 0xb8d47, `ev_89f7bf17…`), and then the text itself, read in full.
  * **What the original has behind the species-tree engine built in the twelfth round: nothing.** `native_bytes_find` finds no
    `astral`, `quartet`, `coalescent`, `newick`, `phylogen`, `phylo`, `consensus`, `posterior`, `gpu` or `workstation` (ten
    searches, 0 hits each: `ev_743ffbf0…`, `ev_43fa4620…`, `ev_095cae83…`, `ev_ccba0451…`, `ev_2723faee…`, `ev_1b5b97a6…`,
    `ev_e5b9bddc…`, `ev_9a3f3809…`, `ev_94471653…`, `ev_924c2752…`), and `species` only four times, all in GenBank qualifier
    lists (0x3070f7, 0x308df2, 0x309811, 0x30b2e2; `ev_9b67ccbd…`, `ev_48f9d69e…`, `ev_397f6686…`). `tree` occurs 735 times
    (`ev_bd20b923…`), every one the Tk `ttk::treeview` that shows the feature table (`treeview_*`, `features_to_tree_view`).
    So the ASTRAL engine has no values, formulas or addresses to follow: it is an addition (section 9) and nothing in it changes.
  * **What did change the build.** The Enzymes menu of the original (L1160–1200) set beside the build's again: **List Enzymes…**
    (L1184 @0xfe05) and **Quick Lists ▸ Site Usage, Unique Sites, Absent Sites** (L1186–1189 @0xfe92, 0xff1a, 0xffa7) were
    missing. They are now built from `list_enzymes_dialog`, `list_enzymes_window`, `list_enzymes_text`, `enz_site_usage_filter`
    and `enz_count` (section 8, "List Enzymes and Quick Lists"). Running the original's own `findenzymes` behind them, with
    its own `initialize_methylation`, `not_in_list` and `toggle_methylation` (the earlier checks had stubs for those), showed
    **a rule the earlier port of the site finder missed**: the 21 enzymes with more than six bases are found by a
    regular-expression search that never reports a site overlapping the one before, and `find_sites` now does the same
    (section 2). Of the 72 random cases, case 9 was the first to differ (MauBI in `CGCGCGCGCG`: one site, not two). A second
    difference came out only with a stress file of 400 cases made of sites repeated at their own period (so that they
    overlap): the original keeps one reverse Type IIS site at the very left end of a linear sequence that cuts a base outside
    it, and `find_sites` now keeps it too, as a site with no cut (section 2). With both rules built, 2,400 more random cases
    (500 plain and 1,900 with such runs; the driver's fourth argument) give the same lists and site lists as the original.
  * **A correction to the eleventh request's list of what is missing.** `enz_digest_partial_dialog` (L17587 @0xc5edc) cannot be
    reached: its menu entry is commented out in the original (L1180 and L59624, `#$w.menubar.enzymes add command -label
    "Partial Digestion..."`), as is the "Do Partial Digests" action of the Enzyme Selector (L17730). It is not built, for that
    reason. Still missing from that menu: Text Map (L1197), the Selection Only switch (L1170 @0xf893), Print and Export as PDF,
    Insert Piece, LabArchives.

### Procedure addresses

| Procedure | Line | Address | | Procedure | Line | Address |
|---|---|---|---|---|---|---|
| `create_window` | L950 | 0xc840 | | `enz_cut_offset` | L18867 | 0xd3d81 |
| `selection_manager` | L11062 | 0x83600 | | `enz_graphic_map` | L20210 | 0xe4e64 |
| `genbank_text` | L11900 | 0x8be2c | | `feature_scan` | L38953 | 0x1c4cce |
| `atgcfilter` | L9708 | 0x76bfa | | `feature_from_found` | L40244 | 0x1d2ff5 |
| `clip_cut_copy` | L12187 | 0x8e9b0 | | `read_feature_library` | L41749 | 0x1e2f8f |
| `clip_paste` | L12274 | 0x8f756 | | `text_to_regexp` | L41808 | 0x1e3ca9 |
| `paste_file` | L12465 | 0x912c2 | | `findenzymes` | L42931 | 0x1f09a2 |
| `get_tags_exact` | L12530 | 0x91e69 | | `initialize_methylation` | L43126 | 0x1f235d |
| `put_tags` | L12847 | 0x947aa | | `primer_scan` | L45305 | 0x210b1d |
| `do_find` | L13739 | 0x9d758 | | `align` | L45551 | 0x2132e9 |
| `insert_enzyme_fill_cascade` | L13374 | 0x99750 | | `NW_affine_align` | L48137 | 0x231136 |
| `insert_enzyme_site` | L13406 | 0x99b3c | | `gibson_find_ends` | L52632 | 0x2629bd |
| `revcom` (later definition wins) | L14451 | 0xa49b0 | | `gibson_find_chains` | L52675 | 0x263188 |
| `gen_translation_map` | L14477 | 0xa4c86 | | `textarea_insert` | L60821 | 0x2c24c5 |
| `translate` | L14518 | 0xa51b0 | | `split_tags_after_insertion` | L60837 | 0x2c2763 |
| `Tm` | L14730 | 0xa72d3 | | `textarea_delete` | L60871 | 0x2c2c1e |
| `GC%` | L14782 | 0xa7a6b | | `read_enzymes` | L64461 | 0x2ea4a6 |
| `extinction` | L14792 | 0xa7bc5 | | `info(graphic_map_properties)` | L67476 | 0x314248 |
| `orf_list` | L15100 | 0xaa4a3 | | `enz_digest` | L17086 | 0xbfeb4 |
| `change_case` | L15162 | 0xaad53 | | `swap` | L15259 | 0xab8c6 |
| `set_origin` | L15267 | 0xaba02 | | `toggle_circ` | L15322 | 0xac08b |
| `linearize` | L15351 | 0xac494 | | `rev_com_window` | L15359 | 0xac54a |

More procedures read (PCR, saving, fonts, undo):

| Procedure | Line | Address | | Procedure | Line | Address |
|---|---|---|---|---|---|---|
| `save_file` | L11737 | 0x8a413 | | `pcr_dialog_scan` | L56729 | 0x297522 |
| `clip_copy_special` | L12354 | 0x90351 | | `pcr_dialog_precompute` | L56880 | 0x29863b |
| `change_app_font_size` | L9280 | 0x72478 | | `pcr_dialog_filter` | L56952 | 0x298e61 |
| `register_undo_separator` | L61140 | 0x2c4f60 | | `pcr_dialog_set_dist` | L56996 | 0x299218 |
| `pcr_dialog` | L55236 | 0x28656f | | `search_with_mm3` | L57959 | 0x2a43a7 |
| `pcr_dialog_do_pcr` | L56473 | 0x2944ee | | | | |

Read in the last pass (Translate window, Copy Special):

| Procedure | Line | Address | | Procedure | Line | Address |
|---|---|---|---|---|---|---|
| `translation_window` | L36460 | 0x1a8803 | | `bankit_text` | L12136 | 0x8e450 |
| …its molecular-weight line | L36616 | 0x1aaa57 | | `feature_coords_text` | L39970 | 0x1d05cf |
| `clip_copy_special` | L12354 | 0x90351 | | `com` | L14468 | 0xa4b5a |

Read after that (ORF stepping, Find, Jump):

| Procedure | Line | Address | | Procedure | Line | Address |
|---|---|---|---|---|---|---|
| `next_ORF` | L15025 | 0xa9b56 | | `do_find_mm` | L13981 | 0x9fdda |
| `orf_list` | L15100 | 0xaa4a3 | | `find_dialog` | L13552 | 0x9b20f |
| `select_from_to` | L13097 | 0x96c8f | | `jump_to` | L13169 | 0x976e3 |
| `check_char` | L13345 | 0x99474 | | | | |

Read in the PCR pass (scan, product, copied features):

| Procedure | Line | Address | | Procedure | Line | Address |
|---|---|---|---|---|---|---|
| `pcr_dialog_do_pcr` | L56473 | 0x2944ee | | `pcr_dialog_precompute` | L56880 | 0x29863b |
| `pcr_dialog_scan` | L56729 | 0x297522 | | `search_with_mm3` | L57959 | 0x2a43a7 |
| `get_tags2` | L12623 | 0x92a8f | | `get_tags_get_coords` | L12755 | 0x93d6c |
| `bp2ix` | L61660 | | | `pcr_dialog_filter` | L56952 | 0x298e61 |

Read in the trace pass (Sanger trace window):

| Procedure | Line | Address | | Procedure | Line | Address |
|---|---|---|---|---|---|---|
| `read_abi_info` | L30572 | 0x1684b9 | | `abi_scale` | L31809 | 0x175ed0 |
| `new_abi_window` | L30863 | 0x16b63a | | `abi_x_scale` | L31823 | 0x1760a8 |
| `abi_draw_phred` | L31524 | 0x173324 | | `abi_file_info` | L31840 | 0x176398 |
| `abi_get_gap` | L31664 | 0x174ddd | | `abi_save_scf` | L31933 | 0x177560 |
| `abi_get_seq` | L31691 | 0x175151 | | `abi_revcom` | L31733 | 0x1754b2 |
| `mix_color` | L27844 | 0x14aaa3 | | | | |

Read in the enzyme-list pass (List Enzymes, Quick Lists, and what they stand on):

| Procedure | Line | Address | | Procedure | Line | Address |
|---|---|---|---|---|---|---|
| `list_enzymes_dialog` | L16305 | 0xb7316 | | `enz_count` | L15950 | 0xb4790 |
| `list_enzymes_window` | L16378 | 0xb82e9 | | `not_in_list` | L43171 | 0x1f2d98 |
| `list_enzymes_text` | L16463 | 0xb8d47 | | `toggle_methylation` | L43206 | 0x1f3111 |
| `enz_site_usage_filter` | L15931 | 0xb4584 | | `ix2bp` (later definition wins) | L62219 | 0x2cf36d |
| menu: List Enzymes… | L1184 | 0xfe05 | | menu: Quick Lists | L1186–1189 | 0xfe92 |
| menu: Selection Only | L1170 | 0xf893 | | start values `enz_list_*` | L66960–66962 | 0x303969 |

Values that sit in the settings block rather than in a procedure are listed by line with their byte offsets in the
appendix at the end; `mechanics.json` carries the same sources.

## 1. Sequence maths (`src/seq/dna.rs`)

| Rule | Value / formula | Address |
|---|---|---|
| Reverse complement | A↔T C↔G B↔V D↔H K↔M R↔Y, S W N `*` self; case kept | `revcom` @0xa49b0 |
| Translation | Codons T,C,A,G order; any non-ACGT base in a codon → `X`; trailing partial codon dropped; U→T | `translate` @0xa51b0, `gen_translation_map` @0xa4c86 |
| Melting temperature | SantaLucia 1998 nearest neighbour. ΔH×10 / ΔS×10 per dinucleotide: AA/TT 79/222, AT 72/204, TA 72/213, CA/TG 85/227, GT/AC 84/224, CT/AG 78/210, GA/TC 82/222, CG 106/272, GC 98/244, GG/CC 80/199 | `Tm` @0xa72d3 |
| Tm ends | AT end: ΔH −23, ΔS −41; GC end: ΔH −1, ΔS +28 (each end) | `Tm` @0xa72d3 (L14745–14749) |
| Tm result | `ΔH·100 / (ΔS·0.1 − 1.987·ln(DNA nM/1e9)) − 273.15 + 12·log10(salt mM/1000)`, clamped ≥ 0, rounded to 0.1; `---` above 100 bases; non-ACGT counted as T | `Tm` @0xa72d3 (L14771) |
| Tm defaults | DNA 250 nM, salt 50 mM | `Tm` @0xa72d3 signature; `ApE_Defaults.txt` primer_dna / primer_salt |
| GC % | every letter except A T N B D H K M R V W Y `*` counts as G/C | `GC%` @0xa7a6b |
| Extinction coefficient | single bases A 15.4 C 7.4 G 11.5 T 8.7; pairs AA 13.7 AC 10.6 AG 12.5 AT 11.4 CA 10.6 CC 7.3 CG 9.0 CT 7.6 GA 12.6 GC 8.8 GG 10.8 GT 10.0 TA 11.7 TC 8.1 TG 9.5 TT 8.4; `e = single(first) + Σ over neighbours (2·pair − single(left))` in L/mmol/cm; an empty sequence or a base outside ACGT is an error. Test vectors from the original: ATGGCGTTTGAACGGCTCAA 194.3, GATCAAGGAA 111.3, AC 21.2. **The original defines this but never calls it**; HelixLoom shows it in its own Tm and GC Calculator (Tools menu; the original's Tools → Calculator is another thing, an arithmetic window, section 8) | `extinction` @0xa7bc5 |
| ORF search | starts = every `ATG` (or the codon after any stop in "Stop+1" mode); first in-frame TAA/TAG/TGA ends it; kept when `length − 1 > min` | `orf_list` @0xaa4a3 |
| ORF defaults | min 90 bp, start rule Met, strand Top; min menu 0 30 60 90 180 300 600 900 1200 1500 | L14928, L67003, L67004, L1214 |
| Selection read-out | `Start`, `Length`, `End` each followed by `<frame>` = (n−1) mod 3 (length: n mod 3); ORF text like `M>85*/<85`; Tm and `%GC` rounded | `selection_manager` @0x83600 |
| Typed keys | only `ACGTNacgtn*` | `ApE_Defaults.txt` allowedkeys |
| Find | IUPAC patterns, both strands (find_revcom 1), case-insensitive (find_casens 0), wrap on (find_wrap 1), optional mismatches | `do_find` @0x9d758, defaults L66966 |

## 2. Restriction enzymes (`src/seq/enzymes.rs`)

* File format "tkDNAview enzymes file type 2": line 4 names, line 5 patterns, `^` = top-strand cut, `NAME(a/b)` = cut offsets past the site. 116 enzymes in the default set (all embedded). `read_enzymes` @0x2ea4a6.
* Cut offsets `[top, bottom]` from the site start: with `^` at index c → `[c, len−c]`; with `(a/b)` → `[len+a, len+b]`. Reverse-strand hits: `[start+len−bottom, start+len−top]`. `enz_cut_offset` @0xd3d81.
* Sites are found on both strands (palindromes once), circular sequences scan the first 63 bases again and keep starts < length; linear sequences drop sites whose cuts fall outside the ends. `findenzymes` @0x1f09a2.
* **Long recognition sequences do not overlap.** The 21 enzymes of the default file with more than **6** bases (`longnamelist`,
  line 2 of the file: AbsI AscI BsaBI FseI MauBI MreI NotI PacI PaqCI PmeI PshAI SapI SbfI SfiI SgfI SgrAI SgrDI SrfI SwaI XcmI
  XmnI) are found with `regexp -all -inline -indices`, once for the forward pattern and once for the reverse (`findenzymes`
  L42960–42998, @0x1f09a2), which goes on after the end of the site it has found: **a site that overlaps the one before it on
  the same strand is not reported**, the two strands' lists are then sorted together by start. The other 95 go through the hash
  scan (L43001–43044, "all sites <=6bp long"), which looks at every start. `find_sites` does the same (`HASH_SCAN_LIMIT`).
  Example: MauBI (CGCGCGCG) in `CGCGCGCGCG` has one site, at the first base; NotI (GCGGCCGC) in `GCGGCCGCGGCCGC`, which has two
  that share two bases, has one. Earlier builds found both; the difference came out when the lists of section 8 were checked
  against the original's `findenzymes`.
* **Sites that cut beyond the ends (linear sequences).** For an enzyme whose cuts reach past its site (BbsI, BsaI, BsmBI,
  PaqCI, SapI), `findenzymes` L43077–43109 takes out a forward site when its number (from 1) is above `length − reach + 1`
  (exact: the cut would fall after the end) and a reverse site when its number is below `reach − length` (**one too few: the
  reverse site whose top-strand cut falls one base before the start stays in the list**: for BsaI a GAGACC at bases 5 to 10).
  The digest then rejects it ("reject sites that CUT outside the range", L17128–17134), so such a site is counted by the
  Enzyme Selector and the lists and cuts nowhere. HelixLoom keeps it the same way (`Site::cuts_inside` is false: counted, no
  cut, no place on the map). Circular sequences drop nothing. DpnI and MalI have no site at all while Dam/Dcm methylation is off (`initialize_methylation`
  L43126 empties them, `toggle_methylation` L43206 brings them back): the build already did that, and the 72 cases of the lists
  run the original's own two procedures with the switch on and off.
* Dam/Dcm blocking (default on, `methylated_default 1`): a site is hidden when `site+d1` or `site+d2` is the start of a GATC (Dam) / CCWGG (Dcm) site. Dam: BclI 1 1, BsaBI 0 6, BspEI −2 4, BspHI −2 4, ClaI −1 3, MboI 0 0, NruI −2 4, TaqI −2 2, XbaI −2 4. Dcm: Acc65I −3 4, ApaI −3 4, BsaI −3 4, MscI −2 3, StuI −2 3. DpnI only cuts methylated DNA. `initialize_methylation` @0x1f235d, table in `Default_Enzymes.txt` line 20.
* Fragment table columns: Size, site1, site2, Mass % (`round(max(1, 100·size/total))`), sorted by size, largest first. `enz_digest` @0xbfeb4.
* **Difference (deliberate):** the original subtracts recognition-site *starts* to get sizes, so two-enzyme digests with different cut offsets can be a few bases off. HelixLoom measures between true top-strand cuts. One-enzyme digests are identical.
* Gel picture: 210 px tall, lane width 50 px, spacing 0.1 of a lane, band 3 px (+1 per repeated ladder band), y = `30 + k·log10(longest) − k·log10(size)` with `k = 210/(log10(longest) − log10(shortest))`, shortest ≥ 10; band shade = `max(10, 100·mass/maxmass)` percent ink. `enz_digest` @0xbfeb4 (gel part L17338–17395, @0xc27d4–0xc3276).
* Ladders: `1Kb plus`, `lambda Hind` (DNA_ladders.txt) and the NEB files; default is the first.
* **Insert Recognition Site** (Enzymes menu): see section 8.

## 3. Features (`src/seq/features.rs`, `src/seq/io.rs`)

* Library line = tab separated: name, sequence (IUPAC; text before `<` and after `>` is context, `#`/`+` are gaps), type, fwd colour, rev colour, graphic format, search type (0 DNA, 1–2 protein), notes. `read_feature_library` @0x1e2f8f, `text_to_regexp` @0x1e3ca9. HelixLoom matches the core sequence only and skips gap/protein entries.
* Scan: both strands, circular sequences searched twice around, duplicates removed. `feature_scan` @0x1c4cce. A forward feature uses the forward colour, a reverse one the reverse colour.
* Default colours: forward cyan, reverse green; per-type table in `feature_default_type_colors` (L67026, @0x304261), e.g. promoter #346ee0, CDS #e9d024, terminator #9d1b1c, misc_feature #7eff74.
* GenBank output: `LOCUS` line `%-16s %11s bp %-3s%-6s     %-8s %-3s %-11s`, the file comment as `COMMENT` lines (long lines wrapped at 67 like the definition; HelixLoom keeps them in `Doc::comment`), `COMMENT ApEinfo:methylated:1`, `FEATURES` with `/label`, `/ApEinfo_fwdcolor`, `/ApEinfo_revcolor`, `ORIGIN` 60 bases per line in blocks of 10. `genbank_text` @0x8be2c.
* FASTA output (`save_file` FASTA branch, L11854–11860, @0x8a413): the header is `>` plus the file's name without its
  extension, and the **whole sequence follows on one line** (no wrapping). HelixLoom does the same. "Copy as FASTA"
  (`clip_copy_special` @0x90351, L12392) puts `>title first-last` and the sequence on one line on the clipboard; that
  is built now with the rest of the Copy Special menu (section 8).

## 4. Alignment and primers

* Align Sequences: affine-gap Needleman–Wunsch, match +1, mismatch −0.1, first gap −1, further gap −0.2 (`mm_penalty`, `gap_penalty`, `gap_extension_penalty`); line width 100; block size 15. `NW_affine_align` @0x231136, defaults L67040–67045.
* Self-complementarity of a primer: align the primer to its own reverse complement (+1 / −1 / −2); report (score at the 3' end, best score). Limits 10 and 3. `align` @0x2132e9. Ported exactly and tested against 8 primers.
* Primer search defaults: length 20–25, Tm 55–60 (unrounded Tm accepted in (min−0.5, max+0.5)), GC 45–60 %, GC clamp 1–2 (consecutive G/C at the 3' end), salt 50 mM, DNA 250 nM. `primer_scan` @0x210b1d, `ApE_Defaults.txt`.
* *Estimated:* the original's primer Tm is computed incrementally with a slightly different end correction; HelixLoom uses the `Tm` proc for every candidate.
* **PCR (the original's PCR dialog, followed rule for rule and checked by running the original's own code, section 8).**
  The dialog defaults are `pcr_min_match` **17** (adjustable 5–30), `pcr_max_mm` **0** (0–10), `pcr_min_3_match` **4**
  (0–4) and "unique in template" shown first (`pcr_dialog` @0x28656f, L55412–55418).
  * **Where a primer anneals.** `pcr_dialog_scan` @0x297522 runs `search_with_mm3` @0x2a43a7 on the primer's reverse
    complement (so it is read from its **3' end**) in both directions. `pcr_dialog_precompute` @0x29863b lists the
    template in upper case and its reverse complement, and for a circle repeats the first **200** bases of each at the
    end, so a primer longer than that cannot anneal across the origin. With no mismatches allowed the first **4** bases
    from the 3' end must match exactly, whatever `pcr_min_3_match` says; with mismatches allowed it is
    `pcr_min_3_match`. A mismatch is counted, the walk ends at the first one over `pcr_max_mm`, and a template position
    past the end counts as a mismatch. The stretch is kept when the length walked (before trimming) reaches
    `pcr_min_match` and is then cut back to end on a matching base, so it can come out a base or two shorter than
    `pcr_min_match` when it ended in tolerated mismatches. A degenerate primer (anything but ACGT) or
    `pcr_min_3_match` 0 takes the slow route, whose quirks are kept: it keeps a stretch one base shorter than the
    minimum, reports the match length as an index (one short), trims only to a base the primer has as itself, and on a
    circle lists a site at the origin twice (starts run from 0 to the template length). Degenerate codes match every base
    they stand for (B→BYSKCGT, D→DRWKAGT, H→HMWYACT, K→KGT, M→MAC, N→all, R→RAG, S→SGC, V→VMSRACG, W→WAT, Y→YCT). A
    row of the table is a direction (`>>>` forward, `<<<` reverse), a location, the match length and the mismatches:
    forward, the top-strand position just after the primer's 3' end; reverse, where the annealed part starts.
  * **The product.** `pcr_dialog_do_pcr` @0x2944ee first resets each annealed length to the bases that agree from the
    3' end (case ignored) and warns "Warning: First forward primer mismatch N bases from 3' end." (or "reverse") when
    N < **15**. (Its two "primer sequence is shorter than the primer binding feature" warnings compare with the length of
    the word `fwd_dna_binding`, 15, so they can never fire; they are not built.) A linear template refuses a pair whose
    forward 3' end is not strictly left of the reverse start, so **primers that touch are refused**: "Primers don't point
    towards each other." On a circle a reverse site left of the forward one is turned round the origin, and when the
    product would be longer than the template ("round the world") the overlap is copied (an end before the start means
    more than once round) and up to four warnings follow, in this order: "Warning: Reverse primer amplifies part of
    forward primer. Primer dimers are very likely." / "Warning: Reverse primer overlaps forward primer and extends less
    than **10**bp from forward primer. Linear amplification is likely." / "Warning: Forwrd primer amplifies part of reverse
    primer. Primer dimers are very likely." (the original's spelling) / the 10bp warning again for the other end. The
    product is the forward primer's tail, the template from the forward stretch to the reverse stretch, and the
    reverse-complemented tail of the reverse primer, in the primers' and the template's own case.
  * **The new window** is named `PCR`, the template's name (**16** characters), the forward primer's name (**9**) and the
    reverse primer's (**9**), and is linear. Its comment is `PCR:`, `template: <name>`, the warnings, `<primer> <sequence>
    Tm1: x°C Tm2: y°C.` (the Tm of the exactly matching 3' part and of the whole primer, `Tm` with its own defaults;
    `Tm: x°C.` when everything matches; `---` over 100 bases; 0.0 for no bases), the same for the reverse primer and
    `Product length: N`. It carries the template's features inside the copied stretch, shifted by the forward tail
    (`get_tags2` @0x92a8f, `get_tags_get_coords` @0x93d6c: a feature is cut to the stretch; one that only touches it
    leaves an empty range) and two `primer_bind` features: the forward one over the whole forward primer in `#14c0bd`,
    the reverse one over the whole reverse primer in `#4ec02b` (`feature_default_type_colors`, L67026), each with the
    qualifier `PCR_conditions` = `primer sequence:<primer>`. The comment is saved as GenBank `COMMENT` lines before
    `ApEinfo:methylated` (`genbank_text` @0x8be2c, L11955–11960).
  * **HelixLoom:** `src/seq/pcr.rs` (`find_sites`, `make_product`, `pcr_products`), `Doc::comment`, the PCR window. Checked
    against 280 scans, 107 primer pairs (65 products and 42 refusals) and 240 feature-copy cases that were run through
    the original's own Tcl (`src/seq/pcr/golden.txt`, `coords.txt`). The earlier build found different sites for 50 of
    those 280 primers and different products for 13 of the 84 pairs it could be compared on, and had no overlap
    warnings, no comment, no primer features and no copied features.
  * **Not ported:** the primer database and tree view, the "unique in template" filter, the features of a product that
    runs through the origin (`get_tags2`'s round-the-origin branch), and `do_pcr` leaving out a primer feature the
    template already has under that name.

## 5. Assembly

* The original has **no read assembler or contig builder** (searched for contig/assemble/consensus). Its assemblers simulate cloning: Gibson (`gibson_find_ends` @0x2629bd, `gibson_find_chains` @0x263188), Golden Gate, restriction–ligation, recombination.
* Gibson rules reused as the "Gibson-style" preset: exact end overlaps of 12–100 bases (`gib_min_overlap 12`, `gib_max_overlap 100`), 0 non-homologous bases, at most 10 fragments, circular product required, fragments used in either orientation. L67101–67106.
* The read assembler is new (the request): overlap detection by shared 12-mers, ungapped verification (length and identity), both strands, containment, greedy best-overlap-first joining, optional circular closure, majority-vote consensus (ties → IUPAC code), coverage. Defaults for reads (**estimated**, no original value): min overlap 20, min identity 90 %, no maximum.

## 6. Look and feel

* Text area 75 characters wide at the start (`default_textarea_width`, L66907 @0x3031da), DNA font Courier **12** (`dnafontsize` in `ApE_Defaults.txt`). **A row is as wide as the window allows, not 75:** `textarea_width_changed` (L62344 @0x2d0413, bound to `<Configure>` at L2067 @0x1d0b4) sets the width at every resize to `int((winfo width textarea − 4) / character width)` and writes the ruler again (a ten of `    *` and the number for as long as the string is shorter than the row). HelixLoom does the same in a bundled monospace font at 12 pt (`row_columns` in `src/ui/editor.rs`; 75 only for a view that has no size yet; at least 10, which the original does not require). An earlier build capped rows at 75; that was wrong. View → Font Size offers the original's choices 8 9 10 11 12 13 14 16 18 20 24 28 32 36 40 (menu code at AppMain.tcl L1382, @0x13a0e; applied by `change_app_font_size` @0x72478).
* Undo: `register_undo_separator` @0x2c4f60 only appends to the undo list and clears redo; the original never trims it. **Deliberate difference:** HelixLoom keeps whole copies of the document, so it keeps the last 100 steps.
* Colours: selection `#6060d0` with white text (`text_select2_*`), find highlight gray75 forward / gray50 reverse, enzyme highlight red. The text area's dark background is `#1a1d21` (`text_bg_color,dark`). **Deliberate difference (platform redesign, section 9):** everything around the editor (rail, pages, cards, buttons, menus, window background) now takes its colours, rounding and spacing from one palette in `src/ui/design.rs`; the original's `#334` window background (`bg_color,dark`, AppMain.tcl:66869) is no longer used.
* Toolbar icons: two-tone line icons (black + blue); HelixLoom draws its own.
* Info strip: Mouse, Sequence, Start / insert@, Length, End, ORF, Tm, %GC, circular/linear, Dam/Dcm. Built inside `create_window` @0xc840 (L1562–1692).
* Circular map: radius 100, 0 at the top running clockwise, title at (0,−10), length at (0,+10), ticks inside, labels outside. `info(graphic_map_properties)` L67476 @0x314248.
* Default file format GenBank (`default_file_format`), extension `.gbk`.
* **Platform shell (added; deliberate difference).** The window is now a rail of eight workspaces (Home, Sequences, Signals, Alignment, Phylogeny, Trees, Structure, Pipelines; Go menu, ⌥⌘1 to ⌥⌘8) beside one page (`src/ui/platform/`). **Sequences is the original editor** with its menus, tabs, toolbar and info strip as described above; the other pages share `src/ui/design.rs`. A file opened from File → Open, dropped on the window or named on the command line goes where its name says (`src/ui/platform/route.rs`): `.ab1`/`.abi` to Signals, `.pdb`/`.ent` to Structure, `.nwk`/`.newick`/`.tre`/`.tree`/`.treefile` to Trees, anything else to Sequences, or into the sequence list when Alignment or Phylogeny is showing. Work that was started (an alignment, an estimation, a prediction, a pipeline run) goes on and finishes whichever page is showing.

## 7. Not built in this version

The trace window's editing and linking (change or add a base, link a trace to a sequence window, search a trace, print or
save as EPS; the rest of the window is built, section 8), the Molecular Reaction Calculator's lookup of a compound's weight by name and its recipe multiplier (both ask PubChem), NCBI BLAST / Entrez download, LabArchives, printing / PDF export, protein-feature search, Golden Gate / recombination / restriction–ligation wizards, partial digests (the original's menu entry is commented out, L1180, so nothing reaches its dialog), text map, RNA documents (the original converts U↔T on paste), per-feature arrow shapes (`arrow_data`, width, offset — HelixLoom draws one arrow shape), multi-range selections, the Enzymes menu's Selection Only switch (L1170; the Digest window has its own check box), Insert Piece.

## 8. Editing rules read in the deep pass (all built)

| Rule | What the original does | Address | In HelixLoom |
|---|---|---|---|
| Change case | Takes the selection (or the whole text) out, puts the upper / lower / swapped text back, then puts every tag back at the same place, so **features stay**; selection and caret are restored | `change_case` @0xaad53 | `Doc::overwrite` (same length, features untouched) |
| Swap letters | Only A B C D G H K M N R S T V W Y, both cases, are swapped; any other letter (U, E…) stays | `swap` @0xab8c6 | `dna::swap_case` |
| Reverse-complement | Nothing selected, or everything selected: the whole sequence turns. With nothing selected the caret goes to `length − caret`; with a selection it stays. A partial selection turns only that stretch | `rev_com_window` @0xac54a | `Doc::reverse_complement`, `reverse_complement_range` |
| …features in the stretch | Mirrored: new range `last − end … last − start`, strand flag flipped, forward/reverse colours swapped. A feature crossing the edge is **split**: the outside part stays as it was, the inside part becomes a new mirrored feature | `put_tags` @0x947aa, `get_tags_exact` @0x91e69 | `Feature::mirror`, `inside`, `remove_range` |
| Set origin | Text before the caret moves to the end; features move with their bases; one crossing the new origin becomes two ranges; caret at the start does nothing | `set_origin` @0xaba02 | `Doc::set_origin` |
| Circular ↔ linear, Linearize | `toggle_circ` flips the topology (Set Origin and Linearize are disabled when linear); `linearize` = set origin at the caret, then toggle | `toggle_circ` @0xac08b, `linearize` @0xac494 | existing actions |
| Insert Recognition Site | Menu groups of **15** enzymes labelled `first-last`; a palindromic enzyme has one entry, others `(Forward)` and `(Reverse)`. The flat recognition sequence is written with ambiguity codes fixed (N W R D H M V → A, S Y B → C, K → G), reverse-complemented for Reverse. A selection is replaced. A `protein_bind` feature named after the enzyme, colour `#C90000` on both strands, is added over the site (strand reverse for the Reverse entry), and the site is selected. A feature the site lands inside grows with it | `insert_enzyme_fill_cascade` @0x99750, `insert_enzyme_site` @0x99b3c, `split_tags_after_insertion` @0x2c2763 | `Enzyme::site_to_insert`, `Doc::insert_site`, Enzymes menu |
| Copy / cut | Text and its features go to the clipboard together; Copy Rev-Com stores the reversed text with mirrored features | `clip_cut_copy` @0x8e9b0 | `clipboard::copy`, `Doc::features_in` |
| Paste | Spaces, digits and line breaks removed, `*` → N. Letters other than ACGTN (either case) trigger a question — remove them, make them N, or cancel (message shows at most 100 of them, then `...`). A selection is replaced; the pasted text ends up selected. When the pasted text is exactly the copied text (or its reverse complement for Paste Rev-Com) the copied features are put back, mirrored for Rev-Com | `clip_paste` @0x8f756, `atgcfilter` @0x76bfa | `clipboard::paste`, `insert`, `Doc::paste` |
| Pasted records | Text that looks like GenBank (`LOCUS … bp … ORIGIN … //`) or FASTA (`>` then a line of bases) is treated as a file, not inserted. **Difference:** the original asks for a file name and saves it; HelixLoom opens it as a new unsaved document | `clip_paste` @0x8f756, `paste_file` @0x912c2 | `clipboard::looks_like_record`, `App::open_text` |
| Typing and deleting inside features | Both neighbouring bases carrying the same feature tag extend it over inserted text; deleting every base of a feature removes it | `textarea_insert` @0x2c24c5, `textarea_delete` @0x2c2c1e | `Feature::after_insert`, `after_delete` |
| Duplicate Selection / Duplicate Full Sequence (File menu L1077–1078, ⌘D and ⇧⌘D) | A new window from the selected stretch (a **linear** document with the features inside it, cut to it) or, with nothing selected or for the full sequence, from all of it (the circle, every feature, the Dam/Dcm methylation, the file's header kept). The new window is called `New_DNA` and its comment starts with `name from a to b` (the file's name without extension, or the window's name; counted from 1) and a line break, then the old comment. A selection across the origin is joined from both parts | `duplicate_window` L7750 @0x5c69b | `Doc::duplicate`, `App::duplicate_current` (no ⌘D keys: the menu and the palette). **Difference:** no empty line is left after the first line when the old comment is empty |
| New features from found highlighting | With no highlight: message "There is no text highlighted from the find dialog." Otherwise one dialog (default name `New Feature`, type `misc_feature`, that type's colours, qualifiers) and then one feature per highlighted match: forward matches first, then reverse ones, each on its own strand | `feature_from_found` @0x1d2ff5 | Features → New Features From Found Highlighting |

### Translate window (`translation_window` @0x1a8803, L36460), built from the original's own text

* **Dialog.** "Translate:" offers **All**; **Selection** only when something is selected (and then it is the first
  choice shown); and every feature whose type is CDS (any case), by name. Choosing a CDS turns off "Reverse
  Complement" and "Copy Highlight". Code: 1 Letter / 3 Letter; Codon Spacing; Reverse Complement; Line Width; Line
  Numbers None / Left / Right / Both; DNA None / Above / Below. **Defaults** (`ApE_Defaults.txt` L46–54): 3 letters,
  spacing off, reverse complement off, width **20**, numbers none, DNA none. An empty width box means **200000**
  residues per line (a width of 0 would loop forever in the original; here 0 means the same as empty).
* **What is translated.** All: bases 1…n. Selection: its first…last base. With Reverse Complement the stretch is
  reverse-complemented and the first base number moves up by `n mod 3`; without it the last base number moves down by
  `n mod 3` (a trailing partial codon is left out either way). A CDS: its bases joined and reverse-complemented when it
  is on the minus strand, numbered 1…its length. Any codon with a non-ACGT base reads X (`Xxx` in 3 letters), U reads as T.
* **Header.** The file's name (and its place); `From %1$s to %2$s.` (start to end, or end to start when reverse
  complemented) or the CDS's name; `Translation %1$s a.a. MW=%2$s` for the residues **up to the first stop**, the mass
  being the sum of the residue masses (G 57.05, A 71.09, S 87.08, T 101.11, C 103.15, V 99.14, L and I 113.16, M
  131.19, P 97.12, F 147.18, Y 163.18, W 186.21, D 115.09, E 129.12, N 114.11, Q 128.14, H 137.14, K 128.17, R
  156.19, X and `*` 0) plus **18.02**, rounded to 0.01 and printed the way Tcl prints a float (`492.6`, `100.0`);
  then `???(%s extra codons after stop)` when codons follow the stop, otherwise an empty line.
* **Lines.** Residues are joined by `spacing` blanks, where spacing is 1 with Codon Spacing, 2 when a one-letter protein
  sits over a DNA line, 3 for both, otherwise 0; with Codon Spacing the DNA is written codon by codon with a blank after
  each; a one-letter protein over a DNA line gets a blank before and after so each letter stands under its codon's
  middle base. One residue takes `1 or 3 + spacing` characters; a line holds Line Width residues. Left numbers are
  `%7u` and a blank, right numbers a blank and `%-7u`: on a DNA line the number of its first and last base (counting
  backwards when reverse complemented; the last line ends at the last translated base), on the protein line the number
  of its first residue and of its last. A DNA line goes above or below the protein line, and a blank line follows each
  group when a DNA line is shown.
* **Not built:** the date and time line, clicking a line to jump to its bases, "Copy Highlight" (feature colours carried
  onto the DNA line) and the index line (its checkbox is commented out in the original).
* HelixLoom: `src/seq/translation.rs`, `src/ui/windows/orfs.rs`. **Checked** on 21 cases taken from the original's
  own code (all options, a selection, a CDS, reverse complement, empty and very short input), every one identical.

### Copy Special (`clip_copy_special` @0x90351, L12354; menu L1117–1139)

* The text is the selection, or the whole sequence when nothing is selected. Items, in the original's order: **Copy All
  as GenBank** (Cmd+B); **Copy Uppercase** and **… Rev-Com** (only the letters A–Z are kept, so lower case and
  everything else is dropped); **Copy Translated**, **Copy Uppercase Translated** (lower case dropped first), and both
  **Rev-Com** (reverse-complemented before translating); **Copy as FASTA**; **Copy as RNA** (T→U, t→u); **Copy Reverse**
  and **Copy Complement** (both marked "NOT A GOOD IDEA"); **Copy Features as NCBI Bankit table**.
* **FASTA:** `>name first-last` where *name* is the file's name without its extension and the two numbers are the count
  of bases before the selection and the count up to its end (so a selection of bases 11–20 gives `10-20`), then the bases
  on one line. With nothing selected the original fails and copies nothing; here the whole sequence is used (`0-n`).
* **Bankit table** (`bankit_text` @0x8e450, `feature_coords_text` @0x1d05cf): per feature a line `first<TAB>last<TAB>type`
  (1-based), a line `first<TAB>last` for each further segment, and `<TAB><TAB><TAB>label<TAB>name`. A feature on the
  minus strand lists its segments in reverse order with each pair written end first; one across the origin starts with
  the part at the end of the sequence. The original lists features in the order Tk returns its tags; HelixLoom uses the
  feature list order. "Copy All But Selected" is commented out in the original and not built.
* **The complement map** (`com`/`revcom` @0xa4b5a / @0xa49b0): A↔T, C↔G, B↔V, D↔H, K↔M, R↔Y, S, W, N and `*` stay, case is
  kept, and **every other letter, U included, is left as it is**. HelixLoom's `complement` used to turn U into A; it now
  follows the original, which also changes Rev-Com on text that contains U.
* HelixLoom: `src/seq/copyspecial.rs`, `src/ui/clipboard.rs`, Edit → Copy Special. **Checked** on 45 cases from the
  original's own code (every conversion, with and without Rev-Com, on mixed-case, ambiguity-code, RNA and empty text).

### ORFs → Find Next / Find Previous (`next_ORF` @0xa9b56 L15025, `orf_list` @0xaa4a3 L15100), built from the original's own text

* **The list** (`orf_list`): when the sequence is circular it is doubled; then for the top strand and (if asked) for the
  reverse complement **of the doubled text**: the starts are every `ATG`, or with "Stop+1" the base after every stop
  codon (`TGA`, `TAA`, `TAG`; the stop's place plus 3); only starts at or before the sequence length `n` count. For a
  start `b` the frame is `b mod 3`; the stop is the first stop at or after `b` in that frame, `c`; the ORF is kept when
  `c + 2 − b > min`; the entry is `[b mod n, (c+2) mod n]` on the top strand and `[(n−c−3) mod n, (n−b−1) mod n]` on the
  bottom one (first and last base on the top strand, both included).
* **Quirks kept.** Every start is its own entry (an ATG inside another ORF with the same stop is listed too; a circular
  sequence lists an ATG at base 0 twice). `c` is remembered from one start to the next and from the top strand into the
  bottom one: a start with no stop in its frame is judged by the stop looked at last, of any frame. On the top strand
  such an ORF is then listed as `[b, n]` (to the end); on the bottom strand it is listed as if that stop were its own,
  with no frame check. A sequence with no stop codon at all stops the original with a Tcl error; here a message says so.
* **The step.** The search position is the selection's first base plus one for Next and minus one for Previous
  (counting round the sequence), or the caret when nothing is selected. The entries are taken in order of start. Next
  selects the first entry whose start is at or after the position (the first entry when there is none); Previous the
  last entry whose start is before it (the last entry when there is none), so with a selection Previous skips an ORF that
  starts one base before it. The selection is bases `start` to `end` included; an ORF across the origin is selected only
  up to the end of the sequence.
* **Defaults** (section 2): minimum 90, start Met, strand Top. HelixLoom: `src/seq/nextorf.rs`, used by ORFs → Find
  Next / Previous; the ORF list window keeps its own list. **Checked** on 240 cases from the original's own code
  (both strands, circular and linear, both start rules, minimum 0 / 30 / 90, caret and selection, both directions).

### Find (`do_find_mm` @0x9fdda L13981; `find_dialog` @0x9b20f; settings L66920 and L66966–66973), built from the original's own text

* **Dialog.** Pattern box; "case sensitive" (off); "also find rev-com of string" (on); "as amino acid" (not built);
  "as literal" (off); "Allow mismatches:" 0 to **25** (0); "Minimum 3' bases:" All or **10** to **60** in steps of 5 (All);
  Find Next, Find Prev, Highlight All, Close. The wrap checkbox is commented out: Find always wraps. An empty box does nothing.
* **Pattern.** Only the letters a b c d g h k m n r s t u v w y, `*` and the separators `|` `;` `,` are kept; anything else
  typed is dropped. Patterns between separators are tried one after the other at each place; with the reverse complement
  each is followed by its reverse complement; "Minimum 3' bases" keeps only the last N bases of each. Without "case
  sensitive" pattern and text are both upper-cased. A letter stands for these bases (and the lower-case letter for the
  lower-case ones): B = BYSKCGT, D = DRWKAGT, H = HMWYACT, K = KGT, M = MAC, N = NBDHVKMRYWSACGT, R = RAG, S = SGC, V =
  VMSRACG, W = WAT, Y = YCT; A, C, G, T, U stand for themselves (so **U does not match T**); "as literal" switches the
  expansion off. A text base matches when it is among the bases its pattern letter stands for; mismatches are counted up
  to the allowed number; a base missing at the end of the text counts as a mismatch; an empty alternative (`GAATTC|`)
  matches everywhere with length 0.
* **Find Next** starts at the selection's first base plus one, or at the caret, and reads the text from there to the
  end, then the text before it (always on a circular sequence; on a linear one because wrap is on), the join marked by a
  tab that matches nothing but counts as a mismatch. The first place where any pattern matches wins. Message `Found
  forward|reverse match at N.` (N from 1), or `Found no match.`.
* **Find Prev** starts at the selection's first base, or at the caret, and reads backwards. **Quirk:** on a linear
  sequence the base just before the start is never tried (the search starts one base further back and wraps to the end),
  so with `AAAA` and the caret at 3, Find Prev of `AA` selects the match at 1 and never the one at 2.
* **Highlight All** tags every place for every pattern (forward tags, and reverse tags for the reverse complements);
  the message `Found N matches.` counts every pattern that matched at every place. A circular sequence repeats its first
  (longest pattern − 1) bases at the end, so matches there are counted twice; with mismatches allowed a pattern may hang
  over the end of a linear sequence.
* HelixLoom: `src/seq/findmm.rs`, the Find window in `src/ui/windows/find.rs`. **Checked** on 300 cases from the original's
  own procedure (random sequences with ambiguity codes and lower case, random patterns with separators, odd letters and
  empty parts, all settings, all three actions): every one identical. The earlier build searched with its own rules
  (no several-pattern syntax, no "Minimum 3' bases", Find Prev reaching the base before the caret).

### Jump To (`jump_to` @0x976e3 L13169) and Select From-To (`select_from_to` @0x96c8f L13097)

* The number boxes take digits only (`check_char` @0x99474). **Jump To** starts with the caret's position; a number `n`
  puts the caret in front of base `n` (position `n − 1`), and a box that is empty, 0 or 1 goes to the start.
* **Select From-To** starts with the selection's first base and last base, or the caret's base and an empty "to". On OK
  leading zeros are dropped (so 0 means empty); numbers the wrong way round are swapped when both are filled in; an
  empty "from" is 1 and an empty "to" is the end of the sequence. HelixLoom: `src/seq/goto.rs`. **Checked** on 28 cases
  from the original's own lines.

### PCR (`pcr_dialog_do_pcr` @0x2944ee L56473, `search_with_mm3` @0x2a43a7 L57959), built from the original's own text

The rules are in section 4. How they were checked: the original's `search_with_mm3` (L57959–58045),
`pcr_dialog_precompute` (L56880–56947), `pcr_dialog_do_pcr` (L56473–56724), `Tm` (L14730–14776) and
`get_tags_get_coords` (L12755–12813) were cut out of `AppMain.tcl` unchanged and run with the system Tcl 8.5 against
stand-ins for the text widget (positions are numbers, and a position at or before the start means the start, as
`bp2ix` @L61660 does), the message boxes (always OK, every message recorded), `create_window` (records the new
document) and `wm title` ("pDemo.gb"). Random templates of 210–309 bases (linear or circular, with a repeat or a
lower-case stretch now and then) got primers cut from them with tails, mismatches, a degenerate base and lower case,
and reverse primers placed to overlap, touch or surround the forward ones, under varied slider settings. The
original's scan and its product for every forward/reverse pair are the expected output of `src/seq/pcr/tests.rs`:
**280 scans, 107 pairs (65 products, 42 refusals; 15 products longer than the template; 35 with warnings; all four
circular overlap warnings) and 240 feature-copy cases**. Replacing the overlap distance 10 by 11 makes the test fail.

What reading and running it changed compared with the earlier build: touching primers on a line are refused; a circle
makes "round the world" products with the four overlap warnings; a stretch is judged before it is trimmed; a template
position past the end counts as a mismatch; a circle repeats only 200 bases; with no mismatches allowed 4 exact 3'
bases are needed whatever the slider says; the slow route's quirks (a stretch one base short, a site at the origin
listed twice); and the new window with its comment, name and primer features. **Evidence this round** (the seventh
request to decompile): `core_target_open` on `AppMain.tcl` → "a file we don't recognise yet" (`t_4467f482ee05`);
`native_strings_find` → "Unknown file magic"; reopened as raw arm64 (Ghidra `AARCH64:LE:64:v8A`):
`native_strings_find` for `annotate` → 0 strings (`ev_58b0c664…`), `native_functions_list` → the single `entry`
(`ev_01903f84…`), and `native_function_pseudocode` → **function_not_found** for `pcr_dialog_do_pcr` 0x2944ee,
`translation_window` 0x1a8803 and `pcr_dialog_filter` 0x298e61. `native_bytes_find` located `pcr_dialog_do_pcr` 0x2944ee
(`ev_4aef6188…`), `pcr_dialog_filter` 0x298e61 (`ev_fc00b2dd…`), `search_with_mm3` 0x2a43a7 (`ev_d0546989…`),
`pcr_dialog_precompute` 0x29863b (`ev_cc898699…`), `pcr_dialog_scan` 0x297522 (`ev_7e564aff…`), `get_tags2` 0x92a8f
(`ev_3bbc62ce…`) and `get_tags_get_coords` 0x93d6c (`ev_06de2700…`); their text was read and run as described above.

### Sanger trace window (`read_abi_info` @0x1684b9 L30572, `new_abi_window` @0x16b63a L30863 and the `abi_*` procedures), built from the original's own text

Found with `^proc abi_` and read in full. `read_abi_info` and `abi_save_scf` were cut out of the file unchanged (only the
writer's file dialog and message box, which need a screen, were replaced) and run with the system Tcl 8.5.9; `mix_color`,
`abi_scale`, `abi_x_scale` and the drawing code were read and their rules written in Rust. Evidence for the addresses:
`read_abi_info` `ev_5c462b76…`, `new_abi_window` `ev_33578936…`, `abi_draw_phred` `ev_45d2fc37…`, `abi_save_scf`
`ev_a13bdd53…`, `abi_revcom` `ev_b65ae2f4…`, `abi_scale` `ev_339f1a32…`, `abi_x_scale` `ev_04fd17c0…`, `abi_file_info`
`ev_f90a71ab…`, `abi_get_gap` `ev_24e0e5e6…`, `abi_get_seq` `ev_3be14efc…`, `mix_color` `ev_66cb2682…`.

| Rule | What the original does | Address (line) | In HelixLoom |
|---|---|---|---|
| File kind | The first four bytes decide: `ABIF`, or `.scf`; anything else is "Wrong Filetype- ABIF not found". Its list of file endings it opens includes `.abi .ab1 .scf` | `read_abi_info` L30607, L30748, L30846; L67221 | `abif::read`; `.ab1 .abi .scf` go to Signals |
| ABIF directory | The entry count is at byte 18 and the directory's place at byte 26 (the header's own entry); 28 bytes per entry; a directory that runs past the end is refused ("Data Past End of File") | L30611–30620 | `Abif::open` |
| Tags used | `DATA` number 9 and up are the signal of tracks 0 to 3 (1 to 8 are never read: "raw data traces- not read at this time"); `FWO_` is the order of the tracks, **A C G T when the file has none**; `PBAS`, `PLOC` and `PCON` **number 1** are the bases, their places (16-bit signed) and their qualities (the original calls them "edited"; number 2 is read into a set it never shows); `SMPL` (Pascal string); `S/N%` four signed 16-bit numbers; `RUND` and `RUNT` numbers 1 and 2 are the start and end of the run (a year in two bytes, month, day; hour, minute, second); `MODL`; every other string or number tag goes into the SCF comments | L30626–30730 | `abif::parse`: set 1 shown (**earlier: set 2**), default order ACGT (**earlier: GATC**), signal strengths and run times kept. Deliberate differences: the raw `DATA1`–`DATA4` are used when `DATA9`–`DATA12` are missing (the original would show empty traces), set 2 when there is no set 1, and places are read unsigned (a signed read breaks traces over 32,767 scans) |
| Base probability | `PCON` is the phred value of the called base; the per-letter probabilities are 0 except the called letter's | L30678–30692, L30741–30747 | `Trace::quality`, and written that way into SCF |
| Run date and time | Texts "m/d/y" and "h:mm:ss" (minutes and seconds padded to two digits); **"01/01/00" and "00:00:00" when the file gives none** | L30586–30589, L30664, L30668 | `Stamp::date_text`, `time_text`, `Default` |
| SCF reading | Header: after the magic eight numbers (samples, their place, bases, two clip points, the bases' place, the comments' size and place), the version (4 characters), the sample size, the code set and the size and place of private data. Version ≥ 3.00 (compared as a number): the four channels one after another, each sample being the change of the change (two running sums), one-byte samples as signed bytes and any other size as signed 16-bit; then the bases' places (4 bytes each), four arrays of probabilities (A, C, G, T) and the letters. Version 2: the four channels together for each scan, then a 12-byte record per base (place, four probabilities, letter, three spare) straight after the samples (the bases' place is not used). The called base's probability is that of its own letter, 0 for any other. Channels are always A, C, G, T | L30748–30845 | `scf::parse`; **equal to its reader on the original's own export and on two version 2 files (one- and two-byte samples)** |
| SCF export | "Export as SCF" writes the whole trace, not the picked bases. Header of 128 bytes: `.scf`, samples, 128, bases, 0, 0, the bases' place (128 + 8 per sample), the comments' size, the comments' place (+ 12 per base), `3.00`, sample size **2**, code set **2**, 0, the end of the comments, 72 spare bytes. Then each channel in the order A, C, G, T as the change of the change in 16 bits; the places; four probability arrays (the quality under the called base's letter only: N and other letters get 0); the letters; **3 spare bytes per base**; the comments | `abi_save_scf` L31933–32011 | `scf::to_scf`; **byte for byte equal to the original's file for a made-up trace** |
| SCF comments | `STRT=` date, `STOP=` date, `SIGN=` with a letter and a strength for each track in the file's order, then `MODL=` and the other tags, joined by newlines. **The list is built with `concat`, so a run's date and time become two lines** ("STRT=2/19/2009", "1:19:30"); the reader here accepts both layouts | L30737–30739 | `scf::comments_text` (the instrument only of the other tags, where the original lists every string and number tag of the file) |
| Colours | A `green`, C `blue`, T `red`, G the foreground colour (black, white in dark mode), N, **every ambiguity code and the gap `magenta`** (Tk's pure colours) | L67217 | `design::trace_color`. **Deliberate difference:** deeper shades of green, blue, red and magenta so the letters can be read on white; the assignment is the original's (the earlier build used grey for the rest) |
| Traces and letters | One polyline per track through every scan (x is the scan), the span from the tallest to the lowest point of all four **fitted between 35 points under the top and 10 above the bottom**; the bases as letters centred on their places under a **number above every tenth base** (counting from 1) | `new_abi_window` L31078–31111, L31098 | `chromatogram.rs`: one scan is one point at scale 1 and the tallest signal fills the trace area (**earlier: the 90th percentile of the peak heights, times 0.8**) |
| Scales | Vertical (`abi_scale`): the slider's 10 to 50 become **6 − v/10, from 5 to 1**, applied about the baseline, 1 when a window opens. Horizontal (`abi_x_scale`): 5 to 50 become **5.5 − v/10, from 5 to 0.5**, 1 when a window opens | L31809–31835, L30868–30869 | Two sliders, 1 to 5 and 0.5 to 5, both 1 at first (**earlier: 0.3 to 4 and 0.3 to 6, with 1.4 across**) |
| Quality graph ("Toggle Phred") | A filled step graph behind the traces: each base's step runs from halfway to the base before to halfway to the base after (the first from the start of the trace, the last to its end), its **height being the quality out of 100** of the way from the baseline up to the bottom of the letters; the fill is `#8080FF` mixed into the page at 20 % and its edge at 35 % (`mix_color` weights its first colour); on when a window opens | `abi_draw_phred` L31524–31562, L31567–31576 | `chromatogram::quality_steps` and a toggle (**earlier: a small bar under each base**) |
| Reverse-complement | The picture is **mirrored** about the middle of the trace; every letter is complemented with the IUPAC pairs (A↔T B↔V C↔G D↔H K↔M R↔Y, S W N `*` as they are, case kept); the track tagged A is drawn in T's colour, C in G's, G in C's and T in A's; the sequence kept for copying is the reverse complement | `abi_revcom` L31733–31759 | `chromatogram::mirrored` and `shown_base` (**earlier: only the output sequence was reversed**). Deliberate difference: the original leaves each tenth-base number on the same base, so the numbers run backwards after a flip; here they count from the left |
| Selection | A click picks the gap between two bases nearest the pointer (halfway between their places; the left edge before the first, the right edge after the last); shift-click moves the nearer end; the bases taken are those whose place lies between the ends | `abi_get_gap` L31664, `abi_get_seq` L31691 | Bases are picked by dragging (a deliberate difference); the picked range is whole bases |
| File info | "Sample:", "Scaling:" (the four strengths with their letters in the file's order), "Run began:", "Run ended:", the moments as `%l:%M%p %a,%b %d, %Y`; for a file with none of it the defaults above | `abi_file_info` L31840–31861 | The File info section of the page, listing only what the file gives (and a space after the comma) |

**What was run.** `dump_read.tcl` (ephemeral, in `/tmp/kv/trace`) sources the cut-out `read_abi_info` and prints what it
read. On the three real files (a 310, a 3100 and a 3730, from Biopython's test suite) it read the sample names `D11F`,
`16S_S2_1387R` and `226032_C-ME-18_pCAGseqF`, the order `GATC`, 868, 795 and 1,165 bases, signal strengths
`93,134,134,52`, `1142,579,708,704` and `500,388,300,461`, and runs of 2009-02-19, 2010-01-27 and 2009-12-12; in all three
the base calls of set 1 and set 2 are identical, and equal to Biopython's. Every value it read, the four signals
included, equals this program's (`cargo test real_files_read_as_the_original_reads_them -- --ignored`, the folder named
by `HELIXLOOM_AB1_DIR`). For SCF a generator made a trace of 600 scans and 40 bases (one an N); the original's
`abi_save_scf` wrote it (5,490 bytes, `src/seq/scf/golden.scf`), `to_scf` writes the same bytes, and the original's reader
on that file gives the arrays in `golden.read`, which `scf::parse` equals; two SCF version 2 files (`golden_v2_*.scf`),
made by a small script, were read both ways with the same result.

### Calculator and Molecular Reaction Calculator (`calc` @0x18d3a6 L34184, `mol_calc` @0x18dde0 L34255), built from the original's own text

Tools menu, L1293 and L1294 @0x1245a (⌘1 and ⌘2; HelixLoom's workspaces keep ⌘⌥1 to ⌘⌥9, which are read first because a key without Alt would also take them). Procedures: `calc_update` L34204 @0x18d6ad; `calc_molar_reaction` L34777
@0x195d0d; `calc_mm_reaction` L35131 @0x199d47; `calc_molarity` L35237 @0x19aec7; `calc_set_units` L35463 @0x19dae4;
`calc_stock_dilution_update` L35540 @0x19e731; `calc_output` L35616 @0x19f21e; the reactions `info(calc_molar_raction_dict)`
L67338 @0x31222b and `info(calc_mm_raction_dict)` L67339 @0x3127df; the first values `info(calc_defaults)` L67340 @0x312850.
HelixLoom: `src/seq/molcalc/` (the arithmetic; `format.rs` writes numbers as Tcl does) and `src/ui/windows/{calculator,molcalc}.rs`.

* **Calculator** (`calc`, `calc_update`). A line of arithmetic goes to Tcl's `expr` after **a point is put after every run
  of digits, so that every number is real** (`1/3` is 0.3333333333333333). Not after digits that follow a point, precede
  one, follow `E`, `E-`, `E+` or `0x`, precede `x`, or touch `<<`, `>>`, `^`, `|` or `&` (whole numbers are wanted there).
  The result is written as Tcl writes a double, with a trailing `.0` dropped (`6.0` is `6`); a line `expr` refuses gives
  `------` (also an empty line once a key was typed; nothing is worked out before that). `$pi` is `acos(-1)` (L34205) and
  `$e` is `exp(1)` (L34206). Zero to a negative power through `**` is a mistake (`0**-1`), while the function `pow(0,-1)` is
  `Inf`. **Quirk kept:** a hex number with a letter before its last run of digits (`0xD07`)
  gets a point and is a mistake. HelixLoom (`expr.rs`) reads numbers (also `0x`), `$pi`, `$e`, the operators of Tcl in its
  order (`**` right to left, binding less tightly than a sign: `-2**2` is 4), `?:`, and the functions `abs acos asin atan
  atan2 bool ceil cos cosh double entier exp floor fmod hypot int isqrt log log10 max min pow round sin sinh sqrt tan tanh
  wide`; whole numbers are 128-bit. **Not built:** commands in `[…]`, quotes, `eq ne in ni`, `rand`; whole numbers beyond
  128 bits (Tcl's go on); and a number typed as the chosen branch of `?:` in an expression that is not all constants,
  which Tcl gives back as typed (`4.`, `2.25e-4`). **Limits of HelixLoom's own** (estimated, section 9): a line over 1,000
  bytes, or nested more than 100 deep (brackets, call arguments, signs, powers, conditions), gives `------`, so that a long
  run of brackets pasted into the window cannot use up the stack.
* **Molar ratio reaction** (`calc_molar_reaction`). Up to 10 fragments (name, size in kb, ng/µl, ratio); a ratio left
  empty is 1 once there are a size and a concentration. A kilobase of double-stranded DNA is **617.96 ng per pmol and the
  ends add 36.04**: `grams = size × 617.96 + 36.04`. For the final amount in **fmole** (the total, for all fragments)
  volume = `final/1000 × ratio/Σratio × grams / conc` and fmoles = `1000 × volume × conc / grams`; in **nM** the same
  times the reaction volume; in **ng** and **ng/ul** the volumes are first `ratio/Σratio × size / conc`, then scaled so
  that the ng add up to the final amount (times the reaction volume for ng/ul). The enzyme mix is `1/X` of the reaction
  and the DNA mix the rest (rounded to hundredths); `X` of 1 or less is DNA only. **Water** (ticked at the start) makes
  up the DNA mix volume; with too little DNA, or without water, the volumes are scaled to fill it and the line says
  `-Input DNA concentrations are too low (n %)` (and the label "Final DNA amount/concentration" turns red, `conditions_frame.l1`
  L34970-34978) or `(n %)` above 100; water under 5 % of the DNA is left out. If
  the smallest volume is under the **smallest pipette volume** everything (volumes, fmoles, the DNA mix, the water) is
  scaled up so that it is not. Volumes and amounts are written with two decimals. **The reactions** (conditions only, volume
  in µl, final amount, units, smallest pipette, enzyme X): HiFi 2-3 fragment 4, 40 fmole, 0.3, 2; HiFi 4-6 fragment 4,
  60 fmole, 0.3, 2; Gibson 2-3 fragment 4, 40 fmole, 0.3, 1.333; Gibson 4-6 fragment 4, 60 fmole, 0.3, 1.333; SapTrap 2.5,
  30 nM, 0.3, 1.25; Golden Gate 10, 6 nM, 0.3, 5; T4 Ligation 10, 10 ng/ul, 0.3, 5; NEB Quick Ligation 20, 10 ng/ul,
  0.3, 2 (each with a note, kept). The original opens on HiFi 2-3 fragment, and **never loads the fragments stored with a
  reaction** because "Preserve Values Between Reactions" is on and nothing turns it off, so choosing a reaction changes the
  conditions only. The enzyme stock menu: 10x (1:9), 5x (1:4), 2x (1:1), 1.33x (3:1) = 1.333, 1.25x (4:1), DNA Only = 1, or
  a number (starting at 10). Copy writes the text of `calc_output` (tab 2).
* **Master mix multiplier** (`calc_mm_reaction`). Phusion: primer .25, primer .25, enzyme .5, dNTP 1, water 38, buffer 10
  (50.0 µl in all). Typing a component's volume makes the single reaction volume `round(100 × Σ) / 100`; typing the single
  volume scales every component to it; each component of the mix is `volume × reactions × (1 + extra/100)` for an extra in
  %, or `× (1 + extra / (single volume × reactions))` for an extra in µl, rounded to hundredths; the total is their sum.
  **Deliberate difference:** with an extra in µl the original divides whole numbers when the three are typed as whole
  numbers, so 20 µl extra for 12 reactions of 50 µl is 20 / 600 = 0 and nothing is added; HelixLoom divides real numbers.
* **Molarity** (`calc_molarity`, `calc_set_units`). `first unit /volume = second unit /volume`, written with `%.3g`; the
  units are g mg ug ng pg, mole mmole umole nmole pmole (and `%(wt/vol)` = 10 g/L for the first), volumes /L /mL /uL. A
  molecular weight is needed between mass and moles: typed in g/mol, or for **nucleotides ssDNA** A 313.23, C 289.20,
  G 329.23, T 304.21 and 18.02 for the ends (without a sequence 308.97 × length + 18.02), **base pairs dsDNA** AT 617.44,
  GC 618.22 and 36.04 (617.96 × length + 36.04), **nucleotides RNA** A 329.2, C 305.2, G 345.2, U 306.2 and 18.02
  (321.47 × length + 18.02); a typed sequence is cleaned to ACGT (U for RNA) and its length is the length. The second line
  gives what a volume (or an amount) comes to, with the prefix `calc_set_units` picks (`floor(ln x / ln 10)`: p below −9, n
  below −6, u below −3, m below 0): `1 mg , 1 umole` on the first screen (1 mmole/L = 1 g/L at 1000 g/mol, 1 mL).
  **Deliberate difference:** for moles per litre alone the original divides by the volume where every other case
  multiplies (1 mmole/L in 0.5 L is "2 mmole") and its other two branches test the amount where they mean its unit, so
  they never run; HelixLoom multiplies, and divides the moles by the concentration for a volume.
* **Stock dilution** (`calc_stock_dilution_update`). `stock volume = final volume × final conc × 1000^(steps between the
  units) / stock conc`, units in the order M mM uM nM pM X; water is the rest; a volume under 1 is written in the next
  smaller unit (L, mL, uL) to hundredths; a final concentration above the stock's is `Error`. A stock in X keeps the final
  in X and the other way round.
* **Not built:** the molecular weight of a compound found by name (asks PubChem), the recipe multiplier (reads the clipboard
  and asks PubChem), saving a reaction of your own, dropping gel bands into the table, the QR code, the µM unit of the
  final amount (in the formulas but not in the menu).
* **Checked** with the original's own procedures, cut out unchanged and run with `tclsh` 8.5.9 on random inputs, widgets stood
  in for (`golden/*.txt` hold inputs and what the procedures printed): 250 molar ratio reactions, 150 master mixes (118
  compared, 32 are the whole-number extra above), 300 conversions (those of moles per litre alone excluded), 150 dilutions,
  40 prefixes, 40 `%.3g`, and 3,000 lines of arithmetic through `calc_update` (the few that are bigger than 128 bits or
  typed numbers are left out). Hand-worked: two fragments of 3.5 and 1.2 kb at 80 and 45 ng/µl in 3:1 (40 fmole, 4 µl, 2x)
  need 0.30 and 0.57 µl.

### Silent Sites and Find Diagnostic Site (`silent_sites` @0x1a0b2c L35791, `do_diagnostic_site` @0x1a3e5b L36062), built from the original's own text

Enzymes menu, L1194 @0x101e1 and L1195. Other procedures: `coordinate_translation` L35922 @0x1a24be;
`add_daignostic_site_dialog` L35993 @0x1a32b3; `degen_transl` L64442 @0x2ea304; `triplet_transl` L64424 @0x2e9fae;
`seq_len_mod` L58532 @0x2aa1ab; the start values `info(enz_make_diag_enz)` L66963 and `info(enz_make_diag_mm)` L66964.
HelixLoom: `src/seq/{silent,diagnostic,sitesearch}.rs` and `src/ui/windows/sitetools.rs` (Enzymes menu, command palette,
`--view silent` and `--view diagnostic`).

* **Silent Sites.** The DNA searched is the selection or all of it (the original also takes a selection across the origin;
  HelixLoom has no such selection); the enzymes are the ones chosen in the Enzyme Selector, or **all of them when none is
  chosen** (a message box in the original, a line in the window here). Each enzyme's flat recognition sequence (`flatpat`: the
  `N`s at both ends cut off, `read_enzymes` L64461) is read in its three frames with two `N` on each side (`degen_transl`);
  each codon becomes the set of residues it can code for in the **standard code** (`triplet_transl` has a table of its own; all
  21 residues is "any"), and the three lists of sets are looked for in the translation of the DNA (in the genetic code
  chosen) with `regexp -all`, so **matches do not overlap**. The frames are tried in the order 2, 1, 0 (the bases of the first
  codon that come before the site). A place whose site start already is a site of that enzyme anywhere in the sequence (both
  strands, as the Enzyme Selector counts them, so Dam and Dcm blocking applies) is dropped. The rest are sorted by the start of
  the site, ties keeping the order enzyme, frame, position.
* **The changed bases** (`coordinate_translation`). The site, padded with `N` (the frame in front, to whole codons behind), is
  laid over the codons one by one: a codon of plain bases is taken as it is; `NNN` keeps the DNA's codon; otherwise each
  letter keeps the DNA's base when it allows it and takes the site's letter when not. When that codes for another residue, the
  first codon in the original's list for the residue (L35946; **its stop list holds TGG, which codes for Trp: kept**) that the
  site's letters allow is taken, its ambiguity codes set back to the DNA's bases where they allow. A residue with no codon in
  the list (an ambiguous codon of the DNA is `X`) gives `xxx`. An ambiguous result is spelled out with every base it allows
  (`S` is `GC`, not `CG`: the order of the original's table), **the first letter changing fastest**.
* **The window.** The date, the name, `From a to b.` (from 1), `Silent Sites`, a blank line, then for each place: its
  position (the first base of the first codon, through `seq_len_mod`) and the DNA's codons; a line for each variant, a changed
  base in red (`<…>` in the copied text); the enzyme and its site cut into threes (three spaces after the site, so the line
  ends in spaces); the residues three spaces apart; a blank line; `No sites Found` when there is none. A click selects the codons
  (in the original it goes to them through `alignlink2`).
* **Find Diagnostic Site.** Options: **selected** or **absent** enzymes (absent: no site in the whole sequence) and **1, 2 or 3
  mismatches**; at the start, selected and 1. Enzymes that need methylated DNA are left out (the 16 names of L36099, which the
  original lists in other loops too). Each enzyme's flat sequence is tried on the forward strand and, unless it is its own
  reverse complement, on the reverse; at every place of the DNA (selection or all; upper case) the letters its pattern letter
  does not allow are counted (an ambiguity code allows itself, the codes inside it and the bases it covers: `B` is `BYSKCGT`;
  a plain base only itself), and a place with at least 1 and at most the chosen number is listed, sorted by position (stable).
  The window shows the position from 1, the DNA, the DNA with the pattern's letter written in red where the DNA does not fit
  it, the enzyme and the pattern. **Quirk:** the heading counts the stretch from 0 (`From 0 to 156.`) while the positions under
  it count from 1; HelixLoom's heading counts from 1.
* **Checked** with the original's own procedures, cut out unchanged and run with `tclsh` 8.5.9 on random inputs with a
  stand-in for the text window (`silent_sites`, `coordinate_translation`, `degen_transl`, `triplet_transl`, `translate`,
  `read_enzymes` on Default_Enzymes.txt, `findenzymes`, `do_diagnostic_site`, `enz_site_usage_filter`, `enz_count`): 80 silent
  cases (429 places, 23 with several variants; linear and circular, with a selection, with lower case and an `N` in the DNA)
  and 30 diagnostic cases (selected and absent, 1 to 3 mismatches). HelixLoom writes **the same window text** in every case
  (`src/seq/silent/golden.txt`, `src/seq/diagnostic/golden.txt`), and the sites it finds for all 116 enzymes in the 110
  sequences equal the lists `findenzymes` made (the `sites` lines of the same files).
* **Differences kept.** A row you can click instead of text (the whole row selects its bases; searching again after such a
  click goes on with the stretch the list was made from, while a stretch you select yourself is the one searched); the copied
  text has `HelixLoom` where the original writes the date; the chosen enzymes keep the order of the enzyme list (the original
  sorts the chosen names with `lsort -unique`, L16209, which is name order too, as the list is); a selection across the
  origin is not searched. **Own limits (estimated):** 100,000 places, 1,024 variants of a place, 1,000 rows drawn (Copy as
  text has all of them). A search runs in the window's thread without a way to stop it: measured in a release build with all
  116 enzymes, 1 Mb of random DNA takes 0.2 s (it stops at the 100,000 places) and 5 Mb of one repeated base takes 1 to 3 s,
  after the sequence's own site list (2 to 3 s for 5 Mb, which the Enzyme Selector needs as well).
* **Letters the original does not keep.** `read_enzymes` keeps only `A C G T N W S Y R K M B D H V` of a pattern (a `U` goes)
  before it cuts the `N`s off the ends, as `flat_pattern` does. A `*` or `?` in the DNA matches whatever the site says
  (`string match "*$orig_char*"`), so `substitute` keeps it.

### List Enzymes and Quick Lists (`list_enzymes_text` @0xb8d47 L16463, `list_enzymes_window` @0xb82e9 L16378, `list_enzymes_dialog` @0xb7316 L16305), built from the original's own text

Enzymes menu, L1184 @0xfe05 and L1186–1189 (`findenzymes $w; list_enzymes_window $w <format> <list>`). Other procedures:
`enz_site_usage_filter` L15931 @0xb4584, `enz_count` L15950 @0xb4790; the start values `info(enz_list_enz)` "all",
`info(enz_list_mode_number)` 0 and `info(enz_list_columns)` 6 (L66960–66962). HelixLoom: `src/seq/enzlist.rs` (the lists),
`src/ui/windows/enzlist.rs` (dialog and window; Enzymes menu, command palette, `--view enzlist` and `--view enzlist-dialog`).

* **The dialog** (List Enzymes…): *Format* (Name · Name, #Sites · Name, Site, #Sites, Positions · Site, Name, numbered 0 to 3),
  *Columns* (1 to 10, turned off for the last two formats), *Enzymes* (All · Selected Enzymes Only · Absent Enzymes · Unique
  Enzymes) and *DNA* (All · Selected DNA only); OK or Return lists, Cancel does nothing. The format, the columns and the
  enzymes keep what was chosen last; DNA starts from the Selection Only switch (all). **Quick Lists** are the same code with
  the choices fixed: Site Usage = format 1, all; Unique Sites = format 0, unique; Absent Sites = format 0, absent; their columns
  are the dialog's (6 at the start); their DNA is the Selection Only switch's, which HelixLoom does not have (all of it).
* **Which enzymes.** *All*: every enzyme of the file in the file's order. *Selected*: the ones chosen in the Enzyme Selector
  (`info(enz_currently_selected)`, which is kept sorted by name); none chosen lists nothing. *Absent* and *Unique*: the enzymes
  whose number of sites (`enz_count`) is 0 or exactly 1 on the DNA looked at. With a selection `enz_count` counts the sites
  numbered `start + 1` to `end` (`start` bases come before the selection, `end` is the last selected base).
* **The text** (`list_enzymes_text`). Each enzyme's sites are the ones `findenzymes` found (numbered from 1, a site that
  matches on both strands twice), kept when **`start` ≤ site ≤ `end`** with a selection, all of them without one. **Quirk: this
  is one base earlier than `enz_count` counts**, so with a selection a site that begins on the base just before it is listed
  although the enzyme was picked as absent or unique without it. Every field is padded on the right (`%-*s`) and none is cut
  (XcmI's pattern is 16 letters, wider than its field of 15).
  * Format 0 (Name): sorted by name (`lsort`, ASCII order, case counted), rows of *Columns* cells, each the name padded to 15.
  * Format 1 (Name, #Sites): the same, each cell the name padded to 12 and the number of sites padded to 8.
  * Format 2 (Name, Site, #Sites, Positions): `lsort -integer` by the number of sites, fewest first, **equal numbers keep the
    order they came in** (the stable sort); a line is the name (15), the pattern as the enzyme file writes it (`G^AATTC`,
    `GGTCTC(1/5)`; 15), the number (7) and the sites separated by single spaces (7), then a line break.
  * Format 3 (Site, Name): sorted by the pattern, the pattern (15) and the name (15), then a line break.
  * **Quirk of the columns:** the rows are `ceil(enzymes / columns)`; the line break belongs to the format of the last column
    of a row, so **a row that is not full has none**: a list whose last row is partly empty does not end with a line break.
* **The window** (`list_enzymes_window`). Its name is the sequence's name and the list's: `All Enzymes`, `Selected Sites`,
  `Absent Sites`, `Unique Sites`. Its text starts with the date, the sequence's name and a space, `From <start> to <end>. `
  (counted from 0: `From 0 to 2724.` for all of the DNA; with a selection of bases 11 to 20, `From 10 to 20.`), a blank line, the
  list's name, a blank line, and then the text. The original makes a text window 120 characters wide and sizes it to the longest
  line (`textbox_geometry`, L16817); it is read-only and has the output menu.
* **Checked** with the original's own procedures, cut out unchanged with `tests/tcl/cut_procs.tcl` and run with `tclsh` 8.5.9
  on random inputs (`tests/tcl/enzyme_list_driver.tcl`): `list_enzymes_window`, `list_enzymes_text`, `enz_site_usage_filter`,
  `enz_count`, `findenzymes`, `initialize_methylation`, `not_in_list`, `toggle_methylation`, `read_enzymes` on
  Default_Enzymes.txt, with a stand-in for the text window only. 72 cases: sequences of 52 to 238 bases with 4 to 17 sites of
  random enzymes put in (some twice, on either strand), 34 linear and 38 circular, Dam/Dcm methylation off in 42 and on in 30,
  51 with a selection (20 of them asking for the list of all of the DNA) and 21 without one (some of those asking for the
  selection, which then falls back to all), each of the four formats in 13 cases or more, each of the four lists in 18, 1 to 10
  columns, 0 to 13 enzymes chosen in the Selector. HelixLoom writes **the same window text** in every
  case (`src/seq/enzlist/golden.txt`), and the sites it finds for all 116 enzymes equal the original's lists in number and
  order. The pattern HelixLoom prints for each enzyme equals what the enzyme file gives (`enz_pattern`). The same two
  comparisons on files made with other seeds are an ignored test (`HELIXLOOM_ENZLIST_GOLDEN=file cargo test enzlist --
  --ignored`): 500 more cases, and 1,900 made with the driver's fourth argument 1, which puts in a site repeated at its own
  period so that sites overlap as far as they can (this is what showed the two site-finding rules of section 2): all equal.
* **Differences kept.** One result window that the next list replaces (the original opens a new one each time), with a *List
  again* button and a note when the sequence has changed since; the copied text has `HelixLoom` where the original writes the
  date; a window opens as wide as its longest line when the screen allows, and keeps its place and size while the same kind
  of list is made again; the text is read-only and *Copy as text* takes all of it. **Not built:** a selection across the
  origin (the original reads a second range, L16406–16413, L16483–16485 and in `enz_count`, L15960–15963; HelixLoom has no such
  selection), the Selection Only switch that Quick Lists read, the Enzyme Selector's own List menu (L15509–15512: the four
  formats for the selected enzymes, which List Enzymes… with *Selected Enzymes Only* gives) and its Text button (L15531,
  L15676), and the output menu of the window (save, print).

## 9. Added by the person during the build (not in the original)

The person asked for primer3, MEGA 12.1, BEAST2, PhyloCNN, MrBayes and IQ-TREE to be rebuilt as parts of the platform,
for more tool families (MAFFT / MUSCLE5 alignment, Minimap2 / SAMtools / GATK, fastp / FastQC, ETE / DendroPy,
ESMFold / OpenFold / NGLView, Nextflow / Snakemake) and for the heavy computing to run on GPUs through Rust tensor
libraries (candle, burn, tch-rs, dfdx), and then on a workstation with several GPUs that a Mac, Windows or Linux
computer is the interface to. None of these exist in ApE, so none has an original address. Their status is
kept in the Craft Anything unknowns list (scope question, still open: which engine first, how deep, and whether
clean-room re-implementation is acceptable under each project's licence).

### Built (each re-implemented from reading the project's rules; no code copied)

* **primer3 → Tools → Primer3 Pair Design** (`src/seq/primer3.rs`, `src/ui/windows/primer3.rs`). Defaults, the melting
  temperature and the penalties were read from primer3's `libprimer3.cc` (`pr_set_default_global_args_1/_2`, `p_obj_fn`,
  `obj_fn`), `oligotm.c` (`oligotm`, `divalent_to_monovalent`, `symmetry`) and `dpal.c` (the dimer alignment).
  Rules: oligos 18–27 bases (optimum 20), Tm 57–63 °C (optimum 60), GC 20–80 %, no run longer than 5 of one base,
  own dimer scores ≤ 8 (any) and ≤ 3 (3' end), pair dimer scores ≤ 8 and ≤ 3, product 100–300, 5 pairs returned. An
  oligo's penalty is |Tm − 60| + |length − 20| and a pair's is the sum of its two; a target must lie inside the product
  and primers may not overlap it. Tm: SantaLucia 1998 nearest-neighbour table, terminal corrections (AT ΔH −23, ΔS −41;
  GC ΔH −1, ΔS +28), symmetry correction +14, salt term `0.368·(n − 1)·ln([mono + 120·√(div − dNTP)]/1000)`, strand
  factor 4e9 (1e9 for a self-complementary oligo); defaults DNA 50 nM, monovalent 50 mM, divalent 1.5 mM, dNTP 0.6 mM.
  Dimer alignment (`dpal` defaults): match +1, mismatch −1, a one-base gap −2 and nothing longer; "any" is the best
  local score, "end" the best alignment that starts at either sequence's start and finishes on the 3' base.
  **Checked against primer3 itself** (built from its source and run): the melting temperature of 10 oligos to 1e-8 °C;
  the five best pairs, with positions, penalties, product sizes and dimer scores, identical on two templates, with and
  without a target. Not ported: the thermodynamic dimer and hairpin checks (thal), template mispriming, internal
  oligos, excluded regions, several product-size ranges, circular templates.
* **MAFFT / MUSCLE → centre-star multiple alignment** (`src/seq/msa.rs`). The sequence sharing most 6-mers is the
  centre; every other sequence is aligned to it with the pairwise aligner (section 4, free ends) and the alignments
  are merged so a gap opened for one sequence stays open for all. A small stand-in that suits clearly related
  sequences; MAFFT's and MUSCLE's own algorithms are not used. While testing it a **bug in the pairwise aligner** was
  found and fixed: when a gap in the second sequence was extended, the traceback decoded the wrong state and wrote a
  garbled alignment (the score was right); Align Sequences benefits too, and a test now guards it.
* **MEGA / IQ-TREE → Tools → Phylogenetic Tree** (`src/seq/tree.rs`, `src/seq/likelihood.rs`,
  `src/ui/windows/tree.rs`): the open sequences are aligned, then a tree is built by
  neighbour joining (Saitou–Nei), UPGMA, or maximum likelihood, drawn with branch lengths to scale and as Newick text.
  Distances: p-distance, Jukes–Cantor `−3/4·ln(1 − 4p/3)` or Kimura 2-parameter
  `−1/2·ln(1 − 2P − Q) − 1/4·ln(1 − 2Q)` with gaps and ambiguity codes skipped (pairwise deletion). Bootstrap support:
  the alignment's columns are resampled with replacement (default 100 times, a fixed seed so results repeat) and each
  inner branch scores the share of resampled trees that contain its split. Maximum likelihood: Felsenstein pruning for
  JC69, K80, F81 and HKY85 on unique site patterns, branch lengths and kappa by golden-section search, the model
  chosen by BIC, then nearest-neighbour-interchange hill climbing from the neighbour-joining tree; no rate variation
  between sites, no partitions, no ultrafast bootstrap. All defaults here are **estimated** (Kimura 2-parameter, 100
  replicates). The four models are fitted side by side on all CPU cores. Tests: an additive matrix is reproduced
  exactly, distance formulas match hand-computed values, the likelihood of a small tree equals an independent
  brute-force sum (computed separately), a wrong topology is repaired by the search, data simulated under a known tree
  gives back its split and branch lengths.

* **MrBayes / BEAST2 → Tools → Phylogenetic Tree → Bayesian** (`src/seq/likelihood/bayes.rs`). A Metropolis–Hastings
  sampler over the topology (swaps of subtrees across an inner branch), every branch length (multiplier moves) and the
  transition/transversion ratio under HKY85 with the alignment's base frequencies. Priors: all topologies equally
  likely; each branch length exponential with mean 0.1 (MrBayes' default `exp(10)`, see `src/command.c` L11035 of
  MrBayes); kappa uniform in log. Two independent runs (one per CPU core, the second starting from a scrambled tree)
  are compared by the **average standard deviation of split frequencies** (below 0.01 = the runs agree, MrBayes'
  guideline); the result is the majority-rule consensus with posterior probabilities as branch labels. Defaults are
  **estimated** (10,000 generations, a tree every 10, the first 25 % discarded). **Checked the way a sampler should be:**
  with no data at all it returns its priors (the three topologies of four sequences a third each, branch lengths
  averaging 0.1), and on simulated data it finds the true split with > 95 % support and the true inner branch length.
  Not included: Metropolis coupling, relaxed clocks, tip dating, the coalescent and the other BEAST2/MrBayes models,
  XML input.
* **fastp / FastQC → Tools → Read Quality and Trimming** (`src/seq/fastq.rs`, `src/ui/windows/reads.rs`). FASTQ
  (Phred+33) is read; the window shows reads, bases, lengths, GC %, mean quality, Q20 and Q30 shares and the mean
  quality along the reads, before and after. Cleaning follows fastp's source (`options.h`, `filter.cpp`,
  `adaptertrimmer.cpp`): a base below Q**15** is low quality; a read with more than **40 %** low-quality bases, more
  than **5** N, or shorter than **15** bases is dropped (in that order of checks); an adapter is found by sliding it
  along the read from 4, 3 or 2 bases before its start (adapters of at least 16, 12 or 8 bases), needing **4** bases of
  overlap and allowing one mismatch per **8** bases compared, and the read is cut where it starts (a read that begins
  inside the adapter is dropped); optional sliding-window cuts (window **4**, mean quality **20**) from the front, the
  tail or "right" exactly as `trimAndCut` does them. fastp has no default adapter (it detects one); the Illumina
  TruSeq read-1 adapter is offered as an **estimated** default. **Checked with** fastp's own adapter self-test. Not
  ported: paired reads, adapter detection, the one-gap adapter fallback, poly-G/X trimming, the low-complexity filter.
  Cleaned reads can be saved as FASTQ or handed to the Contig Assembler.
* **PhyloCNN → Tools → PhyloCNN Epidemic Estimates** (`src/seq/newick.rs`, `src/seq/phylocnn.rs`,
  `src/ui/windows/phylocnn.rs`, the trained networks in `assets/phylocnn/`). Perez & Gascuel (2025),
  github.com/manolofperez/phyloCNN. A rooted binary tree in Newick is encoded as one row of **18** numbers per node,
  root first, level by level: distance to the deepest tip, distance from the root, branch length (0 for the root),
  tips below, distance to the grandparent, branch length and tip count of the minor child (fewer tips; on a tie the
  second child) and of the major child, distance to the great-grandparent, and for each child's own two children their
  distance from the node and their tip counts; −1 where a relative does not exist. First every branch length is
  divided by the mean branch length (the root counting as 0) so the mean becomes **1**; a root with one child is
  dropped, and a tree that is not binary is refused. The network input is **500** rows × 19 features × 2 channels:
  tips and inner nodes each sorted by distance from the root, the inner nodes followed by a copy of the last one, both
  zero-padded to 500 rows, with the sampling probability as the 19th feature of every real row. Network: a grouped
  convolution (1×19 window, 32 filters, 2 groups, so tips and inner nodes are first treated apart), batch
  normalisation, two 1×1 convolutions of 32 filters each followed by batch normalisation (all ELU), the average over
  the 500 rows, dense layers of 64, 32, 16 and 8 units, then 2 / 3 / 4 outputs for BD / BDEI / BDSS or three
  probabilities (softmax) for the choice between them. Outputs: R0, infectious period (times the rescale factor, so in
  the tree's own time units), incubation period (BDEI, also times the rescale factor), X_SS and f_SS (BDSS). The
  training outbreaks (their README) had R0 1–5, infectious period 1–10, incubation 0.2–5, X_SS 3–10, f_SS 0.05–0.2,
  **200**–500 tips and sampling probability **0.01**–1; the HIV set uses priors the authors chose for the Zurich data
  and was trained only for model choice and BDSS. The window warns when a tree has fewer than 200 tips.
  **Checked by running the authors' own software.** (1) Their encoder (`Encoding/PhyloCNN_Encoding_PhyloDyn.py`,
  ete3) on the Zurich HIV tree agrees with this program on all 399 × 18 values to 1e-7 (float32 rounding), and its
  rescale factor, 5.174108653956917, is the one their CI notebook states. (2) All six trained networks were rebuilt
  layer by layer in Keras 3 from their JSON and H5 files and run on exactly the input this program builds: the
  outputs agree to within 4e-6 (relative; the smallest probability, 4e-4, within 7e-5). (3) Their encoder and Keras on
  the Zurich tree give 1.40196, 9.1154, 19.732 and 0.12139 for R0, infectious period, X_SS and f_SS, this program
  1.40044, 9.1208, 19.707 and 0.12167: within 0.25 %. The two channels are paired row by row, so the order of
  equal-age tips (an unstable sort in their code, a stable one here) moves the answer by about a tenth of a percent.
  **Not reproduced:** the numbers printed in their `CI_HIV.ipynb` (1.412436, 11.32054, 20.390572, 0.129327). Their own
  encoder and the shipped model files give the values above with or without the rescale row that notebook leaves in
  the input, so the printed output must come from other files. Their training notebooks drop that row; so does this
  program. Not ported: training, the diversification (BiSSE) networks, the one-generation and no-context networks,
  and the confidence intervals and temperature calibration of their other notebooks. **Licence:** the repository has
  no licence file, so its trained weights (converted to a plain binary file) and the Zurich tree (Rasmussen et al.
  2017) are used here only because the person asked for this rebuild; ask the authors before sharing the program.
* **ETE Toolkit / DendroPy → Tools → Tree Tools (Newick)** (`src/seq/newick.rs`, `src/seq/treetools.rs`,
  `src/ui/windows/treetools.rs`). Newick is read (several trees, quoted names, `[...]` comments; a missing length is 0
  for the root and 1 elsewhere, as ETE reads it) and written back (lengths to twelve decimals; characters Newick uses
  become `_`). **Midpoint root:** the two tips farthest apart are found by two sweeps over the unrooted tree, the root
  goes exactly halfway between them, a root that only joined two branches disappears, and branch labels such as support
  values stay with their branches. **Ladderize:** smaller clades first at every node, equal sizes keep their order
  (ETE's `ladderize`). **Describe:** tips, inner nodes, total length, height, binary or not. **Distance table:** the
  path length between every two tips as tab-separated text. **Robinson–Foulds:** the splits found in only one of two
  trees on the same tips, ignoring the roots; the most it could be is the splits of both added together. **Checked
  against the real tools** on five random ten-tip trees: the ladderized trees, the tip-to-tip distances and the
  Robinson–Foulds distances (12, 10, 6 and 8, each of at most 14) are ETE3's, and the midpoint-rooted trees are
  DendroPy's (every length to 1e-9). Not built: pruning or keeping chosen tips, collapsing weak branches, drawing and
  editing trees, NHX annotations, dating, and the rest of both toolkits.

* **GPU (candle) → PhyloCNN on many trees** (`src/seq/phylocnn/gpu.rs`, `src/seq/phylocnn.rs`, the PhyloCNN window).
  The six networks run through **candle** 0.9 (a Rust tensor library, a default cargo feature `gpu`; `cuda` is a second
  feature, not tried here): the grouped convolution, batch normalisation, ELU, average and dense layers as tensor
  operations, 128 trees at a time, on the Apple GPU through Metal (or CUDA when built for it) or on candle's own CPU
  code. Burn, tch-rs, dfdx and rusttorch are not used: one library was enough to try. The PhyloCNN window now takes several
  trees (a posterior sample), estimates all of them and gives mean, median and the middle 95 % of each estimate.
  **Checked on the real GPU of this Mac:** all six networks give the plain-Rust answers on the Apple GPU to 1e-3 and on
  candle's CPU path to 1e-4 (`cargo test the_gpu_gives -- --ignored`). **Measured honestly** (release build, 1000
  trees of 200 tips, BDSS network): the plain-Rust code on every CPU core takes **0.20 s**; the GPU takes **1.04 s**
  after a **13.5 s one-time build of its shader programs**; the two agree to 3.6e-6. For networks this small the GPU
  does **not** beat the CPU, which is why the window offers it as a choice, with CPU the default, and a single tree always
  runs on the CPU. Everything else that is heavy (multiple alignment, the four model fits, Bayesian runs) runs on all CPU
  cores through `src/seq/par.rs`; those are loops over small matrices that would not gain from a GPU either.
* **MAFFT, IQ-TREE, primer3, Snakemake, Nextflow → Pipelines workspace** (`src/pipeline/`, `src/ui/platform/pipelines.rs`;
  Tools → Pipelines opens it, and it replaced the earlier floating Workflow window).
  The open sequences (or six made-up isolates, for the example) go through **Align**, **Tree** and **Primers** in one run. Each step calls the installed program when
  there is one (found in PATH and in Homebrew's and MacPorts' folders, because a program started from the Finder has a
  short PATH) and the built-in engine when there is not: `mafft --auto --quiet`, IQ-TREE 3 with ModelFinder (`-m MFP`,
  `-seed 1`, optionally `-B 1000`, the least IQ-TREE accepts), and `primer3_core` on a Boulder-IO input made from the
  Primer3 window's settings; the runs happen on a worker thread so the window keeps drawing. **Save as Snakemake and
  Nextflow files** writes `sequences.fasta`, `primer3_input.txt`, a `Snakefile` and a `main.nf` with the same commands.
  The page shows MAFFT, IQ-TREE, primer3, Snakemake and Nextflow as found or not (Home also lists fastp and curl; fastp is
  not called: the read cleaning is built in). Each result is a card that opens in its own workspace: the alignment in
  Alignment, the tree in Trees. Primers are picked on the first open sequence around **its own** selection (the earlier
  window applied the selection of whichever sequence was showing to the first one, a slip that is fixed).
  **Checked with the real programs on this Mac** (MAFFT 7.526, IQ-TREE 3.1.4, primer3 2.6.1 built from its source,
  Snakemake 9.27, Nextflow 26.04): both exported workflows run to the end and write byte-identical alignment, tree and
  primer files; MAFFT's alignment of the demo sequences is the same as the built-in one; on the demo sequences IQ-TREE's
  tree and the built-in search's tree differ by one split (Robinson–Foulds 2 of 6) but have **the same likelihood under the
  built-in models** (−1129.516), so the search is not worse; primer3 and the built-in rules pick the same pairs (same
  primers, product sizes, Tm to 0.05, penalty to 0.01) once primer3 is told to use its score-based dimer limits
  (`PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0`): **without that tag primer3 2.6 uses thermodynamic alignment, with other
  limits, and picks different primers than the settings in the window describe**, so the tag is always written.

* **Platform workspaces: Signals, Alignment, Phylogeny, Trees** (`src/ui/platform/`, built when the person asked for the
  whole interface to be redesigned around raw sequence signals, multiple alignment, phylogenetic estimation, tree
  drawing and 3D protein structure; nothing here is in ApE). Home says what each workspace does, which programs were
  found (MAFFT, IQ-TREE, primer3, fastp, curl, Snakemake, Nextflow) and offers two routes ("From signal to tree", "From
  gene to fold"); every workspace starts from an example (`--view signals`, `alignment`, `phylogeny`, `trees`,
  `trees-circular`, `structure`, `pipelines`, `home`). Layout, colours and spacing come from `src/ui/design.rs` (the
  numbers are **estimated**, listed in `mechanics.json`).
  **Signals** (`src/seq/abif.rs`, `src/seq/scf.rs`, `chromatogram.rs`, `signals.rs`): the original has a whole Sanger
  trace window, found late (section 8, "Sanger trace window"), and this workspace follows its rules. ABIF (`.ab1`,
  `.abi`) and SCF (`.scf`, versions 2 and 3) files are read; the chromatogram is drawn as the original draws it (the
  traces fitted so the tallest peak fills the trace area, which starts 35 points under the top and ends 10 above the
  bottom as in the original (`new_abi_window` L31110–31111), scales across 0.5 to 5 and up 1 to 5, a quality graph behind
  them, a number above every tenth base, A green, C blue, G black, T red and magenta for the rest, reverse-complement
  mirroring the picture); File info gives the sample, the signal strengths and the run times; **Export as SCF** writes
  the whole trace as the original does. What is not in the original: **trimming** is Richard Mott's modified algorithm
  as phred and Biopython use it — every base scores `cutoff − P(error)` with `P = 10^(−Q/10)` and the stretch with the
  highest total is kept; the default cutoff is Q **13** (a 1 in 20 chance of being wrong, **estimated**) and the trimmed
  ends are shaded; bases are picked by dragging; the clean sequence goes to Sequences, to Alignment, to the clipboard
  or to a FASTA file. A second tab holds the FASTQ read summary and cleaning of Tools → Read Quality and Trimming (the
  same code, section above). **Checked** against the original's own reader and writer on real and made-up files
  (section 8) and against Biopython 1.88's `Bio.SeqIO.AbiIO` on the three real files (`cargo test real_files -- --ignored`
  with `HELIXLOOM_AB1_DIR` naming the files). The example trace is **synthetic** (computed from the demo plasmid, not
  recorded). Not built: the original's editing of a trace (change or add a base), linking a trace to a sequence window,
  searching a trace, printing and saving as EPS; mixed-base calling.
  **Alignment** (`alignment.rs`, `alnview.rs`, `src/seq/alnstats.rs`): sequences are gathered (the open ones, a FASTA or
  GenBank file, pasted text, a trimmed trace) and aligned by the built-in centre-star method or by MAFFT, on a worker
  thread. The viewer keeps the names and the header still while the letters scroll and paints only what is visible; it
  colours by nucleotide or amino-acid type, only the letters that differ from the consensus, or not at all, and shows a
  conservation bar and the consensus for every column (the commonest symbol; a column that is half gaps or more gives a
  gap; a gap never counts as conserved; a tie goes to the symbol that sorts first), the identity (mean conservation),
  the share of gaps, and the numbers of variable and parsimony-informative columns. DNA and protein are told apart by
  the letters. **Phylogeny** (`phylogeny.rs`; the same `estimate()` as Tools → Phylogenetic Tree): neighbour joining,
  UPGMA, maximum likelihood (JC69, K80, F81, HKY85 by BIC) or Bayesian sampling (HKY85, two runs), exactly as in the
  paragraphs above, on a worker thread; the tree becomes the one Trees shows, and the page hands on to PhyloCNN and Tree
  Tools. **Trees** (`trees.rs`, `treeview.rs`, `src/seq/treelayout.rs`): phylogram, cladogram and circular layouts are
  worked out once as lines and arcs, and the viewer and the **SVG export** draw the same pieces, so the saved picture is
  the picture shown; support values, a scale bar of a round length (1, 2 or 5 times a power of ten), tip names thinned
  out when rows are closer than a line of text, colours by the start of the tip names (up to eight groups), zoom 0.25 to
  60 times, a clicked clade lists its tips, midpoint root and ladderize (Tree Tools' code), Newick opened, pasted,
  copied and saved. Tested: small trees laid out and checked by hand in each layout, every branch and join drawn once,
  odd trees (one tip, no lengths) without failing, round scale-bar lengths, and an SVG that holds the strokes, the
  (escaped) names and a scale bar.
* **Structure workspace: protein 3D structure** (`src/seq/pdb.rs`, `secstruct.rs`, `fold.rs`, `src/pipeline/predict.rs`,
  `src/ui/platform/structure.rs`, `viewer3d.rs`). PDB files are read from their fixed columns (first model, alternate
  location blank or A, waters skipped, `HELIX` and `SHEET` records kept); a file with no `EXPDTA` record whose
  temperature factors all lie between 0 and 100 (or 0 and 1) is taken for a prediction and those factors for its
  confidence (pLDDT). **Secondary structure** is assigned from the backbone by the rules of DSSP (Kabsch and Sander,
  Biopolymers 22:2577, 1983), written from the published definition: a hydrogen bond from the C=O of one residue to the
  N-H of another when `0.084·332·(1/r(ON) + 1/r(CH) − 1/r(OH) − 1/r(CN))` is below −0.5 kcal/mol, the H placed 1.0 Å from
  the N opposite the previous residue's C=O, only between residues at least 3 apart and with alpha carbons within 9 Å,
  proline never a donor; turns of 3, 4 and 5 make helices (alpha before 3-10 before pi) and parallel or antiparallel
  bridges make strands. **Checked against PyDSSP 0.9.1** (three states) on six real structures, share of residues in
  agreement: 1UBQ 100 %, 3GB1 100, 1TEN 100, 2LZM 100, 1LYZ 95, 1CRN 91 (46 residues, so one boundary residue is 2 %)
  (`cargo test real_structures -- --ignored` with `HELIXLOOM_PDB_DIR`). **Prediction** sends the sequence (20 standard
  amino acids, 10 to 400 of them) to the **ESM Metagenomic Atlas** (`POST https://api.esmatlas.com/foldSequence/v1/pdb/`,
  ESMFold, Lin et al., Science 2023) with the `curl` program, so no network library is linked in (180 seconds allowed),
  and shows the PDB it answers with, coloured by pLDDT in the bands AlphaFold draws (above 90, 70 to 90, 50 to 70,
  below 50). **The sequence leaves the computer when Predict is pressed; the page says so beside the button.** Checked with
  one live request (hen lysozyme, 129 residues: mean pLDDT 95, drawn on the Structure page) and, for the
  client's error handling (an HTTP error, an answer without atoms, a missing or failing `curl`), against a stand-in
  server run locally. The viewer draws a ribbon (a Catmull-Rom curve through the alpha carbons, wider for helices and
  strands), the backbone, or every bond (found with a grid; up to 40,000 atoms), as depth-sorted lines that fade with
  distance, so no graphics-card code is needed; turn by dragging, zoom with the wheel, move with shift or the right
  button, click a residue (in the picture or in the sequence strip) to pick it. "Use the open sequence" translates the
  selection, or else the longest ORF, with the genetic code chosen in the editor. The example is ubiquitin (PDB 1UBQ,
  1.8 Å). Not built: OpenFold or an ESMFold inside the app (a workstation can fold on its GPUs with a program of its own,
  see the Workstation entry), side-chain and surface drawings, structure superposition, saving the picture.

* **Workstation: the heavy jobs on another computer** (`src/remote/`, `src/ui/platform/workstation.rs`,
  `tests/remote.rs`, `DEPLOY.md`, `deploy/`). Asked for by the person: a workstation with several GPUs does the heavy
  work, and the person's Mac, Windows or Linux computer is only the interface. Nothing here is in the original or in any
  program the person named; it is new work, written from the person's request and the standard library alone (no new
  crate). `helixloom serve` takes jobs over the network and computes them. The app's **Workstation** page connects (an
  address and a token, kept in a file that only the person can read), shows the workstation's GPUs, the programs it found
  and its jobs, and, when jobs are set to go there, Alignment, Phylogeny, Pipelines, PhyloCNN and
  Structure prediction send their work there, show what the workstation has (its MAFFT, its GPU) instead of this
  computer's, and show what comes back. A result waits on the workstation (7 days unless told otherwise), also after the
  app has been closed, until the person opens it from the page. **Where a job runs is the person's choice** (asked for
  in a later request: "either my workstation or my local computers, and if not, let the app choose how"): **This
  computer**, **Workstation** or **Automatic**, kept with the connection (`compute=` in the settings file; the earlier
  `use=1` is read as Workstation), set on the Workstation page and beside the start button of every page that starts a
  job. Every analysis asks `Platform::place(&Weight)`, where a `Weight` counts what makes its job slow (sequences,
  letters, bootstrap resamples that the method uses, Bayesian generations, trees, a GPU on the workstation); Automatic
  sends a job there when it passes any limit in `src/remote/weight.rs` (`HEAVY_*`: 100 sequences, 1,000,000 letters, 200 resamples,
  100,000 generations, 100 trees, or a GPU), which are **estimates**, and keeps it here otherwise. The page shows
  which way the job as set up would go and offers that computer's programs. With no workstation connected only This
  computer can be chosen. A new analysis (motif and domain finding, epitope and paratope prediction are planned) plugs into
  the same choice with a kind of job, an arm of `execute`, a `Weight` and a page (`DEPLOY.md`, section 10). The messages are plain HTTP/1.1 over TCP; a body is a
  container of named sections (`@name length`, the bytes, a line break) that is checked whole when it is read; numbers
  are written in the form that reads back as exactly the same number, so a result is the same on both ends. **Security:**
  a token of at least 32 characters (the one `helixloom token` makes comes from `/dev/urandom`), compared in constant time
  and read from the header before any body is accepted, never taken from the command line and never handed to a program
  the server runs; a job is one of six kinds, checked when it arrives, and nothing in a request is used as a path, a
  command or a program (the folding program comes only from the server's own `--fold-command`); limits on what a request
  may hold (5,000 sequences, 64 million letters, 20,000 trees, so that one job's memory stays within a few gigabytes), on
  the header (all of it within 5 seconds, however slowly it is sent), the upload (512 MB, and no more bytes of uploads in
  memory at once than two of the largest), the connections (64, and 8 from one address), the queue (1000 jobs), the jobs
  kept (2,000; a job's request is deleted when it is over) and the time of a job (24 hours); the server listens on
  127.0.0.1 unless told otherwise. The app is as careful with a server's replies: a size and a time for each kind, text
  cut to one short line, job numbers checked, and trees nested more than 2,000 levels deep refused (the readers go down a
  tree by recursion). **There is no encryption**, and a computer that can reach the port can make it busy without the
  token: `DEPLOY.md` says to use an SSH tunnel or a private network with a firewall, and never the open internet.
  **GPUs:** `--gpus auto` lists them with `nvidia-smi` (or takes the Apple GPU), or the person names them; each
  GPU runs one job at a time; a job is started as a process of its own with `CUDA_VISIBLE_DEVICES` set to the one GPU it
  was given; a job behind one that waits for a GPU may use the processor; a job that is stopped, or runs too long, is
  killed with everything it started. The jobs call the functions the app calls for itself (the tree estimators were
  moved from the tree window to `src/seq/estimate.rs` so that both use them). **Tried here, on a Mac with the server and
  the app on one computer:** unit tests of every part; `tests/remote.rs` starts the real program and speaks plain HTTP to
  it (jobs run in their own processes and come back; the token and the limits; results kept across a restart; a running
  job stopped and its process gone; three jobs sharing two pretend GPUs, with a stand-in program that reports the GPU it
  was shown); whole-window tests send an alignment, a tree, a pipeline, PhyloCNN estimates and a structure to a server on
  this computer and show the answers. **Not tried:** a real NVIDIA GPU or the `cuda` build, several real GPUs, the systemd
  unit and the ESMFold example in `deploy/`, the server on Windows, the app on Windows and Linux, an SSH tunnel or
  Tailscale, and uploads of hundreds of megabytes.

* **HelixLoom: the name, the interface, the models folder, the web build** (November 2026, asked for by the person; none of it
  is in ApE and none has an original address). The program is now called **HelixLoom** (it was ApECraft; its saved
  workstation settings are still read from the folder of the former name). The interface was redesigned and is no longer
  modelled on ApE's: a design system (`src/ui/design.rs`, `docs/DESIGN.md`), a top bar with the menus, a command palette on
  Cmd+K (`commands.rs`, `palette.rs`; Annotate moved to Cmd+Shift+K), toasts, a shortcuts sheet, a task-first Home and a mark
  drawn by code (`brand.rs`). The PhyloCNN networks and the Zurich HIV tree belong to their authors, whose repository has
  no licence: they are no longer built into the program but read from a models folder (`src/models.rs`,
  `docs/MODELS.md`), and the window's example is a made-up tree where that file is absent. The program builds for the browser
  (`docs/WEB.md`). What to add next, for the experience and for the pipelines, is in `docs/ROADMAP.md`. The code lives in the
  GitHub repository linhduongtuan/HelixLoom.

* **ASTRAL / ASTER → Tools → Species Tree from Gene Trees** (`src/seq/astral.rs` and `src/seq/astral/`, the job in
  `src/remote/species.rs`, the window `src/ui/windows/species.rs`, `helixloom astral`; November 2026, asked for by the person
  with "reimplement and rebuild ASTER and ASTRAL"; none of it is in ApE). Written from the papers and the programs' tutorials
  (Mirarab et al. 2014; Mirarab and Warnow 2015; Zhang et al. 2018; Sayyari and Mirarab 2016) and **not from the programs'
  code**: ASTRAL is Apache-2.0, ASTER is AGPL-3.0 and its source was not read. How it works (`docs/SPECIES.md` is the
  plain-words version):
  **Objective.** The species tree T (unrooted, binary) that maximizes S(T), the sum over the gene trees g of the number of
  four-species sets that g and T resolve the same way; a gene tree may lack species and may have polytomies, and a quartet it
  does not resolve counts for nothing. The normalized score is S(T) divided by the number of quartets the gene trees resolve.
  **Counting without listing quartets.** A resolved quartet ab|cd has two nodes that join a pair (the node where a and b
  meet with c and d together on one other branch, and the one where c and d meet); counting at those nodes counts each
  quartet twice. For a node of T with clusters A and B below it and the rest C above, `W(A,B)` is the sum over the nodes of
  the gene trees, for each branch k of the node, of (pairs of two species of one cluster on two other branches) times (the
  species of the other two clusters on k), over the three choices of the cluster that holds the pair; S(T) is half the sum of
  W over the inner nodes of T (`weight.rs`).
  **Search.** A dynamic programme over a set X of clusters (a cluster is a side of a split that does not hold a fixed root
  species, so a tree is read hanging from that leaf): best(V) is the highest best(A) + best(B) + W(A,B) over the ways to split V
  in two clusters A and B of X. X holds the single species and the whole; the splits of the gene trees, with the species a
  gene lacks put back (on the branch above the smallest group that holds their nearest relatives in a guide tree, which is
  neighbour joining on the mean number of branches between two species; its own splits are in X too); with `extra` 1 (the
  default, as ASTRAL's `-p 1`) every union of two of the branches of a polytomy of at most 8 branches (the children and the
  rest of the tree, since the tree is unrooted), also for the polytomies of the greedy consensus of the completed trees (the
  most frequent splits first, each kept when it fits those kept); with `exact`, all 2^(n−1) − 1 groups, for at most 12
  species. A polish then swaps neighbouring groups (nearest-neighbour interchanges, made in the tree itself) while S(T)
  rises, at most 30 sweeps and no more work than 2·10^10 word operations in all (a sweep costs nodes × 6 weights × the gene
  nodes × the words of a set), so fewer sweeps or none on very large inputs (`search.rs`, `complete.rs`). The sets, their order and every tie are decided without hash order, so the same gene
  trees give the same tree, branch order and numbers every time (a test runs it 24 times; it found a hash-set order that had
  made two runs list two branches the other way round). Putting the missing species back keeps the sets of species below
  every node up to date instead of working them out again for each species put back, and shares the genes among the cores (a
  test checks that the trees are the ones the plain way gives); a search with more than 1,000,000 groups to choose from is
  refused as soon as that is known, because the search looks at pairs of them.
  **What a branch carries** (`species.rs`, `special.rs`; the branch has L and R on one end, S and O on the other, O including
  the root species). For each gene tree the shares of its quartets around the branch that it resolves as LR|SO, RS|LO and
  RO|LS (each divided by the number it resolves, so a gene adds 1 to the three, or 0 when it resolves none): their sums over
  the genes are f1, f2, f3; EN = f1 + f2 + f3 (the effective number of genes); QC = |L||R||S||O| quartets; the quartet
  supports are q_i = f_i / EN. The local posterior probability of topology i (Sayyari and Mirarab 2016) is proportional to
  `3/2 · 2^−(EN − f_i) · B(f_i + 1, EN − f_i + 1) · (1 − I_{1/3}(f_i + 1, EN − f_i + 1))`, the three scaled to add up to 1
  (a uniform prior on the probability p ∈ [1/3, 1] that a gene has the branch's topology, B the beta function, I the
  regularized incomplete beta function, by a continued fraction in log space). The length in coalescent units is
  `−ln(1.5 · (1 − f1 / (EN + 2λ)))`, 0 when that is not positive, infinite when f1 = EN + 2λ, with **λ = 0.5** as ASTRAL's
  default (`-c`; λ = 0 is the maximum likelihood length). The Newick text follows ASTRAL's `-t` 0, 1, 2, 3 (default), 4 and 8
  and writes numbers as ASTRAL does (two decimals, or Java's shortest form for the full annotation); the root branch of a
  tree hung from a species is annotated once.
  **Reading.** Trees are unrooted; a tree with fewer than max(4, `-m`) leaves is left out and counted; a species twice in a tree
  is refused (gene families need ASTRAL-Pro, not built); at least four species must occur; a binary tree to score may be
  rooted or not and must have exactly the species of the gene trees; the Newick reader refuses trees nested more than 2,000
  levels deep, and the search, the counting and the Newick writer do not recurse, so a species tree as deep as it is wide
  (100,000 species in a ladder is tested) cannot run a thread out of stack.
  **Checked against ASTRAL 5.7.8**, run as a black box and never read: on 12 data sets of 6 to 80 species and 30 to 200 genes
  simulated under the multispecies coalescent (a simulator written for this project; some with polytomies and missing species),
  for 4 species trees each, the quartet score, the normalized score and, on every inner branch, f1, f2 and f3, EN, QC, q1, pp1
  and the length equal ASTRAL's (`golden/`; pp1 to 1e-6, because ASTRAL's own value of a posterior as small as 1.4e-5 is 8e-8 off
  the formula worked out to 40 digits, which this program matches to 1e-12); the best score of the default search is ASTRAL's
  or higher on all 12, and the exact search finds the optimum of ASTRAL's `-x` on the 6 sets small enough; on eight more
  generated sets (25, 40 and 60 species with 200, 400 and 600 genes, and 30, 50, 80, 15 and 24 species with 300, 200, 150, 500
  and 120 genes, the last five with missing species and, in four of them, polytomies) the best scores equal ASTRAL's to the last
  digit. **Checked against brute force**: the score equals counting
  every quartet, f1–f3 equal counting every quartet around the branch, and the exact search equals trying all rooted binary
  trees, for up to 7 species. Speed on those larger sets: 0.2, 1.8 and 7.9 s here, ASTRAL 1, 12 and 42 s.
  **Estimated, not read from anywhere** (`mechanics.json`): the polytomy limit 8 and the 1,500-species limit of the guide tree
  (above it the search runs without one), the 30 polish sweeps, the 12-species limit of the exact search, the 50,000 leaves above
  which Automatic sends a search to the workstation, the rows the table of the window lists (300), the 8 GB of sets of species
  that the search may keep (about three per tree node, a bit for every species; above it the trees are refused with what to
  do) and the 20,000,000 leaves a request may hold. The quartet counts are 64-bit, and trees on so many species and genes that
  they could overflow (more than about 18,000 species in 2,000 genes) are refused, not counted wrongly.
  **On the workstation:** a job kind `species`, run by the same function; a request holds at most 20,000 gene trees and
  20,000,000 leaves in all (counted by commas before any tree is read) and a result is read up to 8 MB; the compute choice
  (this computer, the workstation, automatic) decides where, as for every kind. The app does not trust a reply: a species name
  may be at most 800 bytes, the branch table at most as many rows as species, each row's species in order and once, and the
  table at most 1,000,000 entries and 64 MB of copied names; a table that is larger than that (every branch of a ladder-shaped
  tree lists the species below it) is cut by the server to the rows that fit, with the real number of branches kept in the
  result, and a reply that breaks the bounds is refused. Species names are written in the Newick in single quotes when they
  hold a space, a bracket, a comma, a colon, a semicolon or a quote, so that the text reads back as the same names; a branch
  length that is infinite or not a number (a prior of 0 makes them) is written as ASTRAL writes it and is left off when the
  tree is drawn on the Trees page.
  **Not built:** weighted ASTRAL (gene trees with lengths or supports as weights), ASTRAL-Pro3 (several copies per species),
  ASTRAL-IV's placement heuristic and CASTLES-II lengths in substitution units, CASTER-site and CASTER-pair (from alignments),
  WASTER (from reads), SISTER and the D* test, `-a` (several individuals to one species), bootstrapping over gene trees and
  the polytomy test (`-t 10`), the csv outputs.

### Not built

BEAST2's dated/coalescent analyses, MEGA's model selection for distances and its alignment editor, IQ-TREE's site-rate
models and partitions (IQ-TREE itself does them when called), mapping and variant tools (Minimap2 / SAMtools / GATK),
ETE / DendroPy beyond rooting, sorting, measuring and comparing trees, OpenFold and a locally run ESMFold inside the app
(the prediction asks the ESM Atlas service over the internet, or a folding program that the owner of a workstation names
with `--fold-command`; `deploy/esmfold_predict.py` is an example of one, not run here), NGLView's drawing styles beyond
ribbon, backbone and bonds, fastp as an external step, and a workflow editor (Pipelines runs one fixed chain). The
workstation is not used for the contig assembler, read cleaning, the primer3 window, chromatograms or tree drawing.

**Licences:** primer3, IQ-TREE and MrBayes are GPL, BEAST2 is LGPL. No code from them was copied or ported; primer3's
numbers, formulas and rules were read from its source and re-written, and its behaviour was checked by running it. ASTRAL
(Apache-2.0) was only run, to compare answers; ASTER (AGPL-3.0, which would reach anything copied from it) was not read, and
its programs other than the ASTRAL core are listed above as not built. Any engine added later is re-implemented the same way
and credited in `CREDITS.md`.

## Appendix: settings and inner lines cited in `mechanics.json` (line = byte address)

Menu bar: L1073 0xd8e2, L1074 0xd949, L1079 0xdb5e, L1081 0xdbea, L1082 0xdc57, L1146 0xf0f9, L1147 0xf165,
L1148 0xf1f2, L1150 0xf28a, L1151 0xf304, L1158 0xf5ab, L1161 0xf699, L1168 0xf7e0, L1173 0xf95e, L1179 0xfba0,
L1214 0x108cb, L1220 0x10a01, L1273 0x11d2f.
Inside procedures: L11085 0x838e4, L11921 0x8c286, L11954 0x8c7f5, L12119 0x8e133, L12120 0x8e2a0, L14732 0xa7324,
L14745 0xa768e, L14746 0xa769f, L14748 0xa76bb, L14749 0xa76cb, L14771 0xa7937, L14772 0xa7975, L14784 0xa7a9f,
L14928 0xa8ecb, L15124 0xaa838, L17090 0xbff3e, L17091 0xbff53, L17092 0xbff6b, L17093 0xbff80, L17191 0xc0ca7,
L17196 0xc0de4, L17327 0xc26d8, L17339 0xc280a, L20300 0xe611f, L42956 0x1f0bcb, L43160 0x1f2c31, L45309 0x210bb8,
L45576 0x2135d3, L48167 0x2314cb, L52652 0x262dbf.
Settings block (`set info(...)`): L66869 0x302c5c, L66881 0x302e56, L66907 0x3031da, L66910 0x30323f, L66920 0x3033be,
L66924 0x30342c, L66965 0x303a08, L66966 0x303a22, L66967 0x303a3a, L67003 0x303f98, L67004 0x303fb8, L67024 0x30421a,
L67026 0x304261, L67036 0x30483e, L67040 0x3048bd, L67041 0x3048d7, L67042 0x3048f0, L67043 0x304915, L67045 0x304943,
L67055 0x304aa3, L67101 0x304f2f, L67102 0x304f4c, L67103 0x304f6a, L67105 0x304fa3, L67106 0x304fc0, L67497 0x314c52.
