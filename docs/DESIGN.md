# The HelixLoom interface: design system, version 1

HelixLoom's interface is no longer modelled on the program it started from. This document says what it looks like and why,
so that every new page can follow it. The code is in `src/ui/design.rs` (colours, shapes, shared widgets), `brand.rs` (the
mark), `topbar.rs`, `palette.rs`, `commands.rs`, `toasts.rs`, `sheet.rs` and `platform/` (the pages).

## Principles

1. **Task first.** People arrive with a question ("align these", "is this trace good?"), not with a menu in mind. Home asks
   "What will you work on today?", and the command palette answers a few letters of any task.
2. **One place to find anything.** Cmd+K searches every command, so no one has to remember which menu holds a thing.
3. **Calm, and honest about state.** Quiet surfaces; colour is kept for meaning (the accent for what to press, green, amber
   and red for how things went). What is running, where it runs and what finished is always visible.
4. **Same on every page.** One header, one set of spacings and shapes, the same buttons, so that what is learned on one page
   works on all.
5. **Light and dark are equals.** Both are designed, not inverted; the default follows the system.

## Colour

Two themes, the same roles. The values are in `look()`.

| Role | Light | Dark | Used for |
|------|-------|------|----------|
| `bg` | `#F5F6FA` | `#0B0E15` | the page behind cards |
| `surface` | `#FFFFFF` | `#131825` | cards, windows, the top bar |
| `surface_alt` | `#F0F2F8` | `#181E2D` | fields, quiet buttons |
| `border` | `#E1E5EF` | `#252E42` | outlines |
| `text` | `#131826` | `#E9EDF6` | main text |
| `dim` | `#5A6376` | `#9AA4B9` | secondary text |
| `accent` | `#5046E5` (loom indigo) | `#8F8BFF` | the colour of the program: the main button, the open workspace, links |
| `accent2` | `#0F9D94` (helix teal) | `#3DD1C5` | the second colour of the mark; what is new or alive |
| `good` / `warn` / `bad` | `#1B9A59` / `#B26A00` / `#D13D3D` | `#4CC38A` / `#E0A030` / `#F07178` | outcomes |

Text on a surface and on the accent keeps a luminance difference above 120 and 60 (a test checks it for both themes). The
colours of bases, residues and traces are data colours and stay as they are (`base_fill`, `residue_fill`, `trace_color`).

## Shape, space and depth

* Corners: cards and windows 12, buttons 8, plain widgets 4 (smaller, so a check box stays a square).
* Space: a page keeps an 18 point margin; a card 16; items 8 apart. Controls are 30 high.
* Depth: cards have a one-point shadow, windows and menus a soft larger one. Nothing else casts a shadow.
* Type: titles 22 to 32, section heads 15 to 19, body 13.5, helper text 12. (A bundled font pair with real weights is on the
  roadmap.)

## The pieces

* **Top bar** (`topbar.rs`). The mark and name (a click goes Home), the menus, the search box (Cmd+K), a pill that says where
  heavy jobs go (this computer, workstation, automatic) with a dot for the connection, a spinner with a count while analyses
  run, the shortcuts sheet and the theme button (system, light, dark). On a narrow window the words fall away and only the
  icons stay.
* **Rail.** One button per workspace, with a bar and a tint on the one that is open.
* **Command palette** (`palette.rs`, `commands.rs`, `answers.rs`). Type, move with the arrows, Enter to run, Esc to close. Every command has
  a group, words to find it by, and its key; the keys written there are checked against the keys the program reads. It also answers
  from the work in hand: a position (`go to 1200`, `jump 1200`, `base 1200` or the number alone) goes to that base; `find GAATTC`
  searches the sequence in view; a few letters (two or more) of the name of an open sequence, of a feature of the sequence in view
  (or its kind) or of a recent file list them to switch to, select or open. The answers that are asked for lead the list, the
  commands that match follow, and at most four names close it (three sequences, five features, three recent files at most of each).
* **Toasts** (`toasts.rs`). Short notes at the bottom right that fade by themselves. They are made from the status line, so
  every message the program already said is now seen; and one says when an analysis has finished.
* **Cards, chips, buttons** (`design.rs`). `card`, `chip`, `primary_button` (the one thing a page is for), `ghost_button`
  (everything else), `icon_tile`, `keycap`, `section`, `page_header`.
* **The frame of a page** (`design::page`, `page_header`, `controls_panel`, `results_area`). Every workspace that has
  controls and results (Signals, Alignment, Phylogeny, Trees, Structure, Pipelines, Workstation) is built in one frame, so
  that the eye finds the same thing at the same place whichever page is open: the header; under it, at the left, **one card
  of controls, 320 points wide on every page** (`CONTROLS_WIDTH`; the pages used to choose their own, from 270 to 340, and the
  results jumped sideways from page to page), which scrolls by itself when the window is short; and beside it the results,
  from the same place on. The rules a page keeps: its controls are one card, drawn by `controls_panel`, and the first thing
  in it is a heading (`section`); a page that has modes (the Sanger traces and the reads of Signals, the fixed route and the
  pipeline file of Pipelines) has their switch in the header, where its actions are, not above the controls; *Save figure*
  is among the header's actions, pale while there is nothing to save, so that the header does not move; and the button that
  runs the page's analysis (*Align*, *Estimate tree*, *Predict structure*, *Run pipeline*, *Connect*) is the filled button
  of the controls, paler while it cannot be pressed. An analysis that may run on the workstation ends its controls with
  `run_button` (`platform/workstation/compute.rs`): the line that says where it will run and, right under it, the button,
  with a wheel beside it while it runs, so the two are always found together (they were several controls apart on two
  pages, and the other way round on one). In the results, a picture that fills its frame (an alignment, a tree, a structure,
  a trace) is in `view_frame`, a card with a title is `titled_card`, and cards are `CARD_GAP` apart.
  **In a small window** (a tab of a browser on a small screen, a window beside another, or the appearance's larger sizes,
  which leave the window fewer points) the frame stays as it is, and three things give way: the line of purpose is cut
  short where the header's actions leave it no room (the pointer shows it whole), and the actions go on a row of their own
  under the title, which is never cut, when it has none beside them, so nothing of a header is drawn over anything else; a card that stands at
  the right of a view (what a tree is, on Trees) goes under the view when the view would be left less than 360 points; and a
  row of choices that wraps holds single widgets only, never `add_enabled_ui`, which does not wrap with the others (that hid
  the structure altogether in a small window). The rail squeezes its buttons and then scrolls.
  `ui/platform/same_look.rs` holds the frame for all seven pages, with and without their examples, in the usual window and
  in ones of 853 by 533 and 640 by 400 points. A column of controls
  is often longer than the window, so a test that clicks a button by its words (`Screen::click`) first turns the wheel over
  the column when the button is out of sight, as a person would: adding a control must not break the tests of the ones
  under it.
* **How this was made** (`ui/provenance.rs`). The section under a result that says what made it and what has changed since, with
  *Copy the record* and *Run again with these settings*; a page that makes a result keeps an `Asked` when a run starts and turns it
  into a `Provenance` when the answer comes (the Protein Domains window, `ui/windows/domains.rs`, is the shortest to copy). The rules
  every page keeps: a `described` function gives the inputs and the settings both for the record and for what they are now, so that
  the two cannot be written differently; the record is `None` whenever the result it is of is replaced, cleared or fetched from
  elsewhere; and *Run again* puts back the settings that are kept beside the record before it starts the run.
* **One selection on every page** (`ui/platform/linked.rs`). What is selected on one page is selected on every page that holds the
  same thing: today a column of the alignment (`AlignmentPage::select_column`, the line of `alignment/column.rs`) and the residue of the
  structure that was lined up with it (`Res::column`, set by *Colour by conservation*; `structure/link.rs`). The rules a new pair keeps:
  a page changes the selection only through `Platform` (`select_column`, `residue_picked`), never the other page's fields; the link
  exists only while what it rests on is current (the columns go with the colours when the alignment or its score changes), and a page
  that is not linked is left as it is; a selection that comes from another page is scrolled to once, when the page is next shown (the
  ORF list does the same for its rows: once, so that it can be scrolled away); what is selected is said in words as well as marked, since
  a colour is not read by everyone; and the line that says it is `design::selection_line`: one row whatever it says, with its links at
  the right, so that nothing moves under the pointer (the alignment's column and the Signals page's picked bases use it).
* **Figures** (`ui/figure.rs`). A page gives its view a figure with two lines: `figure::button(ui)` among its buttons (before the
  view is drawn), and where the view is drawn `if let Some(format) = figure::wanted() { message = figure::save(format, "name",
  size, |ui| { the_view::show(ui, &view); }) }`, with the data it is about to draw. The view's function is called once more in a
  window of its own that nobody sees, so the figure cannot differ from the screen; a view that scrolls says where it is scrolled
  to (`alnview::scrolled`), and one that shows a hint for the pointer leaves it out when `figure::is_figure`.
* **The report** (`report.rs` makes the page, `ui/platform/report.rs` gathers its parts). A result is in the report when it
  carries the record of how it was made (`Provenance`): to add one, give `report_parts` a `Section` with the record, the check
  that `checks::of_result` finds for it, the figure drawn by its view (`figure::draw`, in the light colours), and its tables.
  Everything that comes from a file goes through `report::escaped`; the figure's own text is escaped where it is made.
* **Programs HelixLoom calls** (`ui/windows/programs.rs`). A card for each program of the registry (`pipeline/tools.rs`), as
  wide as the list: its name with a chip (*found* in green, *not found* in amber), what it is for, where it is and its
  version or how to install it by hand, its licence. *Install…* shows the plan in a card of its own with the one primary
  button of the window (*Run it*), *Copy the command* and *Cancel*; nothing runs before that button, and while a plan waits
  or runs no other can be asked for. A program added to the registry appears here, in `helixloom tools` and in the
  diagnostics without more work.
* **Diagnostics** (`diagnostics.rs` makes the text, `ui/diagnostics.rs` gathers the facts). What goes in must be safe to
  paste where others read it: no sequence, no name of an opened file, no address, no token, no text that a fault or another
  program wrote freely (a fault gives its place in the code, a program its version only when it looks like one), and every
  path through `tilde`. A new line in the text is a new thing to ask this of.
* **How each analysis was checked** (`ui/checks.rs`, the table in `checks/table.rs`). The page of the Help menu: a chip for each
  grade (green for a reference program, the accent for a second implementation, grey for known answers, amber for not checked;
  `Tone::Warn` is new for the last), a filter, and a row for each analysis that a click opens. A new analysis gets a line in the
  table when it is built, with what was *not* checked written as plainly as what was.
* **Windows** (`design::window`). Tool windows open inside the working area, under the bar and beside the rail, each shifted a
  little so that several open at once do not hide each other. A window that was moved opens there next time (`places.rs`,
  `places.txt` in the settings folder); every tool window goes through `design::window`, which gives it its id and its place.

## Rules for a new page

1. Start with `page_header` (title, one line of purpose, the page's actions on the right).
2. One `primary_button` per page or card; other buttons are `ghost_button` or plain.
3. Put the controls in `controls_panel` (one card, a heading first, `run_button` for the analysis) and the results in
   `results_area` (`view_frame` round a picture, `titled_card` for a card); say what to do in an empty state. Modes and
   *Save figure* go in the header.
4. Say where a job runs (the compute picker) and report its end through the status line, which becomes a toast.
5. Add the page's actions to `commands.rs` so that the palette finds them, and a key if it deserves one.
6. Test it in both themes at 1400, 1100 and 780 points wide (`topbar.rs` shows how).

## Next

Phase 2 (every workspace to this look, linked views and an inspector) is built in part: the frame of a page, the empty
states, the inspector and two linked pairs are there (`ROADMAP.md`, part 4, says what is not). Phase 3: docking panels.
Phase 4: the font pair, motion and an illustration set.
