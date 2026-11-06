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
* **Command palette** (`palette.rs`, `commands.rs`). Type, move with the arrows, Enter to run, Esc to close. Every command has
  a group, words to find it by, and its key; the keys written there are checked against the keys the program reads.
* **Toasts** (`toasts.rs`). Short notes at the bottom right that fade by themselves. They are made from the status line, so
  every message the program already said is now seen; and one says when an analysis has finished.
* **Cards, chips, buttons** (`design.rs`). `card`, `chip`, `primary_button` (the one thing a page is for), `ghost_button`
  (everything else), `icon_tile`, `keycap`, `section`, `page_header`.
* **Windows** (`design::window`). Tool windows open inside the working area, under the bar and beside the rail, each shifted a
  little so that several open at once do not hide each other.

## Rules for a new page

1. Start with `page_header` (title, one line of purpose, the page's actions on the right).
2. One `primary_button` per page or card; other buttons are `ghost_button` or plain.
3. Put controls in a card on the left or the top, results on the right or below; say what to do in an empty state.
4. Say where a job runs (the compute picker) and report its end through the status line, which becomes a toast.
5. Add the page's actions to `commands.rs` so that the palette finds them, and a key if it deserves one.
6. Test it in both themes at 1400, 1100 and 780 points wide (`topbar.rs` shows how).

## Next

Phase 2: every workspace to this look, linked views and an inspector. Phase 3: docking panels. Phase 4: the font pair, motion
and an illustration set. See `ROADMAP.md`, part 4.
