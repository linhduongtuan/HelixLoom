# Trees: drawing, annotations and editing, after FigTree and TreeGraph 2

Asked on 20 November 2026: to improve the tree drawing with what FigTree and TreeGraph 2 do. Both are published under the GNU General
Public Licence, so **neither one's source code was read**: the features were taken from what their documentation, tutorials, change
logs and the TreeGraph 2 paper (Stöver and Müller, *BMC Bioinformatics* 2010) say a person sees and does, and written here from those
descriptions. Where the documentation does not give a rule, the rule is this program's own and is said so below.

The **Trees** page does all of it; the example tree and any tree the Phylogeny page estimates work as before.

## Opening trees, with their annotations

*Open Newick or NEXUS…* (or paste) reads Newick and NEXUS, with the `translate` table of a NEXUS file. What a program writes in square
brackets after a node, `[&posterior=0.98,height_95%_HPD={0.5,1.5},rate=1.2]` (the way BEAST's TreeAnnotator, MrBayes and FigTree write
it), is kept: a number, a range in braces, a word in quotes. `[&R]`, `[&U]` and NHX comments are no annotation. A tree from BEAST shows its
`posterior` beside the inner nodes at once.

*Save NEXUS…* writes the tree with those annotations again; *Save Newick…* writes it plain.

An annotation belongs to the branch above the node. When the root is moved the branch carries what is written on it.

Annotations are kept for the session, not in the results the program keeps between sessions (those keep the tree itself).

## Layouts

Phylogram, cladogram and circular as before, and **Radial**: the unrooted layout of Felsenstein's equal-angle method. Each branch leaves its
parent in the middle of a share of the circle as large as the share of the tips below it; branch lengths are to scale and the root is in the
middle. (FigTree has a layout of that name; its drawing was not compared.)

## What the annotations look like

*Annotations* on the left:

* **Colour branches by** any annotation. A number is coloured along a gradient from its least to its greatest value (violet, green, yellow),
  with a scale in the corner; anything else gets a colour for each different value, with a list of them.
* **Widen branches by** a number: from one to four times as wide.
* **Node shapes by** an annotation: a disc on every node that has it, larger for a greater number.
* **Bars across nodes**: an annotation that is a range of heights (a 95 % HPD interval) is drawn as a bar across the node on a phylogram. The
  range is taken as a height *before the tips*, which is how BEAST writes it.
* **Beside the inner nodes**: any number of annotations, and the node's own label, joined with slashes (`0.98/92`). Whole numbers are written
  whole, others to two decimals.
* **Colour tip names like their branches**.
* Hovering over a node lists everything written on it.

## An axis

*Axis*, for the phylogram: ticks at round numbers under the tree, as *distance from the root*, as *time before the tips* (0 at the tips), or as
*dates* (the tips at a date you give, the units per branch length given). Lines can run up from the ticks. This is FigTree's time scale and
scale axis in the three forms the documentation describes; the details (defaults, how it rounds) are this program's.

## Clades

Click a node, then:

* **Collapse to a triangle** (FigTree's cartoon): the clade becomes one triangle, as tall as its number of tips needs (more tips, taller, up
  to a limit), labelled with its tips; the card still lists every tip's name. *Expand* puts it back.
* **Highlight** with one of eight colours: a block of colour behind the clade on a rectangular tree, the clade's lines in that colour on a
  round one.
* **Rotate** (turn the branches of the node round), **Root here**, **Show only this clade**.

*Find tips* marks the tips whose names hold some text.

## Editing

* **Root here** puts the root on the picked node's branch, in the middle of it. The branch's length is shared between the two new branches,
  both keep what was written on the branch, every branch between the old root and the node is turned round with what is written on it, and an
  old root that is left with one branch is dissolved. (FigTree and TreeGraph 2 do not say what happens to the length; this is the usual
  choice.) *Root at midpoint* keeps node labels but not other annotations.
* **Order: smaller first / larger first** (FigTree's order nodes, TreeGraph 2's ladderize): at every node the branches are put in order of the
  number of tips they lead to. *Ladderize* is the same as the first.
* **Collapse** branches whose value is below a number (TreeGraph 2's collapse nodes by support): any annotation, or the node's label; the
  branches of a taken-out node go to its parent and its length is added to theirs. The root and the tips stay.
* **Rename tips** from a table pasted as lines of *old name, a tab or comma, new name*.

## Support from other trees

*Support from other trees* (TreeGraph 2's add support values): open a file of trees (Newick or NEXUS), or paste them, and give the value a
name.

* **Share of the trees that have the branch** (bootstrap replicates, posterior samples): each branch of the tree on the page gets the percent
  of the other trees that have it as `name`, and the percent of the trees that have the branch most often in conflict with it as
  `name_conflict`. A branch is the same in two trees when it separates the same tips; where the root is does not matter, and a branch
  counts once in a tree.
* **The value on the same branch of one tree**: the number written on the matching branch of another analysis' tree (its label, or an
  annotation you name).

Only tips that every tree has are compared (at least four); the other tips are left out. A conflict is looked for among the 5,000 branches
most of the other trees have. The value is added as an annotation and shown beside the nodes; several values (`bootstrap/posterior`) can be
shown together. A branch the other trees do not have gets 0 for the share, and no value for the written one.

## How it was checked

`seq::treeedit`: after rooting on the branch above every kind of node the set of branches, their total lengths and their values are those of
the tree before; shares and conflicts of replicates are counted on paper; `seq::treelayout`: the lengths and angles of a radial tree, the rows
and numbers after a collapsed clade; `seq::treestyle`: the gradient, the colours, the ticks; `ui::platform::trees`: the page drawn without a
window and read as text and shapes (every layout draws every tip; a collapsed clade is a triangle that says its tips; the dates are on the
axis). They are listed on **Help ▸ Checks**. Nothing was compared with FigTree or TreeGraph 2 themselves, and the page has not been looked at on
a screen by a person: ask for that before relying on how it looks.

## Not done (from what the two programs have)

From FigTree: the fish-eye zoom; the polar layout's root angle, angle range and root length; per-taxon colours written in a NEXUS taxon block and
FigTree's own `figtree` block; importing a tab-delimited file of annotations for tips and inferring the ancestors' values from them; legends
for node shapes; node bars of other kinds than a height range; hiding a branch; selection modes (node, clade, taxa) and a filter by label; the
printing dialog. From TreeGraph 2: several values from different analyses as separate text labels with their own place, format and symbols
(here they are joined with slashes), pie charts, brackets and braces with text, the document's own file format (XTG), PhyloXML and NeXML,
editing a tree by cutting, copying and pasting nodes, sorting the tips by a list, importing a table as node data, calculating a column, a
synchronised selection between open trees, and a command line for the batch. *Save figure* writes what is on the screen as SVG, PDF or PNG
(with the annotations); the older *Save SVG…* writes the plain tree.
