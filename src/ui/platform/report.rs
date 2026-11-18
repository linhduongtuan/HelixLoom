//! The report of the results that are on the pages (`crate::report`, `docs/ROADMAP.md` 1.7): File → Save a Report of the Results….
//! A result is in the report when it carries the record of how it was made: the alignment of the Alignment page, the tree that
//! the Phylogeny page estimated, the species tree, the motifs, the protein domains and the PhyloCNN estimates of their
//! windows, a predicted structure, and a run of the Pipelines page (its fixed route or a pipeline file). What has a drawing is
//! drawn for the report by the view that draws it in the window, in the light colours (a report is a page to read and to
//! print), whatever the window wears.

use super::pipelines::flow::FlowPage;
use super::pipelines::PipelinesPage;
use super::structure::StructurePage;
use super::{alnview, treeview};
use crate::checks;
use crate::pipeline::flow::{Status, Value};
use crate::provenance::Provenance;
use crate::report::{self, Section, Table};
use crate::seq::newick::NwTree;
use crate::seq::tree::Node;
use crate::seq::treelayout::{lay_out, Layout};
use crate::ui::app::App;
use crate::ui::design::Look;
use crate::ui::figure::{self, Looks};
use crate::ui::windows::domains::{self, DomainsWin};
use crate::ui::windows::motifs::{self, MotifWin};
use crate::ui::windows::phylocnn::PhyloCnnWin;
use eframe::egui::{vec2, Vec2};

pub const FILE: &str = "report.html";
pub const TITLE: &str = "HelixLoom results";
pub const NOTHING: &str = "There is nothing to report yet: a report is of the results that say how they were made, and no page or window has one. Run an analysis first.";
/// The most rows of a table in a report: what is left out is counted in a note, and is in the window to copy.
const MOST_TABLE_ROWS: usize = 200;
/// What the line under a tree says when a scale bar is drawn in it.
pub const BAR: &str = "; the bar is the scale of the branches";
pub const FILE_STEPS_CHECKED: &str = "The check above is of the running of a pipeline file: its steps, its cache and its export. The methods of the steps (the aligner, the tree search, the species tree) each have a check of their own on the page of checks (Help).";
pub const STEPS_CHECKED: &str = "A run of the fixed route is several methods, each with a check of its own: the page of checks (Help) has the aligner, the tree search and the primer design that ran.";
/// The room a structure is drawn in for a report.
const STRUCTURE_ROOM: Vec2 = vec2(720.0, 540.0);
/// The most rows and columns of an alignment that a report draws: a page is for reading, and the whole alignment is in its file.
const MOST_ROWS: usize = 60;
const MOST_COLUMNS: usize = 120;
/// The width of a tree in a report, the room of one tip, and the heights a tree can have.
const TREE_WIDTH: f32 = 900.0;
const TIP_ROOM: f32 = 20.0;
const TREE_HEIGHTS: (f32, f32) = (240.0, 6000.0);
pub const CROWDED: &str = " The tips are closer together than their names are tall: for a tree of so many, save a figure of a part of it from the Trees page.";
pub const NOT_DRAWN: &str =
    "The tree could not be drawn from its Newick text, which is in the table below.";

/// The height a tree of so many tips wants in a report.
fn tree_height(tips: usize) -> f32 {
    tips as f32 * TIP_ROOM + 80.0
}

/// What the line under a tree adds when the tree has more tips than a drawing has room for.
fn crowded(tips: usize) -> &'static str {
    if tree_height(tips) > TREE_HEIGHTS.1 {
        CROWDED
    } else {
        ""
    }
}

/// A tree drawn for a report: a row for each tip, the names beside them. With it, what the line under it says of the bar:
/// that it is the scale of the branches, when one is drawn (a tree whose branches have no length has none).
fn tree_figure(look: &Look, looks: Looks, tree: &Node) -> (String, &'static str) {
    let laid = lay_out(tree, Layout::Phylogram);
    let options = treeview::Options::default();
    let height = tree_height(tree.leaf_count()).clamp(TREE_HEIGHTS.0, TREE_HEIGHTS.1);
    let svg = figure::draw(looks, vec2(TREE_WIDTH, height), 1.0, |ui| {
        treeview::show(
            ui,
            look,
            &laid,
            &options,
            &mut treeview::Camera::default(),
            &mut None,
            false,
        );
    })
    .svg();
    let bar = match treeview::scale_length(&laid, &options) {
        Some(_) => BAR,
        None => "",
    };
    (svg, bar)
}

/// A tree as a table of one cell: its Newick text, to copy from the page.
fn newick_table(newick: &str) -> Table {
    Table {
        caption: "The tree in Newick format".to_string(),
        header: vec!["Newick".to_string()],
        rows: vec![vec![newick.trim().to_string()]],
    }
}

/// The first tree of a Newick text, as the views draw it. A length that is infinite or not a number (a prior of 0 makes
/// some) is left out: the reader does not take it.
fn tree_of(newick: &str) -> Option<Node> {
    let drawable = newick.replace(":Infinity", "").replace(":NaN", "");
    NwTree::parse_all(&drawable)
        .ok()
        .and_then(|mut trees| (!trees.is_empty()).then(|| trees.remove(0).to_node()))
}

/// The rows of a table, cut to what a report shows, with a note that says so when they were cut.
fn cut(rows: Vec<Vec<String>>, what: &str, notes: &mut Vec<String>) -> Vec<Vec<String>> {
    if rows.len() > MOST_TABLE_ROWS {
        notes.push(format!(
            "The table has the first {MOST_TABLE_ROWS} of {} {what}.",
            rows.len()
        ));
    }
    rows.into_iter().take(MOST_TABLE_ROWS).collect()
}

fn words(cells: &[&str]) -> Vec<String> {
    cells.iter().map(|c| c.to_string()).collect()
}

/// The motifs of the Motif Finder: their logos, and a row for each.
fn motif_part(w: &MotifWin, look: &Look, looks: Looks) -> Option<Section> {
    let (made, run) = (w.made.as_ref()?, w.run.as_ref()?);
    let mut section = part("Motifs", made);
    let svg = figure::draw(looks, motifs::logos_room(&run.found), 1.0, |ui| {
        motifs::draw_logos(ui, &run.found, look);
    })
    .svg();
    section.figure = Some((
        svg,
        format!(
            "The logo of each motif found in {} sequences: in a column a letter is as tall as its share of what the column tells, two bits at most.",
            run.seqs.len()
        ),
    ));
    section.tables.push(Table {
        caption: "The motifs".to_string(),
        header: words(&[
            "Motif",
            "Likeliest letters",
            "Width",
            "Sites",
            "Information (bits)",
            "Gain (bits)",
        ]),
        rows: run
            .found
            .iter()
            .enumerate()
            .map(|(k, f)| {
                vec![
                    (k + 1).to_string(),
                    f.motif.consensus(),
                    f.motif.probs.len().to_string(),
                    f.sites.len().to_string(),
                    format!("{:.1}", f.motif.information()),
                    format!("{:.1}", f.gain),
                ]
            })
            .collect(),
    });
    Some(section)
}

/// The domains of the Protein Domains window: the sequences with their domains along them, and a row for each domain.
fn domains_part(w: &DomainsWin, look: &Look, looks: Looks) -> Option<Section> {
    let ((made, _), result) = (w.made.as_ref()?, w.result.as_ref()?);
    let mut section = part("Protein domains", made);
    let svg = figure::draw(looks, domains::tracks_room(result), 1.0, |ui| {
        domains::draw_tracks(ui, look, result);
    })
    .svg();
    let (drawn, all) = (domains::drawn(result), result.names.len());
    // there are eight colours: with more families two share one, and a box too narrow for its name is told by the table
    let colours = "A family has the same colour wherever it is found; there are eight colours, so two families can share one, and the table names every domain.";
    let caption = if drawn == all {
        format!("The domains along each of the {all} sequences. {colours}")
    } else {
        format!("The domains along the first {drawn} of {all} sequences. {colours}")
    };
    section.figure = Some((svg, caption));
    let rows = result
        .hits
        .iter()
        .map(|hit| {
            vec![
                result.names.get(hit.sequence).cloned().unwrap_or_default(),
                hit.family.clone(),
                hit.accession.clone(),
                format!("{}–{}", hit.ali_from, hit.ali_to),
                format!("{:.1}", hit.score),
                format!("{:.1e}", hit.i_evalue),
                hit.description.clone(),
            ]
        })
        .collect();
    section.tables.push(Table {
        caption: "The domains".to_string(),
        header: words(&[
            "Sequence",
            "Domain",
            "Accession",
            "From–To",
            "Score",
            "E-value",
            "What it is",
        ]),
        rows: cut(rows, "domains", &mut section.notes),
    });
    Some(section)
}

/// A predicted structure: drawn as the Structure page shows it, with what it is made of.
fn structure_part(page: &StructurePage, look: &Look, looks: Looks) -> Option<Section> {
    let (made, loaded) = (page.made.as_ref()?, page.loaded.as_ref()?);
    let mut section = part("Predicted structure", made);
    let svg = figure::draw(looks, STRUCTURE_ROOM, 1.0, |ui| page.draw_still(ui, look)).svg();
    section.figure = Some((
        svg,
        "The predicted structure, turned and coloured as it is on the Structure page.".to_string(),
    ));
    let s = loaded.fold.summary();
    section.tables.push(Table {
        caption: "The structure".to_string(),
        header: words(&[
            "Residues",
            "Chains",
            "Atoms",
            "Helix (%)",
            "Strand (%)",
            "Mean confidence (pLDDT)",
        ]),
        rows: vec![vec![
            s.residues.to_string(),
            s.chains.to_string(),
            s.atoms.to_string(),
            format!("{:.0}", s.helix),
            format!("{:.0}", s.strand),
            s.confidence
                .map_or("not given".to_string(), |c| format!("{c:.0}")),
        ]],
    });
    Some(section)
}

/// The estimates of the PhyloCNN window, as its table has them, with what the window says of them.
fn phylocnn_part(w: &PhyloCnnWin) -> Option<Section> {
    let (made, _) = w.made.as_ref()?;
    if w.table.is_empty() {
        return None;
    }
    let mut section = part("Epidemic estimates (PhyloCNN)", made);
    if !w.message.is_empty() {
        section.notes.push(w.message.clone());
    }
    section.tables.push(Table {
        caption: "The estimates; times are in the units of the branch lengths".to_string(),
        header: words(&w.header),
        rows: w.table.clone(),
    });
    Some(section)
}

/// A tree of a run for a part: its drawing, or the note that says why there is none, and its Newick text.
fn add_tree(section: &mut Section, newick: &str, what: &str, look: &Look, looks: Looks) {
    match tree_of(newick) {
        Some(tree) => {
            let (svg, bar) = tree_figure(look, looks, &tree);
            section.figure = Some((
                svg,
                format!(
                    "{what} of {} tips{bar}.{}",
                    tree.leaf_count(),
                    crowded(tree.leaf_count())
                ),
            ));
        }
        None => section.notes.push(NOT_DRAWN.to_string()),
    }
    section.tables.push(newick_table(newick));
}

/// A run of the Pipelines page's fixed route: what ran, and what each step gave.
fn pipeline_part(page: &PipelinesPage, look: &Look, looks: Looks) -> Option<Section> {
    let ((made, _), done) = (page.made.as_ref()?, page.outcome.as_ref()?);
    let mut section = part("Pipeline run", made);
    section.notes.push(STEPS_CHECKED.to_string());
    section.tables.push(Table {
        caption: "What ran".to_string(),
        header: words(&["Step"]),
        rows: done.log.iter().map(|line| vec![line.clone()]).collect(),
    });
    if let Some(a) = &done.aligned {
        section.tables.push(Table {
            caption: "The alignment".to_string(),
            header: words(&["Sequences", "Columns", "Made by"]),
            rows: vec![vec![
                a.rows.len().to_string(),
                a.columns().to_string(),
                a.engine.clone(),
            ]],
        });
    }
    if let Some(t) = &done.tree {
        section.tables.push(Table {
            caption: "The tree".to_string(),
            header: words(&["Made by", "Model", "ln L"]),
            rows: vec![vec![
                t.engine.clone(),
                t.model.clone(),
                t.ln_likelihood.map_or(String::new(), |l| format!("{l:.2}")),
            ]],
        });
        add_tree(&mut section, &t.newick, "The tree", look, looks);
    }
    if let Some(p) = &done.primers {
        if p.pairs.is_empty() {
            section
                .notes
                .push("No primer pair fits the primer settings.".to_string());
        }
        section.tables.push(Table {
            caption: format!("Primer pairs, by {}", p.engine),
            header: words(&[
                "#",
                "Left primer",
                "Tm",
                "Right primer",
                "Tm",
                "Product",
                "Penalty",
            ]),
            rows: p
                .pairs
                .iter()
                .enumerate()
                .map(|(k, pair)| {
                    vec![
                        (k + 1).to_string(),
                        pair.left.clone(),
                        format!("{:.1}", pair.left_tm),
                        pair.right.clone(),
                        format!("{:.1}", pair.right_tm),
                        pair.product_size.to_string(),
                        format!("{:.3}", pair.penalty),
                    ]
                })
                .collect(),
        });
    }
    Some(section)
}

/// A run of a pipeline file: what each step did on each file, and the tree it gave.
fn flow_part(page: &FlowPage, look: &Look, looks: Looks) -> Option<Section> {
    let ((made, _), report) = (page.made.as_ref()?, page.report.as_ref()?);
    let mut section = part("Pipeline file run", made);
    section.notes.push(FILE_STEPS_CHECKED.to_string());
    if !report.left_out.is_empty() {
        section.notes.push(format!(
            "Left out of what was gathered: {}.",
            report.left_out.join(", ")
        ));
    }
    section.notes.extend(report.warnings.iter().cloned());
    let rows = report
        .items
        .iter()
        .flat_map(|item| item.steps.iter().map(|step| (item.name.as_str(), step)))
        .chain(report.gathered.iter().map(|step| ("all files", step)))
        .map(|(file, step)| {
            let happened = match &step.status {
                Status::Ran => format!("ran in {:.2} s", step.seconds),
                Status::Cached => "from the cache".to_string(),
                Status::Failed(why) => format!("failed: {why}"),
            };
            vec![
                file.to_string(),
                step.step.name().to_string(),
                happened,
                step.note.clone(),
            ]
        })
        .collect();
    section.tables.push(Table {
        caption: "What each step did".to_string(),
        // a note is what made the answer (an aligner, a tree builder with its model) or what the step found
        header: words(&["File", "Step", "What happened", "Note"]),
        rows: cut(rows, "steps", &mut section.notes),
    });
    if let Some(Value::Tree(newick)) = &report.output {
        add_tree(&mut section, newick, "The species tree", look, looks);
    }
    Some(section)
}

fn part(title: &str, made: &Provenance) -> Section {
    Section {
        title: title.to_string(),
        made: made.clone(),
        checked: checks::of_result(&made.analysis, &made.parameters),
        notes: Vec::new(),
        figure: None,
        tables: Vec::new(),
    }
}

impl App {
    /// The results on the pages that carry a record of how they were made, each as a part of a report.
    pub(in crate::ui) fn report_parts(&self) -> Vec<Section> {
        let looks = Looks {
            dark: false,
            appearance: self.prefs.appearance,
        };
        let look = self.prefs.appearance.look(false);
        let mut parts = Vec::new();

        let page = &self.platform.alignment;
        if let (Some((made, _)), Some(view)) = (&page.made, page.view(&look)) {
            let (rows, columns) = (view.rows.len(), view.columns.len());
            let (shown_rows, shown_columns) = (rows.min(MOST_ROWS), columns.min(MOST_COLUMNS));
            let room = alnview::room(shown_rows, shown_columns);
            let svg = figure::draw(looks, room, 1.0, |ui| {
                alnview::show(ui, &view);
            })
            .svg();
            let caption = if (shown_rows, shown_columns) == (rows, columns) {
                format!("The alignment: {rows} sequences, {columns} columns, with how conserved each column is above it.")
            } else {
                format!(
                    "The first {shown_columns} of {columns} columns and the first {shown_rows} of {rows} sequences of the alignment; the whole alignment is in its FASTA file."
                )
            };
            let mut section = part("Alignment", made);
            section.figure = Some((svg, caption));
            parts.push(section);
        }

        if let (Some((made, _)), Some(tree)) = (&self.platform.phylogeny.made, &self.tree_win.tree)
        {
            let mut section = part("Estimated tree", made);
            let (svg, bar) = tree_figure(&look, looks, tree);
            section.figure = Some((
                svg,
                format!(
                    "The tree of {} sequences{bar}.{}",
                    tree.leaf_count(),
                    crowded(tree.leaf_count())
                ),
            ));
            section.tables.push(newick_table(&tree.newick()));
            parts.push(section);
        }

        let window = &self.species_win;
        if let (Some((made, _)), Some(result)) = (&window.made, &window.result) {
            let mut section = part("Species tree", made);
            match tree_of(&result.newick) {
                Some(tree) => {
                    section.figure = Some((
                        tree_figure(&look, looks, &tree).0,
                        format!(
                            "The species tree from {} gene trees.{}",
                            result.genes,
                            crowded(tree.leaf_count())
                        ),
                    ));
                }
                // a part with no figure says why
                None => section.notes.push(NOT_DRAWN.to_string()),
            }
            section.tables.push(Table {
                caption: "Scores".to_string(),
                header: vec![
                    "Gene trees".into(),
                    "Quartet score".into(),
                    "Share of the quartets".into(),
                ],
                rows: vec![vec![
                    result.genes.to_string(),
                    result.score.to_string(),
                    format!("{:.4}", result.normalized),
                ]],
            });
            section.tables.push(newick_table(&result.newick));
            parts.push(section);
        }
        let pipelines = &self.platform.pipelines;
        parts.extend(
            [
                motif_part(&self.motif_win, &look, looks),
                domains_part(&self.domains_win, &look, looks),
                structure_part(&self.platform.structure, &look, looks),
                phylocnn_part(&self.phylocnn_win),
                pipeline_part(pipelines, &look, looks),
                flow_part(&pipelines.flow, &look, looks),
            ]
            .into_iter()
            .flatten(),
        );
        parts
    }

    /// Save the report of the results on the pages, or say that there is none to report.
    pub(in crate::ui) fn save_report(&mut self) {
        let parts = self.report_parts();
        if parts.is_empty() {
            self.say(NOTHING);
            return;
        }
        let page = report::html(TITLE, &crate::provenance::today(), &parts);
        let said = figure::keep(FILE, page.as_bytes());
        if !said.is_empty() {
            self.say(said);
        }
    }
}

#[cfg(test)]
mod tests;
