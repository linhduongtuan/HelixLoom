use super::*;
use crate::ui::app::{Action, Startup};
use crate::ui::commands;
use crate::ui::platform::Workspace;
use crate::ui::testing::Screen;
use std::time::{Duration, Instant};

/// Draw frames until `done` says so (an analysis runs on a thread of its own).
fn until(screen: &mut Screen, app: &mut App, mut done: impl FnMut(&App) -> bool) {
    let end = Instant::now() + Duration::from_secs(120);
    loop {
        screen.frame(app, vec![]);
        if done(app) {
            return;
        }
        assert!(Instant::now() < end, "gave up waiting");
        std::thread::sleep(Duration::from_millis(10));
    }
}

/// The report that the command saves, as text.
fn saved_report(app: &mut App, screen: &Screen) -> String {
    let _ = figure::take_saved();
    app.run(Action::Report, screen.ctx());
    let mut saved = figure::take_saved();
    assert_eq!(saved.len(), 1, "one file");
    let (name, bytes) = saved.remove(0);
    assert_eq!(name, FILE);
    assert_eq!(app.status, format!("Saved {FILE}."));
    String::from_utf8(bytes).expect("a page is text")
}

#[test]
fn with_no_result_that_says_how_it_was_made_there_is_nothing_to_report_and_nothing_is_saved() {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    screen.settle(&mut app);
    let _ = figure::take_saved();
    app.run(Action::Report, screen.ctx());
    assert!(figure::take_saved().is_empty());
    assert_eq!(app.status, NOTHING);
    assert!(app.report_parts().is_empty());
    // the command is in the File menu's group of the palette, and needs no sequence in view
    let command = commands::all()
        .iter()
        .find(|c| matches!(c.action, Action::Report))
        .expect("a command");
    assert_eq!(
        (command.group, command.title),
        ("File", "Save a report of the results…")
    );
    assert!(!commands::acts_on_the_sequence(&Action::Report));
}

#[test]
fn an_alignment_is_reported_with_its_methods_its_figure_its_inputs_and_how_its_method_was_checked()
{
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    app.load_example(Workspace::Alignment);
    until(&mut screen, &mut app, |a| {
        a.platform.alignment.made.is_some() && !a.platform.alignment.is_running()
    });
    let page = saved_report(&mut app, &screen);
    let (made, _) = app.platform.alignment.made.clone().expect("a record");
    let rows = app
        .platform
        .alignment
        .result
        .as_ref()
        .expect("an alignment")
        .rows
        .clone();

    assert!(page.contains(": 1 result.</p>"), "{page}");
    assert!(page.contains("<h2>1. Alignment</h2>"));
    assert!(page.contains(&report::escaped(&report::methods(&made))));
    // the check of the built-in aligner, with what was not checked
    assert!(page.contains("<span class=\"grade\">second implementation</span>"));
    assert!(page.contains("never scored on a benchmark set"));
    // the figure is in the page: every row by its name, in the light colours whatever the window wears
    assert!(page.contains("<figure>\n<svg xmlns="));
    for (name, _) in &rows {
        assert!(page.contains(&format!(">{name}</text>")), "{name}");
    }
    let [r, g, b, _] = app.prefs.appearance.look(false).bg.to_srgba_unmultiplied();
    assert!(page.contains(&format!("fill=\"#{r:02x}{g:02x}{b:02x}\"/>\n")));
    assert!(page.contains("<figcaption>The "), "a line under the figure");
    // every input with its size and the whole of its checksum, and the programs at the end
    for input in &made.inputs {
        assert!(page.contains(&format!(
            "<td>{}</td><td>{}</td><td class=\"mono\">{}</td>",
            report::escaped(&input.name),
            input.size,
            input.digest
        )));
    }
    let listed = &page[page.find("<h2>Programs</h2>").expect("the programs")..];
    assert_eq!(
        listed.matches("<li>").count(),
        made.programs.len(),
        "{listed}"
    );
    for program in &made.programs {
        assert!(
            listed.contains(&format!("<li>{}</li>", report::escaped(program))),
            "{program} in {listed}"
        );
    }
    assert!(made
        .programs
        .contains(&format!("HelixLoom {}", env!("CARGO_PKG_VERSION"))));

    // the same report in the dark window: the figure stays light
    app.run(
        Action::Theme(crate::ui::state::ThemeMode::Dark),
        screen.ctx(),
    );
    screen.settle(&mut app);
    let dark = saved_report(&mut app, &screen);
    assert!(dark.contains(&format!("fill=\"#{r:02x}{g:02x}{b:02x}\"/>\n")));
}

#[test]
fn an_estimated_tree_is_reported_with_its_drawing_its_newick_text_and_the_check_of_its_method() {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    app.load_example(Workspace::Phylogeny);
    app.workspace = Workspace::Phylogeny;
    until(&mut screen, &mut app, |a| {
        a.platform.phylogeny.made.is_some() && !a.platform.phylogeny.is_running()
    });
    let page = saved_report(&mut app, &screen);
    let (made, _) = app.platform.phylogeny.made.clone().expect("a record");
    let tree = app.tree_win.tree.clone().expect("the tree");
    assert!(page.contains("Estimated tree</h2>"), "{page}");
    assert!(page.contains(&report::escaped(&report::methods(&made))));
    let check =
        checks::of_result(&made.analysis, &made.parameters).expect("the check of its method");
    assert!(page.contains(&report::escaped(check.against)));
    assert!(page.contains(&format!(
        "<h3>The tree in Newick format</h3>\n<table>\n<tr><th>Newick</th></tr>\n<tr><td>{}</td></tr>",
        report::escaped(tree.newick().trim())
    )));
    assert!(page.contains(&format!(
        "<figcaption>The tree of {} sequences",
        tree.leaf_count()
    )));
    // a tree that was changed since (built in the tree window, re-rooted) carries no record and is not reported as estimated
    app.platform.phylogeny.made = None;
    assert!(app
        .report_parts()
        .iter()
        .all(|part| part.title != "Estimated tree"));
}

/// The figure and the line under it of the alignment of three sequences of `length` bases, aligned on the page.
fn figure_of_three(length: usize) -> (App, String, String) {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    for k in 0..3 {
        let bases: Vec<u8> = (0..length).map(|i| b"ACGT"[(i + k) % 4]).collect();
        app.platform.alignment.add(&format!("row {k}"), &bases);
    }
    app.workspace = Workspace::Alignment;
    screen.settle(&mut app);
    screen.click(&mut app, "Align");
    until(&mut screen, &mut app, |a| {
        a.platform.alignment.made.is_some() && !a.platform.alignment.is_running()
    });
    let parts = app.report_parts();
    let (svg, caption) = parts[0].figure.clone().expect("a figure");
    (app, svg, caption)
}

#[test]
fn a_short_alignment_is_drawn_whole() {
    let (app, svg, caption) = figure_of_three(40);
    let columns = app
        .platform
        .alignment
        .result
        .as_ref()
        .expect("aligned")
        .rows[0]
        .1
        .len();
    assert!(columns < 120);
    assert_eq!(
        caption,
        format!("The alignment: 3 sequences, {columns} columns, with how conserved each column is above it.")
    );
    assert!(svg.contains(">row 2</text>"), "{svg}");
}

#[test]
fn a_long_alignment_is_drawn_by_its_first_columns_and_the_caption_says_so() {
    let (app, svg, caption) = figure_of_three(400);
    let columns = app
        .platform
        .alignment
        .result
        .as_ref()
        .expect("aligned")
        .rows[0]
        .1
        .len();
    assert!(columns >= 400);
    assert_eq!(
        caption,
        format!("The first 120 of {columns} columns and the first 3 of 3 sequences of the alignment; the whole alignment is in its FASTA file.")
    );
    // as wide as 120 columns and the names, not as the whole alignment
    let width: f32 = svg
        .split(" width=\"")
        .nth(1)
        .and_then(|rest| rest.split('"').next())
        .and_then(|w| w.parse().ok())
        .expect("the width");
    assert!(width > 1000.0 && width < 1300.0, "{width}");
}

#[test]
fn a_name_from_a_file_that_is_written_like_markup_adds_nothing_to_the_page_or_to_its_figure() {
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    let sly = "<script>alert(1)</script>&\"q\"";
    for (k, name) in [sly, "plain one", "plain two"].iter().enumerate() {
        let bases: Vec<u8> = (0..60).map(|i| b"ACGT"[(i + k) % 4]).collect();
        app.platform.alignment.add(name, &bases);
    }
    app.workspace = Workspace::Alignment;
    screen.settle(&mut app);
    screen.click(&mut app, "Align");
    until(&mut screen, &mut app, |a| {
        a.platform.alignment.made.is_some() && !a.platform.alignment.is_running()
    });
    let name = app
        .platform
        .alignment
        .result
        .as_ref()
        .expect("aligned")
        .rows[0]
        .0
        .clone();
    assert!(
        name.contains("<script>"),
        "the page keeps the name as it is in the file: {name}"
    );
    let page = saved_report(&mut app, &screen);
    assert!(
        !page.contains("<script>") && !page.contains("</script>"),
        "{page}"
    );
    // in the figure and in the table of inputs it is text
    let written = report::escaped(&name);
    assert!(page.matches(&written).count() >= 1, "{written} in {page}");
    assert!(page.contains("&lt;script&gt;alert(1)&lt;/script&gt;&amp;"));
    assert!(page.contains("</text>"));
}

#[test]
fn a_tree_of_more_tips_than_a_drawing_has_room_for_says_so_under_it() {
    assert_eq!(crowded(100), "");
    assert_eq!(crowded(295), "");
    // 296 tips take the 6,000 points a figure can have, and one more does not fit
    assert_eq!(crowded(296), "");
    assert_eq!(crowded(297), CROWDED);
    assert_eq!(tree_height(10), 280.0);
}

/// A record as a page makes it, of an analysis that ran here.
fn record(
    analysis: &'static str,
    (inputs, parameters): (Vec<crate::provenance::Input>, Vec<(String, String)>),
) -> Provenance {
    crate::provenance::Asked {
        analysis,
        inputs,
        parameters,
        ran_on: "this computer",
    }
    .made(vec!["a program 1.0".to_string()], 1.5)
}

/// An app whose windows and pages hold a result of every tool that keeps a record, each as its window would hold it. The
/// motifs and the fixed route are run; the others (which need programs or files that are not everywhere) are put in.
fn with_the_results_of_the_tools() -> (App, Screen) {
    use crate::pipeline::domains::Cutoff;
    use crate::remote::domains::{DomainsJob, DomainsResult};
    use crate::seq::hmmer::DomainHit;
    use crate::seq::phylocnn::{Priors, Task};
    use crate::ui::platform::structure;
    use crate::ui::windows::phylocnn;
    let mut app = App::new(Startup::default());
    let mut screen = Screen::new();
    app.motif_win.show_example();

    let job = DomainsJob {
        queries: vec![("kinase one".to_string(), b"MKV".repeat(40))],
        cutoff: Cutoff::Curated,
    };
    app.domains_win.made = Some((
        record("Protein domain search", domains::described(&job)),
        (true, 1e-3),
    ));
    app.domains_win.result = Some(std::sync::Arc::new(DomainsResult {
        names: vec!["kinase one".to_string()],
        lengths: vec![120],
        hits: vec![DomainHit {
            sequence: 0,
            family: "Pkinase".to_string(),
            accession: "PF00069.30".to_string(),
            family_length: 264,
            score: 180.5,
            i_evalue: 1e-50,
            c_evalue: 1e-52,
            hmm_from: 1,
            hmm_to: 100,
            ali_from: 10,
            ali_to: 110,
            env_from: 8,
            env_to: 112,
            accuracy: 0.95,
            description: "Protein kinase <domain>".to_string(),
        }],
        program: "HMMER 3.4".to_string(),
        database: "Pfam-A.hmm".to_string(),
        seconds: 1.5,
    }));

    app.platform.structure.load_example();
    app.platform.structure.made = Some(record(
        "Structure prediction",
        structure::described(b"MQIFVKTLTG", None),
    ));

    let kept = (0.4, Priors::Simulated, Task::Bd);
    let w = &mut app.phylocnn_win;
    w.header = vec!["", "estimate"];
    w.table = vec![vec!["R0".to_string(), "2.500".to_string()]];
    w.message = "3 tips. Read this as a rough guide.".to_string();
    w.made = Some((
        record(
            "PhyloCNN estimates",
            phylocnn::described("((A:1,B:1):2,C:3);", kept),
        ),
        kept,
    ));

    app.show_example_flow_run();
    let flow = &mut app.platform.pipelines.flow;
    let settings =
        super::super::pipelines::flow::settings_of(&flow.pipeline, "example-loci: 4 files");
    let files = vec![crate::provenance::Input::of("locus1", b"ACGT")];
    flow.made = Some((
        record("Pipeline file", (files, settings)),
        flow.pipeline.clone(),
    ));

    app.platform.pipelines.use_installed = false;
    app.load_example(Workspace::Pipelines);
    until(&mut screen, &mut app, |a| {
        a.platform.pipelines.made.is_some() && !a.platform.pipelines.is_running()
    });
    (app, screen)
}

#[test]
fn the_results_of_the_tools_are_each_a_part_of_the_report_with_their_tables_and_drawings() {
    let (mut app, screen) = with_the_results_of_the_tools();
    let page = saved_report(&mut app, &screen);
    assert!(page.contains(": 6 results.</p>"), "{page}");
    // each part, in its order, from its heading to the next
    let titles = [
        "Motifs",
        "Protein domains",
        "Predicted structure",
        "Epidemic estimates (PhyloCNN)",
        "Pipeline run",
        "Pipeline file run",
    ];
    let part = |k: usize| -> &str {
        let from = page
            .find(&format!("<h2>{}. {}</h2>", k + 1, titles[k]))
            .unwrap_or_else(|| panic!("no part {}", titles[k]));
        let rest = &page[from..];
        &rest[..rest.find("</section>").expect("the end of the part")]
    };
    let grade = |analysis: &str| {
        let check = checks::of_result(analysis, &[]).expect("a check");
        format!("<span class=\"grade\">{}</span>", check.standard.label())
    };

    let motifs = part(0);
    let run = app.motif_win.run.as_ref().expect("motifs");
    assert!(motifs.contains(&grade("Motif search")), "{motifs}");
    assert!(motifs.contains("<figure>\n<svg xmlns="));
    assert!(motifs.contains("<h3>The motifs</h3>"));
    let consensus = run.found[0].motif.consensus();
    assert!(
        motifs.contains(&format!("<td>{consensus}</td>")),
        "{motifs}"
    );
    assert!(
        motifs.contains(&format!(">Motif 1: {consensus}</text>")),
        "the logo is under its line"
    );
    assert!(motifs.contains("<td>width</td><td>10</td>"));

    let domains = part(1);
    assert!(domains.contains(&grade("Protein domain search")));
    assert!(domains.contains("<figure>\n<svg xmlns="));
    assert!(domains.contains("each of the 1 sequences"), "{domains}");
    assert!(
        domains.contains("two families can share one, and the table names every domain"),
        "{domains}"
    );
    assert!(domains.contains("<td>Pkinase</td><td>PF00069.30</td><td>10–110</td>"));
    // what a file says is text on the page, not markup
    assert!(domains.contains("<td>Protein kinase &lt;domain&gt;</td>"));
    assert!(!domains.contains("<domain>"));

    let structure = part(2);
    assert!(structure.contains(&grade("Structure prediction")));
    assert!(structure.contains("<span class=\"grade\">not checked</span>"));
    assert!(structure.contains("<figure>\n<svg xmlns="));
    assert!(structure.contains("<h3>The structure</h3>"));
    assert!(
        structure.contains("<tr><td>76</td><td>1</td>"),
        "{structure}"
    );

    let estimates = part(3);
    assert!(estimates.contains(&grade("PhyloCNN estimates")));
    assert!(estimates.contains("<tr><td>R0</td><td>2.500</td></tr>"));
    assert!(estimates.contains("<p class=\"dim\">3 tips. Read this as a rough guide.</p>"));
    assert!(!estimates.contains("<figure>"), "a table has no drawing");

    let route = part(4);
    assert!(route.contains(&report::escaped(STEPS_CHECKED)), "{route}");
    assert!(route.contains("<h3>What ran</h3>"));
    assert!(route.contains("<h3>The alignment</h3>"));
    assert!(route.contains("<figcaption>The tree of 6 tips;"), "{route}");
    assert!(route.contains("<h3>The tree in Newick format</h3>"));
    assert!(route.contains("<h3>Primer pairs, by "));
    assert!(
        route.contains("<td>primer settings</td>"),
        "an input of the run"
    );

    let file = part(5);
    assert!(file.contains(&grade("Pipeline file")));
    assert!(
        file.contains(&report::escaped(FILE_STEPS_CHECKED)),
        "{file}"
    );
    assert!(
        file.contains("<th>What happened</th><th>Note</th>"),
        "{file}"
    );
    assert!(file.contains("<h3>What each step did</h3>"));
    assert!(
        file.contains("<td>locus1</td><td>tree</td><td>ran in "),
        "{file}"
    );
    assert!(file.contains("<td>all files</td><td>species</td>"));
    assert!(
        file.contains("<figcaption>The species tree of 6 tips;"),
        "{file}"
    );

    // a figure that is cut to a part of its room names its cuts, and a name is looked up in the whole page: every cut that
    // a figure uses is one of its own, and no two figures of the page have a cut of one name
    let mut all_names: Vec<&str> = Vec::new();
    let figures: Vec<&str> = page
        .split("<figure>")
        .skip(1)
        .filter_map(|rest| rest.split("</figure>").next())
        .collect();
    assert_eq!(figures.len(), 5);
    for figure in &figures {
        let own: Vec<&str> = figure
            .split("<clipPath id=\"")
            .skip(1)
            .filter_map(|rest| rest.split('"').next())
            .collect();
        let used: Vec<&str> = figure
            .split("clip-path=\"url(#")
            .skip(1)
            .filter_map(|rest| rest.split(')').next())
            .collect();
        assert!(
            used.iter().all(|name| own.contains(name)),
            "{used:?} of {own:?}"
        );
        all_names.extend(own);
    }
    assert!(
        all_names.len() >= 2,
        "the logos and the domains are cut: {all_names:?}"
    );
    let mut once = all_names.clone();
    once.sort_unstable();
    once.dedup();
    assert_eq!(once.len(), all_names.len(), "{all_names:?}");

    // a result whose record is gone (another structure was opened, say) is no longer a part
    app.platform.structure.load_example();
    app.phylocnn_win.table.clear();
    assert_eq!(app.report_parts().len(), 4);
}

#[test]
fn a_tree_of_a_run_that_cannot_be_read_has_a_note_and_its_text_and_no_drawing() {
    let app = App::new(Startup::default());
    let look = app.prefs.appearance.look(false);
    let looks = Looks {
        dark: false,
        appearance: app.prefs.appearance,
    };
    let made = record("Pipeline file", (Vec::new(), Vec::new()));
    let mut section = part("Pipeline file run", &made);
    add_tree(&mut section, "((A,B),", "The species tree", &look, looks);
    assert!(section.figure.is_none());
    assert_eq!(section.notes, [NOT_DRAWN]);
    assert_eq!(section.tables, [newick_table("((A,B),")]);
    // a length that is infinite or not a number (a prior of 0 makes some) is left out, and the tree is drawn all the same
    let endless = "((A:Infinity,B:0.2):NaN,(C:0.1,D:0.3):0.05);";
    assert_eq!(tree_of(endless).map(|tree| tree.leaf_count()), Some(4));
    let mut section = part("Pipeline file run", &made);
    add_tree(&mut section, endless, "The species tree", &look, looks);
    let (drawn, _) = section.figure.clone().expect("a drawing");
    assert!(drawn.contains(">A</text>") && section.notes.is_empty());
    assert!(!drawn.contains("NaN") && !drawn.contains("inf"), "{drawn}");
    assert_eq!(
        section.tables,
        [newick_table(endless)],
        "the text is as it was given"
    );
    // a tree with lengths has a bar for their scale, and the line under it says so
    let with_lengths = "((A:0.1,B:0.2):0.05,(C:0.1,D:0.3):0.05);";
    let mut section = part("Pipeline file run", &made);
    add_tree(&mut section, with_lengths, "The species tree", &look, looks);
    let (svg, caption) = section.figure.expect("a drawing");
    assert!(svg.contains(">D</text>"));
    assert_eq!(
        caption,
        "The species tree of 4 tips; the bar is the scale of the branches."
    );
    assert!(section.notes.is_empty());
    let laid = lay_out(&tree_of(with_lengths).expect("a tree"), Layout::Phylogram);
    assert!(treeview::scale_length(&laid, &treeview::Options::default()).is_some());
    // a tree whose branches have no length to measure has no bar, and the line does not speak of one
    let flat = "((A:0,B:0):0,(C:0,D:0):0);";
    let laid = lay_out(&tree_of(flat).expect("a tree"), Layout::Phylogram);
    assert_eq!(
        treeview::scale_length(&laid, &treeview::Options::default()),
        None
    );
    let mut section = part("Pipeline file run", &made);
    add_tree(&mut section, flat, "The species tree", &look, looks);
    let (_, caption) = section.figure.expect("a drawing");
    assert_eq!(caption, "The species tree of 4 tips.");
}

#[test]
fn a_table_longer_than_a_report_shows_is_cut_and_says_so() {
    let mut notes = Vec::new();
    let rows: Vec<Vec<String>> = (0..MOST_TABLE_ROWS + 7)
        .map(|k| vec![k.to_string()])
        .collect();
    let kept = cut(rows.clone(), "domains", &mut notes);
    assert_eq!(kept.len(), MOST_TABLE_ROWS);
    assert_eq!(kept[MOST_TABLE_ROWS - 1], rows[MOST_TABLE_ROWS - 1]);
    assert_eq!(
        notes,
        [format!(
            "The table has the first {MOST_TABLE_ROWS} of {} domains.",
            MOST_TABLE_ROWS + 7
        )]
    );
    let all = cut(rows[..MOST_TABLE_ROWS].to_vec(), "domains", &mut notes);
    assert_eq!(
        (all.len(), notes.len()),
        (MOST_TABLE_ROWS, 1),
        "no new note"
    );
}
