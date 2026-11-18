use super::*;
use crate::provenance::Input;

fn made(
    analysis: &str,
    programs: &[&str],
    parameters: &[(&str, &str)],
    inputs: usize,
) -> Provenance {
    Provenance {
        analysis: analysis.to_string(),
        inputs: (0..inputs)
            .map(|k| Input {
                name: format!("seq <{k}>"),
                digest: format!("abc{k}"),
                size: 100 + k,
            })
            .collect(),
        parameters: parameters
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect(),
        programs: programs.iter().map(|p| p.to_string()).collect(),
        ran_on: "this computer".to_string(),
        date: "2026-11-13".to_string(),
        seconds: 1.234,
    }
}

fn part(title: &str, made: Provenance) -> Section {
    Section {
        title: title.to_string(),
        made,
        checked: None,
        notes: Vec::new(),
        figure: None,
        tables: Vec::new(),
    }
}

#[test]
fn the_methods_say_what_made_a_result_when_where_with_which_settings_and_from_how_many_inputs() {
    let record = made(
        "Alignment",
        &["MAFFT v7.526", "HelixLoom 0.2.0"],
        &[("engine", "MAFFT"), ("mode", "auto")],
        6,
    );
    assert_eq!(
        methods(&record),
        "Alignment was made with MAFFT v7.526 and HelixLoom 0.2.0 on 2026-11-13, on this computer, in 1.23 s. \
         Settings: engine: MAFFT; mode: auto. \
         It was made from 6 inputs, which are listed below with their sizes and their checksums."
    );
    // one program, no setting, one input and none
    let plain = made("Species tree", &["HelixLoom 0.2.0"], &[], 1);
    assert_eq!(
        methods(&plain),
        "Species tree was made with HelixLoom 0.2.0 on 2026-11-13, on this computer, in 1.23 s. \
         It has no settings. \
         It was made from 1 input, which is listed below with its size and its checksum."
    );
    assert!(methods(&made("X", &["a", "b", "c"], &[], 0))
        .starts_with("X was made with a, b and c on 2026-11-13, on this computer, in 1.23 s. It has no settings. It had no inputs."));
}

#[test]
fn the_programs_of_all_the_parts_are_listed_once_each_in_the_order_they_come() {
    let parts = [
        part(
            "One",
            made("A", &["MAFFT v7.526", "HelixLoom 0.2.0"], &[], 0),
        ),
        part(
            "Two",
            made("B", &["IQ-TREE 3.1.4", "HelixLoom 0.2.0"], &[], 0),
        ),
        part("Three", made("C", &["MAFFT v7.526"], &[], 0)),
    ];
    assert_eq!(
        programs(&parts),
        ["MAFFT v7.526", "HelixLoom 0.2.0", "IQ-TREE 3.1.4"]
    );
    assert!(programs(&[]).is_empty());
}

#[test]
fn a_page_has_a_part_for_each_result_with_its_methods_its_settings_its_inputs_and_the_programs_at_the_end(
) {
    let mut first = part(
        "The alignment",
        made(
            "Alignment",
            &["HelixLoom 0.2.0"],
            &[("engine", "built in")],
            2,
        ),
    );
    first.figure = Some((
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"10\" height=\"10\"></svg>\n".to_string(),
        "The first 2 rows.".to_string(),
    ));
    first.checked = crate::checks::CHECKS
        .iter()
        .find(|c| c.analysis == "Alignment (built in)");
    first.tables = vec![Table {
        caption: "Scores".to_string(),
        header: vec!["Row".to_string(), "Score".to_string()],
        rows: vec![vec!["a".to_string(), "1".to_string()]],
    }];
    let second = part(
        "The tree",
        made("Tree estimation", &["HelixLoom 0.2.0"], &[], 0),
    );
    let page = html("Results of a project", "2026-11-14", &[first, second]);

    assert!(page.starts_with("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<title>Results of a project</title>"));
    assert!(page.ends_with("</body>\n</html>\n"));
    assert!(page.contains(&format!(
        "<p class=\"dim\">Written by HelixLoom {} on 2026-11-14: 2 results.</p>",
        env!("CARGO_PKG_VERSION")
    )));
    // a list of the parts, and each part under its number
    assert!(page.contains("<li><a href=\"#result1\">The alignment</a></li>\n<li><a href=\"#result2\">The tree</a></li>"));
    assert!(page.contains("<section id=\"result1\">\n<h2>1. The alignment</h2>\n<h3>Methods</h3>\n<p>Alignment was made with HelixLoom 0.2.0 on 2026-11-13"));
    assert!(page.contains("<section id=\"result2\">\n<h2>2. The tree</h2>"));
    // how its method was checked, from the table of checks
    assert!(page.contains("<h3>How this method was checked</h3>\n<p><span class=\"grade\">second implementation</span> The plain aligner"));
    assert!(page.contains(
        "<p><strong>Not checked:</strong> This shows that the fast code is the plain method"
    ));
    // the figure in the page itself, the settings, the part's own table, the inputs with their checksums
    assert!(page.contains("<figure>\n<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"10\" height=\"10\"></svg>\n<figcaption>The first 2 rows.</figcaption>\n</figure>"));
    assert!(page.contains("<h3>Settings</h3>\n<table>\n<tr><th>Setting</th><th>Value</th></tr>\n<tr><td>engine</td><td>built in</td></tr>\n</table>"));
    assert!(page.contains("<h3>Scores</h3>\n<table>\n<tr><th>Row</th><th>Score</th></tr>\n<tr><td>a</td><td>1</td></tr>\n</table>"));
    assert!(
        page.contains("<tr><td>seq &lt;1&gt;</td><td>101</td><td class=\"mono\">abc1</td></tr>")
    );
    // a part with no settings, no inputs and nothing in the table of checks says so, and has no empty tables
    let tree = &page[page
        .find("<section id=\"result2\">")
        .expect("the second part")..];
    assert!(tree.contains("The table of checks has nothing for this result."));
    assert!(
        !tree.contains("<h3>Settings</h3>")
            && !tree.contains("<h3>Inputs</h3>")
            && !tree.contains("<figure>")
    );
    assert!(page.contains("<h2>Programs</h2>\n<ul>\n<li>HelixLoom 0.2.0</li>\n</ul>"));
}

#[test]
fn one_result_is_said_in_the_singular_without_a_list_and_what_a_person_typed_cannot_add_to_the_page(
) {
    let mut sly = part(
        "<script>alert(1)</script>",
        made("A & B", &["prog \"x\" <1>"], &[("k<", "v>")], 0),
    );
    sly.figure = Some((String::new(), "cap<tion>".to_string()));
    sly.notes = vec!["a <b>note</b>".to_string()];
    let page = html("T<i>tle", "2026-11-14", &[sly]);
    assert!(page.contains(": 1 result.</p>"));
    assert!(!page.contains("<ol>"));
    assert!(!page.contains("<script>") && !page.contains("<i>") && !page.contains("<tion>"));
    assert!(page.contains("<title>T&lt;i&gt;tle</title>"));
    assert!(page.contains("<h2>1. &lt;script&gt;alert(1)&lt;/script&gt;</h2>"));
    assert!(page.contains("A &amp; B was made with prog &quot;x&quot; &lt;1&gt;"));
    assert!(page.contains("<tr><td>k&lt;</td><td>v&gt;</td></tr>"));
    assert!(page.contains("<figcaption>cap&lt;tion&gt;</figcaption>"));
    assert!(page.contains("<p class=\"dim\">a &lt;b&gt;note&lt;/b&gt;</p>"));
    assert_eq!(escaped("a<b>&\"c\""), "a&lt;b&gt;&amp;&quot;c&quot;");
}
