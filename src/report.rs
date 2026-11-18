//! A report of results (`docs/ROADMAP.md` 1.7): one HTML file that stands by itself, with a part for each result. A part has the
//! methods in sentences written from the record of how the result was made (`provenance`), how the method was checked
//! (`checks`), the figure as it is drawn in the window (SVG, in the file itself), the settings, and the inputs with their sizes
//! and checksums; the programs are listed at the end with their versions. Nothing here draws or reads the window: `ui/report.rs`
//! gathers the parts.

use crate::checks::Check;
use crate::provenance::Provenance;
use std::fmt::Write;

/// A table of a part: a line above it, the names of its columns and its rows.
#[derive(Clone, Debug, PartialEq)]
pub struct Table {
    pub caption: String,
    pub header: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

/// One result in a report.
#[derive(Clone, Debug)]
pub struct Section {
    pub title: String,
    pub made: Provenance,
    /// How the method that made it was checked, when the table of checks has it.
    pub checked: Option<&'static Check>,
    /// What a reader must be told about this part (a figure that could not be drawn), a line each.
    pub notes: Vec<String>,
    /// The figure as SVG, and the line under it.
    pub figure: Option<(String, String)>,
    pub tables: Vec<Table>,
}

/// What must stand for a character in the text of a page.
pub fn escaped(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// `one, two and three`
pub fn listed(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [one] => one.clone(),
        [most @ .., last] => format!("{} and {last}", most.join(", ")),
    }
}

/// The methods of a result in sentences, from its record: what made it, when and where, with which settings, from how many
/// inputs.
pub fn methods(made: &Provenance) -> String {
    let settings = if made.parameters.is_empty() {
        "It has no settings.".to_string()
    } else {
        let each: Vec<String> = made
            .parameters
            .iter()
            .map(|(k, v)| format!("{k}: {v}"))
            .collect();
        format!("Settings: {}.", each.join("; "))
    };
    let inputs = match made.inputs.len() {
        0 => "It had no inputs.".to_string(),
        1 => "It was made from 1 input, which is listed below with its size and its checksum.".to_string(),
        n => format!("It was made from {n} inputs, which are listed below with their sizes and their checksums."),
    };
    format!(
        "{} was made with {} on {}, on {}, in {:.2} s. {settings} {inputs}",
        made.analysis,
        listed(&made.programs),
        made.date,
        made.ran_on,
        made.seconds
    )
}

/// The programs of every part, each once, in the order they first come.
pub fn programs(sections: &[Section]) -> Vec<String> {
    let mut all: Vec<String> = Vec::new();
    for program in sections.iter().flat_map(|s| &s.made.programs) {
        if !all.contains(program) {
            all.push(program.clone());
        }
    }
    all
}

const STYLE: &str = "\
body { font: 15px/1.5 system-ui, -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif; color: #1a1f2e; max-width: 62rem; \
margin: 2rem auto; padding: 0 1.25rem; }
h1 { font-size: 1.7rem; margin-bottom: 0.2rem; }
h2 { font-size: 1.3rem; margin-top: 2.4rem; border-top: 1px solid #d6dae3; padding-top: 1.2rem; }
h3 { font-size: 1rem; margin: 1.2rem 0 0.3rem; }
p { margin: 0.4rem 0; }
.dim { color: #5b6478; }
figure { margin: 1rem 0; overflow-x: auto; }
figure svg { box-sizing: border-box; max-width: 100%; height: auto; border: 1px solid #d6dae3; border-radius: 6px; }
figcaption { color: #5b6478; font-size: 0.9rem; }
table { border-collapse: collapse; margin: 0.4rem 0; font-size: 0.9rem; }
th, td { border: 1px solid #d6dae3; padding: 0.2rem 0.6rem; text-align: left; vertical-align: top; }
th { background: #f1f3f8; }
td.mono, code { font-family: 'JetBrains Mono', Menlo, Consolas, monospace; font-size: 0.85em; word-break: break-all; }
.grade { display: inline-block; padding: 0 0.5rem; border-radius: 0.7rem; background: #e9ecf5; font-size: 0.85rem; }
";

fn table(out: &mut String, t: &Table) {
    let _ = writeln!(out, "<h3>{}</h3>\n<table>", escaped(&t.caption));
    let cells = |tag: &str, row: &[String]| {
        row.iter()
            .map(|c| format!("<{tag}>{}</{tag}>", escaped(c)))
            .collect::<String>()
    };
    let _ = writeln!(out, "<tr>{}</tr>", cells("th", &t.header));
    for row in &t.rows {
        let _ = writeln!(out, "<tr>{}</tr>", cells("td", row));
    }
    out.push_str("</table>\n");
}

fn section(out: &mut String, number: usize, s: &Section) {
    let _ = writeln!(
        out,
        "<section id=\"result{number}\">\n<h2>{number}. {}</h2>",
        escaped(&s.title)
    );
    let _ = writeln!(
        out,
        "<h3>Methods</h3>\n<p>{}</p>",
        escaped(&methods(&s.made))
    );
    match s.checked {
        Some(check) => {
            let _ = writeln!(
                out,
                "<h3>How this method was checked</h3>\n<p><span class=\"grade\">{}</span> {} {}</p>\n\
                 <p><strong>Not checked:</strong> {}</p>",
                escaped(check.standard.label()),
                escaped(check.against),
                escaped(check.result),
                escaped(check.limits)
            );
        }
        None => out.push_str(
            "<h3>How this method was checked</h3>\n<p class=\"dim\">The table of checks has nothing for this result.</p>\n",
        ),
    }
    for note in &s.notes {
        let _ = writeln!(out, "<p class=\"dim\">{}</p>", escaped(note));
    }
    if let Some((svg, caption)) = &s.figure {
        // the figure is this program's own SVG, whose text is escaped where it is made
        let _ = writeln!(
            out,
            "<figure>\n{svg}<figcaption>{}</figcaption>\n</figure>",
            escaped(caption)
        );
    }
    if !s.made.parameters.is_empty() {
        table(
            out,
            &Table {
                caption: "Settings".into(),
                header: vec!["Setting".into(), "Value".into()],
                rows: s
                    .made
                    .parameters
                    .iter()
                    .map(|(k, v)| vec![k.clone(), v.clone()])
                    .collect(),
            },
        );
    }
    for t in &s.tables {
        table(out, t);
    }
    if !s.made.inputs.is_empty() {
        let _ = writeln!(
            out,
            "<h3>Inputs</h3>\n<table>\n<tr><th>Name</th><th>Bytes</th><th>Checksum</th></tr>"
        );
        for i in &s.made.inputs {
            let _ = writeln!(
                out,
                "<tr><td>{}</td><td>{}</td><td class=\"mono\">{}</td></tr>",
                escaped(&i.name),
                i.size,
                escaped(&i.digest)
            );
        }
        out.push_str("</table>\n");
    }
    out.push_str("</section>\n");
}

/// The report as a page: a heading, a part for each result in the order given, and the programs. `date` is the day the report
/// is written (YYYY-MM-DD).
pub fn html(title: &str, date: &str, sections: &[Section]) -> String {
    let mut out = format!(
        "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<title>{0}</title>\n<style>\n{STYLE}</style>\n\
         </head>\n<body>\n<h1>{0}</h1>\n<p class=\"dim\">Written by HelixLoom {1} on {2}: {3}.</p>\n",
        escaped(title),
        env!("CARGO_PKG_VERSION"),
        escaped(date),
        match sections.len() {
            1 => "1 result".to_string(),
            n => format!("{n} results"),
        }
    );
    if sections.len() > 1 {
        out.push_str("<ol>\n");
        for (k, s) in sections.iter().enumerate() {
            let _ = writeln!(
                out,
                "<li><a href=\"#result{}\">{}</a></li>",
                k + 1,
                escaped(&s.title)
            );
        }
        out.push_str("</ol>\n");
    }
    for (k, s) in sections.iter().enumerate() {
        section(&mut out, k + 1, s);
    }
    out.push_str("<h2>Programs</h2>\n<ul>\n");
    for program in programs(sections) {
        let _ = writeln!(out, "<li>{}</li>", escaped(&program));
    }
    out.push_str(
        "</ul>\n<p class=\"dim\">A checksum is of the bytes of an input: the same input gives the same checksum, so a result can be \
         traced to exactly what it was made from.</p>\n</body>\n</html>\n",
    );
    out
}

#[cfg(test)]
mod tests;
