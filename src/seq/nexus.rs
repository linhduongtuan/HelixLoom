//! NEXUS tree files: the `trees` block with its `translate` table and the trees after it, as BEAST's TreeAnnotator, MrBayes, PAUP* and
//! FigTree write them. Only trees are read; every other block is passed over. The annotations inside the trees (`[&posterior=0.99]`) are
//! kept by `NwTree`'s reader.
//!
//! Written from the NEXUS format as Maddison, Swofford and Maddison describe it (*Systematic Biology* 1997) and from files that these
//! programs write.

use super::newick::NwTree;

/// Whether the text is a NEXUS file: its first word is `#NEXUS`, in any case.
pub fn is_nexus(text: &str) -> bool {
    let t = text.trim_start_matches('\u{feff}').trim_start();
    t.get(..6).is_some_and(|w| w.eq_ignore_ascii_case("#nexus"))
}

/// The text without the comments in square brackets that begin it.
fn after_comments(mut text: &str) -> &str {
    text = text.trim_start();
    while let Some(rest) = text.strip_prefix('[') {
        text = rest
            .find(']')
            .map_or("", |end| rest[end + 1..].trim_start());
    }
    text
}

/// The statements of a block, cut at the semicolons that are neither inside a quoted name nor inside a comment in square brackets.
fn statements(block: &str) -> Vec<&str> {
    let (mut out, mut from, mut quoted, mut depth) = (Vec::new(), 0, false, 0usize);
    for (i, c) in block.char_indices() {
        match c {
            '\'' if depth == 0 => quoted = !quoted,
            '[' if !quoted => depth += 1,
            ']' if !quoted => depth = depth.saturating_sub(1),
            ';' if !quoted && depth == 0 => {
                out.push(after_comments(&block[from..i]));
                from = i + 1;
            }
            _ => {}
        }
    }
    let rest = after_comments(&block[from..]).trim();
    if !rest.is_empty() {
        out.push(rest);
    }
    out
}

/// The words of a statement, cut at spaces that are not inside a quoted name; a quoted name loses its quotes (a doubled quote inside is one).
fn words(text: &str) -> Vec<String> {
    let (mut out, mut word, mut quoted, mut chars) =
        (Vec::new(), String::new(), false, text.chars().peekable());
    while let Some(c) = chars.next() {
        match c {
            '\'' if quoted && chars.peek() == Some(&'\'') => {
                word.push('\'');
                chars.next();
            }
            '\'' => quoted = !quoted,
            c if c.is_whitespace() && !quoted => {
                if !word.is_empty() {
                    out.push(std::mem::take(&mut word));
                }
            }
            c => word.push(c),
        }
    }
    if !word.is_empty() {
        out.push(word);
    }
    out
}

/// The block named `name` (the text between `begin name;` and its `end;`), found without regard to case.
fn block<'a>(text: &'a str, name: &str) -> Option<&'a str> {
    let lower = text.to_ascii_lowercase();
    let start = lower.find(&format!("begin {name};"))? + format!("begin {name};").len();
    let rest = &lower[start..];
    let end = ["\nend;", " end;", "\tend;", "\nendblock;"]
        .iter()
        .filter_map(|e| rest.find(e))
        .min()
        .unwrap_or(rest.len());
    Some(&text[start..start + end])
}

/// The trees of a NEXUS file, each with its name, their tip names put through the `translate` table when there is one.
pub fn parse_trees(text: &str) -> Result<Vec<(String, NwTree)>, String> {
    if !is_nexus(text) {
        return Err("This is not a NEXUS file: its first word should be #NEXUS.".into());
    }
    let body = block(text, "trees")
        .ok_or("This NEXUS file has no trees block (begin trees; ... end;).")?;
    let mut table: Vec<(String, String)> = Vec::new();
    let mut found = Vec::new();
    for statement in statements(body) {
        let lead = statement
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_ascii_lowercase();
        match lead.as_str() {
            "translate" => {
                for pair in statement["translate".len()..].split(',') {
                    if let [key, name, ..] = words(pair).as_slice() {
                        table.push((key.clone(), name.clone()));
                    }
                }
            }
            "tree" | "utree" => {
                let rest = statement[lead.len()..].trim_start();
                let (name, tree) = rest.split_once('=').ok_or_else(|| {
                    format!(
                        "A tree statement needs an equals sign: {}",
                        rest.chars().take(40).collect::<String>()
                    )
                })?;
                let mut parsed = NwTree::parse_all(&format!("{};", tree.trim()))?.remove(0);
                for node in parsed.nodes.iter_mut().filter(|n| n.children.is_empty()) {
                    if let Some((_, to)) = table.iter().find(|(k, _)| *k == node.name) {
                        node.name.clone_from(to);
                    }
                }
                found.push((
                    name.trim()
                        .trim_start_matches('*')
                        .trim()
                        .trim_matches('\'')
                        .to_string(),
                    parsed,
                ));
            }
            _ => {}
        }
    }
    if found.is_empty() {
        return Err("The trees block holds no tree statement (tree NAME = (...);).".into());
    }
    Ok(found)
}

/// The trees of a NEXUS or a Newick text.
pub fn parse_any(text: &str) -> Result<Vec<NwTree>, String> {
    if is_nexus(text) {
        parse_trees(text).map(|trees| trees.into_iter().map(|(_, t)| t).collect())
    } else {
        NwTree::parse_all(text)
    }
}

/// A NEXUS file of one tree, with the annotations when `annotated`. A name with anything but letters, digits and `_.-` is put in quotes.
pub fn write(tree: &NwTree, name: &str, annotated: bool) -> String {
    let body = if annotated {
        tree.newick_annotated()
    } else {
        tree.newick()
    };
    format!(
        "#NEXUS\n\nbegin trees;\n\ttree {} = [&R] {body}\nend;\n",
        quote(name)
    )
}

fn quote(name: &str) -> String {
    if !name.is_empty()
        && name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || "_.-".contains(c))
    {
        name.to_string()
    } else {
        format!("'{}'", name.replace('\'', "''"))
    }
}

#[cfg(test)]
mod tests;
