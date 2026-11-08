//! List Enzymes and the Quick Lists (Enzymes menu): the enzymes of a sequence as a text, in four formats. The rules are the
//! original's: `list_enzymes_text` (AppMain.tcl L16463 @0xb8d47) makes the text, `list_enzymes_window` (L16378 @0xb82e9)
//! writes the lines above it, `enz_site_usage_filter` (L15931 @0xb4584) and `enz_count` (L15950 @0xb4790) pick the absent and
//! the unique enzymes, and `list_enzymes_dialog` (L16305 @0xb7316) holds the choices and what they start as (L66960-66962).

use std::collections::BTreeSet;

/// How many columns the first two formats start with, and the most the dialog offers (`info(enz_list_columns)` is 6; the
/// menu holds 1 to 10).
pub const DEFAULT_COLUMNS: usize = 6;
pub const MAX_COLUMNS: usize = 10;

/// The dialog's four formats (`modelist`, L16315), numbered 0 to 3 in the original.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    /// Names alone, in columns, sorted by name.
    Names,
    /// Names and how many sites, in columns, sorted by name.
    Counts,
    /// One enzyme to a line: name, site, how many sites, where; the fewest sites first.
    Positions,
    /// One enzyme to a line: site and name, sorted by site.
    BySite,
}

impl Format {
    pub const ALL: [Format; 4] = [
        Format::Names,
        Format::Counts,
        Format::Positions,
        Format::BySite,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Format::Names => "Name",
            Format::Counts => "Name, #Sites",
            Format::Positions => "Name, Site, #Sites, Positions",
            Format::BySite => "Site, Name",
        }
    }

    /// The dialog turns the Columns choice off for the two formats that write a line to an enzyme (L16320).
    pub fn uses_columns(self) -> bool {
        matches!(self, Format::Names | Format::Counts)
    }
}

/// Which enzymes are listed (`info(enz_list_enz)`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Which {
    All,
    Selected,
    Absent,
    Unique,
}

impl Which {
    pub const ALL: [Which; 4] = [Which::All, Which::Selected, Which::Absent, Which::Unique];

    /// What the dialog calls it (L16334-16337).
    pub fn label(self) -> &'static str {
        match self {
            Which::All => "All",
            Which::Selected => "Selected Enzymes Only",
            Which::Absent => "Absent Enzymes",
            Which::Unique => "Unique Enzymes",
        }
    }

    /// The name of the list in the window's heading (L16384-16397).
    pub fn title(self) -> &'static str {
        match self {
            Which::All => "All Enzymes",
            Which::Selected => "Selected Sites",
            Which::Absent => "Absent Sites",
            Which::Unique => "Unique Sites",
        }
    }
}

/// The bases looked at when only the selection is: `start` bases come before it and `end` bases up to its last, which is how
/// the heading gives it ("From 10 to 20" for bases 11 to 20).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

/// One enzyme: its name, its pattern as the enzyme file writes it (`G^AATTC`), and where its sites start (numbered from 1, in
/// the order found: two for a site that matches on both strands).
#[derive(Clone, Debug)]
pub struct Entry {
    pub name: String,
    pub pattern: String,
    pub sites: Vec<usize>,
}

impl Entry {
    /// The sites `enz_count` counts when picking absent and unique enzymes: the ones on the selected bases, `start + 1` to `end`.
    pub fn count_in(&self, span: Option<Span>) -> usize {
        let Some(s) = span else {
            return self.sites.len();
        };
        self.sites
            .iter()
            .filter(|&&p| p > s.start && p <= s.end)
            .count()
    }

    /// The sites the text lists, which `list_enzymes_text` takes from `start` (not `start + 1`) to `end`: a site that begins
    /// on the base just before the selection is listed too.
    fn listed(&self, span: Option<Span>) -> Vec<usize> {
        let Some(s) = span else {
            return self.sites.clone();
        };
        let within = |&&p: &&usize| p >= s.start && p <= s.end;
        self.sites.iter().filter(within).copied().collect()
    }
}

/// The enzymes a list is about (`list_enzymes_text`, L16488-16501): all of them in the order of the enzyme file; the selected
/// ones in the order of the Enzyme Selector's list (sorted); or those with no site, or exactly one, on the bases looked at.
pub fn choose<'a>(
    entries: &'a [Entry],
    which: Which,
    selected: &BTreeSet<String>,
    span: Option<Span>,
) -> Vec<&'a Entry> {
    match which {
        Which::All => entries.iter().collect(),
        Which::Selected => {
            let named = |n: &String| entries.iter().find(|e| &e.name == n);
            selected.iter().filter_map(named).collect()
        }
        Which::Absent => entries.iter().filter(|e| e.count_in(span) == 0).collect(),
        Which::Unique => entries.iter().filter(|e| e.count_in(span) == 1).collect(),
    }
}

/// Cells in rows of `columns`; a row ends with a newline only when it is full, so the last row of a list that does not fill it
/// has none (the newline belongs to the last column's format, L16533 and L16549).
fn in_columns(cells: &[String], columns: usize) -> String {
    let columns = columns.clamp(1, MAX_COLUMNS);
    let mut out = String::new();
    for row in cells.chunks(columns) {
        for (j, cell) in row.iter().enumerate() {
            out.push_str(cell);
            if j + 1 == columns {
                out.push('\n');
            }
        }
    }
    out
}

/// The text of a list (`list_enzymes_text`, L16502-16575). Every field is padded on the right and none is cut.
pub fn text(chosen: &[&Entry], format: Format, columns: usize, span: Option<Span>) -> String {
    let mut rows: Vec<(&Entry, Vec<usize>)> = chosen.iter().map(|&e| (e, e.listed(span))).collect();
    match format {
        Format::Names | Format::Counts => {
            rows.sort_by(|a, b| a.0.name.cmp(&b.0.name));
            let cell = |(e, sites): &(&Entry, Vec<usize>)| match format {
                Format::Names => format!("{:<15}", e.name),
                _ => format!("{:<12}{:<8}", e.name, sites.len()),
            };
            let cells: Vec<String> = rows.iter().map(cell).collect();
            in_columns(&cells, columns)
        }
        Format::Positions => {
            rows.sort_by_key(|(_, sites)| sites.len());
            let line = |(e, sites): &(&Entry, Vec<usize>)| {
                let at: Vec<String> = sites.iter().map(usize::to_string).collect();
                format!(
                    "{:<15}{:<15}{:<7}{:<7}\n",
                    e.name,
                    e.pattern,
                    sites.len(),
                    at.join(" ")
                )
            };
            rows.iter().map(line).collect()
        }
        Format::BySite => {
            rows.sort_by(|a, b| a.0.pattern.cmp(&b.0.pattern));
            let line = |(e, _): &(&Entry, Vec<usize>)| format!("{:<15}{:<15}\n", e.pattern, e.name);
            rows.iter().map(line).collect()
        }
    }
}

/// What the original's result window holds (`list_enzymes_window`, L16419-16443): the date, the sequence's name, the bases
/// looked at (`len` of them when no selection is), the name of the list, then the text. The original leaves a space after the
/// name and after the bases.
pub fn report(
    date: &str,
    name: &str,
    which: Which,
    span: Option<Span>,
    len: usize,
    body: &str,
) -> String {
    let (start, end) = span.map_or((0, len), |s| (s.start, s.end));
    format!(
        "{date}\n{name} \nFrom {start} to {end}. \n\n{}\n\n{body}",
        which.title()
    )
}

#[cfg(test)]
mod tests;
