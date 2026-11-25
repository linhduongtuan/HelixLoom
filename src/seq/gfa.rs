//! Inspect the segments and links of a GFA assembly graph. Keep one entry per segment name, including names first seen on a link.

#[cfg(not(target_arch = "wasm32"))]
pub mod cli;
mod record;

use record::Record;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::{self, BufRead, Write};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Summary {
    pub segments: u64,
    /// Number of L records, including repeated records.
    pub links: u64,
    /// Connected components when link orientation is ignored.
    pub components: u64,
    pub isolated: u64,
    /// Sum of declared or explicit segment lengths, counting each S record once.
    pub known_bases: u64,
    pub unknown_lengths: u64,
}

impl Summary {
    const METRICS: [&'static str; 6] = [
        "segments",
        "links",
        "components",
        "isolated_segments",
        "known_segment_bases",
        "segments_without_length",
    ];

    /// The same six rows are used by the CLI, pipeline cache and workstation.
    pub fn to_tsv(&self) -> String {
        let values = [
            self.segments,
            self.links,
            self.components,
            self.isolated,
            self.known_bases,
            self.unknown_lengths,
        ];
        let mut out = String::new();
        for (key, value) in Self::METRICS.into_iter().zip(values) {
            writeln!(out, "{key}\t{value}").expect("writing to a String");
        }
        out
    }

    pub fn from_tsv(text: &str) -> Result<Self, String> {
        let mut values = [0; Self::METRICS.len()];
        let mut seen = [false; Self::METRICS.len()];
        for (row, line) in text.lines().enumerate() {
            let (key, value) = line.split_once('\t').ok_or_else(|| {
                format!(
                    "GFA summary line {} needs a metric and count separated by a tab",
                    row + 1
                )
            })?;
            let index = Self::METRICS
                .iter()
                .position(|&name| name == key)
                .ok_or_else(|| format!("unknown GFA summary metric {key:?}"))?;
            if seen[index] {
                return Err(format!("duplicate GFA summary metric {key:?}"));
            }
            values[index] = value
                .parse::<u64>()
                .map_err(|_| format!("invalid count for GFA summary metric {key:?}"))?;
            seen[index] = true;
        }
        if let Some(index) = seen.iter().position(|present| !present) {
            return Err(format!(
                "missing GFA summary metric {:?}",
                Self::METRICS[index]
            ));
        }
        let [segments, links, components, isolated, known_bases, unknown_lengths] = values;
        if isolated > components || components > segments {
            return Err("inconsistent GFA summary: isolated segments must be at most components, and components at most segments".into());
        }
        if unknown_lengths > segments || (unknown_lengths == segments && known_bases != 0) {
            return Err("inconsistent GFA summary: unknown lengths and known bases disagree with the segment count".into());
        }
        // Every component containing linked segments needs a link (possibly a self-link).
        // Connecting n segments into c components also requires at least n - c links.
        let linked_components = components - isolated;
        if (linked_components == 0 && links != 0)
            || links < linked_components.max(segments - components)
        {
            return Err("inconsistent GFA summary: the link count cannot give these components and isolated segments".into());
        }
        Ok(Self {
            segments,
            links,
            components,
            isolated,
            known_bases,
            unknown_lengths,
        })
    }
}

struct Segment {
    parent: usize,
    first_seen: usize,
    declared_on: Option<usize>,
    length: Option<u64>,
    /// Largest CIGAR consumption at either end of a link, and the line on which it occurs.
    /// A link with an unspecified overlap still makes the segment incident to a link.
    overlap: Option<(u64, usize)>,
}

/// Indexed by the same IDs as the graph, allocated only for the per-segment report.
#[derive(Default)]
struct TableFields {
    depth: Option<f64>,
    kmer_count: Option<u64>,
    link_records: usize,
}

#[derive(Default)]
struct Graph {
    names: HashMap<String, usize>,
    segments: Vec<Segment>,
    summary: Summary,
}

impl Graph {
    fn id(&mut self, name: &str, line: usize) -> usize {
        if let Some(&id) = self.names.get(name) {
            return id;
        }
        let id = self.segments.len();
        self.names.insert(name.to_owned(), id);
        self.segments.push(Segment {
            parent: id,
            first_seen: line,
            declared_on: None,
            length: None,
            overlap: None,
        });
        id
    }

    fn endpoint(&mut self, name: &str, span: Option<u64>, line: usize) -> usize {
        let id = self.id(name, line);
        let segment = &mut self.segments[id];
        let consumed = span.unwrap_or(0);
        if segment
            .overlap
            .is_none_or(|(largest, _)| consumed > largest)
        {
            segment.overlap = Some((consumed, line));
        }
        id
    }

    /// Follow parent links, shortening them on the way to the root.
    fn root(&mut self, mut id: usize) -> usize {
        while self.segments[id].parent != id {
            let parent = self.segments[id].parent;
            let grandparent = self.segments[parent].parent;
            self.segments[id].parent = grandparent;
            id = grandparent;
        }
        id
    }

    fn validate_links(&self) -> Result<(), String> {
        // Break ties by the order names appeared, so hash iteration never changes the diagnostic.
        let failure = self
            .names
            .iter()
            .filter_map(|(name, &id)| {
                let segment = &self.segments[id];
                if segment.declared_on.is_none() {
                    return Some((segment.first_seen, id, name, None));
                }
                let length = segment.length?;
                let (span, line) = segment.overlap?;
                (span > length).then_some((line, id, name, Some((span, length))))
            })
            .min_by_key(|&(line, id, _, _)| (line, id));
        match failure {
            Some((line, _, name, None)) => Err(format!("line {line}: link refers to undefined segment '{name}'")),
            Some((line, _, name, Some((span, length)))) => Err(format!("line {line}: overlap consumes {span} bases of segment '{name}', which has length {length}")),
            None => Ok(()),
        }
    }

    fn write_segments(&self, fields: &[TableFields], mut out: impl Write) -> io::Result<()> {
        // First appearance can be on a link. Only the table needs the S records in input order.
        let mut rows: Vec<_> = self.names.iter().collect();
        rows.sort_unstable_by_key(|&(_, &id)| self.segments[id].declared_on);
        writeln!(out, "segment\tlength_bp\tdepth\tkmer_count\tlink_records")?;
        for (name, &id) in rows {
            let row = &fields[id];
            writeln!(
                out,
                "{name}\t{}\t{}\t{}\t{}",
                integer_field(self.segments[id].length),
                row.depth.map_or_else(|| ".".to_owned(), depth_text),
                integer_field(row.kmer_count),
                row.link_records
            )?;
        }
        Ok(())
    }
}

fn integer_field(value: Option<u64>) -> String {
    value.map_or_else(|| ".".to_owned(), |value| value.to_string())
}

/// Rust's shortest round-trip digits, with scientific notation for tiny or large values.
fn depth_text(depth: f64) -> String {
    if depth == 0.0 {
        return "0".to_owned();
    }
    if (0.0001..1e16).contains(&depth) {
        depth.to_string()
    } else {
        format!("{depth:e}")
    }
}

/// Links may precede their segments. Segment sequences are used only while reading their record.
fn read_graph(
    mut reader: impl BufRead,
    mut table: Option<&mut Vec<TableFields>>,
) -> Result<Graph, String> {
    let mut graph = Graph::default();
    let mut header_seen = false;
    let mut buffer = String::new();
    for line in 1.. {
        buffer.clear();
        let bytes = reader
            .read_line(&mut buffer)
            .map_err(|e| format!("line {line}: {e}"))?;
        if bytes == 0 {
            break;
        }
        let text = buffer.strip_suffix('\n').unwrap_or(&buffer);
        let text = text.strip_suffix('\r').unwrap_or(text);
        match record::parse(text).map_err(|e| format!("line {line}: {e}"))? {
            Record::Skip => (),
            Record::Header => header_seen = true,
            Record::Segment {
                name,
                length,
                depth,
                kmer_count,
            } => {
                let id = graph.id(name, line);
                let segment = &mut graph.segments[id];
                if let Some(first) = segment.declared_on {
                    return Err(format!(
                        "line {line}: duplicate segment '{name}' (first declared on line {first})"
                    ));
                }
                segment.declared_on = Some(line);
                segment.length = length;
                if let Some(fields) = table.as_deref_mut() {
                    fields.resize_with(graph.segments.len(), TableFields::default);
                    fields[id].depth = depth;
                    fields[id].kmer_count = kmer_count;
                }
                if let Some(length) = length {
                    graph.summary.known_bases = graph
                        .summary
                        .known_bases
                        .checked_add(length)
                        .ok_or_else(|| {
                            format!("line {line}: sum of segment lengths exceeds {}", u64::MAX)
                        })?;
                } else {
                    graph.summary.unknown_lengths += 1;
                }
            }
            Record::Link {
                from,
                to,
                from_span,
                to_span,
            } => {
                let a = graph.endpoint(from, from_span, line);
                let b = graph.endpoint(to, to_span, line);
                if let Some(fields) = table.as_deref_mut() {
                    fields.resize_with(graph.segments.len(), TableFields::default);
                    fields[a].link_records += 1;
                    if a != b {
                        fields[b].link_records += 1;
                    }
                }
                let a = graph.root(a);
                let b = graph.root(b);
                graph.segments[a].parent = b;
                graph.summary.links += 1;
            }
        }
    }
    if graph.segments.is_empty() && !header_seen {
        return Err("no GFA header or graph records found".into());
    }
    graph.validate_links()?;
    graph.summary.segments = graph.segments.len() as u64;
    for (id, segment) in graph.segments.iter().enumerate() {
        graph.summary.isolated += u64::from(segment.overlap.is_none());
        graph.summary.components += u64::from(id == segment.parent);
    }
    Ok(graph)
}

/// Inspect the whole graph and report its counts and known bases.
pub fn inspect(reader: impl BufRead) -> Result<Summary, String> {
    Ok(read_graph(reader, None)?.summary)
}

/// Write per-segment metadata in declaration order after validating the entire input.
pub fn segment_table(reader: impl BufRead, out: impl Write) -> Result<(), String> {
    let mut fields = Vec::new();
    read_graph(reader, Some(&mut fields))?
        .write_segments(&fields, out)
        .map_err(|e| format!("standard output: {e}"))
}

#[cfg(test)]
mod tests;
