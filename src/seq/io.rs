//! Reading and writing sequence files: GenBank (in the original's style, with its colour
//! qualifiers), FASTA and bare sequence text.

use super::doc::Doc;
use super::features::{color_from_name, type_colors, Feature};

/// Label qualifiers in the order the original prefers them.
const LABEL_KEYS: &[&str] = &[
    "label",
    "ApEinfo_label",
    "gene",
    "locus_tag",
    "product",
    "standard_name",
    "note",
];

pub fn parse_any(text: &str, default_name: &str) -> Result<Vec<Doc>, String> {
    let trimmed = text.trim_start_matches('\u{feff}').trim_start();
    let docs = if trimmed.starts_with("LOCUS") {
        parse_genbank(trimmed)
    } else if trimmed.starts_with('>') {
        parse_fasta(trimmed)
    } else {
        parse_bare(trimmed, default_name)
    };
    if docs.is_empty() {
        Err("This file doesn't look like a GenBank, FASTA or plain sequence file.".to_string())
    } else {
        Ok(docs)
    }
}

fn parse_bare(text: &str, name: &str) -> Vec<Doc> {
    let letters: Vec<u8> = text.bytes().filter(|c| c.is_ascii_alphabetic()).collect();
    let visible = text
        .bytes()
        .filter(|c| !c.is_ascii_whitespace() && !c.is_ascii_digit())
        .count();
    let nucleotides = letters
        .iter()
        .filter(|c| b"ACGTUNacgtun".contains(c))
        .count();
    // mostly letters, and mostly A/C/G/T/U/N: otherwise it is some other kind of text
    if letters.is_empty()
        || letters.len() * 100 < visible * 95
        || nucleotides * 100 < letters.len() * 85
    {
        return Vec::new();
    }
    vec![Doc::new(name, letters)]
}

pub fn parse_fasta(text: &str) -> Vec<Doc> {
    let mut docs: Vec<Doc> = Vec::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            let (name, rest) = header
                .trim()
                .split_once(char::is_whitespace)
                .unwrap_or((header.trim(), ""));
            let mut d = Doc::new(if name.is_empty() { "Untitled" } else { name }, Vec::new());
            d.definition = rest.trim().to_string();
            docs.push(d);
        } else if let Some(d) = docs.last_mut() {
            d.seq
                .extend(line.bytes().filter(|c| c.is_ascii_alphabetic()));
        }
    }
    docs.retain(|d| !d.seq.is_empty());
    docs
}

/// FASTA text the way the original's Save writes it (`save_file`, AppMain.tcl 0x8a413): a header
/// with the file's name without its extension, then the whole sequence on one line.
pub fn to_fasta(doc: &Doc, file_stem: &str) -> String {
    format!(">{file_stem}\n{}\n", String::from_utf8_lossy(&doc.seq))
}

// ---------------------------------------------------------------------------
// GenBank reading
// ---------------------------------------------------------------------------

pub fn parse_genbank(text: &str) -> Vec<Doc> {
    let mut records: Vec<Vec<&str>> = vec![Vec::new()];
    for line in text.lines() {
        if line.trim_end() == "//" {
            records.push(Vec::new());
        } else if let Some(r) = records.last_mut() {
            r.push(line);
        }
    }
    records.iter().filter_map(|r| parse_record(r)).collect()
}

enum Section {
    Header,
    Features,
    Origin,
}

fn parse_record(lines: &[&str]) -> Option<Doc> {
    let locus = lines.iter().find(|l| l.starts_with("LOCUS"))?;
    let name = locus.split_whitespace().nth(1).unwrap_or("Untitled");
    let mut doc = Doc::new(name, Vec::new());
    doc.circular = locus.to_ascii_lowercase().contains("circular");
    let mut section = Section::Header;
    let mut raw: Vec<RawFeature> = Vec::new();
    let mut in_definition = false;
    // while a COMMENT line goes on: is it one to keep (true) or one of ApE's own `ApEinfo:` lines?
    let mut keeping_comment: Option<bool> = None;
    for &line in lines {
        match section {
            Section::Header => {
                if line.starts_with("FEATURES") {
                    section = Section::Features;
                } else if line.starts_with("ORIGIN") {
                    section = Section::Origin;
                } else if let Some(rest) = line.strip_prefix("DEFINITION") {
                    doc.definition = rest.trim().to_string();
                    (in_definition, keeping_comment) = (true, None);
                } else if in_definition && line.starts_with("            ") {
                    doc.definition.push(' ');
                    doc.definition.push_str(line.trim());
                } else if let Some(rest) = line.strip_prefix("COMMENT") {
                    in_definition = false;
                    keeping_comment = comment_line(&mut doc, rest.trim(), None);
                } else if keeping_comment.is_some() && line.starts_with("            ") {
                    keeping_comment = comment_line(&mut doc, line.trim(), keeping_comment);
                } else {
                    (in_definition, keeping_comment) = (false, None);
                }
            }
            Section::Features => {
                if line.starts_with("ORIGIN") {
                    section = Section::Origin;
                } else {
                    feature_line(&mut raw, line);
                }
            }
            Section::Origin => doc
                .seq
                .extend(line.bytes().filter(|c| c.is_ascii_alphabetic())),
        }
    }
    doc.features = raw.into_iter().filter_map(|r| r.finish()).collect();
    let n = doc.seq.len();
    doc.features
        .retain(|f| f.spans.iter().all(|&(_, e)| e <= n));
    Some(doc)
}

/// One line of a COMMENT block. ApE's own `ApEinfo:` lines are read for what they say (the
/// methylation flag) and left out of the comment; any other line is added to it, a continuation
/// (`after` is `Some`) joined to the line before. Returns what the next continuation does.
fn comment_line(doc: &mut Doc, text: &str, after: Option<bool>) -> Option<bool> {
    let continued = after.is_some();
    let ours = if continued {
        after == Some(false)
    } else {
        text.starts_with("ApEinfo:")
    };
    if ours {
        if let Some(i) = text.find("ApEinfo:methylated:") {
            doc.methylated = text[i..].ends_with('1');
        }
        return Some(false);
    }
    if text.is_empty() {
        return Some(true);
    }
    if continued && !doc.comment.is_empty() {
        doc.comment.push(' ');
    } else if !doc.comment.is_empty() {
        doc.comment.push('\n');
    }
    doc.comment.push_str(text);
    Some(true)
}

#[derive(Default)]
struct RawFeature {
    key: String,
    location: String,
    quals: Vec<(String, String)>,
    open_quote: bool,
}

fn feature_line(raw: &mut Vec<RawFeature>, line: &str) {
    if line.trim().is_empty() || line.starts_with("FEATURES") {
        return;
    }
    let indent = line.len() - line.trim_start().len();
    if indent < 20 {
        let (key, loc) = line
            .trim_start()
            .split_once(char::is_whitespace)
            .unwrap_or((line.trim(), ""));
        raw.push(RawFeature {
            key: key.to_string(),
            location: loc.trim().to_string(),
            ..Default::default()
        });
        return;
    }
    let Some(cur) = raw.last_mut() else { return };
    let text = line.trim();
    if cur.open_quote {
        if let Some((_, v)) = cur.quals.last_mut() {
            v.push(' ');
            v.push_str(text.trim_end_matches('"'));
        }
        cur.open_quote = !text.ends_with('"');
    } else if let Some(q) = text.strip_prefix('/') {
        let (k, v) = q.split_once('=').unwrap_or((q, ""));
        let quoted = v.starts_with('"');
        let value = v.trim_matches('"').to_string();
        cur.open_quote = quoted && !(v.len() > 1 && v.ends_with('"'));
        cur.quals.push((k.to_string(), value));
    } else {
        cur.location.push_str(text);
    }
}

impl RawFeature {
    fn finish(self) -> Option<Feature> {
        let ranges = parse_location(&self.location)?;
        let reverse = ranges.first().is_some_and(|r| r.2);
        let mut spans: Vec<(usize, usize)> = ranges.iter().map(|&(s, e, _)| (s, e)).collect();
        spans.sort_unstable();
        let qual = |k: &str| {
            self.quals
                .iter()
                .find(|(q, _)| q == k)
                .map(|(_, v)| v.clone())
        };
        let name = LABEL_KEYS
            .iter()
            .find_map(|k| qual(k))
            .unwrap_or_else(|| self.key.clone());
        let (def_f, def_r) = type_colors(&self.key);
        let mut f = Feature::new(name, &self.key, 0, 1, reverse);
        f.spans = spans;
        f.fwd_color = qual("ApEinfo_fwdcolor")
            .and_then(|c| color_from_name(&c))
            .unwrap_or(def_f);
        f.rev_color = qual("ApEinfo_revcolor")
            .and_then(|c| color_from_name(&c))
            .unwrap_or(def_r);
        f.qualifiers = self
            .quals
            .into_iter()
            .filter(|(k, _)| {
                !k.starts_with("ApEinfo") && !matches!(k.as_str(), "label" | "locus_tag")
            })
            .collect();
        Some(f)
    }
}

/// Ranges as (start0, end_exclusive, on_minus_strand). `None` for remote or malformed locations.
fn parse_location(s: &str) -> Option<Vec<(usize, usize, bool)>> {
    let s = s.trim();
    if let Some(inner) = s
        .strip_prefix("complement(")
        .and_then(|r| r.strip_suffix(')'))
    {
        return Some(
            parse_location(inner)?
                .into_iter()
                .map(|(a, b, c)| (a, b, !c))
                .collect(),
        );
    }
    for op in ["join(", "order("] {
        if let Some(inner) = s.strip_prefix(op).and_then(|r| r.strip_suffix(')')) {
            let mut out = Vec::new();
            for part in split_top_level(inner) {
                out.extend(parse_location(part)?);
            }
            return Some(out);
        }
    }
    let clean: String = s.chars().filter(|c| !matches!(c, '<' | '>')).collect();
    let (a, b) = match clean.split_once("..") {
        Some((a, b)) => (a.parse::<usize>().ok()?, b.parse::<usize>().ok()?),
        None => {
            let a = clean.split('^').next()?.parse::<usize>().ok()?;
            (a, a)
        }
    };
    (a >= 1 && b >= a).then(|| vec![(a - 1, b, false)])
}

fn split_top_level(s: &str) -> Vec<&str> {
    let (mut depth, mut start, mut parts) = (0i32, 0, Vec::new());
    for (i, c) in s.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => depth -= 1,
            ',' if depth == 0 => {
                parts.push(&s[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    parts.push(&s[start..]);
    parts
}

// ---------------------------------------------------------------------------
// GenBank writing
// ---------------------------------------------------------------------------

pub fn to_genbank(doc: &Doc) -> String {
    let name: String = doc
        .name
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .take(16)
        .collect();
    let mut out = format!(
        "LOCUS       {:<16} {:>11} bp {:<3}{:<6}     {:<8} {:<3} {:<11}\n",
        name,
        doc.seq.len(),
        "",
        "DNA",
        if doc.circular { "circular" } else { "linear" },
        "",
        today()
    );
    let definition = if doc.definition.is_empty() {
        "."
    } else {
        doc.definition.as_str()
    };
    push_wrapped(&mut out, "DEFINITION  ", "            ", definition, 67);
    for line in doc.comment.lines() {
        push_wrapped(&mut out, "COMMENT     ", "            ", line, 67);
    }
    out.push_str(&format!(
        "COMMENT     ApEinfo:methylated:{}\n",
        u8::from(doc.methylated)
    ));
    if !doc.features.is_empty() {
        out.push_str("FEATURES             Location/Qualifiers\n");
    }
    for f in &doc.features {
        out.push_str(&format!("     {:<16}{}\n", f.kind, f.location_text()));
        let hex = |c: [u8; 3]| format!("#{:02x}{:02x}{:02x}", c[0], c[1], c[2]);
        let mut quals: Vec<(String, String)> = f.qualifiers.clone();
        quals.push(("label".into(), f.name.clone()));
        quals.push(("ApEinfo_fwdcolor".into(), hex(f.fwd_color)));
        quals.push(("ApEinfo_revcolor".into(), hex(f.rev_color)));
        for (k, v) in quals {
            push_wrapped(
                &mut out,
                "                     ",
                "                     ",
                &format!("/{k}=\"{v}\""),
                58,
            );
        }
    }
    out.push_str("ORIGIN\n");
    for (i, line) in doc.seq.chunks(60).enumerate() {
        let groups: Vec<String> = line
            .chunks(10)
            .map(|g| String::from_utf8_lossy(g).into_owned())
            .collect();
        out.push_str(&format!("{:>9} {}\n", i * 60 + 1, groups.join(" ")));
    }
    out.push_str("//\n");
    out
}

fn push_wrapped(out: &mut String, first: &str, rest: &str, text: &str, width: usize) {
    let mut prefix = first;
    let mut remaining = text;
    while remaining.len() > width {
        let split = remaining[..=width]
            .rfind(' ')
            .filter(|&i| i > 0)
            .unwrap_or(width);
        out.push_str(prefix);
        out.push_str(remaining[..split].trim_end());
        out.push('\n');
        remaining = remaining[split..].trim_start();
        prefix = rest;
    }
    out.push_str(prefix);
    out.push_str(remaining);
    out.push('\n');
}

/// Today's date as `05-NOV-2026` (UTC), without a date library.
fn today() -> String {
    let secs = web_time::SystemTime::now()
        .duration_since(web_time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    let z = (secs / 86_400) as i64 + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    const MONTHS: [&str; 12] = [
        "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
    ];
    format!("{day:02}-{}-{year}", MONTHS[(month - 1) as usize])
}

#[cfg(test)]
mod tests {
    use super::*;

    // A real multi-line literal: a trailing backslash would strip the indentation GenBank relies on.
    const NCBI: &str = r##"LOCUS       TEST1                     60 bp    DNA     circular SYN 01-JAN-2020
DEFINITION  A test plasmid.
FEATURES             Location/Qualifiers
     source          1..60
                     /organism="synthetic construct"
     CDS             complement(join(5..10,20..28))
                     /gene="abc"
                     /note="first line
                     second line"
     promoter        <35..>40
                     /ApEinfo_fwdcolor="#112233"
ORIGIN
        1 gaattcaaaa aaaaaaaaaa gggggggggg cccccccccc tttttttttt acgtacgtac
//
"##;

    #[test]
    fn reads_locus_features_and_sequence() {
        let docs = parse_any(NCBI, "x").unwrap();
        assert_eq!(docs.len(), 1);
        let d = &docs[0];
        assert_eq!(
            (d.name.as_str(), d.circular, d.seq.len()),
            ("TEST1", true, 60)
        );
        assert_eq!(d.definition, "A test plasmid.");
        assert_eq!(d.features.len(), 3);
        let cds = &d.features[1];
        assert_eq!((cds.name.as_str(), cds.reverse), ("abc", true));
        assert_eq!(cds.spans, vec![(4, 10), (19, 28)]);
        assert!(cds
            .qualifiers
            .iter()
            .any(|(k, v)| k == "note" && v == "first line second line"));
        let promoter = &d.features[2];
        assert_eq!(promoter.spans, vec![(34, 40)]);
        assert_eq!(promoter.fwd_color, [0x11, 0x22, 0x33]);
        assert_eq!(&d.seq[..6], b"gaattc");
    }

    #[test]
    fn genbank_round_trips_through_the_writer() {
        let mut d = parse_any(NCBI, "x").unwrap().remove(0);
        d.methylated = false;
        let text = to_genbank(&d);
        assert!(text.starts_with("LOCUS       TEST1"));
        assert!(text.contains("COMMENT     ApEinfo:methylated:0"));
        assert!(text.contains("     CDS             complement(join(5..10,20..28))"));
        let back = parse_any(&text, "x").unwrap().remove(0);
        assert_eq!(back.seq, d.seq);
        assert_eq!(back.circular, d.circular);
        assert!(!back.methylated);
        assert_eq!(back.features.len(), d.features.len());
        for (a, b) in back.features.iter().zip(&d.features) {
            assert_eq!(
                (&a.name, &a.kind, &a.spans, a.reverse),
                (&b.name, &b.kind, &b.spans, b.reverse)
            );
            assert_eq!((a.fwd_color, a.rev_color), (b.fwd_color, b.rev_color));
        }
        assert!(back.features[1]
            .qualifiers
            .iter()
            .any(|(k, v)| k == "note" && v.contains("second line")));
    }

    #[test]
    fn a_file_comment_is_kept_and_apes_own_lines_are_not_part_of_it() {
        let mut d = Doc::new("p", b"acgtacgtac".to_vec());
        d.methylated = false;
        let long = "Forward GGATCCGATTACAGGCATCCAGTCAGTCGATTACAGGCATCC Tm1: 52.0°C Tm2: 61.6°C.";
        d.comment = format!("PCR:\ntemplate: pDemo\n{long}\nProduct length: 10");
        let text = to_genbank(&d);
        assert!(text
            .contains("COMMENT     PCR:\nCOMMENT     template: pDemo\nCOMMENT     Forward GGATCC"));
        // the long line is wrapped like a definition, and the methylation line follows the comment
        assert!(text.contains("\n            "));
        assert!(text.find("Product length").unwrap() < text.find("ApEinfo:methylated:0").unwrap());
        let back = parse_any(&text, "x").unwrap().remove(0);
        assert_eq!(back.comment, d.comment);
        assert!(!back.methylated);
        // a file without a comment has none, whatever ApE's own lines say
        let plain = parse_any(NCBI, "x").unwrap().remove(0);
        assert_eq!(plain.comment, "");
        let mut flagged = Doc::new("q", b"acgt".to_vec());
        flagged.methylated = true;
        let again = parse_any(&to_genbank(&flagged), "x").unwrap().remove(0);
        assert_eq!((again.comment.as_str(), again.methylated), ("", true));
    }

    #[test]
    fn sequence_lines_are_numbered_in_blocks_of_sixty() {
        let mut d = Doc::new("seq", vec![b'a'; 130]);
        d.circular = false;
        let text = to_genbank(&d);
        assert!(text.contains(
            "\n        1 aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa aaaaaaaaaa\n"
        ));
        assert!(text.contains("\n      121 aaaaaaaaaa\n"));
        assert!(text.contains("linear"));
    }

    #[test]
    fn reads_multi_record_fasta_and_bare_sequences() {
        let docs = parse_any(">r1 first read\nACGT\nacgt\n>r2\nGGGG\n", "f").unwrap();
        assert_eq!(docs.len(), 2);
        assert_eq!(
            (docs[0].name.as_str(), docs[0].definition.as_str()),
            ("r1", "first read")
        );
        assert_eq!(docs[0].seq, b"ACGTacgt");
        let bare = parse_any("1 acgtacgt\n9 ggcc", "raw").unwrap();
        assert_eq!((bare[0].name.as_str(), bare[0].seq.len()), ("raw", 12));
        assert!(parse_any("this is not dna at all, sorry!", "bad").is_err());
        // saved like the original: file name as the header, the whole sequence on one line
        assert_eq!(to_fasta(&docs[0], "reads"), ">reads\nACGTacgt\n");
        let long = Doc::new("long", b"ACGT".repeat(100));
        assert_eq!(to_fasta(&long, "x").lines().count(), 2);
    }

    #[test]
    fn wrapping_keeps_long_values_on_one_logical_line() {
        let mut d = Doc::new("w", b"ACGTACGTAC".to_vec());
        let mut f = Feature::new("long", "misc_feature", 0, 5, false);
        f.qualifiers.push(("note".into(), "word ".repeat(40)));
        d.features.push(f);
        let back = parse_any(&to_genbank(&d), "w").unwrap().remove(0);
        let note = &back.features[0]
            .qualifiers
            .iter()
            .find(|(k, _)| k == "note")
            .unwrap()
            .1;
        assert_eq!(note.split_whitespace().count(), 40);
    }

    #[test]
    fn the_date_helper_formats_a_real_date() {
        let t = today();
        assert_eq!(t.len(), 11);
        assert!(t.as_bytes()[2] == b'-' && t.as_bytes()[6] == b'-');
    }
}
