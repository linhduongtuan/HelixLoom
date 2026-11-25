//! The header, segment and link records of GFA 1. Other record types are reported to the caller.

use std::collections::HashSet;

#[derive(Debug, PartialEq)]
pub(super) enum Record<'a> {
    Skip,
    Header,
    Segment {
        name: &'a str,
        length: Option<u64>,
        depth: Option<f64>,
        kmer_count: Option<u64>,
    },
    Link {
        from: &'a str,
        to: &'a str,
        from_span: Option<u64>,
        to_span: Option<u64>,
    },
}

fn name(field: &str) -> Result<&str, String> {
    if field.is_empty()
        || field.starts_with(['*', '='])
        || field.contains("+,")
        || field.contains("-,")
        || !field.bytes().all(|b| b.is_ascii_graphic())
    {
        Err(format!("invalid segment name {field:?}"))
    } else {
        Ok(field)
    }
}

fn orientation(field: &str) -> Result<(), String> {
    match field {
        "+" | "-" => Ok(()),
        _ => Err(format!(
            "expected + or - for the orientation, got {field:?}"
        )),
    }
}

/// GFA uses the From segment as the SAM reference and the To segment as the query.
fn overlap(mut cigar: &str) -> Result<(u64, u64), String> {
    let (mut from, mut to) = (0u64, 0u64);
    if cigar.is_empty() {
        return Err("missing link CIGAR".into());
    }
    while !cigar.is_empty() {
        let op_at = cigar
            .bytes()
            .position(|b| !b.is_ascii_digit())
            .ok_or("missing operation in link CIGAR")?;
        let count: u64 = cigar[..op_at]
            .parse()
            .map_err(|_| "invalid count in link CIGAR")?;
        let op = cigar.as_bytes()[op_at];
        if !b"MIDNSHPX=".contains(&op) {
            return Err(format!("invalid link CIGAR operation {:?}", char::from(op)));
        }
        if b"MDNX=".contains(&op) {
            from = from
                .checked_add(count)
                .ok_or("From span overflows in link CIGAR")?;
        }
        if b"MISX=".contains(&op) {
            to = to
                .checked_add(count)
                .ok_or("To span overflows in link CIGAR")?;
        }
        cigar = &cigar[op_at + 1..];
    }
    Ok((from, to))
}

#[derive(Default)]
struct SegmentTags {
    length: Option<u64>,
    depth: Option<f64>,
    kmer_count: Option<u64>,
}

/// Check the envelope and uniqueness of every optional field. Interpret VN on headers and LN, KC
/// and DP on segments. Other tag values are opaque (including JSON, arrays and external URIs).
fn tags<'a>(record: &str, fields: impl Iterator<Item = &'a str>) -> Result<SegmentTags, String> {
    let mut seen = HashSet::new();
    let mut segment = SegmentTags::default();
    for field in fields {
        let mut parts = field.splitn(3, ':');
        let tag = parts.next().unwrap_or_default();
        let kind = parts.next().unwrap_or_default();
        let value = parts.next().unwrap_or_default();
        let key = tag.as_bytes();
        if key.len() != 2
            || !key[0].is_ascii_alphabetic()
            || !key[1].is_ascii_alphanumeric()
            || !matches!(kind, "A" | "i" | "f" | "Z" | "J" | "H" | "B")
            || value.is_empty()
            || !value.bytes().all(|b| (b' '..=b'~').contains(&b))
        {
            return Err("invalid optional field (expected TAG:TYPE:VALUE)".into());
        }
        if !seen.insert(tag) {
            return Err(format!("duplicate {tag} tag"));
        }
        match (record, tag) {
            ("H", "VN") if kind != "Z" || !matches!(value, "1.0" | "1.1" | "1.2") => {
                return Err(format!(
                    "expected VN:Z with GFA version 1.0, 1.1 or 1.2, got {field}"
                ));
            }
            ("S", "LN" | "KC") => {
                let count = value
                    .parse::<u64>()
                    .ok()
                    .filter(|_| kind == "i")
                    .ok_or_else(|| format!("{tag} must be a nonnegative integer (type i)"))?;
                if tag == "LN" {
                    segment.length = Some(count);
                } else {
                    segment.kmer_count = Some(count);
                }
            }
            ("S", "DP") => {
                const INVALID_DEPTH: &str = "DP must be a finite nonnegative depth (type f)";
                let mantissa = value.split(['e', 'E']).next().unwrap_or_default();
                if kind != "f"
                    || !mantissa.ends_with(|c: char| c.is_ascii_digit())
                    || (mantissa.starts_with('-')
                        && mantissa.bytes().any(|b| matches!(b, b'1'..=b'9')))
                    || !value.parse::<f32>().is_ok_and(f32::is_finite)
                {
                    return Err(INVALID_DEPTH.into());
                }
                // GFA's f tag has single-precision range. Read the original decimal for the report.
                segment.depth = Some(value.parse().map_err(|_| INVALID_DEPTH)?);
            }
            _ => {}
        }
    }
    Ok(segment)
}

pub(super) fn parse(line: &str) -> Result<Record<'_>, String> {
    if !line.is_ascii() {
        return Err("GFA records must be ASCII".into());
    }
    if line.is_empty() || line.starts_with('#') {
        return Ok(Record::Skip);
    }
    let mut fields = line.split('\t');
    let record = fields.next().unwrap_or_default();
    match record {
        "H" => {
            tags(record, fields)?;
            Ok(Record::Header)
        }
        "S" => {
            let id = name(fields.next().ok_or("missing segment name")?)?;
            let sequence = fields
                .next()
                .ok_or("missing segment sequence (use * if unknown)")?;
            let metadata = tags(record, fields)?;
            let length = if sequence == "*" {
                metadata.length
            } else {
                if sequence.is_empty()
                    || !sequence
                        .bytes()
                        .all(|b| b.is_ascii_alphabetic() || b == b'=' || b == b'.')
                {
                    return Err("invalid segment sequence".into());
                }
                let actual = sequence.len() as u64;
                if metadata.length.is_some_and(|len| len != actual) {
                    return Err(format!(
                        "LN disagrees with the sequence length of {id} ({actual})"
                    ));
                }
                Some(actual)
            };
            Ok(Record::Segment {
                name: id,
                length,
                depth: metadata.depth,
                kmer_count: metadata.kmer_count,
            })
        }
        "L" => {
            let from = name(fields.next().ok_or("missing From segment")?)?;
            orientation(fields.next().ok_or("missing From orientation")?)?;
            let to = name(fields.next().ok_or("missing To segment")?)?;
            orientation(fields.next().ok_or("missing To orientation")?)?;
            let cigar = fields.next().ok_or("missing overlap (use * if unknown)")?;
            let (from_span, to_span) = if cigar == "*" {
                (None, None)
            } else {
                let (reference, query) = overlap(cigar)?;
                (Some(reference), Some(query))
            };
            tags(record, fields)?;
            Ok(Record::Link {
                from,
                to,
                from_span,
                to_span,
            })
        }
        _ => Err(format!(
            "unsupported GFA record {record:?} (expected H, S or L)"
        )),
    }
}

#[cfg(test)]
mod tests;
