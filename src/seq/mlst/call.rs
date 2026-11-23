//! The calls for a scheme, from what the search found at each locus, and what follows from them. The notation and the scoring are those
//! of mlst's README: `n` an exact allele; `~n` a novel allele, all of it present and alike enough, the closest being n; `n?` a partial
//! match to n; `-` missing; `n,m` more than one allele. A scheme's score, out of 100, is 90/N for each exact locus, 63/N for each novel,
//! 18/N for each partial (N the number of loci), and 10 more when the alleles are a known sequence type.
//!
//! What the README does not say, and this program chooses: the thresholds (`min_identity` 95, `min_coverage` 10, as percents); that
//! a locus with several places is called by the best kind among them (exact, then novel, then partial) and reports all the places of that
//! kind; that a locus with more than one allele scores as exact; the order in which schemes are tried when two score alike (the order given).

use super::scheme::Scheme;
use super::search::{find, Contigs, Found};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Params {
    /// How alike (percent) the stretch found must be to an allele for the call to be novel or partial: mlst's `--minid`.
    pub min_identity: f64,
    /// How much of an allele (percent) a partial match must cover: `--mincov`.
    pub min_coverage: f64,
    /// The score (1 to 100) a scheme needs to be chosen when schemes are tried: `--minscore`.
    pub min_score: f64,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            min_identity: 95.0,
            min_coverage: 10.0,
            min_score: 50.0,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Call {
    /// One or more alleles that are there exactly, by id.
    Exact(Vec<String>),
    /// Alleles that are not known, with the closest known allele and the bases found, for each place.
    Novel(Vec<(String, Vec<u8>)>),
    /// Part of an allele.
    Partial(String),
    Missing,
}

impl Call {
    /// The call as it is written: `n`, `n,m`, `~n`, `n?`, `-`.
    pub fn text(&self) -> String {
        match self {
            Call::Exact(ids) => ids.join(","),
            Call::Novel(found) => found
                .iter()
                .map(|(id, _)| format!("~{id}"))
                .collect::<Vec<_>>()
                .join(","),
            Call::Partial(id) => format!("{id}?"),
            Call::Missing => "-".to_string(),
        }
    }

    fn points(&self) -> f64 {
        match self {
            Call::Exact(_) => 90.0,
            Call::Novel(_) => 63.0,
            Call::Partial(_) => 18.0,
            Call::Missing => 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    /// A known sequence type.
    Perfect,
    /// Exact alleles that are in no profile.
    Novel,
    /// Not one allele matches.
    None,
    Mixed,
    Missing,
    /// None of those, and a score under 70.
    Bad,
    Ok,
}

impl Status {
    pub fn word(self) -> &'static str {
        match self {
            Status::Perfect => "PERFECT",
            Status::Novel => "NOVEL",
            Status::None => "NONE",
            Status::Mixed => "MIXED",
            Status::Missing => "MISSING",
            Status::Bad => "BAD",
            Status::Ok => "OK",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Typing {
    pub scheme: String,
    /// The call for each locus, in the order of the scheme.
    pub calls: Vec<(String, Call)>,
    pub sequence_type: Option<String>,
    pub score: f64,
    pub status: Status,
}

fn id_order(a: &str, b: &str) -> std::cmp::Ordering {
    match (a.parse::<u64>(), b.parse::<u64>()) {
        (Ok(x), Ok(y)) => x.cmp(&y),
        _ => a.cmp(b),
    }
}

/// The call for one locus from what was found at its places.
fn call_locus(found: &[Found], p: &Params) -> Call {
    let alike = |f: &&Found| f.identity() * 100.0 >= p.min_identity;
    let mut exact: Vec<String> = found.iter().filter_map(|f| f.exact.clone()).collect();
    if !exact.is_empty() {
        exact.sort_by(|a, b| id_order(a, b));
        exact.dedup();
        return Call::Exact(exact);
    }
    let novel: Vec<(String, Vec<u8>)> = found
        .iter()
        .filter(alike)
        .filter(|f| f.covered == f.allele_len)
        .map(|f| (f.allele.clone(), f.sequence.clone()))
        .collect();
    if !novel.is_empty() {
        return Call::Novel(novel);
    }
    found
        .iter()
        .filter(alike)
        .filter(|f| f.coverage() * 100.0 >= p.min_coverage)
        .max_by(|a, b| a.matches.cmp(&b.matches))
        .map_or(Call::Missing, |f| Call::Partial(f.allele.clone()))
}

/// The scheme's typing of the genome.
pub fn type_contigs(scheme: &Scheme, contigs: &Contigs, p: &Params) -> Typing {
    let calls: Vec<(String, Call)> = scheme
        .loci
        .iter()
        .map(|l| (l.name.clone(), call_locus(&find(contigs, l), p)))
        .collect();
    summarize(scheme, calls)
}

/// The sequence type, the score and the status that the calls come to.
fn summarize(scheme: &Scheme, calls: Vec<(String, Call)>) -> Typing {
    let single: Option<Vec<&str>> = calls
        .iter()
        .map(|(_, c)| match c {
            Call::Exact(ids) if ids.len() == 1 => Some(ids[0].as_str()),
            _ => None,
        })
        .collect();
    let sequence_type = single
        .and_then(|ids| scheme.sequence_type(&ids))
        .map(str::to_string);
    let n = calls.len() as f64;
    // the points are added before they are divided, so that a score that comes to a whole number is one (70, or 100)
    let score = calls.iter().map(|(_, c)| c.points()).sum::<f64>() / n
        + if sequence_type.is_some() { 10.0 } else { 0.0 };
    let all_exact = calls
        .iter()
        .all(|(_, c)| matches!(c, Call::Exact(ids) if ids.len() == 1));
    let status = if sequence_type.is_some() {
        Status::Perfect
    } else if all_exact {
        Status::Novel
    } else if calls.iter().all(|(_, c)| *c == Call::Missing) {
        Status::None
    } else if calls
        .iter()
        .any(|(_, c)| matches!(c, Call::Exact(ids) if ids.len() > 1))
    {
        Status::Mixed
    } else if calls.iter().any(|(_, c)| *c == Call::Missing) {
        Status::Missing
    } else if score < 70.0 {
        Status::Bad
    } else {
        Status::Ok
    };
    Typing {
        scheme: scheme.name.clone(),
        calls,
        sequence_type,
        score,
        status,
    }
}

#[cfg(test)]
mod tests;
