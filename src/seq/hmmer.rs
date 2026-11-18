//! Protein domains as HMMER 3 reports them: `hmmscan` searches protein sequences against a database of profile HMMs (Pfam) and, with
//! `--domtblout`, writes one line for each domain it finds (`docs/DOMAINS.md`). This reads that table. The search itself is HMMER's,
//! run as a program (`pipeline::domains`); nothing of its code is here.
//!
//! The table has 22 whitespace-separated fields and then the description of the family, as the HMMER 3.4 User's Guide gives them:
//! target name and accession (for hmmscan the *profile*), its length, query name and accession, its length, the E-value and score of
//! the whole comparison with the bias, the number of this domain and of all, the conditional and independent E-values, the
//! domain's score and bias, its start and end on the profile, on the sequence (the alignment) and on the sequence again (the
//! envelope), and the mean posterior probability of the aligned residues. Lines that start with `#` are comments.

use super::io::parse_fasta;

/// The most sequences and residues a search takes, and the longest sequence.
pub const MAX_SEQUENCES: usize = 10_000;
pub const MAX_LENGTH: usize = 100_000;
pub const MAX_RESIDUES: usize = 20_000_000;

/// The names and residues of the protein sequences in `text`: FASTA (a name is the first word of its header), or when there is no
/// `>` one sequence for each line. Letters only, in capitals (spaces, digits and `*` are left out).
pub fn read_proteins(text: &str) -> Result<Vec<(String, Vec<u8>)>, String> {
    let clean = |s: &[u8]| -> Vec<u8> {
        s.iter()
            .filter(|b| b.is_ascii_alphabetic())
            .map(u8::to_ascii_uppercase)
            .collect()
    };
    let found: Vec<(String, Vec<u8>)> = if text.contains('>') {
        parse_fasta(text)
            .into_iter()
            .map(|d| (d.name, clean(&d.seq)))
            .collect()
    } else {
        text.lines()
            .map(|l| clean(l.as_bytes()))
            .filter(|s| !s.is_empty())
            .enumerate()
            .map(|(k, s)| (format!("sequence {}", k + 1), s))
            .collect()
    };
    if found.is_empty() {
        return Err("Give at least one protein sequence (FASTA, or one on each line).".into());
    }
    if found.len() > MAX_SEQUENCES {
        return Err(format!(
            "{} sequences are more than the {MAX_SEQUENCES} that can be searched at once.",
            found.len()
        ));
    }
    if let Some((name, s)) = found.iter().find(|(_, s)| s.len() > MAX_LENGTH) {
        return Err(format!(
            "{name} has {} residues; the longest that is searched is {MAX_LENGTH}.",
            s.len()
        ));
    }
    let residues: usize = found.iter().map(|(_, s)| s.len()).sum();
    if residues > MAX_RESIDUES {
        return Err(format!(
            "{residues} residues are more than the {MAX_RESIDUES} that can be searched at once."
        ));
    }
    Ok(found)
}

/// One domain found in a sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct DomainHit {
    /// Which of the query sequences it is in (the order they were given in).
    pub sequence: usize,
    /// The family's name, as Pfam gives it (`Pkinase`), and its accession (`PF00069.30`; empty when there is none).
    pub family: String,
    pub accession: String,
    /// The length of the family's profile, in residues.
    pub family_length: usize,
    /// The bit score of the domain, and the expected number of such domains in the whole database by chance: the independent
    /// E-value (a stringent measure) and the conditional one (a permissive one).
    pub score: f64,
    pub i_evalue: f64,
    pub c_evalue: f64,
    /// Where the domain lies on the profile and on the sequence, counted from 1 and both ends included: the alignment, and the
    /// envelope that may reach beyond it.
    pub hmm_from: usize,
    pub hmm_to: usize,
    pub ali_from: usize,
    pub ali_to: usize,
    pub env_from: usize,
    pub env_to: usize,
    /// The mean posterior probability of the aligned residues, from 0 to 1.
    pub accuracy: f64,
    /// The family's description line.
    pub description: String,
}

impl DomainHit {
    /// Whether the envelopes of the two domains share a residue.
    pub fn overlaps(&self, other: &DomainHit) -> bool {
        self.sequence == other.sequence
            && self.env_from <= other.env_to
            && other.env_from <= self.env_to
    }
}

/// The 22 fields before the description.
const FIELDS: usize = 22;

/// A number in a field, or why the line is no good.
fn field<T: std::str::FromStr>(
    fields: &[&str],
    k: usize,
    what: &str,
    line: usize,
) -> Result<T, String> {
    fields[k - 1].parse().map_err(|_| {
        format!(
            "Line {line} of the domain table: {what} (field {k}) is {:?}, not a number.",
            fields[k - 1]
        )
    })
}

/// The domains in the text of a `--domtblout` table. `queries` are the names the sequences were given to the program under and
/// their lengths, in order; a line about a name that is not among them, or about a sequence of another length, is an error, and so
/// is a coordinate that is not inside its sequence or its profile, so that nothing drawn from the table can run off an end.
pub fn parse_domtblout(text: &str, queries: &[(&str, usize)]) -> Result<Vec<DomainHit>, String> {
    let mut hits = Vec::new();
    for (k, row) in text.lines().enumerate() {
        let line = k + 1;
        let row = row.trim();
        if row.is_empty() || row.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = row.split_whitespace().collect();
        if f.len() < FIELDS {
            return Err(format!(
                "Line {line} of the domain table has {} fields; it needs at least {FIELDS}.",
                f.len()
            ));
        }
        let family_length: usize = field(&f, 3, "the length of the profile", line)?;
        let query_length: usize = field(&f, 6, "the length of the sequence", line)?;
        let Some(sequence) = queries.iter().position(|(name, _)| *name == f[3]) else {
            return Err(format!(
                "Line {line} of the domain table is about {}, a sequence that was not searched.",
                f[3]
            ));
        };
        if queries[sequence].1 != query_length {
            return Err(format!(
                "Line {line} of the domain table gives {} a length of {query_length}, not {}.",
                f[3], queries[sequence].1
            ));
        }
        let hit = DomainHit {
            sequence,
            family: f[0].to_string(),
            accession: if f[1] == "-" {
                String::new()
            } else {
                f[1].to_string()
            },
            family_length,
            c_evalue: field(&f, 12, "the conditional E-value", line)?,
            i_evalue: field(&f, 13, "the independent E-value", line)?,
            score: field(&f, 14, "the domain score", line)?,
            hmm_from: field(&f, 16, "the start on the profile", line)?,
            hmm_to: field(&f, 17, "the end on the profile", line)?,
            ali_from: field(&f, 18, "the start of the alignment", line)?,
            ali_to: field(&f, 19, "the end of the alignment", line)?,
            env_from: field(&f, 20, "the start of the envelope", line)?,
            env_to: field(&f, 21, "the end of the envelope", line)?,
            accuracy: field(&f, 22, "the posterior probability", line)?,
            description: f[FIELDS..].join(" "),
        };
        check(&hit, query_length, line)?;
        hits.push(hit);
    }
    Ok(hits)
}

/// Whether the numbers of a hit make sense for a sequence of `query_length` residues.
pub(crate) fn check(h: &DomainHit, query_length: usize, line: usize) -> Result<(), String> {
    let inside = |from: usize, to: usize, length: usize| 1 <= from && from <= to && to <= length;
    let numbers_ok = [h.score, h.i_evalue, h.c_evalue, h.accuracy]
        .iter()
        .all(|x| x.is_finite())
        && h.i_evalue >= 0.0
        && h.c_evalue >= 0.0
        && (0.0..=1.0).contains(&h.accuracy);
    if !numbers_ok {
        return Err(format!(
            "Line {line} of the domain table has a score, an E-value or a probability that is out of range."
        ));
    }
    if !inside(h.hmm_from, h.hmm_to, h.family_length) {
        return Err(format!(
            "Line {line} of the domain table puts {} at {}-{} on a profile of {} residues.",
            h.family, h.hmm_from, h.hmm_to, h.family_length
        ));
    }
    if !inside(h.ali_from, h.ali_to, query_length) || !inside(h.env_from, h.env_to, query_length) {
        return Err(format!(
            "Line {line} of the domain table puts {} at {}-{} (envelope {}-{}) on a sequence of {query_length} residues.",
            h.family, h.ali_from, h.ali_to, h.env_from, h.env_to
        ));
    }
    Ok(())
}

/// The domains of one sequence in the order they lie along it (by the start of the alignment, then the end).
pub fn along(hits: &[DomainHit], sequence: usize) -> Vec<&DomainHit> {
    let mut mine: Vec<&DomainHit> = hits.iter().filter(|h| h.sequence == sequence).collect();
    mine.sort_by_key(|h| (h.ali_from, h.ali_to));
    mine
}

/// The domains in a text for people and spreadsheets: a line each with a header, tab-separated.
pub fn to_table(hits: &[DomainHit], names: &[String]) -> String {
    let mut text = String::from(
        "sequence\tfamily\taccession\tfrom\tto\tscore\ti_evalue\tc_evalue\tprofile_from\tprofile_to\tprofile_length\tdescription\n",
    );
    for h in hits {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{:e}\t{:e}\t{}\t{}\t{}\t{}\n",
            names[h.sequence].replace(['\t', '\n'], " "),
            h.family,
            h.accession,
            h.ali_from,
            h.ali_to,
            h.score,
            h.i_evalue,
            h.c_evalue,
            h.hmm_from,
            h.hmm_to,
            h.family_length,
            h.description.replace(['\t', '\n'], " ")
        ));
    }
    text
}

#[cfg(test)]
mod tests;
