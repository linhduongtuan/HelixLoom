//! Connectors: fetch a record from a public database by its accession, keep it in a local cache, and remember where it came
//! from and under what licence. NCBI (nucleotide and protein), UniProt, the PDB and the AlphaFold DB.
//!
//! The accession is checked against the form each database gives its own before it goes into an address, so no text that a
//! person typed (or pasted) can add anything to a request. The web is reached by `curl`, as the structure prediction does, so
//! no network library is linked in; the tests use a fetcher that answers from a table and touch no network.

pub mod curl;

use crate::pipeline::flow::digest;
use std::path::{Path, PathBuf};
use web_time::{SystemTime, UNIX_EPOCH};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    /// NCBI's nucleotide records, as GenBank files.
    NcbiNucleotide,
    /// NCBI's protein records, as FASTA.
    NcbiProtein,
    /// UniProtKB, as FASTA.
    UniProt,
    /// The Protein Data Bank (RCSB), as a PDB file.
    Pdb,
    /// The AlphaFold DB's predicted structure of a UniProt protein, as a PDB file.
    AlphaFold,
}

/// What a fetched file is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    GenBank,
    Fasta,
    Pdb,
}

impl Source {
    pub const ALL: [Source; 5] = [
        Source::NcbiNucleotide,
        Source::NcbiProtein,
        Source::UniProt,
        Source::Pdb,
        Source::AlphaFold,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Source::NcbiNucleotide => "NCBI nucleotide",
            Source::NcbiProtein => "NCBI protein",
            Source::UniProt => "UniProt",
            Source::Pdb => "PDB",
            Source::AlphaFold => "AlphaFold DB",
        }
    }

    /// A short name safe in file names.
    fn key(self) -> &'static str {
        match self {
            Source::NcbiNucleotide => "ncbi-nucleotide",
            Source::NcbiProtein => "ncbi-protein",
            Source::UniProt => "uniprot",
            Source::Pdb => "pdb",
            Source::AlphaFold => "alphafold",
        }
    }

    pub fn format(self) -> Format {
        match self {
            Source::NcbiNucleotide => Format::GenBank,
            Source::NcbiProtein | Source::UniProt => Format::Fasta,
            Source::Pdb | Source::AlphaFold => Format::Pdb,
        }
    }

    /// What the database says of the use of its data, in a line, with where it says so. Written from the databases' own policy
    /// pages as they stood when this was written: check them before relying on it for a publication.
    pub fn licence(self) -> &'static str {
        match self {
            Source::NcbiNucleotide | Source::NcbiProtein => {
                "NCBI data are public, but a submitter may claim rights in a record; see https://www.ncbi.nlm.nih.gov/home/about/policies/"
            }
            Source::UniProt => "CC BY 4.0; see https://www.uniprot.org/help/license",
            Source::Pdb => "CC0 1.0 (no rights reserved); see https://www.rcsb.org/pages/policies",
            Source::AlphaFold => "CC BY 4.0; see https://alphafold.ebi.ac.uk/faq",
        }
    }

    /// What an accession of this database looks like, for the box where it is typed.
    pub fn example(self) -> &'static str {
        match self {
            Source::NcbiNucleotide => "NC_000913.3 or U00096.3",
            Source::NcbiProtein => "NP_414542.1",
            Source::UniProt | Source::AlphaFold => "P69905",
            Source::Pdb => "1A3N",
        }
    }
}

/// A record asked for: the database and an accession that has passed `check`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Request {
    pub source: Source,
    pub id: String,
}

impl Request {
    /// Check `text` against the form of the accessions of `source`; the request, with the accession in the database's own
    /// letters (upper case, no surrounding space).
    pub fn check(source: Source, text: &str) -> Result<Request, String> {
        let id = text.trim().to_ascii_uppercase();
        let ok = match source {
            Source::NcbiNucleotide | Source::NcbiProtein => ncbi_accession(&id),
            // UniProt names an isoform by a dash and a number after the accession: P69905-2
            Source::UniProt => {
                let (base, isoform) = id
                    .split_once('-')
                    .map_or((id.as_str(), None), |(b, n)| (b, Some(n)));
                uniprot_accession(base)
                    && isoform.is_none_or(|n| {
                        (1..=3).contains(&n.len()) && n.bytes().all(|b| b.is_ascii_digit())
                    })
            }
            Source::AlphaFold => uniprot_accession(&id),
            Source::Pdb => {
                id.len() == 4
                    && id.as_bytes()[0].is_ascii_digit()
                    && id.bytes().all(|b| b.is_ascii_alphanumeric())
            }
        };
        if ok {
            Ok(Request { source, id })
        } else {
            Err(format!(
                "\"{}\" is not an accession of {} (they look like {}).",
                text.trim().chars().take(40).collect::<String>(),
                source.label(),
                source.example()
            ))
        }
    }

    /// The address of the record (for the AlphaFold DB, of the entry that says where the file is).
    pub fn url(&self) -> String {
        let id = &self.id;
        match self.source {
            Source::NcbiNucleotide => format!(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={id}&rettype=gb&retmode=text&tool=HelixLoom"
            ),
            Source::NcbiProtein => format!(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={id}&rettype=fasta&retmode=text&tool=HelixLoom"
            ),
            Source::UniProt => format!("https://rest.uniprot.org/uniprotkb/{id}.fasta"),
            Source::Pdb => format!("https://files.rcsb.org/download/{id}.pdb"),
            Source::AlphaFold => format!("https://alphafold.ebi.ac.uk/api/prediction/{id}"),
        }
    }
}

/// NCBI accessions: letters (and an underscore for RefSeq) then digits, and a version: NC_000913.3, U00096.3, NP_414542.1,
/// AAA12345.1; and RefSeq's prokaryote genomes and contigs, which have NZ_ in front of a GenBank accession: NZ_CP009072.1.
fn ncbi_accession(id: &str) -> bool {
    let id = id.strip_prefix("NZ_").unwrap_or(id);
    let (body, version) = id.split_once('.').map_or((id, None), |(b, v)| (b, Some(v)));
    let letters = body.bytes().take_while(u8::is_ascii_alphabetic).count();
    let rest = &body[letters..];
    let rest = rest.strip_prefix('_').unwrap_or(rest);
    (1..=6).contains(&letters)
        && rest.len() >= 5
        && rest.len() <= 12
        && rest.bytes().all(|b| b.is_ascii_digit())
        && version
            .is_none_or(|v| !v.is_empty() && v.len() <= 3 && v.bytes().all(|b| b.is_ascii_digit()))
}

/// UniProt's own pattern: O, P or Q then a digit and three letters or digits and a digit (P69905), or the longer form that
/// starts with another letter (A0A024R161).
fn uniprot_accession(id: &str) -> bool {
    let b = id.as_bytes();
    let alnum = |c: u8| c.is_ascii_uppercase() || c.is_ascii_digit();
    match b.len() {
        6 => {
            let first_ok = matches!(b[0], b'O' | b'P' | b'Q')
                && b[1].is_ascii_digit()
                && b[2..5].iter().all(|&c| alnum(c))
                && b[5].is_ascii_digit();
            let other = b[0].is_ascii_uppercase()
                && !matches!(b[0], b'O' | b'P' | b'Q')
                && b[1].is_ascii_digit()
                && b[2].is_ascii_uppercase()
                && alnum(b[3])
                && alnum(b[4])
                && b[5].is_ascii_digit();
            first_ok || other
        }
        10 => {
            b[0].is_ascii_uppercase()
                && b[1].is_ascii_digit()
                && b[2].is_ascii_uppercase()
                && alnum(b[3])
                && alnum(b[4])
                && b[5].is_ascii_digit()
                && b[6].is_ascii_uppercase()
                && alnum(b[7])
                && alnum(b[8])
                && b[9].is_ascii_digit()
        }
        _ => false,
    }
}

// --- where it came from -------------------------------------------------------------------------------------------------------

/// Where a file came from, written beside it and carried into the document that is opened from it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Provenance {
    pub source: String,
    pub id: String,
    /// The address the file itself was fetched from.
    pub url: String,
    /// The day it was fetched, as YYYY-MM-DD (UTC).
    pub fetched: String,
    pub licence: String,
    /// Of the file, so that a copy that has been changed on the disk is known.
    pub digest: String,
}

impl Provenance {
    pub fn to_text(&self) -> String {
        format!(
            "helixloom-provenance 1\nsource: {}\nid: {}\nurl: {}\nfetched: {}\nlicence: {}\ndigest: {}\n",
            self.source, self.id, self.url, self.fetched, self.licence, self.digest
        )
    }

    pub fn from_text(text: &str) -> Option<Provenance> {
        let mut lines = text.lines();
        if lines.next()? != "helixloom-provenance 1" {
            return None;
        }
        let mut get = |key: &str| {
            lines
                .next()
                .and_then(|l| l.strip_prefix(key))
                .map(str::to_string)
        };
        Some(Provenance {
            source: get("source: ")?,
            id: get("id: ")?,
            url: get("url: ")?,
            fetched: get("fetched: ")?,
            licence: get("licence: ")?,
            digest: get("digest: ")?,
        })
    }

    /// The line that goes into the comment of a document made from the file (and so into the GenBank file saved from it).
    pub fn comment(&self) -> String {
        format!(
            "Fetched from {} {} on {} ({}). Licence: {}",
            self.source, self.id, self.fetched, self.url, self.licence
        )
    }
}

/// The date of `unix_seconds` as YYYY-MM-DD in UTC (the days-to-date arithmetic of the proleptic Gregorian calendar).
pub fn utc_date(unix_seconds: u64) -> String {
    let days = (unix_seconds / 86_400) as i64 + 719_468;
    let era = days.div_euclid(146_097);
    let day_of_era = days.rem_euclid(146_097);
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let m = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * m + 2) / 5 + 1;
    let month = if m < 10 { m + 3 } else { m - 9 };
    let year = year_of_era + era * 400 + i64::from(month <= 2);
    format!("{year:04}-{month:02}-{day:02}")
}

fn today() -> String {
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    utc_date(seconds)
}

// --- fetching -----------------------------------------------------------------------------------------------------------------

/// What an address answered.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Reply {
    pub status: u16,
    pub body: String,
}

/// Whatever gets an address: `curl` in the program, a table in the tests.
pub trait Fetcher {
    fn get(&self, url: &str) -> Result<Reply, String>;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fetched {
    pub request: Request,
    pub body: String,
    pub provenance: Provenance,
    /// The file was already in the local cache: nothing was asked of the database.
    pub from_cache: bool,
    /// Why a good record was not kept in the cache (the folder cannot be written): the record is still there to use.
    pub warning: Option<String>,
}

/// A folder of fetched files, each with its provenance beside it.
pub struct Store {
    dir: PathBuf,
}

impl Store {
    pub fn open(dir: &Path) -> Result<Store, String> {
        std::fs::create_dir_all(dir)
            .map_err(|e| format!("Could not make the folder {}: {e}", dir.display()))?;
        Ok(Store {
            dir: dir.to_path_buf(),
        })
    }

    fn files(&self, r: &Request) -> (PathBuf, PathBuf) {
        // the accession has passed `Request::check`: letters, digits, a dot and an underscore, so it is safe in a name
        let stem = format!("{}-{}", r.source.key(), r.id);
        (
            self.dir.join(format!("{stem}.body")),
            self.dir.join(format!("{stem}.provenance")),
        )
    }

    /// The file kept for `r`, if it is whole and unchanged.
    pub fn get(&self, r: &Request) -> Option<Fetched> {
        let (body_path, provenance_path) = self.files(r);
        let body = std::fs::read_to_string(body_path).ok()?;
        let provenance = Provenance::from_text(&std::fs::read_to_string(provenance_path).ok()?)?;
        (provenance.id == r.id && provenance.digest == digest(body.as_bytes())).then(|| Fetched {
            request: r.clone(),
            body,
            provenance,
            from_cache: true,
            warning: None,
        })
    }

    pub fn put(&self, f: &Fetched) -> Result<(), String> {
        let (body_path, provenance_path) = self.files(&f.request);
        std::fs::write(&body_path, &f.body)
            .and_then(|()| std::fs::write(&provenance_path, f.provenance.to_text()))
            .map_err(|e| format!("Could not keep the file in {}: {e}", self.dir.display()))
    }
}

impl Fetched {
    /// The sequences of a GenBank or FASTA file as documents, each with where it came from in its comment (which a GenBank
    /// file saved from it keeps).
    pub fn documents(&self) -> Result<Vec<crate::seq::doc::Doc>, String> {
        if self.request.source.format() == Format::Pdb {
            return Err("A structure file holds no sequence document.".into());
        }
        let mut docs = crate::seq::io::parse_any(&self.body, &self.request.id)?;
        let line = self.provenance.comment();
        for doc in &mut docs {
            if !doc.comment.is_empty() {
                doc.comment.push('\n');
            }
            doc.comment.push_str(&line);
        }
        Ok(docs)
    }
}

/// What the database's answer must hold for it to be the record and not a page that says there is none.
fn check_body(format: Format, body: &str) -> bool {
    match format {
        Format::GenBank => body.trim_start().starts_with("LOCUS"),
        Format::Fasta => body.trim_start().starts_with('>'),
        Format::Pdb => body
            .lines()
            .any(|l| l.starts_with("ATOM") || l.starts_with("HETATM")),
    }
}

/// Why an answer is not the record, in words for the person.
fn not_the_record(source: Source, id: &str, reply: &Reply) -> String {
    match reply.status {
        200 => {
            let quote: String = reply.body.trim().chars().take(80).collect();
            format!(
                "{} answered, but not with {}'s record: {quote}",
                source.label(),
                id
            )
        }
        404 if source == Source::Pdb => format!(
            "The PDB has no file in PDB format for {id}: there is no such entry, or it is a large structure that only comes in the mmCIF format."
        ),
        400 | 404 => format!("{} has no record called {id}.", source.label()),
        403 => format!(
            "{} refused the request (it said 403): some databases block a network that has sent them too many requests.",
            source.label()
        ),
        429 => format!(
            "{} says there have been too many requests; wait a moment and try again.",
            source.label()
        ),
        s if s >= 500 => format!(
            "{} is not answering just now (it said {s}); try again later.",
            source.label()
        ),
        s => format!("{} answered {s} for {id}.", source.label()),
    }
}

/// The address of the file of an AlphaFold DB entry, from the entry's JSON: it must be a file of the AlphaFold DB's own server.
fn alphafold_file(json: &str) -> Option<String> {
    let at = json.find("\"pdbUrl\"")? + "\"pdbUrl\"".len();
    let rest = json[at..]
        .trim_start()
        .strip_prefix(':')?
        .trim_start()
        .strip_prefix('"')?;
    let url = &rest[..rest.find('"')?];
    let safe = url
        .bytes()
        .all(|b| b.is_ascii_alphanumeric() || b"-_./:".contains(&b));
    (safe
        && url.starts_with("https://alphafold.ebi.ac.uk/files/")
        && url.ends_with(".pdb")
        && !url.contains(".."))
    .then(|| url.to_string())
}

/// The model versions of the AlphaFold DB's files, the newest first: tried when its entry service does not answer. Version 6 was
/// the newest when this was written; 4 and 5 answer 404 now.
const ALPHAFOLD_VERSIONS: [&str; 3] = ["v6", "v5", "v4"];

/// The record `request` asks for: from the store when it is there (and `refresh` is not asked), else from the database, then kept.
pub fn fetch(
    request: &Request,
    fetcher: &dyn Fetcher,
    store: Option<&Store>,
    refresh: bool,
) -> Result<Fetched, String> {
    if !refresh {
        if let Some(found) = store.and_then(|s| s.get(request)) {
            return Ok(found);
        }
    }
    let source = request.source;
    let ask = |url: &str| -> Result<Reply, String> {
        fetcher
            .get(url)
            .map_err(|why| format!("Could not reach {}: {why}", source.label()))
    };
    let mut url = request.url();
    let mut reply = ask(&url)?;
    if source == Source::AlphaFold {
        match reply.status {
            200 => {
                // the entry says where the file is
                url = alphafold_file(&reply.body).ok_or(format!(
                    "The AlphaFold DB has an entry for {} but names no structure file for it.",
                    request.id
                ))?;
                reply = ask(&url)?;
            }
            // there is no such entry: no file to look for
            400 | 404 => {}
            // the entry service did not answer (it refuses some networks), but the files are named by a pattern: the newest
            // model versions are tried, and the address that gave the file is the one that is kept
            _ => {
                for version in ALPHAFOLD_VERSIONS {
                    let file = format!(
                        "https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_{version}.pdb",
                        request.id
                    );
                    let tried = ask(&file)?;
                    if tried.status == 200 {
                        (url, reply) = (file, tried);
                        break;
                    }
                }
            }
        }
    }
    if reply.status != 200 || !check_body(source.format(), &reply.body) {
        return Err(not_the_record(source, &request.id, &reply));
    }
    let mut fetched = Fetched {
        request: request.clone(),
        provenance: Provenance {
            source: source.label().to_string(),
            id: request.id.clone(),
            url,
            fetched: today(),
            licence: source.licence().to_string(),
            digest: digest(reply.body.as_bytes()),
        },
        body: reply.body,
        from_cache: false,
        warning: None,
    };
    // a record that cannot be kept is still the record
    if let Some(Err(why)) = store.map(|s| s.put(&fetched)) {
        fetched.warning = Some(why);
    }
    Ok(fetched)
}

/// For the tests of the windows: a fetcher that answers from a table, from any thread.
#[cfg(test)]
pub(crate) struct Table {
    pub replies: Vec<(String, Reply)>,
    pub asked: std::sync::Mutex<Vec<String>>,
}

#[cfg(test)]
impl Fetcher for Table {
    fn get(&self, url: &str) -> Result<Reply, String> {
        self.asked
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .push(url.to_string());
        Ok(self
            .replies
            .iter()
            .find(|(u, _)| u == url)
            .map(|(_, r)| r.clone())
            .unwrap_or(Reply {
                status: 404,
                body: String::new(),
            }))
    }
}

#[cfg(test)]
mod tests;
